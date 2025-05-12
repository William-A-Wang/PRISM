
import numpy as np
import pandas as pd
from numba import njit
from joblib import Parallel, delayed
from design_primers import PrimerSetBadnessFast



@njit
def generate_initial_solution_numba(badness_components, block_indices):
    """
    Construct an initial primer selection by choosing, for each block,
    the primer with the highest badness score on both the forward and
    reverse strands.
    """
    num_blocks = np.max(block_indices) + 1 
    selected_left = np.empty(num_blocks, dtype=np.int32)
    selected_right = np.empty(num_blocks, dtype=np.int32)
    selected_count = 0

    for block in range(num_blocks):
        block_mask = (block_indices == block)
        if np.any(block_mask):
            block_indices_local = np.where(block_mask)[0] 

            # Select forward primer with maximal badness
            left_badness = badness_components[0, block_indices_local]  
            best_left_idx = np.argmax(left_badness) 
            selected_left[selected_count] = block_indices_local[best_left_idx]  

            # Select reverse primer with maximal badness
            right_badness = badness_components[1, block_indices_local]  
            best_right_idx = np.argmax(right_badness)  
            selected_right[selected_count] = block_indices_local[best_right_idx] 

            selected_count += 1 

    return selected_left[:selected_count], selected_right[:selected_count]  



def calculate_objective_value_numba(indices, badness_components, primer_sequences):
    """
    Evaluate the objective function：f1 - f2
    """
    if len(indices) == 0:
        return 0.0

    # Sum of badness across selected primers
    f1_sum = np.sum(badness_components[0, indices]) + np.sum(badness_components[1, indices])

    # Retrieve sequences for thermodynamic evaluation
    selected_left_primers = np.array([primer_sequences[0][i] for i in indices]) 
    selected_right_primers = np.array([primer_sequences[1][i] for i in indices]) 

    # total_badness
    total_badness, _ = PrimerSetBadnessFast(selected_left_primers, selected_right_primers)

    # f2
    f2 = total_badness / len(indices) if len(indices) > 0 else 0

    return f1_sum - f2 


def process_block_candidates(i, remove_idx, block_candidates, current_left, current_right,
                             is_left, badness_components, primer_sequences, current_value):
    """
    Explore all candidate replacements within a single block and
    identify the swap that maximally improves the objective.
    """
    best_result = {
        'delta': -np.inf,
        'remove_idx': -1,
        'add_idx': -1,
        'is_left': is_left,
        'new_value': np.inf
    }

    for candidate_idx in block_candidates:
        if (is_left and candidate_idx not in current_left) or (not is_left and candidate_idx not in current_right):
            if is_left:
                temp_left = current_left.copy()
                temp_left[i] = candidate_idx
                combined = np.concatenate((temp_left, current_right))
            else:
                temp_right = current_right.copy()
                temp_right[i] = candidate_idx
                combined = np.concatenate((current_left, temp_right))

            new_value = calculate_objective_value_numba(
                combined,
                badness_components,
                primer_sequences
            )

            current_delta = current_value - new_value
            if current_delta > best_result['delta']:
                best_result.update({
                    'delta': current_delta,
                    'remove_idx': i,
                    'add_idx': candidate_idx,
                    'is_left': is_left,
                    'new_value': new_value
                })

    return best_result


def fast_local_search_numba(selected_left, selected_right, badness_components, block_indices,
                             primer_sequences, epsilon, n_jobs=-1):
    """
    Perform a local search with parallel candidate evaluation
    to iteratively improve the primer selection.
    """
    current_left = selected_left.copy()
    current_right = selected_right.copy()
    k = len(current_left)

    if k == 0:
        return current_left, current_right

    delta = np.inf

    current_value = calculate_objective_value_numba(
        np.concatenate((current_left, current_right)),
        badness_components,
        primer_sequences
    )

    while delta >= epsilon / k:
        parallel_tasks = []

        
        for i, remove_idx in enumerate(current_left):
            block_id = block_indices[remove_idx]
            block_candidates = np.where(block_indices == block_id)[0]
            parallel_tasks.append(delayed(process_block_candidates)(
                i, remove_idx, block_candidates,
                current_left, current_right,
                True, badness_components, primer_sequences,
                current_value
            ))

        
        for i, remove_idx in enumerate(current_right):
            block_id = block_indices[remove_idx]
            block_candidates = np.where(block_indices == block_id)[0]
            parallel_tasks.append(delayed(process_block_candidates)(
                i, remove_idx, block_candidates,
                current_left, current_right,
                False, badness_components, primer_sequences,
                current_value
            ))

        
        results = Parallel(n_jobs=n_jobs)(parallel_tasks)

        best_result = max(results, key=lambda x: x['delta'])

        if best_result['delta'] > 0:
            delta = best_result['delta']
            if delta >= epsilon / k:
                if best_result['is_left']:
                    current_left[best_result['remove_idx']] = best_result['add_idx']
                else:
                    current_right[best_result['remove_idx']] = best_result['add_idx']
                current_value = best_result['new_value']
        else:
            delta = 0

    return current_left, current_right



def approximation_algorithm(selected_left, selected_right, badness_components, block_indices, primer_sequences, epsilon):
    """
    Approximation algorithm adapted for new variables:
    - selected_left & selected_right: Initial selections
    - badness_components: New badness representation
    - block_indices: Block assignments
    - primer_sequences: Corresponding primer sequences
    - epsilon: Threshold
    """

    # **Step 1: First Local Search**
    S_opt1_left, S_opt1_right = fast_local_search_numba(
        selected_left, selected_right, badness_components, block_indices, primer_sequences, epsilon
    )
    
    # **Step 2: Generate S_R (Remaining Set)**
    total_primers = badness_components.shape[1] 
    S = np.arange(total_primers, dtype=np.int32)

    S_opt1_set = set(S_opt1_left) | set(S_opt1_right) 
    S_R = np.array([p for p in S if p not in S_opt1_set], dtype=np.int32)

    # **Step 3: Randomly Generate S_opt_R**
    S_opt_R_left = []
    S_opt_R_right = []
    unique_blocks = np.unique(block_indices)

    for block in unique_blocks:
        block_mask = (block_indices == block)
        block_candidates = S_R[np.isin(S_R, np.where(block_mask)[0])]

        if len(block_candidates) > 0:
            selected_left = np.random.choice(block_candidates)
            selected_right = np.random.choice(block_candidates)

            S_opt_R_left.append(selected_left)
            S_opt_R_right.append(selected_right)

    S_opt_R_left = np.array(S_opt_R_left, dtype=np.int32)
    S_opt_R_right = np.array(S_opt_R_right, dtype=np.int32)
    


    # **Step 4: Second Local Search**
    S_opt2_left, S_opt2_right = fast_local_search_numba(
        S_opt_R_left, S_opt_R_right, badness_components, block_indices, primer_sequences, epsilon
    )
  

    # **Step 5: Compute Objective Values**
    f_S_opt1 = calculate_objective_value_numba(S_opt1_left, badness_components, primer_sequences) + \
               calculate_objective_value_numba(S_opt1_right, badness_components, primer_sequences)

    f_S_opt2 = calculate_objective_value_numba(S_opt2_left, badness_components, primer_sequences) + \
               calculate_objective_value_numba(S_opt2_right, badness_components, primer_sequences)

    # **Step 6: Return the Better Set**
    return (S_opt1_left, S_opt1_right) if f_S_opt1 >= f_S_opt2 else (S_opt2_left, S_opt2_right)



