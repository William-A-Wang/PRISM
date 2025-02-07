from Bio import SeqIO
import os
from io import StringIO
import random
import primer3
import numpy as np
import pandas as pd
from numba import njit, prange
import matplotlib.pyplot as plt
from typing import List, Tuple, Set, Dict

# data process

def sliding_window_regions(sequence: str, window_sizes: List[int], overlaps: int = 250) -> Dict[int, List[Dict]]:
    all_regions = {}
    sequence_length = len(sequence)

    for window_size in window_sizes:
        step_size = window_size - overlaps
        start = 0
        regions = []

        while start < sequence_length:
            end = min(start + window_size, sequence_length)
           
            is_last_window = (end == sequence_length)

            region = {
                'id': f'size_{window_size}_region_{start}_{end}',
                'sequence': sequence[start:end],
                'start': start,
                'end': end,
                'length': end - start
            }
            regions.append(region)

           
            if is_last_window:
                break

            start += step_size  

        all_regions[window_size] = regions

    return all_regions

def design_primers(regions: dict) -> Dict[str, List[dict]]:

    primer_results = {}

    global_primer3_settings = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_NUM_RETURN': 5,
        'PRIMER_MIN_SIZE': 20,
        'PRIMER_MAX_SIZE': 26,
        'PRIMER_MIN_TM': 52.0,
        'PRIMER_MAX_TM': 62.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0
    }

    for window_size, region_list in regions.items():
        for region in region_list:
            seq_args = {
                'SEQUENCE_ID': region['id'],
                'SEQUENCE_TEMPLATE': region['sequence'],
                'SEQUENCE_INCLUDED_REGION': [0, len(region['sequence'])]
            }

            try:
                result = primer3.bindings.designPrimers(seq_args, global_primer3_settings)
                primer_pairs = []

                num_pairs = result.get('PRIMER_PAIR_NUM_RETURNED', 0)

                for i in range(num_pairs):
                    primer_info = {
                        'pair_id': f"{region['id']}_pair_{i}",
                        'left_primer': {
                            'sequence': result[f'PRIMER_LEFT_{i}_SEQUENCE'],
                            'tm': result[f'PRIMER_LEFT_{i}_TM'],
                            'gc_percent': result[f'PRIMER_LEFT_{i}_GC_PERCENT'],
                            'position': result[f'PRIMER_LEFT_{i}'][0],
                            'length': result[f'PRIMER_LEFT_{i}'][1]
                        },
                        'right_primer': {
                            'sequence': result[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                            'tm': result[f'PRIMER_RIGHT_{i}_TM'],
                            'gc_percent': result[f'PRIMER_RIGHT_{i}_GC_PERCENT'],
                            'position': result[f'PRIMER_RIGHT_{i}'][0],
                            'length': result[f'PRIMER_RIGHT_{i}'][1]
                        },
                        'product_size': result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                    }
                    primer_pairs.append(primer_info)

                primer_results[region['id']] = primer_pairs


                print(f"\n=== Primer Design for {region['id']} ===")
                if num_pairs > 0:
                    for primer in primer_pairs:
                        print(f"  Primer Pair {primer['pair_id']}:")
                        print(f"    Left: {primer['left_primer']['sequence']} (Tm: {primer['left_primer']['tm']:.1f}°C, GC: {primer['left_primer']['gc_percent']:.1f}%)")
                        print(f"    Right: {primer['right_primer']['sequence']} (Tm: {primer['right_primer']['tm']:.1f}°C, GC: {primer['right_primer']['gc_percent']:.1f}%)")
                        print(f"    Product Size: {primer['product_size']} bp")
                else:
                    print("  No primers found.")

            except Exception as e:
                print(f"Error in designing primers for {region['id']}: {str(e)}")

    return primer_results

def clean_sequence(sequence: str) -> str:
    """Clean the fasta sequence"""
    if sequence.startswith('>'):
        sequence = '\n'.join(sequence.split('\n')[1:])
    sequence = ''.join(sequence.split())
    return sequence

def process_consensus_sequence(sequence: str,
                             window_sizes: int = 1000,
                             overlaps: int = 250) -> dict:


    clean_seq = clean_sequence(sequence)
    print(f"Cleaned sequence length: {len(clean_seq)} bp")

    print(f"\nSliding window segmentation with overlap {overlaps}bp...")
    regions = sliding_window_regions(clean_seq, window_sizes, overlaps)
    print(regions)

    print("\nDesigning primers using Primer3...")
    primer_results = design_primers(regions)


    print("\nPrimer design statistics:")
    for window_size, region_list in regions.items():
        print(f"\nWindow size {window_size}bp:")
        success_count = sum(1 for region in region_list if len(primer_results.get(region['id'], [])) > 0)
        print(f"Regions with successful primer designs: {success_count}/{len(region_list)}")


        positions = [(r['start'], r['end']) for r in region_list]
        gaps = [positions[i+1][0] - positions[i][1] for i in range(len(positions)-1)]
        if gaps:
            avg_gap = sum(gaps) / len(gaps)
            min_gap_found = min(gaps)
            print(f"Average gap between regions: {avg_gap:.1f}bp")
            print(f"Minimum gap between regions: {min_gap_found}bp")

    return primer_results

def read_consensus_sequence(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]

        consensus_sequence = "".join([line.strip() for line in lines])

    return consensus_sequence

# badness define

@njit
def is_complementary_numba(base1: int, base2: int) -> bool:
    """Check if two nucleotides are complementary using ASCII values"""
    return ((base1 == 65 and base2 == 84) or  # A-T
            (base1 == 84 and base2 == 65) or  # T-A
            (base1 == 67 and base2 == 71) or  # C-G
            (base1 == 71 and base2 == 67))    # G-C

@njit
def is_same_region(region1: np.array, region2: np.array) -> bool:
    """Check if two regions represent the same complementary sequence
    (either directly or reversed)"""
    # Check if regions are identical in both directions
    direct_match = (region1[0] == region2[0] and
                   region1[1] == region2[1] and
                   region1[2] == region2[2] and
                   region1[3] == region2[3])

    reverse_match = (region1[0] == region2[2] and
                    region1[1] == region2[3] and
                    region1[2] == region2[0] and
                    region1[3] == region2[1])

    return direct_match or reverse_match

@njit
def find_complementary_regions_numba(primer1: np.array, primer2: np.array, min_len: int = 4):
    """Find all complementary regions between two primers using Numba.
    Returns a 2D array where each row contains [start1, end1, start2, end2, length]
    """
    m, n = len(primer1), len(primer2)
    dp = np.zeros((m + 1, n + 1), dtype=np.int32)

    # Pre-allocate maximum possible size for results
    max_possible_results = m * n
    all_results = np.zeros((max_possible_results, 5), dtype=np.int32)
    result_count = 0

    # First pass: find all complementary regions
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if is_complementary_numba(primer1[i-1], primer2[j-1]):
                dp[i, j] = dp[i-1, j-1] + 1
                if dp[i, j] >= min_len:
                    # Check if this region is already found (in either direction)
                    new_result = np.array([
                        i - dp[i, j] + 1,  # start1
                        i,             # end1
                        j - dp[i, j] + 1,  # start2
                        j,             # end2
                        dp[i, j]       # length
                    ])

                    # Check if this region is already in results
                    is_duplicate = False
                    for k in range(result_count):
                        if is_same_region(new_result, all_results[k]):
                            is_duplicate = True
                            break

                    if not is_duplicate:
                        all_results[result_count] = new_result
                        result_count += 1
            else:
                dp[i, j] = 0

    # If no results found, return empty array
    if result_count == 0:
        return np.zeros((0, 5), dtype=np.int32)

    # Get only the valid results
    valid_results = all_results[:result_count]

    # Sort by length (descending) using bubble sort
    for i in range(result_count):
        for j in range(0, result_count - i - 1):
            if (valid_results[j][4] < valid_results[j + 1][4] or  # Compare lengths
                (valid_results[j][4] == valid_results[j + 1][4] and  # If lengths equal
                 valid_results[j][0] > valid_results[j + 1][0])):    # Compare start positions
                # Swap rows
                valid_results[j], valid_results[j + 1] = valid_results[j + 1].copy(), valid_results[j].copy()

    # Filter overlapping regions
    used_regions = np.zeros((m, n), dtype=np.bool_)
    filtered_results = np.zeros((result_count, 5), dtype=np.int32)
    filtered_count = 0

    for i in range(result_count):
        result = valid_results[i]
        start1, end1, start2, end2 = int(result[0]), int(result[1]), int(result[2]), int(result[3])

        # Check for overlap
        is_overlap = False
        for x in range(start1 - 1, end1 - 1):
            for y in range(start2 - 1, end2 - 1):
                if used_regions[x, y]:
                    is_overlap = True
                    break
            if is_overlap:
                break

        if not is_overlap:
            # Mark region as used
            for x in range(start1 - 1, end1 - 1):
                for y in range(start2 - 1, end2 - 1):
                    used_regions[x, y] = True

            # Add to filtered results
            filtered_results[filtered_count] = result
            filtered_count += 1

    return filtered_results[:filtered_count]

@njit
def calculate_gc_pairs_numba(primer1_seq: np.array, primer2_seq: np.array) -> int:
    """Calculate number of G-C pairs using Numba"""
    num_gc = 0
    for i in range(len(primer1_seq)):
        if ((primer1_seq[i] == 71 and primer2_seq[i] == 67) or    # G-C
            (primer1_seq[i] == 67 and primer2_seq[i] == 71)):     # C-G
            num_gc += 1
    return num_gc

@njit
def calculate_badness_numba(result, sequence_length: int, primer1: np.array, primer2: np.array) -> float:
    """Calculate badness value using Numba"""
    start1, end1, start2, end2, length = result

    d1 = sequence_length - end1
    d2 = sequence_length - end2

    primer1_seq = primer1[start1 - 1:end1]
    primer2_seq = primer2[start2 - 1:end2]

    num_gc = calculate_gc_pairs_numba(primer1_seq, primer2_seq)

    return (2.0**length * 2.0**num_gc) / ((d1 + 1.0) * (d2 + 1.0))

@njit
def calculate_pair_badness_numba(primer1: np.array, primer2: np.array) -> float:
    """Calculate total badness between two primers using Numba"""
    sequence_length = len(primer1)
    results = find_complementary_regions_numba(primer1, primer2)
    total_badness = 0.0

    for result in results:
        badness_value = calculate_badness_numba(result, sequence_length, primer1, primer2)
        total_badness += badness_value

    return total_badness


def create_badness_matrix(primer_results):

    #  DataFrame
    primer_df = pd.DataFrame([
        {
            "PAIR_ID": pair["pair_id"],
            "LEFT_SEQUENCE": pair["left_primer"]["sequence"],
            "RIGHT_SEQUENCE": pair["right_primer"]["sequence"]
        }
        for pairs in primer_results.values()
        for pair in pairs
    ])

    num_sets = len(primer_df)
    badness_matrix = np.zeros((num_sets, num_sets))

    #  ASCII
    primer_sequences = {}
    for _, row in primer_df.iterrows():
        primer_sequences[f"{row['PAIR_ID']}_LEFT"] = np.array([ord(c) for c in row['LEFT_SEQUENCE']], dtype=np.int32)
        primer_sequences[f"{row['PAIR_ID']}_RIGHT"] = np.array([ord(c) for c in row['RIGHT_SEQUENCE']], dtype=np.int32)

    print(f"\nCreating {num_sets}x{num_sets} badness matrix...")

    #  badness_matrix
    for i in range(num_sets):
        row = primer_df.iloc[i]
        set1_left = primer_sequences[f"{row['PAIR_ID']}_LEFT"]
        set1_right = primer_sequences[f"{row['PAIR_ID']}_RIGHT"]

        for j in range(num_sets):
            row2 = primer_df.iloc[j]
            set2_left = primer_sequences[f"{row2['PAIR_ID']}_LEFT"]
            set2_right = primer_sequences[f"{row2['PAIR_ID']}_RIGHT"]

            # badness
            total_badness = (
                calculate_pair_badness_numba(set1_left, set2_left) +   # Left1-Left2
                calculate_pair_badness_numba(set1_left, set2_right) +  # Left1-Right2
                calculate_pair_badness_numba(set1_right, set2_left) +  # Right1-Left2
                calculate_pair_badness_numba(set1_right, set2_right) + # Right1-Right2
                calculate_pair_badness_numba(set1_left, set1_left) +   # Left1-Left1
                calculate_pair_badness_numba(set1_right, set1_right) + # Right1-Right1
                calculate_pair_badness_numba(set2_left, set2_left) +   # Left2-Left2
                calculate_pair_badness_numba(set2_right, set2_right)   # Right2-Right2
            )
            if i == j:
                total_badness = (calculate_pair_badness_numba(set1_left, set2_left) +
                                calculate_pair_badness_numba(set1_right, set2_right)
                )

            badness_matrix[i, j] = total_badness

        if (i + 1) % 10 == 0:
            print(f"Progress: {((i + 1) / num_sets * 100):.1f}%")

    return badness_matrix


# generate S_init

def compute_badness_for_blocks(primer_results: Dict[str, List[Dict]]) -> Dict[str, float]:
    block_badness = {}

    for block_id, primers in primer_results.items():
        if len(primers) < 2:
            print(f"Warning: Block {block_id} has less than 2 primers, skipping badness calculation.")
            block_badness[block_id] = 0
            continue

        block_primer_results = {block_id: primers}
        badness_matrix = create_badness_matrix(block_primer_results)

        total_badness = np.sum(badness_matrix) / 2
        block_badness[block_id] = total_badness

        print(f"Block: {block_id}, Badness Total: {total_badness}")

    return block_badness

def adjust_blocks(primer_results: Dict[str, List[Dict]], block_badness: Dict[str, float], sequence: str, extend_size=100, threshold=600, min_block_size=200):
   
    sequence_length = len(sequence)  

    # mean
    mean_badness = np.mean(list(block_badness.values()))

    # max
    max_badness_block = max(block_badness, key=block_badness.get)
    max_badness_value = block_badness[max_badness_block]
    diff = max_badness_value - mean_badness

    print(f"\n🔍 Mean Badness: {mean_badness:.2f}")
    print(f"🔥 Block with Highest Badness: {max_badness_block} → {max_badness_value:.2f} (Diff: {diff:.2f})")

    # adjust
    if diff < threshold:
        print(f"✅ Badness is within acceptable range. No expansion needed.")
        return primer_results  

  
    parts = max_badness_block.split("_")
    target_start, target_end = int(parts[-2]), int(parts[-1])

    # new
    new_end = min(target_end + extend_size, sequence_length)  # sequence_length
    new_target_block = f"size_1000_region_{target_start}_{new_end}"
    print(f"\n🔄 Expanding {max_badness_block} → {new_target_block}")

    updated_primer_results = {}

    # 
    block_ids = sorted(primer_results.keys(), key=lambda x: int(x.split("_")[-2]))
    last_block = block_ids[-1]  

    for block_id in block_ids:
        start, end = map(int, block_id.split("_")[-2:])

        if block_id == max_badness_block:
           
            updated_primer_results[new_target_block] = primer_results[block_id]

        elif start > target_start:  #shift
            new_start = start + extend_size
            new_end = min(end + extend_size, sequence_length)  
            new_block_id = f"size_1000_region_{new_start}_{new_end}"
            updated_primer_results[new_block_id] = primer_results[block_id]
            print(f"📌 Moving {block_id} → {new_block_id}")

        else:
            
            updated_primer_results[block_id] = primer_results[block_id]

   
    sorted_primer_results = dict(sorted(updated_primer_results.items(), key=lambda x: int(x[0].split("_")[-2])))


    return sorted_primer_results 


def convert_primer_results_to_regions(primer_results: Dict[str, List[Dict]], sequence: str) -> Dict[int, List[Dict]]:
    
    regions = {}

    for block_id in primer_results.keys():
        parts = block_id.split("_")

        try:
            start, end = int(parts[-2]), int(parts[-1])  
        except ValueError:
            print(f"⚠ Warning: Invalid region format: {block_id}, skipping...")
            continue  

        window_size = end - start  
        region_data = {
            "id": block_id,
            "sequence": sequence[start:end]  
        }

        if window_size not in regions:
            regions[window_size] = []
        regions[window_size].append(region_data)

    return regions


def iterative_primer_optimization(sequence: str, regions: Dict[int, List[Dict]], max_iterations=20, extend_size=100, threshold=600):
    
    iteration = 1
    primer_results = design_primers(regions) 

    while iteration <= max_iterations:
        print(f"\n🔄 Iteration {iteration}: Computing Badness...")

        # badness
        block_badness = compute_badness_for_blocks(primer_results)

        # block
        new_primer_results = adjust_blocks(primer_results, block_badness, sequence, extend_size, threshold)


    
        if new_primer_results == primer_results:
            print(f"✅ No further adjustments needed. Stopping at iteration {iteration}.")
            break

       
        new_regions = convert_primer_results_to_regions(new_primer_results, sequence)

        print(f"\n⚡ Generating New Primers with Primer3...")
        primer_results = design_primers(new_regions)

        # renew regions
        regions = new_regions

        iteration += 1

    
    primer_results = dict(sorted(primer_results.items(), key=lambda x: int(x[0].split("_")[-2])))

    return primer_results

def generate_S_init(final_result: Dict[str, List[Dict]], df_badness: pd.DataFrame) -> List[str]:

    S_init = []

    for block_id, primer_list in final_result.items():
        max_f1_value = -np.inf
        best_pair_id = None

      
        for primer in primer_list:
            pair_id = primer["pair_id"]

            if pair_id in df_badness.index:
                f1_value = df_badness.loc[pair_id, :].max()

                print(f"Block: {block_id}, Pair ID: {pair_id}, f1 Value: {f1_value}")

                if f1_value > max_f1_value:
                    max_f1_value = f1_value
                    best_pair_id = pair_id

        if best_pair_id is not None:
            S_init.append(best_pair_id)

    return S_init

# obj

def objective_function(S_init, df_badness):
   
    result = []

  
    available_pairs = set(df_badness.index)


    S_init_filtered = [pair for pair in S_init if pair in available_pairs]

    
    if not S_init_filtered:
        print("⚠️ Warning: None of S_init exists in df_badness!")
        return 0

   
    for col in df_badness.columns:
        max_value = df_badness.loc[S_init_filtered, col].max()
        result.append(max_value)
    f1 = sum(result)

    
    f2 = df_badness.loc[S_init_filtered, S_init_filtered].sum().sum() / len(df_badness)

    return f1 - f2

#local search

def fast_local_search(S, S_init, epsilon, df_badness):
   
    k = len(S_init)
    delta = float('inf')

    print(f"Initial S_init: {S_init}")
    print(f"Candidate S: {S}")

  
    available_pairs = set(df_badness.index)
    S_filtered = [p for p in S if p in available_pairs]
    S_init_filtered = [p for p in S_init if p in available_pairs]

    if not S_filtered or not S_init_filtered:
        print("⚠️ Warning: No valid primers in S or S_init!")
        return S_init

    while delta >= epsilon / k:
        delta = -float('inf')
        best_remove_pair = None
        best_add_pair = None

        # block_id replacement
        for remove_pair_id in S_init_filtered:
            block_id = "_".join(remove_pair_id.split("_")[:4]) 
            block_candidates = [p for p in S_filtered if p.startswith(block_id)]  

            for add_pair_id in block_candidates:
                if remove_pair_id == add_pair_id:
                    continue

               
                g_pa = objective_function(S_init_filtered, df_badness)

                
                new_S_opt = S_init_filtered.copy()
                new_S_opt.remove(remove_pair_id)
                new_S_opt.append(add_pair_id)

             
                g_pb = objective_function(new_S_opt, df_badness)

                
                current_delta = g_pa - g_pb

             
                if current_delta > delta:
                    delta = current_delta
                    best_remove_pair = remove_pair_id
                    best_add_pair = add_pair_id

        # delta
        if delta >= epsilon / k and best_remove_pair and best_add_pair:
            S_init_filtered.remove(best_remove_pair)
            S_init_filtered.append(best_add_pair)
            k = len(S_init_filtered)

            print(f"🔄 Replaced {best_remove_pair} with {best_add_pair} in block {block_id}")
            print(f"Updated S_init: {S_init_filtered}")

    return S_init_filtered


def approximation_algorithm(S, S_init, epsilon, df_badness):
  
    S_opt1 = fast_local_search(S, S_init, epsilon, df_badness)
    print(f"S_opt1_output: {S_opt1}")

    # `S_R`
    S_R = [p for p in S if p not in S_opt1]

    # random
    S_opt_R = []
    for region in set([pair.split("_pair_")[0] for pair in S_R]):  
        region_primers = [p for p in S_R if p.startswith(region)]
        if region_primers:
            S_opt_R.append(random.choice(region_primers))

    print(f"S_opt2_random_input: {S_opt_R}")

   
    S_opt2 = fast_local_search(S_R, S_opt_R, epsilon, df_badness)
    print(f"S_opt2_output: {S_opt2}")

    f_S_opt1 = objective_function(S_opt1, df_badness)
    f_S_opt2 = objective_function(S_opt2, df_badness)

    return S_opt1 if f_S_opt1 >= f_S_opt2 else S_opt2

def main():
    # data
    file_path = "/content/drive/MyDrive/ConsensusSequence.fasta" # 
    sequence = read_consensus_sequence(file_path)

    # regions
    window_size = 1000 
    overlap = 250  
    regions = sliding_window_regions(sequence, [window_size], overlap)

    # primer3
    primer_results = design_primers(regions)

    # iterative_primer_optimization
    regions = convert_primer_results_to_regions(primer_results, sequence)
    max_iterations = 20  
    extend_size = 100    
    threshold = 600      
    print("\n🔍 Running iterative primer optimization...")
    optimized_primers = iterative_primer_optimization(sequence, regions)

    # df_badness
    print("\n🔍 Computing badness matrix...")
    df_badness = pd.DataFrame(create_badness_matrix(optimized_primers),
                              index=[pair["pair_id"] for pairs in optimized_primers.values() for pair in pairs],
                              columns=[pair["pair_id"] for pairs in optimized_primers.values() for pair in pairs])
    
    # S_init
    print("\n🔍 Generating initial S_init...")
    S_init = generate_S_init(optimized_primers, df_badness)
    
    # S
    S = list(df_badness.index)

    # epsilon
    epsilon = 0.1

    # approximation_algorithm
    print("\n⚡ Running Approximation Algorithm...")
    optimized_S = approximation_algorithm(S, S_init, epsilon, df_badness)

    # final
    print("\n🎯 Final Optimized Primers:")
    print(optimized_S)



if __name__ == "__main__":
    main()
