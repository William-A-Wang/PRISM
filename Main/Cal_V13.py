import pandas as pd
import numpy as np
import random

def extract_primers_with_unique_pair_id(file_path):
    primers_data = []
    current_id = None
    primer_left = None
    primer_right = None
    pair_counter = 0  # Counter for pair IDs within each SEQUENCE_ID

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()

            # Detect a new SEQUENCE_ID in the format SEQUENCE_ID=region_*
            if line.startswith("SEQUENCE_ID="):
                current_id = line.split("=")[-1].strip()
                pair_counter = 0  # Reset the pair counter for a new region

            # Extract LEFT and RIGHT sequences
            elif line.startswith("PRIMER_LEFT_") and "_SEQUENCE=" in line:
                primer_left = line.split("=")[-1].strip()
            elif line.startswith("PRIMER_RIGHT_") and "_SEQUENCE=" in line:
                primer_right = line.split("=")[-1].strip()

            # If both LEFT and RIGHT sequences are found, add them as a pair
            if primer_left and primer_right:
                pair_counter += 1  # Increment pair ID
                unique_pair_id = f"{current_id}_{pair_counter}"  # Create a unique ID
                primers_data.append({
                    "PAIR_ID": unique_pair_id,
                    "LEFT_SEQUENCE": primer_left,
                    "RIGHT_SEQUENCE": primer_right
                })
                primer_left = None
                primer_right = None

    # Convert the collected data into a DataFrame
    primer_df = pd.DataFrame(primers_data)
    return primer_df

# Example usage
input_file_path = "/content/drive/MyDrive/Colab Notebooks/primer3_output.txt"  # Replace with your file path
primer_tag = extract_primers_with_unique_pair_id(input_file_path)



### Badness
def is_complementary(seq1, seq2):
    return (seq1 == "A" and seq2 == "T") or \
           (seq1 == "T" and seq2 == "A") or \
           (seq1 == "C" and seq2 == "G") or \
           (seq1 == "G" and seq2 == "C")

def find_all_complementary_primers(primer1, primer2, min_len=4):
    m, n = len(primer1), len(primer2)
    dp = [[0] * (n+1) for _ in range(m+1)]
    results = []

    for i in range(1, m+1):
        for j in range(1, n+1):
            if is_complementary(primer1[i-1], primer2[j-1]):
                dp[i][j] = dp[i-1][j-1] + 1
                if dp[i][j] >= min_len:
                    start_position_primer1 = i - dp[i][j]
                    end_position_primer1 = i
                    start_position_primer2 = j - dp[i][j]
                    end_position_primer2 = j
                    results.append({
                        "Primer1_start": start_position_primer1,
                        "Primer1_end": end_position_primer1,
                        "Primer2_start": start_position_primer2,
                        "Primer2_end": end_position_primer2,
                        "Length": dp[i][j]
                    })
            else:
                dp[i][j] = 0

    results.sort(key=lambda x: (x['Primer1_start'], -x['Length']))

    filtered_results = []
    for i, result in enumerate(results):
        is_subset = False
        for j in range(i):
            if (results[j]['Primer1_start'] <= result['Primer1_start'] and
                results[j]['Primer1_end'] >= result['Primer1_end'] and
                results[j]['Primer2_start'] <= result['Primer2_start'] and
                results[j]['Primer2_end'] >= result['Primer2_end']):
                is_subset = True
                break
        if not is_subset:
            filtered_results.append(result)

    return filtered_results


### Badness calculate
def calculate_gc_pairs(primer1_seq, primer2_seq):
    num_gc = 0
    for b1, b2 in zip(primer1_seq, primer2_seq):
        if (b1 == 'G' and b2 == 'C') or (b1 == 'C' and b2 == 'G'):
            num_gc += 1
    return num_gc

def badness(result, sequence_length, primer1, primer2):
    # result
    length = result['Length']
    primer1_end = result['Primer1_end']
    primer2_end = result['Primer2_end']
    #d1 d2
    d1 = sequence_length - primer1_end
    d2 = sequence_length - primer2_end

    primer1_seq = primer1[result['Primer1_start']:result['Primer1_end'] + 1]
    primer2_seq = primer2[result['Primer2_start']:result['Primer2_end'] + 1]

    # numGC
    num_gc = calculate_gc_pairs(primer1_seq, primer2_seq)

    # badness
    badness_value = (2**length * 2**num_gc) / ((d1 + 1) * (d2 + 1))

    return badness_value

def New_Matrix(primer_df):
    num_primers = len(primer_df)
    badness_matrix = np.zeros((num_primers, num_primers))

    for i in range(num_primers):
        primer1 = primer_df.loc[i, 'Sequence']
        for j in range(num_primers):
            primer2 = primer_df.loc[j, 'Sequence']
            sequence_length = len(primer1)
            results = find_all_complementary_primers(primer1, primer2)
            total_badness = 0
            if results:
                for result in results:
                    badness_value = badness(result, sequence_length, primer1, primer2)
                    total_badness += badness_value
            badness_matrix[i, j] = total_badness

    return badness_matrix






###  yingshe
def map_pair_id_to_indices(primer_tag):
    pair_to_indices = {}
    for i, row in primer_tag.iterrows():
        left_index = 2 * i
        right_index = 2 * i + 1
        pair_to_indices[row['PAIR_ID']] = {left_index, right_index}
    return pair_to_indices


def map_pair_id_to_indices(primer_tag):
    if 'PAIR_ID' not in primer_tag.columns:
        raise ValueError("Input DataFrame must contain 'PAIR_ID' column.")

    left_indices = 2 * primer_tag.index
    right_indices = 2 * primer_tag.index + 1
    pair_to_indices = {
        pair_id: {left, right}
        for pair_id, left, right in zip(primer_tag['PAIR_ID'], left_indices, right_indices)
    }
    return pair_to_indices



pair_to_indices = map_pair_id_to_indices(primer_tag)

for pair_id, indices in pair_to_indices.items():
    print(f"{pair_id}: {indices}")


def generate_s_init_and_mapping(primer_tag, pair_to_indices):


    primer_tag['Region'] = primer_tag['PAIR_ID'].apply(lambda x: x.split('_')[1])


    S_init = set()
    S_init_m = []


    for region, group in primer_tag.groupby('Region'):
        random_pair_id = random.choice(group['PAIR_ID'].tolist())
        indices = pair_to_indices[random_pair_id]
        S_init.update(indices)
        S_init_m.append((random_pair_id, indices))

    return list(S_init), S_init_m



def generate_s_init_and_mapping(primer_tag, pair_to_indices):

    primer_tag['Region'] = primer_tag['PAIR_ID'].apply(lambda x: x.split('_')[1])


    S_init = set()
    S_init_m = []


    for region, group in primer_tag.groupby('Region'):
        # select
        random_pair_id = random.choice(group['PAIR_ID'].tolist())
        indices = pair_to_indices[random_pair_id]

        S_init.update(indices)
        S_init_m.append((random_pair_id, indices))


    return list(S_init), S_init_m


S_init, S_init_m = generate_s_init_and_mapping(primer_tag, pair_to_indices)

print("S_init:", S_init)
print("S_init_m (Fixed Mapping):")
for pair_id, indices in S_init_m:
    print(f"{pair_id}: {indices}")


### obj
def objective_function(S_init, badness_matrix):
    # function_1
    result = []
    for i in range(badness_matrix.shape[1]):
        max_value = -np.inf
        for s in S_init:
            row_index = s  #
            value = badness_matrix[row_index, i]  #
            if value > max_value:
                max_value = value
        result.append(max_value)
    f1 = sum(result)  # function_1

    # function_2
    total_sum = 0
    n = len(badness_matrix)

    for i in S_init:
        row_index_i = i  #
        for j in S_init:
            col_index_j = j  #
            total_sum += badness_matrix[row_index_i, col_index_j]

    f2 = (1 / n) * total_sum  # function_2

    return f1 - f2



def fast_local_search(S, S_init, epsilon, badness_df, S_init_m, pair_to_indices):

    k = len(S_init)
    delta = float('inf')

    print(f"Initial S_init: {S_init}")
    print(f"S_init_m: {S_init_m}")
    print(f"Candidate S: {S}")

    while delta >= epsilon / k:
        delta = -float('inf')
        best_remove_pair = None
        best_add_pair = None

        # Iterate over current pairs in S_init_m
        for remove_pair_id, remove_indices in S_init_m:
            for add_pair_id, add_indices in pair_to_indices.items():  # 动态从 pair_to_indices 中获取候选对
                if remove_pair_id == add_pair_id:
                    continue  # Skip replacing the same pair

                # Check if add_indices conflict with S_init
                if any(idx in S_init for idx in add_indices):
                    # print(f"Skipping Replacement: Add {add_indices} conflicts with S_init {S_init}")
                    continue  # Skip if any index of add_indices is already in S_init

                # Calculate the objective function for the current state
                g_pa = objective_function(S_init, badness_df)

                # Simulate replacing the pair in S_init
                new_S_opt = S_init.copy()
                new_S_opt = [idx for idx in new_S_opt if idx not in remove_indices]  # Remove pair
                new_S_opt.extend(add_indices)  # Add new pair

                # Calculate the objective function after replacement
                g_pb = objective_function(new_S_opt, badness_df)

                # Compute the delta
                current_delta = g_pa - g_pb

                # Debug: Print replacement attempt
                # print(f"Trying Replacement: Remove {remove_indices}, Add {add_indices}")
                # print(f"g_pa: {g_pa}, g_pb: {g_pb}, Current Delta: {current_delta}")

                # Find the largest delta
                if current_delta > delta:
                    delta = current_delta
                    best_remove_pair = remove_indices
                    best_add_pair = add_indices

                    # Debug: Print best replacement found
                    # print(f"Found Better Delta: {current_delta}")
                    # print(f"Best Remove Pair: {remove_indices}")
                    # print(f"Best Add Pair: {add_indices}")

        # If an improvement was found, apply the best replacement
        if delta >= epsilon / k and best_remove_pair and best_add_pair:
            S_init = [idx for idx in S_init if idx not in best_remove_pair]  # Remove the pair
            S_init.extend(best_add_pair)  # Add the new pair
            k = len(S_init)  # Update the size of S_init

            # Debug: Print updated S_init
            # print(f"Executing Replacement: Remove {best_remove_pair}, Add {best_add_pair}")
            # print(f"Updated S_init: {S_init}")

    return S_init
#### orignal version of appro
def approximation_algorithm(S, S_init, epsilon, badness_df, S_init_m, pair_to_indices):


    # Step 1: S_opt1
    S_opt1 = fast_local_search(S, S_init, epsilon, badness_df, S_init_m, pair_to_indices)
    print(f"S_opt1_output: {S_opt1}")

    # Step 2: Generate Second Input
    #
    S_R = [p for p in S if p not in S_opt1]

    S_R_pair_ids = [pair_id for pair_id, indices in pair_to_indices.items()
                    if indices.issubset(S_R)]

    S_R_pair = [
        (f"region_{primer_tag.loc[primer_tag['PAIR_ID'] == pair_id, 'Region'].values[0]}_{pair_id}", indices)
        for pair_id, indices in pair_to_indices.items()
        if indices.issubset(S_R)
    ]

    primer_tag['Region'] = primer_tag['PAIR_ID'].apply(lambda x: x.split('_')[1])


    S_opt_R = []
    for region, group in primer_tag[primer_tag['PAIR_ID'].isin(S_R_pair_ids)].groupby('Region'):

        if not group.empty:
            random_pair_id = random.choice(group['PAIR_ID'].tolist())
            S_opt_R.append(random_pair_id)

    S_opt2_m = [(pair_id, pair_to_indices[pair_id]) for pair_id in S_opt_R]


    S_opt_R_indices = set()
    for pair_id, indices in S_opt2_m:
        S_opt_R_indices.update(indices)
    print(f"S_opt2_random_input: {S_opt_R_indices}")

    # Step 3: S_opt2
    S_opt2 = fast_local_search(S_R, S_opt_R_indices, epsilon, badness_df, S_opt2_m, S_R_pair)
    print(f"S_opt2_output: {S_opt2}")

    # Step 4: Compare Objective Function
    f_S_opt1 = objective_function(S_opt1, badness_df)
    f_S_opt2 = objective_function(S_opt2, badness_df)

    if f_S_opt1 >= f_S_opt2:
        return S_opt1
    else:
        return S_opt2



print(f"S_init_output:{S_init}")
S = list(range(len(badness_df)))
epsilon = 0.01
final_S_opt = approximation_algorithm(S, S_init, epsilon, badness_df, S_init_m, pair_to_indices)
print(f"final select:{final_S_opt}")
