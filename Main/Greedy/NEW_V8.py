import pandas as pd
import numpy as np

#Extract the data
def extract_primers(file_path):
    primers = []
    
    with open(file_path, 'r') as f:
        primer_left = None
        primer_right = None
        for line in f:
            line = line.strip()
            if line.startswith("PRIMER_LEFT_") and "_SEQUENCE=" in line:
                primer_left = line.split("=")[-1].strip()
            elif line.startswith("PRIMER_RIGHT_") and "_SEQUENCE=" in line:
                primer_right = line.split("=")[-1].strip()
            #merge the left and right 
            if primer_left and primer_right:
                primers.append(primer_left)  
                primers.append(primer_right)  
                primer_left = None
                primer_right = None

    return primers

input_file_path = "D:\\CISEProgram\\USDA\\Greedy\\primer3_output.txt" #Your path
primers = extract_primers(input_file_path)
# print(f"primers:")
# print(extract_primers(input_file_path))

primer_df = pd.DataFrame({
    "Sequence": primers
})
# print(primer_df)

######上面是数据处理
# Badness 

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


# Badness calculate
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

# Matrix 
def Matrix(primer_df):
    num_primers = len(primer_df)  
    badness_matrix = np.zeros((num_primers, num_primers)) 

    for i in range(num_primers):
        primer1 = primer_df.loc[i, 'Sequence'] 

        for j in range(num_primers):  
            primer2 = primer_df.loc[j, 'Sequence'][::-1]  

            sequence_length = len(primer1)  

            results = find_all_complementary_primers(primer1, primer2)

            total_badness = 0  

            if results:
                for result in results:
                    
                    badness_value = badness(result, sequence_length, primer1, primer2)
                    total_badness += badness_value  

            badness_matrix[i, j] = total_badness

    # Matrix turn to dataFrame
    badness_df = pd.DataFrame(badness_matrix, columns=primer_df['Sequence'], index=primer_df['Sequence'])
    return badness_df

# Test
badness_df = Matrix(primer_df)
# print(badness_df)

import random

def greedy_primer_optimization(badness_df, delta_thres):
    S_opt = [] 
    S = list(range(len(badness_df)))  
    delta_min = np.inf
    
    first_index = random.choice(S)
    S_opt.append(first_index)

    S.remove(first_index)

    S_random = np.random.permutation(S)
    for new_index in S_random:
        delta_temp = 0  
        for selected in S_opt:
            delta_temp += badness_df.iloc[selected, new_index]  
            delta_temp += badness_df.iloc[new_index, selected] 
        
        delta_temp += badness_df.iloc[new_index, new_index]  
        
      
        if delta_temp < delta_min and delta_temp > delta_thres:
            delta_min = delta_temp
            S_opt.append(new_index)  
        
        if delta_min <= delta_thres:
            break
    
    return S_opt  


# Test
delta_thres = 50
S_opt = greedy_primer_optimization(badness_df, delta_thres)
S_opt = [{int(i)} for i in S_opt]
print("Final selected S_opt:", S_opt)

#### Local search

# # function1
# def calculate_function_1(S_opt, badness_df):
    
#     result = []  
#     for i in range(badness_df.shape[1]):  
#         max_value = -np.inf  
#         for s in S_opt:  
#             row_index = list(s)[0]  
#             value = badness_df.iloc[row_index, i]  
#             if value > max_value:
#                 max_value = value  
#         result.append(max_value)  
#     return sum(result)  

# # Test
# function_1_result = calculate_function_1(S_opt, badness_df)
# print(f"Result: {function_1_result}")   


# # function2
# def calculate_function_2(S_opt, badness_df):
  
#     total_sum = 0
#     n = len(badness_df)  
    
#     # All S_opt 
#     for i in S_opt:
#         row_index_i = list(i)[0]  
#         for j in S_opt:
#             col_index_j = list(j)[0]  
#             total_sum += badness_df.iloc[row_index_i, col_index_j]  
    
#     # 1/n
#     function_2_value = (1 / n) * total_sum
#     return function_2_value

# # Test
# function_2_result = calculate_function_2(S_opt, badness_df)
# print(f"Result: {function_2_result}")

# effi = function_1_result - function_2_result
# print(effi)

# def fast_local_search(S, S_opt, epsilon, badness_df):
#     k = len(S_opt)  #
#     delta = float('inf')  

#     while delta >= epsilon / k:
#         delta = -float('inf')  # reset
#         best_pa = None
#         best_pb = None

#         # 
#         for pa in S_opt:
#             for pb in S:
#                 if {pb} not in S_opt:  
#                     # g(pa) = f1(pa) - f2(pa) 和 g(pb) = f1(pb) - f2(pb)
#                     g_pa = calculate_function_1(S_opt, badness_df) - calculate_function_2(S_opt, badness_df) 
#                     new_S_opt = S_opt.copy()
#                     new_S_opt.remove(pa)  
#                     new_S_opt.append({pb}) 
#                     g_pb = calculate_function_1(new_S_opt, badness_df) - calculate_function_2(new_S_opt, badness_df) 
                    
#                     # delta
#                     current_delta = g_pa - g_pb

#                     # best pa and pb
#                     if current_delta > delta:
#                         delta = current_delta
#                         best_pa = pa
#                         best_pb = pb

#        
#         if delta >= epsilon / k:
#             S_opt.remove(best_pa)
#             S_opt.append({best_pb})
#             k = len(S_opt)  
#     return S_opt

# # 
# S = list(range(len(badness_df)))  
# epsilon = 0.01  
# optimized_S_opt = fast_local_search(S, S_opt, epsilon, badness_df)
# print(f"NEW S_opt: {optimized_S_opt}")























