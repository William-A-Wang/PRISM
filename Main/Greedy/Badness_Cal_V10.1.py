import pandas as pd
import numpy as np
import random

### Extract the data
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

### Matrix 
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


### Greedy
def greedy_primer_optimization(badness_df, delta_thres, debug=False):
    S_opt = []  
    S = list(range(len(badness_df))) 
    # first
    first_index = random.choice(S)
    S_opt.append(first_index)
    S.remove(first_index)
    if debug:
        print("Starting with first selected primer:", first_index)

    # Until delta_min >= delta_thres
    while True:
        delta_min = np.inf  
        selected_s = None  

        # First Loop
        for new_index in S:
            
            delta_temp = badness_df.iloc[new_index, new_index]  

            for selected in S_opt:
                delta_temp += badness_df.iloc[selected, new_index]  
                delta_temp += badness_df.iloc[new_index, selected]  

            if delta_temp < delta_min:
                delta_min = delta_temp
                selected_s = new_index

       
        if delta_min >= delta_thres:
            if debug:
                print(f"Terminating loop as delta_min ({delta_min}) >= delta_thres ({delta_thres})")
            break

        if selected_s is not None:
            S_opt.append(selected_s)
            S.remove(selected_s)  

        if debug:
            print(f"Selected primer {selected_s} with delta_min = {delta_min}")

    return S_opt  


### obj
def objective_function(S_init, badness_df):
    # function_1
    result = []  
    for i in range(badness_df.shape[1]):  
        max_value = -np.inf  
        for s in S_init:  
            row_index = list(s)[0]  
            value = badness_df.iloc[row_index, i]  
            if value > max_value:
                max_value = value  
        result.append(max_value)  
    f1 = sum(result)  # function_1 

    # function_2
    total_sum = 0
    n = len(badness_df)  
    
    for i in S_init:
        row_index_i = list(i)[0]  
        for j in S_init:
            col_index_j = list(j)[0]  
            total_sum += badness_df.iloc[row_index_i, col_index_j]  
    
    f2 = (1 / n) * total_sum  # function_2 

    return f1 - f2


### local search
def fast_local_search(S, S_init, epsilon, badness_df):

    k = len(S_init)  
    delta = float('inf')  

    while delta >= epsilon / k:
        delta = -float('inf')  
        best_pa = None
        best_pb = None

        for pa in S_init:
           
            for pb in S:
                if {pb} not in S_init:  
                    
                    g_pa = objective_function(S_init, badness_df)  
                    
                    new_S_opt = S_init.copy()
                    new_S_opt.remove(pa)
                    new_S_opt.append({pb})  # 
                    
                    g_pb = objective_function(new_S_opt, badness_df) 

                    current_delta = g_pa - g_pb

                    # find biggest delta
                    if current_delta > delta:
                        delta = current_delta
                        best_pa = pa
                        best_pb = pb

        if delta >= epsilon / k:
            S_init.remove(best_pa)  
            S_init.append({best_pb})  
            k = len(S_init)  

    return S_init  

### Approximation

def approximation_algorithm(S, S_init, epsilon, badness_df):
    # S_opt1
    S_opt1 = fast_local_search(S, S_init, epsilon, badness_df)

    # Second input
    S_R = [p for p in S if p not in S_opt1] 
    K = len(S_opt1)  # 
    S_opt_R = random.sample(S, K)

    # S_opt2
    S_opt2 = fast_local_search(S_R, S_opt_R, epsilon, badness_df)

    # Compare
    f_S_opt1 = objective_function(S_opt1, badness_df)
    f_S_opt2 = objective_function(S_opt2, badness_df)

    if f_S_opt1 >= f_S_opt2:
        return S_opt1
    else:
        return S_opt2



def main():
    # 1
    input_file_path = "D:\\CISEProgram\\USDA\\Greedy\\primer3_output.txt" 
    primers = extract_primers(input_file_path)
    
    # 2
    primer_df = pd.DataFrame({
        "Sequence": primers
    })

    badness_df = Matrix(primer_df)
    
    # 3
    delta_thres = 10  
    S_init = greedy_primer_optimization(badness_df, delta_thres, debug=True)
    S_init = [{int(i)} for i in S_init]
    print("Final selected S_opt:", S_init)
    
    # 4&5
    S = list(range(len(badness_df)))  
    epsilon = 0.01  
    final_S_opt = approximation_algorithm(S, S_init, epsilon, badness_df)
    
    
    
    print(f"Final Optimized Set: {final_S_opt}")

if __name__ == "__main__":
    main()
