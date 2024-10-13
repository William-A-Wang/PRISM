import pandas as pd
#1:
badness_matrix_file = 'D:\\CISEProgram\\USDA\\Greedy\\badness_matrix.csv'
badness_matrix = pd.read_csv(badness_matrix_file, index_col=0)

#2: 
primer_file_path = 'D:\\CISEProgram\\USDA\\Greedy\\E_output.txt'
primer_sequences = {}


with open(primer_file_path, 'r') as file:
    for index, line in enumerate(file, start=1):
        primer_name = f"Primer{index}"
        primer_sequence = line.strip() 
        primer_sequences[primer_name] = primer_sequence

#3: 
delta_thres = 10.0  

#4: 
primers = badness_matrix.columns.tolist()
S_opt = [] 
delta_min = float('inf')

#5: 
while len(primers) > 0:
    delta_min = float('inf')
    selected_primer = None

    for t in primers:
        delta_L = 0
        if len(S_opt) > 0:
            for q in S_opt:
                delta_L += badness_matrix.loc[t, q]
        
        delta_L += badness_matrix.loc[t, t] 

        ####
        if delta_L < delta_min or delta_L == 0:
            delta_min = delta_L
            selected_primer = t

    #6 
    if delta_min <= delta_thres and selected_primer is not None:
        S_opt.append(selected_primer)
        primers.remove(selected_primer)
    else:
        break


print("Optimized Primer Set with Sequences:")
for primer in S_opt:
    sequence = primer_sequences.get(primer, 'Unknown Sequence') 
    print(f"{primer}: {sequence}")
