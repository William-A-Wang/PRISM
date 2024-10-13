import numpy as np
import pandas as pd
#read
def load_primers(file_path):
    primers = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            parts = line.strip().split()
            if len(parts) == 2:
                primers.append(parts[1])
    return primers

#reverse
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# Badness计算函数
def calculate_badness(p1, p2):
    min_len = 4  # min4
    badness = 0
    p2_rc = reverse_complement(p2)

    # main
    for i in range(len(p1) - min_len + 1):
        for j in range(len(p2_rc) - min_len + 1):
            
            length = 0
            while (i + length < len(p1)) and (j + length < len(p2_rc)) and (p1[i + length] == p2_rc[j + length]):
                length += 1
                if length >= min_len:  
                    # d1 and d2, to 3'
                    d1 = len(p1) - (i + length)
                    d2 = len(p2_rc) - (j + length)
                    
                    numGC = sum(1 for k in range(length) if p1[i + k] in "GC")
                    
                    badness_contribution = (2 ** length) * (2 ** numGC) / ((d1 + 1) * (d2 + 1))
                    badness += badness_contribution

    return badness

primers = load_primers("D:\\CISEProgram\\USDA\Greedy\\E_output2.txt")
n = len(primers)
badness_matrix = np.zeros((n, n))


for i in range(n):
    for j in range(i, n):  # i == j
        p1 = primers[i]
        p2 = primers[j]
        badness_value = calculate_badness(p1, p2)
       
        badness_matrix[i, j] = badness_value
        badness_matrix[j, i] = badness_value 

print("Badness Matrix:")
print(badness_matrix)


# save
df = pd.DataFrame(badness_matrix, index=[f'Primer{i+1}' for i in range(n)], columns=[f'Primer{i+1}' for i in range(n)])

df.to_csv('D:\\CISEProgram\\USDA\Greedy\\badness_matrix2.csv')

print("Badness matrix has been saved to 'badness_matrix2.csv'")

np.save('D:\\CISEProgram\\USDA\Greedy\\badness_matrix2.npy', badness_matrix)

print("Badness matrix has been saved to 'badness_matrix2.npy'")
