import numpy as np
import pandas as pd
# 读取引物列表的函数
def load_primers(file_path):
    primers = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # 提取并存储每条引物
            parts = line.strip().split()
            if len(parts) == 2:
                primers.append(parts[1])
    return primers

# 反向互补序列计算
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# Badness计算函数
def calculate_badness(p1, p2):
    min_len = 4  # 最少4个连续碱基反向互补才能计算badness
    badness = 0

    # 反向互补的p2序列，如果p1 == p2，那么p2_rc是p1自身的反向互补
    p2_rc = reverse_complement(p2)

    # 遍历引物p1，p2之间的每个可能的反向互补序列
    for i in range(len(p1) - min_len + 1):
        for j in range(len(p2_rc) - min_len + 1):
            # 从i, j位置开始，找到最小长度为min_len的反向互补序列
            length = 0
            while (i + length < len(p1)) and (j + length < len(p2_rc)) and (p1[i + length] == p2_rc[j + length]):
                length += 1
                if length >= min_len:  # 找到长度大于等于min_len的反向互补序列
                    # 计算d1 和 d2, 到3'末端的距离
                    d1 = len(p1) - (i + length)
                    d2 = len(p2_rc) - (j + length)

                    # 计算numGC，统计反向互补序列中的G/C碱基对数量
                    numGC = sum(1 for k in range(length) if p1[i + k] in "GC")

                    # 根据公式计算当前反向互补序列的badness贡献
                    badness_contribution = (2 ** length) * (2 ** numGC) / ((d1 + 1) * (d2 + 1))
                    badness += badness_contribution

    return badness

# 加载引物并计算每对引物的Badness
primers = load_primers("D:\\CISEProgram\\USDA\Greedy\\E_output2.txt")

# 初始化n x n的badness矩阵
n = len(primers)
badness_matrix = np.zeros((n, n))

# 计算并存储每对引物的badness值，包括自身与自身的badness
for i in range(n):
    for j in range(i, n):  # 这里包括i == j的情况
        p1 = primers[i]
        p2 = primers[j]
        badness_value = calculate_badness(p1, p2)
        # 存储到矩阵中
        badness_matrix[i, j] = badness_value
        badness_matrix[j, i] = badness_value  # 矩阵是对称的

# 打印badness矩阵
print("Badness Matrix:")
print(badness_matrix)


# 使用 pandas 将矩阵转换为 DataFrame
df = pd.DataFrame(badness_matrix, index=[f'Primer{i+1}' for i in range(n)], columns=[f'Primer{i+1}' for i in range(n)])

# 将矩阵保存为 CSV 文件
df.to_csv('D:\\CISEProgram\\USDA\Greedy\\badness_matrix2.csv')

print("Badness matrix has been saved to 'badness_matrix2.csv'")

np.save('D:\\CISEProgram\\USDA\Greedy\\badness_matrix2.npy', badness_matrix)

print("Badness matrix has been saved to 'badness_matrix2.npy'")