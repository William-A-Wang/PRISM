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
            
            # extract
            if primer_left and primer_right:
                primers.append(primer_left)  
                primers.append(primer_right)  
                primer_left = None
                primer_right = None

    return primers

def save_primers_to_file(primers, output_file):
    with open(output_file, 'w') as f:
        for i, primer in enumerate(primers, start=1):
            f.write(f"Primer{i}: {primer}\n")

#save
input_file_path = "D:\\CISEProgram\\USDA\\USDA_Project\\primer3_output.txt"  
output_file_path = "D:\\CISEProgram\\USDA\Greedy\\E_output.txt"  

primers = extract_primers(input_file_path)
save_primers_to_file(primers, output_file_path)

print(f"Primer save {output_file_path}")
