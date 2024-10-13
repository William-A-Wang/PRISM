import subprocess
import os

def run_programs():
    python_venv = "D:\\CISEProgram\\USDA\\.venv\\Scripts\\python.exe"  

    subprocess.run([python_venv, "D:\\CISEProgram\\USDA\\Greedy\\Extract_V3.py"], check=True)
    print("Program 1 finished, E-output.txt generated.")

    if not os.path.exists("D:\\CISEProgram\\USDA\\Greedy\\E_output.txt"):
        print("Error: E_output.txt not found!")
        return

    subprocess.run([python_venv, "D:\\CISEProgram\\USDA\\Greedy\\Badness_Cal_V6.py"], check=True)
    print("Program 2 finished, badness_matrix.csv generated.")

    if not os.path.exists("D:\\CISEProgram\\USDA\\Greedy\\badness_matrix.csv"):
        print("Error: badness_matrix.csv not found!")
        return

    subprocess.run([python_venv, "D:\\CISEProgram\\USDA\\Greedy\\Greedy_V5.py"], check=True)
    print("Program 3 finished, final processing done.")

if __name__ == "__main__":
    run_programs()
