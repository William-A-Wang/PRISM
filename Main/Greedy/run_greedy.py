import subprocess

def run_programs():
  
    subprocess.run(["python3", "Extract_V3.py"])
    print("Program 1 finished, E-output.txt generated.")

    subprocess.run(["python3", "Badness_Cal_V6.py"])
    print("Program 2 finished, badness_matrix.csv generated.")

    subprocess.run(["python3", "Greedy_V5.py"])
    print("Program 3 finished, final processing done.")

    try:
        with open("optimized_primers.txt", "r") as file:
            content = file.read()
            print("Content of optimized_primers.txt:")
            print(content)
    except FileNotFoundError:
        print("Error: optimized_primers.txt not found.")

if __name__ == "__main__":
    run_programs()
