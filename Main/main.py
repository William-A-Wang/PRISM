import argparse
from tqdm import tqdm
import pandas as pd

from data_io import read_consensus_sequence, sliding_window_regions
from design_primers import design_primers, PrimerSetBadnessFast
from iterative import convert_primer_results_to_regions, iterative_primer_optimization
from badness_utils import compute_badness, process_badness_mapping
from optimization import generate_initial_solution_numba, approximation_algorithm
import warnings

warnings.filterwarnings("ignore", message="Function deprecated please use")

def parse_args():
    """
    Parse command-line arguments:
      -i / --input       Path to the input FASTA file
      -w / --window-size Window size for slicing (default: 250)
      -o / --output-csv  Path to save the optimized primers CSV
    """
    parser = argparse.ArgumentParser(description="Run the primer optimization pipeline")
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to the input FASTA file, e.g. D:\\data\\seq.fna"
    )
    parser.add_argument(
        "-w", "--window-size", type=int, default=250,
        help="Window size for region slicing (default: 250)"
    )
    parser.add_argument(
        "-o", "--output-csv", default="optimized_primers.csv",
        help="Output CSV for optimized primers (default: optimized_primers.csv)"
    )
    return parser.parse_args()

def main():
    args = parse_args()

    steps = [
        "Load sequence",
        "Generate regions",
        "Design primers with Primer3",
        "Iterative optimization",
        "Compute initial badness",
        "Generate initial solution",
        "Approximation optimization",
        "Save results"
    ]
    pbar = tqdm(total=len(steps), desc="Design and Optimization Progress")

    # Step 1: Load sequence from FASTA file
    sequence = read_consensus_sequence(args.input)
    pbar.set_description(steps[0])
    pbar.update(1)

    # Step 2: Generate overlapping regions
    regions = sliding_window_regions(sequence, [args.window_size], overlaps=0)
    pbar.set_description(steps[1])
    pbar.update(1)

    # Step 3: Design primers for each region
    primer_results = design_primers(regions)
    pbar.set_description(steps[2])
    pbar.update(1)

    # Step 4: Perform iterative optimization
    regions = convert_primer_results_to_regions(primer_results, sequence)
    optimized = iterative_primer_optimization(
        sequence,
        regions,
        max_iterations=50,
        extend_size=100,
        k=2.0
    )
    pbar.set_description(steps[3])
    pbar.update(1)

    # Step 5: Compute total badness and mapping
    total_badness, mapping = compute_badness(optimized)
    pbar.set_description(steps[4])
    pbar.update(1)

    # Step 6: Generate initial solution with numba
    badness_components, block_indices, primer_sequences = process_badness_mapping(mapping, optimized)
    sel_left, sel_right = generate_initial_solution_numba(badness_components, block_indices)
    pbar.set_description(steps[5])
    pbar.update(1)

    # Step 7: Run approximation algorithm
    final_left, final_right = approximation_algorithm(
        sel_left,
        sel_right,
        badness_components,
        block_indices,
        primer_sequences,
        epsilon=0.1
    )
    pbar.set_description(steps[6])
    pbar.update(1)

    # Step 8: Map final indices back to sequences and save CSV
    final_left_seqs  = [primer_sequences[0][i] for i in final_left]
    final_right_seqs = [primer_sequences[1][i] for i in final_right]
    df = pd.DataFrame({
        "Primer ID":   [f"Primer{i+1}" for i in range(len(final_left_seqs))],
        "Left Primer":  final_left_seqs,
        "Right Primer": final_right_seqs
    })
    df.to_csv(args.output_csv, index=False)
    pbar.set_description(steps[7])
    pbar.update(1)

    pbar.close()
    tqdm.write(f"All steps completed. Optimized primers saved to '{args.output_csv}'")

if __name__ == "__main__":
    main()
