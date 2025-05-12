from typing import List, Dict, Tuple
import numpy as np
from design_primers import design_primers
from badness_utils import compute_badness, compute_badness_for_blocks

def adjust_blocks(primer_results: Dict[str, List[Dict]], block_badness: Dict[str, float], sequence: str,
                  extend_size=100, k=2.0, min_block_size=200):
    """
    Expand the region with the highest badness if it exceeds the threshold
    defined by mean + k * std, then shift or merge adjacent blocks as needed.
    """

    sequence_length = len(sequence)  
    # compute threshold
    mean_badness = np.mean(list(block_badness.values()))
    std_badness = np.std(list(block_badness.values()))
    threshold = mean_badness + std_badness * k

    max_badness_block = max(block_badness, key=block_badness.get)
    max_badness_value = block_badness[max_badness_block]
    diff = max_badness_value - mean_badness

    if max_badness_value > threshold:
        pass
    else:
        return primer_results
    

    # parse start/end of the block to expand
    parts = max_badness_block.split("_")
    target_start, target_end = int(parts[-2]), int(parts[-1])

    # compute new end position
    new_end = min(target_end + extend_size, sequence_length)
    new_target_block = f"size_1000_region_{target_start}_{new_end}"

    # sort blocks by their start coordinate
    block_ids = sorted(primer_results.keys(), key=lambda x: int(x.split("_")[-2]))

    updated_primer_results = {}

    for i, block_id in enumerate(block_ids):
        start, end = map(int, block_id.split("_")[-2:])

        if block_id == max_badness_block:
            # replace worst block with its expanded version
            updated_primer_results[new_target_block] = primer_results[block_id]

        elif start > target_start: 
            # shift downstream blocks forward
            new_start = start + extend_size
            new_end = min(end + extend_size, sequence_length) 
            new_block_id = f"size_1000_region_{new_start}_{new_end}"
            updated_primer_results[new_block_id] = primer_results[block_id]

        else:
            # keep upstream blocks unchanged
            updated_primer_results[block_id] = primer_results[block_id]

        # after each insertion, check if the last block needs merging
        sorted_blocks = sorted(updated_primer_results.keys(), key=lambda x: int(x.split("_")[-2]))
        last_block = sorted_blocks[-1]

        if last_block == max_badness_block:

            previous_block = sorted_blocks[-2]
            prev_start, prev_end = map(int, previous_block.split("_")[-2:])
            last_start, last_end = map(int, last_block.split("_")[-2:])

            new_merged_block = f"size_1000_region_{prev_start}_{last_end}"
            
            updated_primer_results[new_merged_block] = updated_primer_results.get(previous_block, []) + updated_primer_results.get(last_block, [])
            del updated_primer_results[last_block]


    # ensure the final block is not too small
    sorted_blocks = sorted(updated_primer_results.keys(), key=lambda x: int(x.split("_")[-2]))
    last_block = sorted_blocks[-1]
    last_start, last_end = map(int, last_block.split("_")[-2:])
    last_block_size = last_end - last_start

    if last_block_size < min_block_size:

        previous_block = sorted_blocks[-2]
        prev_start, prev_end = map(int, previous_block.split("_")[-2:])

        new_merged_block = f"size_1000_region_{prev_start}_{last_end}"
        updated_primer_results[new_merged_block] = updated_primer_results.get(previous_block, []) + updated_primer_results.get(last_block, [])

        if last_block in updated_primer_results:
            del updated_primer_results[last_block]
        
    
    sorted_primer_results = dict(sorted(updated_primer_results.items(), key=lambda x: int(x[0].split("_")[-2])))

    return sorted_primer_results  

def convert_primer_results_to_regions(primer_results: Dict[str, List[Dict]], sequence: str) -> Dict[int, List[Dict]]:
    """
    Convert the primer_results mapping into a regions dict;
    each key is window size, and values are dicts with 'id' and 'sequence'.
    """
    regions = {}

    for block_id in primer_results.keys():
        parts = block_id.split("_")

        try:
            start, end = int(parts[-2]), int(parts[-1]) 
        except ValueError:
            # skip any improperly formatted block IDs
            continue  

        # compute window size
        window_size = end - start 
        region_data = {
            "id": block_id,
            "sequence": sequence[start:end] 
        }

        if window_size not in regions:
            regions[window_size] = []
        regions[window_size].append(region_data)

    return regions

def iterative_primer_optimization(sequence: str, regions: Dict[int, List[Dict]], max_iterations=50, extend_size=100, k=2.0):
    """
    Perform iterative primer design:
    1. design primers for current regions
    2. compute badness and expand worst block if needed
    3. repeat until convergence or max_iterations reached
    """
    
    iteration = 1
    primer_results = design_primers(regions) 

    while iteration <= max_iterations:

        # compute total badness and per-pair mapping
        total_badness, badness_mapping = compute_badness(primer_results)

        # aggregate scores per block
        region_badness_totals, max_badness_region, max_badness_value, diff = compute_badness_for_blocks(badness_mapping)

        block_badness = region_badness_totals

        # attempt to expand the worst block
        new_primer_results = adjust_blocks(primer_results, block_badness, sequence, extend_size, k=k)

        if new_primer_results == primer_results:
            
            break

        # prepare for next round
        new_regions = convert_primer_results_to_regions(new_primer_results, sequence)
        primer_results = design_primers(new_regions)

        regions = new_regions
        iteration += 1

    # ensure final ordering by start position
    primer_results = dict(sorted(primer_results.items(), key=lambda x: int(x[0].split("_")[-2])))

    return primer_results