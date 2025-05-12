from design_primers import PrimerSetBadnessFast
from collections import defaultdict
import numpy as np


def extract_primers(primer_results):
    """
    Flatten primer_results into lists of forward and reverse sequences,
    and generate a paired identifier list indicating strand.
    """
    all_fP = []
    all_rP = []
    pair_id_list = [] 

    for region, primer_pairs in primer_results.items():
        for pair in primer_pairs:
            left_seq = pair['left_primer']['sequence'].lower()
            right_seq = pair['right_primer']['sequence'].lower()
            pair_id = pair['pair_id']  

            all_fP.append(left_seq)
            pair_id_list.append(f"{pair_id}_left")  

            all_rP.append(right_seq)
            pair_id_list.append(f"{pair_id}_right") 

    return all_fP, all_rP, pair_id_list



def compute_badness(primer_results):
    """
    Compute total badness and map each primer identifier to its badness component.

    Calls PrimerSetBadnessFast on the flattened primer lists.
    """
    all_fP, all_rP, pair_id_list = extract_primers(primer_results)

    total_badness, badness_components = PrimerSetBadnessFast(all_fP, all_rP)

    badness_mapping = {}

    for i in range(len(pair_id_list)):
        badness_mapping[pair_id_list[i]] = (
            badness_components[0][i//2] if "_left" in pair_id_list[i] else badness_components[1][i//2]
        )

    return total_badness, badness_mapping



def compute_badness_for_blocks(badness_mapping):
    """
    Aggregate primer badness by block, then compute summary statistics.

    Each block key is derived from the first five underscore-delimited fields
    of the primer identifier.
    """
    
    region_badness_totals = defaultdict(float)

    for pair_id, badness_value in badness_mapping.items():
    
        region = "_".join(pair_id.split("_")[:5])  
      
        region_badness_totals[region] += badness_value
 
    region_badness_totals = dict(region_badness_totals)


    total_badness_sum = sum(region_badness_totals.values())
    region_count = len(region_badness_totals)
    mean_badness = total_badness_sum / region_count

    max_badness_region = max(region_badness_totals, key=region_badness_totals.get)
    max_badness_value = region_badness_totals[max_badness_region]
    diff = max_badness_value - mean_badness

    return region_badness_totals, max_badness_region, max_badness_value, diff


def process_badness_mapping(badness_mapping, optimized_primers):
    """
    Translate badness_map and optimized primer data into numerical arrays
    for downstream optimization routines.

    Matches each '<pair_id>_left' or '_right' entry with its sequence
    in optimized_primers, assigns block indices, and constructs arrays.
    """
    left_primers = []
    right_primers = []
    left_sequences = []
    right_sequences = []
    block_map = {}
    block_indices = []

    for pair_id, badness_value in badness_mapping.items():
       
        block_key = "_".join(pair_id.split("_")[:5])  

        if block_key not in block_map:
            block_map[block_key] = len(block_map)

        block_idx = block_map[block_key]
       
        base_pair_id = pair_id.replace("_left", "").replace("_right", "").strip()
        primer_seq = None

        for primer_data in optimized_primers.get(block_key, []):
            if primer_data["pair_id"].strip() == base_pair_id:
                primer_seq = primer_data
                break 

        if primer_seq is None:
            print(f"No match found for {pair_id} (Base: {base_pair_id}) in optimized_primers[{block_key}]")

        if "_left" in pair_id:
            left_primers.append(badness_value)
            left_sequences.append(primer_seq["left_primer"]["sequence"] if primer_seq else "N/A")
            block_indices.append(block_idx)
        elif "_right" in pair_id:
            right_primers.append(badness_value)
            right_sequences.append(primer_seq["right_primer"]["sequence"] if primer_seq else "N/A")

    assert len(left_primers) == len(right_primers), "Left and right primer counts do not match"
    assert len(left_sequences) == len(right_sequences), "Left and right sequence counts do not match"

    badness_components = np.array([left_primers, right_primers])
    block_indices = np.array(block_indices, dtype=np.int32)

    return badness_components, block_indices, [left_sequences, right_sequences]