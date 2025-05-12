from typing import List, Tuple, Set, Dict


# data process
def read_consensus_sequence(file_path):

    """Read FASTA"""

    with open(file_path, 'r') as file:
        lines = file.readlines()[1:]

        consensus_sequence = "".join([line.strip() for line in lines])

    return consensus_sequence

def sliding_window_regions(sequence: str, window_sizes: List[int], overlaps: int = 250) -> Dict[int, List[Dict]]:

    """ 
    Generate a list of intervals with or without overlaps for sequences with different window_sizes.
    Each region_dict contains id, sequence, start, end, length.
    """

    all_regions = {}
    sequence_length = len(sequence)

    for window_size in window_sizes:
        step_size = window_size - overlaps
        start = 0
        regions = []

        while start < sequence_length:
            end = min(start + window_size, sequence_length)

            is_last_window = (end == sequence_length)

            region = {
                'id': f'size_{window_size}_region_{start}_{end}',
                'sequence': sequence[start:end],
                'start': start,
                'end': end,
                'length': end - start
            }
            regions.append(region)


            if is_last_window:
                break

            start += step_size

        all_regions[window_size] = regions

    return all_regions