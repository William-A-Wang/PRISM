"""
PRISM: PCR primer design and optimization pipeline package

This package is distributed under the GNU General Public License v3.0.
See the LICENSE file at the project root for more details.
"""

# Expose key functions at package level
from .data_io import read_consensus_sequence, sliding_window_regions
from .design_primers import design_primers, PrimerSetBadnessFast
from .badness_utils import compute_badness, compute_badness_for_blocks, process_badness_mapping
from .iterative import convert_primer_results_to_regions, iterative_primer_optimization
from .optimization import generate_initial_solution_numba, approximation_algorithm

__all__ = [
    "read_consensus_sequence", "sliding_window_regions",
    "design_primers", "PrimerSetBadnessFast",
    "compute_badness", "compute_badness_for_blocks", "process_badness_mapping",
    "convert_primer_results_to_regions", "iterative_primer_optimization",
    "generate_initial_solution_numba", "approximation_algorithm",
]
