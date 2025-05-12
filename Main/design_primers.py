from typing import List, Dict
from collections import Counter
import primer3
import numpy as np

def design_primers(regions: dict) -> Dict[str, List[dict]]:
    """
    Design PCR primer pairs for each region using Primer3.

    Parameters
    ----------
    regions : dict
        Mapping from window size to a list of region dictionaries.
        Each region dict must contain:
          - 'id': unique identifier for the region
          - 'sequence': DNA sequence string for primer design

    Returns
    -------
    primer_results : dict
        Mapping from region ID to a list of primer-pair dictionaries.
        Each primer-pair dict contains:
          - 'pair_id': unique identifier for this primer pair
          - 'left_primer': dict with sequence, tm, gc_percent, position, length
          - 'right_primer': dict with sequence, tm, gc_percent, position, length
          - 'product_size': expected amplicon length
    """

    primer_results = {}

    global_primer3_settings = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_NUM_RETURN': 5,
        'PRIMER_MIN_SIZE': 20,
        'PRIMER_MAX_SIZE': 26,
        'PRIMER_MIN_TM': 52.0,
        'PRIMER_MAX_TM': 62.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0
    }

    for window_size, region_list in regions.items():
        for region in region_list:
            seq_args = {
                'SEQUENCE_ID': region['id'],
                'SEQUENCE_TEMPLATE': region['sequence'],
                'SEQUENCE_INCLUDED_REGION': [0, len(region['sequence'])]
            }

            try:
                result = primer3.bindings.designPrimers(seq_args, global_primer3_settings)
                primer_pairs = []

                num_pairs = result.get('PRIMER_PAIR_NUM_RETURNED', 0)

                for i in range(num_pairs):
                    primer_info = {
                        'pair_id': f"{region['id']}_pair_{i}",
                        'left_primer': {
                            'sequence': result[f'PRIMER_LEFT_{i}_SEQUENCE'],
                            'tm': result[f'PRIMER_LEFT_{i}_TM'],
                            'gc_percent': result[f'PRIMER_LEFT_{i}_GC_PERCENT'],
                            'position': result[f'PRIMER_LEFT_{i}'][0],
                            'length': result[f'PRIMER_LEFT_{i}'][1]
                        },
                        'right_primer': {
                            'sequence': result[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                            'tm': result[f'PRIMER_RIGHT_{i}_TM'],
                            'gc_percent': result[f'PRIMER_RIGHT_{i}_GC_PERCENT'],
                            'position': result[f'PRIMER_RIGHT_{i}'][0],
                            'length': result[f'PRIMER_RIGHT_{i}'][1]
                        },
                        'product_size': result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
                    }
                    primer_pairs.append(primer_info)

                primer_results[region['id']] = primer_pairs

            except Exception as e:
                pass

    return primer_results

def revcomp(seq):
    """
    Portions of this module are adapted from the Olivar project:
       https://github.com/treangenlab/Olivar
     Olivar is licensed under the GNU GPL-3.0:
       https://github.com/treangenlab/Olivar/blob/main/LICENSE
    
     This file is therefore distributed under the terms of the GNU General Public
     License version 3.0.  You should have received a copy of the license along
     with this program.  If not, see <https://www.gnu.org/licenses/>.    


    """

    '''
    give the reverse complement of input sequence
    base & number conversion:
        {'A':0, 'T':1, 'C':2, 'G':3}
    input:
        string or array sequence
    output:
        reverse complement of input
    '''
    if isinstance(seq, str):
        complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N',
        'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n',
        '0':'1', '1':'0', '2':'3', '3':'2', '4':'4'}
        try:
            bases = [complement[base] for base in seq]
        except KeyError:
            raise ValueError(f"Base(s) other than 'A', 'T', 'C', 'G' is found in '{seq[:10]}...', ambiguous bases are not accepted.")
        bases.reverse()
        return ''.join(bases)
    elif isinstance(seq, list):
        complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N',
        'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n',
        '0':'1', '1':'0', '2':'3', '3':'2', '4':'4',
        0:1, 1:0, 2:3, 3:2, 4:4}
        try:
            bases = [complement[base] for base in seq]
        except KeyError:
            raise ValueError(f"Base(s) other than 'A', 'T', 'C', 'G' is found in '{seq[:10]}...', ambiguous bases are not accepted.")
        bases.reverse()
        return bases
    else:
        raise ValueError('Only string or list is accepted for reverse complement.')


GC_LETTER = ['C', 'G', 'c', 'g']

def PrimerSetBadnessFast(all_fP: list, all_rP=[], existing=[]):
    '''
    Input:
        all_fP: list of strings
        all_rP: list of strings
        existing: list of strings
    Output:
        Badness: total Loss of the primer set (all fP and rP)
        Badness_component: Badness of each primer
    '''
    endhash4 = Counter()
    endhash5 = Counter()
    endhash6 = Counter()
    middlehash7 = Counter()
    middlehash8 = Counter()

    PENALTY_OFFSET = 0 # if=1, then d=1 gives 2/3 socre, if =0, then d=1 gives 1/2 score, if =-0.5, then d=1 gives 1/3 score
    """
    In this section we made changes to make it fit our badness calculation
    """
    END4 = 16
    END5 = 32
    END6 = 64
    MIDDLE7 = 128
    MIDDLE8 = 256

    # set all sequences as lowercase
    all_fP = [p.lower() for p in all_fP]
    all_rP = [p.lower() for p in all_rP]
    existing = [p.lower() for p in existing]

    # Set up end hash tables (check end 4 nt to middle binding)
    for p in all_fP+all_rP+existing:
        endhash4[p[-4:]] += 1

        endhash5[p[-5:]] += 1
        endhash6[p[-6:]] += 1


    # Set up middle hash table
    # Middlehash penalizes based on closeness to 3' end, Badness = 2 / (distance to 3' + 2);
    # So absolute last 7 is worth 1, 1nt away is worth 0.67, 2nt away is worth 0.5, etc.
    for p in all_fP+all_rP+existing:
        l = len(p)
        for j in range(l-6):
            middlehash7[p[j:j+7]] += (PENALTY_OFFSET + 1) / (l - j - 6 + PENALTY_OFFSET)
        for j in range(l-7):
            middlehash8[p[j:j+8]] += (PENALTY_OFFSET + 1) / (l - j - 7 + PENALTY_OFFSET)

    # Run through each sequence's reverse complement to add up badness
    Badness = 0
    Badness_component = []
    for one_side_primer in [all_fP, all_rP]:
        one_side_badness = []
        for p in one_side_primer:
            p_badness = 0
            c = revcomp(p)
            l = len(c)
            for j in range(l-3):
                k = c[j:j+4]
                try:
                    endscore4 = endhash4[k]
                    numGC = len([b for b in k if b in GC_LETTER])
                    p_badness += (endscore4 * END4 * (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET)) * (2**numGC)

                except KeyError:
                    pass
            for j in range(l-4):
                k = c[j:j+5]
                try:
                    endscore5 = endhash5[k]
                    numGC = len([b for b in k if b in GC_LETTER])
                    p_badness += (endscore5 * END5 * (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET)) * (2**numGC)

                except KeyError:
                    pass
            for j in range(l-5):
                k = c[j:j+6]
                try:
                    endscore6 = endhash6[k]
                    numGC = len([b for b in k if b in GC_LETTER])
                    p_badness += (endscore6 * END6 * (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET)) * (2**numGC)

                except KeyError:
                    pass
            for j in range(l-6):
                k = c[j:j+7]
                try:
                    midscore7 = middlehash7[k]
                    numGC = len([b for b in k if b in GC_LETTER])
                    p_badness += (midscore7 * MIDDLE7 * (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET)) * (2**numGC)

                except KeyError:
                    pass
            for j in range(l-7):
                k = c[j:j+8]
                try:
                    midscore8 = middlehash8[k]
                    numGC = len([b for b in k if b in GC_LETTER])
                    p_badness += (midscore8 * MIDDLE8 * (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET)) * (2**numGC)

                except KeyError:
                    pass
            Badness += p_badness
            one_side_badness.append(p_badness)
        Badness_component.append(one_side_badness)

    return Badness, Badness_component