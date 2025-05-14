# PRISM - Primer Design through Submodular Function Estimation

## User Guide  
**PRISM** (PRImer Selection through Submodular Maximization) is an open-source software tool for automated design of multiplex PCR primers, developed to support viral genome sequencing, pathogen surveillance, and other high-throughput molecular biology applications. Released under the **GNU General Public License (GPL)**, PRISM is freely available and encourages community use, extension, and collaboration.

PRISM introduces a principled optimization framework that formulates the primer design task as a **constrained submodular maximization problem**, balancing two competing objectives: maximizing genome coverage and minimizing undesired primer-primer interactions, quantified using the **Badness value**. This formulation allows PRISM to leverage a fast local search algorithm with a **provable constant-factor approximation guarantee**, making it the first primer design method to combine rigorous theoretical guarantees with practical scalability.

In extensive evaluations on viral genome datasets such as foot-and-mouth disease virus (FMDV) and Zika virus, PRISM consistently outperforms leading tools—including **PrimalScheme**, **Olivar**, and **primerJinn**—achieving significantly lower Badness values, tighter distributions of primer quality, and robust genome coverage with low memory usage and runtime.

<p align="center">
  <img src="https://img.shields.io/pypi/v/prism-bio.svg?color=blue" alt="PyPI">
  <img src="https://img.shields.io/github/license/William-A-Wang/PRISM.svg" alt="License">
</p>


---

## Installation


> **Requirements**
> 
>  – Python ≥ 3.9.  
>  – Linux, macOS or Windows
> 
>  Runtime dependencies (`primer3‑py`, `numpy`, `pandas`, `tqdm`, `numba`, `joblib`) will be installed automatically.




### Installation from PyPI using pip
Recommend create and activate a virtual environment:
```bash
python3 -m venv /path/to/prism
source /path/to/prism/bin/activate
```
Or using conda:
```bash
conda create -n prism python=3.9 -y
conda activate prism
```

Please use the following command to install：

```bash
pip install prism-bio
```


Once the installation is complete, use the following commands to check PRISM's version:
```bash
prism --version
```
---

## General usage
### Input files
PRISM requires a **reference genome in FASTA format** as input. This file should contain one or more nucleotide sequences in standard FASTA format. For optimal results in multiplex primer design, the input should represent a high-quality reference or consensus sequence derived from your target population or viral strain.

---

### Command-line interface
PRISM accepts three main sub-commands, `input` (required), `block` (optional, with a default size of 250 bp), `extend` (optional, with a default size of 100 bp) and `output` (optional, specifies file name and output path):

### CLI options

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input` | **required** | Input FASTA/FA file |
| `-b`, `--block-size` | 250 bp | Block size for region slicing |
| `-e`, `--extend-block-size` | 100 bp | Block extension size during optimisation |
| `-o`, `--output-csv` | `optimized_primers.csv` | Output file (CSV) |

Full help:

```bash
prism --help
or
prism -h
```

---

### Example usage

```bash
# Run PRISM on a Zika reference, 300 bp blocks, 100 bp block extend size, export CSV
prism \
  -i data/NC_012532.1.fna \
  -b 250 \
  -e 100 \
  -o output/optimized_primers.csv
```

You can find this example file `NC_012532.1.fna` in `example input`. This is the Zika virus sequence from the [NCBI dataset](https://www.ncbi.nlm.nih.gov/nuccore/NC_012532.1?report=fasta).


---

## Generate files

The process ends with a .csv file in the form of the following table:

  | Primer ID | Left Primer | Right Primer |
  |-----------|------------|--------------|
  | Primer 1  | GTGTGA…    | CGTAGC… |
  | Primer 2  | GCGTAC…    | TAGCCA… |
  | Primer…   | ………………    | ……………… |

This file will provide the designed primer serial numbers as well as the sequence

---


## Developer Setup

```bash
git clone https://github.com/William-A-Wang/PRISM.git
python -m venv .venv && . .venv/Scripts/activate
pip install -e .[dev]
```
Using this approach, you can deploy code locally while performing development-related activities


---

## External resources

We used Olivar's badness calculation function in our code, but modified some of the parameters.

* [Olivar](https://github.com/treangenlab/Olivar)

We also used primer3 as the initial primer generation.
* [Primer3](https://pypi.org/project/primer3-py/)
---

## License

PRISM is released under the **GNU GPL v3.0**.  

See the [LICENSE](LICENSE) file for details.

---

## Citing PRISM
If you use PRISM in an academic setting, please cite:

    @misc{WangPRISM25,
    author      = { Ao Wang and 
                    Yixin Chen and
                    Aaron Hong and
                    Adam Rivers and
                    Alan Kuhnle and
                    Christina Boucher},
    title       = {Primer Design through Submodular Function Estimation},
    note        = {In submission}
    year        = {2025},
    }


