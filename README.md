# PRISM - Primer Design through Submodular Function Estimation

## User Guide  
PRISM is a **multiplex‑PCR primer design framework** that formulates primer selection as a constrained submodular‑maximization problem, jointly optimizing genome coverage while reducing primer‑dimerization risk. It employs a provably bounded local‑search algorithm to output primer sets with high coverage and low dimer risk.

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
>  Runtime dependencies (`primer3‑py`, `numpy`, `pandas`, `tqdm`, `numba`, `joblib`) are installed automatically.




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
The input sequence must be provided in reference FASTA format.

---

### Command-line interface
PRISM accepts three main sub-commands, `input` (required), `windows` (optional, with a default size of 250 bp), and `output` (optional, specifies file name and output path):

### CLI options

| Flag | Default | Description |
|------|---------|-------------|
| `-i`, `--input` | **required** | Input FASTA/FA file |
| `-w`, `--window-size` | 250 | Window size for region slicing |
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
# Run PRISM on a Zika reference, 300 bp windows, export CSV
prism \
  -i data/NC_012532.1.fna \
  -w 300 \
  -o output/optimized_primers.csv
```

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

## License

PRISM is released under the **GNU GPL v3.0**.  

See the [LICENSE](LICENSE) file for details.

---

## Citing PRISM
If you use PRISM, please cite:
> Wang A. *et al.* Primer Design through Submodular Function
Estimation. (*Being submitted,coming soon...*)


