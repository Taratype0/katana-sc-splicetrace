# katana-sc-motiftrace (`v2.2.0`)

[![CI](https://github.com/Taratype0/katana-sc-motiftrace/actions/workflows/ci.yml/badge.svg)](https://github.com/Taratype0/katana-sc-motiftrace/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Standalone pure-Python package and CLI for **Smart-seq3 / HT Smart-seq3 / Smart-seq3xpress single-cell RBP splicing program discovery**, **mature-mRNA motif consequence tracing**, **target susceptibility ranking**, and **11-channel support-precision continuous PSI modeling**, bundling the **Curated Splicing & Epitranscriptome Library `v3.0` (`169` RBPs / `621` PWMs)** offline out-of-the-box.

> **Supported sequencing scope**: Full-length single-cell RNA-seq protocols capturing internal splice-junction reads (`Smart-seq3`, `HT Smart-seq3`, `Smart-seq3xpress`). Not intended for 3'/5'-end tag-only droplet data.

---

## Key Capabilities
- **Donor-Stratified Precision Meta-Contrast & `max-T` FWER**: Eliminates pseudoreplication across cells (`0.0%` empirical null FWER on `GSE270928` Smart-seq3 donor permutations).
- **11-Channel Sequence + Baseline PSI Modeling**: Outperforms sequence-only baselines (`-16.1%` held-out MSE vs `CASREL` architecture).
- **Bundled `169`-RBP / `621`-PWM Library `v3.0` + Extensible Plugin Engine**: Includes `YTHDC1` (`UGGACU`), `RBM20`, `RBM10`, `SRRM4`, `U2AF1`, `--species human|mouse`, and zero-code `.meme` / `.motif` / `.jaspar` / `.tsv` / `.json.gz` expansion.

---

## Installation

```bash
# Install directly from GitHub
pip install git+https://github.com/Taratype0/katana-sc-motiftrace.git

# Or install pre-built wheel from GitHub Releases
pip install https://github.com/Taratype0/katana-sc-motiftrace/releases/download/v2.2.0/katana_sc_motiftrace-2.2.0-py3-none-any.whl
```

---

## Standalone CLI Usage

```bash
# Run single-cell Smart-seq3 RBP program discovery & motif consequence tracing
katana-sc-motiftrace \
  --events EVENTS.tsv \
  --psi PSI.tsv.gz \
  --support SUPPORT.tsv.gz \
  --metadata CELL_METADATA.tsv \
  --fasta hg38.fa \
  --species human \
  --numerator Stim \
  --denominator Ctrl \
  --out ./sc_motiftrace_out

# List bundled RBPs (human or mouse symbols) or install an external motif pack
katana-sc-motiftrace --list-rbps --species mouse --filter-rbp Yth
katana-sc-motiftrace --install-pack ./custom_rbps.meme --pack-name lab_pack
```

---

## Katana Standalone Package Suite
| Package | Scope | Repository |
| :--- | :--- | :--- |
| **`katana-motiftrace`** | Bulk RBP regulator discovery & mature-mRNA motif tracing | [Taratype0/katana-motiftrace](https://github.com/Taratype0/katana-motiftrace) |
| **`katana-cosplice`** | Bulk count-debiased ExcessBioVar + Signed TOM co-splicing networks | [Taratype0/katana-cosplice](https://github.com/Taratype0/katana-cosplice) |
| **`katana-sc-motiftrace`** | Smart-seq3 single-cell RBP discovery, motif tracing & PSI modeling | [Taratype0/katana-sc-motiftrace](https://github.com/Taratype0/katana-sc-motiftrace) |
| **`katana-sc-cosplice`** | Smart-seq3 single-cell Hybrid-Manifold co-splicing module discovery | [Taratype0/katana-sc-cosplice](https://github.com/Taratype0/katana-sc-cosplice) |
