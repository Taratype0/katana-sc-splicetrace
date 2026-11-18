# katana-sc-splicetrace (`v2.4.1`) — KatanaAS Suite

[![CI](https://github.com/Taratype0/katana-sc-splicetrace/actions/workflows/ci.yml/badge.svg)](https://github.com/Taratype0/katana-sc-splicetrace/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Part of the **KatanaAS** 4-package Bulk & Single-Cell (Smart-seq3) alternative splicing suite. Supports **two first-class invocation modes**:
1. **Standalone Python Package & CLI** (`pip install` -> Python API / `katana-sc-splicetrace` / `sc-splicetrace` / `python -m katana_sc_splicetrace`)
2. **Nextflow `nf-core` DSL2 Module & Subworkflow** (`process KATANA_SCSPLICETRACE` with native `--emit-versions-yml`)

---

## Mode 1: Standalone Package Installation & Invocation

### Installation
```bash
# Direct one-line wheel install from GitHub Release
pip install https://github.com/Taratype0/katana-sc-splicetrace/releases/download/v2.4.1/katana_sc_splicetrace-2.4.1-py3-none-any.whl

# Or install from source with server Numba acceleration & Matplotlib visualization extras
pip install "git+https://github.com/Taratype0/katana-sc-splicetrace.git#egg=katana-sc-splicetrace[server,viz]"
```

### CLI Usage (`katana-sc-splicetrace`, `sc-splicetrace`, or `python -m katana_sc_splicetrace`)
```bash
katana-sc-splicetrace \
  --events EVENTS.tsv \
  --psi PSI_CELL_BY_EVENT.tsv.gz \
  --support SUPPORT_CELL_BY_EVENT.tsv.gz \
  --metadata CELL_METADATA.tsv \
  --fasta GRCh38.fa \
  --gtf gencode.v46.annotation.gtf.gz \
  --splicenet-modules ./sc_splicenet_out/sc_cosplice_modules_v21.tsv \
  --species human --numerator Stim --denominator Ctrl --threads 16 \
  --output-dir ./sc_splicetrace_out
```

### Python API Usage
```python
import katana_sc_splicetrace as kscst

kscst.configure_server_threads(16)
lib = kscst.load_bundled_human_rbp_library(species='human')
rna_df, cov_df = kscst.extract_sc_event_sequences_v2(events_df, 'GRCh38.fa')
tensor, ev_ids, rbps, regions = kscst.build_sc_directional_motif_tensor_v2(events_df, rna_df, lib)
meta_c, _, _ = kscst.compute_sc_stratified_meta_contrast_v2(
    psi_df, cell_meta_df, events_df.set_index('event_id')['gene_id'],
    support=sup_df, numerator_group='Stim', denominator_group='Ctrl'
)
res = kscst.discover_sc_rbp_programs_v2(meta_c, tensor, rbps, regions, cov_df)
```

---

## Mode 2: Nextflow `nf-core` DSL2 Integration

Every KatanaAS CLI natively implements `--emit-versions-yml versions.yml --process-name "${task.process}"` and SLURM multi-core `--threads ${task.cpus}`.

```groovy
include { KATANA_SCSPLICETRACE } from './modules/local/katana_scsplicetrace/main'
```

---

## KatanaAS 4-Package Suite Architecture

| Package | Version | Modality | Responsibility |
| :--- | :---: | :---: | :--- |
| **[katana-splicetrace](https://github.com/Taratype0/katana-splicetrace)** | `v2.4.3` | Bulk RNA-seq | RBP discovery (`169` RBPs / `621` PWMs), mature-mRNA motifs, `SpliceNet` bridge & `GTF/UniProt` ORF/NMD + 5 plots |
| **[katana-splicenet](https://github.com/Taratype0/katana-splicenet)** | `v2.3.2` | Bulk RNA-seq | Pure-Python `ExcessBioVar` + Adaptive Scale-Free Signed TOM co-splicing network & GO enrichment |
| **[katana-sc-splicetrace](https://github.com/Taratype0/katana-sc-splicetrace)** | `v2.4.1` | Smart-seq3 SC | Single-cell RBP discovery, mature-mRNA motifs, `SC-SpliceNet` bridge & `GTF/UniProt` ORF/NMD + 5 plots |
| **[katana-sc-splicenet](https://github.com/Taratype0/katana-sc-splicenet)** | `v2.2.2` | Smart-seq3 SC | Hybrid-Manifold `ExcessBioVar` + Empirical-Bayes Fisher-z shrunk Signed TOM single-cell modules |
