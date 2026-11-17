"""Unit tests for katana-sc-splicetrace v2.4.0 (exports, gene_name UniProt mapping, GTF/NMD, & CLI)."""
from __future__ import annotations

from pathlib import Path
import pandas as pd
import katana_sc_splicetrace as ksc
from katana_sc_splicetrace.cli import main as sc_cli_main


def test_sc_exports_uniprot_and_cli() -> None:
    assert ksc.__version__ == "2.4.0"
    for sym in ksc.__all__:
        assert hasattr(ksc, sym), f"Missing export {sym}"
    lib_h = ksc.load_bundled_human_rbp_library(species="human")
    assert lib_h.n_rbps >= 168
    assert "YTHDC1" in lib_h.rbp_to_models
    dom_h = ksc.load_bundled_uniprot_domains("human")
    assert len(dom_h) >= 16000
    assert "ptprc" in dom_h and "pdcd4" in dom_h
    assert sc_cli_main(["--list-rbps", "--filter-rbp", "CELF2"]) == 0


def test_sc_gene_name_column_uniprot_resolution(tmp_path: Path) -> None:
    # Verify single-cell EVENTS.tsv with `gene_name` (no `gene_symbol`!) resolves UniProt domains cleanly!
    ex1 = "ATG" + "GCT" * 40
    alt = "GCT" * 18
    ex3 = "GCT" * 40 + "TGA"
    seq = "N" * 80 + ex1 + "GTA" + "T" * 77 + alt + "GTA" + "T" * 77 + ex3 + "N" * 80
    fa = tmp_path / "sc_syn.fa"
    fa.write_text(f">chr1\n{seq}\n", encoding="utf-8", newline="\n")
    (tmp_path / "sc_syn.fa.fai").write_text(f"chr1\t{len(seq)}\t6\t{len(seq)}\t{len(seq)+1}\n", encoding="utf-8", newline="\n")

    e1_0, e1_1 = 80, 80 + len(ex1)
    al_0, al_1 = e1_1 + 80, e1_1 + 80 + len(alt)
    e3_0, e3_1 = al_1 + 80, al_1 + 80 + len(ex3)
    gtf = tmp_path / "sc_syn.gtf"
    gtf.write_text(
        f'chr1\ttest\texon\t{e1_0+1}\t{e1_1}\t.\t+\t.\tgene_id "ENSG00000150593"; transcript_id "ENST_P1"; gene_name "PDCD4"; transcript_type "protein_coding";\n'
        f'chr1\ttest\tCDS\t{e1_0+1}\t{e1_1}\t.\t+\t0\tgene_id "ENSG00000150593"; transcript_id "ENST_P1"; gene_name "PDCD4"; transcript_type "protein_coding";\n'
        f'chr1\ttest\texon\t{e3_0+1}\t{e3_1}\t.\t+\t.\tgene_id "ENSG00000150593"; transcript_id "ENST_P1"; gene_name "PDCD4"; transcript_type "protein_coding";\n'
        f'chr1\ttest\tCDS\t{e3_0+1}\t{e3_1-3}\t.\t+\t0\tgene_id "ENSG00000150593"; transcript_id "ENST_P1"; gene_name "PDCD4"; transcript_type "protein_coding";\n',
        encoding="utf-8",
    )
    ev = pd.DataFrame([{
        "event_id": f"SS3SE:chr1:{e1_1}:{al_0}:{al_1}:{e3_0}:+:ENSG00000150593",
        "gene_id": "ENSG00000150593",
        "gene_name": "PDCD4",
        "event_type": "SE",
        "chr": "chr1",
        "strand": "+",
        "exonStart_0base": al_0,
        "exonEnd": al_1,
        "upstreamES": e1_0,
        "upstreamEE": e1_1,
        "downstreamES": e3_0,
        "downstreamEE": e3_1,
        "dpsi": -0.22,
    }])
    p_df = ksc.trace_protein_and_nmd_consequences(ev, gtf, fa, indexed_fasta_cls=ksc.IndexedFasta, species="human")
    assert len(p_df) == 1
    assert p_df.iloc[0]["gene_symbol"] == "PDCD4"
    assert p_df.iloc[0]["uniprot_acc"] != "NA" and p_df.iloc[0]["uniprot_acc"] != "predicted"
