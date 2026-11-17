"""Standalone CLI for katana-sc-splicetrace (2.3.0) with Extensible Motif Database v3.0."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import pandas as pd

from katana_sc_splicetrace.contrast_discovery import (
    compute_sc_splicing_contrast_v2,
    compute_sc_stratified_meta_contrast_v2,
    discover_sc_rbp_programs_v2,
)
from katana_sc_splicetrace.sequence_and_consequence import (
    build_sc_directional_motif_tensor_v2,
    extract_sc_event_sequences_v2,
    trace_sc_mature_motif_consequences_v2,
)
from katana_sc_splicetrace.fast_accel import configure_server_threads
from katana_sc_splicetrace.protein_consequences import trace_protein_and_nmd_consequences
from katana_sc_splicetrace.rna_forms import IndexedFasta
from katana_sc_splicetrace.target_and_psi_models import rank_rbp_target_events_v2

__version__ = "2.3.0"


def write_nfcore_versions_yml(path: Path | str, process_name: str = "KATANA_SCSPLICETRACE") -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(f'"{process_name}":\n    katana-sc-splicetrace: "2.3.0"\n    python: "{platform.python_version()}"\n', encoding="utf-8")


def add_extensible_library_args(parser: argparse.ArgumentParser) -> None:
    grp = parser.add_argument_group("Motif Library & Extensible Database Options")
    grp.add_argument("--species", choices=["human", "mouse"], default="human", help="Target species RBP gene symbol namespace (`human` e.g. YTHDC1 or `mouse` e.g. Ythdc1)")
    grp.add_argument("--extra-motifs", nargs="+", type=Path, default=None, help="One or more external motif files (.meme, .motif, .jaspar, .pfm, .tsv, .csv, .json, .json.gz) to merge on the fly")
    grp.add_argument("--no-v3-expansion", action="store_true", help="Disable v3.0 m6A & ENCODE RBNS/eCLIP expansion (revert to legacy CISBP/ATtRACT 146-RBP core)")
    grp.add_argument("--no-user-packs", action="store_true", help="Ignore user-installed motif packs in ~/.katana/motif_packs")
    grp.add_argument("--list-rbps", action="store_true", help="Print full RBP catalog summary (consensus motifs, mouse symbols, sources) and exit")
    grp.add_argument("--filter-rbp", default=None, help="Optional substring/symbol filter when using --list-rbps (e.g. YTHDC1)")
    grp.add_argument("--import-motifs", type=Path, default=None, help="Convert an external (.meme/.motif/.jaspar/.tsv) file into a portable Katana Motif Pack")
    grp.add_argument("--export-pack", type=Path, default=None, help="Destination .json or .json.gz path for --import-motifs")
    grp.add_argument("--install-pack", type=Path, default=None, help="Install an external motif file or pack into persistent user registry (~/.katana/motif_packs)")
    grp.add_argument("--pack-name", default="custom_user_pack", help="Pack name for --install-pack or --remove-pack")
    grp.add_argument("--remove-pack", default=None, help="Remove an installed user motif pack by name and exit")


def handle_library_management_commands(args: argparse.Namespace) -> tuple[bool, object | None]:
    from katana_sc_splicetrace.library import (
        export_katana_motif_pack,
        get_user_motif_packs_dir,
        load_bundled_human_rbp_library,
        load_curated_human_rbp_library,
        load_external_motif_entries,
    )

    if args.remove_pack:
        pdir = get_user_motif_packs_dir()
        removed = []
        for cand in (pdir / f"{args.remove_pack}.json.gz", pdir / f"{args.remove_pack}.json"):
            if cand.is_file():
                cand.unlink()
                removed.append(str(cand))
        print(json.dumps({"action": "remove_pack", "removed": removed}))
        return True, None

    if args.import_motifs is not None:
        if args.export_pack is None:
            raise ValueError("--export-pack <out.json.gz> is required with --import-motifs")
        ents = load_external_motif_entries(args.import_motifs, default_source=args.pack_name)
        info = export_katana_motif_pack(ents, args.export_pack, pack_name=args.pack_name)
        print(json.dumps({"action": "import_motifs", **info}, indent=2))
        return True, None

    if args.install_pack is not None:
        pdir = get_user_motif_packs_dir()
        pdir.mkdir(parents=True, exist_ok=True)
        dst = pdir / f"{args.pack_name}.json.gz"
        ents = load_external_motif_entries(args.install_pack, default_source=args.pack_name)
        info = export_katana_motif_pack(ents, dst, pack_name=args.pack_name)
        print(json.dumps({"action": "install_pack", **info}, indent=2))
        return True, None

    if getattr(args, "models_json", None) is not None and getattr(args, "assignments_tsv", None) is not None:
        lib = load_curated_human_rbp_library(
            args.models_json,
            args.assignments_tsv,
            include_v3_expansion=not args.no_v3_expansion,
            extra_motif_files=args.extra_motifs,
            include_installed_user_packs=not args.no_user_packs,
            species=args.species,
        )
    else:
        lib = load_bundled_human_rbp_library(
            include_v3_expansion=not args.no_v3_expansion,
            extra_motif_files=args.extra_motifs,
            include_installed_user_packs=not args.no_user_packs,
            species=args.species,
        )

    if args.list_rbps:
        df_sum = lib.to_summary_dataframe()
        if args.filter_rbp:
            q = str(args.filter_rbp).upper()
            df_sum = df_sum[
                df_sum["RBP"].str.upper().str.contains(q)
                | df_sum["mouse_symbol"].str.upper().str.contains(q)
                | df_sum["category"].str.upper().str.contains(q)
            ]
        print(df_sum.to_string(index=False))
        return True, lib

    return False, lib


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="katana-sc-splicetrace",
        description="Katana SC MotifTrace: Smart-seq3 RBP discovery, mature-mRNA motif consequence tracing & extensible motif DB v3.0.",
    )
    p.add_argument("--version", action="version", version=f"katana-sc-splicetrace {__version__}")
    p.add_argument("--events", type=Path, required=False)
    p.add_argument("--psi", type=Path, required=False)
    p.add_argument("--metadata", type=Path, required=False)
    p.add_argument("--fasta", type=Path, required=False)
    p.add_argument("--gtf", type=Path, default=None, help="Optional GTF/GTF.gz for automated single-cell ORF/NMD/domain consequences")
    p.add_argument("--threads", type=int, default=0, help="CPU threads for Numba acceleration (0 = auto)")
    p.add_argument("--support", type=Path, default=None)
    p.add_argument("--transpose-psi", action="store_true")
    p.add_argument("--group-col", default="condition")
    p.add_argument("--numerator", default=None)
    p.add_argument("--denominator", default=None)
    p.add_argument("--donor-col", default="donor_id")
    p.add_argument("--query-rbp", default=None)
    p.add_argument("--models-json", type=Path, default=None)
    p.add_argument("--assignments-tsv", type=Path, default=None)
    p.add_argument("--permutations", type=int, default=1999)
    p.add_argument("--seed", type=int, default=20261117)
    p.add_argument("--out", "--output-dir", dest="out", type=Path, required=False)
    p.add_argument("--emit-versions-yml", type=Path, default=None)
    p.add_argument("--process-name", default="KATANA_SCSPLICETRACE")
    add_extensible_library_args(p)
    args = p.parse_args(argv)

    handled, lib = handle_library_management_commands(args)
    if handled:
        return 0

    if args.emit_versions_yml is not None and args.out is None:
        write_nfcore_versions_yml(args.emit_versions_yml, args.process_name)
        return 0
    if args.events is None or args.fasta is None or args.out is None:
        p.error("--events, --fasta, and --out (or --list-rbps / --import-motifs / --install-pack) are required.")

    configure_server_threads(args.threads)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    events_df = pd.read_csv(args.events, sep="\t")
    rna_df, cov_df = extract_sc_event_sequences_v2(events_df, Path(args.fasta))
    rbp_tensor, ev_ids, rbps, region_names = build_sc_directional_motif_tensor_v2(events_df, rna_df, lib)

    summary: dict = {"version": __version__, "library": lib.name, "n_rbps": lib.n_rbps, "n_models": lib.n_models, "status": "complete"}
    if args.query_rbp is not None:
        resolved_q = lib.resolve_rbp_symbol(args.query_rbp) or args.query_rbp
        t_df = rank_rbp_target_events_v2(resolved_q, ev_ids, events_df["gene_id"].astype(str).to_numpy(), rbp_tensor, rbps, region_names)
        t_df.to_csv(out / f"sc_motiftrace_targets_{resolved_q}.tsv", sep="\t", index=False)
        summary["ranked_targets_for"] = resolved_q

    if args.psi is not None and args.metadata is not None and args.numerator is not None and args.denominator is not None:
        psi_df = pd.read_csv(args.psi, sep="\t", index_col=0)
        sup_df = pd.read_csv(args.support, sep="\t", index_col=0) if args.support else None
        if args.transpose_psi:
            psi_df = psi_df.T
            if sup_df is not None:
                sup_df = sup_df.T
        md = pd.read_csv(args.metadata, sep="\t")
        eg = events_df.set_index("event_id")["gene_id"]
        u_aligned = psi_df.reindex(columns=events_df["event_id"].astype(str))
        s_aligned = sup_df.reindex(columns=events_df["event_id"].astype(str)) if sup_df is not None else None

        d_col = args.donor_col if args.donor_col in md.columns else ("donor" if "donor" in md.columns else None)
        n_donors = md[d_col].nunique() if d_col else 1
        if 2 <= n_donors < 5:
            contrast_df, _, c_info = compute_sc_stratified_meta_contrast_v2(
                u_aligned, md, eg, support=s_aligned, group_col=args.group_col,
                numerator_group=args.numerator, denominator_group=args.denominator, donor_col=d_col, seed=args.seed,
            )
            donor_deltas = None
        else:
            contrast_df, donor_deltas, c_info = compute_sc_splicing_contrast_v2(
                u_aligned, md, eg, support=s_aligned, group_col=args.group_col,
                numerator_group=args.numerator, denominator_group=args.denominator, donor_col=args.donor_col, seed=args.seed,
            )

        conseq_df = trace_sc_mature_motif_consequences_v2(rna_df, contrast_df, lib)
        res = discover_sc_rbp_programs_v2(
            contrast_df, rbp_tensor, rbps, region_names, cov_df,
            donor_delta_df=donor_deltas, mature_consequences_df=conseq_df,
            n_permutations=args.permutations, seed=args.seed,
        )
        res.rbp_ranking.to_csv(out / "sc_motiftrace_rbp_ranking_v21.tsv", sep="\t", index=False)
        res.region_tests.to_csv(out / "sc_motiftrace_region_tests_v21.tsv", sep="\t", index=False)
        res.event_contrasts.to_csv(out / "sc_motiftrace_event_contrasts_v21.tsv", sep="\t", index=False)
        res.mature_consequences.to_csv(out / "sc_motiftrace_mature_consequences_v21.tsv", sep="\t", index=False)
        summary.update({**c_info, **res.diagnostics})

    (out / "sc_motiftrace_summary_v21.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if args.emit_versions_yml is not None:
        write_nfcore_versions_yml(args.emit_versions_yml, args.process_name)
    print(json.dumps(summary))
        if args.gtf is not None and args.events is not None and args.fasta is not None:
        ev_in = pd.read_csv(args.events, sep="\t")
        p_df = trace_protein_and_nmd_consequences(ev_in.head(500), args.gtf, args.fasta, indexed_fasta_cls=IndexedFasta)
        p_df.to_csv(out / "sc_protein_and_nmd_consequences_v23.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
