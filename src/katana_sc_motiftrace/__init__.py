"""Katana Single-Cell MotifTrace 2.2.0 (Smart-seq3) with Extensible Library v3.0."""
from katana_sc_motiftrace.contrast_discovery import (
    SCMotifContrastResultV2,
    compute_sc_splicing_contrast_v2,
    compute_sc_stratified_meta_contrast_v2,
    discover_sc_rbp_programs_v2,
)
from katana_sc_motiftrace.extensible_db import (
    export_katana_motif_pack,
    iupac_to_pwm,
    load_external_motif_entries,
    merge_motif_entries_into_library,
)
from katana_sc_motiftrace.library import CuratedMotifLibrary, load_bundled_human_rbp_library
from katana_sc_motiftrace.sequence_and_consequence import (
    build_sc_directional_motif_tensor_v2,
    extract_sc_event_sequences_v2,
    load_default_curated_library_v2,
    trace_sc_mature_motif_consequences_v2,
)
from katana_sc_motiftrace.target_and_psi_models import (
    estimate_sc_psi_eb_baseline_v2,
    gene_equal_psi_mse_v2,
    predict_sc_psi_precision_v2,
    rank_rbp_target_events_v2,
)

__version__ = "2.2.0"

__all__ = [
    "__version__",
    "CuratedMotifLibrary",
    "SCMotifContrastResultV2",
    "load_bundled_human_rbp_library",
    "load_default_curated_library_v2",
    "load_external_motif_entries",
    "merge_motif_entries_into_library",
    "export_katana_motif_pack",
    "iupac_to_pwm",
    "extract_sc_event_sequences_v2",
    "build_sc_directional_motif_tensor_v2",
    "trace_sc_mature_motif_consequences_v2",
    "compute_sc_splicing_contrast_v2",
    "compute_sc_stratified_meta_contrast_v2",
    "discover_sc_rbp_programs_v2",
    "rank_rbp_target_events_v2",
    "estimate_sc_psi_eb_baseline_v2",
    "predict_sc_psi_precision_v2",
    "gene_equal_psi_mse_v2",
]
