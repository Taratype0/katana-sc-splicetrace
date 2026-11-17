"""Unit tests for katana-sc-splicetrace v2.3.0."""
from __future__ import annotations

import numpy as np
import katana_sc_splicetrace as kscst


def test_sc_splicetrace_v23_and_fast_accel() -> None:
    assert kscst.__version__ == "2.3.0"
    lib_h = kscst.load_bundled_human_rbp_library(species="human")
    assert lib_h.n_rbps >= 165
    assert "YTHDC1" in lib_h.rbp_to_models
    scores = kscst.fast_score_mean_log_odds(["UGGACUUGCAUGACGU"], {k: lib_h.models[k] for k in list(lib_h.models)[:10]})
    assert scores.shape == (1, 10)
    assert np.all(np.isfinite(scores))
