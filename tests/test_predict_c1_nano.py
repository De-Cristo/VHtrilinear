from pathlib import Path

import awkward as ak
import numpy as np

from scripts.predict_c1_nano import (
    build_prediction_paths,
    event_has_required_bosons,
)


def test_build_prediction_paths_follow_process_layout():
    repo_root = Path("/repo")
    zh = build_prediction_paths(repo_root, "zh")
    wh = build_prediction_paths(repo_root, "wh")

    assert zh["model"] == repo_root / "output" / "zh" / "c1_regressor" / "c1_regressor.json"
    assert zh["output"] == repo_root / "output" / "zh" / "nano_c1_predictions.root"
    assert zh["plotdir"] == repo_root / "output" / "zh" / "plots" / "nano_validation"

    assert wh["model"] == repo_root / "output" / "wh" / "c1_regressor" / "c1_regressor.json"
    assert wh["output"] == repo_root / "output" / "wh" / "nano_c1_predictions.root"
    assert wh["plotdir"] == repo_root / "output" / "wh" / "plots" / "nano_validation"


def test_event_has_required_bosons_for_zh():
    pdgid = ak.Array([[23, 25], [24, 25], [-24, 25], [23]])
    mask = event_has_required_bosons(pdgid, "zh")
    assert ak.to_list(mask) == [True, False, False, False]


def test_event_has_required_bosons_for_wh():
    pdgid = ak.Array([[23, 25], [24, 25], [-24, 25], [24], [25]])
    mask = event_has_required_bosons(pdgid, "wh")
    assert ak.to_list(mask) == [False, True, True, False, False]


from scripts.predict_c1_nano import get_vector_boson_mask


def test_get_vector_boson_mask_for_zh_selects_only_z():
    pdgid = ak.Array([[23, 25], [24, 25], [-24, 25]])
    mask = get_vector_boson_mask(pdgid, "zh")
    assert ak.to_list(mask) == [[True, False], [False, False], [False, False]]


def test_get_vector_boson_mask_for_wh_selects_both_charges():
    pdgid = ak.Array([[23, 25], [24, 25], [-24, 25]])
    mask = get_vector_boson_mask(pdgid, "wh")
    assert ak.to_list(mask) == [[False, False], [True, False], [True, False]]


import numpy as np

from scripts.predict_c1_nano import compute_kappa_weights


def test_compute_kappa_weights_uses_process_constants():
    c1 = np.asarray([1.0], dtype=np.float64)

    zh = compute_kappa_weights(c1, "zh", [1])
    wh = compute_kappa_weights(c1, "wh", [1])

    assert "weight_kappa_1" in zh
    assert "weight_kappa_abs_1" in wh
    assert zh["weight_kappa_abs_1"].shape == (1,)
    assert wh["weight_kappa_abs_1"].shape == (1,)


from scripts.predict_c1_nano import get_prediction_plot_metadata


def test_prediction_plot_metadata_tracks_process_labels():
    zh = get_prediction_plot_metadata("zh")
    wh = get_prediction_plot_metadata("wh")

    assert zh["process_label"] == "ZH (13.6 TeV)"
    assert zh["vector_label"] == "Z"
    assert zh["dataset_label"] == "ZH"

    assert wh["process_label"] == "WH (13.6 TeV)"
    assert wh["vector_label"] == "W"
    assert wh["dataset_label"] == "WH"
