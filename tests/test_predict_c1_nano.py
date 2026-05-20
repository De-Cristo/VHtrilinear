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
