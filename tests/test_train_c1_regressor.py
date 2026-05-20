from pathlib import Path

from scripts.train_c1_regressor import (
    build_training_paths,
    get_feature_labels,
)


def test_build_training_paths_defaults_follow_process_layout():
    repo_root = Path("/repo")
    zh = build_training_paths(repo_root, "zh")
    wh = build_training_paths(repo_root, "wh")

    assert zh["lo_file"] == repo_root / "output" / "zh" / "events_lo.root"
    assert zh["rw_file"] == repo_root / "output" / "zh" / "events_rwgt.root"
    assert zh["outdir"] == repo_root / "output" / "zh" / "c1_regressor"

    assert wh["lo_file"] == repo_root / "output" / "wh" / "events_lo.root"
    assert wh["rw_file"] == repo_root / "output" / "wh" / "events_rwgt.root"
    assert wh["outdir"] == repo_root / "output" / "wh" / "c1_regressor"


def test_build_training_paths_preserves_explicit_overrides():
    repo_root = Path("/repo")
    paths = build_training_paths(
        repo_root,
        "wh",
        lo_file="/tmp/custom_lo.root",
        rw_file="/tmp/custom_rw.root",
        outdir="/tmp/model_out",
    )

    assert paths["lo_file"] == Path("/tmp/custom_lo.root")
    assert paths["rw_file"] == Path("/tmp/custom_rw.root")
    assert paths["outdir"] == Path("/tmp/model_out")


def test_feature_labels_switch_between_zh_and_wh():
    zh = get_feature_labels("zh")
    wh = get_feature_labels("wh")

    assert zh["v_pt"] == r"$p_T(Z)$ [GeV]"
    assert zh["vh_m"] == r"$m(ZH)$ [GeV]"
    assert zh["vh_delta_eta"] == r"$\Delta\eta(Z,H)$"

    assert wh["v_pt"] == r"$p_T(W)$ [GeV]"
    assert wh["vh_m"] == r"$m(WH)$ [GeV]"
    assert wh["vh_delta_eta"] == r"$\Delta\eta(W,H)$"
