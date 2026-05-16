from pathlib import Path

from scripts.analyze import build_process_paths, build_wh_merge_plan


def test_build_process_paths_isolates_public_outputs():
    repo_root = Path("/repo")
    zh = build_process_paths(repo_root, "zh")
    wh = build_process_paths(repo_root, "wh")

    assert zh["outdir"] == repo_root / "output" / "zh"
    assert wh["outdir"] == repo_root / "output" / "wh"
    assert zh["lo_root"].name == "events_lo.root"
    assert wh["rw_root"].name == "events_rwgt.root"


def test_build_wh_merge_plan_points_to_hidden_subchannel_roots():
    repo_root = Path("/repo")
    plan = build_wh_merge_plan(repo_root)
    assert plan["lo_inputs"][0][0] == repo_root / "output" / "_wh_internal" / "wh_plus" / "events_lo.root"
    assert plan["lo_inputs"][1][0] == repo_root / "output" / "_wh_internal" / "wh_minus" / "events_lo.root"
    assert plan["rw_output"] == repo_root / "output" / "wh" / "events_rwgt.root"


from scripts.vh_processes import get_public_process


def test_process_labels_are_available_for_plotting():
    assert get_public_process("zh").process_label == "ZH (13.6 TeV)"
    assert get_public_process("wh").process_label == "WH (13.6 TeV)"
