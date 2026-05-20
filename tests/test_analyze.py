from pathlib import Path
import sys
import types

import scripts.analyze as analyze


def test_build_process_paths_isolates_public_outputs():
    repo_root = Path("/repo")
    zh = analyze.build_process_paths(repo_root, "zh")
    wh = analyze.build_process_paths(repo_root, "wh")

    assert zh["outdir"] == repo_root / "output" / "zh"
    assert wh["outdir"] == repo_root / "output" / "wh"
    assert zh["lo_root"].name == "events_lo.root"
    assert wh["rw_root"].name == "events_rwgt.root"


def test_build_wh_merge_plan_points_to_hidden_subchannel_roots():
    repo_root = Path("/repo")
    plan = analyze.build_wh_merge_plan(repo_root)
    assert plan["lo_inputs"][0][0] == repo_root / "output" / "_wh_internal" / "wh_plus" / "events_lo.root"
    assert plan["lo_inputs"][1][0] == repo_root / "output" / "_wh_internal" / "wh_minus" / "events_lo.root"
    assert plan["rw_output"] == repo_root / "output" / "wh" / "events_rwgt.root"


def test_step_process_lhe_to_root_builds_wh_internal_roots_before_merge(monkeypatch, tmp_path):
    calls = []
    merges = []

    def fake_step(lo_lhe, rw_lhe, lo_root, rw_root, ebeam, subchannel_code=None):
        calls.append((Path(lo_lhe), Path(rw_lhe), Path(lo_root), Path(rw_root), ebeam, subchannel_code))

    def fake_merge(inputs, output):
        merges.append((inputs, output))

    monkeypatch.setattr(analyze, "step_lhe_to_root", fake_step)
    monkeypatch.setattr(analyze, "merge_root_files", fake_merge)

    process_spec = get_public_process("wh")
    process_paths = analyze.build_process_paths(tmp_path, "wh")
    analyze.step_process_lhe_to_root(process_spec, tmp_path, process_paths, ebeam=6800.0)

    assert [call[0] for call in calls] == [
        tmp_path / "output" / "_wh_internal" / "wh_plus" / "events.lhe",
        tmp_path / "output" / "_wh_internal" / "wh_minus" / "events.lhe",
    ]
    assert [call[2] for call in calls] == [
        tmp_path / "output" / "_wh_internal" / "wh_plus" / "events_lo.root",
        tmp_path / "output" / "_wh_internal" / "wh_minus" / "events_lo.root",
    ]
    assert [call[5] for call in calls] == [None, None]
    assert merges[0][0] == [
        (tmp_path / "output" / "_wh_internal" / "wh_plus" / "events_lo.root", 0),
        (tmp_path / "output" / "_wh_internal" / "wh_minus" / "events_lo.root", 1),
    ]
    assert merges[0][1] == tmp_path / "output" / "wh" / "events_lo.root"
    assert merges[1][1] == tmp_path / "output" / "wh" / "events_rwgt.root"


def test_step_plots_passes_process_specific_theory_c1(monkeypatch, tmp_path):
    argv_calls = []

    def fake_call(main_fn, argv):
        argv_calls.append(argv)

    fake_module = types.SimpleNamespace(main=lambda: 0)
    monkeypatch.setattr(analyze, "_call_with_argv", fake_call)
    monkeypatch.setitem(sys.modules, "plot_weight_ratio", fake_module)
    monkeypatch.setitem(sys.modules, "plot_C1_vs_pt", fake_module)
    monkeypatch.setitem(sys.modules, "compare_and_C1", fake_module)
    monkeypatch.setitem(sys.modules, "compare_roots", fake_module)
    monkeypatch.setitem(sys.modules, "plot_kappa3", types.SimpleNamespace(process_and_plot=lambda **kwargs: None))

    plotdir = tmp_path / "plots"
    plotdir.mkdir()

    analyze.step_plots(
        lo_root=str(tmp_path / "events_lo.root"),
        rw_root=str(tmp_path / "events_rwgt.root"),
        plotdir=str(plotdir),
        kappas=[1.0],
        feature="h_pt",
        process_spec=get_public_process("wh"),
    )

    theory_argvs = [argv for argv in argv_calls if "--theory-c1" in argv]
    assert len(theory_argvs) == 3
    for argv in theory_argvs:
        idx = argv.index("--theory-c1")
        assert argv[idx + 1] == "1.03"


from scripts.vh_processes import get_public_process


def test_process_labels_are_available_for_plotting():
    assert get_public_process("zh").process_label == "ZH (13.6 TeV)"
    assert get_public_process("wh").process_label == "WH (13.6 TeV)"
