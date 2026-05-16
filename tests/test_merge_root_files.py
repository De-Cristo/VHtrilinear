from pathlib import Path

import numpy as np
import uproot

from scripts.merge_root_files import merge_root_files


def _write_root(path: Path, event_ids, weights):
    with uproot.recreate(path) as fout:
        fout["events"] = {
            "event_id": np.asarray(event_ids, dtype=np.int64),
            "weight": np.asarray(weights, dtype=np.float64),
            "h_pt": np.asarray([10.0, 20.0], dtype=np.float64),
        }


def test_merge_root_files_rewrites_event_ids_and_preserves_source_ids(tmp_path):
    wp = tmp_path / "wp.root"
    wm = tmp_path / "wm.root"
    out = tmp_path / "wh.root"

    _write_root(wp, [0, 1], [1.0, 2.0])
    _write_root(wm, [0, 1], [3.0, 4.0])

    merge_root_files(
        inputs=[(wp, 0), (wm, 1)],
        output=out,
    )

    tree = uproot.open(out)["events"]
    event_id = tree["event_id"].array(library="np").tolist()
    source_event_id = tree["source_event_id"].array(library="np").tolist()
    subchannel = tree["subchannel"].array(library="np").tolist()
    weight = tree["weight"].array(library="np").tolist()

    assert event_id == [0, 1, 2, 3]
    assert source_event_id == [0, 1, 0, 1]
    assert subchannel == [0, 0, 1, 1]
    assert weight == [1.0, 2.0, 3.0, 4.0]


def test_merge_root_files_rejects_schema_mismatch(tmp_path):
    wp = tmp_path / "wp.root"
    wm = tmp_path / "wm.root"
    out = tmp_path / "bad.root"

    with uproot.recreate(wp) as fout:
        fout["events"] = {"event_id": np.asarray([0]), "weight": np.asarray([1.0])}
    with uproot.recreate(wm) as fout:
        fout["events"] = {"event_id": np.asarray([0]), "weight": np.asarray([1.0]), "extra": np.asarray([2.0])}

    try:
        merge_root_files(inputs=[(wp, 0), (wm, 1)], output=out)
    except ValueError as exc:
        assert "schema mismatch" in str(exc).lower()
    else:
        raise AssertionError("Expected schema mismatch")
