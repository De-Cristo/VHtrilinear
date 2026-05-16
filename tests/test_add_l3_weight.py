from pathlib import Path

import numpy as np
import uproot

from scripts.add_l3_weight import main


def _write_root(path: Path, weights_lo, weights_rw):
    with uproot.recreate(path) as fout:
        fout["events"] = {
            "event_id": np.asarray([0, 1], dtype=np.int64),
            "source_event_id": np.asarray([10, 11], dtype=np.int64),
            "subchannel": np.asarray([0, 1], dtype=np.int8),
            "weight": np.asarray(weights_lo if "lo" in path.name else weights_rw, dtype=np.float64),
            "h_pt": np.asarray([10.0, 20.0], dtype=np.float64),
        }


def test_add_l3_weight_preserves_metadata_branches(tmp_path, monkeypatch):
    lo = tmp_path / "events_lo.root"
    rw = tmp_path / "events_rwgt.root"
    _write_root(lo, [1.0, 2.0], [0.0, 0.0])
    _write_root(rw, [0.0, 0.0], [1.1, 2.2])

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("sys.argv", ["add_l3_weight.py", str(lo), str(rw), "--l3", "1.0"])
    main()

    out = uproot.open(tmp_path / "events_l3corr_1p0.root")["events"]
    assert out["source_event_id"].array(library="np").tolist() == [10, 11]
    assert out["subchannel"].array(library="np").tolist() == [0, 1]
