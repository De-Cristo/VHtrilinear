from pathlib import Path

import uproot

from scripts.lhe_to_root import run


WH_LHE = """<LesHouchesEvents>
<event>
4 1 1.0 0 0 0
1 -1 0 0 0 0 0 0 6500 6500 0 0 0
-2 -1 0 0 0 0 0 0 -6500 6500 0 0 0
24 1 0 0 0 0 30 0 40 60 0 0 0
25 1 0 0 0 0 -30 0 -40 140 125 0 0
</event>
</LesHouchesEvents>
"""


def test_lhe_to_root_writes_subchannel_metadata_for_wh(tmp_path):
    lhe = tmp_path / "wh.lhe"
    out = tmp_path / "wh.root"
    lhe.write_text(WH_LHE)

    run(str(lhe), str(out), ebeam=6800.0, subchannel_code=1)

    tree = uproot.open(out)["events"]
    assert tree["event_id"].array(library="np").tolist() == [0]
    assert tree["source_event_id"].array(library="np").tolist() == [0]
    assert tree["subchannel"].array(library="np").tolist() == [1]
