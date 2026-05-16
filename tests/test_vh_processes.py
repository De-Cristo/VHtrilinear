from pathlib import Path

from scripts.vh_processes import (
    PUBLIC_PROCESSES,
    get_public_process,
    get_output_dir,
    get_subchannel_code,
)


def test_public_processes_expose_expected_topology():
    assert set(PUBLIC_PROCESSES) == {"zh", "wh"}
    assert [s.key for s in PUBLIC_PROCESSES["zh"].subchannels] == ["zh"]
    assert [s.key for s in PUBLIC_PROCESSES["wh"].subchannels] == ["wh_plus", "wh_minus"]


def test_output_dirs_are_process_isolated():
    repo_root = Path("/repo")
    assert get_output_dir(repo_root, "zh") == repo_root / "output" / "zh"
    assert get_output_dir(repo_root, "wh") == repo_root / "output" / "wh"


def test_wh_subchannel_codes_are_stable():
    assert get_subchannel_code("wh_plus") == 0
    assert get_subchannel_code("wh_minus") == 1


def test_get_public_process_rejects_unknown_keys():
    try:
        get_public_process("bad")
    except KeyError as exc:
        assert "Unknown process" in str(exc)
    else:
        raise AssertionError("Expected KeyError for unknown process")
