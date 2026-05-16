from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SubchannelSpec:
    key: str
    display_name: str
    proc_card: str
    mc_dir: str
    me_dir: str
    subchannel_code: int
    process_label: str
    vector_label: str


@dataclass(frozen=True)
class PublicProcessSpec:
    key: str
    display_name: str
    output_dirname: str
    run_card: str
    param_card: str
    delta_h: float
    k_ew: float
    process_label: str
    vector_label: str
    subchannels: tuple[SubchannelSpec, ...]


PUBLIC_PROCESSES = {
    "zh": PublicProcessSpec(
        key="zh",
        display_name="ZH",
        output_dirname="zh",
        run_card="cards/zh/run_card.dat",
        param_card="cards/zh/param_card.dat",
        delta_h=-1.536e-3,
        k_ew=0.947,
        process_label="ZH (13.6 TeV)",
        vector_label="Z",
        subchannels=(
            SubchannelSpec(
                key="zh",
                display_name="ZH",
                proc_card="cards/zh/proc_card.dat",
                mc_dir="zh_MC",
                me_dir="zh_ME",
                subchannel_code=0,
                process_label="ZH (13.6 TeV)",
                vector_label="Z",
            ),
        ),
    ),
    "wh": PublicProcessSpec(
        key="wh",
        display_name="WH",
        output_dirname="wh",
        run_card="cards/wh/run_card.dat",
        param_card="cards/wh/param_card.dat",
        delta_h=-1.536e-3,
        k_ew=0.947,
        process_label="WH (13.6 TeV)",
        vector_label="W",
        subchannels=(
            SubchannelSpec(
                key="wh_plus",
                display_name="W+H",
                proc_card="cards/wh/proc_card_wp.dat",
                mc_dir="whp_MC",
                me_dir="whp_ME",
                subchannel_code=0,
                process_label="WH (13.6 TeV)",
                vector_label="W",
            ),
            SubchannelSpec(
                key="wh_minus",
                display_name="W-H",
                proc_card="cards/wh/proc_card_wm.dat",
                mc_dir="whm_MC",
                me_dir="whm_ME",
                subchannel_code=1,
                process_label="WH (13.6 TeV)",
                vector_label="W",
            ),
        ),
    ),
}


def get_public_process(process_key: str) -> PublicProcessSpec:
    try:
        return PUBLIC_PROCESSES[process_key]
    except KeyError as exc:
        raise KeyError(f"Unknown process: {process_key}") from exc


def get_output_dir(repo_root: Path, process_key: str) -> Path:
    spec = get_public_process(process_key)
    return repo_root / "output" / spec.output_dirname


def get_subchannel_code(subchannel_key: str) -> int:
    for process in PUBLIC_PROCESSES.values():
        for subchannel in process.subchannels:
            if subchannel.key == subchannel_key:
                return subchannel.subchannel_code
    raise KeyError(f"Unknown subchannel: {subchannel_key}")
