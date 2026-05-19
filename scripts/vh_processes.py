from pathlib import Path


class SubchannelSpec(object):
    __slots__ = (
        "key",
        "display_name",
        "proc_card",
        "mc_dir",
        "me_dir",
        "subchannel_code",
        "process_label",
        "vector_label",
    )

    def __init__(
        self,
        key,
        display_name,
        proc_card,
        mc_dir,
        me_dir,
        subchannel_code,
        process_label,
        vector_label,
    ):
        self.key = key
        self.display_name = display_name
        self.proc_card = proc_card
        self.mc_dir = mc_dir
        self.me_dir = me_dir
        self.subchannel_code = subchannel_code
        self.process_label = process_label
        self.vector_label = vector_label


class PublicProcessSpec(object):
    __slots__ = (
        "key",
        "display_name",
        "output_dirname",
        "run_card",
        "param_card",
        "delta_h",
        "k_ew",
        "process_label",
        "vector_label",
        "subchannels",
    )

    def __init__(
        self,
        key,
        display_name,
        output_dirname,
        run_card,
        param_card,
        delta_h,
        k_ew,
        process_label,
        vector_label,
        subchannels,
    ):
        self.key = key
        self.display_name = display_name
        self.output_dirname = output_dirname
        self.run_card = run_card
        self.param_card = param_card
        self.delta_h = delta_h
        self.k_ew = k_ew
        self.process_label = process_label
        self.vector_label = vector_label
        self.subchannels = subchannels


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


def get_public_process(process_key):
    try:
        return PUBLIC_PROCESSES[process_key]
    except KeyError as exc:
        raise KeyError(f"Unknown process: {process_key}") from exc


def get_output_dir(repo_root, process_key):
    spec = get_public_process(process_key)
    return repo_root / "output" / spec.output_dirname


def get_subchannel_code(subchannel_key):
    for process in PUBLIC_PROCESSES.values():
        for subchannel in process.subchannels:
            if subchannel.key == subchannel_key:
                return subchannel.subchannel_code
    raise KeyError(f"Unknown subchannel: {subchannel_key}")
