import numpy as np
import uproot


def _load_tree(path):
    tree = uproot.open(path)["events"]
    arrays = {name: tree[name].array(library="np") for name in tree.keys()}
    return arrays


def merge_root_files(inputs, output):
    if not inputs:
        raise ValueError("No input ROOT files provided")

    loaded = [(_load_tree(path), code, path) for path, code in inputs]
    reference_keys = set(loaded[0][0].keys())

    for arrays, _, path in loaded[1:]:
        if set(arrays.keys()) != reference_keys:
            raise ValueError(f"schema mismatch for {path}")

    merged = {key: [] for key in loaded[0][0].keys()}
    merged["source_event_id"] = []
    merged["subchannel"] = []

    next_event_id = 0
    for arrays, code, _ in loaded:
        local_ids = np.asarray(arrays["event_id"], dtype=np.int64)
        count = len(local_ids)

        for key, arr in arrays.items():
            if key == "event_id":
                merged["event_id"].append(np.arange(next_event_id, next_event_id + count, dtype=np.int64))
            else:
                merged[key].append(arr)

        merged["source_event_id"].append(local_ids.astype(np.int64))
        merged["subchannel"].append(np.full(count, code, dtype=np.int8))
        next_event_id += count

    out_arrays = {key: np.concatenate(parts) for key, parts in merged.items()}
    output.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(output) as fout:
        fout["events"] = out_arrays
    return output
