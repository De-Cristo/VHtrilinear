#!/usr/bin/env python3
"""Compute a new weight w_l3 = l3*w_lo + w_lo*c1 (with c1 = w_ew / w_lo) and write a small ROOT file with the new branch.

Usage:
  python scripts/add_l3_weight.py lo_file.root ew_file.root --l3 1.0

By default the script looks for branches named 'event_id', 'w_lo' in the first file and 'event_id', 'w_ew' in the second file. Use
--id-branch, --lo-branch and --ew-branch to change names.

The output file will contain branches: event_id, w_lo, w_ew, w_l3 for the intersection of events present in both input files.
"""
import argparse
import numpy as np
import uproot
import os


def load_weights(fname, id_branch="event_id", weight_branch="w_lo"):
    t = uproot.open(fname)["events"]
    ids = t[id_branch].array(library="np")
    w = t[weight_branch].array(library="np")
    return ids, w


def main():
    p = argparse.ArgumentParser()
    p.add_argument('file_lo')
    p.add_argument('file_ew')
    p.add_argument('--l3', type=float, default=1.0, help='Coefficient l3 used for w_l3')
    p.add_argument('--id-branch', default='event_id')
    p.add_argument('--lo-branch', default='weight')
    p.add_argument('--ew-branch', default='weight')
    args = p.parse_args()

    # derive output filename from file_ew and l3 value; place in same directory
    dname = os.path.dirname(args.file_ew) or '.'
    fname = os.path.basename(args.file_ew)
    if fname.lower().endswith('.root'):
        base = fname[:-5]
    else:
        base = fname
    # replace '_rwgt' suffix with '_l3corr' where present, otherwise append '_l3corr'
    if base.endswith('_rwgt'):
        out_base = base[:-5] + '_l3corr'
    else:
        out_base = base + '_l3corr'
    # format l3 with one decimal to get consistent names (e.g. -10.0 -> '-10.0')
    l3_formatted = "{:.1f}".format(args.l3)
    l3_str = l3_formatted.replace('.', 'p').replace('-', 'm')
    outname = os.path.join(dname, f"{out_base}_{l3_str}.root")

    ids_lo, w_lo = load_weights(args.file_lo, id_branch=args.id_branch, weight_branch=args.lo_branch)
    ids_ew, w_ew = load_weights(args.file_ew, id_branch=args.id_branch, weight_branch=args.ew_branch)

    # build maps
    map_lo = {int(i): float(w) for i, w in zip(ids_lo, w_lo)}
    map_ew = {int(i): float(w) for i, w in zip(ids_ew, w_ew)}

    common = sorted(set(map_lo.keys()).intersection(map_ew.keys()))
    if not common:
        raise SystemExit('No common events found')

    # create arrays aligned to common events
    event_ids = np.array(common, dtype=np.int64)
    w_lo_arr = np.array([map_lo[i] for i in common], dtype=np.float64)
    w_ew_arr = np.array([map_ew[i] for i in common], dtype=np.float64)

    # compute c1 and w_l3
    # c1 = w_ew / w_lo (guard divide by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        c1 = np.where(w_lo_arr != 0.0, w_ew_arr / w_lo_arr, 0.0)
    # user-specified formula: w_l3 = l3*w_lo + w_lo*c1
    # algebraically this equals l3*w_lo + w_ew (since w_lo*c1 == w_ew when w_lo!=0)

    # turn scalar parameters into arrays so broadcasting is explicit and per-event
    n_events = len(event_ids)
    abs_l3_scalar = float(args.l3)
    abs_l3 = np.full(n_events, abs_l3_scalar, dtype=np.float64)
    delta_ZH_scalar = -1.536e-3
    delta_ZH = np.full(n_events, delta_ZH_scalar, dtype=np.float64)

    # ZH_BSM depends on abs_l3 and delta_ZH; compute per-event
    ZH_BSM = 1.0 / (1.0 - (abs_l3**2 - 1.0) * delta_ZH)

    # delta_l3 becomes an array (depends on c1 per-event)
    delta_l3 = (ZH_BSM - 1.0) * (1.0 + delta_ZH) + (ZH_BSM * abs_l3 - 1.0) * c1

    # final per-event weight
    w_l3 = delta_l3 * w_lo_arr + w_lo_arr * c1

    # Also create a separate output (derived from file_lo) containing weights corrected as
    # w_nlo_l3 = w_lo * (1 + c1 + delta_ZH), saved under branch name 'weight'.
    # Derive filename from file_lo
    if args.file_lo.lower().endswith('.root'):
        base_lo = args.file_lo[:-5]
    else:
        base_lo = args.file_lo
    out_lo_name = f"{base_lo}_l3corr.root"

    # compute w_nlo_l3 per-event
    w_nlo_l3 = w_lo_arr * (1.0 + c1 + delta_ZH)

    # update w_l3 to include the LO-corrected term (w_nlo_l3)
    w_l3 = delta_l3 * w_lo_arr + w_nlo_l3

    # read all branches from the first (LO) file so we can carry them over as well
    t1 = uproot.open(args.file_lo)['events']
    branches1 = list(t1.keys())

    # read arrays for all branches from first file
    branch_arrays1 = {}
    for b in branches1:
        try:
            branch_arrays1[b] = t1[b].array(library='np')
        except Exception:
            branch_arrays1[b] = t1[b].array(library='np')

    ids_from_lo = branch_arrays1.get(args.id_branch, ids_lo)
    # build maps for LO branches keyed by event id
    branch_maps1 = {}
    for b, arr in branch_arrays1.items():
        if b == args.id_branch:
            continue
        branch_maps1[b] = {int(i): v for i, v in zip(ids_from_lo, arr)}

    # prepare out1 dict: include event id and all branches from LO file aligned to 'common'
    out1_dict = {}
    out1_dict[args.id_branch] = event_ids
    for b in branches1:
        if b == args.id_branch:
            continue
        if b == args.lo_branch:
            # replace LO weight branch by the newly computed w_nlo_l3
            out1_dict[b] = w_nlo_l3
        else:
            out1_dict[b] = np.array([branch_maps1[b][i] for i in common])

    with uproot.recreate(out_lo_name) as f1:
        f1['events'] = out1_dict
    print(f'Wrote corrected LO file: {out_lo_name} (weight branch contains w_nlo_l3)')
    
    # read all branches from the second (reweighted) file so we can carry them over
    t2 = uproot.open(args.file_ew)["events"]
    branches = list(t2.keys())

    # read arrays for all branches from second file
    branch_arrays = {}
    for b in branches:
        try:
            branch_arrays[b] = t2[b].array(library="np")
        except Exception:
            # fallback: read as numpy (may produce object arrays for jagged branches)
            branch_arrays[b] = t2[b].array(library="np")

    # build maps for each branch keyed by event id
    branch_maps = {}
    ids_from_ew = branch_arrays.get(args.id_branch, ids_ew)
    for b, arr in branch_arrays.items():
        if b == args.id_branch:
            continue
        branch_maps[b] = {int(i): v for i, v in zip(ids_from_ew, arr)}

    # --- compute new NLO weight w_nlo_ew and write separate output file ---
    # w_nlo_ew = ZH_BSM * w_lo * (l3*c1 - c1 + Kew)
    Kew_scalar = 0.947
    Kew = np.full(n_events, Kew_scalar, dtype=np.float64)
    # abs_l3, c1, ZH_BSM and w_lo_arr are defined earlier and are per-event arrays
    w_nlo_ew = ZH_BSM * w_lo_arr * (abs_l3 * c1 - c1 + Kew)

    # prepare output dict for the new NLO EW-weighted file (aligned to 'common')
    out_nloew = {}
    out_nloew[args.id_branch] = event_ids
    for b in branches:
        if b == args.id_branch:
            continue
        if b == args.ew_branch:
            out_nloew[b] = w_nlo_ew
        else:
            out_nloew[b] = np.array([branch_maps[b][i] for i in common])

    # name for the new file
    out_nloew_name = os.path.join(dname, f"{out_base}_nloew_{l3_str}.root")
    with uproot.recreate(out_nloew_name) as f_ne:
        f_ne['events'] = out_nloew
    print(f'Wrote NLO-EW weighted file: {out_nloew_name} (branch "{args.ew_branch}" contains w_nlo_ew)')

    # prepare output dict: include event id and all branches from second file
    out_dict = {}
    out_dict[args.id_branch] = event_ids
    for b in branches:
        if b == args.id_branch:
            continue
        if b == args.ew_branch:
            # replace weight branch by the newly computed w_l3
            out_dict[b] = w_l3
        else:
            # align values to the common event ordering
            out_dict[b] = np.array([branch_maps[b][i] for i in common])

    with uproot.recreate(outname) as f:
        f["events"] = out_dict

    # print summary
    n_total = len(event_ids)
    n_zero_lo = np.count_nonzero(w_lo_arr == 0.0)
    print(f'Wrote {outname} with {n_total} events (events with w_lo==0: {n_zero_lo})')


if __name__ == '__main__':
    raise SystemExit(main())
