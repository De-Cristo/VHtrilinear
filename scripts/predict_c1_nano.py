#!/usr/bin/env python3
"""
predict_c1_nano.py — Predict per-event C1 (NLO EW correction factor) on NanoAOD events
using LHE-level particle information and the trained XGBoost regressor.

Workflow:
  1. Read NanoAOD ROOT file(s) containing LHEPart branches
  2. Identify Z (pdgId=23) and H (pdgId=25) from LHEPart arrays
  3. Compute the 6 truth-level kinematic features used by the C1 regressor:
     h_pt, v_pt, vh_m, cos_theta_star, h_y, vh_delta_eta
  4. Run XGBoost inference to predict C1_NLO per event
  5. Save output ROOT file with predictions + LHE kinematics + reco-level info

Usage:
  python3 scripts/predict_c1_nano.py \\
      --input nanoAOD_temp/*.root \\
      --model output/c1_regressor/c1_regressor.json \\
      --output output/nano_c1_predictions.root \\
      --ebeam 6800

Options:
  --input       Input NanoAOD ROOT file(s)                [required]
  --model       Path to trained XGBoost model (.json)     [required]
  --output      Output ROOT file with predictions         [default: output/nano_c1_predictions.root]
  --ebeam       Beam energy per proton [GeV]              [default: 6800]
  --max-events  Max events to process (0 = all)           [default: 0]
"""
import argparse
import os
import sys
import glob
import math
import numpy as np

try:
    import uproot
    import awkward as ak
    import xgboost as xgb
except ImportError as e:
    print(f"Missing dependency: {e}")
    print("Install with: pip install uproot awkward xgboost numpy")
    sys.exit(1)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica', 'DejaVu Sans'],
        'axes.linewidth': 1.2,
        'font.size': 12,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
    })
    try:
        mpl.rcParams['text.usetex'] = True
    except Exception:
        mpl.rcParams['text.usetex'] = False
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


# ── Feature list: must match training order in train_c1_regressor.py ──
FEATURES = ['h_pt', 'v_pt', 'vh_m', 'cos_theta_star', 'h_y', 'vh_delta_eta']
FEATURE_LABELS = {
    'h_pt': r'$p_T(H)$ [GeV]',
    'v_pt': r'$p_T(Z)$ [GeV]',
    'vh_m': r'$m(ZH)$ [GeV]',
    'cos_theta_star': r'$\cos\theta^*$',
    'h_y': r'$y(H)$',
    'vh_delta_eta': r'$|\Delta\eta(Z, H)|$',
}

EXTRA_NANO_FEATURES = ['vh_pt_gen', 'n_lhe_extra']
FEATURE_LABELS.update({
    'vh_pt_gen': r'$p_T(ZH)_{\rm LHE}$ [GeV]',
    'n_lhe_extra': r'$N_{\rm extra}^{\rm LHE}$',
})

# ── LHEPart branches to read ──
LHE_BRANCHES = ['LHEPart_pt', 'LHEPart_eta', 'LHEPart_phi', 'LHEPart_mass', 'LHEPart_pdgId', 'LHEPart_status']

# ── Reco-level branches to carry over (for future reco-regressor training) ──
# These are read if available, otherwise skipped
RECO_BRANCHES_OPTIONAL = [
    'nJet', 'Jet_pt', 'Jet_eta', 'Jet_phi', 'Jet_mass', 'Jet_btagDeepFlavB',
    'nMuon', 'Muon_pt', 'Muon_eta', 'Muon_phi',
    'nElectron', 'Electron_pt', 'Electron_eta', 'Electron_phi',
    'MET_pt', 'MET_phi',
    'nFatJet', 'FatJet_pt', 'FatJet_eta', 'FatJet_phi', 'FatJet_mass',
    'FatJet_msoftdrop', 'FatJet_particleNet_HbbvsQCD',
]


# ═══════════════════════════════════════════════════════
# Kinematic computations (vectorized, matching lhe_to_root.py)
# ═══════════════════════════════════════════════════════

def compute_px_py_pz_e(pt, eta, phi, mass):
    """Convert (pt, eta, phi, mass) to (px, py, pz, E) arrays."""
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = pt * np.sinh(eta)
    e = np.sqrt(mass**2 + px**2 + py**2 + pz**2)
    return px, py, pz, e


def compute_rapidity(e, pz):
    """Compute rapidity y = 0.5 * ln((E+pz)/(E-pz))."""
    with np.errstate(divide='ignore', invalid='ignore'):
        denom = e - pz
        numer = e + pz
        safe = (denom > 0) & np.isfinite(numer) & np.isfinite(denom)
        y = np.where(safe, 0.5 * np.log(numer / denom), np.sign(pz) * 10.0)
    return y


def compute_pseudorapidity(px, py, pz):
    """Compute pseudorapidity eta."""
    pt_val = np.hypot(px, py)
    p = np.hypot(pt_val, pz)
    with np.errstate(divide='ignore', invalid='ignore'):
        safe = np.abs(pz) < p
        eta = np.where(safe, 0.5 * np.log((p + pz) / (p - pz)), np.sign(pz) * 10.0)
    return eta


def lorentz_boost_z_component(e, px, py, pz, beta_x, beta_y, beta_z):
    """Lorentz boost 4-vectors by velocity beta = (beta_x, beta_y, beta_z).
    All inputs are 1D arrays of same length.
    Returns boosted (E', px', py', pz').

    Uses the standard formula: p' = p + β[(γ-1)(β·p)/β² - γE]
    """
    b2 = beta_x**2 + beta_y**2 + beta_z**2
    gamma = 1.0 / np.sqrt(np.maximum(1.0 - b2, 1e-30))
    bp = beta_x * px + beta_y * py + beta_z * pz  # dot product

    # (gamma - 1) / beta^2
    with np.errstate(divide='ignore', invalid='ignore'):
        factor = np.where(b2 > 0, (gamma - 1.0) / b2, 0.0)

    e_new = gamma * (e - bp)
    # p'_i = p_i + beta_i * [factor * bp - gamma * E]
    px_new = px + factor * bp * beta_x - gamma * beta_x * e
    py_new = py + factor * bp * beta_y - gamma * beta_y * e
    pz_new = pz + factor * bp * beta_z - gamma * beta_z * e

    return e_new, px_new, py_new, pz_new


def compute_cos_theta_star(h_px, h_py, h_pz, h_e, v_px, v_py, v_pz, v_e, ebeam):
    """Compute cos(theta*) of the Higgs in the VH rest frame w.r.t. beam axis.
    Matches the logic in lhe_to_root.py lines 265-288.
    """
    n = len(h_px)

    # VH system 4-vector
    q_e = h_e + v_e
    q_px = h_px + v_px
    q_py = h_py + v_py
    q_pz = h_pz + v_pz

    # Boost velocity
    with np.errstate(divide='ignore', invalid='ignore'):
        beta_x = np.where(q_e != 0, q_px / q_e, 0.0)
        beta_y = np.where(q_e != 0, q_py / q_e, 0.0)
        beta_z = np.where(q_e != 0, q_pz / q_e, 0.0)

    # Boost Higgs to VH rest frame
    h_e_star, h_px_star, h_py_star, h_pz_star = lorentz_boost_z_component(
        h_e, h_px, h_py, h_pz, beta_x, beta_y, beta_z
    )

    # Boost beam to VH rest frame
    beam_e = np.full(n, ebeam)
    beam_px = np.zeros(n)
    beam_py = np.zeros(n)
    beam_pz = np.full(n, ebeam)
    b_e_star, b_px_star, b_py_star, b_pz_star = lorentz_boost_z_component(
        beam_e, beam_px, beam_py, beam_pz, beta_x, beta_y, beta_z
    )

    # cos(theta*) = dot(H*, beam*) / (|H*| * |beam*|)
    norm_h = np.sqrt(h_px_star**2 + h_py_star**2 + h_pz_star**2)
    norm_b = np.sqrt(b_px_star**2 + b_py_star**2 + b_pz_star**2)
    dot = h_px_star * b_px_star + h_py_star * b_py_star + h_pz_star * b_pz_star

    with np.errstate(divide='ignore', invalid='ignore'):
        cos_th = np.where(
            (norm_h > 0) & (norm_b > 0),
            dot / (norm_h * norm_b),
            np.nan
        )
    return np.clip(cos_th, -1.0, 1.0)


# ═══════════════════════════════════════════════════════
# Main processing
# ═══════════════════════════════════════════════════════

def process_nanoaod(input_files, model_path, output_path, ebeam=6800.0, max_events=0):
    """Process NanoAOD files: extract LHE kinematics, predict C1, save output."""

    print(f"\n{'='*60}")
    print(f"  NanoAOD → C1 Prediction Pipeline")
    print(f"{'='*60}")
    print(f"  Input files:  {len(input_files)}")
    print(f"  Model:        {model_path}")
    print(f"  Output:       {output_path}")
    print(f"  Beam energy:  {ebeam} GeV")

    # ── Load model ──
    print(f"\n── Loading XGBoost model ──")
    model = xgb.Booster()
    model.load_model(model_path)
    print(f"  [✓] Model loaded")

    # ── Accumulate events from all input files ──
    all_h_pt = []
    all_v_pt = []
    all_vh_m = []
    all_cos_theta_star = []
    all_h_y = []
    all_vh_delta_eta = []
    all_vh_pt_gen = []
    all_n_lhe_extra = []
    all_event_file_idx = []  # track which file each event came from
    all_event_local_idx = []  # event index within its file

    # Also accumulate reco-level info for future training
    reco_data = {b: [] for b in RECO_BRANCHES_OPTIONAL}
    reco_available = set()

    total_events = 0
    total_zh_events = 0

    for fi, fpath in enumerate(input_files):
        print(f"\n── Processing {os.path.basename(fpath)} ──")

        f = uproot.open(fpath)
        tree = f['Events']
        n_entries = tree.num_entries
        print(f"  Total entries: {n_entries}")

        # Determine how many events to read
        stop = None
        if max_events > 0:
            remaining = max_events - total_events
            if remaining <= 0:
                break
            stop = min(n_entries, remaining)

        # Read LHEPart branches
        available_keys = set(tree.keys())
        missing = [b for b in LHE_BRANCHES if b not in available_keys]
        if missing:
            print(f"  [!] Missing LHE branches: {missing} — skipping file")
            continue

        lhe = tree.arrays(LHE_BRANCHES, library='ak', entry_stop=stop)
        n_read = len(lhe)
        total_events += n_read

        # Find Z (pdgId=23) and H (pdgId=25) per event
        pdgid = lhe['LHEPart_pdgId']
        pt_arr = lhe['LHEPart_pt']
        eta_arr = lhe['LHEPart_eta']
        phi_arr = lhe['LHEPart_phi']
        mass_arr = lhe['LHEPart_mass']

        # Mask for Z and H bosons
        is_z = (pdgid == 23)
        is_h = (pdgid == 25)

        # Check which events have both Z and H
        has_z = ak.any(is_z, axis=1)
        has_h = ak.any(is_h, axis=1)
        has_both = has_z & has_h

        n_zh = int(ak.sum(has_both))
        print(f"  Events with Z+H: {n_zh} / {n_read}")

        if n_zh == 0:
            continue

        # Get the local indices of good events
        good_idx = np.where(ak.to_numpy(has_both))[0]

        # Slice to only good events
        pt_zh = pt_arr[has_both]
        eta_zh = eta_arr[has_both]
        phi_zh = phi_arr[has_both]
        mass_zh = mass_arr[has_both]
        pdgid_zh = pdgid[has_both]

        is_z_zh = (pdgid_zh == 23)
        is_h_zh = (pdgid_zh == 25)

        # Extract first Z and first H per event using ak.firsts + boolean masking
        z_pt = ak.to_numpy(ak.firsts(pt_zh[is_z_zh])).astype(np.float64)
        z_eta = ak.to_numpy(ak.firsts(eta_zh[is_z_zh])).astype(np.float64)
        z_phi = ak.to_numpy(ak.firsts(phi_zh[is_z_zh])).astype(np.float64)
        z_mass = ak.to_numpy(ak.firsts(mass_zh[is_z_zh])).astype(np.float64)

        h_pt_val = ak.to_numpy(ak.firsts(pt_zh[is_h_zh])).astype(np.float64)
        h_eta_val = ak.to_numpy(ak.firsts(eta_zh[is_h_zh])).astype(np.float64)
        h_phi_val = ak.to_numpy(ak.firsts(phi_zh[is_h_zh])).astype(np.float64)
        h_mass_val = ak.to_numpy(ak.firsts(mass_zh[is_h_zh])).astype(np.float64)

        # Convert to (px, py, pz, E)
        z_px, z_py, z_pz, z_e = compute_px_py_pz_e(z_pt, z_eta, z_phi, z_mass)
        h_px, h_py, h_pz, h_e = compute_px_py_pz_e(h_pt_val, h_eta_val, h_phi_val, h_mass_val)

        # Compute derived kinematics
        # vh_m: invariant mass of ZH system
        vh_e = z_e + h_e
        vh_px = z_px + h_px
        vh_py = z_py + h_py
        vh_pz = z_pz + h_pz
        vh_m2 = vh_e**2 - (vh_px**2 + vh_py**2 + vh_pz**2)
        vh_m_val = np.sqrt(np.maximum(vh_m2, 0.0))

        # h_y: Higgs rapidity
        h_y_val = compute_rapidity(h_e, h_pz)

        # vh_delta_eta
        z_pseudoeta = compute_pseudorapidity(z_px, z_py, z_pz)
        h_pseudoeta = compute_pseudorapidity(h_px, h_py, h_pz)
        vh_deta = np.abs(z_pseudoeta - h_pseudoeta)

        # cos_theta_star
        cos_th = compute_cos_theta_star(h_px, h_py, h_pz, h_e, z_px, z_py, z_pz, z_e, ebeam)

        # Extra NANO validation variables
        vh_pt_gen_val = np.sqrt(vh_px**2 + vh_py**2)
        
        is_qg_zh = (abs(pdgid_zh) <= 5) | (pdgid_zh == 21)
        # Assuming initial partons may have exactly 0 Pt (if LHEPart_pt > 0.1, it's radiation)
        n_extra_lhe_zh = ak.sum(is_qg_zh & (pt_zh > 0.1), axis=1)
        n_extra_lhe_np = ak.to_numpy(n_extra_lhe_zh).astype(np.int32)

        # Accumulate
        all_h_pt.append(h_pt_val)
        all_v_pt.append(z_pt)
        all_vh_m.append(vh_m_val)
        all_cos_theta_star.append(cos_th)
        all_h_y.append(h_y_val)
        all_vh_delta_eta.append(vh_deta)
        all_vh_pt_gen.append(vh_pt_gen_val)
        all_n_lhe_extra.append(n_extra_lhe_np)
        all_event_file_idx.append(np.full(n_zh, fi, dtype=np.int32))
        all_event_local_idx.append(good_idx.astype(np.int64))

        total_zh_events += n_zh
        print(f"  [✓] Extracted {n_zh} ZH events")

        # Read reco-level branches if available (for future reco-level training)
        for rb in RECO_BRANCHES_OPTIONAL:
            if rb in available_keys:
                reco_available.add(rb)

    if total_zh_events == 0:
        print("\n  [!] No ZH events found in any input file. Exiting.")
        sys.exit(1)

    # ── Concatenate all events ──
    print(f"\n── Merging {total_zh_events} ZH events from {len(input_files)} files ──")
    h_pt_all = np.concatenate(all_h_pt)
    v_pt_all = np.concatenate(all_v_pt)
    vh_m_all = np.concatenate(all_vh_m)
    cos_th_all = np.concatenate(all_cos_theta_star)
    h_y_all = np.concatenate(all_h_y)
    vh_deta_all = np.concatenate(all_vh_delta_eta)
    vh_pt_gen_all = np.concatenate(all_vh_pt_gen)
    n_lhe_extra_all = np.concatenate(all_n_lhe_extra)
    file_idx_all = np.concatenate(all_event_file_idx)
    local_idx_all = np.concatenate(all_event_local_idx)

    # ── Build feature matrix (same order as training) ──
    X = np.column_stack([h_pt_all, v_pt_all, vh_m_all, cos_th_all, h_y_all, vh_deta_all])
    print(f"  Feature matrix: {X.shape}")

    # Stats
    print(f"\n  Feature statistics:")
    for i, feat in enumerate(FEATURES):
        col = X[:, i]
        print(f"    {feat:20s}: mean={np.nanmean(col):10.4f}  std={np.nanstd(col):10.4f}  "
              f"min={np.nanmin(col):10.4f}  max={np.nanmax(col):10.4f}")

    # ── Run XGBoost inference ──
    print(f"\n── Running C1 prediction ──")
    # Use actual feature names to match training
    dmat = xgb.DMatrix(X, feature_names=FEATURES)
    c1_pred = model.predict(dmat)

    print(f"  Predicted C1 stats:")
    print(f"    Mean:   {np.mean(c1_pred):.4f}%")
    print(f"    Std:    {np.std(c1_pred):.4f}%")
    print(f"    Min:    {np.min(c1_pred):.4f}%")
    print(f"    Max:    {np.max(c1_pred):.4f}%")

    # ── Compute R(kappa) scaling weights ──
    print(f"\n── Computing R(kappa) weights ──")
    target_kappas = [0, 10, 1, 2, 5, -10, -2, -5]
    delta_ZH = -1.536e-3

    def compute_R_kappa_norm(c1_percent, k):
        c_dec = c1_percent / 100.0
        z_h = 1.0 / (1.0 - (k**2 - 1.0) * delta_ZH)
        num = z_h * (1.0 + k * c_dec + delta_ZH)
        den = 1.0 + c_dec + delta_ZH
        return num / den

    def compute_R_kappa_abs(c1_percent, k):
        c_dec = c1_percent / 100.0
        z_h = 1.0 / (1.0 - (k**2 - 1.0) * delta_ZH)
        k_ew = 0.947
        return z_h * (k * c_dec - c_dec + k_ew) / k_ew

    kappa_weights = {}
    for k in target_kappas:
        w_norm = compute_R_kappa_norm(c1_pred, k)
        w_abs = compute_R_kappa_abs(c1_pred, k)
        kappa_weights[f'weight_kappa_{k}'] = w_norm.astype(np.float32)
        kappa_weights[f'weight_kappa_abs_{k}'] = w_abs.astype(np.float32)
        print(f"  [✓] kappa={k} (norm & abs weights)")

    # ── Save output ROOT file ──
    print(f"\n── Saving output ──")
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)

    out_dict = {
        'c1_pred': c1_pred.astype(np.float32),
        'h_pt': h_pt_all.astype(np.float32),
        'v_pt': v_pt_all.astype(np.float32),
        'vh_m': vh_m_all.astype(np.float32),
        'cos_theta_star': cos_th_all.astype(np.float32),
        'h_y': h_y_all.astype(np.float32),
        'vh_delta_eta': vh_deta_all.astype(np.float32),
        'vh_pt_gen': vh_pt_gen_all.astype(np.float32),
        'n_lhe_extra': n_lhe_extra_all.astype(np.int32),
        'file_idx': file_idx_all,
        'event_idx': local_idx_all,
    }
    out_dict.update(kappa_weights)

    with uproot.recreate(output_path) as fout:
        fout['events'] = out_dict

    print(f"  [✓] {output_path}")
    print(f"       {total_zh_events} events, {len(out_dict)} branches:")
    for k in out_dict:
        print(f"         {k}")

    # ── Summary ──
    print(f"\n{'='*60}")
    print(f"  C1 prediction complete!")
    print(f"  Input:  {total_events} total events → {total_zh_events} ZH events")
    print(f"  Output: {output_path}")
    print(f"  Model:  {model_path}")
    print(f"{'='*60}")

    # Return data for validation
    nano_data = {
        'c1_pred': c1_pred,
        'h_pt': h_pt_all, 'v_pt': v_pt_all, 'vh_m': vh_m_all,
        'cos_theta_star': cos_th_all, 'h_y': h_y_all, 'vh_delta_eta': vh_deta_all,
        'vh_pt_gen': vh_pt_gen_all, 'n_lhe_extra': n_lhe_extra_all,
    }
    nano_data.update(kappa_weights)
    return output_path, nano_data


# ═══════════════════════════════════════════════════════
# Validation plots: LO truth C1 vs NanoAOD predicted C1
# ═══════════════════════════════════════════════════════

def load_lo_truth_c1(lo_file, rw_file):
    """Load the true C1 and kinematics from LO+reweighted ROOT files."""
    t_lo = uproot.open(lo_file)['events']
    t_rw = uproot.open(rw_file)['events']

    ids_lo = t_lo['event_id'].array(library='np').astype(int)
    w_lo = t_lo['weight'].array(library='np')
    ids_rw = t_rw['event_id'].array(library='np').astype(int)
    w_rw = t_rw['weight'].array(library='np')

    map_lo = dict(zip(ids_lo, w_lo))
    map_rw = dict(zip(ids_rw, w_rw))
    common = sorted(set(map_lo) & set(map_rw))

    w_lo_arr = np.array([map_lo[i] for i in common])
    w_rw_arr = np.array([map_rw[i] for i in common])
    nonzero = w_lo_arr != 0
    c1_true = np.zeros(len(common))
    c1_true[nonzero] = (w_rw_arr[nonzero] / w_lo_arr[nonzero]) * 100.0

    # Load kinematics from LO file
    feats = {}
    for f in FEATURES:
        arr = t_lo[f].array(library='np')
        feat_map = dict(zip(ids_lo, arr))
        feats[f] = np.array([feat_map[i] for i in common])

    # Filter nonzero
    c1_true = c1_true[nonzero]
    for f in FEATURES:
        feats[f] = feats[f][nonzero]

    return c1_true, feats


def make_validation_plots(nano_data, lo_file, rw_file, plotdir):
    """Generate validation plots comparing LO truth C1 vs NanoAOD predicted C1."""
    if not HAS_MPL:
        print("  [!] matplotlib not available — skipping validation plots")
        return

    os.makedirs(plotdir, exist_ok=True)
    print(f"\n── Generating validation plots ──")

    # Load LO truth
    print(f"  Loading LO truth from {lo_file} + {rw_file}...")
    c1_lo, feats_lo = load_lo_truth_c1(lo_file, rw_file)
    c1_nano = nano_data['c1_pred']

    print(f"  LO events:   {len(c1_lo)}")
    print(f"  Nano events: {len(c1_nano)}")

    # ── Plot 1: Global C1 distribution comparison ──
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(7, 5.5))

    bins = np.linspace(
        min(c1_lo.min(), c1_nano.min()) - 0.1,
        max(c1_lo.max(), c1_nano.max()) + 0.1,
        60
    )
    h_lo, _, _ = ax_top.hist(c1_lo, bins=bins, histtype='step', color='C0',
                              linewidth=1.5, label='LO truth', density=True)
    h_nano, _, _ = ax_top.hist(c1_nano, bins=bins, histtype='step', color='C1',
                                linewidth=1.5, linestyle='--',
                                label='NanoAOD predicted', density=True)
    ax_top.set_ylabel('Density')
    ax_top.legend(frameon=False)
    ax_top.minorticks_on()

    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(h_lo > 0, h_nano / h_lo, np.nan)
    ax_bot.plot(bin_centers, ratio, 'ko-', markersize=2, linewidth=0.8)
    ax_bot.axhline(1.0, color='red', linestyle='--', linewidth=1)
    ax_bot.set_xlabel(r'$C_1$ [\%]')
    ax_bot.set_ylabel('Nano/LO')
    ax_bot.set_ylim(0.5, 1.5)
    ax_bot.minorticks_on()

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.06)
    out = os.path.join(plotdir, 'c1_distribution_validation.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  [✓] {out}")

    # ── Plot 2: Per-feature C1 profile comparisons ──
    for feat in FEATURES:
        label = FEATURE_LABELS.get(feat, feat)
        x_lo = feats_lo[feat]
        x_nano = nano_data[feat]

        # Determine binning
        if feat == 'cos_theta_star':
            fbins = np.linspace(-1, 1, 21)
        elif feat == 'h_y':
            fbins = np.linspace(-4, 4, 21)
        elif feat == 'vh_delta_eta':
            fbins = np.linspace(0, 6, 21)
        else:
            q95 = max(np.percentile(x_lo, 95), np.percentile(x_nano, 95))
            fbins = np.linspace(0, q95, 21)

        bin_c = 0.5 * (fbins[:-1] + fbins[1:])
        nb = len(bin_c)

        # Compute mean C1 in each bin for LO and Nano
        lo_means = np.full(nb, np.nan)
        lo_errs = np.full(nb, np.nan)
        nano_means = np.full(nb, np.nan)
        nano_errs = np.full(nb, np.nan)

        for ib in range(nb):
            # LO
            mask_lo = (x_lo >= fbins[ib]) & (x_lo < fbins[ib + 1])
            if np.sum(mask_lo) > 10:
                lo_means[ib] = np.mean(c1_lo[mask_lo])
                lo_errs[ib] = np.std(c1_lo[mask_lo]) / np.sqrt(np.sum(mask_lo))
            # Nano
            mask_nano = (x_nano >= fbins[ib]) & (x_nano < fbins[ib + 1])
            if np.sum(mask_nano) > 10:
                nano_means[ib] = np.mean(c1_nano[mask_nano])
                nano_errs[ib] = np.std(c1_nano[mask_nano]) / np.sqrt(np.sum(mask_nano))

        fig, (ax_t, ax_b) = plt.subplots(
            2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(7, 5.5))

        # Top: C1 profiles
        ax_t.errorbar(bin_c, lo_means, yerr=lo_errs, fmt='o',
                      color='C0', markersize=4, linewidth=1.2, capsize=2,
                      label='LO truth')
        offset = (fbins[1] - fbins[0]) * 0.05
        ax_t.errorbar(bin_c + offset, nano_means, yerr=nano_errs, fmt='s',
                      color='C1', markersize=4, linewidth=1.2, capsize=2,
                      label='NanoAOD predicted')
        ax_t.set_ylabel(r'$\langle C_1 \rangle$ [\%]')
        ax_t.legend(frameon=False)
        ax_t.minorticks_on()

        # Bottom: ratio
        valid = ~np.isnan(lo_means) & (lo_means != 0)
        ratio_f = np.full_like(lo_means, np.nan)
        ratio_f[valid] = nano_means[valid] / lo_means[valid]
        ax_b.plot(bin_c, ratio_f, 'ko-', markersize=3, linewidth=1)
        ax_b.axhline(1.0, color='red', linestyle='--', linewidth=1)
        ax_b.set_xlabel(label)
        ax_b.set_ylabel('Nano/LO')
        ax_b.set_ylim(0.9, 1.1)
        ax_b.minorticks_on()

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.06)
        out = os.path.join(plotdir, f'c1_profile_validation_{feat}.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  [✓] {out}")

    # ── Plot 3: Feature distribution comparisons ──
    for feat in FEATURES:
        label = FEATURE_LABELS.get(feat, feat)
        x_lo = feats_lo[feat]
        x_nano = nano_data[feat]

        fig, (ax_t, ax_b) = plt.subplots(
            2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(7, 5.5))

        if feat == 'cos_theta_star':
            fbins = np.linspace(-1, 1, 41)
        elif feat == 'h_y':
            fbins = np.linspace(-4, 4, 41)
        elif feat == 'vh_delta_eta':
            fbins = np.linspace(0, 6, 41)
        else:
            q95 = max(np.percentile(x_lo, 97), np.percentile(x_nano, 97))
            fbins = np.linspace(0, q95, 41)

        h_lo_f, _, _ = ax_t.hist(x_lo, bins=fbins, histtype='step', color='C0',
                                  linewidth=1.5, label='LO events', density=True)
        h_nano_f, _, _ = ax_t.hist(x_nano, bins=fbins, histtype='step', color='C1',
                                    linewidth=1.5, linestyle='--',
                                    label='NanoAOD events', density=True)
        ax_t.set_ylabel('Density')
        ax_t.legend(frameon=False)
        ax_t.minorticks_on()

        bin_c = 0.5 * (fbins[:-1] + fbins[1:])
        with np.errstate(divide='ignore', invalid='ignore'):
            r = np.where(h_lo_f > 0, h_nano_f / h_lo_f, np.nan)
        ax_b.plot(bin_c, r, 'ko-', markersize=2, linewidth=0.8)
        ax_b.axhline(1.0, color='red', linestyle='--', linewidth=1)
        ax_b.set_xlabel(label)
        ax_b.set_ylabel('Nano/LO')
        ax_b.set_ylim(0.5, 1.5)
        ax_b.minorticks_on()

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.06)
        out = os.path.join(plotdir, f'feature_validation_{feat}.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  [✓] {out}")

    # ── Plot 4: Combined Feature Density + C1 Profile ──
    for feat in FEATURES + EXTRA_NANO_FEATURES:
        label = FEATURE_LABELS.get(feat, feat)
        is_extra = feat in EXTRA_NANO_FEATURES
        
        x_nano = nano_data[feat]
        x_lo = None if is_extra else feats_lo[feat]

        fig, (ax_t, ax_b) = plt.subplots(
            2, 1, sharex=True, gridspec_kw={'height_ratios': [2, 1]}, figsize=(7, 6))

        # Determine binning
        if feat == 'cos_theta_star':
            fbins = np.linspace(-1, 1, 21)
        elif feat == 'h_y':
            fbins = np.linspace(-4, 4, 21)
        elif feat == 'vh_delta_eta':
            fbins = np.linspace(0, 6, 21)
        elif feat == 'n_lhe_extra':
            fbins = np.linspace(-0.5, 5.5, 7) # integers 0-5
        else:
            q95 = max(np.percentile(x_lo, 95) if not is_extra else 0, np.percentile(x_nano, 95))
            fbins = np.linspace(0, q95, 21)

        # Top: Feature density (Fraction of weighted events)
        if not is_extra:
            ax_t.hist(x_lo, bins=fbins, histtype='stepfilled', color='C0', alpha=0.3,
                      label='LO events (Calculated)', density=True)
            ax_t.hist(x_lo, bins=fbins, histtype='step', color='C0', linewidth=1.5, density=True)
            
        ax_t.hist(x_nano, bins=fbins, histtype='step', color='C1',
                  linewidth=1.5, linestyle='--',
                  label='NanoAOD events (Predicted)', density=True)
        
        ax_t.set_ylabel('Fraction of events\n(Density)')
        ax_t.legend(frameon=False)
        ax_t.minorticks_on()

        # Bottom: C1 Profile
        bin_c = 0.5 * (fbins[:-1] + fbins[1:])
        nb = len(bin_c)
        lo_means = np.full(nb, np.nan)
        nano_means = np.full(nb, np.nan)

        for ib in range(nb):
            if not is_extra:
                mask_lo = (x_lo >= fbins[ib]) & (x_lo < fbins[ib + 1])
                if np.sum(mask_lo) > 10:
                    lo_means[ib] = np.mean(c1_lo[mask_lo])
            
            mask_nano = (x_nano >= fbins[ib]) & (x_nano < fbins[ib + 1])
            if np.sum(mask_nano) > 10:
                nano_means[ib] = np.mean(c1_nano[mask_nano])

        # stairs plot requires matplotlib >= 3.4.0
        # draw as step lines
        if hasattr(ax_b, 'stairs'):
            if not is_extra:
                ax_b.stairs(lo_means, fbins, color='C0', linewidth=1.5, label='Calculated')
            ax_b.stairs(nano_means, fbins, color='C1', linewidth=1.5, linestyle='--', label='Predicted')
        else:
            # Fallback for older matplotlib
            if not is_extra:
                ax_b.step(fbins[:-1], lo_means, where='post', color='C0', linewidth=1.5, label='Calculated')
            ax_b.step(fbins[:-1], nano_means, where='post', color='C1', linewidth=1.5, linestyle='--', label='Predicted')

        ax_b.set_ylabel(r'$C_1$ [\%]')
        ax_b.set_xlabel(label)
        ax_b.set_ylim(0, 3)
        ax_b.legend(frameon=False)
        ax_b.minorticks_on()

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.06)
        out = os.path.join(plotdir, f'combined_validation_{feat}.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  [✓] {out}")

    print(f"  All validation plots saved to {plotdir}/")


def make_kappa_plots(nano_data, plotdir, plot_mode='shape'):
    """
    Plot kinematic distributions weighted by R(kappa) for various kappas.
    plot_mode='shape': density-normalized top panel, density ratio bottom panel
    plot_mode='yield': raw weighted yield top panel, yield ratio bottom panel
    Both panels always use the SAME normalization consistently.
    """
    print(f"\n── Generating Kappa {plot_mode} plots ──")
    outdir = os.path.join(plotdir, f'kappa_{plot_mode}')
    os.makedirs(outdir, exist_ok=True)

    target_kappas = [0, 10, 1, 2, 5, -10, -2, -5]
    
    cmap = plt.get_cmap('tab10')
    color_map = {k: cmap(i%10) for i, k in enumerate(target_kappas)}
    color_map[1] = 'black'

    is_density = (plot_mode == 'shape')

    for feat in FEATURES + EXTRA_NANO_FEATURES:
        label = FEATURE_LABELS.get(feat, feat)
        x_nano = nano_data[feat]

        fig, (ax_t, ax_b) = plt.subplots(
            2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(8, 6))

        if feat == 'cos_theta_star':
            fbins = np.linspace(-1, 1, 41)
        elif feat == 'h_y':
            fbins = np.linspace(-4, 4, 41)
        elif feat == 'vh_delta_eta':
            fbins = np.linspace(0, 6, 41)
        elif feat == 'n_lhe_extra':
            fbins = np.linspace(-0.5, 5.5, 7)
        elif feat == 'vh_m':
            fbins = np.linspace(200, 800, 31)
        elif feat in ['h_pt', 'v_pt', 'vh_pt_gen']:
            fbins = np.linspace(0, 300, 31)
        else:
            q95 = np.percentile(x_nano, 95)
            fbins = np.linspace(0, q95, 41)

        # SM baseline: always use the normalized weight (which equals 1.0 for k=1)
        w_1 = nano_data['weight_kappa_1']
        
        # Compute SM histogram with the chosen normalization
        nominal_hist, _ = np.histogram(x_nano, bins=fbins, weights=w_1, density=is_density)
        sum_w1 = np.sum(w_1)
        
        ax_t.hist(x_nano, bins=fbins, weights=w_1, histtype='step', linewidth=2.0, 
                  color='black', label=r'SM ($\kappa_\lambda=1$)', density=is_density)
        
        for k in target_kappas:
            if k == 1:
                continue
            
            w_k = nano_data[f'weight_kappa_{k}']
            c = color_map[k]
            k_label = r'$\kappa_\lambda = ' + str(k) + '$'
            
            # Compute kappa histogram with the SAME normalization
            k_hist, _ = np.histogram(x_nano, bins=fbins, weights=w_k, density=is_density)
            
            # Draw on top panel
            if hasattr(ax_t, 'stairs'):
                ax_t.stairs(k_hist, fbins, color=c, linewidth=1.2, label=k_label)
            else:
                ax_t.step(fbins[:-1], k_hist, where='post', color=c, linewidth=1.2, label=k_label)
            
            # Ratio: kappa / SM (both computed with same normalization)
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.where(nominal_hist > 0, k_hist / nominal_hist, np.nan)
            
            if hasattr(ax_b, 'stairs'):
                ax_b.stairs(ratio, fbins, color=c, linewidth=1.5)
            else:
                ax_b.step(fbins[:-1], ratio, where='post', color=c, linewidth=1.5)
            
            # In yield mode, show inclusive cross-section ratio as dashed line
            if plot_mode == 'yield':
                incl_ratio = np.sum(w_k) / sum_w1 if sum_w1 > 0 else 1.0
                ax_b.axhline(incl_ratio, color=c, linestyle='--', linewidth=1.0, alpha=0.7)

        # Bottom panel
        ax_b.axhline(1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
        ax_b.set_xlabel(label)
        if plot_mode == 'shape':
            ax_b.set_ylabel(r'Shape $\kappa_\lambda$/SM')
            ax_b.set_ylim(0.9, 1.1)
        else:
            ax_b.set_ylabel(r'Yield $\kappa_\lambda$/SM')
            ax_b.set_ylim(0.65, 1.20)
        ax_b.set_xlim(fbins[0], fbins[-1])
        ax_b.minorticks_on()

        # Top panel
        ax_t.set_ylabel('Density' if is_density else 'Events')
        ax_t.legend(frameon=False, ncol=2, fontsize=10, loc='upper right')
        ax_t.minorticks_on()
        
        ax_t.text(0.0, 1.02, r'\textbf{CMS} \textit{Simulation}', 
                  transform=ax_t.transAxes, fontsize=14, va='bottom', ha='left')
        ax_t.text(1.0, 1.02, r'ZH (13.6 TeV)', 
                  transform=ax_t.transAxes, fontsize=12, va='bottom', ha='right')

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.06)
        out = os.path.join(outdir, f'kappa_{plot_mode}_{feat}.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  [✓] {out}")


def main():
    p = argparse.ArgumentParser(
        description='Predict per-event C1 (NLO EW correction) on NanoAOD events using LHE-level info',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('--input', nargs='+', required=True,
                   help='Input NanoAOD ROOT file(s), supports glob patterns')
    p.add_argument('--model', default='output/c1_regressor/c1_regressor.json',
                   help='Path to trained XGBoost model (.json)')
    p.add_argument('--output', default='output/nano_c1_predictions.root',
                   help='Output ROOT file with predictions')
    p.add_argument('--ebeam', type=float, default=6800.0,
                   help='Beam energy per proton [GeV]')
    p.add_argument('--max-events', type=int, default=0,
                   help='Max events to process (0 = all)')
    p.add_argument('--validate', action='store_true',
                   help='Generate validation plots comparing LO truth vs NanoAOD predicted C1')
    p.add_argument('--lo-file', default='output/events_lo.root',
                   help='LO ROOT file for validation (used with --validate)')
    p.add_argument('--rw-file', default='output/events_rwgt.root',
                   help='Reweighted ROOT file for validation (used with --validate)')
    p.add_argument('--plotdir', default='output/plots/nano_validation',
                   help='Directory for validation plots')
    args = p.parse_args()

    # Expand glob patterns
    input_files = []
    for pattern in args.input:
        expanded = sorted(glob.glob(pattern))
        if not expanded:
            print(f"Warning: no files match '{pattern}'")
        input_files.extend(expanded)

    if not input_files:
        print("Error: no input files found")
        sys.exit(1)

    if not os.path.isfile(args.model):
        print(f"Error: model file not found: {args.model}")
        sys.exit(1)

    output_path, nano_data = process_nanoaod(
        input_files, args.model, args.output, args.ebeam, args.max_events
    )

    if args.validate:
        if not os.path.isfile(args.lo_file):
            print(f"Error: LO file not found: {args.lo_file}")
            sys.exit(1)
        if not os.path.isfile(args.rw_file):
            print(f"Error: reweighted file not found: {args.rw_file}")
            sys.exit(1)
        make_validation_plots(nano_data, args.lo_file, args.rw_file, args.plotdir)
        make_kappa_plots(nano_data, args.plotdir, plot_mode='shape')
        make_kappa_plots(nano_data, args.plotdir, plot_mode='yield')


if __name__ == '__main__':
    main()
