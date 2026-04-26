#!/usr/bin/env python3
"""Simple LHE -> ROOT converter for ZH events.

Reads an LHE file, collects event-level kinematics for Z and H candidates,
and writes a ROOT file with arrays for momenta and weights.

This script is intentionally minimal and uses text parsing of the LHE format.
"""
import argparse
import math
from collections import defaultdict

try:
    import awkward as ak
    import uproot
    import numpy as np
except Exception:
    ak = None
    uproot = None
    np = None


def parse_lhe_events(path):
    """Yield events as dicts with particles and weights.

    Particles: list of dicts with keys pid, status, px, py, pz, e
    We also try to read event weights from <rwgt> blocks if present.
    """
    with open(path, 'r') as f:
        in_event = False
        event_lines = []
        rwgt_block = []
        in_rwgt = False
        for line in f:
            if line.lstrip().startswith('<event>'):
                in_event = True
                event_lines = []
                rwgt_block = []
                in_rwgt = False
                continue
            if line.lstrip().startswith('</event>'):
                # parse event
                ev = parse_single_event('\n'.join(event_lines), '\n'.join(rwgt_block))
                yield ev
                in_event = False
                continue
            if in_event:
                # capture rwgt block if present
                if line.lstrip().startswith('<rwgt>'):
                    in_rwgt = True
                    continue
                if line.lstrip().startswith('</rwgt>'):
                    in_rwgt = False
                    continue
                if in_rwgt:
                    rwgt_block.append(line.rstrip('\n'))
                else:
                    event_lines.append(line.rstrip('\n'))


def parse_single_event(body, rwgt_text):
    """Parse the body of a single LHE <event> block (string).

    Returns dict with particles, weights, and original event weight.
    """
    lines = [l.strip() for l in body.splitlines() if l.strip()]
    if not lines:
        return None
    header = lines[0].split()
    try:
        nup = int(header[0])
    except Exception:
        nup = None
    
    # Extract original event weight from header (3rd column: XWGTUP)
    original_weight = 1.0
    try:
        original_weight = float(header[2])
    except (IndexError, ValueError):
        pass

    particles = []
    for pl in lines[1:1+nup]:
        cols = pl.split()
        if len(cols) < 10:
            continue
        pid = int(cols[0])
        status = int(cols[1])
        # LHE columns: id, status, mother1, mother2, color1, color2, px, py, pz, e, m, ...
        try:
            px = float(cols[6])
            py = float(cols[7])
            pz = float(cols[8])
            e = float(cols[9])
        except Exception:
            # fallback: set zeros
            px = py = pz = e = 0.0
        particles.append({'pid': pid, 'status': status, 'px': px, 'py': py, 'pz': pz, 'e': e})

    # Parse reweighting block if present
    weights = {}
    if rwgt_text:
        for l in rwgt_text.splitlines():
            l = l.strip()
            if not l:
                continue
            # common format: <wgt id="someid"> value </wgt>
            if l.startswith('<wgt') and '>' in l and '</wgt>' in l:
                # extract id and value
                try:
                    idpart = l.split('id="', 1)[1].split('"', 1)[0]
                except Exception:
                    idpart = 'unknown'
                try:
                    val = float(l.split('>')[1].split('<')[0])
                except Exception:
                    val = None
                weights[idpart] = val

    return {'particles': particles, 'weights': weights, 'original_weight': original_weight}


def find_bosons(particles):
    """Look for V (Z pid 23 or W pid 24/-24) and H (id 25) among particles. Return their four-vectors if found."""
    v = None
    h = None
    for p in particles:
        if p['pid'] == 23 or abs(p['pid']) == 24:
            v = p
        if p['pid'] == 25:
            h = p
    return v, h


def fourvec(px, py, pz, e):
    return np.array([e, px, py, pz], dtype=float)


def pt(px, py):
    return math.hypot(px, py)


def eta(px, py, pz, e):
    """Calculate pseudorapidity, handling edge cases robustly."""
    pt_val = math.hypot(px, py)
    p = math.hypot(pt_val, pz)
    
    # Handle zero momentum
    if p == 0:
        return 0.0
    
    # Use the standard formula: eta = -ln(tan(theta/2))
    # where cos(theta) = pz/p
    # Equivalent: eta = 0.5 * ln((p+pz)/(p-pz))
    if abs(pz) >= p:
        # Forward/backward: return large eta with correct sign
        return 10.0 * (1 if pz > 0 else -1)
    
    # Safe calculation
    return 0.5 * math.log((p + pz) / (p - pz))


def phi(px, py):
    return math.atan2(py, px)


def boost(vec, beta):
    """Lorentz boost of 4-vector vec=(E,px,py,pz) by velocity beta (3-vector).
    Returns boosted 4-vector [E', px', py', pz'].
    """
    vec = np.array(vec, dtype=float)
    beta = np.array(beta, dtype=float)
    b2 = float(np.dot(beta, beta))
    if b2 == 0.0:
        return vec.copy()
    if b2 >= 1.0:
        # unphysical: return NaNs
        return np.array([np.nan, np.nan, np.nan, np.nan])
    gamma = 1.0 / math.sqrt(1.0 - b2)
    E = float(vec[0])
    p = vec[1:]
    bp = float(np.dot(beta, p))
    Ep = gamma * (E - bp)
    p_star = p + beta * ( (gamma - 1.0) * bp / b2 - gamma * E )
    return np.array([Ep, p_star[0], p_star[1], p_star[2]])


def run(infile, outfile, ebeam=6800.0):
    if ak is None or uproot is None or np is None:
        print('Missing dependencies: please pip install awkward uproot numpy')
        return 1

    event_ids = []
    v_m = []
    v_pt = []
    v_eta = []
    v_phi = []
    h_m = []
    h_pt = []
    h_eta = []
    h_phi = []
    h_y = []
    vh_m = []
    vh_delta_eta = []
    cos_theta_star = []
    all_weights = []

    event_counter = 0
    for ev in parse_lhe_events(infile):
        if ev is None:
            continue
        v, h = find_bosons(ev['particles'])
        if v is None or h is None:
            # skip events without both, but still increment counter
            event_counter += 1
            continue

        # Store the event index (0-based indexing)
        event_ids.append(event_counter)
        event_counter += 1

        # V (Z or W) kinematics
        v_px, v_py, v_pz, v_e = v['px'], v['py'], v['pz'], v['e']
        v_m2 = v_e*v_e - (v_px*v_px + v_py*v_py + v_pz*v_pz)
        v_m.append(math.sqrt(v_m2) if v_m2 > 0 else 0.0)
        v_pt.append(pt(v_px, v_py))
        v_eta_val = eta(v_px, v_py, v_pz, v_e)
        v_eta.append(v_eta_val)
        v_phi.append(phi(v_px, v_py))

        # H kinematics
        h_px, h_py, h_pz, h_e = h['px'], h['py'], h['pz'], h['e']
        h_m2 = h_e*h_e - (h_px*h_px + h_py*h_py + h_pz*h_pz)
        h_m.append(math.sqrt(h_m2) if h_m2 > 0 else 0.0)
        h_pt.append(pt(h_px, h_py))
        h_eta_val = eta(h_px, h_py, h_pz, h_e)
        h_eta.append(h_eta_val)
        h_phi.append(phi(h_px, h_py))

        # Higgs rapidity y = 0.5 * ln((E + pz)/(E - pz)); handle edge cases
        try:
            denom = h_e - h_pz
            if denom <= 0 or not np.isfinite(h_e + h_pz) or not np.isfinite(denom):
                h_y_val = 10.0 * (1 if h_pz >= 0 else -1)
            else:
                h_y_val = 0.5 * math.log((h_e + h_pz) / denom)
        except Exception:
            h_y_val = 0.0
        h_y.append(h_y_val)

        # delta eta between V and H: |eta_V - eta_H|
        vh_delta_eta.append(abs(v_eta_val - h_eta_val))

        # V+H system mass
        vh_e = v_e + h_e
        vh_px = v_px + h_px
        vh_py = v_py + h_py
        vh_pz = v_pz + h_pz
        vh_m2 = vh_e*vh_e - (vh_px*vh_px + vh_py*vh_py + vh_pz*vh_pz)
        vh_m.append(math.sqrt(vh_m2) if vh_m2 > 0 else 0.0)

        # --- scattering angle in VH rest frame ---
        # build 4-vectors
        pH = np.array([h_e, h_px, h_py, h_pz], dtype=float)
        pV = np.array([v_e, v_px, v_py, v_pz], dtype=float)
        q = pH + pV
        beta = q[1:]/q[0] if q[0] != 0 else np.array([0.0, 0.0, 0.0])

        # boost H and beam to VH rest frame
        pH_star = boost(pH, beta)
        P1 = np.array([ebeam, 0.0, 0.0, ebeam], dtype=float)
        P1_star = boost(P1, beta)

        # Option A: beam-axis angle (cos_theta)
        pH_spat = pH_star[1:]
        norm_pH = np.linalg.norm(pH_spat)
        z_beam_vec = P1_star[1:]
        norm_zbeam = np.linalg.norm(z_beam_vec)
        if norm_pH > 0 and norm_zbeam > 0:
            z_beam = z_beam_vec / norm_zbeam
            cos_th = float(np.dot(pH_spat, z_beam) / norm_pH)
            cos_th = max(-1.0, min(1.0, cos_th))
        else:
            cos_th = np.nan
        cos_theta_star.append(cos_th)

        # Use original event weight from header; fallback to rwgt if available
        if ev['original_weight'] != 1.0:
            w = ev['original_weight']
        elif ev['weights']:
            # pick the first non-null weight from rwgt block
            w = next((v for v in ev['weights'].values() if v is not None), 1.0)
        else:
            w = 1.0
        all_weights.append(w)

    # convert to awkward arrays and write to ROOT
    out = {}
    out['event_id'] = np.array(event_ids, dtype=np.int32)
    out['v_m'] = np.array(v_m)
    out['v_pt'] = np.array(v_pt)
    out['v_eta'] = np.array(v_eta)
    out['v_phi'] = np.array(v_phi)
    out['h_m'] = np.array(h_m)
    out['h_pt'] = np.array(h_pt)
    out['h_eta'] = np.array(h_eta)
    out['h_phi'] = np.array(h_phi)
    out['h_y'] = np.array(h_y)
    out['vh_m'] = np.array(vh_m)
    out['vh_delta_eta'] = np.array(vh_delta_eta)
    out['cos_theta_star'] = np.array(cos_theta_star)
    out['weight'] = np.array(all_weights)

    with uproot.recreate(outfile) as f:
        f['events'] = out

    print(f'Wrote {len(v_m)} events to {outfile}')
    return 0


def main():
    p = argparse.ArgumentParser(description='Convert LHE to ROOT for ZH analysis')
    p.add_argument('infile')
    p.add_argument('outfile')
    p.add_argument('--ebeam', type=float, default=6800.0, help='Beam energy per proton [GeV] used to define beam 4-vectors for angle calculation')
    args = p.parse_args()
    return run(args.infile, args.outfile, ebeam=args.ebeam)


if __name__ == '__main__':
    raise SystemExit(main())
