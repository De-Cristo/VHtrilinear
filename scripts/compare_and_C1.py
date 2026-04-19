#!/usr/bin/env python3
"""Two-panel plot: top = direct pT hist comparison (density); bottom = C1[%] vs pT (differential).

Usage:
  python3 scripts/compare_and_C1.py fileA.root fileB.root --out compare_and_C1.png

File A is denominator, File B is numerator.
"""
import argparse
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica'],
    'axes.linewidth': 1.2,
    'font.size': 12,
})
mpl.rcParams['text.usetex'] = True

# Manual override for bottom-panel (C1) y-range. Set to a float to force that limit, or None to keep automatic.
# Example: C1_YMIN_MANUAL = 0.0; C1_YMAX_MANUAL = 5.0
C1_YMIN_MANUAL = -0.3
C1_YMAX_MANUAL = 2.5


def load_tree(rootfile, varname='h_pt'):
    t = uproot.open(rootfile)["events"]
    eid = t["event_id"].array(library="np")
    w = t["weight"].array(library="np")
    # variable (e.g. 'h_pt' or 'z_pt')
    if varname in t.keys():
        pt = t[varname].array(library="np")
    else:
        raise KeyError(f"Variable '{varname}' not found in {rootfile}")
    return eid, w, pt


def process_and_plot(fileA, fileB, varname, outname, args):
    """Load variable `varname` from both files and produce the two-panel plot saved to outname."""
    try:
        idA, wA, ptA = load_tree(fileA, varname=varname)
    except KeyError as e:
        print(f"Skipping {varname}: {e}")
        return 1
    try:
        idB, wB, ptB = load_tree(fileB, varname=varname)
    except KeyError as e:
        print(f"Skipping {varname}: {e}")
        return 1

    # build maps
    mapA_w = {int(i): float(w) for i, w in zip(idA, wA)}
    mapB_w = {int(i): float(w) for i, w in zip(idB, wB)}
    mapA_pt = {int(i): float(p) for i, p in zip(idA, ptA)}
    mapB_pt = {int(i): float(p) for i, p in zip(idB, ptB)}

    # common events that have pT in A and weights in both
    common = sorted(set(mapA_w.keys()).intersection(mapB_w.keys()).intersection(mapA_pt.keys()))
    if not common:
        print(f'No common events found for variable {varname}')
        return 1

    # gather arrays
    pts = np.array([mapA_pt[i] for i in common])
    weightsA = np.array([mapA_w[i] for i in common])
    weightsB = np.array([mapB_w[i] for i in common])

    # filter out zero denominator weights
    nz = (weightsA != 0)
    pts = pts[nz]
    weightsA = weightsA[nz]
    weightsB = weightsB[nz]

    # compute ratio in percent
    ratios = (weightsB / weightsA) * 100.0

    # bins (use same for hist top and C1 bottom)
    # handle special ranges for different variables
    if varname == 'cos_theta_star':
        bins = np.linspace(-1.0, 1.0, args.nbins + 1)
    elif varname == 'h_y':
        bins = np.linspace(-3.0, 3.0, args.nbins + 1)
    else:
        bins = np.linspace(0.0, args.pt_max, args.nbins + 1)
    bin_centers = 0.5*(bins[:-1] + bins[1:])

    # Top: density histograms for pT (use pts from each file separately for comparison)
    ptsA = np.array([mapA_pt[i] for i in set(mapA_pt.keys()).intersection(set(mapA_w.keys()))])
    if mapB_pt:
        ptsB = np.array([mapB_pt[i] for i in set(mapB_pt.keys()).intersection(set(mapB_w.keys()))])
    else:
        ptsB = pts

    # Compute density histograms for top panel
    try:
        wA_all = np.array([mapA_w[int(i)] for i in set(mapA_pt.keys()).intersection(set(mapA_w.keys()))])
    except Exception:
        wA_all = None
    try:
        wB_all = np.array([mapB_w[int(i)] for i in set(mapB_pt.keys()).intersection(set(mapB_w.keys()))]) if mapB_pt else None
    except Exception:
        wB_all = None

    histA, _ = np.histogram(ptsA, bins=bins, density=True, weights=wA_all)
    histB, _ = np.histogram(ptsB, bins=bins, density=True, weights=wB_all)

    # Bottom: mean ratio per bin
    mean_per_bin = np.full(len(bin_centers), np.nan)
    for ib in range(len(bin_centers)):
        lo = bins[ib]
        hi = bins[ib+1]
        # include right edge on last bin
        if ib == len(bin_centers) - 1:
            mask = (pts >= lo) & (pts <= hi)
        else:
            mask = (pts >= lo) & (pts < hi)
        if np.any(mask):
            # use weightsA to compute weighted average of ratios
            w_mask = weightsA[mask]
            r_mask = ratios[mask]
            # guard against all-zero weights
            if np.sum(w_mask) == 0:
                mean_per_bin[ib] = np.nan
            else:
                mean_per_bin[ib] = np.sum(r_mask * w_mask) / np.sum(w_mask)

    inclusive_mean = np.nanmean(ratios)
    theory_val = 1.19  # percent

    # Plotting
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2.5,1]}, figsize=(8,6))

    # Top: filled step histograms (density)
    ax_top.step(bins[:-1], histA, where='post', color='C0', linewidth=1.5, label='LO QCD')
    ax_top.step(bins[:-1], histB, where='post', color='C1', linewidth=1.5, label=rf'$O(\lambda_{{3}})$')
    # ax_top.fill_between(bin_centers, histA, step='mid', alpha=0.15, color='C0')
    # ax_top.fill_between(bin_centers, histB, step='mid', alpha=0.15, color='C1')
    ax_top.set_ylabel('Density', fontsize=12)
    ax_top.legend(frameon=False, loc='upper right')
    ax_top.set_xlim(bins[0], bins[-1])
    ax_top.set_ylim(bottom=0)
    ax_top.tick_params(direction='in', top=True, right=True)
    for spine in ax_top.spines.values():
        spine.set_linewidth(1.2)

    # Bottom: differential C1 vs pT
    # Plot as green step
    ax_bot.step(bins[:-1], mean_per_bin, where='post', color='green', linewidth=1.6, label='Differential')
    ax_bot.hlines(inclusive_mean, bins[0], bins[-1], color='blue', linewidth=1.8, label='Inclusive (code)')
    ax_bot.hlines(theory_val, bins[0], bins[-1], color='red', linewidth=1.8, linestyle='--', label='Inclusive (theory)')
    # ax_bot.legend(frameon=False, loc='upper left', bbox_to_anchor=(0.02, 0.98))

    ax_bot.set_ylabel(r'C1 [\%]', fontsize=12)
    # use clear math labels for different variables
    xlabel_map = {
        'h_pt': r'$p_{T}(H)$ [GeV]',
        'v_pt': r'$p_{T}(V)$ [GeV]',
        'vh_m': r'$m(VH)$ [GeV]',
        'h_y': r'$y(H)$',
        'vh_delta_eta': r'$\Delta\eta(V,H)$',
        'cos_theta_star': r'$\cos\theta^{*}$',
    }
    ax_bot.set_xlabel(xlabel_map.get(varname, varname), fontsize=12)
    ax_bot.set_xlim(bins[0], bins[-1])
    # set y limits with some margin
    # bottom-panel y-limits: allow manual override via constants at top of file
    if C1_YMIN_MANUAL is not None:
        ymin = float(C1_YMIN_MANUAL)
    else:
        ymin = 0.0

    if C1_YMAX_MANUAL is not None:
        ymax = float(C1_YMAX_MANUAL)
    else:
        vals = []
        if np.any(~np.isnan(mean_per_bin)):
            vals.append(np.nanmax(mean_per_bin[~np.isnan(mean_per_bin)]))
        vals.extend([inclusive_mean, theory_val])
        ymax = max(vals) if vals else theory_val * 1.5

    ax_bot.set_ylim(ymin, ymax * 1.2)

    # draw a light gray dashed horizontal line at C1 = 0 for reference
    ax_bot.axhline(0, color='lightgray', linestyle='--', linewidth=1, zorder=0)

    ax_bot.tick_params(direction='in', top=True, right=True)
    # smaller legend, shifted slightly higher to avoid overlapping the lines
    ax_bot.legend(frameon=False, loc='upper right', bbox_to_anchor=(0.95, 1.02), fontsize=10, handlelength=1.5, labelspacing=0.15, handletextpad=0.4)
    for spine in ax_bot.spines.values():
        spine.set_linewidth(1.2)

    # CMS labels on top-left of figure (use fig coordinates)
    fig.text(0.12, 0.96, r'\textbf{CMS}', fontsize=14, va='top')
    fig.text(0.22, 0.96, r'\textit{Simulation}', fontsize=11, va='top')
    fig.text(0.82, 0.96, r'{ZH (13.6 TeV)}', fontsize=11, va='top')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.08, left=0.10, right=0.96, top=0.92, bottom=0.10)
    plt.savefig(outname, dpi=300, bbox_inches='tight')
    print(f'Wrote {outname}')
    plt.close(fig)
    return 0


def main():
    p = argparse.ArgumentParser()
    p.add_argument('fileA')
    p.add_argument('fileB')
    p.add_argument('--out', default='compare_and_C1.png')
    p.add_argument('--pt-max', type=float, default=300.0)
    p.add_argument('--nbins', type=int, default=30)
    args = p.parse_args()
    
    # list of variables to plot: (varname, output_suffix, pt_max_override, nbins_override)
    variables = [
        ('h_pt', '', args.pt_max, args.nbins),
        ('v_pt', '_vpt', 600.0, args.nbins),
        ('vh_m', '_vhm', 1500.0, 50),
        ('h_y', '_hy', 6.0, 40),
        ('vh_delta_eta', '_vhdeta', 8.0, 40),
        ('cos_theta_star', '_costheta', 2.0, 40),
    ]
    
    for varname, suffix, pt_max_var, nbins_var in variables:
        # construct output name
        if suffix:
            out_var = args.out
            if out_var.lower().endswith('.png'):
                out_var = out_var[:-4] + suffix + '.png'
            else:
                out_var = out_var + suffix
        else:
            out_var = args.out
        
        # override pt_max and nbins for this variable
        args_copy = argparse.Namespace(**vars(args))
        args_copy.pt_max = pt_max_var
        args_copy.nbins = nbins_var
        
        rc = process_and_plot(args.fileA, args.fileB, varname, out_var, args_copy)
        if rc != 0:
            print(f'Warning: failed to plot {varname}')
    
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
