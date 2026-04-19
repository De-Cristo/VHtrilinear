#!/usr/bin/env python3
"""Plot feature comparison for a set of kappa3-weighted ROOT files sharing a common prefix.

Usage:
  python scripts/plot_kappa3.py events_500k_13p6_rwgt.root --feature h_pt --nbins 30 --pt-max 300

Behavior:
- The provided file (file_base) is treated as the central sample (denominator).
- All files in the same directory whose names start with the base name (base without .root)
  are collected (including the base file itself and variants like base_m10p0.root).
- For each file found we compute density histograms (top panel) and the per-bin mean ratio
  (weights_file / weights_base * 100) (bottom panel), plotted in the same style as
  scripts/compare_and_C1.py.

The script expects each file to have a tree named 'events' and branches 'event_id',
the feature (e.g. 'h_pt') and a weight branch named 'weight' by default.
"""
import argparse
import os
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


def find_variant_files(base_path):
    """Return sorted list of files that share the base prefix.
    base_path is a path to the central file (e.g. /path/events_500k_13p6_rwgt.root)
    """
    dname = os.path.dirname(base_path) or '.'
    fname = os.path.basename(base_path)
    if fname.lower().endswith('.root'):
        base = fname[:-5]
    else:
        base = fname
    candidates = []
    for f in os.listdir(dname):
        if not f.lower().endswith('.root'):
            continue
        if f.startswith(base):
            candidates.append(os.path.join(dname, f))
    candidates = sorted(candidates)
    return candidates


def load_tree_arrays(path, id_branch='event_id', weight_branch='weight', feature=None):
    t = uproot.open(path)['events']
    ids = t[id_branch].array(library='np')
    w = t[weight_branch].array(library='np')
    if feature is not None:
        if feature in t.keys():
            feat = t[feature].array(library='np')
        else:
            raise KeyError(f"Feature '{feature}' not found in {path}")
    else:
        feat = None
    return ids, w, feat


def parse_l3_from_name(fname, base_name):
    """Attempt to extract l3 value from filename suffix. If not possible, return filename tail.
    Example: base 'events_..._rwgt' and fname 'events_..._rwgt_m10p0.root' -> '-10.0'
    """
    bn = os.path.basename(fname)
    if bn.lower().endswith('.root'):
        bn = bn[:-5]
    if not bn.startswith(base_name + '_'):
        return os.path.basename(fname)
    suff = bn[len(base_name)+1:]
    # convert m10p0 -> -10.0, 1p0 -> 1.0
    s = suff.replace('m', '-', 1).replace('p', '.', 1)
    return s


def process_and_plot(nlo_base_file, lo_file, feature, id_branch='event_id', weight_branch='weight', nbins=30, pt_max=300.0, outname=None):
    # find all NLO variants from the NLO base file prefix
    nlo_files = find_variant_files(nlo_base_file)
    if not nlo_files:
        raise SystemExit('No files found matching NLO base')

    # derive base_name from nlo_base_file for labeling/output
    nlo_base_fname = os.path.basename(nlo_base_file)
    base_name = nlo_base_fname[:-5] if nlo_base_fname.lower().endswith('.root') else nlo_base_fname

    # ---------- Plot 1: NLO-only comparison (use NLO base as baseline) ----------
    # load base NLO arrays
    ids_base, w_base, feat_base = load_tree_arrays(nlo_base_file, id_branch=id_branch, weight_branch=weight_branch, feature=feature)

    # prepare bins
    pt_max_local = pt_max
    bins = np.linspace(0.0, pt_max_local, nbins + 1)
    bin_centers = 0.5*(bins[:-1] + bins[1:])

    fig1, (ax1_top, ax1_bot) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2.5,1]}, figsize=(8,6))
    colors = plt.rcParams['axes.prop_cycle'].by_key().get('color', None)
    if colors is None:
        colors = ['C%d' % i for i in range(10)]

    # compute base histogram from base NLO file
    hist_base_nlo, _ = np.histogram(feat_base, bins=bins, density=True, weights=w_base)
    # plot base on top
    ax1_top.step(bins[:-1], hist_base_nlo, where='post', color='black', linewidth=1.8, label='NLO base')

    for idx, fpath in enumerate(nlo_files):
        try:
            ids_f, w_f, feat_f = load_tree_arrays(fpath, id_branch=id_branch, weight_branch=weight_branch, feature=feature)
        except KeyError as e:
            print(f"Skipping {fpath}: {e}")
            continue

        # density hist for this NLO file
        hist_f, _ = np.histogram(feat_f, bins=bins, density=True, weights=w_f)

        color = colors[idx % len(colors)]
        label = parse_l3_from_name(fpath, base_name)

        if os.path.abspath(fpath) == os.path.abspath(nlo_base_file):
            ax1_top.step(bins[:-1], hist_f, where='post', color=color, linewidth=1.8, label=f'SM EW NLO')
        else:
            ax1_top.step(bins[:-1], hist_f, where='post', color=color, linewidth=1.2, alpha=0.9, label=label)

        # bottom: ratio hist_f / hist_base_nlo
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio_bins = np.where(hist_base_nlo != 0.0, hist_f / hist_base_nlo, np.nan)
        ax1_bot.step(bins[:-1], ratio_bins, where='post', color=color, linewidth=1.4 if os.path.abspath(fpath)!=os.path.abspath(nlo_base_file) else 2.0)
        # inclusive mean (binwise)
        incl = np.nanmean(ratio_bins)
        if os.path.abspath(fpath) == os.path.abspath(nlo_base_file):
            ax1_bot.hlines(1.0, bins[0], bins[-1], color='black', linewidth=1.4, linestyle='-')
        else:
            ax1_bot.hlines(incl, bins[0], bins[-1], color=color, linewidth=1.2, linestyle='--')

    ax1_top.set_ylabel('Density', fontsize=12)
    ax1_top.set_xlim(bins[0], bins[-1])
    ax1_top.set_ylim(bottom=0)
    ax1_top.tick_params(direction='in', top=True, right=True)
    ax1_top.legend(frameon=False, loc='upper right', ncol=2)

    ax1_bot.set_ylabel('BSM/SM NLO($\lambda_{3})$', fontsize=12)
    if feature == 'h_pt':
        xlabel = r'$p_{T}(H)$ [GeV]'
    elif feature == 'z_pt':
        xlabel = r'$p_{T}(Z)$ [GeV]'
    else:
        xlabel = feature
    ax1_bot.set_xlabel(xlabel, fontsize=12)
    ax1_bot.set_xlim(bins[0], bins[-1])
    ax1_bot.tick_params(direction='in', top=True, right=True)

    fig1.text(0.12, 0.96, r'\textbf{CMS}', fontsize=14, va='top')
    fig1.text(0.22, 0.96, r'\textit{Simulation}', fontsize=11, va='top')
    fig1.text(0.82, 0.96, r'{ZH (13.6 TeV)}', fontsize=11, va='top')
    plt.tight_layout()
    out1 = base_name + '_nlo_only.png' if outname is None else outname.replace('.png', '_nlo_only.png')
    plt.subplots_adjust(hspace=0.08, left=0.10, right=0.96, top=0.92, bottom=0.10)
    plt.savefig(out1, dpi=300, bbox_inches='tight')
    print(f'Wrote {out1}')
    plt.close(fig1)

    # ---------- Plot 2: NLO (file1) vs LO (file2); bottom = (NLO-LO)/LO ----------
    # load LO arrays
    ids_lo, w_lo, feat_lo = load_tree_arrays(lo_file, id_branch=id_branch, weight_branch=weight_branch, feature=feature)
    # compute LO baseline density
    hist_lo_full, _ = np.histogram(feat_lo, bins=bins, density=True, weights=w_lo)

    # load only the provided NLO file (nlo_base_file) and compare to LO
    try:
        ids_nlo, w_nlo, feat_nlo = load_tree_arrays(nlo_base_file, id_branch=id_branch, weight_branch=weight_branch, feature=feature)
    except KeyError as e:
        raise SystemExit(f"Failed to load NLO file {nlo_base_file}: {e}")

    hist_nlo, _ = np.histogram(feat_nlo, bins=bins, density=True, weights=w_nlo)

    fig2, (ax2_top, ax2_bot) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [2.5,1]}, figsize=(8,6))
    # plot LO baseline on top
    ax2_top.step(bins[:-1], hist_lo_full, where='post', color='black', linewidth=1.8, label='LO')

    # plot the single NLO curve (from nlo_base_file)
    color = colors[0] if colors else 'C0'
    ax2_top.step(bins[:-1], hist_nlo, where='post', color=color, linewidth=1.6, alpha=0.9, label='NLO')

    # bottom: relative difference (NLO - LO)/LO
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_bins = np.where(hist_lo_full != 0.0, (hist_nlo - hist_lo_full) / hist_lo_full, np.nan)
    ax2_bot.step(bins[:-1], rel_bins, where='post', color=color, linewidth=1.6)
    incl_rel = np.nanmean(rel_bins)
    # plot mean relative difference
    ax2_bot.hlines(incl_rel, bins[0], bins[-1], color=color, linewidth=1.2, linestyle='--')

    ax2_top.set_ylabel('Density', fontsize=12)
    ax2_top.set_xlim(bins[0], bins[-1])
    ax2_top.set_ylim(bottom=0)
    ax2_top.tick_params(direction='in', top=True, right=True)
    ax2_top.legend(frameon=False, loc='upper right', ncol=2)

    ax2_bot.set_ylabel('(NLO-LO)/LO', fontsize=12)
    ax2_bot.set_xlabel(xlabel, fontsize=12)
    ax2_bot.set_xlim(bins[0], bins[-1])
    ax2_bot.tick_params(direction='in', top=True, right=True)

    fig2.text(0.12, 0.96, r'\textbf{CMS}', fontsize=14, va='top')
    fig2.text(0.22, 0.96, r'\textit{Simulation}', fontsize=11, va='top')
    fig2.text(0.82, 0.96, r'{ZH (13.6 TeV)}', fontsize=11, va='top')
    plt.tight_layout()
    out2 = base_name + '_vs_lo.png' if outname is None else outname.replace('.png', '_vs_lo.png')
    plt.subplots_adjust(hspace=0.08, left=0.10, right=0.96, top=0.92, bottom=0.10)
    plt.savefig(out2, dpi=300, bbox_inches='tight')
    print(f'Wrote {out2}')
    plt.close(fig2)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('nlo_base')
    p.add_argument('lo_file')
    p.add_argument('--feature', default='h_pt')
    p.add_argument('--id-branch', default='event_id')
    p.add_argument('--weight-branch', default='weight')
    p.add_argument('--nbins', type=int, default=30)
    p.add_argument('--pt-max', type=float, default=300.0)
    args = p.parse_args()

    process_and_plot(args.nlo_base, args.lo_file, args.feature, id_branch=args.id_branch, weight_branch=args.weight_branch, nbins=args.nbins, pt_max=args.pt_max)


if __name__ == '__main__':
    raise SystemExit(main())
