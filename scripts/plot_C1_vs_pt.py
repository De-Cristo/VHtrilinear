#!/usr/bin/env python3
"""Plot C1[%] (weight ratio *100) as a function of H pT (differential) with inclusive lines.

Usage:
  python3 scripts/plot_C1_vs_pt.py fileA.root fileB.root --out C1_vs_pt.png

- File A is denominator, File B is numerator.
- Computes per-event ratio (B/A)*100, bins in pT(H), computes mean ratio per bin (differential).
- Plots differential (green step), inclusive code mean (solid blue), inclusive theory (dashed blue).
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


def load_tree(rootfile):
    t = uproot.open(rootfile)["events"]
    eid = t["event_id"].array(library="np")
    w = t["weight"].array(library="np")
    hpt = t["h_pt"].array(library="np")
    return eid, w, hpt


def main():
    p = argparse.ArgumentParser()
    p.add_argument('fileA')
    p.add_argument('fileB')
    p.add_argument('--out', default='C1_vs_pt.png')
    p.add_argument('--pt-max', type=float, default=500.0)
    p.add_argument('--nbins', type=int, default=20)
    args = p.parse_args()

    idA, wA, hptA = load_tree(args.fileA)
    idB, wB, hptB = load_tree(args.fileB)

    # build dicts
    mapA_w = {int(i): float(w) for i, w in zip(idA, wA)}
    mapB_w = {int(i): float(w) for i, w in zip(idB, wB)}
    mapA_pt = {int(i): float(p) for i, p in zip(idA, hptA)}

    # find common events
    common = sorted(set(mapA_w.keys()).intersection(mapB_w.keys()).intersection(mapA_pt.keys()))
    if not common:
        print('No common events')
        return 1

    # compute ratios and pT
    ratios = []
    pts = []
    for i in common:
        wa = mapA_w[i]
        wb = mapB_w[i]
        if wa == 0:
            continue
        ratios.append((wb/wa)*100.0)
        pts.append(mapA_pt[i])
    ratios = np.array(ratios)
    pts = np.array(pts)

    # bins
    bins = np.linspace(0.0, args.pt_max, args.nbins+1)
    bin_centers = 0.5*(bins[:-1] + bins[1:])

    # compute mean per bin
    mean_per_bin = np.zeros(len(bin_centers))
    for ib in range(len(bin_centers)):
        lo = bins[ib]
        hi = bins[ib+1]
        mask = (pts >= lo) & (pts < hi)
        if np.any(mask):
            mean_per_bin[ib] = np.mean(ratios[mask])
        else:
            mean_per_bin[ib] = np.nan

    # inclusive mean (code)
    inclusive_mean = np.nanmean(ratios)

    # plotting
    fig, ax = plt.subplots(figsize=(8,4))
    # differential as step
    ax.step(bins[:-1], mean_per_bin, where='post', color='green', linewidth=1.8, label='Differential')
    # inclusive code
    ax.hlines(inclusive_mean, bins[0], bins[-1], color='blue', linewidth=2.0, label='Inclusive (code)')
    # inclusive theory at 1.19%
    theory = 1.19
    ax.hlines(theory, bins[0], bins[-1], color='blue', linewidth=2.0, linestyle='--', label='Inclusive (theory)')

    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylim(0, max(np.nanmax(mean_per_bin)*1.2, theory*1.2))
    ax.set_xlabel(r'$p_{T}(H)$[GeV]')
    ax.set_ylabel(r'C1 [\%]')

    # CMS label
    ax.text(0.02, 0.95, r'\textbf{CMS}', transform=ax.transAxes, fontsize=14, va='top')
    ax.text(0.12, 0.95, r'\textit{Simulation}', transform=ax.transAxes, fontsize=11, va='top')

    ax.legend(frameon=False, loc='upper right')
    ax.tick_params(direction='in', top=True, right=True)
    ax.minorticks_on()
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches='tight')
    print(f'Wrote {args.out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
