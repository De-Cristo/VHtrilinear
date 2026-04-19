#!/usr/bin/env python3
"""Compare two ROOT files event-by-event using event_id and plot h_pt/z_pt scatter.

Usage:
  python3 scripts/compare_roots.py fileA.root fileB.root --out compare.png

The script matches events by `event_id` branch. It will plot scatter of
(h_pt from A, h_pt from B) and similarly for Z. Points on the diagonal indicate
identical values (within float precision) and matching event ordering.
"""
import argparse
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib as mpl

# Matplotlib styling tuned for HEP paper (PRL) aesthetics
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica'],
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'sans',
    'mathtext.it': 'sans:italic',
    'mathtext.bf': 'sans:bold',
    'axes.linewidth': 1.2,
    'font.size': 12,
    'axes.titlesize': 12,
    'axes.labelsize': 12,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'xtick.major.size': 5,
    'xtick.minor.size': 3,
    'ytick.major.size': 5,
    'ytick.minor.size': 3,
})

mpl.rcParams['text.usetex'] = True

def load_tree(rootfile):
    t = uproot.open(rootfile)["events"]
    # read arrays
    event_id = t["event_id"].array(library="np")
    h_pt = t["h_pt"].array(library="np")
    z_pt = t["z_pt"].array(library="np")
    return event_id, h_pt, z_pt


def main():
    p = argparse.ArgumentParser(description='Compare two ROOT files event-by-event')
    p.add_argument('fileA')
    p.add_argument('fileB')
    args = p.parse_args()

    idA, hA, zA = load_tree(args.fileA)
    idB, hB, zB = load_tree(args.fileB)

    # build dicts for fast lookup
    mapA_h = {int(i): float(v) for i, v in zip(idA, hA)}
    mapB_h = {int(i): float(v) for i, v in zip(idB, hB)}
    mapA_z = {int(i): float(v) for i, v in zip(idA, zA)}
    mapB_z = {int(i): float(v) for i, v in zip(idB, zB)}

    # intersection of event ids
    common = sorted(set(mapA_h.keys()).intersection(mapB_h.keys()))
    if not common:
        print('No common events found (by event_id)')
        return 1

    h_vals_A = np.array([mapA_h[i] for i in common])
    h_vals_B = np.array([mapB_h[i] for i in common])
    z_vals_A = np.array([mapA_z[i] for i in common])
    z_vals_B = np.array([mapB_z[i] for i in common])

    # stats
    print(f'Number of common events: {len(common)}')
    print(f'h_pt differences: mean={np.mean(h_vals_A-h_vals_B):.3g}, max={np.max(np.abs(h_vals_A-h_vals_B)):.3g}')
    print(f'z_pt differences: mean={np.mean(z_vals_A-z_vals_B):.3g}, max={np.max(np.abs(z_vals_A-z_vals_B)):.3g}')

    # helper to add CMS-style labels
    def add_cms_labels(fig):
        fig.text(0.02, 0.96, r'\textbf{CMS}', fontsize=13, va='top')
        fig.text(0.14, 0.96, r'\textit{Simulation}', fontsize=12, va='top')
        fig.text(0.98, 0.96, r'Run 3 (13.6 TeV)', fontsize=11, ha='right', va='top')


    # scatter plot
    fig, axes = plt.subplots(1,2, figsize=(12,5.5))
    
    # H pT scatter
    axes[0].scatter(h_vals_A, h_vals_B, s=12, alpha=0.5, edgecolors='none', color='steelblue')
    mn = min(h_vals_A.min(), h_vals_B.min())
    mx = max(h_vals_A.max(), h_vals_B.max())
    axes[0].plot([mn,mx],[mn,mx],'r-', linewidth=1.5, label='y=x')
    axes[0].set_xlabel(r'$p_{\mathrm{T}}^{H}$ (File A) [GeV]', fontsize=12)
    axes[0].set_ylabel(r'$p_{\mathrm{T}}^{H}$ (File B) [GeV]', fontsize=12)
    axes[0].grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    axes[0].legend(frameon=False)
    axes[0].tick_params(direction='in', which='both', top=True, right=True)
    axes[0].minorticks_on()

    # Z pT scatter
    axes[1].scatter(z_vals_A, z_vals_B, s=12, alpha=0.5, edgecolors='none', color='darkorange')
    mn = min(z_vals_A.min(), z_vals_B.min())
    mx = max(z_vals_A.max(), z_vals_B.max())
    axes[1].plot([mn,mx],[mn,mx],'r-', linewidth=1.5, label='y=x')
    axes[1].set_xlabel(r'$p_{\mathrm{T}}^{Z}$ (File A) [GeV]', fontsize=12)
    axes[1].set_ylabel(r'$p_{\mathrm{T}}^{Z}$ (File B) [GeV]', fontsize=12)
    axes[1].grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    axes[1].legend(frameon=False)
    axes[1].tick_params(direction='in', which='both', top=True, right=True)
    axes[1].minorticks_on()

    plt.tight_layout()
    plt.subplots_adjust(left=0.12, right=0.96, top=0.93, bottom=0.10)
    add_cms_labels(fig)
    scatter_out = 'compare_scatter.png'
    plt.savefig(scatter_out, dpi=300, bbox_inches='tight')
    print(f'Wrote scatter plot to {scatter_out}')

    # histogram + ratio plots
    def hist_ratio(valsA, valsB, name, outname):
        # consistent binning across files with fewer bins for cleaner plots
        bins = np.histogram_bin_edges(np.concatenate([valsA, valsB]), bins=25)
        hA, _ = np.histogram(valsA, bins=bins)
        hB, _ = np.histogram(valsB, bins=bins)

        # Normalize to get density (area = 1)
        bin_width = bins[1] - bins[0]
        densA = hA / (hA.sum() * bin_width)
        densB = hB / (hB.sum() * bin_width)

        # compute ratio
        ratio = np.full_like(hA, np.nan, dtype=float)
        mask = (densB > 0)
        ratio[mask] = densA[mask] / densB[mask]

        # figure layout similar to HEP papers: main panel + ratio panel
        fig, (ax_top, ax_bot) = plt.subplots(2,1,sharex=True, gridspec_kw={"height_ratios": [3,1]}, figsize=(6,7))
        
        # Define colors for the two files
        colors = ['C0', 'C1']  # blue and orange
        
        # top: filled step histograms
        ax_top.hist(valsA, bins=bins, histtype='step', linewidth=1.5, 
                   color=colors[0], label='File A', density=True, fill=False)
        ax_top.hist(valsB, bins=bins, histtype='step', linewidth=1.5, 
                   color=colors[1], label='File B', density=True, fill=False)
        
        # Add subtle filling
        ax_top.hist(valsA, bins=bins, histtype='stepfilled', linewidth=0, 
                   color=colors[0], alpha=0.15, density=True)
        ax_top.hist(valsB, bins=bins, histtype='stepfilled', linewidth=0, 
                   color=colors[1], alpha=0.15, density=True)
        
        ax_top.set_ylabel('Density', fontsize=12)
        ax_top.legend(frameon=True, loc='upper right', framealpha=1.0, edgecolor='black')
        ax_top.set_yscale('linear')
        ax_top.grid(False)
        
        # add minor ticks
        ax_top.minorticks_on()

        # style ticks and spines (all visible for PRL style)
        for ax in (ax_top, ax_bot):
            ax.tick_params(direction='in', which='both', top=True, right=True,
                          length=5, width=1.0)
            ax.tick_params(direction='in', which='minor', top=True, right=True,
                          length=3, width=0.8)
            ax.spines['right'].set_visible(True)
            ax.spines['top'].set_visible(True)
            ax.spines['left'].set_linewidth(1.2)
            ax.spines['right'].set_linewidth(1.2)
            ax.spines['top'].set_linewidth(1.2)
            ax.spines['bottom'].set_linewidth(1.2)

        # bottom: ratio with dashed reference line
        ax_bot.axhline(1.0, color='black', linestyle='--', linewidth=1.0)
        
        # Plot ratio as step histogram
        bin_centers = 0.5*(bins[:-1] + bins[1:])
        ax_bot.step(bin_centers, ratio, where='mid', color='C0', linewidth=1.5)
        
        ax_bot.set_ylabel('Ratio', fontsize=12)
        ax_bot.set_xlabel(r'$p_{\mathrm{T}}$ [GeV]', fontsize=12)
        ax_bot.minorticks_on()
        
        # set ratio y-range around 1 for clarity
        valid = ~np.isnan(ratio)
        if np.any(valid):
            maxdev = np.nanmax(np.abs(ratio[valid] - 1.0))
            pad = max(0.15, 1.3 * maxdev)
            ax_bot.set_ylim(1.0 - pad, 1.0 + pad)
        ax_bot.grid(False)

        plt.subplots_adjust(hspace=0.05, left=0.14, right=0.96, top=0.93, bottom=0.11)
        add_cms_labels(fig)
        plt.savefig(outname, dpi=300, bbox_inches='tight')
        print(f'Wrote histogram+ratio plot to {outname}')

    hist_ratio(h_vals_A, h_vals_B, 'h_pt', 'compare_hpt_hist.png')
    hist_ratio(z_vals_A, z_vals_B, 'z_pt', 'compare_zpt_hist.png')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
