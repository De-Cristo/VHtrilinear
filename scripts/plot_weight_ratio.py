#!/usr/bin/env python3
"""Plot histogram of weight ratios between two ROOT files.

Usage:
  python3 scripts/plot_weight_ratio.py fileA.root fileB.root --out weight_ratio.png

The script matches events by `event_id` branch and plots the distribution of
weight_B / weight_A for each matched event. File A is used as the denominator.
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
    """Load event_id and weight from ROOT file."""
    t = uproot.open(rootfile)["events"]
    event_id = t["event_id"].array(library="np")
    weight = t["weight"].array(library="np")
    return event_id, weight


def main():
    p = argparse.ArgumentParser(description='Plot histogram of weight ratios')
    p.add_argument('fileA', help='First ROOT file (denominator)')
    p.add_argument('fileB', help='Second ROOT file (numerator)')
    p.add_argument('--out', default='weight_ratio.png', help='Output filename')
    p.add_argument('--bins', type=int, default=50, help='Number of bins')
    args = p.parse_args()

    # Load data
    print(f'Loading {args.fileA}...')
    idA, wA = load_tree(args.fileA)
    print(f'Loading {args.fileB}...')
    idB, wB = load_tree(args.fileB)

    # Build dicts for fast lookup
    mapA = {int(i): float(w) for i, w in zip(idA, wA)}
    mapB = {int(i): float(w) for i, w in zip(idB, wB)}

    # Find common events
    common = sorted(set(mapA.keys()).intersection(mapB.keys()))
    if not common:
        print('ERROR: No common events found (by event_id)')
        return 1

    print(f'Number of common events: {len(common)}')

    # Compute weight ratios: weight_B / weight_A
    weights_A = np.array([mapA[i] for i in common])
    weights_B = np.array([mapB[i] for i in common])
    
    # Filter out events where weight_A is zero to avoid division by zero
    nonzero_mask = (weights_A != 0)
    if not np.any(nonzero_mask):
        print('ERROR: All weights in file A are zero!')
        return 1
    
    weights_A_nz = weights_A[nonzero_mask]
    weights_B_nz = weights_B[nonzero_mask]
    # compute ratio and convert to percent (times 100)
    weight_ratios = (weights_B_nz / weights_A_nz) * 100.0

    # Statistics
    print(f'\nWeight ratio statistics:')
    print(f'  Mean (%):   {np.mean(weight_ratios):.4f}')
    print(f'  Median (%): {np.median(weight_ratios):.4f}')
    print(f'  Std (%):    {np.std(weight_ratios):.4f}')
    print(f'  Min (%):    {np.min(weight_ratios):.4f}')
    print(f'  Max (%):    {np.max(weight_ratios):.4f}')

    # Helper to add CMS-style labels
    def add_cms_labels(ax):
        ax.text(0.05, 0.95, r'\textbf{CMS}', fontsize=13, 
                transform=ax.transAxes, va='top')
        ax.text(0.17, 0.95, r'\textit{Simulation}', fontsize=12, 
                transform=ax.transAxes, va='top')
        ax.text(0.95, 0.95, r'Run 3 (13.6 TeV)', fontsize=11, 
                transform=ax.transAxes, ha='right', va='top')

    # Create histogram
    fig, ax = plt.subplots(figsize=(7, 6))
    
    # Plot histogram as filled steps
    counts, bins, patches = ax.hist(weight_ratios, bins=args.bins, 
                                     histtype='stepfilled', linewidth=1.5,
                                     edgecolor='C0', facecolor='C0', alpha=0.3,
                                     label=f'Entries: {len(weight_ratios)}')
    
    # Add outline
    ax.hist(weight_ratios, bins=bins, histtype='step', linewidth=1.5,
            edgecolor='C0')
    
    # Add vertical line at mean (percent)
    mean_ratio = np.mean(weight_ratios)
    # escape percent sign for LaTeX rendering
    ax.axvline(mean_ratio, color='green', linestyle='-.', linewidth=1.5,
               label=f'Mean = {mean_ratio:.3f} \\%')

    # Theoretical reference line at 1.19%
    theory_val = 1.19
    ax.axvline(theory_val, color='purple', linestyle='--', linewidth=1.5,
               label='Theory (13TeV) = 1.19 \\%')
    
    # Labels and styling
    # escape percent sign so it renders when text.usetex=True
    ax.set_xlabel(r'C1 [\%]', fontsize=12)
    ax.set_ylabel(r'Events', fontsize=12)
    ax.legend(frameon=True, loc='upper right', framealpha=1.0, edgecolor='black')
    ax.minorticks_on()
    ax.grid(False)
    
    # Style spines
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    
    ax.tick_params(direction='in', which='both', top=True, right=True,
                   length=5, width=1.0)
    ax.tick_params(direction='in', which='minor', top=True, right=True,
                   length=3, width=0.8)
    
    plt.tight_layout()
    plt.subplots_adjust(left=0.14, right=0.96, top=0.93, bottom=0.11)
    add_cms_labels(ax)
    
    plt.savefig(args.out, dpi=300, bbox_inches='tight')
    print(f'\nWrote histogram to {args.out}')
    
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
