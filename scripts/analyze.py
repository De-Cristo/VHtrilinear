#!/usr/bin/env python3
"""
analyze.py — Unified Stage 4 analysis pipeline for VHtrilinear.

Runs all analysis steps in sequence:
  1. LHE → ROOT conversion (events.lhe → events_lo.root, events_rwgt.root)
  2. BSM weight computation for each κ_λ value (add_l3_weight.py)
  3. All validation & plotting scripts

Usage:
  python3 scripts/analyze.py [options]

Options:
  --lo-lhe       Path to the LO LHE file            [default: output/events.lhe]
  --rw-lhe       Path to the reweighted LHE file    [default: output/events_rwgt.lhe]
  --outdir       Output directory                   [default: output]
  --ebeam        Beam energy per proton [GeV]       [default: 6800.0]
  --kappa        Comma-separated κ_λ values         [default: 0,1,2,5,10,-2,-5,-10]
  --feature      Feature for kappa3 overlay plot    [default: h_pt]
  --skip-lhe     Skip LHE→ROOT conversion (use existing ROOT files)
  --skip-weights Skip BSM weight computation
  --skip-plots   Skip plotting
"""
import argparse
import os
import sys
import glob

# ─────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))


def _call_with_argv(main_fn, argv):
    """Run a main() function with a custom sys.argv, then restore it."""
    old = sys.argv
    sys.argv = argv
    try:
        main_fn()
    except SystemExit:
        pass
    finally:
        sys.argv = old


def _step(msg):
    print(f"\n{'─'*60}")
    print(f"  {msg}")
    print('─'*60)


def _ok(msg):
    print(f"  [✓] {msg}")


def _skip(msg):
    print(f"  [–] {msg} (skipped)")


def _warn(msg):
    print(f"  [!] {msg}")


# ─────────────────────────────────────────────────────
# Step 1: LHE → ROOT
# ─────────────────────────────────────────────────────

def step_lhe_to_root(lo_lhe, rw_lhe, lo_root, rw_root, ebeam):
    # import the converter's run() directly — it exposes a proper API
    sys.path.insert(0, SCRIPTS_DIR)
    import lhe_to_root

    if not os.path.isfile(lo_lhe):
        print(f"  ERROR: LO LHE file not found: {lo_lhe}")
        sys.exit(1)
    if not os.path.isfile(rw_lhe):
        print(f"  ERROR: Reweighted LHE file not found: {rw_lhe}")
        sys.exit(1)

    print(f"  Converting: {lo_lhe} → {lo_root}")
    lhe_to_root.run(lo_lhe, lo_root, ebeam=ebeam)
    _ok(f"LO ROOT: {lo_root}")

    print(f"  Converting: {rw_lhe} → {rw_root}")
    lhe_to_root.run(rw_lhe, rw_root, ebeam=ebeam)
    _ok(f"Reweighted ROOT: {rw_root}")


# ─────────────────────────────────────────────────────
# Step 2: BSM Weights (add_l3_weight.py)
# ─────────────────────────────────────────────────────

def step_bsm_weights(lo_root, rw_root, kappas):
    sys.path.insert(0, SCRIPTS_DIR)
    import add_l3_weight

    for kl in kappas:
        print(f"  Computing weights for κ_λ = {kl} ...")
        _call_with_argv(
            add_l3_weight.main,
            ['add_l3_weight.py', lo_root, rw_root, '--l3', str(kl)]
        )
    _ok(f"BSM weights computed for κ_λ ∈ {{{', '.join(str(k) for k in kappas)}}}")


# ─────────────────────────────────────────────────────
# Step 3: Plots
# ─────────────────────────────────────────────────────

def step_plots(lo_root, rw_root, plotdir, kappas, feature):
    sys.path.insert(0, SCRIPTS_DIR)
    import plot_weight_ratio
    import plot_C1_vs_pt
    import compare_and_C1
    import compare_roots
    import plot_kappa3

    os.makedirs(plotdir, exist_ok=True)

    # 3a. Weight ratio (C1 histogram)
    out = os.path.join(plotdir, 'weight_ratio.png')
    _call_with_argv(
        plot_weight_ratio.main,
        ['plot_weight_ratio.py', lo_root, rw_root, '--out', out]
    )
    _ok(f"weight_ratio.png")

    # 3b. C1 vs pT(H)
    out = os.path.join(plotdir, 'C1_vs_pt.png')
    _call_with_argv(
        plot_C1_vs_pt.main,
        ['plot_C1_vs_pt.py', lo_root, rw_root, '--out', out]
    )
    _ok(f"C1_vs_pt.png")

    # 3c. Two-panel density + C1 for 6 variables
    out = os.path.join(plotdir, 'compare_and_C1.png')
    _call_with_argv(
        compare_and_C1.main,
        ['compare_and_C1.py', lo_root, rw_root, '--out', out]
    )
    _ok(f"compare_and_C1.png (+ variant plots)")

    # 3d. Event-by-event scatter — this script writes its output to CWD, so chdir temporarily
    orig_cwd = os.getcwd()
    os.chdir(plotdir)
    _call_with_argv(
        compare_roots.main,
        ['compare_roots.py', lo_root, rw_root]
    )
    os.chdir(orig_cwd)
    _ok(f"compare_scatter.png, compare_hpt_hist.png, compare_zpt_hist.png")

    # 3e–3g. Multi-κ_λ overlay plots (plot_kappa3.py)
    #
    # plot_kappa3.find_variant_files() discovers files sharing the same prefix.
    # add_l3_weight produces:
    #   - events_lo_l3corr.root       (NLO-corrected LO, prefix: events_lo_l3corr)
    #   - events_l3corr_<κ>.root      (full BSM, prefix: events_l3corr)
    #   - events_l3corr_nloew_<κ>.root (NLO EW only, prefix: events_l3corr_nloew)
    #
    # The NLO-corrected LO file has a DIFFERENT prefix from the BSM variants!
    # Fix: symlink events_l3corr.root → events_lo_l3corr.root so they share
    # prefix "events_l3corr" and find_variant_files discovers everything.
    rw_dir = os.path.dirname(rw_root)
    rw_base = os.path.basename(rw_root)
    rw_stem = rw_base[:-5] if rw_base.lower().endswith('.root') else rw_base
    l3corr_stem = rw_stem[:-5] + '_l3corr' if rw_stem.endswith('_rwgt') else rw_stem + '_l3corr'

    lo_stem = os.path.basename(lo_root)
    lo_stem = lo_stem[:-5] if lo_stem.lower().endswith('.root') else lo_stem
    lo_l3corr = os.path.join(os.path.dirname(lo_root), f"{lo_stem}_l3corr.root")

    # Create symlink: events_l3corr.root → events_lo_l3corr.root
    l3corr_base = os.path.join(rw_dir, f"{l3corr_stem}.root")
    if os.path.isfile(lo_l3corr) and not os.path.exists(l3corr_base):
        os.symlink(os.path.abspath(lo_l3corr), l3corr_base)
        _ok(f"Symlinked {os.path.basename(l3corr_base)} → {os.path.basename(lo_l3corr)}")

    # 3e. Full BSM overlay: uses events_l3corr.root as base, discovers events_l3corr_*.root
    if os.path.exists(l3corr_base):
        orig_cwd2 = os.getcwd()
        os.chdir(plotdir)
        # find_variant_files will match: events_l3corr.root, events_l3corr_0p0.root, etc.
        # but also events_l3corr_nloew_*.root — filter those out isn't needed because
        # the base plot shows all variants which is useful.
        plot_kappa3.process_and_plot(
            nlo_base_file=l3corr_base,
            lo_file=lo_root,
            feature=feature,
            nbins=30,
            pt_max=300.0,
            exclude_pattern='nloew'
        )
        os.chdir(orig_cwd2)
        _ok(f"kappa3 overlay: full BSM variants ({feature})")
    else:
        _warn(f"Skipping BSM kappa3 overlay — base not found: {l3corr_base}")

    # 3f. NLO EW-only overlay: need a base with prefix "events_l3corr_nloew"
    # Create symlink events_l3corr_nloew.root → events_lo_l3corr.root as baseline
    nloew_base = os.path.join(rw_dir, f"{l3corr_stem}_nloew.root")
    if os.path.isfile(lo_l3corr) and not os.path.exists(nloew_base):
        os.symlink(os.path.abspath(lo_l3corr), nloew_base)
        _ok(f"Symlinked {os.path.basename(nloew_base)} → {os.path.basename(lo_l3corr)}")

    if os.path.exists(nloew_base):
        orig_cwd3 = os.getcwd()
        os.chdir(plotdir)
        plot_kappa3.process_and_plot(
            nlo_base_file=nloew_base,
            lo_file=lo_root,
            feature=feature,
            nbins=30,
            pt_max=300.0
        )
        os.chdir(orig_cwd3)
        _ok(f"kappa3 overlay: NLO EW-only variants ({feature})")
    else:
        _warn("Skipping NLO EW kappa3 overlay — no nloew base found")


# ─────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description='VHtrilinear Stage 4: unified analysis pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    p.add_argument('--lo-lhe', default='output/events.lhe',
                   help='Path to the LO LHE file')
    p.add_argument('--rw-lhe', default='output/events_rwgt.lhe',
                   help='Path to the reweighted LHE file')
    p.add_argument('--outdir', default='output',
                   help='Output directory for ROOT files')
    p.add_argument('--ebeam', type=float, default=6800.0,
                   help='Beam energy per proton [GeV]')
    p.add_argument('--kappa', default='0,1,2,5,10,-2,-5,-10',
                   help='Comma-separated list of κ_λ values')
    p.add_argument('--feature', default='h_pt',
                   help='Kinematic feature for kappa3 overlay plot')
    p.add_argument('--skip-lhe', action='store_true',
                   help='Skip LHE→ROOT conversion (use existing ROOT files)')
    p.add_argument('--skip-weights', action='store_true',
                   help='Skip BSM weight computation')
    p.add_argument('--skip-plots', action='store_true',
                   help='Skip all plotting steps')
    args = p.parse_args()

    kappas = [float(k) for k in args.kappa.split(',')]
    os.makedirs(args.outdir, exist_ok=True)
    plotdir = os.path.join(args.outdir, 'plots')

    # Derived ROOT file paths
    lo_root = os.path.join(args.outdir, 'events_lo.root')
    rw_root = os.path.join(args.outdir, 'events_rwgt.root')

    print("=" * 60)
    print("  VHtrilinear Analysis Pipeline")
    print(f"  LO LHE:      {args.lo_lhe}")
    print(f"  RW LHE:      {args.rw_lhe}")
    print(f"  Output dir:  {args.outdir}")
    print(f"  Beam energy: {args.ebeam} GeV")
    print(f"  κ_λ values:  {kappas}")
    print("=" * 60)

    # Step 1: LHE → ROOT
    _step("Step 1: LHE → ROOT conversion")
    if args.skip_lhe:
        _skip("LHE→ROOT")
        if not os.path.isfile(lo_root) or not os.path.isfile(rw_root):
            print("  ERROR: --skip-lhe given but ROOT files don't exist.")
            print(f"    Expected: {lo_root}")
            print(f"    Expected: {rw_root}")
            sys.exit(1)
    else:
        step_lhe_to_root(args.lo_lhe, args.rw_lhe, lo_root, rw_root, args.ebeam)

    # Step 2: BSM weights
    _step("Step 2: BSM weight computation")
    if args.skip_weights:
        _skip("BSM weights")
    else:
        step_bsm_weights(lo_root, rw_root, kappas)

    # Step 3: Plots
    _step("Step 3: Validation & Plotting")
    if args.skip_plots:
        _skip("All plots")
    else:
        step_plots(lo_root, rw_root, plotdir, kappas, args.feature)

    # Summary
    print("\n" + "=" * 60)
    print("  Analysis complete!")
    print(f"  ROOT files → {args.outdir}/")
    roots = glob.glob(os.path.join(args.outdir, '*.root'))
    for r in sorted(roots):
        print(f"    {os.path.basename(r)}")
    print(f"  Plots      → {plotdir}/")
    plots = glob.glob(os.path.join(plotdir, '*.png'))
    for pl in sorted(plots):
        print(f"    {os.path.basename(pl)}")
    print("=" * 60)


if __name__ == '__main__':
    main()
