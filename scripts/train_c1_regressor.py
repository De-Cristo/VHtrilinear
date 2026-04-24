#!/usr/bin/env python3
"""
train_c1_regressor.py — Train an XGBoost regressor to predict per-event C1 (NLO EW correction)
from truth-level ZH kinematics.

C1 = (w_ew / w_lo) * 100  [percent]

Usage:
  python3 scripts/train_c1_regressor.py [options]

Options:
  --lo-file    Path to LO ROOT file            [default: output/events_lo.root]
  --rw-file    Path to reweighted ROOT file    [default: output/events_rwgt.root]
  --outdir     Output directory for model+plots [default: output/c1_regressor]
  --test-frac  Fraction of data for testing     [default: 0.2]
  --n-rounds   XGBoost training rounds          [default: 500]
  --lr         Learning rate                    [default: 0.05]
  --max-depth  Max tree depth                   [default: 6]
  --seed       Random seed                      [default: 42]
"""
import argparse
import os
import sys
import json
import pickle
import numpy as np

try:
    import uproot
    import xgboost as xgb
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except ImportError as e:
    print(f"Missing dependency: {e}")
    print("Install with: pip install uproot xgboost scikit-learn matplotlib numpy")
    sys.exit(1)

# ── Matplotlib HEP style ──
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'DejaVu Sans'],
    'axes.linewidth': 1.2,
    'font.size': 12,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'legend.fontsize': 10,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
})
# Only use LaTeX if available
try:
    mpl.rcParams['text.usetex'] = True
except Exception:
    mpl.rcParams['text.usetex'] = False

# ── Features ──
FEATURES = ['h_pt', 'v_pt', 'vh_m', 'cos_theta_star', 'h_y', 'vh_delta_eta']
FEATURE_LABELS = {
    'h_pt': r'$p_T(H)$ [GeV]',
    'v_pt': r'$p_T(Z)$ [GeV]',
    'vh_m': r'$m(ZH)$ [GeV]',
    'cos_theta_star': r'$\cos\theta^*$',
    'h_y': r'$y(H)$',
    'vh_delta_eta': r'$\Delta\eta(Z,H)$',
}


# ═══════════════════════════════════════════════════════
# Data Loading
# ═══════════════════════════════════════════════════════

def load_data(lo_file, rw_file):
    """Load LO and reweighted ROOT files, compute C1 target, return (X, y, feature_names)."""
    print(f"  Loading {lo_file} ...")
    t_lo = uproot.open(lo_file)['events']
    ids_lo = t_lo['event_id'].array(library='np').astype(int)
    w_lo = t_lo['weight'].array(library='np')

    print(f"  Loading {rw_file} ...")
    t_rw = uproot.open(rw_file)['events']
    ids_rw = t_rw['event_id'].array(library='np').astype(int)
    w_rw = t_rw['weight'].array(library='np')

    # Match events by ID
    map_lo = dict(zip(ids_lo, w_lo))
    map_rw = dict(zip(ids_rw, w_rw))
    common = sorted(set(map_lo.keys()) & set(map_rw.keys()))
    print(f"  Common events: {len(common)}")

    # Compute C1 [%] = (w_ew / w_lo) * 100
    w_lo_arr = np.array([map_lo[i] for i in common])
    w_rw_arr = np.array([map_rw[i] for i in common])
    nonzero = w_lo_arr != 0
    C1 = np.zeros(len(common))
    C1[nonzero] = (w_rw_arr[nonzero] / w_lo_arr[nonzero]) * 100.0

    # Load features from LO file (kinematics identical in both)
    feat_arrays = {}
    for f in FEATURES:
        arr = t_lo[f].array(library='np')
        feat_arrays[f] = dict(zip(ids_lo, arr))

    X = np.column_stack([
        np.array([feat_arrays[f][i] for i in common]) for f in FEATURES
    ])

    # Remove events with w_lo == 0
    mask = nonzero
    X = X[mask]
    C1 = C1[mask]

    print(f"  Final dataset: {X.shape[0]} events × {X.shape[1]} features")
    print(f"  C1 stats: mean={C1.mean():.4f}%, std={C1.std():.4f}%, "
          f"min={C1.min():.4f}%, max={C1.max():.4f}%")

    return X, C1, FEATURES


# ═══════════════════════════════════════════════════════
# Training
# ═══════════════════════════════════════════════════════

def train_model(X_train, y_train, X_val, y_val, params):
    """Train XGBoost regressor with early stopping."""
    dtrain = xgb.DMatrix(X_train, label=y_train, feature_names=FEATURES)
    dval = xgb.DMatrix(X_val, label=y_val, feature_names=FEATURES)

    xgb_params = {
        'objective': 'reg:squarederror',
        'eval_metric': 'rmse',
        'max_depth': params['max_depth'],
        'learning_rate': params['lr'],
        'subsample': 0.8,
        'colsample_bytree': 0.8,
        'min_child_weight': 10,
        'seed': params['seed'],
        'verbosity': 1,
    }

    evals_result = {}
    model = xgb.train(
        xgb_params,
        dtrain,
        num_boost_round=params['n_rounds'],
        evals=[(dtrain, 'train'), (dval, 'val')],
        early_stopping_rounds=50,
        evals_result=evals_result,
        verbose_eval=50,
    )

    return model, evals_result


# ═══════════════════════════════════════════════════════
# Evaluation Plots
# ═══════════════════════════════════════════════════════

def plot_training_curve(evals_result, outdir):
    """Plot training/validation RMSE vs boosting round."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    train_rmse = evals_result['train']['rmse']
    val_rmse = evals_result['val']['rmse']
    rounds = np.arange(1, len(train_rmse) + 1)

    ax.plot(rounds, train_rmse, color='C0', linewidth=1.5, label='Train RMSE')
    ax.plot(rounds, val_rmse, color='C1', linewidth=1.5, label='Validation RMSE')
    ax.set_xlabel('Boosting Round')
    ax.set_ylabel('RMSE [%]')
    ax.legend(frameon=False)
    ax.set_xlim(1, len(train_rmse))
    ax.minorticks_on()
    ax.tick_params(direction='in', top=True, right=True)

    plt.tight_layout()
    out = os.path.join(outdir, 'training_curve.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  [✓] {out}")


def plot_pred_vs_true(y_true, y_pred, outdir, label='test'):
    """2D histogram: predicted C1 vs true C1."""
    fig, ax = plt.subplots(figsize=(6.5, 6))

    h = ax.hist2d(y_true, y_pred, bins=100, cmap='viridis',
                  cmin=1, norm=matplotlib.colors.LogNorm())
    plt.colorbar(h[3], ax=ax, label='Events')

    # diagonal
    lims = [min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())]
    ax.plot(lims, lims, 'r--', linewidth=1.2, label='y = x')

    ax.set_xlabel(r'True $C_1$ [\%]')
    ax.set_ylabel(r'Predicted $C_1$ [\%]')
    ax.legend(frameon=False, loc='upper left')
    ax.set_aspect('equal')
    ax.minorticks_on()
    ax.tick_params(direction='in', top=True, right=True)

    plt.tight_layout()
    out = os.path.join(outdir, f'pred_vs_true_{label}.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  [✓] {out}")


def plot_residuals(y_true, y_pred, outdir, label='test'):
    """Residual (pred - true) distribution + residual vs pT."""
    residual = y_pred - y_true

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    # Left: residual histogram
    ax1.hist(residual, bins=80, histtype='stepfilled', color='C0', alpha=0.3,
             edgecolor='C0', linewidth=1.2)
    ax1.axvline(0, color='red', linestyle='--', linewidth=1)
    mean_r = residual.mean()
    std_r = residual.std()
    ax1.axvline(mean_r, color='green', linestyle='-', linewidth=1.2,
                label=f'Mean = {mean_r:.4f}\\%')
    ax1.set_xlabel(r'Residual ($C_1^{\mathrm{pred}} - C_1^{\mathrm{true}}$) [\%]')
    ax1.set_ylabel('Events')
    ax1.legend(frameon=False)
    ax1.text(0.95, 0.92, f'$\\sigma$ = {std_r:.4f}\\%',
             transform=ax1.transAxes, ha='right', fontsize=11)
    ax1.minorticks_on()
    ax1.tick_params(direction='in', top=True, right=True)

    # Right: residual vs true C1
    ax2.hist2d(y_true, residual, bins=[80, 80], cmap='viridis',
               cmin=1, norm=matplotlib.colors.LogNorm())
    ax2.axhline(0, color='red', linestyle='--', linewidth=1)
    ax2.set_xlabel(r'True $C_1$ [\%]')
    ax2.set_ylabel(r'Residual [\%]')
    ax2.minorticks_on()
    ax2.tick_params(direction='in', top=True, right=True)

    plt.tight_layout()
    out = os.path.join(outdir, f'residuals_{label}.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  [✓] {out}")


def plot_feature_importance(model, outdir):
    """Bar chart of feature importances (gain-based)."""
    importance = model.get_score(importance_type='gain')

    # Sort by importance
    feats = sorted(importance.keys(), key=lambda f: importance[f], reverse=True)
    vals = [importance[f] for f in feats]
    labels = [FEATURE_LABELS.get(f, f) for f in feats]

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.barh(range(len(feats)), vals, color='steelblue', edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(feats)))
    ax.set_yticklabels(labels)
    ax.set_xlabel('Feature Importance (Gain)')
    ax.invert_yaxis()
    ax.minorticks_on()
    ax.tick_params(direction='in', top=True, right=True)

    plt.tight_layout()
    out = os.path.join(outdir, 'feature_importance.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  [✓] {out}")


def plot_c1_profile(X_test, y_true, y_pred, outdir):
    """Plot mean predicted vs true C1 in bins of each feature (profile plot)."""
    for fi, feat in enumerate(FEATURES):
        x_vals = X_test[:, fi]
        label = FEATURE_LABELS.get(feat, feat)

        if feat == 'cos_theta_star':
            bins = np.linspace(-1, 1, 21)
        elif feat == 'h_y':
            bins = np.linspace(-4, 4, 21)
        elif feat == 'vh_delta_eta':
            bins = np.linspace(0, 6, 21)
        else:
            q95 = np.percentile(x_vals, 95)
            bins = np.linspace(0, q95, 21)

        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        true_means = np.full(len(bin_centers), np.nan)
        pred_means = np.full(len(bin_centers), np.nan)
        true_stds = np.full(len(bin_centers), np.nan)
        pred_stds = np.full(len(bin_centers), np.nan)

        for ib in range(len(bin_centers)):
            mask = (x_vals >= bins[ib]) & (x_vals < bins[ib + 1])
            if np.sum(mask) > 10:
                true_means[ib] = np.mean(y_true[mask])
                pred_means[ib] = np.mean(y_pred[mask])
                true_stds[ib] = np.std(y_true[mask]) / np.sqrt(np.sum(mask))
                pred_stds[ib] = np.std(y_pred[mask]) / np.sqrt(np.sum(mask))

        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(7, 5.5))

        # Top: C1 profile
        ax_top.errorbar(bin_centers, true_means, yerr=true_stds, fmt='o',
                        color='C0', markersize=4, linewidth=1.2, capsize=2, label='True')
        ax_top.errorbar(bin_centers + (bins[1] - bins[0]) * 0.05, pred_means, yerr=pred_stds,
                        fmt='s', color='C1', markersize=4, linewidth=1.2, capsize=2, label='Predicted')
        ax_top.set_ylabel(r'$\langle C_1 \rangle$ [\%]')
        ax_top.legend(frameon=False)
        ax_top.minorticks_on()
        ax_top.tick_params(direction='in', top=True, right=True)

        # Bottom: ratio pred/true
        valid = ~np.isnan(true_means) & (true_means != 0)
        ratio = np.full_like(true_means, np.nan)
        ratio[valid] = pred_means[valid] / true_means[valid]
        ax_bot.plot(bin_centers, ratio, 'ko-', markersize=3, linewidth=1)
        ax_bot.axhline(1.0, color='red', linestyle='--', linewidth=1)
        ax_bot.set_xlabel(label)
        ax_bot.set_ylabel('Pred/True')
        ax_bot.set_ylim(0.9, 1.1)
        ax_bot.minorticks_on()
        ax_bot.tick_params(direction='in', top=True, right=True)

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.06)
        out = os.path.join(outdir, f'c1_profile_{feat}.png')
        plt.savefig(out, dpi=200, bbox_inches='tight')
        plt.close()
        print(f"  [✓] {out}")


def plot_c1_distribution(y_true, y_pred, outdir, label='test'):
    """Overlay true vs predicted C1 distributions."""
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(7, 5.5))

    bins = np.linspace(y_true.min() - 0.1, y_true.max() + 0.1, 60)

    h_true, _, _ = ax_top.hist(y_true, bins=bins, histtype='step', color='C0',
                                linewidth=1.5, label='True', density=True)
    h_pred, _, _ = ax_top.hist(y_pred, bins=bins, histtype='step', color='C1',
                                linewidth=1.5, linestyle='--', label='Predicted', density=True)
    ax_top.set_ylabel('Density')
    ax_top.legend(frameon=False)
    ax_top.minorticks_on()
    ax_top.tick_params(direction='in', top=True, right=True)

    # Ratio
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.where(h_true > 0, h_pred / h_true, np.nan)
    ax_bot.plot(bin_centers, ratio, 'ko-', markersize=2, linewidth=0.8)
    ax_bot.axhline(1.0, color='red', linestyle='--', linewidth=1)
    ax_bot.set_xlabel(r'$C_1$ [\%]')
    ax_bot.set_ylabel('Pred/True')
    ax_bot.set_ylim(0.5, 1.5)
    ax_bot.minorticks_on()
    ax_bot.tick_params(direction='in', top=True, right=True)

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.06)
    out = os.path.join(outdir, f'c1_distribution_{label}.png')
    plt.savefig(out, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  [✓] {out}")


# ═══════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(
        description='Train XGBoost regressor to predict per-event C1 from ZH kinematics',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument('--lo-file', default='output/events_lo.root')
    p.add_argument('--rw-file', default='output/events_rwgt.root')
    p.add_argument('--outdir', default='output/c1_regressor')
    p.add_argument('--test-frac', type=float, default=0.2)
    p.add_argument('--n-rounds', type=int, default=500)
    p.add_argument('--lr', type=float, default=0.05)
    p.add_argument('--max-depth', type=int, default=6)
    p.add_argument('--seed', type=int, default=42)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── Step 1: Load data ──
    print("\n" + "=" * 60)
    print("  C1 Regressor Training Pipeline")
    print("=" * 60)
    print("\n── Step 1: Loading data ──")
    X, y, features = load_data(args.lo_file, args.rw_file)

    # ── Step 2: Split ──
    print("\n── Step 2: Train/test split ──")
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_frac, random_state=args.seed
    )
    print(f"  Train: {X_train.shape[0]}  Test: {X_test.shape[0]}")

    # ── Step 3: Train ──
    print("\n── Step 3: Training XGBoost ──")
    params = {
        'n_rounds': args.n_rounds,
        'lr': args.lr,
        'max_depth': args.max_depth,
        'seed': args.seed,
    }
    model, evals_result = train_model(X_train, y_train, X_test, y_test, params)
    print(f"  Best iteration: {model.best_iteration}")

    # ── Step 4: Evaluate ──
    print("\n── Step 4: Evaluation ──")
    dtest = xgb.DMatrix(X_test, feature_names=FEATURES)
    y_pred = model.predict(dtest)

    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    mae = mean_absolute_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    mean_c1 = y_test.mean()

    metrics = {
        'RMSE [%]': float(rmse),
        'MAE [%]': float(mae),
        'R2': float(r2),
        'Mean C1 [%]': float(mean_c1),
        'Relative RMSE': float(rmse / abs(mean_c1)) if mean_c1 != 0 else float('inf'),
        'n_train': int(X_train.shape[0]),
        'n_test': int(X_test.shape[0]),
        'best_iteration': int(model.best_iteration),
        'features': features,
        'hyperparams': params,
    }

    print(f"  RMSE:          {rmse:.4f}%")
    print(f"  MAE:           {mae:.4f}%")
    print(f"  R²:            {r2:.4f}")
    print(f"  Mean C1:       {mean_c1:.4f}%")
    print(f"  Relative RMSE: {rmse/abs(mean_c1)*100:.1f}% of mean C1")

    # Save metrics
    metrics_path = os.path.join(args.outdir, 'metrics.json')
    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"  [✓] {metrics_path}")

    # Save model
    model_path = os.path.join(args.outdir, 'c1_regressor.json')
    model.save_model(model_path)
    print(f"  [✓] {model_path}")

    # ── Step 5: Plots ──
    print("\n── Step 5: Generating evaluation plots ──")
    plot_training_curve(evals_result, args.outdir)
    plot_pred_vs_true(y_test, y_pred, args.outdir)
    plot_residuals(y_test, y_pred, args.outdir)
    plot_feature_importance(model, args.outdir)
    plot_c1_distribution(y_test, y_pred, args.outdir)
    plot_c1_profile(X_test, y_test, y_pred, args.outdir)

    # ── Summary ──
    print("\n" + "=" * 60)
    print("  Training complete!")
    print(f"  Model:   {model_path}")
    print(f"  Metrics: {metrics_path}")
    print(f"  Plots:   {args.outdir}/")
    for f in sorted(os.listdir(args.outdir)):
        if f.endswith('.png'):
            print(f"    {f}")
    print("=" * 60)


if __name__ == '__main__':
    main()
