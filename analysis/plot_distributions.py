#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_distributions.py — 微剂量学单事件与多事件分布图

参考 Sato 2012 式 (20)-(28) (DSMK 显式卷积) 与 Inaniwa 2018 风格:
    f_d(z_d, z_n)  = Σ_k [λ(z_n)^k · e^{-λ(z_n)} / k!] · f_{d,k}(z_d)
    f_{d,k}(z_d)   = f_{d,1} * ... * f_{d,1}     (k 次卷积)
    λ(z_n)         = z_n / z̄_{d,F}

    f_n(z_n, D)    = Σ_k [λ(D)^k · e^{-λ(D)} / k!] · f_{n,k}(z_n)
    λ(D)           = D / z̄_{n,F}

输入:
    --roots        一个或多个 .root 文件路径 (从 single_events 树读取 z_d, z_n, weight)

产物 (按 --roots 顺序逐个生成, 文件名以输入 .root 基名为子目录):
    data/distribution/{filename}/{filename}_f_d1.png       (1) z_d × f_{d,1}(z_d) + f_{d,1}(z_d)
    data/distribution/{filename}/{filename}_f_n1.png       (2) z_n × f_{n,1}(z_n) + f_{n,1}(z_n)
    data/distribution/{filename}/{filename}_f_d_at_zn.png  (4) 3 个 z_n 下 f_d(z_d | z_n)
    data/distribution/{filename}/{filename}_z2f_d_at_zn.png (5) 3 个 z_n 下 z_d² · f_d(z_d | z_n)
    data/distribution/{filename}/{filename}_f_n_at_D.png   (6) 3 个 D 下 f_n(z_n | D)
    data/distribution/{filename}/{filename}_z2f_n_at_D.png (7) 3 个 D 下 z_n² · f_n(z_n | D)

用法:
    conda run -n microtrack python analysis/plot_distributions.py \\
        --roots data/ac225_phy_decay_Membrane.root \\
                 data/lu177_phy_decay_Membrane.root
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless 安全
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm


# ============================================================
# 中文字体探测 (消除 CJK glyph 缺失警告)
# ============================================================
def _setup_cn_font():
    candidates = [
        "Noto Sans CJK SC", "Noto Sans CJK JP", "Noto Sans CJK TC",
        "WenQuanYi Zen Hei", "WenQuanYi Micro Hei",
        "SimHei", "Microsoft YaHei", "PingFang SC", "Hiragino Sans GB",
    ]
    for name in candidates:
        try:
            path = _fm.findfont(name, fallback_to_default=False)
            if path and "DejaVu" not in path:
                matplotlib.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
                matplotlib.rcParams["axes.unicode_minus"] = False
                return name
        except Exception:
            continue
    matplotlib.rcParams["axes.unicode_minus"] = False
    return None
_setup_cn_font()


# ============================================================
# ROOT 读取
# ============================================================
def load_single_events(root_path, tree_name="single_events"):
    """
    从 single_events 树读取 (z_d_Gy, z_n_Gy, weight), 仅保留命中事件 (hitFlag==1)。

    :return: (z_d, z_n, w) 三个 numpy 数组, shape 一致
    """
    import ROOT
    f = ROOT.TFile.Open(str(root_path))
    if not f or f.IsZombie():
        raise RuntimeError(f"无法打开 {root_path}")
    t = f.Get(tree_name)
    if not t:
        raise RuntimeError(f"{root_path} 中没有树 {tree_name}")

    n_entries = int(t.GetEntries())
    if n_entries == 0:
        raise RuntimeError(f"{root_path} 中 {tree_name} 为空")

    # 用 RDataFrame 过滤 hitFlag==1 (single_events 几乎都满足, 但保险起见)
    rdf = ROOT.RDataFrame(str(tree_name), str(root_path))
    hits = rdf.Filter("hitFlag==1")
    cols = hits.AsNumpy(["z_d_Gy", "z_n_Gy", "weight"])
    z_d = np.asarray(cols["z_d_Gy"], dtype=float)
    z_n = np.asarray(cols["z_n_Gy"], dtype=float)
    w   = np.asarray(cols["weight"], dtype=float)
    return z_d, z_n, w


# ============================================================
# 微剂量学量计算
# ============================================================
def compute_basic_microdosimetry(z_d, z_n, w):
    """
    计算 4 个基础量 (与 MKmodel.base.compute_basic_microdosimetry 一致):
        z̄_{d,F} = Σ(w·z_d)  / Σw         (Sato 式 10)
        z̄_{d,D} = Σ(w·z_d²) / Σ(w·z_d)  (Sato 式 11)
        z̄_{n,F} = Σz_n / N               (无加权频率均)
        z̄_{n,D} = Σz_n² / Σz_n            (Sato 式 7)
    """
    Sw   = w.sum()
    Swz  = (w * z_d).sum()
    Swzz = (w * z_d * z_d).sum()
    if Sw <= 0 or Swz <= 0:
        raise RuntimeError("Σw 或 Σ(w·z_d) ≤ 0")
    z_d_F = Swz / Sw
    z_d_D = Swzz / Swz

    Szn = z_n.sum()
    if Szn <= 0:
        raise RuntimeError("Σz_n ≤ 0")
    z_n_F = Szn / len(z_n)
    z_n_D = (z_n * z_n).sum() / Szn

    return z_d_F, z_d_D, z_n_F, z_n_D


# ============================================================
# 单事件分布 → 直方图 (返回 bin centers + 概率密度)
# ============================================================
def single_event_pdf(samples, weights, z_min=None, z_max=None, n_bins=80):
    """
    加权单事件分布 f_{1}(z) 直方图 (概率密度)。
    自动用 log-spaced bins, 适合 log-log 绘图。

    :return: (centers, density, edges, mask)
        mask 标识有效 (center > 0) 的 bin
    """
    if z_min is None:
        # 排除 0 与极端 outlier
        positives = samples[samples > 0]
        z_min = max(positives.min() if positives.size else 1e-6, 1e-4)
    if z_max is None:
        z_max = max(samples.max(), 1.0)

    # log-spaced bins
    edges = np.logspace(np.log10(z_min), np.log10(z_max), n_bins + 1)
    counts, _ = np.histogram(samples, bins=edges, weights=weights, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])

    mask = counts > 0
    return centers[mask], counts[mask], edges, mask


# ============================================================
# 多事件分布: Poisson 卷积 + 单事件加权 bootstrap
# ============================================================
def multi_event_pdf(samples, weights, target, mean_for_pdf, second_moment=None,
                    n_mc=40000, k_max=None, rng=None,
                    use_gaussian_threshold=50.0):
    """
    给定 z = target (Gy), 用 Poisson(λ=target/mean) 卷积单事件分布得到多事件分布样本。
    与 analysis/MKmodel/dsmk.py::survival_n 的 MC 卷积逻辑等价, 但不求均值, 而是返回样本。

    对于 λ 大的情况 (例如 z̄_{n,F} 极小但 D 较大), 直接 MC 求和会爆内存,
    此时改用中心极限定理: sum ≈ N(λ·μ, λ·μ·z̄_{·,D}).

    :param samples:       单事件 z (Gy) 样本数组
    :param weights:       samples 对应的加权
    :param target:        目标 z (Gy): 对 f_d 是 z_n, 对 f_n 是 D
    :param mean_for_pdf:  z̄_{d,F} 或 z̄_{n,F}, 用于 λ = target / mean
    :param second_moment: z̄_{d,D} 或 z̄_{n,D} (剂量均), 用于 CLT 近似时的方差
    :param n_mc:          MC 试验数 (CLT 路径也用同样数)
    :param k_max:         Poisson 卷积最大 k
    :param rng:           np.random.Generator
    :param use_gaussian_threshold: λ > 此值改用高斯近似
    :return:              多事件 z 样本数组
    """
    if rng is None:
        rng = np.random.default_rng()

    lam = target / mean_for_pdf
    if lam < 0:
        raise ValueError(f"λ = {lam} < 0 (target={target}, mean={mean_for_pdf})")

    # —— 大 λ: 用中心极限定理近似 ——
    if lam > use_gaussian_threshold and second_moment is not None and second_moment > 0:
        mu = mean_for_pdf
        # Var[sum] = λ · μ · z̄_{·,D}  (Poisson(k) × 和 k 个 μ 样本, 推导见 docstring)
        var = lam * mu * second_moment
        std = np.sqrt(var)
        # 中心 = λ·μ = target  (按定义严格成立)
        return rng.normal(loc=target, scale=max(std, 1e-12), size=n_mc)

    # —— 中小 λ: 显式 MC 卷积 ——
    if k_max is None:
        if lam < 1e-6:
            k_max = 1
        else:
            k_max = int(np.ceil(lam + 8.0 * max(lam, 1.0)))
    ks = np.arange(0, k_max + 1)
    log_pmf = ks * np.log(lam + 1e-300) - lam - np.array(
        [math.lgamma(k + 1) for k in ks])
    p_k = np.exp(log_pmf)
    p_sum = p_k.sum()
    if p_sum <= 0:
        return np.array([])
    p_k = p_k / p_sum

    sampled_ks = rng.choice(ks, size=n_mc, p=p_k)
    w_norm = weights / weights.sum()

    out = np.zeros(n_mc)
    if k_max == 0:
        return out
    for k_val in range(1, k_max + 1):
        mask = (sampled_ks == k_val)
        n_at_k = int(mask.sum())
        if n_at_k == 0:
            continue
        idx = rng.choice(len(samples), size=(n_at_k, k_val), p=w_norm, replace=True)
        out[mask] = samples[idx].sum(axis=1)
    return out


# ============================================================
# 直方图工具 (用于多事件分布)
# ============================================================
def histogram_pdf(samples, edges):
    counts, _ = np.histogram(samples, bins=edges, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mask = counts > 0
    return centers[mask], counts[mask]


def auto_log_edges(samples, q_low=0.001, q_high=0.9999, n_bins=80):
    """按样本分位数取 log-spaced bin 边界."""
    pos = samples[samples > 0]
    if pos.size == 0:
        return np.logspace(-3, 3, n_bins + 1)
    z_min = np.quantile(pos, q_low)
    z_max = np.quantile(pos, q_high)
    if z_max <= z_min:
        z_max = z_min * 100.0
    return np.logspace(np.log10(z_min), np.log10(z_max), n_bins + 1)


# ============================================================
# 7 张图
# ============================================================
def plot_single_event_f_d(out_dir, fname, z_d, w, z_d_F, z_d_D):
    """
    图 (1): f_{d,1}(z_d) + z_d · f_{d,1}(z_d), 双对数
    """
    centers, density, _, _ = single_event_pdf(z_d, w, n_bins=80)

    fig, ax = plt.subplots(figsize=(8, 5.6), dpi=130)
    ax.loglog(centers, density, 'o-', color='tab:blue', lw=1.5, ms=4,
              label=r'$f_{d,1}(z_d)$  (单事件域比能 PDF)')
    ax.loglog(centers, centers * density, 's-', color='tab:red', lw=1.5, ms=4,
              label=r'$z_d \cdot f_{d,1}(z_d)$  (剂量加权)')

    ax.set_xlabel(r'$z_d$  (Gy)')
    ax.set_ylabel('概率密度 / 剂量密度 (1/Gy)')
    ax.set_title(f'{fname} — 单事件域比能分布  '
                 f'($\\bar z_{{d,F}}$={z_d_F:.3f}, $\\bar z_{{d,D}}$={z_d_D:.3f} Gy)')
    ax.legend(loc='best', framealpha=0.92)
    ax.grid(True, which='both', ls='--', alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_dir / f'{fname}_f_d1.png', dpi=150)
    plt.close(fig)


def plot_single_event_f_n(out_dir, fname, z_n, w, z_n_F, z_n_D):
    """
    图 (2): f_{n,1}(z_n) + z_n · f_{n,1}(z_n), 双对数
    """
    centers, density, _, _ = single_event_pdf(z_n, w, n_bins=80)

    fig, ax = plt.subplots(figsize=(8, 5.6), dpi=130)
    ax.loglog(centers, density, 'o-', color='tab:green', lw=1.5, ms=4,
              label=r'$f_{n,1}(z_n)$  (单事件核比能 PDF)')
    ax.loglog(centers, centers * density, 's-', color='tab:purple', lw=1.5, ms=4,
              label=r'$z_n \cdot f_{n,1}(z_n)$  (剂量加权)')

    ax.set_xlabel(r'$z_n$  (Gy)')
    ax.set_ylabel('概率密度 / 剂量密度 (1/Gy)')
    ax.set_title(f'{fname} — 单事件核比能分布  '
                 f'($\\bar z_{{n,F}}$={z_n_F:.4g}, $\\bar z_{{n,D}}$={z_n_D:.4g} Gy)')
    ax.legend(loc='best', framealpha=0.92)
    ax.grid(True, which='both', ls='--', alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_dir / f'{fname}_f_n1.png', dpi=150)
    plt.close(fig)


def plot_f_d_family(out_dir, fname, z_d, w, z_d_F, z_d_D, z_n_targets, z2=False,
                    rng=None):
    """
    图 (4)/(5): 多事件 f_d(z_d | z_n) 或 z_d²·f_d(z_d | z_n), 3 条 z_n 曲线在一张图
    """
    # bin 范围: 围绕 3 个 z_n 目标值的范围, 而非单事件 z_d 范围
    zn_min, zn_max = min(z_n_targets), max(z_n_targets)
    # 多事件 z_d 的典型宽度 ~ z_n (中心), 选 [zn_min/2, zn_max*3] 涵盖所有曲线
    edges = np.logspace(np.log10(zn_min / 2), np.log10(zn_max * 3), 80)
    colors = ['tab:blue', 'tab:orange', 'tab:red']

    fig, ax = plt.subplots(figsize=(8.5, 6), dpi=130)
    for zn, col in zip(z_n_targets, colors):
        z_d_samples = multi_event_pdf(z_d, w, zn, z_d_F,
                                      second_moment=z_d_D,
                                      n_mc=40000, rng=rng)
        c, pdf = histogram_pdf(z_d_samples, edges)
        if z2:
            y = c * c * pdf
            ylabel = r'$z_d^2 \cdot f_d(z_d | z_n)$  (Gy)'
            tag = r'$z_d^2 \cdot f_d$'
        else:
            y = pdf
            ylabel = r'$f_d(z_d | z_n)$  (1/Gy)'
            tag = r'$f_d$'
        ax.loglog(c, y, 'o-', color=col, lw=1.4, ms=4,
                  label=f'{tag} @ $z_n$={zn:g} Gy')

    ax.set_xlabel(r'$z_d$  (Gy)')
    ax.set_ylabel(ylabel)
    ax.set_title(f'{fname} — 多事件域比能分布 (Poisson 卷积, '
                 f'$\\bar z_{{d,F}}$={z_d_F:.3f} Gy, '
                 f'$\\bar z_{{d,D}}$={z_d_D:.3f} Gy)')
    ax.legend(loc='best', framealpha=0.92)
    ax.grid(True, which='both', ls='--', alpha=0.4)
    fig.tight_layout()
    suffix = 'z2f_d_at_zn' if z2 else 'f_d_at_zn'
    fig.savefig(out_dir / f'{fname}_{suffix}.png', dpi=150)
    plt.close(fig)


def plot_f_n_family(out_dir, fname, z_n, w, z_n_F, z_n_D, D_targets, z2=False,
                    rng=None):
    """
    图 (6)/(7): 多事件 f_n(z_n | D) 或 z_n²·f_n(z_n | D), 3 条 D 曲线在一张图
    """
    D_min, D_max = min(D_targets), max(D_targets)
    edges = np.logspace(np.log10(D_min / 2), np.log10(D_max * 3), 80)
    colors = ['tab:green', 'tab:purple', 'tab:brown']

    fig, ax = plt.subplots(figsize=(8.5, 6), dpi=130)
    for D, col in zip(D_targets, colors):
        z_n_samples = multi_event_pdf(z_n, w, D, z_n_F,
                                      second_moment=z_n_D,
                                      n_mc=40000, rng=rng)
        c, pdf = histogram_pdf(z_n_samples, edges)
        if z2:
            y = c * c * pdf
            ylabel = r'$z_n^2 \cdot f_n(z_n | D)$  (Gy)'
            tag = r'$z_n^2 \cdot f_n$'
        else:
            y = pdf
            ylabel = r'$f_n(z_n | D)$  (1/Gy)'
            tag = r'$f_n$'
        ax.loglog(c, y, 'o-', color=col, lw=1.4, ms=4,
                  label=f'{tag} @ D={D:g} Gy')

    ax.set_xlabel(r'$z_n$  (Gy)')
    ax.set_ylabel(ylabel)
    ax.set_title(f'{fname} — 多事件核比能分布 (Poisson 卷积, '
                 f'$\\bar z_{{n,F}}$={z_n_F:.4g} Gy, '
                 f'$\\bar z_{{n,D}}$={z_n_D:.4g} Gy)')
    ax.legend(loc='best', framealpha=0.92)
    ax.grid(True, which='both', ls='--', alpha=0.4)
    fig.tight_layout()
    suffix = 'z2f_n_at_D' if z2 else 'f_n_at_D'
    fig.savefig(out_dir / f'{fname}_{suffix}.png', dpi=150)
    plt.close(fig)


# ============================================================
# CLI
# ============================================================
def parse_args():
    ap = argparse.ArgumentParser(
        description="微剂量学分布图 (单事件 + Poisson 多事件卷积)")
    ap.add_argument("--roots", nargs='+', required=True,
                    help="一个或多个 .root 文件路径")
    ap.add_argument("--tree", default="single_events",
                    help="树名 (默认 single_events)")
    ap.add_argument("--out-root", default="data/distribution",
                    help="产物根目录 (默认 data/distribution, 内置 {filename} 子目录)")
    ap.add_argument("--zn-targets", type=float, nargs='+', default=[1.0, 11.0, 100.0],
                    help="f_d 系列图的 z_n 取值 (默认 1 11 100 Gy)")
    ap.add_argument("--D-targets", type=float, nargs='+', default=[0.3, 2.0, 11.4],
                    help="f_n 系列图的 D 取值 (默认 0.3 2 11.4 Gy)")
    ap.add_argument("--n-mc", type=int, default=40000,
                    help="多事件 Poisson 卷积 MC 试验数 (默认 40000)")
    ap.add_argument("--seed", type=int, default=20260714,
                    help="随机种子 (默认 20260714)")
    return ap.parse_args()


# ============================================================
# 主流程
# ============================================================
def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    out_root = Path(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    for root_path_str in args.roots:
        root_path = Path(root_path_str)
        if not root_path.exists():
            print(f"[SKIP] {root_path} 不存在")
            continue

        fname = root_path.stem  # 例如 ac225_phy_decay_Membrane
        out_dir = out_root / fname
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n=== [{fname}] 加载 {root_path} ===")
        try:
            z_d, z_n, w = load_single_events(root_path, tree_name=args.tree)
        except Exception as e:
            print(f"[ERROR] {e}")
            continue

        print(f"  命中事件: {len(z_d)}")
        z_d_F, z_d_D, z_n_F, z_n_D = compute_basic_microdosimetry(z_d, z_n, w)
        print(f"  z_d_F={z_d_F:.4f}, z_d_D={z_d_D:.4f}, "
              f"z_n_F={z_n_F:.4g}, z_n_D={z_n_D:.4g}")

        # 1) z_d × f_{d,1}(z_d)
        print("  [1] f_{d,1}(z_d) + z_d·f_{d,1}(z_d) ...")
        plot_single_event_f_d(out_dir, fname, z_d, w, z_d_F, z_d_D)

        # 2) z_n × f_{n,1}(z_n)
        print("  [2] f_{n,1}(z_n) + z_n·f_{n,1}(z_n) ...")
        plot_single_event_f_n(out_dir, fname, z_n, w, z_n_F, z_n_D)

        # 4) f_d(z_d | z_n)  (3 条 z_n 曲线)
        print(f"  [4] f_d(z_d | z_n) @ z_n={args.zn_targets} ...")
        plot_f_d_family(out_dir, fname, z_d, w, z_d_F, z_d_D, args.zn_targets,
                        z2=False, rng=rng)

        # 5) z_d² · f_d(z_d | z_n)
        print(f"  [5] z_d²·f_d(z_d | z_n) @ z_n={args.zn_targets} ...")
        plot_f_d_family(out_dir, fname, z_d, w, z_d_F, z_d_D, args.zn_targets,
                        z2=True, rng=rng)

        # 6) f_n(z_n | D)  (3 条 D 曲线)
        print(f"  [6] f_n(z_n | D) @ D={args.D_targets} ...")
        plot_f_n_family(out_dir, fname, z_n, w, z_n_F, z_n_D, args.D_targets,
                        z2=False, rng=rng)

        # 7) z_n² · f_n(z_n | D)
        print(f"  [7] z_n²·f_n(z_n | D) @ D={args.D_targets} ...")
        plot_f_n_family(out_dir, fname, z_n, w, z_n_F, z_n_D, args.D_targets,
                        z2=True, rng=rng)

        print(f"  [OK] 6 张图 → {out_dir}")

    print("\n[ALL DONE]")


if __name__ == "__main__":
    main()