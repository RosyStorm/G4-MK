#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MKplot.py — Ac-225 vs Lu-177 存活曲线叠加 + RBE 曲线

两幅图 (画在同一张 PNG 上):
    上: 两个模型的存活曲线 S(D) vs D (ModifiedSMK Inaniwa2018 式24)
    下: Ac-225 的 RBE(S) = D_Lu177(S) / D_Ac225(S) vs Ac-225 剂量 D_Ac225

用法:
    cd microtrack
    python analysis/MKplot.py [--z0 66.0] [--dmax 10.0] [--outdir result/]

产物: <outdir>/mks_rbe_fit.png  (单图含两个 panel)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from MKmodel.msmk import ModifiedSMK
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")  # 无显示环境也安全 (headless 服务)
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm

# —— 中文字体探测 (消除 CJK glyph missing 警告) ——
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
# 命令行参数
# ============================================================
def params_args():
    ap = argparse.ArgumentParser(
        description="Ac-225 vs Lu-177 存活曲线 + Ac-225 RBE 曲线 (Modified SMK)")
    ap.add_argument("--tree", default="single_events",
                    help="microtrack 模拟 ntuple 树名 (默认 single_events)")
    ap.add_argument("--root-ac225", default="data/ac225_single_decay_Nucleus.root",
                    help="Ac-225 路线 2 产物 .root 路径")
    ap.add_argument("--root-lu177", default="data/lu177_phy_decay_Nucleus.root",
                    help="Lu-177 路线 2 产物 .root 路径")
    ap.add_argument("--alpha0", type=float, default=0.127,
                    help="参考 α₀ (默认 0.226; 也可从 --lq-csv 读)")
    ap.add_argument("--beta0", type=float, default=0.014,
                    help="参考 β₀ (默认 0.026)")
    ap.add_argument("--z0", type=float, default=144.7,
                    help="饱和参数 z₀ [Gy] (Inaniwa HSG)")
    ap.add_argument("--dmax", type=float, default=None,
                    help="存活/RBE 曲线剂量上限 [Gy]; 默认取两份数据最大值的 1.05 倍")
    ap.add_argument("--outdir", default="result/RBE",
                    help="输出目录 (默认 result/RBE)")
    ap.add_argument("--outname", default="mks_rbe_fit",
                    help="输出 PNG 文件名前缀 (默认 mks_rbe_fit)")
    ap.add_argument("--smin", type=float, default=0.01,
                    help="RBE 计算所用最小存活率 (默认 0.01, 避免数值不稳)")
    ap.add_argument("--smax", type=float, default=0.99,
                    help="RBE 计算所用最大存活率 (默认 0.99, 避免数值不稳)")
    ap.add_argument("--npoints-D", type=int, default=1000,
                    help="剂量网格点数 (默认 400)")
    ap.add_argument("--npoints-S", type=int, default=100,
                    help="RBE 计算所用存活率点数 (默认 100)")
    ap.add_argument("--fig-mode", default="normal", choices=["normal", "log"],
                    help="存活曲线 y 轴模式 (默认 normal)")
    args = ap.parse_args()
    return args


# ============================================================
# RBE 反演: 给定 S 求 D
# ============================================================
def find_dose_at_survival(model, S_target, dose_grid):
    """
    对一个 MKM 模型, 在 dose_grid 上算 S_grid, 然后反查指定存活率 S_target 的剂量。
    survival 是 D 的单调递减函数; 用 np.interp 时把序列翻转以满足单调递增假设。

    :param model:      有 .survival(D) 方法的 MKM 模型
    :param S_target:   目标存活率 (标量)
    :param dose_grid:  shape=(N,) 的剂量数组, 单调递增
    :return:           插值得到的 D (np.interp 在范围外返回端点值)
    """
    S_grid = model.survival(dose_grid)
    # S 随 D 单调递减 → 翻转两端让 np.interp 看到单调递增
    return float(np.interp(S_target, S_grid[::-1], dose_grid[::-1]))


def compute_rbe_curve(model_test, model_ref, dose_grid,
                      s_min=0.01, s_max=0.99, n_points=100,
                      survival_test=None, survival_ref=None):
    """
    对给定剂量网格逐点反查每个存活率下 test 与 reference 模型的剂量,
    返回 RBE 数组及对应的剂量点。

    RBE 定义: 在相同存活率 S 下,
        RBE_test(S) = D_ref(S) / D_test(S)
    这里 test = Ac-225, ref = Lu-177 — Ac-225 是 α 高 LET,
    杀伤更强, 同 S 下所需剂量更小, RBE > 1。

    :param survival_test: 可选预计算的 test 模型 S_grid, 避免重复算
    :param survival_ref:  同上, ref 模型 S_grid
    :return: (D_test_arr, D_ref_arr, RBE_arr) 同形数组, 长度 ≤ n_points
    """
    S_test_grid = (model_test.survival(dose_grid)
                   if survival_test is None else survival_test)
    S_ref_grid  = (model_ref.survival(dose_grid)
                   if survival_ref is None else survival_ref)

    # 求公共 (有定义) 的存活率区间:
    # 存活下界 = max(两个模型各自的最低 S) — 取较高的那个 (更深剂量才会有 S 更低)
    # 存活上界 = min(两个模型各自的最高 S) — 取较低的那个
    s_lo = max(s_min, float(np.max([S_test_grid.min(), S_ref_grid.min()])))
    s_hi = min(s_max, float(np.min([S_test_grid.max(), S_ref_grid.max()])))
    if s_hi <= s_lo:
        return np.array([]), np.array([]), np.array([])

    # S 从上界走到下界 (np.interp 要求 S_targets 在 S_grid[::-1] 范围内)
    S_targets = np.linspace(s_hi, s_lo, n_points)

    D_test_arr, D_ref_arr, RBE_arr = [], [], []
    for S_target in S_targets:
        # 反插 D: 把 S_grid 与 dose_grid 同步翻转 (np.interp 要求 x 单调递增)
        d_t = float(np.interp(S_target, S_test_grid[::-1], dose_grid[::-1]))
        d_r = float(np.interp(S_target, S_ref_grid[::-1],  dose_grid[::-1]))
        if d_t > 0 and np.isfinite(d_t) and np.isfinite(d_r):
            D_test_arr.append(d_t)
            D_ref_arr.append(d_r)
            RBE_arr.append(d_r / d_t)
    return np.array(D_test_arr), np.array(D_ref_arr), np.array(RBE_arr)


# ============================================================
# 主流程
# ============================================================
def main():
    args = params_args()
    root_path_ac225 = Path(args.root_ac225)
    root_path_lu177 = Path(args.root_lu177)
    tree_name = args.tree
    alpha0 = args.alpha0
    beta0 = args.beta0

    # —— 构建两个 ModifiedSMK 模型 ——
    print(f"[MKplot] 加载 {root_path_lu177} (Lu-177) ...")
    msmk_lu177 = ModifiedSMK.from_root(root_path_lu177, tree_name, alpha0, beta0, args.z0)
    print(f"[MKplot] 加载 {root_path_ac225} (Ac-225) ...")
    msmk_ac225 = ModifiedSMK.from_root(root_path_ac225, tree_name, alpha0, beta0, args.z0)

    print(f"RBE_M = {msmk_ac225.alpha_S / msmk_lu177.alpha_S:.3f}  ")

    # —— 自动剂量范围 ——
    dmax = args.dmax
    if dmax is None:
        # 取两份产物中最大 αS / βS 决定的剂量
        # (保守取 6 倍 α₀⁻¹ 让曲线延到 D≈0.05)
        dmax = max(6.0 / msmk_ac225.alpha_S, 6.0 / msmk_lu177.alpha_S)
    print(f"[MKplot] 剂量范围: 0 ~ {dmax:.2f} Gy")

    D_grid = np.linspace(0.0, dmax, args.npoints_D)

    # —— 算存活曲线与 RBE 曲线 ——
    S_lu = msmk_lu177.survival(D_grid)
    S_ac = msmk_ac225.survival(D_grid)
    D_test, D_ref, RBE = compute_rbe_curve(
        msmk_ac225, msmk_lu177, D_grid,
        survival_test=S_ac, survival_ref=S_lu,  # 预计算好的 S_grid, 避免函数内重算
        s_min=args.smin, s_max=args.smax, n_points=args.npoints_S,
    )

    print(f"[MKplot] RBE 有效点: {len(RBE)}  (S ∈ [{args.smax:.3f}, {args.smin:.3f}])")
    if len(RBE) > 0:
        print(f"[MKplot] RBE 范围: {RBE.min():.3f} ~ {RBE.max():.3f}")

    # —— 画双面板图 ——
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 9), dpi=130,
        gridspec_kw={"height_ratios": [2, 1]}, sharex=False,
    )

    # 上: 存活曲线
    ax1.plot(D_grid, S_lu, color="tab:blue",   lw=2.0, label="Lu-177 (β 链)")
    ax1.plot(D_grid, S_ac, color="tab:orange", lw=2.0, label="Ac-225 (α 链)")
    ax1.set_ylabel("存活分数 S(D)")
    ax1.set_title(
        f"ModifiedSMK  (α₀={alpha0}, β₀={beta0}, z₀={args.z0}, "
        f"tree={tree_name}, 路径={msmk_lu177.__class__.__module__})"
    )
    ax1.legend(loc="best", framealpha=0.92)
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    if args.fig_mode == "log":
        ax1.set_yscale("log")
    # 在存活图上加 RBE 标注点 (中间值)
    if len(RBE) >= 5:
        mid = len(RBE) // 2
        ax1.annotate(
            f"RBE≈{RBE[mid]:.2f}",
            xy=(D_test[mid], float(np.interp(D_test[mid], D_grid[::-1], S_ac[::-1]))),
            xytext=(D_test[mid] * 0.5, 0.30),
            fontsize=10, color="tab:red",
            arrowprops=dict(arrowstyle="->", color="gray", lw=1, alpha=0.6),
        )

    # 下: RBE vs D_Ac225
    if len(RBE) > 0:
        ax2.plot(D_test, RBE, color="tab:red", lw=2.0,
                 label="Ac-225 RBE (= D_Lu177(S) / D_Ac225(S))")
        ax2.axhline(1.0, color="gray", ls=":", lw=1, alpha=0.6)
        ax2.set_xlabel("Ac-225 剂量 D_Ac225 (Gy)")
        ax2.set_ylabel("RBE")
        ax2.set_title(f"Ac-225 RBE vs Ac-225 剂量  "
                      f"(同存活率下 Lu-177 剂量 / Ac-225 剂量)")
        ax2.legend(loc="best", framealpha=0.92)
        ax2.grid(True, which="both", ls="--", alpha=0.4)
    else:
        ax2.text(0.5, 0.5, "RBE 区间为空 (两个模型的存活范围无重叠)",
                 ha="center", va="center", transform=ax2.transAxes)

    fig.tight_layout()
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"{args.outname}.png"
    fig.savefig(png_path, dpi=300)
    print(f"[MKplot] PNG → {png_path}")

    # —— 同时写 RBE 数值表到 CSV 便于论文 / 报告使用 ——
    if len(RBE) > 0:
        csv_path = out_dir / f"{args.outname}.csv"
        import csv as csv_mod
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv_mod.writer(f)
            w.writerow(["survival_S", "D_Ac225_Gy", "D_Lu177_Gy", "RBE_Ac225_vs_Lu177"])
            for S_t, d_t, d_r, rbe_v in zip(
                np.linspace(args.smax, args.smin, len(RBE)),
                D_test, D_ref, RBE,
            ):
                w.writerow([f"{S_t:.6f}", f"{d_t:.6f}",
                            f"{d_r:.6f}", f"{rbe_v:.6f}"])
        print(f"[MKplot] CSV → {csv_path}")

    print("[MKplot] 完成。")


if __name__ == "__main__":
    main()
