#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_smk_ac225_overlay.py — 把修正 SMK 存活曲线与 Ac-225 实验数据画在一张图上

核心思路: 用 PC-3 PIP 细胞 X 射线 LQ 拟合得到的 α₀, β₀ 作为参考辐射敏感性,
配合 microtrack 模拟 ntuple 算出的微剂量学量, 按修正 SMK 模型 (Inaniwa 2018)
预测 Ac-225 的存活曲线, 再叠加 data/validation/Ac225_survival.csv 的实验点。

链路 (与 analysis/analyze_mk.C 完全一致):
    1) α₀, β₀ ← analysis/fit_lq_xray.py 的输出 CSV (默认 actual 模式)
    2) z₀      ← 饱和参数 (analyze_mk.C 默认 66.0 Gy, Inaniwa HSG; PC-3 无公开值)
    3) 微剂量学量 ← data/microtrack.root 的 events ntuple (仅命中事件 hitFlag==1):
         z̄_{d,F} = Σ(w·z_d)   / Σw            (式10, 域频率均)
         z̄_{d,D} = Σ(w·z_d²)  / Σ(w·z_d)      (式11, 域剂量均)
         z̄*_{d,D}= Σ(w·z*²)   / Σ(w·z_d)      (式12, 域饱和剂量均, z*=z₀(1-e^{-(z_d/z₀)²}))
         z̄_{n,D} = Σ z_n²     / Σ z_n          (核剂量均, 无权)
    4) SMK 系数 (Inaniwa 式15-16):
         α_S = α₀ + β₀·z̄*_{d,D}
         β_S = β₀·z̄*_{d,D} / z̄_{d,D}
    5) SMK 存活 (Inaniwa 式24):
         S(D) = exp(-α_S D - β_S D²) · [1 + (D·z̄_{n,D}/2)·((α_S+2β_S D)² - 2β_S)]

用法 (从 microtrack/ 根目录运行):
    conda run -n microtrack python analysis/plot_smk_ac225_overlay.py
    conda run -n microtrack python analysis/plot_smk_ac225_overlay.py \\
        --lq-csv result/validation/lq_fit_xray_nominal.csv --z0 89 --dmax 13

输出 (默认 outdir=result/mod-SMK):
    smk_ac225_overlay.png / .csv / .txt
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# PyROOT 仅用于读 microtrack.root 的 events ntuple (与 analyze_mk.C 同源)
import ROOT


# ============================================================
# 中文字体设置 (与 fit_lq_xray.py 一致, 防止图像中文乱码)
# ============================================================
def setup_matplotlib_font():
    """探测系统中可用的中文字体 (优先 Noto CJK)."""
    from matplotlib import font_manager
    matplotlib.rcParams["axes.unicode_minus"] = False
    candidates = [
        "Noto Sans CJK JP", "Noto Sans CJK SC", "Noto Sans CJK TC", "Noto Sans CJK KR",
        "Noto Serif CJK JP", "Noto Serif CJK SC",
        "AR PL UMing CN", "AR PL UMing TW MBE", "AR PL UKai CN",
        "WenQuanYi Zen Hei", "WenQuanYi Micro Hei",
        "SimHei", "Microsoft YaHei", "PingFang SC", "Hiragino Sans GB",
    ]
    for name in candidates:
        try:
            path = font_manager.findfont(name, fallback_to_default=False)
            if path and "DejaVu" not in path:
                matplotlib.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
                print(f"[matplotlib] 使用字体: {name}")
                return name
        except Exception:
            continue
    print("[matplotlib] 未找到中文字体, 使用默认 (图像中中文可能为方块)")
    return None


# ============================================================
# 读取 X 射线 LQ 拟合参数 α₀, β₀ (来自 fit_lq_xray.py 输出)
# ============================================================
def read_lq_params(csv_path):
    """
    从 fit_lq_xray.py 输出的拟合曲线 CSV 读取 α₀, β₀.
    该 CSV 含常量列 alpha, beta (每行重复同一拟合值), 取首行即可.
    """
    df = pd.read_csv(csv_path)
    for col in ("alpha", "beta"):
        if col not in df.columns:
            sys.exit(f"[ERROR] {csv_path} 缺少必需列 '{col}' "
                     f"(应含 dose_Gy,S_LQ,alpha,beta,mode; 请先运行 fit_lq_xray.py)")
    alpha0 = float(df["alpha"].iloc[0])
    beta0 = float(df["beta"].iloc[0])
    mode = str(df["mode"].iloc[0]) if "mode" in df.columns else "?"
    return alpha0, beta0, mode


# ============================================================
# 微剂量学量 (复刻 analyze_mk.C 第 2-3 步的累加逻辑)
# ============================================================
def compute_microdosimetry(root_path, z0):
    """
    用 RDataFrame 读 events ntuple, 仅取命中事件 hitFlag==1,
    计算 SMK 所需的单事件微剂量学量。与 analyze_mk.C 严格一一对应。

    返回 dict: z_d_F, z_d_D, zs_d_D, z_n_F, z_n_D, n_hits, n_entries
    """
    rdf = ROOT.RDataFrame("events", str(root_path))
    n_entries = rdf.Count().GetValue()
    if n_entries == 0:
        sys.exit(f"[ERROR] {root_path} 内 events ntuple 为空")

    hits = rdf.Filter("hitFlag==1")
    n_hits = hits.Count().GetValue()
    if n_hits == 0:
        sys.exit(f"[ERROR] {root_path} 无命中事件 (hitFlag==1), 无法计算微剂量学量")

    cols = hits.AsNumpy(["z_d_Gy", "z_n_Gy", "weight"])
    z_d = np.asarray(cols["z_d_Gy"], dtype=float)
    z_n = np.asarray(cols["z_n_Gy"], dtype=float)
    w = np.asarray(cols["weight"], dtype=float)

    # 域级: 带重要性权 w (校正位点随机放置偏置)
    Sw = w.sum()
    Swz = (w * z_d).sum()
    Swzz = (w * z_d * z_d).sum()
    if Sw <= 0 or Swz <= 0:
        sys.exit("[ERROR] 域级求和异常 (Σw 或 Σ(w·z_d) ≤ 0)")

    z_d_F = Swz / Sw                                  # z̄_{d,F}  式10
    z_d_D = Swzz / Swz                                # z̄_{d,D}  式11
    zs = z0 * (1.0 - np.exp(-(z_d / z0) ** 2))        # z*_{d}    式1 (无 sqrt)
    Swzszs = (w * zs * zs).sum()
    zs_d_D = Swzszs / Swz                             # z̄*_{d,D} 式12 (分母为 z̄_{d,F}, P0 修复 #5)

    # 核级: 无权 (w=1)
    Szn = z_n.sum()
    Sznzn = (z_n * z_n).sum()
    if Szn <= 0:
        sys.exit("[ERROR] 核级求和异常 (Σz_n ≤ 0)")
    z_n_F = Szn / n_hits                              # z̄_{n,F}
    z_n_D = Sznzn / Szn                               # z̄_{n,D}

    return {
        "z_d_F": z_d_F, "z_d_D": z_d_D, "zs_d_D": zs_d_D,
        "z_n_F": z_n_F, "z_n_D": z_n_D,
        "n_hits": int(n_hits), "n_entries": int(n_entries),
    }


# ============================================================
# SMK 系数 (Inaniwa 式15-16)
# ============================================================
def smk_coefficients(alpha0, beta0, z_d_D, zs_d_D):
    """α_S = α₀ + β₀·z̄*_{d,D}  (式15);  β_S = β₀·z̄*_{d,D}/z̄_{d,D}  (式16)."""
    alpha_s = alpha0 + beta0 * zs_d_D
    beta_s = beta0 * zs_d_D / z_d_D
    return alpha_s, beta_s


# ============================================================
# SMK 存活模型 (Inaniwa 式24)
# ============================================================
def smk_survival(D, alpha_s, beta_s, z_n_D):
    """
    S(D) = exp(-α_S D - β_S D²) · [1 + (D·z̄_{n,D}/2)·((α_S+2β_S D)² - 2β_S)]
    方括号为 z_n 随机性修正; 低 LET 时 z̄_{n,D}→小, 退化为 LQ。
    bracket 裁剪到 ≥0, 防止数值上出现负存活 (与 analyze_mk.C 一致)。
    """
    D = np.asarray(D, dtype=float)
    lq = np.exp(-alpha_s * D - beta_s * D ** 2)
    bracket = 1.0 + (D * z_n_D / 2.0) * ((alpha_s + 2.0 * beta_s * D) ** 2 - 2.0 * beta_s)
    bracket = np.clip(bracket, 0.0, None)             # 防御
    return lq * bracket


# ============================================================
# 主流程
# ============================================================
def main():
    ap = argparse.ArgumentParser(
        description="修正 SMK 存活曲线 (α₀/β₀ 取自 X 射线 LQ 拟合) × Ac-225 实验数据叠加图"
    )
    ap.add_argument("--root", default="data/microtrack.root",
                    help="microtrack 模拟 ntuple 路径 (默认 data/microtrack.root)")
    ap.add_argument("--lq-csv", default="result/validation/lq_fit_xray_actual.csv",
                    help="fit_lq_xray.py 输出的拟合曲线 CSV (含 alpha,beta 列); "
                         "默认 actual 模式, 可换 lq_fit_xray_nominal.csv")
    ap.add_argument("--ac225", default="data/validation/Ac225_survival.csv",
                    help="Ac-225 实验存活数据 CSV (默认 data/validation/Ac225_survival.csv)")
    ap.add_argument("--z0", type=float, default=66.0,
                    help="饱和参数 z₀ [Gy] (analyze_mk.C 默认 66.0, Inaniwa HSG; PC-3 无公开值)")
    ap.add_argument("--dmax", type=float, default=None,
                    help="存活曲线剂量上限 [Gy]; 默认取数据最大剂量的 1.05 倍")
    ap.add_argument("--outdir", default="result/validation",
                    help="输出目录 (默认 result/validation, 与 analyze_mk.C 一致)")
    args = ap.parse_args()

    root_path = Path(args.root)
    lq_path = Path(args.lq_csv)
    ac225_path = Path(args.ac225)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for p in (root_path, lq_path, ac225_path):
        if not p.exists():
            sys.exit(f"[ERROR] 找不到输入文件: {p}")

    setup_matplotlib_font()

    # ---------- 1. 读 α₀, β₀ (来自 X 射线 LQ 拟合) ----------
    alpha0, beta0, lq_mode = read_lq_params(lq_path)
    print("=" * 64)
    print("[1] 参考辐射敏感性 (来自 fit_lq_xray.py)")
    print(f"    来源: {lq_path}  (模式: {lq_mode})")
    print(f"    α₀ = {alpha0:.6f} Gy⁻¹   β₀ = {beta0:.6f} Gy⁻²")

    # ---------- 2. 算微剂量学量 ----------
    print("=" * 64)
    print("[2] 微剂量学量 (来自 microtrack.root, 仅命中事件)")
    md = compute_microdosimetry(root_path, args.z0)
    print(f"    命中事件: {md['n_hits']}/{md['n_entries']}  "
          f"(命中率 {100.0*md['n_hits']/md['n_entries']:.1f}%)  z₀ = {args.z0} Gy")
    print(f"    z̄_{{d,F}} = {md['z_d_F']:.5f} Gy   z̄_{{d,D}} = {md['z_d_D']:.5f} Gy")
    print(f"    z̄*_{{d,D}}= {md['zs_d_D']:.5f} Gy   (饱和)")
    print(f"    z̄_{{n,F}} = {md['z_n_F']:.5g} Gy   z̄_{{n,D}} = {md['z_n_D']:.5f} Gy")

    # ---------- 3. SMK 系数 (式15-16) ----------
    alpha_s, beta_s = smk_coefficients(alpha0, beta0, md["z_d_D"], md["zs_d_D"])
    print("=" * 64)
    print("[3] SMK 系数 (Inaniwa 式15-16)")
    print(f"    α_S = α₀ + β₀·z̄*_{{d,D}}            = {alpha_s:.5f} Gy⁻¹")
    print(f"    β_S = β₀·z̄*_{{d,D}}/z̄_{{d,D}}        = {beta_s:.5f} Gy⁻²")

    # ---------- 4. Ac-225 实验数据 ----------
    df_ac = pd.read_csv(ac225_path)
    for col in ("dose_Gy", "survival_fraction"):
        if col not in df_ac.columns:
            sys.exit(f"[ERROR] {ac225_path} 缺少必需列 '{col}'")
    D_obs = df_ac["dose_Gy"].to_numpy(float)
    S_obs = df_ac["survival_fraction"].to_numpy(float)

    # ---------- 5. SMK 存活曲线 (式24) ----------
    Dmax = args.dmax if args.dmax is not None else float(D_obs.max()) * 1.05
    if Dmax <= 0:
        sys.exit(f"[ERROR] --dmax 必须 > 0 (当前 {Dmax})")
    D_grid = np.linspace(0.0, Dmax, 600)
    S_grid = smk_survival(D_grid, alpha_s, beta_s, md["z_n_D"])
    S_grid = np.clip(S_grid, 1e-12, None)             # 对数轴下限保护

    # 在实验剂量点处取模型值, 算 log-RMSE (D>0); 取 log 前下限保护, 防 log(0)→-inf
    pos = D_obs > 0
    S_at = smk_survival(D_obs[pos], alpha_s, beta_s, md["z_n_D"])
    lg_obs = np.log(np.clip(S_obs[pos], 1e-30, None))
    lg_at = np.log(np.clip(S_at, 1e-30, None))
    rmse_log = float(np.sqrt(np.mean((lg_obs - lg_at) ** 2)))
    print("=" * 64)
    print("[4] 模型 vs Ac-225 实验 (log-space)")
    print(f"    D 范围: 0–{D_obs.max():.3f} Gy,  数据点 {pos.sum()} 个")
    print(f"    RMSE(ln S) = {rmse_log:.5f}")

    # ---------- 6. 画图 ----------
    fig, ax = plt.subplots(figsize=(8.0, 5.6), dpi=150)
    ax.scatter(D_obs[pos], S_obs[pos], s=66, c="#1f77b4", edgecolor="black",
               label="Ac-225 实验 (PC-3 PIP)", zorder=5)
    ax.plot(D_grid, S_grid, color="#d62728", lw=2.2,
            label=f"修正 SMK 式(24)")
    ax.axhline(0.1, color="gray", ls=":", lw=1)
    ax.axhline(0.01, color="gray", ls=":", lw=1)

    ax.set_xlabel("Ac-225 吸收剂量 D (Gy)")
    ax.set_ylabel("存活分数 S(D)")
    ax.set_title(f"修正 SMK 存活曲线 vs Ac-225 实验数据 (PC-3 PIP, $z_0$={args.z0} Gy)")
    ax.set_xlim(left=0)
    ax.set_ylim(1e-6, 1.1)
    # ax.set_yscale("log")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(loc="upper right", framealpha=0.92)

    txt = (
        f"$\\alpha_0$ = {alpha0:.4f} Gy$^{{-1}}$   $\\beta_0$ = {beta0:.4f} Gy$^{{-2}}$  ({lq_mode})\n"
        f"$z_0$ = {args.z0} Gy\n"
        f"$z^*_{{d,D}}$ = {md['zs_d_D']:.4f}   $z_{{d,D}}$ = {md['z_d_D']:.4f}\n"
        f"$z_{{n,D}}$ = {md['z_n_D']:.5f}\n"
        f"$\\alpha_S$ = {alpha_s:.4f}   $\\beta_S$ = {beta_s:.5f}\n"
        f"RMSE(ln S) = {rmse_log:.4f}"
    )
    ax.text(0.6, 0.6, txt, transform=ax.transAxes, fontsize=9.5, va="bottom",
            bbox=dict(boxstyle="round", facecolor="lightyellow",
                      edgecolor="gray", alpha=0.92))
    fig.tight_layout()

    png = out_dir / f"smk_ac225_overlay_{lq_mode}.png"
    fig.savefig(png, dpi=150)
    plt.close(fig)
    print("=" * 64)
    print(f"[输出] 图像 → {png}")

    # ---------- 7. 曲线 CSV + 摘要 TXT ----------
    curve_csv = out_dir / f"smk_ac225_{lq_mode}.csv"
    pd.DataFrame({
        "dose_Gy": D_grid,
        "S_SMK": S_grid,
        "alpha_S": np.full_like(D_grid, alpha_s),
        "beta_S": np.full_like(D_grid, beta_s),
        "z_n_D": np.full_like(D_grid, md["z_n_D"]),
    }).to_csv(curve_csv, index=False)

    txt_path = out_dir / f"smk_ac225_overlay_{lq_mode}.txt"
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("修正 SMK 存活曲线 vs Ac-225 实验 (PC-3 PIP)\n")
        f.write("=" * 56 + "\n")
        f.write("模型: S(D)=exp(-α_S D - β_S D²)·[1+(D·z̄_{n,D}/2)·((α_S+2β_S D)²-2β_S)]  (Inaniwa 式24)\n")
        f.write(f"α₀, β₀ 来源: {lq_path} ({lq_mode})\n")
        f.write(f"α₀ = {alpha0:.10g} Gy⁻¹\n")
        f.write(f"β₀ = {beta0:.10g} Gy⁻²\n")
        f.write(f"z₀ = {args.z0} Gy\n")
        f.write("-" * 56 + "\n")
        f.write(f"z̄_{{d,F}} = {md['z_d_F']:.10g} Gy\n")
        f.write(f"z̄_{{d,D}} = {md['z_d_D']:.10g} Gy\n")
        f.write(f"z̄*_{{d,D}}= {md['zs_d_D']:.10g} Gy\n")
        f.write(f"z̄_{{n,F}} = {md['z_n_F']:.10g} Gy\n")
        f.write(f"z̄_{{n,D}} = {md['z_n_D']:.10g} Gy\n")
        f.write(f"α_S = {alpha_s:.10g} Gy⁻¹\n")
        f.write(f"β_S = {beta_s:.10g} Gy⁻²\n")
        f.write(f"RMSE(ln S, Ac-225 点) = {rmse_log:.10g}\n")
        f.write("-" * 56 + "\n")
        f.write("dose_Gy\tobserved_S\tmodel_S\n")
        for d, s in zip(D_obs, S_obs):
            sm = float(smk_survival(np.array([d]), alpha_s, beta_s, md["z_n_D"])[0])
            f.write(f"{d:.10g}\t{s:.10g}\t{sm:.10g}\n")
    print(f"[输出] 曲线 → {curve_csv}")
    print(f"[输出] 摘要 → {txt_path}")
    print("=" * 64)


if __name__ == "__main__":
    main()
