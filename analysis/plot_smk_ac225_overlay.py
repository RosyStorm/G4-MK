#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_smk_ac225_overlay.py — 把修正 SMK 存活曲线与 Ac-225 实验数据画在一张图上

核心思路: 用 PC-3 PIP 细胞 X 射线 LQ 拟合得到的 α₀, β₀ 作为参考辐射敏感性,
配合 microtrack 模拟 ntuple 算出的微剂量学量, 按修正 SMK 模型 (Inaniwa 2018)
预测 Ac-225 的存活曲线, 再叠加 data/validation/Ac225_survival.csv 的实验点。

注: SMK 模型逻辑(饱和 z*、系数 α_S/β_S、存活 S(D)、微剂量学量计算)已抽到
    analysis/smk.py 的 ModifiedSMK 类中, 本脚本只负责 I/O(读 α₀/β₀、实验数据)
    与画图, 通过该类完成建模。

链路 (与 analysis/analyze_mk.C 完全一致):
    1) α₀, β₀ ← analysis/fit_lq_xray.py 的输出 CSV (默认 actual 模式)
    2) z₀      ← 饱和参数 (analyze_mk.C 默认 66.0 Gy, Inaniwa HSG; PC-3 无公开值)
    3) 微剂量学量 + SMK 系数 + S(D) ← ModifiedSMK.from_root(...)
    4) 叠加 Ac-225 实验点, 算 log-RMSE, 出图/CSV/TXT

用法 (从 microtrack/ 根目录运行):
    conda run -n microtrack python analysis/plot_smk_ac225_overlay.py
    conda run -n microtrack python analysis/plot_smk_ac225_overlay.py \\
        --lq-csv result/validation/lq_fit_xray_nominal.csv --z0 89 --dmax 13

输出 (默认 outdir=result/validation/smk-fit):
    smk_ac225_overlay.png / .csv / .txt
"""

import argparse
import sys
from pathlib import Path

# —— 让 analysis/ 进入模块搜索路径, 使 from smk import ... 在任意 cwd 下都成立 ——
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from msmk import ModifiedSMK   # ← SMK 模型类 (analysis/msmk.py)
from mk import MK
from smk import SMK
from dsmk import DSMK

# ============================================================
# 中文字体设置 (防止图像中文乱码)
# ============================================================
def setup_matplotlib_font():
    """
    探测系统中可用的中文字体 (优先 Noto CJK),
    并写入 matplotlib 全局配置, 防止图像中中文乱码.

    Returns:
        str 或 None: 成功匹配到的字体名; 若系统未安装任何中文字体, 返回 None
    """
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

    Returns:
        tuple (alpha0, beta0, mode): α₀(Gy⁻¹), β₀(Gy⁻²), 剂量模式(nominal/actual)
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

def get_z0(beta0, rn, rd):
    '''
    计算得到MK修正的z0
    '''
    return rn**2/(rd*np.sqrt(beta0*(rd**2 + rn**2)))
# ============================================================
# 主流程
# ============================================================
def main():
    """
    脚本入口: 解析参数 -> 读 α₀/β₀ -> 用 ModifiedSMK 建模 ->
    生成存活曲线 -> 叠加 Ac-225 实验点 -> 输出 PNG/CSV/TXT.
    """
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
    ap.add_argument("--z0", type=float, default=40.0,
                    help="饱和参数 z₀ [Gy] (analyze_mk.C 默认 66.0, Inaniwa HSG; PC-3 无公开值)")
    ap.add_argument("--dmax", type=float, default=None,
                    help="存活曲线剂量上限 [Gy]; 默认取数据最大剂量的 1.05 倍")
    ap.add_argument("--outdir", default="result/validation/smk-fit",
                    help="输出目录 (默认 result/validation, 与 analyze_mk.C 一致)")
    ap.add_argument("--fig-mode", default="normal",
                    help="图形模式 (normal: 正常模式, log: 对数模式)")
    args = ap.parse_args()

    # —— 路径准备与存在性校验 ——
    root_path = Path(args.root)
    lq_path = Path(args.lq_csv)
    ac225_path = Path(args.ac225)
    out_dir = Path(args.outdir)
    fig_mode = args.fig_mode
    out_dir.mkdir(parents=True, exist_ok=True)
    for p in (root_path, lq_path, ac225_path):
        if not p.exists():
            sys.exit(f"[ERROR] 找不到输入文件: {p}")

    setup_matplotlib_font()

    # —— 1. 读 α₀, β₀ (来自 X 射线 LQ 拟合) ——
    alpha0, beta0, lq_mode = read_lq_params(lq_path)
    rn = 6.4
    rd = 0.274
    z0 = get_z0(beta0, rn, rd)
    print("=" * 64)
    print("[1] 参考辐射敏感性 (来自 fit_lq_xray.py)")
    print(f"    来源: {lq_path}  (模式: {lq_mode})")
    print(f"    α₀ = {alpha0:.6f} Gy⁻¹   β₀ = {beta0:.6f} Gy⁻²")
    print(f"    计算得到MK修正的 z₀ = {z0:.3f} Gy (由 β₀, r_n={rn} μm, r_d={rd} μm 算出)")
    print(f"    输入饱和参数 z₀ = {args.z0:.3f} Gy")

    # —— 2+3. 用 ModifiedSMK 建模: 从 root 算微剂量学量 + SMK 系数 ——
    try:
        msmk = ModifiedSMK.from_root(root_path, alpha0, beta0, args.z0)
        mk = MK.from_root(root_path, alpha0, beta0)
        smk = SMK.from_root(root_path, alpha0, beta0, args.z0)
        dsmk = DSMK.from_root(root_path, alpha0, beta0, args.z0)
    except ValueError as e:
        sys.exit(f"[ERROR] {e}")
    alpha_s, beta_s = msmk.coefficients()
    print("=" * 64)
    print("[2] 微剂量学量 (来自 microtrack.root, 仅命中事件) + SMK 系数 (Inaniwa 式15-16)")
    print(f"    命中事件: {msmk.n_hits}/{msmk.n_entries}  "
          f"(命中率 {100.0*msmk.n_hits/msmk.n_entries:.1f}%)  z₀ = {args.z0} Gy")
    print(f"    z̄_{{d,F}} = {msmk.z_d_F:.5f} Gy   z̄_{{d,D}} = {msmk.z_d_D:.5f} Gy")
    print(f"    z̄*_{{d,D}}= {msmk.zs_d_D:.5f} Gy   (饱和)")
    print(f"    z̄_{{n,F}} = {msmk.z_n_F:.5g} Gy   z̄_{{n,D}} = {msmk.z_n_D:.5f} Gy")
    print(f"    α_S = α₀ + β₀·z̄*_{{d,D}}            = {alpha_s:.5f} Gy⁻¹")
    print(f"    β_S = β₀·z̄*_{{d,D}}/z̄_{{d,D}}        = {beta_s:.5f} Gy⁻²")

    # —— 4. Ac-225 实验数据 ——
    df_ac = pd.read_csv(ac225_path)
    for col in ("dose_Gy", "survival_fraction"):
        if col not in df_ac.columns:
            sys.exit(f"[ERROR] {ac225_path} 缺少必需列 '{col}'")
    D_obs = df_ac["dose_Gy"].to_numpy(float)
    S_obs = df_ac["survival_fraction"].to_numpy(float)

    # —— 5. SMK 存活曲线 (式24) ——
    Dmax = args.dmax if args.dmax is not None else float(D_obs.max()) * 1.05
    if Dmax <= 0:
        sys.exit(f"[ERROR] --dmax 必须 > 0 (当前 {Dmax})")
    D_grid = np.linspace(0.0, Dmax, 60)
    S_MSMK_grid = np.clip(msmk.survival(D_grid), 1e-12, None)   # 对数轴下限保护
    S_MK_grid = np.clip(mk.survival(D_grid), 1e-12, None)   # 对数轴下限保护
    S_SMK_grid = np.clip(smk.survival(D_grid), 1e-12, None)   # 对数轴下限保护
    S_DSMK_grid = np.clip(dsmk.survival(D_grid), 1e-12, None)   # 对数轴下限保护
    # 在实验剂量点处取模型值, 算 log-RMSE (仅 D>0)
    pos = D_obs > 0
    S_MSMK_at = msmk.survival(D_obs[pos])
    S_MK_at = mk.survival(D_obs[pos])
    S_SMK_at = smk.survival(D_obs[pos])
    S_DSMK_at = dsmk.survival(D_obs[pos])
    lg_obs = np.log(np.clip(S_obs[pos], 1e-30, None))
    lg_MSMK_at = np.log(np.clip(S_MSMK_at, 1e-30, None))
    lg_MK_at = np.log(np.clip(S_MK_at, 1e-30, None))
    lg_SMK_at = np.log(np.clip(S_SMK_at, 1e-30, None))
    lg_DSMK_at = np.log(np.clip(S_DSMK_at, 1e-30, None))
    rmse_MSMK_log = float(np.sqrt(np.mean((lg_obs - lg_MSMK_at) ** 2)))
    rmse_MK_log = float(np.sqrt(np.mean((lg_obs - lg_MK_at) ** 2)))
    rmse_SMK_log = float(np.sqrt(np.mean((lg_obs - lg_SMK_at) ** 2)))
    rmse_DSMK_log = float(np.sqrt(np.mean((lg_obs - lg_DSMK_at) ** 2)))
    print("=" * 64)
    print("[3] 模型 vs Ac-225 实验 (对数空间)")
    print(f"    D 范围: 0–{D_obs.max():.3f} Gy,  数据点 {pos.sum()} 个")
    print(f"    RMSE(ln S) = {rmse_MSMK_log:.5f}")

    # —— 6. 画图 ——
    ymin = np.min(S_obs[pos]) / 2
    fig, ax = plt.subplots(figsize=(8.0, 5.6), dpi=150)
    ax.scatter(D_obs[pos], S_obs[pos], s=66, c="#1f77b4", edgecolor="black",
               label="Ac-225 实验 (PC-3 PIP)", zorder=5)
    ax.plot(D_grid, S_MSMK_grid, color="#d62728", lw=2.2, label="修正 SMK 式")
    ax.plot(D_grid, S_MK_grid, color="#2ca02c", lw=2.2, label="MK 式")
    ax.plot(D_grid, S_SMK_grid, color="#ff7f0e", lw=2.2, label="SMK  式")
    ax.plot(D_grid, S_DSMK_grid, color="#9467bd", lw=2.2, label="DSMK  式")
    ax.axhline(0.1, color="gray", ls=":", lw=1)
    ax.axhline(0.01, color="gray", ls=":", lw=1)

    ax.set_xlabel("Ac-225 吸收剂量 D (Gy)")
    ax.set_ylabel("存活分数 S(D)")
    ax.set_title(f"各模型存活曲线 vs Ac-225 实验数据 (PC-3 PIP, $z_0$={args.z0} Gy)")
    ax.set_xlim(left=0)
    ax.set_ylim(1e-6, 1.1)
    if fig_mode == "log":
        ax.set_yscale("log")
        ax.set_ylim(ymin, 1.1)
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(loc="upper right", framealpha=0.92)

    # txt = (
    #     f"$\\alpha_0$ = {alpha0:.4f} Gy$^{{-1}}$   $\\beta_0$ = {beta0:.4f} Gy$^{{-2}}$  ({lq_mode})\n"
    #     f"$z_0$ = {args.z0} Gy\n"
    #     f"$z^*_{{d,D}}$ = {msmk.zs_d_D:.4f}   $z_{{d,D}}$ = {msmk.z_d_D:.4f}\n"
    #     f"$z_{{n,D}}$ = {msmk.z_n_D:.5f}\n"
    #     f"$\\alpha_S$ = {alpha_s:.4f}   $\\beta_S$ = {beta_s:.5f}\n"
    #     f"RMSE(ln S) = {rmse_MSMK_log:.4f}"
    # )
    # ax.text(0.6, 0.6, txt, transform=ax.transAxes, fontsize=9.5, va="bottom",
    #         bbox=dict(boxstyle="round", facecolor="lightyellow",
    #                   edgecolor="gray", alpha=0.92))
    fig.tight_layout()

    png = out_dir / f"smk_ac225_{lq_mode}_{fig_mode}.png"
    fig.savefig(png, dpi=150)
    plt.close(fig)
    print("=" * 64)
    print(f"[输出] 图像 → {png}")

    # —— 7. 曲线 CSV + 摘要 TXT ——
    curve_csv = out_dir / f"smk_ac225_{lq_mode}_{fig_mode}.csv"
    pd.DataFrame({
        "dose_Gy": D_grid,
        "S_MSMK": S_MSMK_grid,
        "S_SMK": S_SMK_grid,
        "S_DSMK": S_DSMK_grid,
        "S_MK": S_MK_grid,
        "alpha_S": np.full_like(D_grid, alpha_s),
        "beta_S": np.full_like(D_grid, beta_s),
        "z_n_D": np.full_like(D_grid, msmk.z_n_D),
    }).to_csv(curve_csv, index=False)

    txt_path = out_dir / f"smk_ac225_{lq_mode}_{fig_mode}.txt"
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("修正 SMK 存活曲线 vs Ac-225 实验 (PC-3 PIP)\n")
        f.write("=" * 56 + "\n")
        f.write("模型: S(D)=exp(-α_S D - β_S D²)·[1+(D·z̄_{n,D}/2)·((α_S+2β_S D)²-2β_S)]  (Inaniwa 式24)\n")
        f.write(f"α₀, β₀ 来源: {lq_path} ({lq_mode})\n")
        f.write(f"α₀ = {alpha0:.10g} Gy⁻¹\n")
        f.write(f"β₀ = {beta0:.10g} Gy⁻²\n")
        f.write(f"z₀ = {args.z0} Gy\n")
        f.write("-" * 56 + "\n")
        f.write(f"z̄_{{d,F}} = {msmk.z_d_F:.10g} Gy\n")
        f.write(f"z̄_{{d,D}} = {msmk.z_d_D:.10g} Gy\n")
        f.write(f"z̄*_{{d,D}}= {msmk.zs_d_D:.10g} Gy\n")
        f.write(f"z̄_{{n,F}} = {msmk.z_n_F:.10g} Gy\n")
        f.write(f"z̄_{{n,D}} = {msmk.z_n_D:.10g} Gy\n")
        f.write(f"α_S = {alpha_s:.10g} Gy⁻¹\n")
        f.write(f"β_S = {beta_s:.10g} Gy⁻²\n")
        f.write(f"RMSE(ln S, Ac-225 点) = {rmse_MSMK_log:.10g}\n")
        f.write("-" * 56 + "\n")
        f.write("dose_Gy\tobserved_S\tmodel_S\n")
        for d, s in zip(D_obs, S_obs):
            sm = float(msmk.survival(np.array([d]))[0])
            f.write(f"{d:.10g}\t{s:.10g}\t{sm:.10g}\n")
    print(f"[输出] 曲线 → {curve_csv}")
    print(f"[输出] 摘要 → {txt_path}")
    print("=" * 64)


if __name__ == "__main__":
    main()
