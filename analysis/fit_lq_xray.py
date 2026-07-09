#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fit_lq_xray.py — 对 PC-3 PIP 细胞的 X 射线存活分数数据做 LQ 拟合

LQ (Linear-Quadratic) 模型:
    S(D) = exp(-α·D - β·D²)

其中:
    D — 剂量 (Gy)
    S — 存活分数
    α — 线性项系数 (Gy^-1), 低剂量斜率
    β — 二次项系数 (Gy^-2)
    α/β — 表征剂量响应曲线形状, Gy

输入: data/validation/Xray_survival.csv
       (列: nominal_dose_Gy, actual_dose_Gy, survival_fraction)

支持两种剂量模式:
    nominal — 不做 0.9 修正, 直接用实验报告的标称剂量
    actual  — 用 0.9 × nominal 修正后的实际吸收剂量
    both    — 两种都拟合, 各出一张图 (默认)

输出 (默认 outdir=result/validation):
    lq_fit_xray_nominal.csv / .png — nominal 拟合
    lq_fit_xray_actual.csv  / .png — actual 拟合

注: 拟合得到的 α, β 可作为后续 ROOT 脚本 (analyze_mk.C, analyze_dsmk.C)
    中 α₀, β₀ 的近似初始值 (PC-3 PIP)。
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# 注: 不再使用 scipy.optimize.curve_fit — 改为 log-space LSTSQ,
#     与 Inaniwa 2018 / Carlson 2008 等放射生物学文献一致

# ============================================================
# 中文字体设置 (防止图像中中文乱码)
# ============================================================
def setup_matplotlib_font():
    """
    探测系统中可用的中文字体 (优先 Noto CJK, .ttc 多语种文件),
    并写入 matplotlib 全局配置, 防止图像中中文乱码.

    Returns:
        str 或 None: 成功匹配到的字体名; 若系统未安装任何中文字体, 返回 None
    """
    import matplotlib
    from matplotlib import font_manager

    # 关闭 Unicode 负号替换, 避免坐标轴负号显示为方块
    matplotlib.rcParams["axes.unicode_minus"] = False

    # 候选中文字体列表 (按优先级从高到低排列)
    candidates = [
        "Noto Sans CJK JP", "Noto Sans CJK SC", "Noto Sans CJK TC", "Noto Sans CJK KR",
        "Noto Serif CJK JP", "Noto Serif CJK SC",
        "AR PL UMing CN", "AR PL UMing TW MBE", "AR PL UKai CN",
        "WenQuanYi Zen Hei", "WenQuanYi Micro Hei",
        "SimHei", "Microsoft YaHei", "PingFang SC", "Hiragino Sans GB",
    ]

    # —— 逐个尝试候选字体, 命中即设置并返回 ——
    for name in candidates:
        try:
            path = font_manager.findfont(name, fallback_to_default=False)
            if path and "DejaVu" not in path:
                matplotlib.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
                print(f"[matplotlib] 使用字体: {name} ({path})")
                return name
        except Exception:
            continue

    # —— 候选列表均未命中, 尝试从备用路径强制注册 Noto Sans CJK ——
    try:
        font_manager.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
        matplotlib.rcParams["font.sans-serif"] = ["Noto Sans CJK JP", "DejaVu Sans"]
        print("[matplotlib] 强制注册 Noto Sans CJK (备用路径)")
        return "Noto Sans CJK JP"
    except Exception:
        pass

    # —— 系统未安装任何中文字体, 给出警告 ——
    print("[matplotlib] 未找到中文字体, 使用默认 (图像中中文可能为方块)")
    return None


# ============================================================
# LQ 模型
# ============================================================
def lq_model(D, alpha, beta):
    """
    计算 LQ 存活模型 S(D) = exp(-α·D - β·D²).

    Args:
        D: 剂量 (Gy), 标量或数组
        alpha: 线性项系数 α (Gy⁻¹), 低剂量区斜率
        beta: 二次项系数 β (Gy⁻²), 高剂量区弯曲程度

    Returns:
        存活分数 S (无量纲, 取值 0–1)
    """
    return np.exp(-alpha * D - beta * D * D)


def r_squared(y_obs, y_pred):
    """
    计算决定系数 R², 衡量拟合优度.

    Args:
        y_obs: 观测值数组
        y_pred: 模型预测值数组 (与 y_obs 等长)

    Returns:
        R² (1.0 表示完美拟合, 越接近 0 拟合越差)
    """
    ss_res = np.sum((y_obs - y_pred) ** 2)            # 残差平方和
    ss_tot = np.sum((y_obs - np.mean(y_obs)) ** 2)    # 总平方和
    return 1.0 - ss_res / ss_tot


def fit_one_mode(D, S):
    """
    对一组 (D, S) 数据做一次 LQ 拟合.

    采用 **对数空间最小二乘 (log-space LSTSQ)** (与 Inaniwa 2018 / Carlson 2008
    等放射生物学文献一致): 对 y = -ln(S) 拟合 y = α·D + β·D².

    物理依据: LQ 模型假设致死事件是 Poisson 累积过程,
    〈N_lethal〉 = α·D + β·D², S = exp(-〈N_lethal〉),
    因此 -ln(S) 才是线性可加的物理量, 应对它做最小二乘.

    Args:
        D: 剂量数组 (Gy), 可能含零剂量点 (零点会被自动剔除)
        S: 存活分数数组 (无量纲), 与 D 等长

    Returns:
        tuple (alpha, beta, alpha_err, beta_err, r2, rmse_log, n_used):
            alpha, beta — LQ 拟合系数 (Gy⁻¹, Gy⁻²)
            alpha_err, beta_err — 系数标准误
            r2 — 决定系数 (在 -ln(S) 空间)
            rmse_log — 均方根误差 (在 -ln(S) 空间)
            n_used — 实际参与拟合的数据点数 (即 D>0 的点数)
    """
    # —— 数据预处理: 剔除零剂量点 (log(1)=0 不影响拟合但保一致) ——
    mask = D > 0
    D_pos = D[mask]
    S_pos = S[mask]

    # —— 构造对数空间观测向量与设计矩阵 ——
    y = -np.log(S_pos)                          # 对数空间观测量 y = -ln(S)
    X = np.column_stack([D_pos, D_pos ** 2])    # 设计矩阵 X = [D, D²]

    # —— 最小二乘求解 X @ [α, β]^T = y ——
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    alpha, beta = coef

    # —— 参数标准误: σ² = SS_res / (n - p), Cov = σ² (X^T X)^-1 ——
    n, p = X.shape
    y_pred = X @ coef
    ss_res = np.sum((y - y_pred) ** 2)
    sigma2 = ss_res / (n - p)
    try:
        cov = sigma2 * np.linalg.inv(X.T @ X)
        perr = np.sqrt(np.diag(cov))
        alpha_err, beta_err = perr
    except np.linalg.LinAlgError:
        # 设计矩阵奇异 (如数据点不足), 标准误置为 NaN
        alpha_err = beta_err = np.nan

    # —— 拟合优度指标 ——
    r2 = r_squared(y, y_pred)
    rmse_log = float(np.sqrt(ss_res / n))
    return alpha, beta, alpha_err, beta_err, r2, rmse_log, int(mask.sum())


# ============================================================
# 单模式拟合 + 出图
# ============================================================
def run_one(df, dose_col, label, out_dir):
    """
    对指定剂量列做一次完整的 LQ 拟合, 输出 CSV 与 PNG, 并在控制台打印参数.

    Args:
        df: 含存活数据的 DataFrame (需有 dose_col 与 survival_fraction 列)
        dose_col: 本轮使用的剂量列名 (nominal_dose_Gy 或 actual_dose_Gy)
        label: 本轮模式标签 (nominal / actual), 用于输出文件名与图标题
        out_dir: 输出目录 (Path 对象)

    Returns:
        dict: 含 mode/alpha/beta/alpha_err/beta_err/alpha_beta/r2/rmse_log/n_used
    """
    D = df[dose_col].to_numpy()                       # 剂量数组 (Gy)
    S = df["survival_fraction"].to_numpy()            # 存活分数数组

    alpha, beta, alpha_err, beta_err, r2, rmse_log, n_used = fit_one_mode(
        D, S,
    )
    alpha_beta = alpha / beta                         # α/β (Gy), 表征曲线弯曲程度
    ci95_a = 1.96 * alpha_err                         # α 的 95% 置信区间半宽
    ci95_b = 1.96 * beta_err                          # β 的 95% 置信区间半宽

    # —— 控制台输出 ——
    print()
    print("=" * 64)
    print(f"[拟合结果 / {label}]  LQ: S(D) = exp(-α·D - β·D²)   (对数空间 LSTSQ)")
    print("=" * 64)
    print(f"  剂量列: {dose_col}   (数据点 {n_used}/{len(D)})")
    print(f"  α     = {alpha:.6f} ± {alpha_err:.6f} Gy^-1   (95%CI: [{alpha-ci95_a:.6f}, {alpha+ci95_a:.6f}])")
    print(f"  β     = {beta:.6f} ± {beta_err:.6f} Gy^-2   (95%CI: [{beta-ci95_b:.6f}, {beta+ci95_b:.6f}])")
    print(f"  α/β   = {alpha_beta:.4f} Gy")
    print(f"  R²    = {r2:.6f}  (on -ln(S))")
    print(f"  RMSE  = {rmse_log:.6f}  (on -ln(S))")
    print(f"  → ROOT 脚本参数: alpha0={alpha:.4f}, beta0={beta:.4f}")

    # —— 拟合曲线 CSV (稠密采样, 便于后续叠加绘图) ——
    D_dense = np.linspace(0.0, float(D[D > 0].max()) * 1.05, 200)
    S_dense = lq_model(D_dense, alpha, beta)
    fit_csv = out_dir / f"lq_fit_xray_{label}.csv"
    pd.DataFrame({
        "dose_Gy": D_dense,
        "S_LQ": S_dense,
        "alpha": np.full_like(D_dense, alpha),
        "beta": np.full_like(D_dense, beta),
        "mode": [label] * len(D_dense),
    }).to_csv(fit_csv, index=False)
    print(f"  [输出] 拟合曲线 → {fit_csv}")

    # —— 绘图 (对数纵轴) ——
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(D, S, s=70, c="#1f77b4", edgecolor="black",
               label=f"实验数据 ({label} 剂量)", zorder=5)
    ax.plot(D_dense, S_dense, "r-", lw=2.2,
            label=f"LQ 拟合: $\\alpha$={alpha:.4f}, $\\beta$={beta:.4f}")
    # 参考线: 10% 与 1% 存活水平
    ax.axhline(0.1, color="gray", ls=":", lw=1)
    ax.axhline(0.01, color="gray", ls=":", lw=1)
    ax.text(D_dense[-1], 0.11, "S=0.1", ha="right", va="bottom",
            color="gray", fontsize=9)
    ax.text(D_dense[-1], 0.011, "S=0.01", ha="right", va="bottom",
            color="gray", fontsize=9)

    # 模式对应的副标题文字
    title_suffix = (
        "不做 0.9 修正" if label == "nominal"
        else "0.9 × 标称剂量 (实际吸收剂量)" if label == "actual"
        else label
    )
    ax.set_xlabel("剂量 D (Gy)")
    ax.set_ylabel("存活分数 S")
    ax.set_title(f"PC-3 PIP 细胞 X 射线存活分数 — LQ 拟合 ({title_suffix})")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=1e-4)
    ax.set_yscale("log")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend(loc="upper right", framealpha=0.92)

    # 图内参数说明框 (LaTeX 数学标签保留原样)
    txt = (f"$\\alpha$ = {alpha:.4f} ± {alpha_err:.4f} Gy⁻¹\n"
           f"$\\beta$ = {beta:.4f} ± {beta_err:.4f} Gy⁻²\n"
           f"$\\alpha/\\beta$ = {alpha_beta:.2f} Gy\n"
           f"$R^2$ = {r2:.4f} (on -ln S)\n"
           f"RMSE = {rmse_log:.4f} (on -ln S)\n"
           f"模式: {label}")
    # 不要用 family="monospace", 默认等宽字体 (DejaVu Sans Mono) 无中文, 会变方块
    ax.text(0.03, 0.05, txt, transform=ax.transAxes,
            fontsize=10,
            bbox=dict(boxstyle="round", facecolor="lightyellow",
                      edgecolor="gray", alpha=0.9))

    fig.tight_layout()
    fit_png = out_dir / f"lq_fit_xray_{label}.png"
    fig.savefig(fit_png, dpi=140)
    plt.close(fig)
    print(f"  [输出] 拟合图像 → {fit_png}")

    return {
        "mode": label,
        "alpha": alpha, "alpha_err": alpha_err,
        "beta": beta, "beta_err": beta_err,
        "alpha_beta": alpha_beta, "r2": r2,
        "rmse_log": rmse_log, "n_used": n_used,
    }


# ============================================================
# 主流程
# ============================================================
def main():
    """
    脚本入口: 解析命令行参数 -> 读数据 -> 设置字体 -> 按剂量模式拟合 -> 输出对比摘要.
    """
    ap = argparse.ArgumentParser(
        description="PC-3 PIP X 射线存活分数 LQ 拟合 (支持 nominal/actual 两种剂量模式)"
    )
    ap.add_argument(
        "--csv",
        default="data/validation/Xray_survival.csv",
        help="输入 CSV 路径",
    )
    ap.add_argument(
        "--outdir",
        default="result/validation",
        help="输出目录 (默认 result/validation)",
    )
    ap.add_argument(
        "--mode",
        choices=["nominal", "actual", "both"],
        default="both",
        help="剂量模式: nominal=不做0.9修正, actual=用0.9×修正, both=两个都做 (默认)",
    )
    args = ap.parse_args()

    csv_path = Path(args.csv)                    # 输入 CSV 路径
    out_dir = Path(args.outdir)                  # 输出目录
    out_dir.mkdir(parents=True, exist_ok=True)

    # —— 1. 读数据 ——
    if not csv_path.exists():
        sys.exit(f"[ERROR] 找不到输入文件: {csv_path}")
    df = pd.read_csv(csv_path)
    print("=" * 64)
    print(f"[输入] {csv_path}  ({len(df)} 个数据点)")
    print(df.to_string(index=False))

    # —— 2. 字体初始化 ——
    setup_matplotlib_font()

    # —— 3. 按模式跑拟合 ——
    modes = ["nominal", "actual"] if args.mode == "both" else [args.mode]
    results = []
    for mode in modes:
        dose_col = "nominal_dose_Gy" if mode == "nominal" else "actual_dose_Gy"
        res = run_one(df, dose_col, mode, out_dir)
        results.append(res)

    # —— 4. 两种模式都跑的话, 加一个对比摘要 ——
    if len(results) == 2:
        nom, act = results
        print()
        print("=" * 64)
        print("[对比摘要]  0.9 修正前 vs 后  (对数空间 LSTSQ)")
        print("=" * 64)
        print(f"  {'参数':<10} {'nominal (未修正)':<22} {'actual (×0.9)':<22} {'比值':<12}")
        print(f"  {'α':<10} {nom['alpha']:<22.6f} {act['alpha']:<22.6f} "
              f"{act['alpha']/nom['alpha']:<12.4f}  (≈1/0.9 = {1/0.9:.4f})")
        print(f"  {'β':<10} {nom['beta']:<22.6f} {act['beta']:<22.6f} "
              f"{act['beta']/nom['beta']:<12.4f}  (≈1/0.81 = {1/0.81:.4f})")
        print(f"  {'α/β':<10} {nom['alpha_beta']:<22.4f} {act['alpha_beta']:<22.4f} "
              f"{act['alpha_beta']/nom['alpha_beta']:<12.4f}  (≈0.9)")
        print(f"  {'R²':<10} {nom['r2']:<22.6f} {act['r2']:<22.6f}")
        print()
        print("[物理解释] (对数空间 LSTSQ 下)")
        print("  - D_actual = 0.9·D_nominal, 代入 y = -ln(S) = α·D + β·D² 得:")
        print("    α_actual = α_nominal / 0.9")
        print("    β_actual = β_nominal / 0.81")
        print("    (α/β)_actual = 0.9·(α/β)_nominal  ← 注意: 在对数空间下 α/β 不是不变量")
        print()
        print("[推荐]")
        print(f"  → ROOT 脚本请用 actual 模式拟合参数: alpha0={act['alpha']:.4f}, "
              f"beta0={act['beta']:.4f}")

    print("=" * 64)


if __name__ == "__main__":
    main()
