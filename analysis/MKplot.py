#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MKplot.py — Ac-225 vs Lu-177 存活曲线叠加 + RBE 曲线 (批量 + 单文件两模式)

两幅图 (画在同一张 PNG 上):
    上: 两个模型的存活曲线 S(D) vs D (ModifiedSMK Inaniwa2018 式24)
    下: Ac-225 的 RBE(S) = D_Lu177(S) / D_Ac225(S) vs Ac-225 剂量 D_Ac225

模式 1: 单文件 (向后兼容)
    python analysis/MKplot.py --root-ac225 <path> --root-lu177 <path>
模式 2: 批量扫描文件夹
    python analysis/MKplot.py --batch \
        --scan-dirs data/PC-ac225-phy-decay data/V79-ac225-single-decay \
        --branch alpha_events --outdir result/RBE
    (branch 选择 alpha_events 时, lu177 自动切到 beta_events;
     branch 选择 single_events 时, lu177 保持 single_events)

产物:
    单文件模式: <outdir>/<outname>.png + .csv
    批量模式: <outdir>/<cell系>/<input前缀>__<branch>.png + .csv
"""

import sys
import re
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
# 细胞系参数表 (Table 2, mSMKM)
# ============================================================
CELL_PARAMS = {
    "V79":   {"alpha0": 0.127, "beta0": 0.014, "z0": 144.7},
    "HSG":   {"alpha0": 0.183, "beta0": 0.051, "z0":  66.1},
    "Renca": {"alpha0": 0.019, "beta0": 0.019, "z0":  66.0},
    "PC":    {"alpha0": 0.2517, "beta0": 0.03225, "z0": 66},  # V79 同参数
}


def get_cell_params(cell):
    """按细胞系名查 α₀, β₀, z₀; 未找到返回 None."""
    return CELL_PARAMS.get(cell)


# ============================================================
# 配对文件查找
# ============================================================
def find_lu177_for_ac225(ac225_path: Path, scan_dirs):
    """
    在 scan_dirs 中找与 ac225_path 同 compartment 同 cell 系 的 lu177 文件。
    匹配规则: 文件名 ac225_<prefix>_<compartment>_<cell>.root
       -> lu177_phy_decay_<compartment>_<cell>.root
    例: ac225_phy_decay_Nucleus_PC.root   -> lu177_phy_decay_Nucleus_PC.root
       ac225_single_decay_CellExceptNucleus_V79.root -> lu177_phy_decay_CellExceptNucleus_V79.root
    prefix 是固定的源型标记 (phy_decay / single_decay), compartment 可能是复合名.
    """
    name = ac225_path.name
    # 已知 prefix 集合 (源型描述符, 来源 PrimaryGeneratorAction)
    KNOWN_PREFIXES = ("phy_decay", "single_decay", "phy_phy_decay",
                       "single_decay_phy")  # 兜底, 实际只需前两个

    stem = ac225_path.stem  # 不带 .root
    if not stem.startswith("ac225"):
        return None
    rest = stem[len("ac225_"):]  # e.g. 'phy_decay_Nucleus_PC'

    # 尝试去掉已知 prefix
    for prefix in KNOWN_PREFIXES:
        if rest.startswith(prefix + "_"):
            rest = rest[len(prefix) + 1:]  # 去掉 'prefix_'
            break

    # rest 现在应该是 '<compartment>_<cell>', cell 在末尾
    # cell 是单段字母数字 (PC, V79, V19, ...)
    m = re.match(r"^(?P<compartment>.+)_(?P<cell>[A-Za-z0-9]+)$", rest)
    if not m:
        return None
    compartment = m.group("compartment")
    cell = m.group("cell")
    # lu177 命名规则: lu177_phy_decay_<compartment>_<cell>.root
    lu177_name = f"lu177_phy_decay_{compartment}_{cell}.root"
    for d in scan_dirs:
        cand = d / lu177_name
        if cand.exists():
            return cand
    return None


def detect_cell_system(ac225_path: Path):
    """从文件名提取 cell 系 (PC, V79, ...), 用于输出子目录划分."""
    name = ac225_path.name
    # 末尾 _<cell>.root, cell 是字母数字组合 (PC, V79, V19, ...)
    m = re.search(r"_(?P<cell>[A-Za-z0-9]+)\.root$", name)
    if m:
        return m.group("cell")
    return ac225_path.stem.split("_")[-1]


# ============================================================
# RBE 反演: 给定 S 求 D
# ============================================================
def find_dose_at_survival(model, S_target, dose_grid):
    """
    对一个 MKM 模型, 在 dose_grid 上算 S_grid, 然后反查指定存活率 S_target 的剂量。
    survival 是 D 的单调递减函数; 用 np.interp 时把序列翻转以满足单调递增假设。
    """
    S_grid = model.survival(dose_grid)
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
    """
    S_test_grid = (model_test.survival(dose_grid)
                   if survival_test is None else survival_test)
    S_ref_grid  = (model_ref.survival(dose_grid)
                   if survival_ref is None else survival_ref)

    s_lo = max(s_min, float(np.max([S_test_grid.min(), S_ref_grid.min()])))
    s_hi = min(s_max, float(np.min([S_test_grid.max(), S_ref_grid.max()])))
    if s_hi <= s_lo:
        return np.array([]), np.array([]), np.array([])

    S_targets = np.linspace(s_hi, s_lo, n_points)

    D_test_arr, D_ref_arr, RBE_arr = [], [], []
    for S_target in S_targets:
        d_t = float(np.interp(S_target, S_test_grid[::-1], dose_grid[::-1]))
        d_r = float(np.interp(S_target, S_ref_grid[::-1],  dose_grid[::-1]))
        if d_t > 0 and np.isfinite(d_t) and np.isfinite(d_r):
            D_test_arr.append(d_t)
            D_ref_arr.append(d_r)
            RBE_arr.append(d_r / d_t)
    return np.array(D_test_arr), np.array(D_ref_arr), np.array(RBE_arr)


# ============================================================
# 单文件计算 + 绘图
# ============================================================
def compute_and_plot(ac225_path, lu177_path, branch, alpha0, beta0, z0,
                     outdir, outname, npoints_D, npoints_S,
                     smin, smax, dmax_override, fig_mode,
                     label_prefix="", lu177_branch=None):
    """对一对 (Ac225, Lu177) 文件计算存活曲线 + RBE, 写 PNG + CSV.
    lu177_branch: Lu-177 用的 ntuple 树名; None 时与 branch 相同.
    """
    if lu177_branch is None:
        lu177_branch = branch
    print(f"\n{'=' * 70}")
    print(f"[compute] Ac-225: {ac225_path}  (branch={branch})")
    print(f"[compute] Lu-177: {lu177_path}  (branch={lu177_branch})")
    print(f"{'=' * 70}")

    print(f"[MKplot] 加载 {lu177_path} (Lu-177, branch={lu177_branch}) ...")
    msmk_lu177 = ModifiedSMK.from_root(str(lu177_path), lu177_branch, alpha0, beta0, z0)
    print(f"[MKplot] 加载 {ac225_path} (Ac-225, branch={branch}) ...")
    msmk_ac225 = ModifiedSMK.from_root(str(ac225_path), branch, alpha0, beta0, z0)

    print(f"  Lu-177 α_S = {msmk_lu177.alpha_S:.4g} Gy⁻¹")
    print(f"  Ac-225 α_S = {msmk_ac225.alpha_S:.4g} Gy⁻¹")
    print(f"  RBE_M = {msmk_ac225.alpha_S / msmk_lu177.alpha_S:.3f}")

    # 自动剂量范围
    if dmax_override is not None:
        dmax = dmax_override
    else:
        dmax = max(6.0 / msmk_ac225.alpha_S, 6.0 / msmk_lu177.alpha_S)
    print(f"[MKplot] 剂量范围: 0 ~ {dmax:.2f} Gy")

    D_grid = np.linspace(0.0, dmax, npoints_D)
    S_lu = msmk_lu177.survival(D_grid)
    S_ac = msmk_ac225.survival(D_grid)
    D_test, D_ref, RBE = compute_rbe_curve(
        msmk_ac225, msmk_lu177, D_grid,
        survival_test=S_ac, survival_ref=S_lu,
        s_min=smin, s_max=smax, n_points=npoints_S,
    )

    print(f"[MKplot] RBE 有效点: {len(RBE)}  (S ∈ [{smax:.3f}, {smin:.3f}])")
    if len(RBE) > 0:
        print(f"[MKplot] RBE 范围: {RBE.min():.3f} ~ {RBE.max():.3f}")

    # —— 画双面板图 ——
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(9, 9), dpi=130,
        gridspec_kw={"height_ratios": [2, 1]}, sharex=False,
    )

    ax1.plot(D_grid, S_lu, color="tab:blue",   lw=2.0,
             label=f"Lu-177 ({lu177_branch})")
    ax1.plot(D_grid, S_ac, color="tab:orange", lw=2.0,
             label=f"Ac-225 ({branch})")
    ax1.set_ylabel("存活分数 S(D)")
    title_label = label_prefix if label_prefix else ""
    ax1.set_title(
        f"ModifiedSMK{title_label}  (α₀={alpha0}, β₀={beta0}, z₀={z0}, "
        f"Ac-225: {branch} / Lu-177: {lu177_branch})"
    )
    ax1.legend(loc="best", framealpha=0.92)
    ax1.grid(True, which="both", ls="--", alpha=0.4)
    if fig_mode == "log":
        ax1.set_yscale("log")
    if len(RBE) >= 5:
        mid = len(RBE) // 2
        S_at_d = float(np.interp(D_test[mid], D_grid[::-1], S_ac[::-1]))
        ax1.annotate(
            f"RBE≈{RBE[mid]:.2f}",
            xy=(D_test[mid], S_at_d),
            xytext=(D_test[mid] * 0.5, 0.30),
            fontsize=10, color="tab:red",
            arrowprops=dict(arrowstyle="->", color="gray", lw=1, alpha=0.6),
        )

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
    out_dir = Path(outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / f"{outname}.png"
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    print(f"[MKplot] PNG → {png_path}")

    # —— 同时写 RBE 数值表到 CSV ——
    if len(RBE) > 0:
        csv_path = out_dir / f"{outname}.csv"
        import csv as csv_mod
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv_mod.writer(f)
            w.writerow(["survival_S", "D_Ac225_Gy", "D_Lu177_Gy", "RBE_Ac225_vs_Lu177",
                        "ac225_file", "lu177_file", "ac225_branch", "lu177_branch"])
            S_targets_arr = np.linspace(smax, smin, len(RBE))
            for S_t, d_t, d_r, rbe_v in zip(S_targets_arr, D_test, D_ref, RBE):
                w.writerow([f"{S_t:.6f}", f"{d_t:.6f}",
                            f"{d_r:.6f}", f"{rbe_v:.6f}",
                            ac225_path.name, lu177_path.name, branch, lu177_branch])
        print(f"[MKplot] CSV → {csv_path}")

    return len(RBE), (float(RBE.min()) if len(RBE) else None,
                     float(RBE.max()) if len(RBE) else None)


# ============================================================
# 批量模式
# ============================================================
def batch_scan(scan_dirs, branch, outdir, alpha0, beta0, z0,
               npoints_D, npoints_S, smin, smax, dmax_override, fig_mode):
    """
    扫描 scan_dirs 下所有 ac225_*.root 文件, 自动配对 lu177, 计算 RBE.
    不同 cell 系 (PC / V79 / ...) 输出到 outdir/<cell系>/ 子目录.
    """
    print(f"[batch] 扫描目录: {scan_dirs}")
    print(f"[batch] branch: {branch} (lu177 自动切到 beta_events 当 branch=alpha_events)")

    # 配对规则: 当 branch == 'alpha_events', lu177 用 'beta_events'
    lu177_branch = "beta_events" if branch == "alpha_events" else branch
    print(f"[batch] Lu-177 实际 branch: {lu177_branch}")

    # 收集所有 ac225_*.root (排除 lu177)
    ac225_files = []
    for d in scan_dirs:
        d = Path(d)
        if not d.is_dir():
            print(f"[batch] WARN: {d} 不是目录, 跳过")
            continue
        for f in sorted(d.glob("ac225_*.root")):
            ac225_files.append(f)
    if not ac225_files:
        print(f"[batch] ERR: 没找到任何 ac225_*.root")
        return

    print(f"[batch] 找到 {len(ac225_files)} 个 ac225_*.root:")
    for f in ac225_files:
        print(f"  - {f}")

    results = []
    for ac225_path in ac225_files:
        lu177_path = find_lu177_for_ac225(ac225_path, scan_dirs)
        if lu177_path is None:
            print(f"\n[batch] SKIP: {ac225_path.name} (找不到匹配的 lu177)")
            continue

        cell = detect_cell_system(ac225_path)
        cell_outdir = Path(outdir) / cell
        cell_outdir.mkdir(parents=True, exist_ok=True)

        # 输出文件名前缀: 保留输入文件名 (不带 cell 系后缀, 因为 cell 系已体现在子目录)
        # 例: ac225_single_decay_Nucleus_PC.root -> ac225_single_decay_Nucleus
        stem = ac225_path.stem
        # 去掉 cell 系后缀 (最后一个 _<cell>)
        stem_no_cell = re.sub(r"_[A-Za-z0-9]+$", "", stem)
        outname = f"{stem_no_cell}__{branch}"

        try:
            n_pts, rbe_range = compute_and_plot(
                ac225_path=ac225_path,
                lu177_path=lu177_path,
                branch=branch,           # Ac-225 用主 branch (e.g. alpha_events)
                lu177_branch=lu177_branch,  # Lu-177 用映射后的 branch (e.g. beta_events)
                alpha0=alpha0, beta0=beta0, z0=z0,
                outdir=cell_outdir, outname=outname,
                npoints_D=npoints_D, npoints_S=npoints_S,
                smin=smin, smax=smax, dmax_override=dmax_override,
                fig_mode=fig_mode,
                label_prefix=f" — {cell}",
            )
            results.append({
                "cell": cell,
                "ac225": ac225_path.name,
                "lu177": lu177_path.name,
                "out_png": str(cell_outdir / f"{outname}.png"),
                "out_csv": str(cell_outdir / f"{outname}.csv"),
                "n_points": n_pts,
                "rbe_min": rbe_range[0],
                "rbe_max": rbe_range[1],
            })
        except Exception as e:
            print(f"[batch] ERR 处理 {ac225_path.name}: {e}")
            import traceback; traceback.print_exc()
            continue

    # 汇总
    print("\n" + "=" * 70)
    print("[batch] 处理汇总:")
    print("=" * 70)
    print(f"{'Cell':<8} {'Ac-225 file':<45} {'N_RBE':>7} {'RBE 范围':<20}")
    print("-" * 70)
    for r in results:
        rng = f"[{r['rbe_min']:.2f}, {r['rbe_max']:.2f}]" if r['rbe_min'] is not None else "N/A"
        print(f"{r['cell']:<8} {r['ac225'][:43]:<45} {r['n_points']:>7d} {rng:<20}")
    print(f"\n[batch] 共处理 {len(results)} 对 (Ac-225, Lu-177)")
    print(f"[batch] 输出根目录: {outdir}/<cell系>/")
    print(f"[batch] 完成.")


# ============================================================
# 命令行参数
# ============================================================
def params_args():
    ap = argparse.ArgumentParser(
        description="Ac-225 vs Lu-177 存活曲线 + Ac-225 RBE 曲线 (Modified SMK), "
                    "支持单文件或批量扫描")
    # —— 细胞类型: 输入即自动读 α₀/β₀/z₀ (Table 2) ——
    ap.add_argument("--cell", default=None,
                    choices=list(CELL_PARAMS.keys()),
                    help="细胞类型; 指定后自动读 α₀/β₀/z₀ (V79/HSG/Renca). "
                         "也可单独用 --alpha0/--beta0/--z0 覆盖.")
    # 通用参数
    ap.add_argument("--branch", default="single_events",
                    choices=["events", "single_events", "alpha_events", "beta_events"],
                    help="microtrack ntuple 树名 (默认 single_events). "
                         "当 batch 模式下选 alpha_events 时, lu177 自动切到 beta_events.")
    ap.add_argument("--alpha0", type=float, default=None,
                    help="参考 α₀ [Gy⁻¹]; 不指定则用 --cell 查表或默认 0.127")
    ap.add_argument("--beta0", type=float, default=None,
                    help="参考 β₀ [Gy⁻²]; 不指定则用 --cell 查表或默认 0.014")
    ap.add_argument("--z0", type=float, default=None,
                    help="饱和参数 z₀ [Gy]; 不指定则用 --cell 查表或默认 144.7")
    ap.add_argument("--dmax", type=float, default=None,
                    help="存活/RBE 曲线剂量上限 [Gy]; 默认取两份数据最大值的 1.05 倍")
    ap.add_argument("--smin", type=float, default=0.01,
                    help="RBE 计算所用最小存活率 (默认 0.01)")
    ap.add_argument("--smax", type=float, default=0.99,
                    help="RBE 计算所用最大存活率 (默认 0.99)")
    ap.add_argument("--npoints-D", type=int, default=1000,
                    help="剂量网格点数 (默认 1000)")
    ap.add_argument("--npoints-S", type=int, default=100,
                    help="RBE 计算所用存活率点数 (默认 100)")
    ap.add_argument("--fig-mode", default="normal", choices=["normal", "log"],
                    help="存活曲线 y 轴模式 (默认 normal)")

    # 单文件模式参数
    ap.add_argument("--root-ac225", default=None,
                    help="Ac-225 .root 路径 (单文件模式)")
    ap.add_argument("--root-lu177", default=None,
                    help="Lu-177 .root 路径 (单文件模式, 可省略则自动按 compartment 配对)")

    # 批量模式参数
    ap.add_argument("--batch", action="store_true",
                    help="批量扫描模式: 遍历 --scan-dirs 中所有 ac225_*.root, 自动配对 lu177")
    ap.add_argument("--scan-dirs", nargs="+", default=[],
                    help="批量模式扫描目录列表 (e.g. data/PC-ac225-phy-decay data/V79-ac225-single-decay)")
    ap.add_argument("--outdir", default="result/RBE",
                    help="输出根目录 (单文件模式直接放这里, 批量模式放 <outdir>/<cell系>/ 下)")

    args = ap.parse_args()
    return args


# ============================================================
# 主流程
# ============================================================
def main():
    args = params_args()

    # —— 解析细胞参数: --cell 查表 → --alpha0 等覆盖 → 默认值 ——
    defaults = {"alpha0": 0.127, "beta0": 0.014, "z0": 144.7}  # V79 默认
    if args.cell:
        cp = get_cell_params(args.cell)
        if cp:
            defaults.update(cp)
            print(f"[cell] {args.cell}: α₀={cp['alpha0']}, β₀={cp['beta0']}, z₀={cp['z0']}")
    args.alpha0 = args.alpha0 if args.alpha0 is not None else defaults["alpha0"]
    args.beta0  = args.beta0  if args.beta0  is not None else defaults["beta0"]
    args.z0     = args.z0     if args.z0     is not None else defaults["z0"]
    print(f"[params] α₀={args.alpha0}, β₀={args.beta0}, z₀={args.z0}")

    if args.batch:
        if not args.scan_dirs:
            print("[MKplot] ERR: 批量模式必须指定 --scan-dirs")
            return
        batch_scan(
            scan_dirs=[Path(d) for d in args.scan_dirs],
            branch=args.branch,
            outdir=args.outdir,
            alpha0=args.alpha0, beta0=args.beta0, z0=args.z0,
            npoints_D=args.npoints_D, npoints_S=args.npoints_S,
            smin=args.smin, smax=args.smax, dmax_override=args.dmax,
            fig_mode=args.fig_mode,
        )
        return

    # 单文件模式 (向后兼容)
    if args.root_ac225 is None:
        print("[MKplot] ERR: 单文件模式必须指定 --root-ac225")
        return

    ac225_path = Path(args.root_ac225)
    if not ac225_path.exists():
        print(f"[MKplot] ERR: {ac225_path} 不存在")
        return

    # Lu-177 路径
    if args.root_lu177:
        lu177_path = Path(args.root_lu177)
    else:
        # 尝试从 ac225 路径所在目录自动配对
        lu177_path = find_lu177_for_ac225(ac225_path, [ac225_path.parent])
        if lu177_path is None:
            print(f"[MKplot] ERR: 无法自动配对 lu177, 请显式指定 --root-lu177")
            return

    # 单文件模式输出文件名: <outname>.png
    stem = ac225_path.stem
    stem_no_cell = re.sub(r"_[A-Za-z0-9]+$", "", stem)
    outname = f"{stem_no_cell}__{args.branch}"

    # 单文件模式也支持 branch 映射 (alpha_events -> beta_events for Lu-177)
    lu177_branch = "beta_events" if args.branch == "alpha_events" else args.branch

    compute_and_plot(
        ac225_path=ac225_path,
        lu177_path=lu177_path,
        branch=args.branch,
        lu177_branch=lu177_branch,
        alpha0=args.alpha0, beta0=args.beta0, z0=args.z0,
        outdir=args.outdir, outname=outname,
        npoints_D=args.npoints_D, npoints_S=args.npoints_S,
        smin=args.smin, smax=args.smax, dmax_override=args.dmax,
        fig_mode=args.fig_mode,
    )
    print("[MKplot] 完成。")


if __name__ == "__main__":
    main()