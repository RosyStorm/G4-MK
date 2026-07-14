import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from MKmodel.msmk import ModifiedSMK
import argparse


def params_args():
    ap = argparse.ArgumentParser(
        description="修正 SMK 存活曲线 (α₀/β₀ 取自 X 射线 LQ 拟合) × Ac-225 实验数据叠加图"
    )
    ap.add_argument("--root", default="data/lu177_phy_decay_Membrane.root",
                    help="microtrack 模拟 ntuple 路径 (默认 data/ac225_phy_decay_Membrane.root)")
    ap.add_argument("--tree", default="single_events", help="microtrack 模拟 ntuple 树名 (默认 events)")
    ap.add_argument("--lq-csv", default="result/validation/lq_fit_xray_actual.csv",
                    help="fit_lq_xray.py 输出的拟合曲线 CSV (含 alpha,beta 列); "
                         "默认 actual 模式, 可换 lq_fit_xray_nominal.csv")
    ap.add_argument("--ac225", default="data/validation/Ac225_survival.csv",
                    help="Ac-225 实验存活数据 CSV (默认 data/validation/Ac225_survival.csv)")
    ap.add_argument("--z0", type=float, default=66.0, help="饱和参数 z₀ [Gy] (analyze_mk.C 默认 66.0, Inaniwa HSG; PC-3 无公开值)")
    ap.add_argument("--dmax", type=float, default=None, help="存活曲线剂量上限 [Gy]; 默认取数据最大剂量的 1.05 倍")
    ap.add_argument("--outdir", default="result/",
                    help="输出目录 (默认 result/, 与 analyze_mk.C 一致)")
    ap.add_argument("--fig-mode", default="normal", help="图形模式 (normal: 正常模式, log: 对数模式)")
    args = ap.parse_args()
    return args

def main():
    args = params_args()
    root_path_ac225 = Path("data/ac225_phy_decay_Membrane.root")
    root_path_lu177 = Path("data/lu177_phy_decay_Membrane.root")
    tree_name = args.tree
    alpha0 = 0.226
    beta0 = 0.026
    msmk_lu117 = ModifiedSMK.from_root(root_path_lu177, tree_name, alpha0, beta0, args.z0)
    msmk_ac225 = ModifiedSMK.from_root(root_path_ac225, tree_name, alpha0, beta0, args.z0)
    fig, ax = msmk_lu117.plot_survival(ax=None, color="blue", label="Lu-177")
    fig, ax = msmk_ac225.plot_survival(ax=ax, color="orange", label="Ac-225")
    ax.set_title(f"Modified SMK")
    fig.savefig(Path(args.outdir) / "mks_fit.png", dpi=300)

if __name__ == "__main__":
    main()