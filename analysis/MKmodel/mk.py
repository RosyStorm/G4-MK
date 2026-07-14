#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mk.py — MKM 基类 (经典 Microdosimetric Kinetic Model)

本类是 SMK / ModifiedSMK / DSMK 的**父类**。
- 内置经典 Sato 2012 式 (6) LQ 闭式 (`α_S = α₀ + (z̄_{d,D}-z̄_{n,D})·β₀`, `β_S = β₀`)
- 子类覆写 `_compute_alpha_beta()` 与 `survival(D)` 即可获得新模型
- 共用 from_root 工厂、summary / microdosimetry / coefficients 等

参考: Hawkins 1994 (微观剂量学);
       Sato & Furusawa 2012 式 (6)–(11) (经典 MK 闭式)

构造参数 (基础 MK):
    alpha0  — 参考 α₀ (Gy⁻¹) 一般 X 射线 LQ 拟合得到
    beta0   — 参考 β₀ (Gy⁻²)
    z_d_D   — 域剂量均 z̄_{d,D} (Gy)
    z_n_D   — 核剂量均 z̄_{n,D} (Gy)
    z_d_F   — 域频率均 z̄_{d,F} (Gy, 可选, 仅记录)
    z_n_F   — 核频率均 z̄_{n,F} (Gy, 可选, 仅记录)

用法:
    mk = MK(alpha0=0.174, beta0=0.0568, z_d_D=54.1, z_n_D=0.107)
    print(mk)            # α_S, β_S, microdosimetry 一览
    S = mk.survival(np.linspace(0, 10, 200))
"""

import numpy as np


# ============================================================
# 模块级 helper: 从 microtrack.root 抽取基础微剂量学量
# ============================================================
def load_microdosimetry_from_root(root_path, tree_name):
    """
    从 microtrack.root 的指定树读取 (z_d_Gy, z_n_Gy, weight) 三列,
    限定 hitFlag==1。

    :param root_path: ROOT 文件路径
    :param tree_name: 树名 (例如 "events" / "single_events")
    :return: (z_d, z_n, w, n_hits, n_entries)
        z_d       — 域比能数组 (Gy), shape=(n_hits,)
        z_n       — 核比能数组 (Gy), shape=(n_hits,)
        w         — 加权数组, shape=(n_hits,)
        n_hits    — 命中事件数
        n_entries — ntuple 总条目数
    :raises ValueError: 树为空或无命中事件
    """
    import ROOT  # 延迟导入: 避免 import mk 时强制依赖 ROOT

    rdf = ROOT.RDataFrame(str(tree_name), str(root_path))
    n_entries = int(rdf.Count().GetValue())
    if n_entries == 0:
        raise ValueError(f"{root_path} 内 {tree_name} ntuple 为空")

    hits = rdf.Filter("hitFlag==1")
    n_hits = int(hits.Count().GetValue())
    if n_hits == 0:
        raise ValueError(f"{root_path} 无命中事件 (hitFlag==1), 无法算微剂量学量")

    cols = hits.AsNumpy(["z_d_Gy", "z_n_Gy", "weight"])
    z_d = np.asarray(cols["z_d_Gy"], dtype=float)
    z_n = np.asarray(cols["z_n_Gy"], dtype=float)
    w   = np.asarray(cols["weight"], dtype=float)

    return z_d, z_n, w, n_hits, n_entries


def compute_basic_microdosimetry(z_d, z_n, w, n_hits):
    """
    从域比能 / 核比能 / 权 / 命中数算 4 个基础微剂量学量。

    :return: (z_d_F, z_d_D, z_n_F, z_n_D)
    """
    # 域级(带权)
    Sw   = w.sum()                       # Σw
    Swz  = (w * z_d).sum()               # Σ(w·z_d)
    Swzz = (w * z_d * z_d).sum()         # Σ(w·z_d²)
    if Sw <= 0 or Swz <= 0:
        raise ValueError("域级求和异常 (Σw 或 Σ(w·z_d) ≤ 0)")
    z_d_F = Swz / Sw                     # 式 (10)
    z_d_D = Swzz / Swz                   # 式 (11)

    # 核级(无权)
    Szn = z_n.sum()
    if Szn <= 0:
        raise ValueError("核级求和异常 (Σz_n ≤ 0)")
    z_n_F = Szn / n_hits                 # z̄_{n,F}
    z_n_D = (z_n * z_n).sum() / Szn      # 式 (7)

    return z_d_F, z_d_D, z_n_F, z_n_D


# ============================================================
# 父类: 经典 MK (Hawkins 1994 / Sato 2012 闭式)
# ============================================================
class MK:
    """
    经典 MKM 基类 (Hawkins 1994 / Sato 2012 闭式)。

    子类 (SMK / ModifiedSMK / DSMK) 通过覆写
      - `_compute_alpha_beta()` 设置 α_S / β_S
      - `survival(D, ...)`         提供存活曲线公式
    即可实现新模型。父类统一处理:
      - 字段存储与校验
      - from_root 工厂
      - 微剂量学量 summary / microdosimetry
      - 共享 saturation / coefficients / __repr__
    """

    # ----- 构造 -----
    def __init__(self, alpha0, beta0, z_d_D, z_n_D,
                 z_d_F=None, z_n_F=None):
        """
        :param alpha0:  参考线性项 α₀ (Gy⁻¹)
        :param beta0:   参考二次项 β₀ (Gy⁻²)
        :param z_d_D:   域剂量均 z̄_{d,D} (Gy), 必须 > 0
        :param z_n_D:   核剂量均 z̄_{n,D} (Gy)
        :param z_d_F:   域频率均 z̄_{d,F} (Gy, 可选)
        :param z_n_F:   核频率均 z̄_{n,F} (Gy, 可选)
        :raises ValueError: z_d_D ≤ 0
        """
        if z_d_D <= 0:
            raise ValueError(f"z_d_D 必须 > 0 (β_S 分母), 当前 {z_d_D}")

        # 细胞系 / 参考辐射
        self.alpha0 = float(alpha0)
        self.beta0  = float(beta0)

        # 单事件基础微剂量学量
        self.z_d_F = None if z_d_F is None else float(z_d_F)
        self.z_d_D = float(z_d_D)
        self.z_n_F = None if z_n_F is None else float(z_n_F)
        self.z_n_D = float(z_n_D)

        # 统计来源 (from_root 时填)
        self.n_hits    = None
        self.n_entries = None

        # 调用子类钩子, 设置 α_S / β_S
        self._compute_alpha_beta()

    # ----- 子类钩子 -----
    def _compute_alpha_beta(self):
        """
        由子类覆写: 计算并设置 self.alpha_S / self.beta_S。
        基类 (经典 MK) 实现 Sato 2012 式 (6):
            α_S = α₀ + (z̄_{d,D} - z̄_{n,D}) · β₀
            β_S = β₀
        """
        self.alpha_S = self.alpha0 + (self.z_d_D - self.z_n_D) * self.beta0
        self.beta_S  = self.beta0

    # ----- 存活曲线 (默认: 经典 LQ 闭式) -----
    def survival(self, D):
        """
        经典 LQ 闭式存活分数:
            S(D) = exp(-α_S·D - β_S·D²)
        子类根据需要覆写 (例如 ModifiedSMK 加 z_n bracket 修正,
        SMK 用 MC, DSMK 用双 MC)。

        :param D: 剂量 (Gy), 标量或数组
        :return: 存活分数 (无量纲, 与 D 同形)
        """
        D = np.asarray(D, dtype=float)
        return np.exp(-self.alpha_S * D - self.beta_S * D ** 2)

    # ----- 绘图 (子类可直接调用: 与子类的 survival() 一致) -----
    def plot_survival(self, D=None, ax=None, mode="normal", **kwargs):
        """
        绘制存活曲线 S(D) vs D (Gy), 使用 self.survival(D)。
        子类覆盖 survival() 后, 调用 plot_survival() 仍能反映子类的模型公式。

        :param D:     剂量数组 (Gy); 默认 np.linspace(0, 10, 200)
        :param ax:    matplotlib Axes; 默认 plt.subplots() 新建
        :param mode:  'normal' 或 'log'
        :param kwargs: 传给 ax.plot (color/label/linestyle 等)
        :return: (fig, ax)
        """
        import matplotlib.pyplot as plt

        if D is None:
            D = np.linspace(0, 10, 200)
        S = self.survival(D)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        ax.plot(D, S, **kwargs)
        if mode == "log":
            ax.set_yscale("log")
        ax.set_xlabel("Dose D (Gy)")
        ax.set_ylabel("Survival fraction S(D)")
        ax.grid(True, which="both", ls="--", lw=0.5)

        # 仅在已有 label 时显示 legend
        if ax.get_legend_handles_labels()[0]:
            ax.legend()
        return fig, ax

    # ----- from_root 工厂 -----
    @classmethod
    def from_root(cls, root_path, tree_name, alpha0, beta0):
        """
        从 microtrack.root 读取命中事件, 计算基础微剂量学量, 构造 cls 实例。
        子类如需饱和 / 域卷积等额外步骤, 请覆写本方法 (例: SMK 加 z0 算 z̄*_{d,D})。

        :param root_path, tree_name: ROOT 文件与树名
        :param alpha0, beta0:        参考辐射 LQ 系数
        :return: cls 实例 (附带 .n_hits / .n_entries)
        """
        z_d, z_n, w, n_hits, n_entries = load_microdosimetry_from_root(
            root_path, tree_name)
        z_d_F, z_d_D, z_n_F, z_n_D = compute_basic_microdosimetry(
            z_d, z_n, w, n_hits)

        obj = cls(alpha0, beta0, z_d_D, z_n_D,
                  z_d_F=z_d_F, z_n_F=z_n_F)
        obj.n_hits    = n_hits
        obj.n_entries = n_entries
        return obj

    # ----- 静态: 饱和修正 (Sato 式 16 sqrt 形式) -----
    @staticmethod
    def saturation(z_d, z0):
        """
        饱和比能 z'_d = z₀ · √(1 − exp(−(z_d/z₀)²))。
        Sato 2012 式 (16), 不带 sqrt 仅作别名, 子类与本类一致地调用。
        """
        z_d = np.asarray(z_d, dtype=float)
        return z0 * np.sqrt(np.maximum(1.0 - np.exp(-(z_d / z0) ** 2), 0.0))

    # ----- 取值辅助 (共用) -----
    def coefficients(self):
        """返回 (α_S, β_S)"""
        return self.alpha_S, self.beta_S

    def microdosimetry(self):
        """返回基础微剂量学量 dict"""
        return {
            "z_d_F": self.z_d_F, "z_d_D": self.z_d_D,
            "z_n_F": self.z_n_F, "z_n_D": self.z_n_D,
            "n_hits": self.n_hits, "n_entries": self.n_entries,
        }

    def summary(self):
        """返回多行 summary 字符串 (子类可覆写扩展饱和项)"""
        lines = [
            f"α₀ = {self.alpha0:.6g} Gy⁻¹   β₀ = {self.beta0:.6g} Gy⁻²",
            f"z̄_{{d,F}} = {self._fmt(self.z_d_F)} Gy   z̄_{{d,D}} = {self._fmt(self.z_d_D)} Gy",
            f"z̄_{{n,F}} = {self._fmt(self.z_n_F)} Gy   z̄_{{n,D}} = {self._fmt(self.z_n_D)} Gy",
            f"α_S = {self.alpha_S:.6g} Gy⁻¹",
            f"β_S = {self.beta_S:.6g} Gy⁻²",
        ]
        return "\n".join(lines)

    @staticmethod
    def _fmt(x):
        return f"{x:.6g}" if x is not None else "N/A"

    def __repr__(self):
        return (f"MK(α₀={self.alpha0:.4f}, β₀={self.beta0:.4f} → "
                f"α_S={self.alpha_S:.4f}, β_S={self.beta_S:.5f}; "
                f"z̄_{{d,D}}={self._fmt(self.z_d_D)}, "
                f"z̄_{{n,D}}={self._fmt(self.z_n_D)})")


# ============================================================
# 命令行自检
# ============================================================
if __name__ == "__main__":
    print("== MK 基类自检 (HSG 量级示例微剂量学量) ==")
    mk = MK(alpha0=0.174, beta0=0.0568, z_d_D=54.1, z_n_D=0.107)
    print(mk)
    print("-" * 50)
    print(mk.summary())
    print("-" * 50)
    D = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
    S = mk.survival(D)
    print("D (Gy):", D.tolist())
    print("S(D)  :", np.round(S, 5).tolist())
