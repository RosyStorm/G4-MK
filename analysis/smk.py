#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
smk.py — SMK (Stochastic Microdosimetric Kinetic) 模型类
参考: sato0(2012)

模型(SMK, 单随机)
微剂量学量 (单事件, 仅命中事件 hitFlag==1, 域级带重要性权 w):
    z̄_{d,F} = Σ(w·z_d)/Σw                   (式10, 域频率均)
    z̄_{d,D} = Σ(w·z_d²)/Σ(w·z_d)            (式11, 域剂量均)
    z̄*_{d,D}= Σ(w·z*²)/Σ(w·z_d)             (式12, 域饱和剂量均)
    z̄_{n,D} = Σz_n²/Σz_n                    (核剂量均, 无权)

用法:
    # 1) 直接从 microtrack.root 构建模型(自动算微剂量学量)
    from analysis.smk import SMK
    smk = SMK.from_root("data/microtrack.root", alpha0=0.174, beta0=0.0568, z0=66.0)
    print(smk)                       # → α_S, β_S 等
    S = smk.survival(np.linspace(0, 10, 200))

    # 2) 或已有微剂量学量, 直接构造
    smk = SMK(alpha0=0.174, beta0=0.0568, z0=66.0,
              z_d_D=54.1, zs_d_D=32.1, z_n_D=0.107)
    smk.alpha_S, smk.beta_S          # SMK 系数
"""

import numpy as np


class SMK:
    """
    修正 SMK 模型 (Inaniwa 2018 闭式 式24)。

    构造后即可:
        - 读 .alpha_S / .beta_S / .z_n_D 等属性
        - 调 .survival(D) 得存活曲线
        - 调 .microdosimetry() / .summary() 取完整结果
    """

    # ============================================================
    # 构造
    # ============================================================
    def __init__(self, alpha0, beta0, z0, z_d_D, zs_d_D, zs_d_F, z_n_D,
                 z_d_F=None, z_n_F=None, z_n_samples=None):
        """
        :param alpha0:  参考线性项 α₀ (Gy⁻¹, 一般 X 射线 LQ 拟合)
        :param beta0:   参考二次项 β₀ (Gy⁻²)
        :param z0:      饱和参数 z₀ (Gy)
        :param z_d_D:   域剂量均 z̄_{d,D} (Gy)
        :param zs_d_D:  域饱和剂量均 z̄*_{d,D} (Gy)
        :param zs_d_F:  域饱和频率均 z̄*_{d,F} ()
        :param z_n_D:   核剂量均 z̄_{n,D} (Gy)
        :param z_d_F:   域频率均 z̄_{d,F} (Gy, 可选, 仅记录/诊断)
        :param z_n_F:   核频率均 z̄_{n,F} (Gy, 可选, 仅记录/诊断)

        注: z_d_D 必须 > 0 (β_S 分母); 否则 ValueError。
        """
        if z_d_D <= 0:
            raise ValueError(f"z_d_D 必须 > 0 (β_S 分母), 当前 {z_d_D}")

        # —— 细胞系 / 饱和参数 ——
        self.alpha0 = float(alpha0)        # α₀ (Gy⁻¹)
        self.beta0 = float(beta0)          # β₀ (Gy⁻²)
        self.z0 = float(z0)                # z₀ (Gy)

        # —— 单事件微剂量学量 ——
        self.z_d_F = None if z_d_F is None else float(z_d_F)   # z̄_{d,F}
        self.z_d_D = float(z_d_D)                               # z̄_{d,D}
        self.zs_d_D = float(zs_d_D)                             # z̄*_{d,D}
        self.zs_d_F = float(zs_d_F)                             # z̄*_{d,F}
        self.z_n_F = None if z_n_F is None else float(z_n_F)    # z̄_{n,F}
        self.z_n_D = float(z_n_D)                               # z̄_{n,D}

        # —— SMK 系数 (式15-16) ——
        self.alpha_S = self.alpha0 * self.zs_d_F / self.z_d_F + self.beta0 * self.zs_d_D        # α_S (Gy⁻¹)
        self.beta_S = self.beta0 * self.zs_d_D / self.z_d_D         # β_S (Gy⁻²)

        # —— f_{n,1} 单事件样本(MC survival 卷积抽样用, from_root 时填) ——
        self.z_n_samples = None if z_n_samples is None else np.asarray(z_n_samples, dtype=float)

        # —— 统计来源(可选, from_root 时填) ——
        self.n_hits = None       # 命中事件数
        self.n_entries = None    # ntuple 总条目数

    # ============================================================
    # 静态: 饱和修正 (式1)
    # ============================================================
    @staticmethod
    def saturation(z_d, z0):
        """
        饱和比能 z*_d = z₀·(1 - e^{-(z_d/z₀)²})  (式1, 无 sqrt, 标量或数组)
        """
        z_d = np.asarray(z_d, dtype=float)
        return z0 * np.sqrt((1.0 - np.exp(-(z_d / z0) ** 2)))

    # ============================================================
    # 备选构造: 从 microtrack.root events ntuple
    # ============================================================
    @classmethod
    def from_root(cls, root_path, alpha0, beta0, z0):
        """
        从 microtrack.root 的 events ntuple 算微剂量学量并构建模型。
        仅用命中事件 (hitFlag==1), 与 analyze_mk.C / plot_smk_ac225_overlay.py 一致。

        :param root_path: microtrack.root 路径 (含 events ntuple)
        :param alpha0, beta0, z0: 细胞系参数 + 饱和参数
        :return: ModifiedSMK 实例 (附带 .n_hits/.n_entries)
        """
        import ROOT  # 延迟导入: 避免 import smk 时强依赖 PyROOT

        rdf = ROOT.RDataFrame("events", str(root_path))
        n_entries = rdf.Count().GetValue()
        if n_entries == 0:
            raise ValueError(f"{root_path} 内 events ntuple 为空")

        hits = rdf.Filter("hitFlag==1")
        n_hits = hits.Count().GetValue()
        if n_hits == 0:
            raise ValueError(f"{root_path} 无命中事件 (hitFlag==1), 无法算微剂量学量")

        cols = hits.AsNumpy(["z_d_Gy", "z_n_Gy", "weight"])
        z_d = np.asarray(cols["z_d_Gy"], dtype=float)
        z_n = np.asarray(cols["z_n_Gy"], dtype=float)
        w = np.asarray(cols["weight"], dtype=float)

        # —— 域级(带权) ——
        Sw = w.sum()                       # Σw
        Swz = (w * z_d).sum()              # Σ(w·z_d)
        Swzz = (w * z_d * z_d).sum()       # Σ(w·z_d²)
        if Sw <= 0 or Swz <= 0:
            raise ValueError("域级求和异常 (Σw 或 Σ(w·z_d) ≤ 0)")
        z_d_F = Swz / Sw                   # 式10
        z_d_D = Swzz / Swz                 # 式11
        zs = cls.saturation(z_d, z0)       # 式1
        zs_d_D = (w * zs * zs).sum() / Swz  # 式12
        zs_d_F = (w * zs).sum() / Sw  

        # —— 核级(无权) ——
        Szn = z_n.sum()                    # Σz_n
        if Szn <= 0:
            raise ValueError("核级求和异常 (Σz_n ≤ 0)")
        z_n_F = Szn / n_hits               # z̄_{n,F}
        z_n_D = (z_n * z_n).sum() / Szn    # z̄_{n,D}

        obj = cls(alpha0, beta0, z0, z_d_D, zs_d_D, zs_d_F,  z_n_D,
                  z_d_F=z_d_F, z_n_F=z_n_F, z_n_samples=z_n)
        obj.n_hits = int(n_hits)
        obj.n_entries = int(n_entries)
        return obj

    # ============================================================
    # 存活曲线 (式24)
    # ============================================================
    def survival(self, D, N=5000):
        """
        SMK 存活分数 S(D) —— Monte Carlo 单随机积分(核比能 z_n 随机; 域已并入 α_S/β_S)。

        逻辑(Sato SMK, 核 z_n 单随机):
            对剂量 D, 每个细胞核命中次数 k ~ Poisson(D/z̄_{n,F});
            k 次命中下核比能 z_n ~ f_{n,k} = f_{n,1}^{*k}  (用 k 个 f_{n,1} 样本求和实现卷积);
            单次抽样存活 S(z_n) = exp(-α_S·z_n - β_S·z_n²);
            S(D) = N 次抽样下 S(z_n) 的均值 = ∫ S(z_n) f_n(z_n|D) dz_n。

        numpy 向量化: 仅对 D 数组外层循环; 内层 N 次试验用 bincount 分组求和,
        不逐试验/逐样本 for。

        :param D: 剂量 (Gy), 标量或数组
        :param N: MC 试验次数(每次代表一个细胞核); 越大越准但越慢
        :return: 存活分数(无量纲, 与 D 同形)
        """
        # —— 前置检查: MC 需要 f_{n,1} 单事件样本与 z̄_{n,F} ——
        if self.z_n_samples is None:
            raise ValueError("MC survival 需 f_{n,1} 单事件样本; 请用 from_root(...) 构造, "
                             "或构造时传 z_n_samples")
        if self.z_n_F is None or self.z_n_F <= 0:
            raise ValueError("z_n_F 缺失或 ≤0, 无法算 Poisson 均值")

        D = np.asarray(D, dtype=float)
        scal = (D.ndim == 0)                       # 记录是否标量输入
        D_flat = np.atleast_1d(D).astype(float)
        samples = self.z_n_samples                 # f_{n,1} 单事件样本
        M = samples.size
        out = np.empty(D_flat.size, dtype=float)

        for i in range(D_flat.size):
            Di = D_flat[i]
            if Di <= 0.0:
                out[i] = 1.0                       # D=0 → k=0 → z_n=0 → S=1
                continue

            mean_k = Di / self.z_n_F               # Poisson 均值(每核命中次数期望)
            k_arr = np.random.poisson(mean_k, size=N)      # (N,) 各试验命中次数
            total = int(k_arr.sum())
            if total == 0:
                out[i] = 1.0
                continue

            # —— 一次性抽 total 个 f_{n,1} 样本(向量化, 不逐次) ——
            idx = np.random.randint(0, M, size=total)
            draws = samples[idx]

            # —— 按 k_arr 分组求和 → 每次试验的 z_n; bincount 法, k=0 的试验自动得 0 ——
            seg = np.repeat(np.arange(N), k_arr)
            z_n = np.bincount(seg, weights=draws, minlength=N)   # (N,) 每试验的核比能

            # —— 单次抽样存活 + 求均值 ——
            S_trial = np.exp(- self.alpha_S * z_n - self.beta_S * z_n * z_n)
            out[i] = float(S_trial.mean())

        return float(out[0]) if scal else out.reshape(D.shape)

    # ============================================================
    # 取值辅助
    # ============================================================
    def coefficients(self):
        """返回 SMK 系数 (α_S, β_S)"""
        return self.alpha_S, self.beta_S

    def microdosimetry(self):
        """返回单事件微剂量学量 dict (z_d_F/z_d_D/zs_d_D/z_n_F/z_n_D)"""
        return {
            "z_d_F": self.z_d_F, "z_d_D": self.z_d_D, "zs_d_D": self.zs_d_D,
            "z_n_F": self.z_n_F, "z_n_D": self.z_n_D,
            "n_hits": self.n_hits, "n_entries": self.n_entries,
        }

    def summary(self):
        """
        返回多行字符串: 细胞系参数 + 微剂量学量 + SMK 系数。
        便于脚本打印或写摘要文件。
        """
        lines = [
            f"α₀ = {self.alpha0:.6g} Gy⁻¹   β₀ = {self.beta0:.6g} Gy⁻²   z₀ = {self.z0:.4g} Gy",
            f"z̄_{{d,F}} = {self._fmt(self.z_d_F)} Gy   z̄_{{d,D}} = {self.z_d_D:.6g} Gy",
            f"z̄*_{{d,D}} = {self.zs_d_D:.6g} Gy   (饱和)",
            f"z̄_{{n,F}} = {self._fmt(self.z_n_F)} Gy   z̄_{{n,D}} = {self.z_n_D:.6g} Gy",
            f"α_S = α₀ + β₀·z̄*_{{d,D}}     = {self.alpha_S:.6g} Gy⁻¹",
            f"β_S = β₀·z̄*_{{d,D}}/z̄_{{d,D}} = {self.beta_S:.6g} Gy⁻²",
        ]
        return "\n".join(lines)

    @staticmethod
    def _fmt(x):
        return f"{x:.6g}" if x is not None else "N/A"

    def __repr__(self):
        return (f"SMK(α₀={self.alpha0:.4f}, β₀={self.beta0:.4f}, z₀={self.z0:.2f} → "
                f"α_S={self.alpha_S:.4f}, β_S={self.beta_S:.5f}; "
                f"z̄*_{{d,D}}={self.zs_d_D:.3f}, z̄_{{n,D}}={self.z_n_D:.4g})")


# ============================================================
# 命令行自检 (不依赖 root 文件, 用示例微剂量学量验证模型)
# ============================================================
if __name__ == "__main__":
    print("== SMK 自检 (示例微剂量学量) ==")
    # 用 HSG 量级示例值
    smk = SMK(alpha0=0.174, beta0=0.0568, z0=66.0,
                      z_d_D=54.1, zs_d_D=32.1, z_n_D=0.107)
    print(smk)
    print("-" * 50)
    print(smk.summary())
    print("-" * 50)
    D = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
    S = smk.survival(D)
    print("D (Gy):", D.tolist())
    print("S(D)  :", np.round(S, 5).tolist())
