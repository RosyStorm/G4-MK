#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
smk.py — SMK (Sato 2012 单随机) 模型类

继承自 MK 基类。
参考: Sato & Furusawa, Radiat. Res. 178 (2012) 341–356 式 (33)。
关键公式:
    α_S = α₀ · z̄*_{d,F}/z̄_{d,F} + β₀ · z̄*_{d,D}   (式 33 第一项)
    β_S = β₀ · z̄*_{d,D} / z̄_{d,D}                 (式 33 第二项)
    S(z_n) = exp(-α_S·z_n - β_S·z_n²)
    S(D) = ⟨ S(z_n) ⟩_{z_n ~ Poisson 卷积}        (核 MC 单随机)

用法:
    smk = SMK.from_root("data/microtrack.root", alpha0=0.174, beta0=0.0568, z0=66.0)
    print(smk)
    S = smk.survival(np.linspace(0, 10, 200))
"""

import numpy as np

from .mk import MK


class SMK(MK):
    """
    SMK 模型 (Sato 2012 式 33): 核 z_n 单随机, 域 z_d 已并入 α_S/β_S。
    """

    # ----- 构造 -----
    def __init__(self, alpha0, beta0, z0, z_d_D, zs_d_D, zs_d_F, z_n_D,
                 z_d_F=None, z_n_F=None, z_n_samples=None):
        """
        :param alpha0:   参考线性项 α₀ (Gy⁻¹)
        :param beta0:    参考二次项 β₀ (Gy⁻²)
        :param z0:       饱和参数 z₀ (Gy)
        :param z_d_D:    域剂量均 z̄_{d,D} (Gy), 必须 > 0
        :param zs_d_D:   域饱和剂量均 z̄*_{d,D} (Gy)
        :param zs_d_F:   域饱和频率均 z̄*_{d,F} (Gy)
        :param z_n_D:    核剂量均 z̄_{n,D} (Gy)
        :param z_d_F:    域频率均 z̄_{d,F} (Gy, 可选)
        :param z_n_F:    核频率均 z̄_{n,F} (Gy, 可选)
        :param z_n_samples: 核比能样本 f_{n,1} (MC survival 用)
        :raises ValueError: z_d_F 或 z_d_D ≤ 0
        """
        if z_d_D <= 0:
            raise ValueError(f"z_d_D 必须 > 0 (β_S 分母), 当前 {z_d_D}")

        # SMK 专属字段
        self.z0          = float(z0)
        self.zs_d_D      = float(zs_d_D)
        self.zs_d_F      = float(zs_d_F)
        self.z_n_samples = (None if z_n_samples is None
                            else np.asarray(z_n_samples, dtype=float))

        # 父类: alpha0/beta0/z_d_F/z_d_D/z_n_F/z_n_D 存储 + _compute_alpha_beta
        super().__init__(alpha0, beta0, z_d_D, z_n_D,
                         z_d_F=z_d_F, z_n_F=z_n_F)

    # ----- 子类覆写: 系数公式 -----
    def _compute_alpha_beta(self):
        """
        Sato 2012 式 (33) — α_S 保留 z̄*_{d,F}/z̄_{d,F} 因子 (与 msmk 的近似 (13) 不同).
            α_S = α₀ · z̄*_{d,F}/z̄_{d,F} + β₀ · z̄*_{d,D}
            β_S = β₀ · z̄*_{d,D} / z̄_{d,D}
        """
        if self.z_d_F is None or self.z_d_F <= 0:
            raise ValueError("SMK 需要 z_d_F > 0 (α_S 第一项分母)")
        self.alpha_S = (self.alpha0 * self.zs_d_F / self.z_d_F
                        + self.beta0 * self.zs_d_D)
        self.beta_S = self.beta0 * self.zs_d_D / self.z_d_D

    # ----- 子类覆写: 存活曲线 (核 z_n MC) -----
    def survival(self, D, N=5000):
        """
        SMK Monte Carlo 单随机 S(D):
            核命中数 k ~ Poisson(D / z̄_{n,F})
            核比能  z_n = k 个 f_{n,1} 样本求和
            单次存活 S(z_n) = exp(-α_S·z_n - β_S·z_n²)
            S(D) = N 次试验均值
        """
        if self.z_n_samples is None:
            raise ValueError("MC survival 需 z_n_samples; 请用 from_root(...) 或传 z_n_samples")
        if self.z_n_F is None or self.z_n_F <= 0:
            raise ValueError("z_n_F 缺失或 ≤0")

        D = np.asarray(D, dtype=float)
        scal = (D.ndim == 0)
        D_flat = np.atleast_1d(D).astype(float)
        samples = self.z_n_samples
        M = samples.size
        out = np.empty(D_flat.size, dtype=float)

        for i in range(D_flat.size):
            Di = D_flat[i]
            if Di <= 0.0:
                out[i] = 1.0
                continue

            mean_k = Di / self.z_n_F
            k_arr = np.random.poisson(mean_k, size=N)
            total = int(k_arr.sum())
            if total == 0:
                out[i] = 1.0
                continue

            idx = np.random.randint(0, M, size=total)
            draws = samples[idx]
            seg = np.repeat(np.arange(N), k_arr)
            z_n = np.bincount(seg, weights=draws, minlength=N)

            S_trial = np.exp(-self.alpha_S * z_n - self.beta_S * z_n * z_n)
            out[i] = float(S_trial.mean())

        return float(out[0]) if scal else out.reshape(D.shape)

    # ----- from_root 工厂 -----
    @classmethod
    def from_root(cls, root_path, tree_name, alpha0, beta0, z0):
        """
        SMK 专属 from_root: 基类基础上额外算饱和微剂量学量 z̄*_{d,D} / z̄*_{d,F}。
        """
        from .mk import (load_microdosimetry_from_root,
                          compute_basic_microdosimetry)

        z_d, z_n, w, n_hits, n_entries = load_microdosimetry_from_root(
            root_path, tree_name)
        z_d_F, z_d_D, z_n_F, z_n_D = compute_basic_microdosimetry(
            z_d, z_n, w, n_hits)

        # 饱和: 与基类相同的 saturation (z₀·√(1−exp(−x²)))
        zs = cls.saturation(z_d, z0)
        zs_d_D = (w * zs * zs).sum() / (w * z_d).sum()
        zs_d_F = (w * zs).sum() / w.sum()

        obj = cls(alpha0, beta0, z0, z_d_D, zs_d_D, zs_d_F, z_n_D,
                  z_d_F=z_d_F, z_n_F=z_n_F, z_n_samples=z_n)
        obj.n_hits    = n_hits
        obj.n_entries = n_entries
        return obj

    # ----- 扩展 summary / __repr__ (含饱和项) -----
    def summary(self):
        lines = [
            f"α₀ = {self.alpha0:.6g} Gy⁻¹   β₀ = {self.beta0:.6g} Gy⁻²   z₀ = {self.z0:.4g} Gy",
            f"z̄_{{d,F}} = {self._fmt(self.z_d_F)} Gy   z̄_{{d,D}} = {self._fmt(self.z_d_D)} Gy",
            f"z̄*_{{d,D}} = {self._fmt(self.zs_d_D)} Gy   (饱和)",
            f"z̄*_{{d,F}} = {self._fmt(self.zs_d_F)} Gy   (饱和)",
            f"z̄_{{n,F}} = {self._fmt(self.z_n_F)} Gy   z̄_{{n,D}} = {self._fmt(self.z_n_D)} Gy",
            f"α_S = α₀·z̄*_{{d,F}}/z̄_{{d,F}} + β₀·z̄*_{{d,D}} = {self.alpha_S:.6g} Gy⁻¹",
            f"β_S = β₀·z̄*_{{d,D}}/z̄_{{d,D}} = {self.beta_S:.6g} Gy⁻²",
        ]
        return "\n".join(lines)

    def microdosimetry(self):
        d = super().microdosimetry()
        d.update({"z0": self.z0, "zs_d_F": self.zs_d_F, "zs_d_D": self.zs_d_D})
        return d

    def __repr__(self):
        return (f"SMK(α₀={self.alpha0:.4f}, β₀={self.beta0:.4f}, z₀={self.z0:.2f} → "
                f"α_S={self.alpha_S:.4f}, β_S={self.beta_S:.5f}; "
                f"z̄*_{{d,D}}={self._fmt(self.zs_d_D)}, "
                f"z̄_{{n,D}}={self._fmt(self.z_n_D)})")


# ============================================================
# 命令行自检
# ============================================================
if __name__ == "__main__":
    print("== SMK 自检 (HSG 量级示例微剂量学量) ==")
    smk = SMK(alpha0=0.174, beta0=0.0568, z0=66.0,
              z_d_D=54.1, zs_d_D=32.1, zs_d_F=15.0, z_n_D=0.107)
    print(smk)
    print("-" * 50)
    print(smk.summary())
    print("-" * 50)
    D = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
    S = smk.survival(D, N=2000)
    print("D (Gy):", D.tolist())
    print("S(D)  :", np.round(S, 5).tolist())
