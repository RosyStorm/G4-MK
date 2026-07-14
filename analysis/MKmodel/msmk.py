#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
msmk.py — 修正 SMK (ModifiedSMK, Inaniwa & Kanematsu 2018 闭式) 模型类

继承自 MK 基类。
参考: Inaniwa & Kanematsu, Phys. Med. Biol. 63 (2018) 095011。

关键公式:
    α_S = α₀ + β₀ · z̄*_{d,D}                       (式 15, 用近似 z̄*_{d,F}/z̄_{d,F} ≈ 1)
    β_S = β₀ · z̄*_{d,D} / z̄_{d,D}                 (式 16)
    S(D) = exp(-α_S·D - β_S·D²) · {1 + (D·z̄_{n,D}/2)·((α_S+2β_S·D)² - 2β_S)}   (式 24, 闭式)

用法:
    smk = ModifiedSMK.from_root("data/microtrack.root", alpha0=0.174, beta0=0.0568, z0=66.0)
    print(smk)
    S = smk.survival(np.linspace(0, 10, 200))
"""

import numpy as np

from .mk import MK


class ModifiedSMK(MK):
    """
    修正 SMK 模型 (Inaniwa 2018 闭式 式 24)。
    """

    # ----- 构造 -----
    def __init__(self, alpha0, beta0, z0, z_d_D, zs_d_D, z_n_D,
                 z_d_F=None, z_n_F=None):
        """
        :param alpha0:  参考 α₀ (Gy⁻¹)
        :param beta0:   参考 β₀ (Gy⁻²)
        :param z0:      饱和参数 z₀ (Gy)
        :param z_d_D:   域剂量均 z̄_{d,D} (Gy), 必须 > 0
        :param zs_d_D:  域饱和剂量均 z̄*_{d,D} (Gy)
        :param z_n_D:   核剂量均 z̄_{n,D} (Gy)
        :param z_d_F:   域频率均 z̄_{d,F} (Gy, 可选)
        :param z_n_F:   核频率均 z̄_{n,F} (Gy, 可选)
        :raises ValueError: z_d_D ≤ 0
        """
        if z_d_D <= 0:
            raise ValueError(f"z_d_D 必须 > 0 (β_S 分母), 当前 {z_d_D}")

        # MSMK 专属字段
        self.z0     = float(z0)
        self.zs_d_D = float(zs_d_D)

        super().__init__(alpha0, beta0, z_d_D, z_n_D,
                         z_d_F=z_d_F, z_n_F=z_n_F)

    # ----- 子类覆写: 系数公式 -----
    def _compute_alpha_beta(self):
        """
        Inaniwa 2018 式 (15)(16):
            α_S = α₀ + β₀ · z̄*_{d,D}
            β_S = β₀ · z̄*_{d,D} / z̄_{d,D}
        """
        self.alpha_S = self.alpha0 + self.beta0 * self.zs_d_D
        self.beta_S  = self.beta0 * self.zs_d_D / self.z_d_D

    # ----- 子类覆写: 存活曲线 (闭式 LQ + z_n bracket 修正) -----
    def survival(self, D):
        """
        Inaniwa 2018 式 (24) 闭式:
            S(D) = exp(-α_S·D - β_S·D²) · [1 + (D·z̄_{n,D}/2)·((α_S+2β_S·D)² - 2β_S)]
        bracket 裁剪 ≥0, 防数值噪声导致负存活。
        """
        D = np.asarray(D, dtype=float)
        lq = np.exp(-self.alpha_S * D - self.beta_S * D ** 2)
        bracket = (1.0
                   + (D * self.z_n_D / 2.0)
                     * ((self.alpha_S + 2.0 * self.beta_S * D) ** 2
                        - 2.0 * self.beta_S))
        bracket = np.clip(bracket, 0.0, None)
        return lq * bracket

    # ----- from_root 工厂 -----
    @classmethod
    def from_root(cls, root_path, tree_name, alpha0, beta0, z0):
        from .mk import (load_microdosimetry_from_root,
                          compute_basic_microdosimetry)

        z_d, z_n, w, n_hits, n_entries = load_microdosimetry_from_root(
            root_path, tree_name)
        z_d_F, z_d_D, z_n_F, z_n_D = compute_basic_microdosimetry(
            z_d, z_n, w, n_hits)
        zs_d_D = (w * cls.saturation(z_d, z0) * cls.saturation(z_d, z0)).sum() \
                 / (w * z_d).sum()

        obj = cls(alpha0, beta0, z0, z_d_D, zs_d_D, z_n_D,
                  z_d_F=z_d_F, z_n_F=z_n_F)
        obj.n_hits    = n_hits
        obj.n_entries = n_entries
        return obj

    # ----- 扩展 summary / __repr__ (含饱和项) -----
    def summary(self):
        lines = [
            f"α₀ = {self.alpha0:.6g} Gy⁻¹   β₀ = {self.beta0:.6g} Gy⁻²   z₀ = {self.z0:.4g} Gy",
            f"z̄_{{d,F}} = {self._fmt(self.z_d_F)} Gy   z̄_{{d,D}} = {self._fmt(self.z_d_D)} Gy",
            f"z̄*_{{d,D}} = {self._fmt(self.zs_d_D)} Gy   (饱和)",
            f"z̄_{{n,F}} = {self._fmt(self.z_n_F)} Gy   z̄_{{n,D}} = {self._fmt(self.z_n_D)} Gy",
            f"α_S = α₀ + β₀·z̄*_{{d,D}} = {self.alpha_S:.6g} Gy⁻¹",
            f"β_S = β₀·z̄*_{{d,D}}/z̄_{{d,D}} = {self.beta_S:.6g} Gy⁻²",
        ]
        return "\n".join(lines)

    def microdosimetry(self):
        d = super().microdosimetry()
        d.update({"z0": self.z0, "zs_d_D": self.zs_d_D})
        return d

    def __repr__(self):
        return (f"ModifiedSMK(α₀={self.alpha0:.4f}, β₀={self.beta0:.4f}, "
                f"z₀={self.z0:.2f} → "
                f"α_S={self.alpha_S:.4f}, β_S={self.beta_S:.5f}; "
                f"z̄*_{{d,D}}={self._fmt(self.zs_d_D)}, "
                f"z̄_{{n,D}}={self._fmt(self.z_n_D)})")


# ============================================================
# 命令行自检
# ============================================================
if __name__ == "__main__":
    print("== ModifiedSMK 自检 (HSG 量级示例微剂量学量) ==")
    smk = ModifiedSMK(alpha0=0.174, beta0=0.0568, z0=66.0,
                      z_d_D=54.1, zs_d_D=32.1, z_n_D=0.107)
    print(smk)
    print("-" * 50)
    print(smk.summary())
    print("-" * 50)
    D = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
    S = smk.survival(D)
    print("D (Gy):", D.tolist())
    print("S(D)  :", np.round(S, 5).tolist())
