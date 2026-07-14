#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dsmk.py — DSMK (Sato 2012 双随机) 模型类

继承自 MK 基类。
参考: Sato & Furusawa, Radiat. Res. 178 (2012) 341–356。

关键公式:
    z'_d = z₀ · √(1 − exp(−(z_d/z₀)²))              (式 16, 饱和)
    α_S = α₀, β_S = β₀                             (饱和放在 z'_d 里, Inaniwa 风格)
    S_n(z_n) = exp(-α₀·⟨z'_d⟩ - β₀·⟨z'_d²⟩)         (式 19, 单 z_n 下域卷积)
    S(D)    = ⟨ S_n(z_n) ⟩_{z_n ~ Poisson(D/z̄_{n,F})}   (双 MC)

用法:
    dsmk = DSMK.from_root("data/microtrack.root", alpha0=0.174, beta0=0.0568, z0=66.0)
    print(dsmk)
    S = dsmk.survival(np.linspace(0, 10, 200))
"""

import numpy as np

from .mk import MK


class DSMK(MK):
    """
    DSMK 模型 (Sato 2012 双随机 MC, 域 + 核都做 Poisson 卷积)。
    """

    # ----- 构造 -----
    def __init__(self, alpha0, beta0, z0, z_d_D, zs_d_D, zs_d_F, z_n_D,
                 z_d_F=None, z_n_F=None, z_n_samples=None, z_d_samples=None):
        """
        :param alpha0:    参考 α₀ (Gy⁻¹)
        :param beta0:     参考 β₀ (Gy⁻²)
        :param z0:        饱和参数 z₀ (Gy)
        :param z_d_D:     域剂量均 z̄_{d,D} (Gy), 必须 > 0
        :param zs_d_D:    域饱和剂量均 z̄*_{d,D} (Gy)
        :param zs_d_F:    域饱和频率均 z̄*_{d,F} (Gy)
        :param z_n_D:     核剂量均 z̄_{n,D} (Gy)
        :param z_d_F:     域频率均 z̄_{d,F} (Gy, 可选)
        :param z_n_F:     核频率均 z̄_{n,F} (Gy, 可选)
        :param z_n_samples: 核比能 f_{n,1} 样本 (MC 卷积)
        :param z_d_samples: 域比能 f_{d,1} 样本 (MC 卷积)
        :raises ValueError: z_d_D ≤ 0
        """
        if z_d_D <= 0:
            raise ValueError(f"z_d_D 必须 > 0 (β_S 分母), 当前 {z_d_D}")

        self.z0          = float(z0)
        self.zs_d_D      = float(zs_d_D)
        self.zs_d_F      = float(zs_d_F)
        self.z_n_samples = (None if z_n_samples is None
                            else np.asarray(z_n_samples, dtype=float))
        self.z_d_samples = (None if z_d_samples is None
                            else np.asarray(z_d_samples, dtype=float))

        super().__init__(alpha0, beta0, z_d_D, z_n_D,
                         z_d_F=z_d_F, z_n_F=z_n_F)

    # ----- 子类覆写: 系数公式 (DSMK 特殊: 用 α₀/β₀, 饱和在 z'_d 里) -----
    def _compute_alpha_beta(self):
        """
        DSMK 不使用 α_S/β_S, 直接 α₀/β₀, 饱和放在 z'_d 里 (Sato 式 16)。
        alpha_S/beta_S 字段保留以兼容父类接口。
        """
        self.alpha_S = self.alpha0
        self.beta_S  = self.beta0

    # ----- 子类覆写: 存活曲线 (核 + 域 双 MC) -----
    def survival(self, D, N=2000):
        """
        DSMK Monte Carlo 双随机 S(D):
            外层 (核):
                k_n ~ Poisson(D / z̄_{n,F});
                z_n = k_n 个 f_{n,1} 样本求和;
                S_n(z_n) 由 self.survival_n 算;
                S(D) = ⟨S_n⟩(核 MC).
            内层 survival_n(z_n): 域 MC:
                k_d ~ Poisson(z_n / z̄_{d,F});
                z_d = k_d 个 f_{d,1} 样本求和;
                z'_d = saturation(z_d);
                S_n = exp(-α₀·⟨z'_d⟩ - β₀·⟨z'_d²⟩).
        """
        if self.z_n_samples is None:
            raise ValueError("MC survival 需 z_n_samples; 请用 from_root(...) 或传 z_n_samples")
        if self.z_n_F is None or self.z_n_F <= 0:
            raise ValueError("z_n_F 缺失或 ≤0")

        D = np.asarray(D, dtype=float)
        scal = (D.ndim == 0)
        D_flat = np.atleast_1d(D).astype(float)
        samples_n = self.z_n_samples
        M_n = samples_n.size
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

            idx = np.random.randint(0, M_n, size=total)
            draws = samples_n[idx]
            seg = np.repeat(np.arange(N), k_arr)
            z_n = np.bincount(seg, weights=draws, minlength=N)

            S_trial = self.survival_n(z_n)
            out[i] = float(S_trial.mean())

        return float(out[0]) if scal else out.reshape(D.shape)

    def survival_n(self, z_n, N=2000):
        """
        DSMK 域级存活 S_n(z_n) (Sato 式 19) — Monte Carlo 域卷积:
            给定 z_n, k_d ~ Poisson(z_n / z̄_{d,F});
            z_d = k_d 个 f_{d,1} 样本求和;
            z'_d = saturation(z_d);   (式 16)
            S_n = exp(-α₀·⟨z'_d⟩ - β₀·⟨z'_d²⟩).   (式 19)
        """
        if self.z_d_samples is None:
            raise ValueError("survival_n 需 z_d_samples; 请用 from_root(...) 或传 z_d_samples")
        if self.z_d_F is None or self.z_d_F <= 0:
            raise ValueError("z_d_F 缺失或 ≤0")

        z_n = np.asarray(z_n, dtype=float)
        scal = (z_n.ndim == 0)
        z_n_flat = np.atleast_1d(z_n).astype(float)
        samples_d = self.z_d_samples
        M_d = samples_d.size
        out = np.empty(z_n_flat.size, dtype=float)

        for i in range(z_n_flat.size):
            zni = z_n_flat[i]
            if zni <= 0.0:
                out[i] = 1.0
                continue

            mean_k = zni / self.z_d_F
            k_arr = np.random.poisson(mean_k, size=N)
            total = int(k_arr.sum())
            if total == 0:
                out[i] = 1.0
                continue

            idx = np.random.randint(0, M_d, size=total)
            draws = samples_d[idx]
            seg = np.repeat(np.arange(N), k_arr)
            zd_multi = np.bincount(seg, weights=draws, minlength=N)

            zd_s = self.saturation(zd_multi, self.z0)
            zd2_s = zd_s * zd_s
            zds_mean  = float(zd_s.mean())
            zds2_mean = float(zd2_s.mean())

            out[i] = float(np.exp(-self.alpha0 * zds_mean
                                   - self.beta0 * zds2_mean))

        return float(out[0]) if scal else out.reshape(z_n.shape)

    # ----- from_root 工厂 -----
    @classmethod
    def from_root(cls, root_path, tree_name, alpha0, beta0, z0):
        from .mk import (load_microdosimetry_from_root,
                          compute_basic_microdosimetry)

        z_d, z_n, w, n_hits, n_entries = load_microdosimetry_from_root(
            root_path, tree_name)
        z_d_F, z_d_D, z_n_F, z_n_D = compute_basic_microdosimetry(
            z_d, z_n, w, n_hits)

        zs = cls.saturation(z_d, z0)
        zs_d_D = (w * zs * zs).sum() / (w * z_d).sum()
        zs_d_F = (w * zs).sum() / w.sum()

        obj = cls(alpha0, beta0, z0, z_d_D, zs_d_D, zs_d_F, z_n_D,
                  z_d_F=z_d_F, z_n_F=z_n_F,
                  z_n_samples=z_n, z_d_samples=z_d)
        obj.n_hits    = n_hits
        obj.n_entries = n_entries
        return obj

    # ----- 扩展 summary / __repr__ -----
    def summary(self):
        lines = [
            f"α₀ = {self.alpha0:.6g} Gy⁻¹   β₀ = {self.beta0:.6g} Gy⁻²   z₀ = {self.z0:.4g} Gy",
            f"z̄_{{d,F}} = {self._fmt(self.z_d_F)} Gy   z̄_{{d,D}} = {self._fmt(self.z_d_D)} Gy",
            f"z̄*_{{d,D}} = {self._fmt(self.zs_d_D)} Gy   (饱和)",
            f"z̄*_{{d,F}} = {self._fmt(self.zs_d_F)} Gy   (饱和)",
            f"z̄_{{n,F}} = {self._fmt(self.z_n_F)} Gy   z̄_{{n,D}} = {self._fmt(self.z_n_D)} Gy",
            f"α_S = α₀  (DSMK 用 α₀/β₀, 饱和在 z'_d 里)",
            f"β_S = β₀",
        ]
        return "\n".join(lines)

    def microdosimetry(self):
        d = super().microdosimetry()
        d.update({"z0": self.z0, "zs_d_F": self.zs_d_F, "zs_d_D": self.zs_d_D})
        return d

    def __repr__(self):
        return (f"DSMK(α₀={self.alpha0:.4f}, β₀={self.beta0:.4f}, z₀={self.z0:.2f} → "
                f"α_S=α₀={self.alpha0:.4f}, β_S=β₀={self.beta0:.5f}; "
                f"z̄*_{{d,D}}={self._fmt(self.zs_d_D)}, "
                f"z̄_{{n,D}}={self._fmt(self.z_n_D)})")


# ============================================================
# 命令行自检
# ============================================================
if __name__ == "__main__":
    print("== DSMK 自检 (HSG 量级示例微剂量学量) ==")
    dsmk = DSMK(alpha0=0.174, beta0=0.0568, z0=66.0,
                z_d_D=54.1, zs_d_D=32.1, zs_d_F=15.0, z_n_D=0.107)
    print(dsmk)
    print("-" * 50)
    print(dsmk.summary())
    print("-" * 50)
    D = np.array([0.0, 0.5, 1.0, 2.0, 5.0])
    S = dsmk.survival(D, N=2000)
    print("D (Gy):", D.tolist())
    print("S(D)  :", np.round(S, 5).tolist())
