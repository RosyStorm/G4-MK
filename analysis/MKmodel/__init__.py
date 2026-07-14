"""
MKmodel — 微剂量学动力学模型 (MKM/SMK/ModifiedSMK/DSMK) 包

层次:
    MK (基类, 经典 LQ 闭式)
    ├── SMK               (Sato 2012 单随机, 矩近似 MC)
    ├── ModifiedSMK       (Inaniwa 2018 闭式)
    └── DSMK              (Sato 2012 双随机 MC)

公开类:
    MK, SMK, ModifiedSMK, DSMK

公开 helper:
    load_microdosimetry_from_root   从 microtrack.root 抽 (z_d, z_n, w, n_hits, n_entries)
    compute_basic_microdosimetry    从样本算 z_d_F, z_d_D, z_n_F, z_n_D
"""

from .mk import MK, load_microdosimetry_from_root, compute_basic_microdosimetry
from .smk import SMK
from .msmk import ModifiedSMK
from .dsmk import DSMK

__all__ = [
    "MK",
    "SMK",
    "ModifiedSMK",
    "DSMK",
    "load_microdosimetry_from_root",
    "compute_basic_microdosimetry",
]
