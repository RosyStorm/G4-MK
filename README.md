# Ac-225 细胞存活的微剂量学模拟（Geant4-DNA + MK/SMK/DSMK）

本项目基于 Geant4-DNA 的 `microtrack` 示例改造，端到端模拟 **Ac-225 α 内照射对细胞存活率的影响**：
**Geant4-DNA 径迹结构模拟**（同心球细胞几何 + Ac-225 α 源 + 双位点打分）→ **逐事件 ntuple** →
**ROOT 后处理**（微剂量学量 + MK/SMK/DSMK 存活曲线 + S 值 + 能量平衡）。

模型参考：
- **DSMK / SMK**：Sato & Furusawa, Radiat. Res. 178 (2012) 341–356
- **修正 SMK（闭式）**：Inaniwa & Kanematsu, Phys. Med. Biol. 63 (2018) 095011
- 原始示例：Baratto-Roldan et al., Front. Phys. 9 (2021) 726787

---

## 1. 当前状态

| 模块 | 状态 | 关键结果 |
|---|---|---|
| 几何（同心球细胞） | ✅ | World→Cell→Nucleus，SD 置核 |
| Ac-225 α 源（路线1：手动4α） | ✅ | 5.83/6.36/7.07/8.40 MeV，各向同性，一α一事件 |
| 源区室（核/质/膜/胞外） | ✅ | `/source/compartment` 可切换 |
| α 物理输运 | ✅ | 射程对 NIST ASTAR；电荷交换修正 |
| 双位点打分 + 逐事件 ntuple | ✅ | 域 $z_d$ + 核 $z_n$ 配对，hit/miss 都记 |
| 速度优化 | ✅ | 出核即杀 + 16 线程，**51× 加速**（0.047 s/事件） |
| ROOT 后处理（MK/SMK/DSMK） | ✅ | 显式卷积 + 矩验证；S(D) 曲线 |
| 全局校验 | ✅ | 能量守恒 1.000；S(N←N)=0.21 Gy/dec；分布效应定量 |

---

## 2. 目录结构

```
microtrack/
├── microtrack.cc              主程序
├── CMakeLists.txt             CMake（GLOB 自动收集 src/*.cc）
├── README.md
│
├── include/ + src/            头文件与实现（成对）
│   ├── ActionInitialization        注册 PGA/RunAction/SteppingAction
│   ├── DetectorConstruction   ★    几何：World→Cell→Nucleus 同心球；SD 挂 Nucleus
│   ├── DetectorMessenger           /mygeom/ 命令（半径 + kill 开关）
│   ├── PhysicsList                 G4EmDNAPhysics_option2（可切换）
│   ├── PrimaryGeneratorAction ★    源：proton | ac225(膜/核/质/胞外) | alpha(单能)
│   ├── PrimaryGeneratorMessenger   /source/ /beam/ 命令
│   ├── RunAction              ★    直方图 + events ntuple 定义 + 统计打印
│   ├── TrackerSD              ★    核内打分：域 z_d + 核 z_n + 全局 edep；EndOfEvent 填 ntuple
│   ├── TrackerHit                  hit（位置+edep）
│   ├── SteppingAction        ★(新) 出核/出细胞即杀（加速）+ 全局 edep 累加 + α 链全链跟踪
│   └── EventAction           ★(新) 逐事件：α 射程 + 全局 edep + 边界步核内 edep 累加器
│
├── macro/                     G4 模拟输入（.mac 宏，按用途分类）
│   ├── run_validation.mac         质子基线（同心球，/source/type proton）
│   ├── run_ac225.mac              Ac-225 膜面源（默认几何，加速全开）
│   ├── run_hsg.mac          ★     HSG 一致几何（r_d=0.274, r_n=6.2 µm，Sato Table1）+ Ac-225
│   ├── run_alpha_range.mac        α 射程验证（单能，kill 关）
│   ├── run.mac                    原始示例宏（旧盒形，过时）
│   └── init_vis.mac / vis.mac     可视化
│
├── analysis/                  ROOT 后处理脚本（按模型/任务分类）
│   ├── analyze_dsmk.C       ★(主) DSMK 模型：显式卷积(Sato 式20-28) + 矩验证 + S(D) + 对比 mod-SMK/MK
│   ├── analyze_mk.C               mod-SMK 闭式(Inaniwa 式24) + 经典 MK
│   ├── task6_svalue_balance.C     6.1 S(N←N) + 6.2 能量平衡
│   ├── task6_compartments.C       6.3 四区室 z̄_{n,D} + DSMK S(D) 对比
│   ├── task6_summary.C            6.3 四区室每衰变核剂量 S(N←source)
│   └── processMicrotrack.C        直方图轴标题（旧）
│
├── data/                      G4 模拟输出（.root 文件，运行时生成）
│   └── microtrack*.root          microtrack.root (默认) + 四区室 (Nuc/Cyt/Mem/Ext)
│
├── result/                    ROOT 后处理输出（图片 + root）
│   ├── DSMK/                     analyze_dsmk.C 输出
│   ├── mod-SMK/                  analyze_mk.C 输出
│   └── compartment/              task6_compartments.C 输出
│
├── build/                     编译产物（cmake .. && make）
├── ref/                       PDF 文献（Sato 2012, Inaniwa 2018 等，不入库）
└── History                    项目历史
```

**核心打分链**：`PrimaryGeneratorAction`(源) → `TrackerSD`(核内 edep→z_d/z_n) → `SteppingAction`(kill加速+全局edep+边界补加) → `RunAction`(ntuple/直方图 → `data/microtrack.root`) → `analysis/analyze_dsmk.C`(后处理)。

---

## 3. 几何与加速

同心球，全为液态水：
```
World(水盒, 半边长=R_cell+maxRange) ─ Cell(球, R_cell) ─ Nucleus(球, R_n, 挂SD)
```
域 $r_d$ 不建实体，由 `TrackerSD` 在核内随机放球采样（加权 $w=N_{sel}/N_{int}$ 校正偏置）。

**加速**（`SteppingAction`）：粒子出核（`killAtNucleus`，kill 半径=R_n）或出细胞（`killOutsideCell`，=R_cell）且**向外运动**时直接 kill——核打分无需核外输运，砍掉 ~90% 计算。`pos·dir>0` 判定保留胞外源向内射入的 α（交叉照射不受影响）。**射程验证须关 kill**。

---

## 4. 物理

默认 `G4EmDNAPhysics_option2`（DNA 径迹结构）。可切换 `dna_opt1..8 / liv / penelope / em_standard_opt4`。

---

## 5. 编译与运行 ⭐（重点）

### 依赖
- Geant4 11.4.1（含 UI/Vis/Analysis + RadioactiveDecay 数据）
- CMake + C++17
- ROOT（在 conda env `microtrack` 内）—— 后处理用

### 编译（out-of-source）
```bash
cd microtrack
mkdir -p build && cd build && cmake .. && make -j4
# 增量: make -C build -j4
```
> ⚠ 本机 ~8 GB 内存、无 swap，`make` 与多线程不要开太高，避免 OOM。

### 运行（批处理，从 microtrack/ 根目录）
```bash
# 1) 质子基线（验证几何/打分）
./build/microtrack macro/run_validation.mac

# 2) Ac-225 膜面源（HSG 一致几何, 主配置）
./build/microtrack macro/run_hsg.mac

# 3) 切换源区室（改宏里 /source/compartment）：
sed 's|^/source/compartment.*|/source/compartment Nucleus|' macro/run_hsg.mac > /tmp/run_N.mac
./build/microtrack /tmp/run_N.mac    # Cytoplasm / Membrane / Extracellular 同理

# 4) α 射程验证（kill 必须关, 否则射程截断）
./build/microtrack macro/run_alpha_range.mac
```
交互模式（可视化）：`./build/microtrack`（加载 init_vis.mac）。

输出：`data/microtrack.root`（直方图 + events ntuple，多线程自动合并）。

### 后处理（ROOT，conda `microtrack` 环境，从 microtrack/ 根目录）
```bash
# 主分析：DSMK 存活曲线 + 微剂量学量 + 卷积矩验证（默认参数 Sato DSMK HSG）
conda run -n microtrack root -b -q analysis/analyze_dsmk.C
# 换细胞系参数:
conda run -n microtrack root -b -q 'analysis/analyze_dsmk.C("data/microtrack.root",0.156,0.0607,89.0,15)'
#   参数: (文件, α0[Gy⁻¹], β0[Gy⁻²], z0[Gy], Dmax[Gy])
#   Sato DSMK HSG 默认; Inaniwa HSG: 0.174,0.0568,66.0

# mod-SMK 闭式 + 经典 MK（Inaniwa 式24）
conda run -n microtrack root -b -q analysis/analyze_mk.C

# 任务6 校验:
conda run -n microtrack root -b -q analysis/task6_summary.C          # 四区室每衰变核剂量 S(N←src)
conda run -n microtrack root -b -q analysis/task6_svalue_balance.C   # S(N←N) + 能量平衡
conda run -n microtrack root -b -q analysis/task6_compartments.C     # 四区室 z̄_{n,D} + S(D) 对比图
```
产物写 `result/DSMK/`、`result/mod-SMK/`、`result/compartment/`（图 png/pdf + root）。

---

## 6. 参数调整（UI 命令）

几何/物理参数须在 `/run/initialize` **之前**；源类型/区室在 **之后**（随 PGA 构造）。

### 几何（/mygeom/...）
| 命令 | 默认 | 含义 |
|---|---|---|
| `/mygeom/cellRadius <L>` | 10 um | 细胞半径 R_cell |
| `/mygeom/nucleusRadius <L>` | 8 um | 核半径 R_n（HSG 用 6.2） |
| `/mygeom/siteRadius <L>` | 0.5 um | 域半径 r_d（HSG 用 0.274） |
| `/mygeom/maxRange <L>` | 8.25 um | 细胞外水层（世界半边长=R_cell+maxRange） |
| `/mygeom/material <name>` | G4_WATER | 材料 |
| `/mygeom/killOutsideCell <bool>` | true | 出细胞且向外→kill（产线加速；射程验证须关） |
| `/mygeom/killAtNucleus <bool>` | true | kill 半径取 R_n（更快，对核打分无偏） |

### 源（/source/...、/beam/...）
| 命令 | 含义 |
|---|---|
| `/source/type proton\|ac225\|alpha` | proton=基线；ac225=Ac-225 4α（默认）；alpha=单能α（/gun/energy 给） |
| `/source/compartment Nucleus\|Cytoplasm\|Membrane\|Extracellular` | ac225 源位置（均匀采样） |
| `/beam/position/Z0 <L>` | proton 模式源 Z 位置 |

### 物理 + 运行控制
| 命令 | 含义 |
|---|---|
| `/physics/addPhysics <name>` | 物理列表 |
| `/run/setCutForAGivenParticle <part> <L>` | 产物截断 |
| `/run/numberOfThreads N` | 线程数（建议 16） |
| `/run/initialize` → `/run/beamOn N` | 初始化后束流 |
| `/random/setSeeds a b` | 固定种子 |

---

## 7. 输出

### data/microtrack.root 直方图
| | 含义 |
|---|---|
| H1[0-2] | 域授予能 $\varepsilon_d$（→ $\bar\varepsilon_{d,F/D}$） |
| H1[3-5] | 线性能量(lineal energy) $y$（→ $\bar y_{F/D}$） |
| H1[6-8] | 域比能 $z_d$（→ $\bar z_{d,F/D}$） |
| H1[9-11] | 命中数 $N_{sel}/N_{site}/N_{int}$ |
| H1[12-13] | 核边界动能 $T_{in}/T_{out}$（⚠ α 下 fGeomBoundary 不稳，仅诊断） |
| H1[14] | alphaRange（α 投影射程） |
| H1[15-17] | 核比能 $z_n$（→ $\bar z_{n,F/D}$） |
| **events ntuple** | 逐事件 13 列（见下） |

### events ntuple（任务4.2，每事件一行，hit/miss 都记）
`eventID, alphaE_MeV, edep_d_keV, z_d_Gy, edep_n_keV, z_n_Gy, weight, nHsel, nHsite, nHint, hitFlag, compartment, edep_total_keV`

后处理由此算：边缘谱、联合 $f(z_d,z_n)$、条件分布、多事件卷积、SMK/DSMK $S(D)$、S 值、能量平衡。

### 加权抽样原理
位点"绕随机 hit 放置"引入命中密度偏置，权 $w=N_{sel}/N_{int}$ 校正为无偏单事件谱。剂量均由加权直方图技巧：`FillH1(7,z,z·w).mean()` = $\bar z_D$。

---

## 8. MK / SMK / DSMK 模型（后处理）

- **DSMK**（`analysis/analyze_dsmk.C`，最严格）：Sato 式7/19，饱和 $z'_d=z_0\sqrt{1-e^{-(z_d/z_0)^2}}$（式16，有 sqrt）。多事件分布用复合 Poisson 卷积（式20-28，MC 实现），**矩验证**：$\langle z_d\rangle=z_n$、$\langle z_d^2\rangle=\bar z_{d,D}z_n+z_n^2$、$\langle z_n\rangle=D$、$\langle z_n^2\rangle=\bar z_{n,D}D+D^2$。
- **mod-SMK**（`analysis/analyze_mk.C`，闭式）：Inaniwa 式24，饱和 $z^*=z_0(1-e^{-(z/z_0)^2})$（式1，无 sqrt），只需单事件矩。
- **经典 MK**：$\ln S=-(\alpha_0+\beta_0\bar z^*_{d,D})D-\beta_0 D^2$。

**HSG 参数**（Sato2012 Table1）：DSMK α0=0.156, β0=0.0607, r_d=0.274, r_n=6.2, z0=89。
Inaniwa2018：α0=0.174, β0=0.0568, r_d=0.28, R_n=8.1, z0=66。

---

## 9. 验证结果摘要（任务6）

| 校验 | 结果 |
|---|---|
| 能量守恒（**killOutsideCell=false, killAtNucleus=false**，即完全输运） | `<edep_total>/<alphaE>` = **1.000** ✓（生产配置 kill=ON 时，<edep_total>/<alphaE><1，因出核 α 被 kill 未沉积剩余能量） |
| S(N←N) Ac-225（r_n=6.2µm） | **0.21 Gy/dec**（自洽；MIRDcell 精确对标需软件） |
| 卷积正确性 | 复合 Poisson 矩关系全部吻合 ✓ |
| 每-Gy 存活（HSG, DSMK） | D@S=0.1 ≈ 1.0 Gy；RBE₁₀ ≈ 5（α 文献 3–5） |
| 分布效应 S(N←src)/dec | N 0.21 > Cy 0.092 > Mem 0.067 > Ext ~0 → **核内化比膜结合致命 ~3×** |

---

## 10. 引用

Geant4-DNA：Med. Phys. 51 (2024) 5873；Med. Phys. 45 (2018) e722；Phys. Med. 31 (2015) 861；Med. Phys. 37 (2010) 4692。
示例：Baratto-Roldan et al., Front. Phys. 9 (2021) 726787。
DSMK/SMK：Sato & Furusawa, Radiat. Res. 178 (2012) 341。
修正 SMK：Inaniwa & Kanematsu, Phys. Med. Biol. 63 (2018) 095011。
MIRDcell：Goddu et al., MIRD Cellular S-Values (1997)；MIRD Pamphlet 25。
