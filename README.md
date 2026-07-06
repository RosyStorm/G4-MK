# Ac-225 细胞存活的微剂量学模拟（Geant4-DNA）

本项目基于 Geant4-DNA 的 `microtrack` 示例改造，用于研究 **Ac-225 α 内照射对细胞存活率的影响**。
它属于整体研究的 **Geant4 模拟域**：用 Geant4-DNA 在液态水细胞几何中模拟带电粒子径迹，产出**单事件比能谱** $f_{d,1}(z_d)$、$f_{n,1}(z_n)$ 等微剂量学量，供后续 **MK / 修正 SMK 模型**（参考 Inaniwa & Kanematsu, Phys. Med. Biol. 63 (2018) 095011）计算细胞存活率。

> 参考：原始示例论文 Baratto-Roldan et al., Front. Phys. 9 (2021) 726787。

---

## 1. 当前状态

| 模块 | 状态 | 说明 |
|---|---|---|
| 几何（同心球细胞） | ✅ 完成（任务 1.1） | World→Cell→Nucleus，SD 置核 |
| 基线复现 | ✅ 完成（任务 0.1） | 质子 20 MeV，$y_D$、$z_D$ 与文献 <10% 一致 |
| Ac-225 α 源 | ⏳ 待实现（任务 2.1） | 当前源仍为质子 |
| 源区室配置（核/质/膜/胞外） | ⏳ 待实现（任务 2.2） | |
| α 物理输运验证 | ⏳ 待验证（任务 3.1） | |
| 双位点打分 + 逐事件 ntuple | ⏳ 待实现（任务 4.x） | 当前仅域级 $f(z_d)$ |
| ROOT 后处理 / MK 计算 | ⏳ 待实现（任务 5.x） | |

---

## 2. 目录结构与文件说明

```
microtrack/
├── microtrack.cc              主程序：创建 RunManager，注册用户初始化类，启动 UI
├── CMakeLists.txt             CMake 构建配置（依赖 Geant4 + UI/Vis）
├── README.md                  本文件
│
├── include/  与  src/         头文件与实现（成对）
│   ├── ActionInitialization        注册 PrimaryGeneratorAction / RunAction（含多线程 master/build）
│   ├── DetectorConstruction        ★ 几何：同心球细胞 World(水盒)→Cell(球)→Nucleus(球)；SD 挂到 Nucleus
│   ├── DetectorMessenger           几何 UI 命令（/mygeom/...）
│   ├── PhysicsList                 物理列表，默认 G4EmDNAPhysics_option2（可切换 DNA 各选项/Livermore/Penelope）
│   ├── PhysicsListMessenger        物理列表 UI 命令（/physics/addPhysics）
│   ├── PrimaryGeneratorAction      粒子枪（当前质子；将改造为 Ac-225 α）
│   ├── PrimaryGeneratorMessenger   粒子枪 UI 命令（/beam/position/Z0）
│   ├── RunAction                   ★ 分析管理器：定义直方图、开/写 root 文件、打印统计
│   ├── TrackerSD                   ★ 敏感探测器：记录核内能量沉积 hit；EndOfEvent 做加权抽样与微剂量学计算
│   └── TrackerHit                  单条 hit（位置 + 能量沉积），用于可视化
│
├── run.mac                    原始示例宏（旧盒形几何，已过时，保留作历史参考）
├── run_validation.mac         ★ 当前验证宏（同心球几何 + 质子源）
├── init_vis.mac               可视化初始化（交互模式）
├── vis.mac                    可视化绘制设置
└── processMicrotrack.C        ROOT 宏：为各直方图设置坐标轴标题
```

**核心三处（打分相关，后续任务会重点改动）**：
- `DetectorConstruction` — 几何与 SD 挂载
- `TrackerSD` — 能量沉积记录 + 微剂量学量计算（加权抽样）
- `RunAction` — 直方图定义与统计输出

---

## 3. 几何模型

同心球结构，全部材料暂为液态水（核成分后续可细化）：

```
World   盒形水，半边长 = R_cell + maxRange   （细胞外加 maxRange 水层以维持次级电子平衡）
 └─ Cell     球，半径 R_cell（细胞膜边界）   中心在原点
     └─ Nucleus   球，半径 R_n（细胞核）     中心在原点，挂敏感探测器
```

域（domain）级位点不建实体，由 `TrackerSD` 在核内**随机放置**半径 $r_d$ 的球进行采样（加权抽样法，见原始示例论文）。

---

## 4. 物理过程

默认 `G4EmDNAPhysics_option2`（Geant4-DNA 径迹结构，适用 DNA 尺度）。可通过 UI 切换：
`dna_opt1..dna_opt8`、`liv`(Livermore)、`penelope`、`em_standard_opt4`。

---

## 5. 编译与运行

### 依赖
- Geant4（本机 11.4.1，已含 UI/Vis/Analysis；RadioactiveDecay 数据已就绪）
- CMake、编译器（支持 C++17）

### 编译（out-of-source）
```bash
cd microtrack
mkdir -p build && cd build
cmake ..
make -j4
```

### 运行
批处理模式：
```bash
./build/microtrack run_validation.mac      # 当前验证宏（同心球 + 质子）
```
交互模式（带可视化）：
```bash
./build/microtrack                          # 进入 UI，自动加载 init_vis.mac
```

多线程：在宏中 `/run/numberOfThreads N` 设置（默认按 `argv[2]` 或宏）。

> ⚠ 本机内存较小（~8 GB、无 swap），`make` 与多线程运行都不要把并行度开太高，避免 OOM。

---

## 6. 参数调整（UI 命令）

所有命令均可写在宏文件里，或交互模式下输入。需在 `/run/initialize` 之前设置几何/物理类参数。

### 几何（/mygeom/...）
| 命令 | 默认 | 含义 |
|---|---|---|
| `/mygeom/material <name>` | `G4_WATER` | 世界/细胞材料（NIST 名） |
| `/mygeom/maxRange <L>` | 8.25 um | 次级电子最大射程；决定细胞外水层厚度（世界半边长 = R_cell + maxRange） |
| `/mygeom/cellRadius <L>` | 10 um | 细胞（膜边界）半径 $R_{cell}$ |
| `/mygeom/nucleusRadius <L>` | 8 um | 细胞核半径 $R_n$（域采样区，须 >0 且 < R_cell） |
| `/mygeom/siteRadius <L>` | 0.5 um | 域(site)半径 $r_d$（须 ≤ R_n） |

### 物理（/physics/...）
| 命令 | 含义 |
|---|---|
| `/physics/addPhysics <name>` | 选择物理列表（见第 4 节） |
| `/run/setCutForAGivenParticle <part> <L>` | 产物射程截断（gamma/e-/e+/proton） |

### 粒子源（/gun/...、/beam/...）
| 命令 | 含义 |
|---|---|
| `/gun/particle <name>` | 粒子种类（当前 `proton`，后续 Ac-225 α） |
| `/gun/energy <E>` | 动能（外照射时可含 dE/dx 补偿） |
| `/beam/position/Z0 <L>` | 源沿 Z 坐标位置 |

### 运行控制
| 命令 | 含义 |
|---|---|
| `/run/numberOfThreads N` | 线程数 |
| `/run/initialize` | 初始化（几何/物理，须在设参数后、束流前） |
| `/run/beamOn N` | 模拟 N 个事件 |
| `/random/setSeeds a b` | 固定随机种子（可复现） |

---

## 7. 输出

运行产出 `microtrack.root`（默认 ROOT 格式，多线程自动合并），含直方图：

| 直方图 | 含义 |
|---|---|
| H1[0-2] | 单事件授予能 $\varepsilon$ 的 $f$、$\varepsilon f$、$\varepsilon^2 f$（→ $\bar\varepsilon_F$、$\bar\varepsilon_D$） |
| H1[3-5] | 线性能量(lineal energy) $y$ 的 $f$、$yf$、$y^2f$（→ $\bar y_F$、$\bar y_D$） |
| H1[6-8] | 比能 $z$ 的 $f$、$zf$、$z^2f$（→ $\bar z_F$、$\bar z_D$） |
| H1[9-11] | 每事件命中数 $N_{sel}$、$N_{site}$、$N_{int}$ |
| H1[12-13] | 核边界入射/出射动能 $T_{in}$、$T_{out}$ |
| H2[0] | $N_{site}$ vs $\varepsilon$ |

运行结束会在终端打印各直方图的频率均/剂量均（$\bar\varepsilon_F$、$\bar\varepsilon_D$、$\bar y_F$、$\bar y_D$、$\bar z_F$、$\bar z_D$）。

### 加权抽样原理（关键）
位点"绕一个随机 hit 随机放置"会引入命中密度偏置，重要性采样权重 $w = N_{sel}/N_{int}$ 将其校正为无偏的单事件谱。剂量均由加权直方图技巧得到：`FillH1(7, z, z·w)` 的 `mean()` 恰好等于 $\bar z_D = \int z^2 f\,dz / \int z f\,dz$。

---

## 8. 后续路线图

1. **任务 2.1** 质子枪 → Ac-225 α 源（4 个 α 各向同性，一个 α = 一个事件）
2. **任务 2.2** 源区室配置（核内/胞质/膜上/胞外）
3. **任务 3.1** α 在 DNA 物理中的射程/输运验证
4. **任务 4.x** 双位点打分（域 $f(z_d)$ + 核 $f(z_n)$）+ 逐事件配对 ntuple + 饱和修正
5. **任务 5.x** ROOT 后处理：$\bar z_{d,D}$、$\bar z^*_{d,D}$、$\bar z_{n,D}$、S 值、SMK 存活曲线
6. **任务 6.x** 验证：S 值对 MIRDcell、能量平衡、定位对比

---

## 9. 引用

使用本代码或其结果时，请引用 Geant4-DNA 协作组出版物：
- Med. Phys. 51 (2024) 5873–5889
- Med. Phys. 45 (2018) e722-e739
- Phys. Med. 31 (2015) 861-874
- Med. Phys. 37 (2010) 4692-4708
- Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178

以及示例参考论文：Baratto-Roldan et al., Front. Phys. 9 (2021) 726787。
SMK 模型：Inaniwa & Kanematsu, Phys. Med. Biol. 63 (2018) 095011。

Geant4-DNA 文档：http://geant4-dna.org
