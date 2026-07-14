//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file PrimaryGeneratorAction.cc
/// \brief PrimaryGeneratorAction 类的实现：初级粒子生成器
///
/// 支持三种源类型(/source/type)：
///   - proton : 原质子枪(基线对比)，位置/方向由成员决定，能量由 /gun/energy 给
///   - alpha  : 单能 α 验证模式(任务3.1)，位置按区室、方向各向同性
///   - ac225  : Ac-225 α 源，从指定区室(/source/compartment)均匀随机点各向同性
///              发射一个 α，能量按 Ac-225 衰变链抽样(一个 α = 一个事件，契合 MK 单事件定义)
///   区室: Nucleus(核内) | Cytoplasm(质内) | Membrane(膜面,默认) | Extracellular(胞外)

#include "PrimaryGeneratorAction.hh"

#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"

#include <cmath>

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
{
  /// 构造函数：创建交互命令对象与粒子枪，设置默认粒子(质子/α)。

  // 创建本类的交互命令对象
  fGunMessenger = std::make_unique<PrimaryGeneratorMessenger>(this);

  // 创建粒子枪（每事件发射 1 个粒子）
  G4int nOfParticles = 1;
  fParticleGun = std::make_unique<G4ParticleGun>(nOfParticles);

  // —— 设置默认粒子定义 ——
  fParticle = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fAlpha = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
  // 注: Ac-225 离子(Z=89,A=225)不在构造函数里取——此时 MT worker 上 GenericIon
  //     尚未就绪会触发 PART105 并返回 nullptr。改为在 GeneratePrimaries 里延迟取(任务7.1)。
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  /// 事件开始时生成初级粒子。
  /// 根据 fSourceType 选择源模式：ac225 / alpha / proton。
  /// @param anEvent 当前事件指针

  if (fSourceType == "ac225") {
    // —— Ac-225 α 源：区室抽样位置 + 各向同性 + 抽样能量 ——
    G4double ekin = SampleAc225AlphaEnergy();
    G4ThreeVector pos = SampleSourcePosition();
    G4ThreeVector dir = SampleIsotropicDirection();

    fParticleGun->SetParticleDefinition(fAlpha);
    fParticleGun->SetParticleEnergy(ekin);
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(dir);

    // 诊断：打印前若干事件的抽样参数(验证区室位置 + 各向同性 + 能量)
    if (anEvent->GetEventID() < 10) {
      G4cout << "[Ac225] event " << anEvent->GetEventID()
             << "  comp=" << fCompartment
             << "  Ekin=" << ekin / MeV << " MeV"
             << "  |pos|=" << pos.mag() / um << " um"
             << "  pos(um)=" << pos / um
             << "  dir=" << dir << G4endl;
    }
  }
  else if (fSourceType == "alpha") {
    // —— 单能 α 验证模式(任务3.1)：能量由 /gun/energy 给，位置按区室，各向同性 ——
    G4ThreeVector pos = SampleSourcePosition();
    G4ThreeVector dir = SampleIsotropicDirection();

    fParticleGun->SetParticleDefinition(fAlpha);
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleMomentumDirection(dir);

    if (anEvent->GetEventID() < 5) {
      G4cout << "[alpha] event " << anEvent->GetEventID()
             << "  comp=" << fCompartment
             << "  |pos|=" << pos.mag() / um << " um"
             << "  dir=" << dir << G4endl;
    }
  }
  else if (fSourceType == "ac225_decay") {
    // —— 路线2(任务7.1): 完整衰变链 ——
    // 静止的 Ac-225 离子置于源点，由 G4RadioactiveDecay 自动衰变(4α+β+γ+反冲核)。
    // 能量=0(原子静止)，方向无意义；位置按区室抽样(与路线1一致，便于对比)。
    // 延迟取 Ac-225 离子定义(构造函数里取会因 GenericIon 未就绪而失败)
    if (!fAc225) {
      fAc225 = G4IonTable::GetIonTable()->GetIon(89, 225, 0.0);
    }
    if (!fAc225) {
      G4ExceptionDescription msg;
      msg << "Ac-225 离子(Z=89,A=225)无法创建 — 检查 PhysicsList 是否构造了离子。";
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "PGA001",
                  FatalException, msg);
    }
    G4ThreeVector pos = SampleSourcePosition();
    fParticleGun->SetParticleDefinition(fAc225);
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleEnergy(0.0);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

    if (anEvent->GetEventID() < 5) {
      G4cout << "[Ac225-decay] event " << anEvent->GetEventID()
             << "  comp=" << fCompartment
             << "  |pos|=" << pos.mag() / um << " um"
             << "  (Ac-225 at rest → full chain 4α+β+γ+recoil)" << G4endl;
    }
  }
  else if (fSourceType == "lu177_decay") {
    // —— 路线2(任务8): Lu-177 β⁻ 衰变链 ———
    // 静止 Lu-177 离子置于源点, G4RadioactiveDecay 自动跑 β⁻ (max 497 keV) +
    // 208 keV (10%) / 113 keV (6.2%) 退激 γ + 反冲核。
    // 几何保持 Ac-225 路线2 一致 (maxRange=100 µm) 以便直接对照 ———
    // 但要注意: β⁻ 在水中射程 ~0.5 mm = 500 µm, 因此核外 β 会被选择性 kill (Z≤2)
    // 截断, 即只能看到"核内+领域内"的 β 贡献。这是有意的: 第一版先确定核内
    // 剂量学下限, 再考虑扩展 maxRange 到 1000 µm 跑完整射程版。
    if (!fLu177) {
      fLu177 = G4IonTable::GetIonTable()->GetIon(71, 177, 0.0);
    }
    if (!fLu177) {
      G4ExceptionDescription msg;
      msg << "Lu-177 离子(Z=71,A=177)无法创建 — 检查 PhysicsList 是否构造了离子。";
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries", "PGA003",
                  FatalException, msg);
    }
    G4ThreeVector pos = SampleSourcePosition();
    fParticleGun->SetParticleDefinition(fLu177);
    fParticleGun->SetParticlePosition(pos);
    fParticleGun->SetParticleEnergy(0.0);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

    if (anEvent->GetEventID() < 5) {
      G4cout << "[Lu177-decay] event " << anEvent->GetEventID()
             << "  comp=" << fCompartment
             << "  |pos|=" << pos.mag() / um << " um"
             << "  (Lu-177 at rest → β⁻ + 208/113 keV γ + recoil)" << G4endl;
    }
  }
  else {
    // —— proton 模式：保持原有基线行为 ——
    fParticleGun->SetParticlePosition(G4ThreeVector(fX0, fY0, fZ0));
    fParticleGun->SetParticleMomentumDirection(
      G4ThreeVector(fMomentumX, fMomentumY, fMomentumZ));
  }

  // 生成初级顶点（在事件开始时调用）
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double PrimaryGeneratorAction::SampleAc225AlphaEnergy() const
{
  /// 按 Ac-225 衰变链抽样一个 α 的动能。
  /// @return 抽样得到的 α 动能(MeV)
  ///
  /// Ac-225 衰变链 4 个 α（能量据 NNDC，实现时已核对）：
  ///   1) 225Ac -> 221Fr : 5.83 MeV  (100%)
  ///   2) 221Fr -> 217At : 6.36 MeV  (100%)
  ///   3) 217At -> 213Bi : 7.07 MeV  (100%)
  ///   4) 经 213Bi: 97.84% -> 213Po -> 209Pb (8.40 MeV)
  ///                 2.16% -> 209Tl 直接 alpha (5.87 MeV)
  /// 一次衰变产生 4 个 α(各槽位各 1 个)，故 4 槽等概率抽样即可重建谱形。

  // 4 个 α 槽位的特征能量
  const G4double eAc = 5.83 * MeV;  // 225Ac -> 221Fr
  const G4double eFr = 6.36 * MeV;  // 221Fr -> 217At
  const G4double eAt = 7.07 * MeV;  // 217At -> 213Bi
  const G4double ePo = 8.40 * MeV;  // 213Po -> 209Pb
  const G4double eBi = 5.87 * MeV;  // 209Tl 直接 α

  // 4 槽等概率抽样（0..3）
  G4int slot = static_cast<G4int>(G4UniformRand() * 4.0);  // 0..3 等概率
  if (slot == 0) return eAc;
  if (slot == 1) return eFr;
  if (slot == 2) return eAt;
  return (G4UniformRand() < 0.9784) ? ePo : eBi;  // slot==3, 含 Bi-213 分支(分支比 97.84%)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::SampleIsotropicDirection() const
{
  /// 抽样一个各向同性单位方向向量。
  /// @return 各向同性单位方向（cosθ 在 [-1,1] 均匀，φ 在 [0,2π) 均匀）

  // 球面均匀采样：极角余弦与方位角均均匀分布
  G4double cosTheta = 2. * G4UniformRand() - 1.;               // cosθ ∈ [-1,1]
  G4double sinTheta = std::sqrt((1. - cosTheta) * (1. + cosTheta));
  G4double phi = CLHEP::twopi * G4UniformRand();               // 方位角 φ ∈ [0,2π)
  return G4ThreeVector(sinTheta * std::cos(phi), sinTheta * std::sin(phi),
                       cosTheta);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector
PrimaryGeneratorAction::SampleSourcePosition() const
{
  /// 按 fCompartment 在对应区室内均匀抽样源点位置。
  /// @return 源点位置（核内/质内/膜面/胞外，取决于当前区室设置）

  // 从探测器构造获取几何尺寸
  const auto* det = static_cast<const DetectorConstruction*>(
    G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double Rn = det ? det->GetNucleusRadius() : 8. * um;     // 核半径(默认 8 um)
  G4double Rc = det ? det->GetCellRadius() : 10. * um;       // 细胞半径(默认 10 um)
  G4double Rworld = Rc + (det ? det->GetMaxRange() : 49. * um);  // 世界半边长

  // 按区室选择对应的抽样方式
  if (fCompartment == "Nucleus")           return SampleInSphere(Rn);
  if (fCompartment == "Cytoplasm")         return SampleInShell(Rn, Rc);
  if (fCompartment == "Extracellular")     return SampleInBoxMinusSphere(Rc, Rworld);
  // CellExceptNucleus: 细胞体均匀, 但排除 Nucleus (即 Cytoplasm ∪ Membrane, 不含核)
  if (fCompartment == "CellExceptNucleus") {
    G4ThreeVector p;
    do {
      p = SampleInSphere(Rc);
    } while (p.mag() < Rn);   // 拒绝落入核内的点
    return p;
  }
  // WholeCell: 整个 Rc 球内均匀 (Nucleus ∪ Cytoplasm ∪ Membrane), 等价于 SampleInSphere(Rc)
  if (fCompartment == "WholeCell") {
    return SampleInSphere(Rc);
  }
  // 默认 / Membrane：细胞膜面
  return SampleOnSphere(Rc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::SampleInSphere(G4double R) const
{
  /// 在半径 R 的球内均匀抽样一点。
  /// @param R 球半径
  /// @return 球内均匀分布的随机点

  // 球内均匀：r ~ R * u^(1/3)，方向各向同性
  G4double r = R * std::pow(G4UniformRand(), 1. / 3.);
  return r * SampleIsotropicDirection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector
PrimaryGeneratorAction::SampleInShell(G4double Rin, G4double Rout) const
{
  /// 在球壳 [Rin, Rout] 内均匀抽样一点。
  /// @param Rin 球壳内半径
  /// @param Rout 球壳外半径
  /// @return 球壳内均匀分布的随机点

  // 球壳内均匀：r^3 在 [Rin^3, Rout^3] 均匀
  G4double rin3 = Rin * Rin * Rin;
  G4double rout3 = Rout * Rout * Rout;
  G4double r = std::pow(rin3 + G4UniformRand() * (rout3 - rin3), 1. / 3.);
  return r * SampleIsotropicDirection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector PrimaryGeneratorAction::SampleOnSphere(G4double R) const
{
  /// 在半径 R 的球面上均匀抽样一点。
  /// @param R 球半径
  /// @return 球面上均匀分布的随机点

  // 球面均匀 = R * 各向同性单位方向
  return R * SampleIsotropicDirection();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreeVector
PrimaryGeneratorAction::SampleInBoxMinusSphere(G4double Rc, G4double wh) const
{
  /// 在半边长 wh 的盒内、排除半径 Rc 的球（胞外区）中均匀抽样一点。
  /// 采用拒绝采样法。
  /// @param Rc 排除球的半径（细胞半径）
  /// @param wh 盒的半边长
  /// @return 盒内且球外的随机点

  // 拒绝采样：在盒内均匀取点，若落在排除球内则重试
  for (G4int i = 0; i < 1000; ++i) {
    G4double x = (2. * G4UniformRand() - 1.) * wh;
    G4double y = (2. * G4UniformRand() - 1.) * wh;
    G4double z = (2. * G4UniformRand() - 1.) * wh;
    G4ThreeVector p(x, y, z);
    if (p.mag() > Rc) return p;
  }
  return G4ThreeVector(wh, 0., 0.);  // fallback(几乎不会到达)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
