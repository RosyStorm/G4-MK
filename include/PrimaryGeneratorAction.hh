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
/// \file PrimaryGeneratorAction.hh
/// \brief PrimaryGeneratorAction 类的定义：初级粒子产生动作

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include <memory>

class PrimaryGeneratorMessenger;
class G4ParticleGun;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// @brief 初级粒子产生动作类
///
/// 支持两种源类型（/source/type）：
///   - proton：原质子枪（基线对比），位置/方向由成员决定，能量由 /gun/energy 给
///   - ac225 ：Ac-225 α 源，从指定区室（/source/compartment）均匀随机点各向同性
///             发射一个 α，能量按 Ac-225 衰变链抽样（一个 α = 一个事件，契合 MK 单事件定义）
///   - alpha ：单能 α 验证模式（任务3.1），能量由 /gun/energy 给
///   - ac225_decay       ：路线2 完整 4α+β+γ+反冲链 (G4RadioactiveDecay)
///   - ac225_single_decay：路线2 单次 Ac-225 → Fr-221 + α (5.83 MeV), 动量守恒反平行,
///                         不跑后续链
///   - lu177_decay       ：路线2 静止 Lu-177 β⁻ + 208/113 keV γ + 反冲
///   - carbon            ：单能 C-12 碳离子验证模式，能量由 /gun/energy 给
///                         LET 参考值：33.7 keV/μm（对应约 150 MeV/u，总动能约 1800 MeV）
///   - am241_decay       ：Am-241 特征 α 能谱直接抽样
///                         5.486 MeV(85.2%)/5.443 MeV(12.8%)，不模拟伴随 γ 或反冲核
///
/// 区室取值：Nucleus(核内) | Cytoplasm(质内) | Membrane(膜面,默认) | Extracellular(胞外) | WholeCell | CellExceptNucleus
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    // ===== 构造与析构 =====
    PrimaryGeneratorAction();                  // 创建枪与 Messenger，设定默认粒子
    virtual ~PrimaryGeneratorAction() override;  // 析构（枪与 Messenger 由智能指针释放）

    // ===== Geant4 强制重载接口 =====
    /// 为每个事件产生初级粒子。
    /// @param anEvent 当前事件指针
    void GeneratePrimaries(G4Event*) override;

    // ===== 查询接口 (GET) =====
    /// 取粒子枪指针（供其它动作查询粒子状态）。
    /// @return 粒子枪原始指针
    G4ParticleGun* GetParticleGun() { return fParticleGun.get(); }

    // ===== 设置接口 (SET)：由 PrimaryGeneratorMessenger 的 UI 命令调用 =====
    void SetPositionZ(G4double z) { fZ0 = z; }                 // 设置源点 Z 坐标（proton 模式）
    void SetSourceType(const G4String& t) { fSourceType = t; } // 设置源类型(proton|ac225|alpha)
    G4String GetSourceType() const { return fSourceType; }     // 取源类型
    void SetCompartment(const G4String& c) { fCompartment = c; }  // 设置源区室(ac225 模式)
    G4String GetCompartment() const { return fCompartment; }      // 取源区室
    void SetSSD(G4double v) { fSSD = v; }                        // 设置源到细胞表面距离(am241)
    G4double GetSSD() const { return fSSD; }                      // 取 SSD

  private:
    // ===== α 源辅助方法 =====
    G4double SampleAc225AlphaEnergy() const;                       // 抽样 Ac-225 衰变链 α 动能
    G4double SampleAm241AlphaEnergy() const;                       // 按 Am-241 主 α 分支比抽样动能
    G4ThreeVector SampleIsotropicDirection() const;                // 各向同性单位方向
    G4ThreeVector SampleConeDirection(G4double halfAngle) const;   // 以 (0,0,-1) 为轴、半角 halfAngle 的锥内均匀抽样方向
    G4ThreeVector SampleSourcePosition() const;                    // 按 fCompartment 抽样源点位置
    G4ThreeVector SampleInSphere(G4double R) const;                // 球内均匀抽样
    G4ThreeVector SampleInShell(G4double Rin, G4double Rout) const;  // 球壳内均匀抽样
    G4ThreeVector SampleOnSphere(G4double R) const;                // 球面均匀抽样
    G4ThreeVector SampleInBoxMinusSphere(G4double Rc, G4double wh) const;  // 盒内排除球(胞外)抽样

    // ===== 枪与交互命令 =====
    std::unique_ptr<G4ParticleGun> fParticleGun;              // 粒子枪
    std::unique_ptr<PrimaryGeneratorMessenger> fGunMessenger;  // UI 命令交互对象

    // ===== 源属性参数 =====
    G4ParticleDefinition* fParticle = nullptr;   // 质子（基线对比用）
    G4ParticleDefinition* fAlpha = nullptr;      // α 粒子
    G4ParticleDefinition* fAc225 = nullptr;      // Ac-225 离子（路线2 完整衰变链，任务7.1）
    G4ParticleDefinition* fLu177 = nullptr;      // Lu-177 离子（路线2 β⁻ 衰变链，任务8）
    G4ParticleDefinition* fFr221 = nullptr;      // Fr-221 反冲核(任务X 单次衰变, 路线2 单步)
    G4ParticleDefinition* fCarbon = nullptr;     // C-12 碳离子（单能模式，LET 参考值 33.7 keV/μm）
    G4String fSourceType = "ac225";              // proton | ac225 | alpha | ac225_decay | ac225_single_decay | lu177_decay | carbon | am241_decay
    G4String fCompartment = "Membrane";          // Nucleus | Cytoplasm | Membrane | Extracellular
    G4double fX0 = 0.;                           // proton 模式源点 X 坐标
    G4double fY0 = 0.;                           // proton 模式源点 Y 坐标
    G4double fZ0 = -10.2 * um;                   // proton 模式源点 Z 坐标
    G4double fSSD = 9.8 * mm;                     // Am-241 源到细胞表面距离(Source-to-Surface Distance)
    G4double fMomentumX = 0.;                    // proton 模式动量方向 X 分量
    G4double fMomentumY = 0.;                    // proton 模式动量方向 Y 分量
    G4double fMomentumZ = 1.;                    // proton 模式动量方向 Z 分量
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
