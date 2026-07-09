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
/// \file DetectorConstruction.hh
/// \brief DetectorConstruction 类的定义：构建细胞微剂量学几何

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <memory>

class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// @brief 探测器构建类
///
/// 负责定义细胞微剂量学模拟的几何结构：同心球嵌套的水盒-细胞-细胞核，
/// 并将敏感探测器(SD)挂载到细胞核。运行时参数可通过 DetectorMessenger 的 UI 命令修改。
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    // ===== 构造与析构 =====
    DetectorConstruction();             // 创建探测器并注册交互命令
    ~DetectorConstruction() override;   // 析构（Messenger 由智能指针自动释放）

    // ===== Geant4 强制重载接口 =====
    G4VPhysicalVolume* Construct() override;   // 构建几何树
    void ConstructSDandField() override;       // 创建并挂载敏感探测器

    // ===== 设置接口 (SET)：由 DetectorMessenger 的 UI 命令调用 =====
    void SetMaterial(const G4String& name);                                  // 设置介质材料
    void SetMaxRange(const G4double& range) { fMaxRange = range; }           // 次级电子最大射程
    void SetCellRadius(const G4double& r) { fCellRadius = r; }               // 细胞半径 R_cell
    void SetNucleusRadius(const G4double& r) { fNucleusRadius = r; }         // 细胞核半径 R_n
    void SetSiteRadius(const G4double& siteRadius) { fSiteRadius = siteRadius; } // 域半径 r_d
    void SetKillOutsideCell(G4bool v) { fKillOutsideCell = v; }              // 是否杀死出细胞的粒子
    void SetKillAtNucleus(G4bool v) { fKillAtNucleus = v; }                  // 杀死半径取 R_n 还是 R_cell
    void PrintParameters(G4VPhysicalVolume*) const;   // 打印体积参数
    void CheckConsistency();                          // 校验几何参数一致性

    // ===== 查询接口 (GET) =====
    G4String GetMaterial() const   // 当前材料名称（未定义时返回 "undefined"）
    {
      return fMat ? fMat->GetName() : G4String("undefined");
    }
    G4double GetMaxRange() const { return fMaxRange; }           // 次级电子最大射程
    G4double GetCellRadius() const { return fCellRadius; }       // 细胞半径
    G4double GetNucleusRadius() const { return fNucleusRadius; } // 细胞核半径
    G4double GetSiteRadius() const { return fSiteRadius; }       // 域半径
    G4bool GetKillOutsideCell() const { return fKillOutsideCell; }
    G4bool GetKillAtNucleus() const { return fKillAtNucleus; }

  private:
    // ===== 私有构建方法 =====
    void DefineMaterials();                                   // 定义介质材料
    G4VPhysicalVolume* DefineWorld();                        // 构建世界体积
    G4VPhysicalVolume* DefineCell(G4VPhysicalVolume* mother);      // 构建细胞体积
    G4VPhysicalVolume* DefineNucleus(G4VPhysicalVolume* mother);   // 构建细胞核体积

    std::unique_ptr<DetectorMessenger> fDetectorMessenger;  // 交互命令对象（智能指针管理）

    // ===== 几何与材料参数 =====
    G4Material* fMat = nullptr;        // 介质指针（默认 G4_WATER）
    G4double fMaxRange = 8.25 * um;    // 次级电子最大射程（决定细胞外水层厚度）
    G4double fCellRadius = 10. * um;   // 细胞(膜)半径 R_cell
    G4double fNucleusRadius = 8. * um; // 细胞核半径 R_n
    G4double fSiteRadius = 0.5 * um;   // 域(site)半径 r_d
    G4bool fKillOutsideCell = true;    // 加速：出细胞且向外运动的粒子 kill 掉(产线默认开,射程验证需关)
    G4bool fKillAtNucleus = true;      // kill 半径取 R_n(默认,更快) 还是 R_cell
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
