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
/// \file DetectorMessenger.hh
/// \brief DetectorMessenger 类的定义：探测器几何的 UI 命令交互

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

#include <memory>
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// @brief 探测器交互命令类
///
/// 定义 /mygeom/ 命令目录下的 UI 命令，使用户可在运行时（宏文件或交互式会话）
/// 修改探测器几何参数：材料、细胞/细胞核/域半径、次级电子最大射程、以及
/// 出细胞 kill 加速开关。命令转发给 DetectorConstruction 的 SET 方法。
class DetectorMessenger : public G4UImessenger
{
  public:
    DetectorMessenger(DetectorConstruction*);   // 注册 /mygeom/ 命令
    ~DetectorMessenger() override;              // 析构（命令对象由智能指针释放）

    /// 响应 UI 命令赋值，转发给 DetectorConstruction 对应的 SET 方法。
    /// @param command 触发的 UI 命令指针
    /// @param newValue 用户输入的新值字符串
    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    // ===== 命令目标与命令对象 =====
    DetectorConstruction* fDetector = nullptr;                          // 关联的探测器构建对象
    std::unique_ptr<G4UIdirectory> fDetectorDir;                        // /mygeom/ 命令目录
    std::unique_ptr<G4UIcmdWithADoubleAndUnit> fMaxRangeCmd;            // 次级电子最大射程
    std::unique_ptr<G4UIcmdWithAString> fMatNameCmd;                    // 介质材料名称
    std::unique_ptr<G4UIcmdWithADoubleAndUnit> fCellRadiusCmd;          // 细胞半径 R_cell
    std::unique_ptr<G4UIcmdWithADoubleAndUnit> fNucleusRadiusCmd;       // 细胞核半径 R_n
    std::unique_ptr<G4UIcmdWithADoubleAndUnit> fSiteRadiusCmd;          // 域半径 r_d
    std::unique_ptr<G4UIcmdWithABool> fKillOutsideCellCmd;              // 出细胞向外粒子 kill 开关
    std::unique_ptr<G4UIcmdWithABool> fKillAtNucleusCmd;                // kill 半径取 R_n 还是 R_cell
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
