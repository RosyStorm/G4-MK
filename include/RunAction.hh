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
/// \file RunAction.hh
/// \brief RunAction 类的定义：逐运行动作

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Run;

/// @brief 逐运行动作类
///
/// 在构造时创建分析管理器并定义全部直方图（能量沉积 ε、线能 y、比能 z、
/// 核比能 z_n 等）与逐事件 ntuple；运行开始时打开输出文件，运行结束时
/// 打印直方图统计量（频率平均、剂量平均）并写出/关闭文件。
class RunAction : public G4UserRunAction
{
  public:
    // ===== 构造与析构 =====
    RunAction();   // 创建分析管理器、定义直方图与 ntuple
    ~RunAction() override;  // 默认析构

    // ===== Geant4 强制重载接口 =====
    void BeginOfRunAction(const G4Run*) override;  // 运行开始：打开输出文件
    void EndOfRunAction(const G4Run*) override;    // 运行结束：打印统计量并写出文件

    // ===== 静态接口（由 microtrack.cc 在 ApplyCommand 之前从宏文件里抽出后注入）=====
    // 背景：master 线程上 GetUserPrimaryGeneratorAction() 返回 nullptr（MT 模式下
    //       PrimaryGeneratorAction 由 Build() 在 worker 线程里创建）。
    //       因此 BeginOfRunAction 不能依赖 PGA 拿源类型，必须靠主线程从 argv[1] 解析宏
    //       文件静态注入。同源类型由 microtrack.cc 在 ApplyCommand() 之前一次性写入。
    static void SetRunMeta(G4String sourceType, G4String compartment);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
