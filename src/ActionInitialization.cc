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
/// \file ActionInitialization.cc
/// \brief ActionInitialization 类的实现：注册用户动作（User Action）
///
/// 负责在主线程（BuildForMaster）和工作线程（Build）中创建并注册
/// 各类用户动作对象：初级粒子产生、运行级、事件级、步级动作。

#include "ActionInitialization.hh"

#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::BuildForMaster() const
{
  /// 主线程入口：仅注册运行级动作（RunAction）。
  /// 多线程模式下主线程不处理事件，故无需注册事件级/步级动作。

  // —— 注册运行级动作（负责汇总、写入输出文件）——
  SetUserAction(new RunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitialization::Build() const
{
  /// 工作线程入口：注册完整的用户动作链。
  /// 注册顺序：初级粒子产生 -> 运行级 -> 事件级 -> 步级。

  // —— 初级粒子产生动作（定义每个事件的初始粒子）——
  SetUserAction(new PrimaryGeneratorAction);

  // —— 运行级动作（负责打开/写入分析输出）——
  SetUserAction(new RunAction);

  // —— 事件级动作（事件内累积量、事件末分析）——
  auto* eventAction = new EventAction;
  SetUserAction(eventAction);

  // —— 径迹级动作（任务X: 每条 track 启动时立即分类, 防 α 链衍生被拆）——
  SetUserAction(new TrackingAction(eventAction));

  // —— 步级动作（每步能量沉积、粒子杀死判定），依赖事件级动作 ——
  SetUserAction(new SteppingAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
