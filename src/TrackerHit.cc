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
/// \file TrackerHit.cc
/// \brief TrackerHit 类的实现：敏感探测器命中类
///
/// 记录细胞核内的能量沉积命中（Hit），包含径迹 ID、能量沉积等字段。
/// 提供 Draw()（可视化）与 Print()（打印）方法，供可视化与回放调用。

#include "TrackerHit.hh"

#include "G4Circle.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 线程局部的命中对象分配器（Geant4 内存池，提升命中创建/销毁效率）
G4ThreadLocal G4Allocator<TrackerHit>* TrackerHitAllocator = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerHit::TrackerHit() : G4VHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerHit::~TrackerHit() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerHit::TrackerHit(const TrackerHit& right) : G4VHit()
{
  /// 拷贝构造函数：从另一个命中对象复制字段。
  /// @param right 被拷贝的命中对象

  fTrackID = right.fTrackID;
  fEdep = right.fEdep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const TrackerHit& TrackerHit::operator=(const TrackerHit& right)
{
  /// 赋值运算符：从另一个命中对象复制字段。
  /// @param right 被赋值的命中对象
  /// @return 本对象的引用

  fTrackID = right.fTrackID;
  fEdep = right.fEdep;

  return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerHit::operator==(const TrackerHit& right) const
{
  /// 相等判定：基于对象地址比较（仅当为同一实例时才视为相等）。
  /// @param right 待比较的命中对象
  /// @return 两对象为同一实例时返回 true，否则 false

  return this == &right;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerHit::Draw()
{
  /// 可视化绘制本命中：以半透明红色实心圆点标记命中位置。
  /// 由 Geant4 可视化管理器在事件回放时调用。

  // —— 获取可视化管理器实例 ——
  G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();
  if (visManager) {
    // —— 构造命中位置的圆形标记 ——
    G4Circle circle(GetPosition());
    circle.SetScreenSize(5.);                 // 屏幕尺寸（像素）
    circle.SetFillStyle(G4Circle::filled);    // 填充样式：实心

    // —— 颜色与可视化属性（半透明红色：RGBA = 1,0,0,0.5）——
    G4Colour colour(1.0, 0.0, 0.0, 0.5);
    G4VisAttributes visAttributes(colour);
    circle.SetVisAttributes(visAttributes);

    // —— 绘制标记 ——
    visManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerHit::Print()
{
  /// 打印本命中的信息（径迹 ID 与能量沉积）。
  /// 能量沉积带最佳单位输出（G4BestUnit）。

  G4cout << "  径迹ID: " << fTrackID
         << "  能量沉积: " << G4BestUnit(fEdep, "Energy") << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
