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
/// \file TrackerSD.hh
/// \brief TrackerSD 类的定义：细胞核敏感探测器

#ifndef TrackerSD_h
#define TrackerSD_h 1

#include "G4VSensitiveDetector.hh"

#include "TrackerHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// @brief 敏感探测器类
///
/// 挂载于细胞核逻辑体积，负责逐事件的命中收集与微剂量学打分：
/// 在 Initialize 中创建命中集合；在 ProcessHits 中为每个有能量沉积的 step
/// 生成一个 TrackerHit；在 EndOfEvent 中执行域(site)随机抽样、计算线能 y、
/// 比能 z 与核比能 z_n，并填充直方图与逐事件 ntuple。
class TrackerSD : public G4VSensitiveDetector
{
  public:
    // ===== 构造与析构 =====
    /// 构造，设定探测器名称与命中集合名称。
    /// @param name 敏感探测器名称（对应细胞核体积）
    /// @param hitsCollectionName 命中集合名称
    /// @param depthIdx 触发深度索引（保留参数）
    TrackerSD(const G4String& name, const G4String& hitsCollectionName,
              G4int depthIdx = 0);
    ~TrackerSD() override;  // 析构

    // ===== 基类重载方法 =====
    /// 事件开始：创建命中集合并注册到本事件的命中容器。
    /// @param hitCollection 本事件的命中容器
    void Initialize(G4HCofThisEvent* hitCollection) override;
    /// 逐步命中处理：为有能量沉积的 step 生成命中。
    /// @param step 当前粒子步
    /// @param history 可触历史（未使用）
    /// @return 本步生成了命中返回 true
    G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
    /// 事件结束：域随机抽样、计算微剂量学量并填充直方图与 ntuple。
    /// @param hitCollection 本事件的命中容器（未使用）
    void EndOfEvent(G4HCofThisEvent* hitCollection) override;

    // ===== 设置接口 (SET) =====
    inline void SetParID(const G4int& parID) { fParID = parID; }  // 设置粒子 ID 过滤值

    // ===== 查询接口 (GET) =====
    inline G4int GetParID() const { return fParID; }  // 取粒子 ID 过滤值

  private:
    // ===== 命中集合与过滤参数 =====
    TrackerHitColl* fHitsCollection = nullptr;  // 本事件的命中集合
    G4int fParID = 0;                           // 粒子 ID（过滤用）
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
