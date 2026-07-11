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
/// \file TrackerHit.hh
/// \brief TrackerHit 类的定义：单步能量沉积命中

#ifndef TrackerHit_h
#define TrackerHit_h 1

#include "G4ParticleDefinition.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4VHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// @brief 敏感探测器命中类
///
/// 记录细胞核内单个 step 的能量沉积信息：径迹 ID、命中位置、能量沉积、
/// 以及产生该沉积的粒子定义。支持 Geant4 的自定义分配器以提高命中对象
/// 的创建/销毁效率（多线程下每个线程独立分配器）。
class TrackerHit : public G4VHit
{
  public:
    // ===== 构造与析构 =====
    TrackerHit();                  // 默认构造
    TrackerHit(const TrackerHit&); // 拷贝构造
    ~TrackerHit() override;        // 析构

    // ===== 运算符 =====
    const TrackerHit& operator=(const TrackerHit& right);  // 赋值运算符
    G4bool operator==(const TrackerHit& right) const;      // 相等判定（按地址）

    inline void* operator new(size_t);   // 自定义 new（走 Geant4 分配器）
    inline void operator delete(void*);  // 自定义 delete（走 Geant4 分配器）

    // ===== 基类重载方法 =====
    void Draw() override;   // 可视化绘制命中（红色圆点）
    void Print() override;  // 打印命中信息（径迹 ID 与能量沉积）

    // ===== 设置接口 (SET) =====
    inline void SetTrackID(const G4int& track) { fTrackID = track; }  // 设置径迹 ID
    inline void SetPosition(const G4ThreeVector& pos) { fPos = pos; }  // 设置命中位置
    inline void SetEdep(const G4double& edep) { fEdep = edep; }       // 设置能量沉积
    inline void SetPartDef(const G4ParticleDefinition* partDef)       // 设置粒子定义
    {
      fPartDef = partDef;
    }
    inline void SetEventParticleID(G4int eid) { fEventParticleID = eid; }  // 设置单事件粒子 ID

    // ===== 查询接口 (GET) =====
    inline G4int GetTrackID() const { return fTrackID; }                        // 取径迹 ID
    inline G4ThreeVector GetPosition() const { return fPos; }                   // 取命中位置
    inline G4double GetEdep() const { return fEdep; }                           // 取能量沉积
    inline const G4ParticleDefinition* GetPartDef() const { return fPartDef; }  // 取粒子定义
    inline G4int GetEventParticleID() const { return fEventParticleID; }        // 取单事件粒子 ID

  private:
    // ===== 命中数据成员 =====
    G4int fTrackID = 0;                             // 径迹 ID
    G4ThreeVector fPos = G4ThreeVector();           // 命中位置
    G4double fEdep = 0.;                            // 能量沉积（内部单位）
    const G4ParticleDefinition* fPartDef = nullptr;  // 粒子定义
    G4int fEventParticleID = -1;                    // 单事件粒子 ID(按粒子分组打分, 路线2)
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<TrackerHit> TrackerHitColl;  // TrackerHit 命中集合类型别名

extern G4ThreadLocal G4Allocator<TrackerHit>* TrackerHitAllocator;  // 线程局部分配器指针

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TrackerHit::operator new(size_t)
{
  // 懒初始化线程局部分配器，并从其中分配单个命中对象
  if (!TrackerHitAllocator) TrackerHitAllocator = new G4Allocator<TrackerHit>;
  return (void*)TrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TrackerHit::operator delete(void* hit)
{
  // 将命中对象归还给线程局部分配器
  TrackerHitAllocator->FreeSingle((TrackerHit*)hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
