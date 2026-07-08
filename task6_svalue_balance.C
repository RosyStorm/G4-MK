// task6_svalue_balance.C  ——  任务 6.1 (S值) + 6.2 (能量平衡)
// 6.1: S(N←N) from microtrack_Nuc.root (源在核内), 对比解析 + 文献量级
// 6.2: 核能量平衡 from microtrack_membrane.root: <KinE_in>-<KinE_out> vs <edep_n>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <cstdio>
void task6_svalue_balance(){
  const double Rn_um=6.2;                 // 核半径 um (HSG, Sato Table1)
  double Rn=Rn_um*1e-6;                   // m
  double m_n=(4./3.)*TMath::Pi()*Rn*Rn*Rn*1000.;  // kg (水 1000 kg/m³)
  printf("核质量 m_n (R_n=%.1f um) = %.4g kg\n\n", Rn_um, m_n);

  // ===== 6.1 S(N←N) =====
  TFile *fN=TFile::Open("microtrack_Nuc.root");
  if(fN && !fN->IsZombie()){
    TTree *tN=(TTree*)fN->Get("events");
    TH1D hN("hN","",2000,0,5e6); tN->Project("hN","edep_n_keV");
    double edep_keV=hN.GetMean();          // 每个α事件的核内沉积均值 (N源: 每事件都命中)
    double Salpha = edep_keV*1.602e-16/m_n;          // Gy per α
    double Sdecay = 4.*Salpha;                        // Gy per decay (Ac-225 释放4α)
    long nall=tN->GetEntries();
    printf("===== 6.1 S(N←N) [源均匀在核内] =====\n");
    printf("  事件=%ld  <edep_n>=%.2f keV/α\n", nall, edep_keV);
    printf("  S(N←N) = %.4f Gy/α  →  %.4f Gy/decay (= %.1f mGy/decay, ×4α)\n",
           Salpha, Sdecay, Sdecay*1e3);
    // 解析交叉校验: 均匀各向同性源在球内的平均弦长 = 3R/4 (几何)
    //   平均LET ≈ Eα(均)/range; Ac-225 α 均~6.9MeV, range~60um → ~115 keV/um
    double meanPath=0.75*Rn_um;            // 3R/4 um
    double LETapprox=115.;                  // keV/um (6.9MeV α 均LET, ASTAR)
    double edep_analytic_keV = LETapprox*meanPath;
    double Sanalytic = 4.*edep_analytic_keV*1.602e-16/m_n;
    printf("  解析估计(LET×3Rn/4×4α): %.4f Gy/decay  (G4/解析=%.2f)\n",
           Sanalytic, Sdecay/Sanalytic);
    printf("  文献量级(MIRDcell/Goddu, R_n~6um Ac-225): ~0.3-0.5 Gy/decay\n");
    fN->Close();
  } else printf("!! 无 microtrack_Nuc.root\n");

  // ===== 6.2 核能量平衡 =====
  printf("\n===== 6.2 核能量平衡 [膜面源, H1(12/13) vs edep_n] =====\n");
  TFile *fM=TFile::Open("microtrack_membrane.root");
  if(fM && !fM->IsZombie()){
    TH1 *hKin=(TH1*)fM->Get("KinE_in"), *hKout=(TH1*)fM->Get("KinE_out");
    TTree *tM=(TTree*)fM->Get("events");
    TH1D hE("hE","",2000,0,5e6); tM->Project("hE","edep_n_keV","hitFlag==1");
    double Kin=hKin->GetMean(), Kout=hKout->GetMean();   // MeV
    double edep_MeV=hE.GetMean()/1e3;
    printf("  <KinE_in>  = %.4f MeV  (初级α进核)\n", Kin);
    printf("  <KinE_out> = %.4f MeV  (初级α出核)\n", Kout);
    printf("  <KinE_in>-<KinE_out> = %.4f MeV  (初级在核内的总能损)\n", Kin-Kout);
    printf("  <edep_n>               = %.4f MeV  (核内总沉积, 含δ电子)\n", edep_MeV);
    printf("  差值 = %.4f MeV (= 在核内产生但逃出核的δ电子动能, 高LET下预期非零)\n",
           (Kin-Kout)-edep_MeV);
    printf("  比值 (KinE损/edep_n) = %.2f  (应略>1, δ逃逸所致)\n", (Kin-Kout)/edep_MeV);
    fM->Close();
  } else printf("!! 无 microtrack_membrane.root\n");
}
