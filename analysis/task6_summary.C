/// \file task6_summary.C
/// \brief 任务6 汇总：四区室每衰变核剂量 S(N←source)（分布效应主因）
///
/// 逐一打开核/胞质/膜/胞外四个源位置的 microtrack ROOT 文件，按
/// S(N←source) = 4 × 命中率 × <edep_n>/m_n 计算每衰变交给细胞核的剂量，
/// 并以核源为基准给出相对值，揭示源位置分布对核剂量的决定性影响。
///
/// 依赖文件：data/microtrack_{Nuc,Cyt,membrane,Ext}.root

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <cstdio>
/// 计算四区室每衰变核剂量 S(N←source) 并打印对比表。
/// z̄_{n,D}、命中率、<edep_n>/hit 已由 task6_compartments.C 给出，
/// 本函数补充各源位置相对核源的 S 值，量化分布效应。
void task6_summary(){
  // z̄_{n,D} 与命中率 + edep_n/hit 已由 task6_compartments.C 给出; 这里补 S(N←source)
  const char* files[]={"data/microtrack_Nuc.root","data/microtrack_Cyt.root",
                       "data/microtrack_membrane.root","data/microtrack_Ext.root"};
  const char* labels[]={"Nucleus","Cytoplasm","Membrane","Extracell"};
  double Rn=6.2e-6, m_n=(4./3.)*TMath::Pi()*Rn*Rn*Rn*1000.;  // 核半径(m), 核质量(kg, 水密度1000kg/m³)
  printf("\n===== 6.3 每衰变核剂量 S(N<-source) [= 4×命中率×<edep_n>/m_n] =====\n");
  printf("%-12s %8s %12s %14s %12s\n","区室","命中率","<edep_n>/hit","S(N<-src) Gy/dec","相对N");
  double SNuc=0;                                            // 核源 S 值(基准), 首个区室填入
  for(int i=0;i<4;i++){
    TFile *f=TFile::Open(files[i]); if(!f||f->IsZombie()){printf("%-12s 缺文件\n",labels[i]);continue;}
    TTree *t=(TTree*)f->Get("events");
    long N=t->GetEntries();                  // 该区室文件的总事件数
    TH1D h("h","",2000,0,5e6); t->Project("h","edep_n_keV","hitFlag==1");
    long nh=h.GetEntries();                  // 命中事件数 (hitFlag==1 投影后的条目数)
    double hitRate=(double)nh/N;             // 命中率 = 命中事件/总事件
    double edep=h.GetMean();                 // keV/hit (仅命中事件)
    double Ssrc = 4.*hitRate*edep*1.602e-16/m_n;  // Gy/decay (4α×命中×edep/m)
    if(i==0) SNuc=Ssrc;
    printf("%-12s %7.1f%% %10.1f keV %14.4f %11.2fx\n",
           labels[i], hitRate*100, edep, Ssrc, SNuc>0?Ssrc/SNuc:0);
    f->Close();
  }
  printf("\n解读: N源每衰变给核的剂量最高(100%%命中+满能量); 膜面源仅~19%%命中且α需穿胞质,\n");
  printf("      故每衰变核剂量远低于N源 → 每衰变杀伤 N >> Cy > Mem > Ext。这是分布效应主因。\n");
  printf("      每-Gy 存活率各区室相近(同为Ac-225 α品质), 差异主要在剂量交付(S值)。\n");
}
