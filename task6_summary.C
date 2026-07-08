// task6_summary.C  ——  任务6 汇总: 4区室 S(N←source) 每衰变核剂量 (=分布效应主因)
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <cstdio>
void task6_summary(){
  // z̄_{n,D} 与命中率 + edep_n/hit 已由 task6_compartments.C 给出; 这里补 S(N←source)
  const char* files[]={"microtrack_Nuc.root","microtrack_Cyt.root",
                       "microtrack_membrane.root","microtrack_Ext.root"};
  const char* labels[]={"Nucleus","Cytoplasm","Membrane","Extracell"};
  double Rn=6.2e-6, m_n=(4./3.)*TMath::Pi()*Rn*Rn*Rn*1000.;
  printf("\n===== 6.3 每衰变核剂量 S(N<-source) [= 4×命中率×<edep_n>/m_n] =====\n");
  printf("%-12s %8s %12s %14s %12s\n","区室","命中率","<edep_n>/hit","S(N<-src) Gy/dec","相对N");
  double SNuc=0;
  for(int i=0;i<4;i++){
    TFile *f=TFile::Open(files[i]); if(!f||f->IsZombie()){printf("%-12s 缺文件\n",labels[i]);continue;}
    TTree *t=(TTree*)f->Get("events");
    long N=t->GetEntries();
    TH1D h("h","",2000,0,5e6); t->Project("h","edep_n_keV","hitFlag==1");
    long nh=h.GetEntries();
    double hitRate=(double)nh/N;
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
