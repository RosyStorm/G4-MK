/// \file task6_compartments.C
/// \brief 任务6.3：四区室 z̄_{n,D} 与 DSMK S(D) 对比
///
/// 读取核/胞质/膜/胞外四个源位置的 ROOT 文件，对每个区室计算单事件矩
/// 与 DSMK 存活曲线（显式卷积），输出 z̄_{n,D} 柱状图与 S(D) 四曲线。
///
/// 输出: z̄_{n,D} 柱状图 + S(D) 四曲线 → result/compartment/comparison.*

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TString.h>
#include <TSystem.h>
#include <vector>
#include <cstdio>
using namespace std;

/// 线性插值（等距节点）。越界时钳位到端点值。
/// \param xv 自变量数组（须等距、单调递增）
/// \param yv 函数值数组（与 xv 等长）
/// \param x  待插值的自变量
/// \return x 处的线性插值结果；x 超出范围时返回最近端点的 yv
double interp(const vector<double>& xv,const vector<double>& yv,double x){
  int n=xv.size(); if(x<=xv[0])return yv[0]; if(x>=xv[n-1])return yv[n-1];
  int i=(int)((x-xv[0])/(xv[1]-xv[0])); if(i<0)i=0; if(i>=n-1)i=n-2;
  return yv[i]+(yv[i+1]-yv[i])/(xv[i+1]-xv[i])*(x-xv[i]);
}

/// 对一棵 events 树用 DSMK（显式卷积）计算存活曲线 S(D)，同时返回单事件矩与命中数。
/// 流程：先填充域/核单事件谱并累加矩，再预计算 <z'_d>(zn) 网格，
/// 最后对每个剂量 D 做复合 Poisson 抽样得到 S(D)。
/// \param t     events 树（含 z_d_Gy, z_n_Gy, weight, hitFlag 分支）
/// \param a0    细胞系线性项系数 α0 (Gy^-1)
/// \param b0    细胞系二次项系数 β0 (Gy^-2)
/// \param z0    域饱和参数 z0 (Gy)
/// \param Dmax  存活曲线剂量上限 (Gy)
/// \param zdF   [out] 域频率平均比能 z̄_{d,F} (Gy)
/// \param znD   [out] 核剂量平均比能 z̄_{n,D} (Gy)
/// \param nh    [out] 命中事件数
/// \return DSMK 存活曲线 TGraph（横轴 D，纵轴 S）
TGraph* dsmkSD(TTree* t,double a0,double b0,double z0,double Dmax,
               double& zdF,double& znD,long& nh){
  double zd=0,zn=0,w=0; int hf=0;
  t->SetBranchAddress("z_d_Gy",&zd); t->SetBranchAddress("z_n_Gy",&zn);
  t->SetBranchAddress("weight",&w);  t->SetBranchAddress("hitFlag",&hf);
  TH1D fd1("fd1","",3000,0,150), fn1("fn1","",3000,0,5);
  double Sw=0,Swz=0,Swzz=0,Szn=0,Sznzn=0; nh=0; long Nall=t->GetEntries();
  for(long i=0;i<Nall;i++){ t->GetEntry(i); if(!hf)continue; nh++;
    fd1.Fill(zd,w); fn1.Fill(zn); Sw+=w; Swz+=w*zd; Swzz+=w*zd*zd; Szn+=zn; Sznzn+=zn*zn; }
  zdF=Swz/Sw; znD=Sznzn/Szn; double znF=fn1.GetMean();
  TRandom3 rnd(12345);
  // 预计算 <z'_d>(zn) 网格 (参数对齐验证版 analyze_dsmk.C)
  int NZN=160; double znMax=50; vector<double> zgv(NZN+1),zb(NZN+1),zb2(NZN+1);
  int Msub=6000;
  for(int iz=0;iz<=NZN;iz++){ double znb=iz*znMax/NZN; zgv[iz]=znb; double s1=0,s2=0;
    for(int m=0;m<Msub;m++){ int k=rnd.Poisson(znb/zdF); double zds=0;
      for(int ki=0;ki<k;ki++) zds+=fd1.GetRandom(&rnd);
      double zp=z0*sqrt(1-exp(-(zds/z0)*(zds/z0))); s1+=zp; s2+=zp*zp; }
    zb[iz]=s1/Msub; zb2[iz]=s2/Msub; }
  TGraph *g=new TGraph();
  int Ntrial=40000;
  for(double D=0;D<=Dmax+1e-9;D+=0.25){ double S=0;
    for(int it=0;it<Ntrial;it++){ int Nn=rnd.Poisson(D/znF); double zns=0;
      for(int n=0;n<Nn;n++) zns+=fn1.GetRandom(&rnd);
      S+=exp(-a0*interp(zgv,zb,zns)-b0*interp(zgv,zb2,zns)); }
    g->SetPoint(g->GetN(),D,S/Ntrial); }
  return g;
}

/// 主入口：对四个区室分别调用 dsmkSD 计算 z̄_{n,D} 与 S(D)，
/// 打印对比表并绘制 S(D) 四曲线与 z̄_{n,D} 柱状图到 result/compartment/。
void task6_compartments(){
  const char* files[]={"data/microtrack_Nuc.root","data/microtrack_Cyt.root",
                       "data/microtrack_membrane.root","data/microtrack_Ext.root"};
  const char* labels[]={"Nucleus","Cytoplasm","Membrane","Extracellular"};
  int cols[]={kBlack,kBlue,kRed,kGreen+2};           // 四区室曲线颜色
  double a0=0.156,b0=0.0607,z0=89,Dmax=12;            // DSMK 参数: α0, β0, z0, 剂量上限
  double znDv[4],zdFv[4],hitR[4]; long nhv[4]; double Ntot[4];  // 各区室: z̄_{n,D}, z̄_{d,F}, 命中率%, 命中数, 总事件
  vector<TGraph*> gs;
  printf("\n===== 6.3 四区室辐射品质对比 =====\n");
  printf("%-13s %8s %10s %10s %10s\n","区室","命中率%","z̄_{d,F}","z̄_{n,D}","D@S=0.1");
  for(int i=0;i<4;i++){
    TFile *f=TFile::Open(files[i]);
    if(!f||f->IsZombie()){ printf("  !! 缺 %s\n",files[i]); znDv[i]=0; gs.push_back(nullptr); continue; }
    TTree *t=(TTree*)f->Get("events"); Ntot[i]=t->GetEntries();
    TGraph *g=dsmkSD(t,a0,b0,z0,Dmax,zdFv[i],znDv[i],nhv[i]);
    hitR[i]=100.0*nhv[i]/Ntot[i]; gs.push_back(g);
    // 找 S=0.1 对应 D
    double D10=0; for(int k=0;k<g->GetN();k++){ double D,y; g->GetPoint(k,D,y); if(y<=0.1){D10=D;break;} }
    printf("%-13s %8.1f %10.2f %10.4g %10.2f\n",labels[i],hitR[i],zdFv[i],znDv[i],D10);
    f->Close();
  }

  // ===== 画 S(D) 四曲线 =====
  TString dir="./result/compartment/"; if(gSystem->AccessPathName(dir))gSystem->mkdir(dir,kTRUE);
  TCanvas *cS=new TCanvas("cS","compartment S(D)",760,560);
  gPad->SetLogy();
  gPad->SetGrid();
  TLegend *lg=new TLegend(0.65,0.75,0.92,0.92);
  for(int i=0;i<4;i++){ if(!gs[i])continue;
    gs[i]->SetLineColor(cols[i]);
    gs[i]->SetLineWidth(2);
    gs[i]->SetTitle("Ac-225/HSG compartment effect (DSMK);Absorbed dose D [Gy];S(D)");
    gs[i]->GetXaxis()->SetRangeUser(0,Dmax);
    gs[i]->GetYaxis()->SetRangeUser(1e-5,1.5);
    gs[i]->Draw(i?"L SAME":"AL");
    lg->AddEntry(gs[i],Form("%s (hit %.0f%%)",labels[i],hitR[i]),"l"); }
  lg->SetFillStyle(0);
  lg->Draw();
  cS->SaveAs(dir+"survival_compartment.png");
  cS->SaveAs(dir+"survival_compartment.pdf");

  // ===== z̄_{n,D} 柱状图 =====
  TCanvas *cZ=new TCanvas("cZ","znD",500,560);
  TH1D *hZ=new TH1D("hZ",";source compartment;z_{n,D} [Gy] (per event)",4,0,4);
  for(int i=0;i<4;i++){ hZ->SetBinContent(i+1,znDv[i]); hZ->GetXaxis()->SetBinLabel(i+1,labels[i]); }
  hZ->SetFillColor(kOrange-3);
  hZ->Draw("bar");
  hZ->SetMinimum(0);
  cZ->SaveAs(dir+"znD_compartment.png");
  printf("\n[输出] %s survival_compartment.png + znD_compartment.png\n",dir.Data());
}
