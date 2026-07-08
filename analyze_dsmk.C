// analyze_dsmk.C  ——  DSMK 模型 (Sato2012)，卷积显式实现 + 矩关系验证
//
// 严格按 Sato2012 式(20)(21)(26)(27) 用复合 Poisson 卷积构造多事件分布:
//   f_{d,k}(zd) = f_{d,1} * f_{d,k-1}            (式20, k 重卷积)
//   f_d(zd,zn) = Σ_k Poisson(k(zn),k) f_{d,k}(zd),  k(zn)=zn/z̄_{d,F}  (式21,22)
//   f_{n,k}(zn) = f_{n,1} * f_{n,k-1}            (式26)
//   f_n(zn,D) = Σ_k Poisson(k(D),k) f_{n,k}(zn),   k(D)=D/z̄_{n,F}     (式27,28)
// 卷积用 MC 实现(Sato Fig3 同算法: 抽 k 个单事件样本求和 = 一次多事件抽样)。
// 验证: <zn>=D, <zn²>=z̄_{n,D}D+D²  (复合 Poisson 矩关系, 证明卷积正确)
//
// DSMK S(D) = ∫ Sn(zn) f_n(zn,D) dzn  (式7)
// Sn(zn)   = exp(-α0<z'_d>(zn) - β0<z'_d²>(zn))  (式19), z'_d=z0√(1-e^{-(zd/z0)²}) (式16)
//
// 默认参数 = Sato2012 DSMK HSG (Table1): α0=0.156 β0=0.0607 r_d=0.274 r_n=6.2 z0=89

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TMath.h>
#include <TParameter.h>
#include <TRandom3.h>
#include <TString.h>
#include <TSystem.h>
#include <vector>
#include <cstdio>
using namespace std;

double interp(const vector<double>& xv, const vector<double>& yv, double x){
  int n=xv.size();
  if(x<=xv[0]) return yv[0]; if(x>=xv[n-1]) return yv[n-1];
  int i=(int)((x-xv[0])/(xv[1]-xv[0])); if(i<0)i=0; if(i>=n-1)i=n-2;
  return yv[i]+(yv[i+1]-yv[i])/(xv[i+1]-xv[i])*(x-xv[i]);
}

void analyze_dsmk(const char* fname = "microtrack.root",
                  double alpha0 = 0.156, double beta0 = 0.0607, double z0 = 89.0,
                  double Dmax = 15.0)
{
  TFile *f=TFile::Open(fname);
  if(!f||f->IsZombie()){Error("analyze_dsmk","打不开 %s",fname);return;}
  TTree *t=dynamic_cast<TTree*>(f->Get("events"));
  if(!t){Error("analyze_dsmk","无 events ntuple");f->Close();return;}
  double zd=0,zn=0,w=0; int hf=0;
  t->SetBranchAddress("z_d_Gy",&zd); t->SetBranchAddress("z_n_Gy",&zn);
  t->SetBranchAddress("weight",&w);  t->SetBranchAddress("hitFlag",&hf);

  // ============ 1. 单事件谱 f_{d,1}, f_{n,1} (per-event) ============
  TH1D *fd1=new TH1D("fd1",";z_{d} [Gy];f_{d,1}",3000,0,150);
  TH1D *fn1=new TH1D("fn1",";z_{n} [Gy];f_{n,1}",3000,0,5);
  double Sw=0,Swz=0,Swzz=0,Szn=0,Sznzn=0; long nh=0; const long Nall=t->GetEntries();
  for(long i=0;i<Nall;i++){ t->GetEntry(i); if(!hf)continue; nh++;
    fd1->Fill(zd,w); fn1->Fill(zn);
    Sw+=w; Swz+=w*zd; Swzz+=w*zd*zd; Szn+=zn; Sznzn+=zn*zn; }
  if(nh==0){Error("analyze_dsmk","无命中事件");f->Close();return;}
  double zdF=Swz/Sw, zdD=Swzz/Swz;   // 单事件 z̄_{d,F}, z̄_{d,D}  (式5)
  double znF=fn1->GetMean(), znD=Sznzn/Szn;  // 单事件 z̄_{n,F}, z̄_{n,D}
  printf("\n===== 单事件 (per-event) 谱与矩 [式5/13, 喂卷积参数 k=zn/zdF] =====\n");
  printf("命中=%ld/%ld (%.1f%%)  z̄_{d,F}=%.3f z̄_{d,D}=%.3f  z̄_{n,F}=%.4g z̄_{n,D}=%.4g Gy\n",
         nh,Nall,100.0*nh/Nall,zdF,zdD,znF,znD);

  TRandom3 rnd(12345);
  // 卷积核心: 从 f1 抽 k~Poisson(kmean) 个样本求和 → 一次多事件抽样 (Sato Fig3)
  auto convolveStep=[&](TH1* f1,double kmean,TH1* outH,int Ns){
    for(int s=0;s<Ns;s++){ int k=rnd.Poisson(kmean); double z=0;
      for(int i=0;i<k;i++) z+=f1->GetRandom(&rnd); outH->Fill(z); } };

  // ============ 2. 显式构造 f_d(zd, zn*) 几个代表 zn, 验证域卷积 ============
  printf("\n===== 域多事件谱 f_d(zd,zn) 卷积验证 [式20-22] =====\n");
  printf("  %8s %10s %10s | %10s %12s | %8s\n","zn","<zd>","(=zn?)","<zd²>","znD·zn+zn²","匹配?");
  TFile *fout=new TFile("result/DSMK/convolved_distributions.root","RECREATE");
  // (假设 result/DSMK 已由脚本外创建; 若无则用当前目录)
  for(double znb : {1.0, 5.0, 10.0}){
    TH1D *hdz=new TH1D(Form("fd_zn%g",znb),Form("f_d(z_d,z_n=%g Gy);z_d [Gy]",znb),2000,0,200);
    convolveStep(fd1, znb/zdF, hdz, 200000);
    double m1=hdz->GetMean(), m2=hdz->GetMean()*hdz->GetMean()+hdz->GetRMS()*hdz->GetRMS();
    double pred2=zdD*znb+znb*znb;  // 复合Poisson: <zd²>=z̄_{d,D}·zn+zn²
    printf("  %8.2f %10.3f %10s | %10.3f %12.3f | %s\n",
           znb,m1,(fabs(m1-znb)<0.05*znb+0.01?"OK":"??"),m2,pred2,
           (fabs(m2-pred2)<0.1*pred2?"OK":"??"));
    hdz->Write();
  }

  // ============ 3. 预计算 <z'_d>(zn), <z'_d²>(zn) 在 zn 网格 (用于 DSMK 式19) ============
  // 注: 必须用精确 zds 直接累加 z'_d, 不能用直方图 bin 中心(否则 zd=0 被归到 bin 中心
  //     产生虚假 z'_d, 使 S(0)<1)。zds=0 → z'_d=0 严格成立。
  int NZN=160; double znGridMax=50;
  vector<double> zgv(NZN+1), zsbar(NZN+1), zs2bar(NZN+1);
  int Msub=6000;
  for(int iz=0; iz<=NZN; iz++){
    double znb=iz*znGridMax/NZN; zgv[iz]=znb;
    double s1=0,s2=0;
    for(int m=0;m<Msub;m++){
      int k=rnd.Poisson(znb/zdF); double zds=0;
      for(int ki=0;ki<k;ki++) zds+=fd1->GetRandom(&rnd);
      double zpr=z0*sqrt(1-exp(-(zds/z0)*(zds/z0)));   // z'_d 式16 (sqrt 饱和), 精确 zds
      s1+=zpr; s2+=zpr*zpr;
    }
    zsbar[iz]=s1/Msub; zs2bar[iz]=s2/Msub;
  }

  // ============ 4. 显式构造 f_n(zn, D*) 几个代表 D, 验证核卷积 + 给 DSMK 用 ============
  printf("\n===== 核多事件谱 f_n(zn,D) 卷积验证 [式26-28] =====\n");
  printf("  %8s %10s %10s | %12s %12s | %8s\n","D","<zn>","(=D?)","<zn²>","znD·D+D²","匹配?");
  for(double Dv : {1.0, 2.0, 5.0, 10.0}){
    TH1D *hnD=new TH1D(Form("fn_D%g",Dv),Form("f_n(z_n,D=%g Gy);z_n [Gy]",Dv),4000,0,60);
    convolveStep(fn1, Dv/znF, hnD, 300000);
    double m1=hnD->GetMean(), m2=hnD->GetMean()*hnD->GetMean()+hnD->GetRMS()*hnD->GetRMS();
    double pred2=znD*Dv+Dv*Dv;
    printf("  %8.2f %10.3f %10s | %12.4f %12.4f | %s\n",
           Dv,m1,(fabs(m1-Dv)<0.03*Dv+0.01?"OK":"??"),m2,pred2,
           (fabs(m2-pred2)<0.1*pred2?"OK":"??"));
    hnD->Write();
  }

  // ============ 5. DSMK S(D): 显式卷积 f_n(zn,D) → S(D)=∫Sn(zn) f_n dzn (式7) ============
  TGraph *gDSMK=new TGraph(); gDSMK->SetName("S_DSMK");
  int Ntrial=40000;
  for(double D=0.0; D<=Dmax+1e-9; D+=0.25){
    double Ssum=0;
    for(int it=0; it<Ntrial; it++){
      int Nn=rnd.Poisson(D/znF);            // 抽核事件数
      double zns=0; for(int n=0;n<Nn;n++) zns+=fn1->GetRandom(&rnd);  // 多事件 zn
      double Sn=exp(-alpha0*interp(zgv,zsbar,zns)-beta0*interp(zgv,zs2bar,zns)); // 式19
      Ssum+=Sn;
    }
    gDSMK->SetPoint(gDSMK->GetN(), D, Ssum/Ntrial);
  }

  // ============ 6. mod-SMK 闭式 + 经典 MK (对比; 它们按设计不卷积) ============
  double zsD_in=z0*(1-exp(-(zdD/z0)*(zdD/z0))); // 仅示意, mod-SMK 需积分形式
  // mod-SMK 用 Inaniwa 式24 (单事件 z̄*_{d,D}), 此处重新积分得 zsD:
  double sZs=0,sZs2=0,sZw=0; // Σw zs, Σw zs², Σw
  for(int b=1;b<=fd1->GetNbinsX();b++){ double zv=fd1->GetBinCenter(b),c=fd1->GetBinContent(b);
    if(c<=0)continue; double zs=z0*(1-exp(-(zv/z0)*(zv/z0))); sZs+=c*zs; sZs2+=c*zs*zs; sZw+=c; }
  double zsD = (sZs>0)? sZs2/sZs : 0;  // z̄*_{d,D} (Inaniwa 无sqrt, 式12)
  double aS=alpha0+beta0*zsD, bS=beta0*zsD/zdD, aM=alpha0+beta0*zsD, bM=beta0;
  TGraph *gMSMK=new TGraph(); gMSMK->SetName("S_mSMK");
  TGraph *gMK  =new TGraph(); gMK->SetName("S_MK");
  for(double D=0; D<=Dmax+1e-9; D+=0.25){
    double corr=1+(D*znD/2.0)*((aS+2*bS*D)*(aS+2*bS*D)-2*bS); if(corr<0)corr=0;
    gMSMK->SetPoint(gMSMK->GetN(),D,exp(-aS*D-bS*D*D)*corr);
    gMK->SetPoint(gMK->GetN(),D,exp(-aM*D-bM*D*D));
  }
  printf("\n模型系数: DSMK(数值卷积) | mod-SMK α=%.3f β=%.4f | MK α=%.3f β=%.4f\n",aS,bS,aM,bM);

  printf("\n===== S(D) =====\n  %6s %10s %10s %10s\n","D","DSMK","modSMK","MK");
  for(int i=0;i<gDSMK->GetN();i++){ double D,yd,ym,yk;
    gDSMK->GetPoint(i,D,yd);gMSMK->GetPoint(i,D,ym);gMK->GetPoint(i,D,yk);
    if(D<1e-6||(int)(D+1e-6)==(int)D) printf("  %6.2f %10.4g %10.4g %10.4g\n",D,yd,ym,yk); }

  // ============ 7. 画图 + 存 ============
  TString outDir="./result/DSMK/";
  if(gSystem->AccessPathName(outDir)) gSystem->mkdir(outDir,kTRUE);
  TCanvas *c=new TCanvas("cDSMK","DSMK",760,560); gPad->SetLogy(); gPad->SetGrid();
  gDSMK->SetLineColor(kBlack); gDSMK->SetLineWidth(2);
  gDSMK->GetXaxis()->SetTitle("Absorbed dose D [Gy]"); gDSMK->GetYaxis()->SetTitle("S(D)");
  gDSMK->GetYaxis()->SetRangeUser(1e-6,1.5);
  gDSMK->SetTitle(Form("Ac-225/HSG DSMK (explicit convolution); z0=%.0f Gy",z0));
  gDSMK->Draw("AL");
  gMSMK->SetLineColor(kRed); gMSMK->SetLineWidth(2); gMSMK->Draw("L SAME");
  gMK->SetLineColor(kBlue); gMK->SetLineWidth(2); gMK->SetLineStyle(2); gMK->Draw("L SAME");
  TLegend *lg=new TLegend(0.60,0.75,0.92,0.92);
  lg->AddEntry(gDSMK,"DSMK (Sato, explicit conv.)","l");
  lg->AddEntry(gMSMK,Form("mod-SMK #alpha=%.2f",aS),"l");
  lg->AddEntry(gMK,Form("MK #alpha=%.2f #beta=%.3f",aM,bM),"l");
  lg->SetFillStyle(0); lg->Draw();
  c->SaveAs(outDir+"survival_dsmk.png"); c->SaveAs(outDir+"survival_dsmk.pdf");

  gDSMK->Write(); gMSMK->Write(); gMK->Write();
  fout->WriteObject(new TParameter<double>("z_d_F",zdF),"z_d_F");
  fout->WriteObject(new TParameter<double>("z_d_D",zdD),"z_d_D");
  fout->WriteObject(new TParameter<double>("z_n_F",znF),"z_n_F");
  fout->WriteObject(new TParameter<double>("z_n_D",znD),"z_n_D");
  fout->Close();
  printf("\n[输出] %s 内: S(D) 图 + convolved_distributions.root(f_d/f_n 显式直方图)\n",outDir.Data());
  f->Close();
}
