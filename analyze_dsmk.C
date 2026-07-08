// analyze_dsmk.C  ——  任务5增强: DSMK 模型 (Sato2012) + 对比 modified-SMK (Inaniwa2018)
//
// DSMK (double-stochastic MK): 严格数值积分 Sato2012 eq(7)(19), 饱和 z'_d=z0*sqrt(1-e^{-(zd/z0)^2}) (式16, 有sqrt)
//   - 单事件谱 f_{d,1}(zd), f_{n,1}(zn) 直接取自 ntuple
//   - 多事件分布用复合 Poisson 卷积 (式20-22, 26-28), MC 实现
// modified-SMK: Inaniwa2018 闭式 式(24), 饱和 z*=z0(1-e^{-(z/z0)^2}) (式1, 无sqrt)
//
// 默认参数 = Sato2012 DSMK HSG (Table1): α0=0.156 β0=0.0607 r_d=0.274 r_n=6.2 z0=89
// 用法: conda run -n microtrack root -b -q analyze_dsmk.C
//       conda run -n microtrack root -b -q 'analyze_dsmk.C("microtrack.root",0.156,0.0607,89.0)'

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
#include <vector>
#include <cstdio>
#include <TString.h>
#include <TSystem.h>
using namespace std;

// 线性插值查表(zn -> 值)
double interp(const vector<double>& xv, const vector<double>& yv, double x){
  int n=xv.size();
  if(x<=xv[0]) return yv[0];
  if(x>=xv[n-1]) return yv[n-1];
  int i=(int)((x-xv[0])/(xv[1]-xv[0]));
  if(i<0)i=0; if(i>=n-1)i=n-2;
  double dy=(yv[i+1]-yv[i])/(xv[i+1]-xv[i]);
  return yv[i]+dy*(x-xv[i]);
}

void analyze_dsmk(const char* fname = "microtrack.root",
                  double alpha0 = 0.156,   // Sato2012 DSMK HSG α0, Gy^-1
                  double beta0  = 0.0607,  // β0, Gy^-2
                  double z0     = 89.0,    // 饱和参数 z0, Gy
                  double Dmax   = 15.0)
{
  TFile *f=TFile::Open(fname);
  if(!f||f->IsZombie()){Error("analyze_dsmk","打不开 %s",fname);return;}
  TTree *t=dynamic_cast<TTree*>(f->Get("events"));
  if(!t){Error("analyze_dsmk","无 events ntuple");f->Close();return;}

  double zd=0,zn=0,w=0; int hf=0;
  t->SetBranchAddress("z_d_Gy",&zd); t->SetBranchAddress("z_n_Gy",&zn);
  t->SetBranchAddress("weight",&w);  t->SetBranchAddress("hitFlag",&hf);

  // ---- 1. 单事件谱(仅命中事件) ----
  TH1D *hd=new TH1D("hd","",3000,0,150);   // f_{d,1}(zd), 用重要性权 w 填充
  TH1D *hn=new TH1D("hn","",3000,0,5);     // f_{n,1}(zn), 无权
  double Sw=0,Swz=0,Swzz=0,Swzs=0,Swzszs=0,Szn=0,Sznzn=0; long nh=0;
  const long Nall=t->GetEntries();
  for(long i=0;i<Nall;i++){
    t->GetEntry(i); if(!hf) continue; nh++;
    hd->Fill(zd,w); hn->Fill(zn);
    // 同时累加微剂量学矩
    double zs_in = z0*(1-exp(-(zd/z0)*(zd/z0)));   // Inaniwa 式1 (无sqrt) 用于 mod-SMK
    Sw+=w; Swz+=w*zd; Swzz+=w*zd*zd; Swzs+=w*zs_in; Swzszs+=w*zs_in*zs_in;
    Szn+=zn; Sznzn+=zn*zn;
  }
  if(nh==0){Error("analyze_dsmk","无命中事件");f->Close();return;}
  double zdF=Swz/Sw, zdD=Swzz/Swz;            // z̄_{d,F}, z̄_{d,D}
  double zsF=Swzs/Sw, zsD=Swzszs/Swzs;        // z̄*_{d,F}, z̄*_{d,D} (Inaniwa 无sqrt)
  double znF=hn->GetMean(), znD=Sznzn/Szn;     // z̄_{n,F}, z̄_{n,D}

  printf("\n========= 微剂量学量 (HSG 一致几何) =========\n");
  printf("文件 %s: 总事件=%ld 命中=%ld (%.1f%%)\n",fname,Nall,nh,100.0*nh/Nall);
  printf("z̄_{d,F}=%.3f  z̄_{d,D}=%.3f Gy\n",zdF,zdD);
  printf("z̄*_{d,F}=%.3f z̄*_{d,D}=%.3f Gy (Inaniwa sat, z0=%.0f)\n",zsF,zsD,z0);
  printf("z̄_{n,F}=%.4g  z̄_{n,D}=%.4g Gy\n",znF,znD);
  printf("=============================================\n");

  TRandom3 rnd(12345);

  // ---- 2. DSMK: 预计算 <z'_d>(zn), <z'_d²>(zn) 在 zn 网格上 (Sato 式16 sqrt 饱和) ----
  int NZN=120; double znGridMax=40;
  vector<double> zgv(NZN+1), zsbar(NZN+1), zs2bar(NZN+1);
  int Msub=3000;
  for(int iz=0; iz<=NZN; iz++){
    double znb=iz*znGridMax/NZN; zgv[iz]=znb;
    double kmean=(zdF>1e-9)? znb/zdF : 0.0;
    double s1=0,s2=0;
    for(int m=0;m<Msub;m++){
      int k = rnd.Poisson(kmean);
      double zds=0; for(int ki=0;ki<k;ki++) zds+=hd->GetRandom(&rnd);
      double zpr = z0*sqrt(1-exp(-(zds/z0)*(zds/z0)));   // z'_d (Sato 式16)
      s1+=zpr; s2+=zpr*zpr;
    }
    zsbar[iz]=s1/Msub; zs2bar[iz]=s2/Msub;
  }

  // ---- 3. DSMK S(D): MC over 细胞试验 ----
  TGraph *gDSMK=new TGraph(); gDSMK->SetName("S_DSMK");
  int Ntrial=30000;
  double kDmax = Dmax/znF;
  for(double D=0.0; D<=Dmax+1e-9; D+=0.25){
    double Ssum=0;
    int Nav = (kDmax<50)? 0 : 1;  // 提示: k(D) 大时仍用 MC(GetRandom 快)
    for(int it=0; it<Ntrial; it++){
      int Nn = rnd.Poisson(D/znF);
      double zns=0; for(int n=0;n<Nn;n++) zns+=hn->GetRandom(&rnd);
      double zs_b = interp(zgv,zsbar,zns);
      double zs2_b= interp(zgv,zs2bar,zns);
      double Sn = exp(-alpha0*zs_b - beta0*zs2_b);   // Sato 式19
      Ssum += Sn;
    }
    gDSMK->SetPoint(gDSMK->GetN(), D, Ssum/Ntrial);
  }

  // ---- 4. modified-SMK 闭式 (Inaniwa 式24) + 经典 MK, 同参数对比 ----
  TGraph *gMSMK=new TGraph(); gMSMK->SetName("S_mSMK");
  TGraph *gMK  =new TGraph(); gMK->SetName("S_MK");
  double aS=alpha0+beta0*zsD, bS=beta0*zsD/zdD;       // Inaniwa 式15-16
  double aM=alpha0+beta0*zsD, bM=beta0;               // 经典 MK(z*-based)
  for(double D=0; D<=Dmax+1e-9; D+=0.25){
    double lq=exp(-aS*D-bS*D*D);
    double corr=1+(D*znD/2.0)*((aS+2*bS*D)*(aS+2*bS*D)-2*bS);
    if(corr<0)corr=0;
    gMSMK->SetPoint(gMSMK->GetN(),D,lq*corr);
    gMK->SetPoint(gMK->GetN(),D,exp(-aM*D-bM*D*D));
  }
  printf("\n模型系数(同 Sato DSMK 参数): α0=%.3f β0=%.4f z0=%.0f\n",alpha0,beta0,z0);
  printf("  DSMK        : 数值积分(式7,19)\n");
  printf("  mod-SMK     : α=%.3f β=%.4f (式15-16+24)\n",aS,bS);
  printf("  经典MK(z*)  : α=%.3f β=%.4f\n",aM,bM);

  printf("\n===== S(D) =====\n  %6s %10s %10s %10s\n","D","DSMK","modSMK","MK");
  for(int i=0;i<gDSMK->GetN();i++){
    double D,yd,ym,yk; gDSMK->GetPoint(i,D,yd);gMSMK->GetPoint(i,D,ym);gMK->GetPoint(i,D,yk);
    if(D<1e-6||(int)(D+1e-6)==(int)D) printf("  %6.2f %10.4g %10.4g %10.4g\n",D,yd,ym,yk);
  }

  // ---- 输出目录检查 ----
  TString outDir = "./result/DSMK/";
  if (gSystem->AccessPathName(outDir))
  {
      printf("[目录] %s 不存在, 创建中...\n", outDir.Data());
      gSystem->mkdir(outDir, kTRUE); // kTRUE 递归创建多级目录
  }

  TString outFile = outDir + "survival_dsmk.root";
  TString outPng  = outDir + "survival_dsmk.png";
  TString outPdf  = outDir + "survival_dsmk.pdf";
  // ---- 5. 画图 ----
  TCanvas *c=new TCanvas("cDSMK","DSMK survival",760,560);
  gPad->SetLogy(); gPad->SetGrid();
  gDSMK->SetLineColor(kBlack); gDSMK->SetLineWidth(2);
  gDSMK->GetXaxis()->SetTitle("Absorbed dose D [Gy]");
  gDSMK->GetYaxis()->SetTitle("S(D)");
  gDSMK->GetYaxis()->SetRangeUser(1e-6,1.5);
  gDSMK->SetTitle(Form("Ac-225 / HSG (r_{d}=0.274, r_{n}=6.2 #mum); z0=%.0f Gy",z0));
  gDSMK->Draw("AL");
  gMSMK->SetLineColor(kRed); gMSMK->SetLineWidth(2); gMSMK->Draw("L SAME");
  gMK->SetLineColor(kBlue); gMK->SetLineWidth(2); gMK->SetLineStyle(2); gMK->Draw("L SAME");
  TLegend *lg=new TLegend(0.60,0.75,0.92,0.92);
  lg->AddEntry(gDSMK,"DSMK (Sato2012)","l");
  lg->AddEntry(gMSMK,Form("mod-SMK #alpha=%.2f",aS),"l");
  lg->AddEntry(gMK,Form("MK #alpha=%.2f #beta=%.3f",aM,bM),"l");
  lg->SetFillStyle(0); lg->Draw();
  c->SaveAs(outPng); c->SaveAs(outPdf);
  printf("\n[图] %s / %s 已保存\n", outPng.Data(), outPdf.Data());

  TFile *fout=new TFile(outFile,"RECREATE");
  gDSMK->Write(); gMSMK->Write(); gMK->Write();
  fout->WriteObject(new TParameter<double>("z_d_F",zdF),"z_d_F");
  fout->WriteObject(new TParameter<double>("z_d_D",zdD),"z_d_D");
  fout->WriteObject(new TParameter<double>("z_n_F",znF),"z_n_F");
  fout->WriteObject(new TParameter<double>("z_n_D",znD),"z_n_D");
  fout->Close();
  printf("[结果] %s 已保存\n", outFile.Data());
  f->Close();
}
