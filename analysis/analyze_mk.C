// analyze_mk.C  ——  任务5: ROOT 后处理
// 读 events ntuple → 算微剂量学量 → MK / 修正 SMK 存活曲线 S(D)
// 参考: Inaniwa & Kanematsu, Phys. Med. Biol. 63 (2018) 095011
//
// 用法(从 microtrack/ 根目录运行):
//   conda run -n microtrack root -b -q analysis/analyze_mk.C
//   conda run -n microtrack root -b -q 'analysis/analyze_mk.C("data/microtrack.root",0.174,0.0568,66.0,15)'
//
// 输入列(任务4.2 ntuple): z_d_Gy, z_n_Gy, weight, hitFlag, alphaE_MeV, ...
// 输出: 终端打印各量 + S(D) 表 + result/mod-SMK/survival.*(图与结果)

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TMath.h>
#include <TParameter.h>
#include <TString.h>
#include <cstdio>
#include <TSystem.h>

/// 任务5 后处理：读 events ntuple 计算微剂量学量，给出 MK 与修正 SMK 存活曲线 S(D)。
/// 模型参考 Inaniwa & Kanematsu, Phys. Med. Biol. 63 (2018) 095011。
/// \param fname  输入 ROOT 文件路径（默认 "data/microtrack.root"）
/// \param alpha0 细胞系线性项系数 α0 (Gy^-1)，默认 0.174 (HSG, Inaniwa2018 Table1)
/// \param beta0  细胞系二次项系数 β0 (Gy^-2)，默认 0.0568 (HSG)
/// \param z0     域饱和参数 z0 (Gy)，默认 66.0 (HSG)
/// \param Dmax   存活曲线剂量上限 (Gy)，默认 5.0
void analyze_mk(const char* fname = "data/microtrack.root",
                double alpha0 = 0.174,   // 细胞系 α0 (HSG), Gy^-1  —— Inaniwa2018 Table1
                double beta0  = 0.0568,  // 细胞系 β0 (HSG), Gy^-2
                double z0     = 66.0,    // 饱和参数 z0 (HSG), Gy
                double Dmax   = 5.0)    // 存活曲线剂量上限, Gy
{
  // ===== 1. 读 ntuple =====
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) { Error("analyze_mk", "打不开 %s", fname); return; }
  TTree *t = dynamic_cast<TTree*>(f->Get("events"));
  if (!t) { Error("analyze_mk", "%s 内无 events ntuple", fname); f->Close(); return; }

  double zd=0, zn=0, w=0, ae=0; int hf=0;
  t->SetBranchAddress("z_d_Gy",   &zd);
  t->SetBranchAddress("z_n_Gy",   &zn);
  t->SetBranchAddress("weight",   &w);
  t->SetBranchAddress("hitFlag",  &hf);
  t->SetBranchAddress("alphaE_MeV",&ae);

  // ===== 2. 累加(仅命中事件 hitFlag==1) =====
  // 域级用重要性权 w(校正位点随机放置偏置)；核级无权(w=1)
  double Sw=0, Swz=0, Swzz=0, Swzs=0, Swzszs=0;  // w, w·z, w·z², w·z*, w·z*²
  double Szn=0, Sznzn=0;                          // z_n, z_n²
  long nh=0;
  const long Nall = t->GetEntries();
  for (long i=0;i<Nall;i++){
    t->GetEntry(i);
    if(!hf) continue;                       // miss 事件不计入单事件谱
    nh++;
    double zs = z0*(1.0 - TMath::Exp(-(zd/z0)*(zd/z0)));   // z*_d, 式(1)
    Sw += w; Swz += w*zd; Swzz += w*zd*zd;
    Swzs += w*zs; Swzszs += w*zs*zs;
    Szn += zn; Sznzn += zn*zn;
  }
  if (nh==0){ Error("analyze_mk","无命中事件, 无法计算"); f->Close(); return; }
  if (Swz<=0 || Swzs<=0 || Szn<=0){ Error("analyze_mk","求和异常(分母<=0)"); f->Close(); return; }

  // ===== 3. 微剂量学量 =====
  double zF  = Swz / Sw;          // z̄_{d,F}  频率均
  double zD  = Swzz / Swz;        // z̄_{d,D}  剂量均 (式11: ∫z²f / ∫zf)
  double zsF = Swzs / Sw;         // z̄*_{d,F} (式10)
  // P0 修复 #5: Inaniwa 式 (12) 分母应为 z̄_{d,F} (z_d 的频率均), 不是 z̄*_{d,F} (z*_d 的频率均)
  // 原代码 Swzszs/Swzs = ∫(z*)²f / ∫z*f 是 z* 的"自身剂量均", 与 Inaniwa 公式不符.
  // 正确公式: z̄*_{d,D} = ∫(z*)²f / z̄_{d,F}, 见 Inaniwa 2018 式(12) 与 Sato 2012 式(13).
  double zsD = (Swz > 0) ? Swzszs / Swz : 0;  // z̄*_{d,D} (Inaniwa 式12)
  double znF = Szn / nh;          // z̄_{n,F}
  double znD = Sznzn / Szn;       // z̄_{n,D}

  printf("\n================ 微剂量学量 (任务5) ================\n");
  printf("文件: %s   总事件=%ld  命中事件=%ld  命中率=%.1f%%\n",
         fname, Nall, nh, 100.0*nh/Nall);
  printf("细胞系参数: α0=%.3f β0=%.4f z0=%.1f (HSG)\n", alpha0, beta0, z0);
  printf("----------------------------------------------------\n");
  printf(" z̄_{d,F} = %.4f Gy      z̄_{d,D} = %.4f Gy\n", zF, zD);
  printf(" z̄*_{d,F}= %.4f Gy      z̄*_{d,D}= %.4f Gy   [饱和]\n", zsF, zsD);
  printf(" z̄_{n,F} = %.5g Gy      z̄_{n,D} = %.5g Gy\n", znF, znD);
  printf("----------------------------------------------------\n");

  // ===== 4. MK / SMK 系数 =====
  // 修正 SMK: 式(15)(16)
  double aS = alpha0 + beta0*zsD;
  double bS = beta0 * zsD / zD;
  // 经典 MK: α=α0+β0·z̄*_{d,D}, β=β0(常数)
  double aM = alpha0 + beta0*zsD;
  double bM = beta0;
  printf(" α_SMK=%.4f  β_SMK=%.4f   (修正SMK, 式15-16)\n", aS, bS);
  printf(" α_MK =%.4f  β_MK =%.4f   (经典MK)\n", aM, bM);
  printf("====================================================\n\n");

  // ===== 5. 存活曲线 S(D) =====
  // SMK 式(24): S = exp(-α_S D - β_S D²) · [1 + (D·z̄_{n,D}/2)·((α_S+2β_S D)² - 2β_S)]
  //   方括号项为 z_n 随机性修正(高LET下使 S 高于 MK)。低 LET 时 z̄_{n,D}→小, 项→1, 退化为 LQ。
  TGraph *gSMK = new TGraph(); gSMK->SetName("S_SMK");
  TGraph *gMK  = new TGraph(); gMK ->SetName("S_MK");
  double Dstep = 0.01;
  for (double D=0.0; D<=Dmax+1e-9; D+=Dstep){
    double lq = TMath::Exp(-aS*D - bS*D*D);
    double bracket = 1.0 + (D*znD/2.0)*((aS+2.0*bS*D)*(aS+2.0*bS*D) - 2.0*bS);
    if (bracket < 0.0) bracket = 0.0;          // 防御
    double Ssmk = lq*bracket;
    gSMK->SetPoint(gSMK->GetN(), D, Ssmk);
    gMK ->SetPoint(gMK->GetN(),  D, TMath::Exp(-aM*D - bM*D*D));
  }

  // S(D) 表(每 1 Gy 一行)
  printf("===== 存活曲线 S(D) =====\n");
  printf(" %6s  %10s  %10s\n","D(Gy)","S_MK","S_SMK");
  for (int i=0;i<gSMK->GetN();i++){
    double D,yM,yS; gMK->GetPoint(i,D,yM); gSMK->GetPoint(i,D,yS);
    if (D < 1e-6 || (int)(D+1e-6) == (int)D)
      printf(" %6.2f  %10.4g  %10.4g\n", D, yM, yS);
  }
  printf("\n");


  // ===== 存结果（路径准备） =====
  TString outDir = "./result/mod-SMK/";
  TString outFile = outDir + "survival.root";
  TString outPng  = outDir + "survival.png";
  TString outPdf  = outDir + "survival.pdf";

  if (gSystem->AccessPathName(outDir))
  {
      printf("[目录] %s 不存在, 创建中...\n", outDir.Data());
      gSystem->mkdir(outDir, kTRUE); // kTRUE 递归创建多级目录
  }
  // ===== 6. 画图 =====
  TCanvas *c = new TCanvas("cSMK","SMK survival",720,520);
  gPad->SetLogy();
  gPad->SetGrid();
  gSMK->SetLineColor(kRed);
  gSMK->SetLineWidth(2);
  gSMK->GetXaxis()->SetTitle("Absorbed dose D [Gy]");
  gSMK->GetYaxis()->SetTitle("Cell survival S(D)");
  gSMK->GetYaxis()->SetRangeUser(1e-6, 1.5);
  gSMK->SetTitle(Form("Ac-225 survival; #alpha_{SMK}=%.3f #beta_{SMK}=%.3f;z0=%.0f Gy",aS,bS,z0));
  gSMK->Draw("AL");
  gMK->SetLineColor(kBlue);
  gMK->SetLineWidth(2);
  gMK->SetLineStyle(2);
  gMK->Draw("L SAME");
  TLegend *lg = new TLegend(0.65,0.78,0.92,0.92);
  lg->AddEntry(gSMK,Form("SMK: #alpha=%.3f #beta=%.3f",aS,bS),"l");
  lg->AddEntry(gMK ,Form("MK:  #alpha=%.3f #beta=%.3f",aM,bM),"l");
  lg->SetFillStyle(0);
  lg->Draw();
  c->SaveAs(outPng);
  c->SaveAs(outPdf);
  printf("[图] %s / %s 已保存\n", outPng.Data(), outPdf.Data());

  // ===== 7. 存结果 =====
  TFile *fout = new TFile(outFile,"RECREATE");
  gSMK->Write();
  gMK->Write();
  fout->WriteObject(new TParameter<double>("z_d_F",zF),"z_d_F");
  fout->WriteObject(new TParameter<double>("z_d_D",zD),"z_d_D");
  fout->WriteObject(new TParameter<double>("zs_d_F",zsF),"zs_d_F");
  fout->WriteObject(new TParameter<double>("zs_d_D",zsD),"zs_d_D");
  fout->WriteObject(new TParameter<double>("z_n_F",znF),"z_n_F");
  fout->WriteObject(new TParameter<double>("z_n_D",znD),"z_n_D");
  fout->WriteObject(new TParameter<double>("alpha_SMK",aS),"alpha_SMK");
  fout->WriteObject(new TParameter<double>("beta_SMK",bS),"beta_SMK");
  fout->Close();
  printf("[结果] %s 已保存(图 + 各 TParameter)\n", outFile.Data());

  f->Close();
}
