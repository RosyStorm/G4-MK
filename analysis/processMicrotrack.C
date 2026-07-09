/// \file processMicrotrack.C
/// \brief 为 microtrack.root 内的已知直方图补写 X/Y 轴标题
///
/// 后处理脚本：打开任务4产出的 ROOT 文件，根据直方图名称为其补写
/// X/Y 轴标题（线能 #varepsilon_{1}、比能 z_{1}、线能 y、f(y) 等微剂量学量），
/// 并以覆盖方式写回文件。
///
/// 用法（从 microtrack/ 根目录运行）：
///   conda run -n microtrack root -b -q analysis/processMicrotrack.C

#include <map>
#include <string>
#include <utility>
#include <vector>

#include <TError.h>
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>

/// 为指定 ROOT 文件内的已知直方图补写 X/Y 轴标题。
/// 以 UPDATE 模式打开文件，覆盖写回修改后的直方图。
/// \param filename ROOT 文件路径（默认 "data/microtrack.root"）
void processMicrotrack(const char* filename = "data/microtrack.root")
{
  // 以更新模式打开 ROOT 文件（允许写回修改后的直方图）
  TFile* f = new TFile(filename, "UPDATE");
  if (!f || f->IsZombie()) {
    Error("processMicrotrack", "Cannot open file %s", filename);
    return;
  }

  // 映射表：直方图名称 -> (X 轴标题, Y 轴标题)
  std::map<std::string, std::pair<const char*, const char*>> h1Axis = {
    {"fe", {"#varepsilon_{1} (keV)", "f(#varepsilon_{1})"}},
    {"efe", {"#varepsilon_{1} (keV)", "#varepsilon_{1} f(#varepsilon_{1})"}},
    {"e2fe", {"#varepsilon_{1} (keV)", "#varepsilon^{2}_{1} f(#varepsilon_{1})"}},
    {"fy", {"y (keV/#mum)", "f(y)"}},
    {"yfy", {"y (keV/#mum)", "y f(y)"}},
    {"y2fy", {"y (keV/#mum)", "y^{2} f(y)"}},
    {"fz", {"z_{1} (Gy)", "f(z_{1})"}},
    {"zfz", {"z_{1} (Gy)", "z_{1} f(z_{1})"}},
    {"z2fz", {"z_{1} (Gy)", "z^{2}_{1} f(z_{1})"}},
    {"Nsel", {"N_{sel}", "Counts"}},
    {"Nsite", {"N_{site}", "Counts"}},
    {"Nint", {"N_{int}", "Counts"}},
    {"KinE_in", {"T_{in} (MeV)", "Counts"}},
    {"KinE_out", {"T_{out} (MeV)", "Counts"}},
  };

  // ===== 1. 补写轴标题（仅对已存在的直方图） =====
  for (const auto& kv : h1Axis) {
    if (TH1* h = (TH1*)f->Get(kv.first.c_str())) {
      h->GetXaxis()->SetTitle(kv.second.first);
      h->GetYaxis()->SetTitle(kv.second.second);
      h->Write("", TObject::kOverwrite);  // 以覆盖方式写回，持久化本次更新
    }
  }

  f->Close();
}

