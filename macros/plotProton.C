// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

// plotProton.C: figures from triveni-proton output.
//
//   root -l -b -q 'macros/plotProton.C("output/proton", "output/figures")'
//
// proton_maps:  thickness T_p(x, y) of the same quark configurations built as
//               hot spots, Y junction and Delta (rows = configurations).
// proton_stats: string lengths, junction position, single-proton eps2, mean profile.

#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace
{

const char* kShapeName[3] = {"hotspot", "Y", "Delta"};
const char* kShapeLabel[3] = {"hot spots", "Y junction", "#Delta"};

int shapeColor(int k)
{
  static const int colors[3] = {TColor::GetColor("#5790fc"), TColor::GetColor("#e42536"), TColor::GetColor("#f89c20")};
  return colors[k];
}

void save(TCanvas& c, const char* outDir, const char* name)
{
  c.SaveAs(Form("%s/%s.pdf", outDir, name));
  c.SaveAs(Form("%s/%s.png", outDir, name));
}

} // namespace

void plotProton(const char* stem = "output/proton", const char* outDir = "output/figures")
{
  gSystem->mkdir(outDir, true);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(64);

  // ---------------------------------------------------------------- maps
  TTree maps("maps", "maps");
  maps.ReadFile(Form("%s_maps.txt", stem));
  TTree quarks("quarks", "quarks");
  quarks.ReadFile(Form("%s_quarks.txt", stem));

  const int nConfig = static_cast<int>(maps.GetMaximum("config")) + 1;
  const int n = static_cast<int>(std::lround(std::sqrt(static_cast<double>(maps.GetEntries()) / nConfig)));
  const double xmin = maps.GetMinimum("x");
  const double xmax = maps.GetMaximum("x");
  const double step = (xmax - xmin) / (n - 1);
  const double lo = xmin - 0.5 * step;
  const double hi = xmax + 0.5 * step;

  std::vector<TH2D*> hist(3 * nConfig);
  for (int c = 0; c < nConfig; ++c) {
    for (int k = 0; k < 3; ++k) {
      hist[3 * c + k] = new TH2D(Form("map_%d_%d", c, k), Form("%s;x (fm);y (fm)", kShapeLabel[k]), n, lo, hi, n, lo, hi);
    }
  }
  Int_t config = 0;
  Double_t x = 0, y = 0, t[3] = {0, 0, 0};
  maps.SetBranchAddress("config", &config);
  maps.SetBranchAddress("x", &x);
  maps.SetBranchAddress("y", &y);
  for (int k = 0; k < 3; ++k) {
    maps.SetBranchAddress(kShapeName[k], &t[k]);
  }
  for (Long64_t i = 0; i < maps.GetEntries(); ++i) {
    maps.GetEntry(i);
    for (int k = 0; k < 3; ++k) {
      hist[3 * config + k]->Fill(x, y, t[k]);
    }
  }

  // Quark (open circles) and junction (star) positions per configuration and shape.
  std::vector<TGraph*> quarkGraph(3 * nConfig);
  std::vector<TGraph*> junctionGraph(3 * nConfig);
  for (auto& g : quarkGraph) {
    g = new TGraph();
  }
  for (auto& g : junctionGraph) {
    g = new TGraph();
  }
  Int_t qConfig = 0, qShape = 0, item = 0;
  Double_t qx = 0, qy = 0;
  quarks.SetBranchAddress("config", &qConfig);
  quarks.SetBranchAddress("shape", &qShape);
  quarks.SetBranchAddress("item", &item);
  quarks.SetBranchAddress("x", &qx);
  quarks.SetBranchAddress("y", &qy);
  for (Long64_t i = 0; i < quarks.GetEntries(); ++i) {
    quarks.GetEntry(i);
    TGraph* g = (item < 3 ? quarkGraph : junctionGraph)[3 * qConfig + qShape];
    g->SetPoint(g->GetN(), qx, qy);
  }

  TCanvas cMaps("cMaps", "cMaps", 3 * 340, nConfig * 310);
  cMaps.Divide(3, nConfig, 0.001, 0.001);
  for (int c = 0; c < nConfig; ++c) {
    double rowMax = 0.0;
    for (int k = 0; k < 3; ++k) {
      rowMax = std::max(rowMax, hist[3 * c + k]->GetMaximum());
    }
    for (int k = 0; k < 3; ++k) {
      auto* pad = cMaps.cd(3 * c + k + 1);
      pad->SetRightMargin(0.15);
      pad->SetLeftMargin(0.13);
      pad->SetFrameFillColor(TColor::GetColorPalette(0)); // cells beyond 5 widths are exactly 0
      TH2D* h = hist[3 * c + k];
      h->SetMaximum(rowMax);
      h->SetMinimum(0.0);
      h->GetZaxis()->SetTitle("T_{p} (fm^{-2})");
      h->Draw("COLZ");
      TGraph* gq = quarkGraph[3 * c + k];
      gq->SetMarkerStyle(24);
      gq->SetMarkerSize(1.1);
      gq->SetMarkerColor(kWhite);
      gq->Draw("P SAME");
      if (k == 1) {
        TGraph* gj = junctionGraph[3 * c + k];
        gj->SetMarkerStyle(29);
        gj->SetMarkerSize(1.6);
        gj->SetMarkerColor(kRed + 1);
        gj->Draw("P SAME");
      }
    }
  }
  save(cMaps, outDir, "proton_maps");

  // ---------------------------------------------------------------- statistics
  TTree stats("stats", "stats");
  stats.ReadFile(Form("%s_stats.txt", stem));
  stats.SetEstimate(stats.GetEntries() + 1);
  const double nStats = static_cast<double>(stats.GetEntries());

  TCanvas cStats("cStats", "cStats", 1150, 950);
  cStats.Divide(2, 2);
  for (int i = 1; i <= 4; ++i) {
    cStats.cd(i)->SetLeftMargin(0.14);
  }

  cStats.cd(1);
  auto* hLY = new TH1D("hLY", ";string length (fm);probability", 80, 0.0, 5.0);
  auto* hLD = new TH1D("hLD", ";string length (fm);probability", 80, 0.0, 5.0);
  stats.Draw("Lsteiner>>hLY", "", "goff");
  stats.Draw("perimeter>>hLD", "", "goff");
  hLY->Scale(1.0 / nStats);
  hLD->Scale(1.0 / nStats);
  hLY->SetLineColor(shapeColor(1));
  hLD->SetLineColor(shapeColor(2));
  hLY->SetLineWidth(2);
  hLD->SetLineWidth(2);
  hLY->SetMaximum(1.15 * std::max(hLY->GetMaximum(), hLD->GetMaximum()));
  hLY->Draw("HIST");
  hLD->Draw("HIST SAME");
  auto* legL = new TLegend(0.55, 0.72, 0.88, 0.88);
  legL->SetBorderSize(0);
  legL->AddEntry(hLY, Form("Y: L_{min}, #LTL#GT = %.2f fm", hLY->GetMean()), "l");
  legL->AddEntry(hLD, Form("#Delta: perimeter, #LTL#GT = %.2f fm", hLD->GetMean()), "l");
  legL->Draw();

  cStats.cd(2);
  auto* hRatio = new TH1D("hRatio", ";L_{Y} / L_{#Delta};probability", 80, 0.49, 0.59);
  stats.Draw("Lsteiner/perimeter>>hRatio", "", "goff");
  hRatio->Scale(1.0 / nStats);
  hRatio->SetFillColorAlpha(shapeColor(1), 0.35);
  hRatio->SetLineColor(shapeColor(1));
  hRatio->Draw("HIST");
  const double equilateral = 1.0 / std::sqrt(3.0);
  auto* line = new TLine(equilateral, 0.0, equilateral, hRatio->GetMaximum());
  line->SetLineStyle(2);
  line->Draw();
  const double atQuark = stats.GetEntries("atQuark == 1") / nStats;
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.042);
  latex.DrawLatex(0.15, 0.84, "dashed: equilateral, 1/#sqrt{3}");
  latex.DrawLatex(0.15, 0.78, Form("junction on a quark (#geq 120#circ): %.1f%%", 100.0 * atQuark));

  cStats.cd(3);
  TH1D* hEps[3];
  double epsMax = 0.0;
  for (int k = 0; k < 3; ++k) {
    hEps[k] = new TH1D(Form("hEps2_%d", k), ";single-proton #varepsilon_{2};probability", 50, 0.0, 1.0);
    stats.Draw(Form("eps2_%s>>hEps2_%d", kShapeName[k], k), "", "goff");
    hEps[k]->Scale(1.0 / nStats);
    hEps[k]->SetLineColor(shapeColor(k));
    hEps[k]->SetLineWidth(2);
    epsMax = std::max(epsMax, hEps[k]->GetMaximum());
  }
  auto* legE = new TLegend(0.15, 0.70, 0.55, 0.88);
  legE->SetBorderSize(0);
  for (int k = 0; k < 3; ++k) {
    hEps[k]->SetMaximum(1.2 * epsMax);
    hEps[k]->Draw(k == 0 ? "HIST" : "HIST SAME");
    legE->AddEntry(hEps[k], Form("%s, #LT#varepsilon_{2}#GT = %.3f", kShapeLabel[k], hEps[k]->GetMean()), "l");
  }
  legE->Draw();

  cStats.cd(4);
  TTree profile("profile", "profile");
  profile.ReadFile(Form("%s_profile.txt", stem));
  auto* frame = new TH1D("frame", ";r (fm);#LTT_{p}#GT(r) (fm^{-2})", 10, 0.0, 1.5);
  frame->SetMaximum(1.1 * std::max(profile.GetMaximum("Y"), profile.GetMaximum("Delta")));
  frame->SetMinimum(0.0);
  frame->Draw();
  auto* legP = new TLegend(0.50, 0.66, 0.88, 0.88);
  legP->SetBorderSize(0);
  for (int k = 0; k < 3; ++k) {
    profile.Draw(Form("r:%s", kShapeName[k]), "", "goff");
    auto* g = new TGraph(profile.GetSelectedRows(), profile.GetV1(), profile.GetV2());
    g->SetLineColor(shapeColor(k));
    g->SetLineWidth(2);
    g->Draw("L SAME");
    legP->AddEntry(g, kShapeLabel[k], "l");
  }
  profile.Draw("r:hotspotAnalytic", "", "goff");
  auto* ga = new TGraph(profile.GetSelectedRows(), profile.GetV1(), profile.GetV2());
  ga->SetLineStyle(2);
  ga->SetLineWidth(2);
  ga->Draw("L SAME");
  legP->AddEntry(ga, "hot spots, analytic", "l");
  legP->Draw();

  save(cStats, outDir, "proton_stats");
}
