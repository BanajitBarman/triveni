// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

// plotPP.C: pp figures for the hot-spot, Y and Delta proton shapes.
//
//   root -l -b -q 'macros/plotPP.C("output/pp200", "output/figures", "#sqrt{s} = 200 GeV")'
//
// Reads <prefix>_{hotspot,Y,Delta}.txt from triveni-pp.
// pp_eccentricity: multiplicity distribution, eps2{2} and eps3{2} vs S/<S>, b of inelastic events.
// pp_events:       entropy density of the first events with the Y proton.

#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TProfile.h>
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

const int kShapeMarker[3] = {20, 21, 33};

void save(TCanvas& c, const char* outDir, const char* name)
{
  c.SaveAs(Form("%s/%s.pdf", outDir, name));
  c.SaveAs(Form("%s/%s.png", outDir, name));
}

/// eps_n{2} = sqrt(<eps_n^2>) per bin, with the error propagated from the profile.
TGraphErrors* rmsGraph(const TProfile* p)
{
  auto* g = new TGraphErrors();
  for (int i = 1; i <= p->GetNbinsX(); ++i) {
    if (p->GetBinEntries(i) < 20) {
      continue;
    }
    const double m = p->GetBinContent(i);
    const int point = g->GetN();
    g->SetPoint(point, p->GetBinCenter(i), std::sqrt(m));
    g->SetPointError(point, 0.0, p->GetBinError(i) / (2.0 * std::sqrt(m)));
  }
  return g;
}

} // namespace

void plotPP(const char* prefix = "output/pp200", const char* outDir = "output/figures", const char* label = "")
{
  gSystem->mkdir(outDir, true);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  gStyle->SetNumberContours(64);

  const std::vector<double> edges = {0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0};
  TH1D* hMult[3];
  TH1D* hB[3];
  TGraphErrors* gEps[2][3];
  int nEvents[3] = {0, 0, 0};

  for (int k = 0; k < 3; ++k) {
    TTree tree(Form("pp_%s", kShapeName[k]), "pp");
    tree.ReadFile(Form("%s_%s.txt", prefix, kShapeName[k]));
    tree.SetEstimate(tree.GetEntries() + 1);
    nEvents[k] = static_cast<int>(tree.GetEntries());
    tree.Draw("S", "", "goff");
    double meanS = 0.0;
    for (Long64_t i = 0; i < tree.GetSelectedRows(); ++i) {
      meanS += tree.GetV1()[i];
    }
    meanS /= tree.GetSelectedRows();

    hMult[k] = new TH1D(Form("hMult_%d", k), ";S / #LTS#GT;probability density", 60, 0.0, 6.0);
    tree.Draw(Form("S/%.10g>>hMult_%d", meanS, k), "", "goff");
    hMult[k]->Scale(1.0 / (nEvents[k] * hMult[k]->GetBinWidth(1)));

    hB[k] = new TH1D(Form("hB_%d", k), ";b (fm);dP/db of inelastic events (fm^{-1})", 70, 0.0, 3.5);
    tree.Draw(Form("b>>hB_%d", k), "", "goff");
    hB[k]->Scale(1.0 / (nEvents[k] * hB[k]->GetBinWidth(1)));

    for (int n = 2; n <= 3; ++n) {
      auto* prof = new TProfile(Form("pEps%d_%d", n, k), "", static_cast<int>(edges.size()) - 1, edges.data());
      tree.Draw(Form("eps%d*eps%d:S/%.10g>>pEps%d_%d", n, n, meanS, n, k), "", "goff prof");
      gEps[n - 2][k] = rmsGraph(prof);
    }
  }

  TCanvas c("cPP", "cPP", 1150, 950);
  c.Divide(2, 2);
  for (int i = 1; i <= 4; ++i) {
    c.cd(i)->SetLeftMargin(0.14);
  }
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.045);

  c.cd(1)->SetLogy();
  auto* leg = new TLegend(0.55, 0.68, 0.88, 0.88);
  leg->SetBorderSize(0);
  for (int k = 0; k < 3; ++k) {
    hMult[k]->SetLineColor(shapeColor(k));
    hMult[k]->SetLineWidth(2);
    hMult[k]->SetMinimum(1e-4);
    hMult[k]->Draw(k == 0 ? "HIST" : "HIST SAME");
    leg->AddEntry(hMult[k], Form("%s (%d ev.)", kShapeLabel[k], nEvents[k]), "l");
  }
  leg->Draw();
  latex.DrawLatex(0.18, 0.20, Form("pp %s", label));

  for (int n = 2; n <= 3; ++n) {
    c.cd(n);
    auto* frame = new TH1D(Form("frame%d", n), Form(";S / #LTS#GT;#varepsilon_{%d}{2}", n), 10, 0.0, 6.0);
    double top = 0.0;
    for (int k = 0; k < 3; ++k) {
      for (int i = 0; i < gEps[n - 2][k]->GetN(); ++i) {
        top = std::max(top, gEps[n - 2][k]->GetY()[i]);
      }
    }
    frame->SetMinimum(0.0);
    frame->SetMaximum(1.25 * top);
    frame->Draw();
    for (int k = 0; k < 3; ++k) {
      TGraphErrors* g = gEps[n - 2][k];
      g->SetMarkerStyle(kShapeMarker[k]);
      g->SetMarkerColor(shapeColor(k));
      g->SetLineColor(shapeColor(k));
      g->SetMarkerSize(1.2);
      g->Draw("PL SAME");
    }
    leg->Draw();
  }

  c.cd(4);
  double bTop = 0.0;
  for (int k = 0; k < 3; ++k) {
    bTop = std::max(bTop, hB[k]->GetMaximum());
  }
  for (int k = 0; k < 3; ++k) {
    hB[k]->SetLineColor(shapeColor(k));
    hB[k]->SetLineWidth(2);
    hB[k]->SetMaximum(1.2 * bTop);
    hB[k]->Draw(k == 0 ? "HIST" : "HIST SAME");
  }
  leg->Draw();
  save(c, outDir, "pp_eccentricity");

  // Event displays for the Y proton.
  TTree maps("maps", "maps");
  maps.ReadFile(Form("%s_Y_maps.txt", prefix));
  TTree quarks("quarks", "quarks");
  quarks.ReadFile(Form("%s_Y_quarks.txt", prefix));
  const int nMaps = static_cast<int>(maps.GetMaximum("event")) + 1;
  const int n = static_cast<int>(std::lround(std::sqrt(static_cast<double>(maps.GetEntries()) / nMaps)));
  const double xmin = maps.GetMinimum("x");
  const double xmax = maps.GetMaximum("x");
  const double step = (xmax - xmin) / (n - 1);

  TCanvas cEv("cEv", "cEv", 330 * nMaps, 330);
  cEv.Divide(nMaps, 1, 0.001, 0.001);
  for (int e = 0; e < nMaps; ++e) {
    auto* h = new TH2D(Form("ev_%d", e), Form("event %d;x (fm);y (fm)", e), n, xmin - 0.5 * step, xmax + 0.5 * step, n,
                       xmin - 0.5 * step, xmax + 0.5 * step);
    maps.Project(Form("ev_%d", e), "y:x", Form("s * (event == %d)", e));
    auto* pad = cEv.cd(e + 1);
    pad->SetRightMargin(0.15);
    pad->SetLeftMargin(0.14);
    pad->SetFrameFillColor(TColor::GetColorPalette(0));
    h->GetXaxis()->SetRangeUser(-2.2, 2.2);
    h->GetYaxis()->SetRangeUser(-2.2, 2.2);
    h->GetZaxis()->SetTitle("s (fm^{-2})");
    h->Draw("COLZ");
    for (int item = 0; item < 4; item += 3) {
      quarks.Draw("y:x", Form("event == %d && item %s 3", e, item == 3 ? "==" : "<"), "goff");
      auto* g = new TGraph(quarks.GetSelectedRows(), quarks.GetV2(), quarks.GetV1());
      g->SetMarkerStyle(item == 3 ? 29 : 24);
      g->SetMarkerSize(item == 3 ? 1.5 : 1.0);
      g->SetMarkerColor(item == 3 ? kRed + 1 : kWhite);
      g->Draw("P SAME");
    }
  }
  save(cEv, outDir, "pp_events");
}
