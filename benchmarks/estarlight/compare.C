// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Banajit Barman

// compare.C: TRIVENI vs eSTARlight, coherent J/psi at the EIC (18 GeV e on 110 GeV/u Au, 275 GeV p).
//
//   root -l -b -q 'benchmarks/estarlight/compare.C("output/bench", "output/figures")'
//
// Reads from <dir>:
//   triveni_{eAu,ep,ep_w89}_{q2,t,fit}.txt          (triveni-exclusive)
//   estarlight_{eAu,ep,ep_w89,eAu_q2test}_{q2,t}.txt (analyze_estarlight.py)
// bench_sigma.pdf/png: gamma p fit; e+p and e+Au dsigma/dQ2; Au/p ratio vs Q2.
// bench_t.pdf/png:     e+Au dsigma/d|t| in three Q2 bins.

#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace
{

int colour(const char* hex) { return TColor::GetColor(hex); }

struct Bins {
  std::vector<double> lo, hi, sigma, err;
};

Bins readQ2(const std::string& path, bool hasErr)
{
  TTree t("q2", "q2");
  t.ReadFile(path.c_str());
  Double_t lo, hi, s, e = 0;
  t.SetBranchAddress("q2lo", &lo);
  t.SetBranchAddress("q2hi", &hi);
  t.SetBranchAddress("sigma", &s);
  if (hasErr) {
    t.SetBranchAddress("err", &e);
  }
  Bins b;
  for (Long64_t i = 0; i < t.GetEntries(); ++i) {
    t.GetEntry(i);
    b.lo.push_back(lo);
    b.hi.push_back(hi);
    b.sigma.push_back(s);
    b.err.push_back(e);
  }
  return b;
}

/// dsigma/dQ2 = sigma_bin / width at the log-centre of each bin, skipping the first (kinematic-minimum) bin.
TGraphErrors* dsdq2(const Bins& b)
{
  auto* g = new TGraphErrors();
  for (std::size_t i = 1; i < b.lo.size(); ++i) {
    const double w = b.hi[i] - b.lo[i];
    const int p = g->GetN();
    g->SetPoint(p, std::sqrt(b.lo[i] * b.hi[i]), b.sigma[i] / w);
    g->SetPointError(p, 0.0, b.err[i] / w);
  }
  return g;
}

TGraphErrors* ratioGraph(const Bins& a, const Bins& p)
{
  auto* g = new TGraphErrors();
  for (std::size_t i = 1; i < a.lo.size(); ++i) {
    const double r = a.sigma[i] / p.sigma[i];
    const double e = r * std::hypot(a.err[i] / a.sigma[i], p.err[i] / p.sigma[i]);
    const int k = g->GetN();
    g->SetPoint(k, std::sqrt(a.lo[i] * a.hi[i]), r);
    g->SetPointError(k, 0.0, e);
  }
  return g;
}

void style(TGraph* g, int col, int marker, int lineStyle = 1)
{
  g->SetLineColor(col);
  g->SetMarkerColor(col);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
  g->SetLineStyle(lineStyle);
}

TH1D* frame(const char* name, const char* title, double x0, double x1, double y0, double y1)
{
  auto* h = new TH1D(name, title, 10, x0, x1);
  h->SetMinimum(y0);
  h->SetMaximum(y1);
  h->SetStats(false);
  return h;
}

} // namespace

void compare(const char* dir = "output/bench", const char* outDir = "output/figures")
{
  gSystem->mkdir(outDir, true);
  gStyle->SetOptStat(0);
  const std::string d = dir;
  const int cTriveni = colour("#2F6FDB");   // accent: TRIVENI
  const int cEstar = colour("#555555");     // eSTARlight as released
  const int cTest = colour("#E0453A");      // eSTARlight with the sqrt(Q2-factor) test

  const Bins tAu = readQ2(d + "/triveni_eAu_q2.txt", false);
  const Bins tP = readQ2(d + "/triveni_ep_q2.txt", false);
  const Bins tP89 = readQ2(d + "/triveni_ep_w89_q2.txt", false);
  const Bins eAu = readQ2(d + "/estarlight_eAu_q2.txt", true);
  const Bins eP = readQ2(d + "/estarlight_ep_q2.txt", true);
  const Bins eP89 = readQ2(d + "/estarlight_ep_w89_q2.txt", true);
  const bool haveTest = !gSystem->AccessPathName((d + "/estarlight_eAu_q2test_q2.txt").c_str());
  Bins eAuTest;
  if (haveTest) {
    eAuTest = readQ2(d + "/estarlight_eAu_q2test_q2.txt", true);
  }

  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.042);

  TCanvas c("c", "c", 1150, 950);
  c.Divide(2, 2);
  for (int i = 1; i <= 4; ++i) {
    c.cd(i)->SetLeftMargin(0.15);
    c.cd(i)->SetBottomMargin(0.12);
  }

  // (a) gamma p -> J/psi p: fit
  {
    auto* pad = c.cd(1);
    pad->SetLogx();
    TTree fit("fit", "fit");
    fit.ReadFile((d + "/triveni_ep_fit.txt").c_str());
    fit.Draw("W:estarlight:triveni", "", "goff");
    auto* gE = new TGraph(fit.GetSelectedRows(), fit.GetV1(), fit.GetV2());
    auto* gT = new TGraph(fit.GetSelectedRows(), fit.GetV1(), fit.GetV3());
    frame("fa", ";W_{#gammap} (GeV);#sigma(#gamma p #rightarrow J/#psi p) (nb)", 20, 250, 0, 160)->Draw();
    style(gE, cEstar, 1, 1);
    style(gT, cTriveni, 20);
    gE->Draw("L SAME");
    gT->Draw("P SAME");
    auto* leg = new TLegend(0.18, 0.72, 0.62, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(gE, "eSTARlight input (HERA fit)", "l");
    leg->AddEntry(gT, "TRIVENI dipole model, fitted", "p");
    leg->Draw();
    label.DrawLatex(0.18, 0.92, "(a) common e+p anchor at Q^{2} = 0");
  }

  // (b) e+p and (c) e+Au: dsigma/dQ2, W > 40 GeV
  auto drawDsdq2 = [&](int padIndex, const Bins& tri, const Bins& est, const Bins* test, const char* title, double y0, double y1) {
    auto* pad = c.cd(padIndex);
    pad->SetLogx();
    pad->SetLogy();
    frame(Form("fq%d", padIndex), ";Q^{2} (GeV^{2});d#sigma/dQ^{2} (nb/GeV^{2})", 5e-4, 50, y0, y1)->Draw();
    auto* gE = dsdq2(est);
    auto* gT = dsdq2(tri);
    style(gE, cEstar, 20);
    style(gT, cTriveni, 21);
    gT->SetMarkerStyle(1);
    gT->Draw("L SAME");
    gE->Draw("P SAME");
    auto* leg = new TLegend(0.18, 0.18, 0.62, 0.36);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(gE, "eSTARlight", "p");
    if (test) {
      auto* gX = dsdq2(*test);
      style(gX, cTest, 24);
      gX->Draw("P SAME");
      leg->AddEntry(gX, "eSTARlight, #sqrt{Q^{2} factor} on #sigma_{VN}", "p");
    }
    leg->AddEntry(gT, "TRIVENI", "l");
    leg->Draw();
    label.DrawLatex(0.18, 0.92, title);
  };
  drawDsdq2(2, tP, eP, nullptr, "(b) e+p #rightarrow e+J/#psi+p, W > 40 GeV", 1e-3, 1e3);
  drawDsdq2(3, tAu, eAu, haveTest ? &eAuTest : nullptr, "(c) e+Au #rightarrow e+J/#psi+Au, W > 40 GeV", 1e-2, 1e6);

  // (d) Au/p ratio in the common W window 40-89 GeV
  {
    auto* pad = c.cd(4);
    pad->SetLogx();
    frame("fr", ";Q^{2} (GeV^{2});#sigma(e+Au) / #sigma(e+p),  40 < W < 89 GeV", 5e-4, 50, 0, 480)->Draw();
    auto* gE = ratioGraph(eAu, eP89);
    auto* gT = ratioGraph(tAu, tP89);
    style(gE, cEstar, 20);
    style(gT, cTriveni, 21);
    gT->Draw("LP SAME");
    gE->Draw("P SAME");
    auto* leg = new TLegend(0.18, 0.18, 0.62, 0.36);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(gE, "eSTARlight", "p");
    if (haveTest) {
      auto* gX = ratioGraph(eAuTest, eP89);
      style(gX, cTest, 24);
      gX->Draw("P SAME");
      leg->AddEntry(gX, "eSTARlight, #sqrt{Q^{2} factor} test", "p");
    }
    leg->AddEntry(gT, "TRIVENI", "lp");
    leg->Draw();
    label.DrawLatex(0.18, 0.92, "(d) nuclear/proton ratio");
  }
  c.SaveAs(Form("%s/bench_sigma.pdf", outDir));
  c.SaveAs(Form("%s/bench_sigma.png", outDir));

  // dsigma/d|t| for e+Au in three Q2 bins
  TTree tt("tt", "tt");
  tt.ReadFile((d + "/triveni_eAu_t.txt").c_str());
  TTree et("et", "et");
  et.ReadFile((d + "/estarlight_eAu_t.txt").c_str());
  TTree xt("xt", "xt");
  if (haveTest) {
    xt.ReadFile((d + "/estarlight_eAu_q2test_t.txt").c_str());
  }
  tt.SetEstimate(tt.GetEntries() + 1);
  et.SetEstimate(et.GetEntries() + 1);
  xt.SetEstimate(xt.GetEntries() + 1);
  const int bins[3] = {2, 4, 6};
  const char* names[3] = {"0.01 < Q^{2} < 0.1 GeV^{2}", "1 < Q^{2} < 3 GeV^{2}", "10 < Q^{2} < 30 GeV^{2}"};
  TCanvas ct("ct", "ct", 1500, 520);
  ct.Divide(3, 1);
  for (int p = 0; p < 3; ++p) {
    auto* pad = ct.cd(p + 1);
    pad->SetLogy();
    pad->SetLeftMargin(0.16);
    pad->SetBottomMargin(0.12);
    tt.Draw("absT:dsdt", Form("bin==%d && absT<=0.12", bins[p]), "goff");
    auto* gT = new TGraph(tt.GetSelectedRows(), tt.GetV1(), tt.GetV2());
    et.Draw("absT:dsdt:err", Form("bin==%d && dsdt>0", bins[p]), "goff");
    std::vector<double> ex(et.GetSelectedRows(), 0.0);
    auto* gE = new TGraphErrors(et.GetSelectedRows(), et.GetV1(), et.GetV2(), ex.data(), et.GetV3());
    double top = 0;
    for (int i = 0; i < gT->GetN(); ++i) {
      top = std::max(top, gT->GetY()[i]);
    }
    for (int i = 0; i < gE->GetN(); ++i) {
      top = std::max(top, gE->GetY()[i]);
    }
    frame(Form("ft%d", p), ";|t| (GeV^{2});d#sigma/d|t| (nb/GeV^{2})", 0, 0.12, top * 1e-5, top * 3)->Draw();
    style(gT, cTriveni, 1);
    style(gE, cEstar, 20);
    gE->SetMarkerSize(0.7);
    gT->Draw("L SAME");
    gE->Draw("P SAME");
    auto* leg = new TLegend(0.45, 0.70, 0.88, 0.86);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(gE, "eSTARlight", "p");
    if (haveTest) {
      xt.Draw("absT:dsdt", Form("bin==%d && dsdt>0", bins[p]), "goff");
      auto* gX = new TGraph(xt.GetSelectedRows(), xt.GetV1(), xt.GetV2());
      style(gX, cTest, 24);
      gX->SetMarkerSize(0.6);
      gX->Draw("P SAME");
      leg->AddEntry(gX, "eSTARlight, test", "p");
    }
    leg->AddEntry(gT, "TRIVENI", "l");
    leg->Draw();
    label.DrawLatex(0.18, 0.92, Form("e+Au, W > 40 GeV, %s", names[p]));
  }
  ct.SaveAs(Form("%s/bench_t.pdf", outDir));
  ct.SaveAs(Form("%s/bench_t.png", outDir));
}
