#include "TRandom3.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooLegendre.h"
#include "RooExponential.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TString.h>
#include <TSystemDirectory.h>
#include "Riostream.h"
#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <TLine.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTimeStamp.h>
#include <TRandom3.h>
#include <TFitResult.h>

#include "NicerRooPlot.c"
#include "Nicer1D.c"
#include "Nicer2D.c"
#include "Nice1D.c"
#include "Nice1Db.c"

using namespace RooFit;

void sig_pdf(){
    
  gROOT->Reset();
  
  gROOT->SetStyle("Bold");
  gStyle->SetCanvasColor(0);
  gStyle->SetLabelColor(1);
  gStyle->SetLabelColor(1,"Y");
  gStyle->SetHistLineColor(1); 
  gStyle->SetHistLineWidth(1); 
  gStyle->SetNdivisions(505);
  gStyle->SetNdivisions(505,"Y");
  //gROOT->Macro("setcolor2.c");
  gStyle->SetHistFillColor(999);
  gROOT->SetStyle("Plain");  // white as bg
  gStyle->SetOptStat("111111");
  gStyle->SetFrameLineWidth(1);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleBorderSize(0);
    
  //Definitions
  Double_t smallBetween1 = .1;
  Double_t smallBetween2 = .1;
  Double_t smallBetween3 = .1;
  Double_t smallBetween4 = .1;
  
  //Double_t small = .00001;
  
  TLine TLine;
  TLatex *t = new TLatex();
  t->SetTextSize(0.03);
  t->SetTextFont(42);
  t->SetTextAlign(12);
  t->SetTextColor(1);
  t->SetTextFont(12);
   
  TCanvas *C1;
  
  //TPad *SmallC;
  //TGaxis *XAxis,*YAxis;
  //TLatex XTitle,YTitle,XtopTitle;
  
  //TLegend *legend;

  ifstream in;
  ofstream out;
  TString Line = ""; 
  TString cleg = "";

  TTimeStamp * time_st = new TTimeStamp();
  double_t timeseed = time_st->GetNanoSec();
  gRandom->SetSeed(timeseed);
  
  TString path = "./";

  int bin_range = 1000;
  double x_min = 0.5;
  double x_max = 1.5;

  TF1 * fct_width = new TF1("fct_width", "pol1", x_min, x_max);
  fct_width->SetParameter(0, 4);
  fct_width->SetParameter(1, 2);

  TGraph * gr_width = new TGraph();

  for (int i = 0; i < 40; i ++) {
    double mass = x_min + (x_max - x_min) / 40. * i + (x_max - x_min) / 40. / 2.;
    double width = fct_width->Eval(mass);
    gr_width->SetPoint(i, mass, width);
  }

  TH2F * DrawW = new TH2F("DrawW", ";#font[42]{#it{m}_{X} [GeV/#it{c}^{2}]};#font[42]{Width MeV/#it{c}^{2}};", 100, x_min, x_max, 100, 0., 9.99);
  Nicer2D(DrawW, 0.05, 42, 505, 1.2, 1.4, 1.2);
  
  smallBetween1 = .15;
  smallBetween2 = .05;
  smallBetween3 = .05;
  smallBetween4 = .15;
    
  cleg = "Signal_PDF";
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  DrawW->Draw();
  gr_width->Draw("CS");
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".png");

  int nevt = 1e3;
  double m = 1.0;
  double w = fct_width->Eval(m) * 1e-3;
  double x_min_fit = m - 5.0 * w;
  double x_max_fit = m + 5.0 * w;
  
  TF1 * fct_sig = new TF1("fct_sig", "gaus", x_min, x_max);
  fct_sig->SetParameter(0, 1e2);
  fct_sig->SetParameter(1, m);
  fct_sig->SetParameter(2, w);
  
  TH1F * h_m_sig = new TH1F("h_m_sig", ";#font[42]{#it{m}_{X} [GeV/#it{c}^{2}]};#font[42]{Events / MeV/#it{c}^{2}}", bin_range, x_min, x_max);
  h_m_sig->Sumw2();
  h_m_sig->FillRandom("fct_sig", nevt);
  Nicer1D(h_m_sig, 0.05, 42, 506, 1.1, 1.4);

  cleg ="Signal";
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  h_m_sig->GetXaxis()->SetRangeUser(x_min_fit , x_max_fit);
  h_m_sig->SetMarkerStyle(20);
  h_m_sig->SetMarkerSize(1.2);
  h_m_sig->SetMarkerColor(1);
  h_m_sig->SetLineColor(1);
  h_m_sig->Draw();
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".png");

  TString TITLE = "#font[42]{#it{m}_{X} [GeV/#it{c}^{2}]}";
  RooRealVar x("x", TITLE, x_min_fit, x_max_fit);
  
  RooRealVar mean("#it{m}", "mass", m, m, m);
  RooRealVar width("w", "width", w, w, w);
  RooGaussian gau_pdf("gau","gau", x, mean, width);
  
  RooRealVar nsig("N_{sig}", "N_{sig}", nevt, 0.0, 2.0 * nevt);
  RooAddPdf model_sig("model_sig", "Signal model", RooArgList(gau_pdf), RooArgList(nsig));
  RooDataHist data_sig("data_sig","dataset_sig", x, h_m_sig);

  RooPlot *xframe_sig = x.frame(Range(x_min_fit, x_max_fit),Title("#font[42]{}"));
  data_sig.plotOn(xframe_sig, XErrorSize(0), DataError(RooAbsData::SumW2)) ; 
  model_sig.fitTo(data_sig);
  model_sig.paramOn(xframe_sig,Format(("NEU"),AutoPrecision(1)),Layout(0.65,0.90,0.90));
  xframe_sig->getAttText()->SetTextSize(0.03);
  model_sig.plotOn(xframe_sig, Components(gau_pdf), LineColor(kGreen),LineStyle(kDashed));
  model_sig.plotOn(xframe_sig);
  
  RooHist* hresid_sig = xframe_sig->residHist() ;
  RooHist* hpull_sig  = xframe_sig->pullHist() ;
  RooPlot* xframe_sig_res = x.frame(Title("#font[42]{Residual}")) ;
  xframe_sig_res->addPlotable(hresid_sig,"P") ;
  RooPlot* xframe_sig_pull = x.frame(Title("#font[42]{}")) ;
  xframe_sig_pull->addPlotable(hpull_sig,"P") ;

  cleg = "Fit-signal";
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  C1->Divide(1,2);
  C1->GetPad(1)->SetPad(0, 0.3, 1., 1.0);
  C1->GetPad(2)->SetPad(0, 0.0, 1., 0.3);
  C1->cd(1);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  NicerRooPlot(xframe_sig, 0.05, 42, 505, 1.3, 1.3);
  xframe_sig->Draw("EZ");
  C1->cd(2);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  xframe_sig_pull->SetMaximum(4.99);
  xframe_sig_pull->SetMinimum(-4.99);
  NicerRooPlot(xframe_sig_pull, 0.1, 42, 505, 1.3, 0.3);
  xframe_sig_pull->GetXaxis()->CenterTitle(kTRUE);
  xframe_sig_pull->GetYaxis()->CenterTitle(kTRUE);
  xframe_sig_pull->SetYTitle("#font[42]{Pull (#sigma)}");
  xframe_sig_pull->SetXTitle("#font[42]{}");
  xframe_sig_pull->Draw("EZ");
  TLine.SetLineColor(2);
  TLine.SetLineWidth(3);
  TLine.SetLineStyle(2);
  TLine.DrawLine(x_min_fit, 0, x_max_fit, 0);
  TLine.SetLineStyle(3);
  TLine.DrawLine(x_min_fit, -3, x_max_fit, -3);
  TLine.DrawLine(x_min_fit, 3, x_max_fit, 3);
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".png");

  RooMCStudy* mcstudy_sig = new RooMCStudy(model_sig, x, Binned(kTRUE) ,Silence(), Extended(),FitOptions(Save(kTRUE),PrintEvalErrors(0))) ;
  mcstudy_sig->generateAndFit(1000) ;
  
  RooPlot* frame_mc_nsig_1 = mcstudy_sig->plotParam(nsig, Bins(50));
  frame_mc_nsig_1->SetXTitle("#font[42]{N_{sig}}");                                                                                                                                                    
  RooPlot* frame_mc_nsig_2 = mcstudy_sig->plotError(nsig, Bins(50));
  frame_mc_nsig_2->SetXTitle("#font[42]{Error}"); 
  RooPlot* frame_mc_nsig_3 = mcstudy_sig->plotPull(nsig, -10,  10 , 500, kTRUE);     //mean,-3,3,40,kTRUE
  frame_mc_nsig_3->SetXTitle("#font[42]{Pull}");
  RooPlot* frame_mc_nsig_4 = mcstudy_sig->plotNLL(Bins(50)) ;
  frame_mc_nsig_4->SetXTitle("#font[42]{-log(Likelihood)}");
  
  cleg = "nsig_mcstudy";
  C1 = new TCanvas(cleg, cleg, 10, 10, 1600, 400);
  C1->Divide(4,1);
  C1->cd(1);
  gPad->SetLeftMargin(smallBetween1);                                                                                                                                                                  
  gPad->SetRightMargin(smallBetween2);                                                                                                                                                                 
  gPad->SetTopMargin(smallBetween3);                                                                                                                                                                   
  gPad->SetBottomMargin(smallBetween4);                                                                                                                                                                
  NicerRooPlot(frame_mc_nsig_1, 0.05, 42, 504, 1.1, 1.6);
  frame_mc_nsig_1->Draw();
  C1->cd(2);                                                                                                                                                                                           
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);                                                                                                                                                                
  NicerRooPlot(frame_mc_nsig_2, 0.05, 42, 504, 1.1, 1.6);
  frame_mc_nsig_2->Draw();
  C1->cd(3);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  NicerRooPlot(frame_mc_nsig_3, 0.05, 42, 504, 1.1, 1.6);
  frame_mc_nsig_3->Draw();
  C1->cd(4);                                                                                                                                                                                           
  gPad->SetLeftMargin(smallBetween1);                                                                                                                                                                  
  gPad->SetRightMargin(smallBetween2);                                                                                                                                                                 
  gPad->SetTopMargin(smallBetween3);                                                                                                                                                                   
  gPad->SetBottomMargin(smallBetween4);                                                                                                                                                                
  NicerRooPlot(frame_mc_nsig_4,0.05, 42, 504, 1.1, 1.6);
  frame_mc_nsig_4->Draw() ;
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".png");
}
