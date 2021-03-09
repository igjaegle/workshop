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


#include "RooStats/HybridCalculatorOriginal.h"
#include "RooStats/HybridPlot.h"
#include "RooStats/HypoTestInverterOriginal.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include "NicerRooPlot.c"
#include "Nicer1D.c"
#include "Nice1D.c"
#include "Nice1Db.c"

using namespace RooFit;
using namespace RooStats ;

void limit(){
    
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
  
  TF1 * fct_bkg = new TF1("fct_bkg", "pol2", x_min, x_max);
  fct_bkg->SetParameter(0, 100);
  fct_bkg->SetParameter(1, 1);
  fct_bkg->SetParameter(2, -10);
    
  TH1F * h_m_bkg1 = new TH1F("h_m_bkg1", ";#font[42]{#it{m}_{X} [GeV/#it{c}^{2}]};#font[42]{Events / MeV/#it{c}^{2}}", bin_range, x_min, x_max);
  TH1F * h_m_bkg2 = new TH1F("h_m_bkg2", ";#font[42]{#it{m}_{X} [GeV/#it{c}^{2}]};#font[42]{Events / MeV/#it{c}^{2}}", bin_range, x_min, x_max);
  h_m_bkg1->Sumw2();
  h_m_bkg2->Sumw2();
  Nicer1D(h_m_bkg1, 0.05, 42, 506, 1.1, 1.4);
  Nicer1D(h_m_bkg2, 0.05, 42, 506, 1.1, 1.4);

  int bkg_level = 1;
  for (int i = 0; i < bkg_level; i ++) h_m_bkg1->Fill(fct_bkg->GetRandom());

  bkg_level = 1e5;
  for (int i = 0; i < bkg_level; i ++) h_m_bkg2->Fill(fct_bkg->GetRandom());

  smallBetween1 = .15;
  smallBetween2 = .05;
  smallBetween3 = .05;
  smallBetween4 = .15;
    
  cleg = "Count_experiment";
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  h_m_bkg1->SetMinimum(0);
  h_m_bkg1->SetMarkerStyle(20);
  h_m_bkg1->SetMarkerSize(1.2);
  h_m_bkg1->SetMarkerColor(1);
  h_m_bkg1->SetLineColor(1);
  h_m_bkg1->Draw();
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".png");

  cleg = "Shape_experiment";
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  h_m_bkg2->SetMinimum(0);
  h_m_bkg2->SetMarkerStyle(20);
  h_m_bkg2->SetMarkerSize(1.2);
  h_m_bkg2->SetMarkerColor(1);
  h_m_bkg2->SetLineColor(1);
  h_m_bkg2->Draw();
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".png");

  TF1 * fct_width = new TF1("fct_width", "pol1", x_min, x_max);
  fct_width->SetParameter(0, 4);
  fct_width->SetParameter(1, 2);

  double m = 1.0;
  double w = fct_width->Eval(m) * 1e-3;
  double x_min_fit = m - 45.0 * w;
  double x_max_fit = m + 45.0 * w;
  
  TFitResultPtr r;
  TF1 * fct_bkg_pdf = new TF1("fct_bkg_pdf", "pol2", x_min_fit, x_max_fit);
  r = h_m_bkg2->Fit(fct_bkg_pdf, "RMBQS0");

  double nChi2 = r->Chi2() / (double) r->Ndf();
  cout <<"nChi2 " << nChi2 << endl;
  
  cleg = "Shape_experiment_fitted_by_ROOT";
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  h_m_bkg2->SetMinimum(0);
  h_m_bkg2->SetMarkerStyle(20);
  h_m_bkg2->SetMarkerSize(1.2);
  h_m_bkg2->SetMarkerColor(1);
  h_m_bkg2->SetLineColor(1);
  h_m_bkg2->Draw();
  fct_bkg_pdf->SetLineColor(2);
  fct_bkg_pdf->Draw("same");
  t->DrawLatex(0.6, ((double) h_m_bkg2->GetMaximum()) * 0.3, Form("#font[42]{#chi^{2}/Ndf = %0.3f/%d = %0.3f}", r->Chi2(),  r->Ndf(), nChi2));
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".png");

  TString TITLE = "#font[42]{#it{m}_{X} [GeV/#it{c}^{2}]}";
  RooRealVar x("x", TITLE, x_min_fit, x_max_fit);
  
  RooRealVar c0("c_{0}", "c_{0}", 10, -1000., 1000.);
  RooRealVar c1("c_{1}", "c_{1}", 10, -1000., 1000.);
  RooPolynomial bkg_PDF("bkg_PDF", "bkg_PDF", x, RooArgSet(c0, c1));
  RooDataHist bkg_pdf("bkg_pdf","bkg_pdf", x, h_m_bkg2); 
  RooRealVar nbkg("N_{bkg}", "N_{bkg}", 1.0, 0, 2.0 * h_m_bkg2->GetEntries());   
  RooAddPdf model_bkg("model_bkg", "Background model", RooArgList(bkg_PDF), RooArgList(nbkg));
  model_bkg.fitTo(bkg_pdf);

  RooPlot *xframe_bkg = x.frame(Range(x_min_fit, x_max_fit),Title("#font[42]{}"));
  bkg_pdf.plotOn(xframe_bkg, XErrorSize(0), DataError(RooAbsData::SumW2)) ; 
  //xframe_bkg->getAttText()->SetTextSize(0.03);
  model_bkg.paramOn(xframe_bkg,Format(("NEU"), AutoPrecision(1)), Layout(0.75,1.0,1.0));
  model_bkg.plotOn(xframe_bkg, Components(bkg_PDF), LineColor(kBlue),LineStyle(kDashed));
  model_bkg.plotOn(xframe_bkg, LineColor(kRed));
  
  RooHist * hresid_bkg = xframe_bkg->residHist() ;
  RooHist * hpull_bkg  = xframe_bkg->pullHist() ;
  RooPlot * xframe_bkg_res = x.frame(Title("#font[42]{Residual}")) ;
  xframe_bkg_res->addPlotable(hresid_bkg,"P") ;
  RooPlot * xframe_bkg_pull = x.frame(Title("#font[42]{}")) ;
  xframe_bkg_pull->addPlotable(hpull_bkg,"P") ;
  
  cleg = "Shape_experiment_fitted_by_ROOFIT";
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  C1->Divide(1,2);
  C1->GetPad(1)->SetPad(0, 0.3, 1., 1.0);
  C1->GetPad(2)->SetPad(0, 0.0, 1., 0.3);
  C1->cd(1);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  NicerRooPlot(xframe_bkg, 0.05, 42, 505, 1.3, 1.3);
  xframe_bkg->Draw("EZ");
  t->DrawLatex(0.8, ((double) h_m_bkg2->GetMaximum()) * 0.3, TString::Format("#font[42]{#chi^{2} = %.3f}", xframe_bkg->chiSquare(2)));
  C1->cd(2);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  xframe_bkg_pull->SetMaximum(4.99);
  xframe_bkg_pull->SetMinimum(-4.99);
  NicerRooPlot(xframe_bkg_pull, 0.1, 42, 505, 1.3, 0.6);
  xframe_bkg_pull->SetXTitle("#font[42]{}");
  xframe_bkg_pull->SetYTitle("#font[42]{Pull (#sigma)}");
  xframe_bkg_pull->Draw("EZ");
  TLine.SetLineColor(2);
  TLine.SetLineWidth(3);
  TLine.SetLineStyle(2);
  TLine.DrawLine(x_min_fit, 0, x_max_fit, 0);
  TLine.SetLineStyle(1);
  TLine.DrawLine(x_min_fit, -2, x_max_fit, -2);
  TLine.DrawLine(x_min_fit, 2, x_max_fit, 2);
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".png");

  RooRealVar mean("#it{m}", "mass", m, m, m);
  RooRealVar width("w", "width", w, w, w);
  RooGaussian sig_PDF("gau","gau", x, mean, width);
  RooRealVar nsig("N_{sig}", "N_{sig}", 1, 0.0, 2.0 * nbkg.getValV());
  
  //Combine signal and background
  RooAddPdf model("model", "model", RooArgList(sig_PDF,bkg_PDF), RooArgList(nsig,nbkg));
  
  RooDataHist data("data", "dataset", x, h_m_bkg2);
  
  RooAbsReal * nll = model.createNLL(data,Extended(kTRUE));
  RooMinuit rm(*nll);
  rm.setVerbose(kFALSE);
  rm.migrad();
  rm.minos(); 

  RooPlot* xframe = x.frame(Range(x_min_fit, x_max_fit),Title("#font[42]{}"));
  data.plotOn(xframe, XErrorSize(0), DataError(RooAbsData::SumW2)) ; 

  //Fit model to data and check output
  RooFitResult* fitres = model.fitTo(data, Extended(kTRUE), Save()); //used with nsig
  bool failed = false;
  if(fitres->status() != 0) {
    ofstream status_file;
    cleg = "CONVERGENCE_FAILED";
    status_file.open(cleg, ios::app);
    status_file << "CONVERGENCE FAILED. STATUS = "<<fitres->status() <<endl; 
    status_file.close();
    failed = true;
  }
  model.plotOn(xframe, LineColor(kRed));
  model.plotOn(xframe, Components(bkg_PDF), LineColor(kRed), LineStyle(kDashed)) ;
  model.plotOn(xframe, Components(sig_PDF), LineColor(kGreen), LineStyle(kDashed)) ;
  model.plotOn(xframe, Name("model"), DrawOption("L"));
    
  // Construct a histogram with the residuals of the data w.r.t. the curve
  RooHist* hresid = xframe->residHist() ;
  // Construct a histogram with the pulls of the data w.r.t the curve
  RooHist* hpull = xframe->pullHist() ;
    
  // Create a new frame to draw the residual distribution and add the distribution to the frame
  RooPlot* frame_red = x.frame(Range(x_min_fit, x_max_fit),Title("#font[42]{}")) ;
  frame_red->addPlotable(hresid,"P") ;
  
  // Create a new frame to draw the pull distribution and add the distribution to the frame
  RooPlot* frame_pul = x.frame(Range(x_min_fit, x_max_fit),Title("#font[42]{}")) ;
  frame_pul->addPlotable(hpull, "P") ;
  
  model.paramOn(xframe,Format(("NEU"),AutoPrecision(0)),Layout(0.75,1.0,1.0)) ;
  xframe->getAttText()->SetTextFont(42);
  xframe->getAttText()->SetTextSize(0.03);
  xframe->getAttFill()->SetFillStyle(0);
  xframe->getAttLine()->SetLineColor(0);
  
  smallBetween1 = .15;
  smallBetween2 = .05;
  smallBetween3 = .05;
  smallBetween4 = .15;
  
  cleg = "Mass_Scan";
    
  cout <<"CLEG " << cleg << endl;
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  C1->Divide(1,2);
  C1->GetPad(1)->SetPad(0, 0.3, 1., 1.0);
  C1->GetPad(2)->SetPad(0, 0.0, 1., 0.3);
  C1->cd(1);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  NicerRooPlot(xframe, 0.05, 42, 505, 1.3, 1.3);
  xframe->SetMaximum(h_m_bkg2->GetMaximum() * 1.8);
  xframe->Draw("EZ");
  
  t->SetTextSize(0.04);
    
  t->DrawLatex(0.8, ((double) h_m_bkg2->GetMaximum()) * 0.3, TString::Format("#font[42]{#chi^{2}/ndf = %.2f}", xframe->chiSquare()));
  t->DrawLatex(0.8, ((double) h_m_bkg2->GetMaximum()) * 0.2, TString::Format("#font[42]{#it{m}_{A'} = %.2f GeV/#it{c}^{2}}", m));
  
  C1->cd(2);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  frame_pul->SetMaximum(4.99);
  frame_pul->SetMinimum(-4.99);
  NicerRooPlot(frame_pul, 0.1, 42, 505, 1.3, 1.3);
  frame_pul->GetXaxis()->CenterTitle(kTRUE);
  frame_pul->GetYaxis()->CenterTitle(kTRUE);
  frame_pul->SetYTitle("#font[42]{Pull}");
  frame_pul->SetXTitle("#font[42]{}");
  frame_pul->Draw("EZ");
  TLine.SetLineColor(2);
  TLine.SetLineWidth(3);
  TLine.SetLineStyle(2);
  TLine.DrawLine(x_min_fit, 0, x_max_fit, 0);
  TLine.SetLineStyle(1);
  TLine.DrawLine(x_min_fit, -2, x_max_fit, -2);
  TLine.DrawLine(x_min_fit, 2, x_max_fit, 2);    
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".png");

  ////////////////
  //SIGNIFIGANCE//
  ////////////////
  smallBetween1 = .185;
  smallBetween2 = .05;
  smallBetween3 = .05;
  smallBetween4 = .185;
  //HybridCalculatorOriginal
  double significance_HCO = 0;
  RooExtendPdf bkg_ext_pdf("bkg_ext_pdf", "", bkg_PDF, nbkg);
  //PDF for CLs testing
  HybridCalculatorOriginal myHybridCalc(data, model, bkg_ext_pdf, 0, 0);
  myHybridCalc.UseNuisance(false);
  myHybridCalc.SetTestStatistic(1);
  myHybridCalc.SetNumberOfToys(2000);
  myHybridCalc.SetGenerateBinned(true);
  //default is 1000
  HybridResult* myHybridResult = myHybridCalc.GetHypoTest();
  //Plot
  HybridPlot* myHybridPlot = myHybridResult->GetPlot("myHybridPlot","#font[41]{Results with H.C.O.}", nbkg.getError());	
  
  cleg = "HybridCalculatorOriginal";
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  myHybridPlot->Draw();
  C1->Print(path + cleg + ".pdf");  
  C1->Print(path + cleg + ".png");  
  
  significance_HCO = myHybridResult->Significance();
  cout<<"significance_HCO = "<<significance_HCO<<endl;
    
  //ProfileLikelihoodCalculator
  RooArgSet POI_PLC_SIGNIFIGANCE(nsig);
  ProfileLikelihoodCalculator plc(data, model, POI_PLC_SIGNIFIGANCE);	
  RooArgSet nullparams;
  nullparams.addClone(nsig);
  ((RooRealVar *) (nullparams.first()))->setVal(0);
  plc.SetNullParameters(nullparams);
    
  HypoTestResult* plcResult = plc.GetHypoTest();
  double significance_PLC = 0;
  significance_PLC = (float) plcResult->Significance();
  
  cleg = "ProfileLikelihoodCalculator";
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  LikelihoodInterval* plInt_significance = plc.GetInterval();
  LikelihoodIntervalPlot lplot(plInt_significance);
  lplot.SetMaximum(10.);
  lplot.Draw();
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".root");
  
  ///////////////
  //UPPER LIMIT//
  ///////////////
  
  //Assigning the Parameter of Interest
  nsig.setRange(0, 500);
  /////////////////////nsig.setRange(-200, 200);
  //RESETTING OF RANGE OF nsig!
  RooArgSet POI_Bayes(nsig);
  RooArgSet POI_PLC_UL(nsig);

  RooWorkspace* ws = new RooWorkspace("w",true);
  ws->factory("SUM::pdf(s[0,15]*Uniform(x[0,1]),b[1,0,2]*Uniform(x))");
  ws->factory("Gaussian::prior_b(b,1,1)");
  ws->factory("Uniform::priorPOI(s)");
  
  RooAbsPdf* priorPOI = ws->pdf("priorPOI");
  
  //First Call to BayesianCalculator
  //Invoke Bayes calculator
  BayesianCalculator bcalc(data, model, POI_Bayes, *priorPOI);
  //90% CL
  bcalc.SetConfidenceLevel(0.90);
  //One-sided confidence interval
  bcalc.SetLeftSideTailFraction(0);							
  //bcalc.SetIntegrationType("TOYMC");
  //bcalc.SetNumIters(10000);
  SimpleInterval* interval = bcalc.GetInterval();
  cout << "90% CL interval: [ " << interval->LowerLimit() << " - " << interval->UpperLimit() << " ] or 90% CL limits\n";
  //Upper Limit
  double UL_Bayes = interval->UpperLimit();
  //Plot Likelihood and integral region. 
  RooPlot * plot = bcalc.GetPosteriorPlot();
  //bcalc.SetScanOfPosterior(50);
  cleg = "BayesianCalculator";	
  C1 = new TCanvas(cleg, cleg, 10, 10, 800, 800);
  gPad->SetLeftMargin(smallBetween1);
  gPad->SetRightMargin(smallBetween2);
  gPad->SetTopMargin(smallBetween3);
  gPad->SetBottomMargin(smallBetween4);
  
  //plot->SetMaximum(1.9);
  plot->Draw();
  C1->Print(path + cleg + ".pdf");
  C1->Print(path + cleg + ".png");
    
  //UPPER LIMIT
  ProfileLikelihoodCalculator plc_UL(data, model, POI_PLC_UL);	//Call to PLC with new range on nsig
  plc_UL.SetConfidenceLevel(0.90); 
  LikelihoodInterval* plInt_UL = plc.GetInterval();
  double UL_PLC = plInt_UL->UpperLimit(nsig);
  cout << " UL_Bayes " << UL_Bayes << " UL_PLC " << UL_PLC << endl;
  
}
