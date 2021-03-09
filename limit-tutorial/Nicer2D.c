void Nicer2D(TH2F *MyHisto, 
	     Double_t Size, 
	     Int_t TFONT, 
	     Int_t NDIV,
	     Double_t OffsetTitleX,
	     Double_t OffsetTitleY,
	     Double_t OffsetTitleZ)
{
  MyHisto->SetLabelSize(Size,"X");
  MyHisto->SetLabelSize(Size,"Y");
  MyHisto->SetLabelSize(Size,"Z");
  MyHisto->SetLabelFont(TFONT,"X");
  MyHisto->SetLabelFont(TFONT,"Y");
  MyHisto->SetLabelFont(TFONT,"Z");
  MyHisto->SetTitleSize(Size,"X");
  MyHisto->SetTitleSize(Size,"Y");
  MyHisto->SetTitleSize(Size,"Z");
  MyHisto->SetNdivisions(NDIV,"X");
  MyHisto->SetNdivisions(NDIV,"Y");
  MyHisto->SetNdivisions(NDIV,"Z");
  MyHisto->SetTitleOffset(OffsetTitleX,"X");
  MyHisto->SetTitleOffset(OffsetTitleY,"Y");
  MyHisto->SetTitleOffset(OffsetTitleZ,"Z");
  MyHisto->GetXaxis()->CenterTitle(kTRUE);
  MyHisto->GetYaxis()->CenterTitle(kTRUE);
  MyHisto->GetZaxis()->CenterTitle(kTRUE);
}


