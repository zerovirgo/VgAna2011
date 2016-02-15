#ifndef EvtWeightFunc_h
#define EvtWeightFunc_h
Double_t getPhoIDWeightFromPtVtx(float pt, float nVtx, float scEta, TH2* hEB, TH2* hEE){
	Int_t iPt(0);
	if ( pt > 20 && pt <= 30 ) iPt = 1 ;
	if ( pt > 30 && pt <= 40 ) iPt = 2 ;
	if ( pt > 40 && pt <= 50 ) iPt = 3 ;
	if ( pt > 50             ) iPt = 4 ;
//	if ( pt > 60             ) iPt = 5 ;
	Int_t iVtx(0);

	if ( nVtx <=   1                ) iVtx = 1 ;
	if ( nVtx ==   2  && nVtx ==  3 ) iVtx = 2 ;
	if ( nVtx ==   4  && nVtx ==  5 ) iVtx = 3 ;
	if ( nVtx ==   6  && nVtx ==  7 ) iVtx = 4 ;
	if ( nVtx ==   8  && nVtx ==  9 ) iVtx = 5 ;
	if ( nVtx >=  10                ) iVtx = 6 ;
	if ( iPt == 0 || iVtx == 0 ) return 1 ;
	if ( fabs(scEta) < 1.4442 ) return hEB->GetBinContent(iVtx,iPt);
	else if ( fabs(scEta) < 2.5    ) return hEE->GetBinContent(iVtx,iPt);
	else return 1. ;
}
Double_t getPhoIDWeightFromEta(float scEta){
	double weight[3] = { 1.038, 0.993 , 1.037 } ;
	int iSide = -1 ;
	if ( scEta < -1.566 && scEta > -2.5) {
		iSide = 0 ;
	}
	else if ( scEta > -1.4442 && scEta < 1.4442 ) {
		iSide = 1 ;
	}
	else if ( scEta > 1.566 && scEta < 2.5 ) {
		iSide = 2 ;
	}
	else return 1. ;
	return weight[iSide] ;
}
Double_t getEleTrgWeightFromTrgIndex(float scEta, Int_t trg_index){
	int iSide = -1 ;
	if (fabs(scEta) < 1.4442 ) iSide = 0 ;
	else iSide = 1 ;
	Double_t weight_ = 1. ;

	Int_t trg_array[] = { 29 , // HLT_Ele8_CaloIdL_CaloIsoVL_v*
	                      30   // HLT_Ele17_CaloIdL_CaloIsoVL_v*
							} ;
	const size_t nBin = sizeof(trg_array) / sizeof(trg_array[0]);
	double rhoeff[][2] = { // {EB} , {EE}
		{ 0.9970 , 0.9979 } ,
		{ 0.9916 , 0.9974 } 
	};
	for (size_t iBin = 0; iBin < nBin; iBin++){
		if (trg_index == trg_array[iBin] ){
			weight_ *= rhoeff[iSide][iBin];
		}
	}
	return weight_ ;
}
Double_t getEleIDWeightFromPt(float scEta, float pt){
	double xBin[] = { 20., 30., 40., 50., 60., 80., 100., 200. } ;
	int iSide = -1 ;
	if (fabs(scEta) < 1.4442 ) iSide = 0 ;
	else iSide = 1 ;
	Double_t weight_ = 1. ;
	const size_t nBin = sizeof(xBin) / sizeof(xBin[0]) - 1 ;
	double rhoeff[][7] = { // {EB} , {EE}
		{ 1.0044 , 1.0092 , 1.0236 , 1.0165 , 0.9913 , 0.9782 , 0.9865 } , 
		{ 1.0091 , 1.0333 , 1.0469 , 1.0440 , 1.0166 , 0.9733 , 0.9522 }   
				};
	for (size_t iBin = 0; iBin < nBin; iBin++){
		if ( pt > xBin[iBin] && pt < xBin[iBin+1] ){
			weight_ = rhoeff[iSide][iBin] ;
		}
	}
	return weight_ ;
}

Double_t getEleRecoWeightFromPt(float scEta, float pt){
	double xBin[] = { 20., 30., 40., 50., 60., 80., 100., 200. } ;
	int iSide = -1 ;
	size_t nBin = sizeof(xBin) / sizeof(xBin[0]) - 1 ;
	if (fabs(scEta) < 1.4442 ) iSide = 0 ;
	else iSide = 1 ;
	Double_t weight_ = 1. ;
	double rhoeff[][7] = { // {EB} , {EE}
		{ 0.9949 , 1.0002, 0.9981, 1.0005, 0.9985, 1.0008, 0.9947 } ,
		{ 0.9983 , 0.9967, 1.0020, 1.0033, 1.0045, 0.9997, 1.0017 }
				};
	for (size_t iBin = 0; iBin < nBin; iBin++){
		if ( pt > xBin[iBin] && pt < xBin[iBin+1] ){
			weight_ = rhoeff[iSide][iBin] ;
		}
	}
	return weight_ ;
}

Double_t getPhoBkgWeightFromPt(float scEta , float pt){
	double xBin[] = { 15. , 20., 25., 30., 35., 40., 60.,999.} ;
	size_t nBin = sizeof(xBin) / sizeof(xBin[0]) - 1 ;
	int iSide = -1 ;
	if (fabs(scEta) < 1.4442 ) iSide = 0 ;
	else iSide = 1 ;
	Double_t weight_ = 1. ;
	double bkgRatio[][7] = { // {EB} , {EE}
		{ 0.9906 , 1.3320 , 0.6230 , 1.5353 , 1.4073 , 0.6413 , 2.7310 } ,
		{ 1.3093 , 0.9576 , 0.8817 , 1.5353 , 1.8328 , 2.0501 , 1.1843 }
				};
	for (size_t iBin = 0; iBin < nBin; iBin++){
		if ( pt > xBin[iBin] && pt <= xBin[iBin+1] ){
			weight_ = bkgRatio[iSide][iBin] ;
		}
	}
	return weight_ ;
}

Double_t  getErrorFromHistoAndArray(TH1* hDist, vector<Double_t> error , Double_t  xBin[]){
	Double_t error_(0.) ;
//	if (hDist->GetNbinsX() != error.size() ){
//		cout << "Bin not consistent !!" << endl ;
//		return error_ ;
//	}
//	
	
	for (size_t i = 0; i < error.size(); i++){
		if ( hDist->Integral() > 0. ){
		double ratio = hDist->Integral(hDist->FindBin(xBin[i]) , hDist->FindBin(xBin[i+1]) -1 ) / hDist->Integral() ;
		error_ += pow(error[i] * ratio ,2);
		}
		else error_ += pow(error[i],2) ;
//		if ( error.size() == 3 ){
//			cout << "Bin in ratio = " << hDist->Integral(hDist->FindBin(xBin[i]) , hDist->FindBin(xBin[i+1]) -1 ) /
//			hDist->Integral() << " error2 = "<<  pow(error[i],2)  << endl ;
//		}
	}
	return sqrt(error_);
}
Double_t getJetWeightFromPhoEt(Float_t scEta, Float_t pt){
	double xBin[] = { 15 , 500 , 1000} ;
	int iSide = -1 ;
	size_t nBin = sizeof(xBin) / sizeof(xBin[0]) - 1 ;
	if (fabs(scEta) < 1.4442 ) iSide = 0 ;
	else iSide = 1 ;
	Double_t weight_ = 1. ;
	double rhoeff[][2] = { // {EB} , {EE}
		{ 1.256 , 1.256 } ,
		{ 1.442 , 1.442 }
				};
	for (size_t iBin = 0; iBin < nBin; iBin++){
		if ( pt > xBin[iBin] && pt < xBin[iBin+1] ){
			weight_ = rhoeff[iSide][iBin] ;
		}
	}
	return weight_ ;
}
Double_t getScaleFactorFromRatio(Float_t scEta){
	if ( fabs(scEta) < 1.4442 ) return 1.093 ;
	else if ( fabs(scEta) > 1.566 && fabs(scEta) < 2.5 ) return 0.76065 ;
	else return 1. ;
}

Double_t getScaleFactorFrom2DH(Double_t x_ , Double_t y_ , TH2* hin, bool print = false)
{
	Int_t xBins = hin->GetXaxis()->GetNbins() ; 
	Int_t yBins = hin->GetYaxis()->GetNbins() ;
	Double_t scale = 1. ;
	// fin value in the 2D map
	for (int iBin = 0; iBin < xBins; iBin++)
	{
		for (int jBin = 0; jBin < yBins; jBin++)
		{
			if (    x_ >  hin->GetXaxis()->GetBinLowEdge(iBin+1) && 
					x_ <= hin->GetXaxis()->GetBinUpEdge(iBin+1)  &&
					y_ >  hin->GetYaxis()->GetBinLowEdge(jBin+1) && 
					y_ <= hin->GetYaxis()->GetBinUpEdge(jBin+1)    )
			{
				scale = hin->GetBinContent(iBin+1,jBin+1);
				if (print)
				{
					cout << "x " << x_ << " ix "<< iBin << " y " << y_ << " iy " << jBin << " scale " << scale << endl ;
				}
				
			}
		}
	}

	// if x out of range ...
	if (x_ > hin->GetXaxis()->GetBinUpEdge(xBins))
	{
		for (int jBin = 0; jBin < yBins; jBin++) 
			if (    
					y_ >  hin->GetYaxis()->GetBinLowEdge(jBin+1) && 
					y_ <= hin->GetYaxis()->GetBinUpEdge(jBin+1)    )
				scale = hin->GetBinContent(xBins,jBin+1);
	}
	if (x_ < hin->GetXaxis()->GetBinLowEdge(1))
	{
		for (int jBin = 0; jBin < yBins; jBin++) 
			if (    
					y_ >  hin->GetYaxis()->GetBinLowEdge(jBin+1) && 
					y_ <= hin->GetYaxis()->GetBinUpEdge(jBin+1)    )
				scale = hin->GetBinContent(1,jBin+1);
	}

	// if y out of range ...
	if (y_ > hin->GetYaxis()->GetBinUpEdge(yBins))
	{
		for (int iBin = 0; iBin < xBins; iBin++)
			if (    x_ >  hin->GetXaxis()->GetBinLowEdge(iBin+1) && 
					x_ <= hin->GetXaxis()->GetBinUpEdge(iBin+1)  )
				scale = hin->GetBinContent(iBin+1,yBins);
	}
	if (y_ < hin->GetYaxis()->GetBinLowEdge(1))
	{
		for (int iBin = 0; iBin < xBins; iBin++)
			if (    x_ >  hin->GetXaxis()->GetBinLowEdge(iBin+1) && 
					x_ <= hin->GetXaxis()->GetBinUpEdge(iBin+1)  )
				scale = hin->GetBinContent(iBin+1,1);
	}

	// if both x and y out of range of the map
	if (y_ > hin->GetYaxis()->GetBinUpEdge(yBins) &&
			x_ > hin->GetXaxis()->GetBinUpEdge(xBins))
		scale = hin->GetBinContent(xBins,yBins);
	return scale ;
}

Double_t getScaleFactorFrom1DH(Double_t x_, TH1* hin)
{
	Int_t xBins = hin->GetXaxis()->GetNbins() ; 
	Double_t scale = 1. ;
	for (int iBin = 0; iBin < xBins; iBin++)
	{
		if (    x_ >  hin->GetXaxis()->GetBinLowEdge(iBin+1) && 
				x_ <= hin->GetXaxis()->GetBinUpEdge(iBin+1) )
		{
			scale = hin->GetBinContent(iBin+1);
		}
	}
		if (    x_ >  hin->GetXaxis()->GetBinLowEdge(xBins) )
		{
			scale = hin->GetBinContent(xBins);
		}
		if (    x_ <=  hin->GetXaxis()->GetBinLowEdge(1) )
		{
			scale = hin->GetBinContent(1);
		}
	return scale ;
}
Double_t getElectronOfflineScaleFactorWP85(double pt, double eta, bool uncert){
  double sf_electrons_RECO[2][4]={{1.002, 1.001, 1.009, 1.008}, {1.004, 0.984, 0.993, 0.944}};
  //double sf_electrons_ID_err[2][4]={{0.005, 0.005, 0.005, 0.005}, {0.005, 0.005, 0.007, 0.010}};
  
  double sf_electrons_IDISO[2][2]={{1.004, 1.043}, {0.988, 1.013}};
  //double sf_electrons_ISO_err[2][2]={{0.005, 0.005}, {0.005, 0.006}};

  //double sf_electrons_TRG[1][2]={{0.999, 0.999}};
  //double sf_electrons_TRG_err[1][2]={{0.01, 0.01}};

  // defaults  
  double w(1);
  //double w_err(1); 

  int idx=-1;
  int idy=-1;

  int isox=-1;
  int isoy=-1;
    
  // -------- pt ------------------------------
  if (pt<=50) {
    idx  =0;
    isox =0;
   } else if (pt>50) {
    idx  =1;
    isox =1;
  }
  
  // -------- eta -----------------------------
  if (eta<=0.8) {
    idy  =0;
    isoy =0;
  } else if (TMath::Abs(eta)>0.8 && TMath::Abs(eta)<=1.44) {   
    idy  =1;
    isoy =0;
  } else if (TMath::Abs(eta)>1.57 && TMath::Abs(eta)<=1.6) {
    idy  =2;
    isoy =0;
  } else if (TMath::Abs(eta)>1.6 && TMath::Abs(eta)<=2.0) {
    idy  =2;
    isoy =1;
  } else if (TMath::Abs(eta)>2.0 && TMath::Abs(eta)<=2.5) {
    idy =3;
    isoy =1;
  }

  // --------------- final SF computation (and error) ------------- 
  if (idx>=0 && idy>=0 && isox>=0 && isoy>=0){
    w = sf_electrons_RECO[idx][idy]*sf_electrons_IDISO[isox][isoy];
    //w_err = TMath::Sqrt( TMath::Power(sf_electrons_ID_err[idx][idy]*sf_electrons_ISO[isox][isoy]*sf_electrons_TRG[0][trgy],2)+ 
    //			 TMath::Power(sf_electrons_ISO_err[isox][isoy]*sf_electrons_ID[idx][idy]*sf_electrons_TRG[0][trgy],2)+
    //			 TMath::Power(sf_electrons_TRG_err[0][trgy]*sf_electrons_ID[idx][idy]*sf_electrons_ISO[isox][isoy],2));// errors to be propagated
  }
  return w;
}


#endif
