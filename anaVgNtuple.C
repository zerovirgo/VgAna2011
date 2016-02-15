//////////////////////////////////////////////////////////////////////////////////
// Use VgNtuple vgamma_42x-v11 to perform analysis
//
// v3-03 12:09 29/Nov/2011 Written by Kuan-Hsin CHEN
//
// v2-02
// Fill histogram
// can use for flow histograms
// v2-03
// Run on MC 
// check single electron trigger leg
// v2-04
// Can use Template and store tree
// v2-05
// Do Ratio method form Anthony
// Draw Flow chart
// check option
// v3-01
// to VgKitV11 
// add nJets Histogram
// add HLT WP function
// v2-02
// fix some bugs
// v2-03
// fix pile-up rewieghting problem
// v2-05
// ElectronID eleReco phoID efficiency 
// v2-06
// Add Energy Scale , nPU shift 
// v2-07
// include pile up stury, nVtx preparing
// v2-08
// update
// v3-01 update for VgKitV14
// v3-02 
// S6 PU reweighting for Fall11 Samples
// Can read multiple input in 1 run
// 2D weight function to Header
// v3-03 
// 3D-Reweighting 
// New Electron Correction from Yurii M.
// Photon timming histograms
// v3-04
// Energy Smering for EM object in MC
// - Please notice that : func of Muon channel analysis and 
//   ShowEvent and relative functions are too old and not useful
// v3-05
// Add more Z information
//
//
// contect : zerovirgo@gamil.com , kschen@cern.ch
//
//////////////////////////////////////////////////////////////////////////////////

#define anaVgNtuple_cxx
#include "anaVgNtuple.h"
#include "Header/EvtWeightFunc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TH2D.h>
#include <TString.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <time.h>
#include <utility>

inline Bool_t anaVgNtuple::CleanEvent(){
		nGoodVtx = 0 ;
		for (int iV = 0; iV < nVtx; iV++){
		if (vtxNDF[iV] >= 4 && fabs(vtx[iV][2]) <= 24 && fabs(vtxD0[iV]) <= 2)  nGoodVtx++ ;
		}
		if (nGoodVtx<1) return false ;
		if (IsVtxGood ==0 ) return false ;
		if (IsTracksGood !=1 ) return false ;
		return true ;
}
//Bool_t anaVgNtuple::GenZllCut()
//{
//	vector<Int_t> lep ;
//	Int_t         iPho    =  -1 ;
//	Float_t       gammaPt = 0.0 ;
//	for (int iMC = 0; iMC < nMC; iMC++){
//		if ( fabs(mcPID[iMC]) == 11 ){
//			if ( fabs(mcMomPID[iMC]) == 11 && mcGMomPID[iMC] == 23 ) lep.push_back(iMC);
//			if ( mcMomPID[iMC] == 23 )                               lep.push_back(iMC);
//		}
//		if ( mcPID[iMC] == 22 ){
//			if ( fabs(mcMomPID[iMC]) == 22 && mcGMomPID[iMC] == 11 && mcPt[iMC] > gammaPt  )
//			{
//				lep.push_back(iMC);
//				gammaPt = mcPt[iMC];
//			}
//			if ( fabs(mcMomPID[iMC]) == 22 && fabs(mcGMomPID[iMC]) < 7 && mcPt[iMC] > gammaPt  )
//			{
//				lep.push_back(iMC);
//				gammaPt = mcPt[iMC];
//			}
//		}
//	}
//}

Int_t anaVgNtuple::phoCleanLevel(Int_t lep1_index, Int_t lep2_index ,Int_t pho_index){
	if ( doRmVJetsPho &&
			( ProcessTag.Contains("ZJet") || ProcessTag.Contains("WJet") ) )
	{
		if ( IsISR(pho_index) || IsFSR(pho_index)         ) 
		{
			return 0 ;
		}
	}
		// Photon kinamatic cut
		if (        phoEt_[pho_index]    <  PhoEtCut     ) return 0 ;
		if ( fabs(phoSCEta[pho_index])   >  2.5          ) return 0 ;
		if ( fabs(phoSCEta[pho_index])   >  1.4442 
		&&   fabs(phoSCEta[pho_index])   <  1.566        ) return 0 ;
		if ( IsSpike(pho_index,"pho")                    ) return 0 ;
		if ( phohasPixelSeed[pho_index]  == 1            ) return 0 ;
		// Pixel Seed Veto
		if ( phoHoverE[pho_index]        > PhoHoverPreCut) return 1 ;
		if (doEle){
		if ( phohasPixelSeed[pho_index]  == 1            ) return 1 ;
		// DeltaR cut
			if ( DeltaR(eleEta[lep1_index],elePhi[lep1_index],phoEta[pho_index],phoPhi[pho_index]) < 0.7 ) return 2 ;
			if ( DeltaR(eleEta[lep2_index],elePhi[lep2_index],phoEta[pho_index],phoPhi[pho_index]) < 0.7 ) return 2 ;
		}
		if (doMu){
			if ( DeltaR(muEta[lep1_index],muPhi[lep1_index],phoEta[pho_index],phoPhi[pho_index])   < 0.7 ) return 2 ;
			if ( DeltaR(muEta[lep2_index],muPhi[lep2_index],phoEta[pho_index],phoPhi[pho_index])   < 0.7 ) return 2 ;
		}
		// full photon clean selection
		return 3 ;
}
Int_t anaVgNtuple::PhoSelection(
		Int_t lep1_index , 
		Int_t lep2_index  
		){
// function to select photon , return photon_index, if no photon candidate found, return -1
	Int_t pho_index = -1 ;
	float leadingPt = PhoEtCut ;
	Int_t nPassPho = 0;
	for (int iPho = 0; iPho < nPho; iPho++){
		if ( doRmVJetsPho &&
				( ProcessTag.Contains("ZJet") || ProcessTag.Contains("WJet") ) )
		{
			if ( IsISR(iPho) || IsFSR(iPho)         ) 
			{
				//if  (mcPhoMatcher(iPho)) cout << "matched" << endl ;
				//else                     cout << "empty"   << endl ;
				//if  (IsISR(iPho)) cout << "ISR" << endl ;
				//if  (IsFSR(iPho)) cout << "FSR" << endl ;
				continue ;
			}
			//else
			//{
			//	cout << "Not ISR nor FSR" << endl ;
			//}

		}
		if ( phoCleanLevel(lep1_index,
	       	lep2_index,iPho)     < 3                ) continue ;
		if ( doPhoStandard   &&  !phoPassID(iPho,3) ) continue ;
		if ( doRatioA        &&  !IsFakeable(iPho)  ) continue ;
		if ( doRatioA                               ) {
			FillHisto(lep1_index,lep2_index,iPho);
		}
		nPassPho++;
		if ( leadingPt <  phoEt_[iPho]){
			leadingPt = phoEt_[iPho] ;
			pho_index = iPho ;
		}
	}
	if (nPassPho>= 2 && isData && doShowList) cout << "event " << event << " has " << nPassPho << " GoodPhotons " << endl ;
	return pho_index ;
}
Int_t anaVgNtuple::PhoSelectLevel(
		Int_t lep1_index , 
		Int_t lep2_index  
		){
	Int_t phoLV(0);
	for (int iPho = 0; iPho < nPho; iPho++){
		Int_t preLevel = 0 ;
		// Pixel Seed Veto
		if ( phohasPixelSeed[iPho]  == 1            ) continue;
		preLevel++ ; phoLV = ( phoLV < preLevel ) ? preLevel : phoLV ;
		// Photon kinematic
		if (        phoEt_[iPho]    <  PhoEtCut     ) continue;
		if ( fabs(phoSCEta[iPho])   >  2.5          ) continue;
		if ( fabs(phoSCEta[iPho])   >  1.4442 
		&&   fabs(phoSCEta[iPho])   <  1.566        ) continue;
		if ( phoHoverE[iPho]        > PhoHoverPreCut) continue;
		if ( IsSpike(iPho,"pho")                    ) continue ;
		preLevel++ ; phoLV = ( phoLV < preLevel ) ? preLevel : phoLV ;
		// DeltaR
		if (doEle){
			if ( DeltaR(eleEta[lep1_index],elePhi[lep1_index],phoEta[iPho],phoPhi[iPho]) < 0.7 ) continue ;
			if ( DeltaR(eleEta[lep2_index],elePhi[lep2_index],phoEta[iPho],phoPhi[iPho]) < 0.7 ) continue ;
		}
		if (doMu){
			if ( DeltaR(muEta[lep1_index],muPhi[lep1_index],phoEta[iPho],phoPhi[iPho])   < 0.7 ) continue ;
			if ( DeltaR(muEta[lep2_index],muPhi[lep2_index],phoEta[iPho],phoPhi[iPho])   < 0.7 ) continue ;
		}
		preLevel++ ; phoLV = ( phoLV < preLevel ) ? preLevel : phoLV ;
		// PhotonID
		preLevel+= phoPassIDLevel(iPho,3);
		phoLV = ( phoLV < preLevel ) ? preLevel : phoLV ;
	}
	return phoLV ;
}

Int_t  anaVgNtuple::getNPassJets(){
	Int_t _nPassJets = 0  ; 
//	cout << "Number of candidate Jets " << nJet <<  endl ;
	for (int iJet=0; iJet < nJet; iJet++){
		if (jetNeutralHadronEnergyFraction[iJet] > 0.99) continue;
		if (jetNeutralEmEnergyFraction[iJet] > 0.99) continue;
		if (jetNConstituents[iJet] < 2) continue;
		if (fabs(jetEta[iJet]) < 2.5 &&
				jetChargedHadronEnergyFraction[iJet] <= 0) continue;
		if (fabs(jetEta[iJet]) < 2.5 &&
				jetChargedHadronMultiplicity[iJet] == 0) continue;
		if (fabs(jetEta[iJet]) < 2.5 &&
				jetChargedHadronEnergyFraction[iJet] > 1.99) continue;
		bool overlapFlag = false ;
		for (int iPho = 0; iPho < nPho; iPho++){
			if ( DeltaR(phoEta[iPho],phoPhi[iPho],jetEta[iJet],jetPhi[iJet]) < 0.5 ) overlapFlag = true ;
		}
		for (int iEle = 0; iEle < nEle; iEle++){
			if ( DeltaR(eleEta[iEle],elePhi[iEle],jetEta[iJet],jetPhi[iJet]) < 0.5 ) overlapFlag = true ;
		}
		for (int iMu = 0; iMu < nMu; iMu++){
			if ( DeltaR(muEta[iMu],muPhi[iMu],jetEta[iJet],jetPhi[iJet]) < 0.5 ) overlapFlag = true ;
		}
		if (overlapFlag) continue ;
		_nPassJets++;
	}
//	cout << "Number of Passed Jets " << _nPassJets <<  endl ;
	if (doHistoPassPho) histo[19]->Fill(_nPassJets,EvtWeight);
	return _nPassJets ;
}

Bool_t anaVgNtuple::IsSpike(Int_t index , TString option){
	Bool_t isSpike = false ;
	if (option == "ele"){
		if (fabs(eleSCEta[index]) < 1.5 && eleSigmaIEtaIEta[index] < 0.001) isSpike = true ;
		if (fabs(eleSCEta[index]) < 1.5 && eleSigmaIPhiIPhi[index] < 0.001) isSpike = true ;
//		if (eleSeverity[index] == 3 || eleSeverity[index] == 4 || eleSeverity[index] == 5) isSpike = true ;
	}
	if (option == "pho"){
		if (fabs(phoSCEta[index]) < 1.5 && phoSigmaIEtaIEta[index] < 0.001) isSpike = true ;
		if (fabs(phoSCEta[index]) < 1.5 && phoSigmaIPhiIPhi[index] < 0.001) isSpike = true ;
//		if (phoSeverity[index] == 3 || phoSeverity[index] == 4 || phoSeverity[index] == 5) isSpike = true ;
	}
	return isSpike ;
}

void anaVgNtuple::ShowEvent(
	Long64_t entry, 
	Int_t index_ele1, 
	Int_t index_ele2,
	Int_t index_pho
){
	if (!fChain) return ;
	fChain->GetEntry(entry);

	// electron energy and pt correction
	doEMCorrection   = isData ;
	if (doEMCorrection){
		eleEn_ = newEleEn();
		elePt_ = newElePt(eleEn_) ;
	} else {
		eleEn_ = vector<Float_t>(eleEn, eleEn + nEle);
		elePt_ = vector<Float_t>(elePt, elePt + nEle);
	}
	// photon energy and et correction
	if (doEMCorrection){
		phoE_  = newPhoE();
		phoEt_ = newPhoEt(phoE_) ;
	} else {
		phoEt_ = vector<Float_t>(phoEt, phoEt + nPho);
		phoE_ =  vector<Float_t>(phoE , phoE  + nPho);
	}
	printEvent(index_ele1,index_ele2,index_pho);
}

void anaVgNtuple::addScale(TString  objName){
	if ( objName == "Ele" || objName == "ele")
	{
		unsigned int nnEle = nEle ;
		if ( eleEn_.size() != nnEle ) return ;
		if ( eleScaleSysEB == 0 && eleScaleSysEE == 0 ) return ;
		for (unsigned int iEle = 0; iEle < eleEn_.size() ; iEle++ )
		{
			if (fabs(eleSCEta[iEle]) < 1.4442)
			{
				eleEn_[iEle] *= 1.0 + eleScaleSysEB ;
				elePt_[iEle] *= 1.0 + eleScaleSysEB ;
			}
			else
			{
				eleEn_[iEle] *= 1.0 + eleScaleSysEE ;
				elePt_[iEle] *= 1.0 + eleScaleSysEE ;
			}
		}
	}

	if ( objName == "Pho" || objName == "pho")
	{
		unsigned int nnPho = nPho ;
		if ( phoE_.size() != nnPho ) return ;
		if ( phoScaleSysEB == 0 && phoScaleSysEE == 0 ) return ;
		for (unsigned int iPho = 0; iPho < phoE_.size() ; iPho++ )
		{
			if (fabs(phoSCEta[iPho]) < 1.4442)
			{
				phoE_ [iPho] *= 1.0 + phoScaleSysEB ;
				phoEt_[iPho] *= 1.0 + phoScaleSysEB ;
			}
			else
			{
				phoE_ [iPho] *= 1.0 + phoScaleSysEE ;
				phoEt_[iPho] *= 1.0 + phoScaleSysEE ;
			}
		}
	}
	return ;
}

void anaVgNtuple::addResolution(TString  objName){
	if ( objName == "Ele" || objName == "ele")
	{
		unsigned int nnEle = nEle ;
		if ( eleEn_.size() != nnEle ) return ;
		if ( eleResoSysEB == 0 && eleResoSysEE == 0 ) return ;
		for (unsigned int iEle = 0; iEle < eleEn_.size() ; iEle++ )
		{
			float myScale = 1.0 ;
			if (fabs(eleSCEta[iEle]) < 1.4442)
			{
				myScale = myRandom->Gaus(1.,eleResoSysEB);
				eleEn_[iEle] *= myScale  ;
				elePt_[iEle] *= myScale  ;
			}
			else
			{
				myScale = myRandom->Gaus(1.,eleResoSysEE);
				eleEn_[iEle] *= myScale ;
				elePt_[iEle] *= myScale ;
			}
		}
	}

	if ( objName == "Pho" || objName == "pho")
	{
		unsigned int nnPho = nPho ;
		if ( phoE_.size() != nnPho ) return ;
		if ( phoResoSysEB == 0 && phoResoSysEE == 0 ) return ;
		for (unsigned int iPho = 0; iPho < phoE_.size() ; iPho++ )
		{
			float myScale = 1.0 ;
			if (fabs(phoSCEta[iPho]) < 1.4442)
			{
				myScale = myRandom->Gaus(1.,phoResoSysEB);
				phoE_ [iPho] *= myScale  ;
				phoEt_[iPho] *= myScale  ;
			}
			else
			{
				myScale = myRandom->Gaus(1.,phoResoSysEE);
				phoE_ [iPho] *= myScale ;
				phoEt_[iPho] *= myScale ;
			}
		}
	}
	return ;
}

vector<Double_t> anaVgNtuple::generate_flat10_weights(TH1D* data_npu_estimated){
	// see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
	const double npu_probs[25] = {0.0698146584, 0.0698146584, 0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584 /* <-- 10*/,
		0.0630151648,0.0526654164,0.0402754482,0.0292988928,0.0194384503,0.0122016783,0.007207042,0.004003637,0.0020278322,
		0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 /* <-- 24 */};
	vector<double> result(25);
	double s = 0.0;
	for(int npu=0; npu<25; ++npu){
		double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
		result[npu] = npu_estimated / npu_probs[npu];
		s += npu_estimated;
	}
	// normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
	for(int npu=0; npu<25; ++npu){
		result[npu] /= s;
	}
	return result;
}

Double_t anaVgNtuple::PU3D_Weight(Int_t pm1 , Int_t p0 , Int_t p1 , TH3* hin){
	pm1 = TMath::Min(49,pm1);
	p0  = TMath::Min(49,p0 );
	p1  = TMath::Min(49,p1 );
	Double_t Weight = 1. ;
	Weight = hin->GetBinContent(pm1,p0,p1);
	return Weight ;
}

vector<Double_t> anaVgNtuple::generate_S3_weights(TH1D* data_npu_estimated){
	// see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
	const double npu_probs[50] = {0.0723615, 0.0755922, 0.0659539,
		0.0713979, 0.0657199, 0.0716397, 0.0778073, 0.0731223, 0.0664508, 0.073128,
		0.0667021, 0.0592924, 0.0501524, 0.033843, 0.0291496, 0.0192694, 0.0126059,
		0.00865009, 0.00222387, 0.00271679, 0.00197434, 0.000246604, 0.000001, 0.000001, 0.000001 ,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0} ;
	vector<double> result(50);
	double s = 0.0;
	for(int npu=0; npu<50; ++npu){
		double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
		result[npu] = npu_estimated / npu_probs[npu];
		s += npu_estimated;
	}
	// normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
	double sum = 0. ;
	for(int npu=0; npu<50; ++npu){
		result[npu] /= s;
		sum += result[npu] ;
	}
	return result;
}

vector<Double_t> anaVgNtuple::generate_S6_weights(TH1D* data_npu_estimated)
{
	const double npu_probs[50] =
	// Sample In-time
	{ 0.00892372, 0.0192435, 0.0315688, 0.0436898, 0.0536932, 0.0601263, 0.0626846, 0.0631723, 0.0610671, 0.0579129, 0.053774, 0.0500593, 0.046286, 0.0428236, 0.0396703, 0.0367312, 0.033888, 0.0311444, 0.0285041, 0.0260051, 0.0235356, 0.0210846, 0.0187871, 0.0164829, 0.0143669, 0.0229481, 0.00894086, 0.00750105, 0.00617546, 0.00511003, 0.00409912, 0.0032948, 0.00264224, 0.00207933, 0.00162502, 0.00123888, 0.000948434, 0.000731111, 0.000551436, 0.000413313, 0.000303918, 0.000222048, 0.000164269, 0.000119871, 8.64774e-05, 6.25912e-05, 4.54835e-05, 3.17211e-05, 2.35341e-05, 1.39972e-05 };
	vector<double> result(50);
	double s = 0.0;
	for(int npu=0; npu<50; ++npu){
		double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
		result[npu] = npu_estimated / npu_probs[npu];
		s += npu_estimated;
	}
	// normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
	for(int npu=0; npu<50; ++npu){
		result[npu] /= s;
	}
	return result;
}

vector<Double_t> anaVgNtuple::generate_S4_weights(TH1D* data_npu_estimated){
	const double npu_probs[4][50] = { 
	// Sample In-time
	{ 0.116601, 0.0643305, 0.0688894, 0.0690878, 0.0690538, 0.0682825, 0.0667944, 0.0653175, 0.0619362, 0.0582307, 0.0528067, 0.0472947, 0.0409742, 0.03453, 0.0283149, 0.0226965, 0.0177862, 0.0135091, 0.0100665, 0.00735822, 0.00524764, 0.0037156, 0.00251507, 0.00167913, 0.00111451, 0.000726347, 0.000454571, 0.000294632, 0.000176598, 0.000108227, 6.675e-05, 4.00841e-05, 2.19752e-05, 1.34182e-05, 8.84125e-06, 4.06527e-06, 2.87127e-06, 1.73414e-06, 7.95997e-07, 6.53854e-07, 3.12713e-07, 8.52854e-08, 1.42142e-07, 1.0e-10, 1.0e-10, 1.0e-10, 2.84285e-08, 1.0e-10, 1.0e-10, 1.0e-10 } ,
	// Sample Average
	{ 0.0910746, 0.0662995, 0.0691138, 0.0697803, 0.0697461, 0.069043, 0.0683954, 0.0681454, 0.0677981, 0.0657211, 0.0614862, 0.0550286, 0.0468644, 0.0378157, 0.0291675, 0.0215018, 0.0153515, 0.0105599, 0.00695455, 0.00440207, 0.00264464, 0.00147601, 0.000778612, 0.000399654, 0.000202144, 0.000110666, 6.01797e-05, 3.68981e-05, 2.14339e-05, 1.09159e-05, 5.31583e-06, 2.38786e-06, 1.10865e-06, 4.5483e-07, 8.52807e-08, 5.68538e-08, 2.84269e-08, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 1.0e-10, 0, 0, 0, 0, 0, 0, 0 } ,
	// True In-time
	{ 0.14551, 0.0644453, 0.0696412, 0.0700311, 0.0694257, 0.0685655, 0.0670929, 0.0646049, 0.0609383, 0.0564597, 0.0508014, 0.0445226, 0.0378796, 0.0314746, 0.0254139, 0.0200091, 0.0154191, 0.0116242, 0.00846857, 0.00614328, 0.00426355, 0.00300632, 0.00203485, 0.00133045, 0.000893794, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } ,
	// True Average
	{ 0.104109, 0.0703573, 0.0698445, 0.0698254, 0.0697054, 0.0697907, 0.0696751, 0.0694486, 0.0680332, 0.0651044, 0.0598036, 0.0527395, 0.0439513, 0.0352202, 0.0266714, 0.019411, 0.0133974, 0.00898536, 0.0057516, 0.00351493, 0.00212087, 0.00122891, 0.00070592, 0.000384744, 0.000219377, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 } };


//const double npu_probs[4][25] = { // from DYJetsToLL_TuneZ2_M50_Madgraph_Summer11
//		{0.116782, 0.0644221, 0.0690056, 0.0696675, 0.0691216, 0.0683445, 0.0668447, 0.0653128, 0.0620426, 0.0584028, 0.0529344, 0.0473904, 0.0410362, 0.0345557, 0.0283562, 0.022728, 0.0177843, 0.0135265, 0.0100753, 0.0073685, 0.00525266, 0.00371981, 0.00252111, 0.00168566, 0.00111913 } ;
//	//const double npu_probs[25] = { 0.104109, 0.0703573, 0.0698445, 0.0698254, 0.0697054, 0.0697907, 0.0696751, 0.0694486, 0.0680332, 0.0651044, 0.0598036, 0.0527395, 0.0439513, 0.0352202, 0.0266714, 0.019411, 0.0133974, 0.00898536, 0.0057516, 0.00351493, 0.00212087, 0.00122891, 0.00070592, 0.000384744, 0.000219377 },
//	} ;
	vector<double> result(25);
	double s = 0.0;
	for(int npu=0; npu<25; ++npu){
		double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
		result[npu] = npu_estimated / npu_probs[PUOption][npu];
		s += npu_estimated;
	}
	// normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
	for(int npu=0; npu<25; ++npu){
		result[npu] /= s;
	}
	return result;
}

void anaVgNtuple::printMC(Int_t pdgID){
	if (isData) return ;
	bool matched = false ;
	if (pdgID != 0 )
	{
		for (int iMC = 0; iMC < nMC; iMC++) if (abs(mcPID[iMC]) == abs(pdgID)) matched = true ;
	}
	else
	{
		matched = true ;
	}
	if (!matched) return ;
	cout << "|  pdgID |   MomID |  GMomID |      Pt |     Eta |     Phi |     Mass | " << endl ;
	for (int iMC = 0; iMC < nMC; iMC++)
	{
		if ( pdgID == 0 || abs(mcPID[iMC]) == abs(pdgID) ){
			cout << fixed ;
			cout << "|" ;
			cout << setw(7)                    << mcPID[iMC]      << " |" ;
			cout << setw(8)                    << mcMomPID[iMC]   << " |" ;
			cout << setw(8)                    << mcGMomPID[iMC]  << " |" ;
			cout << setw(8) << setprecision(3) << mcPt[iMC]       << " |" ;
			cout << setw(8) << setprecision(3) << mcEta[iMC]      << " |" ;
			cout << setw(8) << setprecision(3) << mcPhi[iMC]      << " |" ;
			cout << setw(8) << setprecision(3) << mcMass[iMC]     << " |" ;
			cout << endl ;
		}
	}
	return ;
}

void anaVgNtuple::printEvent(int index_e1, int index_e2, int index_pho){
	float cIso1 = 100. ;
	float ecalIso1(0.);
	float trkIso1(0.) ;
	float hcalIso1(0.) ;
	float cIso2 = 100. ;
	float ecalIso2(0.);
	float trkIso2(0.) ;
	float hcalIso2(0.) ;
//	cout << "PassedZeeGamma: run " << run << " event " << event << " rho " << rho << endl ;
	cout << "+++++++++ Show Event Informations +++++++++++++++++++" << endl ;
	cout << "check Information : run " << run << " event " << event << " rho " << rho << " MET " << pfMET << endl ;
	cout << "nGVtx = " << nGoodVtx << endl ;
	cout << "IsVtxGood = " << IsVtxGood << endl ;
	cout << "IsTracksGood = " << IsTracksGood << endl ;
	int lepTrg_index = -1 ;
	cout << "PassHLT      = " << passHLT(lepTrg_index) ;
	//	if (HLTIndex[121] >= 0){ // HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1
	//		if ( HLT[HLTIndex[121]] != 1 ) cout << "Event not passHLT v1" << endl ;
	//		lepTrg_index = 20 ;
	//	}
	//	else
	//	if (HLTIndex[122] >= 0){ // HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2
	//		if ( HLT[HLTIndex[122]] != 1 ) cout << "Event not passHLT v2" << endl ;
	//		lepTrg_index = 21 ;
	//	}
	//	else
	//	if (HLTIndex[123] >= 0){ // HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3
	//		if ( HLT[HLTIndex[123]] != 1 ) cout << "Event not passHLT v3" << endl ;
	//		lepTrg_index = 22 ;
	//	}
	//	else
	//	if (HLTIndex[124] >= 0){ // HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4
	//		if ( HLT[HLTIndex[124]] != 1 ) cout << "Event not passHLT v4" << endl ;
	//		lepTrg_index = 23 ;
	//	}
	//	else
	//	if (HLTIndex[125] >= 0){ // HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5
	//		if ( HLT[HLTIndex[125]] != 1 ) cout << "Event not passHLT v5" << endl ;
	//		lepTrg_index = 24 ;
	//	}
	//	else
	//	if (HLTIndex[192] >= 0){ // HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6
	//		if ( HLT[HLTIndex[192]] != 1 ) cout << "Event not passHLT v6" << endl ;
	//		lepTrg_index = 28 ;
	//	}
		cout << "ele1Trg = " << eleTrg[index_e1][lepTrg_index] << endl ;
		cout << "ele2Trg = " << eleTrg[index_e2][lepTrg_index] << endl ;
	// calculate cIso
	if ( fabs(eleSCEta[index_e1]) < 1.4442 ){
		ecalIso1 = TMath::Max(eleIsoEcalDR03[index_e1]-1.,0.) ; 
		trkIso1  = eleIsoTrkDR03[index_e1] ;
		hcalIso1 = eleIsoHcalSolidDR03[index_e1] ;
		cIso1 = ( eleIsoTrkDR03[index_e1] + TMath::Max(eleIsoEcalDR03[index_e1]-1.,0.) + eleIsoHcalSolidDR03[index_e1] - rho*TMath::Pi()*0.3*0.3) / elePt_[index_e1];
	}
	else {
		cIso1 = ( eleIsoTrkDR03[index_e1] + eleIsoEcalDR03[index_e1] + eleIsoHcalSolidDR03[index_e1] - rho*TMath::Pi()*0.3*0.3) / elePt_[index_e1];
		ecalIso1 = eleIsoEcalDR03[index_e1] ; 
		trkIso1  = eleIsoTrkDR03[index_e1] ;
		hcalIso1 = eleIsoHcalSolidDR03[index_e1] ;
	}
	if ( fabs(eleSCEta[index_e2]) < 1.4442 ){
		ecalIso2 = TMath::Max(eleIsoEcalDR03[index_e2]-1.,0.) ; 
		trkIso2  = eleIsoTrkDR03[index_e2] ;
		hcalIso2 = eleIsoHcalSolidDR03[index_e2] ;
		cIso2 = ( eleIsoTrkDR03[index_e2] + TMath::Max(eleIsoEcalDR03[index_e2]-1.,0.) + eleIsoHcalSolidDR03[index_e2] - rho*TMath::Pi()*0.3*0.3) / elePt_[index_e2];
	}
	else {
		cIso2 = ( eleIsoTrkDR03[index_e2] + eleIsoEcalDR03[index_e2] + eleIsoHcalSolidDR03[index_e2] - rho*TMath::Pi()*0.3*0.3) / elePt_[index_e2];
		ecalIso2 = eleIsoEcalDR03[index_e2] ; 
		trkIso2  = eleIsoTrkDR03[index_e2] ;
		hcalIso2 = eleIsoHcalSolidDR03[index_e2] ;
	}


//	cout << "PassedZeeGamma: run " << run << " event " << event << " rho " << rho << endl ;
	cout << "electron :index,pt0,pt1,eta,phi,En,trkIso,ecalIso,hcalIso,cIso,scEta\nsceta,sieie,hoe  " << endl ;
	cout << index_e1 << " , " << elePt[index_e1] << " , " << elePt_[index_e1] << " , " << eleEta[index_e1] ;
	cout << " , " << elePhi[index_e1] << " , " << elePt_[index_e1]/elePt[index_e1]*eleEn[index_e1] << " , " << trkIso1 ;
	cout << " , " << ecalIso1 << " , " << hcalIso1 << " , " << cIso1 << " , " << eleSCEta[index_e1] << endl ;
	cout << eleSCEta[index_e1] << " , " << eleSigmaIEtaIEta[index_e1] << " , " << eleHoverE[index_e1] << endl ;
	cout << endl ;
	cout << index_e2 << " , " << elePt[index_e2] << " , " << elePt_[index_e2] << " , " << eleEta[index_e2] ;
	cout << " , " << elePhi[index_e2] << " , " << elePt_[index_e2]/elePt[index_e2]*eleEn[index_e2] << " , " << trkIso2 ;
	cout << " , " << ecalIso2 << " , " << hcalIso2 << " , " << cIso2 << " , " << eleSCEta[index_e2] << endl ;
	cout << eleSCEta[index_e2] << " , " << eleSigmaIEtaIEta[index_e2] << " , " << eleHoverE[index_e2] << endl ;
	cout << endl ;

	cout << "electron :index,dPhiIn,dEtaIn,mhit,Dist,Dcot,charge  " << endl ;
	cout << index_e1 << " , " << eledPhiAtVtx[index_e1] << " , " << eledEtaAtVtx[index_e1] << " , " << eleConvMissinghit[index_e1] ;
	cout << " , " << eleConvDist[index_e1] << " , " << eleConvDcot[index_e1] << " , " << eleCharge[index_e1] << endl ;
	cout << endl ;
	cout << index_e2 << " , " << eledPhiAtVtx[index_e2] << " , " << eledEtaAtVtx[index_e2] << " , " << eleConvMissinghit[index_e2] ;
	cout << " , " << eleConvDist[index_e2] << " , " << eleConvDcot[index_e2] << " , " << eleCharge[index_e1] << endl ;
	cout << endl ;

	cout << "Electron Pass " << endl ;
	cout << "TrackIso " << eleIsoTrkDR03[index_e1] << " EcalIso " << eleIsoEcalDR03[index_e1] << " HcalSolidIso ";
	cout << eleIsoHcalSolidDR03[index_e1] << " passID = " << elePassID(index_e1,3) << endl ;
	cout << "TrackIso " << eleIsoTrkDR03[index_e2] << " EcalIso " << eleIsoEcalDR03[index_e2] << " HcalSolidIso ";
	cout << eleIsoHcalSolidDR03[index_e2] << " passID = " << elePassID(index_e2,3) <<  endl ;
	TLorentzVector ve1, ve2, vZ ;
	ve1.SetPtEtaPhiM(elePt_[index_e1],eleEta[index_e1],elePhi[index_e1],0.000511);
	ve2.SetPtEtaPhiM(elePt_[index_e2],eleEta[index_e2],elePhi[index_e2],0.000511);
	vZ = ve1 + ve2 ;
	
	cout << "ZMass = " << vZ.M() << endl ;

	cout << "Photon : index, pt, corrpt , eta, phi, dr1, dr2 , scEta" << endl ;
	TLorentzVector vp ;
	vp.SetPtEtaPhiM(phoEt_[index_pho],phoEta[index_pho],phoPhi[index_pho],0);

	double ecalIso = phoEcalIsoDR04[index_pho] ;
	double hcalIso = phoHcalIsoDR04[index_pho] ;
	double trkIso =  phoTrkIsoHollowDR04[index_pho] ;
	cout << index_pho << " , " << vp.Pt() << " , " << phoEt[index_pho] << " , " << vp.Eta() << " , " << vp.Phi() ;
	cout << " , " << vp.DeltaR(ve1) << " , " << vp.DeltaR(ve2) << " , " << phoSCEta[index_pho] << endl ;
	cout << "photon Id Variables : sieie , hoe, ecalIso      , trkIso     , hcalIso     " << endl ;
	cout << ": ecalIsoForCut, trkIsForCut, hcalIsForCut,pixelSeed overlap" << endl ;
	cout << phoSigmaIEtaIEta[index_pho] << ", " << phoHoverE[index_pho] << " , " << phoEcalIsoDR04[index_pho] << " , ";
	cout << phoHcalIsoDR04[index_pho] << " , " << phoTrkIsoHollowDR04[index_pho] << " , " << phoOverlap[index_pho] << endl ;
	if ( fabs(phoSCEta[index_pho]) < 1.4442 ){
				trkIso  -= 0.167 * rho ;
				ecalIso -= 0.183 * rho ;
				hcalIso -= 0.062 * rho ;
	}
	else if ( fabs(phoSCEta[index_pho]) > 1.566 && fabs(phoSCEta[index_pho]) < 2.5 ){
				trkIso  -=  0.032 * rho ;
				ecalIso -=  0.090 * rho ;
				hcalIso -=  0.180 * rho ;
	}
				ecalIso = ecalIso - 0.006 * phoEt_[index_pho] ;
				hcalIso = hcalIso - 0.0025 * phoEt_[index_pho] ;
				trkIso  = trkIso  - 0.001 * phoEt_[index_pho] ;
	cout << ecalIso << " " << trkIso << " " << hcalIso << " " << phohasPixelSeed[index_pho] << endl ;
	cout << "passID = " << phoPassID(index_pho,3) << endl ; 
	cout << "=================================================================" << endl ;

}
void anaVgNtuple::printElectron(Int_t iEle){
	if (iEle < 0) cout << "!!!!!! Electron Index = " << iEle << " is not allowed !!!!!!!!" << endl ;
	if (iEle < 0) return ;
		cout << "==================== " ;
		cout << "Now print information of electron " << setw(2) << iEle ;
		cout << " ====================" ;
		cout << endl ;

	// Print 4 vector of this electron
		cout << setw(13) << "pt"          << " |" ;
		cout << setw(13) << "eta"         << " |" ;
		cout << setw(13) << "phi"         << " |" ;
		cout << setw(13) << "energy"      << " |" ;
	if (doEMCorrection) 
	{
		cout << setw(13) << "corr. pt"    << " |" ;
		cout << setw(13) << "corr. en"    << " |" ;
	}
	cout << endl ;
	cout << fixed ;
	cout << setw(13) << setprecision(3) << elePt[iEle]  << " |" ;
	cout << setw(13) << setprecision(3) << eleEta[iEle] << " |" ;
	cout << setw(13) << setprecision(3) << elePhi[iEle] << " |" ;
	cout << setw(13) << setprecision(3) << eleEn[iEle]  << " |" ;
	if (doEMCorrection) 
	{
	cout << setw(13) << setprecision(3) << elePt_[iEle] << " |" ;
	cout << setw(13) << setprecision(3) << eleEn_[iEle] << " |" ;
	}
	cout << endl ;

	// Print trigger information
	cout << "Trigger Matching : " << endl ;
	TString TrgName[] = 
	{
		"HLT_Photon10_L1R" ,
		"HLT_Photon15_Cleaned_L1R" ,
		"HLT_Ele15_LW_L1R" ,
		"HLT_Ele15_SW_L1R" ,
		"HLT_Ele15_SW_CaloEleId_L1R" ,
		"HLT_Ele17_SW_CaloEleId_L1R" ,
		"HLT_Ele17_SW_TightEleId_L1R" ,
		"HLT_Ele17_SW_TighterEleIdIsol_L1R" ,
		"HLT_Ele17_SW_TighterEleIdIsol_L1R_v2" ,
		"HLT_Ele17_SW_TighterEleIdIsol_L1R_v3" ,
		"HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1" ,
		"HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2" ,
		"HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3" ,
		"HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1" ,
		"HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2" ,
		"HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3" ,
		"HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4" ,
		"HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1" ,
		"HLT_Ele25_WP80_PFMT40_v1" ,
		"HLT_Ele27_WP70_PFMT40_PFMHT20_v1" ,
		"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1" ,
		"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2" ,
		"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3" ,
		"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4" ,
		"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5" ,
		"HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5" ,
		"HLT_Ele27_WP80_PFMT50_v1" ,
		"HLT_Ele32_WP70_PFMT50_v1" ,
		"HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6" ,
		"HLT_Ele8_CaloIdL_CaloIsoVL_v*" ,
		"HLT_Ele17_CaloIdL_CaloIsoVL_v*" 
	};
		cout << "              " << setw(82) << "TriggerName |" << " Matched |" <<endl ;
	for (int iTrg = 0; iTrg < 31; iTrg++){
		cout << "TriggerInfo : " << setw(80) << TrgName[iTrg] << " |" << setw(10) ;
		if (eleTrg[iEle][iTrg] == 1) 
			cout << "YES |" << endl ;
		else 
			cout << "NO |" << endl ;
	}
}
Float_t anaVgNtuple::getCorrection_YM(Float_t en, Float_t absEta, Float_t r9, 
		   Int_t run, Bool_t data)
{
  double eta = absEta;
  double mz = 91.19;
  // list of correction coefficiencints for EB
  // |eta| < 1.0 r9 > 0.94, |eta| < 1.0 r9 < 0.94, 4th r9 > 0.94, 4th r9 < 0.94
  float scaleDataCorrEB[4] = { 0.000, -0.355, -0.210,  1.480}; 
  float scaleMCCorrEB[4]   = { 0.190, -0.635, -0.758,  0.445};
  float smearMCCorrEB[4]   = { 0.950,  0.000,  0.322,  1.750}; 

  // list of correction coefficients for EE
  float scaleDataCorrEE[4] = { 0.00,   2.40,   0.00,   1.05}; 
  float scaleMCCorrEE[4]   = {-1.70,   0.58,  -1.57,  -1.25};
  float smearMCCorrEE[4]   = { 2.48,   2.56,   3.32,   2.60}; 

  // Corrections for data, EB
  if ( data && eta < 1.5 ) {
    // run ranges between interventions
    bool run1  = run < 161016;      // very small
    bool run2  = run >= 161016 && run < 162762; // technical stop
    bool run3  = run >= 162762 && run < 163869; // 
    bool run4  = run >= 163869 && run < 165204; // technical stop
    bool run5  = run >= 165204 && run < 166380; // 
    bool run6  = run >= 166380 && run < 166839; // noisy period
    bool run7  = run >= 166839 && run < 167284;
    bool run8  = run >= 167284 && run < 167674; // noisy period
    bool run9  = run >= 167674 && run < 167898;

    bool run10 = run >= 167898 && run < 176100;
    bool run11 = run >= 176100 && run < 176400;
    bool run12 = run >= 176400 && run < 176600;
    bool run13 = run >= 176600 && run < 177000;
    bool run14 = run >= 177000 && run < 177500;
    bool run15 = run >= 177500 && run < 178000;
    bool run16 = run >= 178000 && run < 178300;
    bool run17 = run >= 178300 && run < 178600;
    bool run18 = run >= 178600 && run < 179100;
    bool run19 = run >= 179100 && run < 179600;
    bool run20 = run >= 179600;
    

    int iRun = 0;
    if ( run1  ) iRun = 0;
    if ( run2  ) iRun = 1;
    if ( run3  ) iRun = 2;
    if ( run4  ) iRun = 3;
    if ( run5  ) iRun = 4;
    if ( run6  ) iRun = 5;
    if ( run7  ) iRun = 6;
    if ( run8  ) iRun = 7;
    if ( run9  ) iRun = 8;
    if ( run10 ) iRun = 9;
    if ( run11 ) iRun = 10;
    if ( run12 ) iRun = 11;
    if ( run13 ) iRun = 12;
    if ( run14 ) iRun = 13;
    if ( run15 ) iRun = 14;
    if ( run16 ) iRun = 15;
    if ( run17 ) iRun = 16;
    if ( run18 ) iRun = 17;
    if ( run19 ) iRun = 18;
    if ( run20 ) iRun = 19;

    float corr[20] = {1 - 0.130/mz,  // 1
		      1 + 0.040/mz,  // 2
		      1 - 0.200/mz,  // 3
		      1 - 0.000/mz,  // 4
		      1 + 0.112/mz,  // 5
		      1 + 0.279/mz,  // 6
		      1 + 0.373/mz,  // 7
		      1 + 0.120/mz,  // 8
		      1 + 0.611/mz,  // 9
		      1 + 0.206/mz,  // 10
		      1 + 0.580/mz,  // 11
		      1 + 0.528/mz,  // 12
		      1 + 0.390/mz,  // 13
		      1 + 0.593/mz,  // 14
		      1 + 0.485/mz,  // 15
		      1 + 0.639/mz,  // 16
		      1 + 0.617/mz,  // 17
		      1 + 0.829/mz,  // 18
		      1 + 0.676/mz,  // 19
		      1 + 1.171/mz   // 20
    }; 

    // data scaling factors
    en *= corr[iRun];
    if ( eta < 1.0 ) {
      if ( r9 <= 0.94 ) en *= 1 + scaleDataCorrEB[1]/mz; // DONE!
    } else {
      // 4th module
      en *= 1 + scaleDataCorrEB[2]/mz; // r9 > 0.94; DONE!
      if ( r9 <= 0.94 ) 
	// additional correction for r9 < 0.94 DONE!
	en *= 1 + scaleDataCorrEB[3]/mz; 
    }
  }    

  // MC ***************************************************************
  // Corrections for MC, EB
  if ( !data && eta < 1.5 ) {
    // energy scale
    en *= 1 + scaleMCCorrEB[0]/mz;
    en *= myRandom->Gaus(1, smearMCCorrEB[0]/mz); 

    if ( eta < 1.0 ) {
      if ( r9 <= 0.94 ) {
	en *= 1 + scaleMCCorrEB[1]/mz; // 0.14
	en *= myRandom->Gaus(1, smearMCCorrEB[1]/mz); // extra 0.35% smearing 
      }
    } else {
      // 4th module
      en *= 1 + scaleMCCorrEB[2]/mz; 
      en *= myRandom->Gaus(1, smearMCCorrEB[2]/mz); // extra 1.8% ~0.80 GeV
      if ( r9 <= 0.94 ) {
	en *= 1 + scaleMCCorrEB[3]/mz;
	en *= myRandom->Gaus(1, smearMCCorrEB[3]/mz); // 1.37
      }
    }

  }

  // EE corrections for data
  if ( data && eta > 1.5 ) {

    // Use Paolo's run ranges
    bool run1 = run <= 167913;
    bool run2 = run > 167913 && run <= 172619;
    bool run3 = run > 172619 && run <= 173692;
    bool run4 = run > 173692 && run <= 177139;
    bool run5 = run > 177139 && run <= 178421;
    bool run6 = run > 178421;

    float eeCorr1[6] = {1 - 1.28/mz,
			1 - 2.90/mz,
			1 - 1.65/mz,
			1 - 3.70/mz,
			1 - 4.14/mz,
			1 - 4.14/mz
    };

    float eeCorr2[6] = {1 + 3.00/mz,
			1 - 1.85/mz,
			1 - 1.84/mz,
			1 - 1.33/mz,
			1 - 2.05/mz,
                        1 - 1.62/mz};
    int iRun = 0;
    if ( run1 ) iRun = 0;
    if ( run2 ) iRun = 1;
    if ( run3 ) iRun = 2;
    if ( run4 ) iRun = 3;
    if ( run5 ) iRun = 4;
    if ( run6 ) iRun = 5;

    if ( eta < 2.0 )  {
      en *= eeCorr1[iRun];
      if ( r9 <= 0.94 ) en *= 1 + scaleDataCorrEE[1]/mz;
    } else {
      en *= eeCorr2[iRun];
      if ( r9 <= 0.94 ) en *= 1 + scaleDataCorrEE[3]/mz;
    }
  }
  // MC
  if ( !data && eta > 1.5 ) {
    // scales
    if ( eta <  2.0 && r9 >  0.94 ) en *= 1 + scaleMCCorrEE[0]/mz;
    if ( eta <  2.0 && r9 <= 0.94 ) en *= 1 + scaleMCCorrEE[1]/mz;
    if ( eta >= 2.0 && r9 >  0.94 ) en *= 1 + scaleMCCorrEE[2]/mz;
    if ( eta >= 2.0 && r9 <= 0.94 ) en *= 1 + scaleMCCorrEE[3]/mz;

    // resolutions
    if ( eta <  2.0 && r9 >  0.94 ) en *= myRandom->Gaus(1, smearMCCorrEE[0]/mz);
    if ( eta <  2.0 && r9 <= 0.94 ) en *= myRandom->Gaus(1, smearMCCorrEE[1]/mz);
    if ( eta >= 2.0 && r9 >  0.94 ) en *= myRandom->Gaus(1, smearMCCorrEE[2]/mz);
    if ( eta >= 2.0 && r9 <= 0.94 ) en *= myRandom->Gaus(1, smearMCCorrEE[3]/mz);
  }

  return en;
}

Float_t getCorrection_FSR(Float_t OutValue_, Float_t Pt_, 
Float_t scEta , Float_t R9_,
Int_t option =  0)
{
	// outvalue : could be energy or pt
	//option - 0, defaule , 1, + uncertainty, -1 - uncertainty
	Int_t cat = 0; 
	if (fabs(scEta) < 1.44442 ){
		if ( R9_ > 0.94 ) cat = 0 ;
		else cat = 1 ;
	} else {
		if ( R9_ > 0.95 ) cat = 2 ;
		else cat = 3 ;
	}
	// category : 
	// 0 : EB (R9>0.94) 1 : EB (R9<0.94)
	// 2 : EE (R9>0.94) 3 : EE (R9<0.94)
	Float_t scale[][4] = // [ptBin][category]
	{
		{ 1.049 , 1.027 , 1.048 , 1.038 },
		{ 1.003 , 0.988 , 1.029 , 1.009 },
		{ 1.006 , 0.983 , 1.003 , 0.991 },
		{ 1.007 , 0.993 , 1.010 , 0.977 },
		{ 0.997 , 0.981 , 1.034 , 0.984 },
		{ 0.993 , 0.991 , 1.005 , 0.993 }
	};
	Float_t error[][4] = // [run][category]
	{
		{ 0.114 , 0.141 , 0.104 , 0.133 },
		{ 0.110 , 0.126 , 0.111 , 0.129 },
		{ 0.102 , 0.100 , 0.114 , 0.113 },
		{ 0.083 , 0.085 , 0.079 , 0.095 },
		{ 0.071 , 0.071 , 0.071 , 0.069 },
		{ 0.066 , 0.065 , 0.063 , 0.069 }
	};
	//Float_t scale[][4] = // [ptBin][category]
	//{
	//	{ 1.048 , 1.027 , 1.055 , 1.030 },
	//	{ 1.017 , 0.989 , 1.011 , 0.990 },
	//	{ 1.008 , 0.985 , 1.003 , 0.981 },
	//	{ 1.006 , 0.991 , 0.998 , 0.976 },
	//	{ 0.999 , 0.985 , 1.014 , 0.975 },
	//	{ 0.998 , 0.980 , 0.988 , 0.985 }
	//};
	//Float_t error[][4] = // [run][category]
	//{
	//	{ 0.119 , 0.135 , 0.121 , 0.144 },
	//	{ 0.113 , 0.127 , 0.122 , 0.137 },
	//	{ 0.097 , 0.101 , 0.112 , 0.116 },
	//	{ 0.080 , 0.086 , 0.088 , 0.094 },
	//	{ 0.067 , 0.070 , 0.072 , 0.073 },
	//	{ 0.067 , 0.063 , 0.070 , 0.074 }
	//};
	Int_t iPt = -1 ;
	if (Pt_ > 10. && Pt_ < 12. ) iPt = 0 ;
	if (Pt_ > 12. && Pt_ < 15. ) iPt = 1 ;
	if (Pt_ > 15. && Pt_ < 20. ) iPt = 2 ;
	if (Pt_ > 20. && Pt_ < 25. ) iPt = 3 ;
	if (Pt_ > 25. && Pt_ < 30. ) iPt = 4 ;
	if (Pt_ > 30.              ) iPt = 5 ;
	if ( iPt == -1 ) return OutValue_ ;
	else return OutValue_ * ( scale[iPt][cat] + option * error[iPt][cat] ) ;
}

void anaVgNtuple::EcorrectionYM(TString flavor)
{
	float YMcorr[10] = {
		1 - 0.360/91.19, 
		1 - 0.140/91.19, 
		1 - 0.008/91.19, 
		1 - 0.010/91.19,
		1 + 0.036/91.19,
		1 + 0.046/91.19,
		1 + 0.160/91.19,
		1 + 0.081/91.19,
		1 + 0.435/91.19,
		1 + 0.571/91.19
	};
	//float YMcorrMC = 1 + 0.284/91.19;
	//float YMsmear = 0.75 / 91.19;
	//float YMsmearExtra = 2.6/91.19;  // this is for basket 4

	float YMEEcorr = 1 - 1.85/91.19;
	//float YMEEcorrMC = 1 - 1.62/91.19; // mc scale is imperfect
	//float YMEEsmear[2] = {1.8 / 91.19,
	//2.3 / 91.19};

	// run ranges between interventions
	int iRun = 0;
	if      ( run <  161016 && run != 0     ) iRun = 0 ; // very small
	else if ( run >= 161016 && run < 162762 ) iRun = 1 ; // technical stop
	else if ( run >= 162762 && run < 163869 ) iRun = 2 ; // 
	else if ( run >= 163869 && run < 165204 ) iRun = 3 ; // technical stop
	else if ( run >= 165204 && run < 166380 ) iRun = 4 ; // 
	else if ( run >= 166380 && run < 166839 ) iRun = 5 ; // noisy period
	else if ( run >= 166839 && run < 167284 ) iRun = 6 ; 
	else if ( run >= 167284 && run < 167674 ) iRun = 7 ; // noisy period
	else if ( run >= 167674 && run < 167898 ) iRun = 8 ; 
	else if ( run >= 167898                 ) iRun = 9 ;

	if ( flavor.Contains("Ele") || flavor.Contains("ele") )
	{
	elePt_.clear(); 
	eleEn_.clear();
		for (int i=0; i<nEle; ++i) 
		{  
			Float_t Escale1 = 1;
			//Float_t Esmear1= 1;
			if(fabs(eleEta[i])<1.4442)
				Escale1 = YMcorr[iRun];
			else 
				Escale1 = YMEEcorr;
			elePt_.push_back(elePt[i]*Escale1) ;
			eleEn_.push_back(eleEn[i]*Escale1) ;
		}
	}
}

Float_t getCorrection(Float_t En_, 
Float_t scEta , Float_t R9_, Int_t run_,
Int_t option =  0){
	//option - 0, defaule , 1, + uncertainty, -1 - uncertainty
	Float_t myEn = En_ ;
	Int_t cat = 0; 
	if (fabs(scEta) < 1.44442 ){
		if ( R9_ > 0.94 ) cat = 0 ;
		else cat = 1 ;
	} else {
		if ( R9_ > 0.94 ) cat = 2 ;
		else cat = 3 ;
	}
	// categgory : 
	// 0 : EB (R9>0.94) 1 : EB (R9<0.94)
	// 2 : EB (R9>0.94) 3 : EB (R9<0.94)
	Float_t scale[][4] = // [run][class]
	{
		{ 0.47  , -0.25 , -0.58 , 0.10  },
		{ 0.07  , -0.49 , -2.49 , -0.62 },
		{ -0.03 , -0.67 , -3.76 , -1.33 },
		{ -0.11 , -0.63 , -4.50 , -1.78 },
		{ -0.14 , -0.74 , -5.61 , -2.73 },
		{ 0.00  , 0.00  , -0.11 , -0.06 }
	};
	Float_t error[][4] = // [run][category]
	{
		{ 0.05  , 0.04 , 0.19 , 0.16 },
		{ 0.07  , 0.04 , 0.22 , 0.18 },
		{ 0.06  , 0.04 , 0.19 , 0.15 },
		{ 0.06  , 0.04 , 0.20 , 0.15 },
		{ 0.05  , 0.03 , 0.18 , 0.13 },
		{ 0.04  , 0.02 , 0.09 , 0.07 }
	};
	Int_t iRun = -1 ;
	if      ( run_ >= 160431 && run_ <= 163869 ) iRun = 0 ;
	else if ( run_ >= 165071 && run_ <= 165970 ) iRun = 1 ;
	else if ( run_ >= 165971 && run_ <= 166502 ) iRun = 2 ;
	else if ( run_ >= 166503 && run_ <= 166861 ) iRun = 3 ;
	else if ( run_ >= 166862 && run_ <= 167784 ) iRun = 4 ;
	else if ( run_ >= 166863                   ) iRun = 4 ;
//	if ( run_ >= 160431 && run_ <= 167784 ) iRun = 0 ;
	if ( iRun == -1 ) return En_ ;
	myEn = En_ * ( 1. - ( 0.01 * ( scale[iRun][cat] + option * error[iRun][cat] ) ) ) ;
//	myEn = En_ * ( 1. - ( 0.01 *  scale[iRun][cat]  ) ) ;
//	cout << "before " << En_ << " after " << myEn << " run " << run_ << " iRun " << iRun << " class " << cat << endl ;
	return myEn ;
}
vector<Float_t> anaVgNtuple::newEleEn(){
	vector<Float_t> _eleEn ;
	float corrEn(0.);
	for (int i = 0; i < nEle; i++){
		Float_t r9 = eleE3x3[i] / eleSCRawEn[i] ;
//		corrEn = getCorrection(eleEn[i],eleSCEta[i],r9,run,EnCorrShift);
//		corrEn = getCorrection_FSR(eleEn[i],elePt[i],eleSCEta[i],r9);
		corrEn = getCorrection_YM(eleEn[i],fabs(eleSCEta[i]),r9,run,isData);
		_eleEn.push_back(corrEn);
	}
	return _eleEn ;
}
vector<Float_t> anaVgNtuple::newPhoE(){
	vector<Float_t> _phoEn ;
	float corrEn(0.);
	for (int i = 0; i < nPho; i++){
//		corrEn = getCorrection(phoE[i],phoSCEta[i],phoR9[i],run);
//		corrEn = getCorrection_FSR(phoE[i],phoEt[i],phoSCEta[i],phoR9[i],FSRCorrShift);
		corrEn = getCorrection_YM(phoE[i],fabs(phoSCEta[i]),phoR9[i],run,isData);
		_phoEn.push_back(corrEn);
	}
	return _phoEn ;
}
vector<Float_t> anaVgNtuple::newElePt(vector<Float_t> EnCorr){
	float theta ;
	vector<Float_t> _elePt ;
	unsigned int nnEle = nEle ;
	for (int i = 0; i < nEle; i++){
		if ( nnEle == EnCorr.size() ){
			theta = 2. * atan(exp(-eleEta[i]));
			_elePt.push_back( fabs(sin(theta)) * EnCorr[i] );
		}
		else {
			_elePt.push_back(elePt[i]) ;
		}
	}
	return _elePt ;
}
vector<Float_t> anaVgNtuple::newPhoEt(vector<Float_t> EnCorr){
	float theta ;
	vector<Float_t> _phoPt ;
	unsigned int nnPho = nPho ;
	for (int i = 0; i < nPho; i++){
		if ( nnPho == EnCorr.size() ){
			theta = 2. * atan(exp(-phoEta[i]));
			_phoPt.push_back( fabs(sin(theta)) * EnCorr[i] );
		}
		else {
			_phoPt.push_back(phoEt[i]) ;
		}
	}
	return _phoPt ;
}
vector<int> anaVgNtuple::mcMuon(){
	vector<int> iMC ;
	iMC.clear();
	// find electron decay form W or Z
	for (int i = 0; i < nMC; i++){
		if ( fabs(mcPID[i]) == 13 && mcPt[i] > 10. && ( abs(mcGMomPID[i]) == 23 || abs(mcMomPID[i]) == 23 ) ) {
			// remove decay from tau
			if (abs(mcGMomPID[i]) == 12 || abs(mcGMomPID[i]) == 14 || abs(mcGMomPID[i]) == 16 ) continue ;
				if (abs(mcMomPID[i]) == 12 || abs(mcMomPID[i]) == 14 || abs(mcMomPID[i]) == 16 ) continue ;
			iMC.push_back(i);
		}
		else
			if ( fabs(mcPID[i]) == 13 && mcPt[i] > 10. && ( abs(mcGMomPID[i]) == 24 || abs(mcMomPID[i]) == 24 ) ) {
				// remove decay from tau
				if (abs(mcGMomPID[i]) == 12 || abs(mcGMomPID[i]) == 14 || abs(mcGMomPID[i]) == 16 ) continue ;
				if (abs(mcMomPID[i]) == 12 || abs(mcMomPID[i]) == 14 || abs(mcMomPID[i]) == 16 ) continue ;
				iMC.push_back(i);
			}
	}
	return iMC ;
}
Bool_t anaVgNtuple::mcMuMatcher(Int_t particleIndex, TString particleName){
	if (isData) return false ;
	// do mc match for spacific object
	float eta(0), phi(0) ;
	if (particleName.Contains("jet") ||
		particleName.Contains("Jet"))
	{
		eta = jetEta[particleIndex] ;
		phi = jetPhi[particleIndex] ;
	}
	else
	if (particleName.Contains("ele") || 
		particleName.Contains("Ele") )
	{
		eta = eleEta[particleIndex] ;
		phi = elePhi[particleIndex] ;
	}
	else
	if (particleName.Contains("mu") || 
		particleName.Contains("Mu") )
	{
		eta =  muEta[particleIndex] ;
		phi =  muPhi[particleIndex] ;
	}
	else
	if (particleName.Contains("pho") || 
		particleName.Contains("Pho") )
	{
		eta = phoEta[particleIndex] ;
		phi = phoPhi[particleIndex] ;
	}
	else return false ;
	vector<int> vecMu = mcMuon();
	for (vector<int>::iterator it = vecMu.begin(); 
			it != vecMu.end(); it++)
	{
		if (DeltaR(eta,phi,mcEta[*it],mcPhi[*it]) < 0.15 ) return true ;
	}
	return false ;
}
vector<int> anaVgNtuple::mcElectron(){
	vector<int> iMC ;
	iMC.clear();
	// find electron decay form W or Z
	for (int i = 0; i < nMC; i++){
		if ( fabs(mcPID[i]) == 11 && mcPt[i] > 10. && ( abs(mcGMomPID[i]) == 23 || abs(mcMomPID[i]) == 23 ) ) {
			// remove decay from tau
			if (abs(mcGMomPID[i]) == 12 || abs(mcGMomPID[i]) == 14 || abs(mcGMomPID[i]) == 16 ) continue ;
			if (abs(mcMomPID[i]) == 12 || abs(mcMomPID[i]) == 14 || abs(mcMomPID[i]) == 16 ) continue ;
			iMC.push_back(i);
		}
		else
			if ( fabs(mcPID[i]) == 11 && mcPt[i] > 10. && ( abs(mcGMomPID[i]) == 24 || abs(mcMomPID[i]) == 24 ) ) {
				// remove decay from tau
				if (abs(mcGMomPID[i]) == 12 || abs(mcGMomPID[i]) == 14 || abs(mcGMomPID[i]) == 16 ) continue ;
				if (abs(mcMomPID[i]) == 12 || abs(mcMomPID[i]) == 14 || abs(mcMomPID[i]) == 16 ) continue ;
				iMC.push_back(i);
			}
		else
			if ( fabs(mcPID[i]) == 11 && mcPt[i] > 10. && ( abs(mcGMomPID[i]) == 11 && abs(mcGMomPID[i]) == 11 ) ) {
				// remove decay from tau
				iMC.push_back(i);
			}
	}
	return iMC ;
}
Bool_t anaVgNtuple::mcEleMatcher(Int_t particleIndex, TString particleName){
	if (isData) return false ;
	// do mc match for spacific object
	float eta(0), phi(0) ;
	if (particleName.Contains("jet") ||
		particleName.Contains("Jet"))
	{
		eta = jetEta[particleIndex] ;
		phi = jetPhi[particleIndex] ;
	}
	else
	if (particleName.Contains("ele") || 
		particleName.Contains("Ele") )
	{
		eta = eleEta[particleIndex] ;
		phi = elePhi[particleIndex] ;
	}
	else
	if (particleName.Contains("mu") || 
		particleName.Contains("Mu") )
	{
		eta =  muEta[particleIndex] ;
		phi =  muPhi[particleIndex] ;
	}
	else
	if (particleName.Contains("pho") || 
		particleName.Contains("Pho") )
	{
		eta = phoEta[particleIndex] ;
		phi = phoPhi[particleIndex] ;
	}
	else return false ;
	vector<int> vecEle = mcElectron();
	for (vector<int>::iterator it = vecEle.begin(); 
			it != vecEle.end(); it++)
	{
		if (DeltaR(eta,phi,mcEta[*it],mcPhi[*it]) < 0.15 ) return true ;
	}
	return false ;
}
vector<int> anaVgNtuple::mcPhoton(){
	vector<int> iMC ;
	iMC.clear();
	//for (int i = 0; i < nMC; i++){
	//	bool mcFlag = false ;
	//	if (doShowMCPho && mcPID[i] == 22 ){
	//		cout << "MCPho : index " << i << " Mother=" << mcMomPID[i] << " GMomther=" << mcGMomPID[i] << endl ;
	//	}
	//	if ( mcPID[i] == 22 && abs(mcMomPID[i]) <= 22 ) {
	//		mcFlag = true ;
	//	}
	//	if ( mcPID[i] == 22 && abs(mcGMomPID[i]) <= 24 ) {
	//		mcFlag = true ;
	//	}
	//	if ( mcPID[i] == 22 && abs(mcMomPID[i]) <= 22 ) {
	//		mcFlag = true ;
	//	}
	//	if (mcFlag)
	//		iMC.push_back(i);
	//}
	//return iMC ;

	for (int i = 0; i < nMC; i++){
		if (doShowMCPho && mcPID[i] == 22 ){
			cout << "MCPho : index " << i << " Mother=" << mcMomPID[i] << " GMomther=" << mcGMomPID[i] << endl ;
		}
		if ( mcPID[i] == 22 && mcPt[i] > 10. && abs(mcMomPID[i]) <= 22 ) {
			iMC.push_back(i);
		}
	}
	return iMC ;
}

Bool_t anaVgNtuple::mcPhoMatcher(Int_t iPho){
	if (isData) return false ;
	vector<int> phoMC =  mcPhoton();
	if (phoMC.size() == 0) return false ;
	for (size_t iMC = 0; iMC < phoMC.size(); iMC++){
		if ( DeltaR(mcEta[phoMC[iMC]],mcPhi[phoMC[iMC]],phoEta[iPho],phoPhi[iPho] ) < 0.30 )
			return true ;
	}
	return false ;
}

unsigned int anaVgNtuple::eleTrgWP(Int_t iEle){
	unsigned int value(0);
	unsigned int I(1);
	Float_t ecalIso = ( fabs(eleSCEta[iEle]) < 1.4442 ) ?
		 TMath::Max(eleIsoEcalDR03[iEle]-1.,0.) :
		 eleIsoEcalDR03[iEle] ;
	Float_t hcalIso = eleIsoHcalDR03[iEle] ;
	Float_t trkIso  = eleIsoTrkDR03[iEle] ;
	Int_t   side    = -1 ;
	if   ( fabs(eleSCEta[iEle]) < 1.4442 ) side = 0 ;
	else if ( fabs(eleSCEta[iEle]) > 1.566  && fabs(eleSCEta[iEle]) < 2.5 ) side = 1 ;
	if (side == -1) return value ;
	Float_t cutHoverE[5][2] = {
		{0.15,0.10} ,
		{0.15,0.10} ,
		{0.10,0.10} ,
		{0.15,0.075} ,
		{0.05,0.05} 
	};
	Float_t cutSihih[5][2] = {
		{0.024,0.040} ,
		{0.014,0.035} ,
		{0.014,0.035} ,
		{0.011,0.031} ,
		{0.011,0.031}
	};
	Float_t cutEcalIso[5][2] = {
		{0.2,0.2},
		{0.2,0.2},
		{0.2,0.2},
		{0.125,0.075},
		{0.125,0.075}
	};
	Float_t cutHcalIso[5][2] = {
		{0.2,0.2},
		{0.2,0.2},
		{0.2,0.2},
		{0.125,0.075},
		{0.125,0.075}
	};
	Float_t cutDEtaIn[5][2] = {
		{ 0.01  , 0.01  },
		{ 0.01  , 0.01  },
		{ 0.01  , 0.01  },
		{ 0.008 , 0.008 },
		{ 0.008 , 0.008 }
	};
	Float_t cutDPhiIn[5][2] = {
		{ 0.15  , 0.15  },
		{ 0.15  , 0.15  },
		{ 0.15  , 0.15  },
		{ 0.07  , 0.05  },
		{ 0.07  , 0.05  }
	};
	Float_t cutTrkIso[5][2] = {
		{ 0.20  , 0.20  },
		{ 0.20  , 0.20  },
		{ 0.20  , 0.20  },
		{ 0.125 , 0.075 },
		{ 0.125 , 0.075 }
	};
	// CaloId
	if ( eleHoverE[iEle]        < cutHoverE[0][side] &&
		 eleSigmaIEtaIEta[iEle] < cutSihih[0][side] ) value = value | (I<<_eTrg_CaloIdVL) ;
	if ( eleHoverE[iEle]        < cutHoverE[1][side] &&
		 eleSigmaIEtaIEta[iEle] < cutSihih[1][side] ) value = value | (I<<_eTrg_CaloIdL ) ;
	if ( eleHoverE[iEle]        < cutHoverE[2][side] &&
		 eleSigmaIEtaIEta[iEle] < cutSihih[2][side] ) value = value | (I<<_eTrg_CaloIdXL  ) ;
	if ( eleHoverE[iEle]        < cutHoverE[3][side] &&
		 eleSigmaIEtaIEta[iEle] < cutSihih[3][side] ) value = value | (I<<_eTrg_CaloIdT   ) ;
	if ( eleHoverE[iEle]        < cutHoverE[4][side] &&
		 eleSigmaIEtaIEta[iEle] < cutSihih[4][side] ) value = value | (I<<_eTrg_CaloIdVT  ) ;
	// CaloIso
	if ( ecalIso / elePt[iEle]  < cutEcalIso[0][side] && 
		 hcalIso / elePt[iEle]  < cutHcalIso[0][side]  ) value = value | (I<<_eTrg_CaloIsoVL ) ;
	if ( ecalIso / elePt[iEle]  < cutEcalIso[1][side] && 
		 hcalIso / elePt[iEle]  < cutHcalIso[1][side]  ) value = value | (I<<_eTrg_CaloIsoL  ) ;
	if ( ecalIso / elePt[iEle]  < cutEcalIso[2][side] && 
		 hcalIso / elePt[iEle]  < cutHcalIso[2][side]  ) value = value | (I<<_eTrg_CaloIsoXL ) ;
	if ( ecalIso / elePt[iEle]  < cutEcalIso[3][side] && 
		 hcalIso / elePt[iEle]  < cutHcalIso[3][side]  ) value = value | (I<<_eTrg_CaloIsoT  ) ;
	if ( ecalIso / elePt[iEle]  < cutEcalIso[4][side] && 
		 hcalIso / elePt[iEle]  < cutHcalIso[4][side]  ) value = value | (I<<_eTrg_CaloIsoVT ) ;
	// TrkId
	if ( eledEtaAtVtx[iEle]     < cutDEtaIn[0][side] && 
		 eledPhiAtVtx[iEle]     < cutDPhiIn[0][side]  ) value = value | (I<<_eTrg_TrkIdVL ) ;
	if ( eledEtaAtVtx[iEle]     < cutDEtaIn[1][side] && 
		 eledPhiAtVtx[iEle]     < cutDPhiIn[1][side]  ) value = value | (I<<_eTrg_TrkIdL  ) ;
	if ( eledEtaAtVtx[iEle]     < cutDEtaIn[2][side] && 
		 eledPhiAtVtx[iEle]     < cutDPhiIn[2][side]  ) value = value | (I<<_eTrg_TrkIdXL ) ;
	if ( eledEtaAtVtx[iEle]     < cutDEtaIn[3][side] && 
		 eledPhiAtVtx[iEle]     < cutDPhiIn[3][side]  ) value = value | (I<<_eTrg_TrkIdT  ) ;
	if ( eledEtaAtVtx[iEle]     < cutDEtaIn[4][side] && 
		 eledPhiAtVtx[iEle]     < cutDPhiIn[4][side]  ) value = value | (I<<_eTrg_TrkIdVT ) ;
	// TrkIso
	if ( trkIso  / elePt[iEle]  < cutTrkIso[0][side]  ) value = value | (I<<_eTrg_TrkIsoVL) ;
	if ( trkIso  / elePt[iEle]  < cutTrkIso[1][side]  ) value = value | (I<<_eTrg_TrkIsoL ) ;
	if ( trkIso  / elePt[iEle]  < cutTrkIso[2][side]  ) value = value | (I<<_eTrg_TrkIsoXL) ;
	if ( trkIso  / elePt[iEle]  < cutTrkIso[3][side]  ) value = value | (I<<_eTrg_TrkIsoT ) ;
	if ( trkIso  / elePt[iEle]  < cutTrkIso[4][side]  ) value = value | (I<<_eTrg_TrkIsoVT) ;
	return value; 
}


Bool_t anaVgNtuple::passHLTMu(Int_t &lepTrg_index){
	if (isData){
		if (HLTIndex[107] >= 0){
			lepTrg_index = 28 ;
			if ( HLT[HLTIndex[107]] == 1 ) return true ;
		}
		else
		if (HLTIndex[108] >= 0){
			lepTrg_index = 29 ;
			if ( HLT[HLTIndex[108]] == 1 ) return true ;
		}
		else
		if (HLTIndex[110] >= 0){
			lepTrg_index = 31 ;
			if ( HLT[HLTIndex[110]] == 1 ) return true ;
		}
		else {
			lepTrg_index = -1 ;
			return false ;
		}
	} else {
		if (HLTIndex[107] >= 0){
			lepTrg_index = 28 ;
			if ( HLT[HLTIndex[107]] == 1 ) return true ;
		}
		else {
			lepTrg_index = -1 ;
			return false ;
		}
	}
	return false ;
}

Bool_t anaVgNtuple::passHLT(Int_t &lepTrg_index){
	t_HLTprescale = -1 ;
	if (isData){
		if (run >= 160404 && run <= 167913 ){// HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v
			lepTrg_index  = 18 ;
			t_HLTprescale = HLTprescale[HLTIndex[83]] ;
			if (HLT[HLTIndex[83]] == 1) return true ;
			else return false ;
		}
		if (run >= 167039  ){ // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v
			lepTrg_index = 21 ;
			t_HLTprescale = HLTprescale[HLTIndex[86]] ;
			if (HLT[HLTIndex[86]] == 1) return true ;
			else return false ;
		}
	}
	else if (!isData && HLTIndex[83] >= 0 ){
		lepTrg_index = 18 ;
		t_HLTprescale = HLTprescale[HLTIndex[83]] ;
		if (HLT[HLTIndex[83]] == 1) return true ;
		else return false ;
	}
	else if (!isData && HLTIndex[86] >= 0 ){
		lepTrg_index = 21 ;
		t_HLTprescale = HLTprescale[HLTIndex[86]] ;
		if (HLT[HLTIndex[86]] == 1) return true ;
		else return false ;
	}
	return false ;
}
Bool_t anaVgNtuple::checkHLT(Int_t iLep, Int_t lepTrg_index){
	Int_t trgLeg1(0) ;
	Int_t trgLeg2(0) ;
	if ( 
		lepTrg_index == 20 ||
		lepTrg_index == 21 ||
		lepTrg_index == 22 ||
		lepTrg_index == 23 ||
		lepTrg_index == 24 ||
		lepTrg_index == 25 ||
		lepTrg_index == 28
	){
		trgLeg1 = 29 ; // HLT_Ele8_CaloIdL_CaloIsoVL_v
		trgLeg2 = 30 ; // HLT_Ele17_CaloIdL_CaloIsoVL_v
	} else {
		return false ;
	}
	if ( doEle && eleTrg[iLep][lepTrg_index] != 1 ) return false ;
	if ( doMu  &&  muTrg[iLep][lepTrg_index] != 1 ) return false ;
	if ( doEle && ( eleTrg[iLep][lepTrg_index] ||
					eleTrg[iLep][lepTrg_index])   ) return true ;
	if ( doMu  && (  muTrg[iLep][lepTrg_index] ||
					 muTrg[iLep][lepTrg_index])   ) return true ;
	return false ;
}

Bool_t anaVgNtuple::muClean(Int_t iMu){
	if ( fabs(muEta[iMu]) > 2.4 ) return false ;
	if ( muType[iMu] != 14 ) return false ;
	if ( muPt[iMu] < LepPtCut ) return false ;
	return true ;
}

Bool_t anaVgNtuple::muPassID(Int_t iMu){
	if (iMu >= nMu) return false ;
	return (
		fabs(muPVD0[iMu])              < 0.02 &&
		fabs(muPVDz[iMu])              < 0.1  &&
		muChi2NDF[iMu]                 < 10.0 &&
		muStations[iMu]                > 1    &&
		muNumberOfValidPixelHits[iMu]  > 0    &&
		muNumberOfValidTrkHits[iMu]    > 10   &&
		muNumberOfValidMuonHits[iMu]   > 0    &&
		( muIsoEcal[iMu] + muIsoHcal[iMu] + muIsoTrk[iMu] - rho 
		* 0.3 * 0.3 * TMath::Pi() )     < 0.1 * muPt[iMu] 
	) ;
}

Bool_t anaVgNtuple::DiMuPair(Int_t & mu1_index, Int_t& mu2_index , Int_t trg_index){
	mu1_index = -1 ;
	mu2_index = -1 ;
	vector<Int_t> muIarray ;
	for (int iMu = 0; iMu < nMu; iMu++){
		if (muPt[iMu] < LepPtCut )                         continue ;
		if (!muClean(iMu))                                 continue ;
		if (!muPassID(iMu))                                continue ;
		if (trg_index != -1 && muTrg[iMu][trg_index] != 1) continue ;
		muIarray.push_back(iMu);
	}
	if ( muIarray.size() < 2 ) return false ;
	mu1_index = muIarray.at(0) ;
	mu2_index = muIarray.at(1) ;
	return true ;
}

vector<int> anaVgNtuple::ElectronPassLevel(Int_t index_eID, Int_t index_trg){
	vector<int> eleIDLevel ;

	for (int iEle = 0; iEle < nEle; iEle++){
		int Level = 1 ;
		if ( elePt_[iEle] < LepPtCut ) Level = 0 ;
		if ( fabs(eleSCEta[iEle]) > 2.5    ) Level = 0 ;
		if ( fabs(eleSCEta[iEle]) > 1.4442
				&& fabs(eleSCEta[iEle]) < 1.566  ) Level = 0 ;
		if (IsSpike(iEle,TString("ele"))) Level = 0 ;

		if (Level == 0 ){
			eleIDLevel.push_back(Level) ;
			continue ;
		}
		Level = 0 ;
		if ( elePassID(iEle,index_eID) ){
			Level++ ;
			if ( eleTrg[iEle][index_trg] == 1 || noHLT ){
				Level++ ;
			}
		}
		eleIDLevel.push_back(Level) ;

	}// loop over electrons
	return eleIDLevel ;
}

Bool_t anaVgNtuple::DiElePair(
		Int_t &lep1_index , 
		Int_t &lep2_index, 
		vector<int> eID_Level){

	lep1_index = -1 ; lep2_index = -1 ;
	float leading1(0), leading2(0);
	for (int iEle = 0; iEle < nEle; iEle++){
		if ( eID_Level[iEle] < 2 ) continue ;
		if ( leading1 < elePt_[iEle] && leading2 < elePt_[iEle] ){
			leading2 = leading1 ;
			lep2_index = lep1_index ;
			leading1 = elePt_[iEle] ;
			lep1_index = iEle ;
		} else if ( leading2 < elePt_[iEle] && leading1 > elePt_[iEle] ) {
			lep2_index = iEle ;
			leading2 = elePt_[iEle] ;
		}
	}
	if (lep1_index != -1 && lep2_index != -1 && lep1_index != lep2_index) return true;
	else return false ;
}

float anaVgNtuple::DeltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2){
	float dEta = eta1 - eta2 ;
	float dPhi = phi1 - phi2 ;
	if (dPhi >  TMath::Pi()) dPhi -= 2*TMath::Pi();
	if (dPhi < -TMath::Pi()) dPhi += 2*TMath::Pi();
	return sqrt(pow(dEta,2)+pow(dPhi,2))  ;
}

Bool_t anaVgNtuple::IsISR(Int_t iPho){
	if (doShowMCPho){
		cout << "ISR check index " << iPho << 
		" mother " << phoGenMomPID[iPho] <<
		" grandMother " << phoGenGMomPID[iPho] << endl ;
	}
	if (isData) return false ;
	if ( phoGenMomPID[iPho]==21 ||
			TMath::Abs(phoGenMomPID[iPho])<7) return true ;
	if ( phoGenMomPID[iPho] == 22 ){
		if ( phoGenGMomPID[iPho]==21 ||
				TMath::Abs(phoGenGMomPID[iPho])<7) return true ;
	}
	return false ;
}
Bool_t anaVgNtuple::checkOption(){
	if ( doEle && doMu ){
		cout << "Error !! Electron and muon channel cannot run at the same time" << endl ;
		return false ;
	}
	bool optionPass = true ;
	if (doPhoStandard && doRatioA ) optionPass = false ;
	if (doPhoStandard && doRatioZK ) optionPass = false ;
	Int_t nProcess = 0 ;
	if ( doTemplate ) nProcess++ ;
	if ( doRatioA   ) nProcess++ ;
	if ( doRatioZK      ) nProcess++ ;
	if ( doHistoPassPho ) nProcess++ ;
	if ( nProcess != 1 ) optionPass = false ;
	if ( !optionPass ) {
		cout << "option not allowed !! " << endl ;
	}
	return optionPass ;
}
Bool_t anaVgNtuple::IsFSR(Int_t iPho){
	if (isData) return false ;
	if( (TMath::Abs(phoGenMomPID[iPho])==11  ||
				TMath::Abs(phoGenMomPID[iPho])==13  ||
				TMath::Abs(phoGenMomPID[iPho])==15) &&
			(TMath::Abs(phoGenGMomPID[iPho])==11 ||
			 TMath::Abs(phoGenGMomPID[iPho])==13 ||
			 TMath::Abs(phoGenGMomPID[iPho])==15))   return true ;
	if( (TMath::Abs(phoGenMomPID[iPho])==11 ||
				TMath::Abs(phoGenMomPID[iPho])==13 ||
				TMath::Abs(phoGenMomPID[iPho])==15) &&
			(TMath::Abs(phoGenGMomPID[iPho])==23 ||
			 TMath::Abs(phoGenGMomPID[iPho])==24  )) return true ;
	if ( phoGenMomPID[iPho]==22 &&
			( TMath::Abs(phoGenGMomPID[iPho]==11) ||
			  TMath::Abs(phoGenGMomPID[iPho]==13) ||
			  TMath::Abs(phoGenGMomPID[iPho]==15)
			)
	   ) return true ;
	if ( phoGenMomPID[iPho] == 22 ){
		if ( phoGenGMomPID[iPho]==23 ||
				TMath::Abs(phoGenGMomPID[iPho]) == 24) return true ;
	}
	return false ;
}

Bool_t anaVgNtuple::elePassID(int index_e , int index_eID , unsigned int CutSet){
	Bool_t passID = false ;
	if ( CutSet  == 0 ) return true ;
	if ( fabs(elePVD0[index_e]) > 0.02 || fabs(elePVDz[index_e]) > 0.1 ) return passID ;
	//index_eID 0: wp
	int BE =  ( fabs(eleSCEta[index_e]) < 1.4442) ? 0 : 1 ;
	if (eleConvMissinghit[index_e] != 0) return passID ;
	//ID [WP60,WP70,WP80,WP85,WP90,WP95][Dist,Dcot,cIso,SIEIE,dPhiIn,dEtaIn]
	Double_t cutEB[6][6] = { 
		{0.02,0.02, 0.016 ,0.01  ,0.02  ,0.004},
		{0.02,0.02, 0.03  ,0.01  ,0.02  ,0.004},
		{0.02,0.02, 0.04  ,0.01  ,0.027 ,0.005},
		{0.02,0.02, 0.053 ,0.01  ,0.039 ,0.005},
		{0.,0., 0.085 ,0.01  ,0.071 ,0.007},
		{0.,0., 0.15  ,0.012 ,0.8   ,0.007}
	};
	Double_t cutEE[6][6] = { 
		{0.02,0.02, 0.008 ,0.031 ,0.021 ,0.004},
		{0.02,0.02, 0.016 ,0.031 ,0.021 ,0.005},
		{0.02,0.02, 0.033 ,0.031 ,0.021 ,0.006},
		{0.02,0.02, 0.042 ,0.031 ,0.028 ,0.007},
		{0.,0., 0.051 ,0.031 ,0.047 ,0.011},
		{0.,0., 0.1   ,0.031 ,0.7   ,0.011}
	};
	float cIso = 100. ;
	// calculate cIso
	if ( fabs(eleSCEta[index_e]) < 1.4442 ){
		cIso = ( eleIsoTrkDR03[index_e] + TMath::Max(eleIsoEcalDR03[index_e]-1.,0.) + eleIsoHcalSolidDR03[index_e] - rho*TMath::Pi()*0.3*0.3) / elePt_[index_e];
	}
	else cIso = ( eleIsoTrkDR03[index_e] + eleIsoEcalDR03[index_e] + eleIsoHcalSolidDR03[index_e] - rho*TMath::Pi()*0.3*0.3) / elePt_[index_e];
	// ID 
	if (BE == 0){
		passID = ( (fabs(eleConvDist[index_e]) > cutEB[index_eID][0] || 
					fabs(eleConvDcot[index_e]) > cutEB[index_eID][1]) &&
				cIso < cutEB[index_eID][2] &&
				eleSigmaIEtaIEta[index_e] < cutEB[index_eID][3] &&
				fabs(eledPhiAtVtx[index_e]) < cutEB[index_eID][4] &&
				fabs(eledEtaAtVtx[index_e]) < cutEB[index_eID][5] ) ;
	} else if (BE == 1){
		passID = ( (fabs(eleConvDist[index_e]) > cutEE[index_eID][0] || 
					fabs(eleConvDcot[index_e]) > cutEE[index_eID][1]) &&
				cIso < cutEE[index_eID][2] &&
				eleSigmaIEtaIEta[index_e] < cutEE[index_eID][3] &&
				fabs(eledPhiAtVtx[index_e]) < cutEE[index_eID][4] &&
				fabs(eledEtaAtVtx[index_e]) < cutEE[index_eID][5] ) ;
	}
	return passID ;
}
Bool_t anaVgNtuple::phoPassID(int index_g , int index_gID){
	Int_t requireLevel[] = { 5, 4, 3, 5 } ;
	if (index_gID > int(sizeof(requireLevel)/sizeof(requireLevel[0])) ) return false ;
	return ( phoPassIDLevel(index_g,index_gID) == requireLevel[index_gID] ) ;
}
Int_t anaVgNtuple::phoPassIDLevel(int index_g , int index_gID){
	int flag = 0 ;
	int BE =  ( fabs(phoSCEta[index_g]) < 1.4442) ? 0 : 1 ;

	double ecalIso = phoEcalIsoDR04[index_g] ;
	double hcalIso = phoHcalIsoDR04[index_g] ;
	double trkIso =  phoTrkIsoHollowDR04[index_g] ;
	double hovere =  phoHoverE[index_g] ;
	Double_t PhoIDScale[][3] = {
		{0.006,0.0025,0.001},
		{0.006,0.0025,0.001},
		{0.006,0.0025,0.001},
		{0.006,0.0025,0.001}
	};
	// photonID [scaled EcalIso, scaled HcalIso, 
	// scaled TrkIsoHollowd, h/e, sigmaIEtaIEta, 
	// ,etaWidth, phiWidth]
	Double_t cutEB[][7] = {
		{4.2,2.2,2.0,0.05,0.013 , 999., 999.} ,
		{2.0,2.0,1.5,0.02,0.01 , 999., 999.} ,
		{1.59, 1.5, 0.834 , 0.0196,0.01 , 999., 999.} ,
		{4.2,2.2,2.0,0.05,0.011 , 999., 999.} 
		};
	Double_t cutEE[][7] = {
		{4.2,2.2,2.0,0.05 ,0.03 ,999., 999. } ,
		{2.0,2.0,1.5,0.02 ,0.028 ,999., 999. } ,
		{0.832, 1.25, 0.887, 0.0195 ,0.028 ,999., 999. } ,
		{4.2,2.2,2.0,0.05 ,0.03 ,999., 999. } 
		};

	switch(index_gID){
		case 3: // RS-Ming
			if (BE == 0){
				trkIso  -= 0.0167 * rho ;
				ecalIso -= 0.183 * rho ;
				hcalIso -= 0.062 * rho ;
				if ( ecalIso - PhoIDScale[index_gID][0]*phoEt_[index_g] < cutEB[index_gID][0] 
				     || noPhoEcalIso                                                          ) flag++ ; else return flag ;
				if ( hcalIso - PhoIDScale[index_gID][1]*phoEt_[index_g] < cutEB[index_gID][1] 
				     || noPhoHcalIso                                                          ) flag++ ; else return flag ;
				if ( trkIso -  PhoIDScale[index_gID][2]*phoEt_[index_g] < cutEB[index_gID][2] 
				     || noPhoTrkIso                                                           ) flag++ ; else return flag ;
				if ( hovere                                             < cutEB[index_gID][3] 
				     || noPhoHoverE                                                           ) flag++ ; else return flag ;
				if ( phoSigmaIEtaIEta[index_g]                          < cutEB[index_gID][4] 
				     || noPhoSihih || doTemplate                                              ) flag++ ; else return flag ;
				return flag ;
			} else if (BE == 1){
				trkIso  -=  0.032 * rho ;
				ecalIso -=  0.090 * rho ;
				hcalIso -=  0.180 * rho ;
				if ( ecalIso - PhoIDScale[index_gID][0]*phoEt_[index_g] < cutEE[index_gID][0]
				     || noPhoEcalIso                                                          ) flag++ ; else return flag ;
				if ( hcalIso - PhoIDScale[index_gID][1]*phoEt_[index_g] < cutEE[index_gID][1]
				     || noPhoHcalIso                                                          ) flag++ ; else return flag ;
				if ( trkIso -  PhoIDScale[index_gID][2]*phoEt_[index_g] < cutEE[index_gID][2]
				     || noPhoTrkIso                                                           ) flag++ ; else return flag ;
				if ( hovere                                             < cutEE[index_gID][3]
				     || noPhoHoverE                                                           ) flag++ ; else return flag ;
				if ( phoSigmaIEtaIEta[index_g]                          < cutEE[index_gID][4]
				     || noPhoSihih || doTemplate                                              ) flag++ ; else return flag ;
				return flag ;
			}
			break ;
		case 2: // RomaIso
			if (BE == 0){
				trkIso -= 0.548 * rho ;
				ecalIso -= 0.299 * rho ;
				hcalIso -= 0.245 * rho ;
				hovere -= 0.001 * rho ;
				flag =          ( ecalIso - PhoIDScale[index_gID][0]*phoEt_[index_g] < cutEB[index_gID][0] ) ;
				flag = (flag) ? ( hcalIso - PhoIDScale[index_gID][1]*phoEt_[index_g] < cutEB[index_gID][1] ) : false ;
				flag = (flag) ? ( trkIso -  PhoIDScale[index_gID][2]*phoEt_[index_g] < cutEB[index_gID][2] ) : false ;
				flag = (flag) ? ( hovere                                            < cutEB[index_gID][3] ) : false ; 
				return flag ;
			} else if (BE == 1){
				trkIso -= 0.525 * rho ;
				ecalIso -= 0.192 * rho ;
				hcalIso -= 0.275 * rho ;
				hovere -= 0.001 * rho ;
				flag =          ( ecalIso - PhoIDScale[index_gID][0]*phoEt_[index_g] < cutEE[index_gID][0] ) ;
				flag = (flag) ? ( hcalIso - PhoIDScale[index_gID][1]*phoEt_[index_g] < cutEE[index_gID][1] ) : false ;
				flag = (flag) ? ( trkIso -  PhoIDScale[index_gID][2]*phoEt_[index_g] < cutEE[index_gID][2] ) : false ;
				flag = (flag) ? ( hovere                                            < cutEE[index_gID][3] ) : false ;
				return flag ;
			}
			break ;
		case 1: // MIT HWW
			if (BE == 0){
				if ( rho > -0.5 ){
					double fEffAreaEcalEB(0.162);
					double fEffAreaHcalEB(0.042);
					double fEffAreaTrkEB(0.317);
					ecalIso -= rho * fEffAreaEcalEB ;
					hcalIso -= rho * fEffAreaHcalEB ;
					trkIso  -= rho * fEffAreaTrkEB  ;
				}
				flag =           ( ecalIso - PhoIDScale[index_gID][0]*phoEt_[index_g] < cutEB[index_gID][0] ) ;
				flag = (flag) ?  ( hcalIso - PhoIDScale[index_gID][1]*phoEt_[index_g] < cutEB[index_gID][1] ) : false ;
				flag = (flag) ?  ( trkIso -  PhoIDScale[index_gID][2]*phoEt_[index_g] < cutEB[index_gID][2] ) : false ;
				flag = (flag) ?  ( phoHoverE[index_g]                                < cutEB[index_gID][3] ) : false ;
				flag = (flag) ?  ( phoSigmaIEtaIEta[index_g]                         < cutEB[index_gID][4] ) : false ;
			} else if (BE == 1){
				if ( rho > -0.5 ){
					double fEffAreaEcalEE(0.071);
					double fEffAreaHcalEE(0.095);
					double fEffAreaTrkEE(0.269);
					ecalIso -= rho * fEffAreaEcalEE ;
					hcalIso -= rho * fEffAreaHcalEE ;
					trkIso  -= rho * fEffAreaTrkEE  ;
				}
				flag =          ( ecalIso - PhoIDScale[index_gID][0]*phoEt_[index_g] < cutEE[index_gID][0] ) ;
				flag = (flag) ? ( hcalIso - PhoIDScale[index_gID][1]*phoEt_[index_g] < cutEE[index_gID][1] ) : false ;
				flag = (flag) ? ( trkIso -  PhoIDScale[index_gID][2]*phoEt_[index_g] < cutEE[index_gID][2] ) : false ;
				flag = (flag) ? ( phoHoverE[index_g]                                 < cutEE[index_gID][3] ) : false ;
				flag = (flag) ? ( phoSigmaIEtaIEta[index_g] 						 < cutEE[index_gID][4] ) : false ;
			}
			break ;
		case 0: // pID tight
			if (BE == 0){
				flag =          ( ecalIso - PhoIDScale[index_gID][0]*phoEt_[index_g] < cutEB[index_gID][0] ) ;
				flag = (flag) ? ( hcalIso - PhoIDScale[index_gID][1]*phoEt_[index_g] < cutEB[index_gID][1] ) : false ;
				flag = (flag) ? ( trkIso  - PhoIDScale[index_gID][2]*phoEt_[index_g] < cutEB[index_gID][2] ) : false ;
				flag = (flag) ? ( phoHoverE[index_g]                                 < cutEB[index_gID][3] ) : false ;
				flag = (flag) ? ( phoSigmaIEtaIEta[index_g]                          < cutEB[index_gID][4] ) : false ;
			} else if (BE == 1){
				flag =          ( ecalIso - PhoIDScale[index_gID][0]*phoEt_[index_g] < cutEE[index_gID][0] ) ;
				flag = (flag) ? ( hcalIso - PhoIDScale[index_gID][1]*phoEt_[index_g] < cutEE[index_gID][1] ) : false ;
				flag = (flag) ? ( trkIso  - PhoIDScale[index_gID][2]*phoEt_[index_g] < cutEE[index_gID][2] ) : false ;
				flag = (flag) ? ( phoHoverE[index_g]                                < cutEE[index_gID][3] ) : false ;
				flag = (flag) ? ( phoSigmaIEtaIEta[index_g]                         < cutEE[index_gID][4] ) : false ;
			}
			break ;
	}
	return flag ;
}

Bool_t anaVgNtuple::IsFakeable(int iPho){
	Bool_t flag = false ;
	if (phohasPixelSeed[iPho]) return false ;
	if ( fabs( phoSCEta[iPho] ) < 1.4442 ){
		flag = ( phoHoverE[iPho] < 0.05 &&
				phoSigmaIEtaIEta[iPho] < 0.014 &&
				phoTrkIsoHollowDR04[iPho] <  TMath::Min( 5*(3.5 + 0.001 * phoEt_[iPho] + 0.0167*rho), 0.2* phoEt_[iPho]) &&
				phoEcalIsoDR04[iPho] < TMath::Min( 5*(4.2+0.006* phoEt_[iPho] + 0.183*rho ), 0.2* phoEt_[iPho]) &&
				phoHcalIsoDR04[iPho] <  TMath::Min( 5*(2.2 + 0.0025 * phoEt_[iPho] + 0.062*rho), 0.2* phoEt_[iPho]) && 
				(
				 phoSigmaIEtaIEta[iPho] > 0.011 || 
				 phoTrkIsoHollowDR04[iPho] > (3.5 + 0.001 * phoEt_[iPho] + 0.0167*rho) ||
				 phoEcalIsoDR04[iPho] > ( 4.2+ 0.006 * phoEt_[iPho] + 0.183*rho) ||
				 phoHcalIsoDR04[iPho] > (2.2 + 0.0025 * phoEt_[iPho] + 0.062*rho)
				)
			   ) ;
	} else 
	if ( fabs( phoSCEta[iPho]) < 2.5 && fabs(phoSCEta[iPho]) > 1.566 ){
		flag = (  phoHoverE[iPho] < 0.05 &&
				phoSigmaIEtaIEta[iPho] <0.035 &&
				phoTrkIsoHollowDR04[iPho] <  TMath::Min( 5*(3.5 + 0.001 * phoEt_[iPho] + 0.032*rho), 0.2* phoEt_[iPho]) &&
				phoEcalIsoDR04[iPho] < TMath::Min( 5*(4.2+0.006* phoEt_[iPho] + 0.090*rho ), 0.2* phoEt_[iPho]) &&
				phoHcalIsoDR04[iPho] <  TMath::Min( 5*(2.2+0.0025* phoEt_[iPho] + 0.180*rho), 0.2* phoEt_[iPho]) && 
				(
				 phoSigmaIEtaIEta[iPho] >0.030 || 
				 phoTrkIsoHollowDR04[iPho] > (3.5 + 0.001 * phoEt_[iPho] + 0.032*rho) ||
				 phoEcalIsoDR04[iPho] > ( 4.2+ 0.006 * phoEt_[iPho] + 0.090*rho) ||
				 phoHcalIsoDR04[iPho] > (2.2 + 0.0025 * phoEt_[iPho] + 0.180*rho)
				)
			   ) ;
	}
	return flag ;
}
void anaVgNtuple::Loop(Int_t SampleIndex)
{
	Int_t nCandidate   = 0 ;
	Float_t nGenPass   = 0.  ;
	Float_t nPUCandidate = 0 ;
	Int_t nSample      = 0 ;
	Float_t cSample    = 0 ; 
	Int_t nClean       = 0 ;
	Int_t nPassHLT     = 0 ;
	Int_t nPassDiLep   = 0 ;
	Int_t nZll         = 0 ; // count number of Zll pass mass cuts
	Int_t nMatched     = 0 ; // count number of matched events
	Double_t nYield    = 0 ;
	Double_t nZllw     = 0 ; // count number of Zll pass mass cuts

	// Initialize Input options
	PUShift          = OO.call_int  (string("PUShift")        );
	PUOption         = OO.call_int  (string("PUOption")       );
	EnCorrShift      = OO.call_int  (string("EnCorrShift")    );
	FSRCorrShift     = OO.call_int  (string("FSRCorrShift")   );
	LepPtCut         = OO.call_float(string("LepPtCut")       );
	PhoHoverPreCut   = OO.call_float(string("PhoHoverPreCut") );
	PhoEtCut         = OO.call_float(string("PhoEtCut")       );
	ZMassCutL        = OO.call_float(string("ZMassCutL")      );
	ZMassCutU        = OO.call_float(string("ZMassCutU")      );
	eleScaleSysEB    = OO.call_float(string("eleScaleSysEB")  );
	eleScaleSysEE    = OO.call_float(string("eleScaleSysEE")  );
	phoScaleSysEB    = OO.call_float(string("phoScaleSysEB")  );
	phoScaleSysEE    = OO.call_float(string("phoScaleSysEE")  );
	eleResoSysEB     = OO.call_float(string("eleResoSysEB")   );
	eleResoSysEE     = OO.call_float(string("eleResoSysEE")   );
	phoResoSysEB     = OO.call_float(string("phoResoSysEB")   );
	phoResoSysEE     = OO.call_float(string("phoResoSysEE")   );
	doEle            = OO.call_bool (string("doEle")          );
	doEleTrgCheck    = OO.call_bool (string("doEleTrgCheck")  );
	doGenInfo1       = OO.call_bool (string("doGenInfo1")     );
	doMu             = OO.call_bool (string("doMu")           );
	doEMCorrection   = OO.call_bool (string("doEMCorrection") );
	doDebug          = OO.call_bool (string("doDebug")        );
	doShowList       = OO.call_bool (string("doShowList")     );
	doShowMCPho      = OO.call_bool (string("doShowMCPho")    );
	doJet            = OO.call_bool (string("doJet")          );
	doCleanEvt       = OO.call_bool (string("doCleanEvt")     );
	doPUWeight       = OO.call_bool (string("doPUWeight")     );
	doEleIDWeight    = OO.call_bool (string("doEleIDWeight")  );
	doEleRecoWeight  = OO.call_bool (string("doEleRecoWeight"));
	doEleTrgWeight   = OO.call_bool (string("doEleTrgWeight") );
	doPhoIDWeight    = OO.call_bool (string("doPhoIDWeight")  );
	doPhoPXWeight    = OO.call_bool (string("doPhoPXWeight")  );
	doRmVJetsPho     = OO.call_bool (string("doRmVJetsPho")   );
	doZJetsCorr      = OO.call_bool (string("doZJetsCorr")    );
	doPhoStandard    = OO.call_bool (string("doPhoStandard")  );
	doATGC           = OO.call_bool (string("doATGC")         );
	doTemplate       = OO.call_bool (string("doTemplate")     );
	doRatioA         = OO.call_bool (string("doRatioA")       );
	doRatioZK        = OO.call_bool (string("doRatioZK")      );
	doHistoPassPho   = OO.call_bool (string("doHistoPassPho") );
	SampleWeight     = OO.call_float(string("SampleWeight")   );
	Luminosity       = OO.call_float(string("Luminosity")     );
	noHLT            = OO.call_bool (string("noHLT")          );
	MinRun           = OO.call_int  (string("MinRun")         );
	MaxRun           = OO.call_int  (string("MaxRun")         );
	rSeed            = OO.call_int  (string("rSeed")          );
	myRandom         = new  TRandom3(rSeed);
	myRandom->SetSeed(rSeed);
	//gRandom          ->SetSeed(rSeed);
	noPhoEcalIso     = false ; // prepare for n-1 cut
	noPhoHcalIso     = false ;
	noPhoTrkIso      = false ;
	noPhoHoverE      = false ;
	noPhoSihih       = false ;
	if (OO.call_bool(string("ShowOptions"))) OO.show_options();
	// PhotonID Scale Factor
	//TFile *fEBEffSF = new TFile("File/EBEffsSF.root","READ");
	//TFile *fEEEffSF = new TFile("File/EEEffsSF.root","READ");
	//TH2F  *hEBEffSF = (TH2F*) fEBEffSF->Get("BarrelSF");
	//TH2F  *hEEEffSF = (TH2F*) fEEEffSF->Get("EndcapSF");
	// Pile-Up Re-weighting
	vector<double> puWeight ;
	TH3D*          puWeight3D = new TH3D();
	if (ProcessTag != "Data"){ // flollowing samples are using S4 weight
			if ( OO.call_string(string("PileUpFile")) != "None" ) {
				PileUpFile = OO.call_string(string("PileUpFile")) ;
			} else { // if not apply , use default one
				PileUpFile = "/afs/cern.ch/user/k/kschen/public/forVgamma/puDist_160404-173692.root" ;
			}
			if ( OO.call_string(string("PileUpFile3D")) != "None" ) {
				PileUpFile3D = OO.call_string(string("PileUpFile3D")) ;
			} else { // if not apply , use default one
				PileUpFile3D = "File/73p5mb_pudist_2011A.root" ;
			}
			TFile *fPU   = new TFile(PileUpFile  ,"READ");
			TFile *fPU3D = new TFile(PileUpFile3D,"READ");
		// S6 PU weight
		if		(	ProcessTag.Contains("Zg")     ||
					ProcessTag.Contains("ZJet")   ||
					ProcessTag.Contains("Zgtau")  ||
					ProcessTag.Contains("TTb")    ||
					ProcessTag.Contains("DiB") 
					)
		{
			puWeight3D = (TH3D*) fPU3D->Get("WHist");
			//for (int i = 1; i < 30; i++)
			//for (int j = 1; j < 30; j++)
			//for (int k = 1; k < 30; k++)
			//	cout <<"input File = " <<  i << " " << j << " " << k << " " << puWeight3D->GetBinContent(i,j,k) << endl ;
			puWeight   = generate_S6_weights((TH1D*) fPU->Get("pileup"));
			cout << "Use S6 PU weight " << endl ;
		}
		else if (	ProcessTag.Contains("W")   || 
				 	ProcessTag.Contains("Z")   ||
				 	ProcessTag.Contains("PJet")  || 
				 	ProcessTag.Contains("QCD")   
					)
		{ // S4 PU wieght
			puWeight = generate_S4_weights((TH1D*) fPU->Get("pileup"));
			cout << "Use S4 PU weight " << endl ;
		}
		else // QCD and Gamma+Jets are using S3 weight
		{
			puWeight = generate_S3_weights((TH1D*) fPU->Get("pileup"));
			cout << "Use S3 PU weight " << endl ;
		}
	fPU->Close();
	}

	// Initialize Histograms for ScaleFactor
	TFile* fSFeID ;
	TFile* fSFeRECO ; 
	TFile* fSFeHLT ;
	TFile* fSFpID ;
	TFile* fSFpPX ;
	TFile* fSFfakePho ;
	TH2D* SFeID_EB = new TH2D();
	TH2D* SFeID_EE = new TH2D();
	TH2D* SFeRECO_EB = new TH2D();
	TH2D* SFeRECO_EE = new TH2D();
	TH2D* SFeHLT_EB = new TH2D();
	TH2D* SFeHLT_EE = new TH2D();
	TH2D* SFpID_EB = new TH2D();
	TH2D* SFpID_EE = new TH2D();
	TH1D* SFpPX    = new TH1D();
	TH1D* SFfakePho_EB = new TH1D();
	TH1D* SFfakePho_EE = new TH1D();
	if (doEleIDWeight){ // x : pt , y : nGVtx
		fSFeID      = new TFile(OO.call_string(string("eID_SF_Path")).data(),"READ");
		if (!fSFeID) return ;
		SFeID_EB    = (TH2D*) fSFeID->Get("eleSF_EB");
		SFeID_EE    = (TH2D*) fSFeID->Get("eleSF_EE");
		if (!SFeID_EB)  return ;
		if (!SFeID_EE)  return ;
	}
	if (doEleRecoWeight) { // x : nGVtx , y : pt
		fSFeRECO    = new TFile(OO.call_string(string("eRECO_SF_Path")).data(),"READ");
		SFeRECO_EB  = (TH2D*) fSFeRECO->Get("BarrelSF");
		SFeRECO_EE  = (TH2D*) fSFeRECO->Get("EndcapSF");
		if (!SFeRECO_EB) return ;
		if (!SFeRECO_EE) return ;
	}
	if (doEleTrgWeight) { // x : pt , y : nGVtx
		fSFeHLT     = new TFile(OO.call_string(string("eHLT_SF_Path")).data(),"READ");
		SFeHLT_EB   = (TH2D*) fSFeHLT->Get("eleSF_EB");
		SFeHLT_EE   = (TH2D*) fSFeHLT->Get("eleSF_EE");
		if (!SFeHLT_EB) return ;
		if (!SFeHLT_EE) return ;
	}
	if (doPhoIDWeight) { // x : nGVtx , y : pt
		fSFpID     = new TFile(OO.call_string(string("pID_SF_Path")).data(),"READ");
		SFpID_EB   = (TH2D*) fSFpID->Get("BarrelSF");
		SFpID_EE   = (TH2D*) fSFpID->Get("EndcapSF");
		if (!SFpID_EB) return ;
		if (!SFpID_EE) return ;
	}
	if (doPhoPXWeight) { // x : phoSCEta
		fSFpPX     = new TFile(OO.call_string(string("pPX_SF_Path")).data(),"READ");
		SFpPX      = (TH1D*) fSFpPX->Get("SF");
		if (!SFpPX   ) return ;
	}
	if (doZJetsCorr)
	{
		fSFfakePho     = new TFile(OO.call_string(string("fPho_SF_Path")).data(),"READ");
		SFfakePho_EB   = (TH1D*) fSFfakePho->Get("BarrelSF");
		SFfakePho_EE   = (TH1D*) fSFfakePho->Get("EndcapSF");
		if (!SFfakePho_EB) return ;
		if (!SFfakePho_EE) return ;
	}

	// File create:
	TFile* fout = new TFile("NotUsed.root","RECREATE");
	TDirectory *folder = (TDirectory*) fout->Get(Form("Sample_%d",SampleIndex));
	if ( !OO.call_bool(string("useTDirectory")) 
	     && TString(OO.call_string(string("SaveFileName"))) != "None"
		 && !doATGC)
	{
		fout               = new TFile(OO.call_string(string("SaveFileName")).data(),"RECREATE");
	}
	else if (TString(OO.call_string(string("SaveFileName"))) != "None" &&
	         doATGC                                                ){
		string filePath_ = OO.call_string(string("FilePath"));
		string fileATGC = filePath_.substr(filePath_.find_last_of("/\\")+1);
		if (TString(OO.call_string(string("SaveFileName"))).Contains("2011A")){
			fileATGC    = string("2011A") + fileATGC ;
		}
		if (TString(OO.call_string(string("SaveFileName"))).Contains("2011B")){
			fileATGC    = string("2011B") + fileATGC ;
		}
		fout               = new TFile((string("File/Output_")+fileATGC).data(),"RECREATE");
		if (TString(OO.call_string(string("SaveFileName"))).Contains("Data")){
		fout               = new TFile(TString(OO.call_string(string("SaveFileName"))),"RECREATE");
		}
	}
	else if (TString(OO.call_string(string("SaveFileName"))) != "None")
	{
		fout               = new TFile(OO.call_string(string("SaveFileName")).data(),"UPDATE");
		fout->cd();
	}

	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "Create Folder if use TDir" << endl ;
	if (OO.call_bool(string("useTDirectory")))
	{
	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "try to find a folder" << endl ;
		folder = (TDirectory*) fout->Get(Form("Sample_%d",SampleIndex));
	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "try end find a folder" << endl ;
		if (folder){
	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "after try check a folder" << endl ;
			if (doDebug||OO.call_bool(string("doCheckFile")))
			cout << "found this directory, will be removed!!" << endl ;
			fout->rmdir(Form("Sample_%d",SampleIndex));
		} 
	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "mkdir" << endl ;
		folder = fout->mkdir(Form("Sample_%d",SampleIndex), "Histogram folder");
	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "cd" << endl ;
		folder->cd();
	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "endCD" << endl ;
	}

	// Initialize Histograms and Trees
	InitHisto();
	InitTree();
	Int_t eventPassLevel = 0 ;


	// time counting
	clock_t start, end ;
	start = clock();
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		// if (Cut(ientry) < 0) continue;

		// Initialize variables
		PileUpWeight = 1. ;
		LepIDWeight  = 1. ;
		PhoIDWeight  = 1. ;
		JetBkgWeight = 1. ;
		EvtWeight = SampleWeight ;
		if ( doATGC ) EvtWeight *= Luminosity ;
		if (doEMCorrection){
			phoE_  = newPhoE();
			phoEt_ = newPhoEt(phoE_) ;
		} else {
			phoE_ =  vector<Float_t>(phoE , phoE  + nPho);
			phoEt_ = vector<Float_t>(phoEt, phoEt + nPho);
		}
		// Initialize electron pt and energy
		if (doEMCorrection){
			eleEn_ = newEleEn();
			elePt_ = newElePt(eleEn_) ;
			//EcorrectionYM(TString("ele"));
		} else {
			eleEn_ = vector<Float_t>(eleEn, eleEn + nEle);
			elePt_ = vector<Float_t>(elePt, elePt + nEle);
		}
		
		// energy scale syst. if needed
		addScale(TString("ele"));
		addScale(TString("pho"));
		// energy resolution syst. if needed
		addResolution(TString("ele"));
		addResolution(TString("pho"));

		//if ( run != 166033 || event != 137575277) continue ;
		//cout << "run " << run << " event " << event << endl ;
		//cout << "entry = " << jentry << endl ;
		//cout << "nEle = " << nPho << endl ;
		//cout << "nPho = " << nEle << endl ;
		//printEvent(0,1,2);
		//if (event != 773163598) continue ;
		//if (event == 773163598) printEvent(0,1,2);

		eventPassLevel = 0 ;
		bool accept = true ;

		if (doGenInfo1)
		{
			if (!FillHistoPre()) continue ;
		}
		nGenPass++ ;
		// remove ISR and FSR for Z+Jet and W+Jet
		if ( (ProcessTag == "ZJet" || ProcessTag == "WJet") && !isData && doRmVJetsPho){
			vector<int> realPho = mcPhoton();
			if (realPho.size()!= 0) continue ;
		}
		if ( isData && run == 170722 ) continue ;
		if ( isData && run > MaxRun ) continue ;
		if ( isData && run < MinRun ) continue ;
//		if (doDebug && isData){
//			if (jentry != 29632 ) continue ;
//			if (jentry >  29632 ) break    ;
//		}

		//if (!isData)
		//{
		//	bool evtFlag = (nPU[1] >= 0 && nPU[1] < 5)  ;
		//	if (!evtFlag) continue ;
		//}

		// PU reweighting MC only
		Int_t numberPU(0) ;
		if (!isData && doPUWeight)
		{
			if	( !doATGC	&&
					(
					 ProcessTag.Contains("Zg")     ||
					 ProcessTag.Contains("ZJet")   ||
					 ProcessTag.Contains("Zgtau")  ||
					 ProcessTag.Contains("TTb")    ||
					 ProcessTag.Contains("DiB") 
					)
				)
			{
				if (doDebug) cout << "Using 3D Weight" << endl ;
				PileUpWeight = PU3D_Weight(nPU[0]+1+PUShift,nPU[1]+1+PUShift,nPU[2]+1+PUShift,puWeight3D);
				if (doDebug )cout << "PUWEIGHT = " << PileUpWeight << " " << puWeight.at(TMath::Min(nPU[1],49)) << endl ;
				//PileUpWeight = puWeight3D->GetBinContent(nPU[0],nPU[1],nPU[2]);
				//cout << "PUweight = " << PileUpWeight << endl ;
				EvtWeight   *= PileUpWeight;
			}

			else
			{
				for (int iBX = 0; iBX < nBX; iBX++) 
					if ( BXPU[iBX] == 0 ) numberPU = nPU[iBX] ;
				if (PUOption == 2 || PUOption == 3)
				{ // average reweighting
					numberPU = 0 ;
					for (int iBX = 0; iBX < nBX; iBX++) numberPU += nPU[iBX] ;
					numberPU = ( numberPU / float(nBX) ) + 0.5 ;
				}
				numberPU += PUShift ;
				if (numberPU<0) numberPU = 0 ;
				PileUpWeight = ( numberPU >= int(puWeight.size()) ) ? 0. : puWeight.at(numberPU) ;
				if (doDebug)
					cout << "numberPU1D = " <<  numberPU << " weight = " << PileUpWeight << endl ;
				EvtWeight *= PileUpWeight ;
			}
		}
		cSample+= PileUpWeight ;
		nSample++;
		// clean event
		if (doDebug) cout << "BeforeCut : run " << run << " event " << event << " entry " << jentry << endl ;
		if (doCleanEvt && !CleanEvent() ) continue ;
		nClean++ ;
		if (doDebug) cout << "Pass Clean Event " << endl ;

		// HLT
		Int_t lepTrg_index(-1);
		if ( doEle && !passHLT(lepTrg_index)   ) accept = false ;
		if ( doMu  && !passHLTMu(lepTrg_index) ) accept = false ;
		if ( noHLT )                             accept = true  ;
		if ( accept ) nPassHLT++ ;
		if (doDebug) cout << "Pass HLT " << endl ;

		// Check Electron Working Point
		if (doEleTrgCheck){
			FillTreeEleTrg();
			continue ;
		}
		
		TLorentzVector vl1, vl2, vPho , vZ ,vZg;
		Double_t MassZ(0) , MassZg(0) ;
		Int_t lep1_index = -1 , lep2_index = -1 ;
		if (doEle && accept){
			// Identify electron : 0 not pass, 1 passID, 2 passID and trigger
			vector<int> eleIDLevel = ElectronPassLevel(3,lepTrg_index);
			if (doDebug && doEle) {
				for (unsigned int iEle = 0; iEle < eleIDLevel.size(); iEle++) cout << "Electron " << iEle << " level = " << eleIDLevel.at(iEle) << endl ;
			}

			// find leading 2 electrons pass selection
			if ( !DiElePair(lep1_index,lep2_index,eleIDLevel) ) accept = false ;
			if ( accept ) nPassDiLep++ ;
			if ( accept ) eventPassLevel++;
			if (doDebug && accept) cout << "Pass Di-electron " << endl ;
			if (accept){
				vl1.SetPtEtaPhiE(elePt_[lep1_index],eleEta[lep1_index],elePhi[lep1_index],eleEn_[lep1_index]);
				vl2.SetPtEtaPhiE(elePt_[lep2_index],eleEta[lep2_index],elePhi[lep2_index],eleEn_[lep2_index]);
				vZ = vl1 + vl2 ;
				MassZ = vZ.M() ;
				if ( vZ.M() < ZMassCutL || vZ.M() > ZMassCutU ) accept = false ;
				if (accept){
					nZll++;
					if (doDebug&&accept) cout << "Pass Mass " << endl ;
					if (accept) eventPassLevel++;
					if (!isData && ProcessTag.Contains("Z")){
						if (doEleIDWeight || doEleRecoWeight || doEleTrgWeight  )
						{
							
							if ( fabs(eleSCEta[lep1_index] < 1.4442 ) )
							{
								if (doEleIDWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(float(nGoodVtx),elePt_[lep1_index],SFeID_EB);
								if (doEleRecoWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(float(nGoodVtx),elePt_[lep1_index],SFeRECO_EB);
								if (doEleTrgWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(elePt_[lep1_index],float(nGoodVtx),SFeHLT_EB);
							}
							else
							{
								if (doEleIDWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(float(nGoodVtx),elePt_[lep1_index],SFeID_EE);
								if (doEleRecoWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(float(nGoodVtx),elePt_[lep1_index],SFeRECO_EE);
								if (doEleTrgWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(elePt_[lep1_index],float(nGoodVtx),SFeHLT_EE);
							}
							if ( fabs(eleSCEta[lep2_index] < 1.4442 ) )
							{
								if (doEleIDWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(float(nGoodVtx),elePt_[lep2_index],SFeID_EB);
								if (doEleRecoWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(float(nGoodVtx),elePt_[lep2_index],SFeRECO_EB);
								if (doEleTrgWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(elePt_[lep2_index],float(nGoodVtx),SFeHLT_EB);
							}
							else
							{
								if (doEleIDWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(float(nGoodVtx),elePt_[lep2_index],SFeID_EE);
								if (doEleRecoWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(float(nGoodVtx),elePt_[lep2_index],SFeRECO_EE);
								if (doEleTrgWeight)
									LepIDWeight *=  getScaleFactorFrom2DH(elePt_[lep1_index],float(nGoodVtx),SFeHLT_EE);
							}
						}
						//if (doEleRecoWeight ) LepIDWeight *= getEleRecoWeightFromPt(eleSCEta[lep1_index],elePt_[lep1_index]);
						//if (doEleRecoWeight ) LepIDWeight *= getEleRecoWeightFromPt(eleSCEta[lep2_index],elePt_[lep2_index]);
						//if (eleTrg[lep1_index][29] == 1 )
						//	if (doEleTrgWeight  ) LepIDWeight *= getEleTrgWeightFromTrgIndex(eleSCEta[lep1_index], 29);
						//if (eleTrg[lep1_index][30] == 1 )
						//	if (doEleTrgWeight  ) LepIDWeight *= getEleTrgWeightFromTrgIndex(eleSCEta[lep1_index], 30);
						//if (eleTrg[lep2_index][29] == 1 )
						//	if (doEleTrgWeight  ) LepIDWeight *= getEleTrgWeightFromTrgIndex(eleSCEta[lep2_index], 29);
						//if (eleTrg[lep2_index][30] == 1 )
						//	if (doEleTrgWeight  ) LepIDWeight *= getEleTrgWeightFromTrgIndex(eleSCEta[lep2_index], 30);
						EvtWeight *= LepIDWeight ;
					}
					FillHistoZ(lep1_index,lep2_index) ;
					nZllw+=EvtWeight ;
				}
			}
		}

		// find 2 muons pass selection
		if (doMu&&accept){
		if (doDebug) cout << "Start DiMuon selection " << endl ;
			if ( !DiMuPair(lep1_index,lep2_index,lepTrg_index) ) accept = false ;
				if (accept) eventPassLevel++;
				if (accept){
			float muEn1 = sqrt(muPt[lep1_index]*muPt[lep1_index] + muPz[lep1_index]*muPz[lep1_index]);
			float muEn2 = sqrt(muPt[lep2_index]*muPt[lep2_index] + muPz[lep2_index]*muPz[lep2_index]);
				vl1.SetPtEtaPhiE(muPt[lep1_index],muEta[lep1_index],muPhi[lep1_index],muEn1);
				vl2.SetPtEtaPhiE(muPt[lep2_index],muEta[lep2_index],muPhi[lep2_index],muEn2);
				vZ = vl1 + vl2 ;
				MassZ = vZ.M() ;
				if (doDebug) cout << "pass di muon " << lep1_index << " " << lep2_index <<  " Zmass = " << MassZ << endl ;
				if ( vZ.M() < ZMassCutL || vZ.M() > ZMassCutU ) accept = false ;
				if (accept) eventPassLevel++;
				if (doDebug) cout << "pass Zmass " << endl ;
				}
		}

		// Fill EventFlow
		if (doHistoPassPho){
			if (accept) eventPassLevel += PhoSelectLevel(lep1_index,lep2_index);
			for (int iF = 0; iF < eventPassLevel; iF++) histo[100]->Fill(float(0.5+iF),EvtWeight);
			if (accept) // fill number of vertices and NPU Passing Z selection
			{
				histo[ 21]->Fill(nGoodVtx , EvtWeight);
				histo[ 41]->Fill(nVtx     , EvtWeight);
				if (!isData && doPUWeight){
					histo[ 62]->Fill(numberPU, SampleWeight);
					histo[ 63]->Fill(numberPU, SampleWeight*PileUpWeight);
				}
			}
		}
		if (!accept) continue ;

		// Photon Selection
		
		Int_t pho_index(-1) ;
		pho_index = PhoSelection(lep1_index,lep2_index);

		if (doATGC && pho_index!=-1){
			vPho.SetPtEtaPhiE(phoEt_[pho_index],phoEta[pho_index],phoPhi[pho_index],phoE_[pho_index]);
			vZg = vZ + vPho ;
			if (vZ.M()+vZg.M()<185) pho_index = -1 ;
		}
		if ( pho_index != -1 ){
			if (doDebug) cout << "Pass Photon selection " << endl ;
			nCandidate++ ;
			nPUCandidate+=PileUpWeight ;
			// photon selection correction
			if (!isData && ProcessTag == "Zg" ){
				if (doPhoIDWeight)
				{
//					PhoIDWeight = getPhoIDWeightFromEta(phoSCEta[pho_index]) ;
//					PhoIDWeight = getPhoIDWeightFromPtVtx(phoEt_[pho_index],nGoodVtx,phoSCEta[pho_index],hEBEffSF,hEEEffSF) ;
					if ( fabs(phoSCEta[pho_index] < 1.4442 ) )
					{
						PhoIDWeight *=  getScaleFactorFrom2DH(float(nGoodVtx),elePt_[pho_index],SFpID_EB);
					}
					else
					{
						PhoIDWeight *=  getScaleFactorFrom2DH(float(nGoodVtx),elePt_[pho_index],SFpID_EE);
					}
				}
				if (doPhoPXWeight)
				{
					PhoIDWeight *= getScaleFactorFrom1DH(fabs(phoSCEta[pho_index]),SFpPX);
				}
				EvtWeight *= PhoIDWeight ;
			}
			if (doZJetsCorr){
				if (!mcMuMatcher(pho_index,"pho" ) 
						&& !mcEleMatcher(pho_index,"pho") 
						&& !mcPhoMatcher(pho_index)     )
				{
					if (fabs(phoSCEta[pho_index]) < 1.4442)
					{
						EvtWeight *= getScaleFactorFrom1DH(phoEt_[pho_index],SFfakePho_EB);
					}
					else
					{
						EvtWeight *= getScaleFactorFrom1DH(phoEt_[pho_index],SFfakePho_EE);
					}
				}
			}
			nYield+=EvtWeight ;
			if (doShowList){
				vPho.SetPtEtaPhiE(phoEt_[pho_index],phoEta[pho_index],phoPhi[pho_index],phoE_[pho_index]);
				vZg    = vZ + vPho ;
				MassZg = vZg.M();
				cout<<"#Zgamma:"<<setw(5)<<nCandidate<<setw(8);
				cout<<"run:"<<setw(7)<<run<<setw(10)<<"event:"<<setw(12);
				cout<<event<<setw(12)<<"lumis:"<<setw(5)<<lumis ;
				cout << endl ;
				//cout<<setw(7)<<"MassZ:"<< setw(10) <<MassZ ;
				//cout<<setw(8)<<"MassZg:"<< setw(10) <<MassZg ;
				//cout<<" pho:ele1:ele2 " << pho_index << ":" << lep1_index << ":" << lep2_index ;
				//cout<<" prescale: " << t_HLTprescale << endl;
				if (doEle && doDebug ) { printElectron(lep1_index); printElectron(lep2_index);}
				
			}
			if (doDebug){
				cout << "Entry = " << jentry << " ele1 = " << lep1_index <<  " ele2 = " << lep2_index << " pho = " << pho_index <<endl ;
				if (doDebug && doEle) printEvent(lep1_index,lep2_index,pho_index);
			}
			if (doHistoPassPho||doTemplate)
			{
				//cout << "GenMomPID = " << phoGenMomPID[pho_index] << " GenGMomPID = " << phoGenGMomPID[pho_index] << endl ;
				FillHisto(lep1_index, lep2_index, pho_index);
				FillTree(lep1_index, lep2_index, pho_index);
				getNPassJets();
				if (!isData && doPUWeight){
					histo[ 64]->Fill(numberPU, SampleWeight);
					histo[ 65]->Fill(numberPU, EvtWeight);
				}
			}
			if (doATGC){
				if (ProcessTag.Contains("ZJet")){
					EvtWeight *= getScaleFactorFromRatio(phoSCEta[pho_index]);
				}
				FillTree(lep1_index, lep2_index, pho_index);
			}
			// electron trigger check:
			//if (doEleTrgWeight){
			//	bool lepTrgLeg1 = checkHLT(lep1_index,lepTrg_index);
			//	bool lepTrgLeg2 = checkHLT(lep2_index,lepTrg_index);
			//	if (lepTrgLeg1&&lepTrgLeg2){
			//		nMatched++;
			//	} 
			//	else 
			//	if (!lepTrgLeg1)
			//	{
			//		cout << "Event " << event << " run " << run << "fail match : " << endl ;
			//		printElectron(lep1_index);
			//	}
			//	else
			//	if (!lepTrgLeg2)
			//	{
			//		cout << "Event " << event << " run " << run << "fail match : " << endl ;
			//		printElectron(lep2_index);
			//	}
			//}
		}
		
	}// loop over events
	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "Endloop" << endl ;
	cout << "Procss: " << ProcessTag <<"_"<< SampleIndex << endl ;
	cout << "Sample Coutnt  = " << nSample      << endl ;
	cout << "Sample PUC     = " << cSample      << endl ;
	cout << "Pass vtx trk   = " << nClean       << endl ;
	cout << "Pass   HLT     = " << nPassHLT     << " eff = " << float(nPassHLT)/nSample << endl ;
	cout << "Passed DiLep   = " << nPassDiLep   << " eff = " << float(nPassDiLep)/float(nPassHLT) << endl ;
	cout << "Passed Zs      = " << nZll         << " eff = " << float(nZll)/float(nPassDiLep) << endl ;
	cout << "Passed ZsW     = " << nZllw * Luminosity        << endl ;
	cout << "Passed events  = " << nCandidate   << " eff = " << float(nCandidate)/float(nZll) 
	                            << " " << float(nCandidate)/float(nSample) << endl ;
	cout << "Passed PU evt  = " << nPUCandidate << endl ;
	cout << "Passed Yields  = " << nYield       << endl ;
	if (doEleTrgWeight)
	cout << "Matched events = " << nMatched << endl ;

	// Write histogram and tree to TFile then close file
	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "Start to write histograms" << endl ;
	WrHisto();
	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "Start to write trees" << endl ;
	WrTree();
//	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "Start to remove histograms" << endl ;
//	RmHisto();
	if (doDebug||OO.call_bool(string("doCheckFile"))) cout << "Start to close file" << endl ;
	if (fout) fout->Close();

	end = clock();
	double dif = (double)( end - start) / CLOCKS_PER_SEC ;
	cout << "Execute time = " << dif << " s" << endl ;
}
