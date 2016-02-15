//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 22 14:52:01 2011 by ROOT version 5.27/06b
// from TChain VgAnalyzerKit/EventTree/
//////////////////////////////////////////////////////////

#ifndef anaVgNtuple_h
#define anaVgNtuple_h

#include "Header/useoption.h"
#include <TROOT.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH3D.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>

   const Int_t kMaxmuValidFraction = 1;
   const Int_t kMaxTree            = 10  ;
   const Int_t kMaxHiso1D          = 200 ;
   const Int_t kMaxHiso2D          = 20  ;

class anaVgNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TRandom3       *myRandom ; // Random Seed

   // Bitwise cut
   enum { // electron trigger working points
     _eTrg_CaloIdVL           , // 00
     _eTrg_CaloIdL            , // 01
     _eTrg_CaloIdXL           , // 02
     _eTrg_CaloIdT            , // 03
     _eTrg_CaloIdVT           , // 04
     _eTrg_CaloIsoVL          , // 05
     _eTrg_CaloIsoL           , // 06
     _eTrg_CaloIsoXL          , // 07
     _eTrg_CaloIsoT           , // 08
     _eTrg_CaloIsoVT          , // 09
     _eTrg_TrkIdVL            , // 10
     _eTrg_TrkIdL             , // 11
     _eTrg_TrkIdXL            , // 12
     _eTrg_TrkIdT             , // 13
     _eTrg_TrkIdVT            , // 14
     _eTrg_TrkIsoVL           , // 15
     _eTrg_TrkIsoL            , // 16
     _eTrg_TrkIsoXL           , // 17
     _eTrg_TrkIsoT            , // 18
     _eTrg_TrkIsoVT           , // 19
	 N_ELE_TRG_WP_CATOGORY      // 20
   };
   // electronIDCut
   enum {
   	_ele_PVtx    , // 00
	_ele_Conv    , // 01
	_ele_Iso     , // 02
	_ele_ID      , // 03
	N_ELE_ISOID    // 04
   };



   // Initial cut value and options
   useoption      OO              ;
   Float_t        PhoHoverPreCut  ;
   Float_t        LepPtCut        ;
   Float_t        PhoEtCut        ;
   Float_t        ZMassCutL       ;
   Float_t        ZMassCutU       ;
   Float_t        eleScaleSysEB   ;
   Float_t        eleScaleSysEE   ;
   Float_t        phoScaleSysEB   ;
   Float_t        phoScaleSysEE   ;
   Float_t        eleResoSysEB    ;
   Float_t        eleResoSysEE    ;
   Float_t        phoResoSysEB    ;
   Float_t        phoResoSysEE    ;
   Bool_t         doEle           ;
   Bool_t         doMu            ;
   Bool_t         doEleTrgCheck   ;
   Bool_t         doGenInfo1      ;
   Bool_t         doEMCorrection  ;
   Bool_t         doDebug         ;
   Bool_t         doShowList      ;
   Bool_t         doShowMCPho     ;
   Bool_t         doJet           ;
   Bool_t         doCleanEvt      ;
   Bool_t         doPUWeight      ;
   Bool_t         doEleIDWeight   ;
   Bool_t         doEleRecoWeight ;
   Bool_t         doEleTrgWeight  ;
   Bool_t         doPhoIDWeight   ;
   Bool_t         doPhoPXWeight   ;
   Bool_t         doZJetsCorr     ;
   Bool_t         doATGC          ;
   Bool_t         doTemplate      ;
   Bool_t         doPhoStandard   ;
   Bool_t         doRatioA        ;
   Bool_t         doRatioZK       ;
   Bool_t         doHistoPassPho  ;
   Bool_t         doRmVJetsPho    ;
   Bool_t         noPhoEcalIso    ;
   Bool_t         noPhoHcalIso    ;
   Bool_t         noPhoTrkIso     ;
   Bool_t         noPhoHoverE     ;
   Bool_t         noPhoSihih      ;
   Bool_t         noHLT           ;
   Bool_t         useTDirectory   ;
   Int_t          PUShift         ;
   Int_t          PUOption        ;
   Int_t          EnCorrShift     ;
   Int_t          FSRCorrShift    ;
   TString        ProcessTag      ;
   TString        SaveFileName    ;
   Double_t       SampleWeight    ; 
   Double_t       PileUpWeight    ; 
   Double_t       LepIDWeight     ; 
   Double_t       PhoIDWeight     ; 
   Double_t       JetBkgWeight    ; 
   Double_t       EvtWeight       ;
   Double_t       Luminosity      ;
   Int_t          MaxRun          ;
   Int_t          MinRun          ;
   Int_t          rSeed           ;
   TString        PileUpFile      ;
   TString        PileUpFile3D    ;
   TString        eID_SF_Path     ;
   TString        eHLT_SF_Path    ;
   TString        eRECO_SF_Path   ;
   TString        pID_SF_Path     ;
   TString        pPX_SF_Path     ;
   TString        fPho_SF_Path    ;
   vector<string> InputPath       ;

   // Values for prepare for tree
   Int_t    t_HLTprescale      ;
   Double_t t_phoEt            ; // photon
   Double_t t_phoEta           ; 
   Double_t t_phoSCEta         ; 
   Double_t t_phoPhi           ;
   Double_t t_phoE             ;
   Double_t t_phoSigmaIEtaIEta ;
   Double_t t_phoSCE           ;
   Bool_t   t_phoInEB          ;
   Bool_t   t_phoInEE          ;
   Bool_t   t_phoIsISR         ;
   Bool_t   t_phoIsFSR         ;
   Bool_t   t_phoIsEle         ;
   Bool_t   t_phoIsPho         ;
   Double_t t_ele1Pt           ; // electorn1
   Double_t t_ele1Eta          ;
   Double_t t_ele1SCEta        ;
   Double_t t_ele1Phi          ;
   Double_t t_ele1En           ;
   Bool_t   t_ele1InEB         ;
   Bool_t   t_ele1InEE         ;
   Double_t t_ele2Pt           ; // electorn2
   Double_t t_ele2Eta          ;
   Double_t t_ele2SCEta        ;
   Double_t t_ele2Phi          ;
   Double_t t_ele2En           ;
   Bool_t   t_ele2InEB         ;
   Bool_t   t_ele2InEE         ;
   Double_t t_mu1Pt            ; // muon1 
   Double_t t_mu1Eta           ;
   Double_t t_mu1Phi           ;
   Double_t t_mu1En            ;
   Double_t t_mu2Pt            ; // muon2 
   Double_t t_mu2Eta           ;
   Double_t t_mu2Phi           ;
   Double_t t_mu2En            ;
   Double_t t_massZ            ; // Z mass
   Double_t t_massZg           ; // Zgamma mass
   Double_t t_massZZg          ; // Zgamma+Z mass
   Double_t t_ZPt              ;
   Double_t t_ZgST             ;
   Bool_t   t_isHWW_WP95       ;
   Bool_t   t_isHWW_WP90       ;
   Bool_t   t_isHWW_WP85       ;
   Bool_t   t_isHWW_WP70       ;
   Bool_t   t_isHWW_WP80       ;
   Bool_t   t_isHWW_WP60       ;
   Bool_t   t_eleTrgLeg1Match  ;
   Bool_t   t_eleTrgLeg2Match  ;
   Bool_t   t_CaloIdVL         ; 
   Bool_t   t_CaloIdL          ; 
   Bool_t   t_CaloIdXL         ; 
   Bool_t   t_CaloIdT          ; 
   Bool_t   t_CaloIdVT         ; 
   Bool_t   t_CaloIsoVL        ; 
   Bool_t   t_CaloIsoL         ; 
   Bool_t   t_CaloIsoXL        ; 
   Bool_t   t_CaloIsoT         ; 
   Bool_t   t_CaloIsoVT        ; 
   Bool_t   t_TrkIdVL          ; 
   Bool_t   t_TrkIdL           ; 
   Bool_t   t_TrkIdXL          ; 
   Bool_t   t_TrkIdT           ; 
   Bool_t   t_TrkIdVT          ; 
   Bool_t   t_TrkIsoVL         ; 
   Bool_t   t_TrkIsoL          ; 
   Bool_t   t_TrkIsoXL         ; 
   Bool_t   t_TrkIsoT          ; 
   Bool_t   t_TrkIsoVT         ; 


   // Declaration of variables
   vector<Float_t> phoEt_ ;
   vector<Float_t> phoE_ ;
   vector<Float_t> eleEn_ ;
   vector<Float_t> elePt_ ;

   //histograms and trees to be stored;
   TTree* tree_[kMaxTree] ;
   TH1D* histo[kMaxHiso1D] ;
   TH2D* his2D[kMaxHiso2D] ;

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           orbit;
   Int_t           bx;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           ttbit0;
   Int_t           nHLT;
   Int_t           HLT[500];   //[nHLT]
   Int_t           HLTIndex[125];
   Int_t           HLTprescale[500];   //[nHLT]
   Int_t           nHFTowersP;
   Int_t           nHFTowersN;
   Int_t           nVtx;
   Float_t         vtx[70][3];   //[nVtx]
   Int_t           vtxNTrk[70];   //[nVtx]
   Float_t         vtxNDF[70];   //[nVtx]
   Float_t         vtxD0[70];   //[nVtx]
   Int_t           nGoodVtx;
   Int_t           IsVtxGood;
   Int_t           nTrk;
   Int_t           nGoodTrk;
   Int_t           IsTracksGood;
   Float_t         rho;
   Float_t         sigma;
   Float_t         rhoNeutral;
   Float_t         pdf[7];
   Float_t         pthat;
   Float_t         processID;
   Int_t           nBX;
   Int_t           nPU[3];   //[nBX]
   Int_t           BXPU[3];   //[nBX]
   Int_t           nMC;
   Int_t           mcPID[70];   //[nMC]
   Float_t         mcVtx[70][3];   //[nMC]
   Float_t         mcPt[70];   //[nMC]
   Float_t         mcMass[70];   //[nMC]
   Float_t         mcEta[70];   //[nMC]
   Float_t         mcPhi[70];   //[nMC]
   Int_t           mcGMomPID[70];   //[nMC]
   Int_t           mcMomPID[70];   //[nMC]
   Float_t         mcMomPt[70];   //[nMC]
   Float_t         mcMomMass[70];   //[nMC]
   Float_t         mcMomEta[70];   //[nMC]
   Float_t         mcMomPhi[70];   //[nMC]
   Int_t           mcIndex[70];   //[nMC]
   Int_t           mcDecayType[70];   //[nMC]
   Float_t         mcIsoDR03[70];   //[nMC]
   Float_t         mcCalIsoDR03[70];   //[nMC]
   Float_t         mcTrkIsoDR03[70];   //[nMC]
   Float_t         mcIsoDR04[70];   //[nMC]
   Float_t         mcCalIsoDR04[70];   //[nMC]
   Float_t         mcTrkIsoDR04[70];   //[nMC]
   Float_t         genMET;
   Float_t         genMETx;
   Float_t         genMETy;
   Float_t         genMETPhi;
   Float_t         MET;
   Float_t         METx;
   Float_t         METy;
   Float_t         METPhi;
   Float_t         METsumEt;
   Float_t         uncorrMET[3];
   Float_t         uncorrMETPhi[3];
   Float_t         uncorrMETSumEt[3];
   Float_t         tcMET;
   Float_t         tcMETx;
   Float_t         tcMETy;
   Float_t         tcMETPhi;
   Float_t         tcMETsumEt;
   Float_t         tcMETmEtSig;
   Float_t         tcMETSig;
   Float_t         pfMET;
   Float_t         pfMETx;
   Float_t         pfMETy;
   Float_t         pfMETPhi;
   Float_t         pfMETsumEt;
   Float_t         pfMETmEtSig;
   Float_t         pfMETSig;
   Float_t         TypeIpfMET;
   Float_t         TypeIpfMETx;
   Float_t         TypeIpfMETy;
   Float_t         TypeIpfMETPhi;
   Float_t         TypeIpfMETsumEt;
   Float_t         TypeIpfMETmEtSig;
   Float_t         TypeIpfMETSig;
   Float_t         TypeIpIIpfMET;
   Float_t         TypeIpIIpfMETx;
   Float_t         TypeIpIIpfMETy;
   Float_t         TypeIpIIpfMETPhi;
   Float_t         TypeIpIIpfMETsumEt;
   Float_t         TypeIpIIpfMETmEtSig;
   Float_t         TypeIpIIpfMETSig;
   Float_t         SmearedpfMET;
   Float_t         SmearedpfMETx;
   Float_t         SmearedpfMETy;
   Float_t         SmearedpfMETPhi;
   Float_t         SmearedpfMETsumEt;
   Float_t         SmearedpfMETmEtSig;
   Float_t         SmearedpfMETSig;
   Float_t         SmearedTypeIpfMET;
   Float_t         SmearedTypeIpfMETx;
   Float_t         SmearedTypeIpfMETy;
   Float_t         SmearedTypeIpfMETPhi;
   Float_t         SmearedTypeIpfMETsumEt;
   Float_t         SmearedTypeIpfMETmEtSig;
   Float_t         SmearedTypeIpfMETSig;
   Int_t           npfCharged;
   Float_t         pfChargedSumPt;
   Int_t           npfChargedHadron;
   Float_t         pfChargedHadronSumPt;
   Int_t           npfLepton;
   Float_t         pfLeptonSumPt;
   Int_t           npfNeutral;
   Float_t         pfNeutralSumPt;
   Int_t           npfNeutralHadron;
   Float_t         pfNeutralHadronSumPt;
   Int_t           npfPhoton;
   Float_t         pfPhotonSumPt;
   Int_t           nEle;
   Int_t           eleTrg[80][31];   //[nEle]
   Int_t           eleID[80][30];   //[nEle]
   Float_t         eleIDLH[80];   //[nEle]
   Int_t           eleClass[80];   //[nEle]
   Int_t           eleCharge[80];   //[nEle]
   Float_t         eleEn[80];   //[nEle]
   Float_t         eleSCRawEn[80];   //[nEle]
   Float_t         eleESEn[80];   //[nEle]
   Float_t         eleSCEn[80];   //[nEle]
   Float_t         elePt[80];   //[nEle]
   Float_t         elePz[80];   //[nEle]
   Float_t         eleEta[80];   //[nEle]
   Float_t         elePhi[80];   //[nEle]
   Float_t         eleSCEta[80];   //[nEle]
   Float_t         eleSCPhi[80];   //[nEle]
   Float_t         eleSCEtaWidth[80];   //[nEle]
   Float_t         eleSCPhiWidth[80];   //[nEle]
   Float_t         eleVtx[80][3];   //[nEle]
   Float_t         eleCaloPos[80][3];   //[nEle]
   Float_t         eleSCPos[80][3];   //[nEle]
   Float_t         eleHoverE[80];   //[nEle]
   Float_t         eleEoverP[80];   //[nEle]
   Float_t         elePin[80];   //[nEle]
   Float_t         elePout[80];   //[nEle]
   Float_t         eleBrem[80];   //[nEle]
   Int_t           elenBrem[80];   //[nEle]
   Float_t         eledEtaAtVtx[80];   //[nEle]
   Float_t         eledPhiAtVtx[80];   //[nEle]
   Float_t         eleSigmaEtaEta[80];   //[nEle]
   Float_t         eleSigmaIEtaIEta[80];   //[nEle]
   Float_t         eleSigmaIEtaIPhi[80];   //[nEle]
   Float_t         eleSigmaIPhiIPhi[80];   //[nEle]
   Float_t         eleE3x3[80];   //[nEle]
   Float_t         eleSeedTime[80];   //[nEle]
   Float_t         eleSeedEnergy[80];   //[nEle]
   Int_t           eleRecoFlag[80];   //[nEle]
   Int_t           eleSeverity[80];   //[nEle]
   Int_t           eleGenIndex[80];   //[nEle]
   Int_t           eleGenGMomPID[80];   //[nEle]
   Int_t           eleGenMomPID[80];   //[nEle]
   Float_t         eleGenMomPt[80];   //[nEle]
   Float_t         eleIsoTrkDR03[80];   //[nEle]
   Float_t         eleIsoEcalDR03[80];   //[nEle]
   Float_t         eleIsoHcalDR03[80];   //[nEle]
   Float_t         eleIsoHcalSolidDR03[80];   //[nEle]
   Float_t         eleIsoTrkDR04[80];   //[nEle]
   Float_t         eleIsoEcalDR04[80];   //[nEle]
   Float_t         eleIsoHcalDR04[80];   //[nEle]
   Float_t         eleIsoHcalSolidDR04[80];   //[nEle]
   Float_t         eleConvDist[80];   //[nEle]
   Float_t         eleConvDcot[80];   //[nEle]
   Float_t         eleConvRadius[80];   //[nEle]
   Int_t           eleConvFlag[80];   //[nEle]
   Int_t           eleConvMissinghit[80];   //[nEle]
   Float_t         eleESRatio[80];   //[nEle]
   Float_t         eleESProfileFront[80][123];   //[nEle]
   Float_t         eleESProfileRear[80][123];   //[nEle]
   Float_t         elePV2D[80];   //[nEle]
   Float_t         elePV3D[80];   //[nEle]
   Float_t         eleBS2D[80];   //[nEle]
   Float_t         eleBS3D[80];   //[nEle]
   Float_t         elePVD0[80];   //[nEle]
   Float_t         elePVDz[80];   //[nEle]
   Int_t           nPho;
   Int_t           phoTrg[150][14];   //[nPho]
   Bool_t          phoIsPhoton[150];   //[nPho]
   Float_t         phoE[150];   //[nPho]
   Float_t         phoEt[150];   //[nPho]
   Float_t         phoPz[150];   //[nPho]
   Float_t         phoEta[150];   //[nPho]
   Float_t         phoPhi[150];   //[nPho]
   Float_t         phoR9[150];   //[nPho]
   Float_t         phoTrkIsoSolidDR03[150];   //[nPho]
   Float_t         phoTrkIsoHollowDR03[150];   //[nPho]
   Int_t           phoNTrkSolidDR03[150];   //[nPho]
   Int_t           phoNTrkHollowDR03[150];   //[nPho]
   Float_t         phoEcalIsoDR03[150];   //[nPho]
   Float_t         phoHcalIsoDR03[150];   //[nPho]
   Float_t         phoHcalIsoSolidDR03[150];   //[nPho]
   Float_t         phoTrkIsoSolidDR04[150];   //[nPho]
   Float_t         phoTrkIsoHollowDR04[150];   //[nPho]
   Int_t           phoNTrkSolidDR04[150];   //[nPho]
   Int_t           phoNTrkHollowDR04[150];   //[nPho]
   Float_t         phoEcalIsoDR04[150];   //[nPho]
   Float_t         phoHcalIsoDR04[150];   //[nPho]
   Float_t         phoHcalIsoSolidDR04[150];   //[nPho]
   Float_t         phoEtVtx[150][150];   //[nPho]
   Float_t         phoEtaVtx[150][150];   //[nPho]
   Float_t         phoPhiVtx[150][150];   //[nPho]
   Float_t         phoTrkIsoSolidDR03Vtx[150][150];   //[nPho]
   Float_t         phoTrkIsoHollowDR03Vtx[150][150];   //[nPho]
   Float_t         phoTrkIsoSolidDR04Vtx[150][150];   //[nPho]
   Float_t         phoTrkIsoHollowDR04Vtx[150][150];   //[nPho]
   Float_t         phoHoverE[150];   //[nPho]
   Float_t         phoSigmaEtaEta[150];   //[nPho]
   Float_t         phoSigmaIEtaIEta[150];   //[nPho]
   Float_t         phoSigmaIEtaIPhi[150];   //[nPho]
   Float_t         phoSigmaIPhiIPhi[150];   //[nPho]
   Float_t         phoE3x3[150];   //[nPho]
   Float_t         phoE5x5[150];   //[nPho]
   Float_t         phoSeedTime[150];   //[nPho]
   Float_t         phoSeedEnergy[150];   //[nPho]
   Int_t           phoRecoFlag[150];   //[nPho]
   Int_t           phoSeverity[150];   //[nPho]
   Int_t           phoPos[150];   //[nPho]
   Int_t           phoGenIndex[150];   //[nPho]
   Int_t           phoGenGMomPID[150];   //[nPho]
   Int_t           phoGenMomPID[150];   //[nPho]
   Float_t         phoGenMomPt[150];   //[nPho]
   Float_t         phoSCE[150];   //[nPho]
   Float_t         phoESE[150];   //[nPho]
   Float_t         phoSCEt[150];   //[nPho]
   Float_t         phoSCEta[150];   //[nPho]
   Float_t         phoSCPhi[150];   //[nPho]
   Float_t         phoSCEtaWidth[150];   //[nPho]
   Float_t         phoSCPhiWidth[150];   //[nPho]
   Float_t         phoVtx[150][3];   //[nPho]
   Float_t         phoVtxD0[150];   //[nPho]
   Int_t           phoOverlap[150];   //[nPho]
   Int_t           phohasPixelSeed[150];   //[nPho]
   Int_t           phoIsConv[150];   //[nPho]
   Float_t         phoESRatio[150];   //[nPho]
   Float_t         phoESProfileFront[150][123];   //[nPho]
   Float_t         phoESProfileRear[150][123];   //[nPho]
   Int_t           phoNTracks[150];   //[nPho]
   Float_t         phoConvPairInvariantMass[150];   //[nPho]
   Float_t         phoConvPairCotThetaSeparation[150];   //[nPho]
   Float_t         phoConvPairMomentumEta[150];   //[nPho]
   Float_t         phoConvPairMomentumPhi[150];   //[nPho]
   Float_t         phoConvPairMomentumX[150];   //[nPho]
   Float_t         phoConvPairMomentumY[150];   //[nPho]
   Float_t         phoConvPairMomentumZ[150];   //[nPho]
   Float_t         phoConvDistOfMinimumApproach[150];   //[nPho]
   Float_t         phoConvDPhiTracksAtVtx[150];   //[nPho]
   Float_t         phoConvDPhiTracksAtEcal[150];   //[nPho]
   Float_t         phoConvDEtaTracksAtEcal[150];   //[nPho]
   Float_t         phoConvVtxValid[150];   //[nPho]
   Float_t         phoConvVtxEta[150];   //[nPho]
   Float_t         phoConvVtxPhi[150];   //[nPho]
   Float_t         phoConvVtxR[150];   //[nPho]
   Float_t         phoConvVtxX[150];   //[nPho]
   Float_t         phoConvVtxY[150];   //[nPho]
   Float_t         phoConvVtxZ[150];   //[nPho]
   Float_t         phoConvVtxChi2[150];   //[nPho]
   Float_t         phoConvVtxNdof[150];   //[nPho]
   Float_t         phoConvChi2Prob[150];   //[nPho]
   Float_t         phoConvEoverP[150];   //[nPho]
   Int_t           phoNxtal[150];   //[nPho]
   Float_t         phoXtalTime[150][200];   //[nPho]
   Float_t         phoXtalEnergy[150][200];   //[nPho]
   Int_t           phoXtalZ[150][200];   //[nPho]
   Int_t           phoXtalX[150][200];   //[nPho]
   Int_t           phoXtalY[150][200];   //[nPho]
   Int_t           phoXtalEta[150][200];   //[nPho]
   Int_t           phoXtalPhi[150][200];   //[nPho]
   Float_t         pho5x5Time[150][25];   //[nPho]
   Float_t         pho5x5Energy[150][25];   //[nPho]
   Int_t           pho5x5Z[150][25];   //[nPho]
   Int_t           pho5x5X[150][25];   //[nPho]
   Int_t           pho5x5Y[150][25];   //[nPho]
   Int_t           pho5x5Eta[150][25];   //[nPho]
   Int_t           pho5x5Phi[150][25];   //[nPho]
   Int_t           nMu;
   Int_t           muTrg[180][16];   //[nMu]
   Float_t         muEta[180];   //[nMu]
   Float_t         muPhi[180];   //[nMu]
   Int_t           muCharge[180];   //[nMu]
   Float_t         muPt[180];   //[nMu]
   Float_t         muPz[180];   //[nMu]
   Int_t           muGenIndex[180];   //[nMu]
   Int_t           muGenGMomPID[180];   //[nMu]
   Int_t           muGenMomPID[180];   //[nMu]
   Float_t         muGenMomPt[180];   //[nMu]
   Float_t         muIsoTrk[180];   //[nMu]
   Float_t         muIsoCalo[180];   //[nMu]
   Float_t         muIsoEcal[180];   //[nMu]
   Float_t         muIsoHcal[180];   //[nMu]
   Float_t         muChi2NDF[180];   //[nMu]
   Float_t         muEmVeto[180];   //[nMu]
   Float_t         muHadVeto[180];   //[nMu]
   Int_t           muType[180];   //[nMu]
   Bool_t          muID[180][6];   //[nMu]
   Float_t         muD0[180];   //[nMu]
   Float_t         muDz[180];   //[nMu]
   Float_t         muPVD0[180];   //[nMu]
   Float_t         muPVDz[180];   //[nMu]
   Float_t         muValidFraction[180];   //[nMu]
   Float_t         muTrkdPt[180];   //[nMu]
   Int_t           muNumberOfHits[180];   //[nMu]
   Int_t           muNumberOfValidHits[180];   //[nMu]
   Int_t           muNumberOfInactiveHits[180];   //[nMu]
   Int_t           muNumberOfValidTrkHits[180];   //[nMu]
   Int_t           muNumberOfValidPixelHits[180];   //[nMu]
   Int_t           muNumberOfValidMuonHits[180];   //[nMu]
   Int_t           muStations[180];   //[nMu]
   Int_t           muChambers[180];   //[nMu]
   Float_t         muPV2D[180];   //[nMu]
   Float_t         muPV3D[180];   //[nMu]
   Float_t         muBS2D[180];   //[nMu]
   Float_t         muBS3D[180];   //[nMu]
   Float_t         muVtx[180][3];   //[nMu]
   Int_t           nJet;
   Int_t           jetTrg[180][23];   //[nJet]
   Float_t         jetEn[180];   //[nJet]
   Float_t         jetPt[180];   //[nJet]
   Float_t         jetEta[180];   //[nJet]
   Float_t         jetPhi[180];   //[nJet]
   Float_t         jetMass[180];   //[nJet]
   Float_t         jetEt[180];   //[nJet]
   Int_t           jetpartonFlavour[180];   //[nJet]
   Float_t         jetRawPt[180];   //[nJet]
   Float_t         jetRawEn[180];   //[nJet]
   Float_t         jetCharge[180];   //[nJet]
   Float_t         jetNeutralEmEnergy[180];   //[nJet]
   Float_t         jetNeutralEmEnergyFraction[180];   //[nJet]
   Float_t         jetNeutralHadronEnergy[180];   //[nJet]
   Float_t         jetNeutralHadronEnergyFraction[180];   //[nJet]
   Int_t           jetNConstituents[180];   //[nJet]
   Float_t         jetChargedEmEnergy[180];   //[nJet]
   Float_t         jetChargedEmEnergyFraction[180];   //[nJet]
   Float_t         jetChargedHadronEnergy[180];   //[nJet]
   Float_t         jetChargedHadronEnergyFraction[180];   //[nJet]
   Int_t           jetChargedHadronMultiplicity[180];   //[nJet]
   Float_t         jetChargedMuEnergy[180];   //[nJet]
   Float_t         jetChargedMuEnergyFraction[180];   //[nJet]
   Double_t        jetJVAlpha[180];   //[nJet]
   Double_t        jetJVBeta[180];   //[nJet]
   Int_t           jetGenJetIndex[180];   //[nJet]
   Float_t         jetGenJetEn[180];   //[nJet]
   Float_t         jetGenJetPt[180];   //[nJet]
   Float_t         jetGenJetEta[180];   //[nJet]
   Float_t         jetGenJetPhi[180];   //[nJet]
   Float_t         jetGenJetMass[180];   //[nJet]
   Int_t           jetGenPartonID[180];   //[nJet]
   Int_t           jetGenPartonMomID[180];   //[nJet]
   Int_t           nZee;
   Float_t         ZeeMass[80];   //[nZee]
   Float_t         ZeePt[80];   //[nZee]
   Float_t         ZeeEta[80];   //[nZee]
   Float_t         ZeePhi[80];   //[nZee]
   Int_t           ZeeLeg1Index[80];   //[nZee]
   Int_t           ZeeLeg2Index[80];   //[nZee]
   Int_t           nZmumu;
   Float_t         ZmumuMass[80];   //[nZmumu]
   Float_t         ZmumuPt[80];   //[nZmumu]
   Float_t         ZmumuEta[80];   //[nZmumu]
   Float_t         ZmumuPhi[80];   //[nZmumu]
   Int_t           ZmumuLeg1Index[80];   //[nZmumu]
   Int_t           ZmumuLeg2Index[80];   //[nZmumu]
   Int_t           nWenu;
   Float_t         WenuMassTCaloMET[80];   //[nWenu]
   Float_t         WenuEtCaloMET[80];   //[nWenu]
   Float_t         WenuACopCaloMET[80];   //[nWenu]
   Float_t         WenuMassTTcMET[80];   //[nWenu]
   Float_t         WenuEtTcMET[80];   //[nWenu]
   Float_t         WenuACopTcMET[80];   //[nWenu]
   Float_t         WenuMassTPfMET[80];   //[nWenu]
   Float_t         WenuEtPfMET[80];   //[nWenu]
   Float_t         WenuACopPfMET[80];   //[nWenu]
   Int_t           WenuEleIndex[80];   //[nWenu]
   Int_t           nWmunu;
   Float_t         WmunuMassTCaloMET[80];   //[nWmunu]
   Float_t         WmunuEtCaloMET[80];   //[nWmunu]
   Float_t         WmunuACopCaloMET[80];   //[nWmunu]
   Float_t         WmunuMassTTcMET[80];   //[nWmunu]
   Float_t         WmunuEtTcMET[80];   //[nWmunu]
   Float_t         WmunuACopTcMET[80];   //[nWmunu]
   Float_t         WmunuMassTPfMET[80];   //[nWmunu]
   Float_t         WmunuEtPfMET[80];   //[nWmunu]
   Float_t         WmunuACopPfMET[80];   //[nWmunu]
   Int_t           WmunuMuIndex[80];   //[nWmunu]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_ttbit0;   //!
   TBranch        *b_nHLT;   //!
   TBranch        *b_HLT;   //!
   TBranch        *b_HLTIndex;   //!
   TBranch        *b_HLTprescale;   //!
   TBranch        *b_nHFTowersP;   //!
   TBranch        *b_nHFTowersN;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vtxNTrk;   //!
   TBranch        *b_vtxNDF;   //!
   TBranch        *b_vtxD0;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_IsVtxGood;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_nGoodTrk;   //!
   TBranch        *b_IsTracksGood;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_sigma;   //!
   TBranch        *b_rhoNeutral;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_nBX;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_BXPU;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcDecayType;   //!
   TBranch        *b_mcIsoDR03;   //!
   TBranch        *b_mcCalIsoDR03;   //!
   TBranch        *b_mcTrkIsoDR03;   //!
   TBranch        *b_mcIsoDR04;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcTrkIsoDR04;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETx;   //!
   TBranch        *b_genMETy;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_METx;   //!
   TBranch        *b_METy;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_METsumEt;   //!
   TBranch        *b_uncorrMET;   //!
   TBranch        *b_uncorrMETPhi;   //!
   TBranch        *b_uncorrMETSumEt;   //!
   TBranch        *b_tcMET;   //!
   TBranch        *b_tcMETx;   //!
   TBranch        *b_tcMETy;   //!
   TBranch        *b_tcMETPhi;   //!
   TBranch        *b_tcMETsumEt;   //!
   TBranch        *b_tcMETmEtSig;   //!
   TBranch        *b_tcMETSig;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETx;   //!
   TBranch        *b_pfMETy;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMETsumEt;   //!
   TBranch        *b_pfMETmEtSig;   //!
   TBranch        *b_pfMETSig;   //!
   TBranch        *b_TypeIpfMET;   //!
   TBranch        *b_TypeIpfMETx;   //!
   TBranch        *b_TypeIpfMETy;   //!
   TBranch        *b_TypeIpfMETPhi;   //!
   TBranch        *b_TypeIpfMETsumEt;   //!
   TBranch        *b_TypeIpfMETmEtSig;   //!
   TBranch        *b_TypeIpfMETSig;   //!
   TBranch        *b_TypeIpIIpfMET;   //!
   TBranch        *b_TypeIpIIpfMETx;   //!
   TBranch        *b_TypeIpIIpfMETy;   //!
   TBranch        *b_TypeIpIIpfMETPhi;   //!
   TBranch        *b_TypeIpIIpfMETsumEt;   //!
   TBranch        *b_TypeIpIIpfMETmEtSig;   //!
   TBranch        *b_TypeIpIIpfMETSig;   //!
   TBranch        *b_SmearedpfMET;   //!
   TBranch        *b_SmearedpfMETx;   //!
   TBranch        *b_SmearedpfMETy;   //!
   TBranch        *b_SmearedpfMETPhi;   //!
   TBranch        *b_SmearedpfMETsumEt;   //!
   TBranch        *b_SmearedpfMETmEtSig;   //!
   TBranch        *b_SmearedpfMETSig;   //!
   TBranch        *b_SmearedTypeIpfMET;   //!
   TBranch        *b_SmearedTypeIpfMETx;   //!
   TBranch        *b_SmearedTypeIpfMETy;   //!
   TBranch        *b_SmearedTypeIpfMETPhi;   //!
   TBranch        *b_SmearedTypeIpfMETsumEt;   //!
   TBranch        *b_SmearedTypeIpfMETmEtSig;   //!
   TBranch        *b_SmearedTypeIpfMETSig;   //!
   TBranch        *b_npfCharged;   //!
   TBranch        *b_pfChargedSumPt;   //!
   TBranch        *b_npfChargedHadron;   //!
   TBranch        *b_pfChargedHadronSumPt;   //!
   TBranch        *b_npfLepton;   //!
   TBranch        *b_pfLeptonSumPt;   //!
   TBranch        *b_npfNeutral;   //!
   TBranch        *b_pfNeutralSumPt;   //!
   TBranch        *b_npfNeutralHadron;   //!
   TBranch        *b_pfNeutralHadronSumPt;   //!
   TBranch        *b_npfPhoton;   //!
   TBranch        *b_pfPhotonSumPt;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleTrg;   //!
   TBranch        *b_eleID;   //!
   TBranch        *b_eleIDLH;   //!
   TBranch        *b_eleClass;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_elePz;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleVtx;   //!
   TBranch        *b_eleCaloPos;   //!
   TBranch        *b_eleSCPos;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_elePin;   //!
   TBranch        *b_elePout;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_elenBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaEtaEta;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIPhi;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleE3x3;   //!
   TBranch        *b_eleSeedTime;   //!
   TBranch        *b_eleSeedEnergy;   //!
   TBranch        *b_eleRecoFlag;   //!
   TBranch        *b_eleSeverity;   //!
   TBranch        *b_eleGenIndex;   //!
   TBranch        *b_eleGenGMomPID;   //!
   TBranch        *b_eleGenMomPID;   //!
   TBranch        *b_eleGenMomPt;   //!
   TBranch        *b_eleIsoTrkDR03;   //!
   TBranch        *b_eleIsoEcalDR03;   //!
   TBranch        *b_eleIsoHcalDR03;   //!
   TBranch        *b_eleIsoHcalSolidDR03;   //!
   TBranch        *b_eleIsoTrkDR04;   //!
   TBranch        *b_eleIsoEcalDR04;   //!
   TBranch        *b_eleIsoHcalDR04;   //!
   TBranch        *b_eleIsoHcalSolidDR04;   //!
   TBranch        *b_eleConvDist;   //!
   TBranch        *b_eleConvDcot;   //!
   TBranch        *b_eleConvRadius;   //!
   TBranch        *b_eleConvFlag;   //!
   TBranch        *b_eleConvMissinghit;   //!
   TBranch        *b_eleESRatio;   //!
   TBranch        *b_eleESProfileFront;   //!
   TBranch        *b_eleESProfileRear;   //!
   TBranch        *b_elePV2D;   //!
   TBranch        *b_elePV3D;   //!
   TBranch        *b_eleBS2D;   //!
   TBranch        *b_eleBS3D;   //!
   TBranch        *b_elePVD0;   //!
   TBranch        *b_elePVDz;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoTrg;   //!
   TBranch        *b_phoIsPhoton;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoPz;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoTrkIsoSolidDR03;   //!
   TBranch        *b_phoTrkIsoHollowDR03;   //!
   TBranch        *b_phoNTrkSolidDR03;   //!
   TBranch        *b_phoNTrkHollowDR03;   //!
   TBranch        *b_phoEcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoDR03;   //!
   TBranch        *b_phoHcalIsoSolidDR03;   //!
   TBranch        *b_phoTrkIsoSolidDR04;   //!
   TBranch        *b_phoTrkIsoHollowDR04;   //!
   TBranch        *b_phoNTrkSolidDR04;   //!
   TBranch        *b_phoNTrkHollowDR04;   //!
   TBranch        *b_phoEcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoDR04;   //!
   TBranch        *b_phoHcalIsoSolidDR04;   //!
   TBranch        *b_phoEtVtx;   //!
   TBranch        *b_phoEtaVtx;   //!
   TBranch        *b_phoPhiVtx;   //!
   TBranch        *b_phoTrkIsoSolidDR03Vtx;   //!
   TBranch        *b_phoTrkIsoHollowDR03Vtx;   //!
   TBranch        *b_phoTrkIsoSolidDR04Vtx;   //!
   TBranch        *b_phoTrkIsoHollowDR04Vtx;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoSigmaEtaEta;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedEnergy;   //!
   TBranch        *b_phoRecoFlag;   //!
   TBranch        *b_phoSeverity;   //!
   TBranch        *b_phoPos;   //!
   TBranch        *b_phoGenIndex;   //!
   TBranch        *b_phoGenGMomPID;   //!
   TBranch        *b_phoGenMomPID;   //!
   TBranch        *b_phoGenMomPt;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoESE;   //!
   TBranch        *b_phoSCEt;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoVtx;   //!
   TBranch        *b_phoVtxD0;   //!
   TBranch        *b_phoOverlap;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoIsConv;   //!
   TBranch        *b_phoESRatio;   //!
   TBranch        *b_phoESProfileFront;   //!
   TBranch        *b_phoESProfileRear;   //!
   TBranch        *b_phoNTracks;   //!
   TBranch        *b_phoConvPairInvariantMass;   //!
   TBranch        *b_phoConvPairCotThetaSeparation;   //!
   TBranch        *b_phoConvPairMomentumEta;   //!
   TBranch        *b_phoConvPairMomentumPhi;   //!
   TBranch        *b_phoConvPairMomentumX;   //!
   TBranch        *b_phoConvPairMomentumY;   //!
   TBranch        *b_phoConvPairMomentumZ;   //!
   TBranch        *b_phoConvDistOfMinimumApproach;   //!
   TBranch        *b_phoConvDPhiTracksAtVtx;   //!
   TBranch        *b_phoConvDPhiTracksAtEcal;   //!
   TBranch        *b_phoConvDEtaTracksAtEcal;   //!
   TBranch        *b_phoConvVtxValid;   //!
   TBranch        *b_phoConvVtxEta;   //!
   TBranch        *b_phoConvVtxPhi;   //!
   TBranch        *b_phoConvVtxR;   //!
   TBranch        *b_phoConvVtxX;   //!
   TBranch        *b_phoConvVtxY;   //!
   TBranch        *b_phoConvVtxZ;   //!
   TBranch        *b_phoConvVtxChi2;   //!
   TBranch        *b_phoConvVtxNdof;   //!
   TBranch        *b_phoConvChi2Prob;   //!
   TBranch        *b_phoConvEoverP;   //!
   TBranch        *b_phoNxtal;   //!
   TBranch        *b_phoXtalTime;   //!
   TBranch        *b_phoXtalEnergy;   //!
   TBranch        *b_phoXtalZ;   //!
   TBranch        *b_phoXtalX;   //!
   TBranch        *b_phoXtalY;   //!
   TBranch        *b_phoXtalEta;   //!
   TBranch        *b_phoXtalPhi;   //!
   TBranch        *b_pho5x5Time;   //!
   TBranch        *b_pho5x5Energy;   //!
   TBranch        *b_pho5x5Z;   //!
   TBranch        *b_pho5x5X;   //!
   TBranch        *b_pho5x5Y;   //!
   TBranch        *b_pho5x5Eta;   //!
   TBranch        *b_pho5x5Phi;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muTrg;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muPz;   //!
   TBranch        *b_muGenIndex;   //!
   TBranch        *b_muGenGMomPID;   //!
   TBranch        *b_muGenMomPID;   //!
   TBranch        *b_muGenMomPt;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muIsoCalo;   //!
   TBranch        *b_muIsoEcal;   //!
   TBranch        *b_muIsoHcal;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muEmVeto;   //!
   TBranch        *b_muHadVeto;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muID;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muPVD0;   //!
   TBranch        *b_muPVDz;   //!
   TBranch        *b_muValidFraction;   //!
   TBranch        *b_muTrkdPt;   //!
   TBranch        *b_muNumberOfHits;   //!
   TBranch        *b_muNumberOfValidHits;   //!
   TBranch        *b_muNumberOfInactiveHits;   //!
   TBranch        *b_muNumberOfValidTrkHits;   //!
   TBranch        *b_muNumberOfValidPixelHits;   //!
   TBranch        *b_muNumberOfValidMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muChambers;   //!
   TBranch        *b_muPV2D;   //!
   TBranch        *b_muPV3D;   //!
   TBranch        *b_muBS2D;   //!
   TBranch        *b_muBS3D;   //!
   TBranch        *b_muVtx;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetTrg;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetEt;   //!
   TBranch        *b_jetpartonFlavour;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetCharge;   //!
   TBranch        *b_jetNeutralEmEnergy;   //!
   TBranch        *b_jetNeutralEmEnergyFraction;   //!
   TBranch        *b_jetNeutralHadronEnergy;   //!
   TBranch        *b_jetNeutralHadronEnergyFraction;   //!
   TBranch        *b_jetNConstituents;   //!
   TBranch        *b_jetChargedEmEnergy;   //!
   TBranch        *b_jetChargedEmEnergyFraction;   //!
   TBranch        *b_jetChargedHadronEnergy;   //!
   TBranch        *b_jetChargedHadronEnergyFraction;   //!
   TBranch        *b_jetChargedHadronMultiplicity;   //!
   TBranch        *b_jetChargedMuEnergy;   //!
   TBranch        *b_jetChargedMuEnergyFraction;   //!
   TBranch        *b_jetJVAlpha;   //!
   TBranch        *b_jetJVBeta;   //!
   TBranch        *b_jetGenJetIndex;   //!
   TBranch        *b_jetGenJetEn;   //!
   TBranch        *b_jetGenJetPt;   //!
   TBranch        *b_jetGenJetEta;   //!
   TBranch        *b_jetGenJetPhi;   //!
   TBranch        *b_jetGenJetMass;   //!
   TBranch        *b_jetGenPartonID;   //!
   TBranch        *b_jetGenPartonMomID;   //!
   TBranch        *b_nZee;   //!
   TBranch        *b_ZeeMass;   //!
   TBranch        *b_ZeePt;   //!
   TBranch        *b_ZeeEta;   //!
   TBranch        *b_ZeePhi;   //!
   TBranch        *b_ZeeLeg1Index;   //!
   TBranch        *b_ZeeLeg2Index;   //!
   TBranch        *b_nZmumu;   //!
   TBranch        *b_ZmumuMass;   //!
   TBranch        *b_ZmumuPt;   //!
   TBranch        *b_ZmumuEta;   //!
   TBranch        *b_ZmumuPhi;   //!
   TBranch        *b_ZmumuLeg1Index;   //!
   TBranch        *b_ZmumuLeg2Index;   //!
   TBranch        *b_nWenu;   //!
   TBranch        *b_WenuMassTCaloMET;   //!
   TBranch        *b_WenuEtCaloMET;   //!
   TBranch        *b_WenuACopCaloMET;   //!
   TBranch        *b_WenuMassTTcMET;   //!
   TBranch        *b_WenuEtTcMET;   //!
   TBranch        *b_WenuACopTcMET;   //!
   TBranch        *b_WenuMassTPfMET;   //!
   TBranch        *b_WenuEtPfMET;   //!
   TBranch        *b_WenuACopPfMET;   //!
   TBranch        *b_WenuEleIndex;   //!
   TBranch        *b_nWmunu;   //!
   TBranch        *b_WmunuMassTCaloMET;   //!
   TBranch        *b_WmunuEtCaloMET;   //!
   TBranch        *b_WmunuACopCaloMET;   //!
   TBranch        *b_WmunuMassTTcMET;   //!
   TBranch        *b_WmunuEtTcMET;   //!
   TBranch        *b_WmunuACopTcMET;   //!
   TBranch        *b_WmunuMassTPfMET;   //!
   TBranch        *b_WmunuEtPfMET;   //!
   TBranch        *b_WmunuACopPfMET;   //!
   TBranch        *b_WmunuMuIndex;   //!

   anaVgNtuple();
   anaVgNtuple(TString pTag_ , vector<string> op_, TString FilePath_ = "");
   virtual ~anaVgNtuple();
   virtual Int_t           Cut(Long64_t entry);
   virtual Int_t           GetEntry(Long64_t entry);
   virtual Long64_t        LoadTree(Long64_t entry);
   virtual void            Init(TTree *tree, Bool_t runData = false);
   virtual void            Loop(Int_t SampleIndex = 0);

   virtual Bool_t          checkOption();
   virtual Bool_t          CleanEvent();
   virtual vector<int>     mcElectron();
   virtual Bool_t          mcEleMatcher(Int_t pIndex, TString pName);
   virtual vector<int>     mcMuon();
   virtual Bool_t          mcMuMatcher(Int_t pIndex, TString pName);
   virtual vector<int>     mcPhoton();
   virtual Bool_t          mcPhoMatcher(Int_t iPho);
   virtual void            addScale(TString objName);
   virtual void            addResolution(TString objName);
   virtual Bool_t          IsSpike(Int_t index , TString option);
   virtual Int_t           getNPassJets();
//   virtual vector<int>     mcElectron();
//   virtual Bool_t          mcEleMatcher(Int_t iEle);
   virtual Bool_t          passHLT(Int_t &trgIndex);
   virtual Bool_t          passHLTMu(Int_t &trgIndex);
   virtual Bool_t          checkHLT(Int_t iLep, Int_t lepTrg_index);
   virtual Bool_t          muClean(Int_t iMu);
   virtual Bool_t          muPassID(Int_t iMu);
   virtual Bool_t          IsISR(Int_t iPho);
   virtual Bool_t          IsFSR(Int_t iPho);
   virtual Bool_t          DiMuPair(Int_t & iMu1, Int_t& iMu2 , Int_t trg_index);
   virtual vector<int>     ElectronPassLevel(Int_t index_eID, Int_t index_trg);
   virtual Bool_t          DiElePair(Int_t &ele1_index , Int_t &ele2_index, vector<int> eID_Level) ;
   virtual Bool_t          elePassID(Int_t index_e , Int_t index_eID, unsigned int CutSet = 15);
   virtual unsigned int    eleTrgWP(Int_t iEle);
   virtual Int_t           phoCleanLevel(Int_t lep1_index, Int_t lep2_index ,Int_t iPho);
   virtual Int_t           phoPassIDLevel(Int_t index_g , Int_t index_gID);
   virtual Bool_t          phoPassID(Int_t index_g , Int_t index_gID);
   virtual void            printMC(Int_t pdgID = 0);
   virtual void            printEvent(int index_e1, int index_e2, int index_pho) ;
   virtual void            printElectron(int iEle) ;
   virtual Bool_t          Notify();
   virtual void            Show(Long64_t entry = -1);
   virtual void            ShowEvent(Long64_t entry, Int_t iEle1, Int_t iEle2, Int_t iPho);
   virtual Bool_t          IsFakeable(int iPho);
   virtual Int_t           PhoSelection(Int_t lep1, Int_t lep2);
   virtual Int_t           PhoSelectLevel(Int_t lep1, Int_t lep2);
   virtual vector<Float_t> newEleEn();
   virtual vector<Float_t> newElePt(vector<Float_t> EnCorr);
   virtual vector<Float_t> newPhoE();
   virtual vector<Float_t> newPhoEt(vector<Float_t> EnCorr);
   virtual void            EcorrectionYM(TString flavor); // energy correction from Yurii M.
   virtual Float_t         getCorrection_YM(Float_t en, Float_t absEta, Float_t r9, 
				   							Int_t run, Bool_t data = true); // energy correction from Yurii M.
   virtual Float_t         DeltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2);
   virtual vector<Double_t>  generate_S3_weights(TH1D* data_npu_estimated);
   virtual vector<Double_t>  generate_S4_weights(TH1D* data_npu_estimated);
   virtual vector<Double_t>  generate_S6_weights(TH1D* data_npu_estimated);
   virtual vector<Double_t>  generate_flat10_weights(TH1D* data_npu_estimated);
   virtual Double_t          PU3D_Weight(Int_t pm1 , Int_t p0 , Int_t p1 , TH3* hin);
   virtual void              InitHisto();
   virtual Bool_t            FillHistoPre();
   virtual void              FillHistoZ(Int_t lep1_index , Int_t lep2_index);
   virtual void              FillHisto(Int_t lep1_index , Int_t lep2_index, Int_t pho_index);
   virtual void              WrHisto();
   virtual void              RmHisto();
   virtual void              InitTree();
   virtual void              FillTree(Int_t lep1_index , Int_t lep2_index, Int_t pho_index);
   virtual void              FillTreeEleTrg();
   virtual void              WrTree();
   //testlineline2efjaiore
};

#endif

#ifdef anaVgNtuple_cxx
anaVgNtuple::anaVgNtuple(){
	TChain * Chain = new TChain("VgAnalyzerKit/EventTree","");
	Chain->Add("/data2/kschen/Vgamma/Zgamma2011/Skim/VgKitV11/VgKitV11/DoubleElectron_May10ReReco_JSON_v3.root");
    Init(Chain,true);
}
anaVgNtuple::anaVgNtuple(TString pTag_ , vector<string> op_, TString FilePath_)
{
	OO.set_option(op_);
	ProcessTag = pTag_ ;
	InputPath = OO.call_vector_string(string("InputPath"));
	
	if (FilePath_ != ""){
		TString myOption = "FilePath" ;
		myOption = myOption + " " + FilePath_ ;
		OO.add_option(string(myOption.Data()));
	}
	TChain * Chain = new TChain("VgAnalyzerKit/EventTree","");
	if (InputPath.size() != 0)
	{
		for (vector<string>::iterator it = InputPath.begin() ; it != InputPath.end(); it++)
		{
			Chain->Add(TString(*it));
		}
	}
	else if ( ProcessTag == "Data" && OO.call_bool(string("doEle")) ) 
	{
		//Chain->Add("/data3/ncuhep/data/428p5_vgamma_42x_v14/DoubleElectron_Run2011A-May10ReReco-v1/DoubleElectron_Run2011A-May10ReReco-v1.root");
		//Chain->Add("/data4/kschen/Skim/VgKitV14/test/DoubleElectron_Run2011A-May10ReReco-v1.root");

		Chain->Add("/data4/kschen/Skim/VgKitV14/test/DoubleElectron_Run2011A-May10ReReco-v1.root");
		Chain->Add("/data4/kschen/Skim/VgKitV14/test/DoubleElectron_Run2011A-PromptReco-v4.root");
		Chain->Add("/data4/kschen/Skim/VgKitV14/test/DoubleElectron_Run2011A-05Aug2011-v1_excTrg_addCSC.root");
		Chain->Add("/data4/kschen/Skim/VgKitV14/test/DoubleElectron_Run2011A-03Oct2011-v1.root");
		Chain->Add("/data4/kschen/Skim/VgKitV14/test/DoubleElectron_Run2011B-PromptReco-v1_run175860_178078.root");
		Chain->Add("/data4/kschen/Skim/VgKitV14/test/DoubleElectron_Run2011B-PromptReco-v1_run178098_178677.root");
		Chain->Add("/data4/kschen/Skim/VgKitV14/test/DoubleElectron_Run2011B-PromptReco-v1_run178703_179431.root");
		Chain->Add("/data4/kschen/Skim/VgKitV14/test/DoubleElectron_Run2011B-PromptReco-v1_run179434_180252.root");
		//Chain->Add("/data4/kschen/Skim/VgKitV14/test/");
		//Chain->Add("/data4/kschen/Skim/VgKitV14/test/");
		//Chain->Add("/data4/kschen/Skim/VgKitV14/test/");
		//Chain->Add("/data4/kschen/Skim/VgKitV14/test/");
		//Chain->Add("/data4/kschen/Skim/VgKitV14/test/");


		//Chain->Add("/data3/ncuhep/data/428p5_vgamma_42x_v14/DoubleElectron_Run2011B-PromptReco-v1_run178098_178677/DoubleElectron_Run2011B-PromptReco-v1_run178098_178677.root");

		//Chain->Add("/data3/ncuhep/data/425_vgamma_42x_v11/DoubleElectron_Run2011A_May10ReRecov1_JSONv3/DoubleElectron_May10ReReco_JSON_v3.root");
//		Chain->Add("/scratch/kschen/DataProduction/VgKitV9/CMSSW_4_2_5/src/CRABjob/2011A_PromptReco172798_173243/Sample/DoubleElectron_Run2011A-PromptReco-v6_run172798_173243/DoubleElectron_Run2011A-PromptReco-v6_run172798_173243.root");

		//Chain->Add("/data3/ncuhep/data/428p5_vgamma_42x_v14/DoubleElectron_Run2011B-PromptReco-v1_run178098_178677/DoubleElectron_Run2011B-PromptReco-v1_run178098_178677.root");

	}
	else
	if ( ProcessTag == "Data" && OO.call_bool(string("doMu")) ) 
	{
		Chain->Add("/data2/kschen/Vgamma/Zgamma2011/EventSelection/MC/Skim/MuonSkim/DoubleMu_Run2011A_May10ReRecov1.root") ;
		Chain->Add("/data2/kschen/Vgamma/Zgamma2011/EventSelection/MC/Skim/MuonSkim/DoubleMu_Run2011A_PromptRecov4_Run165088to167913.root") ;
	}
	else {
		Chain->Add(FilePath_);
	}

   Init(Chain,ProcessTag.Contains("Data"));
}

anaVgNtuple::~anaVgNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t anaVgNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t anaVgNtuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void anaVgNtuple::Init(TTree *tree, Bool_t runData)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchStatus("*",0);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
//   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
//   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
//   fChain->SetBranchAddress("ttbit0", &ttbit0, &b_ttbit0);
   fChain->SetBranchAddress("nHLT", &nHLT, &b_nHLT);
   fChain->SetBranchAddress("HLT", HLT, &b_HLT);
   fChain->SetBranchAddress("HLTIndex", HLTIndex, &b_HLTIndex);
   fChain->SetBranchAddress("HLTprescale", HLTprescale, &b_HLTprescale);
//   fChain->SetBranchAddress("nHFTowersP", &nHFTowersP, &b_nHFTowersP);
//   fChain->SetBranchAddress("nHFTowersN", &nHFTowersN, &b_nHFTowersN);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("vtxNTrk", vtxNTrk, &b_vtxNTrk);
   fChain->SetBranchAddress("vtxNDF", vtxNDF, &b_vtxNDF);
   fChain->SetBranchAddress("vtxD0", vtxD0, &b_vtxD0);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("IsVtxGood", &IsVtxGood, &b_IsVtxGood);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("nGoodTrk", &nGoodTrk, &b_nGoodTrk);
   fChain->SetBranchAddress("IsTracksGood", &IsTracksGood, &b_IsTracksGood);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("sigma", &sigma, &b_sigma);
   if (!runData) fChain->SetBranchAddress("pdf", pdf, &b_pdf);
   if (!runData) fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   if (!runData) fChain->SetBranchAddress("processID", &processID, &b_processID);
   if (!runData) fChain->SetBranchAddress("nBX", &nBX, &b_nBX);
   if (!runData) fChain->SetBranchAddress("nPU", nPU, &b_nPU);
   if (!runData) fChain->SetBranchAddress("BXPU", BXPU, &b_BXPU);
   if (!runData) fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   if (!runData) fChain->SetBranchAddress("mcPID", mcPID, &b_mcPID);
   if (!runData) fChain->SetBranchAddress("mcPt", mcPt, &b_mcPt);
   if (!runData) fChain->SetBranchAddress("mcMass", mcMass, &b_mcMass);
   if (!runData) fChain->SetBranchAddress("mcEta", mcEta, &b_mcEta);
   if (!runData) fChain->SetBranchAddress("mcPhi", mcPhi, &b_mcPhi);
   if (!runData) fChain->SetBranchAddress("mcGMomPID", mcGMomPID, &b_mcGMomPID);
   if (!runData) fChain->SetBranchAddress("mcMomPID", mcMomPID, &b_mcMomPID);
   if (!runData) fChain->SetBranchAddress("mcMomPt", mcMomPt, &b_mcMomPt);
   if (!runData) fChain->SetBranchAddress("mcMomMass", mcMomMass, &b_mcMomMass);
   if (!runData) fChain->SetBranchAddress("mcMomEta", mcMomEta, &b_mcMomEta);
   if (!runData) fChain->SetBranchAddress("mcMomPhi", mcMomPhi, &b_mcMomPhi);
   if (!runData) fChain->SetBranchAddress("mcIndex", mcIndex, &b_mcIndex);
   if (!runData) fChain->SetBranchAddress("mcDecayType", mcDecayType, &b_mcDecayType);
   if (!runData) fChain->SetBranchAddress("mcIsoDR03", mcIsoDR03, &b_mcIsoDR03);
   if (!runData) fChain->SetBranchAddress("mcCalIsoDR03", mcCalIsoDR03, &b_mcCalIsoDR03);
   if (!runData) fChain->SetBranchAddress("mcTrkIsoDR03", mcTrkIsoDR03, &b_mcTrkIsoDR03);
   if (!runData) fChain->SetBranchAddress("mcIsoDR04", mcIsoDR04, &b_mcIsoDR04);
   if (!runData) fChain->SetBranchAddress("mcCalIsoDR04", mcCalIsoDR04, &b_mcCalIsoDR04);
   if (!runData) fChain->SetBranchAddress("mcTrkIsoDR04", mcTrkIsoDR04, &b_mcTrkIsoDR04);
   if (!runData) fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
//   if (!runData) fChain->SetBranchAddress("genMETx", &genMETx, &b_genMETx);
//   if (!runData) fChain->SetBranchAddress("genMETy", &genMETy, &b_genMETy);
   if (!runData) fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
//   fChain->SetBranchAddress("METx", &METx, &b_METx);
//   fChain->SetBranchAddress("METy", &METy, &b_METy);
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
//   fChain->SetBranchAddress("METsumEt", &METsumEt, &b_METsumEt);
//   fChain->SetBranchAddress("uncorrMET", uncorrMET, &b_uncorrMET);
//   fChain->SetBranchAddress("uncorrMETPhi", uncorrMETPhi, &b_uncorrMETPhi);
//   fChain->SetBranchAddress("uncorrMETSumEt", uncorrMETSumEt, &b_uncorrMETSumEt);
   fChain->SetBranchAddress("tcMET", &tcMET, &b_tcMET);
//   fChain->SetBranchAddress("tcMETx", &tcMETx, &b_tcMETx);
//   fChain->SetBranchAddress("tcMETy", &tcMETy, &b_tcMETy);
   fChain->SetBranchAddress("tcMETPhi", &tcMETPhi, &b_tcMETPhi);
//   fChain->SetBranchAddress("tcMETsumEt", &tcMETsumEt, &b_tcMETsumEt);
//   fChain->SetBranchAddress("tcMETmEtSig", &tcMETmEtSig, &b_tcMETmEtSig);
//   fChain->SetBranchAddress("tcMETSig", &tcMETSig, &b_tcMETSig);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
//   fChain->SetBranchAddress("pfMETx", &pfMETx, &b_pfMETx);
//   fChain->SetBranchAddress("pfMETy", &pfMETy, &b_pfMETy);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMETsumEt", &pfMETsumEt, &b_pfMETsumEt);
   fChain->SetBranchAddress("pfMETmEtSig", &pfMETmEtSig, &b_pfMETmEtSig);
   fChain->SetBranchAddress("pfMETSig", &pfMETSig, &b_pfMETSig);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleTrg", eleTrg, &b_eleTrg);
   fChain->SetBranchAddress("eleID", eleID, &b_eleID);
   fChain->SetBranchAddress("eleIDLH", eleIDLH, &b_eleIDLH);
   fChain->SetBranchAddress("eleClass", eleClass, &b_eleClass);
   fChain->SetBranchAddress("eleCharge", eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleEn", eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleSCRawEn", eleSCRawEn, &b_eleSCRawEn);
//   fChain->SetBranchAddress("eleESEn", eleESEn, &b_eleESEn);
   fChain->SetBranchAddress("eleSCEn", eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("elePt", elePt, &b_elePt);
   fChain->SetBranchAddress("elePz", elePz, &b_elePz);
   fChain->SetBranchAddress("eleEta", eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleSCEta", eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCEtaWidth", eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleVtx", eleVtx, &b_eleVtx);
//   fChain->SetBranchAddress("eleCaloPos", eleCaloPos, &b_eleCaloPos);
//   fChain->SetBranchAddress("eleSCPos", eleSCPos, &b_eleSCPos);
   fChain->SetBranchAddress("eleHoverE", eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("elePin", elePin, &b_elePin);
   fChain->SetBranchAddress("elePout", elePout, &b_elePout);
   fChain->SetBranchAddress("eleBrem", eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("elenBrem", elenBrem, &b_elenBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaEtaEta", eleSigmaEtaEta, &b_eleSigmaEtaEta);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIPhi", eleSigmaIEtaIPhi, &b_eleSigmaIEtaIPhi);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
//   fChain->SetBranchAddress("eleE2overE9", eleE2overE9, &b_eleE2overE9);
   fChain->SetBranchAddress("eleE3x3", eleE3x3, &b_eleE3x3);
//   fChain->SetBranchAddress("eleSeedTime", eleSeedTime, &b_eleSeedTime);
//   fChain->SetBranchAddress("eleRecoFlag", eleRecoFlag, &b_eleRecoFlag);
//   fChain->SetBranchAddress("eleSeverity", eleSeverity, &b_eleSeverity);
   if (!runData) fChain->SetBranchAddress("eleGenIndex", eleGenIndex, &b_eleGenIndex);
   if (!runData) fChain->SetBranchAddress("eleGenGMomPID", eleGenGMomPID, &b_eleGenGMomPID);
   if (!runData) fChain->SetBranchAddress("eleGenMomPID", eleGenMomPID, &b_eleGenMomPID);
   if (!runData) fChain->SetBranchAddress("eleGenMomPt", eleGenMomPt, &b_eleGenMomPt);
   fChain->SetBranchAddress("eleIsoTrkDR03", eleIsoTrkDR03, &b_eleIsoTrkDR03);
   fChain->SetBranchAddress("eleIsoEcalDR03", eleIsoEcalDR03, &b_eleIsoEcalDR03);
   fChain->SetBranchAddress("eleIsoHcalDR03", eleIsoHcalDR03, &b_eleIsoHcalDR03);
   fChain->SetBranchAddress("eleIsoHcalSolidDR03", eleIsoHcalSolidDR03, &b_eleIsoHcalSolidDR03);
//   fChain->SetBranchAddress("eleIsoTrkDR04", eleIsoTrkDR04, &b_eleIsoTrkDR04);
//   fChain->SetBranchAddress("eleIsoEcalDR04", eleIsoEcalDR04, &b_eleIsoEcalDR04);
//   fChain->SetBranchAddress("eleIsoHcalDR04", eleIsoHcalDR04, &b_eleIsoHcalDR04);
//   fChain->SetBranchAddress("eleIsoHcalSolidDR04", eleIsoHcalSolidDR04, &b_eleIsoHcalSolidDR04);
   fChain->SetBranchAddress("eleConvDist", eleConvDist, &b_eleConvDist);
   fChain->SetBranchAddress("eleConvDcot", eleConvDcot, &b_eleConvDcot);
   fChain->SetBranchAddress("eleConvRadius", eleConvRadius, &b_eleConvRadius);
   fChain->SetBranchAddress("eleConvFlag", eleConvFlag, &b_eleConvFlag);
   fChain->SetBranchAddress("eleConvMissinghit", eleConvMissinghit, &b_eleConvMissinghit);
   fChain->SetBranchAddress("eleESRatio", eleESRatio, &b_eleESRatio);
//   fChain->SetBranchAddress("eleESProfileFront", eleESProfileFront, &b_eleESProfileFront);
//   fChain->SetBranchAddress("eleESProfileRear", eleESProfileRear, &b_eleESProfileRear);
//   fChain->SetBranchAddress("elePV2D", elePV2D, &b_elePV2D);
//   fChain->SetBranchAddress("elePV3D", elePV3D, &b_elePV3D);
//   fChain->SetBranchAddress("eleBS2D", eleBS2D, &b_eleBS2D);
//   fChain->SetBranchAddress("eleBS3D", eleBS3D, &b_eleBS3D);
   fChain->SetBranchAddress("elePVD0", elePVD0, &b_elePVD0);
   fChain->SetBranchAddress("elePVDz", elePVDz, &b_elePVDz);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoTrg", phoTrg, &b_phoTrg);
   fChain->SetBranchAddress("phoIsPhoton", phoIsPhoton, &b_phoIsPhoton);
   fChain->SetBranchAddress("phoE", phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoPz", phoPz, &b_phoPz);
   fChain->SetBranchAddress("phoEta", phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoR9", phoR9, &b_phoR9);
//   fChain->SetBranchAddress("phoTrkIsoSolidDR03", phoTrkIsoSolidDR03, &b_phoTrkIsoSolidDR03);
//   fChain->SetBranchAddress("phoTrkIsoHollowDR03", phoTrkIsoHollowDR03, &b_phoTrkIsoHollowDR03);
//   fChain->SetBranchAddress("phoNTrkSolidDR03", phoNTrkSolidDR03, &b_phoNTrkSolidDR03);
//   fChain->SetBranchAddress("phoNTrkHollowDR03", phoNTrkHollowDR03, &b_phoNTrkHollowDR03);
//   fChain->SetBranchAddress("phoEcalIsoDR03", phoEcalIsoDR03, &b_phoEcalIsoDR03);
//   fChain->SetBranchAddress("phoHcalIsoDR03", phoHcalIsoDR03, &b_phoHcalIsoDR03);
//   fChain->SetBranchAddress("phoHcalIsoSolidDR03", phoHcalIsoSolidDR03, &b_phoHcalIsoSolidDR03);
   fChain->SetBranchAddress("phoTrkIsoSolidDR04", phoTrkIsoSolidDR04, &b_phoTrkIsoSolidDR04);
   fChain->SetBranchAddress("phoTrkIsoHollowDR04", phoTrkIsoHollowDR04, &b_phoTrkIsoHollowDR04);
   fChain->SetBranchAddress("phoNTrkSolidDR04", phoNTrkSolidDR04, &b_phoNTrkSolidDR04);
   fChain->SetBranchAddress("phoNTrkHollowDR04", phoNTrkHollowDR04, &b_phoNTrkHollowDR04);
   fChain->SetBranchAddress("phoEcalIsoDR04", phoEcalIsoDR04, &b_phoEcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoDR04", phoHcalIsoDR04, &b_phoHcalIsoDR04);
   fChain->SetBranchAddress("phoHcalIsoSolidDR04", phoHcalIsoSolidDR04, &b_phoHcalIsoSolidDR04);
   fChain->SetBranchAddress("phoHoverE", phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoSigmaEtaEta", phoSigmaEtaEta, &b_phoSigmaEtaEta);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
//   fChain->SetBranchAddress("phoE2overE9", phoE2overE9, &b_phoE2overE9);
   fChain->SetBranchAddress("phoE3x3", phoE3x3, &b_phoE3x3);
   fChain->SetBranchAddress("phoSeedTime", phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoRecoFlag", phoRecoFlag, &b_phoRecoFlag);
   fChain->SetBranchAddress("phoSeverity", phoSeverity, &b_phoSeverity);
//   fChain->SetBranchAddress("phoPos", phoPos, &b_phoPos);
//   fChain->SetBranchAddress("phoRoundness", phoRoundness, &b_phoRoundness);
//   fChain->SetBranchAddress("phoAngle", phoAngle, &b_phoAngle);
   if (!runData) fChain->SetBranchAddress("phoGenIndex", phoGenIndex, &b_phoGenIndex);
   if (!runData) fChain->SetBranchAddress("phoGenGMomPID", phoGenGMomPID, &b_phoGenGMomPID);
   if (!runData) fChain->SetBranchAddress("phoGenMomPID", phoGenMomPID, &b_phoGenMomPID);
   if (!runData) fChain->SetBranchAddress("phoGenMomPt", phoGenMomPt, &b_phoGenMomPt);
   fChain->SetBranchAddress("phoSCE", phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCEt", phoSCEt, &b_phoSCEt);
   fChain->SetBranchAddress("phoSCEta", phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoOverlap", phoOverlap, &b_phoOverlap);
   fChain->SetBranchAddress("phohasPixelSeed", phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoIsConv", phoIsConv, &b_phoIsConv);
//   fChain->SetBranchAddress("phoPi0Disc", phoPi0Disc, &b_phoPi0Disc);
//   fChain->SetBranchAddress("phoESRatio", phoESRatio, &b_phoESRatio);
//   fChain->SetBranchAddress("phoESProfileFront", phoESProfileFront, &b_phoESProfileFront);
//   fChain->SetBranchAddress("phoESProfileRear", phoESProfileRear, &b_phoESProfileRear);
   fChain->SetBranchAddress("phoNTracks", phoNTracks, &b_phoNTracks);
   fChain->SetBranchAddress("phoConvPairInvariantMass", phoConvPairInvariantMass, &b_phoConvPairInvariantMass);
   fChain->SetBranchAddress("phoConvPairCotThetaSeparation", phoConvPairCotThetaSeparation, &b_phoConvPairCotThetaSeparation);
   fChain->SetBranchAddress("phoConvPairMomentumEta", phoConvPairMomentumEta, &b_phoConvPairMomentumEta);
   fChain->SetBranchAddress("phoConvPairMomentumPhi", phoConvPairMomentumPhi, &b_phoConvPairMomentumPhi);
   fChain->SetBranchAddress("phoConvPairMomentumX", phoConvPairMomentumX, &b_phoConvPairMomentumX);
   fChain->SetBranchAddress("phoConvPairMomentumY", phoConvPairMomentumY, &b_phoConvPairMomentumY);
   fChain->SetBranchAddress("phoConvPairMomentumZ", phoConvPairMomentumZ, &b_phoConvPairMomentumZ);
   fChain->SetBranchAddress("phoConvDistOfMinimumApproach", phoConvDistOfMinimumApproach, &b_phoConvDistOfMinimumApproach);
   fChain->SetBranchAddress("phoConvDPhiTracksAtVtx", phoConvDPhiTracksAtVtx, &b_phoConvDPhiTracksAtVtx);
   fChain->SetBranchAddress("phoConvDPhiTracksAtEcal", phoConvDPhiTracksAtEcal, &b_phoConvDPhiTracksAtEcal);
   fChain->SetBranchAddress("phoConvDEtaTracksAtEcal", phoConvDEtaTracksAtEcal, &b_phoConvDEtaTracksAtEcal);
   fChain->SetBranchAddress("phoConvVtxValid", phoConvVtxValid, &b_phoConvVtxValid);
   fChain->SetBranchAddress("phoConvVtxEta", phoConvVtxEta, &b_phoConvVtxEta);
   fChain->SetBranchAddress("phoConvVtxPhi", phoConvVtxPhi, &b_phoConvVtxPhi);
   fChain->SetBranchAddress("phoConvVtxR", phoConvVtxR, &b_phoConvVtxR);
   fChain->SetBranchAddress("phoConvVtxX", phoConvVtxX, &b_phoConvVtxX);
   fChain->SetBranchAddress("phoConvVtxY", phoConvVtxY, &b_phoConvVtxY);
   fChain->SetBranchAddress("phoConvVtxZ", phoConvVtxZ, &b_phoConvVtxZ);
   fChain->SetBranchAddress("phoConvVtxChi2", phoConvVtxChi2, &b_phoConvVtxChi2);
   fChain->SetBranchAddress("phoConvVtxNdof", phoConvVtxNdof, &b_phoConvVtxNdof);
   fChain->SetBranchAddress("phoConvChi2Prob", phoConvChi2Prob, &b_phoConvChi2Prob);
   fChain->SetBranchAddress("phoConvEoverP", phoConvEoverP, &b_phoConvEoverP);
   fChain->SetBranchAddress("phoNxtal", phoNxtal, &b_phoNxtal);
   fChain->SetBranchAddress("phoXtalTime", phoXtalTime, &b_phoXtalTime);
   fChain->SetBranchAddress("phoXtalEnergy", phoXtalEnergy, &b_phoXtalEnergy);
   fChain->SetBranchAddress("phoXtalZ", phoXtalZ, &b_phoXtalZ);
   fChain->SetBranchAddress("phoXtalX", phoXtalX, &b_phoXtalX);
   fChain->SetBranchAddress("phoXtalY", phoXtalY, &b_phoXtalY);
   fChain->SetBranchAddress("phoXtalEta", phoXtalEta, &b_phoXtalEta);
   fChain->SetBranchAddress("phoXtalPhi", phoXtalPhi, &b_phoXtalPhi);
   fChain->SetBranchAddress("pho5x5Time", pho5x5Time, &b_pho5x5Time);
   fChain->SetBranchAddress("pho5x5Energy", pho5x5Energy, &b_pho5x5Energy);
   fChain->SetBranchAddress("pho5x5Z", pho5x5Z, &b_pho5x5Z);
   fChain->SetBranchAddress("pho5x5X", pho5x5X, &b_pho5x5X);
   fChain->SetBranchAddress("pho5x5Y", pho5x5Y, &b_pho5x5Y);
   fChain->SetBranchAddress("pho5x5Eta", pho5x5Eta, &b_pho5x5Eta);
   fChain->SetBranchAddress("pho5x5Phi", pho5x5Phi, &b_pho5x5Phi);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muTrg", muTrg, &b_muTrg);
   fChain->SetBranchAddress("muEta", muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", muCharge, &b_muCharge);
   fChain->SetBranchAddress("muPt", muPt, &b_muPt);
   fChain->SetBranchAddress("muPz", muPz, &b_muPz);
   if (!runData) fChain->SetBranchAddress("muGenIndex", muGenIndex, &b_muGenIndex);
   if (!runData) fChain->SetBranchAddress("muGenGMomPID", muGenGMomPID, &b_muGenGMomPID);
   if (!runData) fChain->SetBranchAddress("muGenMomPID", muGenMomPID, &b_muGenMomPID);
   if (!runData) fChain->SetBranchAddress("muGenMomPt", muGenMomPt, &b_muGenMomPt);
   fChain->SetBranchAddress("muIsoTrk", muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muIsoCalo", muIsoCalo, &b_muIsoCalo);
   fChain->SetBranchAddress("muIsoEcal", muIsoEcal, &b_muIsoEcal);
   fChain->SetBranchAddress("muIsoHcal", muIsoHcal, &b_muIsoHcal);
   fChain->SetBranchAddress("muChi2NDF", muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muEmVeto", muEmVeto, &b_muEmVeto);
   fChain->SetBranchAddress("muHadVeto", muHadVeto, &b_muHadVeto);
   fChain->SetBranchAddress("muType", muType, &b_muType);
   fChain->SetBranchAddress("muID", muID, &b_muID);
   fChain->SetBranchAddress("muD0", muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", muDz, &b_muDz);
   fChain->SetBranchAddress("muPVD0", muPVD0, &b_muPVD0);
   fChain->SetBranchAddress("muPVDz", muPVDz, &b_muPVDz);
   fChain->SetBranchAddress("muValidFraction", muValidFraction, &b_muValidFraction);
   fChain->SetBranchAddress("muTrkdPt", muTrkdPt, &b_muTrkdPt);
   fChain->SetBranchAddress("muNumberOfHits", muNumberOfHits, &b_muNumberOfHits);
   fChain->SetBranchAddress("muNumberOfValidHits", muNumberOfValidHits, &b_muNumberOfValidHits);
   fChain->SetBranchAddress("muNumberOfInactiveHits", muNumberOfInactiveHits, &b_muNumberOfInactiveHits);
   fChain->SetBranchAddress("muNumberOfValidTrkHits", muNumberOfValidTrkHits, &b_muNumberOfValidTrkHits);
   fChain->SetBranchAddress("muNumberOfValidPixelHits", muNumberOfValidPixelHits, &b_muNumberOfValidPixelHits);
   fChain->SetBranchAddress("muNumberOfValidMuonHits", muNumberOfValidMuonHits, &b_muNumberOfValidMuonHits);
   fChain->SetBranchAddress("muStations", muStations, &b_muStations);
   fChain->SetBranchAddress("muChambers", muChambers, &b_muChambers);
   fChain->SetBranchAddress("muPV2D", muPV2D, &b_muPV2D);
   fChain->SetBranchAddress("muPV3D", muPV3D, &b_muPV3D);
   fChain->SetBranchAddress("muBS2D", muBS2D, &b_muBS2D);
   fChain->SetBranchAddress("muBS3D", muBS3D, &b_muBS3D);
   fChain->SetBranchAddress("muVtx", muVtx, &b_muVtx);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetTrg", jetTrg, &b_jetTrg);
   fChain->SetBranchAddress("jetEn", jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetPt", jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetMass", jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetEt", jetEt, &b_jetEt);
   if (!runData) fChain->SetBranchAddress("jetpartonFlavour", jetpartonFlavour, &b_jetpartonFlavour);
   fChain->SetBranchAddress("jetRawPt", jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetCharge", jetCharge, &b_jetCharge);
   fChain->SetBranchAddress("jetNeutralEmEnergy", jetNeutralEmEnergy, &b_jetNeutralEmEnergy);
   fChain->SetBranchAddress("jetNeutralEmEnergyFraction", jetNeutralEmEnergyFraction, &b_jetNeutralEmEnergyFraction);
   fChain->SetBranchAddress("jetNeutralHadronEnergy", jetNeutralHadronEnergy, &b_jetNeutralHadronEnergy);
   fChain->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction, &b_jetNeutralHadronEnergyFraction);
   fChain->SetBranchAddress("jetNConstituents", jetNConstituents, &b_jetNConstituents);
   fChain->SetBranchAddress("jetChargedEmEnergy", jetChargedEmEnergy, &b_jetChargedEmEnergy);
   fChain->SetBranchAddress("jetChargedEmEnergyFraction", jetChargedEmEnergyFraction, &b_jetChargedEmEnergyFraction);
   fChain->SetBranchAddress("jetChargedHadronEnergy", jetChargedHadronEnergy, &b_jetChargedHadronEnergy);
   fChain->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction, &b_jetChargedHadronEnergyFraction);
   fChain->SetBranchAddress("jetChargedHadronMultiplicity", jetChargedHadronMultiplicity, &b_jetChargedHadronMultiplicity);
   fChain->SetBranchAddress("jetChargedMuEnergy", jetChargedMuEnergy, &b_jetChargedMuEnergy);
   fChain->SetBranchAddress("jetChargedMuEnergyFraction", jetChargedMuEnergyFraction, &b_jetChargedMuEnergyFraction);
   if (!runData) fChain->SetBranchAddress("jetGenJetIndex", jetGenJetIndex, &b_jetGenJetIndex);
   if (!runData) fChain->SetBranchAddress("jetGenJetEn", jetGenJetEn, &b_jetGenJetEn);
   if (!runData) fChain->SetBranchAddress("jetGenJetPt", jetGenJetPt, &b_jetGenJetPt);
   if (!runData) fChain->SetBranchAddress("jetGenJetEta", jetGenJetEta, &b_jetGenJetEta);
   if (!runData) fChain->SetBranchAddress("jetGenJetPhi", jetGenJetPhi, &b_jetGenJetPhi);
   if (!runData) fChain->SetBranchAddress("jetGenJetMass", jetGenJetMass, &b_jetGenJetMass);
   if (!runData) fChain->SetBranchAddress("jetGenPartonID", jetGenPartonID, &b_jetGenPartonID);
   if (!runData) fChain->SetBranchAddress("jetGenPartonMomID", jetGenPartonMomID, &b_jetGenPartonMomID);
   fChain->SetBranchAddress("nZee", &nZee, &b_nZee);
   fChain->SetBranchAddress("ZeeMass", ZeeMass, &b_ZeeMass);
//   fChain->SetBranchAddress("ZeePt", ZeePt, &b_ZeePt);
//   fChain->SetBranchAddress("ZeeEta", ZeeEta, &b_ZeeEta);
//   fChain->SetBranchAddress("ZeePhi", ZeePhi, &b_ZeePhi);
//   fChain->SetBranchAddress("ZeeLeg1Index", ZeeLeg1Index, &b_ZeeLeg1Index);
//   fChain->SetBranchAddress("ZeeLeg2Index", ZeeLeg2Index, &b_ZeeLeg2Index);
   fChain->SetBranchAddress("nZmumu", &nZmumu, &b_nZmumu);
   fChain->SetBranchAddress("ZmumuMass", ZmumuMass, &b_ZmumuMass);
//   fChain->SetBranchAddress("ZmumuPt", ZmumuPt, &b_ZmumuPt);
//   fChain->SetBranchAddress("ZmumuEta", ZmumuEta, &b_ZmumuEta);
//   fChain->SetBranchAddress("ZmumuPhi", ZmumuPhi, &b_ZmumuPhi);
//   fChain->SetBranchAddress("ZmumuLeg1Index", ZmumuLeg1Index, &b_ZmumuLeg1Index);
//   fChain->SetBranchAddress("ZmumuLeg2Index", ZmumuLeg2Index, &b_ZmumuLeg2Index);
//   fChain->SetBranchAddress("nWenu", &nWenu, &b_nWenu);
//   fChain->SetBranchAddress("WenuMassTCaloMET", WenuMassTCaloMET, &b_WenuMassTCaloMET);
//   fChain->SetBranchAddress("WenuEtCaloMET", WenuEtCaloMET, &b_WenuEtCaloMET);
//   fChain->SetBranchAddress("WenuACopCaloMET", WenuACopCaloMET, &b_WenuACopCaloMET);
//   fChain->SetBranchAddress("WenuMassTTcMET", WenuMassTTcMET, &b_WenuMassTTcMET);
//   fChain->SetBranchAddress("WenuEtTcMET", WenuEtTcMET, &b_WenuEtTcMET);
//   fChain->SetBranchAddress("WenuACopTcMET", WenuACopTcMET, &b_WenuACopTcMET);
//   fChain->SetBranchAddress("WenuMassTPfMET", WenuMassTPfMET, &b_WenuMassTPfMET);
//   fChain->SetBranchAddress("WenuEtPfMET", WenuEtPfMET, &b_WenuEtPfMET);
//   fChain->SetBranchAddress("WenuACopPfMET", WenuACopPfMET, &b_WenuACopPfMET);
//   fChain->SetBranchAddress("WenuEleIndex", WenuEleIndex, &b_WenuEleIndex);
//   fChain->SetBranchAddress("nWmunu", &nWmunu, &b_nWmunu);
//   fChain->SetBranchAddress("WmunuMassTCaloMET", WmunuMassTCaloMET, &b_WmunuMassTCaloMET);
//   fChain->SetBranchAddress("WmunuEtCaloMET", WmunuEtCaloMET, &b_WmunuEtCaloMET);
//   fChain->SetBranchAddress("WmunuACopCaloMET", WmunuACopCaloMET, &b_WmunuACopCaloMET);
//   fChain->SetBranchAddress("WmunuMassTTcMET", WmunuMassTTcMET, &b_WmunuMassTTcMET);
//   fChain->SetBranchAddress("WmunuEtTcMET", WmunuEtTcMET, &b_WmunuEtTcMET);
//   fChain->SetBranchAddress("WmunuACopTcMET", WmunuACopTcMET, &b_WmunuACopTcMET);
//   fChain->SetBranchAddress("WmunuMassTPfMET", WmunuMassTPfMET, &b_WmunuMassTPfMET);
//   fChain->SetBranchAddress("WmunuEtPfMET", WmunuEtPfMET, &b_WmunuEtPfMET);
//   fChain->SetBranchAddress("WmunuACopPfMET", WmunuACopPfMET, &b_WmunuACopPfMET);
//   fChain->SetBranchAddress("WmunuMuIndex", WmunuMuIndex, &b_WmunuMuIndex);
   Notify();
}

Bool_t anaVgNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void anaVgNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t anaVgNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void anaVgNtuple::InitTree(){
	tree_[0] = new TTree(OO.call_string(string("TreeName")).data(),"");

	// Initialize variables
	t_phoEt            = 0 ; 
	t_phoEta           = 0 ; 
	t_phoSCEta         = 0 ; 
	t_phoPhi           = 0 ; 
	t_phoE             = 0 ; 
	t_phoSigmaIEtaIEta = 0 ; 
	t_phoSCE           = 0 ;
	t_phoInEE          = false ; 
	t_phoInEB          = false ; 
	t_phoIsISR         = false ; 
	t_phoIsFSR         = false ; 
	t_phoIsEle         = false ; 
	t_phoIsPho         = false ; 
	t_ele1Pt           = 0 ; 
	t_ele1Eta          = 0 ; 
	t_ele1SCEta        = 0 ; 
	t_ele1Phi          = 0 ; 
	t_ele1En           = 0 ; 
	t_ele1InEE         = false ; 
	t_ele1InEB         = false ; 
	t_ele2Pt           = 0 ; 
	t_ele2Eta          = 0 ; 
	t_ele2SCEta        = 0 ; 
	t_ele2Phi          = 0 ; 
	t_ele2En           = 0 ; 
	t_ele2InEE         = false ; 
	t_ele2InEB         = false ; 
	t_mu1Pt            = 0 ; 
	t_mu1Eta           = 0 ; 
	t_mu1Phi           = 0 ; 
	t_mu1En            = 0 ; 
	t_mu2Pt            = 0 ; 
	t_mu2Eta           = 0 ; 
	t_mu2Phi           = 0 ; 
	t_mu2En            = 0 ; 
	t_massZ            = 0 ; 
	t_massZg           = 0 ; 
    t_massZZg          = 0 ;
    t_ZPt              = 0 ;
    t_ZgST             = 0 ;
	t_eleTrgLeg1Match  = 0 ;
	t_eleTrgLeg2Match  = 0 ;
	t_CaloIdVL         = 0 ; 
	t_CaloIdL          = 0 ; 
	t_CaloIdXL         = 0 ; 
	t_CaloIdT          = 0 ; 
	t_CaloIdVT         = 0 ; 
	t_CaloIsoVL        = 0 ; 
	t_CaloIsoL         = 0 ; 
	t_CaloIsoXL        = 0 ; 
	t_CaloIsoT         = 0 ; 
	t_CaloIsoVT        = 0 ; 
	t_TrkIdVL          = 0 ; 
	t_TrkIdL           = 0 ; 
	t_TrkIdXL          = 0 ; 
	t_TrkIdT           = 0 ; 
	t_TrkIdVT          = 0 ; 
	t_TrkIsoVL         = 0 ; 
	t_TrkIsoL          = 0 ; 
	t_TrkIsoXL         = 0 ; 
	t_TrkIsoT          = 0 ; 
	t_TrkIsoVT         = 0 ;

	if (doHistoPassPho){
	tree_[0] = new TTree("photonTreePass","store information pass all event selection");
	tree_[0]->Branch("run"             , &run                , "run/I"             );
	tree_[0]->Branch("event"           , &event              , "event/I"           );
	tree_[0]->Branch("nGoodVtx"        , &nGoodVtx           , "nGoodVtx/I"        );
	tree_[0]->Branch("nVtx"            , &nVtx               , "nVtx/I"            );
	tree_[0]->Branch("rho"             , &rho                , "rho/D"             );
	if (!ProcessTag.Contains("Data"))
	tree_[0]->Branch("nPU"             , &nPU[1]             , "nPU/I"             );
	tree_[0]->Branch("phoEt"           , &t_phoEt            , "phoEt/D"           );
	tree_[0]->Branch("phoEta"          , &t_phoEta           , "phoEta/D"          );
	tree_[0]->Branch("phoSCEta"        , &t_phoSCEta         , "phoSCEta/D"        );
	tree_[0]->Branch("phoPhi"          , &t_phoPhi           , "phoPhi/D"          );
	tree_[0]->Branch("phoE"            , &t_phoE             , "phoE/D"            );
	tree_[0]->Branch("phoSCE"          , &t_phoSCE           , "phoSCE/D"          );
	tree_[0]->Branch("phoSigmaIEtaIEta", &t_phoSigmaIEtaIEta , "phoSigmaIEtaIEta/D");
	tree_[0]->Branch("phoInEE"         , &t_phoInEE          , "phoInEE/O"         );
	tree_[0]->Branch("phoInEB"         , &t_phoInEB          , "phoInEB/O"         );
	tree_[0]->Branch("phoIsISR"        , &t_phoIsISR         , "phoIsISR/O"        );
	tree_[0]->Branch("phoIsFSR"        , &t_phoIsFSR         , "phoIsFSR/O"        );
	tree_[0]->Branch("phoIsEle"        , &t_phoIsEle         , "phoIsEle/O"        );
	tree_[0]->Branch("phoIsPho"        , &t_phoIsPho         , "phoIsPho/O"        );
	if (doEle){
	tree_[0]->Branch("ele1Pt"          , &t_ele1Pt           , "ele1Pt/D"          );
	tree_[0]->Branch("ele1Eta"         , &t_ele1Eta          , "ele1Eta/D"         );
	tree_[0]->Branch("ele1SCEta"       , &t_ele1SCEta        , "ele1SCEta/D"       );
	tree_[0]->Branch("ele1Phi"         , &t_ele1Phi          , "ele1Phi/D"         );
	tree_[0]->Branch("ele1En"          , &t_ele1En           , "ele1En/D"          );
	tree_[0]->Branch("ele1InEE"        , &t_ele1InEE         , "ele1InEE/O"        );
	tree_[0]->Branch("ele1InEB"        , &t_ele1InEB         , "ele1InEB/O"        );
	tree_[0]->Branch("ele2Pt"          , &t_ele2Pt           , "ele2Pt/D"          );
	tree_[0]->Branch("ele2Eta"         , &t_ele2Eta          , "ele2Eta/D"         );
	tree_[0]->Branch("ele2SCEta"       , &t_ele2SCEta        , "ele2SCEta/D"       );
	tree_[0]->Branch("ele2Phi"         , &t_ele2Phi          , "ele2Phi/D"         );
	tree_[0]->Branch("ele2En"          , &t_ele2En           , "ele2En/D"          );
	tree_[0]->Branch("ele2InEE"        , &t_ele2InEE         , "ele2InEE/O"        );
	tree_[0]->Branch("ele2InEB"        , &t_ele2InEB         , "ele2InEB/O"        );
	}
	if (doMu){
	tree_[0]->Branch("mu1Pt"           , &t_mu1Pt            , "mu1Pt/D"           );
	tree_[0]->Branch("mu1Eta"          , &t_mu1Eta           , "mu1Eta/D"          );
	tree_[0]->Branch("mu1Phi"          , &t_mu1Phi           , "mu1Phi/D"          );
	tree_[0]->Branch("mu1En"           , &t_mu1En            , "mu1En/D"           );
	tree_[0]->Branch("mu2Pt"           , &t_mu2Pt            , "mu2Pt/D"           );
	tree_[0]->Branch("mu2Eta"          , &t_mu2Eta           , "mu2Eta/D"          );
	tree_[0]->Branch("mu2Phi"          , &t_mu2Phi           , "mu2Phi/D"          );
	tree_[0]->Branch("mu2En"           , &t_mu2En            , "mu2En/D"           );
	}
	tree_[0]->Branch("massZ"           , &t_massZ            , "massZ/D"           );
	tree_[0]->Branch("massZg"          , &t_massZg           , "massZg/D"          );
	tree_[0]->Branch("SampleWeight"    , &SampleWeight       , "SampleWeight/D"    );
	tree_[0]->Branch("PileUpWeight"    , &PileUpWeight       , "PileUpWeight/D"    );
	tree_[0]->Branch("LepIDWeight"     , &LepIDWeight        , "LepIDWeight/D"     );
	tree_[0]->Branch("PhoIDWeight"     , &PhoIDWeight        , "PhoIDWeight/D"     );
	tree_[0]->Branch("JetBkgWeight"    , &JetBkgWeight       , "JetBkgWeight/D"    );
	tree_[0]->Branch("EvtWeight"       , &EvtWeight          , "EvtWeight/D"       );
	}
	if (doTemplate){
	tree_[1] = new TTree("DataTree","Selected pho for template");
	if (!ProcessTag.Contains("Data"))
	tree_[1] = new TTree("MCTree","Selected pho for template");
	tree_[1]->Branch("Pho_Pt"           , &t_phoEt              , "Pho_Pt/D");
	tree_[1]->Branch("Pho_Eta"          , &t_phoSCEta           , "Pho_Eta/D");
	tree_[1]->Branch("Pho_SigmaIEtaIEta", &t_phoSigmaIEtaIEta   , "Pho_SigmaIEtaIEta/D");
	tree_[1]->Branch("weight"           , &EvtWeight            , "weight/D");
	}
	if (doATGC){
	tree_[2] = new TTree("photonTree","");
	tree_[2]->Branch("photonEt"         , &t_phoEt              , "photonEt/D"     );
	tree_[2]->Branch("massZ"            , &t_massZ              , "massZ/D"        );
	tree_[2]->Branch("massZg"           , &t_massZg             , "massZg/D"       );
	tree_[2]->Branch("massZZg"          , &t_massZZg            , "massZZg/D"      );
	tree_[2]->Branch("ZPt"              , &t_ZPt                , "ZPt/D"          );
	tree_[2]->Branch("ZgST"             , &t_ZgST               , "ZgST/D"         );
	tree_[2]->Branch("weight"           , &EvtWeight            , "weight/D"       );
	}
	if (doEleTrgCheck){
	tree_[3] = new TTree("electronTree","electron tree information");
	tree_[3]->Branch("EvtWeight"      , &EvtWeight            ,"EvtWeight/F"       );
	tree_[3]->Branch("nVtx"           , &nVtx                 ,"nVtx/I"            );
	tree_[3]->Branch("nGoodVtx"       , &nGoodVtx             ,"nGoodVtx/I"        );
	tree_[3]->Branch("elePt"          , &t_ele1Pt             , "elePt/D"          );
	tree_[3]->Branch("eleEta"         , &t_ele1Eta            , "eleEta/D"         );
	tree_[3]->Branch("eleSCEta"       , &t_ele1SCEta          , "eleSCEta/D"       );
	tree_[3]->Branch("elePhi"         , &t_ele1Phi            , "elePhi/D"         );
	tree_[3]->Branch("eleEn"          , &t_ele1En             , "eleEn/D"          );
	tree_[3]->Branch("eleInEE"        , &t_ele1InEE           , "eleInEE/O"        );
	tree_[3]->Branch("eleInEB"        , &t_ele1InEB           , "eleInEB/O"        );
	tree_[3]->Branch("isHWW_WP95"     , &t_isHWW_WP95         , "isHWW_WP95/O"     );
	tree_[3]->Branch("isHWW_WP90"     , &t_isHWW_WP90         , "isHWW_WP90/O"     );
	tree_[3]->Branch("isHWW_WP85"     , &t_isHWW_WP85         , "isHWW_WP85/O"     );
	tree_[3]->Branch("isHWW_WP80"     , &t_isHWW_WP80         , "isHWW_WP80/O"     );
	tree_[3]->Branch("isHWW_WP70"     , &t_isHWW_WP70         , "isHWW_WP70/O"     );
	tree_[3]->Branch("isHWW_WP60"     , &t_isHWW_WP60         , "isHWW_WP60/O"     );
	tree_[3]->Branch("CaloIdVL"       , &t_CaloIdVL           , "CaloIdVL/O"       );
	tree_[3]->Branch("CaloIdL"        , &t_CaloIdL            , "CaloIdL/O"        );
	tree_[3]->Branch("CaloIdXL"       , &t_CaloIdXL           , "CaloIdXL/O"       );
	tree_[3]->Branch("CaloIdT"        , &t_CaloIdT            , "CaloIdT/O"        );
	tree_[3]->Branch("CaloIdVT"       , &t_CaloIdVT           , "CaloIdVT/O"       );
	tree_[3]->Branch("CaloIsoVL"      , &t_CaloIsoVL          , "CaloIsoVL/O"      );
	tree_[3]->Branch("CaloIsoL"       , &t_CaloIsoL           , "CaloIsoL/O"       );
	tree_[3]->Branch("CaloIsoXL"      , &t_CaloIsoXL          , "CaloIsoXL/O"      );
	tree_[3]->Branch("CaloIsoT"       , &t_CaloIsoT           , "CaloIsoT/O"       );
	tree_[3]->Branch("CaloIsoVT"      , &t_CaloIsoVT          , "CaloIsoVT/O"      );
	tree_[3]->Branch("TrkIdVL"        , &t_TrkIdVL            , "TrkIdVL/O"        );
	tree_[3]->Branch("TrkIdL"         , &t_TrkIdL             , "TrkIdL/O"         );
	tree_[3]->Branch("TrkIdXL"        , &t_TrkIdXL            , "TrkIdXL/O"        );
	tree_[3]->Branch("TrkIdT"         , &t_TrkIdT             , "TrkIdT/O"         );
	tree_[3]->Branch("TrkIdVT"        , &t_TrkIdVT            , "TrkIdVT/O"        );
	tree_[3]->Branch("TrkIsoVL"       , &t_TrkIsoVL           , "TrkIsoVL/O"       );
	tree_[3]->Branch("TrkIsoL"        , &t_TrkIsoL            , "TrkIsoL/O"        );
	tree_[3]->Branch("TrkIsoXL"       , &t_TrkIsoXL           , "TrkIsoXL/O"       );
	tree_[3]->Branch("TrkIsoT"        , &t_TrkIsoT            , "TrkIsoT/O"        );
	tree_[3]->Branch("TrkIsoVT"       , &t_TrkIsoVT           , "TrkIsoVT/O"       );
	tree_[3]->Branch("eleTrgLeg1Match", &t_eleTrgLeg1Match    , "eleTrgLeg1Match/O");
	tree_[3]->Branch("eleTrgLeg2Match", &t_eleTrgLeg2Match    , "eleTrgLeg2Match/O");
	}
}


void anaVgNtuple::InitHisto(){
	
	if (doHistoPassPho){
		histo[  0] = new TH1D("hpEt"    ,"hpEt"        , 500 , 0, 500);
		histo[  1] = new TH1D("hpEt_EB" ,"hpEt_EB"     , 500 , 0, 500);
		histo[  2] = new TH1D("hpEt_EE" ,"hpEt_EE"     , 500 , 0, 500);
		histo[  3] = new TH1D("hpEta"   ,"hpEta"       , 50, -2.5, 2.5 );
		histo[  4] = new TH1D("hpPhi"   ,"hpPhi"       , 50, -3.15 , 3.15);
		histo[  5] = new TH1D("hpPhi_EB","hpPhi_EB"    , 50,-3.15 , 3.15);
		histo[  6] = new TH1D("hpPhi_EE","hpPhi_EE"    , 50,-3.15 , 3.15);
		histo[  7] = new TH1D("hVtx"    ,"Vtx"         , 30 , -0.5 , 29.5 );
		histo[  8] = new TH1D("hGVtx"   ,"hGVtx"       , 30 , -0.5 , 29.5 );
		histo[  9] = new TH1D("hMll"    ,"hMll"        ,1000,0,1000 );
		histo[ 10] = new TH1D("hMllg"   ,"hMllg"       ,1000,0,1000);
		histo[ 11] = new TH1D("lep1Pt"  ,"lep1Pt"      , 500 , 0, 500);
		histo[ 12] = new TH1D("lep1Eta" ,"lep1Eta"     , 50, -2.5, 2.5 );
		histo[ 13] = new TH1D("lep1Phi" ,"lep1Phi"     , 50,-3.15 , 3.15);
		histo[ 14] = new TH1D("lep1DR"  ,"lep1DR"      , 65, 0 , 6.5 );
		histo[ 15] = new TH1D("lep2Pt"  ,"lep2Pt"      , 500 , 0, 500);
		histo[ 16] = new TH1D("lep2Eta" ,"lep2Eta"     , 50, -2.5, 2.5);
		histo[ 17] = new TH1D("lep2Phi" ,"lep2Phi"     , 50,-3.15 , 3.15);
		histo[ 18] = new TH1D("lep2DR"  ,"lep2DR"      , 65, 0 , 6.5 );

		histo[ 19] = new TH1D("nJets"   ,"nJets"       , 26, -0.5 , 25.5 );

		histo[ 20] = new TH1D("nGVtx1"   ,"nGVtx1"     , 26, -0.5 , 25.5 );
		histo[ 21] = new TH1D("nGVtx2"   ,"nGVtx2"     , 26, -0.5 , 25.5 );
		histo[ 22] = new TH1D("nGVtx3"   ,"nGVtx3"     , 26, -0.5 , 25.5 );
		histo[ 23] = new TH1D("nGVtx4"   ,"nGVtx4"     , 26, -0.5 , 25.5 );
		histo[ 24] = new TH1D("nGVtx5"   ,"nGVtx5"     , 26, -0.5 , 25.5 );
		histo[ 25] = new TH1D("nGVtx6"   ,"nGVtx6"     , 26, -0.5 , 25.5 );
		histo[ 26] = new TH1D("nGVtx7"   ,"nGVtx7"     , 26, -0.5 , 25.5 );
		histo[ 27] = new TH1D("nGVtx8"   ,"nGVtx8"     , 26, -0.5 , 25.5 );
		histo[ 28] = new TH1D("nGVtx9"   ,"nGVtx9"     , 26, -0.5 , 25.5 );
		histo[ 29] = new TH1D("nGVtx10"  ,"nGVtx10"    , 26, -0.5 , 25.5 );
		histo[ 30] = new TH1D("nGVtx11"  ,"nGVtx11"    , 26, -0.5 , 25.5 );
		histo[ 31] = new TH1D("nGVtx12"  ,"nGVtx12"    , 26, -0.5 , 25.5 );

		histo[ 40] = new TH1D("nVtx1"    ,"nVtx1"      , 26, -0.5 , 25.5 );
		histo[ 41] = new TH1D("nVtx2"    ,"nVtx2"      , 26, -0.5 , 25.5 );
		histo[ 42] = new TH1D("nVtx3"    ,"nVtx3"      , 26, -0.5 , 25.5 );
		histo[ 43] = new TH1D("nVtx4"    ,"nVtx4"      , 26, -0.5 , 25.5 );
		histo[ 44] = new TH1D("nVtx5"    ,"nVtx5"      , 26, -0.5 , 25.5 );
		histo[ 45] = new TH1D("nVtx6"    ,"nVtx6"      , 26, -0.5 , 25.5 );
		histo[ 46] = new TH1D("nVtx7"    ,"nVtx7"      , 26, -0.5 , 25.5 );
		histo[ 47] = new TH1D("nVtx8"    ,"nVtx8"      , 26, -0.5 , 25.5 );
		histo[ 48] = new TH1D("nVtx9"    ,"nVtx9"      , 26, -0.5 , 25.5 );
		histo[ 49] = new TH1D("nVtx10"   ,"nVtx10"     , 26, -0.5 , 25.5 );
		histo[ 50] = new TH1D("nVtx11"   ,"nVtx11"     , 26, -0.5 , 25.5 );
		histo[ 51] = new TH1D("nVtx12"   ,"nVtx12"     , 26, -0.5 , 25.5 );

		histo[ 60] = new TH1D("hphoR9"         ,"hphoR9"         , 105, 0. , 1.05  );
		histo[ 61] = new TH1D("EmSCEn"         ,"phoE-SCEn"      , 100, -15. , 5.0 );
		histo[ 62] = new TH1D("hNPU_Z"         ,"hNPU_Z"         , 50 , -0.5 ,49.5 );
		histo[ 63] = new TH1D("hReNPU_Z"       ,"hReNPU_Z"       , 50 , -0.5 ,49.5 );
		histo[ 64] = new TH1D("hNPU_Full"      ,"hNPU_Full"      , 50 , -0.5 ,49.5 );
		histo[ 65] = new TH1D("hReNPU_Full"    ,"hReNPU_Full"    , 50 , -0.5 ,49.5 );

		histo[ 71] = new TH1D("hp5x5TimeHR9EB"  ,""     , 240., -60 , 60. );
		histo[ 72] = new TH1D("hpXtalTimeLR9EB" ,""     , 240., -60 , 60. );
		histo[ 73] = new TH1D("hp5x5TimeHR9EE"  ,""     , 240., -60 , 60. );
		histo[ 74] = new TH1D("hpXtalTimeLR9EE" ,""     , 240., -60 , 60. );

		histo[ 80] = new TH1D("ZmassAll"        ,"ZmassNoPho"         , 100, 50 , 150. );
		histo[ 81] = new TH1D("ZmassBB"         ,"Barrel-Barrel"      , 100, 50 , 150. );
		histo[ 82] = new TH1D("ZmassBE"         ,"Barrel-Endcap"      , 100, 50 , 150. );
		histo[ 83] = new TH1D("ZmassEE"         ,"Endcap-Endcap"      , 100, 50 , 150. );
		histo[ 84] = new TH1D("ZlepPt1"         ,"lepton1 Pt"         , 200, 20 , 220. );
		histo[ 85] = new TH1D("ZlepPt2"         ,"lepton2 Pt"         , 200, 20 , 220. );
		histo[ 86] = new TH1D("ZlepPt1EB"       ,"lepton1 Pt in EB"   , 200, 20 , 220. );
		histo[ 87] = new TH1D("ZlepPt2EB"       ,"lepton2 Pt in EB"   , 200, 20 , 220. );
		histo[ 88] = new TH1D("ZlepPt1EE"       ,"lepton1 Pt in EE"   , 200, 20 , 220. );
		histo[ 89] = new TH1D("ZlepPt2EE"       ,"lepton2 Pt in EE"   , 200, 20 , 220. );
		histo[ 90] = new TH1D("ZmassBB_R9up"    ,"EB_EB_R9>0.94"      , 100, 50 , 150. );

		histo[100] = new TH1D("hPass"           ,"Sequence of passed events",10,0,10);

		histo[110] = new TH1D("nGoodVtx_pre"    ,"nGVtx noSelection"    , 50 , -0.5,  49.5);
		histo[111] = new TH1D("elePt_pre"       ,"elePt_pre noSelection", 200, -0.5, 199.5);
		histo[112] = new TH1D("eleEta_pre"      ,"elePt_pre noSelection", 50 , -2.5, 2.5  );
		histo[113] = new TH1D("phoPt_pre"       ,"elePt_pre noSelection", 200, -0.5, 199.5);
		histo[114] = new TH1D("phoEta_pre"      ,"elePt_pre noSelection", 50 , -2.5, 2.5  );
		histo[115] = new TH1D("lepDR_pre"       ,"lepDR noSelection"    , 65 , 0   , 6.5  );

		// ISR selectioin
		histo[120] = new TH1D("MllgISR"         ,"MllgISR"              , 500,    0, 500  );
		histo[121] = new TH1D("nGVtxISR"        ,"nGgVtxISR"            , 50,  -0.5, 49.5 );
		// FSR selectioin
		histo[122] = new TH1D("MllgFSR"         ,"MllgFSR"              , 500,    0, 500  );
		histo[123] = new TH1D("nGVtxFSR"        ,"nGgVtxFSR"            , 50 , -0.5, 49.5 );
	}
	//	histo[ 19] = new TH1D("","");
	//	histo[ 20] = new TH1D("","");
	if (doHistoPassPho){
		his2D[0] = new TH2D("MeeVsMllg","MeeVsMllg",300,0,300,300,0,300);
		his2D[1] = new TH2D("MeeVsMllgISR","MeeVsMllgISR",300,0,300,300,0,300);
		his2D[2] = new TH2D("MeeVsMllgFSR","MeeVsMllgFSR",300,0,300,300,0,300);
		if (ProcessTag != "Data")
		{
			his2D[3] = new TH2D("nGVtxVSnITPU","nGVtxVSnITPU",25,0.5,50.5,50,-0.5,49.5);
			his2D[4] = new TH2D("nGVtxVSnRho","nGVtxVSnRho"  ,25,0.5,50.5,20,0,20);
		}
	}
	if (doRatioA){
		histo[101] = new TH1D("hTightEB",   "hTightEB"   ,500,0,500);
		histo[102] = new TH1D("hTightEE",   "hTightEE"   ,500,0,500);
		histo[103] = new TH1D("hFakeableEB","hFakeableEB",500,0,500);
		histo[104] = new TH1D("hFakeableEE","hFakeableEE",500,0,500);
	}
	if (doTemplate){
//		double xBin[]={0,15,20,25,30,35,40,60,500};
//		int   nxBin = sizeof(xBin)/sizeof(xBin[0])-1;
		his2D[11] = new TH2D("PhoPtSigmaIEtaIEta_Barrel","PhoPtSigmaIEtaIEta_Barrel",500,0,500,9000,0.0005,0.0905);
		his2D[12] = new TH2D("PhoPtSigmaIEtaIEta_Endcap","PhoPtSigmaIEtaIEta_Endcap",500,0,500,9000,0.0005,0.0905);
	}
}
Bool_t anaVgNtuple::FillHistoPre(){
	Int_t ele1Flag = -1 ;
	Int_t ele2Flag = -1 ;
	Int_t phoFlag  = -1 ;
	for (int iMC = 0 ; iMC < nMC; iMC++){
		// true electron selection
		if ( fabs(mcPID[iMC]) == 11 ){
			if (fabs(mcMomPID[iMC]) <= 23 || fabs(mcGMomPID[iMC]) == 23 )
				if ( fabs(mcEta[iMC]) < 2.5 && mcPt[iMC] > 10. ){
					if (ele1Flag != -1 ) ele2Flag = iMC ;
					else ele1Flag = iMC ;
				}
		}
		if ( fabs(mcPID[iMC]) == 22 ){
			if (fabs(mcMomPID[iMC]) <= 22 || fabs(mcGMomPID[iMC]) <= 23 ){
				
				if ( fabs(mcEta[iMC]) < 2.5 && mcPt[iMC] > 10. ){
					phoFlag = iMC ;
				}
			}
		}
	}
	if (ele1Flag==-1||ele2Flag==-1||phoFlag==-1) return false ;
//	cout << mcEta[ele1Flag] << " " << mcPhi[ele1Flag] << " " << mcEta[phoFlag] << 
//	" " << mcPhi[phoFlag] << endl ;
//	if (DeltaR(mcEta[ele1Flag],mcPhi[ele1Flag],mcEta[phoFlag],mcPhi[phoFlag]) < 0.65) cout << DeltaR(mcEta[ele1Flag],mcPhi[ele1Flag],mcEta[phoFlag],mcPhi[phoFlag]) << endl ;
//	if (DeltaR(mcEta[ele2Flag],mcPhi[ele2Flag],mcEta[phoFlag],mcPhi[phoFlag]) < 0.65) cout << DeltaR(mcEta[ele1Flag],mcPhi[ele1Flag],mcEta[phoFlag],mcPhi[phoFlag]) << endl ;
	if (DeltaR(mcEta[ele2Flag],mcPhi[ele2Flag],mcEta[phoFlag],mcPhi[phoFlag]) < 0.65) return false ;
	if (DeltaR(mcEta[ele1Flag],mcPhi[ele1Flag],mcEta[phoFlag],mcPhi[phoFlag]) < 0.65) return false ;
	TLorentzVector vele1, vele2, vZ, vZg ;
	vele1.SetPtEtaPhiM(mcPt[ele1Flag],mcEta[ele1Flag],mcPhi[ele1Flag],mcMass[ele1Flag]);
	vele2.SetPtEtaPhiM(mcPt[ele2Flag],mcEta[ele2Flag],mcPhi[ele2Flag],mcMass[ele2Flag]);
	vZ  = vele1 + vele2 ;
	if ( vZ.M() < 50 ) return false ;
	histo[110]->Fill(nGoodVtx);
	vector<Int_t> listIndexEle ;
	Int_t pho_index = -1 ;
	for (Int_t iEle = 0; iEle < nEle; iEle++)
	{
		if (DeltaR(mcEta[ele1Flag],mcPhi[ele1Flag],eleEta[iEle],elePhi[iEle]) < 0.3 )
			listIndexEle.push_back(iEle);
		if (DeltaR(mcEta[ele2Flag],mcPhi[ele2Flag],eleEta[iEle],elePhi[iEle]) < 0.3 )
			listIndexEle.push_back(iEle);
	}
	for (Int_t iPho = 0; iPho < nPho; iPho++)
	{
		if (DeltaR(mcEta[phoFlag],mcPhi[phoFlag],phoEta[iPho],phoPhi[iPho]) < 0.3 )
			pho_index = iPho ;
	}
	if ( pho_index != -1 )
	{
		histo[113]->Fill(phoEt[pho_index]);
		histo[114]->Fill(phoEta[pho_index]);
		for (unsigned int iEle = 0; iEle < listIndexEle.size(); iEle++){
			Float_t eEta = eleEta[listIndexEle[iEle]] ;
			Float_t ePhi = elePhi[listIndexEle[iEle]] ;
			Float_t ePt  = elePt [listIndexEle[iEle]] ;
			histo[111]->Fill(ePt);
			histo[112]->Fill(eEta);
			histo[115]->Fill(DeltaR(eEta,ePhi,phoEta[pho_index],phoPhi[pho_index]));
		}
	}
	return true ;
}

void anaVgNtuple::FillHistoZ(Int_t lep1_index , Int_t lep2_index)
{
	if (doHistoPassPho){
		t_ele1SCEta               = eleSCEta[lep1_index] ;
		t_ele2SCEta               = eleSCEta[lep2_index] ;
		t_ele1InEB                =  ( fabs(t_ele1SCEta) < 1.4442 ) ;
		t_ele1InEE                =  ( fabs(t_ele1SCEta) < 2.5 && fabs(t_ele1SCEta) > 1.566 ) ;
		t_ele2InEB                =  ( fabs(t_ele2SCEta) < 1.4442 ) ;
		t_ele2InEE                =  ( fabs(t_ele2SCEta) < 2.5 && fabs(t_ele2SCEta) > 1.566 ) ;
		TLorentzVector vLep1, vLep2, vPho, vZ ;
		vLep1.SetPtEtaPhiE(elePt_[lep1_index],eleEta[lep1_index],elePhi[lep1_index],eleEn_[lep1_index]);
		vLep2.SetPtEtaPhiE(elePt_[lep2_index],eleEta[lep2_index],elePhi[lep2_index],eleEn_[lep2_index]);
		vZ   = vLep1 + vLep2 ;
		histo[84]->Fill(vLep1.Pt());
		histo[85]->Fill(vLep2.Pt());
		histo[80]->Fill(vZ.M(),EvtWeight);
		if ( t_ele1InEB && t_ele2InEB )
		{
			histo[81]->Fill(vZ.M(),EvtWeight);
			if ( eleE3x3[lep1_index] / eleSCRawEn[lep1_index] > 0.94 &&
			     eleE3x3[lep2_index] / eleSCRawEn[lep2_index] > 0.94  )
			histo[90]->Fill(vZ.M(),EvtWeight);
		}
		else if ( t_ele1InEE && t_ele2InEE )
		{
			histo[83]->Fill(vZ.M(),EvtWeight);
		}
		else
		{
			histo[82]->Fill(vZ.M(),EvtWeight);
		}
		if ( t_ele1InEB ) histo[86]->Fill(vLep1.Pt());
		if ( t_ele2InEB ) histo[87]->Fill(vLep2.Pt());
		if ( t_ele1InEE ) histo[88]->Fill(vLep1.Pt());
		if ( t_ele2InEE ) histo[89]->Fill(vLep2.Pt());
	}
}

void anaVgNtuple::FillHisto(Int_t lep1_index , Int_t lep2_index , Int_t pho_index){
	if (doHistoPassPho){
		TLorentzVector vLep1, vLep2, vPho, vZ , vZg ;
		if (doMu){
			float muEn1 = sqrt(muPt[lep1_index]*muPt[lep1_index] + muPz[lep1_index]*muPz[lep1_index]);
			float muEn2 = sqrt(muPt[lep2_index]*muPt[lep2_index] + muPz[lep2_index]*muPz[lep2_index]);
			vLep1.SetPtEtaPhiE(muPt[lep1_index],muEta[lep1_index],muPhi[lep1_index],muEn1);
			vLep2.SetPtEtaPhiE(muPt[lep2_index],muEta[lep2_index],muPhi[lep2_index],muEn2);
		}
		if (doEle){
			vLep1.SetPtEtaPhiE(elePt_[lep1_index],eleEta[lep1_index],elePhi[lep1_index],eleEn_[lep1_index]);
			vLep2.SetPtEtaPhiE(elePt_[lep2_index],eleEta[lep2_index],elePhi[lep2_index],eleEn_[lep2_index]);
		}
		vPho.SetPtEtaPhiE(phoEt_[pho_index],phoEta[pho_index],phoPhi[pho_index],phoE_[pho_index]);
		vZ = vLep1 + vLep2 ;
		vZg = vZ + vPho ;
		// Fill hitograms
		float DR1(0) ;
		float DR2(0) ;
		if (doMu){
			DR1 = DeltaR(muEta[lep1_index],muPhi[lep1_index],phoEta[pho_index],phoPhi[pho_index]);
			DR2 = DeltaR(muEta[lep2_index],muPhi[lep2_index],phoEta[pho_index],phoPhi[pho_index]);
		}
		if (doEle){
			DR1 = DeltaR(eleEta[lep1_index],elePhi[lep1_index],phoEta[pho_index],phoPhi[pho_index]);
			DR2 = DeltaR(eleEta[lep2_index],elePhi[lep2_index],phoEta[pho_index],phoPhi[pho_index]);
		}

		int isEndcap = ( fabs( phoSCEta[pho_index] ) < 1.4442 ) ? 0 : 1 ;
		                      histo[  0]->Fill(phoEt_[pho_index]     ,EvtWeight);
		if (isEndcap == 0 )   histo[  1]->Fill(phoEt_[pho_index]     ,EvtWeight);
		if (isEndcap == 1 )   histo[  2]->Fill(phoEt_[pho_index]     ,EvtWeight);
		                      histo[  3]->Fill(phoEta[pho_index]     ,EvtWeight);
		                      histo[  4]->Fill(phoPhi[pho_index]     ,EvtWeight);
		if (isEndcap == 0 )   histo[  5]->Fill(phoPhi[pho_index]     ,EvtWeight);
		if (isEndcap == 1 )   histo[  6]->Fill(phoPhi[pho_index]     ,EvtWeight);
		                      histo[  7]->Fill(nVtx                  ,EvtWeight);
		                      histo[  8]->Fill(nGoodVtx              ,EvtWeight);
		                      histo[  9]->Fill(vZ.M()                ,EvtWeight);
		                      histo[ 10]->Fill(vZg.M()               ,EvtWeight);
		                      histo[ 11]->Fill(vLep1.Pt()            ,EvtWeight);
		                      histo[ 12]->Fill(vLep1.Eta()           ,EvtWeight);
		                      histo[ 13]->Fill(vLep1.Phi()           ,EvtWeight);
		                      histo[ 14]->Fill(DR1                   ,EvtWeight);
		                      histo[ 15]->Fill(vLep2.Pt()            ,EvtWeight);
		                      histo[ 16]->Fill(vLep2.Eta()           ,EvtWeight);
		                      histo[ 17]->Fill(vLep2.Phi()           ,EvtWeight);
		                      histo[ 18]->Fill(DR2                   ,EvtWeight);
							  histo[ 60]->Fill(phoR9[pho_index]      ,EvtWeight);
							  histo[ 61]->Fill(phoE[pho_index]-
							  				      phoSCE[pho_index]  ,EvtWeight);

		                      his2D[0]->Fill(vZ.M(),vZg.M(),EvtWeight);
		if (IsISR(pho_index)) his2D[1]->Fill(vZ.M(),vZg.M(),EvtWeight);
		if (IsFSR(pho_index)) his2D[2]->Fill(vZ.M(),vZg.M(),EvtWeight);
		if (vZ.M()+vZg.M() > 185.){
			histo[120]->Fill(vZg.M(),EvtWeight);
			histo[121]->Fill(nGoodVtx,EvtWeight);
		}
		else
		{
			histo[122]->Fill(vZg.M(),EvtWeight);
			histo[123]->Fill(nGoodVtx,EvtWeight);
		}
		// Photon 5x5 Time
		// For barrel
		if (fabs(phoSCEta[pho_index]) < 1.4442) {
			// if R9 > 0.94, photon energy is from 5x5
			if (phoR9[pho_index] > 0.94) {
				for (int a=0; a<25; a++) {
					if (pho5x5Energy[pho_index][a] < 2) continue;
					histo[71]->Fill(pho5x5Time[pho_index][a],EvtWeight);
				}
			}
			// if R9 < 0.94, photon energy is from SC
			else {
				for (int a=0; a<phoNxtal[pho_index]; a++) {
					if (phoXtalEnergy[pho_index][a] < 2) continue;
					histo[72]->Fill(phoXtalTime[pho_index][a],EvtWeight);
				}
			}
		}
		// For endcap
		else {
			// if R9 > 0.95, photon energy is from 5x5
			if (phoR9[pho_index] > 0.95) {
				for (int a=0; a<25; a++) {
					if (pho5x5Energy[pho_index][a] < 2) continue;
					histo[73]->Fill(pho5x5Time[pho_index][a]);
				}
			}
			// if R9 < 0.95, photon energy is from 5x5
			else {
				for (int a=0; a<phoNxtal[pho_index]; a++) {
					if (phoXtalEnergy[pho_index][a] < 2) continue;
					histo[74]->Fill(phoXtalTime[pho_index][a]);
				}
			}
		}

	}

	if (doRatioA){
		int isEndcap = ( fabs( phoSCEta[pho_index] ) < 1.4442 ) ? 0 : 1 ;
//		if (isEndcap == 0)     histo[101]->Fill(phoEt_[pho_index],EvtWeight);
//		if (isEndcap == 1)     histo[102]->Fill(phoEt_[pho_index],EvtWeight);
		if (isEndcap == 0)     histo[103]->Fill(phoEt_[pho_index],EvtWeight);
		if (isEndcap == 1)     histo[104]->Fill(phoEt_[pho_index],EvtWeight);
	}
	if (doTemplate){
		int isEndcap = ( fabs( phoSCEta[pho_index] ) < 1.4442 ) ? 0 : 1 ;
		if (isEndcap == 0)     his2D[11] ->Fill(phoEt_[pho_index],phoSigmaIEtaIEta[pho_index],EvtWeight);
		if (isEndcap == 1)     his2D[12] ->Fill(phoEt_[pho_index],phoSigmaIEtaIEta[pho_index],EvtWeight);
	}
}
void anaVgNtuple::FillTree(Int_t lep1_index, Int_t lep2_index, Int_t pho_index){

		TLorentzVector vLep1, vLep2, vPho, vZ , vZg ;
	// Set Photon Variables
		t_phoEt                   = phoEt_[pho_index] ;
		t_phoEta                  = phoEta[pho_index] ;
		t_phoSCEta                = phoSCEta[pho_index] ;
		t_phoPhi                  = phoPhi[pho_index] ;
		t_phoE                    = phoE_[pho_index] ;
		t_phoSCE                  = phoSCE[pho_index] ;
		t_phoSigmaIEtaIEta        = phoSigmaIEtaIEta[pho_index] ;
		t_phoInEB                 = ( fabs(t_phoSCEta) < 1.4442 ) ;
		t_phoInEE                 = ( fabs(t_phoSCEta) < 2.5 && fabs(t_phoSCEta) > 1.566 ) ;
		t_phoIsISR                = IsISR(pho_index);
		t_phoIsFSR                = IsFSR(pho_index);
   		t_phoIsEle                = mcEleMatcher(pho_index,TString("pho"));
   		t_phoIsPho                = mcPhoMatcher(pho_index);
	if (doHistoPassPho || doATGC){
		vPho.SetPtEtaPhiE(t_phoEt,t_phoEta,t_phoPhi,t_phoE);
	// Set Electron Variables
		if (doEle && lep1_index != -1 && lep2_index != -1 ){
		t_ele1Pt                  =  elePt_[lep1_index]   ; // electorn1
		t_ele1Eta                 =  eleEta[lep1_index]  ;
		t_ele1SCEta               =  eleSCEta[lep1_index];
		t_ele1Phi                 =  elePhi[lep1_index]  ;
		t_ele1En                  =  eleEn_[lep1_index]   ;
		t_ele1InEB                =  ( fabs(t_ele1SCEta) < 1.4442 ) ;
		t_ele1InEE                =  ( fabs(t_ele1SCEta) < 2.5 && fabs(t_ele1SCEta) > 1.566 ) ;
		t_ele2Pt                  =  elePt_[lep2_index]   ; // electorn2
		t_ele2Eta                 =  eleEta[lep2_index]  ;
		t_ele2SCEta               =  eleSCEta[lep2_index];
		t_ele2Phi                 =  elePhi[lep2_index]  ;
		t_ele2En                  =  eleEn_[lep2_index]   ;
		t_ele2InEB                =  ( fabs(t_ele2SCEta) < 1.4442 ) ;
		t_ele2InEE                =  ( fabs(t_ele2SCEta) < 2.5 && fabs(t_ele2SCEta) > 1.566 ) ;
		vLep1.SetPtEtaPhiE(t_ele1Pt,t_ele1Eta,t_ele1Phi,t_ele1En);
		vLep2.SetPtEtaPhiE(t_ele2Pt,t_ele2Eta,t_ele2Phi,t_ele2En);
		vZ                        = vLep1 + vLep2 ;
		vZg                       = vZ + vPho ;
		t_massZ                   = vZ.M();
		t_massZg                  = vZg.M();
		t_massZZg                 = t_massZ + t_massZg ;
		t_ZPt                     = vZ.Pt();
		t_ZgST                    = sqrt( t_ZPt*t_ZPt+t_phoEt*t_phoEt ) ;
		}
	// Set Muon Variables
		if (doMu  && lep1_index != -1 && lep2_index != -1 )
		{
		t_mu1Pt                   =  muPt[lep1_index]   ; // electorn1
		t_mu1Eta                  =  muEta[lep1_index]  ;
		t_mu1Phi                  =  muPhi[lep1_index]  ;
		t_mu1En                   =  sqrt(muPt[lep1_index]*muPt[lep1_index] + muPz[lep1_index]*muPz[lep1_index]);
		t_mu2Pt                   =  muPt[lep2_index]   ; // electorn2
		t_mu2Eta                  =  muEta[lep2_index]  ;
		t_mu2Phi                  =  muPhi[lep2_index]  ;
		t_mu2En                   =  sqrt(muPt[lep2_index]*muPt[lep2_index] + muPz[lep2_index]*muPz[lep2_index]);
		vLep1.SetPtEtaPhiE(t_mu1Pt,t_mu1Eta,t_mu1Phi,t_mu1En);
		vLep2.SetPtEtaPhiE(t_mu2Pt,t_mu2Eta,t_mu2Phi,t_mu2En);
		vZ                        = vLep1 + vLep2 ;
		vZg                       = vZ + vPho ;
		t_massZ                   = vZ.M();
		t_massZg                  = vZg.M();
		t_massZZg                 = t_massZ + t_massZg ;
		t_ZPt                     = vZ.Pt();
		t_ZgST                    = sqrt( t_ZPt*t_ZPt+t_phoEt*t_phoEt ) ;
		}
	}
	if ( doHistoPassPho ) tree_[0]->Fill();
	if (doTemplate){
		tree_[1]->Fill();
	}
	if (doATGC){
		tree_[2]->Fill();
	}
}

void anaVgNtuple::FillTreeEleTrg(){
	Bool_t accept = true ;
	Int_t myHLTIndex =   -1 ;
	Int_t eleTrg_index = -1 ;
	if (!CleanEvent()) accept = false ;
	if (isData){
		if (run >= 160431 && run <= 161016){ // HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1
			myHLTIndex   = 111 ;
			eleTrg_index = 10  ;
		} 
		else 
		if (run >= 162762 && run <= 163261){ // HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2
			myHLTIndex   = 112 ;
			eleTrg_index = 11  ;
		} 
		else 
		if (run >= 163270 && run <= 163869){ // HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3
			myHLTIndex   = 113 ;
			eleTrg_index = 12  ;
		}
		else 
		if (run >= 165088 && run <= 165633){ // HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3
			myHLTIndex   = 116 ;
			eleTrg_index = 15;
		}
		else 
		if (run >= 165970 && run <= 166967){ // HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1
			myHLTIndex   = 118 ;
			eleTrg_index = 17  ;
		}
		if (run >= 167039 && run <= 167913){ // HLT_Ele27_WP80_PFMT50_v1
			myHLTIndex   = 190 ;
			eleTrg_index = 27  ;
		}
		if (run >= 170826 && run <= 173198){ // HLT_Ele32_WP70_PFMT50_v3
			myHLTIndex   = 219 ;
			eleTrg_index = 28  ;
		}
		if (run >= 173236 && run <= 173236){ // HLT_Ele32_WP70_PFMT50_v4
			myHLTIndex   = 220 ;
			eleTrg_index = 28  ;
		}
	} else {
			myHLTIndex   = 112 ;
			eleTrg_index = 11  ;
	}
	if (myHLTIndex  == -1 ) return ;
	if (HLT[HLTIndex[myHLTIndex]] != 1 ) accept = false ;
	if (!accept) return ;
//	cout << "myHLTIndex = " << myHLTIndex << " eleTrg_index " << eleTrg_index << endl ; 
	// HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*
	unsigned int leg1 = (1<<_eTrg_CaloIdT) | (1<<_eTrg_TrkIdVL) | (1<<_eTrg_CaloIsoVL) | (1<<_eTrg_TrkIsoVL) ;
	unsigned int leg2 = (1<<_eTrg_CaloIdT) | (1<<_eTrg_TrkIdVL) | (1<<_eTrg_CaloIsoVL) | (1<<_eTrg_TrkIsoVL) ;
//	cout << "leg1 = " << leg1 << " leg2 =  " << leg2 << endl ;
	vector<Int_t> nMatched ;
	for (int iEle = 0; iEle < nEle; iEle++){
		if (eleTrg[iEle][eleTrg_index]==1) nMatched.push_back(iEle) ;
	}
	if (nMatched.size()<=1) return ;
	for (int iEle = 0; iEle < nEle; iEle++){
		if (iEle == nMatched.at(0)) continue ;
		if ( fabs(eleSCEta[iEle]) > 2.5 ) continue ;
		if ( fabs(eleSCEta[iEle]) < 1.566 && fabs(eleSCEta[iEle]) > 1.4442) continue ;
		t_ele1Pt                  =  elePt_[iEle]   ;
		t_ele1Eta                 =  eleEta[iEle]   ;
		t_ele1SCEta               =  eleSCEta[iEle] ;
		t_ele1Phi                 =  elePhi[iEle]   ;
		t_ele1En                  =  eleEn_[iEle]   ;
		t_ele1InEB                =  ( fabs(t_ele1SCEta) < 1.4442 ) ;
		t_ele1InEE                =  ( fabs(t_ele1SCEta) < 2.5 && fabs(t_ele1SCEta) > 1.566 ) ;
		t_isHWW_WP95              = elePassID(iEle,5)  ;
		t_isHWW_WP90              = elePassID(iEle,4)  ;
		t_isHWW_WP85              = elePassID(iEle,3)  ;
		t_isHWW_WP80              = elePassID(iEle,2)  ;
		t_isHWW_WP70              = elePassID(iEle,1)  ;
		t_isHWW_WP60              = elePassID(iEle,0)  ;
		t_CaloIdVL         = (eleTrgWP(iEle) & (1<<_eTrg_CaloIdVL ) ) != 0 ; 
		t_CaloIdL          = (eleTrgWP(iEle) & (1<<_eTrg_CaloIdL  ) ) != 0 ; 
		t_CaloIdXL         = (eleTrgWP(iEle) & (1<<_eTrg_CaloIdXL ) ) != 0 ; 
		t_CaloIdT          = (eleTrgWP(iEle) & (1<<_eTrg_CaloIdT  ) ) != 0 ; 
		t_CaloIdVT         = (eleTrgWP(iEle) & (1<<_eTrg_CaloIdVT ) ) != 0 ; 
		t_CaloIsoVL        = (eleTrgWP(iEle) & (1<<_eTrg_CaloIsoVL) ) != 0 ; 
		t_CaloIsoL         = (eleTrgWP(iEle) & (1<<_eTrg_CaloIsoL ) ) != 0 ; 
		t_CaloIsoXL        = (eleTrgWP(iEle) & (1<<_eTrg_CaloIsoXL) ) != 0 ; 
		t_CaloIsoT         = (eleTrgWP(iEle) & (1<<_eTrg_CaloIsoT ) ) != 0 ; 
		t_CaloIsoVT        = (eleTrgWP(iEle) & (1<<_eTrg_CaloIsoVT) ) != 0 ; 
		t_TrkIdVL          = (eleTrgWP(iEle) & (1<<_eTrg_TrkIdVL  ) ) != 0 ; 
		t_TrkIdL           = (eleTrgWP(iEle) & (1<<_eTrg_TrkIdL   ) ) != 0 ; 
		t_TrkIdXL          = (eleTrgWP(iEle) & (1<<_eTrg_TrkIdXL  ) ) != 0 ; 
		t_TrkIdT           = (eleTrgWP(iEle) & (1<<_eTrg_TrkIdT   ) ) != 0 ; 
		t_TrkIdVT          = (eleTrgWP(iEle) & (1<<_eTrg_TrkIdVT  ) ) != 0 ; 
		t_TrkIsoVL         = (eleTrgWP(iEle) & (1<<_eTrg_TrkIsoVL ) ) != 0 ; 
		t_TrkIsoL          = (eleTrgWP(iEle) & (1<<_eTrg_TrkIsoL  ) ) != 0 ; 
		t_TrkIsoXL         = (eleTrgWP(iEle) & (1<<_eTrg_TrkIsoXL ) ) != 0 ; 
		t_TrkIsoT          = (eleTrgWP(iEle) & (1<<_eTrg_TrkIsoT  ) ) != 0 ; 
		t_TrkIsoVT         = (eleTrgWP(iEle) & (1<<_eTrg_TrkIsoVT ) ) != 0 ; 
		t_eleTrgLeg1Match         = ( elePt_[iEle] > 17. && ( eleTrgWP(iEle) & leg1 ) == leg1 ) ;
		t_eleTrgLeg2Match         = ( elePt_[iEle] >  8. && ( eleTrgWP(iEle) & leg2 ) == leg2 ) ;
		tree_[3]->Fill();
	}
}

void anaVgNtuple::WrHisto(){
	if (doHistoPassPho){
		histo[  0]->Write();
		histo[  1]->Write();
		histo[  2]->Write();
		histo[  3]->Write();
		histo[  4]->Write();
		histo[  5]->Write();
		histo[  6]->Write();
		histo[  7]->Write();
		histo[  8]->Write();
		histo[  9]->Write();
		histo[ 10]->Write();
		histo[ 11]->Write();
		histo[ 12]->Write();
		histo[ 13]->Write();
		histo[ 14]->Write();
		histo[ 15]->Write();
		histo[ 16]->Write();
		histo[ 17]->Write();
		histo[ 18]->Write();
		histo[ 19]->Write();

		histo[ 21]->Write();
		histo[ 41]->Write();

		histo[ 60]->Write();
		histo[ 61]->Write();

		histo[ 71]->Write();
		histo[ 72]->Write();
		histo[ 73]->Write();
		histo[ 74]->Write();

		histo[ 80]->Write();
		histo[ 81]->Write();
		histo[ 82]->Write();
		histo[ 83]->Write();
		histo[ 84]->Write();
		histo[ 85]->Write();
		histo[ 90]->Write();

		histo[100]->Write();
		histo[120]->Write();
		histo[121]->Write();
		histo[122]->Write();
		histo[123]->Write();
	if (doPUWeight){
		histo[ 62]->Write();
		histo[ 63]->Write();
		histo[ 64]->Write();
		histo[ 65]->Write();
	}

	}
	if (doHistoPassPho){
		his2D[0]->Write();
		if (!isData) his2D[1]->Write();
		if (!isData) his2D[2]->Write();
	}
	if (doRatioA){	
		histo[103]->Write();
		histo[104]->Write();
	}
	if (doTemplate    ){
		his2D[11] ->Write();
		his2D[12] ->Write();
	}
	if (doGenInfo1){
		histo[110]->Write();
		histo[111]->Write();
		histo[112]->Write();
		histo[113]->Write();
		histo[114]->Write();
		histo[115]->Write();
	}
}
void anaVgNtuple::WrTree(){
	if (OO.call_bool(string("doCheckFile")) ) cout << "WriteTree1 " << endl ;
	if (doHistoPassPho) tree_[0]->Write() ;
	if (doHistoPassPho) tree_[0]->Clear() ;
	if (OO.call_bool(string("doCheckFile")) ) cout << "WriteTree2 " << endl ;
	if (doTemplate    ) tree_[1]->Write() ;
	if (doTemplate    ) tree_[1]->Clear() ;
	if (OO.call_bool(string("doCheckFile")) ) cout << "WriteTree3 " << endl ;
	if (doATGC        ) tree_[2]->Write() ;
	if (doATGC        ) tree_[2]->Clear() ;
	if (OO.call_bool(string("doCheckFile")) ) cout << "WriteTree4 " << endl ;
	if (doEleTrgCheck ) tree_[3]->Write() ;
	if (doEleTrgCheck ) tree_[3]->Clear() ;
}
void anaVgNtuple::RmHisto(){
	for (int i = 0; i < kMaxHiso1D; i++){
		if (histo[i]){
		if (OO.call_bool(string("doCheckFile"))) cout << "Removing " << i << " 1D histograms" << endl ;
			histo[i]->Delete();
		}
	}
	for (int i = 0; i < kMaxHiso2D; i++){
		if (his2D[i]){
		if (OO.call_bool(string("doCheckFile"))) cout << "Removing " << i << " 2D histograms" << endl ;
			his2D[i]->Delete();
		}
	}
}
#endif // #ifdef anaVgNtuple_cxx
