//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 11 14:06:15 2022 by ROOT version 6.20/04
// from TTree EventTree/Event data
// found on file: flatPtPi0_HiForestAOD_31.root
//////////////////////////////////////////////////////////

#ifndef MVANtuple_h
#define MVANtuple_h
#define MVANtuple_cxx

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "Math/GenVector/PxPyPzE4D.h"

class MVANtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxCastorTower_p4 = 4;

   // Declaration of leaf types
   UInt_t          run;
   ULong64_t       event;
   UInt_t          lumis;
   Bool_t          isData;
   Float_t         rho;
   Int_t           nPUInfo;
   vector<int>     *nPU;
   vector<int>     *puBX;
   vector<float>   *puTrue;
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<int>     *mcStatus;
   vector<float>   *mcVtx_x;
   vector<float>   *mcVtx_y;
   vector<float>   *mcVtx_z;
   vector<float>   *mcPt;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<float>   *mcMass;
   vector<int>     *mcParentage;
   vector<int>     *mcMomPID;
   vector<float>   *mcMomPt;
   vector<float>   *mcMomEta;
   vector<float>   *mcMomPhi;
   vector<float>   *mcMomMass;
   vector<int>     *mcGMomPID;
   vector<int>     *mcIndex;
   vector<float>   *mcCalIsoDR03;
   vector<float>   *mcCalIsoDR04;
   vector<float>   *mcTrkIsoDR03;
   vector<float>   *mcTrkIsoDR04;
   Int_t           nSC;
   vector<float>   *scE;
   vector<float>   *scRawE;
   vector<float>   *scEta;
   vector<float>   *scPhi;
   Int_t           nEle;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<int>     *eleSCPixCharge;
   vector<int>     *eleCtfCharge;
   vector<float>   *eleEn;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *eleIP3D;
   vector<float>   *eleD0Err;
   vector<float>   *eleDzErr;
   vector<float>   *eleIP3DErr;
   vector<float>   *eleTrkPt;
   vector<float>   *eleTrkEta;
   vector<float>   *eleTrkPhi;
   vector<int>     *eleTrkCharge;
   vector<float>   *eleTrkPtErr;
   vector<float>   *eleTrkChi2;
   vector<float>   *eleTrkNdof;
   vector<float>   *eleTrkNormalizedChi2;
   vector<int>     *eleTrkValidHits;
   vector<int>     *eleTrkLayers;
   vector<float>   *elePt;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleSCEn;
   vector<float>   *eleESEn;
   vector<float>   *eleSCEta;
   vector<float>   *eleSCPhi;
   vector<float>   *eleSCRawEn;
   vector<float>   *eleSCEtaWidth;
   vector<float>   *eleSCPhiWidth;
   vector<float>   *eleHoverE;
   vector<float>   *eleHoverEBc;
   vector<float>   *eleEoverP;
   vector<float>   *eleEoverPInv;
   vector<float>   *eleEcalE;
   vector<float>   *elePAtVtx;
   vector<float>   *elePAtSC;
   vector<float>   *elePAtCluster;
   vector<float>   *elePAtSeed;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   vector<float>   *eledEtaSeedAtVtx;
   vector<float>   *eleSigmaIEtaIEta;
   vector<float>   *eleSigmaIEtaIEta_2012;
   vector<float>   *eleSigmaIPhiIPhi;
   vector<int>     *eleConvVeto;
   vector<int>     *eleMissHits;
   vector<float>   *eleESEffSigmaRR;
   vector<float>   *elePFChIso;
   vector<float>   *elePFPhoIso;
   vector<float>   *elePFNeuIso;
   vector<float>   *elePFPUIso;
   vector<float>   *elePFChIso03;
   vector<float>   *elePFPhoIso03;
   vector<float>   *elePFNeuIso03;
   vector<float>   *elePFChIso04;
   vector<float>   *elePFPhoIso04;
   vector<float>   *elePFNeuIso04;
   vector<float>   *elePFRelIsoWithEA;
   vector<float>   *elePFRelIsoWithDBeta;
   vector<float>   *eleEffAreaTimesRho;
   vector<float>   *eleR9;
   vector<float>   *eleE3x3;
   vector<float>   *eleE5x5;
   vector<float>   *eleR9Full5x5;
   vector<float>   *eleE3x3Full5x5;
   vector<float>   *eleE5x5Full5x5;
   vector<int>     *NClusters;
   vector<int>     *NEcalClusters;
   vector<float>   *eleSeedEn;
   vector<float>   *eleSeedEta;
   vector<float>   *eleSeedPhi;
   vector<float>   *eleSeedCryEta;
   vector<float>   *eleSeedCryPhi;
   vector<float>   *eleSeedCryIeta;
   vector<float>   *eleSeedCryIphi;
   vector<float>   *eleBC1E;
   vector<float>   *eleBC1Eta;
   vector<float>   *eleBC2E;
   vector<float>   *eleBC2Eta;
   Int_t           nPho;
   vector<float>   *phoE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoEcorrStdEcal;
   vector<float>   *phoEcorrPhoEcal;
   vector<float>   *phoEcorrRegr1;
   vector<float>   *phoEcorrRegr2;
   vector<float>   *phoEcorrErrStdEcal;
   vector<float>   *phoEcorrErrPhoEcal;
   vector<float>   *phoEcorrErrRegr1;
   vector<float>   *phoEcorrErrRegr2;
   vector<float>   *phoSCE;
   vector<float>   *phoSCEt;
   vector<float>   *phoSCRawE;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<float>   *phoSCBrem;
   vector<int>     *phoSCnHits;
   vector<unsigned int> *phoSCflags;
   vector<int>     *phoSCinClean;
   vector<int>     *phoSCinUnClean;
   vector<int>     *phoSCnBC;
   vector<float>   *phoESEn;
   vector<float>   *phoPSCE;
   vector<float>   *phoPSCRawE;
   vector<float>   *phoPSCEta;
   vector<float>   *phoPSCPhi;
   vector<float>   *phoPSCEtaWidth;
   vector<float>   *phoPSCPhiWidth;
   vector<float>   *phoPSCBrem;
   vector<int>     *phoPSCnHits;
   vector<unsigned int> *phoPSCflags;
   vector<int>     *phoPSCinClean;
   vector<int>     *phoPSCinUnClean;
   vector<int>     *phoPSCnBC;
   vector<float>   *phoPESEn;
   vector<int>     *phoIsPFPhoton;
   vector<int>     *phoIsStandardPhoton;
   vector<int>     *phoHasPixelSeed;
   vector<int>     *phoHasConversionTracks;
   vector<float>   *phoHadTowerOverEm;
   vector<float>   *phoHoverE;
   vector<int>     *phoHoverEValid;
   vector<float>   *phoSigmaIEtaIEta;
   vector<float>   *phoR9;
   vector<float>   *phoE1x5;
   vector<float>   *phoE2x5;
   vector<float>   *phoE3x3;
   vector<float>   *phoE5x5;
   vector<float>   *phoMaxEnergyXtal;
   vector<float>   *phoSigmaEtaEta;
   vector<float>   *phoSigmaIEtaIEta_2012;
   vector<float>   *phoR9_2012;
   vector<float>   *phoE1x5_2012;
   vector<float>   *phoE2x5_2012;
   vector<float>   *phoE3x3_2012;
   vector<float>   *phoE5x5_2012;
   vector<float>   *phoMaxEnergyXtal_2012;
   vector<float>   *phoSigmaEtaEta_2012;
   vector<float>   *phoHadTowerOverEm1;
   vector<float>   *phoHadTowerOverEm2;
   vector<float>   *phoHoverE1;
   vector<float>   *phoHoverE2;
   vector<float>   *phoSigmaIEtaIPhi;
   vector<float>   *phoSigmaIPhiIPhi;
   vector<float>   *phoR1x5;
   vector<float>   *phoR2x5;
   vector<float>   *phoE2nd;
   vector<float>   *phoETop;
   vector<float>   *phoEBottom;
   vector<float>   *phoELeft;
   vector<float>   *phoERight;
   vector<float>   *phoE1x3;
   vector<float>   *phoE2x2;
   vector<float>   *phoE2x5Max;
   vector<float>   *phoE2x5Top;
   vector<float>   *phoE2x5Bottom;
   vector<float>   *phoE2x5Left;
   vector<float>   *phoE2x5Right;
   vector<float>   *phoSigmaIEtaIPhi_2012;
   vector<float>   *phoSigmaIPhiIPhi_2012;
   vector<float>   *phoR1x5_2012;
   vector<float>   *phoR2x5_2012;
   vector<float>   *phoE2nd_2012;
   vector<float>   *phoETop_2012;
   vector<float>   *phoEBottom_2012;
   vector<float>   *phoELeft_2012;
   vector<float>   *phoERight_2012;
   vector<float>   *phoE1x3_2012;
   vector<float>   *phoE2x2_2012;
   vector<float>   *phoE2x5Max_2012;
   vector<float>   *phoE2x5Top_2012;
   vector<float>   *phoE2x5Bottom_2012;
   vector<float>   *phoE2x5Left_2012;
   vector<float>   *phoE2x5Right_2012;
   vector<float>   *phoBC1E;
   vector<float>   *phoBC1Ecorr;
   vector<float>   *phoBC1Eta;
   vector<float>   *phoBC1Phi;
   vector<int>     *phoBC1size;
   vector<unsigned int> *phoBC1flags;
   vector<int>     *phoBC1inClean;
   vector<int>     *phoBC1inUnClean;
   vector<unsigned int> *phoBC1rawID;
   vector<float>   *pho_ecalClusterIsoR2;
   vector<float>   *pho_ecalClusterIsoR3;
   vector<float>   *pho_ecalClusterIsoR4;
   vector<float>   *pho_ecalClusterIsoR5;
   vector<float>   *pho_hcalRechitIsoR1;
   vector<float>   *pho_hcalRechitIsoR2;
   vector<float>   *pho_hcalRechitIsoR3;
   vector<float>   *pho_hcalRechitIsoR4;
   vector<float>   *pho_hcalRechitIsoR5;
   vector<float>   *pho_trackIsoR1PtCut20;
   vector<float>   *pho_trackIsoR2PtCut20;
   vector<float>   *pho_trackIsoR3PtCut20;
   vector<float>   *pho_trackIsoR4PtCut20;
   vector<float>   *pho_trackIsoR5PtCut20;
   vector<float>   *pho_swissCrx;
   vector<float>   *pho_seedTime;
   vector<int>     *pho_genMatchedIndex;
   vector<float>   *pfcIso1;
   vector<float>   *pfcIso2;
   vector<float>   *pfcIso3;
   vector<float>   *pfcIso4;
   vector<float>   *pfcIso5;
   vector<float>   *pfpIso1;
   vector<float>   *pfpIso2;
   vector<float>   *pfpIso3;
   vector<float>   *pfpIso4;
   vector<float>   *pfpIso5;
   vector<float>   *pfnIso1;
   vector<float>   *pfnIso2;
   vector<float>   *pfnIso3;
   vector<float>   *pfnIso4;
   vector<float>   *pfnIso5;
   vector<float>   *pfpIso1subSC;
   vector<float>   *pfpIso2subSC;
   vector<float>   *pfpIso3subSC;
   vector<float>   *pfpIso4subSC;
   vector<float>   *pfpIso5subSC;
   vector<float>   *pfcIso1subUE;
   vector<float>   *pfcIso2subUE;
   vector<float>   *pfcIso3subUE;
   vector<float>   *pfcIso4subUE;
   vector<float>   *pfcIso5subUE;
   vector<float>   *pfpIso1subUE;
   vector<float>   *pfpIso2subUE;
   vector<float>   *pfpIso3subUE;
   vector<float>   *pfpIso4subUE;
   vector<float>   *pfpIso5subUE;
   vector<float>   *pfnIso1subUE;
   vector<float>   *pfnIso2subUE;
   vector<float>   *pfnIso3subUE;
   vector<float>   *pfnIso4subUE;
   vector<float>   *pfnIso5subUE;
   vector<float>   *pfpIso1subSCsubUE;
   vector<float>   *pfpIso2subSCsubUE;
   vector<float>   *pfpIso3subSCsubUE;
   vector<float>   *pfpIso4subSCsubUE;
   vector<float>   *pfpIso5subSCsubUE;
   vector<float>   *pfcIso1pTgt1p0subUE;
   vector<float>   *pfcIso2pTgt1p0subUE;
   vector<float>   *pfcIso3pTgt1p0subUE;
   vector<float>   *pfcIso4pTgt1p0subUE;
   vector<float>   *pfcIso5pTgt1p0subUE;
   vector<float>   *pfcIso1pTgt2p0subUE;
   vector<float>   *pfcIso2pTgt2p0subUE;
   vector<float>   *pfcIso3pTgt2p0subUE;
   vector<float>   *pfcIso4pTgt2p0subUE;
   vector<float>   *pfcIso5pTgt2p0subUE;
   vector<float>   *pfcIso1pTgt3p0subUE;
   vector<float>   *pfcIso2pTgt3p0subUE;
   vector<float>   *pfcIso3pTgt3p0subUE;
   vector<float>   *pfcIso4pTgt3p0subUE;
   vector<float>   *pfcIso5pTgt3p0subUE;
   Int_t           nMu;
   vector<float>   *muPt;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<int>     *muCharge;
   vector<int>     *muType;
   vector<int>     *muIsGood;
   vector<int>     *muIsGlobal;
   vector<int>     *muIsTracker;
   vector<int>     *muIsPF;
   vector<int>     *muIsSTA;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muIP3D;
   vector<float>   *muD0Err;
   vector<float>   *muDzErr;
   vector<float>   *muIP3DErr;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnerD0;
   vector<float>   *muInnerDz;
   vector<float>   *muInnerD0Err;
   vector<float>   *muInnerDzErr;
   vector<float>   *muInnerPt;
   vector<float>   *muInnerPtErr;
   vector<float>   *muInnerEta;
   vector<int>     *muTrkLayers;
   vector<int>     *muPixelLayers;
   vector<int>     *muPixelHits;
   vector<int>     *muMuonHits;
   vector<int>     *muTrkQuality;
   vector<int>     *muStations;
   vector<float>   *muIsoTrk;
   vector<float>   *muPFChIso;
   vector<float>   *muPFPhoIso;
   vector<float>   *muPFNeuIso;
   vector<float>   *muPFPUIso;
   vector<int>     *muIDSoft;
   vector<int>     *muIDLoose;
   vector<int>     *muIDMedium;
   vector<int>     *muIDMediumPrompt;
   vector<int>     *muIDTight;
   vector<int>     *muIDGlobalHighPt;
   vector<int>     *muIDTrkHighPt;
   vector<int>     *muIDInTime;
   Int_t           nTrk;
   vector<float>   *trkPt;
   vector<float>   *trkP;
   vector<float>   *trkEta;
   vector<float>   *trkPhi;
   vector<int>     *trkcharge;
   vector<float>   *trkvx;
   vector<float>   *trkvy;
   vector<float>   *trkvz;
   vector<float>   *trknormchi2;
   vector<float>   *trkchi2;
   vector<float>   *trkd0;
   vector<float>   *trkdxy;
   vector<float>   *trkdz;
   vector<float>   *trkdxyBS;
   vector<float>   *trkdzBS;
   vector<float>   *trkdxyError;
   vector<float>   *trkdzError;
   vector<int>     *trkValidHits;
   vector<int>     *trkMissHits;
   vector<int>     *trkPurity;
   Int_t           nDisplacedTracks;
   Int_t           nTower;
   vector<float>   *CaloTower_hadE;
   vector<float>   *CaloTower_emE;
   vector<float>   *CaloTower_e;
   vector<float>   *CaloTower_et;
   vector<float>   *CaloTower_eta;
   vector<float>   *CaloTower_phi;
   Int_t           nCastorTower;
   vector<float>   *CastorTower_hadE;
   vector<float>   *CastorTower_emE;
   Int_t           CastorTower_p4_;
   Double_t        CastorTower_p4_fCoordinates_fX[kMaxCastorTower_p4];   //[CastorTower_p4_]
   Double_t        CastorTower_p4_fCoordinates_fY[kMaxCastorTower_p4];   //[CastorTower_p4_]
   Double_t        CastorTower_p4_fCoordinates_fZ[kMaxCastorTower_p4];   //[CastorTower_p4_]
   Double_t        CastorTower_p4_fCoordinates_fT[kMaxCastorTower_p4];   //[CastorTower_p4_]
   vector<int>     *CastorTower_NrecHits;
   Int_t           nTrackerHits;
   Int_t           nPixelClusters;
   Int_t           nPixelRecHits;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcVtx_x;   //!
   TBranch        *b_mcVtx_y;   //!
   TBranch        *b_mcVtx_z;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcIndex;   //!
   TBranch        *b_mcCalIsoDR03;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcTrkIsoDR03;   //!
   TBranch        *b_mcTrkIsoDR04;   //!
   TBranch        *b_nSC;   //!
   TBranch        *b_scE;   //!
   TBranch        *b_scRawE;   //!
   TBranch        *b_scEta;   //!
   TBranch        *b_scPhi;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleSCPixCharge;   //!
   TBranch        *b_eleCtfCharge;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleIP3D;   //!
   TBranch        *b_eleD0Err;   //!
   TBranch        *b_eleDzErr;   //!
   TBranch        *b_eleIP3DErr;   //!
   TBranch        *b_eleTrkPt;   //!
   TBranch        *b_eleTrkEta;   //!
   TBranch        *b_eleTrkPhi;   //!
   TBranch        *b_eleTrkCharge;   //!
   TBranch        *b_eleTrkPtErr;   //!
   TBranch        *b_eleTrkChi2;   //!
   TBranch        *b_eleTrkNdof;   //!
   TBranch        *b_eleTrkNormalizedChi2;   //!
   TBranch        *b_eleTrkValidHits;   //!
   TBranch        *b_eleTrkLayers;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleESEn;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleHoverEBc;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleEcalE;   //!
   TBranch        *b_elePAtVtx;   //!
   TBranch        *b_elePAtSC;   //!
   TBranch        *b_elePAtCluster;   //!
   TBranch        *b_elePAtSeed;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eledEtaSeedAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEta;   //!
   TBranch        *b_eleSigmaIEtaIEta_2012;   //!
   TBranch        *b_eleSigmaIPhiIPhi;   //!
   TBranch        *b_eleConvVeto;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_elePFChIso03;   //!
   TBranch        *b_elePFPhoIso03;   //!
   TBranch        *b_elePFNeuIso03;   //!
   TBranch        *b_elePFChIso04;   //!
   TBranch        *b_elePFPhoIso04;   //!
   TBranch        *b_elePFNeuIso04;   //!
   TBranch        *b_elePFRelIsoWithEA;   //!
   TBranch        *b_elePFRelIsoWithDBeta;   //!
   TBranch        *b_eleEffAreaTimesRho;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleE3x3;   //!
   TBranch        *b_eleE5x5;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleE3x3Full5x5;   //!
   TBranch        *b_eleE5x5Full5x5;   //!
   TBranch        *b_NClusters;   //!
   TBranch        *b_NEcalClusters;   //!
   TBranch        *b_eleSeedEn;   //!
   TBranch        *b_eleSeedEta;   //!
   TBranch        *b_eleSeedPhi;   //!
   TBranch        *b_eleSeedCryEta;   //!
   TBranch        *b_eleSeedCryPhi;   //!
   TBranch        *b_eleSeedCryIeta;   //!
   TBranch        *b_eleSeedCryIphi;   //!
   TBranch        *b_eleBC1E;   //!
   TBranch        *b_eleBC1Eta;   //!
   TBranch        *b_eleBC2E;   //!
   TBranch        *b_eleBC2Eta;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEcorrStdEcal;   //!
   TBranch        *b_phoEcorrPhoEcal;   //!
   TBranch        *b_phoEcorrRegr1;   //!
   TBranch        *b_phoEcorrRegr2;   //!
   TBranch        *b_phoEcorrErrStdEcal;   //!
   TBranch        *b_phoEcorrErrPhoEcal;   //!
   TBranch        *b_phoEcorrErrRegr1;   //!
   TBranch        *b_phoEcorrErrRegr2;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCEt;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phoSCnHits;   //!
   TBranch        *b_phoSCflags;   //!
   TBranch        *b_phoSCinClean;   //!
   TBranch        *b_phoSCinUnClean;   //!
   TBranch        *b_phoSCnBC;   //!
   TBranch        *b_phoESEn;   //!
   TBranch        *b_phoPSCE;   //!
   TBranch        *b_phoPSCRawE;   //!
   TBranch        *b_phoPSCEta;   //!
   TBranch        *b_phoPSCPhi;   //!
   TBranch        *b_phoPSCEtaWidth;   //!
   TBranch        *b_phoPSCPhiWidth;   //!
   TBranch        *b_phoPSCBrem;   //!
   TBranch        *b_phoPSCnHits;   //!
   TBranch        *b_phoPSCflags;   //!
   TBranch        *b_phoPSCinClean;   //!
   TBranch        *b_phoPSCinUnClean;   //!
   TBranch        *b_phoPSCnBC;   //!
   TBranch        *b_phoPESEn;   //!
   TBranch        *b_phoIsPFPhoton;   //!
   TBranch        *b_phoIsStandardPhoton;   //!
   TBranch        *b_phoHasPixelSeed;   //!
   TBranch        *b_phoHasConversionTracks;   //!
   TBranch        *b_phoHadTowerOverEm;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoHoverEValid;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoE1x5;   //!
   TBranch        *b_phoE2x5;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoMaxEnergyXtal;   //!
   TBranch        *b_phoSigmaEtaEta;   //!
   TBranch        *b_phoSigmaIEtaIEta_2012;   //!
   TBranch        *b_phoR9_2012;   //!
   TBranch        *b_phoE1x5_2012;   //!
   TBranch        *b_phoE2x5_2012;   //!
   TBranch        *b_phoE3x3_2012;   //!
   TBranch        *b_phoE5x5_2012;   //!
   TBranch        *b_phoMaxEnergyXtal_2012;   //!
   TBranch        *b_phoSigmaEtaEta_2012;   //!
   TBranch        *b_phoHadTowerOverEm1;   //!
   TBranch        *b_phoHadTowerOverEm2;   //!
   TBranch        *b_phoHoverE1;   //!
   TBranch        *b_phoHoverE2;   //!
   TBranch        *b_phoSigmaIEtaIPhi;   //!
   TBranch        *b_phoSigmaIPhiIPhi;   //!
   TBranch        *b_phoR1x5;   //!
   TBranch        *b_phoR2x5;   //!
   TBranch        *b_phoE2nd;   //!
   TBranch        *b_phoETop;   //!
   TBranch        *b_phoEBottom;   //!
   TBranch        *b_phoELeft;   //!
   TBranch        *b_phoERight;   //!
   TBranch        *b_phoE1x3;   //!
   TBranch        *b_phoE2x2;   //!
   TBranch        *b_phoE2x5Max;   //!
   TBranch        *b_phoE2x5Top;   //!
   TBranch        *b_phoE2x5Bottom;   //!
   TBranch        *b_phoE2x5Left;   //!
   TBranch        *b_phoE2x5Right;   //!
   TBranch        *b_phoSigmaIEtaIPhi_2012;   //!
   TBranch        *b_phoSigmaIPhiIPhi_2012;   //!
   TBranch        *b_phoR1x5_2012;   //!
   TBranch        *b_phoR2x5_2012;   //!
   TBranch        *b_phoE2nd_2012;   //!
   TBranch        *b_phoETop_2012;   //!
   TBranch        *b_phoEBottom_2012;   //!
   TBranch        *b_phoELeft_2012;   //!
   TBranch        *b_phoERight_2012;   //!
   TBranch        *b_phoE1x3_2012;   //!
   TBranch        *b_phoE2x2_2012;   //!
   TBranch        *b_phoE2x5Max_2012;   //!
   TBranch        *b_phoE2x5Top_2012;   //!
   TBranch        *b_phoE2x5Bottom_2012;   //!
   TBranch        *b_phoE2x5Left_2012;   //!
   TBranch        *b_phoE2x5Right_2012;   //!
   TBranch        *b_phoBC1E;   //!
   TBranch        *b_phoBC1Ecorr;   //!
   TBranch        *b_phoBC1Eta;   //!
   TBranch        *b_phoBC1Phi;   //!
   TBranch        *b_phoBC1size;   //!
   TBranch        *b_phoBC1flags;   //!
   TBranch        *b_phoBC1inClean;   //!
   TBranch        *b_phoBC1inUnClean;   //!
   TBranch        *b_phoBC1rawID;   //!
   TBranch        *b_pho_ecalClusterIsoR2;   //!
   TBranch        *b_pho_ecalClusterIsoR3;   //!
   TBranch        *b_pho_ecalClusterIsoR4;   //!
   TBranch        *b_pho_ecalClusterIsoR5;   //!
   TBranch        *b_pho_hcalRechitIsoR1;   //!
   TBranch        *b_pho_hcalRechitIsoR2;   //!
   TBranch        *b_pho_hcalRechitIsoR3;   //!
   TBranch        *b_pho_hcalRechitIsoR4;   //!
   TBranch        *b_pho_hcalRechitIsoR5;   //!
   TBranch        *b_pho_trackIsoR1PtCut20;   //!
   TBranch        *b_pho_trackIsoR2PtCut20;   //!
   TBranch        *b_pho_trackIsoR3PtCut20;   //!
   TBranch        *b_pho_trackIsoR4PtCut20;   //!
   TBranch        *b_pho_trackIsoR5PtCut20;   //!
   TBranch        *b_pho_swissCrx;   //!
   TBranch        *b_pho_seedTime;   //!
   TBranch        *b_pho_genMatchedIndex;   //!
   TBranch        *b_pfcIso1;   //!
   TBranch        *b_pfcIso2;   //!
   TBranch        *b_pfcIso3;   //!
   TBranch        *b_pfcIso4;   //!
   TBranch        *b_pfcIso5;   //!
   TBranch        *b_pfpIso1;   //!
   TBranch        *b_pfpIso2;   //!
   TBranch        *b_pfpIso3;   //!
   TBranch        *b_pfpIso4;   //!
   TBranch        *b_pfpIso5;   //!
   TBranch        *b_pfnIso1;   //!
   TBranch        *b_pfnIso2;   //!
   TBranch        *b_pfnIso3;   //!
   TBranch        *b_pfnIso4;   //!
   TBranch        *b_pfnIso5;   //!
   TBranch        *b_pfpIso1subSC;   //!
   TBranch        *b_pfpIso2subSC;   //!
   TBranch        *b_pfpIso3subSC;   //!
   TBranch        *b_pfpIso4subSC;   //!
   TBranch        *b_pfpIso5subSC;   //!
   TBranch        *b_pfcIso1subUE;   //!
   TBranch        *b_pfcIso2subUE;   //!
   TBranch        *b_pfcIso3subUE;   //!
   TBranch        *b_pfcIso4subUE;   //!
   TBranch        *b_pfcIso5subUE;   //!
   TBranch        *b_pfpIso1subUE;   //!
   TBranch        *b_pfpIso2subUE;   //!
   TBranch        *b_pfpIso3subUE;   //!
   TBranch        *b_pfpIso4subUE;   //!
   TBranch        *b_pfpIso5subUE;   //!
   TBranch        *b_pfnIso1subUE;   //!
   TBranch        *b_pfnIso2subUE;   //!
   TBranch        *b_pfnIso3subUE;   //!
   TBranch        *b_pfnIso4subUE;   //!
   TBranch        *b_pfnIso5subUE;   //!
   TBranch        *b_pfpIso1subSCsubUE;   //!
   TBranch        *b_pfpIso2subSCsubUE;   //!
   TBranch        *b_pfpIso3subSCsubUE;   //!
   TBranch        *b_pfpIso4subSCsubUE;   //!
   TBranch        *b_pfpIso5subSCsubUE;   //!
   TBranch        *b_pfcIso1pTgt1p0subUE;   //!
   TBranch        *b_pfcIso2pTgt1p0subUE;   //!
   TBranch        *b_pfcIso3pTgt1p0subUE;   //!
   TBranch        *b_pfcIso4pTgt1p0subUE;   //!
   TBranch        *b_pfcIso5pTgt1p0subUE;   //!
   TBranch        *b_pfcIso1pTgt2p0subUE;   //!
   TBranch        *b_pfcIso2pTgt2p0subUE;   //!
   TBranch        *b_pfcIso3pTgt2p0subUE;   //!
   TBranch        *b_pfcIso4pTgt2p0subUE;   //!
   TBranch        *b_pfcIso5pTgt2p0subUE;   //!
   TBranch        *b_pfcIso1pTgt3p0subUE;   //!
   TBranch        *b_pfcIso2pTgt3p0subUE;   //!
   TBranch        *b_pfcIso3pTgt3p0subUE;   //!
   TBranch        *b_pfcIso4pTgt3p0subUE;   //!
   TBranch        *b_pfcIso5pTgt3p0subUE;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIsGood;   //!
   TBranch        *b_muIsGlobal;   //!
   TBranch        *b_muIsTracker;   //!
   TBranch        *b_muIsPF;   //!
   TBranch        *b_muIsSTA;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muIP3D;   //!
   TBranch        *b_muD0Err;   //!
   TBranch        *b_muDzErr;   //!
   TBranch        *b_muIP3DErr;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muInnerD0Err;   //!
   TBranch        *b_muInnerDzErr;   //!
   TBranch        *b_muInnerPt;   //!
   TBranch        *b_muInnerPtErr;   //!
   TBranch        *b_muInnerEta;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!
   TBranch        *b_muIDSoft;   //!
   TBranch        *b_muIDLoose;   //!
   TBranch        *b_muIDMedium;   //!
   TBranch        *b_muIDMediumPrompt;   //!
   TBranch        *b_muIDTight;   //!
   TBranch        *b_muIDGlobalHighPt;   //!
   TBranch        *b_muIDTrkHighPt;   //!
   TBranch        *b_muIDInTime;   //!
   TBranch        *b_nTrk;   //!
   TBranch        *b_trkPt;   //!
   TBranch        *b_trkP;   //!
   TBranch        *b_trkEta;   //!
   TBranch        *b_trkPhi;   //!
   TBranch        *b_trkcharge;   //!
   TBranch        *b_trkvx;   //!
   TBranch        *b_trkvy;   //!
   TBranch        *b_trkvz;   //!
   TBranch        *b_trknormchi2;   //!
   TBranch        *b_trkchi2;   //!
   TBranch        *b_trkd0;   //!
   TBranch        *b_trkdxy;   //!
   TBranch        *b_trkdz;   //!
   TBranch        *b_trkdxyBS;   //!
   TBranch        *b_trkdzBS;   //!
   TBranch        *b_trkdxyError;   //!
   TBranch        *b_trkdzError;   //!
   TBranch        *b_trkValidHits;   //!
   TBranch        *b_trkMissHits;   //!
   TBranch        *b_trkPurity;   //!
   TBranch        *b_nDisplacedTracks;   //!
   TBranch        *b_nTower;   //!
   TBranch        *b_CaloTower_hadE;   //!
   TBranch        *b_CaloTower_emE;   //!
   TBranch        *b_CaloTower_e;   //!
   TBranch        *b_CaloTower_et;   //!
   TBranch        *b_CaloTower_eta;   //!
   TBranch        *b_CaloTower_phi;   //!
   TBranch        *b_nCastorTower;   //!
   TBranch        *b_CastorTower_hadE;   //!
   TBranch        *b_CastorTower_emE;   //!
   TBranch        *b_CastorTower_p4_;   //!
   TBranch        *b_CastorTower_p4_fCoordinates_fX;   //!
   TBranch        *b_CastorTower_p4_fCoordinates_fY;   //!
   TBranch        *b_CastorTower_p4_fCoordinates_fZ;   //!
   TBranch        *b_CastorTower_p4_fCoordinates_fT;   //!
   TBranch        *b_CastorTower_NrecHits;   //!
   TBranch        *b_nTrackerHits;   //!
   TBranch        *b_nPixelClusters;   //!
   TBranch        *b_nPixelRecHits;   //!

   MVANtuple(TTree *tree=0);
   virtual ~MVANtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MVANtuple_cxx
MVANtuple::MVANtuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     std::cout << "Initializing empty tree !! \n";
   }
   else {
     Init(tree);
   }
}

MVANtuple::~MVANtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MVANtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MVANtuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MVANtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   mcPID = 0;
   mcStatus = 0;
   mcVtx_x = 0;
   mcVtx_y = 0;
   mcVtx_z = 0;
   mcPt = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcMass = 0;
   mcParentage = 0;
   mcMomPID = 0;
   mcMomPt = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcMomMass = 0;
   mcGMomPID = 0;
   mcIndex = 0;
   mcCalIsoDR03 = 0;
   mcCalIsoDR04 = 0;
   mcTrkIsoDR03 = 0;
   mcTrkIsoDR04 = 0;
   scE = 0;
   scRawE = 0;
   scEta = 0;
   scPhi = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleSCPixCharge = 0;
   eleCtfCharge = 0;
   eleEn = 0;
   eleD0 = 0;
   eleDz = 0;
   eleIP3D = 0;
   eleD0Err = 0;
   eleDzErr = 0;
   eleIP3DErr = 0;
   eleTrkPt = 0;
   eleTrkEta = 0;
   eleTrkPhi = 0;
   eleTrkCharge = 0;
   eleTrkPtErr = 0;
   eleTrkChi2 = 0;
   eleTrkNdof = 0;
   eleTrkNormalizedChi2 = 0;
   eleTrkValidHits = 0;
   eleTrkLayers = 0;
   elePt = 0;
   eleEta = 0;
   elePhi = 0;
   eleSCEn = 0;
   eleESEn = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCRawEn = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleHoverE = 0;
   eleHoverEBc = 0;
   eleEoverP = 0;
   eleEoverPInv = 0;
   eleEcalE = 0;
   elePAtVtx = 0;
   elePAtSC = 0;
   elePAtCluster = 0;
   elePAtSeed = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eledEtaSeedAtVtx = 0;
   eleSigmaIEtaIEta = 0;
   eleSigmaIEtaIEta_2012 = 0;
   eleSigmaIPhiIPhi = 0;
   eleConvVeto = 0;
   eleMissHits = 0;
   eleESEffSigmaRR = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   elePFChIso03 = 0;
   elePFPhoIso03 = 0;
   elePFNeuIso03 = 0;
   elePFChIso04 = 0;
   elePFPhoIso04 = 0;
   elePFNeuIso04 = 0;
   elePFRelIsoWithEA = 0;
   elePFRelIsoWithDBeta = 0;
   eleEffAreaTimesRho = 0;
   eleR9 = 0;
   eleE3x3 = 0;
   eleE5x5 = 0;
   eleR9Full5x5 = 0;
   eleE3x3Full5x5 = 0;
   eleE5x5Full5x5 = 0;
   NClusters = 0;
   NEcalClusters = 0;
   eleSeedEn = 0;
   eleSeedEta = 0;
   eleSeedPhi = 0;
   eleSeedCryEta = 0;
   eleSeedCryPhi = 0;
   eleSeedCryIeta = 0;
   eleSeedCryIphi = 0;
   eleBC1E = 0;
   eleBC1Eta = 0;
   eleBC2E = 0;
   eleBC2Eta = 0;
   phoE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoEcorrStdEcal = 0;
   phoEcorrPhoEcal = 0;
   phoEcorrRegr1 = 0;
   phoEcorrRegr2 = 0;
   phoEcorrErrStdEcal = 0;
   phoEcorrErrPhoEcal = 0;
   phoEcorrErrRegr1 = 0;
   phoEcorrErrRegr2 = 0;
   phoSCE = 0;
   phoSCEt = 0;
   phoSCRawE = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phoSCnHits = 0;
   phoSCflags = 0;
   phoSCinClean = 0;
   phoSCinUnClean = 0;
   phoSCnBC = 0;
   phoESEn = 0;
   phoPSCE = 0;
   phoPSCRawE = 0;
   phoPSCEta = 0;
   phoPSCPhi = 0;
   phoPSCEtaWidth = 0;
   phoPSCPhiWidth = 0;
   phoPSCBrem = 0;
   phoPSCnHits = 0;
   phoPSCflags = 0;
   phoPSCinClean = 0;
   phoPSCinUnClean = 0;
   phoPSCnBC = 0;
   phoPESEn = 0;
   phoIsPFPhoton = 0;
   phoIsStandardPhoton = 0;
   phoHasPixelSeed = 0;
   phoHasConversionTracks = 0;
   phoHadTowerOverEm = 0;
   phoHoverE = 0;
   phoHoverEValid = 0;
   phoSigmaIEtaIEta = 0;
   phoR9 = 0;
   phoE1x5 = 0;
   phoE2x5 = 0;
   phoE3x3 = 0;
   phoE5x5 = 0;
   phoMaxEnergyXtal = 0;
   phoSigmaEtaEta = 0;
   phoSigmaIEtaIEta_2012 = 0;
   phoR9_2012 = 0;
   phoE1x5_2012 = 0;
   phoE2x5_2012 = 0;
   phoE3x3_2012 = 0;
   phoE5x5_2012 = 0;
   phoMaxEnergyXtal_2012 = 0;
   phoSigmaEtaEta_2012 = 0;
   phoHadTowerOverEm1 = 0;
   phoHadTowerOverEm2 = 0;
   phoHoverE1 = 0;
   phoHoverE2 = 0;
   phoSigmaIEtaIPhi = 0;
   phoSigmaIPhiIPhi = 0;
   phoR1x5 = 0;
   phoR2x5 = 0;
   phoE2nd = 0;
   phoETop = 0;
   phoEBottom = 0;
   phoELeft = 0;
   phoERight = 0;
   phoE1x3 = 0;
   phoE2x2 = 0;
   phoE2x5Max = 0;
   phoE2x5Top = 0;
   phoE2x5Bottom = 0;
   phoE2x5Left = 0;
   phoE2x5Right = 0;
   phoSigmaIEtaIPhi_2012 = 0;
   phoSigmaIPhiIPhi_2012 = 0;
   phoR1x5_2012 = 0;
   phoR2x5_2012 = 0;
   phoE2nd_2012 = 0;
   phoETop_2012 = 0;
   phoEBottom_2012 = 0;
   phoELeft_2012 = 0;
   phoERight_2012 = 0;
   phoE1x3_2012 = 0;
   phoE2x2_2012 = 0;
   phoE2x5Max_2012 = 0;
   phoE2x5Top_2012 = 0;
   phoE2x5Bottom_2012 = 0;
   phoE2x5Left_2012 = 0;
   phoE2x5Right_2012 = 0;
   phoBC1E = 0;
   phoBC1Ecorr = 0;
   phoBC1Eta = 0;
   phoBC1Phi = 0;
   phoBC1size = 0;
   phoBC1flags = 0;
   phoBC1inClean = 0;
   phoBC1inUnClean = 0;
   phoBC1rawID = 0;
   pho_ecalClusterIsoR2 = 0;
   pho_ecalClusterIsoR3 = 0;
   pho_ecalClusterIsoR4 = 0;
   pho_ecalClusterIsoR5 = 0;
   pho_hcalRechitIsoR1 = 0;
   pho_hcalRechitIsoR2 = 0;
   pho_hcalRechitIsoR3 = 0;
   pho_hcalRechitIsoR4 = 0;
   pho_hcalRechitIsoR5 = 0;
   pho_trackIsoR1PtCut20 = 0;
   pho_trackIsoR2PtCut20 = 0;
   pho_trackIsoR3PtCut20 = 0;
   pho_trackIsoR4PtCut20 = 0;
   pho_trackIsoR5PtCut20 = 0;
   pho_swissCrx = 0;
   pho_seedTime = 0;
   pho_genMatchedIndex = 0;
   pfcIso1 = 0;
   pfcIso2 = 0;
   pfcIso3 = 0;
   pfcIso4 = 0;
   pfcIso5 = 0;
   pfpIso1 = 0;
   pfpIso2 = 0;
   pfpIso3 = 0;
   pfpIso4 = 0;
   pfpIso5 = 0;
   pfnIso1 = 0;
   pfnIso2 = 0;
   pfnIso3 = 0;
   pfnIso4 = 0;
   pfnIso5 = 0;
   pfpIso1subSC = 0;
   pfpIso2subSC = 0;
   pfpIso3subSC = 0;
   pfpIso4subSC = 0;
   pfpIso5subSC = 0;
   pfcIso1subUE = 0;
   pfcIso2subUE = 0;
   pfcIso3subUE = 0;
   pfcIso4subUE = 0;
   pfcIso5subUE = 0;
   pfpIso1subUE = 0;
   pfpIso2subUE = 0;
   pfpIso3subUE = 0;
   pfpIso4subUE = 0;
   pfpIso5subUE = 0;
   pfnIso1subUE = 0;
   pfnIso2subUE = 0;
   pfnIso3subUE = 0;
   pfnIso4subUE = 0;
   pfnIso5subUE = 0;
   pfpIso1subSCsubUE = 0;
   pfpIso2subSCsubUE = 0;
   pfpIso3subSCsubUE = 0;
   pfpIso4subSCsubUE = 0;
   pfpIso5subSCsubUE = 0;
   pfcIso1pTgt1p0subUE = 0;
   pfcIso2pTgt1p0subUE = 0;
   pfcIso3pTgt1p0subUE = 0;
   pfcIso4pTgt1p0subUE = 0;
   pfcIso5pTgt1p0subUE = 0;
   pfcIso1pTgt2p0subUE = 0;
   pfcIso2pTgt2p0subUE = 0;
   pfcIso3pTgt2p0subUE = 0;
   pfcIso4pTgt2p0subUE = 0;
   pfcIso5pTgt2p0subUE = 0;
   pfcIso1pTgt3p0subUE = 0;
   pfcIso2pTgt3p0subUE = 0;
   pfcIso3pTgt3p0subUE = 0;
   pfcIso4pTgt3p0subUE = 0;
   pfcIso5pTgt3p0subUE = 0;
   muPt = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muType = 0;
   muIsGood = 0;
   muIsGlobal = 0;
   muIsTracker = 0;
   muIsPF = 0;
   muIsSTA = 0;
   muD0 = 0;
   muDz = 0;
   muIP3D = 0;
   muD0Err = 0;
   muDzErr = 0;
   muIP3DErr = 0;
   muChi2NDF = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muInnerD0Err = 0;
   muInnerDzErr = 0;
   muInnerPt = 0;
   muInnerPtErr = 0;
   muInnerEta = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muMuonHits = 0;
   muTrkQuality = 0;
   muStations = 0;
   muIsoTrk = 0;
   muPFChIso = 0;
   muPFPhoIso = 0;
   muPFNeuIso = 0;
   muPFPUIso = 0;
   muIDSoft = 0;
   muIDLoose = 0;
   muIDMedium = 0;
   muIDMediumPrompt = 0;
   muIDTight = 0;
   muIDGlobalHighPt = 0;
   muIDTrkHighPt = 0;
   muIDInTime = 0;
   trkPt = 0;
   trkP = 0;
   trkEta = 0;
   trkPhi = 0;
   trkcharge = 0;
   trkvx = 0;
   trkvy = 0;
   trkvz = 0;
   trknormchi2 = 0;
   trkchi2 = 0;
   trkd0 = 0;
   trkdxy = 0;
   trkdz = 0;
   trkdxyBS = 0;
   trkdzBS = 0;
   trkdxyError = 0;
   trkdzError = 0;
   trkValidHits = 0;
   trkMissHits = 0;
   trkPurity = 0;
   CaloTower_hadE = 0;
   CaloTower_emE = 0;
   CaloTower_e = 0;
   CaloTower_et = 0;
   CaloTower_eta = 0;
   CaloTower_phi = 0;
   CastorTower_hadE = 0;
   CastorTower_emE = 0;
   CastorTower_NrecHits = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("mcVtx_x", &mcVtx_x, &b_mcVtx_x);
   fChain->SetBranchAddress("mcVtx_y", &mcVtx_y, &b_mcVtx_y);
   fChain->SetBranchAddress("mcVtx_z", &mcVtx_z, &b_mcVtx_z);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcIndex", &mcIndex, &b_mcIndex);
   fChain->SetBranchAddress("mcCalIsoDR03", &mcCalIsoDR03, &b_mcCalIsoDR03);
   fChain->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04, &b_mcCalIsoDR04);
   fChain->SetBranchAddress("mcTrkIsoDR03", &mcTrkIsoDR03, &b_mcTrkIsoDR03);
   fChain->SetBranchAddress("mcTrkIsoDR04", &mcTrkIsoDR04, &b_mcTrkIsoDR04);
   fChain->SetBranchAddress("nSC", &nSC, &b_nSC);
   fChain->SetBranchAddress("scE", &scE, &b_scE);
   fChain->SetBranchAddress("scRawE", &scRawE, &b_scRawE);
   fChain->SetBranchAddress("scEta", &scEta, &b_scEta);
   fChain->SetBranchAddress("scPhi", &scPhi, &b_scPhi);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleSCPixCharge", &eleSCPixCharge, &b_eleSCPixCharge);
   fChain->SetBranchAddress("eleCtfCharge", &eleCtfCharge, &b_eleCtfCharge);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleIP3D", &eleIP3D, &b_eleIP3D);
   fChain->SetBranchAddress("eleD0Err", &eleD0Err, &b_eleD0Err);
   fChain->SetBranchAddress("eleDzErr", &eleDzErr, &b_eleDzErr);
   fChain->SetBranchAddress("eleIP3DErr", &eleIP3DErr, &b_eleIP3DErr);
   fChain->SetBranchAddress("eleTrkPt", &eleTrkPt, &b_eleTrkPt);
   fChain->SetBranchAddress("eleTrkEta", &eleTrkEta, &b_eleTrkEta);
   fChain->SetBranchAddress("eleTrkPhi", &eleTrkPhi, &b_eleTrkPhi);
   fChain->SetBranchAddress("eleTrkCharge", &eleTrkCharge, &b_eleTrkCharge);
   fChain->SetBranchAddress("eleTrkPtErr", &eleTrkPtErr, &b_eleTrkPtErr);
   fChain->SetBranchAddress("eleTrkChi2", &eleTrkChi2, &b_eleTrkChi2);
   fChain->SetBranchAddress("eleTrkNdof", &eleTrkNdof, &b_eleTrkNdof);
   fChain->SetBranchAddress("eleTrkNormalizedChi2", &eleTrkNormalizedChi2, &b_eleTrkNormalizedChi2);
   fChain->SetBranchAddress("eleTrkValidHits", &eleTrkValidHits, &b_eleTrkValidHits);
   fChain->SetBranchAddress("eleTrkLayers", &eleTrkLayers, &b_eleTrkLayers);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleESEn", &eleESEn, &b_eleESEn);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleHoverEBc", &eleHoverEBc, &b_eleHoverEBc);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleEcalE", &eleEcalE, &b_eleEcalE);
   fChain->SetBranchAddress("elePAtVtx", &elePAtVtx, &b_elePAtVtx);
   fChain->SetBranchAddress("elePAtSC", &elePAtSC, &b_elePAtSC);
   fChain->SetBranchAddress("elePAtCluster", &elePAtCluster, &b_elePAtCluster);
   fChain->SetBranchAddress("elePAtSeed", &elePAtSeed, &b_elePAtSeed);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eledEtaSeedAtVtx", &eledEtaSeedAtVtx, &b_eledEtaSeedAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEta", &eleSigmaIEtaIEta, &b_eleSigmaIEtaIEta);
   fChain->SetBranchAddress("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012, &b_eleSigmaIEtaIEta_2012);
   fChain->SetBranchAddress("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi, &b_eleSigmaIPhiIPhi);
   fChain->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleESEffSigmaRR", &eleESEffSigmaRR, &b_eleESEffSigmaRR);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("elePFChIso03", &elePFChIso03, &b_elePFChIso03);
   fChain->SetBranchAddress("elePFPhoIso03", &elePFPhoIso03, &b_elePFPhoIso03);
   fChain->SetBranchAddress("elePFNeuIso03", &elePFNeuIso03, &b_elePFNeuIso03);
   fChain->SetBranchAddress("elePFChIso04", &elePFChIso04, &b_elePFChIso04);
   fChain->SetBranchAddress("elePFPhoIso04", &elePFPhoIso04, &b_elePFPhoIso04);
   fChain->SetBranchAddress("elePFNeuIso04", &elePFNeuIso04, &b_elePFNeuIso04);
   fChain->SetBranchAddress("elePFRelIsoWithEA", &elePFRelIsoWithEA, &b_elePFRelIsoWithEA);
   fChain->SetBranchAddress("elePFRelIsoWithDBeta", &elePFRelIsoWithDBeta, &b_elePFRelIsoWithDBeta);
   fChain->SetBranchAddress("eleEffAreaTimesRho", &eleEffAreaTimesRho, &b_eleEffAreaTimesRho);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleE3x3", &eleE3x3, &b_eleE3x3);
   fChain->SetBranchAddress("eleE5x5", &eleE5x5, &b_eleE5x5);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleE3x3Full5x5", &eleE3x3Full5x5, &b_eleE3x3Full5x5);
   fChain->SetBranchAddress("eleE5x5Full5x5", &eleE5x5Full5x5, &b_eleE5x5Full5x5);
   fChain->SetBranchAddress("NClusters", &NClusters, &b_NClusters);
   fChain->SetBranchAddress("NEcalClusters", &NEcalClusters, &b_NEcalClusters);
   fChain->SetBranchAddress("eleSeedEn", &eleSeedEn, &b_eleSeedEn);
   fChain->SetBranchAddress("eleSeedEta", &eleSeedEta, &b_eleSeedEta);
   fChain->SetBranchAddress("eleSeedPhi", &eleSeedPhi, &b_eleSeedPhi);
   fChain->SetBranchAddress("eleSeedCryEta", &eleSeedCryEta, &b_eleSeedCryEta);
   fChain->SetBranchAddress("eleSeedCryPhi", &eleSeedCryPhi, &b_eleSeedCryPhi);
   fChain->SetBranchAddress("eleSeedCryIeta", &eleSeedCryIeta, &b_eleSeedCryIeta);
   fChain->SetBranchAddress("eleSeedCryIphi", &eleSeedCryIphi, &b_eleSeedCryIphi);
   fChain->SetBranchAddress("eleBC1E", &eleBC1E, &b_eleBC1E);
   fChain->SetBranchAddress("eleBC1Eta", &eleBC1Eta, &b_eleBC1Eta);
   fChain->SetBranchAddress("eleBC2E", &eleBC2E, &b_eleBC2E);
   fChain->SetBranchAddress("eleBC2Eta", &eleBC2Eta, &b_eleBC2Eta);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEcorrStdEcal", &phoEcorrStdEcal, &b_phoEcorrStdEcal);
   fChain->SetBranchAddress("phoEcorrPhoEcal", &phoEcorrPhoEcal, &b_phoEcorrPhoEcal);
   fChain->SetBranchAddress("phoEcorrRegr1", &phoEcorrRegr1, &b_phoEcorrRegr1);
   fChain->SetBranchAddress("phoEcorrRegr2", &phoEcorrRegr2, &b_phoEcorrRegr2);
   fChain->SetBranchAddress("phoEcorrErrStdEcal", &phoEcorrErrStdEcal, &b_phoEcorrErrStdEcal);
   fChain->SetBranchAddress("phoEcorrErrPhoEcal", &phoEcorrErrPhoEcal, &b_phoEcorrErrPhoEcal);
   fChain->SetBranchAddress("phoEcorrErrRegr1", &phoEcorrErrRegr1, &b_phoEcorrErrRegr1);
   fChain->SetBranchAddress("phoEcorrErrRegr2", &phoEcorrErrRegr2, &b_phoEcorrErrRegr2);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCEt", &phoSCEt, &b_phoSCEt);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phoSCnHits", &phoSCnHits, &b_phoSCnHits);
   fChain->SetBranchAddress("phoSCflags", &phoSCflags, &b_phoSCflags);
   fChain->SetBranchAddress("phoSCinClean", &phoSCinClean, &b_phoSCinClean);
   fChain->SetBranchAddress("phoSCinUnClean", &phoSCinUnClean, &b_phoSCinUnClean);
   fChain->SetBranchAddress("phoSCnBC", &phoSCnBC, &b_phoSCnBC);
   fChain->SetBranchAddress("phoESEn", &phoESEn, &b_phoESEn);
   fChain->SetBranchAddress("phoPSCE", &phoPSCE, &b_phoPSCE);
   fChain->SetBranchAddress("phoPSCRawE", &phoPSCRawE, &b_phoPSCRawE);
   fChain->SetBranchAddress("phoPSCEta", &phoPSCEta, &b_phoPSCEta);
   fChain->SetBranchAddress("phoPSCPhi", &phoPSCPhi, &b_phoPSCPhi);
   fChain->SetBranchAddress("phoPSCEtaWidth", &phoPSCEtaWidth, &b_phoPSCEtaWidth);
   fChain->SetBranchAddress("phoPSCPhiWidth", &phoPSCPhiWidth, &b_phoPSCPhiWidth);
   fChain->SetBranchAddress("phoPSCBrem", &phoPSCBrem, &b_phoPSCBrem);
   fChain->SetBranchAddress("phoPSCnHits", &phoPSCnHits, &b_phoPSCnHits);
   fChain->SetBranchAddress("phoPSCflags", &phoPSCflags, &b_phoPSCflags);
   fChain->SetBranchAddress("phoPSCinClean", &phoPSCinClean, &b_phoPSCinClean);
   fChain->SetBranchAddress("phoPSCinUnClean", &phoPSCinUnClean, &b_phoPSCinUnClean);
   fChain->SetBranchAddress("phoPSCnBC", &phoPSCnBC, &b_phoPSCnBC);
   fChain->SetBranchAddress("phoPESEn", &phoPESEn, &b_phoPESEn);
   fChain->SetBranchAddress("phoIsPFPhoton", &phoIsPFPhoton, &b_phoIsPFPhoton);
   fChain->SetBranchAddress("phoIsStandardPhoton", &phoIsStandardPhoton, &b_phoIsStandardPhoton);
   fChain->SetBranchAddress("phoHasPixelSeed", &phoHasPixelSeed, &b_phoHasPixelSeed);
   fChain->SetBranchAddress("phoHasConversionTracks", &phoHasConversionTracks, &b_phoHasConversionTracks);
   fChain->SetBranchAddress("phoHadTowerOverEm", &phoHadTowerOverEm, &b_phoHadTowerOverEm);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoHoverEValid", &phoHoverEValid, &b_phoHoverEValid);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoE1x5", &phoE1x5, &b_phoE1x5);
   fChain->SetBranchAddress("phoE2x5", &phoE2x5, &b_phoE2x5);
   fChain->SetBranchAddress("phoE3x3", &phoE3x3, &b_phoE3x3);
   fChain->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoMaxEnergyXtal", &phoMaxEnergyXtal, &b_phoMaxEnergyXtal);
   fChain->SetBranchAddress("phoSigmaEtaEta", &phoSigmaEtaEta, &b_phoSigmaEtaEta);
   fChain->SetBranchAddress("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012, &b_phoSigmaIEtaIEta_2012);
   fChain->SetBranchAddress("phoR9_2012", &phoR9_2012, &b_phoR9_2012);
   fChain->SetBranchAddress("phoE1x5_2012", &phoE1x5_2012, &b_phoE1x5_2012);
   fChain->SetBranchAddress("phoE2x5_2012", &phoE2x5_2012, &b_phoE2x5_2012);
   fChain->SetBranchAddress("phoE3x3_2012", &phoE3x3_2012, &b_phoE3x3_2012);
   fChain->SetBranchAddress("phoE5x5_2012", &phoE5x5_2012, &b_phoE5x5_2012);
   fChain->SetBranchAddress("phoMaxEnergyXtal_2012", &phoMaxEnergyXtal_2012, &b_phoMaxEnergyXtal_2012);
   fChain->SetBranchAddress("phoSigmaEtaEta_2012", &phoSigmaEtaEta_2012, &b_phoSigmaEtaEta_2012);
   fChain->SetBranchAddress("phoHadTowerOverEm1", &phoHadTowerOverEm1, &b_phoHadTowerOverEm1);
   fChain->SetBranchAddress("phoHadTowerOverEm2", &phoHadTowerOverEm2, &b_phoHadTowerOverEm2);
   fChain->SetBranchAddress("phoHoverE1", &phoHoverE1, &b_phoHoverE1);
   fChain->SetBranchAddress("phoHoverE2", &phoHoverE2, &b_phoHoverE2);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi, &b_phoSigmaIEtaIPhi);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi, &b_phoSigmaIPhiIPhi);
   fChain->SetBranchAddress("phoR1x5", &phoR1x5, &b_phoR1x5);
   fChain->SetBranchAddress("phoR2x5", &phoR2x5, &b_phoR2x5);
   fChain->SetBranchAddress("phoE2nd", &phoE2nd, &b_phoE2nd);
   fChain->SetBranchAddress("phoETop", &phoETop, &b_phoETop);
   fChain->SetBranchAddress("phoEBottom", &phoEBottom, &b_phoEBottom);
   fChain->SetBranchAddress("phoELeft", &phoELeft, &b_phoELeft);
   fChain->SetBranchAddress("phoERight", &phoERight, &b_phoERight);
   fChain->SetBranchAddress("phoE1x3", &phoE1x3, &b_phoE1x3);
   fChain->SetBranchAddress("phoE2x2", &phoE2x2, &b_phoE2x2);
   fChain->SetBranchAddress("phoE2x5Max", &phoE2x5Max, &b_phoE2x5Max);
   fChain->SetBranchAddress("phoE2x5Top", &phoE2x5Top, &b_phoE2x5Top);
   fChain->SetBranchAddress("phoE2x5Bottom", &phoE2x5Bottom, &b_phoE2x5Bottom);
   fChain->SetBranchAddress("phoE2x5Left", &phoE2x5Left, &b_phoE2x5Left);
   fChain->SetBranchAddress("phoE2x5Right", &phoE2x5Right, &b_phoE2x5Right);
   fChain->SetBranchAddress("phoSigmaIEtaIPhi_2012", &phoSigmaIEtaIPhi_2012, &b_phoSigmaIEtaIPhi_2012);
   fChain->SetBranchAddress("phoSigmaIPhiIPhi_2012", &phoSigmaIPhiIPhi_2012, &b_phoSigmaIPhiIPhi_2012);
   fChain->SetBranchAddress("phoR1x5_2012", &phoR1x5_2012, &b_phoR1x5_2012);
   fChain->SetBranchAddress("phoR2x5_2012", &phoR2x5_2012, &b_phoR2x5_2012);
   fChain->SetBranchAddress("phoE2nd_2012", &phoE2nd_2012, &b_phoE2nd_2012);
   fChain->SetBranchAddress("phoETop_2012", &phoETop_2012, &b_phoETop_2012);
   fChain->SetBranchAddress("phoEBottom_2012", &phoEBottom_2012, &b_phoEBottom_2012);
   fChain->SetBranchAddress("phoELeft_2012", &phoELeft_2012, &b_phoELeft_2012);
   fChain->SetBranchAddress("phoERight_2012", &phoERight_2012, &b_phoERight_2012);
   fChain->SetBranchAddress("phoE1x3_2012", &phoE1x3_2012, &b_phoE1x3_2012);
   fChain->SetBranchAddress("phoE2x2_2012", &phoE2x2_2012, &b_phoE2x2_2012);
   fChain->SetBranchAddress("phoE2x5Max_2012", &phoE2x5Max_2012, &b_phoE2x5Max_2012);
   fChain->SetBranchAddress("phoE2x5Top_2012", &phoE2x5Top_2012, &b_phoE2x5Top_2012);
   fChain->SetBranchAddress("phoE2x5Bottom_2012", &phoE2x5Bottom_2012, &b_phoE2x5Bottom_2012);
   fChain->SetBranchAddress("phoE2x5Left_2012", &phoE2x5Left_2012, &b_phoE2x5Left_2012);
   fChain->SetBranchAddress("phoE2x5Right_2012", &phoE2x5Right_2012, &b_phoE2x5Right_2012);
   fChain->SetBranchAddress("phoBC1E", &phoBC1E, &b_phoBC1E);
   fChain->SetBranchAddress("phoBC1Ecorr", &phoBC1Ecorr, &b_phoBC1Ecorr);
   fChain->SetBranchAddress("phoBC1Eta", &phoBC1Eta, &b_phoBC1Eta);
   fChain->SetBranchAddress("phoBC1Phi", &phoBC1Phi, &b_phoBC1Phi);
   fChain->SetBranchAddress("phoBC1size", &phoBC1size, &b_phoBC1size);
   fChain->SetBranchAddress("phoBC1flags", &phoBC1flags, &b_phoBC1flags);
   fChain->SetBranchAddress("phoBC1inClean", &phoBC1inClean, &b_phoBC1inClean);
   fChain->SetBranchAddress("phoBC1inUnClean", &phoBC1inUnClean, &b_phoBC1inUnClean);
   fChain->SetBranchAddress("phoBC1rawID", &phoBC1rawID, &b_phoBC1rawID);
   fChain->SetBranchAddress("pho_ecalClusterIsoR2", &pho_ecalClusterIsoR2, &b_pho_ecalClusterIsoR2);
   fChain->SetBranchAddress("pho_ecalClusterIsoR3", &pho_ecalClusterIsoR3, &b_pho_ecalClusterIsoR3);
   fChain->SetBranchAddress("pho_ecalClusterIsoR4", &pho_ecalClusterIsoR4, &b_pho_ecalClusterIsoR4);
   fChain->SetBranchAddress("pho_ecalClusterIsoR5", &pho_ecalClusterIsoR5, &b_pho_ecalClusterIsoR5);
   fChain->SetBranchAddress("pho_hcalRechitIsoR1", &pho_hcalRechitIsoR1, &b_pho_hcalRechitIsoR1);
   fChain->SetBranchAddress("pho_hcalRechitIsoR2", &pho_hcalRechitIsoR2, &b_pho_hcalRechitIsoR2);
   fChain->SetBranchAddress("pho_hcalRechitIsoR3", &pho_hcalRechitIsoR3, &b_pho_hcalRechitIsoR3);
   fChain->SetBranchAddress("pho_hcalRechitIsoR4", &pho_hcalRechitIsoR4, &b_pho_hcalRechitIsoR4);
   fChain->SetBranchAddress("pho_hcalRechitIsoR5", &pho_hcalRechitIsoR5, &b_pho_hcalRechitIsoR5);
   fChain->SetBranchAddress("pho_trackIsoR1PtCut20", &pho_trackIsoR1PtCut20, &b_pho_trackIsoR1PtCut20);
   fChain->SetBranchAddress("pho_trackIsoR2PtCut20", &pho_trackIsoR2PtCut20, &b_pho_trackIsoR2PtCut20);
   fChain->SetBranchAddress("pho_trackIsoR3PtCut20", &pho_trackIsoR3PtCut20, &b_pho_trackIsoR3PtCut20);
   fChain->SetBranchAddress("pho_trackIsoR4PtCut20", &pho_trackIsoR4PtCut20, &b_pho_trackIsoR4PtCut20);
   fChain->SetBranchAddress("pho_trackIsoR5PtCut20", &pho_trackIsoR5PtCut20, &b_pho_trackIsoR5PtCut20);
   fChain->SetBranchAddress("pho_swissCrx", &pho_swissCrx, &b_pho_swissCrx);
   fChain->SetBranchAddress("pho_seedTime", &pho_seedTime, &b_pho_seedTime);
   fChain->SetBranchAddress("pho_genMatchedIndex", &pho_genMatchedIndex, &b_pho_genMatchedIndex);
   fChain->SetBranchAddress("pfcIso1", &pfcIso1, &b_pfcIso1);
   fChain->SetBranchAddress("pfcIso2", &pfcIso2, &b_pfcIso2);
   fChain->SetBranchAddress("pfcIso3", &pfcIso3, &b_pfcIso3);
   fChain->SetBranchAddress("pfcIso4", &pfcIso4, &b_pfcIso4);
   fChain->SetBranchAddress("pfcIso5", &pfcIso5, &b_pfcIso5);
   fChain->SetBranchAddress("pfpIso1", &pfpIso1, &b_pfpIso1);
   fChain->SetBranchAddress("pfpIso2", &pfpIso2, &b_pfpIso2);
   fChain->SetBranchAddress("pfpIso3", &pfpIso3, &b_pfpIso3);
   fChain->SetBranchAddress("pfpIso4", &pfpIso4, &b_pfpIso4);
   fChain->SetBranchAddress("pfpIso5", &pfpIso5, &b_pfpIso5);
   fChain->SetBranchAddress("pfnIso1", &pfnIso1, &b_pfnIso1);
   fChain->SetBranchAddress("pfnIso2", &pfnIso2, &b_pfnIso2);
   fChain->SetBranchAddress("pfnIso3", &pfnIso3, &b_pfnIso3);
   fChain->SetBranchAddress("pfnIso4", &pfnIso4, &b_pfnIso4);
   fChain->SetBranchAddress("pfnIso5", &pfnIso5, &b_pfnIso5);
   fChain->SetBranchAddress("pfpIso1subSC", &pfpIso1subSC, &b_pfpIso1subSC);
   fChain->SetBranchAddress("pfpIso2subSC", &pfpIso2subSC, &b_pfpIso2subSC);
   fChain->SetBranchAddress("pfpIso3subSC", &pfpIso3subSC, &b_pfpIso3subSC);
   fChain->SetBranchAddress("pfpIso4subSC", &pfpIso4subSC, &b_pfpIso4subSC);
   fChain->SetBranchAddress("pfpIso5subSC", &pfpIso5subSC, &b_pfpIso5subSC);
   fChain->SetBranchAddress("pfcIso1subUE", &pfcIso1subUE, &b_pfcIso1subUE);
   fChain->SetBranchAddress("pfcIso2subUE", &pfcIso2subUE, &b_pfcIso2subUE);
   fChain->SetBranchAddress("pfcIso3subUE", &pfcIso3subUE, &b_pfcIso3subUE);
   fChain->SetBranchAddress("pfcIso4subUE", &pfcIso4subUE, &b_pfcIso4subUE);
   fChain->SetBranchAddress("pfcIso5subUE", &pfcIso5subUE, &b_pfcIso5subUE);
   fChain->SetBranchAddress("pfpIso1subUE", &pfpIso1subUE, &b_pfpIso1subUE);
   fChain->SetBranchAddress("pfpIso2subUE", &pfpIso2subUE, &b_pfpIso2subUE);
   fChain->SetBranchAddress("pfpIso3subUE", &pfpIso3subUE, &b_pfpIso3subUE);
   fChain->SetBranchAddress("pfpIso4subUE", &pfpIso4subUE, &b_pfpIso4subUE);
   fChain->SetBranchAddress("pfpIso5subUE", &pfpIso5subUE, &b_pfpIso5subUE);
   fChain->SetBranchAddress("pfnIso1subUE", &pfnIso1subUE, &b_pfnIso1subUE);
   fChain->SetBranchAddress("pfnIso2subUE", &pfnIso2subUE, &b_pfnIso2subUE);
   fChain->SetBranchAddress("pfnIso3subUE", &pfnIso3subUE, &b_pfnIso3subUE);
   fChain->SetBranchAddress("pfnIso4subUE", &pfnIso4subUE, &b_pfnIso4subUE);
   fChain->SetBranchAddress("pfnIso5subUE", &pfnIso5subUE, &b_pfnIso5subUE);
   fChain->SetBranchAddress("pfpIso1subSCsubUE", &pfpIso1subSCsubUE, &b_pfpIso1subSCsubUE);
   fChain->SetBranchAddress("pfpIso2subSCsubUE", &pfpIso2subSCsubUE, &b_pfpIso2subSCsubUE);
   fChain->SetBranchAddress("pfpIso3subSCsubUE", &pfpIso3subSCsubUE, &b_pfpIso3subSCsubUE);
   fChain->SetBranchAddress("pfpIso4subSCsubUE", &pfpIso4subSCsubUE, &b_pfpIso4subSCsubUE);
   fChain->SetBranchAddress("pfpIso5subSCsubUE", &pfpIso5subSCsubUE, &b_pfpIso5subSCsubUE);
   fChain->SetBranchAddress("pfcIso1pTgt1p0subUE", &pfcIso1pTgt1p0subUE, &b_pfcIso1pTgt1p0subUE);
   fChain->SetBranchAddress("pfcIso2pTgt1p0subUE", &pfcIso2pTgt1p0subUE, &b_pfcIso2pTgt1p0subUE);
   fChain->SetBranchAddress("pfcIso3pTgt1p0subUE", &pfcIso3pTgt1p0subUE, &b_pfcIso3pTgt1p0subUE);
   fChain->SetBranchAddress("pfcIso4pTgt1p0subUE", &pfcIso4pTgt1p0subUE, &b_pfcIso4pTgt1p0subUE);
   fChain->SetBranchAddress("pfcIso5pTgt1p0subUE", &pfcIso5pTgt1p0subUE, &b_pfcIso5pTgt1p0subUE);
   fChain->SetBranchAddress("pfcIso1pTgt2p0subUE", &pfcIso1pTgt2p0subUE, &b_pfcIso1pTgt2p0subUE);
   fChain->SetBranchAddress("pfcIso2pTgt2p0subUE", &pfcIso2pTgt2p0subUE, &b_pfcIso2pTgt2p0subUE);
   fChain->SetBranchAddress("pfcIso3pTgt2p0subUE", &pfcIso3pTgt2p0subUE, &b_pfcIso3pTgt2p0subUE);
   fChain->SetBranchAddress("pfcIso4pTgt2p0subUE", &pfcIso4pTgt2p0subUE, &b_pfcIso4pTgt2p0subUE);
   fChain->SetBranchAddress("pfcIso5pTgt2p0subUE", &pfcIso5pTgt2p0subUE, &b_pfcIso5pTgt2p0subUE);
   fChain->SetBranchAddress("pfcIso1pTgt3p0subUE", &pfcIso1pTgt3p0subUE, &b_pfcIso1pTgt3p0subUE);
   fChain->SetBranchAddress("pfcIso2pTgt3p0subUE", &pfcIso2pTgt3p0subUE, &b_pfcIso2pTgt3p0subUE);
   fChain->SetBranchAddress("pfcIso3pTgt3p0subUE", &pfcIso3pTgt3p0subUE, &b_pfcIso3pTgt3p0subUE);
   fChain->SetBranchAddress("pfcIso4pTgt3p0subUE", &pfcIso4pTgt3p0subUE, &b_pfcIso4pTgt3p0subUE);
   fChain->SetBranchAddress("pfcIso5pTgt3p0subUE", &pfcIso5pTgt3p0subUE, &b_pfcIso5pTgt3p0subUE);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muIsGood", &muIsGood, &b_muIsGood);
   fChain->SetBranchAddress("muIsGlobal", &muIsGlobal, &b_muIsGlobal);
   fChain->SetBranchAddress("muIsTracker", &muIsTracker, &b_muIsTracker);
   fChain->SetBranchAddress("muIsPF", &muIsPF, &b_muIsPF);
   fChain->SetBranchAddress("muIsSTA", &muIsSTA, &b_muIsSTA);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muIP3D", &muIP3D, &b_muIP3D);
   fChain->SetBranchAddress("muD0Err", &muD0Err, &b_muD0Err);
   fChain->SetBranchAddress("muDzErr", &muDzErr, &b_muDzErr);
   fChain->SetBranchAddress("muIP3DErr", &muIP3DErr, &b_muIP3DErr);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muInnerD0Err", &muInnerD0Err, &b_muInnerD0Err);
   fChain->SetBranchAddress("muInnerDzErr", &muInnerDzErr, &b_muInnerDzErr);
   fChain->SetBranchAddress("muInnerPt", &muInnerPt, &b_muInnerPt);
   fChain->SetBranchAddress("muInnerPtErr", &muInnerPtErr, &b_muInnerPtErr);
   fChain->SetBranchAddress("muInnerEta", &muInnerEta, &b_muInnerEta);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
   fChain->SetBranchAddress("muIDSoft", &muIDSoft, &b_muIDSoft);
   fChain->SetBranchAddress("muIDLoose", &muIDLoose, &b_muIDLoose);
   fChain->SetBranchAddress("muIDMedium", &muIDMedium, &b_muIDMedium);
   fChain->SetBranchAddress("muIDMediumPrompt", &muIDMediumPrompt, &b_muIDMediumPrompt);
   fChain->SetBranchAddress("muIDTight", &muIDTight, &b_muIDTight);
   fChain->SetBranchAddress("muIDGlobalHighPt", &muIDGlobalHighPt, &b_muIDGlobalHighPt);
   fChain->SetBranchAddress("muIDTrkHighPt", &muIDTrkHighPt, &b_muIDTrkHighPt);
   fChain->SetBranchAddress("muIDInTime", &muIDInTime, &b_muIDInTime);
   fChain->SetBranchAddress("nTrk", &nTrk, &b_nTrk);
   fChain->SetBranchAddress("trkPt", &trkPt, &b_trkPt);
   fChain->SetBranchAddress("trkP", &trkP, &b_trkP);
   fChain->SetBranchAddress("trkEta", &trkEta, &b_trkEta);
   fChain->SetBranchAddress("trkPhi", &trkPhi, &b_trkPhi);
   fChain->SetBranchAddress("trkcharge", &trkcharge, &b_trkcharge);
   fChain->SetBranchAddress("trkvx", &trkvx, &b_trkvx);
   fChain->SetBranchAddress("trkvy", &trkvy, &b_trkvy);
   fChain->SetBranchAddress("trkvz", &trkvz, &b_trkvz);
   fChain->SetBranchAddress("trknormchi2", &trknormchi2, &b_trknormchi2);
   fChain->SetBranchAddress("trkchi2", &trkchi2, &b_trkchi2);
   fChain->SetBranchAddress("trkd0", &trkd0, &b_trkd0);
   fChain->SetBranchAddress("trkdxy", &trkdxy, &b_trkdxy);
   fChain->SetBranchAddress("trkdz", &trkdz, &b_trkdz);
   fChain->SetBranchAddress("trkdxyBS", &trkdxyBS, &b_trkdxyBS);
   fChain->SetBranchAddress("trkdzBS", &trkdzBS, &b_trkdzBS);
   fChain->SetBranchAddress("trkdxyError", &trkdxyError, &b_trkdxyError);
   fChain->SetBranchAddress("trkdzError", &trkdzError, &b_trkdzError);
   fChain->SetBranchAddress("trkValidHits", &trkValidHits, &b_trkValidHits);
   fChain->SetBranchAddress("trkMissHits", &trkMissHits, &b_trkMissHits);
   fChain->SetBranchAddress("trkPurity", &trkPurity, &b_trkPurity);
   fChain->SetBranchAddress("nDisplacedTracks", &nDisplacedTracks, &b_nDisplacedTracks);
   fChain->SetBranchAddress("nTower", &nTower, &b_nTower);
   fChain->SetBranchAddress("CaloTower_hadE", &CaloTower_hadE, &b_CaloTower_hadE);
   fChain->SetBranchAddress("CaloTower_emE", &CaloTower_emE, &b_CaloTower_emE);
   fChain->SetBranchAddress("CaloTower_e", &CaloTower_e, &b_CaloTower_e);
   fChain->SetBranchAddress("CaloTower_et", &CaloTower_et, &b_CaloTower_et);
   fChain->SetBranchAddress("CaloTower_eta", &CaloTower_eta, &b_CaloTower_eta);
   fChain->SetBranchAddress("CaloTower_phi", &CaloTower_phi, &b_CaloTower_phi);
   fChain->SetBranchAddress("nCastorTower", &nCastorTower, &b_nCastorTower);
   fChain->SetBranchAddress("CastorTower_hadE", &CastorTower_hadE, &b_CastorTower_hadE);
   fChain->SetBranchAddress("CastorTower_emE", &CastorTower_emE, &b_CastorTower_emE);
   fChain->SetBranchAddress("CastorTower_p4", &CastorTower_p4_, &b_CastorTower_p4_);
   fChain->SetBranchAddress("CastorTower_p4.fCoordinates.fX", CastorTower_p4_fCoordinates_fX, &b_CastorTower_p4_fCoordinates_fX);
   fChain->SetBranchAddress("CastorTower_p4.fCoordinates.fY", CastorTower_p4_fCoordinates_fY, &b_CastorTower_p4_fCoordinates_fY);
   fChain->SetBranchAddress("CastorTower_p4.fCoordinates.fZ", CastorTower_p4_fCoordinates_fZ, &b_CastorTower_p4_fCoordinates_fZ);
   fChain->SetBranchAddress("CastorTower_p4.fCoordinates.fT", CastorTower_p4_fCoordinates_fT, &b_CastorTower_p4_fCoordinates_fT);
   fChain->SetBranchAddress("CastorTower_NrecHits", &CastorTower_NrecHits, &b_CastorTower_NrecHits);
   fChain->SetBranchAddress("nTrackerHits", &nTrackerHits, &b_nTrackerHits);
   fChain->SetBranchAddress("nPixelClusters", &nPixelClusters, &b_nPixelClusters);
   fChain->SetBranchAddress("nPixelRecHits", &nPixelRecHits, &b_nPixelRecHits);
   Notify();
}

Bool_t MVANtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MVANtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MVANtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MVANtuple_cxx
