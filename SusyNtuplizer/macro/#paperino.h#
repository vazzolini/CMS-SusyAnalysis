//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 18 12:37:20 2012 by ROOT version 5.32/00
// from TTree susyTree/SUSY Event
// found on file: /tmp/azzolini/susyEvents.root
//////////////////////////////////////////////////////////

#ifndef paperino_h
#define paperino_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TObject.h>
#include <TVector3.h>
#include <TVector2.h>
#include <TLorentzVector.h>

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxl1Map = 119;
const Int_t kMaxhltMap = 450;
const Int_t kMaxmetMap = 2;
const Int_t kMaxvertices = 26;
const Int_t kMaxtracks = 174;
const Int_t kMaxsuperClusters = 7;
const Int_t kMaxclusters = 41;
const Int_t kMaxmuons = 3;
const Int_t kMaxgeneralTracks = 1;
const Int_t kMaxpu = 1;
const Int_t kMaxsimVertices = 1;
const Int_t kMaxgenParticles = 1;
const Int_t kMaxgridParams = 1;

class paperino {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //susy::Event     *susyEvent;
   UChar_t         isRealData;
   Int_t           runNumber;
   ULong_t         eventNumber;
   Int_t           luminosityBlockNumber;
   Int_t           bunchCrossing;
   Float_t         avgInsRecLumi;
   Float_t         intgRecLumi;
   UChar_t         cosmicFlag;
   Float_t         rho;
   Float_t         rhoBarrel;
   Float_t         rho25;
   Int_t           metFilterBit;
   UInt_t          beamSpot_fUniqueID;
   UInt_t          beamSpot_fBits;
   Double_t        beamSpot_fX;
   Double_t        beamSpot_fY;
   Double_t        beamSpot_fZ;
   Int_t           l1Map_;
   TString         l1Map_first[kMaxl1Map];
   Int_t           l1Map_second_first[kMaxl1Map];   //[l1Map_]
   UChar_t         l1Map_second_second[kMaxl1Map];   //[l1Map_]
   Int_t           hltMap_;
   TString         hltMap_first[kMaxhltMap];
   Int_t           hltMap_second_first[kMaxhltMap];   //[hltMap_]
   UChar_t         hltMap_second_second[kMaxhltMap];   //[hltMap_]
   Int_t           metMap_;
   TString         metMap_first[kMaxmetMap];
   Float_t         metMap_second_sumEt[kMaxmetMap];   //[metMap_]
   Float_t         metMap_second_significance[kMaxmetMap];   //[metMap_]
   UInt_t          metMap_second_mEt_fUniqueID[kMaxmetMap];   //[metMap_]
   UInt_t          metMap_second_mEt_fBits[kMaxmetMap];   //[metMap_]
   Double_t        metMap_second_mEt_fX[kMaxmetMap];   //[metMap_]
   Double_t        metMap_second_mEt_fY[kMaxmetMap];   //[metMap_]
   UInt_t          metMap_second_vertex_fUniqueID[kMaxmetMap];   //[metMap_]
   UInt_t          metMap_second_vertex_fBits[kMaxmetMap];   //[metMap_]
   Double_t        metMap_second_vertex_fX[kMaxmetMap];   //[metMap_]
   Double_t        metMap_second_vertex_fY[kMaxmetMap];   //[metMap_]
   Double_t        metMap_second_vertex_fZ[kMaxmetMap];   //[metMap_]
 //vector<susy::CorrMETData> metMap_second_mEtCorr[kMaxmetMap];
   Int_t           vertices_;
   Float_t         vertices_chi2[kMaxvertices];   //[vertices_]
   Float_t         vertices_ndof[kMaxvertices];   //[vertices_]
   UChar_t         vertices_tracksSize[kMaxvertices];   //[vertices_]
   UInt_t          vertices_position_fUniqueID[kMaxvertices];   //[vertices_]
   UInt_t          vertices_position_fBits[kMaxvertices];   //[vertices_]
   Double_t        vertices_position_fX[kMaxvertices];   //[vertices_]
   Double_t        vertices_position_fY[kMaxvertices];   //[vertices_]
   Double_t        vertices_position_fZ[kMaxvertices];   //[vertices_]
   Int_t           tracks_;
   Int_t           tracks_algorithm[kMaxtracks];   //[tracks_]
   Int_t           tracks_quality[kMaxtracks];   //[tracks_]
   UChar_t         tracks_numberOfValidHits[kMaxtracks];   //[tracks_]
   UChar_t         tracks_numberOfValidTrackerHits[kMaxtracks];   //[tracks_]
   UChar_t         tracks_numberOfValidMuonHits[kMaxtracks];   //[tracks_]
   UChar_t         tracks_numberOfValidPixelHits[kMaxtracks];   //[tracks_]
   UChar_t         tracks_numberOfValidStripHits[kMaxtracks];   //[tracks_]
   Float_t         tracks_chi2[kMaxtracks];   //[tracks_]
   Float_t         tracks_ndof[kMaxtracks];   //[tracks_]
   Float_t         tracks_charge[kMaxtracks];   //[tracks_]
   Float_t         tracks_error[kMaxtracks][5];   //[tracks_]
   UInt_t          tracks_vertex_fUniqueID[kMaxtracks];   //[tracks_]
   UInt_t          tracks_vertex_fBits[kMaxtracks];   //[tracks_]
   Double_t        tracks_vertex_fX[kMaxtracks];   //[tracks_]
   Double_t        tracks_vertex_fY[kMaxtracks];   //[tracks_]
   Double_t        tracks_vertex_fZ[kMaxtracks];   //[tracks_]
   UInt_t          tracks_momentum_fUniqueID[kMaxtracks];   //[tracks_]
   UInt_t          tracks_momentum_fBits[kMaxtracks];   //[tracks_]
   UInt_t          tracks_momentum_fP_fUniqueID[kMaxtracks];   //[tracks_]
   UInt_t          tracks_momentum_fP_fBits[kMaxtracks];   //[tracks_]
   Double_t        tracks_momentum_fP_fX[kMaxtracks];   //[tracks_]
   Double_t        tracks_momentum_fP_fY[kMaxtracks];   //[tracks_]
   Double_t        tracks_momentum_fP_fZ[kMaxtracks];   //[tracks_]
   Double_t        tracks_momentum_fE[kMaxtracks];   //[tracks_]
 //map<TString,TVector3> tracks_extrapolatedPositions[kMaxtracks];
   Int_t           superClusters_;
   Short_t         superClusters_seedClusterIndex[kMaxsuperClusters];   //[superClusters_]
   Float_t         superClusters_energy[kMaxsuperClusters];   //[superClusters_]
   Float_t         superClusters_preshowerEnergy[kMaxsuperClusters];   //[superClusters_]
   Float_t         superClusters_phiWidth[kMaxsuperClusters];   //[superClusters_]
   Float_t         superClusters_etaWidth[kMaxsuperClusters];   //[superClusters_]
   UInt_t          superClusters_position_fUniqueID[kMaxsuperClusters];   //[superClusters_]
   UInt_t          superClusters_position_fBits[kMaxsuperClusters];   //[superClusters_]
   Double_t        superClusters_position_fX[kMaxsuperClusters];   //[superClusters_]
   Double_t        superClusters_position_fY[kMaxsuperClusters];   //[superClusters_]
   Double_t        superClusters_position_fZ[kMaxsuperClusters];   //[superClusters_]
   vector<Int_t>   superClusters_basicClusterIndices[kMaxsuperClusters];
   Int_t           clusters_;
   UChar_t         clusters_nCrystals[kMaxclusters];   //[clusters_]
   Float_t         clusters_energy[kMaxclusters];   //[clusters_]
   UInt_t          clusters_position_fUniqueID[kMaxclusters];   //[clusters_]
   UInt_t          clusters_position_fBits[kMaxclusters];   //[clusters_]
   Double_t        clusters_position_fX[kMaxclusters];   //[clusters_]
   Double_t        clusters_position_fY[kMaxclusters];   //[clusters_]
   Double_t        clusters_position_fZ[kMaxclusters];   //[clusters_]
   Int_t           muons_;
   UChar_t         muons_type[kMaxmuons];   //[muons_]
   UChar_t         muons_nMatches[kMaxmuons];   //[muons_]
   UChar_t         muons_nValidHits[kMaxmuons];   //[muons_]
   UChar_t         muons_nValidTrackerHits[kMaxmuons];   //[muons_]
   UChar_t         muons_nValidMuonHits[kMaxmuons];   //[muons_]
   UChar_t         muons_nChambers[kMaxmuons];   //[muons_]
   UChar_t         muons_timeNDof[kMaxmuons];   //[muons_]
   Char_t          muons_timeDirection[kMaxmuons];   //[muons_]
   Float_t         muons_timeAtIp[kMaxmuons];   //[muons_]
   Float_t         muons_timeAtIpError[kMaxmuons];   //[muons_]
   Float_t         muons_caloCompatibility[kMaxmuons];   //[muons_]
   Float_t         muons_emEnergy[kMaxmuons];   //[muons_]
   Float_t         muons_hadEnergy[kMaxmuons];   //[muons_]
   Float_t         muons_trackIsoR03[kMaxmuons];   //[muons_]
   Float_t         muons_ecalIsoR03[kMaxmuons];   //[muons_]
   Float_t         muons_hcalIsoR03[kMaxmuons];   //[muons_]
   Float_t         muons_trackIsoR05[kMaxmuons];   //[muons_]
   Float_t         muons_ecalIsoR05[kMaxmuons];   //[muons_]
   Float_t         muons_hcalIsoR05[kMaxmuons];   //[muons_]
   Float_t         muons_sumChargedHadronPt03[kMaxmuons];   //[muons_]
   Float_t         muons_sumChargedParticlePt03[kMaxmuons];   //[muons_]
   Float_t         muons_sumNeutralHadronEt03[kMaxmuons];   //[muons_]
   Float_t         muons_sumPhotonEt03[kMaxmuons];   //[muons_]
   Float_t         muons_sumNeutralHadronEtHighThreshold03[kMaxmuons];   //[muons_]
   Float_t         muons_sumPhotonEtHighThreshold03[kMaxmuons];   //[muons_]
   Float_t         muons_sumPUPt03[kMaxmuons];   //[muons_]
   Float_t         muons_sumChargedHadronPt04[kMaxmuons];   //[muons_]
   Float_t         muons_sumChargedParticlePt04[kMaxmuons];   //[muons_]
   Float_t         muons_sumNeutralHadronEt04[kMaxmuons];   //[muons_]
   Float_t         muons_sumPhotonEt04[kMaxmuons];   //[muons_]
   Float_t         muons_sumNeutralHadronEtHighThreshold04[kMaxmuons];   //[muons_]
   Float_t         muons_sumPhotonEtHighThreshold04[kMaxmuons];   //[muons_]
   Float_t         muons_sumPUPt04[kMaxmuons];   //[muons_]
   Short_t         muons_trackIndex[kMaxmuons];   //[muons_]
   Short_t         muons_standAloneTrackIndex[kMaxmuons];   //[muons_]
   Short_t         muons_combinedTrackIndex[kMaxmuons];   //[muons_]
   UInt_t          muons_momentum_fUniqueID[kMaxmuons];   //[muons_]
   UInt_t          muons_momentum_fBits[kMaxmuons];   //[muons_]
   UInt_t          muons_momentum_fP_fUniqueID[kMaxmuons];   //[muons_]
   UInt_t          muons_momentum_fP_fBits[kMaxmuons];   //[muons_]
   Double_t        muons_momentum_fP_fX[kMaxmuons];   //[muons_]
   Double_t        muons_momentum_fP_fY[kMaxmuons];   //[muons_]
   Double_t        muons_momentum_fP_fZ[kMaxmuons];   //[muons_]
   Double_t        muons_momentum_fE[kMaxmuons];   //[muons_]
 //map<TString,UChar_t> muons_idPairs[kMaxmuons];
 //map<TString,susy::ElectronCollection> electrons;
 //map<TString,susy::PhotonCollection> photons;
 //map<TString,susy::CaloJetCollection> caloJets;
 //map<TString,susy::PFJetCollection> pfJets;
 //map<TString,susy::JPTJetCollection> jptJets;
 //map<TString,susy::PFParticleCollection> pfParticles;
   Int_t           generalTracks_;
   Int_t           generalTracks_algorithm[kMaxgeneralTracks];   //[generalTracks_]
   Int_t           generalTracks_quality[kMaxgeneralTracks];   //[generalTracks_]
   UChar_t         generalTracks_numberOfValidHits[kMaxgeneralTracks];   //[generalTracks_]
   UChar_t         generalTracks_numberOfValidTrackerHits[kMaxgeneralTracks];   //[generalTracks_]
   UChar_t         generalTracks_numberOfValidMuonHits[kMaxgeneralTracks];   //[generalTracks_]
   UChar_t         generalTracks_numberOfValidPixelHits[kMaxgeneralTracks];   //[generalTracks_]
   UChar_t         generalTracks_numberOfValidStripHits[kMaxgeneralTracks];   //[generalTracks_]
   Float_t         generalTracks_chi2[kMaxgeneralTracks];   //[generalTracks_]
   Float_t         generalTracks_ndof[kMaxgeneralTracks];   //[generalTracks_]
   Float_t         generalTracks_charge[kMaxgeneralTracks];   //[generalTracks_]
   Float_t         generalTracks_error[kMaxgeneralTracks][5];   //[generalTracks_]
   UInt_t          generalTracks_vertex_fUniqueID[kMaxgeneralTracks];   //[generalTracks_]
   UInt_t          generalTracks_vertex_fBits[kMaxgeneralTracks];   //[generalTracks_]
   Double_t        generalTracks_vertex_fX[kMaxgeneralTracks];   //[generalTracks_]
   Double_t        generalTracks_vertex_fY[kMaxgeneralTracks];   //[generalTracks_]
   Double_t        generalTracks_vertex_fZ[kMaxgeneralTracks];   //[generalTracks_]
   UInt_t          generalTracks_momentum_fUniqueID[kMaxgeneralTracks];   //[generalTracks_]
   UInt_t          generalTracks_momentum_fBits[kMaxgeneralTracks];   //[generalTracks_]
   UInt_t          generalTracks_momentum_fP_fUniqueID[kMaxgeneralTracks];   //[generalTracks_]
   UInt_t          generalTracks_momentum_fP_fBits[kMaxgeneralTracks];   //[generalTracks_]
   Double_t        generalTracks_momentum_fP_fX[kMaxgeneralTracks];   //[generalTracks_]
   Double_t        generalTracks_momentum_fP_fY[kMaxgeneralTracks];   //[generalTracks_]
   Double_t        generalTracks_momentum_fP_fZ[kMaxgeneralTracks];   //[generalTracks_]
   Double_t        generalTracks_momentum_fE[kMaxgeneralTracks];   //[generalTracks_]
 //map<TString,TVector3> generalTracks_extrapolatedPositions[kMaxgeneralTracks];
   Int_t           pu_;
   Int_t           pu_numInteractions[kMaxpu];   //[pu_]
   vector<float>   pu_zPositions[kMaxpu];
   vector<float>   pu_sumPTLowPT[kMaxpu];
   vector<float>   pu_sumPTHighPT[kMaxpu];
   vector<int>     pu_numTracksLowPT[kMaxpu];
   vector<int>     pu_numTracksHighPT[kMaxpu];
   vector<float>   pu_instLumi[kMaxpu];
   vector<unsigned int> pu_dataMixerRun[kMaxpu];
   vector<unsigned int> pu_dataMixerEvt[kMaxpu];
   vector<unsigned int> pu_dataMixerLumiSection[kMaxpu];
   Int_t           pu_BX[kMaxpu];   //[pu_]
   Float_t         pu_trueNumInteractions[kMaxpu];   //[pu_]
   Int_t           simVertices_;
   UInt_t          simVertices_fUniqueID[kMaxsimVertices];   //[simVertices_]
   UInt_t          simVertices_fBits[kMaxsimVertices];   //[simVertices_]
   Double_t        simVertices_fX[kMaxsimVertices];   //[simVertices_]
   Double_t        simVertices_fY[kMaxsimVertices];   //[simVertices_]
   Double_t        simVertices_fZ[kMaxsimVertices];   //[simVertices_]
   Int_t           genParticles_;
   UChar_t         genParticles_status[kMaxgenParticles];   //[genParticles_]
   Int_t           genParticles_motherId[kMaxgenParticles];   //[genParticles_]
   Int_t           genParticles_pdgId[kMaxgenParticles];   //[genParticles_]
   Char_t          genParticles_charge[kMaxgenParticles];   //[genParticles_]
   UInt_t          genParticles_vertex_fUniqueID[kMaxgenParticles];   //[genParticles_]
   UInt_t          genParticles_vertex_fBits[kMaxgenParticles];   //[genParticles_]
   Double_t        genParticles_vertex_fX[kMaxgenParticles];   //[genParticles_]
   Double_t        genParticles_vertex_fY[kMaxgenParticles];   //[genParticles_]
   Double_t        genParticles_vertex_fZ[kMaxgenParticles];   //[genParticles_]
   UInt_t          genParticles_momentum_fUniqueID[kMaxgenParticles];   //[genParticles_]
   UInt_t          genParticles_momentum_fBits[kMaxgenParticles];   //[genParticles_]
   UInt_t          genParticles_momentum_fP_fUniqueID[kMaxgenParticles];   //[genParticles_]
   UInt_t          genParticles_momentum_fP_fBits[kMaxgenParticles];   //[genParticles_]
   Double_t        genParticles_momentum_fP_fX[kMaxgenParticles];   //[genParticles_]
   Double_t        genParticles_momentum_fP_fY[kMaxgenParticles];   //[genParticles_]
   Double_t        genParticles_momentum_fP_fZ[kMaxgenParticles];   //[genParticles_]
   Double_t        genParticles_momentum_fE[kMaxgenParticles];   //[genParticles_]
   Int_t           gridParams_;
   TString         gridParams_first[kMaxgridParams];
   Float_t         gridParams_second[kMaxgridParams];   //[gridParams_]

   // List of branches
   TBranch        *b_susyEvent_isRealData;   //!
   TBranch        *b_susyEvent_runNumber;   //!
   TBranch        *b_susyEvent_eventNumber;   //!
   TBranch        *b_susyEvent_luminosityBlockNumber;   //!
   TBranch        *b_susyEvent_bunchCrossing;   //!
   TBranch        *b_susyEvent_avgInsRecLumi;   //!
   TBranch        *b_susyEvent_intgRecLumi;   //!
   TBranch        *b_susyEvent_cosmicFlag;   //!
   TBranch        *b_susyEvent_rho;   //!
   TBranch        *b_susyEvent_rhoBarrel;   //!
   TBranch        *b_susyEvent_rho25;   //!
   TBranch        *b_susyEvent_metFilterBit;   //!
   TBranch        *b_susyEvent_beamSpot_fUniqueID;   //!
   TBranch        *b_susyEvent_beamSpot_fBits;   //!
   TBranch        *b_susyEvent_beamSpot_fX;   //!
   TBranch        *b_susyEvent_beamSpot_fY;   //!
   TBranch        *b_susyEvent_beamSpot_fZ;   //!
   TBranch        *b_susyEvent_l1Map_;   //!
   TBranch        *b_l1Map_first;   //!
   TBranch        *b_l1Map_second_first;   //!
   TBranch        *b_l1Map_second_second;   //!
   TBranch        *b_susyEvent_hltMap_;   //!
   TBranch        *b_hltMap_first;   //!
   TBranch        *b_hltMap_second_first;   //!
   TBranch        *b_hltMap_second_second;   //!
   TBranch        *b_susyEvent_metMap_;   //!
   TBranch        *b_metMap_first;   //!
   TBranch        *b_metMap_second_sumEt;   //!
   TBranch        *b_metMap_second_significance;   //!
   TBranch        *b_metMap_second_mEt_fUniqueID;   //!
   TBranch        *b_metMap_second_mEt_fBits;   //!
   TBranch        *b_metMap_second_mEt_fX;   //!
   TBranch        *b_metMap_second_mEt_fY;   //!
   TBranch        *b_metMap_second_vertex_fUniqueID;   //!
   TBranch        *b_metMap_second_vertex_fBits;   //!
   TBranch        *b_metMap_second_vertex_fX;   //!
   TBranch        *b_metMap_second_vertex_fY;   //!
   TBranch        *b_metMap_second_vertex_fZ;   //!
   TBranch        *b_susyEvent_vertices_;   //!
   TBranch        *b_vertices_chi2;   //!
   TBranch        *b_vertices_ndof;   //!
   TBranch        *b_vertices_tracksSize;   //!
   TBranch        *b_vertices_position_fUniqueID;   //!
   TBranch        *b_vertices_position_fBits;   //!
   TBranch        *b_vertices_position_fX;   //!
   TBranch        *b_vertices_position_fY;   //!
   TBranch        *b_vertices_position_fZ;   //!
   TBranch        *b_susyEvent_tracks_;   //!
   TBranch        *b_tracks_algorithm;   //!
   TBranch        *b_tracks_quality;   //!
   TBranch        *b_tracks_numberOfValidHits;   //!
   TBranch        *b_tracks_numberOfValidTrackerHits;   //!
   TBranch        *b_tracks_numberOfValidMuonHits;   //!
   TBranch        *b_tracks_numberOfValidPixelHits;   //!
   TBranch        *b_tracks_numberOfValidStripHits;   //!
   TBranch        *b_tracks_chi2;   //!
   TBranch        *b_tracks_ndof;   //!
   TBranch        *b_tracks_charge;   //!
   TBranch        *b_tracks_error;   //!
   TBranch        *b_tracks_vertex_fUniqueID;   //!
   TBranch        *b_tracks_vertex_fBits;   //!
   TBranch        *b_tracks_vertex_fX;   //!
   TBranch        *b_tracks_vertex_fY;   //!
   TBranch        *b_tracks_vertex_fZ;   //!
   TBranch        *b_tracks_momentum_fUniqueID;   //!
   TBranch        *b_tracks_momentum_fBits;   //!
   TBranch        *b_tracks_momentum_fP_fUniqueID;   //!
   TBranch        *b_tracks_momentum_fP_fBits;   //!
   TBranch        *b_tracks_momentum_fP_fX;   //!
   TBranch        *b_tracks_momentum_fP_fY;   //!
   TBranch        *b_tracks_momentum_fP_fZ;   //!
   TBranch        *b_tracks_momentum_fE;   //!
   TBranch        *b_susyEvent_superClusters_;   //!
   TBranch        *b_superClusters_seedClusterIndex;   //!
   TBranch        *b_superClusters_energy;   //!
   TBranch        *b_superClusters_preshowerEnergy;   //!
   TBranch        *b_superClusters_phiWidth;   //!
   TBranch        *b_superClusters_etaWidth;   //!
   TBranch        *b_superClusters_position_fUniqueID;   //!
   TBranch        *b_superClusters_position_fBits;   //!
   TBranch        *b_superClusters_position_fX;   //!
   TBranch        *b_superClusters_position_fY;   //!
   TBranch        *b_superClusters_position_fZ;   //!
   TBranch        *b_superClusters_basicClusterIndices;   //!
   TBranch        *b_susyEvent_clusters_;   //!
   TBranch        *b_clusters_nCrystals;   //!
   TBranch        *b_clusters_energy;   //!
   TBranch        *b_clusters_position_fUniqueID;   //!
   TBranch        *b_clusters_position_fBits;   //!
   TBranch        *b_clusters_position_fX;   //!
   TBranch        *b_clusters_position_fY;   //!
   TBranch        *b_clusters_position_fZ;   //!
   TBranch        *b_susyEvent_muons_;   //!
   TBranch        *b_muons_type;   //!
   TBranch        *b_muons_nMatches;   //!
   TBranch        *b_muons_nValidHits;   //!
   TBranch        *b_muons_nValidTrackerHits;   //!
   TBranch        *b_muons_nValidMuonHits;   //!
   TBranch        *b_muons_nChambers;   //!
   TBranch        *b_muons_timeNDof;   //!
   TBranch        *b_muons_timeDirection;   //!
   TBranch        *b_muons_timeAtIp;   //!
   TBranch        *b_muons_timeAtIpError;   //!
   TBranch        *b_muons_caloCompatibility;   //!
   TBranch        *b_muons_emEnergy;   //!
   TBranch        *b_muons_hadEnergy;   //!
   TBranch        *b_muons_trackIsoR03;   //!
   TBranch        *b_muons_ecalIsoR03;   //!
   TBranch        *b_muons_hcalIsoR03;   //!
   TBranch        *b_muons_trackIsoR05;   //!
   TBranch        *b_muons_ecalIsoR05;   //!
   TBranch        *b_muons_hcalIsoR05;   //!
   TBranch        *b_muons_sumChargedHadronPt03;   //!
   TBranch        *b_muons_sumChargedParticlePt03;   //!
   TBranch        *b_muons_sumNeutralHadronEt03;   //!
   TBranch        *b_muons_sumPhotonEt03;   //!
   TBranch        *b_muons_sumNeutralHadronEtHighThreshold03;   //!
   TBranch        *b_muons_sumPhotonEtHighThreshold03;   //!
   TBranch        *b_muons_sumPUPt03;   //!
   TBranch        *b_muons_sumChargedHadronPt04;   //!
   TBranch        *b_muons_sumChargedParticlePt04;   //!
   TBranch        *b_muons_sumNeutralHadronEt04;   //!
   TBranch        *b_muons_sumPhotonEt04;   //!
   TBranch        *b_muons_sumNeutralHadronEtHighThreshold04;   //!
   TBranch        *b_muons_sumPhotonEtHighThreshold04;   //!
   TBranch        *b_muons_sumPUPt04;   //!
   TBranch        *b_muons_trackIndex;   //!
   TBranch        *b_muons_standAloneTrackIndex;   //!
   TBranch        *b_muons_combinedTrackIndex;   //!
   TBranch        *b_muons_momentum_fUniqueID;   //!
   TBranch        *b_muons_momentum_fBits;   //!
   TBranch        *b_muons_momentum_fP_fUniqueID;   //!
   TBranch        *b_muons_momentum_fP_fBits;   //!
   TBranch        *b_muons_momentum_fP_fX;   //!
   TBranch        *b_muons_momentum_fP_fY;   //!
   TBranch        *b_muons_momentum_fP_fZ;   //!
   TBranch        *b_muons_momentum_fE;   //!
   TBranch        *b_susyEvent_generalTracks_;   //!
   TBranch        *b_generalTracks_algorithm;   //!
   TBranch        *b_generalTracks_quality;   //!
   TBranch        *b_generalTracks_numberOfValidHits;   //!
   TBranch        *b_generalTracks_numberOfValidTrackerHits;   //!
   TBranch        *b_generalTracks_numberOfValidMuonHits;   //!
   TBranch        *b_generalTracks_numberOfValidPixelHits;   //!
   TBranch        *b_generalTracks_numberOfValidStripHits;   //!
   TBranch        *b_generalTracks_chi2;   //!
   TBranch        *b_generalTracks_ndof;   //!
   TBranch        *b_generalTracks_charge;   //!
   TBranch        *b_generalTracks_error;   //!
   TBranch        *b_generalTracks_vertex_fUniqueID;   //!
   TBranch        *b_generalTracks_vertex_fBits;   //!
   TBranch        *b_generalTracks_vertex_fX;   //!
   TBranch        *b_generalTracks_vertex_fY;   //!
   TBranch        *b_generalTracks_vertex_fZ;   //!
   TBranch        *b_generalTracks_momentum_fUniqueID;   //!
   TBranch        *b_generalTracks_momentum_fBits;   //!
   TBranch        *b_generalTracks_momentum_fP_fUniqueID;   //!
   TBranch        *b_generalTracks_momentum_fP_fBits;   //!
   TBranch        *b_generalTracks_momentum_fP_fX;   //!
   TBranch        *b_generalTracks_momentum_fP_fY;   //!
   TBranch        *b_generalTracks_momentum_fP_fZ;   //!
   TBranch        *b_generalTracks_momentum_fE;   //!
   TBranch        *b_susyEvent_pu_;   //!
   TBranch        *b_pu_numInteractions;   //!
   TBranch        *b_pu_zPositions;   //!
   TBranch        *b_pu_sumPTLowPT;   //!
   TBranch        *b_pu_sumPTHighPT;   //!
   TBranch        *b_pu_numTracksLowPT;   //!
   TBranch        *b_pu_numTracksHighPT;   //!
   TBranch        *b_pu_instLumi;   //!
   TBranch        *b_pu_dataMixerRun;   //!
   TBranch        *b_pu_dataMixerEvt;   //!
   TBranch        *b_pu_dataMixerLumiSection;   //!
   TBranch        *b_pu_BX;   //!
   TBranch        *b_pu_trueNumInteractions;   //!
   TBranch        *b_susyEvent_simVertices_;   //!
   TBranch        *b_simVertices_fUniqueID;   //!
   TBranch        *b_simVertices_fBits;   //!
   TBranch        *b_simVertices_fX;   //!
   TBranch        *b_simVertices_fY;   //!
   TBranch        *b_simVertices_fZ;   //!
   TBranch        *b_susyEvent_genParticles_;   //!
   TBranch        *b_genParticles_status;   //!
   TBranch        *b_genParticles_motherId;   //!
   TBranch        *b_genParticles_pdgId;   //!
   TBranch        *b_genParticles_charge;   //!
   TBranch        *b_genParticles_vertex_fUniqueID;   //!
   TBranch        *b_genParticles_vertex_fBits;   //!
   TBranch        *b_genParticles_vertex_fX;   //!
   TBranch        *b_genParticles_vertex_fY;   //!
   TBranch        *b_genParticles_vertex_fZ;   //!
   TBranch        *b_genParticles_momentum_fUniqueID;   //!
   TBranch        *b_genParticles_momentum_fBits;   //!
   TBranch        *b_genParticles_momentum_fP_fUniqueID;   //!
   TBranch        *b_genParticles_momentum_fP_fBits;   //!
   TBranch        *b_genParticles_momentum_fP_fX;   //!
   TBranch        *b_genParticles_momentum_fP_fY;   //!
   TBranch        *b_genParticles_momentum_fP_fZ;   //!
   TBranch        *b_genParticles_momentum_fE;   //!
   TBranch        *b_susyEvent_gridParams_;   //!
   TBranch        *b_gridParams_first;   //!
   TBranch        *b_gridParams_second;   //!

   paperino(TTree *tree=0);
   virtual ~paperino();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef paperino_cxx
paperino::paperino(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/tmp/azzolini/susyEvents.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/tmp/azzolini/susyEvents.root");
      }
      f->GetObject("susyTree",tree);

   }
   Init(tree);
}

paperino::~paperino()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t paperino::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t paperino::LoadTree(Long64_t entry)
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

void paperino::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("isRealData", &isRealData, &b_susyEvent_isRealData);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_susyEvent_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_susyEvent_eventNumber);
   fChain->SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber, &b_susyEvent_luminosityBlockNumber);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_susyEvent_bunchCrossing);
   fChain->SetBranchAddress("avgInsRecLumi", &avgInsRecLumi, &b_susyEvent_avgInsRecLumi);
   fChain->SetBranchAddress("intgRecLumi", &intgRecLumi, &b_susyEvent_intgRecLumi);
   fChain->SetBranchAddress("cosmicFlag", &cosmicFlag, &b_susyEvent_cosmicFlag);
   fChain->SetBranchAddress("rho", &rho, &b_susyEvent_rho);
   fChain->SetBranchAddress("rhoBarrel", &rhoBarrel, &b_susyEvent_rhoBarrel);
   fChain->SetBranchAddress("rho25", &rho25, &b_susyEvent_rho25);
   fChain->SetBranchAddress("metFilterBit", &metFilterBit, &b_susyEvent_metFilterBit);
   fChain->SetBranchAddress("beamSpot.fUniqueID", &beamSpot_fUniqueID, &b_susyEvent_beamSpot_fUniqueID);
   fChain->SetBranchAddress("beamSpot.fBits", &beamSpot_fBits, &b_susyEvent_beamSpot_fBits);
   fChain->SetBranchAddress("beamSpot.fX", &beamSpot_fX, &b_susyEvent_beamSpot_fX);
   fChain->SetBranchAddress("beamSpot.fY", &beamSpot_fY, &b_susyEvent_beamSpot_fY);
   fChain->SetBranchAddress("beamSpot.fZ", &beamSpot_fZ, &b_susyEvent_beamSpot_fZ);
   fChain->SetBranchAddress("l1Map", &l1Map_, &b_susyEvent_l1Map_);
   fChain->SetBranchAddress("l1Map.first", l1Map_first, &b_l1Map_first);
   fChain->SetBranchAddress("l1Map.second.first", l1Map_second_first, &b_l1Map_second_first);
   fChain->SetBranchAddress("l1Map.second.second", l1Map_second_second, &b_l1Map_second_second);
   fChain->SetBranchAddress("hltMap", &hltMap_, &b_susyEvent_hltMap_);
   fChain->SetBranchAddress("hltMap.first", hltMap_first, &b_hltMap_first);
   fChain->SetBranchAddress("hltMap.second.first", hltMap_second_first, &b_hltMap_second_first);
   fChain->SetBranchAddress("hltMap.second.second", hltMap_second_second, &b_hltMap_second_second);
   fChain->SetBranchAddress("metMap", &metMap_, &b_susyEvent_metMap_);
   fChain->SetBranchAddress("metMap.first", metMap_first, &b_metMap_first);
   fChain->SetBranchAddress("metMap.second.sumEt", metMap_second_sumEt, &b_metMap_second_sumEt);
   fChain->SetBranchAddress("metMap.second.significance", metMap_second_significance, &b_metMap_second_significance);
   fChain->SetBranchAddress("metMap.second.mEt.fUniqueID", metMap_second_mEt_fUniqueID, &b_metMap_second_mEt_fUniqueID);
   fChain->SetBranchAddress("metMap.second.mEt.fBits", metMap_second_mEt_fBits, &b_metMap_second_mEt_fBits);
   fChain->SetBranchAddress("metMap.second.mEt.fX", metMap_second_mEt_fX, &b_metMap_second_mEt_fX);
   fChain->SetBranchAddress("metMap.second.mEt.fY", metMap_second_mEt_fY, &b_metMap_second_mEt_fY);
   fChain->SetBranchAddress("metMap.second.vertex.fUniqueID", metMap_second_vertex_fUniqueID, &b_metMap_second_vertex_fUniqueID);
   fChain->SetBranchAddress("metMap.second.vertex.fBits", metMap_second_vertex_fBits, &b_metMap_second_vertex_fBits);
   fChain->SetBranchAddress("metMap.second.vertex.fX", metMap_second_vertex_fX, &b_metMap_second_vertex_fX);
   fChain->SetBranchAddress("metMap.second.vertex.fY", metMap_second_vertex_fY, &b_metMap_second_vertex_fY);
   fChain->SetBranchAddress("metMap.second.vertex.fZ", metMap_second_vertex_fZ, &b_metMap_second_vertex_fZ);
   fChain->SetBranchAddress("vertices", &vertices_, &b_susyEvent_vertices_);
   fChain->SetBranchAddress("vertices.chi2", vertices_chi2, &b_vertices_chi2);
   fChain->SetBranchAddress("vertices.ndof", vertices_ndof, &b_vertices_ndof);
   fChain->SetBranchAddress("vertices.tracksSize", vertices_tracksSize, &b_vertices_tracksSize);
   fChain->SetBranchAddress("vertices.position.fUniqueID", vertices_position_fUniqueID, &b_vertices_position_fUniqueID);
   fChain->SetBranchAddress("vertices.position.fBits", vertices_position_fBits, &b_vertices_position_fBits);
   fChain->SetBranchAddress("vertices.position.fX", vertices_position_fX, &b_vertices_position_fX);
   fChain->SetBranchAddress("vertices.position.fY", vertices_position_fY, &b_vertices_position_fY);
   fChain->SetBranchAddress("vertices.position.fZ", vertices_position_fZ, &b_vertices_position_fZ);
   fChain->SetBranchAddress("tracks", &tracks_, &b_susyEvent_tracks_);
   fChain->SetBranchAddress("tracks.algorithm", tracks_algorithm, &b_tracks_algorithm);
   fChain->SetBranchAddress("tracks.quality", tracks_quality, &b_tracks_quality);
   fChain->SetBranchAddress("tracks.numberOfValidHits", tracks_numberOfValidHits, &b_tracks_numberOfValidHits);
   fChain->SetBranchAddress("tracks.numberOfValidTrackerHits", tracks_numberOfValidTrackerHits, &b_tracks_numberOfValidTrackerHits);
   fChain->SetBranchAddress("tracks.numberOfValidMuonHits", tracks_numberOfValidMuonHits, &b_tracks_numberOfValidMuonHits);
   fChain->SetBranchAddress("tracks.numberOfValidPixelHits", tracks_numberOfValidPixelHits, &b_tracks_numberOfValidPixelHits);
   fChain->SetBranchAddress("tracks.numberOfValidStripHits", tracks_numberOfValidStripHits, &b_tracks_numberOfValidStripHits);
   fChain->SetBranchAddress("tracks.chi2", tracks_chi2, &b_tracks_chi2);
   fChain->SetBranchAddress("tracks.ndof", tracks_ndof, &b_tracks_ndof);
   fChain->SetBranchAddress("tracks.charge", tracks_charge, &b_tracks_charge);
   fChain->SetBranchAddress("tracks.error[5]", tracks_error, &b_tracks_error);
   fChain->SetBranchAddress("tracks.vertex.fUniqueID", tracks_vertex_fUniqueID, &b_tracks_vertex_fUniqueID);
   fChain->SetBranchAddress("tracks.vertex.fBits", tracks_vertex_fBits, &b_tracks_vertex_fBits);
   fChain->SetBranchAddress("tracks.vertex.fX", tracks_vertex_fX, &b_tracks_vertex_fX);
   fChain->SetBranchAddress("tracks.vertex.fY", tracks_vertex_fY, &b_tracks_vertex_fY);
   fChain->SetBranchAddress("tracks.vertex.fZ", tracks_vertex_fZ, &b_tracks_vertex_fZ);
   fChain->SetBranchAddress("tracks.momentum.fUniqueID", tracks_momentum_fUniqueID, &b_tracks_momentum_fUniqueID);
   fChain->SetBranchAddress("tracks.momentum.fBits", tracks_momentum_fBits, &b_tracks_momentum_fBits);
   fChain->SetBranchAddress("tracks.momentum.fP.fUniqueID", tracks_momentum_fP_fUniqueID, &b_tracks_momentum_fP_fUniqueID);
   fChain->SetBranchAddress("tracks.momentum.fP.fBits", tracks_momentum_fP_fBits, &b_tracks_momentum_fP_fBits);
   fChain->SetBranchAddress("tracks.momentum.fP.fX", tracks_momentum_fP_fX, &b_tracks_momentum_fP_fX);
   fChain->SetBranchAddress("tracks.momentum.fP.fY", tracks_momentum_fP_fY, &b_tracks_momentum_fP_fY);
   fChain->SetBranchAddress("tracks.momentum.fP.fZ", tracks_momentum_fP_fZ, &b_tracks_momentum_fP_fZ);
   fChain->SetBranchAddress("tracks.momentum.fE", tracks_momentum_fE, &b_tracks_momentum_fE);
   fChain->SetBranchAddress("superClusters", &superClusters_, &b_susyEvent_superClusters_);
   fChain->SetBranchAddress("superClusters.seedClusterIndex", superClusters_seedClusterIndex, &b_superClusters_seedClusterIndex);
   fChain->SetBranchAddress("superClusters.energy", superClusters_energy, &b_superClusters_energy);
   fChain->SetBranchAddress("superClusters.preshowerEnergy", superClusters_preshowerEnergy, &b_superClusters_preshowerEnergy);
   fChain->SetBranchAddress("superClusters.phiWidth", superClusters_phiWidth, &b_superClusters_phiWidth);
   fChain->SetBranchAddress("superClusters.etaWidth", superClusters_etaWidth, &b_superClusters_etaWidth);
   fChain->SetBranchAddress("superClusters.position.fUniqueID", superClusters_position_fUniqueID, &b_superClusters_position_fUniqueID);
   fChain->SetBranchAddress("superClusters.position.fBits", superClusters_position_fBits, &b_superClusters_position_fBits);
   fChain->SetBranchAddress("superClusters.position.fX", superClusters_position_fX, &b_superClusters_position_fX);
   fChain->SetBranchAddress("superClusters.position.fY", superClusters_position_fY, &b_superClusters_position_fY);
   fChain->SetBranchAddress("superClusters.position.fZ", superClusters_position_fZ, &b_superClusters_position_fZ);
   fChain->SetBranchAddress("superClusters.basicClusterIndices", superClusters_basicClusterIndices, &b_superClusters_basicClusterIndices);
   fChain->SetBranchAddress("clusters", &clusters_, &b_susyEvent_clusters_);
   fChain->SetBranchAddress("clusters.nCrystals", clusters_nCrystals, &b_clusters_nCrystals);
   fChain->SetBranchAddress("clusters.energy", clusters_energy, &b_clusters_energy);
   fChain->SetBranchAddress("clusters.position.fUniqueID", clusters_position_fUniqueID, &b_clusters_position_fUniqueID);
   fChain->SetBranchAddress("clusters.position.fBits", clusters_position_fBits, &b_clusters_position_fBits);
   fChain->SetBranchAddress("clusters.position.fX", clusters_position_fX, &b_clusters_position_fX);
   fChain->SetBranchAddress("clusters.position.fY", clusters_position_fY, &b_clusters_position_fY);
   fChain->SetBranchAddress("clusters.position.fZ", clusters_position_fZ, &b_clusters_position_fZ);
   fChain->SetBranchAddress("muons", &muons_, &b_susyEvent_muons_);
   fChain->SetBranchAddress("muons.type", muons_type, &b_muons_type);
   fChain->SetBranchAddress("muons.nMatches", muons_nMatches, &b_muons_nMatches);
   fChain->SetBranchAddress("muons.nValidHits", muons_nValidHits, &b_muons_nValidHits);
   fChain->SetBranchAddress("muons.nValidTrackerHits", muons_nValidTrackerHits, &b_muons_nValidTrackerHits);
   fChain->SetBranchAddress("muons.nValidMuonHits", muons_nValidMuonHits, &b_muons_nValidMuonHits);
   fChain->SetBranchAddress("muons.nChambers", muons_nChambers, &b_muons_nChambers);
   fChain->SetBranchAddress("muons.timeNDof", muons_timeNDof, &b_muons_timeNDof);
   fChain->SetBranchAddress("muons.timeDirection", muons_timeDirection, &b_muons_timeDirection);
   fChain->SetBranchAddress("muons.timeAtIp", muons_timeAtIp, &b_muons_timeAtIp);
   fChain->SetBranchAddress("muons.timeAtIpError", muons_timeAtIpError, &b_muons_timeAtIpError);
   fChain->SetBranchAddress("muons.caloCompatibility", muons_caloCompatibility, &b_muons_caloCompatibility);
   fChain->SetBranchAddress("muons.emEnergy", muons_emEnergy, &b_muons_emEnergy);
   fChain->SetBranchAddress("muons.hadEnergy", muons_hadEnergy, &b_muons_hadEnergy);
   fChain->SetBranchAddress("muons.trackIsoR03", muons_trackIsoR03, &b_muons_trackIsoR03);
   fChain->SetBranchAddress("muons.ecalIsoR03", muons_ecalIsoR03, &b_muons_ecalIsoR03);
   fChain->SetBranchAddress("muons.hcalIsoR03", muons_hcalIsoR03, &b_muons_hcalIsoR03);
   fChain->SetBranchAddress("muons.trackIsoR05", muons_trackIsoR05, &b_muons_trackIsoR05);
   fChain->SetBranchAddress("muons.ecalIsoR05", muons_ecalIsoR05, &b_muons_ecalIsoR05);
   fChain->SetBranchAddress("muons.hcalIsoR05", muons_hcalIsoR05, &b_muons_hcalIsoR05);
   fChain->SetBranchAddress("muons.sumChargedHadronPt03", muons_sumChargedHadronPt03, &b_muons_sumChargedHadronPt03);
   fChain->SetBranchAddress("muons.sumChargedParticlePt03", muons_sumChargedParticlePt03, &b_muons_sumChargedParticlePt03);
   fChain->SetBranchAddress("muons.sumNeutralHadronEt03", muons_sumNeutralHadronEt03, &b_muons_sumNeutralHadronEt03);
   fChain->SetBranchAddress("muons.sumPhotonEt03", muons_sumPhotonEt03, &b_muons_sumPhotonEt03);
   fChain->SetBranchAddress("muons.sumNeutralHadronEtHighThreshold03", muons_sumNeutralHadronEtHighThreshold03, &b_muons_sumNeutralHadronEtHighThreshold03);
   fChain->SetBranchAddress("muons.sumPhotonEtHighThreshold03", muons_sumPhotonEtHighThreshold03, &b_muons_sumPhotonEtHighThreshold03);
   fChain->SetBranchAddress("muons.sumPUPt03", muons_sumPUPt03, &b_muons_sumPUPt03);
   fChain->SetBranchAddress("muons.sumChargedHadronPt04", muons_sumChargedHadronPt04, &b_muons_sumChargedHadronPt04);
   fChain->SetBranchAddress("muons.sumChargedParticlePt04", muons_sumChargedParticlePt04, &b_muons_sumChargedParticlePt04);
   fChain->SetBranchAddress("muons.sumNeutralHadronEt04", muons_sumNeutralHadronEt04, &b_muons_sumNeutralHadronEt04);
   fChain->SetBranchAddress("muons.sumPhotonEt04", muons_sumPhotonEt04, &b_muons_sumPhotonEt04);
   fChain->SetBranchAddress("muons.sumNeutralHadronEtHighThreshold04", muons_sumNeutralHadronEtHighThreshold04, &b_muons_sumNeutralHadronEtHighThreshold04);
   fChain->SetBranchAddress("muons.sumPhotonEtHighThreshold04", muons_sumPhotonEtHighThreshold04, &b_muons_sumPhotonEtHighThreshold04);
   fChain->SetBranchAddress("muons.sumPUPt04", muons_sumPUPt04, &b_muons_sumPUPt04);
   fChain->SetBranchAddress("muons.trackIndex", muons_trackIndex, &b_muons_trackIndex);
   fChain->SetBranchAddress("muons.standAloneTrackIndex", muons_standAloneTrackIndex, &b_muons_standAloneTrackIndex);
   fChain->SetBranchAddress("muons.combinedTrackIndex", muons_combinedTrackIndex, &b_muons_combinedTrackIndex);
   fChain->SetBranchAddress("muons.momentum.fUniqueID", muons_momentum_fUniqueID, &b_muons_momentum_fUniqueID);
   fChain->SetBranchAddress("muons.momentum.fBits", muons_momentum_fBits, &b_muons_momentum_fBits);
   fChain->SetBranchAddress("muons.momentum.fP.fUniqueID", muons_momentum_fP_fUniqueID, &b_muons_momentum_fP_fUniqueID);
   fChain->SetBranchAddress("muons.momentum.fP.fBits", muons_momentum_fP_fBits, &b_muons_momentum_fP_fBits);
   fChain->SetBranchAddress("muons.momentum.fP.fX", muons_momentum_fP_fX, &b_muons_momentum_fP_fX);
   fChain->SetBranchAddress("muons.momentum.fP.fY", muons_momentum_fP_fY, &b_muons_momentum_fP_fY);
   fChain->SetBranchAddress("muons.momentum.fP.fZ", muons_momentum_fP_fZ, &b_muons_momentum_fP_fZ);
   fChain->SetBranchAddress("muons.momentum.fE", muons_momentum_fE, &b_muons_momentum_fE);
   fChain->SetBranchAddress("generalTracks", &generalTracks_, &b_susyEvent_generalTracks_);
   fChain->SetBranchAddress("generalTracks.algorithm", &generalTracks_algorithm, &b_generalTracks_algorithm);
   fChain->SetBranchAddress("generalTracks.quality", &generalTracks_quality, &b_generalTracks_quality);
   fChain->SetBranchAddress("generalTracks.numberOfValidHits", &generalTracks_numberOfValidHits, &b_generalTracks_numberOfValidHits);
   fChain->SetBranchAddress("generalTracks.numberOfValidTrackerHits", &generalTracks_numberOfValidTrackerHits, &b_generalTracks_numberOfValidTrackerHits);
   fChain->SetBranchAddress("generalTracks.numberOfValidMuonHits", &generalTracks_numberOfValidMuonHits, &b_generalTracks_numberOfValidMuonHits);
   fChain->SetBranchAddress("generalTracks.numberOfValidPixelHits", &generalTracks_numberOfValidPixelHits, &b_generalTracks_numberOfValidPixelHits);
   fChain->SetBranchAddress("generalTracks.numberOfValidStripHits", &generalTracks_numberOfValidStripHits, &b_generalTracks_numberOfValidStripHits);
   fChain->SetBranchAddress("generalTracks.chi2", &generalTracks_chi2, &b_generalTracks_chi2);
   fChain->SetBranchAddress("generalTracks.ndof", &generalTracks_ndof, &b_generalTracks_ndof);
   fChain->SetBranchAddress("generalTracks.charge", &generalTracks_charge, &b_generalTracks_charge);
   fChain->SetBranchAddress("generalTracks.error[5]", &generalTracks_error, &b_generalTracks_error);
   fChain->SetBranchAddress("generalTracks.vertex.fUniqueID", &generalTracks_vertex_fUniqueID, &b_generalTracks_vertex_fUniqueID);
   fChain->SetBranchAddress("generalTracks.vertex.fBits", &generalTracks_vertex_fBits, &b_generalTracks_vertex_fBits);
   fChain->SetBranchAddress("generalTracks.vertex.fX", &generalTracks_vertex_fX, &b_generalTracks_vertex_fX);
   fChain->SetBranchAddress("generalTracks.vertex.fY", &generalTracks_vertex_fY, &b_generalTracks_vertex_fY);
   fChain->SetBranchAddress("generalTracks.vertex.fZ", &generalTracks_vertex_fZ, &b_generalTracks_vertex_fZ);
   fChain->SetBranchAddress("generalTracks.momentum.fUniqueID", &generalTracks_momentum_fUniqueID, &b_generalTracks_momentum_fUniqueID);
   fChain->SetBranchAddress("generalTracks.momentum.fBits", &generalTracks_momentum_fBits, &b_generalTracks_momentum_fBits);
   fChain->SetBranchAddress("generalTracks.momentum.fP.fUniqueID", &generalTracks_momentum_fP_fUniqueID, &b_generalTracks_momentum_fP_fUniqueID);
   fChain->SetBranchAddress("generalTracks.momentum.fP.fBits", &generalTracks_momentum_fP_fBits, &b_generalTracks_momentum_fP_fBits);
   fChain->SetBranchAddress("generalTracks.momentum.fP.fX", &generalTracks_momentum_fP_fX, &b_generalTracks_momentum_fP_fX);
   fChain->SetBranchAddress("generalTracks.momentum.fP.fY", &generalTracks_momentum_fP_fY, &b_generalTracks_momentum_fP_fY);
   fChain->SetBranchAddress("generalTracks.momentum.fP.fZ", &generalTracks_momentum_fP_fZ, &b_generalTracks_momentum_fP_fZ);
   fChain->SetBranchAddress("generalTracks.momentum.fE", &generalTracks_momentum_fE, &b_generalTracks_momentum_fE);
   fChain->SetBranchAddress("pu", &pu_, &b_susyEvent_pu_);
   fChain->SetBranchAddress("pu.numInteractions", &pu_numInteractions, &b_pu_numInteractions);
   fChain->SetBranchAddress("pu.zPositions", &pu_zPositions, &b_pu_zPositions);
   fChain->SetBranchAddress("pu.sumPTLowPT", &pu_sumPTLowPT, &b_pu_sumPTLowPT);
   fChain->SetBranchAddress("pu.sumPTHighPT", &pu_sumPTHighPT, &b_pu_sumPTHighPT);
   fChain->SetBranchAddress("pu.numTracksLowPT", &pu_numTracksLowPT, &b_pu_numTracksLowPT);
   fChain->SetBranchAddress("pu.numTracksHighPT", &pu_numTracksHighPT, &b_pu_numTracksHighPT);
   fChain->SetBranchAddress("pu.instLumi", &pu_instLumi, &b_pu_instLumi);
   fChain->SetBranchAddress("pu.dataMixerRun", &pu_dataMixerRun, &b_pu_dataMixerRun);
   fChain->SetBranchAddress("pu.dataMixerEvt", &pu_dataMixerEvt, &b_pu_dataMixerEvt);
   fChain->SetBranchAddress("pu.dataMixerLumiSection", &pu_dataMixerLumiSection, &b_pu_dataMixerLumiSection);
   fChain->SetBranchAddress("pu.BX", &pu_BX, &b_pu_BX);
   fChain->SetBranchAddress("pu.trueNumInteractions", &pu_trueNumInteractions, &b_pu_trueNumInteractions);
   fChain->SetBranchAddress("simVertices", &simVertices_, &b_susyEvent_simVertices_);
   fChain->SetBranchAddress("simVertices.fUniqueID", &simVertices_fUniqueID, &b_simVertices_fUniqueID);
   fChain->SetBranchAddress("simVertices.fBits", &simVertices_fBits, &b_simVertices_fBits);
   fChain->SetBranchAddress("simVertices.fX", &simVertices_fX, &b_simVertices_fX);
   fChain->SetBranchAddress("simVertices.fY", &simVertices_fY, &b_simVertices_fY);
   fChain->SetBranchAddress("simVertices.fZ", &simVertices_fZ, &b_simVertices_fZ);
   fChain->SetBranchAddress("genParticles", &genParticles_, &b_susyEvent_genParticles_);
   fChain->SetBranchAddress("genParticles.status", &genParticles_status, &b_genParticles_status);
   fChain->SetBranchAddress("genParticles.motherId", &genParticles_motherId, &b_genParticles_motherId);
   fChain->SetBranchAddress("genParticles.pdgId", &genParticles_pdgId, &b_genParticles_pdgId);
   fChain->SetBranchAddress("genParticles.charge", &genParticles_charge, &b_genParticles_charge);
   fChain->SetBranchAddress("genParticles.vertex.fUniqueID", &genParticles_vertex_fUniqueID, &b_genParticles_vertex_fUniqueID);
   fChain->SetBranchAddress("genParticles.vertex.fBits", &genParticles_vertex_fBits, &b_genParticles_vertex_fBits);
   fChain->SetBranchAddress("genParticles.vertex.fX", &genParticles_vertex_fX, &b_genParticles_vertex_fX);
   fChain->SetBranchAddress("genParticles.vertex.fY", &genParticles_vertex_fY, &b_genParticles_vertex_fY);
   fChain->SetBranchAddress("genParticles.vertex.fZ", &genParticles_vertex_fZ, &b_genParticles_vertex_fZ);
   fChain->SetBranchAddress("genParticles.momentum.fUniqueID", &genParticles_momentum_fUniqueID, &b_genParticles_momentum_fUniqueID);
   fChain->SetBranchAddress("genParticles.momentum.fBits", &genParticles_momentum_fBits, &b_genParticles_momentum_fBits);
   fChain->SetBranchAddress("genParticles.momentum.fP.fUniqueID", &genParticles_momentum_fP_fUniqueID, &b_genParticles_momentum_fP_fUniqueID);
   fChain->SetBranchAddress("genParticles.momentum.fP.fBits", &genParticles_momentum_fP_fBits, &b_genParticles_momentum_fP_fBits);
   fChain->SetBranchAddress("genParticles.momentum.fP.fX", &genParticles_momentum_fP_fX, &b_genParticles_momentum_fP_fX);
   fChain->SetBranchAddress("genParticles.momentum.fP.fY", &genParticles_momentum_fP_fY, &b_genParticles_momentum_fP_fY);
   fChain->SetBranchAddress("genParticles.momentum.fP.fZ", &genParticles_momentum_fP_fZ, &b_genParticles_momentum_fP_fZ);
   fChain->SetBranchAddress("genParticles.momentum.fE", &genParticles_momentum_fE, &b_genParticles_momentum_fE);
   fChain->SetBranchAddress("gridParams", &gridParams_, &b_susyEvent_gridParams_);
   fChain->SetBranchAddress("gridParams.first", &gridParams_first, &b_gridParams_first);
   fChain->SetBranchAddress("gridParams.second", &gridParams_second, &b_gridParams_second);
   Notify();
}

Bool_t paperino::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void paperino::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t paperino::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef paperino_cxx
