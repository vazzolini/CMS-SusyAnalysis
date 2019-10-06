// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.cc
// 
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.cc,v 1.12 2012/05/03 19:58:51 dwjang Exp $
//

#define SusyEventAnalyzer_cxx

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>

#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <utility>

#include "SusyEventAnalyzer.h"
#include "SusyEventPrinter.h"

#include "../jec/JetMETObjects/interface/JetCorrectorParameters.h"
#include "../jec/JetMETObjects/interface/FactorizedJetCorrector.h"


template<typename T> bool EtGreater(const T* p1, const T* p2) {
  return (p1->momentum.Et() > p2->momentum.Et());
}


void SusyEventAnalyzer::InitializePerEvent() {

}


bool SusyEventAnalyzer::isSameObject(TLorentzVector& p1, TLorentzVector& p2) {

  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < 0.5) return true;
  return false;
}


float SusyEventAnalyzer::d0correction(TVector3& beamSpot, susy::Track& track) const {

  float d0 = track.d0() - beamSpot.X()*std::sin(track.phi()) + beamSpot.Y()*std::cos(track.phi());
  return d0;
}


bool SusyEventAnalyzer::PassTrigger(TString path) {
  bool pass = false;
  for(susy::TriggerMap::iterator it = event->hltMap.begin(); it != event->hltMap.end(); it++) {
    if(it->first.Contains(path) && (int(it->second.second)) ) {
      pass = true;
      break;
    }
  }
  return pass;
}


bool SusyEventAnalyzer::PassTriggers() {
  bool pass = false;
  for(std::vector<TString>::iterator it = hltNames.begin(); it != hltNames.end(); it++) {
    if(PassTrigger(*it)) {
      pass = true;
      std::cout << "VIRpasstrigger" << pass << std::endl;
      break;
    }
  } 
  return pass;
}



void SusyEventAnalyzer::Loop() {

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  std::cout << "total events in files  : " << nentries << std::endl;

  if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;

  std::cout << "events to be processed : " << processNEvents << std::endl; 
  
  
  if(printLevel > 0) std::cout << "Initialize event counters." << std::endl;
  const int NCNT = 20;
  int nCnt[NCNT];
  for(int i=0; i<NCNT; i++) nCnt[i] = 0;
  
  //// --- VIR -----
  const int NVIR = 20;
  int nVir[NVIR];
  for(int i=0; i<NVIR; i++) nVir[i] = 0;
  
  // event counter
  nCnt[0]++; // total number of events
  nVir[0]++; // total number of events
  
  std::cout << " ----------------- Job Summary prova ----------------- " << std::endl;
  std::cout << " Total events            : " << nCnt[0] << std::endl;
  std::cout << " Total events            : " << nCnt[0] << std::endl;
  std::cout << " Total events            : " << nCnt[0] << " Total events vir        : " << nVir[0] <<std::endl;
  
  
  //// --- VIR uncomment for default -----
  if(printLevel > 0) std::cout <<" Apply trigger selection in the event." << std::endl;
  bool passHLTa = (useTrigger ? PassTriggers() : true); // default
  if(printLevel > 0) std::cout <<" Select which met will be used in the event." << std::endl;
  if(printLevel > 0) std::cout <<"---VIR-- " << "useTrigger: " << useTrigger<< "PassTriggers()"  << PassTriggers() << "passHLT: " << passHLTa << std::endl;
  //   if(printLevel > 0 && useTrigger==1 && PassTriggers()==1)    std::cout << "---VIR--alla fine SI" << "passHLTa: "<< passHLTa << std::endl;
  //// --- VIR uncomment for default -----
  
  
  if(!passHLTa) continue;
  nCnt[1]++;
  nVir[1]++;
  
  std::cout << " HLT passed              : " << nCnt[1] << " (" << nCnt[1]/float(nCnt[0]) << ") wrt total events" << std::endl;
  std::cout << " HLT passed vir          : " << nVir[1] << " (" << nVir[1]/float(nVir[0]) << ") wrt total events" << std::endl;



//  total events - vir" << std::endl;
//   std::cout << " Total events            : " << nCnt[1] << std::endl;
//   std::cout << " Total events            : " << nCnt[1] << " Total events vir        : " << nVir[1] <<std::endl;
  
  int nFiltered = 0;
  TTree* filterTree = 0;

  if(enableFilter) {
    TFile* filterFile = new TFile(filtered_file_name,"RECREATE");
    filterTree = (TTree*) fChain->GetTree()->CloneTree(0);
    filterTree->SetAutoSave();
  }
  
  
  // open hist file and define histograms
  TFile* fout = new TFile("hist-"+ds+".root","RECREATE");
  
  fout->cd();

  TH1F* h_vtxZ = new TH1F("vtxZ","Z position of the primary vertex;Z (cm);Events",100,-50.0,50.0);
  TH1F* h_bsZ = new TH1F("bsZ","Z position of the beam spot;Z (cm);Events",100,-50.0,50.0);
  TH1F* h_met = new TH1F("met","missing transverse energy;#slash{E}_{T} (GeV);Events",200,0.0,1000.0);
  TH1F* h_sumEt = new TH1F("sumEt","Scalar sum of all calorimeter energy;#sigmaE_{T} (GeV);Events",200,0.0,2000.0);

  TH1F* h_EleMomX = new TH1F("EleMomentumX","Electron momentum X(GeV/c);Events",200,0.0,1000.0);  //---VIR--
  TH1F* h_EleMomY = new TH1F("EleMomentumY","Electron momentum Y(Ge/cV);Events",200,0.0,1000.0);  //---VIR--
  TH1F* h_EleMomZ = new TH1F("EleMomentumZ","Electron momentum Z(GeV/c);Events",200,0.0,1000.0);  //---VIR--
  TH1F* h_ElePt = new TH1F("ElePt","Electron Pt (GeV/c);Events",200,0.0,1000.0);                  //---VIR--
  TH1F* h_ElePhi = new TH1F("ElePhi","Phi of the Electron(GeV)",200,-5.0,5.0);                    //---VIR--
  TH1F* h_EleEta = new TH1F("EleEta","Eta of the Electron(GeV)",200,-3.0,3.0);                    //---VIR--
  TH1F* h_ElePtTrigger = new TH1F("ElePtTrigger","Electron Pt (GeV/c)after passedTRIGGER;Events",200,0.0,1000.0);                  //---VIR--

  
  TH1F* h_MuMomX = new TH1F("MuonMomentumX","Muon momentum X(GeV/c)",200,0.0,1000.0);    //---VIR--
  TH1F* h_MuMomY = new TH1F("MuonMomentumY","Muon momentum Y(GeV/c)",200,0.0,1000.0);    //---VIR--
  TH1F* h_MuMomZ = new TH1F("MuonMomentumZ","Muon momentum Z(GeV/c)",200,0.0,1000.0);    //---VIR--
  TH1F* h_MuPt = new TH1F("MuonPt","Muon Pt (GeV/c)",200,0.0,1000.0);                    //---VIR--
  TH1F* h_MuPhi = new TH1F("MuonPhi","Phi of the Muon(GeV)",200,-5.0,5.0);               //---VIR--
  TH1F* h_MuEta = new TH1F("MuonEta","Eta of the Muon(GeV)",200,-3.0,3.0);               //---VIR--
  TH1F* h_MuPtTrigger = new TH1F("MuonPtTrigger","Muon Pt (GeV/c)after passedTRIGGER",200,0.0,1000.0);                    //---VIR--




  TH1F* h_PhoMomX = new TH1F("PhoMomentumX","Photon momentum X(GeV/c);Events",200,0.0,1000.0);   //---VIR--
  TH1F* h_PhoMomY = new TH1F("PhoMomentumY","Photon momentum Y(Ge/cV);Events",200,0.0,1000.0);   //---VIR--
  TH1F* h_PhoMomZ = new TH1F("PhoMomentumZ","Photon momentum Z(GeV/c);Events",200,0.0,1000.0);   //---VIR--
  TH1F* h_PhoPt = new TH1F("PhoPt","Photon Pt (GeV/c);Events",200,0.0,1000.0);                   //---VIR--
  TH1F* h_PhoPhi = new TH1F("PhoPhi","Phi of the Photon(GeV)",200,-5.0,5.0);                     //---VIR--
  TH1F* h_PhoEta = new TH1F("PhoEta","Eta of the Photon(GeV)",200,-3.0,3.0);                     //---VIR--
  TH1F* h_PhoPtTrigger = new TH1F("PhoPtTrigger","Photon Pt (GeV/c)after passedTRIGGER;Events",200,0.0,1000.0);                   //---VIR--




  // to check duplicated events
  std::map<int, std::set<int> > allEvents;

  // start event looping

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {

    if(printLevel > 0) std::cout << "Get the tree contents." << std::endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0)) ) {
      std::cout << int(jentry) << " events processed with run="
		<< event->runNumber << ", event=" << event->eventNumber << std::endl;
    }


    if(printLevel > 0) std::cout << "Initialize any global variables to be reset per event." << std::endl;

    InitializePerEvent();


    if(printLevel > 0) std::cout << "Apply good run list." << std::endl;
    // uncomment this to use the Json file to flag good data (or bad depending on your outlook)    
    // if(!isInJson(event->runNumber,event->luminosityBlockNumber)) continue;

    // uncomment this to print all ntuple variables
    Print(*event);

    if(printLevel > 0) std::cout << "Check duplicated events for data only." << std::endl;
    if(printLevel > 0) std::cout << "isRealData : " << int(event->isRealData) << std::endl;

    bool duplicateEvent = ! (allEvents[event->runNumber].insert(event->eventNumber)).second;
    if(event->isRealData && duplicateEvent) continue; // uncomment
 
    // remove events filtered by optional met filters
    if(event->isRealData) {
      if(!event->passMetFilters()) continue;}

    if(printLevel > 0) std::cout << "Setup object vectors." << std::endl;

    // classify photon objects

    // loose objects have all standard cuts except for isolation
    std::vector<susy::Photon*>   loose_photons;

    // tight objects hava isolation cuts applied on top of loose objects
    std::vector<susy::Photon*>   tight_photons;

    // same as tight except for nPixelSeeds > 0
    std::vector<susy::Photon*>   ele_photons;

    // same as tight except for reversing either trackIso or sigmaIetaIeta
    std::vector<susy::Photon*>   fake_photons;

    std::vector<susy::CaloJet*>  caloJets;
    std::vector<susy::PFJet*>    pfJets;

    if(printLevel > 0) std::cout << "Find primary vertex in the event." << std::endl;

    TVector3* primVtx = 0;
    if(event->vertices.size() > 0) primVtx = &(event->vertices[0].position);

    if(primVtx) h_vtxZ->Fill(primVtx->Z());
    h_bsZ->Fill(event->beamSpot.Z());

    // ----   MUON   --------------------------------------------------
    // -- VIR
    if(printLevel > 0) std::cout << "Find a muon: "<< event->muons.size() << std::endl;   // -- VIR
    
    TLorentzVector* MuonMome = 0;
    if(event->muons.size() > 0) MuonMome = &(event->muons[0].momentum);
    
    if(MuonMome) h_MuMomX->Fill(MuonMome->X());
    if(MuonMome) h_MuMomY->Fill(MuonMome->Y());
    if(MuonMome) h_MuMomZ->Fill(MuonMome->Z());
    if(MuonMome) h_MuPt->Fill(MuonMome->Pt());
    //    if(MuonMome) h_MuE->Fill(MuonMome->E());
    if(MuonMome) h_MuPhi->Fill(MuonMome->Phi());
    if(MuonMome) h_MuEta->Fill(MuonMome->Eta());
    
    if(event->muons.size() > 0) {
      std::cout << "Found "<< event->muons.size()<<
	" muon, with Pt:" << event->muons[0].momentum.Pt() <<
	" and Energy: " <<event->muons[0].momentum.E()<<     
	" and Phi: "<<event->muons[0].momentum.Phi()<<
	" and eta: "<<event->muons[0].momentum.Eta()<< std::endl;}
    
    
    if(printLevel > 0) std::cout << "Apply trigger selection in the event." << std::endl;
    bool passHLT2 = (useTrigger ? PassTriggers() : true); // default
    std::cout << "Look for muon, with Pt e passato HLT " << passHLT2<< std::endl;  
    
    if(!passHLT2)continue;
    std::cout << "Found 1 muon, with Pt e passato HLT" << passHLT2<< std::endl;  	  
    if(MuonMome) h_MuPtTrigger->Fill(MuonMome->Pt());  //---VIR  trigger prova




    //if(MuonMome) h_MuMomPt->Fill(MuonMome->Pt());
    // -- end VIR

    // ----   ELECTRON   -----------------------------------------------
    if(printLevel > 0) std::cout << "Find an electron in the event." << std::endl;
    
    std::map<TString, std::vector<susy::Electron> >::iterator eleMap = event->electrons.find("gsfElectrons");
    
    if(eleMap != event->electrons.end()) {
      
      for(std::vector<susy::Electron>::iterator it = eleMap->second.begin(); 
	  it != eleMap->second.end(); it++) {
	
	
	std::cout << "Found  electron, with Pt:" << it->momentum.Pt() <<
	  " and Energy: " << it->momentum.Et()  <<
	  " and Phi: "    << it->momentum.Phi() <<
	  " and eta: "    << it->momentum.Eta() <<
	  " and pX: "     << it->momentum.X()   << 
	  " and EB:"      << it->isEB() << std::endl;

	
	float elepx = it->momentum.X();
	float elepy = it->momentum.Y();
	float elepz = it->momentum.Z();
	float elept = it->momentum.Pt();
	float elephi = it->momentum.Phi();
	float eleeta = it->momentum.Eta();

	h_EleMomX->Fill(elepx);
	h_EleMomY->Fill(elepy);
	h_EleMomZ->Fill(elepz);
	h_ElePt->Fill(elept);
	h_ElePhi->Fill(elephi);
	h_EleEta->Fill(eleeta);

	if(printLevel > 0) std::cout << "Apply trigger selection in the event." << std::endl;
	bool passHLT = (useTrigger ? PassTriggers() : true); // default
	std::cout << "Look for electron, with Pt e passato HLT " << passHLT<< std::endl;  
	
	if(!passHLT)continue;
	std::cout << "Found 1 electron, with Pt e passato HLT" << passHLT<< std::endl;  	  
	float eleptTrigger = it->momentum.Pt();
	
	h_ElePtTrigger->Fill(eleptTrigger); //---VIR  trigger prova

	// Et cuts, 100 GeV 
	//if(it->momentum.Et() < 100.0) std::cout << "Find an electron E: " << >momentum.Et() << std::endl;
	
	// -----  electron cuts  -------------

	// fiducial region (barrel: |eta|<1.444, endcap 1.566 < |eta| < 3.0)
	// An electron is considered to be within this ECAL acceptance 
	// if its associated SuperCluster (SC) is within the ECAL acceptance
	// isEE
	// isEB
	// isPF
	
	// sigma_ietaieta (in trigger 0.01 for EB, 0.031 for EE)
        // bool sIetaCut = (it->sigmaIetaIeta < 0.01);
	
	// delta phi and delta eta
	// 
	//     Float_t        deltaPhiSuperClusterTrackAtVtx;
	//     Float_t        deltaPhiSeedClusterTrackAtCalo;
	//     Float_t        deltaPhiEleClusterTrackAtCalo;
	
	//     Float_t        deltaEtaSuperClusterTrackAtVtx;
	//     Float_t        deltaEtaSeedClusterTrackAtCalo;
	//     Float_t        deltaEtaEleClusterTrackAtCalo;
	
	// --- /. electron cuts  -------------
      }
    }
    // -- end VIR
    // ---- end  ELECTRON   -----------------------------------------------
    
      

    // ----   PHOTON   -----------------------------------------------
    if(printLevel > 0) std::cout << "Find loose and tight photons in the event." << std::endl;

    std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");

    if(phoMap != event->photons.end()) {

      for(std::vector<susy::Photon>::iterator it = phoMap->second.begin();
	  it != phoMap->second.end(); it++) {

	std::cout << "Found  photon, with Pt:" << it->momentum.Pt() <<
	  " and Energy: " << it->momentum.Et()  <<
	  " and Phi: "    << it->momentum.Phi() <<
	  " and eta: "    << it->momentum.Eta() <<
	  " and pX: "     << it->momentum.X()   << 
	  " and EB:"      << it->isEB() << std::endl;

	//	std::cout << "Look for photon, with Pt e passato HLT " << passHLT1<< std::endl;  


	//float phopt = it->momentum.Pt();
	float phopx = it->momentum.X();
	float phopy = it->momentum.Y();
	float phopz = it->momentum.Z();
	float phopt = it->momentum.Pt();
	float phophi = it->momentum.Phi();
	float phoeta = it->momentum.Eta();
	
	//h_PhoPt->Fill(phopt);
	h_PhoMomX->Fill(phopx);
	h_PhoMomY->Fill(phopy);
	h_PhoMomZ->Fill(phopz);
	h_PhoPt->Fill(phopt);
	h_PhoPhi->Fill(phophi);
	h_PhoEta->Fill(phoeta);
	


	if(printLevel > 0) std::cout << "Apply trigger selection in the event." << std::endl;
	bool passHLT1 = (useTrigger ? PassTriggers() : true); // default
	std::cout << "Look for electron, with Pt e passato HLT " << passHLT1<< std::endl;  

	if(!passHLT1)continue;
	std::cout << "Found 1 electron, with Pt e passato HLT" << passHLT1<< std::endl;  	  
	float phoptTrigger = it->momentum.Pt();
	
	h_PhoPtTrigger->Fill(phoptTrigger); //---VIR  trigger prova





	// fiducial cuts. Look for only barrel now
	if(!it->isEB()) continue;

	// Et cuts, 25 GeV for trailing photons. Will apply tighter for the leading one.
	if(it->momentum.Et() < 25.0) continue;
	std::cout << "it->momentum.Et(): " <<it->momentum.Et() <<std::endl; // ---- VIR

        // optional Spike cleaning
        if(it->r9 > 1.0) continue;

        // H/E (in trigger, 0.15 for EB, 0.10 for EE)
        bool heCut = (it->hadronicOverEm < 0.05);
        
        // sigma_ietaieta (in trigger 0.014 for EB, 0.034 for EE)
        bool sIetaCut = (it->sigmaIetaIeta < 0.013);

        // Ecal Isolation
        bool ecalIsoCut = (it->ecalRecHitSumEtConeDR04 < 4.2 + 0.006 * it->momentum.Et());

        // Hcal Isolation
        bool hcalIsoCut = (it->hcalTowerSumEtConeDR04() < 2.2 + 0.0025 * it->momentum.Et());

        // Track Isolation
        bool trackIsoCut = (it->trkSumPtHollowConeDR04 < 2.0 + 0.001 * it->momentum.Et());

        bool pixelCut = (it->nPixelSeeds == 0);

        // loose & tight ID variables
        bool looseCut = heCut && ecalIsoCut && hcalIsoCut;
        bool tightCut = looseCut && pixelCut && sIetaCut && trackIsoCut;
        bool eleClass  = looseCut && !pixelCut && sIetaCut && trackIsoCut;
        bool fakeClass = looseCut && pixelCut && !(sIetaCut && trackIsoCut);

        if(looseCut) {
          loose_photons.push_back(&*it);
        }
        if(tightCut) {
          tight_photons.push_back(&*it);
        }
        if(eleClass) {
          ele_photons.push_back(&*it);
        }
        if(fakeClass) {
          fake_photons.push_back(&*it);
        }

      }// for photon
    }// else

    // sort photons by Et
    std::sort(loose_photons.begin(),loose_photons.end(),EtGreater<susy::Photon>);
    std::sort(tight_photons.begin(),tight_photons.end(),EtGreater<susy::Photon>);
    std::sort(ele_photons.begin(),ele_photons.end(),EtGreater<susy::Photon>);
    std::sort(fake_photons.begin(),fake_photons.end(),EtGreater<susy::Photon>);

    // ----   end PHOTON   -----------------------------------------------





    // ----   CALOJETS   -----------------------------------------------

    if(printLevel > 0) std::cout << "Find caloJets in the event." << std::endl;
      
    std::map<TString,susy::CaloJetCollection>::iterator caloJets_it = event->caloJets.find("ak5");

    if(caloJets_it != event->caloJets.end()){

      susy::CaloJetCollection& jetColl = caloJets_it->second;

      for(std::vector<susy::CaloJet>::iterator it = jetColl.begin();
	  it != jetColl.end(); it++) {

	std::map<TString,Float_t>::iterator s_it = it->jecScaleFactors.find("L2L3");
	if (s_it == it->jecScaleFactors.end()) {
	  std::cout << "JEC is not available for this jet!!!" << std::endl;
	  continue;
	}
	float scale = s_it->second;

        if(printLevel > 2) std::cout << "CaloJet stored (" << scale << ")" << std::endl;

	TLorentzVector corrP4 = scale * it->momentum;

	if(std::abs(corrP4.Eta()) > 3.0) continue;

	bool same = false;

	for(std::vector<susy::Photon*>::iterator m_it = tight_photons.begin();
	    m_it != tight_photons.end(); m_it++){
	  if(isSameObject(corrP4,(*m_it)->momentum)){
	    same = true;
	    break;
	  }
	}
	if(same) continue;

	//	if(pt < 20) continue;

	caloJets.push_back(&*it);

      }// for jet
    }// else

    std::sort(caloJets.begin(),caloJets.end(),EtGreater<susy::CaloJet>);


    // ----   PFJETS   -------------------------------------------------
    if(printLevel > 0) std::cout << "Find pfJets in the event." << std::endl;
      
    std::map<TString,susy::PFJetCollection>::iterator pfJets_it = event->pfJets.find("ak5");
    if(pfJets_it == event->pfJets.end()){
      if(event->pfJets.size() > 0) std::cout << "JetCollection is not available!!!" << std::endl;
    }
    else {

      susy::PFJetCollection& jetColl = pfJets_it->second;

      for(std::vector<susy::PFJet>::iterator it = jetColl.begin();
	  it != jetColl.end(); it++) {

	std::map<TString,Float_t>::iterator s_it = it->jecScaleFactors.find("L2L3");
	if (s_it == it->jecScaleFactors.end()) {
	  std::cout << "JEC is not available for this jet!!!" << std::endl;
	  continue;
	}
	float scale = s_it->second;

        if(printLevel > 2) std::cout << "PFJet stored (" << scale << ")" << std::endl;

	TLorentzVector corrP4 = scale * it->momentum;

	if(std::abs(corrP4.Eta()) > 3.0) continue;

	bool same = false;

	for(std::vector<susy::Photon*>::iterator m_it = tight_photons.begin();
	    m_it != tight_photons.end(); m_it++){
	  if(isSameObject(corrP4,(*m_it)->momentum)){
	    same = true;
	    break;
	  }
	}
	if(same) continue;

	//	if(pt < 20) continue;

	pfJets.push_back(&*it);

      }// for jet
    }// else

    std::sort(pfJets.begin(),pfJets.end(),EtGreater<susy::PFJet>);

    // ----   TRIGGER   ------------------------------------------------

    //     // --- VIR uncomment for default -----
    //     if(printLevel > 0) std::cout << "Apply trigger selection in the event." << std::endl;
    //     bool passHLT3 = (useTrigger ? PassTriggers() : true); // default
    //     if(printLevel > 0) std::cout << "Select which met will be used in the event." << std::endl;
    //     if(printLevel > 0) std::cout <<"---VIR--" << "useTrigger: " << useTrigger<< "PassTriggers()"  << PassTriggers() << std::endl;
    //     if(printLevel > 0) std::cout <<"---VIR--" << "passHLT: " << passHLT3<< std::endl;
    //     //if(printLevel > 0 && useTrigger==1 && PassTriggers()==1) {
    //     //  std::cout << "---VIR--alla fine SI" << "passHLT: "<< passHLT << std::endl;}
    //     // --- VIR uncomment for default -----
    
    
    std::map<TString, susy::MET>::iterator met_it = event->metMap.find("pfMet");
    if(met_it == event->metMap.end()) {
      std::cout << "MET map is not available!!!" << std::endl;
      continue;
    }
    susy::MET* met = &(met_it->second);

    if(printLevel > 0) {
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "              event summary" << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "loose_photons     : " << loose_photons.size() << std::endl;
      std::cout << "tight_photons     : " << tight_photons.size() << std::endl;
      std::cout << "ele_photons       : " << ele_photons.size() << std::endl;
      std::cout << "fake_photons      : " << fake_photons.size() << std::endl;
      std::cout << "caloJets          : " << caloJets.size() << std::endl;
      std::cout << "pfJets            : " << pfJets.size() << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "met               : " << met->met() << std::endl;
    } 


    if(printLevel > 0) std::cout << "Apply event level cuts from now on..." << std::endl;


    // filter conditions

    if(enableFilter) {
      bool filterThis = (loose_photons.size() > 0);
      if(filterThis) {
	nFiltered++;
	filterTree->Fill();
      }
    }// if(enableFilter)
    
    
    // event counter
    
    nCnt[0]++; // total number of events
    
    if(!passHLTa) continue;
    nCnt[1]++;
    
    if(loose_photons.size() == 0) continue;
    nCnt[2]++;
    
    h_met->Fill(met->met());
    h_sumEt->Fill(met->sumEt);
    // h_EleE->Fill(met->met()); // ---VIR--


    // two photons
    if(tight_photons.size() >= 2) {
      nCnt[3]++;
    }

    // one photon + one electron
    if(tight_photons.size() >= 1 && ele_photons.size() >= 1) {
      nCnt[4]++;
    }


    // two electrons
    if(ele_photons.size() >= 2) {
      nCnt[5]++;
    }

    // one photon + one fake
    if(tight_photons.size() >= 1 && fake_photons.size() >= 1) {
      nCnt[6]++;
    }

    // two fakes
    if(fake_photons.size() >= 2) {
      nCnt[7]++;
    }


    if(met->met() < 50.0) continue;

    nCnt[8]++;

  } // for jentry


  // end of event loop and print summary

  std::cout << " ----------------- Job Summary ----------------- " << std::endl;
  std::cout << " Total events            : " << nCnt[0] << std::endl;
  std::cout << " HLT passed              : " << nCnt[1] << " (" << nCnt[1]/float(nCnt[0]) << ") wrt total events" << std::endl;
  std::cout << " loose_photons > 0       : " << nCnt[2] << " (" << nCnt[2]/float(nCnt[1]) << ") wrt HLT" << std::endl;
  std::cout << " gg events               : " << nCnt[3] << " (" << nCnt[3]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " ge events               : " << nCnt[4] << " (" << nCnt[4]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " ee events               : " << nCnt[5] << " (" << nCnt[5]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " gf events               : " << nCnt[6] << " (" << nCnt[6]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " ff events               : " << nCnt[7] << " (" << nCnt[7]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " met > 50 GeV            : " << nCnt[8] << " (" << nCnt[8]/float(nCnt[1]) << ")" << std::endl;

  if(enableFilter) {
    std::cout << " --------------- Filtered events --------------- " << std::endl;
    std::cout << " filtered events         : " << nFiltered << " (" << nFiltered/float(nCnt[0]) << ")" << std::endl;
  }
  std::cout << " ----------------------------------------------- " << std::endl;

  // close the output file

  fout->cd();
  fout->Write();
  fout->Close();

  if(enableFilter) {
    filterTree->GetCurrentFile()->cd();
    filterTree->GetCurrentFile()->Write();
    filterTree->GetCurrentFile()->Close();
  }

}

