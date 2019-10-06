// Original Author:  Dongwook Jang
// $Id: ana.C,v 1.8 2011/11/01 22:14:51 dwjang Exp $
//
// Jet energy correction is possible at ntuple level.
// $ cd ../jec/JetMETObjects
// $ make
// This will create a shared library in jec/lib
// which is included below as libJetMETObjects.so
//
// Come back to this directory and do
// $ make
// $ root -b -q -l ana.C
// will produce hist_"physics"_"ds".root

// void ana(TString ds="400_420_375", TString physics="binolike") {
// void ana(TString ds="400_1520_375", TString physics="binolike") {
// void ana(TString ds="1500_420_375", TString physics="binolike") {
// void ana(TString ds="1500_1520_375", TString physics="binolike") {

//void ana(TString ds="400_420_375-Mu22_Photon22_CaloIdL_v5", TString physics="binolike") {
//void ana(TString ds="400_1520_375-Mu22_Photon22_CaloIdL_v5", TString physics="binolike") {
//void ana(TString ds="1500_420_375-Mu22_Photon22_CaloIdL_v5", TString physics="binolike") {
//void ana(TString ds="1500_1520_375-Mu22_Photon22_CaloIdL_v5", TString physics="binolike") {

//void ana(TString ds="400_420_375-IsoMu24_v15", TString physics="binolike") {
//void ana(TString ds="400_1520_375-IsoMu24_v15", TString physics="binolike") {
//void ana(TString ds="1500_420_375-IsoMu24_v15", TString physics="binolike") {
//void ana(TString ds="1500_1520_375-IsoMu24_v15", TString physics="binolike") {

//void ana(TString ds="400_420_375-IsoMu15_eta2p1_L1ETM20_v5", TString physics="binolike") {
//void ana(TString ds="400_1520_375-IsoMu15_eta2p1_L1ETM20_v5", TString physics="binolike") {
//void ana(TString ds="1500_420_375-IsoMu15_eta2p1_L1ETM20_v5", TString physics="binolike") {
void ana(TString ds="1500_1520_375-IsoMu15_eta2p1_L1ETM20_v5", TString physics="binolike") {

  gSystem->Load("libSusyEvent.so");

  // Look ../jec/JetMETObjects/README
  gSystem->Load("../jec/lib/libJetMETObjects.so");

  // Printing utility for ntuple variables
  gROOT->LoadMacro("SusyEventPrinter.cc+");

  // Main analysis code
  gROOT->LoadMacro("SusyEventAnalyzer.cc+");

  // chain of inputs
  TChain* chain = new TChain("susyTree");
  // ----- MC ---------------------------------
  //chain->Add("/tmp/azzolini/MC/tree_400_420_375.root");
  //chain->Add("/tmp/azzolini/MC/tree_400_1520_375.root");
  //chain->Add("/tmp/azzolini/MC/tree_1500_420_375.root");
  chain->Add("/tmp/azzolini/MC/tree_1500_1520_375.root");

  // ---- DATA --------------------------------
  //chain->Add("/tmp/azzolini/susyEvents.root");
  //  chain->Add("../susyEvents.root");
  //  chain->Add("dcap:///pnfs/cms/WAX/resilient/lpcpjm/SusyNtuples/cms423v2_v1/Run2011A-May10ReReco-v1/Photon/susyEvent_1_1_dLs.root");
  
  SusyEventAnalyzer* sea = new SusyEventAnalyzer(chain);
  
  // configuration parameters
  // any values given here will replace the default values
  sea->SetDataset(physics+"_"+ds);        // dataset name
  sea->SetPrintInterval(1e4);             // print frequency
  sea->SetPrintLevel(1);                  // print level for event contents --- VIR
  //sea->SetPrintLevel(0);                  // print level for event contents --- default
  //sea->SetUseTrigger(false); // default
  sea->SetUseTrigger(true); // --- VIR --

  //  sea->AddHltName("HLT_Mu22_Photon22_CaloIdL_v5");  // --- VIR ---add HLT trigger path name
  //  sea->AddHltName("HLT_IsoMu24_v15");               // --- VIR ---add HLT trigger path name
      sea->AddHltName("HLT_IsoMu15_eta2p1_L1ETM20_v5"); // --- VIR ---add HLT trigger path name
  
  //  sea->AddHltName("HLT_Photon36_CaloIdL_Photon22_CaloIdL");  // add HLT trigger path name
  //  sea->AddHltName("HLT_Photon32_CaloIdL_Photon26_CaloIdL");  // add HLT trigger path name
  sea->SetFilter(false);               // filter events passing final cuts
  sea->SetProcessNEvents(999999999);   // number of events to be processed --- VIR
  //sea->SetProcessNEvents(200);             // number of events to be processed --- default
  
  // as an example -- add your favorite Json here.  More than one can be "Include"ed
  //  sea->IncludeAJson("Cert_161079-161352_7TeV_PromptReco_Collisions11_JSON_noESpbl_v2.txt");
  //sea->IncludeAJson("anotherJSON.txt");

  TStopwatch ts;

  ts.Start();

  sea->Loop();

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
