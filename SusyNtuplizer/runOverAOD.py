import FWCore.ParameterSet.Config as cms

# change this to 0 if you run on MC files
realData = 1

# change this to 0 if you run on FullSim MC
isFastSim = 1

process = cms.Process("RA3")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 50
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring( )
                            )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections/Configuration/JetCorrectionServices_cff")
process.load("JetMETCorrections/Type1MET/caloMETCorrections_cff")

# SusyNtuplizer options
process.load("SusyAnalysis.SusyNtuplizer.susyNtuplizer_cfi")
process.susyNtuplizer.debugLevel = cms.int32(0)
process.susyNtuplizer.bTagCollectionTags = cms.vstring()
process.susyNtuplizer.storePFJetPartonMatches = cms.bool(False)
if isFastSim:
    process.susyNtuplizer.isFastSim = cms.bool(True)

process.metAnalysisSequence = cms.Sequence(process.producePFMETCorrections*
                                           process.produceCaloMETCorrections)

#Calculate rho restricted to barrel for photon pileup subtraction
process.kt6PFJetsRhoBarrelOnly = process.kt4PFJets.clone(
    src = cms.InputTag('particleFlow'),
    rParam = cms.double(0.6),
    #Eta range of jets to be considered for Rho calculation
    #Should be at most (jet acceptance - jet radius)
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax=cms.double(1.4442),
    #Eta range of ghost jets to be considered for Rho calculation - must be greater than Rho_EtaMax
    Ghost_EtaMax=cms.double(2.5)
    )

#Calculate rho restricted to barrel for photon pileup subtraction
process.kt6PFJetsRho25 = process.kt4PFJets.clone(
    src = cms.InputTag('particleFlow'),
    rParam = cms.double(0.6),
    #Eta range of jets to be considered for Rho calculation
    #Should be at most (jet acceptance - jet radius)
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax=cms.double(2.5)
    )

# HBHENoiseFilterResultProducer
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

# HCAL laser events filter
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.taggingMode = cms.bool(True)

# EcalDeadCellTriggerPrimitiveFilter
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

# EcalDeadCellBoundaryEnergyFilter
process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(True)

# Tracking failure filter
process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.JetSource = cms.InputTag('ak5PFJetsL2L3Residual')
process.trackingFailureFilter.taggingMode = cms.bool(True)

# EE Bad SC Filter
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.taggingMode = cms.bool(True)

#Add up all MET filters
if realData or not isFastSim:
    process.metFiltersSequence = cms.Sequence(
        process.HBHENoiseFilterResultProducer *
        process.hcalLaserEventFilter *
        process.EcalDeadCellTriggerPrimitiveFilter *
        process.EcalDeadCellBoundaryEnergyFilter *
        process.goodVertices *
        process.trackingFailureFilter *
	process.eeBadScFilter
        )
else:
    process.metFiltersSequence = cms.Sequence(
        process.hcalLaserEventFilter *
        process.EcalDeadCellTriggerPrimitiveFilter *
        process.EcalDeadCellBoundaryEnergyFilter *  
        process.goodVertices *
        process.trackingFailureFilter *
	process.eeBadScFilter
        )

# IsoDeposit
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.isoDeposit = cms.Sequence( process.pfParticleSelectionSequence + process.eleIsoSequence + process.phoIsoSequence)

if realData:
    process.source.fileNames = cms.untracked.vstring(
	'/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/190/706/DA8B61A9-BE83-E111-8BCB-001D09F2906A.root'
        )
    process.GlobalTag.globaltag = 'GR_R_52_V9::All'

    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3Residual")
    # JEC for data
    process.jet = cms.Sequence(
        # CaloJets
        process.ak5CaloJetsL2L3Residual * process.ak5CaloJetsL1L2L3Residual *
        # PFJets
        process.ak5PFJetsL2L3Residual * process.ak5PFJetsL1FastL2L3Residual *
        # Barrel only Rho calculation
        process.kt6PFJetsRhoBarrelOnly *
	# Rho25 calculation
	process.kt6PFJetsRho25
        )

else:
    process.source.fileNames = cms.untracked.vstring(
        'dcap:///pnfs/cms/WAX/resilient/lpcpjm/PrivateMC/BinoSignalPoints_5_7_11/reco/1250_1200_225/reco_1250_1200_225_1.root'
        )
    process.GlobalTag.globaltag = 'START52_V11::All'
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3")
    # JEC for MC
    process.jet = cms.Sequence(
        # CaloJets
        process.ak5CaloJetsL2L3 * process.ak5CaloJetsL1L2L3 *
        # PFJets
        process.ak5PFJetsL2L3 * process.ak5PFJetsL1FastL2L3 *
        # Barrel only Rho calculation
        process.kt6PFJetsRhoBarrelOnly *
	# Rho25 calculation
	process.kt6PFJetsRho25
        )
    process.trackingFailureFilter.JetSource = cms.InputTag('ak5PFJetsL2L3')
    if isFastSim:
	process.susyNtuplizer.muonIDCollectionTags = cms.vstring()


process.p = cms.Path(
    process.metAnalysisSequence *
    process.jet *
    process.metFiltersSequence *
    process.isoDeposit *
    process.susyNtuplizer
    )
