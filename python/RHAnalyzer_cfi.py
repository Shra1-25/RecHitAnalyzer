import FWCore.ParameterSet.Config as cms 

fevt = cms.EDAnalyzer('RecHitAnalyzer'
    , reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
    , mode = cms.string("JetLevel")
    # Jet level cfg
    , nJets = cms.int32(-1)
    , minJetPt = cms.double(35.)
    , maxJetEta = cms.double(2.4)
    , z0PVCut  = cms.double(0.1)
    )
    
