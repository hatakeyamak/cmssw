import FWCore.ParameterSet.Config as cms

calotowersAnalyzer = cms.EDAnalyzer("CaloTowersValidation",
    outputFile                 = cms.untracked.string(''),
    CaloTowerCollectionLabel   = cms.untracked.InputTag('towerMaker'),
    PFCandidateCollectionLabel = cms.untracked.InputTag('particleFlow'),
    hcalselector               = cms.untracked.string('all'),
    mc                         = cms.untracked.string('yes')
)
