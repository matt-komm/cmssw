# The following comments couldn't be translated into the new config version:

#keep the logging output to a nice level

import FWCore.ParameterSet.Config as cms

process = cms.Process("RPCSegmentEff")

process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("DQMServices.Components.MEtoEDMConverter_cfi")
process.load("DQMServices.Core.DQM_cfg")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START3X_V18::All'
process.prefer("GlobalTag")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('/store/data/Commissioning10/Cosmics/RECO/v3/000/127/337/122E1F96-F416-DF11-B529-000423D6CAF2.root')
)

process.load("DQM.RPCMonitorDigi.RPCEfficiency_cfi")

process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring('/tmp/carrillo/RPCEfficiency.log')
)


process.FEVT = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *_MEtoEDMConverter_*_*'),
    fileName = cms.untracked.string('/tmp/carrillo/first.root')
)

process.p = cms.Path(process.rpcEfficiency*process.MEtoEDMConverter)
process.outpath = cms.EndPath(process.FEVT)

process.DQM.collectorHost = ''
process.DQM.collectorPort = 9090
process.DQM.debug = False


