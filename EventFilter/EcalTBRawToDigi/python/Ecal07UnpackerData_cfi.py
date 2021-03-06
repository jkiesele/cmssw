import FWCore.ParameterSet.Config as cms

ecalEBunpacker = cms.EDProducer("EcalDCCTB07UnpackingModule",
    fedRawDataCollectionTag = cms.InputTag('rawDataCollector'),
    ccuIDs = cms.untracked.vint32(1, 71, 80, 45),
    statusIDs = cms.untracked.vint32(1, 2, 3, 4),
    positionIDs = cms.untracked.vint32(6, 2, 5, 1),
    # use this for h2 data
    tbName = cms.untracked.string('h2'),
    # use this for h4 data
    #   untracked string     tbName = "h4"
    #   include "EventFilter/EcalTBRawToDigi/data/h4_mapping.cfi"
    produceEEdigi = cms.untracked.bool(True),
    stripIDs = cms.untracked.vint32(1, 2, 3, 4, 5, 
        5, 4, 3, 2, 1, 
        1, 2, 3, 4, 5, 
        5, 4, 3, 2, 1, 
        1, 2, 3, 4, 5, 
        5, 4, 3, 2, 1, 
        1, 2, 3, 4, 5, 
        5, 4, 3, 2, 1, 
        1, 2, 3, 4, 5, 
        5, 4, 3, 2, 1, 
        1, 2, 3, 4, 5, 
        5, 4, 3, 2, 1, 
        1, 2, 3, 4, 5, 
        5, 4, 3, 2, 1, 
        1, 2, 3, 4, 5, 
        5, 4, 3, 2, 1, 
        1, 2, 3, 4, 5, 
        5, 4, 3, 2, 1, 
        1, 2, 3, 4, 5, 
        5, 4, 3, 2, 1),
    towerIDs = cms.untracked.vint32(1, 1, 1, 1, 1, 
        2, 2, 2, 2, 2, 
        1, 1, 1, 1, 1, 
        2, 2, 2, 2, 2, 
        1, 1, 1, 1, 1, 
        2, 2, 2, 2, 2, 
        1, 1, 1, 1, 1, 
        2, 2, 2, 2, 2, 
        1, 1, 1, 1, 1, 
        2, 2, 2, 2, 2, 
        5, 5, 5, 5, 5, 
        6, 6, 6, 6, 6, 
        5, 5, 5, 5, 5, 
        6, 6, 6, 6, 6, 
        5, 5, 5, 5, 5, 
        6, 6, 6, 6, 6, 
        5, 5, 5, 5, 5, 
        6, 6, 6, 6, 6, 
        5, 5, 5, 5, 5, 
        6, 6, 6, 6, 6),
    produceEBdigi = cms.untracked.bool(False),
    ics = cms.untracked.vint32(1, 2, 3, 4, 5, 
        6, 7, 8, 9, 10, 
        21, 22, 23, 24, 25, 
        26, 27, 28, 29, 30, 
        41, 42, 43, 44, 45, 
        46, 47, 48, 49, 50, 
        61, 62, 63, 64, 65, 
        66, 67, 68, 69, 70, 
        81, 82, 83, 84, 85, 
        86, 87, 88, 89, 90, 
        101, 102, 103, 104, 105, 
        106, 107, 108, 109, 110, 
        121, 122, 123, 124, 125, 
        126, 127, 128, 129, 130, 
        141, 142, 143, 144, 145, 
        146, 147, 148, 149, 150, 
        161, 162, 163, 164, 165, 
        166, 167, 168, 169, 170, 
        181, 182, 183, 184, 185, 
        186, 187, 188, 189, 190),
    channelIDs = cms.untracked.vint32(1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5, 
        1, 1, 1, 1, 1, 
        1, 1, 1, 1, 1, 
        2, 2, 2, 2, 2, 
        2, 2, 2, 2, 2, 
        3, 3, 3, 3, 3, 
        3, 3, 3, 3, 3, 
        4, 4, 4, 4, 4, 
        4, 4, 4, 4, 4, 
        5, 5, 5, 5, 5, 
        5, 5, 5, 5, 5)
)


