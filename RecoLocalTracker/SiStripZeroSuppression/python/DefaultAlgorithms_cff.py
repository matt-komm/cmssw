import FWCore.ParameterSet.Config as cms

DefaultAlgorithms = cms.PSet(

    ## Pedestal subtraction ----------------
    PedestalSubtractionFedMode = cms.bool(True),

    ## Baseline finder ---------------------
    ## Supported CMN modes: Median, Percentile, IteratedMedian, TT6, FastLinear
    CommonModeNoiseSubtractionMode = cms.string('Median'),     

    #CutToAvoidSignal = cms.double(3.0), ## for TT6
    
    #Percentile = cms.double(25.0),      ## for Percentile

    CutToAvoidSignal = cms.double(2.0),  ## for IteratedMedian (needed by APVRestorer)
    Iterations = cms.int32(3),           ##

    ## APV restoration ---------------------
    ## Supported inspect modes: BaselineFollower, AbnormalBaseline, Null, BaselineAndSaturation
    APVInspectMode = cms.string("BaselineFollower"),
    ForceNoRestore = cms.bool(False),
    SelfSelectRestoreAlgo = cms.bool(False),
    useRealMeanCM = cms.bool(False),
    DeltaCMThreshold = cms.uint32(20),       # for BaselineFollower inspect
    distortionThreshold = cms.uint32(40),    # " "
    Fraction = cms.double(0.2),              # for AbnormalBaseline inspect
    Deviation = cms.uint32(25),              # " "
    restoreThreshold = cms.double(0.5),      # for Null inspect
    nSaturatedStrip = cms.uint32(2),         # for BaselineAndSaturation inspect

    ## Supported restore modes: Flat, BaselineFollower, IterativeMedian
    APVRestoreMode = cms.string("BaselineFollower"),
    nSigmaNoiseDerTh = cms.uint32(4),        # threshold for rejecting hit strips: nSigma * noise
    consecThreshold = cms.uint32(5),         # minimum length of flat region
    hitStripThreshold = cms.uint32(40),      # height above median when strip is definitely a hit
    nSmooth = cms.uint32(9),                 # for smoothing and local minimum determination (odd number)
    minStripsToFit = cms.uint32(4),          # minimum strips to try spline algo (otherwise default to median)
    
    ## Zero suppression --------------------
    SiStripFedZeroSuppressionMode = cms.uint32(4),
    TruncateInSuppressor = cms.bool(True)
    
    )
