import FWCore.ParameterSet.Config as cms

from RecoTracker.DeDx.dedxTruncated40_cfi import *
from RecoTracker.DeDx.dedxMedian_cfi import *
from RecoTracker.DeDx.dedxHarmonic2_cfi import *
from RecoTracker.DeDx.dedxUnbinned_cfi import *

doAlldEdXEstimators = cms.Sequence(dedxTruncated40 + dedxMedian + dedxHarmonic2)

