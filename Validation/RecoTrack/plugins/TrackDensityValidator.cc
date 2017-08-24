#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "Validation/RecoTrack/interface/MultiTrackValidatorBase.h"
#include "Validation/RecoTrack/interface/MTVHistoProducerAlgoForTracker.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelectorBase.h"
#include "SimTracker/TrackAssociation/interface/ParametersDefinerForTP.h"
#include "CommonTools/Utils/interface/DynArray.h"

#include "Validation/RecoTrack/interface/TrackDensityValidator.h"

TrackDensityValidator::TrackDensityValidator(const edm::ParameterSet& pset)
{
}

void TrackDensityValidator::analyze(edm::Event const&, edm::EventSetup const&)
{
}

void TrackDensityValidator::bookHistograms(DQMStore::IBooker &i, edm::Run const&, edm::EventSetup const&)
{
}
