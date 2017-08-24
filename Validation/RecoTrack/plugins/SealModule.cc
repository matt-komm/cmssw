#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "Validation/RecoTrack/interface/MultiTrackValidator.h"
#include "Validation/RecoTrack/interface/MultiTrackValidatorGenPs.h"

#include "Validation/RecoTrack/interface/SiStripTrackingRecHitsValid.h"
#include "Validation/RecoTrack/interface/SiPixelTrackingRecHitsValid.h"

#include "Validation/RecoTrack/interface/TrackDensityValidator.h"

DEFINE_FWK_MODULE(MultiTrackValidator);
DEFINE_FWK_MODULE(MultiTrackValidatorGenPs);
DEFINE_FWK_MODULE(SiStripTrackingRecHitsValid);
DEFINE_FWK_MODULE(SiPixelTrackingRecHitsValid);

DEFINE_FWK_MODULE(TrackDensityValidator);
