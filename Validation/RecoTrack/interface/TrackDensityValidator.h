#ifndef TrackDensityValidator_h
#define TrackDensityValidator_h


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "Validation/RecoTrack/interface/MultiTrackValidatorBase.h"
#include "Validation/RecoTrack/interface/MTVHistoProducerAlgoForTracker.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelectorBase.h"
#include "SimTracker/TrackAssociation/interface/ParametersDefinerForTP.h"
#include "CommonTools/Utils/interface/DynArray.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "MagneticField/Engine/interface/MagneticField.h" 
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "CommonTools/RecoAlgos/interface/CosmicTrackingParticleSelector.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <DQMServices/Core/interface/DQMStore.h>

#include "FWCore/Framework/interface/ConsumesCollector.h"

class TrackDensityValidator:
    public DQMEDAnalyzer
    
{
    public:
        TrackDensityValidator(const edm::ParameterSet& pset);
        virtual void analyze(edm::Event const&, edm::EventSetup const&);
        virtual void bookHistograms(DQMStore::IBooker &i, edm::Run const&, edm::EventSetup const&);
};

#endif
