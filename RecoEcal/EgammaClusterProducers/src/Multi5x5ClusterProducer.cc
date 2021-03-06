// C/C++ headers
#include <iostream>
#include <vector>
#include <memory>

// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

// Reconstruction Classes
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"

// Geometry
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"

// EgammaCoreTools
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

// Class header file
#include "RecoEcal/EgammaClusterProducers/interface/Multi5x5ClusterProducer.h"


Multi5x5ClusterProducer::Multi5x5ClusterProducer(const edm::ParameterSet& ps)
{
  // The verbosity level
  std::string verbosityString = ps.getParameter<std::string>("VerbosityLevel");
  if      (verbosityString == "DEBUG")   verbosity = Multi5x5ClusterAlgo::pDEBUG;
  else if (verbosityString == "WARNING") verbosity = Multi5x5ClusterAlgo::pWARNING;
  else if (verbosityString == "INFO")    verbosity = Multi5x5ClusterAlgo::pINFO;
  else                                   verbosity = Multi5x5ClusterAlgo::pERROR;

  // Parameters to identify the hit collections
  barrelHitProducer_   = ps.getParameter<std::string>("barrelHitProducer");
  endcapHitProducer_   = ps.getParameter<std::string>("endcapHitProducer");
  barrelHitCollection_ = ps.getParameter<std::string>("barrelHitCollection");
  endcapHitCollection_ = ps.getParameter<std::string>("endcapHitCollection");

  // should cluster algo be run in barrel and endcap?
  doEndcap_ = ps.getParameter<bool>("doEndcap");
  doBarrel_ = ps.getParameter<bool>("doBarrel");

  // The names of the produced cluster collections
  barrelClusterCollection_  = ps.getParameter<std::string>("barrelClusterCollection");
  endcapClusterCollection_  = ps.getParameter<std::string>("endcapClusterCollection");

  // Island algorithm parameters
  double barrelSeedThreshold = ps.getParameter<double>("IslandBarrelSeedThr");
  double endcapSeedThreshold = ps.getParameter<double>("IslandEndcapSeedThr");

  std::vector<int> v_chstatus = ps.getParameter<std::vector<int> >("RecHitFlagToBeExcluded");

  // Parameters for the position calculation:
  edm::ParameterSet posCalcParameters = 
    ps.getParameter<edm::ParameterSet>("posCalcParameters");
  posCalculator_ = PositionCalc(posCalcParameters);

  // Produces a collection of barrel and a collection of endcap clusters
  produces< reco::BasicClusterCollection >(endcapClusterCollection_);
  produces< reco::BasicClusterCollection >(barrelClusterCollection_);

  island_p = new Multi5x5ClusterAlgo(barrelSeedThreshold, endcapSeedThreshold,  v_chstatus, posCalculator_,verbosity);

  nEvt_ = 0;
}


Multi5x5ClusterProducer::~Multi5x5ClusterProducer()
{
  delete island_p;
}


void Multi5x5ClusterProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{

  if (doEndcap_) {
    clusterizeECALPart(evt, es, endcapHitProducer_, endcapHitCollection_, endcapClusterCollection_, reco::CaloID::DET_ECAL_ENDCAP); 
  }
  if (doBarrel_) {
    clusterizeECALPart(evt, es, barrelHitProducer_, barrelHitCollection_, barrelClusterCollection_, reco::CaloID::DET_ECAL_BARREL);
  }

  nEvt_++;
}


const EcalRecHitCollection * Multi5x5ClusterProducer::getCollection(edm::Event& evt,
                                                                  const std::string& hitProducer_,
                                                                  const std::string& hitCollection_)
{
  edm::Handle<EcalRecHitCollection> rhcHandle;
  try
    {
      evt.getByLabel(hitProducer_, hitCollection_, rhcHandle);
      if (!(rhcHandle.isValid())) 
	{
	  std::cout << "could not get a handle on the EcalRecHitCollection!" << std::endl;
	  return 0;
	}
    }
  catch ( cms::Exception& ex ) 
    {
      edm::LogError("Multi5x5ClusterProducerError") << "Error! can't get the product " << hitCollection_.c_str() ;
      return 0;
    }
  return rhcHandle.product();
}


void Multi5x5ClusterProducer::clusterizeECALPart(edm::Event &evt, const edm::EventSetup &es,
                                               const std::string& hitProducer,
                                               const std::string& hitCollection,
                                               const std::string& clusterCollection,
                                               const reco::CaloID::Detectors detector)
{
  // get the hit collection from the event:
  const EcalRecHitCollection *hitCollection_p = getCollection(evt, hitProducer, hitCollection);

  // get the geometry and topology from the event setup:
  edm::ESHandle<CaloGeometry> geoHandle;
  es.get<CaloGeometryRecord>().get(geoHandle);

  const CaloSubdetectorGeometry *geometry_p;
  CaloSubdetectorTopology *topology_p;

  if (detector == reco::CaloID::DET_ECAL_BARREL) 
    {
      geometry_p = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
      topology_p = new EcalBarrelTopology(geoHandle);
    }
  else
    {
      geometry_p = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
      topology_p = new EcalEndcapTopology(geoHandle); 
   }

  const CaloSubdetectorGeometry *geometryES_p;
  geometryES_p = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);

  // Run the clusterization algorithm:
  reco::BasicClusterCollection clusters;
  clusters = island_p->makeClusters(hitCollection_p, geometry_p, topology_p, geometryES_p, detector);

  // create an auto_ptr to a BasicClusterCollection, copy the barrel clusters into it and put in the Event:
  std::auto_ptr< reco::BasicClusterCollection > clusters_p(new reco::BasicClusterCollection);
  clusters_p->assign(clusters.begin(), clusters.end());
  edm::OrphanHandle<reco::BasicClusterCollection> bccHandle;
  if (detector == reco::CaloID::DET_ECAL_BARREL) 
    bccHandle = evt.put(clusters_p, barrelClusterCollection_);
  else
    bccHandle = evt.put(clusters_p, endcapClusterCollection_);

  delete topology_p;
}
