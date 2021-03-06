#ifndef RecoEcal_EgammaClusterProducers_EgammaSCCorrectionMaker_h
#define RecoEcal_EgammaClusterProducers_EgammaSCCorrectionMaker_h

// -*- C++ -*-
//
// Package:    EgammaSCCorrectionMaker
// Class:      EgammaSCCorrectionMaker
// 
/**\class EgammaSCCorrectionMaker EgammaSCCorrectionMaker.cc EgammaSCCorrectionMaker/EgammaSCCorrectionMaker/src/EgammaSCCorrectionMaker.cc

 Description: Producer of corrected SuperClusters

*/
//
// Original Author:  Dave Evans
//         Created:  Thu Apr 13 15:50:17 CEST 2006
// $Id: EgammaSCCorrectionMaker.h,v 1.9 2009/03/25 17:50:24 ferriff Exp $
//
//

#include <memory>
#include <string>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

#include "RecoEcal/EgammaClusterAlgos/interface/EgammaSCEnergyCorrectionAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h" 
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h" 

class EgammaSCCorrectionMaker : public edm::EDProducer {
	
   public:
     explicit EgammaSCCorrectionMaker(const edm::ParameterSet&);
     ~EgammaSCCorrectionMaker();
     virtual void produce(edm::Event&, const edm::EventSetup&);

   private:

     EcalClusterFunctionBaseClass* EnergyCorrection_;

     // the debug level
     EgammaSCEnergyCorrectionAlgo::VerbosityLevel verbosity_;

     // pointer to the correction algo object
     EgammaSCEnergyCorrectionAlgo *energyCorrector_;
    
     // vars for the correction algo
     bool applyEnergyCorrection_;
     //     bool oldEnergyScaleCorrection_;
     double sigmaElectronicNoise_;
     double etThresh_;
     
     // vars to get products
     edm::InputTag rHInputProducer_;
     edm::InputTag sCInputProducer_;

     reco::CaloCluster::AlgoId sCAlgo_;
     std::string outputCollection_;

};
#endif
