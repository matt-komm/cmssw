///
/// Description: Firmware headers
///
/// Implementation:
///    Concrete firmware implementations
///
/// \author: Jim Brooke - University of Bristol
/// Modified: Adam Elwood - ICL

//
//

#ifndef Stage2Layer2JetAlgorithmFirmware_H
#define Stage2Layer2JetAlgorithmFirmware_H

#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2JetAlgorithm.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"

namespace l1t {

  // Imp1 is for v1 and v2
  class Stage2Layer2JetAlgorithmFirmwareImp1 : public Stage2Layer2JetAlgorithm {
  public:
    Stage2Layer2JetAlgorithmFirmwareImp1(CaloParamsHelper* params);
    ~Stage2Layer2JetAlgorithmFirmwareImp1() override;
    void processEvent(const std::vector<CaloTower> & towers,
			      std::vector<Jet> & jets, std::vector<Jet> & alljets) override;

    void create(const std::vector<CaloTower> & towers,
	                      std::vector<Jet> & jets, std::vector<Jet> & alljets, std::string PUSubMethod);

    void accuSort(std::vector<Jet> & jets);

    void calibrate(std::vector<Jet> & jets, int calibThreshold, bool isAllJets);

    double calibFit(double, double*);
    double calibFitErr(double,double*);
    
    void calculateCentralSums(l1t::Jet& jet, const std::vector<l1t::CaloTower> & towers) const;
    void calculateDonutSums(l1t::Jet& jet, const std::vector<l1t::CaloTower> & towers) const;


    int donutPUEstimate(int jetEta, int jetPhi, int size,
                        const std::vector<l1t::CaloTower> & towers);

    int chunkyDonutPUEstimate(Jet & jet, int pos,
                              const std::vector<l1t::CaloTower> & towers,const std::string& PUSubMethodCfg);

  private:

    CaloParamsHelper* const params_;

  };

}

#endif
