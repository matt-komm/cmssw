#ifndef __L1Analysis_L1AnalysisL1UpgradeDataFormat_H__
#define __L1Analysis_L1AnalysisL1UpgradeDataFormat_H__

//-------------------------------------------------------------------------------
// Created 20/04/2010 - E. Conte, A.C. Le Bihan
// 
// 
// Original code : L1TriggerDPG/L1Ntuples/L1UpgradeTreeProducer - Jim Brooke
//-------------------------------------------------------------------------------


#include <vector>
#include <array>

namespace L1Analysis
{

  // copied from DataFormats/L1Trigger/interface/EtSum.h, for use in standalone ROOT macros which use this class.
  enum EtSumType {
      kTotalEt,
      kTotalHt,
      kMissingEt,
      kMissingHt,
      kTotalEtx,
      kTotalEty,
      kTotalHtx,
      kTotalHty,
      kMissingEtHF,
      kTotalEtxHF,
      kTotalEtyHF,
      kMinBiasHFP0,
      kMinBiasHFM0,
      kMinBiasHFP1,
      kMinBiasHFM1,
      kTotalEtHF,
      kTotalEtEm,
      kTotalHtHF,
      kTotalHtxHF,
      kTotalHtyHF,
      kMissingHtHF,
      kTowerCount      
  };
  
  struct CentralCellSums
  {
    short int cell0,cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8;
    CentralCellSums(){}
    CentralCellSums(const std::array<short int,9>& sums):
        cell0(sums[0]),cell1(sums[1]),cell2(sums[2]),cell3(sums[3]),cell4(sums[4]),
        cell5(sums[5]),cell6(sums[6]),cell7(sums[7]),cell8(sums[8])
    {
    }
  };
  
  struct DonutCellSums
  {
    short int cell0,cell1,cell2,cell3,cell4,cell5,cell6,cell7,cell8,cell9;
    short int cell10,cell11,cell12,cell13,cell14,cell15,cell16,cell17,cell18,cell19;
    short int cell20,cell21,cell22,cell23,cell24,cell25,cell26,cell27,cell28,cell29;
    short int cell30,cell31,cell32,cell33,cell34,cell35,cell36,cell37,cell38,cell39;
    DonutCellSums(){}
    DonutCellSums(const std::array<short int,40>& sums):
        cell0(sums[0]),cell1(sums[1]),cell2(sums[2]),cell3(sums[3]),cell4(sums[4]),
        cell5(sums[5]),cell6(sums[6]),cell7(sums[7]),cell8(sums[8]),cell9(sums[9]),
        cell10(sums[10]),cell11(sums[11]),cell12(sums[12]),cell13(sums[13]),cell14(sums[14]),
        cell15(sums[15]),cell16(sums[16]),cell17(sums[17]),cell18(sums[18]),cell19(sums[19]),
        cell20(sums[20]),cell21(sums[21]),cell22(sums[22]),cell23(sums[23]),cell24(sums[24]),
        cell25(sums[25]),cell26(sums[26]),cell27(sums[27]),cell28(sums[28]),cell29(sums[29]),
        cell30(sums[30]),cell31(sums[31]),cell32(sums[32]),cell33(sums[33]),cell34(sums[34]),
        cell35(sums[35]),cell36(sums[36]),cell37(sums[37]),cell38(sums[38]),cell39(sums[39])
    {
    }
  };
  
  struct L1AnalysisL1UpgradeDataFormat
  {
  
    L1AnalysisL1UpgradeDataFormat(){ Reset();};
    ~L1AnalysisL1UpgradeDataFormat(){};
    
    void Reset()
    {
      nEGs = 0;
      egEt.clear();
      egEta.clear();
      egPhi.clear();
      egIEt.clear();
      egIEta.clear();
      egIPhi.clear();
      egIso.clear();
      egBx.clear();
      egTowerIPhi.clear();
      egTowerIEta.clear();
      egRawEt.clear();
      egIsoEt.clear();
      egFootprintEt.clear();
      egNTT.clear();
      egShape.clear();
      egTowerHoE.clear();
      egHwQual.clear();

      nTaus = 0;
      tauEt.clear();
      tauEta.clear();
      tauPhi.clear(); 
      tauIEt.clear();
      tauIEta.clear();
      tauIPhi.clear(); 
      tauIso.clear();
      tauBx.clear();
      tauTowerIPhi.clear();
      tauTowerIEta.clear();
      tauRawEt.clear();
      tauIsoEt.clear();
      tauNTT.clear();
      tauHasEM.clear();
      tauIsMerged.clear();
      tauHwQual.clear();

      nJets = 0;
      jetEt.clear();
      jetEta.clear();
      jetPhi.clear();
      jetIEt.clear();
      jetIEta.clear();
      jetIPhi.clear();
      jetBx.clear();
      jetTowerIPhi.clear();
      jetTowerIEta.clear();
      jetRawEt.clear();
      jetSeedEt.clear();
      jetPUEt.clear();
      jetPUDonutEt0.clear();
      jetPUDonutEt1.clear();
      jetPUDonutEt2.clear();
      jetPUDonutEt3.clear();

      nMuons = 0;
      muonEt.clear();
      muonEta.clear();
      muonPhi.clear();
      muonEtaAtVtx.clear();
      muonPhiAtVtx.clear();
      muonIEt.clear();
      muonIEta.clear();
      muonIPhi.clear();
      muonIEtaAtVtx.clear();
      muonIPhiAtVtx.clear();
      muonIDEta.clear();
      muonIDPhi.clear();
      muonChg.clear();
      muonIso.clear();
      muonQual.clear();
      muonTfMuonIdx.clear();
      muonBx.clear();
      
      nSums = 0;
      sumType.clear();
      sumEt.clear();
      sumPhi.clear();
      sumIEt.clear();
      sumIPhi.clear();
      sumBx.clear();
        
      jetCentralCellSums.clear();
      jetDonutCellSums.clear();
    }
   
    unsigned short int nEGs;
    std::vector<float> egEt;
    std::vector<float> egEta;
    std::vector<float> egPhi;
    std::vector<short int> egIEt;
    std::vector<short int> egIEta;
    std::vector<short int> egIPhi;
    std::vector<short int> egIso;
    std::vector<short int> egBx;
    std::vector<short int> egTowerIPhi;
    std::vector<short int> egTowerIEta;
    std::vector<short int> egRawEt;
    std::vector<short int> egIsoEt;
    std::vector<short int> egFootprintEt;
    std::vector<short int> egNTT;
    std::vector<short int> egShape;
    std::vector<short int> egTowerHoE;
    std::vector<short int> egHwQual;
 
    unsigned short int nTaus;
    std::vector<float> tauEt;
    std::vector<float> tauEta;
    std::vector<float> tauPhi;
    std::vector<short int> tauIEt;
    std::vector<short int> tauIEta;
    std::vector<short int> tauIPhi;
    std::vector<short int> tauIso;
    std::vector<short int> tauBx;
    std::vector<short int> tauTowerIPhi;
    std::vector<short int> tauTowerIEta;
    std::vector<short int> tauRawEt;    
    std::vector<short int> tauIsoEt;
    std::vector<short int> tauNTT;
    std::vector<short int> tauHasEM;
    std::vector<short int> tauIsMerged;
    std::vector<short int> tauHwQual;

    unsigned short int nJets;
    std::vector<float> jetEt;
    std::vector<float> jetEta;
    std::vector<float> jetPhi;
    std::vector<short int> jetIEt;
    std::vector<short int> jetIEta;
    std::vector<short int> jetIPhi;
    std::vector<short int> jetBx;
    std::vector<short int> jetTowerIPhi;
    std::vector<short int> jetTowerIEta;
    std::vector<short int> jetRawEt;    
    std::vector<short int> jetSeedEt;
    std::vector<short int> jetPUEt;
    std::vector<short int> jetPUDonutEt0;
    std::vector<short int> jetPUDonutEt1;
    std::vector<short int> jetPUDonutEt2;
    std::vector<short int> jetPUDonutEt3;
    
    std::vector<CentralCellSums> jetCentralCellSums;
    std::vector<DonutCellSums> jetDonutCellSums;

    unsigned short int nMuons;
    std::vector<float>   muonEt;
    std::vector<float>   muonEta;
    std::vector<float>   muonPhi;
    std::vector<float>   muonEtaAtVtx;
    std::vector<float>   muonPhiAtVtx;
    std::vector<short int>   muonIEt;
    std::vector<short int>   muonIEta;
    std::vector<short int>   muonIPhi;
    std::vector<short int>   muonIEtaAtVtx;
    std::vector<short int>   muonIPhiAtVtx;
    std::vector<short int>   muonIDEta;
    std::vector<short int>   muonIDPhi;
    std::vector<short int>      muonChg;
    std::vector<unsigned short int> muonIso;
    std::vector<unsigned short int> muonQual;
    std::vector<unsigned short int> muonTfMuonIdx;
    std::vector<short int>      muonBx;

    
    unsigned short int nSums;
    std::vector<short int> sumType;
    std::vector<float> sumEt;
    std::vector<float> sumPhi;
    std::vector<short int> sumIEt;
    std::vector<short int> sumIPhi;
    std::vector<float> sumBx;

  }; 
}
#endif


