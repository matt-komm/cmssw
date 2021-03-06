// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Calibration/HcalCalibAlgos/interface/Analyzer_minbias.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtfeWord.h"
#include "FWCore/Framework/interface/Run.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <fstream>
#include <sstream>

#include "CondFormats/HcalObjects/interface/HcalRespCorrs.h"
#include "CondFormats/DataRecord/interface/HcalRespCorrsRcd.h"

using namespace std;
using namespace reco;
//
// constructors and destructor
//
namespace cms{
Analyzer_minbias::Analyzer_minbias(const edm::ParameterSet& iConfig)
{
  // get name of output file with histogramms
  fOutputFileName = iConfig.getUntrackedParameter<string>("HistOutFile"); 
  // get names of modules, producing object collections
  
  hbherecoMB = iConfig.getParameter<edm::InputTag>("hbheInputMB");
  horecoMB   = iConfig.getParameter<edm::InputTag>("hoInputMB");
  hfrecoMB   = iConfig.getParameter<edm::InputTag>("hfInputMB");
  
  hbherecoNoise = iConfig.getParameter<edm::InputTag>("hbheInputNoise");
  horecoNoise   = iConfig.getParameter<edm::InputTag>("hoInputNoise");
  hfrecoNoise   = iConfig.getParameter<edm::InputTag>("hfInputNoise");
  
  theRecalib = iConfig.getParameter<bool>("Recalib"); 
   
//
//
  for(int i=0; i<73; i++)
  {
     for(int j=0; j<43; j++)
     {
        noise_min[i][j] = 0.;
	noise_pl[i][j] = 0.;
     } 
  }
//
//
     
}

Analyzer_minbias::~Analyzer_minbias()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

void Analyzer_minbias::beginRun( const edm::Run& r, const edm::EventSetup& iSetup)
{
  nevent_run = 0;
}
void Analyzer_minbias::endRun( const edm::Run& r, const edm::EventSetup& iSetup)
{
 std::cout<<" Runnumber "<<r.run()<<" Nevents  "<<nevent_run<<std::endl;
}

void Analyzer_minbias::beginJob()
{
   
   hOutputFile   = new TFile( fOutputFileName.c_str(), "RECREATE" ) ;
   
   myTree = new TTree("RecJet","RecJet Tree");
   myTree->Branch("mydet",  &mydet, "mydet/I");
   myTree->Branch("mysubd",  &mysubd, "mysubd/I");
   myTree->Branch("depth",  &depth, "depth/I");
   myTree->Branch("ieta",  &ieta, "ieta/I");
   myTree->Branch("iphi",  &iphi, "iphi/I");
   myTree->Branch("eta",  &eta, "eta/F");
   myTree->Branch("phi",  &phi, "phi/F");
   
   myTree->Branch("mom0_MB",  &mom0_MB, "mom0_MB/F");
   myTree->Branch("mom1_MB",  &mom1_MB, "mom1_MB/F");
   myTree->Branch("mom2_MB",  &mom2_MB, "mom2_MB/F");
   myTree->Branch("mom4_MB",  &mom4_MB, "mom4_MB/F");
   
   myTree->Branch("mom0_Noise",  &mom0_Noise, "mom0_Noise/F");
   myTree->Branch("mom1_Noise",  &mom1_Noise, "mom1_Noise/F");
   myTree->Branch("mom2_Noise",  &mom2_Noise, "mom2_Noise/F");
   myTree->Branch("mom4_Noise",  &mom4_Noise, "mom4_Noise/F");
   
   myTree->Branch("mom0_Diff",  &mom0_Diff, "mom0_Diff/F");
   myTree->Branch("mom1_Diff",  &mom1_Diff, "mom1_Diff/F");
   myTree->Branch("mom2_Diff",  &mom2_Diff, "mom2_Diff/F");

   myTree->Branch("occup",  &occup, "occup/F");
   
   std::cout<<" Before ordering Histos "<<std::endl;
  
   char str0[15];
   char str1[15];

   char str10[15];
   char str11[15];

   int k=0;
   nevent = 0;
// Size of collections

   hHBHEsize_vs_run = new TH2F("hHBHEsize_vs_run","hHBHEsize_vs_run",500,111500.,112000.,6101,-100.5,6000.5);
   hHFsize_vs_run = new TH2F("hHFsize_vs_run","hHFsize_vs_run",500,111500.,112000.,6101,-100.5,6000.5); 

   for(int i=1;i<73;i++){
    for(int j=1;j<43;j++){

       meannoise_pl[i][j] = 0.;
       meannoise_min[i][j] = 0.;

//     for(int l=1;l<5;l++){
        k = i*1000+j;
        sprintf(str0,"mpl%d",k);
        sprintf(str1,"mmin%d",k);

        sprintf(str10,"vpl%d",k);
        sprintf(str11,"vmin%d",k);
//      cout<<" "<<i<<" "<<j<<endl;
    if( j < 30 )
    {
// first order moment
    hCalo1[i][j] = new TH1F(str0, "h0", 320, -10., 10.);
    hCalo2[i][j] = new TH1F(str1, "h1", 320, -10., 10.);

// second order moment
    hCalo1mom2[i][j] = new TH1F(str10, "h10", 320, 0., 20.);
    hCalo2mom2[i][j] = new TH1F(str11, "h11", 320, 0., 20.);
    }
      else
      {
// HF
// first order moment
//   cout<<" "<<i<<" "<<j<<" "<<k<<endl;
   if(j < 40)
   {
    hCalo1[i][j] = new TH1F(str0, "h0", 320, -10., 10.);
    hCalo2[i][j] = new TH1F(str1, "h1", 320, -10., 10.);
//
// second order moment
    hCalo1mom2[i][j] = new TH1F(str10, "h10", 320, 0., 40.);
    hCalo2mom2[i][j] = new TH1F(str11, "h11", 320, 0., 40.);
   }
     else
     {
    hCalo1[i][j] = new TH1F(str0,"h0" , 320, -10., 10.);
    hCalo2[i][j] = new TH1F(str1, "h1", 320, -10., 10.);

// second order moment
    hCalo1mom2[i][j] = new TH1F(str10, "h10", 320, 0., 120.);
    hCalo2mom2[i][j] = new TH1F(str11, "h11", 320, 0., 120.);

     }
    } // HE/HF boundary
//     } // l
    } // j
   } // i


   std::cout<<" After ordering Histos "<<std::endl;

   std::string ccc = "noise_0.dat";

   myout_hcal = new ofstream(ccc.c_str());
   if(!myout_hcal) cout << " Output file not open!!! "<<endl;

//
   for (int i=0; i<5;i++)
   {
    for (int j=0; j<5;j++)
    {
     for (int k=0; k<73;k++)
     {
       for (int l=0; l<43;l++)
       {
        theMBFillDetMapPl0[i][j][k][l] = 0.;
        theMBFillDetMapPl1[i][j][k][l] = 0.;
        theMBFillDetMapPl2[i][j][k][l] = 0.;
        theMBFillDetMapPl4[i][j][k][l] = 0.;

        theMBFillDetMapMin0[i][j][k][l] = 0.;
        theMBFillDetMapMin1[i][j][k][l] = 0.;
        theMBFillDetMapMin2[i][j][k][l] = 0.;
        theMBFillDetMapMin4[i][j][k][l] = 0.;


        theNSFillDetMapPl0[i][j][k][l] = 0.;
        theNSFillDetMapPl1[i][j][k][l] = 0.;
        theNSFillDetMapPl2[i][j][k][l] = 0.;
        theNSFillDetMapPl4[i][j][k][l] = 0.;

        theNSFillDetMapMin0[i][j][k][l] = 0.;
        theNSFillDetMapMin1[i][j][k][l] = 0.;
        theNSFillDetMapMin2[i][j][k][l] = 0.;
        theNSFillDetMapMin4[i][j][k][l] = 0.;

        theDFFillDetMapPl0[i][j][k][l] = 0.;
        theDFFillDetMapPl1[i][j][k][l] = 0.;
        theDFFillDetMapPl2[i][j][k][l] = 0.;
        theDFFillDetMapMin0[i][j][k][l] = 0.;
        theDFFillDetMapMin1[i][j][k][l] = 0.;
        theDFFillDetMapMin2[i][j][k][l] = 0.;
       }
     }  
    }
   }    
    
   return ;
}
//
//  EndJob
//
void Analyzer_minbias::endJob()
{
   int ii=0;
      
   for (int i=1; i<5;i++)
   {
    for (int j=1; j<5;j++)
    {
     for (int k=1; k<73;k++)
     {
       for (int l=1; l<43;l++)
       {
	  if(theMBFillDetMapPl0[i][j][k][l] > 0)
	  { 
            mom0_MB = theMBFillDetMapPl0[i][j][k][l];
            mom1_MB = theMBFillDetMapPl1[i][j][k][l];
            mom2_MB = theMBFillDetMapPl2[i][j][k][l];
            mom4_MB = theMBFillDetMapPl4[i][j][k][l];
            mom0_Noise = theNSFillDetMapPl0[i][j][k][l];
            mom1_Noise = theNSFillDetMapPl1[i][j][k][l];
            mom2_Noise = theNSFillDetMapPl2[i][j][k][l];
            mom4_Noise = theNSFillDetMapPl4[i][j][k][l];
            mom0_Diff = theDFFillDetMapPl0[i][j][k][l];
            mom1_Diff = theDFFillDetMapPl1[i][j][k][l];
            mom2_Diff = theDFFillDetMapPl2[i][j][k][l];
	    
            mysubd = i;
            depth = j;
            ieta = l;
            iphi = k;
            cout<<" Result Plus= "<<mysubd<<" "<<ieta<<" "<<iphi<<" mom0  "<<mom0_MB<<" mom1 "<<mom1_MB<<" mom2 "<<mom2_MB<<endl;
            myTree->Fill();
            ii++;
	   } // Pl > 0
	
         
	  if(theMBFillDetMapMin0[i][j][k][l] > 0)
	  { 
            mom0_MB = theMBFillDetMapMin0[i][j][k][l];
            mom1_MB = theMBFillDetMapMin1[i][j][k][l];
            mom2_MB = theMBFillDetMapMin2[i][j][k][l];
            mom4_MB = theMBFillDetMapMin4[i][j][k][l];
            mom0_Noise = theNSFillDetMapMin0[i][j][k][l];
            mom1_Noise = theNSFillDetMapMin1[i][j][k][l];
            mom2_Noise = theNSFillDetMapMin2[i][j][k][l];
            mom4_Noise = theNSFillDetMapMin4[i][j][k][l];
            mom0_Diff = theDFFillDetMapMin0[i][j][k][l];
            mom1_Diff = theDFFillDetMapMin1[i][j][k][l];
            mom2_Diff = theDFFillDetMapMin2[i][j][k][l];
	    
	    
            mysubd = i;
            depth = j;
            ieta = -1*l;
            iphi = k;
            cout<<" Result Minus= "<<mysubd<<" "<<ieta<<" "<<iphi<<" mom0  "<<mom0_MB<<" mom1 "<<mom1_MB<<" mom2 "<<mom2_MB<<endl;
            myTree->Fill();
            ii++;
	    
	  } // Min>0  
      } // ieta
     } // iphi  
    } // depth
   } //subd    
   
   
   
   cout<<" Number of cells "<<ii<<endl; 
      
   hOutputFile->Write() ;   
   hOutputFile->cd();
   hHBHEsize_vs_run->Write() ;
   hHFsize_vs_run->Write() ;
   for(int i=1;i<73;i++){
    for(int j=1;j<43;j++){
    hCalo1[i][j]->Write();
    hCalo2[i][j]->Write();
    hCalo1mom2[i][j]->Write();
    hCalo2mom2[i][j]->Write();
   }
  }
  
   
   myTree->Write();
   hOutputFile->Close() ;
   
   return ;
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
Analyzer_minbias::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout<<" Start Analyzer_minbias::analyze "<<nevent<<std::endl;
  nevent++; 
  nevent_run++;
  using namespace edm;

  float rnnum = (float)iEvent.run(); 
/*
  std::vector<Provenance const*> theProvenance;
  iEvent.getAllProvenance(theProvenance);

  for( std::vector<Provenance const*>::const_iterator ip = theProvenance.begin();
                                                      ip != theProvenance.end(); ip++)
  {
     cout<<" Print all module/label names "<<(**ip).moduleName()<<" "<<(**ip).moduleLabel()<<
     " "<<(**ip).productInstanceName()<<endl;
  }
*/

   edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
   iEvent.getByLabel("gtDigis", gtRecord);
   
   if (!gtRecord.isValid()) {

//     LogDebug("L1GlobalTriggerRecordProducer")
//       << "\n\n Error: no L1GlobalTriggerReadoutRecord found with input tag "
//       << m_l1GtReadoutRecord
//       << "\n Returning empty L1GlobalTriggerRecord.\n\n"
//       << std::endl;
      cout<<" No L1 trigger record "<<endl;
//     return;
   }


  const HcalRespCorrs* myRecalib=0;
  if( theRecalib ) {
// Radek:   
  edm::ESHandle <HcalRespCorrs> recalibCorrs;
  iSetup.get<HcalRespCorrsRcd>().get("recalibrate",recalibCorrs);
  myRecalib = recalibCorrs.product();
// end
  } // theRecalib

// Noise part for HB HE

     double tmpNSFillDetMapPl1[5][5][73][43]; 
     double tmpNSFillDetMapPl2[5][5][73][43];
     double tmpNSFillDetMapMin1[5][5][73][43]; 
     double tmpNSFillDetMapMin2[5][5][73][43];

   for (int i=0; i<5;i++)
   {
    for (int j=0; j<5;j++)
    {
     for (int k=0; k<73;k++)
     {
       for (int l=0; l<43;l++)
       {
        tmpNSFillDetMapPl1[i][j][k][l] = 0.;
        tmpNSFillDetMapPl2[i][j][k][l] = 0.;
        tmpNSFillDetMapMin1[i][j][k][l] = 0.;
        tmpNSFillDetMapMin2[i][j][k][l] = 0.;
       }
     }  
    }
   }    
   edm::Handle<HBHERecHitCollection> hbheNormal;
   iEvent.getByLabel("hbhereco", hbheNormal);
   if(!hbheNormal.isValid()){  
     cout<<" hbheNormal failed "<<endl;
   } else {
       cout<<" The size of the normal collection "<<hbheNormal->size()<<endl;
   } 


   edm::Handle<HBHERecHitCollection> hbheNS;
   iEvent.getByLabel(hbherecoNoise, hbheNS);


   if(!hbheNS.isValid()){
     LogDebug("") << "HcalCalibAlgos: Error! can't get hbhe product!" << std::endl;
     cout<<" No HBHE MS "<<endl;
     return ;
   }

  
   const HBHERecHitCollection HithbheNS = *(hbheNS.product());
   cout<<" HBHE NS size of collection "<<HithbheNS.size()<<endl;
   hHBHEsize_vs_run->Fill(rnnum,(float)HithbheNS.size());

   if(HithbheNS.size()!= 5184) {
          cout<<" HBHE problem "<<rnnum<<" "<<HithbheNS.size()<<endl;
          return;
   }
   edm::Handle<HBHERecHitCollection> hbheMB;
   iEvent.getByLabel(hbherecoMB, hbheMB);

   if(!hbheMB.isValid()){
     LogDebug("") << "HcalCalibAlgos: Error! can't get hbhe product!" << std::endl;
     cout<<" No HBHE MB"<<endl;
     return ;
   }

  const HBHERecHitCollection HithbheMB = *(hbheMB.product());
  cout<<" HBHE MB size of collection "<<HithbheMB.size()<<endl;
  if(HithbheMB.size()!= 5184) {
     cout<<" HBHE problem "<<rnnum<<" "<<HithbheMB.size()<<endl;
     return;
  }

   edm::Handle<HFRecHitCollection> hfNS;
   iEvent.getByLabel(hfrecoNoise, hfNS);

   if(!hfNS.isValid()){
     LogDebug("") << "HcalCalibAlgos: Error! can't get hbhe product!" << std::endl;
     cout<<" No HF NS "<<endl;
     return ;
   }

   const HFRecHitCollection HithfNS = *(hfNS.product());
   cout<<" HFE NS size of collection "<<HithfNS.size()<<endl;
   hHFsize_vs_run->Fill(rnnum,(float)HithfNS.size());
   if(HithfNS.size()!= 1728) {
          cout<<" HF problem "<<rnnum<<" "<<HithfNS.size()<<endl;
          return;
   }

   edm::Handle<HFRecHitCollection> hfMB;
   iEvent.getByLabel(hfrecoMB, hfMB);

   if(!hfMB.isValid()){
     LogDebug("") << "HcalCalibAlgos: Error! can't get hbhe product!" << std::endl;
     cout<<" No HBHE MB"<<endl;
     return ;
   }

  const HFRecHitCollection HithfMB = *(hfMB.product());
  cout<<" HF MB size of collection "<<HithfMB.size()<<endl;
   if(HithfMB.size()!= 1728) {
          cout<<" HF problem "<<rnnum<<" "<<HithfMB.size()<<endl;
          return;
   }



  for(HBHERecHitCollection::const_iterator hbheItr=HithbheNS.begin(); hbheItr!=HithbheNS.end(); hbheItr++)
        {
// Recalibration of energy
         float icalconst=1.;	 
         DetId mydetid = hbheItr->id().rawId();
	 if( theRecalib ) icalconst=myRecalib->getValues(mydetid)->getValue();
    
	 HBHERecHit aHit(hbheItr->id(),hbheItr->energy()*icalconst,hbheItr->time());
	 
         double energyhit = aHit.energy();
	 
	 DetId id = (*hbheItr).detid(); 
	 HcalDetId hid=HcalDetId(id);
 
         int mysu = ((hid).rawId()>>25)&0x7;
         if( hid.ieta() > 0 ) {
	 theNSFillDetMapPl0[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theNSFillDetMapPl0[mysu][hid.depth()][hid.iphi()][hid.ieta()]+ 1.;
         theNSFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theNSFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()]+energyhit;
	 theNSFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theNSFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow(energyhit,2);
         theNSFillDetMapPl4[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theNSFillDetMapPl4[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow(energyhit,4);
	 
         tmpNSFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = energyhit;
	 tmpNSFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()] = pow(energyhit,2);
	 
	 
         } else {
	 theNSFillDetMapMin0[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theNSFillDetMapMin0[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+ 1.;
         theNSFillDetMapMin1[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theNSFillDetMapMin1[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+energyhit;
	 theNSFillDetMapMin2[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theNSFillDetMapMin2[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+pow(energyhit,2);
         theNSFillDetMapMin4[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theNSFillDetMapMin4[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+pow(energyhit,4);


         tmpNSFillDetMapMin1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = energyhit;
	 tmpNSFillDetMapMin2[mysu][hid.depth()][hid.iphi()][hid.ieta()] = pow(energyhit,2);
	 
	 }  
	 
         if(hid.depth() == 1) {
         if( hid.ieta() > 0 ) {
          hCalo1[hid.iphi()][hid.ieta()]->Fill(energyhit-noise_pl[hid.iphi()][hid.ieta()]);
          hCalo1mom2[hid.iphi()][hid.ieta()]->Fill(pow(energyhit,2));
         } else {
          hCalo2[hid.iphi()][abs(hid.ieta())]->Fill(energyhit-noise_min[hid.iphi()][abs(hid.ieta())]);
          hCalo2mom2[hid.iphi()][abs(hid.ieta())]->Fill(pow(energyhit,2));
         } // eta><0
	 } // depth=1
        } // HBHE_NS


// Signal part for HB HE
     

  for(HBHERecHitCollection::const_iterator hbheItr=HithbheMB.begin(); hbheItr!=HithbheMB.end(); hbheItr++)
        {
// Recalibration of energy
         float icalconst=1.;	 
         DetId mydetid = hbheItr->id().rawId();
	 if( theRecalib ) icalconst=myRecalib->getValues(mydetid)->getValue();
    
	 HBHERecHit aHit(hbheItr->id(),hbheItr->energy()*icalconst,hbheItr->time());
	 
         double energyhit = aHit.energy();
	 
	 DetId id = (*hbheItr).detid(); 
	 HcalDetId hid=HcalDetId(id);
 
         int mysu = ((hid).rawId()>>25)&0x7;
         if( hid.ieta() > 0 ) {
	 theMBFillDetMapPl0[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theMBFillDetMapPl0[mysu][hid.depth()][hid.iphi()][hid.ieta()]+ 1.;
         theMBFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theMBFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()]+energyhit;
	 theMBFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theMBFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow(energyhit,2);
         theMBFillDetMapPl4[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theMBFillDetMapPl4[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow(energyhit,4);
	 float mydiff = energyhit - tmpNSFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()];
	 
	 
         theDFFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theDFFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()]+mydiff;
	 theDFFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theDFFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow(mydiff,2);
         } else {
	 theMBFillDetMapMin0[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theMBFillDetMapMin0[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+ 1.;
         theMBFillDetMapMin1[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theMBFillDetMapMin1[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+energyhit;
	 theMBFillDetMapMin2[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theMBFillDetMapMin2[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+pow(energyhit,2);
         theMBFillDetMapMin4[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theMBFillDetMapMin4[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+pow(energyhit,4);


	 float mydiff = energyhit - tmpNSFillDetMapMin1[mysu][hid.depth()][hid.iphi()][hid.ieta()];
         theDFFillDetMapMin1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theDFFillDetMapMin1[mysu][hid.depth()][hid.iphi()][hid.ieta()]+mydiff;
	 theDFFillDetMapMin2[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theDFFillDetMapMin2[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow(mydiff,2);
	 }  

	 
         if(hid.depth() == 1) {
         if( hid.ieta() > 0 ) {
          hCalo1[hid.iphi()][hid.ieta()]->Fill(energyhit);
          hCalo1mom2[hid.iphi()][hid.ieta()]->Fill(pow(energyhit,2));
         } else {
          hCalo2[hid.iphi()][abs(hid.ieta())]->Fill(energyhit);
          hCalo2mom2[hid.iphi()][abs(hid.ieta())]->Fill(pow(energyhit,2));
         } // eta><0
	 } // depth=1
        } // HBHE_MB
	
// HF
 
  for(HFRecHitCollection::const_iterator hbheItr=HithfNS.begin(); hbheItr!=HithfNS.end(); hbheItr++)
        {
// Recalibration of energy
         float icalconst=1.;	 
         DetId mydetid = hbheItr->id().rawId();
	 if( theRecalib ) icalconst=myRecalib->getValues(mydetid)->getValue();
    
	 HFRecHit aHit(hbheItr->id(),hbheItr->energy()*icalconst,hbheItr->time());
	 
         double energyhit = aHit.energy();
	 
	 DetId id = (*hbheItr).detid(); 
	 HcalDetId hid=HcalDetId(id);
 
         int mysu = ((hid).rawId()>>25)&0x7;
         if( hid.ieta() > 0 ) {
	 theNSFillDetMapPl0[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theNSFillDetMapPl0[mysu][hid.depth()][hid.iphi()][hid.ieta()]+ 1.;
         theNSFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theNSFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()]+energyhit;
	 theNSFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theNSFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow(energyhit,2);
         theNSFillDetMapPl4[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theNSFillDetMapPl4[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow(energyhit,4);
	 
         tmpNSFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = energyhit;
	 tmpNSFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()] = pow(energyhit,2);
	 
	 
         } else {
	 theNSFillDetMapMin0[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theNSFillDetMapMin0[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+ 1.;
         theNSFillDetMapMin1[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theNSFillDetMapMin1[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+energyhit;
	 theNSFillDetMapMin2[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theNSFillDetMapMin2[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+pow(energyhit,2);
         theNSFillDetMapMin4[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theNSFillDetMapMin4[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+pow(energyhit,4);

 
         tmpNSFillDetMapMin1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = energyhit;
	 tmpNSFillDetMapMin2[mysu][hid.depth()][hid.iphi()][hid.ieta()] = pow(energyhit,2);
	 
	 }  
	 
         if(hid.depth() == 1) {
         if( hid.ieta() > 0 ) {
          hCalo1[hid.iphi()][hid.ieta()]->Fill(energyhit-noise_pl[hid.iphi()][hid.ieta()]);
          hCalo1mom2[hid.iphi()][hid.ieta()]->Fill(pow(energyhit,2));
         } else {
          hCalo2[hid.iphi()][abs(hid.ieta())]->Fill(energyhit-noise_min[hid.iphi()][abs(hid.ieta())]);
          hCalo2mom2[hid.iphi()][abs(hid.ieta())]->Fill(pow(energyhit,2));
         } // eta><0
	 } // depth=1
        } // HBHE_NS


// Signal part for HB HE

  for(HFRecHitCollection::const_iterator hbheItr=HithfMB.begin(); hbheItr!=HithfMB.end(); hbheItr++)
        {
// Recalibration of energy
         float icalconst=1.;	 
         DetId mydetid = hbheItr->id().rawId();
	 if( theRecalib ) icalconst=myRecalib->getValues(mydetid)->getValue();
    
	 HFRecHit aHit(hbheItr->id(),hbheItr->energy()*icalconst,hbheItr->time());
	 
         double energyhit = aHit.energy();
	 
	 DetId id = (*hbheItr).detid(); 
	 HcalDetId hid=HcalDetId(id);
 
         int mysu = ((hid).rawId()>>25)&0x7;
         if( hid.ieta() > 0 ) {
	 theMBFillDetMapPl0[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theMBFillDetMapPl0[mysu][hid.depth()][hid.iphi()][hid.ieta()]+ 1.;
         theMBFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theMBFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()]+energyhit;
	 theMBFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theMBFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow(energyhit,2);
         theMBFillDetMapPl4[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theMBFillDetMapPl4[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow(energyhit,4);


         theDFFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theDFFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()]+energyhit-tmpNSFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()];
	 theDFFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()] =
	 theDFFillDetMapPl2[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow((energyhit-tmpNSFillDetMapPl1[mysu][hid.depth()][hid.iphi()][hid.ieta()]),2);
         } else {
	 theMBFillDetMapMin0[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theMBFillDetMapMin0[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+ 1.;
         theMBFillDetMapMin1[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theMBFillDetMapMin1[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+energyhit;
	 theMBFillDetMapMin2[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theMBFillDetMapMin2[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+pow(energyhit,2);
         theMBFillDetMapMin4[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())] = theMBFillDetMapMin4[mysu][hid.depth()][hid.iphi()][abs(hid.ieta())]+pow(energyhit,4);
 
         theDFFillDetMapMin1[mysu][hid.depth()][hid.iphi()][hid.ieta()] = theDFFillDetMapMin1[mysu][hid.depth()][hid.iphi()][hid.ieta()]+energyhit-tmpNSFillDetMapMin1[mysu][hid.depth()][hid.iphi()][hid.ieta()];
	 theDFFillDetMapMin2[mysu][hid.depth()][hid.iphi()][hid.ieta()] =
	 theDFFillDetMapMin2[mysu][hid.depth()][hid.iphi()][hid.ieta()]+pow((energyhit-tmpNSFillDetMapMin1[mysu][hid.depth()][hid.iphi()][hid.ieta()]),2);
	 }  

	 
         if(hid.depth() == 1) {
         if( hid.ieta() > 0 ) {
          hCalo1[hid.iphi()][hid.ieta()]->Fill(energyhit);
          hCalo1mom2[hid.iphi()][hid.ieta()]->Fill(pow(energyhit,2));
         } else {
          hCalo2[hid.iphi()][abs(hid.ieta())]->Fill(energyhit);
          hCalo2mom2[hid.iphi()][abs(hid.ieta())]->Fill(pow(energyhit,2));
         } // eta><0
	 } // depth=1
        } // HF_MB

   std::cout<<" Event is finished "<<std::endl;
}
}
//define this as a plug-in
//DEFINE_FWK_MODULE(Analyzer_minbias)

