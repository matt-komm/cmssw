
/*
    This is the DQM code for UE physics plots
    11/12/2009 Sunil Bansal
*/
#include "DQM/Physics/src/QcdUeDQM.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "CommonTools/RecoAlgos/src/TrackToRefCandidate.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfo.h"
#include "AnalysisDataFormats/TrackInfo/interface/TrackInfoTrackAssociation.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include <TString.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
using namespace std;
using namespace edm;

#define CP(level) \
  if (level>=verbose_)

struct deleter {
  void operator()(TH3F *&h) { delete h; h=0;}
};


//--------------------------------------------------------------------------------------------------
QcdUeDQM::QcdUeDQM(const ParameterSet &parameters) :
  hltResName_(parameters.getUntrackedParameter<string>("hltTrgResults")),
  verbose_(parameters.getUntrackedParameter<int>("verbose",3)),
  tgeo_(0),
  theDbe_(0),
  repSumMap_(0),
  repSummary_(0),
  h2TrigCorr_(0),
  ptMin_(parameters.getParameter<double>("ptMin")),
  minRapidity_(parameters.getParameter<double>("minRapidity")),
  maxRapidity_(parameters.getParameter<double>("maxRapidity")),
  tip_(parameters.getParameter<double>("tip")),
  lip_(parameters.getParameter<double>("lip")),
  diffvtxbs_(parameters.getParameter<double>("diffvtxbs")),
  ptErr_pt_(parameters.getParameter<double>("ptErr_pt")),
  vtxntk_(parameters.getParameter<double>("vtxntk")),
  minHit_(parameters.getParameter<int>("minHit")),
  pxlLayerMinCut_(parameters.getParameter<double>("pxlLayerMinCut")),
  requirePIX1_(parameters.getParameter<bool>("requirePIX1")),
  min3DHit_(parameters.getParameter<int>("min3DHit")),
  maxChi2_(parameters.getParameter<double>("maxChi2")),
  bsuse_(parameters.getParameter<bool>("bsuse")),
  allowTriplets_(parameters.getParameter<bool>("allowTriplets")),
  bsPos_(parameters.getParameter<double>("bsPos")),
  caloJetLabel_(parameters.getUntrackedParameter<edm::InputTag>("caloJetTag")),
  chargedJetLabel_(parameters.getUntrackedParameter<edm::InputTag>("chargedJetTag")),
  trackLabel_(parameters.getUntrackedParameter<edm::InputTag>("trackTag")),
  vtxLabel_(parameters.getUntrackedParameter<edm::InputTag>("vtxTag")),
  bsLabel_(parameters.getParameter<edm::InputTag>("beamSpotTag")) 
{
  // Constructor.
  std::vector<std::string> quality = parameters.getParameter<std::vector<std::string> >("quality");
  for (unsigned int j=0;j<quality.size();j++) quality_.push_back(reco::TrackBase::qualityByName(quality[j])); 
  std::vector<std::string> algorithm = parameters.getParameter<std::vector<std::string> >("algorithm");
  for (unsigned int j=0;j<algorithm.size();j++) algorithm_.push_back(reco::TrackBase::algoByName(algorithm[j])); 

  if (parameters.exists("hltTrgNames"))
    hltTrgNames_ = parameters.getUntrackedParameter<vector<string> >("hltTrgNames");

  if (parameters.exists("hltProcNames"))
     hltProcNames_ = parameters.getUntrackedParameter<vector<string> >("hltProcNames");
  else {
    //     hltProcNames_.push_back("FU");
     hltProcNames_.push_back("HLT");
  }

 
}

//--------------------------------------------------------------------------------------------------
QcdUeDQM::~QcdUeDQM()
{
  // Destructor.


}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::analyze(const Event &iEvent, const EventSetup &iSetup) 
{
  if( ! isHltConfigSuccessful_ ) return;

  // Analyze the given event.
   
   edm::Handle<reco::BeamSpot> beamSpot;
   bool ValidBS_ = iEvent.getByLabel(bsLabel_,beamSpot);
   if(!ValidBS_)return;

   edm::Handle<reco::TrackCollection>tracks ;
   bool ValidTrack_ = iEvent.getByLabel(trackLabel_,tracks);
   if(!ValidTrack_)return;

   edm::Handle<reco::CandidateView> trkJets;
   bool ValidTrackJet_ = iEvent.getByLabel (chargedJetLabel_,trkJets);
   if(!ValidTrackJet_)return;
   
   edm::Handle<reco::CaloJetCollection> calJets;
   bool ValidCaloJet_ = iEvent.getByLabel (caloJetLabel_,calJets);
   if(!ValidCaloJet_)return;
 
   edm::Handle< reco::VertexCollection > vertexColl;
   bool ValidVtxColl_ = iEvent.getByLabel (vtxLabel_, vertexColl);
   if(!ValidVtxColl_)return;

   reco::TrackCollection tracks_sort = *tracks;
   std::sort(tracks_sort.begin(), tracks_sort.end(), PtSorter()); 

  // get tracker geometry
  ESHandle<TrackerGeometry> trackerHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle);
  tgeo_ = trackerHandle.product();
  if (!tgeo_)return;
  selected_.clear(); 
  fillHltBits(iEvent);
  // select good tracks
  if(fillVtxPlots(vertexColl))
  {
  fill1D(hNevts_,1);
  for(reco::TrackCollection::const_iterator Trk = tracks_sort.begin(); Trk != tracks_sort.end(); ++Trk)
   {
    
   if ( trackSelection(*Trk,beamSpot.product(),vtx1,vertexColl->size()) ) selected_.push_back( & * Trk );   
   }
 
      
    fillpTMaxRelated(selected_);
    fillChargedJetSpectra(trkJets);  
    fillCaloJetSpectra(calJets);
    fillUE_with_MaxpTtrack(selected_);
    if(trkJets->size() > 0)fillUE_with_ChargedJets(selected_,trkJets); 
    if(calJets->size()>0)fillUE_with_CaloJets(selected_,calJets);
  
 }

}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::beginJob() 
{
  // Begin job and setup the DQM store.

  theDbe_ = Service<DQMStore>().operator->();
  if (!theDbe_)return;
  
  //  theDbe_->setCurrentFolder("Physics/QcdUe");
  
}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::beginLuminosityBlock(const LuminosityBlock &l, 
                                       const EventSetup &iSetup)
{
  if( ! isHltConfigSuccessful_ ) return;

  // At the moment, nothing needed to be done.
}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::beginRun(const Run &run, const EventSetup &iSetup)
{

 // indicating change of HLT cfg at run boundries
 // for HLTConfigProvider::init()
 bool isHltCfgChange = false;
 isHltConfigSuccessful_ = false; // init

 string teststr;
 for(size_t i=0; i<hltProcNames_.size(); ++i) {
   if (i>0) 
     teststr += ", ";
   teststr += hltProcNames_.at(i);
   if ( hltConfig.init( run, iSetup, hltProcNames_.at(i), isHltCfgChange ) ) {
     isHltConfigSuccessful_ = true;
     hltUsedResName_ = hltResName_;
     if (hltResName_.find(':')==string::npos)
       hltUsedResName_ += "::";
     else 
       hltUsedResName_ += ":";
     hltUsedResName_ += hltProcNames_.at(i);
     break;
   }
 }
 
 if ( ! isHltConfigSuccessful_ )return;

  // setup "Any" bit
  hltTrgBits_.clear();
  hltTrgBits_.push_back(-1);
  hltTrgDeci_.clear();
  hltTrgDeci_.push_back(true);
  hltTrgUsedNames_.clear();
  hltTrgUsedNames_.push_back("Any");

  // figure out relation of trigger name to trigger bit and store used trigger names/bits
  for(size_t i=0;i<hltTrgNames_.size();++i) {
    const string &n1(hltTrgNames_.at(i));
    bool found = 0;
    for(size_t j=0;j<hltConfig.size();++j) {
      const string &n2(hltConfig.triggerName(j));
      if (n2==n1) {
        hltTrgBits_.push_back(j);
        hltTrgUsedNames_.push_back(n1);
        hltTrgDeci_.push_back(false);
        found = 1;
        break;
      }
    }      
    if (!found) {
      CP(2) cout<<"Could not find trigger bit"<<endl ;
    }
  }
 
  // book monitoring histograms
  createHistos();
  isHltConfigSuccessful_ = true;

}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::book1D(std::vector<MonitorElement*> &mes, 
                         const std::string &name, const std::string &title, 
                         int nx, double x1, double x2, bool sumw2, bool sbox)
{
  // Book 1D histos.

  for(size_t i=0;i<hltTrgUsedNames_.size();++i) {
    std::string folderName = "Physics/QcdUe/" + hltTrgUsedNames_.at(i);
    theDbe_->setCurrentFolder(folderName);
    MonitorElement *e = theDbe_->book1D(Form("%s_%s",name.c_str(),hltTrgUsedNames_.at(i).c_str()),
                                        Form("%s: %s",hltTrgUsedNames_.at(i).c_str(), title.c_str()), 
                                        nx, x1, x2);
    TH1 *h1 = e->getTH1();
    if (sumw2) {
      if( 0 == h1->GetSumw2N() ) { // protect against re-summing (would cause exception)
	h1->Sumw2();
      }
    }
    h1->SetStats(sbox);
    mes.push_back(e);
  }
}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::book2D(std::vector<MonitorElement*> &mes, 
                         const std::string &name, const std::string &title, 
                         int nx, double x1, double x2, int ny, double y1, double y2, 
                         bool sumw2, bool sbox)
{
  // Book 2D histos.

  for(size_t i=0;i<hltTrgUsedNames_.size();++i) {
    std::string folderName = "Physics/QcdUe/" + hltTrgUsedNames_.at(i);
    theDbe_->setCurrentFolder(folderName);
    MonitorElement *e = theDbe_->book2D(Form("%s_%s",name.c_str(),hltTrgUsedNames_.at(i).c_str()),
                                        Form("%s: %s",hltTrgUsedNames_.at(i).c_str(), title.c_str()), 
                                        nx, x1, x2, ny, y1, y2);
    TH1 *h1 = e->getTH1();
    if (sumw2) {
      if( 0 == h1->GetSumw2N() ) { // protect against re-summing (would cause exception)
	h1->Sumw2();
      }
    }
    h1->SetStats(sbox);
    mes.push_back(e);
  }
}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::bookProfile(std::vector<MonitorElement*> &mes, 
                         const std::string &name, const std::string &title, 
                         int nx, double x1, double x2,  double y1, double y2, 
                         bool sumw2, bool sbox)
{
  // Book Profile histos.

  for(size_t i=0;i<hltTrgUsedNames_.size();++i) {
    std::string folderName = "Physics/QcdUe/" + hltTrgUsedNames_.at(i);
    theDbe_->setCurrentFolder(folderName);
    MonitorElement *e = theDbe_->bookProfile(Form("%s_%s",name.c_str(),hltTrgUsedNames_.at(i).c_str()),
                                        Form("%s: %s",hltTrgUsedNames_.at(i).c_str(), title.c_str()), 
                                        nx, x1, x2, y1, y2," ");
    mes.push_back(e);
  }
}
//--------------------------------------------------------------------------------------------------
void QcdUeDQM::createHistos()
{
  // Book histograms if needed.


  if (1) {
    theDbe_->setCurrentFolder("Physics/EventInfo/");
    repSumMap_  = theDbe_->book2D("reportSummaryMap","reportSummaryMap",1,0,1,1,0,1);
    repSummary_ = theDbe_->bookFloat("reportSummary");
  }
   
   theDbe_->setCurrentFolder("Physics/QcdUe");

  if (1) {
    const int Nx = hltTrgUsedNames_.size();
    const double x1 = -0.5;
    const double x2 = Nx-0.5;
    h2TrigCorr_ = theDbe_->book2D("h2TriCorr","Trigger bit x vs y;y&&!x;x&&y",Nx,x1,x2,Nx,x1,x2);
    for(size_t i=1;i<=hltTrgUsedNames_.size();++i) {
      h2TrigCorr_->setBinLabel(i,hltTrgUsedNames_.at(i-1),1);
      h2TrigCorr_->setBinLabel(i,hltTrgUsedNames_.at(i-1),2);
    }
    TH1 *h = h2TrigCorr_->getTH1();
    if (h)
      h->SetStats(0);
  }
  book1D(hNevts_,"hNevts","number of events",2,0,2);
  book1D(hNtrackerLayer_,"hNtrackerLayer","number of tracker layers;multiplicity",20,-0.5,19.5 );
  book1D(hNtrackerPixelLayer_,"hNtrackerPixelLayer","number of pixel layers;multiplicity",10,-0.5,9.5 );
  book1D(hNtrackerStripPixelLayer_,"hNtrackerStripPixelLayer","number of strip + pixel layers;multiplicity",30,-0.5,39.5 );
  book1D(hRatioPtErrorPt_,"hRatioPtErrorPt","ratio of pT error and track pT",25,0.,5.);
  book1D(hTrkPt_,"hTrkPt","pT of all tracks",50,0.,50.);
  book1D(hTrkEta_,"hTrkEta","eta of all tracks",40,-4.,4.);
  book1D(hTrkPhi_,"hTrkPhi","phi of all tracks",40,-4.,4.);
  book1D(hRatioDxySigmaDxyBS_,"hRatioDxySigmaDxyBS","ratio of transverse impact parameter and its significance wrt beam spot",60,-10.,10);
  book1D(hRatioDxySigmaDxyPV_,"hRatioDxySigmaDxyPV","ratio of transverse impact parameter and its significance wrt PV",60,-10.,10);
  book1D(hRatioDzSigmaDzBS_,"hRatioDzSigmaDzBS","ratio of longitudinal impact parameter and its significance wrt beam spot",80,-20.,20);
  book1D(hRatioDzSigmaDzPV_,"hRatioDzSigmaDzPV","ratio of longitudinal impact parameter and its significance wrt PV",80,-20.,20);
  book1D(hTrkChi2_,"hTrkChi2","track chi2",30,0.,30);
  book1D(hTrkNdof_,"hTrkNdof","track NDOF",100,0,100);

  book1D(hNgoodTrk_,"hNgoodTrk","number of good tracks",50,-0.5,49.5);

  book1D(hGoodTrkPt500_,"hGoodTrkPt500","pT of all good tracks with pT > 500 MeV",50,0.,50.);
  book1D(hGoodTrkEta500_,"hGoodTrkEta500","eta of all good tracks pT > 500 MeV",40,-4.,4.);
  book1D(hGoodTrkPhi500_,"hGoodTrkPhi500","phi of all good tracks pT > 500 MeV",40,-4.,4.);

  book1D(hGoodTrkPt900_,"hGoodTrkPt900","pT of all good tracks with pT > 900 MeV",50,0.,50.);
  book1D(hGoodTrkEta900_,"hGoodTrkEta900","eta of all good tracks pT > 900 MeV",40,-4.,4.);
  book1D(hGoodTrkPhi900_,"hGoodTrkPhi900","phi of all good tracks pT > 900 MeV",40,-4.,4.);

  book1D(hNvertices_,"hNvertices","number of vertices",5,-0.5,4.5);
  book1D(hVertex_z_,"hVertex_z","z position of vertex; z[cm]",100,-20,20);
  book1D(hVertex_y_,"hVertex_y","y position of vertex; y[cm]",100,-5,5);
  book1D(hVertex_x_,"hVertex_x","x position of vertex; x[cm]",100,-5,5);

  book1D(hBeamSpot_z_,"hBeamSpot_z","z position of beamspot; z[cm]",100,-20,20);
  book1D(hBeamSpot_y_,"hBeamSpot_y","y position of beamspot; y[cm]",50,-10,10);
  book1D(hBeamSpot_x_,"hBeamSpot_x","x position of beamspot; x[cm]",50,-10,10);


  if (1) {
    const int Nx = 25;
    const double x1 = 0.0;
    const double x2 = 50.0;
    book1D(hLeadingTrack_pTSpectrum_,"hLeadingTrack_pTSpectrum","pT spectrum of leading track;pT(GeV/c)",Nx,x1,x2);
    book1D(hLeadingCaloJet_pTSpectrum_,"hLeadingCalo_pTSpectrum","pT spectrum of leading calo jet;pT(GeV/c)",Nx,x1,x2);
    book1D(hLeadingChargedJet_pTSpectrum_,"hLeadingChargedJet_pTSpectrum","pT spectrum of leading track jet;pT(GeV/c)",Nx,x1,x2);
    
  }
  
  if (1) {
    const int Nx = 24;
    const double x1 = -4.;
    const double x2 =  4.;
    book1D(hLeadingTrack_phiSpectrum_,"hLeadingTrack_phiSpectrum","#phi spectrum of leading track;#phi",Nx,x1,x2);
    book1D(hLeadingCaloJet_phiSpectrum_,"hLeadingCaloJet_phiSpectrum","#phi spectrum of leading calo jet;#phi",Nx,x1,x2);
    book1D(hLeadingChargedJet_phiSpectrum_,"hLeadingChargedJet_phiSpectrum","#phi spectrum of leading track jet;#phi",Nx,x1,x2);

  }
  
  if (1) {
    const int Nx = 24;
    const double x1 = -4.;
    const double x2 =  4.;
    book1D(hLeadingTrack_etaSpectrum_,"hLeadingTrack_etaSpectrum","#eta spectrum of leading track;#eta",Nx,x1,x2);
    book1D(hLeadingCaloJet_etaSpectrum_,"hLeadingCaloJet_etaSpectrum","#eta spectrum of leading calo jet;#eta",Nx,x1,x2);
    book1D(hLeadingChargedJet_etaSpectrum_,"hLeadingChargedJet_etaSpectrum","#eta spectrum of leading track jet;#eta",Nx,x1,x2);

  }


if (1) {
    const int Nx = 75;
    const double x1 = 0.0;
    const double x2 = 75.0;
    const double y1 = 0.;
    const double y2 = 10.;
    bookProfile(hdNdEtadPhi_pTMax_Toward500_,"hdNdEtadPhi_pTMax_Toward500", 
                 "Average number of tracks (pT > 500 MeV) in toward region vs leading track pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hdNdEtadPhi_pTMax_Transverse500_,"hdNdEtadPhi_pTMax_Transverse500", 
                 "Average number of tracks (pT > 500 MeV) in transverse region vs leading track pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hdNdEtadPhi_pTMax_Away500_,"hdNdEtadPhi_pTMax_Away500", 
                 "Average number of tracks (pT > 500 MeV) in away region vs leading track pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
 
    bookProfile(hdNdEtadPhi_caloJet_Toward500_,"hdNdEtadPhi_caloJet_Toward500", 
                 "Average number of tracks (pT > 500 MeV) in toward region vs leading calo jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hdNdEtadPhi_caloJet_Transverse500_,"hdNdEtadPhi_caloJet_Transverse500", 
                 "Average number of tracks (pT > 500 MeV) in transverse region vs leading calo jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hdNdEtadPhi_caloJet_Away500_,"hdNdEtadPhi_caloJet_Away500", 
                 "Average number of tracks (pT > 500 MeV) in away region vs leading calo jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);   
  
    bookProfile(hdNdEtadPhi_trackJet_Toward500_,"hdNdEtadPhi_trackJet_Toward500", 
                 "Average number of tracks (pT > 500 MeV) in toward region vs leading track jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2);
    bookProfile(hdNdEtadPhi_trackJet_Transverse500_,"hdNdEtadPhi_trackJet_Transverse500", 
                 "Average number of tracks (pT > 500 MeV) in transverse region vs leading track jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hdNdEtadPhi_trackJet_Away500_,"hdNdEtadPhi_trackJet_Away500", 
                 "Average number of tracks (pT > 500 MeV) in away region vs leading track jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);


 
    bookProfile(hpTSumdEtadPhi_pTMax_Toward500_,"hpTSumdEtadPhi_pTMax_Toward500", 
                 "Average number of tracks (pT > 500 MeV) in toward region vs leading track pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_pTMax_Transverse500_,"hpTSumdEtadPhi_pTMax_Transverse500", 
                 "Average number of tracks (pT > 500 MeV) in transverse region vs leading track pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_pTMax_Away500_,"hpTSumdEtadPhi_pTMax_Away500", 
                 "Average number of tracks (pT > 500 MeV) in away region vs leading track pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
 
    bookProfile(hpTSumdEtadPhi_caloJet_Toward500_,"hpTSumdEtadPhi_caloJet_Toward500", 
                 "Average number of tracks (pT > 500 MeV) in toward region vs leading calo jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_caloJet_Transverse500_,"hpTSumdEtadPhi_caloJet_Transverse500", 
                 "Average number of tracks (pT > 500 MeV) in transverse region vs leading calo jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_caloJet_Away500_,"hpTSumdEtadPhi_caloJet_Away500", 
                 "Average number of tracks (pT > 500 MeV) in away region vs leading calo jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);   
  
    bookProfile(hpTSumdEtadPhi_trackJet_Toward500_,"hpTSumdEtadPhi_trackJet_Toward500", 
                 "Average number of tracks (pT > 500 MeV) in toward region vs leading track jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_trackJet_Transverse500_,"hpTSumdEtadPhi_trackJet_Transverse500", 
                 "Average number of tracks (pT > 500 MeV) in transverse region vs leading track jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_trackJet_Away500_,"hpTSumdEtadPhi_trackJet_Away500", 
                 "Average number of tracks (pT > 500 MeV) in away region vs leading track jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);

   
   bookProfile(hdNdEtadPhi_pTMax_Toward900_,"hdNdEtadPhi_pTMax_Toward900",
                 "Average number of tracks (pT > 900 MeV) in toward region vs leading track pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hdNdEtadPhi_pTMax_Transverse900_,"hdNdEtadPhi_pTMax_Transverse900",
                 "Average number of tracks (pT > 900 MeV) in transverse region vs leading track pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hdNdEtadPhi_pTMax_Away900_,"hdNdEtadPhi_pTMax_Away900",
                 "Average number of tracks (pT > 900 MeV) in away region vs leading track pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);

    bookProfile(hdNdEtadPhi_caloJet_Toward900_,"hdNdEtadPhi_caloJet_Toward900",
                 "Average number of tracks (pT > 900 MeV) in toward region vs leading calo jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hdNdEtadPhi_caloJet_Transverse900_,"hdNdEtadPhi_caloJet_Transverse900",
                 "Average number of tracks (pT > 900 MeV) in transverse region vs leading calo jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hdNdEtadPhi_caloJet_Away900_,"hdNdEtadPhi_caloJet_Away900",
                 "Average number of tracks (pT > 900 MeV) in away region vs leading calo jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);

    bookProfile(hdNdEtadPhi_trackJet_Toward900_,"hdNdEtadPhi_trackJet_Toward900",
                 "Average number of tracks (pT > 900 MeV) in toward region vs leading track jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2);
    bookProfile(hdNdEtadPhi_trackJet_Transverse900_,"hdNdEtadPhi_trackJet_Transverse900",
                 "Average number of tracks (pT > 900 MeV) in transverse region vs leading track jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hdNdEtadPhi_trackJet_Away900_,"hdNdEtadPhi_trackJet_Away900",
                 "Average number of tracks (pT > 900 MeV) in away region vs leading track jet pT;pT(GeV/c);dN/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);



    bookProfile(hpTSumdEtadPhi_pTMax_Toward900_,"hpTSumdEtadPhi_pTMax_Toward900",
                 "Average number of tracks (pT > 900 MeV) in toward region vs leading track pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_pTMax_Transverse900_,"hpTSumdEtadPhi_pTMax_Transverse900",
                 "Average number of tracks (pT > 900 MeV) in transverse region vs leading track pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_pTMax_Away900_,"hpTSumdEtadPhi_pTMax_Away900",
                 "Average number of tracks (pT > 900 MeV) in away region vs leading track pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);

    bookProfile(hpTSumdEtadPhi_caloJet_Toward900_,"hpTSumdEtadPhi_caloJet_Toward900",
                 "Average number of tracks (pT > 900 MeV) in toward region vs leading calo jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_caloJet_Transverse900_,"hpTSumdEtadPhi_caloJet_Transverse900",
                 "Average number of tracks (pT > 900 MeV) in transverse region vs leading calo jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_caloJet_Away900_,"hpTSumdEtadPhi_caloJet_Away900",
                 "Average number of tracks (pT > 900 MeV) in away region vs leading calo jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);

    bookProfile(hpTSumdEtadPhi_trackJet_Toward900_,"hpTSumdEtadPhi_trackJet_Toward900",
                 "Average number of tracks (pT > 900 MeV) in toward region vs leading track jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_trackJet_Transverse900_,"hpTSumdEtadPhi_trackJet_Transverse900",
                 "Average number of tracks (pT > 900 MeV) in transverse region vs leading track jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
    bookProfile(hpTSumdEtadPhi_trackJet_Away900_,"hpTSumdEtadPhi_trackJet_Away900",
                 "Average number of tracks (pT > 900 MeV) in away region vs leading track jet pT;pT(GeV/c);dpTSum/d#eta d#phi",Nx,x1,x2,y1,y2,0,0);
 
   
  }

if (1) {
    const int Nx = 20;
    const double x1 = 0.0;
    const double x2 = 20.0;

        book1D(hChargedJetMulti_,"hChargedJetMulti","Charged jet multiplicity;multiplicities",Nx,x1,x2);
        book1D(hCaloJetMulti_,"hCaloJetMulti","Calo jet multiplicity;multiplicities",Nx,x1,x2);

  }


if (1) {
    const int Nx = 60;
    const double x1 = -180.0;
    const double x2 = 180.0;

        book1D(hdPhi_maxpTTrack_tracks_,"hdPhi_maxpTTrack_tracks","delta phi between leading tracks and other tracks;#Delta#phi(leading track-track)",Nx,x1,x2);
        book1D(hdPhi_caloJet_tracks_,"hdPhi_caloJet_tracks","delta phi between leading calo jet  and tracks;#Delta#phi(leading calo jet-track)",Nx,x1,x2);
        book1D(hdPhi_chargedJet_tracks_,"hdPhi_chargedJet_tracks","delta phi between leading charged jet  and tracks;#Delta#phi(leading charged jet-track)",Nx,x1,x2);

}
            

}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::endJob(void) 
{
}

//--------------------------------------------------------------------------------------------------

void QcdUeDQM::endLuminosityBlock(const LuminosityBlock &l, 
                                     const EventSetup &iSetup)
{
  if( ! isHltConfigSuccessful_ ) return;

  // Update various histograms.

  repSummary_->Fill(1.);
  repSumMap_->Fill(0.5,0.5,1.);

}

//--------------------------------------------------------------------------------------------------

void QcdUeDQM::endRun(const Run &, const EventSetup &)
{
  if( ! isHltConfigSuccessful_ ) return;

  // End run, cleanup. TODO: can this be called several times in DQM???

}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::fill1D(std::vector<TH1F*> &hs, double val, double w)
{
  // Loop over histograms and fill if trigger has fired.

  for(size_t i=0;i<hs.size();++i) {
    if (!hltTrgDeci_.at(i))
      continue;
    hs.at(i)->Fill(val,w);
  }
}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::fill1D(std::vector<MonitorElement*> &mes, double val, double w)
{
  // Loop over histograms and fill if trigger has fired.

  for(size_t i=0;i<mes.size();++i) {
    if (!hltTrgDeci_.at(i))
      continue;
    mes.at(i)->Fill(val,w);
  }
}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::setLabel1D(std::vector<MonitorElement*> &mes)
{
  // Loop over histograms and fill if trigger has fired.
  string cut[5] = {"Nevt","vtx!=bmspt","Zvtx<10cm","pT>1GeV","trackFromVtx"};
  for(size_t i=0;i<mes.size();++i) {
    if (!hltTrgDeci_.at(i))
      continue;
    for(size_t j = 1;j < 6;j++)mes.at(i)->setBinLabel(j,cut[j-1],1);
  }
}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::fill2D(std::vector<TH2F*> &hs, double valx, double valy, double w)
{
  // Loop over histograms and fill if trigger has fired.

  for(size_t i=0;i<hs.size();++i) {
    if (!hltTrgDeci_.at(i))
      continue;
    hs.at(i)->Fill(valx, valy ,w);
  }
}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::fill2D(std::vector<MonitorElement*> &mes, double valx, double valy, double w)
{
  // Loop over histograms and fill if trigger has fired.

  for(size_t i=0;i<mes.size();++i) {
    if (!hltTrgDeci_.at(i))
      continue;
    mes.at(i)->Fill(valx, valy ,w);
  }
}
//--------------------------------------------------------------------------------------------------
void QcdUeDQM::fillProfile(std::vector<TProfile*> &hs, double valx, double valy, double w)
{
  // Loop over histograms and fill if trigger has fired.

  for(size_t i=0;i<hs.size();++i) {
    if (!hltTrgDeci_.at(i))
      continue;
    hs.at(i)->Fill(valx, valy ,w);
  }
}

//--------------------------------------------------------------------------------------------------
void QcdUeDQM::fillProfile(std::vector<MonitorElement*> &mes, double valx, double valy, double w)
{
  // Loop over histograms and fill if trigger has fired.

  for(size_t i=0;i<mes.size();++i) {
    if (!hltTrgDeci_.at(i))
      continue;
   const double y = valy*w; 
    mes.at(i)->Fill(valx, y);
  }
}

//--------------------------------------------------------------------------------------------------
bool QcdUeDQM::trackSelection(const reco::Track &trk, const reco::BeamSpot* bs, const reco::Vertex vtx, int sizevtx )
{
//-------------Fill basic histograms---------


//-------------------------------------------

 bool goodTrk = false;

 if(sizevtx!=1) return 0;    //selection events with only a vertex
 if(vtx.z()-bs->z0()==0) return 0;   //selection events with good vertex
 vtxntk_=vtxntk_-1;
    if(vtx.tracksSize()<= vtxntk_) return 0;    //selection evets with vertex recostructed with at least vtxntk  tracks

 //Fill basic information of all the tracks
  fill1D(hNtrackerLayer_,trk.hitPattern().trackerLayersWithMeasurement());
  fill1D(hNtrackerPixelLayer_,trk.hitPattern().pixelLayersWithMeasurement());
  fill1D(hNtrackerStripPixelLayer_,(trk.hitPattern().pixelLayersWithMeasurement() +  trk.hitPattern().numberOfValidStripLayersWithMonoAndStereo()));
  fill1D(hRatioPtErrorPt_,(trk.ptError()/trk.pt()));
  fill1D(hTrkPt_,trk.pt());
  fill1D(hTrkEta_,trk.eta());  
  fill1D(hTrkPhi_,trk.phi());
  fill1D(hRatioDxySigmaDxyBS_,(trk.dxy(bs->position())/trk.dxyError()));
  fill1D(hRatioDxySigmaDxyPV_,(trk.dxy(vtx.position())/trk.dxyError()));
  fill1D(hRatioDzSigmaDzBS_,(trk.dz(bs->position())/trk.dzError()));
  fill1D(hRatioDzSigmaDzPV_,(trk.dz(vtx.position())/trk.dzError()));
  fill1D(hTrkChi2_,trk.normalizedChi2());
  fill1D(hTrkNdof_,trk.ndof());   

  fill1D(hBeamSpot_x_,bs->x0());
  fill1D(hBeamSpot_y_,bs->y0()); 
  fill1D(hBeamSpot_z_,bs->z0());
 
 //number of layers
 bool layerMinCutbool=false;
    if (trk.hitPattern().trackerLayersWithMeasurement() >= minHit_ ||
                (trk.hitPattern().trackerLayersWithMeasurement()==3 && trk.hitPattern().pixelLayersWithMeasurement()==3 && allowTriplets_))
      layerMinCutbool=true;


  //number of pixel layers
 bool pxlLayerMinCutbool=false;
    if (trk.hitPattern().pixelLayersWithMeasurement() >=pxlLayerMinCut_) pxlLayerMinCutbool=true;


 // cut on the hits in pixel layers
 bool hasPIX1 = false;
    if (requirePIX1_) {
      const reco::HitPattern& p = trk.hitPattern();
      for (int i=0; i<p.numberOfHits(); i++) {
        uint32_t hit = p.getHitPattern(i);
        if (p.validHitFilter(hit) && p.pixelHitFilter(hit) && p.getLayer(hit)==1) hasPIX1 = true;
      }
    }else hasPIX1 = true;
 
 // cut on the pT error
 bool ptErrorbool=false; 
    if (trk.ptError()/trk.pt() < ptErr_pt_ || 
        (trk.hitPattern().trackerLayersWithMeasurement()==3 && trk.hitPattern().pixelLayersWithMeasurement()==3 && allowTriplets_)) ptErrorbool=true; 
 // quality cut
 bool quality_ok = true;  
 if (quality_.size()!=0) {
      quality_ok = false;
      for (unsigned int i = 0; i<quality_.size();++i) {
        if (trk.quality(quality_[i])){
          quality_ok = true;
          break;
        }
      }
    }
 //-----
 bool algo_ok = true;
    if (algorithm_.size()!=0) {
      if (std::find(algorithm_.begin(),algorithm_.end(),trk.algo())==algorithm_.end()) algo_ok = false;
    }
 

 if(bsuse_==1)
      {
    if(hasPIX1 &&  pxlLayerMinCutbool && layerMinCutbool &&  (trk.hitPattern().pixelLayersWithMeasurement() +  trk.hitPattern().numberOfValidStripLayersWithMonoAndStereo()) >= min3DHit_ && ptErrorbool && fabs(trk.pt()) >= ptMin_ &&  trk.eta() >= minRapidity_ && trk.eta() <= maxRapidity_ &&  fabs(trk.dxy(bs->position())/trk.dxyError()) < tip_ &&  fabs(trk.dz(bs->position())/trk.dzError()) < lip_  &&  trk.normalizedChi2()<=maxChi2_ &&  quality_ok &&  algo_ok)goodTrk=true ;
      }

    if(bsuse_==0)
     {
       if(hasPIX1 &&  pxlLayerMinCutbool &&  layerMinCutbool &&  (trk.hitPattern().pixelLayersWithMeasurement() +  trk.hitPattern().numberOfValidStripLayersWithMonoAndStereo())  >= min3DHit_ && ptErrorbool && fabs(trk.pt()) >= ptMin_ && trk.eta() >= minRapidity_ && trk.eta() <= maxRapidity_ && fabs(trk.dxy(vtx.position())/trk.dxyError()) < tip_ && fabs(trk.dz(vtx.position())/trk.dzError()) < lip_  && trk.normalizedChi2()<=maxChi2_ && quality_ok &&  algo_ok)goodTrk=true;
     
     }

  return goodTrk;

 }
//--------------------------------------------------------------------------------------------------
bool  QcdUeDQM::fillVtxPlots( const edm::Handle< reco::VertexCollection > vtxColl)
{
  const reco::VertexCollection theVertices = *(vtxColl.product());
  bool goodVtx = false;
  fill1D(hNvertices_,theVertices.size()); 
    for (reco::VertexCollection::const_iterator vertexIt = theVertices.begin(); vertexIt != theVertices.end(); ++vertexIt) 
      {
	fill1D(hVertex_z_,vertexIt->z());
        fill1D(hVertex_y_,vertexIt->y());
        fill1D(hVertex_x_,vertexIt->x());

        if(fabs(vertexIt->z()) < diffvtxbs_)
         {
         goodVtx = true;
         vtx1=(*vertexIt);
         
         break;
         }
      } // Loop over vertcies
   return goodVtx;
}
//--------------------------------------------------------------------------------------------------
void QcdUeDQM::fillpTMaxRelated(const std::vector<const reco::Track *> &track)
 {
   fill1D(hNgoodTrk_,track.size());
   if(track.size()>0)
   {
   fill1D(hLeadingTrack_pTSpectrum_,track[0]->pt());
   fill1D(hLeadingTrack_phiSpectrum_,track[0]->phi());
   fill1D(hLeadingTrack_etaSpectrum_,track[0]->eta());
   }
     for(size_t i = 0; i < track.size(); i++)
       {
        fill1D(hGoodTrkPt500_,track[i]->pt());
        fill1D(hGoodTrkEta500_,track[i]->eta());
        fill1D(hGoodTrkPhi500_,track[i]->phi());
        if(track[i]->pt() > 0.9)
        {
        fill1D(hGoodTrkPt900_,track[i]->pt());
        fill1D(hGoodTrkEta900_,track[i]->eta());
        fill1D(hGoodTrkPhi900_,track[i]->phi());
        }
        }
		     
 }


void QcdUeDQM::fillChargedJetSpectra(const edm::Handle<reco::CandidateView> trackJets)
{
  fill1D(hChargedJetMulti_,trackJets->size());
  for( reco::CandidateView::const_iterator f  = trackJets->begin();  f != trackJets->end(); f++) 
    {
      if(f != trackJets->begin())continue;
      fill1D(hLeadingChargedJet_pTSpectrum_,f->pt());
      fill1D(hLeadingChargedJet_etaSpectrum_,f->eta());
      fill1D(hLeadingChargedJet_phiSpectrum_,f->phi());
    } 
	
}

void QcdUeDQM::fillCaloJetSpectra(const edm::Handle<reco::CaloJetCollection> caloJets)
{
  fill1D(hCaloJetMulti_,caloJets->size());
   for( reco::CaloJetCollection::const_iterator f  = caloJets->begin();  f != caloJets->end(); f++)
     {
       if(f != caloJets->begin())continue;
       fill1D(hLeadingCaloJet_pTSpectrum_,f->pt()); 
       fill1D(hLeadingCaloJet_etaSpectrum_,f->eta());
       fill1D(hLeadingCaloJet_phiSpectrum_,f->phi());
     }
   
}

/*
 weight for transverse/toward/away region = 0.12
 

*/

void QcdUeDQM::fillUE_with_MaxpTtrack(const std::vector<const reco::Track*>  &track)
{
double w = 0.119;          
//double w = 1.;
double nTrk500_TransReg = 0;
double nTrk500_AwayReg = 0;
double nTrk500_TowardReg = 0;
 
double pTSum500_TransReg = 0;
double pTSum500_AwayReg = 0;
double pTSum500_TowardReg = 0;


double nTrk900_TransReg = 0;
double nTrk900_AwayReg = 0;
double nTrk900_TowardReg = 0;

double pTSum900_TransReg = 0;
double pTSum900_AwayReg = 0;
double pTSum900_TowardReg = 0;
   if(track.size() > 0) 
    {
     if(track[0]->pt() > 1.)
         {
	    for(size_t i = 1; i < track.size();i++)
		 {
                        
                       double dphi = (180./PI)*(deltaPhi(track[0]->phi(),track[i]->phi()));
                       fill1D(hdPhi_maxpTTrack_tracks_,dphi);
		       if(fabs(dphi)>60. && fabs(dphi)<120.)
			 {
			       pTSum500_TransReg =  pTSum500_TransReg + track[i]->pt();     
			       nTrk500_TransReg++;
                               if(track[i]->pt() > 0.9)
                               {
                               pTSum900_TransReg =  pTSum900_TransReg + track[i]->pt();
                               nTrk900_TransReg++;
                               }
			 }            
			
		       if(fabs(dphi)>120. && fabs(dphi)<180.)
			 {
			       pTSum500_AwayReg =  pTSum500_AwayReg + track[i]->pt();   
			       nTrk500_AwayReg++;
                                if(track[i]->pt() > 0.9)
                                {
                                pTSum900_AwayReg =  pTSum900_AwayReg + track[i]->pt();
                                nTrk900_AwayReg++;

                                } 
			 } 
		       
		       if(fabs(dphi)<60.)
			 {
			       pTSum500_TowardReg =  pTSum500_TowardReg + track[i]->pt();
			       nTrk500_TowardReg++;
                               if(track[i]->pt() > 0.9)
                               {
                               pTSum900_TowardReg =  pTSum900_TowardReg + track[i]->pt();
                               nTrk900_TowardReg++;
                               } 
			 }           
		     } // track loop 
		 }// leading track
             // non empty collection
               fillProfile(hdNdEtadPhi_pTMax_Toward500_, track[0]->pt(),nTrk500_TowardReg,w);
               fillProfile(hdNdEtadPhi_pTMax_Transverse500_, track[0]->pt(),nTrk500_TransReg,w);
               fillProfile(hdNdEtadPhi_pTMax_Away500_, track[0]->pt(),nTrk500_AwayReg,w);

               fillProfile(hpTSumdEtadPhi_pTMax_Toward500_,track[0]->pt() ,pTSum500_TowardReg,w);
               fillProfile(hpTSumdEtadPhi_pTMax_Transverse500_,track[0]->pt(),pTSum500_TransReg,w);
               fillProfile(hpTSumdEtadPhi_pTMax_Away500_, track[0]->pt(),pTSum500_AwayReg,w);

               fillProfile(hdNdEtadPhi_pTMax_Toward900_, track[0]->pt(),nTrk900_TowardReg,w);
               fillProfile(hdNdEtadPhi_pTMax_Transverse900_, track[0]->pt(),nTrk900_TransReg,w);
               fillProfile(hdNdEtadPhi_pTMax_Away900_, track[0]->pt(),nTrk900_AwayReg,w);

               fillProfile(hpTSumdEtadPhi_pTMax_Toward900_,track[0]->pt() ,pTSum900_TowardReg,w);
               fillProfile(hpTSumdEtadPhi_pTMax_Transverse900_,track[0]->pt(),pTSum900_TransReg,w);
               fillProfile(hpTSumdEtadPhi_pTMax_Away900_, track[0]->pt(),pTSum900_AwayReg,w);
     }
}

void QcdUeDQM::fillUE_with_ChargedJets(const std::vector<const reco::Track *>  &track, const edm::Handle<reco::CandidateView> &trackJets)
{
double w = 0.119;
double nTrk500_TransReg = 0;
double nTrk500_AwayReg = 0;
double nTrk500_TowardReg = 0;
  
double pTSum500_TransReg = 0;
double pTSum500_AwayReg = 0;
double pTSum500_TowardReg = 0;


double nTrk900_TransReg = 0;
double nTrk900_AwayReg = 0;
double nTrk900_TowardReg = 0;

double pTSum900_TransReg = 0;
double pTSum900_AwayReg = 0;
double pTSum900_TowardReg = 0;

   if(!(trackJets->empty()) && (trackJets->begin())->pt() > 1.)
         {
         double jetPhi = (trackJets->begin())->phi();
           for(size_t i = 0; i < track.size();i++)
		   {
			 double dphi = (180./PI)*(deltaPhi(jetPhi,track[i]->phi()));
			 fill1D(hdPhi_chargedJet_tracks_,dphi);
			 if(fabs(dphi)>60. && fabs(dphi)<120.)
			   {
				 pTSum500_TransReg =  pTSum500_TransReg + track[i]->pt();
				 nTrk500_TransReg++;
                                 if(track[i]->pt() > 0.9)
                                 {
                                 pTSum900_TransReg =  pTSum900_TransReg + track[i]->pt();
                                 nTrk900_TransReg++;
                                 }
			   }
			
			 if(fabs(dphi)>120. && fabs(dphi)<180.)
			   {
				 pTSum500_AwayReg =  pTSum500_AwayReg + track[i]->pt();
				 nTrk500_AwayReg++;
                                  if(track[i]->pt() > 0.9)
                                 {
                                 pTSum900_AwayReg =  pTSum900_AwayReg + track[i]->pt();
                                 nTrk900_AwayReg++;
                                 }
			   }
			 if(fabs(dphi)<60.)
			   {
				 pTSum500_TowardReg =  pTSum500_TowardReg + track[i]->pt();
				 nTrk500_TowardReg++;
                                  if(track[i]->pt() > 0.9)
                                 {
                                 pTSum900_TowardReg =  pTSum900_TowardReg + track[i]->pt();
                                 nTrk900_TowardReg++;
                                 } 
			   }
		       }// tracks loop 

		   }// leading track jet
		 
		 fillProfile(hdNdEtadPhi_trackJet_Toward500_, (trackJets->begin())->pt(),nTrk500_TowardReg,w);
		 fillProfile(hdNdEtadPhi_trackJet_Transverse500_, (trackJets->begin())->pt(),nTrk500_TransReg,w);
		 fillProfile(hdNdEtadPhi_trackJet_Away500_, (trackJets->begin())->pt(),nTrk500_AwayReg,w);
		 
		 fillProfile(hpTSumdEtadPhi_trackJet_Toward500_, (trackJets->begin())->pt(),pTSum500_TowardReg,w);
		 fillProfile(hpTSumdEtadPhi_trackJet_Transverse500_, (trackJets->begin())->pt(),pTSum500_TransReg,w);
		 fillProfile(hpTSumdEtadPhi_trackJet_Away500_, (trackJets->begin())->pt(),pTSum500_AwayReg,w);

                 fillProfile(hdNdEtadPhi_trackJet_Toward900_, (trackJets->begin())->pt(),nTrk900_TowardReg,w);
                 fillProfile(hdNdEtadPhi_trackJet_Transverse900_, (trackJets->begin())->pt(),nTrk900_TransReg,w);
                 fillProfile(hdNdEtadPhi_trackJet_Away900_, (trackJets->begin())->pt(),nTrk900_AwayReg,w);

                 fillProfile(hpTSumdEtadPhi_trackJet_Toward900_, (trackJets->begin())->pt(),pTSum900_TowardReg,w);
                 fillProfile(hpTSumdEtadPhi_trackJet_Transverse900_, (trackJets->begin())->pt(),pTSum900_TransReg,w);
                 fillProfile(hpTSumdEtadPhi_trackJet_Away900_, (trackJets->begin())->pt(),pTSum900_AwayReg,w);  
}

void QcdUeDQM:: fillUE_with_CaloJets(const std::vector<const reco::Track *>  &track, const edm::Handle<reco::CaloJetCollection> &caloJets)
{
double w = 0.119;
double nTrk500_TransReg = 0;
double nTrk500_AwayReg = 0;
double nTrk500_TowardReg = 0;
                 
double pTSum500_TransReg = 0;
double pTSum500_AwayReg = 0;
double pTSum500_TowardReg = 0;

double nTrk900_TransReg = 0;
double nTrk900_AwayReg = 0;
double nTrk900_TowardReg = 0;

double pTSum900_TransReg = 0;
double pTSum900_AwayReg = 0;
double pTSum900_TowardReg = 0;
    if(!(caloJets->empty()) && (caloJets->begin())->pt() > 1.)
          {
            double jetPhi = (caloJets->begin())->phi();
	     for(size_t i = 0; i < track.size();i++)
		   {
			 double dphi = (180./PI)*(deltaPhi(jetPhi,track[i]->phi()));
			 fill1D(hdPhi_caloJet_tracks_,dphi);
			 if(fabs(dphi)>60. && fabs(dphi)<120.)
			   {
				 pTSum500_TransReg =  pTSum500_TransReg + track[i]->pt();
				 nTrk500_TransReg++;
                                 if(track[i]->pt() > 0.9)
                                 {
                                 pTSum900_TransReg =  pTSum900_TransReg + track[i]->pt();
                                 nTrk900_TransReg++;
                                 }
			   }
			 if(fabs(dphi)>120. && fabs(dphi)<180.)
			   {
				 pTSum500_AwayReg =  pTSum500_AwayReg + track[i]->pt();
				 nTrk500_AwayReg++;
                                 if(track[i]->pt() > 0.9)
                                 {
                                 pTSum900_AwayReg =  pTSum900_AwayReg + track[i]->pt();
                                 nTrk900_AwayReg++;
                                 }
			   }
			 if(fabs(dphi)<60.)
			   {
				 pTSum500_TowardReg =  pTSum500_TowardReg + track[i]->pt();
				 nTrk500_TowardReg++;
                                 if(track[i]->pt() > 0.9)
                                 {
                                 pTSum900_TowardReg =  pTSum900_TowardReg + track[i]->pt();
                                 nTrk900_TowardReg++; 
                                 }
			   }
		       }// tracks loop 

		   }// leading calo jet
		 fillProfile(hdNdEtadPhi_caloJet_Toward500_, (caloJets->begin())->pt(),nTrk500_TowardReg,w);
		 fillProfile(hdNdEtadPhi_caloJet_Transverse500_, (caloJets->begin())->pt(),nTrk500_TransReg,w);
		 fillProfile(hdNdEtadPhi_caloJet_Away500_, (caloJets->begin())->pt(),nTrk500_AwayReg,w);
		   
		 fillProfile(hpTSumdEtadPhi_caloJet_Toward500_, (caloJets->begin())->pt(),pTSum500_TowardReg,w);
		 fillProfile(hpTSumdEtadPhi_caloJet_Transverse500_, (caloJets->begin())->pt(),pTSum500_TransReg,w);
                 fillProfile(hpTSumdEtadPhi_caloJet_Away500_, (caloJets->begin())->pt(),pTSum500_AwayReg,w);

                 fillProfile(hdNdEtadPhi_caloJet_Toward900_, (caloJets->begin())->pt(),nTrk900_TowardReg,w);
                 fillProfile(hdNdEtadPhi_caloJet_Transverse900_, (caloJets->begin())->pt(),nTrk900_TransReg,w);
                 fillProfile(hdNdEtadPhi_caloJet_Away900_, (caloJets->begin())->pt(),nTrk900_AwayReg,w);

                 fillProfile(hpTSumdEtadPhi_caloJet_Toward900_, (caloJets->begin())->pt(),pTSum900_TowardReg,w);
                 fillProfile(hpTSumdEtadPhi_caloJet_Transverse900_, (caloJets->begin())->pt(),pTSum900_TransReg,w);
                 fillProfile(hpTSumdEtadPhi_caloJet_Away900_, (caloJets->begin())->pt(),pTSum900_AwayReg,w);

                 
}

void QcdUeDQM::fillHltBits(const Event &iEvent)
{
  // Fill HLT trigger bits.

  Handle<TriggerResults> triggerResultsHLT;
  getProduct(hltUsedResName_, triggerResultsHLT, iEvent);

  for(size_t i=0;i<hltTrgBits_.size();++i) {
    if (hltTrgBits_.at(i)<0) 
      continue; //ignore unknown trigger 
    size_t tbit = hltTrgBits_.at(i);
    if (tbit<triggerResultsHLT->size()) {
      hltTrgDeci_[i] = triggerResultsHLT->accept(tbit);
  }
 }
  // fill correlation histogram
  for(size_t i=0;i<hltTrgBits_.size();++i) {
    if (hltTrgDeci_.at(i))
      h2TrigCorr_->Fill(i,i);
    for(size_t j=i+1;j<hltTrgBits_.size();++j) {
      if (hltTrgDeci_.at(i) && hltTrgDeci_.at(j))
        h2TrigCorr_->Fill(i,j);
      if (hltTrgDeci_.at(i) && !hltTrgDeci_.at(j))
        h2TrigCorr_->Fill(j,i);
    }
  }
}



//--------------------------------------------------------------------------------------------------

