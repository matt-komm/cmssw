#include "Validation/RecoMuon/plugins/MuonTrackValidator.h"
#include "DQMServices/ClientConfig/interface/FitSlicesYTool.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/plugins/CosmicParametersDefinerForTPESProducer.h"

#include "TMath.h"
#include <TF1.h>

using namespace std;
using namespace edm;

void MuonTrackValidator::beginRun(Run const&, EventSetup const& setup) {

  //  dbe_->showDirStructure();

  int j=0;
  for (unsigned int ww=0;ww<associators.size();ww++){
    for (unsigned int www=0;www<label.size();www++){

      dbe_->cd();
      InputTag algo = label[www];
      string dirName=dirName_;
      if (algo.process()!="")
	dirName+=algo.process()+"_";
      if(algo.label()!="")
	dirName+=algo.label()+"_";
      if(algo.instance()!="")
	dirName+=algo.instance()+"_";      
      if (dirName.find("Tracks")<dirName.length()){
	dirName.replace(dirName.find("Tracks"),6,"");
      }
      string assoc= associators[ww];
      if (assoc.find("Track")<assoc.length()){
	assoc.replace(assoc.find("Track"),5,"");
      }
      dirName+=assoc;
      std::replace(dirName.begin(), dirName.end(), ':', '_');
      dbe_->setCurrentFolder(dirName.c_str());

      setUpVectors();

      dbe_->goUp();
      string subDirName = dirName + "/simulation";
      dbe_->setCurrentFolder(subDirName.c_str());
      h_ptSIM.push_back( dbe_->book1D("ptSIM", "generated p_{t}", 5500, 0, 110 ) );
      h_etaSIM.push_back( dbe_->book1D("etaSIM", "generated pseudorapidity", 500, -2.5, 2.5 ) );
      h_tracksSIM.push_back( dbe_->book1D("tracksSIM","number of simulated tracks",200,-0.5,99.5) );
      h_vertposSIM.push_back( dbe_->book1D("vertposSIM","Transverse position of sim vertices",100,0.,120.) );
      
      dbe_->cd();
      dbe_->setCurrentFolder(dirName.c_str());
      h_tracks.push_back( dbe_->book1D("tracks","number of reconstructed tracks",200,-0.5,19.5) );
      h_fakes.push_back( dbe_->book1D("fakes","number of fake reco tracks",20,-0.5,19.5) );
      h_charge.push_back( dbe_->book1D("charge","charge",3,-1.5,1.5) );
      h_hits.push_back( dbe_->book1D("hits", "number of hits per track", nintHit,minHit,maxHit ) );
      h_losthits.push_back( dbe_->book1D("losthits", "number of lost hits per track", nintHit,minHit,maxHit) );
      h_nchi2.push_back( dbe_->book1D("chi2", "normalized #chi^{2}", 200, 0, 20 ) );
      h_nchi2_prob.push_back( dbe_->book1D("chi2_prob", "normalized #chi^{2} probability",100,0,1));

      h_effic.push_back( dbe_->book1D("effic","efficiency vs #eta",nint,min,max) );
      h_efficPt.push_back( dbe_->book1D("efficPt","efficiency vs pT",nintpT,minpT,maxpT) );
      h_effic_vs_hit.push_back( dbe_->book1D("effic_vs_hit","effic vs hit",nintHit,minHit,maxHit) );
      h_effic_vs_phi.push_back( dbe_->book1D("effic_vs_phi","effic vs phi",nintPhi,minPhi,maxPhi) );
      h_effic_vs_dxy.push_back( dbe_->book1D("effic_vs_dxy","effic vs dxy",nintDxy,minDxy,maxDxy) );
      h_effic_vs_dz.push_back( dbe_->book1D("effic_vs_dz","effic vs dz",nintDz,minDz,maxDz) );
      h_effic_vs_vertpos.push_back( dbe_->book1D("effic_vs_vertpos","effic vs vertpos",nintVertpos,minVertpos,maxVertpos) );
      h_effic_vs_zpos.push_back( dbe_->book1D("effic_vs_zpos","effic vs zpos",nintZpos,minZpos,maxZpos) );

      h_fakerate.push_back( dbe_->book1D("fakerate","fake rate vs #eta",nint,min,max) );
      h_fakeratePt.push_back( dbe_->book1D("fakeratePt","fake rate vs pT",nintpT,minpT,maxpT) );
      h_fake_vs_hit.push_back( dbe_->book1D("fakerate_vs_hit","fake rate vs hit",nintHit,minHit,maxHit) );
      h_fake_vs_phi.push_back( dbe_->book1D("fakerate_vs_phi","fake vs phi",nintPhi,minPhi,maxPhi) );
      h_fake_vs_dxy.push_back( dbe_->book1D("fakerate_vs_dxy","fake rate vs dxy",nintDxy,minDxy,maxDxy) );
      h_fake_vs_dz.push_back( dbe_->book1D("fakerate_vs_dz","fake vs dz",nintDz,minDz,maxDz) );


      /// this are needed to calculate efficiency during tha harvesting for the automated validation
      h_recoeta.push_back( dbe_->book1D("num_reco_eta","N of reco track vs eta",nint,min,max) );
      h_assoceta.push_back( dbe_->book1D("num_assoc(simToReco)_eta","N of associated tracks (simToReco) vs eta",nint,min,max) );
      h_assoc2eta.push_back( dbe_->book1D("num_assoc(recoToSim)_eta","N of associated (recoToSim) tracks vs eta",nint,min,max) );
      h_simuleta.push_back( dbe_->book1D("num_simul_eta","N of simulated tracks vs eta",nint,min,max) );
      h_recopT.push_back( dbe_->book1D("num_reco_pT","N of reco track vs pT",nintpT,minpT,maxpT) );
      h_assocpT.push_back( dbe_->book1D("num_assoc(simToReco)_pT","N of associated tracks (simToReco) vs pT",nintpT,minpT,maxpT) );
      h_assoc2pT.push_back( dbe_->book1D("num_assoc(recoToSim)_pT","N of associated (recoToSim) tracks vs pT",nintpT,minpT,maxpT) );
      h_simulpT.push_back( dbe_->book1D("num_simul_pT","N of simulated tracks vs pT",nintpT,minpT,maxpT) );
      //
      h_recohit.push_back( dbe_->book1D("num_reco_hit","N of reco track vs hit",nintHit,minHit,maxHit) );
      h_assochit.push_back( dbe_->book1D("num_assoc(simToReco)_hit","N of associated tracks (simToReco) vs hit",nintHit,minHit,maxHit) );
      h_assoc2hit.push_back( dbe_->book1D("num_assoc(recoToSim)_hit","N of associated (recoToSim) tracks vs hit",nintHit,minHit,maxHit) );
      h_simulhit.push_back( dbe_->book1D("num_simul_hit","N of simulated tracks vs hit",nintHit,minHit,maxHit) );
      //
      h_recophi.push_back( dbe_->book1D("num_reco_phi","N of reco track vs phi",nintPhi,minPhi,maxPhi) );
      h_assocphi.push_back( dbe_->book1D("num_assoc(simToReco)_phi","N of associated tracks (simToReco) vs phi",nintPhi,minPhi,maxPhi) );
      h_assoc2phi.push_back( dbe_->book1D("num_assoc(recoToSim)_phi","N of associated (recoToSim) tracks vs phi",nintPhi,minPhi,maxPhi) );
      h_simulphi.push_back( dbe_->book1D("num_simul_phi","N of simulated tracks vs phi",nintPhi,minPhi,maxPhi) );

      h_recodxy.push_back( dbe_->book1D("num_reco_dxy","N of reco track vs dxy",nintDxy,minDxy,maxDxy) );
      h_assocdxy.push_back( dbe_->book1D("num_assoc(simToReco)_dxy","N of associated tracks (simToReco) vs dxy",nintDxy,minDxy,maxDxy) );
      h_assoc2dxy.push_back( dbe_->book1D("num_assoc(recoToSim)_dxy","N of associated (recoToSim) tracks vs dxy",nintDxy,minDxy,maxDxy) );
      h_simuldxy.push_back( dbe_->book1D("num_simul_dxy","N of simulated tracks vs dxy",nintDxy,minDxy,maxDxy) );
      
      h_recodz.push_back( dbe_->book1D("num_reco_dz","N of reco track vs dz",nintDz,minDz,maxDz) );
      h_assocdz.push_back( dbe_->book1D("num_assoc(simToReco)_dz","N of associated tracks (simToReco) vs dz",nintDz,minDz,maxDz) );
      h_assoc2dz.push_back( dbe_->book1D("num_assoc(recoToSim)_dz","N of associated (recoToSim) tracks vs dz",nintDz,minDz,maxDz) );
      h_simuldz.push_back( dbe_->book1D("num_simul_dz","N of simulated tracks vs dz",nintDz,minDz,maxDz) );

      h_assocvertpos.push_back( dbe_->book1D("num_assoc(simToReco)_vertpos","N of associated tracks (simToReco) vs transverse vert position",nintVertpos,minVertpos,maxVertpos) );
      h_simulvertpos.push_back( dbe_->book1D("num_simul_vertpos","N of simulated tracks vs transverse vert position",nintVertpos,minVertpos,maxVertpos) );

      h_assoczpos.push_back( dbe_->book1D("num_assoc(simToReco)_zpos","N of associated tracks (simToReco) vs z vert position",nintZpos,minZpos,maxZpos) );
      h_simulzpos.push_back( dbe_->book1D("num_simul_zpos","N of simulated tracks vs z vert position",nintZpos,minZpos,maxZpos) );


      /////////////////////////////////

      h_eta.push_back( dbe_->book1D("eta", "pseudorapidity residue", 1000, -0.1, 0.1 ) );
      h_pt.push_back( dbe_->book1D("pullPt", "pull of p_{t}", 100, -10, 10 ) );
      h_pullTheta.push_back( dbe_->book1D("pullTheta","pull of #theta parameter",250,-25,25) );
      h_pullPhi.push_back( dbe_->book1D("pullPhi","pull of #phi parameter",250,-25,25) );
      h_pullDxy.push_back( dbe_->book1D("pullDxy","pull of dxy parameter",250,-25,25) );
      h_pullDz.push_back( dbe_->book1D("pullDz","pull of dz parameter",250,-25,25) );
      h_pullQoverp.push_back( dbe_->book1D("pullQoverp","pull of qoverp parameter",250,-25,25) );
      
      if (associators[ww]=="TrackAssociatorByChi2"){
	h_assochi2.push_back( dbe_->book1D("assocChi2","track association #chi^{2}",1000000,0,100000) );
	h_assochi2_prob.push_back(dbe_->book1D("assocChi2_prob","probability of association #chi^{2}",100,0,1));
      } else if (associators[ww]=="TrackAssociatorByHits"){
	h_assocFraction.push_back( dbe_->book1D("assocFraction","fraction of shared hits",200,0,2) );
	h_assocSharedHit.push_back(dbe_->book1D("assocSharedHit","number of shared hits",20,0,20));
      }

      chi2_vs_nhits.push_back( dbe_->book2D("chi2_vs_nhits","#chi^{2} vs nhits",25,0,25,100,0,10) );
      h_chi2meanhitsh.push_back( dbe_->bookProfile("chi2mean_vs_nhits","mean #chi^{2} vs nhits",25,0,25,100,0,10) );

      etares_vs_eta.push_back( dbe_->book2D("etares_vs_eta","etaresidue vs eta",nint,min,max,200,-0.1,0.1) );
      nrec_vs_nsim.push_back( dbe_->book2D("nrec_vs_nsim","nrec vs nsim",20,-0.5,19.5,20,-0.5,19.5) );

      chi2_vs_eta.push_back( dbe_->book2D("chi2_vs_eta","chi2_vs_eta",nint,min,max, 200, 0, 20 ));
      h_chi2meanh.push_back( dbe_->bookProfile("chi2mean","mean #chi^{2} vs #eta",nint,min,max, 200, 0, 20) );
      chi2_vs_phi.push_back( dbe_->book2D("chi2_vs_phi","#chi^{2} vs #phi",nintPhi,minPhi,maxPhi, 200, 0, 20 ) );
      h_chi2mean_vs_phi.push_back( dbe_->bookProfile("chi2mean_vs_phi","mean of #chi^{2} vs #phi",nintPhi,minPhi,maxPhi, 200, 0, 20) );

      nhits_vs_eta.push_back( dbe_->book2D("nhits_vs_eta","nhits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      nDThits_vs_eta.push_back( dbe_->book2D("nDThits_vs_eta","# DT hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      nCSChits_vs_eta.push_back( dbe_->book2D("nCSChits_vs_eta","# CSC hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      nRPChits_vs_eta.push_back( dbe_->book2D("nRPChits_vs_eta","# RPC hits vs eta",nint,min,max,nintHit,minHit,maxHit) );

      h_DThits_eta.push_back( dbe_->bookProfile("DThits_eta","mean # DT hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      h_CSChits_eta.push_back( dbe_->bookProfile("CSChits_eta","mean # CSC hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      h_RPChits_eta.push_back( dbe_->bookProfile("RPChits_eta","mean # RPC hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      h_hits_eta.push_back( dbe_->bookProfile("hits_eta","mean #hits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      nhits_vs_phi.push_back( dbe_->book2D("nhits_vs_phi","#hits vs #phi",nintPhi,minPhi,maxPhi,nintHit,minHit,maxHit) );
      h_hits_phi.push_back( dbe_->bookProfile("hits_phi","mean #hits vs #phi",nintPhi,minPhi,maxPhi, nintHit,minHit,maxHit) );

      nlosthits_vs_eta.push_back( dbe_->book2D("nlosthits_vs_eta","nlosthits vs eta",nint,min,max,nintHit,minHit,maxHit) );
      h_losthits_eta.push_back( dbe_->bookProfile("losthits_eta","losthits_eta",nint,min,max,nintHit,minHit,maxHit) );

      ptres_vs_eta.push_back(dbe_->book2D("ptres_vs_eta","ptres_vs_eta",nint,min,max, ptRes_nbin, ptRes_rangeMin, ptRes_rangeMax));
      h_ptrmsh.push_back( dbe_->book1D("ptres_vs_eta_Sigma","#sigma(#deltap_{t}/p_{t}) vs #eta",nint,min,max) );

      ptres_vs_phi.push_back( dbe_->book2D("ptres_vs_phi","p_{t} res vs #phi",nintPhi,minPhi,maxPhi, ptRes_nbin, ptRes_rangeMin, ptRes_rangeMax));
      h_ptmeanhPhi.push_back( dbe_->book1D("ptres_vs_phi_Mean","mean of p_{t} resolution vs #phi",nintPhi,minPhi,maxPhi));
      h_ptrmshPhi.push_back( dbe_->book1D("ptres_vs_phi_Sigma","#sigma(#deltap_{t}/p_{t}) vs #phi",nintPhi,minPhi,maxPhi) );

      ptres_vs_pt.push_back(dbe_->book2D("ptres_vs_pt","ptres_vs_pt",nintpT,minpT,maxpT, ptRes_nbin, ptRes_rangeMin, ptRes_rangeMax));
      h_ptmeanhPt.push_back( dbe_->book1D("ptres_vs_pt_Mean","mean of p_{t} resolution vs p_{t}",nintpT,minpT,maxpT));
      h_ptrmshPt.push_back( dbe_->book1D("ptres_vs_pt_Sigma","#sigma(#deltap_{t}/p_{t}) vs pT",nintpT,minpT,maxpT) );

      cotThetares_vs_eta.push_back(dbe_->book2D("cotThetares_vs_eta","cotThetares_vs_eta",nint,min,max,cotThetaRes_nbin, cotThetaRes_rangeMin, cotThetaRes_rangeMax));
      h_cotThetameanh.push_back( dbe_->book1D("cotThetares_vs_eta_Mean","#sigma(cot(#theta)) vs #eta Mean",nint,min,max) );
      h_cotThetarmsh.push_back( dbe_->book1D("cotThetares_vs_eta_Sigma","#sigma(cot(#theta)) vs #eta Sigma",nint,min,max) );


      cotThetares_vs_pt.push_back(dbe_->book2D("cotThetares_vs_pt","cotThetares_vs_pt",nintpT,minpT,maxpT, cotThetaRes_nbin, cotThetaRes_rangeMin, cotThetaRes_rangeMax));      
      h_cotThetameanhPt.push_back( dbe_->book1D("cotThetares_vs_pt_Mean","#sigma(cot(#theta)) vs pT Mean",nintpT,minpT,maxpT) );
      h_cotThetarmshPt.push_back( dbe_->book1D("cotThetares_vs_pt_Sigma","#sigma(cot(#theta)) vs pT Sigma",nintpT,minpT,maxpT) );


      phires_vs_eta.push_back(dbe_->book2D("phires_vs_eta","phires_vs_eta",nint,min,max, phiRes_nbin, phiRes_rangeMin, phiRes_rangeMax));
      h_phimeanh.push_back(dbe_->book1D("phires_vs_eta_Mean","mean of #phi res vs #eta",nint,min,max));
      h_phirmsh.push_back( dbe_->book1D("phires_vs_eta_Sigma","#sigma(#delta#phi) vs #eta",nint,min,max) );

      phires_vs_pt.push_back(dbe_->book2D("phires_vs_pt","phires_vs_pt",nintpT,minpT,maxpT, phiRes_nbin, phiRes_rangeMin, phiRes_rangeMax));
      h_phimeanhPt.push_back(dbe_->book1D("phires_vs_pt_Mean","mean of #phi res vs pT",nintpT,minpT,maxpT));
      h_phirmshPt.push_back( dbe_->book1D("phires_vs_pt_Sigma","#sigma(#delta#phi) vs pT",nintpT,minpT,maxpT) );

      phires_vs_phi.push_back(dbe_->book2D("phires_vs_phi","#phi res vs #phi",nintPhi,minPhi,maxPhi,phiRes_nbin, phiRes_rangeMin, phiRes_rangeMax));
      h_phimeanhPhi.push_back(dbe_->book1D("phires_vs_phi_Mean","mean of #phi res vs #phi",nintPhi,minPhi,maxPhi));
      h_phirmshPhi.push_back( dbe_->book1D("phires_vs_phi_Sigma","#sigma(#delta#phi) vs #phi",nintPhi,minPhi,maxPhi) );

      dxyres_vs_eta.push_back(dbe_->book2D("dxyres_vs_eta","dxyres_vs_eta",nint,min,max,dxyRes_nbin, dxyRes_rangeMin, dxyRes_rangeMax));
      h_dxymeanh.push_back( dbe_->book1D("dxyres_vs_eta_Mean","mean of dxyres vs #eta",nint,min,max) );
      h_dxyrmsh.push_back( dbe_->book1D("dxyres_vs_eta_Sigma","#sigma(#deltadxy) vs #eta",nint,min,max) );

      dxyres_vs_pt.push_back( dbe_->book2D("dxyres_vs_pt","dxyres_vs_pt",nintpT,minpT,maxpT,dxyRes_nbin, dxyRes_rangeMin, dxyRes_rangeMax));
      h_dxymeanhPt.push_back( dbe_->book1D("dxyres_vs_pt_Mean","mean of dxyres vs pT",nintpT,minpT,maxpT) );
      h_dxyrmshPt.push_back( dbe_->book1D("dxyres_vs_pt_Sigma","#sigmadxy vs pT",nintpT,minpT,maxpT) );

      dzres_vs_eta.push_back(dbe_->book2D("dzres_vs_eta","dzres_vs_eta",nint,min,max,dzRes_nbin, dzRes_rangeMin, dzRes_rangeMax));
      h_dzmeanh.push_back( dbe_->book1D("dzres_vs_eta_Mean","mean of dzres vs #eta",nint,min,max) );
      h_dzrmsh.push_back( dbe_->book1D("dzres_vs_eta_Sigma","#sigma(#deltadz) vs #eta",nint,min,max) );

      dzres_vs_pt.push_back(dbe_->book2D("dzres_vs_pt","dzres_vs_pt",nintpT,minpT,maxpT,dzRes_nbin, dzRes_rangeMin, dzRes_rangeMax));
      h_dzmeanhPt.push_back( dbe_->book1D("dzres_vs_pt_Mean","mean of dzres vs pT",nintpT,minpT,maxpT) );
      h_dzrmshPt.push_back( dbe_->book1D("dzres_vs_pt_Sigma","#sigma(#deltadz vs pT",nintpT,minpT,maxpT) );

      ptmean_vs_eta_phi.push_back(dbe_->bookProfile2D("ptmean_vs_eta_phi","mean p_{t} vs #eta and #phi",nintPhi,minPhi,maxPhi,nint,min,max,1000,0,1000));
      phimean_vs_eta_phi.push_back(dbe_->bookProfile2D("phimean_vs_eta_phi","mean #phi vs #eta and #phi",nintPhi,minPhi,maxPhi,nint,min,max,nintPhi,minPhi,maxPhi));

      //pulls of track params vs eta: to be used with fitslicesytool
      dxypull_vs_eta.push_back(dbe_->book2D("dxypull_vs_eta","dxypull_vs_eta",nint,min,max,100,-10,10));
      ptpull_vs_eta.push_back(dbe_->book2D("ptpull_vs_eta","ptpull_vs_eta",nint,min,max,100,-10,10)); 
      dzpull_vs_eta.push_back(dbe_->book2D("dzpull_vs_eta","dzpull_vs_eta",nint,min,max,100,-10,10)); 
      phipull_vs_eta.push_back(dbe_->book2D("phipull_vs_eta","phipull_vs_eta",nint,min,max,100,-10,10)); 
      thetapull_vs_eta.push_back(dbe_->book2D("thetapull_vs_eta","thetapull_vs_eta",nint,min,max,100,-10,10));

      h_dxypulletamean.push_back( dbe_->book1D("h_dxypulleta_Mean","mean of dxy pull vs #eta",nint,min,max) ); 
      h_ptpulletamean.push_back( dbe_->book1D("h_ptpulleta_Mean","mean of p_{t} pull vs #eta",nint,min,max) ); 
      h_dzpulletamean.push_back( dbe_->book1D("h_dzpulleta_Mean","mean of dz pull vs #eta",nint,min,max) ); 
      h_phipulletamean.push_back( dbe_->book1D("h_phipulleta_Mean","mean of #phi pull vs #eta",nint,min,max) ); 
      h_thetapulletamean.push_back( dbe_->book1D("h_thetapulleta_Mean","mean of #theta pull vs #eta",nint,min,max) );
      //      h_ptshiftetamean.push_back( dbe_->book1D("h_ptshifteta_Mean","<#deltapT/pT>[%] vs #eta",nint,min,max) ); 

      h_dxypulleta.push_back( dbe_->book1D("h_dxypulleta_Sigma","#sigma of dxy pull vs #eta",nint,min,max) ); 
      h_ptpulleta.push_back( dbe_->book1D("h_ptpulleta_Sigma","#sigma of p_{t} pull vs #eta",nint,min,max) ); 
      h_dzpulleta.push_back( dbe_->book1D("h_dzpulleta_Sigma","#sigma of dz pull vs #eta",nint,min,max) ); 
      h_phipulleta.push_back( dbe_->book1D("h_phipulleta_Sigma","#sigma of #phi pull vs #eta",nint,min,max) ); 
      h_thetapulleta.push_back( dbe_->book1D("h_thetapulleta_Sigma","#sigma of #theta pull vs #eta",nint,min,max) );
      h_ptshifteta.push_back( dbe_->book1D("ptres_vs_eta_Mean","<#deltapT/pT>[%] vs #eta",nint,min,max) ); 

      //pulls of track params vs phi
      ptpull_vs_phi.push_back(dbe_->book2D("ptpull_vs_phi","p_{t} pull vs #phi",nintPhi,minPhi,maxPhi,100,-10,10)); 
      phipull_vs_phi.push_back(dbe_->book2D("phipull_vs_phi","#phi pull vs #phi",nintPhi,minPhi,maxPhi,100,-10,10)); 
      thetapull_vs_phi.push_back(dbe_->book2D("thetapull_vs_phi","#theta pull vs #phi",nintPhi,minPhi,maxPhi,100,-10,10));

      h_ptpullphimean.push_back( dbe_->book1D("h_ptpullphi_Mean","mean of p_{t} pull vs #phi",nintPhi,minPhi,maxPhi) ); 
      h_phipullphimean.push_back( dbe_->book1D("h_phipullphi_Mean","mean of #phi pull vs #phi",nintPhi,minPhi,maxPhi) );
      h_thetapullphimean.push_back( dbe_->book1D("h_thetapullphi_Mean","mean of #theta pull vs #phi",nintPhi,minPhi,maxPhi) );

      h_ptpullphi.push_back( dbe_->book1D("h_ptpullphi_Sigma","#sigma of p_{t} pull vs #phi",nintPhi,minPhi,maxPhi) ); 
      h_phipullphi.push_back( dbe_->book1D("h_phipullphi_Sigma","#sigma of #phi pull vs #phi",nintPhi,minPhi,maxPhi) );
      h_thetapullphi.push_back( dbe_->book1D("h_thetapullphi_Sigma","#sigma of #theta pull vs #phi",nintPhi,minPhi,maxPhi) );
  
      
      nrecHit_vs_nsimHit_sim2rec.push_back( dbe_->book2D("nrecHit_vs_nsimHit_sim2rec","nrecHit vs nsimHit (Sim2RecAssoc)",nintHit,minHit,maxHit, nintHit,minHit,maxHit ));
      nrecHit_vs_nsimHit_rec2sim.push_back( dbe_->book2D("nrecHit_vs_nsimHit_rec2sim","nrecHit vs nsimHit (Rec2simAssoc)",nintHit,minHit,maxHit, nintHit,minHit,maxHit ));
      
      if (MABH) {
	h_PurityVsQuality.push_back( dbe_->book2D("PurityVsQuality","Purity vs Quality (MABH)",20,0.01,1.01,20,0.01,1.01) );
	h_effic_Quality05.push_back( dbe_->book1D("effic_Q05","efficiency vs #eta (Quality>0.5)",nint,min,max) );
	h_effic_Quality075.push_back( dbe_->book1D("effic_Q075","efficiency vs #eta (Quality>0.75)",nint,min,max) );
	h_efficPt_Quality05.push_back( dbe_->book1D("efficPt_Q05","efficiency vs pT (Quality>0.5)",nintpT,minpT,maxpT) );
	h_efficPt_Quality075.push_back( dbe_->book1D("efficPt_Q075","efficiency vs pT (Quality>0.75)",nintpT,minpT,maxpT) );
	h_effic_vs_phi_Quality05.push_back( dbe_->book1D("effic_vs_phi_Q05","effic vs phi (Quality>0.5)",nintPhi,minPhi,maxPhi) );
	h_effic_vs_phi_Quality075.push_back( dbe_->book1D("effic_vs_phi_Q075","effic vs phi (Quality>0.75)",nintPhi,minPhi,maxPhi) );
	//
	h_assoceta_Quality05.push_back( dbe_->book1D("num_assoc(simToReco)_eta_Q05","N of associated tracks (simToReco) vs eta (Quality>0.5)",nint,min,max) );
	h_assoceta_Quality075.push_back( dbe_->book1D("num_assoc(simToReco)_eta_Q075","N of associated tracks (simToReco) vs eta (Quality>0.75)",nint,min,max) );
	h_assocpT_Quality05.push_back( dbe_->book1D("num_assoc(simToReco)_pT_Q05","N of associated tracks (simToReco) vs pT (Quality>0.5)",nintpT,minpT,maxpT) );
	h_assocpT_Quality075.push_back( dbe_->book1D("num_assoc(simToReco)_pT_Q075","N of associated tracks (simToReco) vs pT (Quality>0.75)",nintpT,minpT,maxpT) );
	h_assocphi_Quality05.push_back( dbe_->book1D("num_assoc(simToReco)_phi_Q05","N of associated tracks (simToReco) vs phi (Quality>0.5)",nintPhi,minPhi,maxPhi) );
	h_assocphi_Quality075.push_back( dbe_->book1D("num_assoc(simToReco)_phi_Q075","N of associated tracks (simToReco) vs phi (Quality>0.75)",nintPhi,minPhi,maxPhi) );
      }

      if(useLogPt){
      BinLogX(dzres_vs_pt[j]->getTH2F());
      BinLogX(h_dzmeanhPt[j]->getTH1F());
      BinLogX(h_dzrmshPt[j]->getTH1F());

      BinLogX(dxyres_vs_pt[j]->getTH2F());
      BinLogX(h_dxymeanhPt[j]->getTH1F());
      BinLogX(h_dxyrmshPt[j]->getTH1F());

      BinLogX(phires_vs_pt[j]->getTH2F());
      BinLogX(h_phimeanhPt[j]->getTH1F());
      BinLogX(h_phirmshPt[j]->getTH1F());

      BinLogX(cotThetares_vs_pt[j]->getTH2F());
      BinLogX(h_cotThetameanhPt[j]->getTH1F());
      BinLogX(h_cotThetarmshPt[j]->getTH1F());

      BinLogX(ptres_vs_pt[j]->getTH2F());
      BinLogX(h_ptmeanhPt[j]->getTH1F());
      BinLogX(h_ptrmshPt[j]->getTH1F());

      BinLogX(h_efficPt[j]->getTH1F());
      BinLogX(h_fakeratePt[j]->getTH1F());
      BinLogX(h_recopT[j]->getTH1F());
      BinLogX(h_assocpT[j]->getTH1F());
      BinLogX(h_assoc2pT[j]->getTH1F());
      BinLogX(h_simulpT[j]->getTH1F());

      if (MABH) 
	{
	  BinLogX(h_efficPt_Quality05[j]->getTH1F());
	  BinLogX(h_efficPt_Quality075[j]->getTH1F());
	  //
	  BinLogX(h_assocpT_Quality05[j]->getTH1F());
	  BinLogX(h_assocpT_Quality075[j]->getTH1F());
	}
      
      j++;
      }
    }
  }
  if (UseAssociators) {
    edm::ESHandle<TrackAssociatorBase> theAssociator;
    for (unsigned int w=0;w<associators.size();w++) {
      setup.get<TrackAssociatorRecord>().get(associators[w],theAssociator);
      associator.push_back( theAssociator.product() );
    }
  }
}

void MuonTrackValidator::analyze(const edm::Event& event, const edm::EventSetup& setup){
  using namespace reco;
  
  edm::LogInfo("MuonTrackValidator") << "\n====================================================" << "\n"
				     << "Analyzing new event" << "\n"
				     << "====================================================\n" << "\n";
  edm::ESHandle<ParametersDefinerForTP> parametersDefinerTP; 
  setup.get<TrackAssociatorRecord>().get(parametersDefiner,parametersDefinerTP);    
  
  edm::Handle<TrackingParticleCollection>  TPCollectionHeff ;
  event.getByLabel(label_tp_effic,TPCollectionHeff);
  const TrackingParticleCollection tPCeff = *(TPCollectionHeff.product());
  
  edm::Handle<TrackingParticleCollection>  TPCollectionHfake ;
  event.getByLabel(label_tp_fake,TPCollectionHfake);
  const TrackingParticleCollection tPCfake = *(TPCollectionHfake.product());
  
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  event.getByLabel(bsSrc,recoBeamSpotHandle);
  reco::BeamSpot bs = *recoBeamSpotHandle;      
  
  int w=0;
  for (unsigned int ww=0;ww<associators.size();ww++){
    for (unsigned int www=0;www<label.size();www++){
      //
      //get collections from the event
      //
      edm::Handle<View<Track> >  trackCollection;

      reco::RecoToSimCollection recSimColl;
      reco::SimToRecoCollection simRecColl;
      unsigned int trackCollectionSize = 0;

      //      if(!event.getByLabel(label[www], trackCollection)&&ignoremissingtkcollection_) continue;
      if(!event.getByLabel(label[www], trackCollection)&&ignoremissingtkcollection_) {

	recSimColl.post_insert();
	simRecColl.post_insert();

      }

      else {

	trackCollectionSize = trackCollection->size();
	//associate tracks
	if(UseAssociators){
	  edm::LogVerbatim("MuonTrackValidator") << "Analyzing " 
						 << label[www].process()<<":"
						 << label[www].label()<<":"
						 << label[www].instance()<<" with "
						 << associators[ww].c_str() <<"\n";
	
	  LogTrace("MuonTrackValidator") << "Calling associateRecoToSim method" << "\n";
	  recSimColl=associator[ww]->associateRecoToSim(trackCollection,
							TPCollectionHfake,
							&event);
	  LogTrace("MuonTrackValidator") << "Calling associateSimToReco method" << "\n";
	  simRecColl=associator[ww]->associateSimToReco(trackCollection,
							TPCollectionHeff, 
							&event);
	}
	else{
	  edm::LogVerbatim("MuonTrackValidator") << "Analyzing " 
						 << label[www].process()<<":"
						 << label[www].label()<<":"
						 << label[www].instance()<<" with "
						 << associatormap.process()<<":"
						 << associatormap.label()<<":"
						 << associatormap.instance()<<"\n";
	
	  Handle<reco::SimToRecoCollection > simtorecoCollectionH;
	  event.getByLabel(associatormap,simtorecoCollectionH);
	  simRecColl= *(simtorecoCollectionH.product()); 
	
	  Handle<reco::RecoToSimCollection > recotosimCollectionH;
	  event.getByLabel(associatormap,recotosimCollectionH);
	  recSimColl= *(recotosimCollectionH.product()); 
	}

      }

      
      //
      //fill simulation histograms
      //compute number of tracks per eta interval
      //
      edm::LogVerbatim("MuonTrackValidator") << "\n# of TrackingParticles: " << tPCeff.size() << "\n";
      int ats = 0;
      int st=0;
      for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){
	bool TP_is_matched = false;
	double quality = 0.;
	bool Quality05  = false;
	bool Quality075 = false;

	TrackingParticleRef tpr(TPCollectionHeff, i);
	TrackingParticle* tp=const_cast<TrackingParticle*>(tpr.get());
	ParticleBase::Vector momentumTP; 
	ParticleBase::Point vertexTP;
	double dxySim = 0;
	double dzSim = 0; 

	//If the TrackingParticle is collison like, get the momentum and vertex at production state
	if(parametersDefiner=="LhcParametersDefinerForTP")
	  {
	    if(! tpSelector(*tp)) continue;
	    momentumTP = tp->momentum();
	    vertexTP = tp->vertex();
	    //Calcualte the impact parameters w.r.t. PCA
	    ParticleBase::Vector momentum = parametersDefinerTP->momentum(event,setup,*tp);
	    ParticleBase::Point vertex = parametersDefinerTP->vertex(event,setup,*tp);
	    dxySim = (-vertex.x()*sin(momentum.phi())+vertex.y()*cos(momentum.phi()));
	    dzSim = vertex.z() - (vertex.x()*momentum.x()+vertex.y()*momentum.y())/sqrt(momentum.perp2()) * momentum.z()/sqrt(momentum.perp2());
	  }
	//If the TrackingParticle is comics, get the momentum and vertex at PCA
	if(parametersDefiner=="CosmicParametersDefinerForTP")
	  {
	    if(! cosmictpSelector(*tp,&bs,event,setup)) continue;	
	    momentumTP = parametersDefinerTP->momentum(event,setup,*tp);
	    vertexTP = parametersDefinerTP->vertex(event,setup,*tp);
	    dxySim = (-vertexTP.x()*sin(momentumTP.phi())+vertexTP.y()*cos(momentumTP.phi()));
	    dzSim = vertexTP.z() - (vertexTP.x()*momentumTP.x()+vertexTP.y()*momentumTP.y())/sqrt(momentumTP.perp2()) * momentumTP.z()/sqrt(momentumTP.perp2());
	  }
	edm::LogVerbatim("MuonTrackValidator") <<"--------------------Selected TrackingParticle #"<<tpr.key();
	st++;

	h_ptSIM[w]->Fill(sqrt(momentumTP.perp2()));
	h_etaSIM[w]->Fill(momentumTP.eta());
	h_vertposSIM[w]->Fill(sqrt(vertexTP.perp2()));
	
	std::vector<std::pair<RefToBase<Track>, double> > rt;
	if(simRecColl.find(tpr) != simRecColl.end()){
	  rt = (std::vector<std::pair<RefToBase<Track>, double> >) simRecColl[tpr];
	  if (rt.size()!=0) {
	    RefToBase<Track> assoc_recoTrack = rt.begin()->first;
	    edm::LogVerbatim("MuonTrackValidator")<<"-----------------------------associated Track #"<<assoc_recoTrack.key();
	    TP_is_matched = true;
	    ats++;
	    quality = rt.begin()->second;
	    edm::LogVerbatim("MuonTrackValidator") << "TrackingParticle #" <<tpr.key()  
						   << " with pt=" << sqrt(momentumTP.perp2()) 
						   << " associated with quality:" << quality <<"\n";
	    if (MABH) {
	      if (quality > 0.75) {
		Quality075 = true;
		Quality05  = true;
	      } 
	      else if (quality > 0.5) {
		Quality05  = true;
	      }
	    }	    
	  }
	}else{
	  edm::LogVerbatim("MuonTrackValidator") 
	    << "TrackingParticle #" << tpr.key()
	    << " with pt,eta,phi: " 
	    << sqrt(momentumTP.perp2()) << " , "
	    << momentumTP.eta() << " , "
	    << momentumTP.phi() << " , "
	    << " NOT associated to any reco::Track" << "\n";
	}
	
	for (unsigned int f=0; f<etaintervals[w].size()-1; f++){
	  if (getEta(momentumTP.eta())>etaintervals[w][f]&&
	      getEta(momentumTP.eta())<etaintervals[w][f+1]) {
	    totSIMeta[w][f]++;
	    if (TP_is_matched) {
	      totASSeta[w][f]++;

	      if (MABH) {
		if (Quality075) {
		  totASSeta_Quality075[w][f]++;
		  totASSeta_Quality05[w][f]++;
		}
		else if (Quality05) {
		  totASSeta_Quality05[w][f]++;
		}
	      }
	    }
	  }
	} // END for (unsigned int f=0; f<etaintervals[w].size()-1; f++){
	
	for (unsigned int f=0; f<phiintervals[w].size()-1; f++){
	  if (momentumTP.phi() > phiintervals[w][f]&&
	      momentumTP.phi() <phiintervals[w][f+1]) {
	    totSIM_phi[w][f]++;
	    if (TP_is_matched) {
	      totASS_phi[w][f]++;
	      
	      if (MABH) {
		if (Quality075) {
		  totASS_phi_Quality075[w][f]++;
		  totASS_phi_Quality05[w][f]++;
		}
		else if (Quality05) {
		  totASS_phi_Quality05[w][f]++;
		}
	      }
	    }
	  }
	} // END for (unsigned int f=0; f<phiintervals[w].size()-1; f++){
	
	
	for (unsigned int f=0; f<pTintervals[w].size()-1; f++){
          if (getPt(sqrt(momentumTP.perp2()))>pTintervals[w][f]&&
              getPt(sqrt(momentumTP.perp2()))<pTintervals[w][f+1]) {
            totSIMpT[w][f]++;
	    if (TP_is_matched) {
	      totASSpT[w][f]++;
	      
	      if (MABH) {
		if (Quality075) {
		  totASSpT_Quality075[w][f]++;
		  totASSpT_Quality05[w][f]++;
		}
		else if (Quality05) { 
		  totASSpT_Quality05[w][f]++;
		}
	      }
	    }
	  }
	} // END for (unsigned int f=0; f<pTintervals[w].size()-1; f++){
	
 	for (unsigned int f=0; f<dxyintervals[w].size()-1; f++){
	  if (dxySim>dxyintervals[w][f]&&
	      dxySim<dxyintervals[w][f+1]) {
	    totSIM_dxy[w][f]++;
	    if (TP_is_matched) {
	      totASS_dxy[w][f]++;
	    }
	  }
	} // END for (unsigned int f=0; f<dxyintervals[w].size()-1; f++){

	for (unsigned int f=0; f<dzintervals[w].size()-1; f++){
	  if (dzSim>dzintervals[w][f]&&
	      dzSim<dzintervals[w][f+1]) {
	    totSIM_dz[w][f]++;
	    if (TP_is_matched) {
 	      totASS_dz[w][f]++;
 	    }
 	  }
 	} // END for (unsigned int f=0; f<dzintervals[w].size()-1; f++){

	for (unsigned int f=0; f<vertposintervals[w].size()-1; f++){
	  if (sqrt(vertexTP.perp2())>vertposintervals[w][f]&&
	      sqrt(vertexTP.perp2())<vertposintervals[w][f+1]) {
	    totSIM_vertpos[w][f]++;
	    if (TP_is_matched) {
	      totASS_vertpos[w][f]++;
	    }
	  }
	} // END for (unsigned int f=0; f<vertposintervals[w].size()-1; f++){

	for (unsigned int f=0; f<zposintervals[w].size()-1; f++){
	  if (vertexTP.z()>zposintervals[w][f]&&
	      vertexTP.z()<zposintervals[w][f+1]) {
	    totSIM_zpos[w][f]++;
	    if (TP_is_matched) {
 	      totASS_zpos[w][f]++;
 	    }
 	  }
 	} // END for (unsigned int f=0; f<zposintervals[w].size()-1; f++){
	
	std::vector<PSimHit> simhits;
	
	if (usetracker && usemuon) {
	  simhits=tp->trackPSimHit();
	} 
	else if (!usetracker && usemuon) {
	  simhits=tp->trackPSimHit(DetId::Muon);
	}
	else if (usetracker && !usemuon) {
	  simhits=tp->trackPSimHit(DetId::Tracker);
	}
	
        int tmp = std::min((int)(simhits.end()-simhits.begin()),int(maxHit-1));
	edm::LogVerbatim("MuonTrackValidator") << "\t N simhits = "<< (int)(simhits.end()-simhits.begin())<<"\n";

	totSIM_hit[w][tmp]++;
	if (TP_is_matched) totASS_hit[w][tmp]++;

	if (TP_is_matched)	 
	  {
	    RefToBase<Track> assoctrack = rt.begin()->first; 
	    nrecHit_vs_nsimHit_sim2rec[w]->Fill( assoctrack->numberOfValidHits(),(int)(simhits.end()-simhits.begin() ));
	  }
      } // End  for (TrackingParticleCollection::size_type i=0; i<tPCeff.size(); i++){
      if (st!=0) h_tracksSIM[w]->Fill(st);
      

      //
      //fill reconstructed track histograms
      // 
      edm::LogVerbatim("MuonTrackValidator") << "\n# of reco::Tracks with "
					 << label[www].process()<<":"
					 << label[www].label()<<":"
					 << label[www].instance()
					 << ": " << trackCollectionSize << "\n";
      int at=0;
      int rT=0;
      for(View<Track>::size_type i=0; i<trackCollectionSize; ++i){
        bool Track_is_matched = false; 
	RefToBase<Track> track(trackCollection, i);
	rT++;

	std::vector<std::pair<TrackingParticleRef, double> > tp;
	TrackingParticleRef tpr;

	// new logic (bidirectional)
	if (BiDirectional_RecoToSim_association) {	  
	  edm::LogVerbatim("MuonTrackValidator")<<"----------------------------------------Track #"<< track.key();

	  if(recSimColl.find(track) != recSimColl.end()) {
	    tp = recSimColl[track];	    
	    if (tp.size() != 0) {
	      tpr = tp.begin()->first;	      
	      // RtS and StR must associate the same pair !
	      if(simRecColl.find(tpr) != simRecColl.end()) {
		std::vector<std::pair<RefToBase<Track>, double> > track_checkback  = simRecColl[tpr];
		RefToBase<Track> assoc_track_checkback;
		assoc_track_checkback = track_checkback.begin()->first;

		if ( assoc_track_checkback.key() == track.key() ) {
		  edm::LogVerbatim("MuonTrackValidator")<<"------------------associated TrackingParticle #"<<tpr.key();
		  Track_is_matched = true;
		  at++;
		  double Purity = tp.begin()->second;
		  double Quality = track_checkback.begin()->second;
		  edm::LogVerbatim("MuonTrackValidator") << "reco::Track #" << track.key() << " with pt=" << track->pt() 
							 << " associated with quality:" << Purity <<"\n";
		  if (MABH) h_PurityVsQuality[w]->Fill(Quality,Purity);
		}
	      }
	    }
	  }

	  if (!Track_is_matched) edm::LogVerbatim("MuonTrackValidator") 
	    << "reco::Track #" << track.key() << " with pt=" << track->pt() << " NOT associated to any TrackingParticle" << "\n";
	}
	// old logic (bugged)
	else {
	  if(recSimColl.find(track) != recSimColl.end()){
	    tp = recSimColl[track];
	    if (tp.size()!=0) {
	      Track_is_matched = true;
	      tpr = tp.begin()->first;
	      at++;
	      edm::LogVerbatim("MuonTrackValidator") << "reco::Track #" << track.key() << " with pt=" << track->pt() 
						     << " associated with quality:" << tp.begin()->second <<"\n";
	    }
	  } else {
	    edm::LogVerbatim("MuonTrackValidator") << "reco::Track #" << track.key() << " with pt=" << track->pt()
						   << " NOT associated to any TrackingParticle" << "\n";		  
	  }
	}
	
	//Compute fake rate vs eta
	for (unsigned int f=0; f<etaintervals[w].size()-1; f++){
	  if (getEta(track->momentum().eta())>etaintervals[w][f]&&
	      getEta(track->momentum().eta())<etaintervals[w][f+1]) {
	    totRECeta[w][f]++; 
	    if (Track_is_matched) {
	      totASS2eta[w][f]++;
	    }		
	  }
	} // End for (unsigned int f=0; f<etaintervals[w].size()-1; f++){

	for (unsigned int f=0; f<phiintervals[w].size()-1; f++){
	  if (track->momentum().phi()>phiintervals[w][f]&&
	      track->momentum().phi()<phiintervals[w][f+1]) {
	    totREC_phi[w][f]++; 
	    if (Track_is_matched) {
	      totASS2_phi[w][f]++;
	    }		
	  }
	} // End for (unsigned int f=0; f<phiintervals[w].size()-1; f++){

	
	for (unsigned int f=0; f<pTintervals[w].size()-1; f++){
	  if (getPt(sqrt(track->momentum().perp2()))>pTintervals[w][f]&&
	      getPt(sqrt(track->momentum().perp2()))<pTintervals[w][f+1]) {
	    totRECpT[w][f]++; 
	    if (Track_is_matched) {
	      totASS2pT[w][f]++;
	    }	      
	  }
	} // End for (unsigned int f=0; f<pTintervals[w].size()-1; f++){

	for (unsigned int f=0; f<dxyintervals[w].size()-1; f++){
	  if (track->dxy(bs.position())>dxyintervals[w][f]&&
	      track->dxy(bs.position())<dxyintervals[w][f+1]) {
	    totREC_dxy[w][f]++; 
	    if (Track_is_matched) {
	      totASS2_dxy[w][f]++;
	    }	      
	  }
	} // End for (unsigned int f=0; f<dxyintervals[w].size()-1; f++){

	for (unsigned int f=0; f<dzintervals[w].size()-1; f++){
	  if (track->dz(bs.position())>dzintervals[w][f]&&
	      track->dz(bs.position())<dzintervals[w][f+1]) {
	    totREC_dz[w][f]++; 
	    if (Track_is_matched) {
	      totASS2_dz[w][f]++;
	    }	      
	  }
	} // End for (unsigned int f=0; f<dzintervals[w].size()-1; f++){

	int tmp = std::min((int)track->found(),int(maxHit-1));
 	totREC_hit[w][tmp]++;
	if (Track_is_matched) totASS2_hit[w][tmp]++;

	edm::LogVerbatim("MuonTrackValidator") << "\t N valid rechits = "<< (int)track->found() <<"\n";

	//Fill other histos
 	try{
	  if (!Track_is_matched) continue;

	  if (associators[ww]=="TrackAssociatorByChi2"){
	    //association chi2
	    double assocChi2 = -tp.begin()->second;//in association map is stored -chi2
	    h_assochi2[www]->Fill(assocChi2);
	    h_assochi2_prob[www]->Fill(TMath::Prob((assocChi2)*5,5));
	  }
	  else if (associators[ww]=="TrackAssociatorByHits"){
	    double fraction = tp.begin()->second;
	    h_assocFraction[www]->Fill(fraction);
	    h_assocSharedHit[www]->Fill(fraction*track->numberOfValidHits());
	  }
    
	  //nchi2 and hits global distributions
	  h_nchi2[w]->Fill(track->normalizedChi2());
	  h_nchi2_prob[w]->Fill(TMath::Prob(track->chi2(),(int)track->ndof()));
	  h_hits[w]->Fill(track->numberOfValidHits());
	  h_losthits[w]->Fill(track->numberOfLostHits());
	  chi2_vs_nhits[w]->Fill(track->numberOfValidHits(),track->normalizedChi2());
	  h_charge[w]->Fill( track->charge() );
	  
	  //Get tracking particle parameters at point of closest approach to the beamline
	  ParticleBase::Vector momentumTP = parametersDefinerTP->momentum(event,setup,*(tpr.get()));
	  ParticleBase::Point vertexTP = parametersDefinerTP->vertex(event,setup,*(tpr.get()));
	  double ptSim = sqrt(momentumTP.perp2());
	  double qoverpSim = tpr->charge()/sqrt(momentumTP.x()*momentumTP.x()+momentumTP.y()*momentumTP.y()+momentumTP.z()*momentumTP.z());
	  double thetaSim = momentumTP.theta();
	  double lambdaSim = M_PI/2-momentumTP.theta();
	  double phiSim    = momentumTP.phi();
	  double dxySim    = (-vertexTP.x()*sin(momentumTP.phi())+vertexTP.y()*cos(momentumTP.phi()));
	  double dzSim     = vertexTP.z() - (vertexTP.x()*momentumTP.x()+vertexTP.y()*momentumTP.y())/sqrt(momentumTP.perp2()) * momentumTP.z()/sqrt(momentumTP.perp2());
	  
	  TrackBase::ParameterVector rParameters = track->parameters();

	  double qoverpRec(0);
	  double qoverpErrorRec(0); 
	  double ptRec(0);
	  double ptErrorRec(0);
	  double lambdaRec(0); 
	  double lambdaErrorRec(0);
	  double phiRec(0);
	  double phiErrorRec(0);


	  //loop to decide whether to take gsfTrack (utilisation of mode-function) or common track
	  const GsfTrack* gsfTrack(0);
	  if(useGsf){
	    gsfTrack = dynamic_cast<const GsfTrack*>(&(*track));
	    if (gsfTrack==0) edm::LogInfo("MuonTrackValidator") << "Trying to access mode for a non-GsfTrack";
	  }
	  
	  if (gsfTrack) {
	    // get values from mode
	    getRecoMomentum(*gsfTrack, ptRec, ptErrorRec, qoverpRec, qoverpErrorRec, 
			    lambdaRec,lambdaErrorRec, phiRec, phiErrorRec); 
	  }
	 
	  else {
	    // get values from track (without mode) 
	    getRecoMomentum(*track, ptRec, ptErrorRec, qoverpRec, qoverpErrorRec, 
			    lambdaRec,lambdaErrorRec, phiRec, phiErrorRec); 
	  }
	 
	  double thetaRec = track->theta();
	  double ptError = ptErrorRec;
	  double ptres = ptRec - ptSim; 
	  double etares = track->eta()-momentumTP.Eta();
	  double dxyRec    = track->dxy(bs.position());
	  double dzRec     = track->dz(bs.position());
	  // eta residue; pt, k, theta, phi, dxy, dz pulls
	  double qoverpPull=(qoverpRec-qoverpSim)/qoverpErrorRec;
	  double thetaPull=(lambdaRec-lambdaSim)/lambdaErrorRec;
	  double phiDiff = phiRec - phiSim;
	  if (abs(phiDiff) > M_PI) {
	    if (phiDiff >0.) phiDiff = phiDiff - 2.*M_PI;
            else phiDiff = phiDiff + 2.*M_PI;
	  }
	  double phiPull=phiDiff/phiErrorRec;
	  double dxyPull=(dxyRec-dxySim)/track->dxyError();
	  double dzPull=(dzRec-dzSim)/track->dzError();

	  double contrib_Qoverp = ((qoverpRec-qoverpSim)/qoverpErrorRec)*
	    ((qoverpRec-qoverpSim)/qoverpErrorRec)/5;
	  double contrib_dxy = ((dxyRec-dxySim)/track->dxyError())*((dxyRec-dxySim)/track->dxyError())/5;
	  double contrib_dz = ((dzRec-dzSim)/track->dzError())*((dzRec-dzSim)/track->dzError())/5;
	  double contrib_theta = ((lambdaRec-lambdaSim)/lambdaErrorRec)*
	    ((lambdaRec-lambdaSim)/lambdaErrorRec)/5;
	  double contrib_phi = (phiDiff/phiErrorRec)*(phiDiff/phiErrorRec)/5;
	  
	  LogTrace("MuonTrackValidator") << "assocChi2=" << tp.begin()->second << "\n"
					 << "" <<  "\n"
					 << "ptREC=" << ptRec << "\n"
					 << "etaREC=" << track->eta() << "\n"
					 << "qoverpREC=" << qoverpRec << "\n"
					 << "dxyREC=" << dxyRec << "\n"
					 << "dzREC=" << dzRec << "\n"
					 << "thetaREC=" << track->theta() << "\n"
					 << "phiREC=" << phiRec << "\n"
					 << "" <<  "\n"
					 << "qoverpError()=" << qoverpErrorRec << "\n"
					 << "dxyError()=" << track->dxyError() << "\n"
					 << "dzError()=" << track->dzError() << "\n"
					 << "thetaError()=" << lambdaErrorRec << "\n"
					 << "phiError()=" << phiErrorRec << "\n"
					 << "" <<  "\n"
					 << "ptSIM=" << ptSim << "\n"
					 << "etaSIM=" << momentumTP.Eta() << "\n"    
					 << "qoverpSIM=" << qoverpSim << "\n"
					 << "dxySIM=" << dxySim << "\n"
					 << "dzSIM=" << dzSim << "\n"
					 << "thetaSIM=" << M_PI/2-lambdaSim << "\n"
					 << "phiSIM=" << phiSim << "\n"
					 << "" << "\n"
					 << "contrib_Qoverp=" << contrib_Qoverp << "\n"
					 << "contrib_dxy=" << contrib_dxy << "\n"
					 << "contrib_dz=" << contrib_dz << "\n"
					 << "contrib_theta=" << contrib_theta << "\n"
					 << "contrib_phi=" << contrib_phi << "\n"
					 << "" << "\n"
					 <<"chi2PULL="<<contrib_Qoverp+contrib_dxy+contrib_dz+contrib_theta+contrib_phi<<"\n";
	  
	  h_pullQoverp[w]->Fill(qoverpPull);
	  h_pullTheta[w]->Fill(thetaPull);
	  h_pullPhi[w]->Fill(phiPull);
	  h_pullDxy[w]->Fill(dxyPull);
	  h_pullDz[w]->Fill(dzPull);


	  h_pt[w]->Fill(ptres/ptError);
	  h_eta[w]->Fill(etares);
	  etares_vs_eta[w]->Fill(getEta(track->eta()),etares);
 

	  //chi2 and #hit vs eta: fill 2D histos
	  chi2_vs_eta[w]->Fill(getEta(track->eta()),track->normalizedChi2());
	  nhits_vs_eta[w]->Fill(getEta(track->eta()),track->numberOfValidHits());
	  nDThits_vs_eta[w]->Fill(getEta(track->eta()),track->hitPattern().numberOfValidMuonDTHits());
	  nCSChits_vs_eta[w]->Fill(getEta(track->eta()),track->hitPattern().numberOfValidMuonCSCHits());
	  nRPChits_vs_eta[w]->Fill(getEta(track->eta()),track->hitPattern().numberOfValidMuonRPCHits());

	  nlosthits_vs_eta[w]->Fill(getEta(track->eta()),track->numberOfLostHits());

	  //resolution of track params: fill 2D histos
	  dxyres_vs_eta[w]->Fill(getEta(track->eta()),dxyRec-dxySim);
	  ptres_vs_eta[w]->Fill(getEta(track->eta()),(ptRec-ptSim)/ptRec);
	  dzres_vs_eta[w]->Fill(getEta(track->eta()),dzRec-dzSim);
	  phires_vs_eta[w]->Fill(getEta(track->eta()),phiDiff);
	  cotThetares_vs_eta[w]->Fill(getEta(track->eta()), cos(thetaRec)/sin(thetaRec) - cos(thetaSim)/sin(thetaSim));
	  
	  //same as before but vs pT
	  dxyres_vs_pt[w]->Fill(getPt(ptRec),dxyRec-dxySim);
	  ptres_vs_pt[w]->Fill(getPt(ptRec),(ptRec-ptSim)/ptRec);
	  dzres_vs_pt[w]->Fill(getPt(ptRec),dzRec-dzSim);
	  phires_vs_pt[w]->Fill(getPt(ptRec),phiDiff);
 	  cotThetares_vs_pt[w]->Fill(getPt(ptRec), cos(thetaRec)/sin(thetaRec) - cos(thetaSim)/sin(thetaSim));
 	   	 
	  //pulls of track params vs eta: fill 2D histos
	  dxypull_vs_eta[w]->Fill(getEta(track->eta()),dxyPull);
	  ptpull_vs_eta[w]->Fill(getEta(track->eta()),ptres/ptError);
	  dzpull_vs_eta[w]->Fill(getEta(track->eta()),dzPull);
	  phipull_vs_eta[w]->Fill(getEta(track->eta()),phiPull);
	  thetapull_vs_eta[w]->Fill(getEta(track->eta()),thetaPull);

	  //plots vs phi
	  nhits_vs_phi[w]->Fill(phiRec,track->numberOfValidHits());
	  chi2_vs_phi[w]->Fill(phiRec,track->normalizedChi2());
	  ptmean_vs_eta_phi[w]->Fill(phiRec,getEta(track->eta()),ptRec);
	  phimean_vs_eta_phi[w]->Fill(phiRec,getEta(track->eta()),phiRec);
	  ptres_vs_phi[w]->Fill(phiRec,(ptRec-ptSim)/ptRec);
	  phires_vs_phi[w]->Fill(phiRec,phiDiff);
	  ptpull_vs_phi[w]->Fill(phiRec,ptres/ptError);
	  phipull_vs_phi[w]->Fill(phiRec,phiPull); 
	  thetapull_vs_phi[w]->Fill(phiRec,thetaPull); 
	  
	  std::vector<PSimHit> simhits;
	  
	  if (usetracker && usemuon) {
	    simhits=tpr.get()->trackPSimHit();
	  } 
	  else if (!usetracker && usemuon) {
	    simhits=tpr.get()->trackPSimHit(DetId::Muon);
	  }
	  else if (usetracker && !usemuon) {
	    simhits=tpr.get()->trackPSimHit(DetId::Tracker);
	  }
	  
	  nrecHit_vs_nsimHit_rec2sim[w]->Fill(track->numberOfValidHits(), (int)(simhits.end()-simhits.begin() ));
	  
	} // End of try{
	catch (cms::Exception e){
	  LogTrace("MuonTrackValidator") << "exception found: " << e.what() << "\n";
	}
      } // End of for(View<Track>::size_type i=0; i<trackCollectionSize; ++i){
      if (at!=0) h_tracks[w]->Fill(at);
      h_fakes[w]->Fill(rT-at);
      edm::LogVerbatim("MuonTrackValidator") << "Total Simulated: " << st << "\n"
					     << "Total Associated (simToReco): " << ats << "\n"
					     << "Total Reconstructed: " << rT << "\n"
					     << "Total Associated (recoToSim): " << at << "\n"
					     << "Total Fakes: " << rT-at << "\n";
      nrec_vs_nsim[w]->Fill(rT,st);
      w++;
    } // End of  for (unsigned int www=0;www<label.size();www++){
  } //END of for (unsigned int ww=0;ww<associators.size();ww++){
}

void MuonTrackValidator::endRun(Run const&, EventSetup const&) {

  int w=0;
  for (unsigned int ww=0;ww<associators.size();ww++){
    for (unsigned int www=0;www<label.size();www++){

      //resolution of track params: get sigma from 2D histos
      if(!skipHistoFit){
      FitSlicesYTool fsyt_dxy(dxyres_vs_eta[w]);
      fsyt_dxy.getFittedSigmaWithError(h_dxyrmsh[w]);
      fsyt_dxy.getFittedMeanWithError(h_dxymeanh[w]);
      FitSlicesYTool fsyt_dxyPt(dxyres_vs_pt[w]);
      fsyt_dxyPt.getFittedSigmaWithError(h_dxyrmshPt[w]);
      fsyt_dxyPt.getFittedMeanWithError(h_dxymeanhPt[w]);
      FitSlicesYTool fsyt_pt(ptres_vs_eta[w]);
      fsyt_pt.getFittedSigmaWithError(h_ptrmsh[w]);
      fsyt_pt.getFittedMeanWithError(h_ptshifteta[w]);      
      FitSlicesYTool fsyt_ptPt(ptres_vs_pt[w]);
      fsyt_ptPt.getFittedSigmaWithError(h_ptrmshPt[w]);
      fsyt_ptPt.getFittedMeanWithError(h_ptmeanhPt[w]);
      FitSlicesYTool fsyt_ptPhi(ptres_vs_phi[w]); 
      fsyt_ptPhi.getFittedSigmaWithError(h_ptrmshPhi[w]);
      fsyt_ptPhi.getFittedMeanWithError(h_ptmeanhPhi[w]);
      FitSlicesYTool fsyt_dz(dzres_vs_eta[w]);
      fsyt_dz.getFittedSigmaWithError(h_dzrmsh[w]);
      fsyt_dz.getFittedMeanWithError(h_dzmeanh[w]);
      FitSlicesYTool fsyt_dzPt(dzres_vs_pt[w]);
      fsyt_dzPt.getFittedSigmaWithError(h_dzrmshPt[w]);
      fsyt_dzPt.getFittedMeanWithError(h_dzmeanhPt[w]);
      FitSlicesYTool fsyt_phi(phires_vs_eta[w]);
      fsyt_phi.getFittedSigmaWithError(h_phirmsh[w]);
      fsyt_phi.getFittedMeanWithError(h_phimeanh[w]);
      FitSlicesYTool fsyt_phiPt(phires_vs_pt[w]);
      fsyt_phiPt.getFittedSigmaWithError(h_phirmshPt[w]);
      fsyt_phiPt.getFittedMeanWithError(h_phimeanhPt[w]);
      FitSlicesYTool fsyt_phiPhi(phires_vs_phi[w]); 
      fsyt_phiPhi.getFittedSigmaWithError(h_phirmshPhi[w]); 
      fsyt_phiPhi.getFittedMeanWithError(h_phimeanhPhi[w]); 
      FitSlicesYTool fsyt_cotTheta(cotThetares_vs_eta[w]);
      fsyt_cotTheta.getFittedSigmaWithError(h_cotThetarmsh[w]);
      fsyt_cotTheta.getFittedMeanWithError(h_cotThetameanh[w]);
      FitSlicesYTool fsyt_cotThetaPt(cotThetares_vs_pt[w]);
      fsyt_cotThetaPt.getFittedSigmaWithError(h_cotThetarmshPt[w]);
      fsyt_cotThetaPt.getFittedMeanWithError(h_cotThetameanhPt[w]);
      }
      //chi2 and #hit vs eta: get mean from 2D histos
      doProfileX(chi2_vs_eta[w],h_chi2meanh[w]);
      doProfileX(nhits_vs_eta[w],h_hits_eta[w]);    
      doProfileX(nDThits_vs_eta[w],h_DThits_eta[w]);    
      doProfileX(nCSChits_vs_eta[w],h_CSChits_eta[w]);    
      doProfileX(nRPChits_vs_eta[w],h_RPChits_eta[w]);    

      doProfileX(nlosthits_vs_eta[w],h_losthits_eta[w]);    
      //vs phi
      doProfileX(chi2_vs_nhits[w],h_chi2meanhitsh[w]); 
      //      doProfileX(ptres_vs_eta[w],h_ptresmean_vs_eta[w]);
      //      doProfileX(phires_vs_eta[w],h_phiresmean_vs_eta[w]);
      doProfileX(chi2_vs_phi[w],h_chi2mean_vs_phi[w]);
      doProfileX(nhits_vs_phi[w],h_hits_phi[w]);
      //       doProfileX(ptres_vs_phi[w],h_ptresmean_vs_phi[w]);
      //       doProfileX(phires_vs_phi[w],h_phiresmean_vs_phi[w]);
      if(!skipHistoFit){      
      //pulls of track params vs eta: get sigma from 2D histos
      FitSlicesYTool fsyt_dxyp(dxypull_vs_eta[w]);
      fsyt_dxyp.getFittedSigmaWithError(h_dxypulleta[w]);
      fsyt_dxyp.getFittedMeanWithError(h_dxypulletamean[w]);
      FitSlicesYTool fsyt_ptp(ptpull_vs_eta[w]);
      fsyt_ptp.getFittedSigmaWithError(h_ptpulleta[w]);
      fsyt_ptp.getFittedMeanWithError(h_ptpulletamean[w]);
      FitSlicesYTool fsyt_dzp(dzpull_vs_eta[w]);
      fsyt_dzp.getFittedSigmaWithError(h_dzpulleta[w]);
      fsyt_dzp.getFittedMeanWithError(h_dzpulletamean[w]);
      FitSlicesYTool fsyt_phip(phipull_vs_eta[w]);
      fsyt_phip.getFittedSigmaWithError(h_phipulleta[w]);
      fsyt_phip.getFittedMeanWithError(h_phipulletamean[w]);
      FitSlicesYTool fsyt_thetap(thetapull_vs_eta[w]);
      fsyt_thetap.getFittedSigmaWithError(h_thetapulleta[w]);
      fsyt_thetap.getFittedMeanWithError(h_thetapulletamean[w]);
      //vs phi
      FitSlicesYTool fsyt_ptpPhi(ptpull_vs_phi[w]);
      fsyt_ptpPhi.getFittedSigmaWithError(h_ptpullphi[w]);
      fsyt_ptpPhi.getFittedMeanWithError(h_ptpullphimean[w]);
      FitSlicesYTool fsyt_phipPhi(phipull_vs_phi[w]);
      fsyt_phipPhi.getFittedSigmaWithError(h_phipullphi[w]);
      fsyt_phipPhi.getFittedMeanWithError(h_phipullphimean[w]);
      FitSlicesYTool fsyt_thetapPhi(thetapull_vs_phi[w]);
      fsyt_thetapPhi.getFittedSigmaWithError(h_thetapullphi[w]);
      fsyt_thetapPhi.getFittedMeanWithError(h_thetapullphimean[w]);
      //effic&fake
      fillPlotFromVectors(h_effic[w],totASSeta[w],totSIMeta[w],"effic");
      fillPlotFromVectors(h_fakerate[w],totASS2eta[w],totRECeta[w],"fakerate");
      fillPlotFromVectors(h_efficPt[w],totASSpT[w],totSIMpT[w],"effic");
      fillPlotFromVectors(h_fakeratePt[w],totASS2pT[w],totRECpT[w],"fakerate");
      fillPlotFromVectors(h_effic_vs_hit[w],totASS_hit[w],totSIM_hit[w],"effic");
      fillPlotFromVectors(h_fake_vs_hit[w],totASS2_hit[w],totREC_hit[w],"fakerate");
      fillPlotFromVectors(h_effic_vs_phi[w],totASS_phi[w],totSIM_phi[w],"effic");
      fillPlotFromVectors(h_fake_vs_phi[w],totASS2_phi[w],totREC_phi[w],"fakerate");
      fillPlotFromVectors(h_effic_vs_dxy[w],totASS_dxy[w],totSIM_dxy[w],"effic");
      fillPlotFromVectors(h_fake_vs_dxy[w],totASS2_dxy[w],totREC_dxy[w],"fakerate");
      fillPlotFromVectors(h_effic_vs_dz[w],totASS_dz[w],totSIM_dz[w],"effic");
      fillPlotFromVectors(h_fake_vs_dz[w],totASS2_dz[w],totREC_dz[w],"fakerate");
      fillPlotFromVectors(h_effic_vs_vertpos[w],totASS_vertpos[w],totSIM_vertpos[w],"effic");
      fillPlotFromVectors(h_effic_vs_zpos[w],totASS_zpos[w],totSIM_zpos[w],"effic");

      if (MABH) 
	{
	  fillPlotFromVectors(h_effic_Quality05[w] ,totASSeta_Quality05[w] ,totSIMeta[w],"effic");
	  fillPlotFromVectors(h_effic_Quality075[w],totASSeta_Quality075[w],totSIMeta[w],"effic");
	  fillPlotFromVectors(h_efficPt_Quality05[w] ,totASSpT_Quality05[w] ,totSIMpT[w],"effic");
	  fillPlotFromVectors(h_efficPt_Quality075[w],totASSpT_Quality075[w],totSIMpT[w],"effic");
	  fillPlotFromVectors(h_effic_vs_phi_Quality05[w] ,totASS_phi_Quality05[w],totSIM_phi[w],"effic");
	  fillPlotFromVectors(h_effic_vs_phi_Quality075[w],totASS_phi_Quality075[w],totSIM_phi[w],"effic");
	}
      }

      fillPlotFromVector(h_recoeta[w],totRECeta[w]);
      fillPlotFromVector(h_simuleta[w],totSIMeta[w]);
      fillPlotFromVector(h_assoceta[w],totASSeta[w]);
      fillPlotFromVector(h_assoc2eta[w],totASS2eta[w]);

      fillPlotFromVector(h_recopT[w],totRECpT[w]);
      fillPlotFromVector(h_simulpT[w],totSIMpT[w]);
      fillPlotFromVector(h_assocpT[w],totASSpT[w]);
      fillPlotFromVector(h_assoc2pT[w],totASS2pT[w]);

      fillPlotFromVector(h_recohit[w],totREC_hit[w]);
      fillPlotFromVector(h_simulhit[w],totSIM_hit[w]);
      fillPlotFromVector(h_assochit[w],totASS_hit[w]);
      fillPlotFromVector(h_assoc2hit[w],totASS2_hit[w]);

      fillPlotFromVector(h_recophi[w],totREC_phi[w]);
      fillPlotFromVector(h_simulphi[w],totSIM_phi[w]);
      fillPlotFromVector(h_assocphi[w],totASS_phi[w]);
      fillPlotFromVector(h_assoc2phi[w],totASS2_phi[w]);

      fillPlotFromVector(h_recodxy[w],totREC_dxy[w]);
      fillPlotFromVector(h_simuldxy[w],totSIM_dxy[w]);
      fillPlotFromVector(h_assocdxy[w],totASS_dxy[w]);
      fillPlotFromVector(h_assoc2dxy[w],totASS2_dxy[w]);

      fillPlotFromVector(h_recodz[w],totREC_dz[w]);
      fillPlotFromVector(h_simuldz[w],totSIM_dz[w]);
      fillPlotFromVector(h_assocdz[w],totASS_dz[w]);
      fillPlotFromVector(h_assoc2dz[w],totASS2_dz[w]);

      fillPlotFromVector(h_simulvertpos[w],totSIM_vertpos[w]);
      fillPlotFromVector(h_assocvertpos[w],totASS_vertpos[w]);

      fillPlotFromVector(h_simulzpos[w],totSIM_zpos[w]);
      fillPlotFromVector(h_assoczpos[w],totASS_zpos[w]);
      
      if (MABH) {
	fillPlotFromVector(h_assoceta_Quality05[w] ,totASSeta_Quality05[w]);
	fillPlotFromVector(h_assoceta_Quality075[w],totASSeta_Quality075[w]);
	fillPlotFromVector(h_assocpT_Quality05[w] ,totASSpT_Quality05[w]);
	fillPlotFromVector(h_assocpT_Quality075[w],totASSpT_Quality075[w]);
	fillPlotFromVector(h_assocphi_Quality05[w] ,totASS_phi_Quality05[w]);
	fillPlotFromVector(h_assocphi_Quality075[w],totASS_phi_Quality075[w]);
      }
      
      w++;
    }
  }
  
  if ( out.size() != 0 && dbe_ ) dbe_->save(out);
}


void 
MuonTrackValidator::getRecoMomentum (const reco::Track& track, double& pt, double& ptError,
				      double& qoverp, double& qoverpError, double& lambda,double& lambdaError,  double& phi, double& phiError ) const {
  pt = track.pt();
  ptError = track.ptError();
  qoverp = track.qoverp();
  qoverpError = track.qoverpError();
  lambda = track.lambda();
  lambdaError = track.lambdaError(); 
  phi = track.phi(); 
  phiError = track.phiError();

}

void 
MuonTrackValidator::getRecoMomentum (const reco::GsfTrack& gsfTrack, double& pt, double& ptError,
				      double& qoverp, double& qoverpError, double& lambda,double& lambdaError,  double& phi, double& phiError  ) const {

  pt = gsfTrack.ptMode();
  ptError = gsfTrack.ptModeError();
  qoverp = gsfTrack.qoverpMode();
  qoverpError = gsfTrack.qoverpModeError();
  lambda = gsfTrack.lambdaMode();
  lambdaError = gsfTrack.lambdaModeError(); 
  phi = gsfTrack.phiMode(); 
  phiError = gsfTrack.phiModeError();

}

