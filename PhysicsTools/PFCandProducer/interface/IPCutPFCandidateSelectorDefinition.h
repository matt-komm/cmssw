#ifndef PhysicsTools_PFCandProducer_IPCutPFCandidateSelectorDefinition
#define PhysicsTools_PFCandProducer_IPCutPFCandidateSelectorDefinition

/**
   \class    pf2pat::IPCutPFCandidateSelectorDefinition IPCutPFCandidateSelectorDefinition.h "PhysicsTools/PFCandProducer/interface/IPCutPFCandidateSelectorDefinition.h"
   \brief    Selects PFCandidates basing on their compatibility with vertex

   \author   Giovanni Petrucciani
   \version  $Id: IPCutPFCandidateSelectorDefinition.h,v 1.0 2010/08/09 12:47:01 gpetrucc Exp $
*/

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PhysicsTools/PFCandProducer/interface/PFCandidateSelectorDefinition.h"
#include <DataFormats/VertexReco/interface/Vertex.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include <DataFormats/GsfTrackReco/interface/GsfTrack.h>

namespace pf2pat {

  struct IPCutPFCandidateSelectorDefinition : public PFCandidateSelectorDefinition {

    IPCutPFCandidateSelectorDefinition ( const edm::ParameterSet & cfg ) :
      vertices_( cfg.getParameter<edm::InputTag> ( "vertices" ) ),
      d0Cut_( cfg.getParameter<double>("d0Cut") ),
      dzCut_( cfg.getParameter<double>("dzCut") ),
      d0SigCut_( cfg.getParameter<double>("d0SigCut") ),
      dzSigCut_( cfg.getParameter<double>("dzSigCut") ) {}
    
    void select( const HandleToCollection & hc, 
		 const edm::Event & e,
		 const edm::EventSetup& s) {
      selected_.clear();
    
      edm::Handle<reco::VertexCollection> vertices;
      e.getByLabel(vertices_, vertices);
      if (vertices->empty()) return;
      const reco::Vertex &vtx = (*vertices)[0];

      unsigned key=0;
      for( collection::const_iterator pfc = hc->begin(); 
	   pfc != hc->end(); ++pfc, ++key) {
	
	bool passing = true;
	const reco::Track *tk = 0;
	if (pfc->gsfTrackRef().isNonnull())   tk = pfc->gsfTrackRef().get();
	else if (pfc->trackRef().isNonnull()) tk = pfc->trackRef().get();

	if (tk != 0) {
	  double d0 =  fabs(tk->dxy(vtx.position()));
	  double dz =  fabs(tk->dz(vtx.position()));
	  double d0e = hypot(tk->dxyError(), hypot(vtx.xError(), vtx.yError()));
	  double dze = hypot(tk->dzError(),  vtx.zError());
	  if (d0Cut_ > 0 && d0 > d0Cut_) passing = false;
	  if (dzCut_ > 0 && dz > dzCut_) passing = false;
	  if (d0SigCut_ > 0 && d0e > 0 && d0/d0e > d0SigCut_) passing = false;
	  if (dzSigCut_ > 0 && dze > 0 && dz/dze > dzSigCut_) passing = false;
	}
	
	if( passing ) {
	  selected_.push_back( reco::PFCandidate(*pfc) );
	  reco::PFCandidatePtr ptrToMother( hc, key );
	  selected_.back().setSourceCandidatePtr( ptrToMother );
	}
      }
    }
    
    private:
    edm::InputTag vertices_;
    double d0Cut_;
    double dzCut_;
    double d0SigCut_;
    double dzSigCut_;
  };
}

#endif
