#include "RecoTracker/TkDetLayers/interface/CompositeTECPetal.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoTracker/TkDetLayers/interface/ForwardDiskSectorBuilderFromWedges.h"
#include "RecoTracker/TkDetLayers/interface/LayerCrossingSide.h"
#include "RecoTracker/TkDetLayers/interface/DetGroupMerger.h"
#include "RecoTracker/TkDetLayers/interface/CompatibleDetToGroupAdder.h"

#include "TrackingTools/DetLayers/interface/DetLayerException.h"
#include "TrackingTools/PatternTools/interface/MeasurementEstimator.h"
#include "TrackingTools/GeomPropagators/interface/HelixForwardPlaneCrossing.h"

#include "RecoTracker/TkDetLayers/interface/TkDetUtil.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include<algorithm>
#include<numeric>
#include<iterator>

using namespace std;

typedef GeometricSearchDet::DetWithState DetWithState;


namespace details {

  struct Mean {
    float operator()(const GeometricSearchDet*  a, const GeometricSearchDet* b) const {
      return 0.5*(b->position().perp()+a->position().perp());
    }
    float operator()(float  a, float b) const {
      return 0.5*(b+a);
    }
  };

  void fillBoundaries(std::vector<const GeometricSearchDet*> const & dets,
		      std::vector<float> & boundaries) {
    boundaries.resize(dets.size());
    std::transform(dets.begin(), dets.end(),  boundaries.begin(),
		   boost::bind(&GlobalPoint::perp,boost::bind(&GeometricSearchDet::position,_1))
		   );
    std::adjacent_difference(boundaries.begin(),boundaries.end(),  boundaries.begin(), Mean());
  }

  int findBin(std::vector<float> const & boundaries, float r) {
    return  
      std::lower_bound(boundaries.begin()+1,boundaries.end(),r)
      -boundaries.begin()-1;
  }

}


CompositeTECPetal::CompositeTECPetal(vector<const TECWedge*>& innerWedges,
				     vector<const TECWedge*>& outerWedges) : 
  theFrontComps(innerWedges.begin(),innerWedges.end()), 
  theBackComps(outerWedges.begin(),outerWedges.end())
{
  theComps.assign(theFrontComps.begin(),theFrontComps.end());
  theComps.insert(theComps.end(),theBackComps.begin(),theBackComps.end());

  details::fillBoundaries( theFrontComps, theFrontBoundaries);
  details::fillBoundaries( theBackComps, theBackBoundaries);


  for(vector<const GeometricSearchDet*>::const_iterator it=theComps.begin();
      it!=theComps.end();it++){  
    theBasicComps.insert(theBasicComps.end(),	
			 (**it).basicComponents().begin(),
			 (**it).basicComponents().end());
  }


  //the Wedge are already R ordered
  //sort( theWedges.begin(), theWedges.end(), DetLessR());
  //sort( theFrontWedges.begin(), theFrontWedges.end(), DetLessR() );
  //sort( theBackWedges.begin(), theBackWedges.end(), DetLessR() );
  vector<const TECWedge*> allWedges;
  allWedges.assign(innerWedges.begin(),innerWedges.end());
  allWedges.insert(allWedges.end(),outerWedges.begin(),outerWedges.end());

  theDiskSector  = ForwardDiskSectorBuilderFromWedges()( allWedges );
  theFrontSector = ForwardDiskSectorBuilderFromWedges()( innerWedges);
  theBackSector  = ForwardDiskSectorBuilderFromWedges()( outerWedges);

  //--------- DEBUG INFO --------------
  LogDebug("TkDetLayers") << "DEBUG INFO for CompositeTECPetal" ;

  for(vector<const GeometricSearchDet*>::const_iterator it=theFrontComps.begin(); 
      it!=theFrontComps.end(); it++){
    LogDebug("TkDetLayers") << "frontWedge phi,z,r: " 
			    << (*it)->surface().position().phi() << " , "
			    << (*it)->surface().position().z() <<   " , "
			    << (*it)->surface().position().perp() ;
  }

  for(vector<const GeometricSearchDet*>::const_iterator it=theBackComps.begin(); 
      it!=theBackComps.end(); it++){
    LogDebug("TkDetLayers") << "backWedge phi,z,r: " 
			    << (*it)->surface().position().phi() << " , "
			    << (*it)->surface().position().z() <<   " , "
			    << (*it)->surface().position().perp() ;
  }
  //----------------------------------- 


}


CompositeTECPetal::~CompositeTECPetal(){
  vector<const GeometricSearchDet*>::const_iterator i;
  for (i=theComps.begin(); i!=theComps.end(); i++) {
    delete *i;
  }
} 

  
pair<bool, TrajectoryStateOnSurface>
CompositeTECPetal::compatible( const TrajectoryStateOnSurface& ts, const Propagator&, 
		  const MeasurementEstimator&) const{
  edm::LogError("TkDetLayers") << "temporary dummy implementation of CompositeTECPetal::compatible()!!" ;
  return pair<bool,TrajectoryStateOnSurface>();
}


void
CompositeTECPetal::groupedCompatibleDetsV( const TrajectoryStateOnSurface& tsos,
					   const Propagator& prop,
					   const MeasurementEstimator& est,
			 std::vector<DetGroup> & result) const {

  vector<DetGroup> closestResult;
  SubLayerCrossings  crossings; 
  crossings = computeCrossings( tsos, prop.propagationDirection());
  if(! crossings.isValid()) return;

  addClosest( tsos, prop, est, crossings.closest(), closestResult); 
  LogDebug("TkDetLayers") << "in TECPetal, closestResult.size(): "<< closestResult.size();

  if (closestResult.empty()){
    vector<DetGroup> nextResult;
    addClosest( tsos, prop, est, crossings.other(), nextResult); 
    LogDebug("TkDetLayers") << "in TECPetal, nextResult.size(): "<< nextResult.size() ;
    if(nextResult.empty())    return;
    
    DetGroupElement nextGel( nextResult.front().front());  
    int crossingSide = LayerCrossingSide().endcapSide( nextGel.trajectoryState(), prop);
    DetGroupMerger::orderAndMergeTwoLevels( closestResult, nextResult, result, 
					    crossings.closestIndex(), crossingSide);
  } else {
    
  DetGroupElement closestGel( closestResult.front().front());  
  float window = computeWindowSize( closestGel.det(), closestGel.trajectoryState(), est); 

  searchNeighbors( tsos, prop, est, crossings.closest(), window,
		   closestResult, false); 

  vector<DetGroup> nextResult;
  searchNeighbors( tsos, prop, est, crossings.other(), window,
		   nextResult, true); 

  int crossingSide = LayerCrossingSide().endcapSide( closestGel.trajectoryState(), prop);
  DetGroupMerger::orderAndMergeTwoLevels( closestResult, nextResult, result, 
				          crossings.closestIndex(), crossingSide);
  }
}

SubLayerCrossings 
CompositeTECPetal::computeCrossings(const TrajectoryStateOnSurface& startingState,
				   PropagationDirection propDir) const
{
  double rho( startingState.transverseCurvature());
  
  HelixPlaneCrossing::PositionType startPos( startingState.globalPosition() );
  HelixPlaneCrossing::DirectionType startDir( startingState.globalMomentum() );
  HelixForwardPlaneCrossing crossing(startPos,startDir,rho,propDir);
  pair<bool,double> frontPath = crossing.pathLength( *theFrontSector);

  if (!frontPath.first) return SubLayerCrossings();

  GlobalPoint gFrontPoint(crossing.position(frontPath.second));
  LogDebug("TkDetLayers") 
    << "in TECPetal,front crossing : r,z,phi: (" 
    << gFrontPoint.perp() << ","
    << gFrontPoint.z() << "," 
    << gFrontPoint.phi() << ")";
  

  int frontIndex = findBin(gFrontPoint.perp(),0);
  float frontDist = fabs( findPosition(frontIndex,0).perp() - gFrontPoint.perp());
  SubLayerCrossing frontSLC( 0, frontIndex, gFrontPoint);



  pair<bool,double> backPath = crossing.pathLength( *theBackSector);

  if (!backPath.first) return SubLayerCrossings();
  

  GlobalPoint gBackPoint( crossing.position(backPath.second));
  LogDebug("TkDetLayers") 
    << "in TECPetal,back crossing r,z,phi: (" 
    << gBackPoint.perp() << ","
    << gBackPoint.z() << "," 
    << gBackPoint.phi() << ")" ;

  int backIndex = findBin(gBackPoint.perp(),1);
  float backDist = fabs( findPosition(backIndex,1).perp() - gBackPoint.perp());
  
  SubLayerCrossing backSLC( 1, backIndex, gBackPoint);
  
  
  // 0ss: frontDisk has index=0, backDisk has index=1
  if (frontDist < backDist) {
    return SubLayerCrossings( frontSLC, backSLC, 0);
  }
  else {
    return SubLayerCrossings( backSLC, frontSLC, 1);
  } 
}

bool CompositeTECPetal::addClosest( const TrajectoryStateOnSurface& tsos,
				    const Propagator& prop,
				    const MeasurementEstimator& est,
				    const SubLayerCrossing& crossing,
				    vector<DetGroup>& result) const
{
  const vector<const GeometricSearchDet*>& sub( subLayer( crossing.subLayerIndex()));
  const GeometricSearchDet* det(sub[crossing.closestDetIndex()]);

  LogDebug("TkDetLayers") 
    << "in TECPetal, adding Wedge at r,z,phi: (" 
    << det->position().perp() << "," 
    << det->position().z() << "," 
    << det->position().phi() << ")" ;
  LogDebug("TkDetLayers") 
    << "wedge comps size: " 
    << det->basicComponents().size();

  return CompatibleDetToGroupAdder::add( *det, tsos, prop, est, result);
}



void 
CompositeTECPetal::searchNeighbors( const TrajectoryStateOnSurface& tsos,
				    const Propagator& prop,
				    const MeasurementEstimator& est,
				    const SubLayerCrossing& crossing,
				    float window, 
				    vector<DetGroup>& result,
				    bool checkClosest) const
{
  GlobalPoint gCrossingPos = crossing.position();

  const vector<const GeometricSearchDet*>& sLayer( subLayer( crossing.subLayerIndex()));
 
  int closestIndex = crossing.closestDetIndex(); 
  int negStartIndex = closestIndex-1;
  int posStartIndex = closestIndex+1;

  if (checkClosest) { // must decide if the closest is on the neg or pos side
    if ( gCrossingPos.perp() < sLayer[closestIndex]->position().perp() ) {
      posStartIndex = closestIndex;
    }
    else {
      negStartIndex = closestIndex;
    }
  }


  //const BinFinderType& binFinder = (crossing.subLayerIndex()==0 ? theFrontBinFinder : theBackBinFinder);
  int theSize = crossing.subLayerIndex()==0 ? theFrontComps.size() : theBackComps.size();
  
  typedef CompatibleDetToGroupAdder Adder;
  for (int idet=negStartIndex; idet >= 0; idet--) {
    //if(idet<0 || idet>= theSize) {edm::LogInfo(TkDetLayers) << "===== error! gone out vector bounds.idet: " << idet ;exit;}
    const GeometricSearchDet & neighborWedge = *sLayer[idet];
    if (!overlap( gCrossingPos, neighborWedge, window)) break;  // --- to check
    if (!Adder::add( neighborWedge, tsos, prop, est, result)) break;
    // maybe also add shallow crossing angle test here???
  }
  for (int idet=posStartIndex; idet <theSize; idet++) {
    //if(idet<0 || idet>= theSize) {edm::LogInfo(TkDetLayers) << "===== error! gone out vector bounds.idet: " << idet ;exit;}
    const GeometricSearchDet & neighborWedge = *sLayer[idet];
    if (!overlap( gCrossingPos, neighborWedge, window)) break;  // ---- to check
    if (!Adder::add( neighborWedge, tsos, prop, est, result)) break;
    // maybe also add shallow crossing angle test here???
  }
}

bool 
CompositeTECPetal::overlap( const GlobalPoint& gpos, const GeometricSearchDet& gsdet, float ymax)
{
  // this method is just a duplication of overlapInR 
  // adapeted for groupedCompatibleDets() needs

  // assume "fixed theta window", i.e. margin in local y = r is changing linearly with z
  float tsRadius = gpos.perp();
  float thetamin = ( max((float)0.,tsRadius-ymax))/(fabs(gpos.z())+10.); // add 10 cm contingency 
  float thetamax = ( tsRadius + ymax)/(fabs(gpos.z())-10.);

  const BoundDiskSector& wedgeSector = static_cast<const BoundDiskSector&>( gsdet.surface());                   
  float wedgeMinZ = fabs( wedgeSector.position().z()) - 0.5*wedgeSector.bounds().thickness();
  float wedgeMaxZ = fabs( wedgeSector.position().z()) + 0.5*wedgeSector.bounds().thickness(); 
  float thetaWedgeMin =  wedgeSector.innerRadius()/ wedgeMaxZ;
  float thetaWedgeMax =  wedgeSector.outerRadius()/ wedgeMinZ;
  
  // do the theta regions overlap ?

  return  !( thetamin > thetaWedgeMax || thetaWedgeMin > thetamax);
  
} 





float CompositeTECPetal::computeWindowSize( const GeomDet* det, 
					    const TrajectoryStateOnSurface& tsos, 
					    const MeasurementEstimator& est)
{
  return est.maximalLocalDisplacement(tsos, det->surface()).y();
}


int CompositeTECPetal::findBin( float R, int diskSectorType) const 
{
  return details::findBin(diskSectorType==0 ? theFrontBoundaries : theBackBoundaries,R);
}




GlobalPoint CompositeTECPetal::findPosition(int index,int diskSectorType) const 
{
  vector<const GeometricSearchDet*> const & diskSector = diskSectorType == 0 ? theFrontComps : theBackComps; 
  return (diskSector[index])->position();
}

