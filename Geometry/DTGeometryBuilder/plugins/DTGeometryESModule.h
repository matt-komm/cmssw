#ifndef DTGeometryBuilder_DTGeometryESModule_h
#define DTGeometryBuilder_DTGeometryESModule_h

/** \class DTGeometryESModule
 * 
 *  ESProducer for DTGeometry in MuonGeometryRecord
 *
 *  $Date: 2009/01/16 11:11:46 $
 *  $Revision: 1.4 $
 *  \author N. Amapane - CERN
 */

#include <FWCore/Framework/interface/ESProducer.h>
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>

#include <string>

class DTGeometryESModule : public edm::ESProducer {
public:
  /// Constructor
  DTGeometryESModule(const edm::ParameterSet & p);

  /// Destructor
  virtual ~DTGeometryESModule();

  /// Produce DTGeometry.
  boost::shared_ptr<DTGeometry> produce(const MuonGeometryRecord& record);

private:  
  void geometryCallback_( const MuonNumberingRecord& record ) ;
  boost::shared_ptr<DTGeometry> _dtGeometry;

  bool applyAlignment_; // Switch to apply alignment corrections
  const std::string alignmentsLabel_;
  const std::string myLabel_;
  bool fromDDD_;
};
#endif






