
/*----------------------------------------------------------------------

Toy EDAnalyzer for testing purposes only.

----------------------------------------------------------------------*/

#include <stdexcept>
#include <iostream>
#include <fstream>
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "CondTools/DT/test/stubs/DTCompactMapDump.h"
#include "CondTools/DT/interface/DTCompactMapWriter.h"

namespace edmtest {

  DTCompactMapDump::DTCompactMapDump(edm::ParameterSet const& p) {
// parameters to setup 
    fileName   = p.getParameter< std::string >( "fileName" );
  }

  DTCompactMapDump::DTCompactMapDump(int i) {
  }

  DTCompactMapDump::~DTCompactMapDump() {
  }

  void DTCompactMapDump::analyze( const edm::Event& e,
                                  const edm::EventSetup& context ) {
  }

  void DTCompactMapDump::endJob() {
    std::ifstream jobDesc( fileName.c_str() );
    DTCompactMapWriter::buildSteering( jobDesc );
  }

  DEFINE_FWK_MODULE(DTCompactMapDump);
}
