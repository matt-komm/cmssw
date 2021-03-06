#ifndef Streamer_StreamedProducts_h
#define Streamer_StreamedProducts_h

/*
  Simple packaging of all the event data that is needed to be serialized
  for transfer.

  The "other stuff in the SendEvent still needs to be
  populated.

  The product is paired with its provenance, and the entire event
  is captured in the SendEvent structure.
 */

#include <vector>
#include "DataFormats/Provenance/interface/BranchDescription.h"
#include "DataFormats/Provenance/interface/BranchID.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Provenance/interface/ParameterSetID.h"
#include "DataFormats/Provenance/interface/ParameterSetBlob.h"
#include "DataFormats/Provenance/interface/EventSelectionID.h"
#include "DataFormats/Provenance/interface/BranchListIndex.h"
#include "DataFormats/Provenance/interface/ProcessHistory.h"
#include "DataFormats/Provenance/interface/ProductStatus.h"
#include "DataFormats/Provenance/interface/BranchIDList.h"

namespace edm {

  class EDProduct;
  // ------------------------------------------

  class StreamedProduct {
  public:
    StreamedProduct() : prod_(0), desc_(0), status_(productstatus::neverCreated()), parents_(0) {}
    explicit StreamedProduct(BranchDescription const& desc) :
      prod_(0), desc_(&desc), status_(productstatus::neverCreated()), parents_(0) {}

    StreamedProduct(EDProduct const* prod,
                    BranchDescription const& desc,
                    ProductStatus status,
                    std::vector<BranchID> const* parents);

    EDProduct const* prod() const {return prod_;}
    BranchDescription const* desc() const {return desc_;}
    BranchID branchID() const {return desc_->branchID();}
    ProductStatus status() const {return status_;}
    std::vector<BranchID> const* parents() const {return parents_;}

   void clear() {
     prod_= 0;
     delete desc_;
     desc_= 0;
     status_ = productstatus::neverCreated();
     delete parents_;
     parents_ = 0;
  }

  private:
    EDProduct const* prod_;
    BranchDescription const* desc_;
    ProductStatus status_;
    std::vector<BranchID> const* parents_;
  };

  // ------------------------------------------

  typedef std::vector<StreamedProduct> SendProds;

  // ------------------------------------------

  class SendEvent {
  public:
    SendEvent() { }
    SendEvent(EventAuxiliary const& aux,
              ProcessHistory const& processHistory,
              EventSelectionIDVector const& eventSelectionIDs,
              BranchListIndexes const& branchListIndexes) :
        aux_(aux),
        processHistory_(processHistory),
        eventSelectionIDs_(eventSelectionIDs),
        branchListIndexes_(branchListIndexes),
        products_() {}
    EventAuxiliary const& aux() const {return aux_;}
    SendProds const& products() const {return products_;}
    ProcessHistory const& processHistory() const {return processHistory_;}
    EventSelectionIDVector const& eventSelectionIDs() const {return eventSelectionIDs_;}
    BranchListIndexes const& branchListIndexes() const {return branchListIndexes_;}
    SendProds & products() {return products_;}
  private:
    EventAuxiliary aux_;
    ProcessHistory processHistory_;
    EventSelectionIDVector eventSelectionIDs_;
    BranchListIndexes branchListIndexes_;
    SendProds products_;

    // other tables necessary for provenance lookup
  };

  typedef std::vector<BranchDescription> SendDescs;

  class SendJobHeader {
  public:
    typedef std::map<ParameterSetID, ParameterSetBlob> ParameterSetMap;
    SendJobHeader() { }
    SendDescs const& descs() const {return descs_;}
    ParameterSetMap const& processParameterSet() const {return processParameterSet_;}
    BranchIDLists const& branchIDLists() const {return branchIDLists_;}
    std::vector<ProcessConfiguration> const& processConfigurations() const {return processConfigurations_;}
    void push_back(BranchDescription const& bd) {descs_.push_back(bd);}
    void setParameterSetMap(ParameterSetMap const& psetMap) {processParameterSet_ = psetMap;}
    void setBranchIDLists(BranchIDLists const& bidlists) {branchIDLists_ = bidlists;}
    void setProcessConfigurations(std::vector<ProcessConfiguration> const& pcs) {processConfigurations_ = pcs;}

  private:
    SendDescs descs_;
    ParameterSetMap processParameterSet_;
    BranchIDLists branchIDLists_;
    std::vector<ProcessConfiguration> processConfigurations_;
    // trigger bit descriptions will be added here and permanent
    //  provenance values
  };


}
#endif

