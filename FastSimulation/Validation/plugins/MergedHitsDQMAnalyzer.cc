#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/MonitorElement.h"



class MergedHitsDQMAnalyzer:
    public DQMEDAnalyzer
{
    public:
        explicit MergedHitsDQMAnalyzer(const edm::ParameterSet& parameterSet)
        {
        }
        
        virtual ~MergedHitsDQMAnalyzer()
        {
        }

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    private:

        virtual void dqmBeginRun (const edm::Run& run, const edm::EventSetup& eventSetup) override
        {
        }
        
        virtual void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override
        {
        }
        
        virtual void bookHistograms(DQMStore::IBooker& booker, const edm::Run& run, const edm::EventSetup& eventSetup) override
        {
        }

  
};


void
MergedHitsDQMAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
{
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}


DEFINE_FWK_MODULE(MergedHitsDQMAnalyzer);

