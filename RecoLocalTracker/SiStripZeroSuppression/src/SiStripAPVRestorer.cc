#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripAPVRestorer.h"

#include <cmath>
#include <iostream>
#include <algorithm>




SiStripAPVRestorer::SiStripAPVRestorer(const edm::ParameterSet& conf):
  quality_cache_id(-1), noise_cache_id(-1),
  ForceNoRestore_(conf.getParameter<bool>("ForceNoRestore")),
  SelfSelectRestoreAlgo_(conf.getParameter<bool>("SelfSelectRestoreAlgo")),
  InspectAlgo_(conf.getParameter<std::string>("APVInspectMode")),
  RestoreAlgo_(conf.getParameter<std::string>("APVRestoreMode")),
  useRealMeanCM_(conf.getParameter<bool>("useRealMeanCM")),
  fraction_(conf.getParameter<double>("Fraction")),
  deviation_(conf.getParameter<uint32_t>("Deviation")),
  restoreThreshold_(conf.getParameter<double>("restoreThreshold")),
  DeltaCMThreshold_(conf.getParameter<uint32_t>("DeltaCMThreshold")),
  nSigmaNoiseDerTh_(conf.getParameter<uint32_t>("nSigmaNoiseDerTh")),
  consecThreshold_(conf.getParameter<uint32_t>("consecThreshold")),
  hitStripThreshold_(conf.getParameter<uint32_t>("hitStripThreshold")),
  nSmooth_(conf.getParameter<uint32_t>("nSmooth")),
  minStripsToFit_(conf.getParameter<uint32_t>("minStripsToFit")),
  distortionThreshold_(conf.getParameter<uint32_t>("distortionThreshold")),
  CutToAvoidSignal_(conf.getParameter<double>("CutToAvoidSignal")),
  nSaturatedStrip_(conf.getParameter<uint32_t>("nSaturatedStrip"))
	
  
{
  apvFlags_.clear();
  median_.clear();
  SmoothedMaps_.clear();
  BaselineMap_.erase(BaselineMap_.begin(), BaselineMap_.end());
}


void SiStripAPVRestorer::init(const edm::EventSetup& es){
  uint32_t n_cache_id = es.get<SiStripNoisesRcd>().cacheIdentifier();
  uint32_t q_cache_id = es.get<SiStripQualityRcd>().cacheIdentifier();
  
  if(n_cache_id != noise_cache_id) {
    es.get<SiStripNoisesRcd>().get( noiseHandle );
    noise_cache_id = n_cache_id;
  } else {
    noise_cache_id = n_cache_id;
  }
  if(q_cache_id != quality_cache_id) {
    es.get<SiStripQualityRcd>().get( qualityHandle );
    quality_cache_id = q_cache_id;
  }else {
    quality_cache_id = q_cache_id;
  }
  
  
}


int16_t SiStripAPVRestorer::inspect( const uint32_t& detId,std::vector<int16_t>& digis, const std::vector< std::pair<short,float> >& vmedians) {
  
  detId_ = detId;
  
  apvFlags_.clear();
  median_.clear();
  SmoothedMaps_.clear();
  BaselineMap_.erase(BaselineMap_.begin(), BaselineMap_.end()); 
    
  for(size_t i=0; i< vmedians.size(); ++i) median_.push_back(vmedians[i].second);
	
  if(InspectAlgo_=="BaselineFollower") return this->BaselineFollowerInspect(digis);
  if(InspectAlgo_=="AbnormalBaseline") return this->AbnormalBaselineInspect(digis);
  if(InspectAlgo_=="Null") return this->NullInspect(digis);
  if(InspectAlgo_=="BaselineAndSaturation") return this->BaselineAndSaturationInspect(digis);
  throw cms::Exception("Unregistered Inspect Algorithm") << "SiStripAPVRestorer possibilities: (Null), (AbnormalBaseline),(BaselineFollower)";
  
}


void SiStripAPVRestorer::restore( std::vector<int16_t>& digis ) {
	
  if(ForceNoRestore_) return;
  
  for( uint16_t APV=0; APV< digis.size()/128; ++APV){
    std::string	algoToUse = *( apvFlags_.begin() + APV );
    
    if ( algoToUse != ""){
      if(!SelfSelectRestoreAlgo_) algoToUse = RestoreAlgo_;
      
      if(algoToUse=="Flat"){
	this->FlatRestore(digis, APV);
      }else if(algoToUse=="BaselineFollower"){
	
	this->BaselineFollowerRestore(digis, APV, median_[APV]);
   //   }else if(algoToUse=="IterativeMedian"){
	//this->IterativeMedian(digis, APV);
      }else{
	throw cms::Exception("Unregistered Restore Algorithm") << "SiStripAPVRestorer possibilities: (Flat), (BaselineFollower)";
      }
    }
  }
  
}


//Inspect method implementation ==========================================
//========================================================================
template<typename T>
inline
int16_t SiStripAPVRestorer::BaselineFollowerInspect(std::vector<T>& digis){
  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId_);
  
  std::vector<T> singleAPVdigi;
  singleAPVdigi.clear();
  
  
  int16_t nAPVflagged = 0;
  
  CMMap::iterator itCMMap;
  if(useRealMeanCM_) itCMMap = MeanCMmap_.find(detId_);
  
  for( uint16_t APV=0; APV< digis.size()/128; ++APV){
    
    int MeanAPVCM = 128;
    if(useRealMeanCM_&&itCMMap!= MeanCMmap_.end()) MeanAPVCM =(itCMMap->second)[APV];
    
    singleAPVdigi.clear(); 
    for(int16_t strip = APV*128; strip < (APV+1)*128; ++strip){
      singleAPVdigi.push_back(digis[strip]); 
      
    }
    
    
    float DeltaCM = median_[APV] -MeanAPVCM; 
    //if(DeltaCM > DeltaCMThreshold_){		      //to be modified when code is extended
    //	apvFlags_.push_back( RestoreAlgo_ );              
    //	nAPVflagged++;
    //}else 
    
    DigiMap smoothedmap;	
    if(DeltaCM < 0 && std::abs(DeltaCM) > DeltaCMThreshold_){
      
      bool isFlat= FlatRegionsFinder(singleAPVdigi,smoothedmap, median_[APV], APV);
      if(!isFlat){
	apvFlags_.push_back( "BaselineFollower" );    //specify any algo to make the restore
	nAPVflagged++;
      }else{
	apvFlags_.push_back( "" );
      }
    } else{
      apvFlags_.push_back( "" );
    }	
    SmoothedMaps_.push_back(smoothedmap);
    
  }
  
  return nAPVflagged;
}


template<typename T>
inline
int16_t SiStripAPVRestorer::BaselineAndSaturationInspect(std::vector<T>& digis){
  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId_);
  
   
  std::vector<T> singleAPVdigi;
  singleAPVdigi.clear();
  
  
  int16_t nAPVflagged = 0;
  
  CMMap::iterator itCMMap;
  if(useRealMeanCM_) itCMMap = MeanCMmap_.find(detId_);
  
  for( uint16_t APV=0; APV< digis.size()/128; ++APV){
    
    int MeanAPVCM = 128;
    if(useRealMeanCM_&&itCMMap!= MeanCMmap_.end()) MeanAPVCM =(itCMMap->second)[APV];
    
    singleAPVdigi.clear();
   
    uint16_t nSatStrip =0;
    for(int16_t strip = APV*128; strip < (APV+1)*128; ++strip){
      singleAPVdigi.push_back(digis[strip]);
      if(digis[strip] >=1023) ++nSatStrip;
    }
    
    float DeltaCM = median_[APV] -MeanAPVCM; 
    
    
    if(DeltaCM < 0 && std::abs(DeltaCM) > DeltaCMThreshold_&&nSatStrip>= nSaturatedStrip_){
      
      apvFlags_.push_back( "RestoreAlgo_" );    //specify any algo to make the restore
      nAPVflagged++;
    } else{
      apvFlags_.push_back( "" );
    }	
    
  }
  
  return nAPVflagged;
}


template<typename T>
inline
int16_t SiStripAPVRestorer::AbnormalBaselineInspect(std::vector<T>& digis){

  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId_);
  
  typename std::vector<T>::iterator fs;
  
  int16_t nAPVflagged=0;
  
  CMMap::iterator itCMMap;
  if(useRealMeanCM_) itCMMap = MeanCMmap_.find(detId_);
  
  
  int devCount = 0, qualityCount = 0, minstrip = 0; 
  for( uint16_t APV=0; APV< digis.size()/128; ++APV){
    int MeanAPVCM = 128;
    if(useRealMeanCM_&&itCMMap!= MeanCMmap_.end()) MeanAPVCM =(itCMMap->second)[APV];
    for (uint16_t istrip=APV*128; istrip<(APV+1)*128; ++istrip){
      fs = digis.begin() + istrip;
      if ( !qualityHandle->IsStripBad(detQualityRange,istrip) ){
	qualityCount++; 
	if ( std::abs((int) *fs - MeanAPVCM) > (int)deviation_ ) devCount++;
	minstrip = std::min((int) *fs, minstrip);
      }
    }
    
    if( devCount > fraction_ * qualityCount ) {
      apvFlags_.push_back( RestoreAlgo_ );      //specify any algo to make the restore
      nAPVflagged++;
    } else {
      apvFlags_.push_back( "" );
    }
    
  }
  
  return nAPVflagged;
  
}



template<typename T>
inline
int16_t SiStripAPVRestorer::NullInspect(std::vector<T>& digis){

  SiStripQuality::Range detQualityRange = qualityHandle->getRange(detId_);

  typename std::vector<T>::iterator fs;

  int16_t nAPVflagged = 0;

  for( uint16_t APV=0; APV< digis.size()/128; ++APV){
    int zeroCount = 0, qualityCount = 0; 
    for (uint16_t istrip=APV*128; istrip<(APV+1)*128; ++istrip){
      fs = digis.begin() + istrip;
      if ( !qualityHandle->IsStripBad(detQualityRange,istrip) ){
        qualityCount++; 
        if ( (int) *fs < 1 ) zeroCount++;
      }
    }
    
    if( zeroCount > restoreThreshold_ * qualityCount ) {
      apvFlags_.push_back( RestoreAlgo_ );     //specify any algo to make the restore
      nAPVflagged++;
    } else {
      apvFlags_.push_back( "" );
    }
    
  }
 
  return nAPVflagged;

}





//Restore method implementation =====================================
//===================================================================
template<typename T>
inline
void SiStripAPVRestorer::BaselineFollowerRestore( std::vector<T>& digis, uint16_t APVn , float median){
  typename std::vector<T>::iterator firstStrip(digis.begin() + APVn*128), lastStrip(firstStrip + 128), actualStrip;
  
  std::vector<int16_t> baseline;
  baseline.clear();
  baseline.insert(baseline.begin(),128, 0);
  
  std::vector<int16_t> adcs;
  adcs.clear();
	
	 
  //============================= Copying only ADCs of one APV =============================    
  for(actualStrip= firstStrip; actualStrip < lastStrip; ++actualStrip ) adcs.push_back(*actualStrip); 	
  
  //============================= Calculate Median =========================================
  //this code was here in order to use the median calculated during the inspect 
  //but actually now we use the one calculated by the cmnsubtractor
  //float median;
  //if(median_.size()) median = median_[APVn];
  //else median = this->median(adcs);
  
  //============================= Find Flat Regions & Interpolating the baseline & subtracting the baseline  =================	
  
  if(SmoothedMaps_.size()){
    this->BaselineFollower(SmoothedMaps_[APVn], baseline, median);	
    
  } else {
    median=0;
    DigiMap  smoothedpoints;
    this->FlatRegionsFinder(adcs,smoothedpoints, median, APVn );
    this->BaselineFollower(smoothedpoints, baseline, median);		
    
  }	
  
  //============================= subtracting the baseline =============================================
  
  for(int16_t itStrip= 0 ; itStrip< 128; ++itStrip){
    //int tempDigi = digis[itStrip+APVn*128] - baseline[itStrip] + median;
    //std::cout << "BaselineFollowerRestore - detId: " << detId_ << " APV: " << APVn << " strip: " << itStrip << " digis: " << digis[itStrip+APVn*128] << " baseline: " << baseline[itStrip] << " median: " << median << " digis after baseline subtraction: " << tempDigi << std::endl;
    digis[itStrip+APVn*128] -= baseline[itStrip] - median;
  }
  
		
  //============================= storing baseline to the map =============================================	
  BaselineMap_.insert(BaselineMap_.end(),  std::pair< uint16_t, std::vector < int16_t> >(APVn, baseline));
  
  
}


template<typename T>
inline
void SiStripAPVRestorer::FlatRestore( std::vector<T>& digis, uint16_t APVn ){
 
  std::vector<int16_t> baseline;
  baseline.clear();
  baseline.insert(baseline.begin(),128, 150);
  baseline[0]=0; baseline[127]=0;
  BaselineMap_.insert(BaselineMap_.end(),  std::pair< uint16_t, std::vector < int16_t> >(APVn, baseline));  
  
  typename std::vector<T>::iterator strip(digis.begin() + APVn*128), lastStrip(strip + 128);
  
  int counter = 0;
  while (strip < lastStrip) {
    *strip = baseline[counter];
    counter++;
    strip++;
  }
  
}



//Baseline calculation implementation =======================================
//===========================================================================

bool inline SiStripAPVRestorer::FlatRegionsFinder(std::vector<int16_t>& adcs, DigiMap& smoothedpoints, float median, uint16_t APVn ){
  SiStripNoises::Range detNoiseRange = noiseHandle->getRange(detId_);
  
  DigiMap consecpoints;
  DigiMapIter itConsecpoints, itSmoothedpoints;
  consecpoints.erase(consecpoints.begin(), consecpoints.end());
  smoothedpoints.erase(smoothedpoints.begin(), smoothedpoints.end());
  
  
  //uint32_t idToLook = 369120382;  //to be removed
  //============================= Height above local minimum ===============================                    
  std::vector<float> adcsLocalMinSubtracted;
  adcsLocalMinSubtracted.clear();
  adcsLocalMinSubtracted.insert(adcsLocalMinSubtracted.begin(), 128,0);
  for(uint32_t istrip=0; istrip<128; ++istrip) {
    float localmin = 999.9;		
    for(uint16_t jstrip=std::max(0,(int)(istrip-nSmooth_/2)); jstrip<std::min(128,(int)(istrip+nSmooth_/2)); ++jstrip) {
      float nextvalue = adcs[jstrip];
      if(nextvalue < localmin) localmin=nextvalue;			
    }
    adcsLocalMinSubtracted[istrip] = adcs[istrip] - localmin;
  }
  
  
  //============================= Find regions with stable slopes ========================
  std::vector<uint16_t> nConsStrip;
  nConsStrip.clear();
  
  //Creating maps with all the neighborhood strip and putting in a nCosntStip vector how many we have
  uint16_t consecStrips=0;
  for(uint32_t istrip=0; istrip<128; ++istrip) {    
    int16_t adc = adcs[istrip]; 
    
    
    if( adcsLocalMinSubtracted[istrip] < nSigmaNoiseDerTh_ * (float)noiseHandle->getNoiseFast(istrip+APVn*128,detNoiseRange)
	&& adc - median < hitStripThreshold_) {
      //&& std::abs(adc - adcsLocalMinSubtracted[istrip]) < hitStripThreshold_) {           //count of many consecutive strips
      consecpoints.insert(consecpoints.end(), std::pair<uint16_t, int16_t >(istrip, adc));
      ++consecStrips;
      
    } else if (consecStrips >0){
      nConsStrip.push_back(consecStrips);
      consecStrips = 0;
    }
    
    
  }     		
  
  //to cope with the last flat region of the APV
  if(consecStrips >0) nConsStrip.push_back(consecStrips);
  
  
  //removing from the map the fist and last points in wide flat regions and erasing from the map too small regions
  itConsecpoints = consecpoints.begin();
  float MinSmoothValue=2000., MaxSmoothValue=0.;
  for(std::vector<uint16_t>::iterator itnConsStrip = nConsStrip.begin(); itnConsStrip < nConsStrip.end(); ++itnConsStrip){
    
    consecStrips = *itnConsStrip;
    if(consecStrips >=consecThreshold_){
      ++itConsecpoints;  //skipping first point
      uint16_t nFirstStrip = itConsecpoints->first;
      uint16_t nLastStrip;
      float smoothValue = 0.0;
      float stripCount =1;
      for(uint16_t n =0; n < consecStrips-2; ++n){
	smoothValue += itConsecpoints->second;
	if(stripCount == consecThreshold_){
	  smoothValue /= (float)stripCount;
	  nLastStrip = nFirstStrip + stripCount -1;				                    
	  smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nFirstStrip, smoothValue));
	  smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nLastStrip, smoothValue));
	  if(smoothValue > MaxSmoothValue) MaxSmoothValue = smoothValue;
	  if(smoothValue < MinSmoothValue) MinSmoothValue = smoothValue;
	  nFirstStrip = nLastStrip+1;
	  smoothValue=0;
	  stripCount=0;
	}
	++stripCount;
	++itConsecpoints;
      }
      ++itConsecpoints;  //and putting the pointer to the new seies of point 
      
      if(smoothValue>0){
	--stripCount;
	smoothValue /= (float)(stripCount);
	nLastStrip = nFirstStrip + stripCount -1;
	smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nFirstStrip, smoothValue));
	smoothedpoints.insert(smoothedpoints.end(), std::pair<uint16_t, int16_t >(nLastStrip, smoothValue));
	if(smoothValue > MaxSmoothValue) MaxSmoothValue = smoothValue;
	if(smoothValue < MinSmoothValue) MinSmoothValue = smoothValue;
      }
    } else{
      for(int n =0; n< consecStrips ; ++n) ++itConsecpoints;
      
    }
  }
  	
  if( (MaxSmoothValue-MinSmoothValue) > distortionThreshold_) return false;
  return true;
}


void inline SiStripAPVRestorer::BaselineFollower(DigiMap& smoothedpoints, std::vector<int16_t>& baseline, float median){
  
  baseline.clear();
  DigiMapIter itSmoothedpoints;
  
  //if not enough points
  if(smoothedpoints.size() < minStripsToFit_){
    baseline.clear();
    baseline.insert(baseline.begin(),128, median);
  } else {
    baseline.insert(baseline.begin(),128, 0);  
    
    DigiMapIter itSmoothedpointsBegin, itSmoothedpointsEnd;
    itSmoothedpointsBegin = smoothedpoints.begin();
    itSmoothedpointsEnd = --(smoothedpoints.end());
    
				
    uint16_t firstStripFlat = itSmoothedpointsBegin->first;
    uint16_t lastStripFlat = itSmoothedpointsEnd->first;
    int16_t firstStipFlatADC= itSmoothedpointsBegin->second;
    int16_t lastStipFlatADC= itSmoothedpointsEnd->second;
    
    //adding here the costant line at the extremities 
    baseline.erase(baseline.begin(), baseline.begin()+firstStripFlat);
    baseline.insert(baseline.begin(), firstStripFlat, firstStipFlatADC);
    
    baseline.erase(baseline.begin()+lastStripFlat, baseline.end());
    baseline.insert(baseline.end(), 128 - lastStripFlat, lastStipFlatADC);
    
    
    //IMPORTANT: the itSmoothedpointsEnd should be at least smaller than smoothedpoints.end() -1
    for(itSmoothedpoints = itSmoothedpointsBegin; itSmoothedpoints != itSmoothedpointsEnd; ++itSmoothedpoints){  
      DigiMapIter itSmoothedpointsNext = itSmoothedpoints;
      ++itSmoothedpointsNext;
      float strip1 = itSmoothedpoints->first;
      float strip2 = itSmoothedpointsNext->first;
      float adc1 = itSmoothedpoints->second;
      float adc2 = itSmoothedpointsNext->second;
      
      baseline[strip1] = adc1;
      baseline[strip2] = adc2;
      float m = (adc2 -adc1)/(strip2 -strip1);
      uint16_t itStrip = strip1 +1;
      float stripadc = adc1 + m; 
      while(itStrip < strip2){
	baseline[itStrip] = stripadc;
	++itStrip;
	stripadc+=m;
      }
      
    }
    
  }
}


//CM subtracion method implementation ================================
//====================================================================
/*
template<typename T>
inline
float SiStripAPVRestorer::median( std::vector<T>& sample) {
  
  typename std::vector<T>::iterator mid = sample.begin() + sample.size()/2;
  std::nth_element(sample.begin(), mid, sample.end());
  if( sample.size() & 1 ) //odd size
    return *mid;
  return ( *std::max_element(sample.begin(), mid) + *mid ) / 2.;
}
*/







//Other method implementation ==============================================
//==========================================================================

void SiStripAPVRestorer::fixAPVsCM(edm::DetSet<SiStripProcessedRawDigi>& cmdigis) {
  
  // cmdigis should be the same size as apvFlags_
  // otherwise something pathological has happened and we do nothing
  if ( cmdigis.size() != apvFlags_.size() ) return;
  
  edm::DetSet<SiStripProcessedRawDigi>::iterator cm_iter = cmdigis.begin();
  std::vector<std::string>::const_iterator apvf_iter = apvFlags_.begin();
  
  // No way to change the adc value of a SiStripProcessedRawDigi
  // so we just extract the values, clear the DetSet, and
  // replace with the proper values.
  
  std::vector<float> cmvalues;
  for( ; cm_iter != cmdigis.end(); ++cm_iter  ) cmvalues.push_back( (*cm_iter).adc() );
  cmdigis.clear();
  
  std::vector<float>::const_iterator cmv_iter = cmvalues.begin();
  while( apvf_iter != apvFlags_.end() ){
    if( *apvf_iter != "") {
      //std::cout << "  apvFlag was " << *apvf_iter << std::endl;
      //std::cout << "  baseline was " << *cmv_iter << std::endl;
      cmdigis.push_back( SiStripProcessedRawDigi( -999.) );
    }
    else
      cmdigis.push_back( SiStripProcessedRawDigi( *cmv_iter ) );
    apvf_iter++;
    cmv_iter++;
  }
}

void SiStripAPVRestorer::LoadMeanCMMap(edm::Event& iEvent){
  if(useRealMeanCM_){  
    edm::Handle< edm::DetSetVector<SiStripProcessedRawDigi> > inputCM;
    iEvent.getByLabel(inputTag_,inputCM);
    this->CreateCMMap(*inputCM);
  }
}

void SiStripAPVRestorer::CreateCMMap(const edm::DetSetVector<SiStripProcessedRawDigi>& Input){
  
  MeanCMmap_.erase(MeanCMmap_.begin(), MeanCMmap_.end());
  
  uint32_t detId_;
  edm::DetSetVector<SiStripProcessedRawDigi>::const_iterator itInput;
  edm::DetSet<SiStripProcessedRawDigi>::const_iterator itCM;
  std::vector<float> MeanCMNValue;
  
  for(itInput = Input.begin(); itInput != Input.end(); ++itInput){
    detId_ = itInput->id;
    MeanCMNValue.clear();
    for(itCM = itInput->begin(); itCM != itInput->end(); ++itCM) MeanCMNValue.push_back(itCM->adc()); 			
    MeanCMmap_.insert(std::pair<uint32_t, std::vector<float> >(detId_,MeanCMNValue));
  }
  
}
