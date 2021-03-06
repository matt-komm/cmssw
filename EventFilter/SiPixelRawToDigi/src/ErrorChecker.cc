#include "EventFilter/SiPixelRawToDigi/interface/ErrorChecker.h"

#include "CondFormats/SiPixelObjects/interface/SiPixelFrameConverter.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/FEDRawData/interface/FEDHeader.h"
#include "DataFormats/FEDRawData/interface/FEDTrailer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <bitset>
#include <sstream>
#include <iostream>

using namespace std;
using namespace edm;
using namespace sipixelobjects;

const int CRC_bits = 1;
const int LINK_bits = 6;
const int ROC_bits  = 5;
const int DCOL_bits = 5;
const int PXID_bits = 8;
const int ADC_bits  = 8;

const int CRC_shift = 2;
const int ADC_shift  = 0;
const int PXID_shift = ADC_shift + ADC_bits;
const int DCOL_shift = PXID_shift + PXID_bits;
const int ROC_shift  = DCOL_shift + DCOL_bits;
const int LINK_shift = ROC_shift + ROC_bits;

const uint32_t dummyDetId = 0xffffffff;

const ErrorChecker::Word64 CRC_mask = ~(~ErrorChecker::Word64(0) << CRC_bits);
const ErrorChecker::Word32 ERROR_mask = ~(~ErrorChecker::Word32(0) << ROC_bits);
const ErrorChecker::Word32 LINK_mask = ~(~ErrorChecker::Word32(0) << LINK_bits);
const ErrorChecker::Word32 ROC_mask  = ~(~ErrorChecker::Word32(0) << ROC_bits);
const ErrorChecker::Word32 DCOL_mask = ~(~ErrorChecker::Word32(0) << DCOL_bits);
const ErrorChecker::Word32 PXID_mask = ~(~ErrorChecker::Word32(0) << PXID_bits);


ErrorChecker::ErrorChecker() {

  includeErrors = false;
}

void ErrorChecker::setErrorStatus(bool ErrorStatus)
{
  includeErrors = ErrorStatus;
}

bool ErrorChecker::checkCRC(bool& errorsInEvent, int fedId, const Word64* trailer, Errors& errors)
{
  int CRC_BIT = (*trailer >> CRC_shift) & CRC_mask;
  if (CRC_BIT == 0) return true;
  errorsInEvent = true;
  if (includeErrors) {
    int errorType = 39;
    SiPixelRawDataError error(*trailer, errorType, fedId);
    errors[dummyDetId].push_back(error);
  }
  return false;
}

bool ErrorChecker::checkHeader(bool& errorsInEvent, int fedId, const Word64* header, Errors& errors)
{
  FEDHeader fedHeader( reinterpret_cast<const unsigned char*>(header));
  if ( !fedHeader.check() ) return false; // throw exception?
  if ( fedHeader.sourceID() != fedId) { 
    LogDebug("PixelDataFormatter::interpretRawData, fedHeader.sourceID() != fedId")
      <<", sourceID = " <<fedHeader.sourceID()
      <<", fedId = "<<fedId<<", errorType = 32"; 
    errorsInEvent = true;
    if (includeErrors) {
      int errorType = 32;
      SiPixelRawDataError error(*header, errorType, fedId);
      errors[dummyDetId].push_back(error);
    }
  }
  return fedHeader.moreHeaders();
}

bool ErrorChecker::checkTrailer(bool& errorsInEvent, int fedId, int nWords, const Word64* trailer, Errors& errors)
{
  FEDTrailer fedTrailer(reinterpret_cast<const unsigned char*>(trailer));
  if ( !fedTrailer.check()) { 
    if(includeErrors) {
      int errorType = 33;
      SiPixelRawDataError error(*trailer, errorType, fedId);
      errors[dummyDetId].push_back(error);
    }
    errorsInEvent = true;
    LogError("PixelDataFormatter::interpretRawData, fedTrailer.check: ")
      <<"fedTrailer.check failed, Fed: " << fedId << ", errorType = 33";
    return false; 
  } 
  if ( fedTrailer.lenght()!= nWords) {
    LogError("PROBLEM in PixelDataFormatter,  fedTrailer.lenght()!= nWords !!")<< " Fed: " << fedId << ", errorType = 34";
    errorsInEvent = true;
    if(includeErrors) {
      int errorType = 34;
      SiPixelRawDataError error(*trailer, errorType, fedId);
      errors[dummyDetId].push_back(error);
    }
  }
  return fedTrailer.moreTrailers();
}

bool ErrorChecker::checkROC(bool& errorsInEvent, int fedId, const SiPixelFrameConverter* converter, Word32& errorWord, Errors& errors)
{
 int errorType = (errorWord >> ROC_shift) & ERROR_mask;

 switch (errorType) {
    case(25) : {
     LogDebug("")<<"  invalid ROC=25 found (errorType=25)";
     errorsInEvent = true;
     break;
   }
   case(26) : {
     //LogDebug("")<<"  gap word found (errorType=26)";
     return false;
   }
   case(27) : {
     //LogDebug("")<<"  dummy word found (errorType=27)";
     return false;
   }
   case(28) : {
     LogDebug("")<<"  error fifo nearly full (errorType=28)";
     errorsInEvent = true;
     break;
   }
   case(29) : {
     LogDebug("")<<"  timeout on a channel (errorType=29)";
     errorsInEvent = true;
     break;
   }
   case(30) : {
     LogDebug("")<<"  TBM error trailer (errorType=30)";
     errorsInEvent = true;
     break;
   }
   case(31) : {
     LogDebug("")<<"  event number error (errorType=31)";
     errorsInEvent = true;
     break;
   }
   default: return true;
 };

 if(includeErrors) {
   SiPixelRawDataError error(errorWord, errorType, fedId);
   uint32_t detId;
   detId = errorDetId(converter, errorType, errorWord);
   errors[detId].push_back(error);
 }
 return false;
}

void ErrorChecker::conversionError(int fedId, const SiPixelFrameConverter* converter, int status, Word32& errorWord, Errors& errors)
{
  switch (status) {
  case(1) : {
    LogDebug("ErrorChecker::conversionError") << " Fed: " << fedId << "  invalid channel Id (errorType=35)";
    if(includeErrors) {
      int errorType = 35;
      SiPixelRawDataError error(errorWord, errorType, fedId);
      uint32_t detId = errorDetId(converter, errorType, errorWord);
      errors[detId].push_back(error);
    }
    break;
  }
  case(2) : {
    LogDebug("ErrorChecker::conversionError")<< " Fed: " << fedId << "  invalid ROC Id (errorType=36)";
    if(includeErrors) {
      int errorType = 36;
      SiPixelRawDataError error(errorWord, errorType, fedId);
      uint32_t detId = errorDetId(converter, errorType, errorWord);
      errors[detId].push_back(error);
    }
    break;
  }
  case(3) : {
    LogDebug("ErrorChecker::conversionError")<< " Fed: " << fedId << "  invalid dcol/pixel value (errorType=37)";
    if(includeErrors) {
      int errorType = 37;
      SiPixelRawDataError error(errorWord, errorType, fedId);
      uint32_t detId = errorDetId(converter, errorType, errorWord);
      errors[detId].push_back(error);
    }
    break;
  }
  case(4) : {
    LogDebug("ErrorChecker::conversionError")<< " Fed: " << fedId << "  dcol/pixel read out of order (errorType=38)";
    if(includeErrors) {
      int errorType = 38;
      SiPixelRawDataError error(errorWord, errorType, fedId);
      uint32_t detId = errorDetId(converter, errorType, errorWord);
      errors[detId].push_back(error);
    }
    break;
  }
  default: LogDebug("ErrorChecker::conversionError")<<"  cabling check returned unexpected result, status = "<< status;
  };
}

// this function finds the detId for an error word, which cannot be processed in word2digi
uint32_t ErrorChecker::errorDetId(const SiPixelFrameConverter* converter, 
    int errorType, const Word32 & word) const
{
  if (!converter) return dummyDetId;

  ElectronicIndex cabling;

  switch (errorType) {
    case  30 : case  31: case  36: {
      // set dummy values for cabling just to get detId from link if in Barrel
      cabling.dcol = 0;
      cabling.pxid = 2;
      cabling.roc  = 1;
      cabling.link = (word >> LINK_shift) & LINK_mask;  

      DetectorIndex detIdx;
      int status = converter->toDetector(cabling, detIdx);
      if (status) break;
      
      if(DetId(detIdx.rawId).subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel)) return detIdx.rawId;
      break;
    }
    case  37 : case  38: {
      cabling.dcol = 0;
      cabling.pxid = 2;
      cabling.roc  = (word >> ROC_shift) & ROC_mask;
      cabling.link = (word >> LINK_shift) & LINK_mask;

      DetectorIndex detIdx;
      int status = converter->toDetector(cabling, detIdx);
      if (status) break;

      return detIdx.rawId;
      break;
    }
  default : break;
  };
  return dummyDetId;
}
