#ifndef QIE_SHAPE_H
#define QIE_SHAPE_H

/** \class QieShape
    
    basic linearization function for HCAL/castor QIE
   $Author: ratnikov
   $Date: 2008/02/15 15:31:57 $
   $Revision: 1.1 $
*/
class QieShape {
 public:
  QieShape (const double fAdcShape [32], const double fAdcBin [32]);
  // center of the nominal linearized QIE bin
  double linearization (int fAdc) const {return mLinearization [fAdc];}
  // width of the nominal linearized bin
  double binSize (int fAdc) const {return mBinSize [fAdc];}
 private:
  double mLinearization [128];
  double mBinSize [128];
};

#endif
