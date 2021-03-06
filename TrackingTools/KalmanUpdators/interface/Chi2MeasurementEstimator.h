#ifndef CommonDet_Chi2MeasurementEstimator_H
#define CommonDet_Chi2MeasurementEstimator_H

/** \class Chi2MeasurementEstimator
 *  A Chi2 Measurement Estimator. 
 *  Computhes the Chi^2 of a TrajectoryState with a RecHit or a 
 *  BoundPlane. The TrajectoryState must have errors.
 *  Works for any RecHit dimension. Ported from ORCA.
 *
 *  $Date: 2007/05/09 13:58:19 $
 *  $Revision: 1.1.2.1 $
 *  \author todorov, cerati
 */

#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimatorBase.h"

class Chi2MeasurementEstimator : public Chi2MeasurementEstimatorBase {
public:

  /** Construct with cuts on chi2 and nSigma.
   *  The cut on Chi2 is used to define the acceptance of RecHits.
   *  The errors of the trajectory state are multiplied by nSigma 
   *  to define acceptance of BoundPlane and maximalLocalDisplacement.
   */
  explicit Chi2MeasurementEstimator(double maxChi2, double nSigma = 3.) : 
    Chi2MeasurementEstimatorBase( maxChi2, nSigma) {}

  virtual std::pair<bool,double> estimate(const TrajectoryStateOnSurface&,
				     const TransientTrackingRecHit&) const;
  template <unsigned int D> std::pair<bool,double> estimate(const TrajectoryStateOnSurface&,
				     const TransientTrackingRecHit&) const;

  virtual Chi2MeasurementEstimator* clone() const {
    return new Chi2MeasurementEstimator(*this);
  }

};

#endif
