#ifndef RKSolver_H
#define RKSolver_H

#include "TrackPropagation/RungeKutta/interface/RKSmallVector.h"
#include "TrackPropagation/RungeKutta/interface/RKDerivative.h"
#include "TrackPropagation/RungeKutta/interface/RKDistance.h"

/// ABC for Runge-Kutta solvers

template <typename T, int N>
class RKSolver {
public:

    typedef T                                   Scalar;
    typedef RKSmallVector<T,N>                  Vector;

    virtual ~RKSolver() {}

/** Advance starting state (startPar,startState) by step.
 *  The accuracy of the result should be better than eps.
 *  The accuracy is computed as the distance (using the "dist" argument)
 *  between different internal estimates of the resulting state.
 *  The "deriv" argument computes the derivatives.
 */
    virtual Vector operator()( Scalar startPar, const Vector& startState,
			       Scalar step, const RKDerivative<T,N>& deriv,
			       const RKDistance<T,N>& dist,
			       Scalar eps) = 0;


};

#endif
