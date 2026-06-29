//$Id: UnscentedKalmanFilter.hpp 1398 2011-04-21 20:39:37Z  $
//------------------------------------------------------------------------------
//                         UnscentedKalmanFilter
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// Copyright (c) 2002-2026 United States Government as represented by the
// Administrator of The National Aeronautics and Space Administration.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// You may not use this file except in compliance with the License. 
// You may obtain a copy of the License at:
// http://www.apache.org/licenses/LICENSE-2.0 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
// express or implied.   See the License for the specific language
// governing permissions and limitations under the License.
//
// Developed jointly by NASA/GSFC and Thinking Systems, Inc. under contract
// number NNG06CA54C
//
// Author: Darrel J. Conway, Thinking Systems, Inc.
// Created: 2009/09/03
// Modified by: Jamie LaPointe, University of Arizona, Spring 2016
/**
 * Definition for a basic extended Kalman filter
 */
//------------------------------------------------------------------------------


#ifndef UnscentedKalmanFilter_hpp
#define UnscentedKalmanFilter_hpp

#include "kalman_defs.hpp"
#include "SeqEstimator.hpp"
#include "GmatState.hpp"   // ukfPriorState
#include "GmatTime.hpp"    // ukfPriorEpochGT

class PropSetup;   // dedicated sigma-point propagator


/**
 * This implementation provides an unscented Kalman filter (UKF).
 *
 * SPIKE 1 SCAFFOLD
 * ----------------
 * This class is presently a clone of the ExtendedKalmanFilter math, registered
 * under its own type name so that "Create UnscentedKalmanFilter" scripts build a
 * working, runnable filter.  Its current numerical behavior is identical to the
 * EKF.  The purpose of Spike 1 is to prove that the sigma-point time update can
 * re-seed and re-propagate the spacecraft state through GMAT's propagation
 * subsystem.  That experiment is driven from PropagateSigmaPoints(), invoked by
 * the overridden TimeUpdate().  Once the propagation mechanism is validated, the
 * EKF-derived measurement update (ComputeObs/ComputeGain via the analytic H) is
 * replaced by the unscented-transform measurement update, and the STM-based
 * covariance map in TimeUpdate() is removed.
 *
 * Target form: Square-Root UKF (Van der Merwe & Wan), which reuses the existing
 * sqrtP / thinQR plumbing inherited from the EKF.
 */
class KALMAN_API UnscentedKalmanFilter : public SeqEstimator
{
public:
   UnscentedKalmanFilter(const std::string name);
   virtual ~UnscentedKalmanFilter();
   UnscentedKalmanFilter(const UnscentedKalmanFilter &ukf);
   UnscentedKalmanFilter&      operator=(const UnscentedKalmanFilter &ukf);

   GmatBase*               Clone() const;
   virtual void            Copy(const GmatBase*);

   // GmatBase parameter access for the UKF tuning constants
   virtual std::string  GetParameterText(const Integer id) const;
   virtual Integer      GetParameterID(const std::string &str) const;
   virtual Gmat::ParameterType
                        GetParameterType(const Integer id) const;
   virtual std::string  GetParameterTypeString(const Integer id) const;
   virtual bool         IsParameterReadOnly(const Integer id) const;
   virtual Real         GetRealParameter(const Integer id) const;
   virtual Real         SetRealParameter(const Integer id, const Real value);

protected:
   /// The measurement Residuals (O-C)
   Rvector                 yi;
   /// The Kalman gain
   Rmatrix                 kalman;

   // Factorization variables
   Rmatrix                 sqrtP;
   Rmatrix                 sqrtPupdate;

   /// Unscented innovation covariance Pyy, carried from the gain step to the
   /// covariance update (P = pBar - K Pyy K').
   Rmatrix                 ukfPyy;

   /// UKF sigma-point spread parameter (small, e.g. 1e-3)
   Real                    alpha;
   /// UKF prior-knowledge parameter (2.0 is optimal for Gaussian)
   Real                    beta;
   /// UKF secondary scaling parameter (often 0 or 3 - n)
   Real                    kappa;

   virtual void            CompleteInitialization();
   virtual void            Estimate();

   virtual void            FilterUpdate();

   virtual void            TimeUpdate();

   virtual std::string     DataFileCovHeader() const;
   virtual void            WriteCovarianceToDataFile();
   virtual void            ReadCovarianceFromDataFile(StringArray header, StringArray restartData,
                                                      UnsignedInt firstStateIndex, IntegerArray stateColumnNum);
   virtual void            UpdateCovarianceNominalValues(RealArray prevEpsilonConversions);

   virtual void            SqrtCovarianceEpsilonConversion(Rmatrix& sqrtCov);
   virtual void            SetCovariance(const Rmatrix& cov);

   virtual void            UpdateCov();

   // ---- Unscented transform machinery ----
   /// Generate the 2n+1 sigma points from (mean, sqrtCov) into sigmaPts
   void                    GenerateSigmaPoints(const Rvector &mean,
                                               const Rmatrix &sqrtCov,
                                               std::vector<Rvector> &sigmaPts) const;
   /// Re-seed and propagate each seed sigma point from startEpoch forward by dt
   /// through sigmaProp, returning the propagated estimation states in propPts.
   void                    PropagateSigmaPoints(const std::vector<Rvector> &seeds,
                                                const GmatTime &startEpoch, Real dt,
                                                std::vector<Rvector> &propPts);
   /// Unscented time update (predict): propagate the prior distribution from
   /// ukfPriorEpochGT to currentEpochGT, producing the predicted mean (written
   /// into the objects), pBar, and sqrtP.  Replaces the EKF STM covariance map.
   void                    UnscentedTimeUpdate();
   /// Refresh sqrtP as the lower-triangular Cholesky factor of a covariance P
   /// (P = sqrtP * sqrtP'), warning if P is not positive definite.
   void                    RefactorCovariance(const Rmatrix &P);
   /// Build sigmaProp: a dedicated, fully-initialized PropSetup wired to the same
   /// spacecraft esm maps to, so sigma points can be re-seeded and stepped.  The
   /// base Estimator's propagators[] are configured but never Initialize()d for
   /// standalone stepping (only RunEstimator's private clones are).
   void                    BuildSigmaPropagator();

   /// UKF weight vectors (sized 2n+1), filled in CompleteInitialization()
   RealArray               wm;   ///< mean weights
   RealArray               wc;   ///< covariance weights
   Real                    lambda;

   /// Dedicated, initialized propagator used to step sigma points (owned)
   PropSetup              *sigmaProp;

   /// UKF prior distribution carried between updates: the filter retains its own
   /// (state, epoch) because by the time TimeUpdate() runs the command has already
   /// advanced the spacecraft and reset prevUpdateEpochGT.  sqrtP (inherited) holds
   /// the matching prior square-root covariance.
   GmatState               ukfPriorState;
   GmatTime                ukfPriorEpochGT;
   bool                    ukfHasPrior;

private:
   void                    SetupMeas();
   void                    ComputeObs(UpdateInfoType &updateStat);
   void                    ComputeGain(UpdateInfoType &updateStat);
   void                    UpdateElements(UpdateInfoType &updateStat);
   void                    AdvanceEpoch();

   void                    UpdateCovarianceSimple();
   void                    UpdateCovarianceJoseph();

   /// Parameter IDs for the UKF tuning constants
   enum
   {
      ALPHA = SeqEstimatorParamCount,
      BETA,
      KAPPA,
      UnscentedKalmanFilterParamCount
   };

   const MeasurementData *calculatedMeas;
   const ObservationData *currentObs;
};

#endif /* UnscentedKalmanFilter_hpp */
