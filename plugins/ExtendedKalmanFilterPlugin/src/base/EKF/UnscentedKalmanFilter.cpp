//$Id: UnscentedKalmanFilter.cpp 1398 2011-04-21 20:39:37Z  $
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
//
// Process Noise added by:
// Author: Jamie J. LaPointe, University of Arizona
// Modifed: 2016/05/09
//
/**
 * A simple extended Kalman filter
 */
//------------------------------------------------------------------------------

#include "UnscentedKalmanFilter.hpp"
#include "EstimatorException.hpp"
#include "MessageInterface.hpp"
#include "StringUtil.hpp"
#include "UtilityException.hpp"
#include "GmatConstants.hpp"     // GmatTimeConstants (Spike 1)
#include "PropSetup.hpp"         // PropSetup (Spike 1)
#include "ODEModel.hpp"          // ODEModel re-seed API (Spike 1)
#include "Propagator.hpp"        // Propagator::Step (Spike 1)
#include "PropagationStateManager.hpp" // sigma-point PropSetup setup (Spike 1)
#include <cmath>
#include <limits>


//#define DEBUG_ESTIMATION
//#define DEBUG_JOSEPH
//#define DEBUG_ESTIMATION_COVARIACE_PROP

// Spike 1: enable to run the sigma-point re-propagation experiment inside the
// time update.  Off by default so UnscentedKalmanFilter runs as a clean EKF
// baseline; turning this on perturbs and restores the live propagator while it
// logs the recombined sigma-point mean for comparison.
#define DEBUG_UKF_SPIKE

//------------------------------------------------------------------------------
// UnscentedKalmanFilter(const std::string name)
//------------------------------------------------------------------------------
/**
 * Default constructor
 *
 * @param name The name of the new instance
 */
//------------------------------------------------------------------------------
UnscentedKalmanFilter::UnscentedKalmanFilter(const std::string name) :
   SeqEstimator  ("UnscentedKalmanFilter", name),
   alpha         (1.0e-3),
   beta          (2.0),
   kappa         (0.0),
   lambda        (0.0),
   spikeCount    (0),
   sigmaProp     (NULL),
   calculatedMeas(0),
   currentObs(0)
{
   objectTypeNames.push_back("UnscentedKalmanFilter");

   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage(" EKF default constructor: stateSize = %o, "
            "measSize = %o\n", stateSize, measSize);
   #endif
}

//------------------------------------------------------------------------------
// ~UnscentedKalmanFilter()
//------------------------------------------------------------------------------
/**
 * Destructor
 */
//------------------------------------------------------------------------------
UnscentedKalmanFilter::~UnscentedKalmanFilter()
{
   // Spike 1: release the dedicated sigma-point propagator (NULL outside the spike)
   if (sigmaProp != NULL)
      delete sigmaProp;
}


//------------------------------------------------------------------------------
// UnscentedKalmanFilter(const UnscentedKalmanFilter & ekf) :
//------------------------------------------------------------------------------
/**
 * Copy constructor
 *
 * @param ekf The instance used to configure this instance
 */
//------------------------------------------------------------------------------
UnscentedKalmanFilter::UnscentedKalmanFilter(const UnscentedKalmanFilter & ukf) :
   SeqEstimator  (ukf),
   alpha         (ukf.alpha),
   beta          (ukf.beta),
   kappa         (ukf.kappa),
   lambda        (ukf.lambda),
   spikeCount    (0),
   sigmaProp     (NULL),
   calculatedMeas(0),
   currentObs(0)
{
}


//------------------------------------------------------------------------------
// UnscentedKalmanFilter& operator=(const UnscentedKalmanFilter &ekf)
//------------------------------------------------------------------------------
/**
 * Assignment operator
 *
 * @param ekf The instance used to configure this instance
 *
 * @return this instance, configured to match ekf.
 */
//------------------------------------------------------------------------------
UnscentedKalmanFilter& UnscentedKalmanFilter::operator=(const UnscentedKalmanFilter &ukf)
{
   if (this != &ukf)
   {
      SeqEstimator::operator=(ukf);
      measSize = ukf.measSize;
      alpha    = ukf.alpha;
      beta     = ukf.beta;
      kappa    = ukf.kappa;
      lambda   = ukf.lambda;
   }

   return *this;
}


//------------------------------------------------------------------------------
// GmatBase parameter access for the UKF tuning constants
//------------------------------------------------------------------------------
std::string UnscentedKalmanFilter::GetParameterText(const Integer id) const
{
   switch (id)
   {
   case ALPHA:  return "Alpha";
   case BETA:   return "Beta";
   case KAPPA:  return "Kappa";
   default:     return SeqEstimator::GetParameterText(id);
   }
}

Integer UnscentedKalmanFilter::GetParameterID(const std::string &str) const
{
   if (str == "Alpha")  return ALPHA;
   if (str == "Beta")   return BETA;
   if (str == "Kappa")  return KAPPA;
   return SeqEstimator::GetParameterID(str);
}

Gmat::ParameterType UnscentedKalmanFilter::GetParameterType(const Integer id) const
{
   if (id >= ALPHA && id < UnscentedKalmanFilterParamCount)
      return Gmat::REAL_TYPE;
   return SeqEstimator::GetParameterType(id);
}

std::string UnscentedKalmanFilter::GetParameterTypeString(const Integer id) const
{
   if (id >= ALPHA && id < UnscentedKalmanFilterParamCount)
      return GmatBase::PARAM_TYPE_STRING[GetParameterType(id)];
   return SeqEstimator::GetParameterTypeString(id);
}

bool UnscentedKalmanFilter::IsParameterReadOnly(const Integer id) const
{
   if (id >= ALPHA && id < UnscentedKalmanFilterParamCount)
      return false;
   return SeqEstimator::IsParameterReadOnly(id);
}

Real UnscentedKalmanFilter::GetRealParameter(const Integer id) const
{
   switch (id)
   {
   case ALPHA:  return alpha;
   case BETA:   return beta;
   case KAPPA:  return kappa;
   default:     return SeqEstimator::GetRealParameter(id);
   }
}

Real UnscentedKalmanFilter::SetRealParameter(const Integer id, const Real value)
{
   switch (id)
   {
   case ALPHA:  alpha = value;  return alpha;
   case BETA:   beta  = value;  return beta;
   case KAPPA:  kappa = value;  return kappa;
   default:     return SeqEstimator::SetRealParameter(id, value);
   }
}


//------------------------------------------------------------------------------
// GmatBase* Clone() const
//------------------------------------------------------------------------------
/**
 * Object cloner
 *
 * @return A clone of this object
 */
//------------------------------------------------------------------------------
GmatBase* UnscentedKalmanFilter::Clone() const
{
   return new UnscentedKalmanFilter(*this);
}

//---------------------------------------------------------------------------
//  void Copy(const GmatBase* orig)
//---------------------------------------------------------------------------
/**
 * Sets this object to match another one.
 *
 * @param orig The original that is being copied.
 */
//---------------------------------------------------------------------------
void  UnscentedKalmanFilter::Copy(const GmatBase* orig)
{
   operator=(*((UnscentedKalmanFilter*)(orig)));
}



//------------------------------------------------------------------------------
// protected methods
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// void CompleteInitialization()
//------------------------------------------------------------------------------
/**
 * Prepares the estimator for a run
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::CompleteInitialization()
{
   SeqEstimator::CompleteInitialization();
   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage(" EKF CompleteInitialization: stateSize = %o, "
            "measSize = %o\n", stateSize, measSize);
   #endif

   Integer size = stateCovariance->GetDimension();
   if (size != (Integer)stateSize)
   {
      throw EstimatorException("In UnscentedKalmanFilter::Estimate(), the "
            "covariance matrix is not sized correctly!!!");
   }

   I = Rmatrix::Identity(stateSize);

   sqrtPupdate.SetSize(stateSize, stateSize);
   if (!sqrtP.IsSized())
   {
      sqrtP.SetSize(stateSize, stateSize);
      cf.Factor(*(stateCovariance->GetCovariance()), sqrtP);
      sqrtP = sqrtP.Transpose();
   }

   currentObs =  measManager.GetObsData();
   if (currentObs == NULL)
   {
      throw EstimatorException("Error: No observation data was used for estimation.\n");
   }

   Rmatrix outSqrtCov = sqrtP;
   SqrtCovarianceEpsilonConversion(outSqrtCov);
   updateStats[0].sqrtCov.SetSize(outSqrtCov.GetNumRows(), outSqrtCov.GetNumColumns());
   updateStats[0].sqrtCov = outSqrtCov;

   prevUpdateEpochGT = currentEpochGT;

   // ---- UKF sigma-point weights (Van der Merwe scaling) ----
   // lambda = alpha^2 (n + kappa) - n
   const Real n = (Real)stateSize;
   lambda = alpha * alpha * (n + kappa) - n;

   const UnsignedInt numSigma = 2 * stateSize + 1;
   wm.assign(numSigma, 0.0);
   wc.assign(numSigma, 0.0);
   wm[0] = lambda / (n + lambda);
   wc[0] = wm[0] + (1.0 - alpha * alpha + beta);
   for (UnsignedInt i = 1; i < numSigma; ++i)
   {
      wm[i] = 0.5 / (n + lambda);
      wc[i] = wm[i];
   }

   spikeCount = 0;

#ifdef DEBUG_UKF_SPIKE
   MessageInterface::ShowMessage("[UKF Spike] CompleteInitialization: n=%d, "
      "alpha=%g, beta=%g, kappa=%g, lambda=%g, %d sigma points\n",
      stateSize, alpha, beta, kappa, lambda, numSigma);

   // Build the dedicated, initialized propagator used to step sigma points.
   BuildSigmaPropagator();
#endif
}


//------------------------------------------------------------------------------
// void Estimate()
//------------------------------------------------------------------------------
/**
 * Implements the time update, compute, and orbit updates for the estimator
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::Estimate()
{
   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("\n\n--UnscentedKalmanFilter::Estimate----\n");
      MessageInterface::ShowMessage("Current covariance:\n");
      for (UnsignedInt i = 0; i < stateSize; ++i)
      {
         for (UnsignedInt j = 0; j < stateSize; ++j)
            MessageInterface::ShowMessage("   %.12le", stateCovariance->GetCovariance()->GetElement(i,j));
         MessageInterface::ShowMessage("\n");
      }
      MessageInterface::ShowMessage("\n");

      MessageInterface::ShowMessage("Current stm:\n");
      for (UnsignedInt i = 0; i < stateSize; ++i)
      {
         for (UnsignedInt j = 0; j < stateSize; ++j)
            MessageInterface::ShowMessage("   %.12lf", (*stm)(i,j));
         MessageInterface::ShowMessage("\n");
      }
      MessageInterface::ShowMessage("\n");

      MessageInterface::ShowMessage("Current State: [ ");
      for (UnsignedInt i = 0; i < stateSize; ++i)
         MessageInterface::ShowMessage(" %.12lf ", (*estimationState)[i]);
      MessageInterface::ShowMessage("\n");
   #endif

   UpdateInfoType updateStat;

   // setup the measurement objects for the rest of this frame of data to use
   SetupMeas();

   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("Time updated matrix \\bar P:\n");
      for (UnsignedInt i = 0; i < stateSize; ++i)
      {
         for (UnsignedInt j = 0; j < stateSize; ++j)
            MessageInterface::ShowMessage("   %.12le", pBar(i,j));
         MessageInterface::ShowMessage("\n");
      }
      MessageInterface::ShowMessage("\n");
   #endif

   // Construct the O-C data and H tilde
   ComputeObs(updateStat);

   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("hTilde:\n");
      for (UnsignedInt i = 0; i < measSize; ++i)
      {
         for (UnsignedInt j = 0; j < stateSize; ++j)
            MessageInterface::ShowMessage("   %.12lf", hTilde[i][j]);
         MessageInterface::ShowMessage("\n");
      }
      MessageInterface::ShowMessage("\n");
   #endif

   // Then the Kalman gain
   ComputeGain(updateStat);

   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("The Kalman gain is: \n");
      for (UnsignedInt i = 0; i < stateSize; ++i)
      {
         for (UnsignedInt j = 0; j < measSize; ++j)
            MessageInterface::ShowMessage("   %.12lf", kalman(i,j));
         MessageInterface::ShowMessage("\n");
      }
      MessageInterface::ShowMessage("\n");
   #endif

   // Finally, update everything
   UpdateElements(updateStat);

   // Plot residuals if set
   if (showAllResiduals)
   {
      PlotResiduals();
   }

   FillUpdateInfo(updateStat);
   updateStat.isObs = true;

   Rmatrix outSqrtCov = sqrtP;
   SqrtCovarianceEpsilonConversion(outSqrtCov);
   updateStat.sqrtCov.SetSize(outSqrtCov.GetNumRows(), outSqrtCov.GetNumColumns());
   updateStat.sqrtCov = outSqrtCov;

   updateStats.push_back(updateStat);
   BuildMeasurementLine(updateStat.measStat);
   WriteDataFile();
   AddMatlabData(updateStat.measStat);
   AddMatlabFilterData(updateStat);


   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("Updated covariance:\n");
      for (UnsignedInt i = 0; i < stateSize; ++i)
      {
         for (UnsignedInt j = 0; j < stateSize; ++j)
            MessageInterface::ShowMessage("   %.12le", stateCovariance->GetCovariance()->GetElement(i,j));
         MessageInterface::ShowMessage("\n");
      }
      MessageInterface::ShowMessage("\n");

      MessageInterface::ShowMessage("Updated State: [ ");
      for (UnsignedInt i = 0; i < stateSize; ++i)
         MessageInterface::ShowMessage(" %.12lf ", (*estimationState)[i]);
      MessageInterface::ShowMessage("\n\n---------------------\n");
   #endif

   if (textFileMode == "Verbose")
      ReportProgress();

   AdvanceEpoch();
}


//------------------------------------------------------------------------------
// void FilterUpdate()
//------------------------------------------------------------------------------
/**
 * This method performs actions common to sequential estimators to update the
 * filter information and add it to the updateStats vector
 */
 //------------------------------------------------------------------------------
void UnscentedKalmanFilter::FilterUpdate()
{
   // Update the STM
   esm.MapObjectsToSTM();
   esm.MapObjectsToVector();

   UpdateProcessNoise();
   TimeUpdate();
   (*(stateCovariance->GetCovariance())) = pBar;
   esm.MapCovariancesToObjects();

   if (currentState != CALCULATING)
   {
      // Get data for the covariance report
      UpdateInfoType updateStat;
      FillUpdateInfo(updateStat);
      updateStat.isObs = false;

      Rmatrix outSqrtCov = sqrtP;
      SqrtCovarianceEpsilonConversion(outSqrtCov);
      updateStat.sqrtCov.SetSize(outSqrtCov.GetNumRows(), outSqrtCov.GetNumColumns());
      updateStat.sqrtCov = outSqrtCov;

      // Don't want to add data to vector if predicting to an epoch that's not an anchor epoch
      if (!dontWriteDataInUpdate)
      {
         if (!(isPredicting && !hasAnchorEpoch))
         {
            WriteDataFile();
            AddMatlabFilterData(updateStat);
            updateStats.push_back(updateStat);
         }
         else if (addPredictToMatlab)
         {
            AddMatlabFilterData(updateStat);
         }
      }

      // Reset the STM
      PrepareForStep();
      esm.MapObjectsToVector();
      PropagationStateManager *psm = propagators[0]->GetPropStateManager();
      psm->MapObjectsToVector();
      resetState = true;
   }
}

//------------------------------------------------------------------------------
// void TimeUpdate()
//------------------------------------------------------------------------------
/**
 * Performs the time update of the state error covariance
 *
 * This method uses Cholesky factorization for covariance
 * This is based on section 5.7 of Brown and Hwang 4e
 */
 //------------------------------------------------------------------------------
void UnscentedKalmanFilter::TimeUpdate()
{
#ifdef DEBUG_ESTIMATION
   MessageInterface::ShowMessage("Performing time update\n");
#endif

#ifdef DEBUG_UKF_SPIKE
   // ---- SPIKE 1 EXPERIMENT ----------------------------------------------
   // Generate sigma points about the current reference state, re-seed and
   // propagate each one forward by a representative test interval, recombine the
   // weighted mean, and restore the reference so the live filter is unperturbed.
   // Run only on the first few time updates to keep the log readable; uses a
   // fixed dt because the goal is to validate the re-propagation mechanism, not
   // to drive the covariance (the EKF STM update below still owns the filter).
   {
      const Integer spikeMaxRuns = 3;
      const Real    spikeTestDt  = 60.0;   // seconds
      if (spikeCount < spikeMaxRuns)
      {
         MessageInterface::ShowMessage(
            "\n[UKF Spike] === experiment %d of %d at epoch %s ===\n",
            spikeCount + 1, spikeMaxRuns, currentEpochGT.ToString().c_str());
         PropagateSigmaPoints(spikeTestDt);
         ++spikeCount;
      }
   }
   // ---- END SPIKE 1 EXPERIMENT ------------------------------------------
#endif

#ifdef DEBUG_ESTIMATION_COVARIACE_PROP
   MessageInterface::ShowMessage("Q = \n");
   for (UnsignedInt i = 0; i < stateSize; ++i)
   {
      for (UnsignedInt j = 0; j < stateSize; ++j)
         MessageInterface::ShowMessage("   %.12le", Q(i, j));
      MessageInterface::ShowMessage("\n");
   }
   MessageInterface::ShowMessage("\n");
#endif

   /// Calculate conversion derivative matrixes
   // Calculate conversion derivative matrix [dX/dS] from Cartesian to Solve-for state
   cart2SolvMatrix = esm.CartToSolveForStateConversionDerivativeMatrix();
   // Calculate conversion derivative matrix [dS/dK] from solve-for state to Keplerian
   solv2KeplMatrix = esm.SolveForStateToKeplConversionDerivativeMatrix();

   Rmatrix dX_dS = cart2SolvMatrixPrev;
   Rmatrix dS_dX = cart2SolvMatrix.Inverse();

   Rmatrix Q_S = dS_dX * Q * dS_dX.Transpose();


#ifdef DEBUG_ESTIMATION_COVARIACE_PROP
   MessageInterface::ShowMessage("dS_dX = \n");
   for (UnsignedInt i = 0; i < stateSize; ++i)
   {
      for (UnsignedInt j = 0; j < stateSize; ++j)
         MessageInterface::ShowMessage("   %.12le", dS_dX(i, j));
      MessageInterface::ShowMessage("\n");
   }
   MessageInterface::ShowMessage("\n");
   MessageInterface::ShowMessage("dX_dS = \n");
   for (UnsignedInt i = 0; i < 6; ++i)
   {
      for (UnsignedInt j = 0; j < 6; ++j)
         MessageInterface::ShowMessage("   %.12le", dX_dS(i, j));
      MessageInterface::ShowMessage("\n");
   }
   MessageInterface::ShowMessage("\n");
   MessageInterface::ShowMessage("stm = \n");
   for (UnsignedInt i = 0; i < 6; ++i)
   {
      for (UnsignedInt j = 0; j < 6; ++j)
         MessageInterface::ShowMessage("   %.12le", (*stm)(i, j));
      MessageInterface::ShowMessage("\n");
   }
   MessageInterface::ShowMessage("\n");
   MessageInterface::ShowMessage("Q_s = \n");
   for (UnsignedInt i = 0; i < stateSize; ++i)
   {
      for (UnsignedInt j = 0; j < stateSize; ++j)
         MessageInterface::ShowMessage("   %.12le", Q_S(i, j));
      MessageInterface::ShowMessage("\n");
   }
   MessageInterface::ShowMessage("\n");
#endif

   Rmatrix stm_S = dS_dX * (*stm) * dX_dS;

   // Update offset from reference trajectory
   if (esm.HasStateOffset())
   {
      GmatState *offsetState = esm.GetStateOffset();
      Rvector xOffset(stateSize);
      xOffset.Set(offsetState->GetState(), stateSize);

      xOffset = (*stm) * xOffset;

      for (UnsignedInt i = 0; i < stateSize; ++i)
         (*offsetState)[i] = xOffset[i];
   }

   // Form perform thinQR decomposition to calculate pBar

   Rmatrix sqrtQ_T(stateSize, stateSize);

   bool hasZeroDiag = false;
   for (UnsignedInt ii = 0U; ii < stateSize; ii++)
   {
      if (Q_S(ii, ii) == 0U)
      {
         hasZeroDiag = true;
         break;
      }
   }

   if (!hasZeroDiag)
   {
      try
      {
         cf.Factor(Q_S, sqrtQ_T);
      }
      catch (UtilityException e)
      {
         throw EstimatorException("The process noise matrix is not positive definite!");
      }
   }
   else
   {
      // Remove all zero rows/columns first
      IntegerArray removedIndexes;
      IntegerArray auxVector;
      Integer numRemoved;
      Rmatrix reducedQ_S = MatrixFactorization::CompressNormalMatrix(Q_S,
         removedIndexes, auxVector, numRemoved);

      Rmatrix reducedSqrtQ_T(stateSize - numRemoved, stateSize - numRemoved);
      try
      {
         cf.Factor(reducedQ_S, reducedSqrtQ_T);
      }
      catch (UtilityException e)
      {
         throw EstimatorException("The process noise matrix is not positive definite!");
      }

      sqrtQ_T = MatrixFactorization::ExpandNormalMatrixInverse(reducedSqrtQ_T,
         auxVector, numRemoved);
   }

   Rmatrix stmP = stm_S * sqrtP;
   Rmatrix sqrtQ = sqrtQ_T.Transpose();

   #ifdef DEBUG_ESTIMATION_COVARIACE_PROP
      MessageInterface::ShowMessage("stmP = \n");
      for (UnsignedInt i = 0; i < stateSize; ++i)
      {
         for (UnsignedInt j = 0; j < stateSize; ++j)
         {
            MessageInterface::ShowMessage("  %.12le", stmP(i, j));
         }
         MessageInterface::ShowMessage("\n");
      }
      MessageInterface::ShowMessage("sqrtP = \n");
      for (UnsignedInt i = 0; i < stateSize; ++i)
      {
         for (UnsignedInt j = 0; j < stateSize; ++j)
         {
            MessageInterface::ShowMessage("  %.12le", sqrtP(i, j));
         }
         MessageInterface::ShowMessage("\n");
      }
      MessageInterface::ShowMessage("sqrtQ = \n");
      for (UnsignedInt i = 0; i < 6; ++i)
      {
         for (UnsignedInt j = 0; j < 6; ++j)
         {
            MessageInterface::ShowMessage("  %.12le", sqrtQ(i, j));
         }
         MessageInterface::ShowMessage("\n");
      }
#endif
   sqrtP = thinQR(stmP, sqrtQ);

   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("sqrtP = \n");
      for (UnsignedInt i = 0; i < stateSize; ++i)
      {
         for (UnsignedInt j = 0; j < stateSize; ++j)
         {
            MessageInterface::ShowMessage("  %.12le", sqrtP(i,j));
         }
         MessageInterface::ShowMessage("\n");
      }
   #endif

   // Warn if covariance is not positive definite
   for (UnsignedInt ii = 0U; ii < stateSize; ii++)
   {
      if (GmatMathUtil::Abs(sqrtP(ii, ii)) < 1e-16)
      {
         MessageInterface::ShowMessage("WARNING The covariance is no longer positive definite! Epoch = %s\n", currentEpochGT.ToString().c_str());
         break;
      }
   }

   pBar = sqrtP * sqrtP.Transpose();

   // make it symmetric!
   Symmetrize(pBar);

#ifdef DEBUG_ESTIMATION_COVARIACE_PROP
   MessageInterface::ShowMessage("pBar = \n");
   for (UnsignedInt i = 0; i < stateSize; ++i)
   {
      for (UnsignedInt j = 0; j < stateSize; ++j)
      {
         MessageInterface::ShowMessage("  %.12le", pBar(i, j));
      }
      MessageInterface::ShowMessage("\n");
   }
#endif
}

//------------------------------------------------------------------------------
// void SetupMeas()
//------------------------------------------------------------------------------
/**
 * This sets up the measurement information for others to use later
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::SetupMeas()
{
   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("Performing measurement setup\n");
   #endif

   modelsToAccess = measManager.GetValidMeasurementList();
   currentObs =  measManager.GetObsData();

   if (modelsToAccess.size() > 0U)
   {
      measCount = measManager.CountFeasibleMeasurements(modelsToAccess[0]);
      calculatedMeas = measManager.GetMeasurement(modelsToAccess[0]);

      // verify media correction to be in acceptable range. It is [0m, 60m] for troposphere correction and [0m, 20m] for ionosphere correction
      ValidateMediaCorrection(calculatedMeas);

      // Make correction for observation value before running data filter
      if ((iterationsTaken == 0) && (currentObs->typeName == "DSN_SeqRange" || currentObs->typeName == "DSN_PNRange"))
      {
         // value correction is only applied for DSN_SeqRange and it is only performed at the first time
         for (Integer index = 0; index < currentObs->value.size(); ++index)
            measManager.GetObsDataObject()->value[index] = ObservationDataCorrection(calculatedMeas->value[index], currentObs->value[index], currentObs->rangeModulo);
      }

      // Get pre-update covariance and symmetrize
      pBar = sqrtP * sqrtP.Transpose();
      Symmetrize(pBar);

      // Size the measurement matricies
      measSize = currentObs->value.size();

      H.SetSize(measSize, stateSize);
      yi.SetSize(measSize);
      kalman.SetSize(stateSize, measSize);
   }

   /// Calculate conversion derivative matrixes
   // Calculate conversion derivative matrix [dX/dS] from Cartesian to Solve-for state
   cart2SolvMatrix = esm.CartToSolveForStateConversionDerivativeMatrix();
   // Calculate conversion derivative matrix [dS/dK] from solve-for state to Keplerian
   solv2KeplMatrix = esm.SolveForStateToKeplConversionDerivativeMatrix();
}


//------------------------------------------------------------------------------
// void ComputeObs(UpdateInfoType &updateStat)
//------------------------------------------------------------------------------
/**
 * Computes the measurement residuals and the H-tilde matrix
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::ComputeObs(UpdateInfoType &updateStat)
{
   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("Computing obs and hTilde\n");
   #endif
   // Compute the O-C, Htilde, and Kalman gain

   // Populate measurement statistics
   FilterMeasurementInfoType measStat;
   CalculateResiduals(measStat);

   // Populate H and y
   if (modelsToAccess.size() > 0)
   {
      // Adjust computed and residual based on value of xOffset
      Rvector xOffset(stateSize);
      xOffset.Set(esm.GetEstimationStateOffset().GetState(), stateSize);

      Rvector H_x = H * xOffset;
      for (UnsignedInt k = 0; k < measStat.residual.size(); ++k)
      {
         measStat.measValue[k] += H_x(k);
         measStat.residual[k] -= H_x(k);
      }

      if (measStat.isCalculated)
      {
         for (UnsignedInt k = 0; k < measStat.residual.size(); ++k)
            yi(k) = measStat.residual[k];
      }

      // get scaled residuals
      Rmatrix R = *(GetMeasurementCovariance()->GetCovariance());

      // Keep this line for when we implement the scaled residual for the entire measurement
      // instead of for each element of the measurement:
      // measStat.scaledResid = GmatMathUtil::Sqrt(yi * (H * pBar * H.Transpose() + R).Inverse() * yi);

      // The element-by-element scaled residual calculation:
      for (UnsignedInt k = 0; k < measStat.residual.size(); ++k)
      {
         Rmatrix Rbar = H * pBar * H.Transpose() + R;
         Real sigmaVal = GmatMathUtil::Sqrt(Rbar(k, k));
         Real scaledResid = measStat.residual[k] / sigmaVal;
         measStat.scaledResid.push_back(scaledResid);
      }

   }  // end of if (modelsToAccess.size() > 0)

   GmatState currentState = esm.GetEstimationMJ2000EqCartesianStateForReport();
   for (UnsignedInt ii = 0; ii < stateSize; ii++)
      measStat.state.push_back(currentState[ii]);

   // Add state offset if not rectified
   if (esm.HasStateOffset())
   {
      GmatState xOffset = *esm.GetStateOffset();
      for (UnsignedInt ii = 0; ii < stateSize; ii++)
      {
         Real conv = GetEpsilonConversion(ii);
         measStat.state[ii] += xOffset[ii] * conv;
      }
   }

   Rmatrix outCov = pBar;
   CovarianceEpsilonConversion(outCov);
   measStat.cov.SetSize(outCov.GetNumRows(), outCov.GetNumColumns());
   measStat.cov = outCov;
   Rmatrix covVNB = GetCovarianceVNB(pBar);
   measStat.covVNB.SetSize(covVNB.GetNumRows(), covVNB.GetNumColumns());
   measStat.covVNB = covVNB;

   Rmatrix outSqrtCov = sqrtP;
   SqrtCovarianceEpsilonConversion(outSqrtCov);
   measStat.sqrtCov.SetSize(outSqrtCov.GetNumRows(), outSqrtCov.GetNumColumns());
   measStat.sqrtCov = outSqrtCov;

   measStats.push_back(measStat);

   updateStat.epoch = currentEpochGT;
   updateStat.isObs = true;
   updateStat.measStat = measStat;

   BuildMeasurementLine(measStat);
   WriteToTextFile();
}


//------------------------------------------------------------------------------
// void ComputeGain(UpdateInfoType &updateStat)
//------------------------------------------------------------------------------
/**
 * Computes the Kalman gain
 *
 * The error estimates used for error bars on the residuals plots are calculated
 * as
 *
 *    sigma = sqrt(H P H' + R)
 *
 * Since the argument of the square root is calculated as part of the Kalman
 * gain calculation, this value is also stored in this method
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::ComputeGain(UpdateInfoType &updateStat)
{
   if (updateStat.measStat.isCalculated)
   {
      #ifdef DEBUG_ESTIMATION
         MessageInterface::ShowMessage("Computing Kalman Gain\n");
      #endif

      // Set up measurement underweighting (Lear's method)
      Real sqrtScale = 1.0;
      Real posCovTraceSqrt = GmatMathUtil::Sqrt(pBar(0, 0) + pBar(1, 1) + pBar(2, 2));
      if (posCovTraceSqrt > deweightThreshold && deweightCoeff > 0)
      {
         sqrtScale = GmatMathUtil::Sqrt(1.0 + deweightCoeff);

         bool handleLeapSecond;
         GmatTime utcMjdEpoch = theTimeConverter->Convert(updateStat.measStat.epoch, TimeSystemConverter::A1MJD, TimeSystemConverter::UTCMJD,
            GmatTimeConstants::JD_JAN_5_1941, &handleLeapSecond);
         std::string utcEpoch = theTimeConverter->ConvertMjdToGregorian(utcMjdEpoch.GetMjd(), handleLeapSecond);

         MessageInterface::ShowMessage("Measurement %d of type %s at %s UTCG was underweighted. (1 sigma pos uncertainty was %g km)\n",
            updateStat.measStat.recNum, updateStat.measStat.type.c_str(), utcEpoch.c_str(), posCovTraceSqrt);
      }

      // Perform thinQR decomposition to calculate K and P
      Rmatrix R = *(GetMeasurementCovariance()->GetCovariance());
      Integer measSize = R.GetNumRows();

      Rmatrix sqrtR_T(measSize, measSize);
      cf.Factor(R, sqrtR_T);

      Rmatrix Spbar = sqrtP;
      Rmatrix Sr = sqrtR_T.Transpose();
      Rmatrix Sw = thinQR(sqrtScale*H*Spbar, Sr);

      #ifdef DEBUG_ESTIMATION
         MessageInterface::ShowMessage("Sw = \n");
         for (UnsignedInt i = 0; i < Sw.GetNumRows(); ++i)
         {
           for (UnsignedInt j = 0; j < Sw.GetNumColumns(); ++j)
           {
             MessageInterface::ShowMessage("  %.12le", Sw(i,j));
           }
           MessageInterface::ShowMessage("\n");
         }
      #endif

      #ifdef DEBUG_ESTIMATION
         MessageInterface::ShowMessage("Calculating the Kalman gain\n");
      #endif

      kalman = Spbar * Spbar.Transpose() * H.Transpose() * (Sw*Sw.Transpose()).Inverse();
      sqrtPupdate = thinQR((I - kalman * H) * Spbar, kalman*Sr);

      updateStat.measStat.kalmanGain.SetSize(kalman.GetNumRows(), kalman.GetNumColumns());
      updateStat.measStat.kalmanGain = kalman;
   }
}


//------------------------------------------------------------------------------
// void UpdateElements(UpdateInfoType &updateStat)
//------------------------------------------------------------------------------
/**
 * Updates the estimation state and covariance matrix
 *
 * Programmers can select the covariance update method at the end of this
 * method.  The resulting covariance is symmetrized before returning.
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::UpdateElements(UpdateInfoType &updateStat)
{
   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("Updating elements\n");
   #endif

   if (updateStat.measStat.editFlag == NORMAL_FLAG)
   {
      dx = kalman * yi;

      if (esm.HasStateOffset())
      {
         GmatState offsetState = esm.GetEstimationStateOffset();
         for (UnsignedInt i = 0; i < stateSize; ++i)
         {
            offsetState[i] += dx[i];
         }
         esm.SetEstimationStateOffset(offsetState);

      }
      else
      {
         // Update the state, covariances, and so forth
         estimationStateS = esm.GetEstimationState();
         for (UnsignedInt i = 0; i < stateSize; ++i)
         {
            estimationStateS[i] += dx[i];
         }

         // Convert estimation state from Keplerian to Cartesian
         esm.SetEstimationState(estimationStateS);                       // update the value of estimation state
         esm.MapVectorToObjects();
         esm.MapCovariancesToObjects();
      }

      #ifdef DEBUG_ESTIMATION
         MessageInterface::ShowMessage("Calculated state change: [");
         for (UnsignedInt i = 0; i < stateSize; ++i)
            MessageInterface::ShowMessage(" %.12lf ", dx[i]);
         MessageInterface::ShowMessage("\n");
      #endif

      // Select the method used to update the covariance here:
      // UpdateCovarianceSimple();
      // UpdateCovarianceJoseph();

      Rmatrix P2 = sqrtPupdate * sqrtPupdate.Transpose();
      sqrtP = sqrtPupdate;

      // Warn if covariance is not positive definite
      for (UnsignedInt ii = 0U; ii < stateSize; ii++)
      {
         if (GmatMathUtil::Abs(sqrtP(ii, ii)) < 1e-16)
         {
            MessageInterface::ShowMessage("WARNING The covariance is no longer positive definite! Epoch = %s\n", currentEpochGT.ToString().c_str());
            break;
         }
      }

      (*(stateCovariance->GetCovariance())) = P2;
   }
   else
      (*(stateCovariance->GetCovariance())) = sqrtP * sqrtP.Transpose();

   Symmetrize(*stateCovariance);
   informationInverse = (*(stateCovariance->GetCovariance()));
   //information = informationInverse.Inverse(COV_INV_TOL);
}


//------------------------------------------------------------------------------
// void UpdateCovarianceSimple()
//------------------------------------------------------------------------------
/**
 * Applies equation (4.7.12) to update the state error covariance matrix
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::UpdateCovarianceSimple()
{
   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("Updating covariance using simple "
            "method\n");
   #endif

   // P = (I - K * H) * Pbar
   (*(stateCovariance->GetCovariance())) = (I - (kalman * H)) * pBar;
}


//------------------------------------------------------------------------------
// void UnscentedKalmanFilter::UpdateCovarianceJoseph()
//------------------------------------------------------------------------------
/**
 * This method updates the state error covariance matrix using the method
 * developed by Bucy and Joseph, as presented in Tapley, Schutz and Born
 * eq (4.7.19)
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::UpdateCovarianceJoseph()
{
   #ifdef DEBUG_ESTIMATION
      MessageInterface::ShowMessage("Updating covariance using Joseph "
            "method\n");
   #endif

   Rmatrix *r = GetMeasurementCovariance()->GetCovariance();

   // P = (I - K * H) * Pbar * (I - K * H)^T + K * R * K^T
   (*(stateCovariance->GetCovariance())) =
         ((I - (kalman * H)) * pBar * (I - (kalman * H)).Transpose()) +
         (kalman * (*r) * kalman.Transpose());

   #ifdef DEBUG_JOSEPH
      for (UnsignedInt i = 0; i < stateSize; ++i)
      {
         for (UnsignedInt j = 0; j < stateSize; ++j)
            MessageInterface::ShowMessage("  %.12lf  ", (*(stateCovariance->GetCovariance()))(i,j));
         MessageInterface::ShowMessage("\n");
      }
      MessageInterface::ShowMessage("\n");

      throw EstimatorException("Intentional debug break!");
   #endif
}

void UnscentedKalmanFilter::AdvanceEpoch()
{
   // Reset the STM
   PrepareForStep();
   esm.MapVectorToObjects();
   esm.MapCovariancesToObjects();
   //PropagationStateManager *psm = propagators[0]->GetPropStateManager();
   //psm->MapObjectsToVector();

   for (int i = 0; i < propagators.size(); i++)
   {
      PropagationStateManager *psm = propagators[i]->GetPropStateManager();
      psm->MapObjectsToVector();
   }
   // Flag that a new current state has been loaded in the objects
   resetState = true;

   // Advance MeasMan to the next measurement and get its epoch
   bool isEndOfTable = measManager.AdvanceObservation();
   if (isEndOfTable)
   {
      currentState = CHECKINGRUN;
   }
   else
   {
      nextMeasurementEpochGT = measManager.GetEpochGT();

      // Check if rectification should begin here
      if (esm.HasStateOffset())
      {
         // Check if next measurement is after the delayed rectification span
         Real elapsedTime = (nextMeasurementEpochGT - estimationEpochGT).GetTimeInSec();

         if (GmatMathUtil::Abs(elapsedTime) > delayRectifySpan)
         {
            // Update the state with the state offset
            GmatState offsetStateS = esm.GetEstimationStateOffset();
            estimationStateS = esm.GetEstimationState();
            for (UnsignedInt i = 0; i < stateSize; ++i)
            {
               estimationStateS[i] += offsetStateS[i];
            }

            // Convert estimation state from Keplerian to Cartesian
            esm.SetEstimationState(estimationStateS);                       // update the value of estimation state
            esm.MapVectorToObjects();
            esm.MapCovariancesToObjects();

            // Zero out the state offset
            GmatState *offsetState = esm.GetStateOffset();
            for (UnsignedInt i = 0; i < stateSize; ++i)
               (*offsetState)[i] = 0.0;

            esm.SetHasStateOffset(false);

            bool handleLeapSecond;
            GmatTime utcMjdEpoch = theTimeConverter->Convert(currentEpochGT, TimeSystemConverter::A1MJD, TimeSystemConverter::UTCMJD,
               GmatTimeConstants::JD_JAN_5_1941, &handleLeapSecond);
            std::string utcEpoch = theTimeConverter->ConvertMjdToGregorian(utcMjdEpoch.GetMjd(), handleLeapSecond);

            MessageInterface::ShowMessage("Exiting DelayRectifyTimeSpan at %s UTCG\n", utcEpoch.c_str());
            WriteDataFile();
         }
      }

      // Reset nextNoiseUpdateGT if it is in the wrong direction
      Real dtNoise = (nextNoiseUpdateGT - currentEpochGT).GetTimeInSec();
      Real dtMeasurement = (nextMeasurementEpochGT - currentEpochGT).GetTimeInSec();

      // If filter has passed the noise epoch
      if (measManager.IsForward() && dtNoise < ESTTIME_ROUNDOFF * GmatTimeConstants::SECS_PER_DAY)
         nextNoiseUpdateGT.AddSeconds(processNoiseStep);
      else if (!measManager.IsForward() && dtNoise > -ESTTIME_ROUNDOFF * GmatTimeConstants::SECS_PER_DAY)
         nextNoiseUpdateGT.SubtractSeconds(processNoiseStep);

      FindTimeStep();

      #ifdef DEBUG_ESTIMATION
         MessageInterface::ShowMessage("UnscentedKalmanFilter::AdvanceEpoch CurrentEpoch = %.12lf, next "
               "epoch = %.12lf, timeStep = %.12lf\n", currentEpochGT.GetMjd(),
               nextMeasurementEpochGT.GetMjd(), timeStep);
      #endif

      // this magical number is from the Batch Estimator in its accumulating state...
      //if (currentEpoch <= (nextMeasurementEpoch+5.0e-12))
      if (nextMeasurementEpochGT >= 5.0e-12)
      {
         currentState = PROPAGATING;
      }
      else
      {
         currentState = CHECKINGRUN;  // Should this just go to FINISHED?
      }
   }
}


//------------------------------------------------------------------------------
// std::string DataFileCovHeader() const
//----------------------------------------------------------------------
/**
 * Write the name of the covariance type in header of the data file
 */
 //------------------------------------------------------------------------------
std::string UnscentedKalmanFilter::DataFileCovHeader() const
{
   return "SqrtCovariance";
}


//------------------------------------------------------------------------------
// void WriteCovarianceToDataFile()
//----------------------------------------------------------------------
/**
 * Writes the covariance to the data file in factorized form
 */
 //------------------------------------------------------------------------------
void UnscentedKalmanFilter::WriteCovarianceToDataFile()
{
   if (!sqrtP.IsSized())
   {
      sqrtP.SetSize(stateSize, stateSize);
      cf.Factor(*(stateCovariance->GetCovariance()), sqrtP);
      sqrtP = sqrtP.Transpose();
   }

   Rmatrix sqrtPOut = sqrtP;

   // Calculate conversion derivative matrix [dX/dS] from Cartesian to Solve-for state
   Rmatrix dX_dS = esm.CartToSolveForStateConversionDerivativeMatrix();
   Rmatrix cov = *(stateCovariance->GetCovariance());
   bool notCartesian = false;

   for (UnsignedInt ii = 0; ii < stateSize; ii++)
   {
      for (UnsignedInt jj = 0; jj < stateSize; jj++)
      {
         Real identityValue = (ii == jj) ? 1.0 : 0.0;

         if (dX_dS(ii, jj) != identityValue)
         {
            notCartesian = true;
            break;
         }
      }

      if (notCartesian)
         break;
   }

   if (notCartesian)
   {
      // If the solve for state is not Cartesian, convert covariance to cartesian
      // and find Cholesky decompsition
      cov = dX_dS * cov * dX_dS.Transpose();
      cf.Factor(cov, sqrtPOut);
      sqrtPOut = sqrtPOut.Transpose();
   }

   // Write covariance lower triangle
   for (UnsignedInt ii = 0; ii < stateSize; ii++)
   {
      Real conv = GetEpsilonConversion(ii);
      for (UnsignedInt jj = 0; jj <= ii; jj++)
      {
         Real value = sqrtPOut(ii, jj) * conv;

         std::string valueStr = GmatStringUtil::RealToString(value, false, true, true);
         dataFile << "," << valueStr;
      }
   }
}


//------------------------------------------------------------------------------
// void ReadCovarianceFromDataFile(StringArray header, StringArray restartData,
//                                 UnsignedInt firstStateIndex, IntegerArray stateColumnNum)
//----------------------------------------------------------------------
/**
 * Reads the covariance from the data file
 *
 * @param header The array of header elements.
 * @param restartData The array of data elements for the restart row.
 * @param firstStateIndex The column index of the first state element in restartData.
 * @param stateColumnNum The maping of the solve for states to the restart file states.
 */
 //------------------------------------------------------------------------------
void UnscentedKalmanFilter::ReadCovarianceFromDataFile(StringArray header, StringArray restartData,
                                                      UnsignedInt firstStateIndex, IntegerArray stateColumnNum)
{
   // Find column of first covariance element
   bool covFound = false;
   UnsignedInt firstCovIndex;
   for (UnsignedInt ii = 0U; ii < header.size(); ii++)
   {
      if (header[ii] == "SqrtCovariance_1_1")
      {
         covFound = true;
         firstCovIndex = ii;
         break;
      }
   }

   if (!covFound)
   {
      // If the EKF factorized covariance isn't found, try the default full covariance
      SeqEstimator::ReadCovarianceFromDataFile(header, restartData, firstStateIndex, stateColumnNum);
      return;
   }

   UnsignedInt fileStateSize = firstCovIndex - firstStateIndex;
   Integer stateAlignedIndex = -1;

   // Indicies in stateColumnNum must be in sequential order for factorized covariance to work
   for (UnsignedInt ii = 0; ii < fileStateSize && ii < stateColumnNum.size(); ii++)
   {
      if (stateColumnNum[ii] - firstStateIndex != ii)
         break;

      stateAlignedIndex = ii;
   }

   // Set covariance
   Rmatrix fileCov(fileStateSize, fileStateSize);
   UnsignedInt index = firstCovIndex;
   // Convert lower triangle array to lower triangle matrix 
   for (UnsignedInt ii = 0; ii < fileStateSize; ii++)
   {
      for (UnsignedInt jj = 0; jj <= ii; jj++)
      {
         Real value;
         GmatStringUtil::ToReal(restartData[index], value);
         fileCov(ii, jj) = value;

         index++;
      }
   }

   // Map file covariance matrix to state order of the covariance in the esm
   Covariance *cov = esm.GetCovariance();

   sqrtP.SetSize(esm.GetStateSize(), esm.GetStateSize());

   for (UnsignedInt ii = 0; ii < esm.GetStateSize(); ii++)
   {
      if (stateColumnNum[ii] >= 0)
      {
         Real conv = GetEpsilonConversion(ii);
         for (UnsignedInt jj = 0; jj <= ii; jj++)
         {
            if (stateColumnNum[jj] >= 0)
            {
               UnsignedInt idx1 = stateColumnNum[ii] - firstStateIndex;
               UnsignedInt idx2 = stateColumnNum[jj] - firstStateIndex;
               Real value = fileCov(idx1, idx2);
               value /= conv;
               sqrtP(ii, jj) = value;
            }
         }
      }
   }

   // If there are additional states being solved for that are not in the warm start file
   if (stateAlignedIndex + 1 < esm.GetStateSize())
   {
      Rmatrix covIn = fileCov * fileCov.Transpose();
      Rmatrix covInAll = (*(cov->GetCovariance()));

      for (UnsignedInt ii = 0; ii < fileStateSize && ii < stateColumnNum.size(); ii++)
      {
         if (stateColumnNum[ii] >= 0)
         {
            Real iiConv = GetEpsilonConversion(ii);
            for (UnsignedInt jj = 0; jj < fileStateSize && jj < stateColumnNum.size(); jj++)
            {
               if (stateColumnNum[jj] >= 0)
               {
                  Real jjConv = GetEpsilonConversion(jj);
                  UnsignedInt idx1 = stateColumnNum[ii] - firstStateIndex;
                  UnsignedInt idx2 = stateColumnNum[jj] - firstStateIndex;
                  Real value = covIn(idx1, idx2);
                  value /= iiConv * jjConv;
                  covInAll(ii, jj) = value;
               }
            }
         }
      }

      // Cholesky decomposition of the full covariance, including states not in the warm start file
      Rmatrix sqrtPFull(esm.GetStateSize(), esm.GetStateSize());
      cf.Factor(covInAll, sqrtPFull);
      sqrtPFull = sqrtPFull.Transpose();

      for (UnsignedInt ii = stateAlignedIndex + 1; ii < esm.GetStateSize(); ii++)
         for (UnsignedInt jj = ii; jj < esm.GetStateSize(); jj++)
            sqrtP(ii, jj) = sqrtPFull(ii, jj);
   }

   (*(cov->GetCovariance())) = sqrtP * sqrtP.Transpose();

   // Calculate conversion derivative matrix [dX/dS] from Cartesian to Solve-for state
   Rmatrix dX_dS = esm.CartToSolveForStateConversionDerivativeMatrix();
   bool notCartesian = false;

   for (UnsignedInt ii = 0; ii < esm.GetStateSize(); ii++)
   {
      for (UnsignedInt jj = 0; jj < esm.GetStateSize(); jj++)
      {
         Real identityValue = (ii == jj) ? 1.0 : 0.0;

         if (dX_dS(ii, jj) != identityValue)
         {
            notCartesian = true;
            break;
         }
      }

      if (notCartesian)
         break;
   }

   if (notCartesian)
   {
      // Convert covariance back to Keplerian
      Rmatrix dS_dX = dX_dS.Inverse();
      Rmatrix covMat = *(cov->GetCovariance());
      covMat = dS_dX * covMat * dS_dX.Transpose();

      for (UnsignedInt ii = 0; ii < esm.GetStateSize(); ii++)
         for (UnsignedInt jj = 0; jj < esm.GetStateSize(); jj++)
            (*cov)(ii, jj) = covMat(ii, jj);

      cf.Factor(covMat, sqrtP);
      sqrtP = sqrtP.Transpose();
   }
}


//------------------------------------------------------------------------------
// void UpdateCovarianceNominalValues(RealArray prevEpsilonConversions)
//----------------------------------------------------------------------
/**
 * Updates the nominal values of solve-fors which use epsilon values
 *
 * @param prevEpsilonConversions The array the previous epsilon conversion values.
 */
 //------------------------------------------------------------------------------
void UnscentedKalmanFilter::UpdateCovarianceNominalValues(RealArray prevEpsilonConversions)
{
   Covariance * cov = esm.GetCovariance();
   bool covChanged = false;

   if (!sqrtP.IsSized())
   {
      sqrtP.SetSize(esm.GetStateSize(), esm.GetStateSize());
      cf.Factor(*(cov->GetCovariance()), sqrtP);
      sqrtP = sqrtP.Transpose();
   }

   for (UnsignedInt ii = 0; ii < esm.GetStateSize(); ii++)
   {
      Real iiConv = GetEpsilonConversion(ii);

      if (iiConv != prevEpsilonConversions[ii])
      {
         for (UnsignedInt jj = 0; jj <= ii; jj++)
         {
            //Real jjConv = GetEpsilonConversion(jj);

            Real value = sqrtP(ii, jj);
            value *= prevEpsilonConversions[ii] / iiConv;
            sqrtP(ii, jj) = value;
            covChanged = true;
         }
      }
   }

   if (covChanged)
      (*(cov->GetCovariance())) = sqrtP * sqrtP.Transpose();
}


//------------------------------------------------------------------------------
// void CovarianceEpsilonConversion(Rmatrix& cov)
//------------------------------------------------------------------------------
/**
 * This method will convert the terms in the factorized covariance that require a
 * conversion from their epsilon value for reporting.
 *
 * @param cov The covariance matrix to convert
 */
 //------------------------------------------------------------------------------
void UnscentedKalmanFilter::SqrtCovarianceEpsilonConversion(Rmatrix& sqrtCov)
{
   for (UnsignedInt ii = 0; ii < stateSize; ii++)
   {
      Real conv = GetEpsilonConversion(ii);
      for (UnsignedInt jj = 0; jj <= ii; jj++)
      {
         sqrtCov(ii, jj) = sqrtCov(ii, jj) * conv;
      }
   }
}


//------------------------------------------------------------------------------
// void SetCovariance(const Rmatrix& cov)
//------------------------------------------------------------------------------
/**
 * Set the covariance for the filter
 *
 * @param cov The covariance to set
 */
 //------------------------------------------------------------------------------
void UnscentedKalmanFilter::SetCovariance(const Rmatrix& cov)
{
   SeqEstimator::SetCovariance(cov);

   sqrtP.SetSize(cov.GetNumRows(), cov.GetNumColumns()); // Zero out the matrix
   cf.Factor(*(esm.GetCovariance()->GetCovariance()), sqrtP);
   sqrtP = sqrtP.Transpose();
}


//------------------------------------------------------------------------------
// void UnscentedKalmanFilter::UpdateCov()
//------------------------------------------------------------------------------
/**
 * Updates the covariance without updating state or publishing
 */
 //------------------------------------------------------------------------------
void UnscentedKalmanFilter::UpdateCov()
{
   esm.MapObjectsToSTM();
   esm.MapObjectsToVector();
   TimeUpdate();
   (*(stateCovariance->GetCovariance())) = pBar;
   esm.MapCovariancesToObjects();
}



//------------------------------------------------------------------------------
// void GenerateSigmaPoints(const Rvector &mean, const Rmatrix &sqrtCov,
//                          std::vector<Rvector> &sigmaPts) const
//------------------------------------------------------------------------------
/**
 * Builds the 2n+1 sigma points for the unscented transform.
 *
 *    X_0     = mean
 *    X_i     = mean + sqrt(n + lambda) * S_i      (i = 1..n)
 *    X_{i+n} = mean - sqrt(n + lambda) * S_i      (i = 1..n)
 *
 * where S_i is the i-th column of the lower-triangular matrix square root of the
 * covariance (sqrtCov here is the stored sqrtP, with P = sqrtP * sqrtP^T).
 *
 * This is pure linear algebra and is independent of GMAT's propagation; it is
 * the low-risk half of the unscented transform.
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::GenerateSigmaPoints(const Rvector &mean,
                                                const Rmatrix &sqrtCov,
                                                std::vector<Rvector> &sigmaPts) const
{
   const UnsignedInt n = stateSize;
   const Real        gamma = GmatMathUtil::Sqrt((Real)n + lambda);

   sigmaPts.clear();
   sigmaPts.reserve(2 * n + 1);

   // Central point
   sigmaPts.push_back(mean);

   // Plus/minus the scaled columns of the covariance square root
   for (UnsignedInt i = 0; i < n; ++i)
   {
      Rvector col(n);
      for (UnsignedInt r = 0; r < n; ++r)
         col(r) = gamma * sqrtCov(r, i);

      sigmaPts.push_back(mean + col);
   }
   for (UnsignedInt i = 0; i < n; ++i)
   {
      Rvector col(n);
      for (UnsignedInt r = 0; r < n; ++r)
         col(r) = gamma * sqrtCov(r, i);

      sigmaPts.push_back(mean - col);
   }
}


//------------------------------------------------------------------------------
// void BuildSigmaPropagator()
//------------------------------------------------------------------------------
/**
 * Builds sigmaProp: a dedicated, fully-initialized PropSetup used to step sigma
 * points.
 *
 * The base Estimator holds propagators[] straight from the script.  Those are
 * configured but never Initialize()d for standalone stepping -- during a run it
 * is RunEstimator that clones them into its own private list and initializes and
 * steps those (RunEstimator::Execute, p[i]->Initialize()).  The estimator's own
 * propagators[0] therefore has no built propagation-state vector; calling
 * UpdateFromSpaceObject()/Step() on it dereferences unallocated state and crashes.
 *
 * Here we clone propagators[0] and prepare it the same way the command does, but
 * critically wire its PropagationStateManager to the *same* spacecraft objects
 * the estimation state manager (esm) maps to.  That makes the sigma-point
 * round-trip consistent: esm.MapVectorToObjects() writes the seed onto those
 * spacecraft, ode->UpdateFromSpaceObject() reads them into this propagator, we
 * Step(), ode->UpdateSpaceObject() writes back, and esm.GetEstimationState()
 * reads the result.  PrepareInternals() builds the state, wires the ODE model to
 * the PropSetup's own PSM (BuildModelFromMap), and initializes the propagator.
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::BuildSigmaPropagator()
{
   if (sigmaProp != NULL)
   {
      delete sigmaProp;
      sigmaProp = NULL;
   }

   if (propagators.empty() || propagators[0] == NULL)
   {
      MessageInterface::ShowMessage("[UKF Spike] BuildSigmaPropagator: no base "
         "propagator to clone; sigma-point experiment will be skipped.\n");
      return;
   }

   // Spacecraft (SpaceObjects) that esm maps the estimation state onto
   ObjectArray stateObjs;
   esm.GetStateObjects(stateObjs, Gmat::SPACEOBJECT);
   if (stateObjs.empty())
   {
      MessageInterface::ShowMessage("[UKF Spike] BuildSigmaPropagator: esm has no "
         "space objects; sigma-point experiment will be skipped.\n");
      return;
   }

   // Clone the configured PropSetup and re-point it at esm's spacecraft
   sigmaProp = (PropSetup*)propagators[0]->Clone();
   sigmaProp->SetPrecisionTimeFlag(true);
   sigmaProp->SetSolarSystem(solarSystem);

   PropagationStateManager *psm = sigmaProp->GetPropStateManager();
   for (UnsignedInt i = 0; i < stateObjs.size(); ++i)
   {
      psm->SetObject(stateObjs[i]);
      sigmaProp->SetRefObject(stateObjs[i], stateObjs[i]->GetType(),
                              stateObjs[i]->GetName());
   }
   psm->SetProperty("CartesianState");
   // Match the command's prop vector (Cartesian + STM) for a faithful round-trip;
   // the STM is unused by the unscented transform and can be dropped later to cut
   // per-step cost once the mechanism is validated.
   psm->SetProperty("STM");

   try
   {
      sigmaProp->PrepareInternals();
      MessageInterface::ShowMessage("[UKF Spike] BuildSigmaPropagator: initialized "
         "sigma-point propagator '%s' over %d space object(s), state dim %d.\n",
         sigmaProp->GetName().c_str(), (Integer)stateObjs.size(),
         (Integer)psm->GetState()->GetSize());
   }
   catch (BaseException &ex)
   {
      MessageInterface::ShowMessage("[UKF Spike] BuildSigmaPropagator: failed to "
         "initialize sigma-point propagator: %s\n", ex.GetFullMessage().c_str());
      delete sigmaProp;
      sigmaProp = NULL;
   }
}


//------------------------------------------------------------------------------
// void PropagateSigmaPoints(Real dt)
//------------------------------------------------------------------------------
/**
 * SPIKE 1 EXPERIMENT
 *
 * Proves that GMAT's propagation subsystem can re-seed an arbitrary state and
 * re-propagate it from inside the estimator, which is the core capability a true
 * UKF time update requires.  The method:
 *
 *   1. Saves the current reference estimation state (so the live filter is left
 *      unperturbed on exit).
 *   2. Generates the 2n+1 sigma points about that reference using sqrtP.
 *   3. For each sigma point: writes it into the estimation state, maps it onto
 *      the propagated objects, re-seeds the ODE model from those objects, steps
 *      the propagator by dt, then reads the propagated state back.
 *   4. Recombines the propagated points (weighted mean) and logs it against the
 *      reference-trajectory propagation for comparison.
 *   5. Restores the saved reference state and re-seeds the ODE model.
 *
 * Known items to resolve during the build/run loop (the reason this is a spike,
 * not the finished time update):
 *   - Coordinate frame of the estimation state (Cartesian solve-for vs Keplerian
 *     vs epsilon-scaled) must be consistent for propagation; this draft seeds in
 *     the native estimation state and will be refined once the round-trip is
 *     confirmed.
 *   - The propagator state vector carries the STM (variational equations); for
 *     the UKF these are unnecessary and should be disabled once validated.
 *   - State restoration must exactly recover the propagator so the live EKF math
 *     that follows is bit-for-bit unaffected.
 */
//------------------------------------------------------------------------------
void UnscentedKalmanFilter::PropagateSigmaPoints(Real dt)
{
   // Use the dedicated, fully-initialized sigma-point propagator.  The base
   // Estimator's propagators[] are configured but never Initialize()d for
   // standalone stepping (only RunEstimator's private clones are), so stepping
   // them directly would integrate an unbuilt state vector and crash.
   if (sigmaProp == NULL)
   {
      MessageInterface::ShowMessage("[UKF Spike] No initialized sigma-point "
         "propagator available; skipping sigma-point propagation experiment.\n");
      return;
   }

   PropSetup    *prop = sigmaProp;
   ODEModel     *ode  = prop->GetODEModel();
   Propagator   *p    = prop->GetPropagator();

   if (ode == NULL || p == NULL)
   {
      MessageInterface::ShowMessage("[UKF Spike] Propagator has no ODE model "
         "(ephemeris propagator?); skipping sigma-point experiment.\n");
      return;
   }

   const UnsignedInt n = stateSize;

   // 1. Save the reference estimation state so we can restore it afterwards
   GmatState refState = esm.GetEstimationState();
   Rvector   refMean(n);
   for (UnsignedInt i = 0; i < n; ++i)
      refMean(i) = refState[i];

   // 2. Generate sigma points about the reference using the current sqrtP
   std::vector<Rvector> sigmaPts;
   GenerateSigmaPoints(refMean, sqrtP, sigmaPts);

   MessageInterface::ShowMessage(
      "\n[UKF Spike] PropagateSigmaPoints: dt=%.6f s, n=%d, %d sigma points\n",
      dt, n, (Integer)sigmaPts.size());

   const GmatEpoch startEpoch = refState.GetEpoch();
   const GmatEpoch newEpoch   = startEpoch + dt / GmatTimeConstants::SECS_PER_DAY;

   // 3. Propagate each sigma point, re-seeding from the saved reference each time
   std::vector<Rvector> propPts;
   propPts.reserve(sigmaPts.size());

   for (UnsignedInt s = 0; s < sigmaPts.size(); ++s)
   {
      // Seed the estimation state with this sigma point
      GmatState seed = refState;
      for (UnsignedInt i = 0; i < n; ++i)
         seed[i] = sigmaPts[s](i);

      esm.SetEstimationState(seed);
      esm.MapVectorToObjects();

      // Re-seed the integrator from the freshly-set objects and step by dt
      ode->UpdateFromSpaceObject();
      p->Step(dt);
      ode->UpdateSpaceObject(newEpoch);

      // Read the propagated estimation state back out
      GmatState afterState = esm.GetEstimationState();
      Rvector   after(n);
      for (UnsignedInt i = 0; i < n; ++i)
         after(i) = afterState[i];
      propPts.push_back(after);
   }

   // 4. Recombine: weighted mean of the propagated sigma points
   Rvector ukfMean(n);
   for (UnsignedInt i = 0; i < n; ++i)
      ukfMean(i) = 0.0;
   for (UnsignedInt s = 0; s < propPts.size(); ++s)
      for (UnsignedInt i = 0; i < n; ++i)
         ukfMean(i) += wm[s] * propPts[s](i);

   MessageInterface::ShowMessage("[UKF Spike] Recombined sigma-point mean (UT):\n   ");
   for (UnsignedInt i = 0; i < n; ++i)
      MessageInterface::ShowMessage(" % .9e", ukfMean(i));
   MessageInterface::ShowMessage("\n[UKF Spike] Central point propagated (X_0):\n   ");
   for (UnsignedInt i = 0; i < n; ++i)
      MessageInterface::ShowMessage(" % .9e", propPts[0](i));
   MessageInterface::ShowMessage("\n[UKF Spike] UT-mean minus central (nonlinearity signature):\n   ");
   for (UnsignedInt i = 0; i < n; ++i)
      MessageInterface::ShowMessage(" % .3e", ukfMean(i) - propPts[0](i));
   MessageInterface::ShowMessage("\n");

   // 5. Restore the reference state so the live EKF update is unperturbed
   esm.SetEstimationState(refState);
   esm.MapVectorToObjects();
   ode->UpdateFromSpaceObject();
   ode->UpdateSpaceObject(startEpoch);
}
