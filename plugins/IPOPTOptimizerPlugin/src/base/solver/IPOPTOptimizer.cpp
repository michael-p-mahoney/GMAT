//------------------------------------------------------------------------------
//                              IPOPTOptimizer
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
//
//------------------------------------------------------------------------------

#include "IPOPTOptimizer.hpp"
#include "MessageInterface.hpp"
#include "SolverException.hpp"
#include "GmatGlobal.hpp"

// IPOPT headers
#include "IpIpoptApplication.hpp"
#include "IpTNLP.hpp"
#include "IpSolveStatistics.hpp"

#include <cassert>
#include <cmath>
#include <sstream>
#include <limits>

//#define DEBUG_IPOPT_OPTIMIZER
//#define DEBUG_IPOPT_STATE_MACHINE
//#define DEBUG_IPOPT_EVAL

// ============================================================================
//  GmatTNLP  —  bridges Ipopt::TNLP callbacks to GMAT's state machine.
//  All heavy lifting is done by signalling IPOPTOptimizer via shared data.
// ============================================================================
class GmatTNLP : public Ipopt::TNLP
{
public:
   GmatTNLP(IPOPTOptimizer *opt) : optimizer(opt) {}

   // ------------------------------------------------------------------
   // Ipopt::TNLP interface
   // ------------------------------------------------------------------
   virtual bool get_nlp_info(Ipopt::Index &n, Ipopt::Index &m,
                              Ipopt::Index &nnz_jac_g,
                              Ipopt::Index &nnz_h_lag,
                              IndexStyleEnum &index_style)
   {
      n           = optimizer->variableCount;
      m           = optimizer->eqConstraintCount + optimizer->ineqConstraintCount;
      nnz_jac_g   = n * m;   // Dense Jacobian (all elements non-zero)
      nnz_h_lag   = 0;        // Not used; L-BFGS Hessian approximation
      index_style = C_STYLE;
      return true;
   }

   virtual bool get_bounds_info(Ipopt::Index n,
                                 Ipopt::Number *x_l, Ipopt::Number *x_u,
                                 Ipopt::Index m,
                                 Ipopt::Number *g_l, Ipopt::Number *g_u)
   {
      const double INF = 2.0e19;

      for (int i = 0; i < n; ++i)
      {
         x_l[i] = optimizer->variableMinimum.at(i);
         x_u[i] = optimizer->variableMaximum.at(i);
         // IPOPT uses large-but-finite values for unbounded variables
         if (x_l[i] <= -9.9e299) x_l[i] = -INF;
         if (x_u[i] >=  9.9e299) x_u[i] =  INF;
      }

      int nEq   = optimizer->eqConstraintCount;
      int nIneq = optimizer->ineqConstraintCount;

      for (int j = 0; j < nEq; ++j)
      {
         g_l[j] = optimizer->eqConstraintDesiredValues[j];
         g_u[j] = optimizer->eqConstraintDesiredValues[j];
      }

      for (int j = 0; j < nIneq; ++j)
      {
         int idx = nEq + j;
         int op  = optimizer->ineqConstraintOp[j];
         double desired = optimizer->ineqConstraintDesiredValues[j];

         if (op == 1)         // constraint >= desired  (achieved >= desired)
         {
            g_l[idx] = desired;
            g_u[idx] = INF;
         }
         else if (op == -1)   // constraint <= desired
         {
            g_l[idx] = -INF;
            g_u[idx] = desired;
         }
         else                 // equality
         {
            g_l[idx] = desired;
            g_u[idx] = desired;
         }
      }
      return true;
   }

   virtual bool get_starting_point(Ipopt::Index n, bool init_x,
                                    Ipopt::Number *x,
                                    bool init_z,
                                    Ipopt::Number * /*z_L*/,
                                    Ipopt::Number * /*z_U*/,
                                    Ipopt::Index /*m*/,
                                    bool init_lambda,
                                    Ipopt::Number * /*lambda*/)
   {
      assert(init_x);
      assert(!init_z);
      assert(!init_lambda);

      for (int i = 0; i < n; ++i)
         x[i] = optimizer->variable.at(i);

      return true;
   }

   // ------------------------------------------------------------------
   // Evaluation entry point — triggers GMAT mission sequence execution.
   //
   // Strategy: when IPOPT gives us a new x, we request a FUNC_AND_GRAD
   // evaluation from the main thread (nominal + perturbations), then
   // return the pre-computed values for all subsequent callbacks at
   // the same x.
   // ------------------------------------------------------------------
   virtual bool eval_f(Ipopt::Index n, const Ipopt::Number *x,
                        bool /*new_x*/, Ipopt::Number &obj_value)
   {
      if (!requestEvalIfNeeded(n, x))
         return false;
      obj_value = cachedCost;
      return true;
   }

   virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number *x,
                             bool /*new_x*/, Ipopt::Number *grad_f)
   {
      if (!requestEvalIfNeeded(n, x))
         return false;
      for (int i = 0; i < n; ++i)
         grad_f[i] = cachedGradient.at(i);
      return true;
   }

   virtual bool eval_g(Ipopt::Index n, const Ipopt::Number *x,
                        bool /*new_x*/, Ipopt::Index m,
                        Ipopt::Number *g)
   {
      if (!requestEvalIfNeeded(n, x))
         return false;
      for (int j = 0; j < m; ++j)
         g[j] = cachedConstraints.at(j);
      return true;
   }

   virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number *x,
                             bool /*new_x*/,
                             Ipopt::Index m,
                             Ipopt::Index /*nele_jac*/,
                             Ipopt::Index *iRow, Ipopt::Index *jCol,
                             Ipopt::Number *values)
   {
      if (values == NULL)
      {
         // Sparsity structure: dense, row-major
         int idx = 0;
         for (int r = 0; r < m; ++r)
            for (int c = 0; c < n; ++c)
            {
               iRow[idx] = r;
               jCol[idx] = c;
               ++idx;
            }
         return true;
      }

      if (!requestEvalIfNeeded(n, x))
         return false;

      // cachedJacobian is row-major [nCon x nVar]
      int nnz = n * m;
      for (int i = 0; i < nnz; ++i)
         values[i] = cachedJacobian.at(i);

      return true;
   }

   // L-BFGS is used; IPOPT will not call eval_h.
   virtual bool eval_h(Ipopt::Index /*n*/, const Ipopt::Number * /*x*/,
                        bool /*new_x*/, Ipopt::Number /*obj_factor*/,
                        Ipopt::Index /*m*/, const Ipopt::Number * /*lambda*/,
                        bool /*new_lambda*/, Ipopt::Index /*nele_hess*/,
                        Ipopt::Index * /*iRow*/, Ipopt::Index * /*jCol*/,
                        Ipopt::Number * /*values*/)
   {
      return false; // Not needed with L-BFGS
   }

   virtual void finalize_solution(Ipopt::SolverReturn status,
                                   Ipopt::Index n,
                                   const Ipopt::Number *x,
                                   const Ipopt::Number * /*z_L*/,
                                   const Ipopt::Number * /*z_U*/,
                                   Ipopt::Index /*m*/,
                                   const Ipopt::Number * /*g*/,
                                   const Ipopt::Number * /*lambda*/,
                                   Ipopt::Number obj_value,
                                   const Ipopt::IpoptData * /*ip_data*/,
                                   Ipopt::IpoptCalculatedQuantities * /*ip_cq*/)
   {
      // Signal the main GMAT thread that IPOPT is done
      {
         std::unique_lock<std::mutex> lock(optimizer->sharedMtx);
         auto &sd = optimizer->shared;
         sd.requestType = IPOPTOptimizer::SharedEval::DONE;
         sd.converged   = (status == Ipopt::SUCCESS ||
                           status == Ipopt::STOP_AT_ACCEPTABLE_POINT);
         sd.ipoptStatus = static_cast<int>(status);
         sd.cost        = obj_value;
         for (int i = 0; i < n; ++i)
            sd.x[i] = x[i];
      }
      optimizer->cvRequest.notify_one();

#ifdef DEBUG_IPOPT_OPTIMIZER
      MessageInterface::ShowMessage(
         "IPOPTOptimizer: finalize_solution called, status=%d, obj=%.12f\n",
         (int)status, obj_value);
#endif
   }

private:
   IPOPTOptimizer *optimizer;

   // Cache for the most recently evaluated point
   std::vector<double> cachedX;
   double              cachedCost     = 0.0;
   std::vector<double> cachedConstraints;
   std::vector<double> cachedGradient;
   std::vector<double> cachedJacobian; // row-major [nCon x nVar]

   bool xEquals(int n, const Ipopt::Number *x) const
   {
      if ((int)cachedX.size() != n) return false;
      for (int i = 0; i < n; ++i)
         if (cachedX[i] != x[i]) return false;
      return true;
   }

   // Post an eval request to the GMAT main thread and wait for results.
   bool requestEvalIfNeeded(int n, const Ipopt::Number *x)
   {
      if (xEquals(n, x))
         return true; // Cache hit — IPOPT is re-querying the same point

      // Post the request
      {
         std::unique_lock<std::mutex> lock(optimizer->sharedMtx);
         auto &sd = optimizer->shared;
         sd.requestType = IPOPTOptimizer::SharedEval::FUNC_AND_GRAD;
         sd.x.assign(x, x + n);
      }
      optimizer->cvRequest.notify_one();

#ifdef DEBUG_IPOPT_EVAL
      MessageInterface::ShowMessage("GmatTNLP: posted FUNC_AND_GRAD request\n");
#endif

      // Wait for GMAT to deliver results
      {
         std::unique_lock<std::mutex> lock(optimizer->sharedMtx);
         optimizer->cvResult.wait(lock, [this] { return optimizer->resultReady; });
         optimizer->resultReady = false;

         // Copy results from shared data into local cache
         auto &sd    = optimizer->shared;
         cachedX     = sd.x;
         cachedCost  = sd.cost;
         cachedConstraints = sd.constraints;
         cachedGradient    = sd.gradient;
         cachedJacobian    = sd.jacobian;
      }

#ifdef DEBUG_IPOPT_EVAL
      MessageInterface::ShowMessage(
         "GmatTNLP: received result, cost=%.12f\n", cachedCost);
#endif
      return true;
   }
};

// ============================================================================
//  IPOPTOptimizer — static data
// ============================================================================
const std::string IPOPTOptimizer::PARAMETER_TEXT[
      IPOPTParamCount - InternalOptimizerParamCount] =
{
   "FeasibilityTolerance",
   "OptimalityTolerance",
   "MaximumIterations",
   "MaximumFunctionEvals",
   "UseCentralDifferences",
   "OutputFileName",
   "PrintLevel"
};

const Gmat::ParameterType IPOPTOptimizer::PARAMETER_TYPE[
      IPOPTParamCount - InternalOptimizerParamCount] =
{
   Gmat::REAL_TYPE,
   Gmat::REAL_TYPE,
   Gmat::INTEGER_TYPE,
   Gmat::INTEGER_TYPE,
   Gmat::BOOLEAN_TYPE,
   Gmat::STRING_TYPE,
   Gmat::INTEGER_TYPE
};

// ============================================================================
//  Constructor / Destructor
// ============================================================================
IPOPTOptimizer::IPOPTOptimizer(const std::string &name) :
   InternalOptimizer    ("IPOPT", name),
   feasibilityTolerance (1.0e-6),
   optimalityTolerance  (1.0e-6),
   maximumIterations    (1000),
   maximumFunctionEvals (10000),
   useCentralDifferences(false),
   outputFileName       ("IPOPT.out"),
   printLevel           (5),
   ipoptThreadStarted   (false),
   pertNumber           (-1),
   lastUnperturbedValue (0.0),
   currentPertState     (0),
   evalPhase            (PHASE_NOMINAL)
{
   objectTypeNames.push_back("IPOPT");
   parameterCount = IPOPTParamCount;
   tolerance      = 1.0e-6;
   maxIterations  = 1000;
   isInitialized  = false;
}

IPOPTOptimizer::~IPOPTOptimizer()
{
   // If thread is still running, we need to signal it to stop
   if (ipoptThreadStarted && ipoptThread.joinable())
   {
      // Signal IPOPT to stop by marking as done with failure
      {
         std::lock_guard<std::mutex> lock(sharedMtx);
         shared.requestType = SharedEval::DONE;
         shared.converged   = false;
      }
      cvResult.notify_all();
      ipoptThread.join();
   }
}

IPOPTOptimizer::IPOPTOptimizer(const IPOPTOptimizer &copy) :
   InternalOptimizer    (copy),
   feasibilityTolerance (copy.feasibilityTolerance),
   optimalityTolerance  (copy.optimalityTolerance),
   maximumIterations    (copy.maximumIterations),
   maximumFunctionEvals (copy.maximumFunctionEvals),
   useCentralDifferences(copy.useCentralDifferences),
   outputFileName       (copy.outputFileName),
   printLevel           (copy.printLevel),
   ipoptThreadStarted   (false),
   pertNumber           (-1),
   lastUnperturbedValue (0.0),
   currentPertState     (0),
   evalPhase            (PHASE_NOMINAL)
{
   isInitialized = false;
}

IPOPTOptimizer& IPOPTOptimizer::operator=(const IPOPTOptimizer &copy)
{
   if (&copy != this)
   {
      InternalOptimizer::operator=(copy);
      feasibilityTolerance  = copy.feasibilityTolerance;
      optimalityTolerance   = copy.optimalityTolerance;
      maximumIterations     = copy.maximumIterations;
      maximumFunctionEvals  = copy.maximumFunctionEvals;
      useCentralDifferences = copy.useCentralDifferences;
      outputFileName        = copy.outputFileName;
      printLevel            = copy.printLevel;
      isInitialized = false;
   }
   return *this;
}

void IPOPTOptimizer::Copy(const GmatBase *obj)
{
   if (obj->IsOfType("IPOPT"))
      (*this) = *((IPOPTOptimizer*)obj);
   else
      throw SolverException("Cannot copy non-IPOPT object to IPOPTOptimizer " +
                            instanceName);
}

GmatBase* IPOPTOptimizer::Clone() const
{
   return new IPOPTOptimizer(*this);
}

// ============================================================================
//  Parameter interface
// ============================================================================
bool IPOPTOptimizer::IsParameterReadOnly(const Integer id) const
{
   if (id == OPTIMIZER_TOLERANCE) return true;
   return InternalOptimizer::IsParameterReadOnly(id);
}

bool IPOPTOptimizer::IsParameterReadOnly(const std::string &label) const
{
   return IsParameterReadOnly(GetParameterID(label));
}

Gmat::ParameterType IPOPTOptimizer::GetParameterType(const Integer id) const
{
   if (id >= InternalOptimizerParamCount && id < IPOPTParamCount)
      return PARAMETER_TYPE[id - InternalOptimizerParamCount];
   return InternalOptimizer::GetParameterType(id);
}

std::string IPOPTOptimizer::GetParameterTypeString(const Integer id) const
{
   return InternalOptimizer::PARAM_TYPE_STRING[GetParameterType(id)];
}

std::string IPOPTOptimizer::GetParameterText(const Integer id) const
{
   if (id >= InternalOptimizerParamCount && id < IPOPTParamCount)
      return PARAMETER_TEXT[id - InternalOptimizerParamCount];
   return InternalOptimizer::GetParameterText(id);
}

Integer IPOPTOptimizer::GetParameterID(const std::string &str) const
{
   for (Integer i = InternalOptimizerParamCount; i < IPOPTParamCount; ++i)
      if (str == PARAMETER_TEXT[i - InternalOptimizerParamCount])
         return i;
   return InternalOptimizer::GetParameterID(str);
}

Real IPOPTOptimizer::GetRealParameter(const Integer id) const
{
   switch (id)
   {
   case FEASIBILITY_TOLERANCE: return feasibilityTolerance;
   case OPTIMALITY_TOLERANCE:  return optimalityTolerance;
   default: break;
   }
   return InternalOptimizer::GetRealParameter(id);
}

Real IPOPTOptimizer::SetRealParameter(const Integer id, const Real value)
{
   switch (id)
   {
   case FEASIBILITY_TOLERANCE: feasibilityTolerance = value; return value;
   case OPTIMALITY_TOLERANCE:  optimalityTolerance  = value; return value;
   default: break;
   }
   return InternalOptimizer::SetRealParameter(id, value);
}

Real IPOPTOptimizer::GetRealParameter(const std::string &label) const
{
   return GetRealParameter(GetParameterID(label));
}

Real IPOPTOptimizer::SetRealParameter(const std::string &label, const Real value)
{
   return SetRealParameter(GetParameterID(label), value);
}

Integer IPOPTOptimizer::GetIntegerParameter(const Integer id) const
{
   switch (id)
   {
   case MAXIMUM_ITERATIONS:    return maximumIterations;
   case MAXIMUM_FUNCTION_EVALS: return maximumFunctionEvals;
   case PRINT_LEVEL:           return printLevel;
   default: break;
   }
   return InternalOptimizer::GetIntegerParameter(id);
}

Integer IPOPTOptimizer::SetIntegerParameter(const Integer id, const Integer value)
{
   switch (id)
   {
   case MAXIMUM_ITERATIONS:     maximumIterations    = value; return value;
   case MAXIMUM_FUNCTION_EVALS: maximumFunctionEvals = value; return value;
   case PRINT_LEVEL:            printLevel           = value; return value;
   default: break;
   }
   return InternalOptimizer::SetIntegerParameter(id, value);
}

bool IPOPTOptimizer::GetBooleanParameter(const Integer id) const
{
   if (id == USE_CENTRAL_DIFFERENCES) return useCentralDifferences;
   return InternalOptimizer::GetBooleanParameter(id);
}

bool IPOPTOptimizer::SetBooleanParameter(const Integer id, const bool value)
{
   if (id == USE_CENTRAL_DIFFERENCES) { useCentralDifferences = value; return value; }
   return InternalOptimizer::SetBooleanParameter(id, value);
}

bool IPOPTOptimizer::GetBooleanParameter(const std::string &label) const
{
   return GetBooleanParameter(GetParameterID(label));
}

bool IPOPTOptimizer::SetBooleanParameter(const std::string &label, const bool value)
{
   return SetBooleanParameter(GetParameterID(label), value);
}

std::string IPOPTOptimizer::GetStringParameter(const Integer id) const
{
   if (id == OUTPUT_FILE_NAME) return outputFileName;
   return InternalOptimizer::GetStringParameter(id);
}

bool IPOPTOptimizer::SetStringParameter(const Integer id, const std::string &value)
{
   if (id == OUTPUT_FILE_NAME) { outputFileName = value; return true; }
   return InternalOptimizer::SetStringParameter(id, value);
}

std::string IPOPTOptimizer::GetStringParameter(const std::string &label) const
{
   return GetStringParameter(GetParameterID(label));
}

bool IPOPTOptimizer::SetStringParameter(const std::string &label,
                                         const std::string &value)
{
   return SetStringParameter(GetParameterID(label), value);
}

const StringArray& IPOPTOptimizer::GetRefObjectNameArray(
      const Gmat::ObjectType type)
{
   dummyArray.clear();
   return dummyArray;
}

// ============================================================================
//  Solver state machine setup
// ============================================================================
Integer IPOPTOptimizer::SetSolverResults(Real *data, const std::string &name,
                                          const std::string &type)
{
   return InternalOptimizer::SetSolverResults(data, name, type);
}

void IPOPTOptimizer::SetResultValue(Integer id, Real value,
                                     const std::string &resultType)
{
#ifdef DEBUG_IPOPT_OPTIMIZER
   MessageInterface::ShowMessage(
      "IPOPTOptimizer::SetResultValue id=%d type=%s value=%.12f phase=%d\n",
      id, resultType.c_str(), value, (int)evalPhase);
#endif

   bool plusEffect = true;
   if (useCentralDifferences && currentPertState == -1)
      plusEffect = false;

   if (resultType == "Objective")
   {
      cost = value;
      gradientCalculator.Achieved(pertNumber, 0,
         (pertNumber < 0 ? 0.0 : (pertNumber < (int)perturbation.size() ?
                                   perturbation[pertNumber] : 0.0)),
         value, plusEffect);
   }
   else
   {
      Integer idToUse;
      if (resultType == "EqConstraint")
      {
         idToUse = id - 1000;
         if (pertNumber < 0)
            eqConstraintValues[idToUse] = value;
      }
      else // IneqConstraint
      {
         idToUse = (id - 2000) + eqConstraintCount;
         if (pertNumber < 0)
            ineqConstraintValues[id - 2000] = value;
      }

      jacobianCalculator.Achieved(pertNumber, idToUse,
         (pertNumber < 0 ? 0.0 : (pertNumber < (int)perturbation.size() ?
                                   perturbation[pertNumber] : 0.0)),
         value, plusEffect);
   }
}

bool IPOPTOptimizer::Initialize()
{
#ifdef DEBUG_IPOPT_OPTIMIZER
   MessageInterface::ShowMessage(
      "IPOPTOptimizer::Initialize() entered for %s; "
      "%d vars, %d constraints\n",
      instanceName.c_str(), registeredVariableCount, registeredComponentCount);
#endif

   bool retval = InternalOptimizer::Initialize();

   if (retval)
   {
      DerivativeModel::derivativeMode diffMode = useCentralDifferences ?
         DerivativeModel::CENTRAL_DIFFERENCE : DerivativeModel::FORWARD_DIFFERENCE;

      gradientCalculator.SetDifferenceMode(diffMode);
      jacobianCalculator.SetDifferenceMode(diffMode);

      retval = gradientCalculator.Initialize(registeredVariableCount);
   }

   if (retval && registeredComponentCount > 0)
      retval = jacobianCalculator.Initialize(registeredVariableCount,
                                             registeredComponentCount);

   // Pre-size accumulators
   gradientAccum.assign(registeredVariableCount, 0.0);
   jacobianAccum.assign(registeredVariableCount * registeredComponentCount, 0.0);

   // Pre-size shared state vectors
   {
      std::lock_guard<std::mutex> lock(sharedMtx);
      shared.x.resize(registeredVariableCount, 0.0);
      shared.constraints.resize(eqConstraintCount + ineqConstraintCount, 0.0);
      shared.gradient.resize(registeredVariableCount, 0.0);
      shared.jacobian.resize(registeredVariableCount *
                             (eqConstraintCount + ineqConstraintCount), 0.0);
      shared.requestType = SharedEval::NONE;
      resultReady        = false;
   }

   ipoptThreadStarted = false;
   pertNumber         = -1;
   currentPertState   = 0;
   evalPhase          = PHASE_NOMINAL;
   isInitialized      = true;

#ifdef DEBUG_IPOPT_OPTIMIZER
   MessageInterface::ShowMessage("IPOPTOptimizer::Initialize() returning %d\n",
                                 retval);
#endif
   return retval;
}

// ============================================================================
//  AdvanceState — main state machine
// ============================================================================
Solver::SolverState IPOPTOptimizer::AdvanceState()
{
#ifdef DEBUG_IPOPT_STATE_MACHINE
   MessageInterface::ShowMessage(
      "IPOPTOptimizer::AdvanceState() state=%d phase=%d pert=%d\n",
      (int)currentState, (int)evalPhase, pertNumber);
#endif

   switch (currentState)
   {
   case INITIALIZING:
      iterationsTaken = 0;
      WriteToTextFile();
      ReportProgress();
      if (!isInitialized)
         Initialize();
      CompleteInitialization();
      break;

   case NOMINAL:
      // GMAT has finished running the mission sequence.
      // Results have been accumulated via SetResultValue().
      // Decide what to do next based on the current eval phase.
      WriteToTextFile();
      ReportProgress();
      if (evalPhase == PHASE_NOMINAL)
      {
         // Nominal done — begin perturbation series
         pertNumber = -1;
         RunPerturbation();
      }
      else if (evalPhase == PHASE_PERTURBING)
      {
         // A perturbation run just finished — advance to next
         RunPerturbation();
      }
      break;

   case PERTURBING:
      // This state is set by RunPerturbation when it needs another run.
      // The Optimize command loop will execute the mission sequence and
      // call back into AdvanceState with NOMINAL after the run.
      RunPerturbation();
      break;

   case CALCULATING:
      CalculateParameters();
      WriteToTextFile();
      ReportProgress();
      break;

   case CHECKINGRUN:
      CheckCompletion();
      break;

   case FINISHED:
      WriteToTextFile();
      ReportProgress();
      RunComplete();
      break;

   default:
      break;
   }

#ifdef DEBUG_IPOPT_STATE_MACHINE
   MessageInterface::ShowMessage(
      "IPOPTOptimizer::AdvanceState() returning state=%d\n",
      (int)currentState);
#endif
   return currentState;
}

bool IPOPTOptimizer::Optimize()
{
   return true;
}

// ============================================================================
//  State machine helpers
// ============================================================================
void IPOPTOptimizer::RunNominal()
{
   // Called when we need to run a nominal (unperturbed) mission sequence.
   // Set variables to the values in variable[], then trigger execution.
   pertNumber   = -1;
   currentState = NOMINAL;
   evalPhase    = PHASE_NOMINAL;
   status       = Gmat::RUN;
   GmatGlobal::Instance()->SetSolverStatus(instanceName, status);
}

void IPOPTOptimizer::RunPerturbation()
{
   if (!useCentralDifferences)
   {
      // Restore previous perturbation if any
      if (pertNumber >= 0 && pertNumber < variableCount)
         variable.at(pertNumber) = lastUnperturbedValue;

      ++pertNumber;

      if (pertNumber == variableCount)
      {
         // All perturbations done — compute gradient/Jacobian
         currentState = CALCULATING;
         evalPhase    = PHASE_DONE;
         pertNumber   = -1;
         return;
      }

      lastUnperturbedValue = variable.at(pertNumber);
      variable.at(pertNumber) += perturbation.at(pertNumber);

      currentState = NOMINAL;
      evalPhase    = PHASE_PERTURBING;
      status       = Gmat::RUN;
      GmatGlobal::Instance()->SetSolverStatus(instanceName, status);
   }
   else
   {
      // Central differences: each variable requires + and - perturbations
      if (pertNumber >= 0)
         variable.at(pertNumber) = lastUnperturbedValue;

      if (currentPertState == -1)
         currentPertState = 0;

      if (currentPertState == 1)
      {
         currentPertState = -1;
      }
      else
      {
         ++pertNumber;
         if (pertNumber == variableCount)
         {
            currentState = CALCULATING;
            evalPhase    = PHASE_DONE;
            pertNumber   = -1;
            currentPertState = 0;
            return;
         }
         currentPertState = 1;
      }

      lastUnperturbedValue  = variable.at(pertNumber);
      variable.at(pertNumber) += currentPertState * perturbation.at(pertNumber);

      currentState = NOMINAL;
      evalPhase    = PHASE_PERTURBING;
      status       = Gmat::RUN;
      GmatGlobal::Instance()->SetSolverStatus(instanceName, status);
   }
}

void IPOPTOptimizer::CalculateParameters()
{
   // Compute gradient and Jacobian from accumulated finite-difference data
   if (objectiveDefined)
      gradientCalculator.Calculate(gradientAccum);
   if (eqConstraintCount + ineqConstraintCount > 0)
      jacobianCalculator.Calculate(jacobianAccum);

#ifdef DEBUG_IPOPT_OPTIMIZER
   MessageInterface::ShowMessage("IPOPTOptimizer: gradient=[");
   for (size_t i = 0; i < gradientAccum.size(); ++i)
      MessageInterface::ShowMessage("%.8f ", gradientAccum[i]);
   MessageInterface::ShowMessage("]\n");
#endif

   // Package results into shared data and signal IPOPT thread
   {
      std::lock_guard<std::mutex> lock(sharedMtx);
      auto &sd         = shared;
      sd.cost          = cost;
      sd.constraints   = std::vector<double>();

      int nEq   = eqConstraintCount;
      int nIneq = ineqConstraintCount;
      for (int j = 0; j < nEq; ++j)
         sd.constraints.push_back(eqConstraintValues[j]);
      for (int j = 0; j < nIneq; ++j)
         sd.constraints.push_back(ineqConstraintValues[j]);

      sd.gradient = std::vector<double>(gradientAccum.begin(),
                                        gradientAccum.end());

      // Jacobian: row-major [nCon x nVar]
      int nCon = nEq + nIneq;
      int nVar = variableCount;
      sd.jacobian.resize(nCon * nVar, 0.0);
      for (int r = 0; r < nCon; ++r)
         for (int c = 0; c < nVar; ++c)
            sd.jacobian[r * nVar + c] = jacobianAccum[c + nVar * r];

      resultReady = true;
   }
   cvResult.notify_one();

   // Transition to CHECKINGRUN to wait for IPOPT's next request
   currentState = CHECKINGRUN;
}

void IPOPTOptimizer::CheckCompletion()
{
   if (!ipoptThreadStarted)
   {
      // First entry: launch the IPOPT background thread
      ipoptThread        = std::thread(RunIPOPT, this);
      ipoptThreadStarted = true;

#ifdef DEBUG_IPOPT_OPTIMIZER
      MessageInterface::ShowMessage(
         "IPOPTOptimizer: IPOPT thread launched\n");
#endif
   }

   // Wait for IPOPT thread to post a request (or signal completion)
   std::unique_lock<std::mutex> lock(sharedMtx);
   cvRequest.wait(lock, [this]
   {
      return shared.requestType != SharedEval::NONE;
   });

   if (shared.requestType == SharedEval::DONE)
   {
      // IPOPT has finished
      bool didConverge = shared.converged;
      lock.unlock();

      if (ipoptThread.joinable())
         ipoptThread.join();

      converged = didConverge;
      status    = didConverge ? Gmat::CONVERGED : Gmat::FAILED;
      GmatGlobal::Instance()->SetSolverStatus(instanceName, status);
      currentState = FINISHED;
      ++iterationsTaken;
      return;
   }

   // IPOPT needs a FUNC_AND_GRAD evaluation at shared.x
   // Set variable[] and kick off a new nominal+perturbation cycle
   for (int i = 0; i < variableCount; ++i)
      variable.at(i) = shared.x[i];
   shared.requestType = SharedEval::NONE;
   resultReady        = false;
   lock.unlock();

   // Reset finite-difference calculators for new x
   gradientCalculator.Initialize(registeredVariableCount);
   if (registeredComponentCount > 0)
      jacobianCalculator.Initialize(registeredVariableCount,
                                    registeredComponentCount);
   gradientAccum.assign(registeredVariableCount, 0.0);
   jacobianAccum.assign(registeredVariableCount * registeredComponentCount, 0.0);

   // Start the nominal run
   pertNumber   = -1;
   evalPhase    = PHASE_NOMINAL;
   currentState = NOMINAL;
   status       = Gmat::RUN;
   GmatGlobal::Instance()->SetSolverStatus(instanceName, status);
   ++iterationsTaken;
}

void IPOPTOptimizer::RunComplete()
{
   hasFired = true;
}

void IPOPTOptimizer::FreeArrays()
{
}

// ============================================================================
//  IPOPT background thread
// ============================================================================
void IPOPTOptimizer::RunIPOPT(IPOPTOptimizer *optimizer)
{
   using namespace Ipopt;

   SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

   // Configure IPOPT
   app->Options()->SetNumericValue("tol",              optimizer->optimalityTolerance);
   app->Options()->SetNumericValue("constr_viol_tol",  optimizer->feasibilityTolerance);
   app->Options()->SetIntegerValue("max_iter",         optimizer->maximumIterations);
   app->Options()->SetIntegerValue("print_level",      optimizer->printLevel);
   app->Options()->SetStringValue ("hessian_approximation", "limited-memory");
   app->Options()->SetStringValue ("derivative_test",  "none");
   app->Options()->SetStringValue ("output_file",      optimizer->outputFileName);

   // Suppress console output if print level is low
   if (optimizer->printLevel <= 0)
      app->Options()->SetIntegerValue("print_level", 0);

   app->Initialize();

#ifdef DEBUG_IPOPT_OPTIMIZER
   MessageInterface::ShowMessage("IPOPTOptimizer: IpoptApplication initialized\n");
#endif

   SmartPtr<TNLP> nlp = new GmatTNLP(optimizer);
   app->OptimizeTNLP(nlp);

   // finalize_solution in GmatTNLP signals the main thread, so we're done.
#ifdef DEBUG_IPOPT_OPTIMIZER
   MessageInterface::ShowMessage("IPOPTOptimizer: IPOPT thread returning\n");
#endif
}

// ============================================================================
//  Reporting
// ============================================================================
void IPOPTOptimizer::WriteToTextFile(SolverState stateToUse)
{
   if (!showProgress) return;
   if (!textFile.is_open())
      OpenSolverTextFile();
   if (!isInitialized) return;

   SolverState s = (stateToUse == UNDEFINED_STATE) ? currentState : stateToUse;
   switch (s)
   {
   case INITIALIZING:
      textFile << "IPOPT Optimizer " << instanceName << " initializing.\n";
      break;
   case NOMINAL:
      if (pertNumber < 0)
         textFile << "Iteration " << iterationsTaken
                  << ": running nominal mission sequence.\n";
      break;
   case FINISHED:
      textFile << "Optimization complete. Converged: "
               << (converged ? "Yes" : "No") << "\n";
      break;
   default:
      break;
   }
   textFile.flush();
}

std::string IPOPTOptimizer::GetProgressString()
{
   std::stringstream ss;
   ss << "IPOPT[" << instanceName << "]: iter=" << iterationsTaken
      << " cost=" << cost;
   return ss.str();
}
