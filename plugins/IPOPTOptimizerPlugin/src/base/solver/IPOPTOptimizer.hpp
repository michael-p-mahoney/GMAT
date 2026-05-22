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

#ifndef IPOPTOptimizer_hpp
#define IPOPTOptimizer_hpp

#include "ipopt_defs.hpp"
#include "InternalOptimizer.hpp"
#include "Gradient.hpp"
#include "Jacobian.hpp"

#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <atomic>
#include <ctime>

/**
 * IPOPT-based NLP optimizer for GMAT scripting.
 *
 * Bridges IPOPT's synchronous callback interface (Ipopt::TNLP) with GMAT's
 * iterative state-machine-based mission sequence execution. IPOPT runs on a
 * background thread; each function/gradient evaluation request triggers GMAT
 * to run the mission sequence on the main thread via condition variables.
 *
 * Usage in a GMAT script:
 *   Create IPOPT myOpt;
 *   myOpt.FeasibilityTolerance = 1e-6;
 *   myOpt.OptimalityTolerance  = 1e-6;
 *   myOpt.MaximumIterations    = 1000;
 *   ...
 *   Optimize myOpt { ... };
 */
class IPOPT_PLUGIN_API IPOPTOptimizer : public InternalOptimizer
{
public:
   IPOPTOptimizer(const std::string &name);
   virtual ~IPOPTOptimizer();
   IPOPTOptimizer(const IPOPTOptimizer &copy);
   IPOPTOptimizer& operator=(const IPOPTOptimizer &copy);

   virtual void         Copy(const GmatBase *obj);
   virtual GmatBase*    Clone() const;

   // Parameter interface
   virtual bool         IsParameterReadOnly(const Integer id) const;
   virtual bool         IsParameterReadOnly(const std::string &label) const;
   virtual Gmat::ParameterType
                        GetParameterType(const Integer id) const;
   virtual std::string  GetParameterTypeString(const Integer id) const;
   virtual std::string  GetParameterText(const Integer id) const;
   virtual Integer      GetParameterID(const std::string &str) const;

   virtual Real         GetRealParameter(const Integer id) const;
   virtual Real         SetRealParameter(const Integer id, const Real value);
   virtual Real         GetRealParameter(const std::string &label) const;
   virtual Real         SetRealParameter(const std::string &label, const Real value);

   virtual Integer      GetIntegerParameter(const Integer id) const;
   virtual Integer      SetIntegerParameter(const Integer id, const Integer value);

   virtual bool         GetBooleanParameter(const Integer id) const;
   virtual bool         SetBooleanParameter(const Integer id, const bool value);
   virtual bool         GetBooleanParameter(const std::string &label) const;
   virtual bool         SetBooleanParameter(const std::string &label, const bool value);

   virtual std::string  GetStringParameter(const Integer id) const;
   virtual bool         SetStringParameter(const Integer id, const std::string &value);
   virtual std::string  GetStringParameter(const std::string &label) const;
   virtual bool         SetStringParameter(const std::string &label, const std::string &value);

   virtual const StringArray& GetRefObjectNameArray(const Gmat::ObjectType type);

   virtual Integer      SetSolverResults(Real *data, const std::string &name,
                                         const std::string &type = "");
   virtual void         SetResultValue(Integer id, Real value,
                                       const std::string &resultType = "");

   virtual bool         Initialize();
   virtual SolverState  AdvanceState();
   virtual bool         Optimize();

   virtual void         WriteToTextFile(SolverState stateToUse = UNDEFINED_STATE);
   virtual std::string  GetProgressString();

   // Parameter IDs
   enum
   {
      FEASIBILITY_TOLERANCE = InternalOptimizerParamCount,
      OPTIMALITY_TOLERANCE,
      MAXIMUM_ITERATIONS,
      MAXIMUM_FUNCTION_EVALS,
      MAXIMUM_RUN_TIME,
      USE_CENTRAL_DIFFERENCES,
      OUTPUT_FILE_NAME,
      PRINT_LEVEL,
      IPOPTParamCount
   };

   static const std::string PARAMETER_TEXT[IPOPTParamCount - InternalOptimizerParamCount];
   static const Gmat::ParameterType PARAMETER_TYPE[IPOPTParamCount - InternalOptimizerParamCount];

   DEFAULT_TO_NO_CLONES
   DEFAULT_TO_NO_REFOBJECTS

   // -----------------------------------------------------------------------
   // Shared state between main thread and IPOPT background thread.
   // All fields protected by sharedMtx.
   // -----------------------------------------------------------------------
   struct SharedEval
   {
      /// Type of evaluation the IPOPT thread is requesting
      enum RequestType { NONE, FUNC_AND_GRAD, DONE };
      RequestType requestType = NONE;

      /// Decision vector IPOPT wants evaluated
      std::vector<double> x;

      /// Results GMAT provides back
      double              cost         = 0.0;
      std::vector<double> constraints;   ///< all constraint values (eq then ineq)
      std::vector<double> gradient;      ///< objective gradient
      std::vector<double> jacobian;      ///< constraint Jacobian, row-major [nCon x nVar]

      /// Convergence info (set when requestType == DONE)
      bool                converged    = false;
      int                 ipoptStatus  = 0;
   };

   std::mutex              sharedMtx;
   std::condition_variable cvRequest;   ///< IPOPT → GMAT: eval ready to be picked up
   std::condition_variable cvResult;    ///< GMAT → IPOPT: result ready

   SharedEval              shared;      ///< protected by sharedMtx
   bool                    resultReady  = false;

   /// Cooperative cancel flag.  Set by the destructor (or any future GUI
   /// abort path) to ask the IPOPT thread to unwind.  intermediate_callback
   /// and requestEvalIfNeeded both return false to IPOPT when this is true,
   /// which makes OptimizeTNLP exit through finalize_solution with
   /// USER_REQUESTED_STOP / NonIpopt_Exception_Thrown.
   std::atomic<bool>       abortRequested{false};

   /// IPOPT solve start time (set in RunIPOPT, read in intermediate_callback
   /// for wall-clock budget enforcement).
   time_t                  nlpStartTime{0};

protected:
   // Optimizer parameters
   Real     feasibilityTolerance;
   Real     optimalityTolerance;
   Integer  maximumIterations;
   Integer  maximumFunctionEvals;
   Real     maximumRunTime;       ///< seconds; 0 = no wall-clock limit
   bool     useCentralDifferences;
   std::string outputFileName;
   Integer  printLevel;

   // Internal state
   bool     ipoptThreadStarted;
   std::thread ipoptThread;

   /// Current perturbation index (-1 = nominal)
   Integer  pertNumber;
   Real     lastUnperturbedValue;
   Integer  currentPertState;    ///< used for central differences (+1 / -1)

   /// Finite-difference calculators (reuse GMAT's existing classes)
   Gradient gradientCalculator;
   Jacobian jacobianCalculator;

   /// Accumulated gradient/Jacobian before packaging for IPOPT
   std::vector<Real> gradientAccum;
   std::vector<Real> jacobianAccum;

   /// Phase of the current eval cycle
   enum EvalPhase { PHASE_NOMINAL, PHASE_PERTURBING, PHASE_DONE };
   EvalPhase evalPhase;

   StringArray dummyArray;

   // State machine helpers
   void RunNominal();
   void RunPerturbation();
   void CalculateParameters();
   void CheckCompletion();
   void RunComplete();
   void FreeArrays();

   // IPOPT thread entry point (static to pass as std::thread target)
   static void RunIPOPT(IPOPTOptimizer *optimizer);

   // GmatTNLP accesses protected base class members (Solver/Optimizer)
   // through an IPOPTOptimizer pointer; friendship is required.
   friend class GmatTNLP;
};

#endif // IPOPTOptimizer_hpp
