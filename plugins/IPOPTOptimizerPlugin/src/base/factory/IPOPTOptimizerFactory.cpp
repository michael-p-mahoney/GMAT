//------------------------------------------------------------------------------
//                          IPOPTOptimizerFactory
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

#include "IPOPTOptimizerFactory.hpp"
#include "IPOPTOptimizer.hpp"

//------------------------------------------------------------------------------
// IPOPTOptimizerFactory()
//------------------------------------------------------------------------------
IPOPTOptimizerFactory::IPOPTOptimizerFactory() :
   Factory(Gmat::SOLVER)
{
   if (creatables.empty())
      creatables.push_back("IPOPT");
}

//------------------------------------------------------------------------------
// ~IPOPTOptimizerFactory()
//------------------------------------------------------------------------------
IPOPTOptimizerFactory::~IPOPTOptimizerFactory()
{
}

//------------------------------------------------------------------------------
// IPOPTOptimizerFactory(const IPOPTOptimizerFactory &copy)
//------------------------------------------------------------------------------
IPOPTOptimizerFactory::IPOPTOptimizerFactory(const IPOPTOptimizerFactory &copy) :
   Factory(copy)
{
   if (creatables.empty())
      creatables.push_back("IPOPT");
}

//------------------------------------------------------------------------------
// IPOPTOptimizerFactory& operator=(const IPOPTOptimizerFactory &fact)
//------------------------------------------------------------------------------
IPOPTOptimizerFactory& IPOPTOptimizerFactory::operator=(
      const IPOPTOptimizerFactory &fact)
{
   if (&fact != this)
      Factory::operator=(fact);
   return *this;
}

//------------------------------------------------------------------------------
// GmatBase* CreateObject(const std::string &ofType, const std::string &withName)
//------------------------------------------------------------------------------
GmatBase* IPOPTOptimizerFactory::CreateObject(const std::string &ofType,
                                               const std::string &withName)
{
   return CreateSolver(ofType, withName);
}

//------------------------------------------------------------------------------
// Solver* CreateSolver(const std::string &ofType, const std::string &withName)
//------------------------------------------------------------------------------
Solver* IPOPTOptimizerFactory::CreateSolver(const std::string &ofType,
                                             const std::string &withName)
{
   if (ofType == "IPOPT")
      return new IPOPTOptimizer(withName);
   return NULL;
}
