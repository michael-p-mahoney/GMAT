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

#ifndef IPOPTOptimizerFactory_hpp
#define IPOPTOptimizerFactory_hpp

#include "ipopt_defs.hpp"
#include "Factory.hpp"
#include "Solver.hpp"

class IPOPT_PLUGIN_API IPOPTOptimizerFactory : public Factory
{
public:
   IPOPTOptimizerFactory();
   virtual ~IPOPTOptimizerFactory();
   IPOPTOptimizerFactory(const IPOPTOptimizerFactory &copy);
   IPOPTOptimizerFactory& operator=(const IPOPTOptimizerFactory &fact);

   virtual GmatBase* CreateObject(const std::string &ofType,
                                  const std::string &withName = "");
   virtual Solver*   CreateSolver(const std::string &ofType,
                                  const std::string &withName = "");
};

#endif // IPOPTOptimizerFactory_hpp
