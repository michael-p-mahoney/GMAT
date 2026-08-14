//------------------------------------------------------------------------------
//                              ipopt_defs
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

#ifndef ipopt_defs_hpp
#define ipopt_defs_hpp

#include "gmatdefs.hpp"

#ifdef _WIN32
   #ifdef IPOPT_EXPORTS
      #define IPOPT_PLUGIN_API __declspec(dllexport)
   #else
      #define IPOPT_PLUGIN_API __declspec(dllimport)
   #endif
#else
   #define IPOPT_PLUGIN_API
#endif

#endif // ipopt_defs_hpp
