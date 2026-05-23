//------------------------------------------------------------------------------
//                            FollowDefs
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool
//
// Copyright (c) 2002-2026 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
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
// Author: Michael Mahoney, NASA/GSFC
// Created: 2026-04-10
//
/**
 * DLL interface definitions for the Follow command plugin
 */
//------------------------------------------------------------------------------

#ifndef FollowDefs_hpp
#define FollowDefs_hpp

#include "gmatdefs.hpp"

#ifdef _WIN32  // Windows
   #ifdef _MSC_VER  // Microsoft Visual C++

      #define WIN32_LEAN_AND_MEAN
      #include <windows.h>
      #define  _USE_MATH_DEFINES

   #endif

   #ifdef _DYNAMICLINK
      #ifdef FOLLOW_EXPORTS
         #define FOLLOW_API __declspec(dllexport)
      #else
         #define FOLLOW_API __declspec(dllimport)
      #endif

      #ifdef EXP_STL
      #    define DECLSPECIFIER __declspec(dllexport)
      #    define EXPIMP_TEMPLATE
      #else
      #    define DECLSPECIFIER __declspec(dllimport)
      #    define EXPIMP_TEMPLATE extern
      #endif

      #define EXPORT_TEMPLATES
   #endif
#endif //  End of OS nits

#ifndef FOLLOW_API
   #define FOLLOW_API
#endif

#endif /* FollowDefs_hpp */
