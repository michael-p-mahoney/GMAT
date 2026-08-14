//------------------------------------------------------------------------------
//                          GmatPluginFunctions
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

#ifndef GmatPluginFunctions_hpp
#define GmatPluginFunctions_hpp

#include "ipopt_defs.hpp"
#include "Factory.hpp"

class MessageReceiver;

extern "C"
{
   IPOPT_PLUGIN_API Integer    GetFactoryCount();
   IPOPT_PLUGIN_API Factory*   GetFactoryPointer(Integer index);
   IPOPT_PLUGIN_API void       SetMessageReceiver(MessageReceiver* mr);
};

#endif // GmatPluginFunctions_hpp
