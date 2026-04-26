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

#include "GmatPluginFunctions.hpp"
#include "IPOPTOptimizerFactory.hpp"
#include "MessageInterface.hpp"

extern "C"
{

//------------------------------------------------------------------------------
// Integer GetFactoryCount()
//------------------------------------------------------------------------------
Integer GetFactoryCount()
{
   return 1;
}

//------------------------------------------------------------------------------
// Factory* GetFactoryPointer(Integer index)
//------------------------------------------------------------------------------
Factory* GetFactoryPointer(Integer index)
{
   if (index == 0)
      return new IPOPTOptimizerFactory();
   return NULL;
}

//------------------------------------------------------------------------------
// void SetMessageReceiver(MessageReceiver* mr)
//------------------------------------------------------------------------------
void SetMessageReceiver(MessageReceiver* mr)
{
   MessageInterface::SetMessageReceiver(mr);
}

} // extern "C"
