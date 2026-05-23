//------------------------------------------------------------------------------
//                            GmatPluginFunctions
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
 * Implementation of FollowCommandPlugin library entry points.
 *
 * The three C-linkage functions below are required by the GMAT plugin loader.
 * The Moderator calls GetFactoryCount() and then GetFactoryPointer(i) for
 * each index to register all factories provided by this plugin.
 */
//------------------------------------------------------------------------------

#include "GmatPluginFunctions.hpp"
#include "MessageInterface.hpp"
#include "FollowCommandFactory.hpp"

extern "C"
{

//------------------------------------------------------------------------------
// Integer GetFactoryCount()
//------------------------------------------------------------------------------
/**
 * Returns the number of factories provided by this plugin (always 1)
 */
//------------------------------------------------------------------------------
Integer GetFactoryCount()
{
   return 1;
}


//------------------------------------------------------------------------------
// Factory* GetFactoryPointer(Integer index)
//------------------------------------------------------------------------------
/**
 * Returns a pointer to the factory at the given index.
 *
 * @param index  Factory index; only 0 is valid for this plugin.
 * @return       Pointer to a new FollowCommandFactory, or NULL for bad index.
 */
//------------------------------------------------------------------------------
Factory* GetFactoryPointer(Integer index)
{
   Factory *factory = NULL;

   switch (index)
   {
      case 0:
         factory = new FollowCommandFactory;
         break;
      default:
         break;
   }

   return factory;
}


//------------------------------------------------------------------------------
// void SetMessageReceiver(MessageReceiver *mr)
//------------------------------------------------------------------------------
/**
 * Wires the GMAT message receiver into this shared library.
 *
 * @param mr  The message receiver provided by the GMAT application.
 *
 * @note This function is deprecated and may not be needed in future builds.
 */
//------------------------------------------------------------------------------
void SetMessageReceiver(MessageReceiver *mr)
{
   MessageInterface::SetMessageReceiver(mr);
}

}; // extern "C"
