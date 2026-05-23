//------------------------------------------------------------------------------
//                         FollowCommandFactory
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
 * Implementation of the factory that creates Follow command objects
 */
//------------------------------------------------------------------------------

#include "FollowCommandFactory.hpp"
#include "Follow.hpp"


//------------------------------------------------------------------------------
// FollowCommandFactory()
//------------------------------------------------------------------------------
/**
 * Constructor — registers the "Follow" type string with the Factory base class
 */
//------------------------------------------------------------------------------
FollowCommandFactory::FollowCommandFactory() :
   Factory(Gmat::COMMAND)
{
   if (creatables.empty())
      creatables.push_back("Follow");
}


//------------------------------------------------------------------------------
// ~FollowCommandFactory()
//------------------------------------------------------------------------------
FollowCommandFactory::~FollowCommandFactory()
{
}


//------------------------------------------------------------------------------
// FollowCommandFactory(const FollowCommandFactory &f)
//------------------------------------------------------------------------------
FollowCommandFactory::FollowCommandFactory(const FollowCommandFactory &f) :
   Factory(f)
{
   if (creatables.empty())
      creatables.push_back("Follow");
}


//------------------------------------------------------------------------------
// FollowCommandFactory& operator=(const FollowCommandFactory &f)
//------------------------------------------------------------------------------
FollowCommandFactory& FollowCommandFactory::operator=(
      const FollowCommandFactory &f)
{
   if (this != &f)
   {
      Factory::operator=(f);
      if (creatables.empty())
         creatables.push_back("Follow");
   }
   return *this;
}


//------------------------------------------------------------------------------
// GmatCommand* CreateCommand(const std::string &ofType,
//                            const std::string &withName)
//------------------------------------------------------------------------------
/**
 * Creates a Follow command if ofType matches "Follow"; returns NULL otherwise.
 */
//------------------------------------------------------------------------------
GmatCommand* FollowCommandFactory::CreateCommand(const std::string &ofType,
                                                 const std::string &withName)
{
   if (ofType == "Follow")
      return new Follow();

   return NULL;
}
