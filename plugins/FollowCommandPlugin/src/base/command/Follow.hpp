//------------------------------------------------------------------------------
//                                Follow
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
 * Definition of the Follow command class.
 *
 * The Follow command keeps one or more spacecraft at a fixed state expressed
 * in a named coordinate system while time advances, until a stopping condition
 * is met.  No numerical force model or integrator is used; the position and
 * velocity of each frozen spacecraft in the specified coordinate system remain
 * constant throughout the command duration.  Frame rotations (e.g. Earth body-
 * fixed rotation) are handled automatically by CoordinateConverter.
 *
 * Script syntax:
 *   Follow CoordSysName(Sat1 [, Sat2 ...]) {Sat1.ElapsedSecs = 3600.0};
 *
 * The stopping condition block follows the same conventions as Propagate.
 * Only time-based stopping conditions (ElapsedSecs, ElapsedDays) are
 * supported in this initial implementation.
 *
 * Note: Following a numerically propagating leader spacecraft ("Chase" mode)
 * is reserved for a future command.
 */
//------------------------------------------------------------------------------

#ifndef Follow_hpp
#define Follow_hpp

#include "FollowDefs.hpp"
#include "GmatCommand.hpp"
#include "CoordinateConverter.hpp"
#include "CoordinateSystem.hpp"
#include "Spacecraft.hpp"
#include "Rvector6.hpp"
#include <vector>
#include <string>

/**
 * Command that freezes spacecraft in a named coordinate system for an elapsed
 * duration or until a scalar stopping condition is met.
 */
class FOLLOW_API Follow : public GmatCommand
{
public:
   Follow();
   virtual ~Follow();
   Follow(const Follow &f);
   Follow& operator=(const Follow &f);

   // GmatBase interface
   virtual GmatBase*          Clone() const override;
   virtual bool               RenameRefObject(const UnsignedInt type,
                                              const std::string &oldName,
                                              const std::string &newName) override;
   virtual const StringArray& GetRefObjectNameArray(const UnsignedInt type) override;

   // Script generation
   virtual const std::string& GetGeneratingString(
                                  Gmat::WriteMode mode = Gmat::SCRIPTING,
                                  const std::string &prefix = "",
                                  const std::string &useName = "") override;

   // GmatCommand interface
   virtual bool               InterpretAction() override;
   virtual bool               Initialize() override;
   virtual bool               Execute() override;
   virtual GmatCommand*       GetNext() override;
   DEFAULT_TO_NO_CLONES
   virtual void               RunComplete() override;

private:
   // -----------------------------------------------------------------------
   // Configuration (populated by InterpretAction, stored for Clone/script)
   // -----------------------------------------------------------------------

   /// Name of the coordinate system in which spacecraft are frozen
   std::string              coordSysName;

   /// Names of the spacecraft being frozen
   StringArray              satNames;

   /// Raw text of the stopping condition block, e.g. "Sat1.ElapsedSecs = 3600"
   std::string              stopCondText;

   /// Step size in seconds used to advance time each Execute() call
   Real                     stepSize;

   // -----------------------------------------------------------------------
   // Parsed stop condition (filled by InterpretAction)
   // -----------------------------------------------------------------------

   /// Stopping condition parameter string (e.g. "Sat1.ElapsedSecs")
   std::string              stopParam;

   /// Stopping condition operator (currently always "=")
   std::string              stopOp;

   /// Stopping condition goal value in seconds (converted from ElapsedDays)
   Real                     stopGoalSecs;

   /// True when the stop condition is elapsed days rather than elapsed seconds
   bool                     stopIsElapsedDays;

   // -----------------------------------------------------------------------
   // Runtime state (filled by Initialize / first Execute)
   // -----------------------------------------------------------------------

   /// Resolved coordinate system pointer
   CoordinateSystem*        fixedCS;

   /// Resolved spacecraft pointers (same order as satNames)
   std::vector<Spacecraft*> sats;

   /// Frozen states in fixedCS, captured on first Execute call
   std::vector<Rvector6>    frozenStates;

   /// Epoch at the start of the command (A1MJD)
   Real                     startEpoch;

   /// Accumulated elapsed time in seconds
   Real                     elapsedSecs;

   /// True while the command is running (drives GetNext re-entrancy)
   bool                     commandRunning;

   /// True once the stop condition has been met
   bool                     stopCondMet;

   /// Reusable coordinate converter
   CoordinateConverter      coordConverter;

   // -----------------------------------------------------------------------
   // Publisher data (for OFI, GroundTrackPlot, and other subscribers)
   // -----------------------------------------------------------------------

   /// Packed data array published each step: [epoch, X, Y, Z, Vx, Vy, Vz, ...]
   Real*                    pubdata;

   /// Number of elements in pubdata (1 + 6 * numSats)
   Integer                  pubdataSize;

   // -----------------------------------------------------------------------
   // Helpers
   // -----------------------------------------------------------------------

   void  CaptureInitialStates();
   void  UpdateSpacecraftStates(Real newEpoch);
   bool  CheckStopCondition() const;
   void  ParseStopCondition(const std::string &text);
   void  SetupPublisher();
   void  PublishState(Real epoch);
};

#endif /* Follow_hpp */
