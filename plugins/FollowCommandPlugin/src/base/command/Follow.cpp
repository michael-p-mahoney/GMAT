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
 * Implementation of the Follow command.
 *
 * See Follow.hpp for full documentation.
 */
//------------------------------------------------------------------------------

#include "Follow.hpp"
#include "CommandException.hpp"
#include "MessageInterface.hpp"
#include "GmatConstants.hpp"   // GmatTimeConstants::SECS_PER_DAY
#include "StringUtil.hpp"
#include "GmatBase.hpp"

//------------------------------------------------------------------------------
// Static helpers
//------------------------------------------------------------------------------

namespace
{
   /// Default time step in seconds per Execute() call
   const Real DEFAULT_STEP_SIZE = 60.0;
}


//------------------------------------------------------------------------------
// Follow()
//------------------------------------------------------------------------------
/**
 * Default constructor
 */
//------------------------------------------------------------------------------
Follow::Follow() :
   GmatCommand         ("Follow"),
   coordSysName        (""),
   stopCondText        (""),
   stepSize            (DEFAULT_STEP_SIZE),
   stopParam           (""),
   stopOp              ("="),
   stopGoalSecs        (0.0),
   stopIsElapsedDays   (false),
   fixedCS             (NULL),
   startEpoch          (0.0),
   elapsedSecs         (0.0),
   commandRunning      (false),
   stopCondMet         (false),
   pubdata             (NULL),
   pubdataSize         (0)
{
}


//------------------------------------------------------------------------------
// ~Follow()
//------------------------------------------------------------------------------
/**
 * Destructor
 */
//------------------------------------------------------------------------------
Follow::~Follow()
{
   // sats and fixedCS are sandbox-owned; do not delete here
   if (pubdata)
   {
      delete[] pubdata;
      pubdata = NULL;
   }
}


//------------------------------------------------------------------------------
// Follow(const Follow &f)
//------------------------------------------------------------------------------
/**
 * Copy constructor
 */
//------------------------------------------------------------------------------
Follow::Follow(const Follow &f) :
   GmatCommand         (f),
   coordSysName        (f.coordSysName),
   satNames            (f.satNames),
   stopCondText        (f.stopCondText),
   stepSize            (f.stepSize),
   stopParam           (f.stopParam),
   stopOp              (f.stopOp),
   stopGoalSecs        (f.stopGoalSecs),
   stopIsElapsedDays   (f.stopIsElapsedDays),
   fixedCS             (NULL),          // resolved fresh on Initialize
   sats                (),              // resolved fresh on Initialize
   frozenStates        (),              // captured fresh on first Execute
   startEpoch          (0.0),
   elapsedSecs         (0.0),
   commandRunning      (false),
   stopCondMet         (false),
   pubdata             (NULL),
   pubdataSize         (0)
{
}


//------------------------------------------------------------------------------
// Follow& operator=(const Follow &f)
//------------------------------------------------------------------------------
/**
 * Assignment operator
 */
//------------------------------------------------------------------------------
Follow& Follow::operator=(const Follow &f)
{
   if (this != &f)
   {
      GmatCommand::operator=(f);
      coordSysName      = f.coordSysName;
      satNames          = f.satNames;
      stopCondText      = f.stopCondText;
      stepSize          = f.stepSize;
      stopParam         = f.stopParam;
      stopOp            = f.stopOp;
      stopGoalSecs      = f.stopGoalSecs;
      stopIsElapsedDays = f.stopIsElapsedDays;
      fixedCS           = NULL;
      sats.clear();
      frozenStates.clear();
      startEpoch        = 0.0;
      elapsedSecs       = 0.0;
      commandRunning    = false;
      stopCondMet       = false;
      if (pubdata)
      {
         delete[] pubdata;
         pubdata = NULL;
      }
      pubdataSize       = 0;
   }
   return *this;
}


//------------------------------------------------------------------------------
// GmatBase* Clone() const
//------------------------------------------------------------------------------
/**
 * Clone this command
 */
//------------------------------------------------------------------------------
GmatBase* Follow::Clone() const
{
   return new Follow(*this);
}


//------------------------------------------------------------------------------
// bool RenameRefObject(const UnsignedInt type, const std::string &oldName,
//                      const std::string &newName)
//------------------------------------------------------------------------------
/**
 * Renames object references if a sandbox object is renamed
 */
//------------------------------------------------------------------------------
bool Follow::RenameRefObject(const UnsignedInt type,
                             const std::string &oldName,
                             const std::string &newName)
{
   // Rename coordinate system reference
   if (type == Gmat::COORDINATE_SYSTEM)
   {
      if (coordSysName == oldName)
         coordSysName = newName;
   }

   // Rename spacecraft references
   if (type == Gmat::SPACECRAFT)
   {
      for (UnsignedInt i = 0; i < satNames.size(); ++i)
         if (satNames[i] == oldName)
            satNames[i] = newName;

      // Also rename inside the stop condition text
      std::string::size_type pos = stopCondText.find(oldName);
      while (pos != std::string::npos)
      {
         stopCondText.replace(pos, oldName.size(), newName);
         pos = stopCondText.find(oldName, pos + newName.size());
      }

      // Keep stopParam in sync
      pos = stopParam.find(oldName);
      if (pos != std::string::npos)
         stopParam.replace(pos, oldName.size(), newName);
   }

   // Pass to base for any generic processing
   return GmatCommand::RenameRefObject(type, oldName, newName);
}


//------------------------------------------------------------------------------
// const StringArray& GetRefObjectNameArray(const UnsignedInt type)
//------------------------------------------------------------------------------
/**
 * Returns names of referenced objects so the Sandbox can resolve them
 */
//------------------------------------------------------------------------------
const StringArray& Follow::GetRefObjectNameArray(const UnsignedInt type)
{
   refObjectNames.clear();

   if (type == Gmat::COORDINATE_SYSTEM || type == Gmat::UNKNOWN_OBJECT)
      refObjectNames.push_back(coordSysName);

   if (type == Gmat::SPACECRAFT || type == Gmat::UNKNOWN_OBJECT)
      for (const auto &name : satNames)
         refObjectNames.push_back(name);

   return refObjectNames;
}


//------------------------------------------------------------------------------
// const std::string& GetGeneratingString(...)
//------------------------------------------------------------------------------
/**
 * Reconstructs the script line that defines this command
 */
//------------------------------------------------------------------------------
const std::string& Follow::GetGeneratingString(Gmat::WriteMode mode,
                                               const std::string &prefix,
                                               const std::string &useName)
{
   std::string gen = prefix + "Follow " + coordSysName + "(";

   for (UnsignedInt i = 0; i < satNames.size(); ++i)
   {
      if (i > 0) gen += ", ";
      gen += satNames[i];
   }
   gen += ")";

   if (!stopCondText.empty())
      gen += " {" + stopCondText + "}";

   gen += ";";

   generatingString = gen;
   return generatingString;
}


//------------------------------------------------------------------------------
// bool InterpretAction()
//------------------------------------------------------------------------------
/**
 * Parses the generating string to extract configuration.
 *
 * Expected form:
 *   Follow CoordSysName(Sat1 [, Sat2 ...]) [{StopCondText}];
 *
 * Examples:
 *   Follow EarthFixed(Sat1) {Sat1.ElapsedSecs = 3600};
 *   Follow EarthMJ2000Eq(Sat1, Sat2) {Sat1.ElapsedDays = 1.0};
 */
//------------------------------------------------------------------------------
bool Follow::InterpretAction()
{
   // generatingString is set by the interpreter before InterpretAction is called
   std::string gs = generatingString;

   // Strip trailing whitespace and semicolons
   std::string::size_type end = gs.find_last_not_of(" \t\n\r;");
   if (end != std::string::npos)
      gs = gs.substr(0, end + 1);

   // Strip leading "Follow" keyword
   std::string::size_type start = gs.find("Follow");
   if (start == std::string::npos)
      throw CommandException("Follow::InterpretAction: missing 'Follow' keyword in: " + gs);
   gs = gs.substr(start + 6);  // skip "Follow"

   // Trim leading whitespace
   start = gs.find_first_not_of(" \t");
   if (start == std::string::npos)
      throw CommandException("Follow: empty command body");
   gs = gs.substr(start);

   // Extract optional stop condition block first: text between { }
   std::string::size_type lbrace = gs.find('{');
   std::string::size_type rbrace = gs.find('}');
   if (lbrace != std::string::npos && rbrace != std::string::npos && rbrace > lbrace)
   {
      stopCondText = gs.substr(lbrace + 1, rbrace - lbrace - 1);
      // Trim whitespace from stop condition text
      std::string::size_type s = stopCondText.find_first_not_of(" \t");
      std::string::size_type e = stopCondText.find_last_not_of(" \t");
      if (s != std::string::npos)
         stopCondText = stopCondText.substr(s, e - s + 1);
      else
         stopCondText.clear();

      // Remove the brace block from gs so remaining parsing is simpler
      gs = gs.substr(0, lbrace);
      std::string::size_type trailEnd = gs.find_last_not_of(" \t");
      if (trailEnd != std::string::npos)
         gs = gs.substr(0, trailEnd + 1);
   }

   // Extract coordinate system name: token before '('
   std::string::size_type lparen = gs.find('(');
   std::string::size_type rparen = gs.find(')');
   if (lparen == std::string::npos || rparen == std::string::npos || rparen < lparen)
      throw CommandException("Follow: expected CoordSysName(SatName) syntax in: " + gs);

   coordSysName = gs.substr(0, lparen);
   // Trim trailing whitespace from coordSysName
   std::string::size_type csEnd = coordSysName.find_last_not_of(" \t");
   if (csEnd != std::string::npos)
      coordSysName = coordSysName.substr(0, csEnd + 1);
   else
      throw CommandException("Follow: coordinate system name is empty");

   // Extract spacecraft names from inside parentheses (comma-separated)
   std::string satList = gs.substr(lparen + 1, rparen - lparen - 1);
   satNames.clear();

   std::string::size_type pos = 0;
   while (pos < satList.size())
   {
      std::string::size_type comma = satList.find(',', pos);
      std::string token = (comma == std::string::npos)
                           ? satList.substr(pos)
                           : satList.substr(pos, comma - pos);

      // Trim whitespace
      std::string::size_type ts = token.find_first_not_of(" \t");
      std::string::size_type te = token.find_last_not_of(" \t");
      if (ts != std::string::npos)
         satNames.push_back(token.substr(ts, te - ts + 1));

      if (comma == std::string::npos) break;
      pos = comma + 1;
   }

   if (satNames.empty())
      throw CommandException("Follow: no spacecraft specified");

   // Parse the stop condition text
   if (!stopCondText.empty())
      ParseStopCondition(stopCondText);

   return true;
}


//------------------------------------------------------------------------------
// void ParseStopCondition(const std::string &text)
//------------------------------------------------------------------------------
/**
 * Parses a stopping condition string of the form "SatName.ElapsedSecs = N"
 * or "SatName.ElapsedDays = N".
 *
 * Only time-based stopping conditions are supported in this implementation.
 * Non-time conditions (e.g. orbital parameter comparisons) are reserved for a
 * future enhancement.
 */
//------------------------------------------------------------------------------
void Follow::ParseStopCondition(const std::string &text)
{
   // Expect: "Sat1.ElapsedSecs = 3600" or "Sat1.ElapsedDays = 1.0"
   std::string::size_type eqPos = text.find('=');
   if (eqPos == std::string::npos)
      throw CommandException("Follow: stop condition must contain '=': " + text);

   std::string lhs = text.substr(0, eqPos);
   std::string rhs = text.substr(eqPos + 1);

   // Trim whitespace
   auto trim = [](const std::string &s) -> std::string {
      std::string::size_type a = s.find_first_not_of(" \t");
      std::string::size_type b = s.find_last_not_of(" \t");
      return (a == std::string::npos) ? "" : s.substr(a, b - a + 1);
   };

   lhs = trim(lhs);
   rhs = trim(rhs);

   if (lhs.empty() || rhs.empty())
      throw CommandException("Follow: malformed stop condition: " + text);

   // Identify ElapsedSecs vs ElapsedDays
   if (lhs.find("ElapsedDays") != std::string::npos)
   {
      stopIsElapsedDays = true;
      stopParam = lhs;
   }
   else if (lhs.find("ElapsedSecs") != std::string::npos)
   {
      stopIsElapsedDays = false;
      stopParam = lhs;
   }
   else
   {
      throw CommandException(
         "Follow: only ElapsedSecs and ElapsedDays stopping conditions are "
         "supported in this version. Got: " + lhs);
   }

   Real goalValue = 0.0;
   try
   {
      goalValue = std::stod(rhs);
   }
   catch (...)
   {
      throw CommandException("Follow: stop condition goal is not a number: " + rhs);
   }

   // Convert to seconds
   stopGoalSecs = stopIsElapsedDays
                    ? goalValue * GmatTimeConstants::SECS_PER_DAY
                    : goalValue;

   if (stopGoalSecs <= 0.0)
      throw CommandException("Follow: stop condition goal must be positive");

   stopOp = "=";
}


//------------------------------------------------------------------------------
// bool Initialize()
//------------------------------------------------------------------------------
/**
 * Resolves object pointers from the sandbox after the mission sequence is
 * fully populated.  Called by the Sandbox before the first Execute().
 */
//------------------------------------------------------------------------------
bool Follow::Initialize()
{
   if (!GmatCommand::Initialize())
      return false;

   // Resolve the coordinate system
   GmatBase *csObj = FindObject(coordSysName);
   if (csObj == NULL)
      throw CommandException(
         "Follow::Initialize: coordinate system '" + coordSysName + "' not found");
   if (csObj->GetType() != Gmat::COORDINATE_SYSTEM)
      throw CommandException(
         "Follow::Initialize: '" + coordSysName + "' is not a CoordinateSystem");
   fixedCS = (CoordinateSystem*) csObj;

   // Resolve spacecraft
   sats.clear();
   for (const auto &name : satNames)
   {
      GmatBase *satObj = FindObject(name);
      if (satObj == NULL)
         throw CommandException(
            "Follow::Initialize: spacecraft '" + name + "' not found");
      if (satObj->GetType() != Gmat::SPACECRAFT)
         throw CommandException(
            "Follow::Initialize: '" + name + "' is not a Spacecraft");
      sats.push_back((Spacecraft*) satObj);
   }

   // Verify stop condition was provided
   if (stopGoalSecs <= 0.0 && !stopCondText.empty())
   {
      // Re-parse in case the object was cloned without the parsed values
      ParseStopCondition(stopCondText);
   }

   if (stopGoalSecs <= 0.0)
      throw CommandException(
         "Follow::Initialize: no valid stopping condition — "
         "provide e.g. {Sat.ElapsedSecs = 3600}");

   // Verify internal coord system is available (injected by Sandbox)
   if (internalCoordSys == NULL)
      throw CommandException(
         "Follow::Initialize: internal coordinate system is null "
         "(Sandbox has not injected it yet)");

   // Reset runtime state
   frozenStates.clear();
   startEpoch     = 0.0;
   elapsedSecs    = 0.0;
   commandRunning = false;
   stopCondMet    = false;

   // Register published data so subscribers (OFI, GroundTrackPlot, etc.)
   // receive spacecraft state updates during Follow execution.
   SetupPublisher();

   return true;
}


//------------------------------------------------------------------------------
// bool Execute()
//------------------------------------------------------------------------------
/**
 * Advances time by one step and holds each spacecraft at its frozen state in
 * the specified coordinate system.
 *
 * On the first call the current spacecraft state is captured in fixedCS.
 * On every subsequent call the epoch advances by stepSize seconds and each
 * spacecraft's internal state is updated by converting the frozen fixedCS
 * state back to the internal (MJ2000Eq) frame at the new epoch.
 */
//------------------------------------------------------------------------------
bool Follow::Execute()
{
   // First call: capture frozen states and record start epoch
   if (!commandRunning)
   {
      CaptureInitialStates();
      commandRunning = true;
      elapsedSecs    = 0.0;
   }

   // Advance time by one step
   elapsedSecs += stepSize;

   // Clamp to the stop goal so we land exactly on the requested epoch
   if (elapsedSecs > stopGoalSecs)
      elapsedSecs = stopGoalSecs;

   Real newEpoch = startEpoch + elapsedSecs / GmatTimeConstants::SECS_PER_DAY;

   // Move each spacecraft to the frozen position in the new epoch's frame
   UpdateSpacecraftStates(newEpoch);

   // Push state to all subscribers (OFI, GroundTrackPlot, ReportFile, etc.)
   PublishState(newEpoch);

   // Evaluate stop condition
   stopCondMet = CheckStopCondition();

   BuildCommandSummary(false);
   return true;
}


//------------------------------------------------------------------------------
// GmatBase* GetNext()
//------------------------------------------------------------------------------
/**
 * Returns 'this' to re-enter Execute() while running; returns next command
 * once the stop condition is satisfied.
 */
//------------------------------------------------------------------------------
GmatCommand* Follow::GetNext()
{
   if (commandRunning && !stopCondMet)
      return this;

   commandRunning = false;
   return next;
}


//------------------------------------------------------------------------------
// void RunComplete()
//------------------------------------------------------------------------------
/**
 * Cleans up runtime state at the end of a mission run.
 */
//------------------------------------------------------------------------------
void Follow::RunComplete()
{
   frozenStates.clear();
   commandRunning = false;
   stopCondMet    = false;
   elapsedSecs    = 0.0;
   startEpoch     = 0.0;
   if (pubdata)
   {
      delete[] pubdata;
      pubdata     = NULL;
      pubdataSize = 0;
   }
   GmatCommand::RunComplete();
}


//------------------------------------------------------------------------------
// void CaptureInitialStates()
//------------------------------------------------------------------------------
/**
 * Captures each spacecraft's current state expressed in fixedCS and stores it
 * in frozenStates.  Also records the start epoch from the first spacecraft.
 *
 * The spacecraft's internal state is always in the internalCoordSys frame
 * (EarthMJ2000Eq).  We convert to fixedCS to get the "frozen" representation.
 */
//------------------------------------------------------------------------------
void Follow::CaptureInitialStates()
{
   frozenStates.clear();
   frozenStates.reserve(sats.size());

   startEpoch = sats[0]->GetRealParameter("A1Epoch");

   for (auto *sat : sats)
   {
      Real epoch = sat->GetRealParameter("A1Epoch");

      // State in the internal CS (EarthMJ2000Eq)
      Rvector6 internalState(sat->GetState().GetState());

      // Convert to fixedCS at the current epoch
      Rvector6 stateInFixedCS;
      if (!coordConverter.Convert(epoch, internalState, internalCoordSys,
                                  stateInFixedCS, fixedCS))
         throw CommandException(
            "Follow::CaptureInitialStates: coordinate conversion failed for " +
            sat->GetName());

      frozenStates.push_back(stateInFixedCS);

      MessageInterface::ShowMessage(
         "Follow: Captured %s state in %s: [%.6f %.6f %.6f %.6f %.6f %.6f]\n",
         sat->GetName().c_str(), coordSysName.c_str(),
         stateInFixedCS[0], stateInFixedCS[1], stateInFixedCS[2],
         stateInFixedCS[3], stateInFixedCS[4], stateInFixedCS[5]);
   }
}


//------------------------------------------------------------------------------
// void UpdateSpacecraftStates(Real newEpoch)
//------------------------------------------------------------------------------
/**
 * Converts each spacecraft's frozen fixedCS state back to internalCoordSys
 * at the new epoch and writes it to the spacecraft state vector.
 *
 * Because the coordinate system transformation is epoch-dependent (e.g. Earth
 * rotation for EarthFixed, or updated Sun/Moon directions for solar frames),
 * the spacecraft's inertial position will change each step even though its
 * position in fixedCS remains constant.
 */
//------------------------------------------------------------------------------
void Follow::UpdateSpacecraftStates(Real newEpoch)
{
   for (UnsignedInt i = 0; i < sats.size(); ++i)
   {
      Spacecraft *sat = sats[i];
      Rvector6 stateInInternalCS;

      if (!coordConverter.Convert(newEpoch, frozenStates[i], fixedCS,
                                  stateInInternalCS, internalCoordSys))
         throw CommandException(
            "Follow::UpdateSpacecraftStates: coordinate conversion failed for " +
            sat->GetName());

      // Write the new internal-frame state
      sat->SetState(stateInInternalCS);

      // Advance the epoch
      sat->SetRealParameter("A1Epoch", newEpoch);
   }
}


//------------------------------------------------------------------------------
// bool CheckStopCondition() const
//------------------------------------------------------------------------------
/**
 * Returns true when the accumulated elapsed time has reached the stop goal.
 */
//------------------------------------------------------------------------------
bool Follow::CheckStopCondition() const
{
   // Use a small tolerance (0.001 s) to avoid floating-point overshoot issues
   return (elapsedSecs >= stopGoalSecs - 1.0e-3);
}


//------------------------------------------------------------------------------
// void SetupPublisher()
//------------------------------------------------------------------------------
/**
 * Registers the spacecraft state elements with the Publisher so that
 * subscribers (OpenFramesInterface, GroundTrackPlot, ReportFile, etc.) receive
 * data during Follow execution.
 *
 * Data layout published each step:
 *   index 0           : epoch (A1MJD)
 *   index 1 + 6*i     : sat[i] X  (EarthMJ2000Eq, km)
 *   index 2 + 6*i     : sat[i] Y
 *   index 3 + 6*i     : sat[i] Z
 *   index 4 + 6*i     : sat[i] Vx (km/s)
 *   index 5 + 6*i     : sat[i] Vy
 *   index 6 + 6*i     : sat[i] Vz
 */
//------------------------------------------------------------------------------
void Follow::SetupPublisher()
{
   if (publisher == NULL)
      return;

   pubdataOwners.clear();
   pubdataElements.clear();

   pubdataOwners.push_back("All");
   pubdataElements.push_back("All.epoch");

   for (const auto &name : satNames)
   {
      for (int k = 0; k < 6; ++k)
         pubdataOwners.push_back(name);

      pubdataElements.push_back(name + ".X");
      pubdataElements.push_back(name + ".Y");
      pubdataElements.push_back(name + ".Z");
      pubdataElements.push_back(name + ".Vx");
      pubdataElements.push_back(name + ".Vy");
      pubdataElements.push_back(name + ".Vz");
   }

   streamID = publisher->RegisterPublishedData(this, streamID,
                                               pubdataOwners, pubdataElements);

   // Allocate or re-allocate pubdata buffer
   if (pubdata)
   {
      delete[] pubdata;
      pubdata = NULL;
   }
   pubdataSize = 1 + 6 * (Integer)satNames.size();
   pubdata = new Real[pubdataSize];
   for (Integer k = 0; k < pubdataSize; ++k)
      pubdata[k] = 0.0;
}


//------------------------------------------------------------------------------
// void PublishState(Real epoch)
//------------------------------------------------------------------------------
/**
 * Packs the current spacecraft states into pubdata and calls
 * publisher->Publish().  The state values come directly from the spacecraft
 * objects, which were already updated by UpdateSpacecraftStates() to reflect
 * the current epoch in the internalCoordSys (EarthMJ2000Eq).
 */
//------------------------------------------------------------------------------
void Follow::PublishState(Real epoch)
{
   if (publisher == NULL || pubdata == NULL || pubdataSize == 0)
      return;

   pubdata[0] = epoch;

   for (UnsignedInt i = 0; i < sats.size(); ++i)
   {
      const Real *state = sats[i]->GetState().GetState();
      pubdata[1 + 6*i]     = state[0];   // X
      pubdata[2 + 6*i]     = state[1];   // Y
      pubdata[3 + 6*i]     = state[2];   // Z
      pubdata[4 + 6*i]     = state[3];   // Vx
      pubdata[5 + 6*i]     = state[4];   // Vy
      pubdata[6 + 6*i]     = state[5];   // Vz
   }

   publisher->Publish(this, streamID, pubdata, pubdataSize);
}
