# GMAT: General Mission Analysis Tool.
#
# CMAKE script file for the GMAT FollowCommand Plugin.
# This file must be installed in plugins/FollowCommandPlugin/src/base/
#
# Author: Michael Mahoney, NASA/GSFC
# Created: 2026-04-10

MESSAGE("==============================")
MESSAGE("GMAT FollowCommand Plugin setup " ${VERSION})

SET(TargetName FollowCommand)

# ====================================================================
# Directories containing source/header files
SET(PLUGIN_DIRS command factory include plugin)

# ====================================================================
# Source files
SET(PLUGIN_SRCS
    command/Follow.cpp
    factory/FollowCommandFactory.cpp
    plugin/GmatPluginFunctions.cpp
)

# ====================================================================
# Setup common plugin definitions, targets, and install rules
_SETUPPLUGIN(${TargetName} "${PLUGIN_DIRS}" "${PLUGIN_SRCS}" plugins)

# Rename the DLL export preprocessor definition to match FollowDefs.hpp
SET_TARGET_PROPERTIES(${TargetName} PROPERTIES DEFINE_SYMBOL "FOLLOW_EXPORTS")
