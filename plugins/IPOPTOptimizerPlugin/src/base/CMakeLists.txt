# GMAT IPOPTOptimizer Plugin
# Modeled after YukonOptimizerPlugin/src/base/CMakeLists.txt

MESSAGE("==============================")
MESSAGE("GMAT IPOPT Optimizer Plugin setup " ${VERSION})

SET(TargetName IPOPTOptimizer)

# Source directories
SET(PLUGIN_DIRS factory include plugin solver)

# Source files
SET(PLUGIN_SRCS
    factory/IPOPTOptimizerFactory.cpp
    plugin/GmatPluginFunctions.cpp
    solver/IPOPTOptimizer.cpp
)

# Set up the plugin (adds include dirs, links GmatBase, etc.)
_SETUPPLUGIN(${TargetName} "${PLUGIN_DIRS}" "${PLUGIN_SRCS}" plugins)

# DLL export symbol
SET_TARGET_PROPERTIES(${TargetName} PROPERTIES DEFINE_SYMBOL "IPOPT_EXPORTS")

# IPOPT include and link
TARGET_INCLUDE_DIRECTORIES(${TargetName} PRIVATE ${IPOPT_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(${TargetName} PRIVATE ${IPOPT_LIBRARY})

# C++11 threads
if(UNIX)
    TARGET_LINK_LIBRARIES(${TargetName} PRIVATE Threads::Threads)
endif()
