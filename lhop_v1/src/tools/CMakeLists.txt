# -*- Mode: CMake; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- *
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(TOOLS TRUE CACHE BOOL "Make command line tools")
SET(BUILD_TEST_TOOL FALSE CACHE BOOL "Make test command line tool")

IF (TOOLS)

ADD_EXECUTABLE(hop1create lay1create/lay1create.cpp)
TARGET_LINK_LIBRARIES(hop1create hop)
CONFIGURE_FILE(lay1create/lay1config.cfg.in ${PROJECT_BINARY_DIR}/lay1config.cfg @ONLY)

ADD_EXECUTABLE(hopncreate layncreate/layncreate.cpp)
TARGET_LINK_LIBRARIES(hopncreate hop)

ADD_EXECUTABLE(hopdisplay laydisplay/laydisplay.cpp)
TARGET_LINK_LIBRARIES(hopdisplay hopstatic)

ADD_EXECUTABLE(hoplearning laylearning/laylearning.cpp)
TARGET_LINK_LIBRARIES(hoplearning hop)

ADD_EXECUTABLE(hopserver layserver/layserver.cpp)
TARGET_LINK_LIBRARIES(hopserver hop)

ENDIF (TOOLS)

IF (BUILD_TEST_TOOL)

ADD_EXECUTABLE(test test/test.cpp)
TARGET_LINK_LIBRARIES(test hopstatic)

ENDIF (BUILD_TEST_TOOL)

ADD_EXECUTABLE(apisample apisample/apisample.cpp)
TARGET_LINK_LIBRARIES(apisample hop)


IF (WIN32)
    IF (TOOLS)
        # hack to get around the "Debug" and "Release" directories cmake tries to add on Windows
        SET_TARGET_PROPERTIES (hop1create hopncreate hopdisplay hoplearning hopserver PROPERTIES RELEASE_PREFIX ${WIN32_INSTALL})
        INSTALL(TARGETS hop1create hopncreate hopdisplay hoplearning hopserver RUNTIME DESTINATION ${WIN32_INSTALL})
    ENDIF (TOOLS)
ELSE (WIN32)
#   INSTALL(TARGETS hop1create hopncreate hopdisplay hoplearning RUNTIME DESTINATION bin)
ENDIF (WIN32)

