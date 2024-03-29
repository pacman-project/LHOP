# -*- Mode: CMake; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- *
PROJECT (libhop)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
GET_FILENAME_COMPONENT(PROJECT_ROOT . ABSOLUTE)

INCLUDE(presets.cmake OPTIONAL)

# On 64bit Linux machines libraries should be compiled using -fPIC
IF(UNIX AND NOT WIN32)
  FIND_PROGRAM(CMAKE_UNAME uname /bin /usr/bin /usr/local/bin )
  
  IF(CMAKE_UNAME)
    EXEC_PROGRAM(uname ARGS -m OUTPUT_VARIABLE CMAKE_SYSTEM_PROCESSOR)
	SET(LINUX_ENABLE_OPTIMIZATIONS TRUE CACHE BOOL "All possible optimizations for linux")
    SET(CMAKE_SYSTEM_PROCESSOR ${CMAKE_SYSTEM_PROCESSOR} CACHE INTERNAL
"processor type (i386 and x86_64)")
    IF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
	ADD_DEFINITIONS(-fPIC)
    ELSE(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
	ADD_DEFINITIONS(-fPIC)
    ENDIF(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
	IF(LINUX_ENABLE_OPTIMIZATIONS)
		ADD_DEFINITIONS(-O3)
		ADD_DEFINITIONS(-march=native)
		ADD_DEFINITIONS(-mtune=native)
		ADD_DEFINITIONS(-fomit-frame-pointer)
		ADD_DEFINITIONS(-msse)
		ADD_DEFINITIONS(-msse2)
	ENDIF()
	ADD_DEFINITIONS(-std=c++0x)
  ENDIF(CMAKE_UNAME)
ENDIF(UNIX AND NOT WIN32)

IF (WIN32)
	IF (CMAKE_GENERATOR MATCHES Win64)
		ADD_DEFINITIONS(-DWIN64)
		SET(TBB_ARCHITECTURE "intel64" CACHE STRING "Processor architecture")
	ELSE (CMAKE_GENERATOR MATCHES Win64)
		SET(TBB_ARCHITECTURE "ia32" CACHE STRING "Processor architecture")
	ENDIF (CMAKE_GENERATOR MATCHES Win64)

    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP" )
ENDIF (WIN32)

IF (NOT WIN32)
# Disable deprecation warnings on GCC
ADD_DEFINITIONS(-Wno-deprecated)
ENDIF (NOT WIN32)

IF (DEBUG_DEFAULT)
        ADD_DEFINITIONS(-DDEBUG_DEFAULT)
ENDIF(DEBUG_DEFAULT)

SET(CMAKE_MODULE_PATH ${PROJECT_ROOT}/cmake)

SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

SET (RESOURCES_DIR ${PROJECT_ROOT}/res)

# add a target to generate API documentation with Doxygen
 FIND_PACKAGE(Doxygen)
 
 IF(DOXYGEN_FOUND)
 CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doxyfile @ONLY)
 ADD_CUSTOM_TARGET(doc
	${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Generating API documentation with Doxygen" VERBATIM)
 ENDIF(DOXYGEN_FOUND)

ADD_SUBDIRECTORY(src)




