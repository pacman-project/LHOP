# -*- Mode: CMake; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- *
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

SET(TOOLSET TRUE CACHE BOOL "Make command line toolset")

FILE(GLOB MAIN_TOOLSET_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
FILE(GLOB MAIN_TOOLSET_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

ADD_LIBRARY(static_main_toolset OBJECT ${MAIN_TOOLSET_SOURCES} ${MAIN_TOOLSET_HEADERS})

###ADD_SUBDIRECTORY(interface) 

IF (TOOLSET)

ADD_EXECUTABLE(hopdisplay laydisplay/laydisplay.cpp 
				$<TARGET_OBJECTS:static_core> 
				$<TARGET_OBJECTS:static_utils> 
				$<TARGET_OBJECTS:static_main_toolset>
				$<TARGET_OBJECTS:static_modules_hoc> 
				$<TARGET_OBJECTS:static_modules_shapes> 
				)

LIST(APPEND LHOP_SOURCES lhop/main_lhop.cpp)

ADD_EXECUTABLE(lhop ${LHOP_SOURCES}
				$<TARGET_OBJECTS:static_core> 
				$<TARGET_OBJECTS:static_utils> 
				$<TARGET_OBJECTS:static_main_toolset>
				$<TARGET_OBJECTS:static_modules_hoc> 
				$<TARGET_OBJECTS:static_modules_shapes> 
				)

ENDIF (TOOLSET)

###ADD_EXECUTABLE(apisample apisample/apisample.cpp)
###TARGET_LINK_LIBRARIES(apisample shared_hop)


IF (WIN32)
    IF (TOOLSET)
        # hack to get around the "Debug" and "Release" directories cmake tries to add on Windows		
        SET_TARGET_PROPERTIES (lhop hopdisplay PROPERTIES RELEASE_PREFIX ${WIN32_INSTALL})
        INSTALL(TARGETS lhop hopdisplay RUNTIME DESTINATION ${WIN32_INSTALL})
    ENDIF (TOOLSET)
ELSE (WIN32)
#   INSTALL(TARGETS hop1create hopncreate hopdisplay RUNTIME DESTINATION bin)
ENDIF (WIN32)

