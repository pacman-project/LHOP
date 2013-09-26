# -*- Mode: CMake; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- *
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(OPENCL_ENABLE FALSE CACHE BOOL "Enable OpenCL support. Set 'use_opencl = true' and appropriate 'opencl_kernel_paths' and 'opencl_devices' in config file to use it.")
SET(OPENMP_ENABLE TRUE CACHE BOOL "Enable OpenMP support")

FIND_PACKAGE(OpenCV REQUIRED)
FIND_PACKAGE(TBB REQUIRED)
FIND_PACKAGE(ZLIB REQUIRED)
FIND_PACKAGE(OpenMP) 

IF (OPENMP_ENABLE)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF (OPENMP_ENABLE)

IF (OPENCL_ENABLE) 
	SET(OPENCL_ENABLE_RESULT_VERIFICATION FALSE CACHE BOOL "Enable support for verifcation of opencl result aginst original results. Set 'opencl_verify_result = true' in config file to use it.")
	
	IF (OPENCL_ENABLE_RESULT_VERIFICATION)
		ADD_DEFINITIONS(-DOPENCL_VERIFICATION_SUPPORT)
	ENDIF (OPENCL_ENABLE_RESULT_VERIFICATION)
	
	ADD_DEFINITIONS(-DOPENCL)
	FIND_PACKAGE(OpenCL REQUIRED)
	LINK_LIBRARIES(${OPENCL_LIBRARIES})
	INCLUDE_DIRECTORIES(AFTER ${OPENCL_INCLUDE_DIR})
ENDIF (OPENCL_ENABLE)

LINK_LIBRARIES(${OpenCV_LIBS} ${TBB_LIBRARIES} ${ZLIB_LIBRARIES})
LINK_DIRECTORIES(${TBB_LIBRARY_DIRS})

INCLUDE_DIRECTORIES(AFTER ${OpenCV_INCLUDE_DIR} ${TBB_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIR} .)

LIST(APPEND LIB_SOURCES utils/structures.cpp utils/utils.cpp utils/ocv.cpp utils/atom.cpp utils/img.cpp utils/misc.cpp utils/streaming.cpp utils/zlib_compress.cpp utils/matching/matching.cpp utils/matching/munkres.cpp utils/hopmath.cpp)

LIST(APPEND LIB_SOURCES layers/layer_1_creators.cpp layers/initialization.cpp layers/layers.cpp layers/layer_1.cpp layers/layer_1_result.cpp layers/part_lib.cpp layers/layer_n_creators.cpp layers/optimization.cpp layers/layer_learning.cpp layers/ocl_layer_1_creators.cpp layers/ocl_layer_n_creators.cpp layers/hoc.cpp)

LIST(APPEND LIB_SOURCES graphs/graph.cpp graphs/img_graph.cpp graphs/graph_utils.cpp)

LIST(APPEND LIB_HEADERS utils/platform.h utils/structures.h utils/matrix.h utils/utils.h utils/convert.h utils/ocv.h utils/atom.h utils/fonts.h utils/img.h utils/misc.h utils/streaming.h utils/zlib_compress.h utils/smemory.h utils/exceptions.h utils/matching/matching.h utils/matching/matrix.h utils/matching/munkres.h utils/hopmath.h utils/mapreduce.h)

LIST(APPEND LIB_HEADERS layers/layer_1_creators.h layers/initialization.h layers/layers.h layers/layer_1.h layers/layer_1_result.h layers/part_lib.h layers/layer_n_creators.h layers/optimization.h layers/layer_learning.h layers/hoc.h)

LIST(APPEND LIB_HEADERS graphs/graph.h graphs/img_graph.h graphs/graph_utils.h)

LIST(APPEND LIB_HEADERS opencl/clContextManager.h)
LIST(APPEND LIB_SOURCES opencl/clContextManager.cpp)
	
IF (OPENCL_ENABLE) 
	LIST(APPEND LIB_HEADERS opencl/cl_structures.h opencl/cl_utils.h opencl/clContextManager.h opencl/kernel_paths.h )
	LIST(APPEND LIB_SOURCES opencl/cl_utils.cpp opencl/clContextManager.cpp)
ENDIF (OPENCL_ENABLE)

IF (WIN32)
    LINK_LIBRARIES(Psapi)
    LIST(APPEND LIB_HEADERS interface/targetver.h)
    # Windows dlls requre special export declarations for external functions
    # ... we do not have them yet.
	# Also create winhop_set_paths.bat files for winhop
	IF (CMAKE_GENERATOR MATCHES Win64)
		CONFIGURE_FILE(winhop_set_paths_x64.bat.in ${PROJECT_BINARY_DIR}/winhop_set_paths_x64.bat @ONLY)
	ELSE (CMAKE_GENERATOR MATCHES Win64)
		CONFIGURE_FILE(winhop_set_paths_x86.bat.in ${PROJECT_BINARY_DIR}/winhop_set_paths_x86.bat @ONLY)
	ENDIF (CMAKE_GENERATOR MATCHES Win64)
#    
#ELSE (WIN32)
ENDIF (WIN32)

ADD_LIBRARY(hopstatic STATIC ${LIB_SOURCES} ${LIB_HEADERS})

SOURCE_GROUP ("Header Files" FILES ${LIB_HEADERS})

ADD_SUBDIRECTORY(interface)
ADD_SUBDIRECTORY(tools)


