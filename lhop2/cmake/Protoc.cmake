# -*- Mode: CMake; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */

FUNCTION(PROTOBUF_GEN_CPP SRCS HDRS)
      if(NOT ARGN)
        message(SEND_ERROR "Error: PROTOBUF_GENERATE() called without any proto files")
        return()
      endif(NOT ARGN)

      set(${SRCS})
      set(${HDRS})
      foreach(FIL ${ARGN})
        get_filename_component(ABS_FIL ${FIL} ABSOLUTE)
        get_filename_component(FIL_WE ${FIL} NAME_WE)
        get_filename_component(FIL_FULL ${FIL} NAME)
        GET_FILENAME_COMPONENT(FILEPATH ${ABS_FIL} PATH )        

        list(APPEND ${SRCS} "${CMAKE_CURRENT_BINARY_DIR}/${FIL_WE}.pb.cc")
        list(APPEND ${HDRS} "${CMAKE_CURRENT_BINARY_DIR}/${FIL_WE}.pb.h")

        add_custom_command(
          OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${FIL_WE}.pb.cc"
                 "${CMAKE_CURRENT_BINARY_DIR}/${FIL_WE}.pb.h"
          COMMAND  ${PROTOBUF_PROTOC_EXECUTABLE}
          ARGS --cpp_out=${CMAKE_CURRENT_BINARY_DIR} ${FIL_FULL}
          DEPENDS ${ABS_FIL}
          COMMENT "Running protocol buffer compiler on ${FIL_FULL} for C++ interface"
          WORKING_DIRECTORY ${FILEPATH}
          VERBATIM )
      ENDFOREACH()

      SET_SOURCE_FILES_PROPERTIES(${${SRCS}} ${${HDRS}} PROPERTIES GENERATED TRUE)
      SET(${SRCS} ${${SRCS}} PARENT_SCOPE)
      SET(${HDRS} ${${HDRS}} PARENT_SCOPE)
ENDFUNCTION()

FUNCTION(PROTOBUF_GEN_JAVA SRCS PACKAGE)
      INCLUDE(Java)

      IF(NOT ARGN)
        MESSAGE(SEND_ERROR "Error: PROTOBUF_GENERATE() called without any proto files")
        RETURN()
      ENDIF(NOT ARGN)

      SET(${SRCS})

      FOREACH(FIL ${ARGN})
        GET_FILENAME_COMPONENT(ABS_FIL ${FIL} ABSOLUTE)
        GET_FILENAME_COMPONENT(FIL_WE ${FIL} NAME_WE)
        GET_FILENAME_COMPONENT(FIL_FULL ${FIL} NAME)
        GET_FILENAME_COMPONENT(FILEPATH ${ABS_FIL} PATH )        

        STRING(REPLACE "." "/" DIR ${PACKAGE})

        CAPITALIZE(${FIL_WE} FIL_WE)
        SET(JAVA_SRC "${CMAKE_CURRENT_BINARY_DIR}/${DIR}/${FIL_WE}.java")

        LIST(APPEND ${SRCS} "${JAVA_SRC}")

        ADD_CUSTOM_COMMAND(
          OUTPUT "${JAVA_SRC}"
          COMMAND  ${PROTOBUF_PROTOC_EXECUTABLE}
          ARGS --java_out=${CMAKE_CURRENT_BINARY_DIR} ${FIL_FULL}
          DEPENDS ${ABS_FIL}
          COMMENT "Running protocol buffer compiler on ${FIL_FULL} for Java interface"
          WORKING_DIRECTORY ${FILEPATH}
          VERBATIM )
      ENDFOREACH()

      SET_SOURCE_FILES_PROPERTIES(${${SRCS}} PROPERTIES GENERATED TRUE)
      SET(${SRCS} ${${SRCS}} PARENT_SCOPE)
       
ENDFUNCTION()

