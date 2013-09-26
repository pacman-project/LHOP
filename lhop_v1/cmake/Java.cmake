# -*- Mode: CMake; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */

FUNCTION(JAVA_FIND_JAR JAR)
    IF(NOT ARGN)
        MESSAGE(SEND_ERROR "Error: At least one possible name of the jar required")
        RETURN()
    ENDIF(NOT ARGN)
    IF (JAVA_SEARCH_DIRS)
        FIND_FILE(${JAR} ${ARGN} PATHS ENV CLASSPATH ENV PATH ENV JAVA_DIRS ${JAVA_SEARCH_DIRS})
    ELSE ()
        FIND_FILE(${JAR} ${ARGN} PATHS ENV CLASSPATH ENV PATH ENV JAVA_DIRS )
    ENDIF()
ENDFUNCTION()

FUNCTION(JAVA_GET_PACKAGE SOURCE_FILE PACKAGE)
    FILE(STRINGS SOURCE_FILE A LIMIT_COUNT 1 REGEX "^package +[a-zA-Z\\._0-9]+;" NO_HEX_CONVERSION)
    STRING(REGEX_MATCH "[a-zA-Z\\._0-9]" "${A}" ${PACKAGE})
ENDFUNCTION()

FUNCTION(JAVA_DEPENDS SRCPATH DSTPATH CLASS_LIST)
    FILE(MAKE_DIRECTORY ${DSTPATH})
    FILE(GLOB_RECURSE JAVA_FILES_FULL ${SRCPATH}/*.java)
    FOREACH(SRCFILE ${JAVA_FILES_FULL})
        GET_FILENAME_COMPONENT(FILENAME ${SRCFILE} NAME_WE )
        GET_FILENAME_COMPONENT(FILEPATH ${SRCFILE} PATH )
        FILE_RELATIVE(${SRCPATH} ${FILEPATH} RELPATH)
        #STRING(REPLACE "${SRCPATH}" "${DSTPATH}" CLASSPATH ${FILEPATH})
        LIST(APPEND CLASSLIST ${DSTPATH}/${RELPATH}/${FILENAME}.class)
    ENDFOREACH(SRCFILE ${JAVA_FILES_FULL} )
    SET(${CLASS_LIST} ${CLASSLIST} PARENT_SCOPE)
ENDFUNCTION()

FUNCTION(JAVA_LIST_CLASSPATH LIST CLASSPATH)
    IF (UNIX)
        STRING(REPLACE ";" ":" B "${LIST}")
    ELSE ()
        SET(B "${LIST}")
    ENDIF ()
    SET(${CLASSPATH} ${B} PARENT_SCOPE)
ENDFUNCTION()


FUNCTION(CAPITALIZE IN OUT)
    STRING(SUBSTRING ${IN} 0 1 A)
    STRING(LENGTH ${IN} L)
    MATH(EXPR LL ${L}-1)
    STRING(SUBSTRING ${IN} 1 ${LL} B)
    STRING(TOUPPER ${A} C)
    SET(${OUT} "${C}${B}" PARENT_SCOPE)
ENDFUNCTION()

FUNCTION(STRING_STARTS_WITH FULL START RESULT)
    STRING(LENGTH ${FULL} FL)
    STRING(LENGTH ${START} SL)
    IF (${FL} LESS ${SL})
        SET(${RESULT} FALSE)
        RETURN()
    ENDIF()
    STRING(SUBSTRING ${FULL} 1 ${SL} B)
    IF (${START} STREQUAL ${B})
        SET(${RESULT} TRUE)
    ELSE()
        SET(${RESULT} FALSE)
    ENDIF()
ENDFUNCTION()

FUNCTION(FILE_RELATIVE PARENT CHILD RELATIVE)
    GET_FILENAME_COMPONENT(PARENT_A ${PARENT} ABSOLUTE)
    GET_FILENAME_COMPONENT(CHILD_A ${CHILD} ABSOLUTE)
    STRING_STARTS_WITH("${CHILD_A}" "${PARENT_A}" RES)
    IF (NOT ${RES})
        MESSAGE(SEND_ERROR "Error: ${CHILD} is not a descendant directory of ${PARENT}")
        RETURN()
    ENDIF()
    STRING(LENGTH ${PARENT_A} PL)
    STRING(LENGTH ${CHILD_A} CL)
    MATH(EXPR LL ${CL}-${PL}-1)
    MATH(EXPR LS ${PL}+1)
    STRING(SUBSTRING ${CHILD_A} ${LS} ${LL} B)
    SET(${RELATIVE} ${B} PARENT_SCOPE)

ENDFUNCTION ()
