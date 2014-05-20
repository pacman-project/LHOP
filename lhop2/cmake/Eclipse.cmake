# -*- Mode: CMake; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */

INCLUDE(Java)

FUNCTION(ECLIPSE_GENERATE_PROJECT DIRECTORY NAME SRC_DIRS BIN_DIR DEPENDENCIES)

    FILE_RELATIVE(${DIRECTORY} ${BIN_DIR} BIN_DIR_R)
    SET(A "<classpathentry kind=\"output\" path=\"${BIN_DIR_R}\"/>\n")

    FOREACH(SRC_DIR ${SRC_DIRS})
        FILE_RELATIVE(${DIRECTORY} ${SRC_DIR} SRC_DIR_R)
        SET(A "${A}\n<classpathentry kind=\"src\" path=\"${SRC_DIR_R}\" />")
    ENDFOREACH()

    FOREACH(DEPENDENCY ${DEPENDENCIES})
        IF (NOT IS_DIRECTORY ${DEPENDENCY})
            SET(A "${A}\n<classpathentry kind=\"lib\" path=\"${DEPENDENCY}\" />")
        ELSE()

        ENDIF()
    ENDFOREACH()

    FILE(WRITE ${DIRECTORY}/.classpath "<?xml version=\"1.0\" encoding=\"UTF-8\"?><classpath><classpathentry kind=\"con\" path=\"org.eclipse.jdt.launching.JRE_CONTAINER\"/>\n${A}\n</classpath>")


    FILE(WRITE ${DIRECTORY}/.project "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<projectDescription>
    <name>${NAME}</name>
    <comment></comment>
    <projects>
    </projects>
    <buildSpec>
        <buildCommand>
			<name>org.eclipse.jdt.core.javabuilder</name>
			<arguments>
			</arguments>
		</buildCommand>
	</buildSpec>
	<natures>
		<nature>org.eclipse.jdt.core.javanature</nature>
	</natures>
</projectDescription>")

ENDFUNCTION()


