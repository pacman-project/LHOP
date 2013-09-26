#ifndef __CL_UTILS__
#define __CL_UTILS__

#ifdef OPENCL

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <CL/opencl.h>
#include "clContextManager.h"

void ocl_set_memory_zero(const best_cmd_queue_info& cmd_queue_info, cl_mem input_mem, cl_int input_size, cl_int unalinged_input_size, cl_int wait_event_count, cl_event* wait_evnt, cl_event* evnt) ;

// util functions (some copied form nvidia SDK)

// *********************************************************************
// Full error handling macro with Cleanup() callback (if supplied)... 
// (Companion Inline Function lower on page)

#define ocl_check_error_ex(a, b, c) __ocl_check_error_ex(a, b, c, __FILE__ , __LINE__) 

// Short version without Cleanup() callback pointer
// Both Input (a) and Reference (b) are specified as args
#define ocl_check_error(a, b) ocl_check_error_ex(a, b, 0) 

//////////////////////////////////////////////////////////////////////////////
//! Loads a Program file and prepends the cPreamble to the code.
//!
//! @return the source string if succeeded, 0 otherwise
//! @param cFilename        program filename
//! @param cPreamble        code that is prepended to the loaded file, typically a set of #defines or a header
//! @param szFinalLength    returned length of the code string
//////////////////////////////////////////////////////////////////////////////
char* ocl_load_prog_source(const char* cFilename, const char* cPreamble, size_t* szFinalLength);

// Helper function to get OpenCL error string from constant
// *********************************************************************
const char* ocl_error_string(cl_int error);

inline void __ocl_check_error_ex(cl_int iSample, cl_int iReference, void (*pCleanup)(int), const char* cFile, const int iLine)
{
    // An error condition is defined by the sample/test value not equal to the reference
    if (iReference != iSample)
    {
        // If the sample/test value isn't equal to the ref, it's an error by defnition, so override 0 sample/test value
        iSample = (iSample == 0) ? -9999 : iSample; 

        // Log the error info
		std::cout << "\n !!! Error # " << iSample << " (" << ocl_error_string(iSample) << ") at line " << iLine << " , in file " << cFile << " !!!" << std::endl << std::endl;

        // Cleanup and exit, or just exit if no cleanup function pointer provided.  Use iSample (error code in this case) as process exit code.
        if (pCleanup != NULL)
        {
            pCleanup(iSample);
        }
        else 
        {
			std::cout << "Press any key to exit ..." << std::endl;
			getchar();
            exit(iSample);
        }
    }
}

#endif 

#endif