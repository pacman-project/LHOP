// dllmain.cpp : Defines the entry point for the DLL application.

#ifdef WIN32
#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
#include <windows.h>
#endif

#include "layers/initialization.h"

#ifdef WIN32

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
					 )
{
	switch (ul_reason_for_call)
	{
	case DLL_PROCESS_ATTACH:
	case DLL_THREAD_ATTACH:
        init_atoms();
        init_streaming();
        break;
	case DLL_THREAD_DETACH:
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}

#else

void __attribute__ ((constructor)) hop_so_load(void);
void __attribute__ ((destructor)) hop_so_unload(void);

#include <execinfo.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>

void error_sig_handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, 2);
  exit(1);
}


// Called when the library is loaded and before dlopen() returns
void hop_so_load(void)
{
	// register with handling for error signals	
	//signal(SIGSEGV, error_sig_handler);

    init_atoms();
	init_streaming();        
}

// Called when the library is unloaded and before dlclose()
// returns
void hop_so_unload(void)
{
    // Add clean-up code
}

#endif


