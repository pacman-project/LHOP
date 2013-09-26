// layncreate.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <stdio.h>
#include "../../interface/hop.h"

#if defined WIN32 | defined WIN64
#include <Windows.h>
#else
#include <dlfcn.h>
#endif

typedef bool (*hop_n_inference_func)(char*, char*, char*);
typedef char* (*hop_time_stamp_func)();

using namespace std;

int main(int argc, char* argv[])
{
    //_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
    //_CrtSetBreakAlloc(20067884);

	hop_n_inference_func ref_hop_n_inference = (hop_n_inference_func)hop_n_inference;
	hop_time_stamp_func ref_hop_time_stamp = (hop_time_stamp_func)hop_time_stamp;
	
	bool quit_on_error = false;
	int argv_offset = 1;
	for (int i = 1; i < argc; i++) {
		const string argv_1(argv[i]);
		if (argv_1.find("--use-dll=") == 0) {			
			
			argv_offset++;
			const string name = argv_1.substr(10);

#if defined WIN32 | defined WIN64
			HINSTANCE dll_handle = LoadLibrary(name.c_str());
#else
		void* dll_handle = dlopen(name.c_str(), RTLD_LAZY);
#endif
			if (dll_handle == nullptr) {
				cout << "Error: Unable to find DLL '" << name.c_str() << "'" << endl;
				return 0;
			}
#if defined WIN32 | defined WIN64
			ref_hop_n_inference = (hop_n_inference_func)GetProcAddress(dll_handle, "hop_n_inference");
			ref_hop_time_stamp = (hop_time_stamp_func)GetProcAddress(dll_handle, "hop_time_stamp");
#else
			ref_hop_n_inference = (hop_n_inference_func)dlsym(dll_handle, "hop_n_inference");
			ref_hop_time_stamp = (hop_time_stamp_func)dlsym(dll_handle, "hop_time_stamp");
#endif
		} else if (argv_1.find("--quit-on-error") == 0) {
			argv_offset++;
			quit_on_error = true;
		}
	}
		
	cout << "layer n creator " << ref_hop_time_stamp() << endl;
	bool ok = true;
	switch (argc - argv_offset) {
		case 2:
			ok = ref_hop_n_inference(argv[argv_offset + 0], argv[argv_offset + 1], "");
			break;
		case 3:
			ok = ref_hop_n_inference(argv[argv_offset + 0], argv[argv_offset + 1], argv[argv_offset + 2]);
			break;
		default:			
			cout << "Usage: hopncreate [--use-dll=<dll_name>] [--quit-on-error] config_file file_pattern." << endl;
	}
	if (ok == false && quit_on_error == false ) {
		cout << "Press any key to exit ..." << endl;
		getchar();
	}
    return ok ? 0 : 1;
}
