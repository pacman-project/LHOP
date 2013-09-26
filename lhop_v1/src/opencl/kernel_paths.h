#pragma once

#ifndef __OPENCL_KERNEL_PATH
#define __OPENCL_KERNEL_PATH

#include <string>

using namespace std;

static const string getDirname(const string file, const string suffix) {
	int found = file.find_last_of("/\\");
	return file.substr(0,found) + suffix;	
}

static const int opencl_program_listing_size = 8;

/**
 * List of all *.cl files i.e. cl_program objects that are created from this file.
 *
 * !!!!!!!!!
 * CAUTION: When adding or deleting !! DO NOT FORGET !!! to update length by changing opencl_program_listing_size.
 * !!!!!!!!!
 *
 * Be careful to use / slash and not \\ or \ since OpenCL preprocessor compile will cause problems
 */
/*static const string opencl_program_listing[] = {
	getDirname(__FILE__, "/kernels/parallel_scan.cl"),
	getDirname(__FILE__, "/kernels/create_layer0.cl"),	
	getDirname(__FILE__, "/kernels/add_layer1.cl"),	
	getDirname(__FILE__, "/kernels/convolution.cl"),
	getDirname(__FILE__, "/kernels/convolution4.cl"),
	getDirname(__FILE__, "/kernels/convolution8.cl"),
	getDirname(__FILE__, "/kernels/convolution16.cl"),
};*/

static const string opencl_program_listing[] = {
	"add_layer1.cl",
	"create_layer0.cl",	
	"parallel_sort.cl",
	"parallel_scan.cl",		
	"convolution.cl",
	"convolution4.cl",
	"convolution8.cl",
	"convolution16.cl",	
};

#endif