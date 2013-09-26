#pragma once

#ifndef __CL_CONTEXT_MANAGER
#define __CL_CONTEXT_MANAGER 1

#include <vector>
#include <string>

using namespace std;

#ifdef OPENCL

#include <iostream>

#include <list>
#include <map>
#include <set>

#include <algorithm>
#include <CL/opencl.h>

#include "kernel_paths.h"

struct best_cmd_queue_info {
	// how many max work-items per work-group can this device_id have
	size_t max_work_group_size;

	// pre-calculated best allocations of max work-items per group based on dimensions (no need for dim1 since it will be same value as max_work_group_size)
	struct {
		int x, y;
		int pow_x, pow_y;
	} best_work_groups_dim2;

	struct {
		int x, y, z;
		int pow_x, pow_y, pow_z;
	} best_work_groups_dim3;

	/**
	 * Some other device info that might be useful for best optimization
	 */
	// does device_id support images  	
	cl_bool has_img_support;

	// type and size of local memory
	cl_device_local_mem_type local_mem_type;
	cl_ulong local_mem_size;

	// number of compute units
	cl_int max_compute_units;
	
	// type of device (CPU, GPU, APU or default)
	cl_device_type device_type;

	unsigned long long max_mem_alocation;

	cl_device_id device_id;
	int context_number;
	cl_command_queue queue;		
};

// main context manager class that must be able to create contextes for all platforms and devices 
// and must be able to return best command queue for user to use
class clContextManager
{
private:
	// set this value to true when you want only one device in one context (i.e. each device will have its own context)
	static const bool seperate_contex_for_devices = true;
	int counter;

	vector<string> kernel_paths;
	string opencl_include_paths;

	struct t_context {
		cl_context c;
		int context_number;
	};

	// list of all contexes
	vector<t_context> contextList;


	// storage for quick access to command queue and kernels

	// mapping of kernel_type number to array of all kernels for each context
	// stored by context_number (application can get kernel list first, and when
	// it gets assigned command queue it can also quickly find correct kernel by using context_number)
	map<string,vector<cl_kernel> > kernelList;

	// array of all queues with saved cl_device and context_number
	vector<best_cmd_queue_info> commandQueueList;

	// list of all program, we must know them so we can destroy them on the end
	list<cl_program> programList;

	/**
	 * Build all programs for selected context and saves them to programsList.
	 * Programs can be loaded from source and build or it can be loaded from binary file if @by_binary == true.
	 * (by default @by_binary is false)
	 */
	void buildProgramsAndKernels(t_context context, cl_device_id* device_list, int device_list_size, bool by_binary = false);
	cl_program compileWithSource(const t_context context, const cl_device_id* device_list, const int device_list_size, const char * programSource, const size_t progLength, const char* program_name);
	cl_program compileWithBinary(const t_context context, const cl_device_id* device_list, const int device_list_size, const char * programSource, const size_t progLength, const char* program_name);

	static string getDeviceName(const cl_device_id device);
	static string getKernelName(const cl_kernel kernel);

	static string createCompileOptFromDeviceCapabilities(const cl_device_id device);
	static void calculateBestWorkGroupSizes(best_cmd_queue_info* cmd_queue);
public:

	bool is_initialized;

	clContextManager(void);
	~clContextManager(void);

	void initializeContexts(vector<string> &kernels_path, vector<string>& allowed_devices, bool force_reinitialize = false);

	const best_cmd_queue_info& getBestCommandQueue(const string& preferred_device = string(), const vector<string>& use_only_devices = vector<string>());

	void addKernelPath(const string &krnl_path) {
		kernel_paths.push_back(krnl_path);
	}

	vector<cl_kernel> getCommonKernel(const string& name) {
		return kernelList[name];
	}

	cl_context getContext(const int context_num) {
		return contextList[context_num].c;
	}

	list<cl_context> getAllContextes() {
		list<cl_context> result;
		for (int i = 0; i < contextList.size(); i++) {
			if (contextList[i].context_number == i)
				result.push_back(contextList[i].c);
		}
	}
};


// make getBestCommandQueue inline function to avoid to many function calls
inline const best_cmd_queue_info& clContextManager::getBestCommandQueue(const string& preferred_device, const vector<string>& use_only_devices)
{
	// currently just get first command queue from list
	if (commandQueueList.empty()) {
		cout << "No command queue found - maybe there is no OpenCL device present?" << endl;
		throw new std::exception();
	} else {

		cl_device_type pref_device_type = CL_DEVICE_TYPE_ALL;

		if (preferred_device.length() > 0) {
			string pref_dev = preferred_device;
			std::transform(pref_dev.begin(), pref_dev.end(), pref_dev.begin(), ::tolower);

			if (pref_dev.compare("gpu") == 0) pref_device_type = CL_DEVICE_TYPE_GPU;
			else if (pref_dev.compare("cpu") == 0) pref_device_type = CL_DEVICE_TYPE_CPU;
			else if (pref_dev.compare("apu") == 0) pref_device_type = CL_DEVICE_TYPE_ACCELERATOR;
		}
		
		set<best_cmd_queue_info*> allowedQueueList;
		
		if (use_only_devices.size() > 0) {
			for (int i = 0; i < use_only_devices.size(); i++) {
				
				string allowed_dev = use_only_devices[i];
				std::transform(allowed_dev.begin(), allowed_dev.end(), allowed_dev.begin(), ::tolower);

				cl_device_type allowed_device_type = CL_DEVICE_TYPE_ALL;

				if (allowed_dev.compare("gpu") == 0) allowed_device_type = CL_DEVICE_TYPE_GPU;
				else if (allowed_dev.compare("cpu") == 0) allowed_device_type = CL_DEVICE_TYPE_CPU;
				else if (allowed_dev.compare("apu") == 0) allowed_device_type = CL_DEVICE_TYPE_ACCELERATOR;

				for (int j = 0; j < commandQueueList.size(); j++) {
					if (allowed_device_type == CL_DEVICE_TYPE_ALL || commandQueueList[i].device_type == allowed_device_type) {
						allowedQueueList.insert(&commandQueueList[i]);
					}
				}
			}
		} else {
			for (int i = 0; i < commandQueueList.size(); i++) {
				allowedQueueList.insert(&commandQueueList[i]);
			}
		}

		best_cmd_queue_info& best_queue = *(*(allowedQueueList.begin()));

		if (pref_device_type != CL_DEVICE_TYPE_ALL && allowedQueueList.size() > 1) {
			// try to find prefered device
			for (set<best_cmd_queue_info*>::iterator it = allowedQueueList.begin(); it != allowedQueueList.end(); it++) {
				if ((*it)->device_type == pref_device_type) {
					best_queue = *(*it);
				}
			}
		}

		return best_queue;
	}
}

#endif // opencl endif

// define static context manager in different class so access to it can be made with more
// normal name 
class OpenCL {
public:
#ifdef OPENCL
	static clContextManager* context_manager;
#endif
	static void initializeContexts(vector<string> &kernels_path, vector<string>& allowed_devices);
};



#endif // file declaration endif