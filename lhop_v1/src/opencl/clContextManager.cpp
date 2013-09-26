
#include "clContextManager.h"

#include "utils/exceptions.h"
#include "utils/utils.h"

#include <vector>
#include <string>

using namespace std;

#ifdef OPENCL

#include "cl_utils.h"


#include <math.h>

#include <sstream>

clContextManager* OpenCL::context_manager = new clContextManager();

#endif

void OpenCL::initializeContexts(vector<string> &kernels_path, vector<string>& allowed_devices) {
#ifdef OPENCL
	context_manager->initializeContexts(kernels_path, allowed_devices);
#endif
}

#ifdef OPENCL
void CL_CALLBACK opencl_pfn_report_error(const char *errinfo, const void *private_info, size_t cb, void *user_data) {
	cout << "opencl error: "  << errinfo << endl;
}

clContextManager::clContextManager(void):
	counter(0), is_initialized(false)
{
}

clContextManager::~clContextManager(void)
{
	cl_int ret = 0;

	// release all kernels
	for (map<string,vector<cl_kernel> >::iterator iter = kernelList.begin(); 
			iter != kernelList.end(); iter++) {
		vector<cl_kernel> all_kernels = (*iter).second;
		for (int i = 0; i < all_kernels.size(); i++) {
			ret = clReleaseKernel(all_kernels[i]);
			ocl_check_error(ret, CL_SUCCESS);
		}
	}

	// release all queues
	for (int i = 0; i < commandQueueList.size(); i++) {
		ret = clReleaseCommandQueue(commandQueueList[i].queue);
		ocl_check_error(ret, CL_SUCCESS);
	}

	// then release all program objects
	for (list<cl_program>::iterator iter = programList.begin(); 
			iter != programList.end(); iter++) {
		ret = clReleaseProgram(*iter);
		ocl_check_error(ret, CL_SUCCESS);
	}

	// and finaly release each context
	for (int i = 0; i < contextList.size(); i++) {
		ret = clReleaseContext(contextList[i].c);
		ocl_check_error(ret, CL_SUCCESS);
	}
}

void clContextManager::initializeContexts(vector<string> &kernels_path, vector<string>& allowed_devices, bool force_reinitialize)
{
	if (is_initialized == true && force_reinitialize == false)
		return;

	const cl_device_type default_device = CL_DEVICE_TYPE_GPU; // CL_DEVICE_TYPE_DEFAULT;

	cl_int ret = 0;

	// find all platforms
	cl_uint platform_count = 0;	
	ret = clGetPlatformIDs(0, nullptr, &platform_count);
	ocl_check_error(ret, CL_SUCCESS);

	cl_platform_id* platform_ids = new cl_platform_id[platform_count];

	ret = clGetPlatformIDs(platform_count, platform_ids, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	if (allowed_devices.size() <= 0) {
		allowed_devices.push_back("default");
	}

	// set kernel paths
	stringstream kernels_include_path;

	for (vector<string>::const_iterator it = kernels_path.begin(); it != kernels_path.end(); it++) {
		kernel_paths.push_back((*it));
		kernels_include_path << "-I " << (*it) << " ";
	}

	// add default folder to end
	string current_path(get_current_path());
	kernel_paths.push_back(current_path);
	kernels_include_path << "-I " << current_path << " ";

	opencl_include_paths = kernels_include_path.str();

	for (vector<string>::iterator it = allowed_devices.begin(); it != allowed_devices.end(); it++) {

		cl_device_type only_device_type = CL_DEVICE_TYPE_ALL;

		std::transform((*it).begin(), (*it).end(), (*it).begin(), ::tolower);

		// split (*it) by '::', left side defines type (gpu, cpu ..) while right side can specify platform name (different drivers)
		size_t sep_position = (*it).find("::");
		
		string device_type = (*it).substr(0, sep_position);
		string platform_filter = (sep_position != string::npos ? (*it).substr(sep_position + 2) : string());

		if (device_type.compare("gpu") == 0) only_device_type = CL_DEVICE_TYPE_GPU;
		else if (device_type.compare("cpu") == 0) only_device_type = CL_DEVICE_TYPE_CPU;
		else if (device_type.compare("apu") == 0) only_device_type = CL_DEVICE_TYPE_ACCELERATOR;
		else if (device_type.compare("default") == 0) only_device_type = default_device;
		else if (device_type.length() > 0) {
			throw custom_libhop_exception(opencl_exception, "Error: Invalid OpenCL device type, allowed only CPU, GPU, APU, default or empty with platform filter set (e.g. '::Intel')");
		}

		// find all devices for each platform
		for (int i = 0; i < platform_count; i++) {

			// if platform_filter is set then filter platforms that do not contain this name
			if (platform_filter.length() > 0) {
				size_t str_name_size = 0;
				ret = clGetPlatformInfo(platform_ids[i], CL_PLATFORM_NAME, 0, nullptr, &str_name_size);
				ocl_check_error(ret, CL_SUCCESS);

				char* c_platform_name = new char[str_name_size];
				ret = clGetPlatformInfo(platform_ids[i], CL_PLATFORM_NAME, str_name_size, c_platform_name, nullptr);
				ocl_check_error(ret, CL_SUCCESS);

				string platform_name(c_platform_name);

				std::transform(platform_name.begin(), platform_name.end(), platform_name.begin(), ::tolower);
				std::transform(platform_filter.begin(), platform_filter.end(), platform_filter.begin(), ::tolower);
				
				bool skip = platform_name.find(platform_filter) == string::npos;

				delete[] c_platform_name;

				if (skip)
					continue;
			}
			// find only devices for device type specified by argument @only_device_type
			cl_uint devices_count = 0;
			ret = clGetDeviceIDs(platform_ids[i], only_device_type, 0, nullptr, &devices_count);
			
			// skip to next platform if no suitable devices found - must check this before checking for error since error may also be CL_DEVICE_NOT_FOUND
			if (devices_count <= 0 && ret == CL_DEVICE_NOT_FOUND)
				continue;

			// now check for errors
			ocl_check_error(ret, CL_SUCCESS);

			// make room for all device ids
			cl_device_id* device_ids = new cl_device_id[devices_count];

			// get device ids
			ret = clGetDeviceIDs(platform_ids[i], only_device_type, devices_count, device_ids, nullptr);
			ocl_check_error(ret, CL_SUCCESS);

			// create context for this platform with all devices in it
			cl_context_properties propreties[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)platform_ids[i], 0};

			int device_count_per_context = seperate_contex_for_devices == true ? 1 : devices_count;
			int context_count = seperate_contex_for_devices == true ? devices_count : 1;

			for (int dev = 0; dev < context_count; dev++) {
				cl_device_id* context_device_ids = device_ids + dev;
				t_context context;
				context.c = clCreateContext(propreties, device_count_per_context, context_device_ids, opencl_pfn_report_error, nullptr, &ret);
				ocl_check_error(ret, CL_SUCCESS);

				// retain context in case if it might be deleted by user
				ret = clRetainContext(context.c);
				ocl_check_error(ret, CL_SUCCESS);

				// create context structure
				context.context_number = counter++;

				// just push context into vector list

				if (this->contextList.size() <= context.context_number)
					this->contextList.resize(context.context_number + 1); 
				this->contextList[context.context_number] = context;

				// also create default command queue for each device
				for (int j = 0; j < device_count_per_context; j++) {
					cl_device_id d_id = context_device_ids[j];
					cl_command_queue queue = clCreateCommandQueue(context.c, context_device_ids[j], CL_QUEUE_PROFILING_ENABLE, &ret);
					ocl_check_error(ret, CL_SUCCESS);
					
					best_cmd_queue_info cmd_queue;
					cmd_queue.context_number = context.context_number;
					cmd_queue.device_id = d_id;
					cmd_queue.queue = queue;

					// also assing some info about device (max work group size, has img support etc)
					ret = clGetDeviceInfo(d_id, CL_DEVICE_IMAGE_SUPPORT, sizeof(cl_bool), &cmd_queue.has_img_support, nullptr);
					ocl_check_error(ret, CL_SUCCESS);
					
					ret = clGetDeviceInfo(d_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &cmd_queue.max_work_group_size, nullptr);
					ocl_check_error(ret, CL_SUCCESS);

					ret = clGetDeviceInfo(d_id, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(cl_device_local_mem_type), &cmd_queue.local_mem_type, nullptr);
					ocl_check_error(ret, CL_SUCCESS);

					ret = clGetDeviceInfo(d_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &cmd_queue.local_mem_size, nullptr);
					ocl_check_error(ret, CL_SUCCESS);

					ret = clGetDeviceInfo(d_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_int), &cmd_queue.max_compute_units, nullptr);
					ocl_check_error(ret, CL_SUCCESS);

					ret = clGetDeviceInfo(d_id, CL_DEVICE_TYPE, sizeof(cl_device_type), &cmd_queue.device_type, nullptr);
					ocl_check_error(ret, CL_SUCCESS);

					ret = clGetDeviceInfo(d_id, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(unsigned long long), &cmd_queue.max_mem_alocation, nullptr);
					ocl_check_error(ret, CL_SUCCESS);
					

					calculateBestWorkGroupSizes(&cmd_queue);

					this->commandQueueList.push_back(cmd_queue);
				}

				cout << "Building OpenCL kernels ... " ;
				// create programs for this context
				buildProgramsAndKernels(context, context_device_ids, device_count_per_context);
				cout << " (done)" << endl;
			}
			delete[] device_ids;
		}
	}
	delete[] platform_ids;

	is_initialized = true;
}

void clContextManager::buildProgramsAndKernels(t_context context, cl_device_id* device_list, int device_list_size, bool by_binary) {
	// for each program :
	for (int i = 0; i < opencl_program_listing_size; i++) {
		// get information where program can be loaded
		//const ocl_program_definition& prog_def = opencl_program_listing[i];
		const string name = opencl_program_listing[i];
		//const char* name = opencl_program_listing[i].c_str();
		// find file, open it and load its content into buffer
		size_t progLength;
		char * programSource = nullptr;
		
		vector<string>::const_iterator kernel_paths_iter = kernel_paths.begin();
		while (programSource == nullptr && kernel_paths_iter != kernel_paths.end()) {
			programSource = ocl_load_prog_source( ((*kernel_paths_iter) + "/" + name).c_str(), "", &progLength);
			kernel_paths_iter++;
		}

		if (programSource == nullptr) {
			stringstream exception_msg;
			exception_msg << "Unable to find OpenCL kernel src file '" << name << "'." << endl;
			exception_msg << "Search paths used: " << endl;
			for (vector<string>::const_iterator it = kernel_paths.begin(); it != kernel_paths.end(); it++) {
				exception_msg << "\t" << ((*it)) << endl;;				
			}
			throw custom_libhop_exception(opencl_exception, exception_msg.str());
			//cout << "Stoped on error. Press any key to exit ... " << endl;
			//getchar();
			//exit(0);
		}

		cl_program prog = nullptr;
		if (by_binary == false) {
			// load it from source and build it			
			prog = compileWithSource(context, device_list, device_list_size, programSource, progLength, name.c_str());
		} else {
			// or just load it from binary if @by_binary is true 
			prog = compileWithBinary(context, device_list, device_list_size, programSource, progLength, name.c_str());
		}

		// program source can be deleted after we create program
		delete[] programSource;

		cl_uint ret;
		// now continue with loading of kernels
		
		// load all kernels in program at once
		cl_uint num_kernels = 0;
		ret = clCreateKernelsInProgram(prog, 0, nullptr, &num_kernels);
		ocl_check_error(ret, CL_SUCCESS);

		// create memory for all kernels
		cl_kernel* prog_kernels = new cl_kernel[num_kernels];
		ret = clCreateKernelsInProgram(prog, num_kernels, prog_kernels, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		// now for each kernel find its name
		// and save it to its mapping 
		for (int j = 0; j < num_kernels; j++) {
			// get kernel name 
			string kernel_name = getKernelName(prog_kernels[j]);

			// add it to mapping
			vector<cl_kernel>& kernel_list = this->kernelList[kernel_name];

			// always resize list, so one new element can be added
			kernel_list.resize(kernel_list.size() + 1);

			// assign kernel to correct position
			kernel_list.at(context.context_number) = prog_kernels[j];
		}

		// we can delete array of kernels, since they were copied to kernel_list
		delete[] prog_kernels;
	}
}

cl_program clContextManager::compileWithSource(const t_context context, const cl_device_id* device_list, const int device_list_size, 
											   const char * programSource, const size_t progLength, const char* program_name) {

	// first load source
	cl_int ret = 0;
	cl_program prog = clCreateProgramWithSource(context.c, 1, (const char **)&programSource, &progLength, &ret);
	ocl_check_error(ret, CL_SUCCESS);

	stringstream compileOptions;
	compileOptions << " -cl-fast-relaxed-math " << opencl_include_paths << " ";

	// if device_list_size is one then also add some device capabilities to compile options
	if (device_list_size == 1) {
		compileOptions << createCompileOptFromDeviceCapabilities(device_list[0]);
	}

	ret = clBuildProgram(prog, 0, nullptr, compileOptions.str().c_str(), nullptr, nullptr);
	
/*
	size_t* binary_sizes = new size_t[device_list_size];
	ret = clGetProgramInfo(prog, CL_PROGRAM_BINARY_SIZES, sizeof(binary_sizes), binary_sizes, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	unsigned char ** prog_binary = new unsigned char *[device_list_size];

	for (int j = 0; j < device_list_size; j++) {

		prog_binary[j] = new unsigned char[binary_sizes[j]];
	}

	ret = clGetProgramInfo(prog, CL_PROGRAM_BINARIES, binary_sizes[0], prog_binary, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	string bin_file;
	bin_file.append(program_name);
	bin_file.append(".bin");
	FILE* fp = fopen(bin_file.c_str(), "wb");

	fwrite(prog_binary[0], sizeof(unsigned char), binary_sizes[0], fp);
	fflush(fp);
	fclose(fp);

	*/
	// check for build errors
	if (ret == CL_BUILD_PROGRAM_FAILURE) {
		// error with building - report to user
		cout << "Error while building OpenCL program '" << program_name << "':" << endl;
		// get build log for all devices
		for (int j = 0; j < device_list_size; j++) {

			// first get size of memory needed for log
			size_t param_value_size_ret;
			ret = clGetProgramBuildInfo(prog, device_list[j], CL_PROGRAM_BUILD_LOG, 0, nullptr, &param_value_size_ret);
			ocl_check_error(ret, CL_SUCCESS);

			// then reserve memory and get log
			char* debug_info = new char[param_value_size_ret];
			ret = clGetProgramBuildInfo(prog, device_list[j], CL_PROGRAM_BUILD_LOG, param_value_size_ret, debug_info, nullptr);
			ocl_check_error(ret, CL_SUCCESS);

			// finally display log
			cout << "Build log for device '" << getDeviceName(device_list[j]) << "' : " << endl; 
			cout << debug_info << endl;

			delete[] debug_info;
		}
		
		throw custom_libhop_exception(opencl_exception, "Error while building/compiling OpenCL program.");
		/*cout << "Error while building/compiling OpenCL program. Do you want to continue with building next program (yes/no)?" << endl;
		string resp;
		cin >> resp;

		//ret = clReleaseProgram(prog);
		ocl_check_error(ret, CL_SUCCESS);


		
		std::transform(resp.begin(), resp.end(), resp.begin(), ::tolower);
		if (resp.compare("y") != 0 && resp.compare("yes") != 0)
			exit(0);
			*/
	} else {
		// if nothing wrong with syntax, then also check for other errors
		ocl_check_error(ret, CL_SUCCESS);
	}
	programList.push_back(prog);

	return prog;
}
cl_program clContextManager::compileWithBinary(const t_context context, const cl_device_id* device_list, const int device_list_size, 
											   const char * programSource, const size_t progLength, const char* program_name) {
	throw custom_libhop_exception(opencl_exception, "Compiling with binary is not yet supported by this app.");
}

string clContextManager::getKernelName(const cl_kernel kernel) {
	cl_int ret = 0;
	// get display name size
	size_t param_value_size_ret;
	ret = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, 0, &param_value_size_ret);
	ocl_check_error(ret, CL_SUCCESS);

	// allocate memory and get name
	char* name = new char[param_value_size_ret];
	ret = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, param_value_size_ret, name, 0);
	ocl_check_error(ret, CL_SUCCESS);

	string result(name);

	delete[] name;
	return result;
}

string clContextManager::getDeviceName(const cl_device_id device) {

	cl_int ret = 0;
	// get display name size
	size_t param_value_size_ret;
	ret = clGetDeviceInfo(device, CL_DEVICE_NAME, 0, 0, &param_value_size_ret);
	ocl_check_error(ret, CL_SUCCESS);

	// allocate memory and get name
	char* name = new char[param_value_size_ret];
	ret = clGetDeviceInfo(device, CL_DEVICE_NAME, param_value_size_ret, name, 0);
	ocl_check_error(ret, CL_SUCCESS);

	string result(name);

	delete[] name;
	return result;
}


string clContextManager::createCompileOptFromDeviceCapabilities(const cl_device_id device) {
	cl_int ret;
	
	stringstream result;

	cl_device_type device_type;
	ret = clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(cl_device_type), &device_type, nullptr);
	ocl_check_error(ret, CL_SUCCESS);
	
	switch (device_type) {
		case CL_DEVICE_TYPE_CPU:
			result << " -D __DEVICE_CPU__ ";
			break;
		case CL_DEVICE_TYPE_GPU:
			result << " -D __DEVICE_GPU__ ";
			break;
		case CL_DEVICE_TYPE_ACCELERATOR:
			result << " -D __DEVICE_APU__ ";
			break;
		case CL_DEVICE_TYPE_DEFAULT:			
		default:
			result << " -D __DEVICE_DEFAULT_TYPE__ ";
			break;
	}

	cl_bool supports_image = false;
	ret = clGetDeviceInfo(device, CL_DEVICE_IMAGE_SUPPORT, sizeof(cl_bool), &supports_image, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	if (supports_image == true && device_type != CL_DEVICE_TYPE_CPU)
		result << " -D __DEVICE_IMAGE_SUPPORTED__ ";

	cl_device_local_mem_type local_mem_type;
	ret = clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(cl_device_local_mem_type), &local_mem_type, nullptr);
	ocl_check_error(ret, CL_SUCCESS);
	
	//switch (CL_LOCAL) {
	//switch (CL_GLOBAL) {
	switch (local_mem_type) {
		case CL_GLOBAL:
			result << " -D __DEVICE_LOCAL_MEM_TYPE__=0 ";
			break;
		case CL_LOCAL:
			result << " -D __DEVICE_LOCAL_MEM_TYPE__=1 ";
			break;
	}

	cl_uint max_compute_unit;
	ret = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_unit, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	result << " -D __DEVICE_MAX_COMPUTE_UNITS__=" << max_compute_unit << " " ;

	return result.str();
}

void clContextManager::calculateBestWorkGroupSizes(best_cmd_queue_info* cmd_queue) {

	// calculate max and best work group size 
	
	// first for dimension 2
	{
		// try using values that are power of two
		frexp(::sqrt((float)cmd_queue->max_work_group_size), &cmd_queue->best_work_groups_dim2.pow_x);
		cmd_queue->best_work_groups_dim2.pow_x--;

		//cmd_queue->best_work_groups_dim2.pow_x = ::log(::sqrt((float)cmd_queue->max_work_group_size))/::log(2.0f);
		cmd_queue->best_work_groups_dim2.x = ::pow(2.0, cmd_queue->best_work_groups_dim2.pow_x);

		// based on x get best y 
		cmd_queue->best_work_groups_dim2.y = cmd_queue->max_work_group_size / cmd_queue->best_work_groups_dim2.x;
		//cmd_queue->best_work_groups_dim2.pow_y = ::log((float)cmd_queue->best_work_groups_dim2.y)/::log(2.0f);
		frexp((float)cmd_queue->best_work_groups_dim2.y, &cmd_queue->best_work_groups_dim2.pow_y);
		cmd_queue->best_work_groups_dim2.pow_y--;
	}

	// then for dimension 3 
	{
		// try using values that are power of two
		//cmd_queue->best_work_groups_dim3.pow_x = ::log( ::pow((float)cmd_queue->max_work_group_size, 1/(float)3) ) /::log(2.0f);
		frexp(::pow((float)cmd_queue->max_work_group_size, 1/(float)3), &cmd_queue->best_work_groups_dim3.pow_x);
		cmd_queue->best_work_groups_dim3.pow_x--;
		cmd_queue->best_work_groups_dim3.x = ::pow(2.0, cmd_queue->best_work_groups_dim3.pow_x);

		// based on x get best y 
		float best_yz = (float)cmd_queue->max_work_group_size / (float)cmd_queue->best_work_groups_dim3.x;
		//cmd_queue->best_work_groups_dim3.pow_y = ::log( ::sqrt(best_yz) ) /::log(2.0f);
		frexp(::sqrt(best_yz), &cmd_queue->best_work_groups_dim3.pow_y);
		cmd_queue->best_work_groups_dim3.pow_y--;
		cmd_queue->best_work_groups_dim3.y = ::pow(2.0, cmd_queue->best_work_groups_dim3.pow_y);

		// based on y get best z
		cmd_queue->best_work_groups_dim3.z = best_yz / cmd_queue->best_work_groups_dim3.y;
		//cmd_queue->best_work_groups_dim3.pow_z = ::log((float)cmd_queue->best_work_groups_dim3.z)/::log(2.0f);
		frexp((float)cmd_queue->best_work_groups_dim3.z, &cmd_queue->best_work_groups_dim3.pow_z);
		cmd_queue->best_work_groups_dim3.pow_z--;
	}
}

#endif // opencl endif
