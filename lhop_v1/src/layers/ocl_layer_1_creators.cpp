
#ifdef OPENCL

// ocl_layer_1_creators
///////////////////////////////////////////////////////////////////////////////
#include "../utils/atom.h"
#include "../graphs/img_graph.h"
#include "../utils/misc.h"
#include "../utils/utils.h"
#include "layer_1_creators.h"

#include "../opencl/cl_utils.h"
#include "../opencl/clContextManager.h"

#include <sstream>

//#define OCL_PROFILE_KERNELS_HOST // profile opencl kernels using HOST CStopWatch and clFinish()
//#define OCL_PROFILE_KERNELS_OPENCL // profile opencl kernels using opencl event profiling
//#define OCL_PROFILE_HOST_COMPLETE // profile each HOST code individually
//#define OCL_PROFILE_HOST_QUICK // profile whole HOST code as one (including waiting on opencl kernels)

int layer1_creator::ocl_create_result(layer1_result* result, const img* im, const img* maskimg) {

	// resize and blur image if neccesary
	int ret = 0;

	// use GPU as prefered type
	const best_cmd_queue_info& cmd_queue_info = OpenCL::context_manager->getBestCommandQueue("gpu", use_opencl_devices);

	cl_context context = OpenCL::context_manager->getContext(cmd_queue_info.context_number);

	cl_command_queue cmd_queue = cmd_queue_info.queue;
	cl_device_id device_id = cmd_queue_info.device_id;

	int context_number = cmd_queue_info.context_number;

	cl_event evnt_markers_start[10];
	cl_event evnt_markers_stop[10];



#if defined OCL_PROFILE_HOST_QUICK || defined OCL_PROFILE_HOST_COMPLETE
clock_3.startTimer();
#endif

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_1.startTimer();
#endif

	result->ocl_shape_nodes.clear();
	result->ocl_edges.clear();

	result->ocl_shape_nodes_coord.clear();
	result->ocl_shape_nodes_inhib_coord.clear();

	result->ocl_shape_nodes_coord_non_zero_count.clear();
	result->ocl_shape_nodes_inhib_coord_non_zero_count.clear();

	result->info.clear();

	result->ocl_shape_nodes.push_back(make_pair<ocl_layer1_data*, int>(nullptr,0));
	result->ocl_edges.push_back(make_pair<ocl_edge_data_ip2*, int>(nullptr,0));
	
	result->ocl_shape_nodes_coord.push_back(make_pair<ocl_layer1_data_coordinates*, int>(nullptr,0));
	result->ocl_shape_nodes_inhib_coord.push_back(make_pair<ocl_layer1_data_coordinates*, int>(nullptr,0));
	
	result->ocl_shape_nodes_coord_non_zero_count.push_back(0);
	result->ocl_shape_nodes_inhib_coord_non_zero_count.push_back(0);

	result->info.push_back(layer_info());

	cl_ulong p_queued, p_submit, p_start, p_end;	

	int memory_copy_type = CL_MEM_COPY_HOST_PTR;

	{

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_2.startTimer();
#endif

		cl_int2 image_size, new_image_size;

		image_size.s[0] = im->width;
		image_size.s[1] = im->height;

		new_image_size.s[0] = result->x_size(0);
		new_image_size.s[1] = result->y_size(0);

		int new_image_full_size = new_image_size.s[0] * new_image_size.s[1];

		int kernel_err_host = 0;
		cl_mem kernel_err = clCreateBuffer(context, CL_MEM_READ_WRITE | memory_copy_type, sizeof(cl_int), &kernel_err_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem made_images = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * ceil((float)(im->size() * make_images_get_count())/16)*16, nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);
	
		cl_mem max_image = clCreateBuffer(context, CL_MEM_READ_WRITE , sizeof(cl_float) * ceil((float)(im->size())/16)*16, nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem max_image_index = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * ceil((float)(im->size())/16)*16, nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_event evnts_set_zero[3];
		ocl_set_memory_zero(cmd_queue_info, made_images, im->size() * make_images_get_count() / 16,0,  0, nullptr, &evnts_set_zero[0]);
		ocl_set_memory_zero(cmd_queue_info, max_image, im->size() / 16, 0, 0, nullptr, &evnts_set_zero[1]);
		ocl_set_memory_zero(cmd_queue_info, max_image_index, im->size() / 16, 0, 0, nullptr, &evnts_set_zero[2]);

		cl_mem sum_2x2 = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) * im->size(), nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem maximum_max_image = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) , nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem maximum_sum_2x2 = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(cl_float) , nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_2.stopTimer();
clock_2.sumElapsedTime(); clock_2.sumElapsedTime(); cout << " buffers part time : " << clock_2.getElapsedTime() << endl;
#endif

		cl_event evnt_convolve;
		ocl_make_images(cmd_queue_info, result, im, made_images, max_image, max_image_index, kernel_err, memory_copy_type, 3, evnts_set_zero, &evnt_convolve);

		cl_event evnt_sum_2x2, evnt_max_sum_2x2, evnt_max_of_max_img;
		ocl_make_sum2x2(cmd_queue_info, im, max_image, sum_2x2, kernel_err, memory_copy_type, 1, &evnt_convolve, &evnt_sum_2x2);
		ocl_get_maximus(cmd_queue_info, im, max_image, maximum_max_image, kernel_err, memory_copy_type, 1, &evnt_convolve, &evnt_max_of_max_img);
		ocl_get_maximus(cmd_queue_info, im, sum_2x2, maximum_sum_2x2, kernel_err, memory_copy_type, 1, &evnt_sum_2x2, &evnt_max_sum_2x2);

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_2.startTimer();
#endif
		
		// get memory needed for coordinates of first layer
		int adjusted_layer1_coords_size = new_image_full_size;

		// align size to 4 
		if (adjusted_layer1_coords_size % 4 != 0)
			adjusted_layer1_coords_size +=  (4 - (adjusted_layer1_coords_size % 4));

		//adjusted_layer1_coords_size = ::pow(2,::ceil(::log((float)adjusted_layer1_coords_size)/::log((float)2)));

		cl_mem layer1_coords_compact_sort_keys = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_coord_sorting_pointer) * adjusted_layer1_coords_size, nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem layer1_coords_compact = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_layer1_data_coordinates) * adjusted_layer1_coords_size, nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		/*{
			cl_int in_all_bytes = (adjusted_layer1_coords_size * sizeof(ocl_layer1_data_coordinates));
			cl_int in_unalinged_bytes = in_all_bytes - (in_all_bytes / sizeof(cl_int16)) * sizeof(cl_int16);

			ocl_set_memory_zero(cmd_queue_info, layer1_coords_compact, in_all_bytes / sizeof(cl_int16), in_unalinged_bytes / sizeof(cl_int), 0, nullptr, &evnts_set_zero[0]);
		}*/
		{
			cl_int in_all_bytes = (adjusted_layer1_coords_size * sizeof(ocl_coord_sorting_pointer));
			cl_int in_unalinged_bytes = in_all_bytes - (in_all_bytes / sizeof(cl_int16)) * sizeof(cl_int16);

			ocl_set_memory_zero(cmd_queue_info, layer1_coords_compact_sort_keys, in_all_bytes / sizeof(cl_int16), in_unalinged_bytes / sizeof(cl_int), 0, nullptr, &evnts_set_zero[0]);
		}
		

		int number_created_results_host = 0;
		cl_mem number_created_result = clCreateBuffer(context, CL_MEM_READ_WRITE | memory_copy_type, sizeof(cl_int) , &number_created_results_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		int number_positions_used_result_host = 0;
		cl_mem number_positions_used_result = clCreateBuffer(context, CL_MEM_READ_WRITE | memory_copy_type, sizeof(cl_int) , &number_positions_used_result_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_2.stopTimer();
clock_2.sumElapsedTime(); clock_2.sumElapsedTime(); cout << " buffers part time : " << clock_2.getElapsedTime() << endl;
#endif

		cl_event evnt_prepare_results;	
		
		cl_event waiting_events[3] = {evnt_max_of_max_img, evnt_max_sum_2x2, evnts_set_zero[0]};
		ocl_prepare_results(cmd_queue_info, result, im, maskimg,
							made_images, max_image, max_image_index, sum_2x2, maximum_max_image, maximum_sum_2x2,
							layer1_coords_compact, layer1_coords_compact_sort_keys, number_created_result, number_positions_used_result, kernel_err, 
							memory_copy_type, 3, waiting_events, &evnt_prepare_results);

#ifdef OCL_PROFILE_KERNELS_OPENCL
ret = clEnqueueMarker(cmd_queue, &evnt_markers_start[1]);
ocl_check_error(ret, CL_SUCCESS);
#endif

		// number_created_result => is number of all elements that will need to be reserved for layer1_result_data 
		ret = clEnqueueReadBuffer(cmd_queue, number_created_result, CL_FALSE, 0, sizeof(cl_int), &number_created_results_host, 1, &evnt_prepare_results, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		// number_positions_used_result => is number of elements in layer1_coords_compact array (i.e. number of different positions used)
		ret = clEnqueueReadBuffer(cmd_queue, number_positions_used_result, CL_FALSE, 0, sizeof(cl_int), &number_positions_used_result_host, 1, &evnt_prepare_results, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clEnqueueReadBuffer(cmd_queue, kernel_err, CL_FALSE, 0, sizeof(cl_int), &kernel_err_host, 1, &evnt_prepare_results, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

#ifdef OCL_PROFILE_KERNELS_OPENCL
ret = clEnqueueMarker(cmd_queue, &evnt_markers_stop[1]);
ocl_check_error(ret, CL_SUCCESS);
#endif
		ret = clFinish(cmd_queue);
		ocl_check_error(ret, CL_SUCCESS);

		clReleaseMemObject(max_image_index);
		clReleaseMemObject(sum_2x2);
		clReleaseMemObject(maximum_max_image);
		clReleaseMemObject(maximum_sum_2x2);

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_2.startTimer();
#endif

		ocl_layer1_data* layer1_result_data_host = nullptr;
		ocl_layer1_data_coordinates* layer1_coords_compact_host = nullptr;
		ocl_layer1_data_coordinates* layer1_coords_inhib_host = new ocl_layer1_data_coordinates[new_image_full_size];

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_2.stopTimer();
clock_2.sumElapsedTime(); clock_2.sumElapsedTime(); cout << " buffers part time : " << clock_2.getElapsedTime() << endl;
#endif
		cl_event envt_sort_coords, envt_set_offsets;
		cl_event evnt_make_results, evnt_sort_results; 

		if (number_created_results_host > 0) {

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_2.startTimer();
#endif

			layer1_coords_compact_host = new ocl_layer1_data_coordinates[number_positions_used_result_host];

			cl_mem layer1_result_data = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_layer1_data) * number_created_results_host, nullptr, &ret);
			ocl_check_error(ret, CL_SUCCESS);

			cl_mem merge_sort_temp_data = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_layer1_data) * number_created_results_host, nullptr, &ret);
			ocl_check_error(ret, CL_SUCCESS);

			cl_mem layer1_coords_inhib = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(ocl_layer1_data_coordinates) * adjusted_layer1_coords_size, nullptr, &ret);
			ocl_check_error(ret, CL_SUCCESS);

			cl_event before_sort_evnt[2];
			cl_uint num_events_to_wait_before_sort = 0;
			// in case adjusted_layer1_coords_size is not big enough for bitonic sort of array of length number_positions_used_result_host 
			// then we need to create new buffer and copy old values to it
			int bitonic_sort_buffer_size = ::pow(2,::ceil(::log((float)number_positions_used_result_host)/::log((float)2)));
			if (bitonic_sort_buffer_size > adjusted_layer1_coords_size) {
				// not big enough buffer; we need to create new buffers
				cl_mem new_keys = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_coord_sorting_pointer) * bitonic_sort_buffer_size, nullptr, &ret);
				ocl_check_error(ret, CL_SUCCESS);

				cl_mem new_compact_coords = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_layer1_data_coordinates) * bitonic_sort_buffer_size, nullptr, &ret);
				ocl_check_error(ret, CL_SUCCESS);
				
				// set memory to zero (only new_keys need to be set to zero)
				cl_event new_keys_set_zero_evnt;
				cl_int in_all_bytes = (bitonic_sort_buffer_size * sizeof(ocl_coord_sorting_pointer));
				cl_int in_unalinged_bytes = in_all_bytes - (in_all_bytes / sizeof(cl_int16)) * sizeof(cl_int16);

				// using ocl_set_memory_zero for whole new memory might be unoptimized since in line below we will copy most of it from original buffer
				// but since this if case is extremly rare it should not cause serioius performance issues
				ocl_set_memory_zero(cmd_queue_info, new_keys, in_all_bytes / sizeof(cl_int16), in_unalinged_bytes / sizeof(cl_int), 0, nullptr, &new_keys_set_zero_evnt);			

				// then we need to copy values from old buffer to new one
				ret = clEnqueueCopyBuffer(cmd_queue_info.queue, layer1_coords_compact_sort_keys, new_keys, 0, 0, sizeof(ocl_coord_sorting_pointer) * number_positions_used_result_host, 1, &new_keys_set_zero_evnt, &before_sort_evnt[0]);
				ret |= clEnqueueCopyBuffer(cmd_queue_info.queue, layer1_coords_compact, new_compact_coords, 0, 0, sizeof(ocl_layer1_data_coordinates) * number_positions_used_result_host, 0, nullptr, &before_sort_evnt[1]);
				ocl_check_error(ret, CL_SUCCESS);

				// release old buffers
				clReleaseMemObject(layer1_coords_compact_sort_keys);
				clReleaseMemObject(layer1_coords_compact);

				// and reset old buffers to point to new one
				layer1_coords_compact_sort_keys = new_keys;
				layer1_coords_compact = new_compact_coords;

				num_events_to_wait_before_sort = 2;
			}
			

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_2.stopTimer();
clock_2.sumElapsedTime(); clock_2.sumElapsedTime(); cout << " buffers part time : " << clock_2.getElapsedTime() << endl;
#endif

			// first sort resultes from ocl_prepare_results (should sort pointers in layer1_coords_compact_sort_keys and their coresponding layer1_coords_compact values)
			// and also update coordinates offset based on sorted results
			int num_events_to_wait_before_update_offset = 0;
			if (number_positions_used_result_host > 1) {
				ocl_sort_compact_coords(cmd_queue_info, im, number_positions_used_result_host, layer1_coords_compact, layer1_coords_compact_sort_keys, kernel_err, memory_copy_type, num_events_to_wait_before_sort, num_events_to_wait_before_sort > 0 ? before_sort_evnt: nullptr, &envt_sort_coords);
				num_events_to_wait_before_update_offset = 1;
			}
			ocl_update_results_offsets(cmd_queue_info, im, number_positions_used_result_host, layer1_coords_compact, kernel_err, memory_copy_type, num_events_to_wait_before_update_offset, num_events_to_wait_before_update_offset > 0 ? &envt_sort_coords : nullptr, &envt_set_offsets, &evnt_markers_start[0], &evnt_markers_stop[0]);

			// finaly make results and sort within each position
			ocl_make_results(cmd_queue_info, result, im, number_positions_used_result_host, made_images, layer1_coords_compact, layer1_coords_compact_sort_keys, layer1_result_data, kernel_err, memory_copy_type, 1, &envt_set_offsets, &evnt_make_results);
			ocl_sort_results(cmd_queue_info, result, im, layer1_result_data, number_positions_used_result_host, layer1_coords_compact, merge_sort_temp_data, kernel_err, memory_copy_type, 1, &evnt_make_results, &evnt_sort_results);
	
			// also make inhibited results
			cl_event evnt_inhib_result;
			//ocl_inhibit_results(cmd_queue_info, im, layer1_result_data, layer1_coords, adjusted_layer1_coords_size, layer1_coords_inhib, kernel_err, memory_copy_type, 1, &evnt_sort_results, &evnt_inhib_result);

#ifdef OCL_PROFILE_KERNELS_OPENCL
ret = clEnqueueMarker(cmd_queue, &evnt_markers_start[2]);
ocl_check_error(ret, CL_SUCCESS);
#endif		

			layer1_result_data_host = new ocl_layer1_data[number_created_results_host];
			
			ret = clEnqueueReadBuffer(cmd_queue, layer1_result_data, CL_FALSE, 0, sizeof(ocl_layer1_data) * number_created_results_host, layer1_result_data_host, 1, &evnt_sort_results, nullptr);
			ocl_check_error(ret, CL_SUCCESS);	

			//ret = clEnqueueReadBuffer(cmd_queue, layer1_coords_inhib, CL_FALSE, 0, sizeof(ocl_layer1_data_coordinates) * new_image_full_size, layer1_coords_inhib_host, 1, &evnt_inhib_result, nullptr);
			//ocl_check_error(ret, CL_SUCCESS);

			ret = clEnqueueReadBuffer(cmd_queue, kernel_err, CL_FALSE, 0, sizeof(cl_int), &kernel_err_host, 1, &evnt_sort_results, nullptr);
			ocl_check_error(ret, CL_SUCCESS);

			ret = clEnqueueReadBuffer(cmd_queue, layer1_coords_compact, CL_FALSE, 0, sizeof(ocl_layer1_data_coordinates) * number_positions_used_result_host, layer1_coords_compact_host, 1, &envt_set_offsets, nullptr);
			ocl_check_error(ret, CL_SUCCESS);

#ifdef OCL_PROFILE_KERNELS_OPENCL
ret = clEnqueueMarker(cmd_queue, &evnt_markers_stop[2]);
ocl_check_error(ret, CL_SUCCESS);
#endif
			ret = clFinish(cmd_queue);
			ocl_check_error(ret, CL_SUCCESS);

			clReleaseMemObject(layer1_result_data);
			clReleaseMemObject(merge_sort_temp_data);
			clReleaseMemObject(layer1_coords_inhib);

		} else {
			memset(layer1_coords_inhib_host, 0, sizeof(ocl_layer1_data_coordinates) * new_image_full_size);
		}

		clReleaseMemObject(max_image);
		clReleaseMemObject(kernel_err);
		clReleaseMemObject(made_images);
		clReleaseMemObject(layer1_coords_compact);
		clReleaseMemObject(layer1_coords_compact_sort_keys);
		clReleaseMemObject(number_created_result);

		if (kernel_err_host > 0) {
			cout << "exception in ocl_layer_1_creators with error number '" << kernel_err_host << "' !!!!!!!!!!!!!!!!!!!!!" << endl;
			getchar();
			throw std::exception();
		}

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_1.startTimer();
#endif

		// we have only created compact form of layer1_coordinates so we need to copy it to full layer1_coord where each position 
		// implicitly represets one image position
		ocl_layer1_data_coordinates* layer1_coords_host = new ocl_layer1_data_coordinates[new_image_full_size];
		memset(layer1_coords_host, 0, sizeof(ocl_layer1_data_coordinates) * new_image_full_size);
		
		for (int i = 0; i < number_positions_used_result_host; i++) {
			// get position for coordinates
			int offset = layer1_coords_compact_host[i].offset;
			int pos = layer1_result_data_host[offset].y * new_image_size.s[0] + layer1_result_data_host[offset].x;

			// copy values
			layer1_coords_host[pos] = layer1_coords_compact_host[i];
		}

		// also count number of non-zero values in layer1_coords_inhib_host and layer1_coords_host (last one is only for double checking)
		int count_non_zero_position = 0;
		int count_non_zero_position_inhib = 0;

		for (int i = 0; i < new_image_full_size; i++) {			
			count_non_zero_position += layer1_coords_host[i].size > 0 ? 1 : 0;
			//count_non_zero_position_inhib += layer1_coords_inhib_host[i].size > 0 ? 1 : 0;
		}

		if (layer1_coords_compact_host != nullptr)
			delete[] layer1_coords_compact_host;

#ifdef OCL_PROFILE_HOST_COMPLETE
clock_1.stopTimer();
#endif

		// now save all parts, edges, coordinates and sizes to layer_1_result
		result->ocl_shape_nodes[0].first = layer1_result_data_host;
		result->ocl_shape_nodes[0].second = number_created_results_host;

		result->ocl_edges[0].first = nullptr;
		result->ocl_edges[0].second = 0;

		result->ocl_shape_nodes_coord[0].first = layer1_coords_host;
		result->ocl_shape_nodes_coord[0].second = new_image_full_size;

		result->ocl_shape_nodes_inhib_coord[0].first = layer1_coords_inhib_host;
		result->ocl_shape_nodes_inhib_coord[0].second = new_image_full_size;

		result->ocl_shape_nodes_coord_non_zero_count[0] = count_non_zero_position;
		result->ocl_shape_nodes_inhib_coord_non_zero_count[0] = count_non_zero_position_inhib;

		result->make_data_from_ocl(0, false);

#if defined OCL_PROFILE_HOST_QUICK || defined OCL_PROFILE_HOST_COMPLETE
clock_3.stopTimer();
#endif

#ifdef OCL_PROFILE_KERNELS_OPENCL
		ret = clGetEventProfilingInfo(evnt_convolve, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &p_queued, nullptr);
		ret |= clGetEventProfilingInfo(evnt_convolve, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &p_submit, nullptr);
		ret |= clGetEventProfilingInfo(evnt_convolve, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &p_start, nullptr);
		ret |= clGetEventProfilingInfo(evnt_convolve, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &p_end, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		cout << "OpenCL convolution only:" << endl;
		cout << "\tqueued : " << p_queued << endl;
		cout << "\tsumbited : " << p_submit -p_queued << endl;
		cout << "\tstarted : " << p_start -p_queued<< endl;
		cout << "\tended : " << p_end -p_queued << " (run time: " << p_end - p_start  << " )" << endl;

		cout << "time for sum_2x2:" << endl;
		ret = clGetEventProfilingInfo(evnt_sum_2x2, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &p_queued, nullptr);
		ret |= clGetEventProfilingInfo(evnt_sum_2x2, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &p_submit, nullptr);
		ret |= clGetEventProfilingInfo(evnt_sum_2x2, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &p_start, nullptr);
		ret |= clGetEventProfilingInfo(evnt_sum_2x2, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &p_end, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		cout << "\tqueued : " << p_queued << endl;
		cout << "\tsumbited : " << p_submit -p_queued << endl;
		cout << "\tstarted : " << p_start -p_queued<< endl;
		cout << "\tended : " << p_end -p_queued<< " (run time: " << p_end - p_start  << " )" << endl;		
		
		ret = clGetEventProfilingInfo(evnt_max_sum_2x2, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &p_queued, nullptr);
		ret |= clGetEventProfilingInfo(evnt_max_sum_2x2, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &p_submit, nullptr);
		ret |= clGetEventProfilingInfo(evnt_max_sum_2x2, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &p_start, nullptr);
		ret |= clGetEventProfilingInfo(evnt_max_sum_2x2, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &p_end, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		cout << "time for get_max:" << endl;
		cout << "\tqueued : " << p_queued << endl;
		cout << "\tsumbited : " << p_submit -p_queued << endl;
		cout << "\tstarted : " << p_start -p_queued<< endl;
		cout << "\tended : " << p_end -p_queued<< " (run time: " << p_end - p_start  << " )" << endl;		

		cout << "time for prepare_results:" << endl;
		ret = clGetEventProfilingInfo(evnt_prepare_results, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &p_queued, nullptr);
		ret |= clGetEventProfilingInfo(evnt_prepare_results, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &p_submit, nullptr);
		ret |= clGetEventProfilingInfo(evnt_prepare_results, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &p_start, nullptr);
		ret |= clGetEventProfilingInfo(evnt_prepare_results, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &p_end, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		cout << "\tqueued : " << p_queued << endl;
		cout << "\tsumbited : " << p_submit -p_queued << endl;
		cout << "\tstarted : " << p_start -p_queued<< endl;
		cout << "\tended : " << p_end -p_queued<< " (run time: " << p_end - p_start  << " )" << endl;		
		
		
		cout << "full time for set_offsets:" << endl;
		ret = clGetEventProfilingInfo(evnt_markers_start[0], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &p_start, nullptr);
		ret |= clGetEventProfilingInfo(evnt_markers_stop[0], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &p_end, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		cout << "\tstarted : " << p_start<< endl;
		cout << "\tended : " << p_end << " (run time: " << p_end - p_start  << " )" << endl;
	

		cout << "time for make_results:" << endl;
		ret = clGetEventProfilingInfo(evnt_make_results, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &p_queued, nullptr);
		ret |= clGetEventProfilingInfo(evnt_make_results, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &p_submit, nullptr);
		ret |= clGetEventProfilingInfo(evnt_make_results, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &p_start, nullptr);
		ret |= clGetEventProfilingInfo(evnt_make_results, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &p_end, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		cout << "\tqueued : " << p_queued << endl;
		cout << "\tsumbited : " << p_submit -p_queued << endl;
		cout << "\tstarted : " << p_start -p_queued<< endl;
		cout << "\tended : " << p_end -p_queued<< " (run time: " << p_end - p_start  << " )" << endl;		

		cout << "time for sotr_results:" << endl;
		ret = clGetEventProfilingInfo(evnt_sort_results, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &p_queued, nullptr);
		ret |= clGetEventProfilingInfo(evnt_sort_results, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &p_submit, nullptr);
		ret |= clGetEventProfilingInfo(evnt_sort_results, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &p_start, nullptr);
		ret |= clGetEventProfilingInfo(evnt_sort_results, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &p_end, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		cout << "\tqueued : " << p_queued << endl;
		cout << "\tsumbited : " << p_submit -p_queued << endl;
		cout << "\tstarted : " << p_start -p_queued<< endl;
		cout << "\tended : " << p_end -p_queued<< " (run time: " << p_end - p_start  << " )" << endl;


		/*cout << "time for inhib_result:" << endl;
		ret = clGetEventProfilingInfo(evnt_inhib_result, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &p_queued, nullptr);
		ret |= clGetEventProfilingInfo(evnt_inhib_result, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &p_submit, nullptr);
		ret |= clGetEventProfilingInfo(evnt_inhib_result, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &p_start, nullptr);
		ret |= clGetEventProfilingInfo(evnt_inhib_result, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &p_end, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		cout << "\tqueued : " << p_queued << endl;
		cout << "\tsumbited : " << p_submit -p_queued << endl;
		cout << "\tstarted : " << p_start -p_queued<< endl;
		cout << "\tended : " << p_end -p_queued<< " (run time: " << p_end - p_start  << " )" << endl;		
		*/

		cout << "time for marker 1:" << endl;
		ret = clGetEventProfilingInfo(evnt_markers_start[1], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &p_start, nullptr);
		ret |= clGetEventProfilingInfo(evnt_markers_stop[1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &p_end, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		cout << "\tstarted : " << p_start<< endl;
		cout << "\tended : " << p_end << " (run time: " << p_end - p_start  << " )" << endl;

		cout << "time for marker 2:" << endl;
		ret = clGetEventProfilingInfo(evnt_markers_start[2], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &p_start, nullptr);
		ret |= clGetEventProfilingInfo(evnt_markers_stop[2], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &p_end, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		cout << "\tstarted : " << p_start<< endl;
		cout << "\tended : " << p_end << " (run time: " << p_end - p_start  << " )" << endl;
		
#endif
		
#if defined OCL_PROFILE_HOST_QUICK || defined OCL_PROFILE_HOST_COMPLETE
		cout << "\t\tALL done in " << clock_3.getElapsedTime() << endl;
		cout << "\t\t  export to org structure took : " << clock_5.getElapsedTime() << endl;
#endif
#ifdef OCL_PROFILE_HOST_COMPLETE
		cout << "\t\tsetting kernel argumenst done in " << clock_4.getElapsedTime() << endl;
		cout << "\t\ttime create memory : " << clock_2.getSumTime() << endl;
		cout << "\t\taditional host operation on end done in " << clock_1.getElapsedTime() << endl;
#endif

#ifdef OPENCL_VERIFICATION_SUPPORT
		// in case binary was build with support for verification of opencl result then use config settings to indicate wheaterh to make verification
		if (opencl_verify_result == true) {
			if (number_created_results_host > 0) {
				layer1_result* ocl_result = result;
				result = new layer1_result();

				delete_images();
				init_result();
				result->new_grid(ocl_result->x_size(0), ocl_result->y_size(0), 0);

				img im_ = *im;
				make_images(im_);
				make_max_image();
				make_result(*maskimg);
				finalize();
				inhibit_result();				
				vector<node*>& s_nodes = result->shape_nodes[0];
				vector<node*>& ocl_s_nodes = ocl_result->shape_nodes[0];

				vector<node*>& s_nodes_inhib = result->shape_nodes_inhib[0];
				vector<node*>& ocl_s_nodes_inhib = ocl_result->shape_nodes_inhib[0];

				FILE* fs_inhib = fopen("results_org_inhib.txt","w");
				for (int i = 0; i < s_nodes_inhib.size(); i++) {
					fprintf(fs_inhib, "%d\t%d\t%f\n", ((layer1_data*)s_nodes_inhib[i]->data)->x, ((layer1_data*)s_nodes_inhib[i]->data)->y, ((layer1_data*)s_nodes_inhib[i]->data)->r(R_RESPONSE));
				}
				fclose(fs_inhib);

				std::sort(ocl_s_nodes_inhib.begin(), ocl_s_nodes_inhib.end(), layer1_data::greater1n);
				fs_inhib = fopen("results_opencl_inhib.txt","w");
				for (int i = 0; i < ocl_s_nodes_inhib.size(); i++) {
					fprintf(fs_inhib, "%d\t%d\t%f\n", ((layer1_data*)ocl_s_nodes_inhib[i]->data)->x, ((layer1_data*)ocl_s_nodes_inhib[i]->data)->y, ((layer1_data*)ocl_s_nodes_inhib[i]->data)->r(R_RESPONSE));
				}
				fclose(fs_inhib);
				int invalid_results = 0;
				node* next = nullptr;
				node* ocl_next = nullptr;

				int ocl_s_nodes_count = 0;
				int s_nodes_count = 0;

				for (int j = 0; j < result->y_size(0); j++) {
					for (int i = 0; i < result->x_size(0); i++) {

						node* org_next__ = result->node_at(i,j,0);
						node* ocl_next__ = ocl_result->node_at(i,j,0);
						node* org_next = org_next__;
						node* ocl_next = ocl_next__;

						while (ocl_next != nullptr || org_next != nullptr) {

							if (ocl_next == nullptr || org_next == nullptr) {
								invalid_results++;
								cout << "invalid number of parts in one location" << endl;
								continue;
							}

							if (org_next->neighbors.size() != ocl_next->neighbors.size()) {
								//invalid_results++;
							}

							layer1_data* org_part = (layer1_data*)(org_next->data);
							layer1_data* ocl_result_part = (layer1_data*)(ocl_next->data);
							
							if (ocl_result_part->x != org_part->x) {
								invalid_results++;
							}

							if (ocl_result_part->y != org_part->y) {
								invalid_results++;
							}
							if (ocl_result_part->z != org_part->z) {
								invalid_results++;
							}
							if (ocl_result_part->m != org_part->m) {
								cout << "diff at x: " << ocl_result_part->x << " y: " << ocl_result_part->y << " in type ocl: " << ocl_result_part->m << " and type org: " << org_part->m << endl;
								invalid_results++;
							}
							if (fabs(ocl_result_part->r(R_RESPONSE) - org_part->r(R_RESPONSE)) > 1e-4) {
								cout << "diff at x: " << ocl_result_part->x << " y: " << ocl_result_part->y << " m: " << ocl_result_part->m << "  in R_response: " << fabs(ocl_result_part->r(R_RESPONSE) - org_part->r(R_RESPONSE)) << endl;
								invalid_results++;
							}

							if (fabs(ocl_result_part->r(G_RESPONSE) - org_part->r(G_RESPONSE)) > 1e-4) {
								cout << "diff at x: " << ocl_result_part->x << " y: " << ocl_result_part->y << " m: " << ocl_result_part->m << " in G_response: " << fabs(ocl_result_part->r(G_RESPONSE) - org_part->r(G_RESPONSE)) << endl;
								invalid_results++;
							}
							if (fabs(ocl_result_part->r(RR_RESPONSE) - org_part->r(RR_RESPONSE)) > 1e-4) {
								invalid_results++;
							}
							ocl_next = ((layer1_data*)ocl_next->data)->next;
							org_next = ((layer1_data*)org_next->data)->next;

						}
					}
				}

				cout << "number invalid results : " << invalid_results << endl;

				delete result;
				result = ocl_result;
			} else {
				cout << "\tskipping validation : no parts created" << endl;
			}
		}
#endif

	}
	return 0;
}

void layer1_creator::ocl_set_offsets_recursively(const best_cmd_queue_info& cmd_queue_info, cl_mem coords, int coords_size, int additional_offsets, int work_group_size_used, 
												cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) {

	cl_context context = OpenCL::context_manager->getContext(cmd_queue_info.context_number);

	cl_command_queue cmd_queue = cmd_queue_info.queue;

	cl_kernel prefix_sum_step2 = OpenCL::context_manager->getCommonKernel("scanExclusiveLocal2")[cmd_queue_info.context_number];
	cl_kernel prefix_sum_step3 = OpenCL::context_manager->getCommonKernel("uniformUpdate")[cmd_queue_info.context_number];

	int ret = 0;

	// offset must be adjusted to size of 4
	if (additional_offsets % 4 != 0)
		additional_offsets += 4 - (additional_offsets % 4);

	cl_mem tmp_buffer = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_layer1_data_coordinates) * additional_offsets * 2, nullptr, &ret);
	ocl_check_error(ret, CL_SUCCESS);

	cl_event evnt_step2_end;

	int best_local_work_size_x = cmd_queue_info.max_work_group_size;
	{
		
		int best_group_size_x = (additional_offsets/4) / best_local_work_size_x + 1;
		int best_global_work_size_x = best_group_size_x * best_local_work_size_x;

		size_t global_work_load[1] = {best_global_work_size_x};
		size_t local_work_load[1] = {best_local_work_size_x};

		
		ret = clSetKernelArg(prefix_sum_step2, 0, sizeof(cl_mem), (void*)&coords);
		ret |= clSetKernelArg(prefix_sum_step2, 1, sizeof(cl_mem), (void*)&tmp_buffer);
		ret |= clSetKernelArg(prefix_sum_step2, 2, sizeof(cl_int) * best_local_work_size_x * 2, nullptr);
		ret |= clSetKernelArg(prefix_sum_step2, 3, sizeof(cl_int), (void*)&additional_offsets);
		ret |= clSetKernelArg(prefix_sum_step2, 4, sizeof(cl_int), (void*)&work_group_size_used);
		ocl_check_error(ret, CL_SUCCESS);
		
		cl_event evnt_step2;
		ret = clEnqueueNDRangeKernel(cmd_queue, prefix_sum_step2, 1, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, &evnt_step2);
		ocl_check_error(ret, CL_SUCCESS);

		// if resutls still bigger then group size then call recursively
		int next_additional_offsets = ceil((float)additional_offsets/(best_local_work_size_x*4)); //always round up

		if (next_additional_offsets > 1) {
			ocl_set_offsets_recursively(cmd_queue_info, tmp_buffer, additional_offsets, next_additional_offsets, best_local_work_size_x * 4, 1, &evnt_step2, &evnt_step2_end);
		} else {
			evnt_step2_end = evnt_step2;
		}

	}
	// enqueue step 3
	{
		
		int best_group_size_x = ceil((float)(coords_size/4) / best_local_work_size_x);
		int best_global_work_size_x = best_group_size_x * best_local_work_size_x;

		size_t global_work_load[1] = {best_global_work_size_x};
		size_t local_work_load[1] = {best_local_work_size_x};

		ret = clSetKernelArg(prefix_sum_step3, 0, sizeof(cl_mem), (void*)&coords);
		ret |= clSetKernelArg(prefix_sum_step3, 1, sizeof(cl_int), (void*)&coords_size);
		ret |= clSetKernelArg(prefix_sum_step3, 2, sizeof(cl_mem), (void*)&tmp_buffer);				
		ocl_check_error(ret, CL_SUCCESS);
		
		ret = clEnqueueNDRangeKernel(cmd_queue, prefix_sum_step3, 1, nullptr, global_work_load, local_work_load, 1, &evnt_step2_end, evnt);
		ocl_check_error(ret, CL_SUCCESS);
	}

	clReleaseMemObject(tmp_buffer);
}

void layer1_creator::ocl_make_images(const best_cmd_queue_info& cmd_queue_info, layer1_result* result, const img* im, cl_mem made_images, cl_mem max_image, cl_mem max_image_index, cl_mem kernel_err,
									 int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) {

	cl_context context = OpenCL::context_manager->getContext(cmd_queue_info.context_number);
	cl_command_queue cmd_queue = cmd_queue_info.queue;

	
#ifdef OCL_PROFILE_HOST_COMPLETE
	clock_1.stopTimer();
	cout << "TIME header and mapping: " << clock_1.getElapsedTime() << endl;
	clock_2.startTimer();
#endif
	// use one buffer for both input objects
	cl_float* input_output_host = new cl_float[im->size() + filter_kernels[0]->size() * filter_kernels.size()];

	int i_pos = 0;
	for (int j = 0; j < im->height; j++) {
		for (int i = 0; i < im->width; i++) {
			input_output_host[i_pos++] = im->at(i,j);
		}
	}

	cl_float* convolution_mask_host = input_output_host + im->size();
	//cl_float* convolution_mask_host = new float[masks[0]->size() * masks.size()];

	i_pos = 0;
	for (int j = 0; j < filter_kernels.size(); j++) {
		img& m = *filter_kernels[j];
		for (int i = 0; i < m.size(); i++) {				
			convolution_mask_host[i_pos++] = m[i];
		}
	}
#ifdef OCL_PROFILE_HOST_COMPLETE
	clock_2.stopTimer();
	cout << "TIME copying input data: " << clock_2.getElapsedTime() << endl;

	clock_2.startTimer();
#endif

	int ret;
	cl_mem input_image = clCreateBuffer(context, CL_MEM_READ_ONLY | memory_copy_type, sizeof(cl_float) * im->size(), input_output_host, &ret);
	ocl_check_error(ret, CL_SUCCESS);

	cl_mem convolution_mask = clCreateBuffer(context, CL_MEM_READ_ONLY | memory_copy_type, sizeof(cl_float) * filter_kernels[0]->size() * filter_kernels.size(), convolution_mask_host, &ret);
	ocl_check_error(ret, CL_SUCCESS);

#ifdef OCL_PROFILE_HOST_COMPLETE
	clock_2.stopTimer();
#endif

#ifdef OCL_PROFILE_KERNELS_HOST
	clock_4.startTimer();
#endif

	// use 16 vector for CPU to fully utilize SSE3 otherwise use vector size 4
	int vector_size = cmd_queue_info.device_type == CL_DEVICE_TYPE_GPU ? 4 : 16;

	stringstream convol_kernel_name;
	convol_kernel_name << "image_convolve_vector" << vector_size << "";

	cl_kernel convolution_kernel = OpenCL::context_manager->getCommonKernel(convol_kernel_name.str())[cmd_queue_info.context_number];
	cl_kernel convolution_img_kernel = OpenCL::context_manager->getCommonKernel("image_convolve_img")[cmd_queue_info.context_number];

	cl_int2 image_size, new_image_size, local_data_size;
	cl_int4 convol_mask_size;

	image_size.s[0] = im->width;
	image_size.s[1] = im->height;

	new_image_size.s[0] = result->x_size(0);
	new_image_size.s[1] = result->y_size(0);

	convol_mask_size.s[0] = filter_kernels[0]->width;
	convol_mask_size.s[1] = filter_kernels[0]->height;
	convol_mask_size.s[2] = filter_kernels.size();

	// set arguments and calculate work-group sizes
	//int best_local_work_size_x = 16;//cmd_queue_info.best_work_groups_dim2.x;
	//int best_local_work_size_y = 32;//cmd_queue_info.best_work_groups_dim2.y;
		
	// get memory usage and best work group size from opencl driver
	size_t max_kernel_work_group_size;
	cl_ulong kernel_local_mem_use;
	
	ret = clGetKernelWorkGroupInfo(convolution_kernel, cmd_queue_info.device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &max_kernel_work_group_size, nullptr);
	ret |= clGetKernelWorkGroupInfo(convolution_kernel, cmd_queue_info.device_id, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(cl_ulong), &kernel_local_mem_use, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	// use integer division in case if even size where convol_mask_size.s[0]-1 would not be ok
	int added_border_size = 2*(convol_mask_size.s[0]/2);

	// use only pre-calculated values by default (this are max values that we can use)
	// take into account number of compute uints that will work in parallel		
	int max_local_work_size_x_from_compute_units = im->width / (cmd_queue_info.max_compute_units * vector_size/2);

	frexp((float)max_local_work_size_x_from_compute_units, &max_local_work_size_x_from_compute_units);
	max_local_work_size_x_from_compute_units--;

	// use either default value (pre-computed max value) or max value based on number of computing units
	//int pow_x = min(cmd_queue_info.best_work_groups_dim2.pow_x, (int)::log((float)max_local_work_size_x_from_compute_units)/::log(2.0f));
	int pow_x = min<int>(cmd_queue_info.best_work_groups_dim2.pow_x, max_local_work_size_x_from_compute_units);
	int pow_y = cmd_queue_info.best_work_groups_dim2.pow_y;

	int best_local_work_size_x = ::pow(2.0, pow_x);
	int best_local_work_size_y = cmd_queue_info.best_work_groups_dim2.y;

	// if will be using local memory then calculate best work-group size
	if (cmd_queue_info.local_mem_type == CL_LOCAL) {

		// calculate how much local memory we will need
		int memory_required = 4 * (best_local_work_size_x * vector_size + added_border_size) * (best_local_work_size_y + added_border_size) + kernel_local_mem_use;
	
		// select best work-group size by decresing size until device can satisfy local memory and work group size requirements for this kernel
		while (memory_required > cmd_queue_info.local_mem_size || best_local_work_size_x * best_local_work_size_y > max_kernel_work_group_size) {
			
			if (best_local_work_size_x >= best_local_work_size_y) {
				pow_x--;
				best_local_work_size_x = ::pow(2.0, pow_x);
			} else {
				pow_y--;
				best_local_work_size_y = ::pow(2.0, pow_y);
			}

			memory_required = 4 * (best_local_work_size_x * vector_size + added_border_size) * (best_local_work_size_y + added_border_size) + kernel_local_mem_use;
		}

		local_data_size.s[0] = best_local_work_size_x * vector_size + added_border_size;
		local_data_size.s[1] = best_local_work_size_y + added_border_size;
	} else {
		// else just use default values and
		// use min local data size since local mem will not be used

		local_data_size.s[0] = 1;
		local_data_size.s[1] = 1;
	}
	
	

	//local_data_size.s[0] = 1;//best_local_work_size_x*vector_size + 2*(convol_mask_size.s[0] / 2); // use integer division in case if even size where convol_mask_size.s[0]-1 would not be ok
	//local_data_size.s[1] = 1;//best_local_work_size_y + 2*(convol_mask_size.s[1] / 2);

	int best_group_size_x = ceil((float)((im->width - 2*(convol_mask_size.s[0]/2))/vector_size) / best_local_work_size_x);
	int best_group_size_y = ceil((float)((im->height - 2*(convol_mask_size.s[1]/2))) / best_local_work_size_y);

	int best_global_work_size_x = best_group_size_x * best_local_work_size_x;
	int best_global_work_size_y = best_group_size_y * best_local_work_size_y;

	size_t global_work_load[2] = {best_global_work_size_x, best_global_work_size_y};
	size_t local_work_load[2] = {best_local_work_size_x, best_local_work_size_y};

	cl_uint save_max = 1;
	cl_uint layer1_creator_type = 1;

	switch (lib_type()) {
		case STRUCT_LIB_ATTR:
			layer1_creator_type = LAYER1_TYPE_STRUCT;
			break;
		case LOG_LIB_ATTR:
			layer1_creator_type = LAYER1_TYPE_LOGGABOR;
			break;
		case APP_LIB_ATTR:
			layer1_creator_type = LAYER1_TYPE_APP;
			break;
		case DOG_LIB_ATTR:
			layer1_creator_type = LAYER1_TYPE_DOG;
//			break;
		default:
			cout << "ERROR: OpenCL does not currently support layer1_creator type '" << lib_type() << "' (struct=0,app=1,dog=2...). " << endl;
			throw std::exception();
			break;
	}

	ret = clSetKernelArg(convolution_kernel, 0, sizeof(cl_mem), (void*)&input_image);
	ret |= clSetKernelArg(convolution_kernel, 1, sizeof(cl_int2), (void*)&image_size);
	ret |= clSetKernelArg(convolution_kernel, 2, sizeof(cl_mem), (void*)&convolution_mask);
	ret |= clSetKernelArg(convolution_kernel, 3, sizeof(cl_int4), (void*)&convol_mask_size);
	ret |= clSetKernelArg(convolution_kernel, 4, sizeof(cl_float) * local_data_size.s[0] * local_data_size.s[1] , nullptr);
	ret |= clSetKernelArg(convolution_kernel, 5, sizeof(cl_int2), (void*)&local_data_size);
	ret |= clSetKernelArg(convolution_kernel, 6, sizeof(cl_uint), (void*)&save_max);
	ret |= clSetKernelArg(convolution_kernel, 7, sizeof(cl_mem), (void*)&max_image);
	ret |= clSetKernelArg(convolution_kernel, 8, sizeof(cl_mem), (void*)&max_image_index);
	ret |= clSetKernelArg(convolution_kernel, 9, sizeof(cl_mem), (void*)&made_images);			
	ret |= clSetKernelArg(convolution_kernel, 10, sizeof(cl_uint), (void*)&layer1_creator_type);
	ret |= clSetKernelArg(convolution_kernel, 11, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);

	ret = clEnqueueNDRangeKernel(cmd_queue, convolution_kernel, 2, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, evnt);
	ocl_check_error(ret, CL_SUCCESS);

	clReleaseMemObject(input_image);
	clReleaseMemObject(convolution_mask);

	delete[] input_output_host;

#ifdef OCL_PROFILE_KERNELS_HOST
	ret = clFinish(cmd_queue);
	clock_4.stopTimer();
#endif

}

void layer1_creator::ocl_make_sum2x2(const best_cmd_queue_info& cmd_queue_info, const img* im, cl_mem max_image, cl_mem sum_2x2, cl_mem kernel_err,
									 int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) {

	cl_kernel sum2x2_kernel = OpenCL::context_manager->getCommonKernel("make_2x2sum")[cmd_queue_info.context_number];

#ifdef OCL_PROFILE_KERNELS_HOST
	clock_4.startTimer();
#endif
	
	int vector_size = 1;

	int num_sums = 2;
	
	int best_local_work_size_x = cmd_queue_info.best_work_groups_dim2.x;
	int best_local_work_size_y = cmd_queue_info.best_work_groups_dim2.y;

	cl_int2 image_size,input_offset,local_data_size;
	
	image_size.s[0] = im->width;
	image_size.s[1] = im->height;

	input_offset.s[0] = 0;
	input_offset.s[1] = 0;//im->height * make_images_get_count();
		
	local_data_size.s[0] = best_local_work_size_x * vector_size + (num_sums - 1);
	local_data_size.s[1] = best_local_work_size_y + (num_sums - 1);

	int best_group_size_x = ceil((float)((im->width - (num_sums - 1))/vector_size) / best_local_work_size_x);
	int best_group_size_y = ceil((float)((im->height - (num_sums - 1))) / best_local_work_size_y);

	int best_global_work_size_x = best_group_size_x * best_local_work_size_x;
	int best_global_work_size_y = best_group_size_y * best_local_work_size_y;

	size_t global_work_load[2] = {best_global_work_size_x, best_global_work_size_y};
	size_t local_work_load[2] = {best_local_work_size_x, best_local_work_size_y};

	int ret;
	ret = clSetKernelArg(sum2x2_kernel, 0, sizeof(cl_mem), (void*)&max_image);
	ret |= clSetKernelArg(sum2x2_kernel, 1, sizeof(cl_int2), (void*)&image_size);
	ret |= clSetKernelArg(sum2x2_kernel, 2, sizeof(cl_int2), (void*)&input_offset);
	ret |= clSetKernelArg(sum2x2_kernel, 3, sizeof(cl_float) * local_data_size.s[0] * local_data_size.s[1] , nullptr);
	ret |= clSetKernelArg(sum2x2_kernel, 4, sizeof(cl_int2), (void*)&local_data_size);
	ret |= clSetKernelArg(sum2x2_kernel, 5, sizeof(cl_mem), (void*)&sum_2x2);
	ret |= clSetKernelArg(sum2x2_kernel, 6, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);
	
	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, sum2x2_kernel, 2, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, evnt);
	ocl_check_error(ret, CL_SUCCESS);

#ifdef OCL_PROFILE_KERNELS_HOST
	ret = clFinish(cmd_queue_info.queue);
	clock_4.stopTimer();
#endif
}
void layer1_creator::ocl_get_maximus(const best_cmd_queue_info& cmd_queue_info, const img* im, cl_mem input, cl_mem output_max, cl_mem kernel_err,
									 int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) {

#ifdef OCL_PROFILE_KERNELS_HOST
	clock_4.startTimer();
#endif

	int ret;
	cl_kernel get_max_kernel = OpenCL::context_manager->getCommonKernel("get_max")[cmd_queue_info.context_number];

	size_t wg_size;
	ret = clGetKernelWorkGroupInfo(get_max_kernel, cmd_queue_info.device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &wg_size, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	cl_int image_full_size = im->size();
	cl_int input_offset = 0;
	cl_int max_local_size = cmd_queue_info.device_type == CL_DEVICE_TYPE_CPU ? 4 : wg_size;//cmd_queue_info.max_work_group_size/2;

	
	ret = clSetKernelArg(get_max_kernel, 0, sizeof(cl_mem), (void*)&input);
	ret |= clSetKernelArg(get_max_kernel, 1, sizeof(cl_int), (void*)&image_full_size);
	ret |= clSetKernelArg(get_max_kernel, 2, sizeof(cl_int), (void*)&input_offset);
	ret |= clSetKernelArg(get_max_kernel, 3, sizeof(cl_float) * max_local_size, nullptr);
	ret |= clSetKernelArg(get_max_kernel, 4, sizeof(cl_int), (void*)&max_local_size);
	ret |= clSetKernelArg(get_max_kernel, 5, sizeof(cl_mem), (void*)&output_max);
	ret |= clSetKernelArg(get_max_kernel, 6, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);

	size_t global_work_load[1] = {max_local_size};
	size_t local_work_load[1] = {max_local_size};

	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, get_max_kernel, 1, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, evnt);
	ocl_check_error(ret, CL_SUCCESS);

#ifdef OCL_PROFILE_KERNELS_HOST			
	ret = clFinish(cmd_queue_info.queue);

	clock_4.stopTimer();
	clock_4.sumElapsedTime(); cout << " maximus kernels part time : " << clock_4.getElapsedTime() << endl;
#endif
}

void layer1_creator::ocl_prepare_results(const best_cmd_queue_info& cmd_queue_info, layer1_result* result, const img* im, const img* maskimg,
									 cl_mem made_images, cl_mem max_image, cl_mem max_image_index, cl_mem sum_2x2, cl_mem maximum_max_image, cl_mem maximum_sum_2x2,
									 cl_mem layer1_coords_compact, cl_mem layer1_coords_compact_sort_keys, cl_mem number_created_result, cl_mem number_positions_used_result, cl_mem kernel_err,
									 int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) {
#ifdef OCL_PROFILE_HOST_COMPLETE
	clock_2.startTimer();
#endif	
	int ret;
	cl_context context = OpenCL::context_manager->getContext(cmd_queue_info.context_number);

	cl_mem input_image_mask;

	if (maskimg != nullptr && maskimg->size() > 0 && false) {
		float* input_image_mask_host = new float[maskimg->width * maskimg->height];

		for (int i = 0; i < maskimg->size(); i++) {
			input_image_mask_host[i] = (*maskimg)[i];
		}

		input_image_mask = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * (maskimg->size()), input_image_mask_host , &ret);

		delete[] input_image_mask_host;
	} else {
		input_image_mask = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * (1), nullptr , &ret);
	}
	

#ifdef OCL_PROFILE_HOST_COMPLETE
	clock_2.stopTimer();
#endif

#ifdef OCL_PROFILE_KERNELS_HOST
	clock_4.startTimer();
#endif
	cl_kernel prepare_results_kernel = OpenCL::context_manager->getCommonKernel("prepare_results")[cmd_queue_info.context_number];

	size_t wg_size;
	ret = clGetKernelWorkGroupInfo(prepare_results_kernel, cmd_queue_info.device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &wg_size, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	int vector_size = 1;
	
	int neighboors_area = 3;

	int best_local_work_size_x = cmd_queue_info.best_work_groups_dim2.x;
	int best_local_work_size_y = wg_size/best_local_work_size_x; //cmd_queue_info.best_work_groups_dim2.y;

	cl_int2 local_data_size;
	local_data_size.s[0] = best_local_work_size_x * vector_size + 2*(neighboors_area/2);
	local_data_size.s[1] = best_local_work_size_y + 2*(neighboors_area/2);

	int best_group_size_x = ceil((float)((im->width )/vector_size) / best_local_work_size_x);
	int best_group_size_y = ceil((float)((im->height )) / best_local_work_size_y);

	int best_global_work_size_x = best_group_size_x * best_local_work_size_x;
	int best_global_work_size_y = best_group_size_y * best_local_work_size_y;

	size_t global_work_load[2] = {best_global_work_size_x, best_global_work_size_y};
	size_t local_work_load[2] = {best_local_work_size_x, best_local_work_size_y};

	cl_int2 image_size, new_image_size, input_image_mask_size;

	image_size.s[0] = im->width;
	image_size.s[1] = im->height;

	new_image_size.s[0] = result->x_size(0);
	new_image_size.s[1] = result->y_size(0);

	input_image_mask_size.s[0] = 0;//maskimg->width;
	input_image_mask_size.s[1] = 0;//maskimg->height;

	float layer1_threshold_float = layer1_threshold;
	int number_images = make_images_get_count();

	float response_percent_float = response_percent;
	
	ret = clSetKernelArg(prepare_results_kernel, 0, sizeof(cl_mem), (void*)&made_images);
	ret |= clSetKernelArg(prepare_results_kernel, 1, sizeof(cl_mem), (void*)&max_image);
	ret |= clSetKernelArg(prepare_results_kernel, 2, sizeof(cl_mem), (void*)&max_image_index);
	ret |= clSetKernelArg(prepare_results_kernel, 3, sizeof(cl_mem), (void*)&sum_2x2);
	ret |= clSetKernelArg(prepare_results_kernel, 4, sizeof(cl_int2), (void*)&image_size);
	ret |= clSetKernelArg(prepare_results_kernel, 5, sizeof(cl_int2), (void*)&new_image_size);
	ret |= clSetKernelArg(prepare_results_kernel, 6, sizeof(cl_mem), (void*)&input_image_mask);
	ret |= clSetKernelArg(prepare_results_kernel, 7, sizeof(cl_int2), (void*)&input_image_mask_size);
	ret |= clSetKernelArg(prepare_results_kernel, 8, sizeof(cl_mem), (void*)&maximum_max_image);
	ret |= clSetKernelArg(prepare_results_kernel, 9, sizeof(cl_mem), (void*)&maximum_sum_2x2);
	ret |= clSetKernelArg(prepare_results_kernel, 10, sizeof(cl_float) * local_data_size.s[0] * local_data_size.s[1] , nullptr);
	ret |= clSetKernelArg(prepare_results_kernel, 11, sizeof(cl_float) * local_data_size.s[0] * local_data_size.s[1] , nullptr);
	ret |= clSetKernelArg(prepare_results_kernel, 12, sizeof(cl_int2), (void*)&local_data_size);			
	ret |= clSetKernelArg(prepare_results_kernel, 13, sizeof(cl_float), (void*)&layer1_threshold_float);
	ret |= clSetKernelArg(prepare_results_kernel, 14, sizeof(cl_int), (void*)&layer1_3x3bound);
	ret |= clSetKernelArg(prepare_results_kernel, 15, sizeof(cl_int), (void*)&number_images);
	ret |= clSetKernelArg(prepare_results_kernel, 16, sizeof(cl_int), (void*)&border);
	ret |= clSetKernelArg(prepare_results_kernel, 17, sizeof(cl_int), (void*)&response_percent_float);
	ret |= clSetKernelArg(prepare_results_kernel, 18, sizeof(cl_mem), (void*)&layer1_coords_compact);
	ret |= clSetKernelArg(prepare_results_kernel, 19, sizeof(cl_mem), (void*)&layer1_coords_compact_sort_keys);
	ret |= clSetKernelArg(prepare_results_kernel, 20, sizeof(cl_mem), (void*)&number_created_result);
	ret |= clSetKernelArg(prepare_results_kernel, 21, sizeof(cl_mem), (void*)&number_positions_used_result);	
	ret |= clSetKernelArg(prepare_results_kernel, 22, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);

	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, prepare_results_kernel, 2, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, evnt);
	ocl_check_error(ret, CL_SUCCESS);

	clReleaseMemObject(input_image_mask);

#ifdef OCL_PROFILE_KERNELS_HOST
	ret = clFinish(cmd_queue_info.queue);
	clock_4.stopTimer();
	clock_4.sumElapsedTime(); cout << " prepare results kernels part time : " << clock_4.getElapsedTime() << endl;
#endif
}


void layer1_creator::ocl_sort_compact_coords(const best_cmd_queue_info& cmd_queue_info, const img* im, const cl_int coords_size_,
												 cl_mem layer1_coords_compact, cl_mem layer1_coords_compact_sort_keys, cl_mem kernel_err,
												 int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) {
	
#ifdef OCL_PROFILE_KERNELS_HOST
clock_4.startTimer();
#endif
	int ret;
	// we have already make sure size of layer1_coords_compact and layer1_coords_compact_sort_keys is adjusted to new coords_size
	int coords_size =  ::pow(2,::ceil(::log((float)coords_size_)/::log((float)2)));
	
	if (cmd_queue_info.local_mem_type == CL_LOCAL)  {

		cl_kernel bitonic_sort_local_kernel = OpenCL::context_manager->getCommonKernel("bitonicSortLocal")[cmd_queue_info.context_number];
		cl_kernel bitonic_merge_local_kernel = OpenCL::context_manager->getCommonKernel("bitonicMergeLocal")[cmd_queue_info.context_number];
		cl_kernel bitonic_merge_global_kernel = OpenCL::context_manager->getCommonKernel("bitonicMergeGlobal")[cmd_queue_info.context_number];

		int best_local_work_size_x = min<int>(cmd_queue_info.max_work_group_size, coords_size/2);
		int best_group_size_x = ceil((float)(coords_size/2) / best_local_work_size_x);
		int best_global_work_size_x = best_group_size_x * best_local_work_size_x;

		size_t global_work_load[1] = {best_global_work_size_x};
		size_t local_work_load[1] = {best_local_work_size_x};

		cl_uint coords_size_uint = coords_size;
		cl_uint best_local_work_size_x_uint = 2*best_local_work_size_x;
		cl_uint dir = 0;

		dir = (dir != 0);

		//Launch bitonicSortLocal
        ret = clSetKernelArg(bitonic_sort_local_kernel, 0,  sizeof(cl_mem), (void *)&layer1_coords_compact_sort_keys);
		ret |= clSetKernelArg(bitonic_sort_local_kernel, 1,  sizeof(cl_mem), (void *)&layer1_coords_compact);
        ret |= clSetKernelArg(bitonic_sort_local_kernel, 2,  sizeof(ocl_coord_sorting_pointer) * best_local_work_size_x_uint, nullptr);
		ret |= clSetKernelArg(bitonic_sort_local_kernel, 3,  sizeof(ocl_layer1_data_coordinates) * best_local_work_size_x_uint, nullptr);
        ret |= clSetKernelArg(bitonic_sort_local_kernel, 4,  sizeof(cl_uint), (void *)&best_local_work_size_x_uint);
        ocl_check_error(ret, CL_SUCCESS);

		int event_a = 0;
		int event_b = 1;
		cl_event evnt_list[2];
        ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, bitonic_sort_local_kernel, 1, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, &evnt_list[event_a]);
        ocl_check_error(ret, CL_SUCCESS);

		
        for(cl_uint size = 2* best_local_work_size_x_uint; size <= coords_size; size <<= 1)
        {
            for(cl_uint stride = size / 2; stride > 0; stride >>= 1)
            {

                if(stride >= best_local_work_size_x_uint)
                {
                    //Launch bitonicMergeGlobal
                    ret = clSetKernelArg(bitonic_merge_global_kernel, 0,  sizeof(cl_mem), (void *)&layer1_coords_compact_sort_keys);
					ret |= clSetKernelArg(bitonic_merge_global_kernel, 1,  sizeof(cl_mem), (void *)&layer1_coords_compact);
                    ret |= clSetKernelArg(bitonic_merge_global_kernel, 2, sizeof(cl_uint), (void *)&coords_size_uint);
                    ret |= clSetKernelArg(bitonic_merge_global_kernel, 3, sizeof(cl_uint), (void *)&size);
                    ret |= clSetKernelArg(bitonic_merge_global_kernel, 4, sizeof(cl_uint), (void *)&stride);
                    ret |= clSetKernelArg(bitonic_merge_global_kernel, 5, sizeof(cl_uint), (void *)&dir);
                    ocl_check_error(ret, CL_SUCCESS);

                    ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, bitonic_merge_global_kernel, 1, nullptr, global_work_load, local_work_load, 1, &evnt_list[event_a], &evnt_list[event_b]);
                    ocl_check_error(ret, CL_SUCCESS);
                }
                else
                {
					
                    //Launch bitonicMergeLocal
                    ret = clSetKernelArg(bitonic_merge_local_kernel, 0, sizeof(cl_mem), (void *)&layer1_coords_compact_sort_keys);
					ret |= clSetKernelArg(bitonic_merge_local_kernel, 1, sizeof(cl_mem), (void *)&layer1_coords_compact);
                    ret |= clSetKernelArg(bitonic_merge_local_kernel, 2, sizeof(ocl_coord_sorting_pointer) * best_local_work_size_x_uint, nullptr);
					ret |= clSetKernelArg(bitonic_merge_local_kernel, 3, sizeof(ocl_layer1_data_coordinates) * best_local_work_size_x_uint, nullptr);
					ret |= clSetKernelArg(bitonic_merge_local_kernel, 4, sizeof(cl_uint), (void *)&best_local_work_size_x_uint);
                    ret |= clSetKernelArg(bitonic_merge_local_kernel, 5, sizeof(cl_uint), (void *)&coords_size_uint);
                    ret |= clSetKernelArg(bitonic_merge_local_kernel, 6, sizeof(cl_uint), (void *)&stride);
                    ret |= clSetKernelArg(bitonic_merge_local_kernel, 7, sizeof(cl_uint), (void *)&size);
                    ret |= clSetKernelArg(bitonic_merge_local_kernel, 8, sizeof(cl_uint), (void *)&dir);
                    ocl_check_error(ret, CL_SUCCESS);

					ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, bitonic_merge_local_kernel, 1, nullptr, global_work_load, local_work_load, 1, &evnt_list[event_a], &evnt_list[event_b]);
                    ocl_check_error(ret, CL_SUCCESS);
                    break;
                }
				event_a = (event_a + 1) % 2;
				event_b = (event_b + 1) % 2;
            }
        }
		*evnt = evnt_list[event_a];

	} else {
		// use only kernel for global merge when we do not have fast LOCAL memory 
		// (i.e. when barrier(CLK_LOCAL_MEM_FENCE) is same as barrier(CLK_GLOBAL_MEM_FENCE) which will cause slow execution)
		cl_kernel bitonic_merge_global_kernel = OpenCL::context_manager->getCommonKernel("bitonicMergeGlobal")[cmd_queue_info.context_number];

		int best_local_work_size_x = min<int>(cmd_queue_info.max_work_group_size, coords_size/2);
		int best_group_size_x = ceil((float)(coords_size/2) / best_local_work_size_x);
		int best_global_work_size_x = best_group_size_x * best_local_work_size_x;

		size_t global_work_load[1] = {best_global_work_size_x};
		size_t local_work_load[1] = {best_local_work_size_x};

		cl_uint coords_size_uint = coords_size;
		cl_uint best_local_work_size_x_uint = 2*best_local_work_size_x;
		cl_uint dir = 0;

		dir = (dir != 0);

		int event_a = 1;
		int event_b = 0;
		cl_event* evnt_list = new cl_event[1 + max<int>(1,wait_event_count)];
		
		cl_uint number_events = wait_event_count;
		if (number_events > 0 && wait_evnt != nullptr) {
			memcpy(evnt_list + 1, wait_evnt, sizeof(cl_event) * number_events);
		}
		
        for(cl_uint size = 1; size <= coords_size; size <<= 1)
        {
            for(cl_uint stride = size / 2; stride > 0; stride >>= 1)
            {

				//Launch bitonicMergeGlobal
                ret = clSetKernelArg(bitonic_merge_global_kernel, 0,  sizeof(cl_mem), (void *)&layer1_coords_compact_sort_keys);
				ret |= clSetKernelArg(bitonic_merge_global_kernel, 1,  sizeof(cl_mem), (void *)&layer1_coords_compact);
                ret |= clSetKernelArg(bitonic_merge_global_kernel, 2, sizeof(cl_uint), (void *)&coords_size_uint);
                ret |= clSetKernelArg(bitonic_merge_global_kernel, 3, sizeof(cl_uint), (void *)&size);
                ret |= clSetKernelArg(bitonic_merge_global_kernel, 4, sizeof(cl_uint), (void *)&stride);
                ret |= clSetKernelArg(bitonic_merge_global_kernel, 5, sizeof(cl_uint), (void *)&dir);
                ocl_check_error(ret, CL_SUCCESS);

				ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, bitonic_merge_global_kernel, 1, nullptr, global_work_load, local_work_load, number_events, number_events > 0 ? &evnt_list[event_a] : nullptr, &evnt_list[event_b]);
                ocl_check_error(ret, CL_SUCCESS);
            
				// each next clEnqueueNDRangeKernel should wait this one
				event_a = (event_a + 1) % 2;
				event_b = (event_b + 1) % 2;
				// since first clEnqueueNDRangeKernel might have have to waited more then one event we need to re-set this number back to 1
				number_events = 1;
            }
        }
		
		*evnt = evnt_list[event_a];
		delete[] evnt_list;
	}

#ifdef OCL_PROFILE_KERNELS_HOST
	ret = clFinish(cmd_queue_info.queue);
	clock_4.stopTimer();
	clock_4.sumElapsedTime(); cout << " sort coordinates kernels part time : " << clock_4.getElapsedTime() << endl;			
#endif	
	/*
	ocl_coord_sorting_pointer* layer1_coords_compact_keys_host = new ocl_coord_sorting_pointer[coords_size];
	ret = clEnqueueReadBuffer(cmd_queue_info.queue, layer1_coords_compact_sort_keys, CL_FALSE, 0, sizeof(ocl_coord_sorting_pointer) * coords_size, layer1_coords_compact_keys_host, 1, evnt, nullptr);

	ocl_layer1_data_coordinates* layer1_coords_compact_host = new ocl_layer1_data_coordinates[coords_size];
	ret = clEnqueueReadBuffer(cmd_queue_info.queue, layer1_coords_compact, CL_FALSE, 0, sizeof(ocl_layer1_data_coordinates) * coords_size, layer1_coords_compact_host, 1, evnt, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	getchar();
	float before = 0;
	for (int i = 0; i < coords_size; i++) {
		float current = layer1_coords_compact_keys_host[i].sort_value;
		ocl_layer1_data_coordinates cc = layer1_coords_compact_host[i];
		//cout << current << endl;
		if (before < current)
			cout << " INVALID ONE at " << i << endl;
		before = current;
	}
	*/
}

void layer1_creator::ocl_update_results_offsets(const best_cmd_queue_info& cmd_queue_info, const img* im, const cl_int coords_size, cl_mem layer1_coords, cl_mem kernel_err,
												int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt,
												cl_event* profile_start, cl_event* profile_stop) {

	int ret;
#ifdef OCL_PROFILE_KERNELS_OPENCL
	ret = clEnqueueMarker(cmd_queue_info.queue, profile_start);
	ocl_check_error(ret, CL_SUCCESS);
#endif
	
	if (cmd_queue_info.local_mem_type == CL_LOCAL)  {
			
		cl_kernel prefix_sum_step1_kernel = OpenCL::context_manager->getCommonKernel("scanExclusiveLocal1")[cmd_queue_info.context_number];			
	
		int best_local_work_size_x = cmd_queue_info.max_work_group_size;
		int best_group_size_x = ceil((float)(coords_size/4) / best_local_work_size_x);
		int best_global_work_size_x = best_group_size_x * best_local_work_size_x;

		size_t global_work_load[1] = {best_global_work_size_x};
		size_t local_work_load[1] = {best_local_work_size_x};

		ret = clSetKernelArg(prefix_sum_step1_kernel, 0, sizeof(cl_mem), (void*)&layer1_coords);
		ret |= clSetKernelArg(prefix_sum_step1_kernel, 1, sizeof(cl_int), (void*)&coords_size);
		ret |= clSetKernelArg(prefix_sum_step1_kernel, 2, sizeof(cl_int) * best_local_work_size_x * 2, nullptr);			
		ret |= clSetKernelArg(prefix_sum_step1_kernel, 3, sizeof(cl_int), (void*)&best_local_work_size_x);
		ocl_check_error(ret, CL_SUCCESS);

		cl_event envt_set_offsets_step1;

		ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, prefix_sum_step1_kernel, 1, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, &envt_set_offsets_step1);
		ocl_check_error(ret, CL_SUCCESS);

		int additional_offsets = ceil((float)coords_size/(best_local_work_size_x*4));

		// if unable to set offsets for all in one turn then adjust them recursively			
		if (additional_offsets > 1) {
			ocl_set_offsets_recursively(cmd_queue_info, layer1_coords, coords_size, additional_offsets, 4 * best_local_work_size_x, 1, &envt_set_offsets_step1, evnt);
		} else {
			*evnt = envt_set_offsets_step1;
		}

	} else {
		
		cl_kernel seq_prefix_sum_kernel = OpenCL::context_manager->getCommonKernel("sequential_prefix_sum")[cmd_queue_info.context_number];

		ret = clSetKernelArg(seq_prefix_sum_kernel, 0, sizeof(cl_mem), (void*)&layer1_coords);
		ret |= clSetKernelArg(seq_prefix_sum_kernel, 1, sizeof(cl_int), (void*)&coords_size);
		ret |= clSetKernelArg(seq_prefix_sum_kernel, 2, sizeof(cl_mem), (void*)&kernel_err);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clEnqueueTask(cmd_queue_info.queue, seq_prefix_sum_kernel, wait_event_count, wait_evnt, evnt);
		ocl_check_error(ret, CL_SUCCESS);
	}

#ifdef OCL_PROFILE_KERNELS_OPENCL
	ret = clEnqueueMarker(cmd_queue_info.queue, profile_stop);
	ocl_check_error(ret, CL_SUCCESS);
#endif

}

void layer1_creator::ocl_make_results(const best_cmd_queue_info& cmd_queue_info, layer1_result* result, const img* im, int number_positions_used,
									  cl_mem made_images, cl_mem layer1_coords_compact, cl_mem layer1_coords_compact_sort_keys, cl_mem layer1_result_data, cl_mem kernel_err,
									  int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) {
	
#ifdef OCL_PROFILE_KERNELS_HOST
clock_4.startTimer();
#endif

	int best_local_work_size_x = cmd_queue_info.max_work_group_size;
	
	int best_group_size_x = ceil((float)(number_positions_used ) / best_local_work_size_x);
	
	int best_global_work_size_x = best_group_size_x * best_local_work_size_x;
	
	size_t global_work_load[1] = {best_global_work_size_x};
	size_t local_work_load[1] = {best_local_work_size_x};

	cl_int2 image_size, new_image_size;
	
	image_size.s[0] = im->width;
	image_size.s[1] = im->height;

	new_image_size.s[0] = result->x_size(0);
	new_image_size.s[1] = result->y_size(0);
	
	int number_images = make_images_get_count();

	int ret;
	cl_kernel make_results_kernel = OpenCL::context_manager->getCommonKernel("make_results")[cmd_queue_info.context_number];

	int normalization_part_index;
	if (normalization_percent <= 0.0 || normalization_percent >= 1.0)
		normalization_part_index = 0;
	else
		normalization_part_index = min<int>(max<int>(number_positions_used * normalization_percent, 0), number_positions_used - 1);

	float power_correction_float = power_correction;
	float position_percent_float = response_percent;

	ret = clSetKernelArg(make_results_kernel, 0, sizeof(cl_mem), (void*)&made_images);	
	ret |= clSetKernelArg(make_results_kernel, 1, sizeof(cl_int2), (void*)&image_size);
	ret |= clSetKernelArg(make_results_kernel, 2, sizeof(cl_int2), (void*)&new_image_size);
	ret |= clSetKernelArg(make_results_kernel, 3, sizeof(cl_int), (void*)&number_images);
	ret |= clSetKernelArg(make_results_kernel, 4, sizeof(cl_mem), (void*)&layer1_coords_compact);
	ret |= clSetKernelArg(make_results_kernel, 5, sizeof(cl_mem), (void*)&layer1_coords_compact_sort_keys);	
	ret |= clSetKernelArg(make_results_kernel, 6, sizeof(cl_int), (void*)&number_positions_used);	
	ret |= clSetKernelArg(make_results_kernel, 7, sizeof(cl_mem), (void*)&layer1_result_data);
	ret |= clSetKernelArg(make_results_kernel, 8, sizeof(cl_int), (void*)&border);
	ret |= clSetKernelArg(make_results_kernel, 9, sizeof(cl_float), (void*)&position_percent_float);
	ret |= clSetKernelArg(make_results_kernel, 10, sizeof(cl_float), (void*)&power_correction_float);	
	ret |= clSetKernelArg(make_results_kernel, 11, sizeof(cl_int), (void*)&normalization_part_index);
	ret |= clSetKernelArg(make_results_kernel, 12, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);
	
	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, make_results_kernel, 1, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, evnt);
	ocl_check_error(ret, CL_SUCCESS);

#ifdef OCL_PROFILE_KERNELS_HOST
	ret = clFinish(cmd_queue_info.queue);
	clock_4.stopTimer();
	clock_4.sumElapsedTime(); cout << " make result kernels part time : " << clock_4.getElapsedTime() << endl;			
#endif
}

void layer1_creator::ocl_sort_results(const best_cmd_queue_info& cmd_queue_info, layer1_result* result, const img* im, cl_mem layer1_result_data, cl_int compact_count, cl_mem layer1_coords_compact, cl_mem merge_sort_temp_data, cl_mem kernel_err,
										int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) 
{

#ifdef OCL_PROFILE_KERNELS_HOST
clock_4.startTimer();
#endif

	int ret;
	cl_context context = OpenCL::context_manager->getContext(cmd_queue_info.context_number);	

	cl_int2 new_image_size;

	new_image_size.s[0] = result->x_size(0);
	new_image_size.s[1] = result->y_size(0);

	cl_kernel merge_sort_within_pos_result_kernel = OpenCL::context_manager->getCommonKernel("merge_sort_within_pos_result")[cmd_queue_info.context_number];

	ret = clSetKernelArg(merge_sort_within_pos_result_kernel, 0, sizeof(cl_mem), (void*)&layer1_result_data);
	ret |= clSetKernelArg(merge_sort_within_pos_result_kernel, 1, sizeof(cl_int), (void*)&compact_count);
	ret |= clSetKernelArg(merge_sort_within_pos_result_kernel, 2, sizeof(cl_mem), (void*)&layer1_coords_compact);
	ret |= clSetKernelArg(merge_sort_within_pos_result_kernel, 3, sizeof(cl_int2), (void*)&new_image_size);
	ret |= clSetKernelArg(merge_sort_within_pos_result_kernel, 4, sizeof(cl_mem), (void*)&merge_sort_temp_data);	
	ret |= clSetKernelArg(merge_sort_within_pos_result_kernel, 5, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);
	
	size_t wg_size;
	ret = clGetKernelWorkGroupInfo(merge_sort_within_pos_result_kernel, cmd_queue_info.device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &wg_size, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	int best_local_work_size_x = wg_size;
	int best_group_size_x = (compact_count) / best_local_work_size_x + 1;
	int best_global_work_size_x = best_group_size_x * best_local_work_size_x;

	size_t global_work_load[1] = {best_global_work_size_x};
	size_t local_work_load[1] = {best_local_work_size_x};

	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, merge_sort_within_pos_result_kernel, 1, nullptr, global_work_load , local_work_load, wait_event_count, wait_evnt, evnt); 
	ocl_check_error(ret, CL_SUCCESS);

#ifdef OCL_PROFILE_KERNELS_HOST
ret = clFinish(cmd_queue_info.queue);
clock_4.stopTimer();
clock_4.sumElapsedTime(); cout << " merge sort kernels part time : " << clock_4.getElapsedTime() << endl;
#endif

}
void layer1_creator::ocl_inhibit_results(const best_cmd_queue_info& cmd_queue_info, layer1_result* result, const img* im, 
											cl_mem layer1_result_data, cl_mem layer1_coords, cl_int adjusted_layer1_coords_size, cl_mem layer1_coords_inhib, cl_mem kernel_err,
											int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) {
	
#ifdef OCL_PROFILE_KERNELS_HOST
	clock_4.startTimer();
#endif

	cl_context context = OpenCL::context_manager->getContext(cmd_queue_info.context_number);
	cl_command_queue cmd_queue = cmd_queue_info.queue;

	cl_kernel inhibit_result_kernel = OpenCL::context_manager->getCommonKernel("inhibit_result")[cmd_queue_info.context_number];

	cl_uint2 new_image_size;

	new_image_size.s[0] = result->x_size(0);
	new_image_size.s[1] = result->y_size(0);

	int ret;
	ret = clSetKernelArg(inhibit_result_kernel, 0, sizeof(cl_mem), (void*)&layer1_result_data);
	ret |= clSetKernelArg(inhibit_result_kernel, 1, sizeof(cl_mem), (void*)&layer1_coords);		
	ret |= clSetKernelArg(inhibit_result_kernel, 2, sizeof(cl_mem), (void*)&layer1_coords_inhib);	
	ret |= clSetKernelArg(inhibit_result_kernel, 3, sizeof(cl_uint2), (void*)&new_image_size);
	ret |= clSetKernelArg(inhibit_result_kernel, 4, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);

	int best_local_work_size_x = cmd_queue_info.max_work_group_size;

	//int best_group_size_x = (adjusted_layer1_coords_size) / best_local_work_size_x + 1;
	int best_group_size_x = (new_image_size.s[0] * new_image_size.s[1]) / best_local_work_size_x + 1;

	int best_global_work_size_x = best_group_size_x * best_local_work_size_x;

	size_t global_work_load[1] = {best_global_work_size_x};
	size_t local_work_load[1] = {best_local_work_size_x};

	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, inhibit_result_kernel, 1, nullptr, global_work_load , local_work_load, wait_event_count, wait_evnt, evnt); 
	ocl_check_error(ret, CL_SUCCESS);

#ifdef OCL_PROFILE_KERNELS_HOST
	ret = clFinish(cmd_queue_info.queue);
	clock_4.stopTimer();
	clock_4.sumElapsedTime(); cout << " inhib result kernels part time : " << clock_4.getElapsedTime() << endl;
#endif
}

// undef all profile and debug values so they do not get mixed up with others cpp files (in case someone includes this file directly)
#undef OCL_PROFILE_KERNELS_HOST 
#undef OCL_PROFILE_KERNELS_OPENCL 
#undef OCL_PROFILE_HOST_COMPLETE 
#undef OCL_PROFILE_HOST_QUICK 

#endif
