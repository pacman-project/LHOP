
#ifdef OPENCL

#include <queue>
#include "layer_1.h"
#include "layer_n_creators.h"
#include "opencl/clContextManager.h"
#include "opencl/cl_utils.h"

#include <limits>

//#define OPENCL_VERIFICATION_SUPPORT
//#define OCL_PROFILE_HOST_QUICK // profile whole HOST code as one (including waiting on opencl kernels)

void layern_creator::ocl_add_layer1(layer1_result* res, int layer) {

	int test_debug_out = 0;
	if (layer < 1) return;

	// k = layer we MUST create
    int k = layer - 1;
	// k1 = layer which we will use to create layer k
    int k1 = k - 1;

	

#ifdef OCL_PROFILE_HOST_QUICK 
CStopWatch clock_1;

clock_1.startTimer();
#endif

	int x_size_k1 = res->x_size(k1);
	int y_size_k1 = res->y_size(k1);
	int xy_size_k1 = x_size_k1 * y_size_k1;

    int x_size_k = int_round(x_size_k1/layer_contraction); // !floor!
    int y_size_k = int_round(y_size_k1/layer_contraction); // !floor!
    int xy_size_k = x_size_k * y_size_k;

    if (ignore_texture) 
		mark_texture(res, k1);

	// this new_grid allcates unnecessary memory since it is only needed due to x_size and y_size
    res->new_grid(x_size_k, y_size_k, k);  
	
	while ((int)res->ocl_shape_nodes.size() <= k) res->ocl_shape_nodes.push_back(make_pair<ocl_layer1_data*, int>(nullptr,0));
	while ((int)res->ocl_edges.size() <= k) res->ocl_edges.push_back(make_pair<ocl_edge_data_ip2*, int>(nullptr,0));
	
	while ((int)res->ocl_shape_nodes_coord.size() <= k) res->ocl_shape_nodes_coord.push_back(make_pair<ocl_layer1_data_coordinates*, int>(nullptr,0));
	while ((int)res->ocl_shape_nodes_inhib_coord.size() <= k) res->ocl_shape_nodes_inhib_coord.push_back(make_pair<ocl_layer1_data_coordinates*, int>(nullptr,0));
	
	while ((int)res->ocl_shape_nodes_coord_non_zero_count.size() <= k) res->ocl_shape_nodes_coord_non_zero_count.push_back(0);
	while ((int)res->ocl_shape_nodes_inhib_coord_non_zero_count.size() <= k) res->ocl_shape_nodes_inhib_coord_non_zero_count.push_back(0);

	while ((int)res->info.size() <= k) res->info.push_back(layer_info());

	// Edges to all nodes (i.e. parts) in the next layer (one level up) which have same 
    // part for its center (in library; class part_lib).
    int to_center = atom("lyrCenterBack"); 

	// Edges to all subparts (one layer down) other then center (in library; class part_lib).
    int to_part = atom("lyrSrc");

	// Edges to forbidden subparts -- parts which are not allowed to be present at specific
    // positions -- (one layer down, in library; class part_lib). 
    int to_forbidden = atom("lyrForbidden").get_index();

    // Edges to all subparts fro, the previous layer (one layer down) that were used in the 
    // inference process of this part (in class layer1_result)
    int to_prev_layer = atom("toPrevLayer").get_index();

	// Edges to all parts in same layer that are similar (in library for class part_data)
    int to_similar = atom("toSimilar").get_index();
    
	int distr_size_x = 5;
	int distr_size_y = 5;

	matrix<double> distr;
    gaussian_mask(distr_size_x, distr_size_y, 2.0, distr);

	// use CPU as prefered type
	const best_cmd_queue_info& cmd_queue_info = OpenCL::context_manager->getBestCommandQueue("cpu", use_opencl_devices);

	cl_context context = OpenCL::context_manager->getContext(cmd_queue_info.context_number);

	// load all memory objects
	
	// all objects from lib for current and next layer
	ocl_part_data* lib_parts_layer_k_host = library->ocl_parts[k].first;
	ocl_part_data_2* lib_edges_layer_k_host = library->ocl_edges[k].first;
	
	ocl_app_data* lib_apps_layer_k_host = library->ocl_apps[k].first;

	ocl_part_data* lib_parts_layer_k1_host = library->ocl_parts[k1].first;
	ocl_part_data_2* lib_edges_layer_k1_host = library->ocl_edges[k1].first;

	// all objects from img for current layer
	ocl_layer1_data* img_s_nodes_layer_k1_host = res->ocl_shape_nodes[k1].first;
	ocl_edge_data_ip2* img_edges_layer_k1_host = res->ocl_edges[k1].first;
	ocl_layer1_data_coordinates* img_coord_layer_k1_host = res->ocl_shape_nodes_coord[k1].first;

	// all new objects from img for new layer
	ocl_layer1_data* new_layer_img_s_nodes_host = nullptr;
	ocl_edge_data_ip2* new_layer_img_edges_host = nullptr;
	ocl_layer1_data_coordinates* new_layer_img_coord_host = new ocl_layer1_data_coordinates[xy_size_k];
	ocl_layer1_data_coordinates* new_layer_img_coord_inhib_host = new ocl_layer1_data_coordinates[xy_size_k];

	// new_layer_img_coord_host and new_layer_img_coord_inhib_host MUST be set to zero since not all positions will be set to valid values
	memset(new_layer_img_coord_host, 0, sizeof(ocl_layer1_data_coordinates) * xy_size_k);
	memset(new_layer_img_coord_inhib_host, 0, sizeof(ocl_layer1_data_coordinates) * xy_size_k);

	cl_int ret;

	int ocl_shape_nodes_k1_count = res->ocl_shape_nodes[k1].second;
	ocl_set_candidate_thresholds(img_s_nodes_layer_k1_host, ocl_shape_nodes_k1_count, k1);

	// generate and count all possible canddiates for this layer
	int layer_candidates_count = 0;
	int max_schur_count = 0;
	ocl_layer_candidate* layer_candidates_host = ocl_generate_candidates(lib_parts_layer_k1_host, lib_edges_layer_k1_host, lib_parts_layer_k_host, img_s_nodes_layer_k1_host, ocl_shape_nodes_k1_count, x_size_k, y_size_k, k, layer_candidates_count, max_schur_count);

	int count_new_parts = 0;
	int count_new_edges = 0;

	// if we have any candidates we can continue
	if (layer_candidates_count > 0 && layer_candidates_host != nullptr) {
		
		// buffers for library parts could be global
		cl_mem lib_parts_layer_k = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, library->ocl_parts[k].second * sizeof(ocl_part_data), lib_parts_layer_k_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);
		cl_mem lib_edges_layer_k = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, library->ocl_edges[k].second * sizeof(ocl_part_data_2), lib_edges_layer_k_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem lib_apps_layer_k = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, library->ocl_apps[k].second * sizeof(ocl_app_data), lib_apps_layer_k_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem lib_parts_layer_k1 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, library->ocl_parts[k1].second * sizeof(ocl_part_data), lib_parts_layer_k1_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);
		cl_mem lib_edges_layer_k1 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, library->ocl_edges[k1].second * sizeof(ocl_part_data_2), lib_edges_layer_k1_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem img_s_nodes_layer_k1;
		if (res->ocl_shape_nodes[k1].second > 0) {
			img_s_nodes_layer_k1 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, res->ocl_shape_nodes[k1].second * sizeof(ocl_layer1_data), img_s_nodes_layer_k1_host, &ret);
			ocl_check_error(ret, CL_SUCCESS);
		}

		cl_mem img_edges_layer_k1;
		if (add_reconstruction_edges && res->ocl_edges[k1].second > 0) {
			img_edges_layer_k1 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, res->ocl_edges[k1].second * sizeof(ocl_edge_data_ip2), img_edges_layer_k1_host, &ret); 		
		} else {
			// just create dummy memory 
			img_edges_layer_k1 = clCreateBuffer(context, CL_MEM_READ_ONLY, 1, nullptr, &ret);
		}
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem img_coord_layer_k1 = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, res->ocl_shape_nodes_coord[k1].second * sizeof(ocl_layer1_data_coordinates), img_coord_layer_k1_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		// arguments common to some kernels
		cl_uint2 img_size_layer_k1, new_layer_img_size, dummy_m_size;

		new_layer_img_size.s[0] = x_size_k;
		new_layer_img_size.s[1] = y_size_k;

		img_size_layer_k1.s[0] = x_size_k1;
		img_size_layer_k1.s[1] = y_size_k1;

		dummy_m_size.s[0] = distr_size_x;
		dummy_m_size.s[1] = distr_size_y;

		cl_event evnt[5];

		cl_float* distr_gaussian_mask_host = new cl_float[distr.size()];
		for (int i = 0; i < distr.size(); i++) {
			distr_gaussian_mask_host[i] = (float)distr[i];
		}

		cl_mem distr_gaussian_mask = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, distr.size() * sizeof(cl_float), distr_gaussian_mask_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);
		
		delete[] distr_gaussian_mask_host;

		cl_int kernel_err_host = 0;
		cl_mem kernel_err = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, 1 * sizeof(cl_int), &kernel_err_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		// adjust to make empty space for bitonic sort (array length must be power of 2)
		//int adjusted_layer_candidates_count = pow(2,::ceil(::log((float)layer_candidates_count)/::log((float)2)));
		int adjusted_layer_candidates_count = layer_candidates_count;

		// also align size to 4 just in case
		if (adjusted_layer_candidates_count % 4 != 0)
			adjusted_layer_candidates_count +=  (4 - (adjusted_layer_candidates_count % 4));

		// make sure memory size can be handled by OpenCL device
		if (adjusted_layer_candidates_count * sizeof(ocl_layer_candidate) > cmd_queue_info.max_mem_alocation) {
			cout << endl << "WARNING: To many candidates of parts (" << adjusted_layer_candidates_count << ") created for OpenCL to be able to handle them. Adjust config values otherwise OpenCL may not be able to compute !!" << endl;
		}

		// create memory objects for all counters
		cl_mem count_new_parts_result = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 1 * sizeof(cl_int), &count_new_parts, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem count_new_edges_result = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 1 * sizeof(cl_int), &count_new_edges , &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_int count_valid_non_zero_pos = 0;
		cl_mem count_valid_non_zero_pos_result = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 1 * sizeof(cl_int), &count_valid_non_zero_pos, &ret);
		ocl_check_error(ret, CL_SUCCESS);


		// make opencl memory objects
		cl_mem layer_candidates = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, layer_candidates_count * sizeof(ocl_layer_candidate), layer_candidates_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);
		
		// create memory object for holding best max_schur products of all candidates
		cl_mem max_schur_memory_obj = clCreateBuffer(context, CL_MEM_READ_WRITE, max_schur_count * sizeof(cl_float), nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);		


		// create memory objects for non-zero position list of all candidates
		cl_mem candidate_non_zero_positions = clCreateBuffer(context, CL_MEM_READ_WRITE, layer_candidates_count * sizeof(ocl_layer1_data_coordinates), nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);
		
		cl_int candidate_non_zero_pos_count_host = 0;
		cl_mem candidate_non_zero_pos_count = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 1 * sizeof(cl_int), &candidate_non_zero_pos_count_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);			


		// create memory objects for compact offsets (i.e. non-zero position of only valid final candidates)
		cl_mem new_layer_s_nodes_offsets_keys = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_coord_sorting_pointer) * adjusted_layer_candidates_count, nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem new_layer_s_nodes_offsets = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_layer1_data_coordinates) * adjusted_layer_candidates_count, nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_int new_layer_s_nodes_offsets_count_host = 0;
		cl_mem new_layer_s_nodes_offsets_count = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 1 * sizeof(cl_int), &new_layer_s_nodes_offsets_count_host, &ret);
		ocl_check_error(ret, CL_SUCCESS);
		

		{
			// we need to set new_layer_s_nodes_offsets_keys to zero since it will be used in sorting
			cl_int in_all_bytes = (adjusted_layer_candidates_count * sizeof(ocl_coord_sorting_pointer));
			cl_int in_unalinged_bytes = in_all_bytes - (in_all_bytes / sizeof(cl_int16)) * sizeof(cl_int16);

			ocl_set_memory_zero(cmd_queue_info, new_layer_s_nodes_offsets_keys, in_all_bytes / sizeof(cl_int16), in_unalinged_bytes / sizeof(cl_int), 0, nullptr, &evnt[3]);
			
			// first eliminate all forbidden parts from candidates
			ocl_eliminate_forbidden(cmd_queue_info, lib_parts_layer_k, lib_edges_layer_k, lib_parts_layer_k1, img_s_nodes_layer_k1, img_coord_layer_k1, layer_candidates, 
									layer_candidates_count, img_size_layer_k1, k, kernel_err, 0, nullptr, &evnt[0]);


			// Next, for all parts still appropriate for processing we need to make matching with original image
			ocl_candidates_matching(cmd_queue_info, lib_parts_layer_k,lib_edges_layer_k, lib_apps_layer_k, lib_parts_layer_k1, img_s_nodes_layer_k1, img_coord_layer_k1, img_edges_layer_k1, distr_gaussian_mask,
										layer_candidates, layer_candidates_count, max_schur_memory_obj, img_size_layer_k1, dummy_m_size, k1,
										candidate_non_zero_positions, candidate_non_zero_pos_count, kernel_err, 1, &evnt[0], &evnt[1]);

			// then from candidates that are stil left valid we select best one among ones with same position and type
			ocl_select_best_candidates(cmd_queue_info, layer_candidates, layer_candidates_count, kernel_err, 1, &evnt[1], &evnt[2]);

			// we also need to create compact offset needed for sorting from remaining valid candidates 
			// (at the same time we also remove excessive candidtes if g_response_threshold_percent is valid number)
			ocl_make_compact_offsets(cmd_queue_info, layer_candidates, layer_candidates_count, candidate_non_zero_positions, candidate_non_zero_pos_count,
								  new_layer_s_nodes_offsets, new_layer_s_nodes_offsets_keys, new_layer_s_nodes_offsets_count,
								  count_new_parts_result, count_new_edges_result, count_valid_non_zero_pos_result, kernel_err, 2, &evnt[2], &evnt[4]);

		}
		// we wait until all have processed and then read all data (using events to synchronize)

		ret = clEnqueueReadBuffer(cmd_queue_info.queue, kernel_err, CL_FALSE, 0, sizeof(cl_int), &kernel_err_host, 1, &evnt[4], nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clEnqueueReadBuffer(cmd_queue_info.queue, count_new_parts_result, CL_FALSE, 0, 1 * sizeof(cl_int), &count_new_parts, 1, &evnt[4], nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clEnqueueReadBuffer(cmd_queue_info.queue, count_new_edges_result, CL_FALSE, 0, 1 * sizeof(cl_int), &count_new_edges, 1, &evnt[4], nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clEnqueueReadBuffer(cmd_queue_info.queue, count_valid_non_zero_pos_result, CL_FALSE, 0, 1 * sizeof(cl_int), &count_valid_non_zero_pos, 1, &evnt[4], nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clFinish(cmd_queue_info.queue);
		ocl_check_error(ret, CL_SUCCESS);


		if (candidate_non_zero_positions != nullptr) clReleaseMemObject(candidate_non_zero_positions);
		if (candidate_non_zero_pos_count != nullptr) clReleaseMemObject(candidate_non_zero_pos_count);
		if (new_layer_s_nodes_offsets_count != nullptr) clReleaseMemObject(new_layer_s_nodes_offsets_count);

		if (count_new_parts_result != nullptr) clReleaseMemObject(count_new_parts_result);
		if (count_new_edges_result != nullptr) clReleaseMemObject(count_new_edges_result);
		if (count_valid_non_zero_pos_result != nullptr) clReleaseMemObject(count_valid_non_zero_pos_result);		

		if (kernel_err_host > 0) {
			cout << "exception in candidates_do_matching_kernel or candidates_select_best_parts_kernel with error number '" << kernel_err_host << "' !!!!!!!!!!!!!!!!!!!!!" << endl;
			getchar();
			throw std::exception();
		}

		if (count_new_parts > 0 && count_new_edges > 0) {

			if (sizeof(ocl_edge_data_ip2) * count_new_edges > cmd_queue_info.max_mem_alocation) {
				cout << endl << "WARNING: To many edges created (" << count_new_edges << ") for OpenCL to be able to handle them. Adjust config values otherwise OpenCL may not be able to compute !!" << endl;
			}
			if (sizeof(ocl_layer1_data) * count_new_parts > cmd_queue_info.max_mem_alocation) {
				cout << endl << "WARNING: To many parts created (" << count_new_parts << ") for OpenCL to be able to handle them. Adjust config values otherwise OpenCL may not be able to compute  !!" << endl;
			}

			// now that we know size of new layer (number of parts, edges and offsets of where to save them) we can create all memory for new layer		

			cl_mem new_layer_img_s_nodes = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_layer1_data) * count_new_parts, nullptr, &ret);
			ocl_check_error(ret, CL_SUCCESS);

			cl_mem new_layer_img_edges = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_edge_data_ip2) * count_new_edges, nullptr, &ret);
			ocl_check_error(ret, CL_SUCCESS);

			cl_mem merge_sort_temp_data = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_layer1_data) * count_new_parts, nullptr, &ret);
			ocl_check_error(ret, CL_SUCCESS);

			cl_event before_sort_evnt[2];
			cl_uint num_events_to_wait_before_sort = 0;
			// in case adjusted_layer1_coords_size is not big enough for bitonic sort of array of length number_positions_used_result_host 
			// then we need to create new buffer and copy old values to it
			int bitonic_sort_buffer_size = ::pow(2,::ceil(::log((float)count_valid_non_zero_pos)/::log((float)2)));

			if (bitonic_sort_buffer_size > adjusted_layer_candidates_count) {
				// not big enough buffer; we need to create new buffers
				cl_mem adjusted_new_layer_s_nodes_offsets_keys = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_coord_sorting_pointer) * bitonic_sort_buffer_size, nullptr, &ret);
				ocl_check_error(ret, CL_SUCCESS);

				cl_mem addjusted_new_layer_s_nodes_offsets = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(ocl_layer1_data_coordinates) * bitonic_sort_buffer_size, nullptr, &ret);
				ocl_check_error(ret, CL_SUCCESS);
				
				// set memory to zero (only new_keys need to be set to zero)
				cl_event new_keys_set_zero_evnt;
				cl_int in_all_bytes = (bitonic_sort_buffer_size * sizeof(ocl_coord_sorting_pointer));
				cl_int in_unalinged_bytes = in_all_bytes - (in_all_bytes / sizeof(cl_int16)) * sizeof(cl_int16);

				// using ocl_set_memory_zero for whole new memory might be unoptimized since in line below we will copy most of it from original buffer
				// but since this if case is extremly rare it should not cause serioius performance issues
				ocl_set_memory_zero(cmd_queue_info, adjusted_new_layer_s_nodes_offsets_keys, in_all_bytes / sizeof(cl_int16), in_unalinged_bytes / sizeof(cl_int), 0, nullptr, &new_keys_set_zero_evnt);			

				// then we need to copy values from old buffer to new one
				ret = clEnqueueCopyBuffer(cmd_queue_info.queue, new_layer_s_nodes_offsets_keys, adjusted_new_layer_s_nodes_offsets_keys, 0, 0, sizeof(ocl_coord_sorting_pointer) * count_valid_non_zero_pos, 1, &new_keys_set_zero_evnt, &before_sort_evnt[0]);
				ret |= clEnqueueCopyBuffer(cmd_queue_info.queue, new_layer_s_nodes_offsets, addjusted_new_layer_s_nodes_offsets, 0, 0, sizeof(ocl_layer1_data_coordinates) * count_valid_non_zero_pos, 0, nullptr, &before_sort_evnt[1]);
				ocl_check_error(ret, CL_SUCCESS);

				// release old buffers
				clReleaseMemObject(new_layer_s_nodes_offsets);
				clReleaseMemObject(new_layer_s_nodes_offsets_keys);

				// and reset old buffers to point to new one
				new_layer_s_nodes_offsets = addjusted_new_layer_s_nodes_offsets;
				new_layer_s_nodes_offsets_keys = adjusted_new_layer_s_nodes_offsets_keys;

				num_events_to_wait_before_sort = 2;				
				//cout << "\t\tadjusting size of array to " << bitonic_sort_buffer_size << endl;
			}
			
			// first sort shape_nodes offsets (i.e. compact coordinates) and update theirs offset values to define memory position for each position
			int num_events_to_wait_before_update_offset = 0;
			if (count_valid_non_zero_pos > 1) {
				//cout << "\t\tsize of count_valid_non_zero_pos : " << count_valid_non_zero_pos << " and size of adjusted_layer_candidates_count: " << adjusted_layer_candidates_count << endl; 
				ocl_sort_compact_coords(cmd_queue_info, count_valid_non_zero_pos, new_layer_s_nodes_offsets, new_layer_s_nodes_offsets_keys, kernel_err, num_events_to_wait_before_sort, num_events_to_wait_before_sort > 0 ? before_sort_evnt : nullptr, &evnt[0]);
				num_events_to_wait_before_update_offset = 1;
			}

			//ocl_update_results_offsets(cmd_queue_info, count_valid_non_zero_pos, new_layer_s_nodes_offsets, kernel_err, 0, nullptr, &evnt[1]);
			ocl_update_results_offsets(cmd_queue_info, count_valid_non_zero_pos, new_layer_s_nodes_offsets, kernel_err, num_events_to_wait_before_update_offset, num_events_to_wait_before_update_offset > 0 ? &evnt[0] : nullptr, &evnt[1]);

			// make resuls (parts) from valid candidtes
			ocl_make_results(cmd_queue_info, lib_parts_layer_k, lib_edges_layer_k, lib_apps_layer_k, lib_parts_layer_k1, img_s_nodes_layer_k1, img_coord_layer_k1, img_edges_layer_k1, distr_gaussian_mask,
							 new_layer_s_nodes_offsets, new_layer_s_nodes_offsets_keys, count_valid_non_zero_pos, layer_candidates, layer_candidates_count, max_schur_memory_obj, img_size_layer_k1, dummy_m_size, k1,
							 new_layer_img_s_nodes, new_layer_img_edges, new_layer_img_size, kernel_err, 1, &evnt[1], &evnt[2]);
			
			// sort parts within each position
			ocl_sort_within_positions(cmd_queue_info, new_layer_img_s_nodes, count_valid_non_zero_pos, new_layer_s_nodes_offsets, new_layer_img_size, merge_sort_temp_data, kernel_err, 1, &evnt[2], &evnt[3]);

			// when merge sort will be done we can read final data to new host memory (using events for synchronization)

			ocl_layer1_data_coordinates* new_layer_s_nodes_offsets_host = new ocl_layer1_data_coordinates[count_valid_non_zero_pos];

			new_layer_img_s_nodes_host = new ocl_layer1_data[count_new_parts];
			new_layer_img_edges_host = new ocl_edge_data_ip2[count_new_edges];			

			ret = clEnqueueReadBuffer(cmd_queue_info.queue, new_layer_s_nodes_offsets, CL_FALSE, 0, sizeof(ocl_layer1_data_coordinates) * count_valid_non_zero_pos, new_layer_s_nodes_offsets_host, 1, &evnt[3], nullptr);
			ocl_check_error(ret, CL_SUCCESS);

			ret = clEnqueueReadBuffer(cmd_queue_info.queue, new_layer_img_s_nodes, CL_FALSE, 0, sizeof(ocl_layer1_data) * count_new_parts, new_layer_img_s_nodes_host, 1, &evnt[3], nullptr);
			ocl_check_error(ret, CL_SUCCESS);

			ret = clEnqueueReadBuffer(cmd_queue_info.queue, new_layer_img_edges, CL_FALSE, 0, sizeof(ocl_edge_data_ip2) * count_new_edges, new_layer_img_edges_host, 1, &evnt[3], nullptr);
			ocl_check_error(ret, CL_SUCCESS);

			ret = clEnqueueReadBuffer(cmd_queue_info.queue, kernel_err, CL_FALSE, 0, sizeof(cl_int), &kernel_err_host, 1, &evnt[3], nullptr);
			ocl_check_error(ret, CL_SUCCESS);

			ret = clFinish(cmd_queue_info.queue);
			ocl_check_error(ret, CL_SUCCESS);

			if (new_layer_img_s_nodes != nullptr) clReleaseMemObject(new_layer_img_s_nodes);
			if (new_layer_img_edges != nullptr) clReleaseMemObject(new_layer_img_edges);
			if (merge_sort_temp_data != nullptr) clReleaseMemObject(merge_sort_temp_data);

			if (kernel_err_host > 0) {
				cout << "exception in candidates_make_results_kernel or merge_sort_within_pos_result_kernel with error number '" << kernel_err_host << "' !!!!!!!!!!!!!!!!!!!!!" << endl;
				getchar();
				kernel_err_host = 0;
				clEnqueueWriteBuffer(cmd_queue_info.queue, kernel_err, CL_TRUE, 0, sizeof(cl_int), &kernel_err_host, 0, nullptr, nullptr);
			}
		

			// before we finish we also need to construct array (img coordinates) of inhibited shape nodes)
			/*{
				cl_kernel inhibit_result_kernel = OpenCL::context_manager->getCommonKernel("inhibit_result")[context_num];

				//memset(new_layer_img_coord_inhib_host, 0, sizeof(ocl_layer1_data_coordinates) * xy_size_k);

				new_layer_img_coord_inhib = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(ocl_layer1_data_coordinates) * xy_size_k, nullptr, &ret);
				ocl_check_error(ret, CL_SUCCESS);

				ret = clSetKernelArg(inhibit_result_kernel, 0, sizeof(cl_mem), (void*)&new_layer_img_s_nodes);
				ret |= clSetKernelArg(inhibit_result_kernel, 1, sizeof(cl_mem), (void*)&new_layer_img_coord);		
				ret |= clSetKernelArg(inhibit_result_kernel, 2, sizeof(cl_mem), (void*)&new_layer_img_coord_inhib);	

				ret |= clSetKernelArg(inhibit_result_kernel, 3, sizeof(cl_uint2), (void*)&new_layer_img_size);		
				ret |= clSetKernelArg(inhibit_result_kernel, 4, sizeof(cl_mem), (void*)&kernel_err);
				ocl_check_error(ret, CL_SUCCESS);

				ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, inhibit_result_kernel, 1, nullptr, global_work_load_by_img , nullptr, 0, nullptr, &evnt[0]); 
				ocl_check_error(ret, CL_SUCCESS);

				ret = clEnqueueReadBuffer(cmd_queue_info.queue, new_layer_img_coord_inhib, CL_FALSE, 0, sizeof(ocl_layer1_data_coordinates) * xy_size_k, new_layer_img_coord_inhib_host, 1, &evnt[0], nullptr);
				ocl_check_error(ret, CL_SUCCESS);

				ret = clEnqueueReadBuffer(cmd_queue_info.queue, kernel_err, CL_FALSE, 0, sizeof(cl_int), &kernel_err_host, 1, &evnt[0], nullptr);
				ocl_check_error(ret, CL_SUCCESS);

				ret = clFinish(cmd_queue_info.queue);
				ocl_check_error(ret, CL_SUCCESS);

				if (kernel_err_host > 0) {
					cout << "exception in inhibit_result_kernel with error number '" << kernel_err_host << "' !!!!!!!!!!!!!!!!!!!!!" << endl;
					getchar();
					throw std::exception();
				}
			}*/
			
		
			// we have only created compact form of layer1_coordinates so we need to copy it to full new_layer_img_coord_host where each position 
			// implicitly represets one image position		
			for (int i = 0; i < count_valid_non_zero_pos; i++) {
				// get position for coordinates
				int offset = new_layer_s_nodes_offsets_host[i].offset;
				int pos = new_layer_img_s_nodes_host[offset].y * x_size_k + new_layer_img_s_nodes_host[offset].x;

				// copy values
				new_layer_img_coord_host[pos] = new_layer_s_nodes_offsets_host[i];
			}

			delete[] new_layer_s_nodes_offsets_host;
		}
		
		// release library memory objects
		if (lib_parts_layer_k != nullptr) clReleaseMemObject(lib_parts_layer_k); 
		if (lib_edges_layer_k != nullptr) clReleaseMemObject(lib_edges_layer_k);

		if (lib_apps_layer_k != nullptr) clReleaseMemObject(lib_apps_layer_k);
		
		if (lib_parts_layer_k1 != nullptr) clReleaseMemObject(lib_parts_layer_k1);
		if (lib_edges_layer_k1 != nullptr) clReleaseMemObject(lib_edges_layer_k1);

		// release memory objects of previous layer
		if (img_s_nodes_layer_k1 != nullptr) clReleaseMemObject(img_s_nodes_layer_k1);
		if (img_edges_layer_k1 != nullptr) clReleaseMemObject(img_edges_layer_k1);
		if (img_coord_layer_k1 != nullptr) clReleaseMemObject(img_coord_layer_k1);
		
		// release memory objects of new layer
		if (new_layer_s_nodes_offsets != nullptr) clReleaseMemObject(new_layer_s_nodes_offsets);
		if (new_layer_s_nodes_offsets_keys != nullptr) clReleaseMemObject(new_layer_s_nodes_offsets_keys);
		if (layer_candidates != nullptr) clReleaseMemObject(layer_candidates);
		if (max_schur_memory_obj != nullptr) clReleaseMemObject(max_schur_memory_obj);
		
		// release other memory objects
		if (distr_gaussian_mask != nullptr) clReleaseMemObject(distr_gaussian_mask);
		if (kernel_err != nullptr) clReleaseMemObject(kernel_err); 

	}
	
	// count number of non-zero values in new_layer_img_coord_host and new_layer_img_coord_inhib_host (last one is only for double checking)
	int count_non_zero_position = 0;
	int count_non_zero_position_inhib = 0;

	if (count_new_parts > 0) {
		for (int i = 0; i < xy_size_k; i++) {
			count_non_zero_position += new_layer_img_coord_host[i].size > 0 ? 1 : 0;
			count_non_zero_position_inhib += new_layer_img_coord_inhib_host[i].size > 0 ? 1 : 0;
		}
	}

	// now save all parts, edges, coordinates and sizes to layer_1_result
	res->ocl_shape_nodes[k].first = new_layer_img_s_nodes_host;
	res->ocl_shape_nodes[k].second = count_new_parts;

	res->ocl_edges[k].first = new_layer_img_edges_host;
	res->ocl_edges[k].second = count_new_edges;

	res->ocl_shape_nodes_coord[k].first = new_layer_img_coord_host;
	res->ocl_shape_nodes_coord[k].second = xy_size_k;

	res->ocl_shape_nodes_inhib_coord[k].first = new_layer_img_coord_inhib_host;
	res->ocl_shape_nodes_inhib_coord[k].second = xy_size_k;

	res->ocl_shape_nodes_coord_non_zero_count[k] = count_non_zero_position;
	res->ocl_shape_nodes_inhib_coord_non_zero_count[k] = count_non_zero_position_inhib;

	// find out how many parts in layer k - 1 were finaly used as center for new part in layer k
	int count_loop = 0;
	for (int i = 0 ; i < layer_candidates_count; i++) {
		if (layer_candidates_host[i].skip == 0) {
			count_loop++;
			unsigned int* attr = &(img_s_nodes_layer_k1_host[layer_candidates_host[i].img_part_offset].attr);
			if (ocl_is_attr_set(*attr, HAS_NEXT_LAYER) == false) {
				ocl_set_attr(attr, HAS_NEXT_LAYER);
				res->inc_covered(k1);
			}
		}
	}

	if (layer_candidates_host != nullptr) delete[] layer_candidates_host;

	double quot = res->ocl_cover_quotient(k1);

	count++;
	quot_sum += quot;
	if (quot < min_cover_quot)
		min_cover_quot = quot;

#ifdef OPENCL_VERIFICATION_SUPPORT
	layer1_result* org_res = nullptr;
	
	if (opencl_verify_result == true && count_new_parts > 0){
		org_res = new layer1_result();

		cloner cl;
		res->copy_to(org_res, cl);
	}
#endif

	res->make_data_from_ocl(k, add_edge_names);

#ifdef OCL_PROFILE_HOST_QUICK
clock_1.stopTimer();
cout << "done in : " << clock_1.getElapsedTime() << endl;
#endif

#ifdef OPENCL_VERIFICATION_SUPPORT
// in case binary was build with support for verification of opencl result then use config settings to indicate wheaterh to make verification
	if (opencl_verify_result == true) {

		if (count_new_parts > 0) {
			
				
			int invalid_sort = 0;
			for (int j = 0; j < res->y_size(k); j++) {
				for (int i = 0; i < res->x_size(k); i++) {

					node* ocl_next = res->node_at(i,j,k);
					float before = FLT_MAX;
					while (ocl_next != nullptr) {

						layer1_data* ocl_result_part = (layer1_data*)(ocl_next->data);
						
						if (before < ocl_result_part->r(G_RESPONSE) ) {
							invalid_sort++;
						}

						before = ocl_result_part->r(G_RESPONSE);

						ocl_next = ((layer1_data*)ocl_next->data)->next;
					}
				}
			}
			cout << endl << "number invalid sorted opencl results : " << invalid_sort << endl;


			add_layer1(org_res, layer);

			vector<node*>& org_s_nodes = org_res->shape_nodes[k];
			vector<node*>& ocl_s_nodes = res->shape_nodes[k];

			int invalid_results = 0;
			int all_nodes = 0;
			for (int j = 0; j < res->y_size(k); j++) {
				for (int i = 0; i < res->x_size(k); i++) {

					node* org_next__ = org_res->node_at(i,j,k);
					node* ocl_next__ = res->node_at(i,j,k);
					node* org_next = org_next__;
					node* ocl_next = ocl_next__;

					map<int,node*> org_map, ocl_map;

					while (ocl_next != nullptr) {						
						ocl_map[((layer1_data*)ocl_next->data)->m] = ocl_next;
						ocl_next = ((layer1_data*)ocl_next->data)->next;
					}
						
					while (org_next != nullptr) {
						org_map[((layer1_data*)org_next->data)->m] = org_next;
						org_next = ((layer1_data*)org_next->data)->next;
					}

					if (org_map.size() != ocl_map.size()) {
						invalid_results++;
						cout << "invalid number of parts in one location (org size: " << org_map.size() << ", ocl size: " << ocl_map.size() << ")" << endl;
					}

					for (map<int,node*>::iterator iter = org_map.begin(); iter != org_map.end(); iter++) {

						int type = iter->first;
						org_next = iter->second;

						all_nodes++;
						layer1_data* org_part = (layer1_data*)(org_next->data);					

						// find same type in ocl
						map<int,node*>::iterator ocl_found = ocl_map.find(type);						

						if (ocl_found == ocl_map.end()) {
							cout << "\t unable to find part with type " << type << " in ocl at location x: " << org_part->x << " y: " << org_part->y << endl;
							continue;
						}

						ocl_next = ocl_found->second;
						layer1_data* ocl_result_part = (layer1_data*)(ocl_next->data);

						if (org_next->neighbors.size() != ocl_next->neighbors.size() && false) {
							invalid_results++;
						} else {
							// for each neighboor in org_next find match in neighboor in ocl_next
							forall_neighbors(org_next, n_iter) {
								int x = ((layer1_data*)n_iter->second.first->data)->x;
								int y = ((layer1_data*)n_iter->second.first->data)->y;
								int m = ((layer1_data*)n_iter->second.first->data)->m;
								int z = ((layer1_data*)n_iter->second.first->data)->z;
								response_map r = ((layer1_data*)n_iter->second.first->data)->r;
								
								//if (n_iter->first == 39)
									//continue;

								bool found = false;
								foreach_neighbor(ocl_next, n_iter->first, m_iter) {
									if (((layer1_data*)m_iter->second.first->data)->x == x && 
										((layer1_data*)m_iter->second.first->data)->y == y && 
										((layer1_data*)m_iter->second.first->data)->m == m && 
										((layer1_data*)m_iter->second.first->data)->z == z &&
										fabs(((layer1_data*)m_iter->second.first->data)->r(R_RESPONSE) - r(R_RESPONSE)) <= 1e-3  && 
										fabs(((layer1_data*)m_iter->second.first->data)->r(G_RESPONSE) - r(G_RESPONSE)) <= 1e-3 ) {
										found = true;
										break;
									}									
								}
								if (found == false) {
									invalid_results++;
								}
							}
						}

						
						
						if (ocl_result_part->x != org_part->x) {
							cout << "diff x in part at: " << ocl_result_part->x << " y: " << ocl_result_part->y << 
                                "and correct part at x: " << org_part->x << " y: " << org_part->y << endl;
							invalid_results++;
						}

						if (ocl_result_part->y != org_part->y) {
							cout << "diff y in part at: " << ocl_result_part->x << " y: " << ocl_result_part->y << 
                                "and correct part at x: " << org_part->x << " y: " << org_part->y << endl;
							invalid_results++;
						}
						if (ocl_result_part->z != org_part->z) {
							cout << "diff z in part at: " << ocl_result_part->x << " y: " << ocl_result_part->y << 
                                "and correct part at x: " << org_part->x << " y: " << org_part->y << endl;
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
							cout << "diff at x: " << ocl_result_part->x << " y: " << ocl_result_part->y << " m: " << ocl_result_part->m << " in RR_response: " << fabs(ocl_result_part->r(RR_RESPONSE) - org_part->r(RR_RESPONSE)) << endl;
							invalid_results++;
						}
						//ocl_next = ((layer1_data*)ocl_next->data)->next;
						//org_next = ((layer1_data*)org_next->data)->next;

					}
				}
			}
			cout << "number invalid results : " << invalid_results << endl;
	/*
			
			{
				static int counter_for_test = 0;
				stringstream ss;
				ss << "results_ocl_" << counter_for_test++ << ".txt";
				cout << "writing to file '" << ss.str() << "'" << endl;
				FILE* fp = fopen(ss.str().c_str(),"w+");

				for (int i = 0; i < res->y_size(k); i++) {
					for (int m = 0; m < res->x_size(k); m++) {		
						node* n = res->node_at(m,i, k);

						while (n != nullptr) {
						
							layer1_data* nd = (layer1_data*)n->data;

							if (nd != nullptr) {

								fprintf(fp, "x: %d, y: %d, z: %d, m: %d \n", nd->x, nd->y, nd->z , nd->m);
								fprintf(fp, "\tR: %.3f, G: %.3f, RR: %.3f\n\n", nd->r.get_response(R_RESPONSE,0), nd->r.get_response(G_RESPONSE,0), nd->r.get_response(RR_RESPONSE,0));
							}

							if (n->neighbors.size() > 0) {

								fprintf(fp, "neighboors: \n");

								forall_neighbors(n, n_iter) {
									int type = n_iter->first;
									layer1_data* p_n = (layer1_data*)n_iter->second.first->data;
									
									fprintf(fp, "\ttype: %d,x: %d, y: %d, z: %d, m: %d, R: %.3f, attr: %d\n", type, p_n->x, nd->y, p_n->z , p_n->m, p_n->r.get_response(R_RESPONSE,0), n_iter->second.first->attr);
								}
							}

							n = nd->next;
						}
					}
				}
				fflush(fp);
				fclose(fp);
			}
			{
				static int counter_for_test = 0;
				stringstream ss;
				ss << "results_org_" << counter_for_test++ << ".txt";
				cout << "writing to file '" << ss.str() << "'" << endl;
				FILE* fp = fopen(ss.str().c_str(),"w+");

				for (int i = 0; i < org_res->y_size(k); i++) {
					for (int m = 0; m < org_res->x_size(k); m++) {		
						node* n = org_res->node_at(m,i, k);

						while (n != nullptr) {
						
							layer1_data* nd = (layer1_data*)n->data;

							if (nd != nullptr) {

								fprintf(fp, "x: %d, y: %d, z: %d, m: %d \n", nd->x, nd->y, nd->z , nd->m);
								fprintf(fp, "\tR: %.3f, G: %.3f, RR: %.3f\n\n", nd->r.get_response(R_RESPONSE,0), nd->r.get_response(G_RESPONSE,0), nd->r.get_response(RR_RESPONSE,0));
							}

							if (n->neighbors.size() > 0) {

								fprintf(fp, "neighboors: \n");

								for (multimap<int, edge_pair>::reverse_iterator n_iter = (n)->neighbors.rbegin();
									n_iter != (n)->neighbors.rend(); ++n_iter) {
									int type = n_iter->first;
									layer1_data* p_n = (layer1_data*)n_iter->second.first->data;
									
									fprintf(fp, "\ttype: %d,x: %d, y: %d, z: %d, m: %d, R: %.3f, attr: %d\n", type, p_n->x, nd->y, p_n->z , p_n->m, p_n->r.get_response(R_RESPONSE,0), n_iter->second.first->attr);
								}
							}

							n = nd->next;
						}
					}
				}
				fflush(fp);
				fclose(fp);
			}
	*/
			
		} else {

		}
		if (org_res != nullptr)
			delete org_res;
	}
#endif
}


void layern_creator::ocl_set_candidate_thresholds(ocl_layer1_data* img_s_nodes_layer_k1_host, const int ocl_shape_nodes_count, const int k1) {

    if (ocl_shape_nodes_count <= 0) 
        return;

	if (!_isnan(candidate_r_threshold_percent)) {
		float max_r = 0;

		// first layer is sorted by R_RESPONSE so we can easly get only first one
        if (k1 == 0)
			max_r = img_s_nodes_layer_k1_host[0].response.R;
        else {
			// other layers are sorted by G_RESPONSE so we need to find best R_RESPONSE by checking all
		    for (int i = 0; i < ocl_shape_nodes_count; i++) {
			    float resp = img_s_nodes_layer_k1_host[i].response.R;
			    if (resp > max_r) max_r = resp;
		    }		
        }
		candidate_r_threshold = max_r * candidate_r_threshold_percent;
    }
    if (!_isnan(candidate_g_threshold_percent)) {
		double max_r = 0;

		// all but first layers are sorted by G_RESPONSE so only get first one in line
        if (k1 == 0) max_r = 1.0;
        else max_r = img_s_nodes_layer_k1_host[0].response.G;
		candidate_g_threshold = max_r * candidate_g_threshold_percent;
    }
}

ocl_layer_candidate* layern_creator::ocl_generate_candidates(ocl_part_data* lib_parts_layer_k1_host,
											ocl_part_data_2* lib_edges_layer_k1_host,
											ocl_part_data* lib_parts_layer_k_host,
											ocl_layer1_data* img_s_nodes_layer_k1_host,
											const int ocl_shape_nodes_count, const int x_size_k, const int y_size_k, const int k,
											int &layer_candidates_count, int &max_schur_count) { 

	
	const int xy_size_k = x_size_k * y_size_k;
	int count_skiped = 0;
	
	ocl_layer1_data_coordinates* new_layer_img_coord_host = new ocl_layer1_data_coordinates[xy_size_k];
	memset(new_layer_img_coord_host, 0, sizeof(ocl_layer1_data_coordinates) * xy_size_k);

	// count how many candidates we will have (we need to count them first since we will be creating array of fixed size)
	layer_candidates_count = 0;	

	cl_float2 candidate_threshold;
	for (int i = 0; i < ocl_shape_nodes_count; i++) {
		// CATION: when changing this loop make sure you also change lower loop where candidates are filled

		ocl_layer1_data* part = &img_s_nodes_layer_k1_host[i];
	
		// skip part if not within threshold value
		if (part->response.R < candidate_r_threshold || part->response.G < candidate_g_threshold ||
			ocl_is_attr_set(part->attr, HAS_NEXT_LAYER) /*|| n->is_attr_set(NEW_ATTR)*/ || (ignore_texture && ocl_is_attr_set(part->attr, TEXTURE_NODE))) {
			count_skiped++;
			continue;
		}
	
		int x_part_loc = int_round(part->x/layer_contraction); // !floor!
		int y_part_loc = int_round(part->y/layer_contraction); // !floor!
        
		// New part coordinates must be within this bounds to be used..
		if (x_part_loc <= 0 || x_part_loc + 1 >= x_size_k || y_part_loc <= 0 || y_part_loc + 1 >= y_size_k) {
			continue;
		}

		// get part type
		int part_type = part->m;

		// count how many candidates we can get for this part
		//layer_candidates_count += lib_parts_layer_k1_host[part_type].edges_loc[OCL_LYR_CENTER_BACK_EDGE].size;
		int offset_to_center_parts = lib_parts_layer_k1_host[part_type].edges_loc[OCL_LYR_CENTER_BACK_EDGE].offset;
		int offset_end_to_center_parts = offset_to_center_parts + lib_parts_layer_k1_host[part_type].edges_loc[OCL_LYR_CENTER_BACK_EDGE].size;

		for (int j = offset_to_center_parts; j < offset_end_to_center_parts; j++) {

			ocl_part_data_2* edge = &lib_edges_layer_k1_host[j];

			if (edge->node.layer != k) {
				cout << "not supported to other layers" << endl;
				throw std::exception();
			}
			// Check whether this part is allowed to be added; allowed_parts is a parameter.
			// Empty value means no restriction.
			int part_type = lib_parts_layer_k_host[edge->node.offset].type;

			if (!allowed_parts.empty() && allowed_parts.find(part_type) == allowed_parts.end())
					continue;

			int cmx = lib_parts_layer_k_host[edge->node.offset].cmx;
			int cmy = lib_parts_layer_k_host[edge->node.offset].cmy;
		
			int new_x = int_round((part->x + cmx)/layer_contraction); // !floor!
			int new_y = int_round((part->y + cmy)/layer_contraction); // !floor!

			int pos = int_round((part->y + cmy)/layer_contraction) * x_size_k + int_round((part->x + cmx)/layer_contraction) ;
			
			// count how many candidates we have for each position
			new_layer_img_coord_host[pos].size++;
			layer_candidates_count++;
		}
	}

	ocl_layer_candidate* layer_candidates_host = nullptr;

	// if we have any candidates we need to create them in buffer
	if (layer_candidates_count > 0) {


		// set initial offsets for candidates of new layer img coordinates that will
		// be used only during construction (actual offsets of final parts are set later)
		int counter = 0;
		for (int i = 0 ; i < xy_size_k;i++) {
			if (new_layer_img_coord_host[i].size > 0) {
				new_layer_img_coord_host[i].offset = counter;
				counter += new_layer_img_coord_host[i].size;
				// set size back to 0 since it was only used as temp
				new_layer_img_coord_host[i].size = 0;
			}
		}

		// make place for all candidates for this layer
		layer_candidates_host = new ocl_layer_candidate[layer_candidates_count];
		
		// then fill memory with candidates ( MUST be same procedure for creating candidates as above)
		max_schur_count = convolution_link_threshold > 0 ? 0 : 1 ; // use size one for dummy cl_mem obejct when no convolution_link_threshold is set
		layer_candidates_count = 0;
		for (int i = 0; i < ocl_shape_nodes_count; i++) {

			ocl_layer1_data* part = &img_s_nodes_layer_k1_host[i];

			// skip part if not within threshold value
			if (part->response.R < candidate_r_threshold || part->response.G < candidate_g_threshold ||
				ocl_is_attr_set(part->attr, HAS_NEXT_LAYER) /*|| n->is_attr_set(NEW_ATTR)*/ || (ignore_texture && ocl_is_attr_set(part->attr, TEXTURE_NODE))){
				continue;
			}
		
			int x_part_loc = int_round(part->x/layer_contraction); // !floor!
			int y_part_loc = int_round(part->y/layer_contraction); // !floor!
	        
			// New part coordinates must be within this bounds to be used..
			if (x_part_loc <= 0 || x_part_loc + 1 >= x_size_k || y_part_loc <= 0 || y_part_loc + 1 >= y_size_k) {
				continue;
			}

			// get part type
			int img_part_type = part->m;
			// get actual part from lib

			int offset_to_center_parts = lib_parts_layer_k1_host[img_part_type].edges_loc[OCL_LYR_CENTER_BACK_EDGE].offset;
			int offset_end_to_center_parts = offset_to_center_parts + lib_parts_layer_k1_host[img_part_type].edges_loc[OCL_LYR_CENTER_BACK_EDGE].size;


			for (int j = offset_to_center_parts; j < offset_end_to_center_parts; j++) {

				ocl_part_data_2* edge = &lib_edges_layer_k1_host[j];

				if (edge->node.layer != k) {
					cout << "not supported to other layers" << endl;
					throw std::exception();
				}
				int part_type = lib_parts_layer_k_host[edge->node.offset].type;

				if (!allowed_parts.empty() && allowed_parts.find(part_type) == allowed_parts.end())
					continue;

				int new_x = int_round((part->x + lib_parts_layer_k_host[edge->node.offset].cmx)/layer_contraction); // !floor!
				int new_y = int_round((part->y + lib_parts_layer_k_host[edge->node.offset].cmy)/layer_contraction); // !floor!

				int coord_pos = new_y * x_size_k + new_x;

				int candidate_pos = new_layer_img_coord_host[coord_pos].offset + new_layer_img_coord_host[coord_pos].size++; // we use size as temp offset

				// we only need to set x,y and part_offset other values are already zeroed by memset
				layer_candidates_host[candidate_pos].img_part_type = img_part_type; // cpartd in original add_layer1()
				layer_candidates_host[candidate_pos].part_type = part_type; // pcd in in original add_layer1()
				layer_candidates_host[candidate_pos].x = part->x;				
				layer_candidates_host[candidate_pos].y = part->y;
				layer_candidates_host[candidate_pos].new_x = new_x;
				layer_candidates_host[candidate_pos].new_y = new_y;
				layer_candidates_host[candidate_pos].lib_part_offset = edge->node.offset;
				layer_candidates_host[candidate_pos].img_part_offset = i;
				layer_candidates_host[candidate_pos].result_edges_prev_lay_size = 0;
				layer_candidates_host[candidate_pos].result_edges_layer0_size = 0;
				layer_candidates_host[candidate_pos].skip = 0;
				layer_candidates_host[candidate_pos].r_sum = 0;
				layer_candidates_host[candidate_pos].realization_ratio = 0;
				layer_candidates_host[candidate_pos].g_prod = 0;
				layer_candidates_host[candidate_pos].max_schur_mem_offset = max_schur_count;
	
				// recount all candidates
				layer_candidates_count++;
				// count number of times max_schur_product will be called (we need to reserve one float for each call)
				if (convolution_link_threshold > 0) 
					max_schur_count += lib_parts_layer_k_host[edge->node.offset].edges_loc[OCL_LYR_SRC_EDGE].size + 1;
			}
		}
	}

	delete[] new_layer_img_coord_host;

	return layer_candidates_host;
}

void layern_creator::ocl_eliminate_forbidden(const best_cmd_queue_info& cmd_queue_info, 
											 cl_mem lib_parts_layer_k, cl_mem lib_edges_layer_k, cl_mem lib_parts_layer_k1, cl_mem img_s_nodes_layer_k1, cl_mem img_coord_layer_k1, 
											 cl_mem layer_candidates, cl_int layer_candidates_count, cl_uint2 img_size_layer_k1, cl_int k,
											 cl_mem kernel_err, cl_int wait_event_count, cl_event* wait_evnt, cl_event* evnt) {
	
	int ret;
	cl_kernel candidates_eliminate_forbidden_kernel = OpenCL::context_manager->getCommonKernel("candidates_eliminate_forbidden")[cmd_queue_info.context_number];

	cl_float c = 1.0f;
	cl_float convolution_threshold_float = (float)convolution_threshold;
	cl_float proj_max_forb_threshold_float = (float)proj_max_forb_threshold;
	cl_float forb_quot_threshold_float = (float)forb_quot_threshold;

	// Now enqueue kernel for processing forbidden parts.
	// each candidate will be processed in its own work-item

	// set arguments	
	ret = clSetKernelArg(candidates_eliminate_forbidden_kernel, 0, sizeof(cl_mem), (void*)&lib_parts_layer_k);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 1, sizeof(cl_mem), (void*)&lib_edges_layer_k);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 2, sizeof(cl_mem), (void*)&lib_parts_layer_k1);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 3, sizeof(cl_mem), (void*)&img_s_nodes_layer_k1);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 4, sizeof(cl_mem), (void*)&img_coord_layer_k1);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 5, sizeof(cl_mem), (void*)&layer_candidates);		
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 6, sizeof(cl_int), (void*)&layer_candidates_count);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 7, sizeof(cl_uint2), (void*)&img_size_layer_k1);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 8, sizeof(cl_float), (void*)&c);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 9, sizeof(cl_int), (void*)&k);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 10, sizeof(cl_float), (void*)&convolution_threshold_float);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 11, sizeof(cl_float), (void*)&proj_max_forb_threshold_float);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 12, sizeof(cl_float), (void*)&forb_quot_threshold_float);
	ret |= clSetKernelArg(candidates_eliminate_forbidden_kernel, 13, sizeof(cl_mem), (void*)&kernel_err);	
	ocl_check_error(ret, CL_SUCCESS);

	// calculate best work-gorup size based on number of candidates
	int best_local_work_size = cmd_queue_info.max_work_group_size;
	int best_group_size = layer_candidates_count / best_local_work_size + 1;
	int best_global_work_size = best_group_size * best_local_work_size;

	size_t global_work_load[1] = {best_global_work_size };
	size_t local_work_load[1] = {best_local_work_size};

	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, candidates_eliminate_forbidden_kernel, 1, nullptr, global_work_load , local_work_load, wait_event_count, wait_evnt, evnt); 
	ocl_check_error(ret, CL_SUCCESS);

	ret = clFinish(cmd_queue_info.queue);
	ocl_check_error(ret, CL_SUCCESS);
	
}

void layern_creator::ocl_candidates_matching(const best_cmd_queue_info& cmd_queue_info, 
											 cl_mem lib_parts_layer_k, cl_mem lib_edges_layer_k, cl_mem lib_apps_layer_k,
											 cl_mem lib_parts_layer_k1, cl_mem img_s_nodes_layer_k1, 
											 cl_mem img_coord_layer_k1, cl_mem img_edges_layer_k1, cl_mem distr_gaussian_mask,
											 cl_mem layer_candidates, cl_int layer_candidates_count, cl_mem max_schur_memory_obj, cl_uint2 img_size_layer_k1, cl_uint2 dummy_m_size, cl_int k1,
											 cl_mem non_zero_positions, cl_mem non_zero_pos_count, cl_mem kernel_err, cl_int wait_event_count, cl_event* wait_evnt, cl_event* evnt) {
	
	int ret;
	cl_kernel candidates_do_matching_kernel = OpenCL::context_manager->getCommonKernel("candidates_do_matching")[cmd_queue_info.context_number];

	cl_float4 response_thershold;
	response_thershold.s[0] = r_response_threshold;
	response_thershold.s[1] = realization_ratio_threshold;
	response_thershold.s[2] = g_response_threshold;

	cl_float c = 1.0f;
	cl_float g_response_var_factor_float = (float)g_response_var_factor;
	cl_float convolution_threshold_float = (float)convolution_threshold;
	cl_float convolution_link_threshold_float = (float)convolution_link_threshold;

	cl_float r_response_pow_float = (float)r_response_pow;
	cl_float g_response_pow_float = (float)g_response_pow;
	int use_manual_thresholds = manual_thresholds ? 1 : 0;

	int add_reconstruction_edges_int = add_reconstruction_edges ? 1 : 0;

	int schur_product_use_distribution = SCHUR_PROD_USE_G_DISTR_RESPONSE;
		
	if (identity_g_response)
		schur_product_use_distribution = SCHUR_PROD_USE_IDENTITY_RESPONSE;
	else if (simple_g_response)
		schur_product_use_distribution = SCHUR_PROD_USE_R_RESPONSE;
	else if (ignore_g_distribution)
		schur_product_use_distribution = SCHUR_PROD_USE_SIMPLE_G_DISTR_RESPONSE;


	// set arguments		
	ret = clSetKernelArg(candidates_do_matching_kernel, 0, sizeof(cl_mem), (void*)&lib_parts_layer_k);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 1, sizeof(cl_mem), (void*)&lib_edges_layer_k);			
	ret |= clSetKernelArg(candidates_do_matching_kernel, 2, sizeof(cl_mem), (void*)&lib_apps_layer_k);			
	ret |= clSetKernelArg(candidates_do_matching_kernel, 3, sizeof(cl_mem), (void*)&lib_parts_layer_k1);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 4, sizeof(cl_mem), (void*)&img_s_nodes_layer_k1);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 5, sizeof(cl_mem), (void*)&img_coord_layer_k1);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 6, sizeof(cl_mem), (void*)&img_edges_layer_k1);			
	ret |= clSetKernelArg(candidates_do_matching_kernel, 7, sizeof(cl_mem), (void*)&distr_gaussian_mask);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 8, sizeof(cl_mem), (void*)&layer_candidates);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 9, sizeof(cl_int), (void*)&layer_candidates_count);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 10, sizeof(cl_mem), (void*)&max_schur_memory_obj);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 11, sizeof(cl_uint2), (void*)&img_size_layer_k1);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 12, sizeof(cl_uint2), (void*)&dummy_m_size);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 13, sizeof(cl_float), (void*)&c);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 14, sizeof(cl_int), (void*)&k1);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 15, sizeof(cl_float), (void*)&convolution_threshold_float);	
	ret |= clSetKernelArg(candidates_do_matching_kernel, 16, sizeof(cl_float), (void*)&convolution_link_threshold_float);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 17, sizeof(cl_int), (void*)&g_response_var_factor_float);		
	ret |= clSetKernelArg(candidates_do_matching_kernel, 18, sizeof(cl_float4), (void*)&response_thershold);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 19, sizeof(cl_int), (void*)&use_manual_thresholds);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 20, sizeof(cl_float), (void*)&r_response_pow_float);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 21, sizeof(cl_float), (void*)&g_response_pow_float);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 22, sizeof(cl_int), (void*)&schur_product_use_distribution);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 23, sizeof(cl_int), (void*)&add_reconstruction_edges_int);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 24, sizeof(cl_mem), (void*)&non_zero_positions);
	ret |= clSetKernelArg(candidates_do_matching_kernel, 25, sizeof(cl_mem), (void*)&non_zero_pos_count);			
	ret |= clSetKernelArg(candidates_do_matching_kernel, 26, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);

	// calculate best work-gorup size based on number of candidates
	int best_local_work_size = cmd_queue_info.max_work_group_size;
	int best_group_size = layer_candidates_count / best_local_work_size + 1;
	int best_global_work_size = best_group_size * best_local_work_size;

	size_t global_work_load[1] = {best_global_work_size };
	size_t local_work_load[1] = {best_local_work_size};	

	// enqueue kernel for matching of parts
	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, candidates_do_matching_kernel, 1, nullptr, global_work_load , nullptr, wait_event_count, wait_evnt, evnt); 
	ocl_check_error(ret, CL_SUCCESS);

	ret = clFinish(cmd_queue_info.queue);
	ocl_check_error(ret, CL_SUCCESS);
}

void layern_creator::ocl_select_best_candidates(const best_cmd_queue_info& cmd_queue_info, cl_mem layer_candidates, cl_int layer_candidates_count,
												cl_mem kernel_err, cl_int wait_event_count, cl_event* wait_evnt, cl_event* evnt) {
	
	int ret;
	cl_kernel candidates_select_best_parts_kernel = OpenCL::context_manager->getCommonKernel("candidates_select_best_parts")[cmd_queue_info.context_number];

	// prepare next kernels while matching is being done
	ret = clSetKernelArg(candidates_select_best_parts_kernel, 0, sizeof(cl_mem), (void*)&layer_candidates);		
	ret |= clSetKernelArg(candidates_select_best_parts_kernel, 1, sizeof(cl_int), (void*)&layer_candidates_count);
	ret |= clSetKernelArg(candidates_select_best_parts_kernel, 2, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);

	
	// calculate best work-gorup size based on number of candidates
	int best_local_work_size = cmd_queue_info.max_work_group_size;
	int best_group_size = layer_candidates_count / best_local_work_size + 1;
	int best_global_work_size = best_group_size * best_local_work_size;

	size_t global_work_load[1] = {best_global_work_size };
	size_t local_work_load[1] = {best_local_work_size};	

	// and enqueue kernel to select best parts for each new location and type 
	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, candidates_select_best_parts_kernel, 1, nullptr, global_work_load , local_work_load, wait_event_count, wait_evnt, evnt); 
	ocl_check_error(ret, CL_SUCCESS);

	ret = clFinish(cmd_queue_info.queue);
	ocl_check_error(ret, CL_SUCCESS);

}

void layern_creator::ocl_make_compact_offsets(const best_cmd_queue_info& cmd_queue_info,
											  cl_mem layer_candidates, cl_int layer_candidates_count, cl_mem candidate_non_zero_positions, cl_mem candidate_non_zero_pos_count,
											  cl_mem compact_offsets, cl_mem compact_offsets_keys, cl_mem compact_offsets_count,
											  cl_mem count_new_parts_result, cl_mem count_new_edges_result, cl_mem count_valid_non_zero_pos_result,
											  cl_mem kernel_err, cl_int wait_event_count, cl_event* wait_evnt, cl_event* evnt) {
	
	int ret;
	cl_kernel make_compact_offsets_kernel = OpenCL::context_manager->getCommonKernel("make_compact_offsets")[cmd_queue_info.context_number]; 

	// remove additional candidates if g_response_threshold_percent is valid
	cl_float g_response_threshold_percent_float = g_response_threshold_percent > 0.0 && g_response_threshold_percent <= 1.0 ? g_response_threshold_percent : 0;

	// prepare next kernels while matching is being done
	ret = clSetKernelArg(make_compact_offsets_kernel, 0, sizeof(cl_mem), (void*)&layer_candidates);		
	ret |= clSetKernelArg(make_compact_offsets_kernel, 1, sizeof(cl_int), (void*)&layer_candidates_count);
	ret |= clSetKernelArg(make_compact_offsets_kernel, 2, sizeof(cl_float), (void*)&g_response_threshold_percent_float);
	ret |= clSetKernelArg(make_compact_offsets_kernel, 3, sizeof(cl_mem), (void*)&candidate_non_zero_positions);
	ret |= clSetKernelArg(make_compact_offsets_kernel, 4, sizeof(cl_mem), (void*)&candidate_non_zero_pos_count);
	ret |= clSetKernelArg(make_compact_offsets_kernel, 5, sizeof(cl_mem), (void*)&compact_offsets);
	ret |= clSetKernelArg(make_compact_offsets_kernel, 6, sizeof(cl_mem), (void*)&compact_offsets_keys);
	ret |= clSetKernelArg(make_compact_offsets_kernel, 7, sizeof(cl_mem), (void*)&compact_offsets_count);
	ret |= clSetKernelArg(make_compact_offsets_kernel, 8, sizeof(cl_mem), (void*)&count_new_parts_result);
	ret |= clSetKernelArg(make_compact_offsets_kernel, 9, sizeof(cl_mem), (void*)&count_new_edges_result);
	ret |= clSetKernelArg(make_compact_offsets_kernel, 10, sizeof(cl_mem), (void*)&count_valid_non_zero_pos_result);
	ret |= clSetKernelArg(make_compact_offsets_kernel, 11, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);

	int best_local_work_size_x = cmd_queue_info.max_work_group_size;
	int best_group_size_x = (layer_candidates_count) / best_local_work_size_x + 1;

	int best_global_work_size_x = best_group_size_x * best_local_work_size_x;

	size_t global_work_load[1] = {best_global_work_size_x};
	size_t local_work_load[1] = {best_local_work_size_x};

	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, make_compact_offsets_kernel, 1, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, evnt); 
	ocl_check_error(ret, CL_SUCCESS);
	
	ret = clFinish(cmd_queue_info.queue);
	ocl_check_error(ret, CL_SUCCESS);
}




void layern_creator::ocl_sort_compact_coords(const best_cmd_queue_info& cmd_queue_info, const cl_int coords_size_, cl_mem layer1_coords_compact, cl_mem layer1_coords_compact_sort_keys, 
												 cl_mem kernel_err, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) {
	
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

		//cout << "\t\tactual coord size in bitonic sort: " << coords_size << endl;

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

	/*
	ocl_coord_sorting_pointer* layer1_coords_compact_keys_host = new ocl_coord_sorting_pointer[coords_size];
	ret = clEnqueueReadBuffer(cmd_queue_info.queue, layer1_coords_compact_sort_keys, CL_FALSE, 0, sizeof(ocl_coord_sorting_pointer) * coords_size, layer1_coords_compact_keys_host, 1, evnt, nullptr);

	ocl_layer1_data_coordinates* layer1_coords_compact_host = new ocl_layer1_data_coordinates[coords_size];
	ret = clEnqueueReadBuffer(cmd_queue_info.queue, layer1_coords_compact, CL_FALSE, 0, sizeof(ocl_layer1_data_coordinates) * coords_size, layer1_coords_compact_host, 1, evnt, nullptr);
	ocl_check_error(ret, CL_SUCCESS);

	ret = clFinish(cmd_queue_info.queue);
	ocl_check_error(ret, CL_SUCCESS);

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

void layern_creator::ocl_update_results_offsets(const best_cmd_queue_info& cmd_queue_info, 
												const cl_int coords_size, cl_mem layer1_coords, cl_mem kernel_err, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt) {

	int ret;
	
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

}


void layern_creator::ocl_set_offsets_recursively(const best_cmd_queue_info& cmd_queue_info, cl_mem coords, int coords_size, int additional_offsets, int work_group_size_used, 
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

void layern_creator::ocl_make_results(const best_cmd_queue_info& cmd_queue_info, 
										 cl_mem lib_parts_layer_k, cl_mem lib_edges_layer_k, cl_mem lib_apps_layer_k,
										 cl_mem lib_parts_layer_k1, cl_mem img_s_nodes_layer_k1, 
										 cl_mem img_coord_layer_k1, cl_mem img_edges_layer_k1, cl_mem distr_gaussian_mask,
										 cl_mem compact_offsets, cl_mem compact_offsets_keys, cl_int compact_offsets_count,
										 cl_mem layer_candidates, cl_int layer_candidates_count, cl_mem max_schur_memory_obj, cl_uint2 img_size_layer_k1, cl_uint2 dummy_m_size, cl_int k1,
										 cl_mem new_layer_img_s_nodes, cl_mem new_layer_img_edges, cl_uint2 new_layer_img_size,
										 cl_mem kernel_err, cl_int wait_event_count, cl_event* wait_evnt, cl_event* evnt) {
	
	int ret;
	cl_kernel candidates_make_results_kernel = OpenCL::context_manager->getCommonKernel("candidates_make_results")[cmd_queue_info.context_number];

	cl_context context = OpenCL::context_manager->getContext(cmd_queue_info.context_number);

	int edges_mutex_host = 0;
	cl_mem edges_mutex = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_int), &edges_mutex_host, &ret);
	ocl_check_error(ret, CL_SUCCESS);

	cl_float c = 1.0f;	
	cl_float g_response_var_factor_float = (float)g_response_var_factor;
	cl_float convolution_threshold_float = (float)convolution_threshold;
	cl_float convolution_link_threshold_float = (float)convolution_link_threshold;

	int schur_product_use_distribution = SCHUR_PROD_USE_G_DISTR_RESPONSE;
		
	if (identity_g_response)
		schur_product_use_distribution = SCHUR_PROD_USE_IDENTITY_RESPONSE;
	else if (simple_g_response)
		schur_product_use_distribution = SCHUR_PROD_USE_R_RESPONSE;
	else if (ignore_g_distribution)
		schur_product_use_distribution = SCHUR_PROD_USE_SIMPLE_G_DISTR_RESPONSE;

	int add_reconstruction_edges_int = add_reconstruction_edges ? 1 : 0; 

	ret = clSetKernelArg(candidates_make_results_kernel, 0, sizeof(cl_mem), (void*)&lib_parts_layer_k);
	ret |= clSetKernelArg(candidates_make_results_kernel, 1, sizeof(cl_mem), (void*)&lib_edges_layer_k);
	ret |= clSetKernelArg(candidates_make_results_kernel, 2, sizeof(cl_mem), (void*)&lib_apps_layer_k);
	ret |= clSetKernelArg(candidates_make_results_kernel, 3, sizeof(cl_mem), (void*)&lib_parts_layer_k1);
	ret |= clSetKernelArg(candidates_make_results_kernel, 4, sizeof(cl_mem), (void*)&img_s_nodes_layer_k1);
	ret |= clSetKernelArg(candidates_make_results_kernel, 5, sizeof(cl_mem), (void*)&img_coord_layer_k1);
	ret |= clSetKernelArg(candidates_make_results_kernel, 6, sizeof(cl_mem), (void*)&img_edges_layer_k1);
	ret |= clSetKernelArg(candidates_make_results_kernel, 7, sizeof(cl_mem), (void*)&distr_gaussian_mask);
	ret |= clSetKernelArg(candidates_make_results_kernel, 8, sizeof(cl_mem), (void*)&compact_offsets);		
	ret |= clSetKernelArg(candidates_make_results_kernel, 9, sizeof(cl_mem), (void*)&compact_offsets_keys);		
	ret |= clSetKernelArg(candidates_make_results_kernel, 10, sizeof(cl_int), (void*)&compact_offsets_count);
	ret |= clSetKernelArg(candidates_make_results_kernel, 11, sizeof(cl_mem), (void*)&layer_candidates);		
	ret |= clSetKernelArg(candidates_make_results_kernel, 12, sizeof(cl_int), (void*)&layer_candidates_count);
	ret |= clSetKernelArg(candidates_make_results_kernel, 13, sizeof(cl_mem), (void*)&max_schur_memory_obj);
	ret |= clSetKernelArg(candidates_make_results_kernel, 14, sizeof(cl_uint2), (void*)&img_size_layer_k1);
	ret |= clSetKernelArg(candidates_make_results_kernel, 15, sizeof(cl_uint2), (void*)&dummy_m_size);
	ret |= clSetKernelArg(candidates_make_results_kernel, 16, sizeof(cl_float), (void*)&c);
	ret |= clSetKernelArg(candidates_make_results_kernel, 17, sizeof(cl_int), (void*)&k1);
	ret |= clSetKernelArg(candidates_make_results_kernel, 18, sizeof(cl_float), (void*)&convolution_threshold_float);	
	ret |= clSetKernelArg(candidates_make_results_kernel, 19, sizeof(cl_float), (void*)&convolution_link_threshold_float);
	ret |= clSetKernelArg(candidates_make_results_kernel, 20, sizeof(cl_float), (void*)&g_response_var_factor_float);
	ret |= clSetKernelArg(candidates_make_results_kernel, 21, sizeof(cl_int), (void*)&schur_product_use_distribution);
	ret |= clSetKernelArg(candidates_make_results_kernel, 22, sizeof(cl_int), (void*)&add_reconstruction_edges_int);
	ret |= clSetKernelArg(candidates_make_results_kernel, 23, sizeof(cl_mem), (void*)&new_layer_img_s_nodes);
	ret |= clSetKernelArg(candidates_make_results_kernel, 24, sizeof(cl_mem), (void*)&new_layer_img_edges);
	ret |= clSetKernelArg(candidates_make_results_kernel, 25, sizeof(cl_mem), (void*)&edges_mutex);		
	ret |= clSetKernelArg(candidates_make_results_kernel, 26, sizeof(cl_uint2), (void*)&new_layer_img_size);		
	ret |= clSetKernelArg(candidates_make_results_kernel, 27, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);
	
	int best_local_work_size_x = cmd_queue_info.max_work_group_size;
	int best_group_size_x = (compact_offsets_count) / best_local_work_size_x + 1;

	int best_global_work_size_x = best_group_size_x * best_local_work_size_x;

	size_t global_work_load[1] = {best_global_work_size_x};
	size_t local_work_load[1] = {best_local_work_size_x};

	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, candidates_make_results_kernel, 1, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, evnt); 
	ocl_check_error(ret, CL_SUCCESS);

	clReleaseMemObject(edges_mutex);

}

void layern_creator::ocl_sort_within_positions(const best_cmd_queue_info& cmd_queue_info, 
												cl_mem new_layer_img_s_nodes, cl_int compact_offsets_count,
												cl_mem new_layer_s_nodes_offsets, cl_uint2 new_layer_img_size, cl_mem merge_sort_temp_data,
												cl_mem kernel_err, cl_int wait_event_count, cl_event* wait_evnt, cl_event* evnt) {

	int ret;
	cl_kernel merge_sort_within_pos_result_kernel = OpenCL::context_manager->getCommonKernel("merge_sort_within_pos_result")[cmd_queue_info.context_number];

	cl_int2 new_layer_img_size_int;
	new_layer_img_size_int.s[0] = new_layer_img_size.s[0];
	new_layer_img_size_int.s[1] = new_layer_img_size.s[1];

	ret = clSetKernelArg(merge_sort_within_pos_result_kernel, 0, sizeof(cl_mem), (void*)&new_layer_img_s_nodes);
	ret |= clSetKernelArg(merge_sort_within_pos_result_kernel, 1, sizeof(cl_int), (void*)&compact_offsets_count);
	ret |= clSetKernelArg(merge_sort_within_pos_result_kernel, 2, sizeof(cl_mem), (void*)&new_layer_s_nodes_offsets);
	ret |= clSetKernelArg(merge_sort_within_pos_result_kernel, 3, sizeof(cl_int2), (void*)&new_layer_img_size_int);
	ret |= clSetKernelArg(merge_sort_within_pos_result_kernel, 4, sizeof(cl_mem), (void*)&merge_sort_temp_data);	
	ret |= clSetKernelArg(merge_sort_within_pos_result_kernel, 5, sizeof(cl_mem), (void*)&kernel_err);
	ocl_check_error(ret, CL_SUCCESS);
	
	int best_local_work_size_x = cmd_queue_info.max_work_group_size;
	int best_group_size_x = (compact_offsets_count) / best_local_work_size_x + 1;

	int best_global_work_size_x = best_group_size_x * best_local_work_size_x;

	size_t global_work_load[1] = {best_global_work_size_x};
	size_t local_work_load[1] = {best_local_work_size_x};

	ret = clEnqueueNDRangeKernel(cmd_queue_info.queue, merge_sort_within_pos_result_kernel, 1, nullptr, global_work_load, local_work_load, wait_event_count, wait_evnt, evnt); 
	ocl_check_error(ret, CL_SUCCESS);

}
//#undef OCL_PROFILE_HOST_QUICK 

#endif