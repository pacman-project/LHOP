
// include files
#include <utils.cl>

typedef struct {
	const __global ocl_layer1_data* layer_k1_img_parts;	
	const __global ocl_layer1_data_coordinates* layer_k1_img_parts_coord;	
	const __global ocl_edge_data_ip2* layer_k1_img_edges;
	
	float thresh; 
	float add_link_threshold;
	uint use_response_type;
	bool check_all_img_parts; // define whether to imitate as schur_product_max_all (check_all_img_parts = true ) or schur_product_max (check_all_img_parts = false)
	uint2 img_size_layer_k1;
	int layer_k1;
	float g_response_var_factor;
	float org_img_part_log_r_response;
	float org_img_part_g_response;
	int debug;	
	bool add_reconstruction_edges;
}  schur_product_in_data;

typedef struct {

	__global ocl_edge_data_ip2* result_prev_layer_data;
	__global ocl_edge_data_ip2* result_layer0_data;	
	
	int result_prev_layer_count;
	int result_layer0_count;
	
}  schur_product_out_data ;

float2 schur_product_max_all(const __private ocl_part_data_2* edge, const __constant ocl_app_data* edge_app_data,
								const uint type, const uint2 center_img_pos,
								const __private schur_product_in_data* in_data, __private schur_product_out_data* out_data);


/**
 * Calculate best R_RESPONSE and G_RESPONSE for each image position.
 * I.e. for each position find part with best R_RESPONSE and use that value 
 * multiplied with candidate_threshold_vf for this specific position only. (same for G_RESPONSE) 
 */
__kernel 
void find_best_candidate_threshold( const __global uint2* position_list, const int size, // list of offsets ands sizes that must be processed (each work-item gets one offset to process
									const __global ocl_layer1_data* layer_k1_img_parts, // pointer to all image part data for layer k - 1
									const float2 candidate_threshold_vf, // x <= R, y <= G
									const int img_width,
									__global float2* result) {

	// must calculate best R and G response for all parts in this position,
	// multiply it by candidate_vf factor and save them to list 
	
	// uint id = get_global_id(0); is faster when using directly

	// some simple error checking (errors should be checked by user app)
	if (get_work_dim() != 1)
		return; // only one dimention

	// skip ones with invalid position
	if (get_global_id(0) < size) {
		
		// get offset of part that needs to be proccessed
		uint offset = position_list[get_global_id(0)].x;
		uint part_count = position_list[get_global_id(0)].y;
		float2 max_responses = (float2)(0.0f,0.0f); // x <= R, y <= G

		uint img_x = layer_k1_img_parts[offset].x;
		uint img_y = layer_k1_img_parts[offset].y;

		// check each part in this position
		for (int i = 0; i < part_count; i++) {
			// get r and g values 
			max_responses = fmax(max_responses, (float2)(layer_k1_img_parts[offset + i].response.R, 
															layer_k1_img_parts[offset + i].response.G));
		}

		result[img_y * img_width + img_x] = max_responses * candidate_threshold_vf;
	}
}

//#define STOP_AT_X 65
//#define STOP_AT_Y 29
//#define STOP_AT_TYPE 189
/**
 * x and y coordinates must be valid for specified layer (this is not checked here, but must be done by user)
 * Use max_schur_procudt with R_RESPONSE to check if this candidate have any forbidden parts and
 * eliminate it if number of forbidden parts is above thershold.
 */
__kernel 
void candidates_eliminate_forbidden( const __global ocl_part_data* layer_k_lib_parts, // pointer to all library part data for layer k
							const __global ocl_part_data_2* layer_k_lib_edges, // pointer to all edges for layer k
							const __global ocl_part_data* layer_k1_lib_parts, // pointer to all library part data for layer k - 1
							const __global ocl_layer1_data* layer_k1_img_parts, // pointer to all image part data for layer k - 1
							const __global ocl_layer1_data_coordinates* layer_k1_img_parts_coord, // pointer to all image part coordinates for layer k - 1
							__global layer_candidate* part_candidates, const int part_candidate_size,
							const uint2 img_size_layer_k1, // size of image for layer k - 1; [width, height]
							const float c, const int layer_k1, const float convolution_thresh, 
							const float proj_max_forb_threshold, const float forb_quot_threshold,
							__global int *error) {
	//*error = 0;
	/**
	 * First just do some simple error checking.
	 * Most of error checking should be done by host.
	 */	
	// CAN ONLY EXECUTE IN 1 dimension
	if (get_work_dim() != 1) {
		*error = 1;
		return; // error
	}
		
	// get edge id number that this work-group is processing
	//uint candidate_i = get_global_id(0); is faster if using directly

	// if using more work-items then there is edges to be processed then just skip invalid ones
	if (get_global_id(0) < part_candidate_size && *error == 0) { 
		
	
		//get info for candidate
		layer_candidate candidate = part_candidates[get_global_id(0)];
	
		uint max_forbidden_count = 0;
		float max_forbidden_sum = 0;

		// get start and end offset
		uint edges_start_offset = layer_k_lib_parts[candidate.lib_part_offset].edges_loc[OCL_LYR_FORBIDDEN_EDGE].offset;
		uint edges_end_offset = edges_start_offset + layer_k_lib_parts[candidate.lib_part_offset].edges_loc[OCL_LYR_FORBIDDEN_EDGE].size;
	
		
		// prepare input data for schur_product_max_all call
		schur_product_in_data schr_in_data;
		
		schr_in_data.layer_k1_img_parts = layer_k1_img_parts;
		schr_in_data.layer_k1_img_parts_coord = layer_k1_img_parts_coord;			
		schr_in_data.layer_k1_img_edges = NULL;
		schr_in_data.thresh = convolution_thresh;
		schr_in_data.add_link_threshold = 0;
		schr_in_data.use_response_type = SCHUR_PROD_USE_VAL_RESPONSE;
		schr_in_data.check_all_img_parts = true; 
		schr_in_data.img_size_layer_k1 = img_size_layer_k1;
		schr_in_data.layer_k1 = layer_k1;
		schr_in_data.g_response_var_factor = 0;
		schr_in_data.org_img_part_log_r_response = 0;
		schr_in_data.org_img_part_g_response = 0;
		schr_in_data.add_reconstruction_edges = false;
		
		schr_in_data.layer_k1_img_parts = layer_k1_img_parts;
	
		for (uint edge_i = edges_start_offset; edge_i < edges_end_offset; edge_i++) {
			
			// each work-item processes one edge and data it is pointing to
			ocl_part_data_2 edge = layer_k_lib_edges[edge_i];
			
			// get type for this work-group
			uint type = layer_k1_lib_parts[edge.node.offset].type;

			// calculate center position where part being processed over whole kernle is
			uint2 center_img_pos = convert_uint2(convert_float2((uint2)(candidate.x, candidate.y) + (uint2)(edge.x, edge.y)) * (float2)c);
			
			#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
			if (part_candidates[get_global_id(0)].new_x ==  STOP_AT_X && part_candidates[get_global_id(0)].new_y ==  STOP_AT_Y && part_candidates[get_global_id(0)].part_type ==  STOP_AT_TYPE)
				schr_in_data.debug = 1337;
			else
				schr_in_data.debug = 0;
			#endif
				
			float2 max_schur = schur_product_max_all(&edge, NULL, type, center_img_pos, &schr_in_data, NULL);
							
			#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
			if (part_candidates[get_global_id(0)].new_x ==  STOP_AT_X && part_candidates[get_global_id(0)].new_y ==  STOP_AT_Y && part_candidates[get_global_id(0)].part_type ==  STOP_AT_TYPE)
				printf("ocl: max_schur result: %f > proj_max_forb_threshold: %f\n", max_schur.y, proj_max_forb_threshold);
			#endif
			// set forbidden if above thershold
			if (max_schur.y > proj_max_forb_threshold) {			
				max_forbidden_count += 1;
				max_forbidden_sum += max_schur.y; 
			}	
		}
						
		// if this one has forbidden parts, then set flag in candidate			
		//if (max_forbidden_count != 0 && max_forbidden_sum/max_forbidden_count >= forb_quot_threshold) {
		// lower if-clause is same as upper but without divide operation
		if (max_forbidden_sum > forb_quot_threshold * max_forbidden_count) {
			part_candidates[get_global_id(0)].skip = 1;
		}

		#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
		if (part_candidates[get_global_id(0)].new_x ==  STOP_AT_X && part_candidates[get_global_id(0)].new_y ==  STOP_AT_Y && part_candidates[get_global_id(0)].part_type ==  STOP_AT_TYPE)
			printf("ocl: skip: %d =>  max_forbidden_sum: %f > forb_quot_threshold * max_forbidden_count: %f \n",part_candidates[get_global_id(0)].skip, max_forbidden_sum, forb_quot_threshold * max_forbidden_count);
		#endif
	}
}

/**
 * x and y coordinates must be valid for specified layer (this is not checked here, but must be done by user)
 * Do matching for each candidate (one candidate - one work-item) using max_schur_product (with G_DISTR_RESPONSE)
 * and eliminate candidates that are above thershold (r_sum, g_prod, realization_ration). Also calculate how many edges 
 * we will use for this candidate. We cannot save them in this fucntions since we do not have memory for them.
 * (We need to know how many edges we will have in the first place - opencl does not support dynamic array so host must allocate memory)
 */
__kernel 
void candidates_do_matching(const __global ocl_part_data* layer_k_lib_parts, // pointer to all library part data for layer k
							const __global ocl_part_data_2* layer_k_lib_edges, // pointer to all edges for layer k
							__constant ocl_app_data* layer_k_lib_apps,
							const __global ocl_part_data* layer_k1_lib_parts, // pointer to all library part data for layer k - 1
							const __global ocl_layer1_data* layer_k1_img_parts, // pointer to all image part data for layer k - 1
							const __global ocl_layer1_data_coordinates* layer_k1_img_parts_coord, // pointer to all image part coordinates for layer k - 1
							const __global ocl_edge_data_ip2* layer_k1_img_edges, // pointer to all image edge data for layer k - 1
							__constant float* gaussian_mask,
							__global layer_candidate* part_candidates, const int part_candidate_size,
							__global float* max_schur_memory_obj,
							const uint2 img_size_layer_k1, // size of image for layer k - 1; [width, height]
							const uint2 dummy_m_size, // size of dummy distribution matrix m; [width, height] (used only in last shur_prod)
							const float c, const int layer_k1, const float convolution_threshold, const float convolution_link_threshold,						
							const float g_response_var_factor,
							const float4 default_threshold, const int manual_thresholds,
							const float r_response_pow, const float g_response_pow,
							const int schur_product_use_distribution, const int add_reconstruction_edges,
							__global ocl_layer1_data_coordinates* non_zero_positions,
							__global int *non_zero_pos_count,
							__global int *error) {
	//*error = 0;
	/**
	 * First just do some simple error checking.
	 * Most of error checking should be done by host.
	 */	
	// CAN ONLY EXECUTE IN 1 dimension
	if (get_work_dim() != 1) {
		*error = 1;
		return; // error
	}

	// get edge id number that this work-group is processing
	//uint candidate_i = get_global_id(0); is faster if using directly

	// if using more work-items then there is edges to be processed then just skip invalid ones
	if (get_global_id(0) < part_candidate_size && *error == 0) {  
	
		//get info for candidate		
		layer_candidate candidate = part_candidates[get_global_id(0)];
		
		float r_sum = 0;
		float g_prod = 1;		
		float realization_ratio = 0;
		float bpn = layer_k_lib_parts[candidate.lib_part_offset].bpcount; // basic part number
		float bpr = 0; // count how many basic parts are realized
		
		float org_img_part_g_response = layer_k1_img_parts[candidate.img_part_offset].response.G; // nd->r(G_RESPONSE);
		float org_img_part_rr_response = layer_k1_img_parts[candidate.img_part_offset].response.RR; // nd->r(RR_RESPONSE);
		float org_img_part_r_response = layer_k1_img_parts[candidate.img_part_offset].response.R; // nd->r(R_RESPONSE);
		
		// get start and end offset
		uint edges_start_offset = layer_k_lib_parts[candidate.lib_part_offset].edges_loc[OCL_LYR_SRC_EDGE].offset;
		uint part_count = layer_k_lib_parts[candidate.lib_part_offset].edges_loc[OCL_LYR_SRC_EDGE].size;
		uint edges_end_offset = edges_start_offset + part_count;
		
		// prepare input data for schur_product_max_all call
		schur_product_in_data schr_in_data;		
		schr_in_data.layer_k1_img_parts = layer_k1_img_parts;
		schr_in_data.layer_k1_img_parts_coord = layer_k1_img_parts_coord;	
		schr_in_data.layer_k1_img_edges = NULL;		
		schr_in_data.thresh = convolution_threshold;		
		schr_in_data.use_response_type = schur_product_use_distribution;
		schr_in_data.check_all_img_parts = true;
		schr_in_data.img_size_layer_k1 = img_size_layer_k1;
		schr_in_data.layer_k1 = layer_k1;
		schr_in_data.g_response_var_factor = g_response_var_factor;
		schr_in_data.org_img_part_log_r_response = log(org_img_part_r_response);
		schr_in_data.org_img_part_g_response = org_img_part_g_response;		
		schr_in_data.add_reconstruction_edges = add_reconstruction_edges;
		
		// ignore convolution_link_thershold in here since we only need number of elements used
		schr_in_data.add_link_threshold = 0; 

		if (add_reconstruction_edges)
			schr_in_data.layer_k1_img_edges = layer_k1_img_edges;
		
		// prepare output data for schur_product_max_all call
		schur_product_out_data schr_out_data;		
		schr_out_data.result_prev_layer_data = NULL;
		schr_out_data.result_prev_layer_count = 0;
		schr_out_data.result_layer0_data = NULL;
		schr_out_data.result_layer0_count = 0;
		
		uint best_max_schur_offset = candidate.max_schur_mem_offset;
		for (uint edge_i = edges_start_offset; edge_i < edges_end_offset; edge_i++) {
			
			// each work-item processes one edge and data it is pointing to
			ocl_part_data_2 edge = layer_k_lib_edges[edge_i];
			
			// get type for this work-group
			uint type = layer_k1_lib_parts[edge.node.offset].type;
			float basic_part_number = layer_k1_lib_parts[edge.node.offset].bpcount;
			
			// calculate center position where part being processed over whole kernle is
			uint2 center_img_pos = convert_uint2(convert_float2((uint2)(candidate.x, candidate.y) + (uint2)(edge.x, edge.y)) * (float2)c);

			#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
			if (candidate.new_x ==  STOP_AT_X && candidate.new_y ==  STOP_AT_Y && candidate.part_type ==  STOP_AT_TYPE)
				schr_in_data.debug = 1337;
			else
				schr_in_data.debug = 0;
			#endif

			float2 max_schur = schur_product_max_all(&edge, &layer_k_lib_apps[edge.app.offset], type, center_img_pos, &schr_in_data, &schr_out_data);
			
			#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
			if (candidate.new_x ==  STOP_AT_X && candidate.new_y ==  STOP_AT_Y && candidate.part_type ==  STOP_AT_TYPE)
				printf("ocl: max_schur result: %f > convolution_threshold: %f\n", max_schur.y, convolution_threshold);
			#endif
			
			if (max_schur.y >= convolution_threshold) {	

				#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
				if (candidate.new_x ==  STOP_AT_X && candidate.new_y ==  STOP_AT_Y && candidate.part_type ==  STOP_AT_TYPE)
					printf("ocl: selected: part loc: %d; x: %d; y: %d; m: %d; r_resp: %f; current r_sum: %f\n",
						convert_int(max_schur.x), layer_k1_img_parts[convert_int(max_schur.x)].x, layer_k1_img_parts[convert_int(max_schur.x)].y, layer_k1_img_parts[convert_int(max_schur.x)].m, layer_k1_img_parts[convert_int(max_schur.x)].response.R, r_sum);
				#endif
								
				// Update r_sum and g_sum (used in the calculation of R_RESPONSE and G_RESPONSE)
                // and bpr (used in the calculation of RR_RESPONSE)
                bpr += basic_part_number * layer_k1_img_parts[convert_int(max_schur.x)].response.RR; // pd->get_basic_part_number(p)*nnd->r(RR_RESPONSE);
				r_sum += layer_k1_img_parts[convert_int(max_schur.x)].response.R; // pd->get_basic_part_number(p)*nnd->r(R_RESPONSE);
				g_prod *= max_schur.y;				
				
				//realization_ratio += 1; (depricated with introduction of bpcount)
				
			}

			if (convolution_link_threshold > 0)
				max_schur_memory_obj[best_max_schur_offset++] = max_schur.y;
			
		}	
		
		#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
		if (candidate.new_x ==  STOP_AT_X && candidate.new_y ==  STOP_AT_Y && candidate.part_type ==  STOP_AT_TYPE)
			printf("ocl: first  bpr: %f; bpn: %f; r_sum: %f; g_prod: %f\n", bpr, bpn, r_sum, g_prod);
		#endif
		
		// Calculate RR_RESPONSE
        bpr += layer_k1_lib_parts[candidate.img_part_type].bpcount * org_img_part_rr_response; // cpartd->get_basic_part_number(cpart) * nd->r(RR_RESPONSE);
		realization_ratio = bpr/bpn;
		
		// Calculate R_RESPONSE (from r_sum) and G_RESPONSE (from g_prod)				
		r_sum += org_img_part_r_response;
		r_sum /= (part_count + 1);
		r_sum = pow(r_sum, r_response_pow);		
		g_prod = pow(g_prod * realization_ratio, g_response_pow);

		#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
		if (candidate.new_x ==  STOP_AT_X && candidate.new_y ==  STOP_AT_Y && candidate.part_type ==  STOP_AT_TYPE)
			printf("ocl:  bpr: %f; bpn: %f; r_sum: %f; g_prod: %f\n", bpr, bpn, r_sum, g_prod);
		#endif
		
		const __global ocl_part_data* part = layer_k_lib_parts + candidate.lib_part_offset;
		
		// use thershold from lib part or default value when using manual thershold or thershold in lib part does not exists (thershold does not exists if is negative) 
		float4 thershold = manual_thresholds != 0 ? default_threshold : (float4)((part->thresh.R >=0 ? part->thresh.R : default_threshold.x),
																					(part->thresh.RR >=0 ? part->thresh.RR : default_threshold.y),
																					(part->thresh.G >=0 ? part->thresh.G : default_threshold.z), 0);
		
		bool this_candidate_is_valid = r_sum >= thershold.x && realization_ratio >= thershold.y && g_prod >= thershold.z ? true : false;
		
		if (this_candidate_is_valid == true) {
			//part_candidates[get_global_id(0)].skip = 0;
			part_candidates[get_global_id(0)].r_sum = r_sum;
			part_candidates[get_global_id(0)].realization_ratio = realization_ratio;
			part_candidates[get_global_id(0)].g_prod = g_prod;
		} else {
			part_candidates[get_global_id(0)].skip = 1;
			part_candidates[get_global_id(0)].r_sum = 0;
			part_candidates[get_global_id(0)].realization_ratio = 0;
			part_candidates[get_global_id(0)].g_prod = 0;
		}
		
		// in case this candidate might be valid, we still need to count subparts from center since they are not included in OCL_LYR_SRC_EDGE 
		if (this_candidate_is_valid) {
			ocl_part_data_2 edge;
			edge.x = 0; edge.y = 0;			

			// copy gaussion mask from global to private			
			edge.distr_size = dummy_m_size.y * dummy_m_size.x;
			for (int i = 0; i < edge.distr_size; i++) {
				edge.distr[i] = gaussian_mask[i];
			}
			
			uint type = layer_k1_img_parts[candidate.img_part_offset].m;

			// only change type of shcur product to use; other values should stay the same (as when they are defined)
			schr_in_data.use_response_type = SCHUR_PROD_USE_R_RESPONSE;

			// we only need this to calculate how many edges we will use
			float2 max_schur = schur_product_max_all(&edge, NULL, type, (uint2)(candidate.x, candidate.y), &schr_in_data, &schr_out_data);

			if (convolution_link_threshold > 0)
				max_schur_memory_obj[best_max_schur_offset++] = max_schur.y;
		}

		part_candidates[get_global_id(0)].result_edges_prev_lay_size = schr_out_data.result_prev_layer_count;
		part_candidates[get_global_id(0)].result_edges_layer0_size = schr_out_data.result_layer0_count;		

		// in case this part candidate is first in its position then we need to set non_zero_positions to this offset
		if (*non_zero_pos_count >= 0) {

			if (get_global_id(0) == 0 || 
					candidate.new_x != part_candidates[get_global_id(0) - 1].new_x || 
					candidate.new_y != part_candidates[get_global_id(0) - 1].new_y) {
				int pos = atomic_inc(non_zero_pos_count);
				non_zero_positions[pos].offset = get_global_id(0);
				non_zero_positions[pos].size = 1;
			}
		}
	}
}

/**
 * Among all candidates in same img position we need to use only one part for each type.
 * For each position and type we select only part with best val().
 * This is done using simple loop where we check all neighbors for their val()
 * and eliminate candidates with lower val().
 */
__kernel
void candidates_select_best_parts(__global layer_candidate* part_candidates, const int part_candidate_size,
									__global int *error) {

	// CAN ONLY EXECUTE IN 1 dimension
	if (get_work_dim() != 1) {
		*error = 1;
		return; // error
	}

	// if using more work-items then there is edges to be processed then just skip invalid ones
	if (get_global_id(0) < part_candidate_size && *error == 0) {
	
		// get info for candidate
		int new_x = part_candidates[get_global_id(0)].new_x;
		int new_y = part_candidates[get_global_id(0)].new_y;
		int part_type = part_candidates[get_global_id(0)].part_type;
		float r_sum = part_candidates[get_global_id(0)].r_sum;
		float val = ocl_layer_candidate_get_val(part_candidates[get_global_id(0)]);
		
		bool stop = part_candidates[get_global_id(0)].skip == 0 ? false : true;
		uint c_offset = get_global_id(0) + 1;
		// this loop is executed for each candidate as its own work-item and does the following:
			// Check for all neighbor parts (that belong to same new position - all with same position are in continuous memory)
			// if they are of same type and if they are better (i.e. val() is higher).
			// Always mark with "skip" flag parts with lower val().
			// Each candidate will check all other candidtes in same position 			
			// therefore there is some owerhead, which can be eliminated by
			// checking only parts saved after this one.			
		while (!stop) {
			// watch for out-of-bounds and candidates on new position
			if (c_offset > part_candidate_size ||
				part_candidates[c_offset].new_x != new_x || part_candidates[c_offset].new_y != new_y ) {
				stop = true;
				continue;
			}			
			// check if is valid and of same type
			if (part_candidates[c_offset].skip == 0 && 
				part_candidates[c_offset].part_type == part_type) {
				
				float current_part_check_value = ocl_layer_candidate_get_val(part_candidates[c_offset]);
				
				if (fabs(current_part_check_value - val) < 1e-6) {
					// in case both val() are the same then use one with greater r_sum
					
					current_part_check_value = part_candidates[c_offset].r_sum;
					
					if (current_part_check_value >= r_sum) {
					
						#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
						if (new_x ==  STOP_AT_X && new_y ==  STOP_AT_Y && part_type ==  STOP_AT_TYPE) 
							printf("org: skipping from r_sum a: r_sum %f; g_prod: %f (with respect to b: r_sum %f; g_prod: %f ) \n", part_candidates[get_global_id(0)].r_sum, part_candidates[get_global_id(0)].g_prod, part_candidates[c_offset].r_sum, part_candidates[c_offset].g_prod);
						#endif
					
						// part that this work-item is assigned to is not good
						part_candidates[get_global_id(0)].skip = 1;
						stop = true;
					} else if (current_part_check_value < r_sum) {
					
						#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
						if (new_x ==  STOP_AT_X && new_y ==  STOP_AT_Y && part_type ==  STOP_AT_TYPE) 
							printf("org: skipping from r_sum b: r_sum %f; g_prod: %f (with respect to a: r_sum %f; g_prod: %f )\n", part_candidates[c_offset].r_sum, part_candidates[c_offset].g_prod,  part_candidates[get_global_id(0)].r_sum, part_candidates[get_global_id(0)].g_prod);
						#endif
						
						// part that we are now checking is not good
						part_candidates[c_offset].skip = 1;
					}				
				} else if (current_part_check_value > val) {
					
					#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
					if (new_x ==  STOP_AT_X && new_y ==  STOP_AT_Y && part_type ==  STOP_AT_TYPE) 
						printf("org: skipping from g_prod a: r_sum %f; g_prod: %f (with respect to b: r_sum %f; g_prod: %f )\n", part_candidates[get_global_id(0)].r_sum, part_candidates[get_global_id(0)].g_prod, part_candidates[c_offset].r_sum, part_candidates[c_offset].g_prod);
					#endif
					
					// part that this work-item is assigned to is not good
					part_candidates[get_global_id(0)].skip = 1;
					stop = true;
				} else if (current_part_check_value < val) { // must be 'if' here in case when both values would be the same

					#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
					if (new_x ==  STOP_AT_X && new_y ==  STOP_AT_Y && part_type ==  STOP_AT_TYPE) 
						printf("org: skipping from g_prod b: r_sum %f; g_prod: %f (with respect to a: r_sum %f; g_prod: %f )\n", part_candidates[c_offset].r_sum, part_candidates[c_offset].g_prod, part_candidates[get_global_id(0)].r_sum, part_candidates[get_global_id(0)].g_prod);				
					#endif
						
					// part that we are now checking is not good
					part_candidates[c_offset].skip = 1;
				} 
			}
			c_offset++;
		}	
	}
}


__kernel
void make_compact_offsets(__global layer_candidate* part_candidates, const int part_candidate_size, const float g_response_threshold_percent,
							const __global ocl_layer1_data_coordinates* non_zero_positions,	const __global int *non_zero_pos_count,
							__global ocl_layer1_data_coordinates* compact_offsets, __global ocl_coord_sorting_pointer* compact_offsets_keys, __global int *compact_offsets_count,
							__global int *count_new_parts, __global int *count_new_edges, __global int *count_valid_non_zero_pos,
							__global int *error) {
	// CAN ONLY EXECUTE IN 1 dimension
	if (get_work_dim() != 1) {
		*error = 1;
		return; // error
	}

	if (get_global_id(0) < *non_zero_pos_count && non_zero_positions[get_global_id(0)].size > 0) {
		
		int offset = non_zero_positions[get_global_id(0)].offset;

		float max_value = 0;

		int current_count_new_parts = 0;
		int current_count_new_edges = 0;		

		int current_pos_x = part_candidates[offset].new_x;
		int current_pos_y = part_candidates[offset].new_y;
		
		int size = -1;
		while (offset + (++size) < part_candidate_size && 
			current_pos_x == part_candidates[offset + size].new_x && 
			current_pos_y == part_candidates[offset + size].new_y) {
			
			if (part_candidates[offset + size].skip == 0) {
				// get best value in this position
				max_value = max(ocl_layer_candidate_get_val(part_candidates[offset + size]), max_value);

				// calculate number of parts and edges
				current_count_new_parts++;
				current_count_new_edges += part_candidates[offset + size].result_edges_prev_lay_size + part_candidates[offset + size].result_edges_layer0_size;

			}
		}

		// in case we need to use g_response_threshold remove some invlaid parts
		if (g_response_threshold_percent > 0) {
			float thresh = max_value * g_response_threshold_percent;
		
			while (--size >= 0) { // need prefix decrese since we have also incresed size in last while loop
				// skip candidate if nd->val() is less then thershold value
				if (part_candidates[offset + size].skip == 0 && ocl_layer_candidate_get_val(part_candidates[offset + size]) < thresh) {
					// set to skip 
					part_candidates[offset + size].skip = 1;

					// and remove its parts from count of count_new_parts and count_new_edges
					current_count_new_parts--;
					current_count_new_edges -= part_candidates[offset + size].result_edges_prev_lay_size + part_candidates[offset + size].result_edges_layer0_size;

				}
			}
		}

		if (current_count_new_parts > 0) {			
			// update count_new_parts, count_new_edges and count_valid_non_zero_pos
			atomic_add(count_new_parts, current_count_new_parts);
			atomic_add(count_new_edges, current_count_new_edges);
			int mem_pos = atomic_inc(count_valid_non_zero_pos);

			// set compact coordinates and its keys

			// set number of parts in this location
			compact_offsets[mem_pos].size = current_count_new_parts;
			compact_offsets[mem_pos].offset = 0;
			// set max value for sorting and offset that indicates where in part_candidates memory does first part of this position start
			compact_offsets_keys[mem_pos].sort_value = max_value;
			compact_offsets_keys[mem_pos].coord_ind = offset;
		}
	}
}


__kernel
void candidates_make_results(const __global ocl_part_data* layer_k_lib_parts,
								const __global ocl_part_data_2* layer_k_lib_edges,
								__constant ocl_app_data* layer_k_lib_apps,
								const __global ocl_part_data* layer_k1_lib_parts, 
								const __global ocl_layer1_data* layer_k1_img_parts,
								const __global ocl_layer1_data_coordinates* layer_k1_img_parts_coord, 
								const __global ocl_edge_data_ip2* layer_k1_img_edges, 
								__constant float* gaussian_mask,
								const __global ocl_layer1_data_coordinates* compact_offsets, 
								const __global ocl_coord_sorting_pointer* compact_offsets_keys,
								const int compact_offsets_count,
								const __global layer_candidate* part_candidates, const int part_candidate_size,			
								const __global float* max_schur_memory_obj,
								const uint2 img_size_layer_k1,
								const uint2 dummy_m_size, // size of dummy distribution matrix m; [width, height] (used only in last shur_prod)
								const float c, const int layer_k1, const float convolution_threshold, const float convolution_link_threshold, 
								const float g_response_var_factor, const int schur_product_use_distribution,
								const int add_reconstruction_edges,
								__global ocl_layer1_data* new_layer_img_s_nodes,
								__global ocl_edge_data_ip2* new_layer_img_edges, __global int* edges_mutex,
								const uint2 new_layer_img_size,
								__global int *error) {
	// CAN ONLY EXECUTE IN 1 dimension
	if (get_work_dim() != 1) {
		*error = 1;
		return; // error
	}
	
	// if using more work-items then there is edges to be processed then just skip invalid ones
	if (get_global_id(0) < compact_offsets_count && *error == 0) {	

		// we saved original position (where saved in part_candidates) to candidtes for this position
		int candidates_position_offset = compact_offsets_keys[get_global_id(0)].coord_ind;
		
		int save_result_offset = compact_offsets[get_global_id(0)].offset;

		// get first candidate in this position (first candidate should always be valid one)
		layer_candidate final_candidate = part_candidates[candidates_position_offset];

		int current_pos_x = final_candidate.new_x;
		int current_pos_y = final_candidate.new_y;

		int offset = 0;
		int save_local_offset = 0;

		// process for each candidate in this position
		bool process_next_candidate = true;
		while (process_next_candidate == true) {
			
			// make parts only from valid candidates
			if (final_candidate.skip == 0) {		
				// for every valid candidate reserve place in new_layer_img_parts and for its edges in new_layer_img_edges
				
				// get pointer to this reserved location
				__global ocl_layer1_data* new_layer1_data = &new_layer_img_s_nodes[save_result_offset + save_local_offset];

				// and set all values exept for edges
				new_layer1_data->x = current_pos_x;
				new_layer1_data->y = current_pos_y;
				new_layer1_data->z = layer_k1 + 1;
				new_layer1_data->m = final_candidate.part_type;
				new_layer1_data->attr = 0;
				
				
				ocl_set_attr(&new_layer1_data->attr, GRID_NODE_ATTR);
				
				new_layer1_data->response.R = final_candidate.r_sum;
				new_layer1_data->response.RR = final_candidate.realization_ratio;
				new_layer1_data->response.G = final_candidate.g_prod;
				
				// set all edges to 0
				for (int j = 0; j < OCL_EDGE_COUNT; j++) {
					new_layer1_data->edges_loc[j].offset = 0;
					new_layer1_data->edges_loc[j].size = 0;
				}		
		
				#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
				if (new_layer1_data->x ==  STOP_AT_X && new_layer1_data->y ==  STOP_AT_Y)
					printf("ocl: final  r_sum: %f; g_prod: %f; type: %d\n", new_layer1_data->response.R, new_layer1_data->response.G, new_layer1_data->m);
				#endif
			
				
				// set attribute for img part from which this new part is created
	//			ocl_set_attr(&layer_k1_img_parts[final_candidate->img_part_offset].attr, HAS_NEXT_LAYER);
				
				// finally we also need to find all connections (i.e. edges back to layer k - 1)
				// but first we need to reserve enough space for all subparts in this part
				int reserved_memory_size = final_candidate.result_edges_prev_lay_size + final_candidate.result_edges_layer0_size;
				int subparts_offset = atomic_add(edges_mutex, reserved_memory_size);

				// we need to save that offset
				new_layer1_data->edges_loc[OCL_TO_PREV_LAYER_EDGE].offset = subparts_offset;
				new_layer1_data->edges_loc[OCL_TO_LAYER_0].offset = subparts_offset + final_candidate.result_edges_prev_lay_size; // layer0 should be seperated from prev_layer edges
				
				__global ocl_edge_data_ip2* subparts_reserved_space = &new_layer_img_edges[subparts_offset];
				
				// get start and end offset
				uint edges_start_offset = layer_k_lib_parts[final_candidate.lib_part_offset].edges_loc[OCL_LYR_SRC_EDGE].offset;
				uint part_count = layer_k_lib_parts[final_candidate.lib_part_offset].edges_loc[OCL_LYR_SRC_EDGE].size;
				uint edges_end_offset = edges_start_offset + part_count;
				
				// prepare input data for schur_product_max_all call
				schur_product_in_data schr_in_data;
				schr_in_data.layer_k1_img_parts = layer_k1_img_parts;			
				schr_in_data.layer_k1_img_parts_coord = layer_k1_img_parts_coord;
				schr_in_data.layer_k1_img_edges = NULL;
				schr_in_data.thresh = convolution_threshold;
				schr_in_data.use_response_type = schur_product_use_distribution;
				schr_in_data.check_all_img_parts = true;
				schr_in_data.img_size_layer_k1 = img_size_layer_k1;
				schr_in_data.layer_k1 = layer_k1;
				schr_in_data.g_response_var_factor = g_response_var_factor;
				schr_in_data.org_img_part_log_r_response = log(layer_k1_img_parts[final_candidate.img_part_offset].response.R);
				schr_in_data.org_img_part_g_response = layer_k1_img_parts[final_candidate.img_part_offset].response.G;
				schr_in_data.add_reconstruction_edges = add_reconstruction_edges;
				
				// use 0 as default actual convolution_link_thershold
				schr_in_data.add_link_threshold = 0;

				if (add_reconstruction_edges) 
					schr_in_data.layer_k1_img_edges = layer_k1_img_edges;
					
				// prepare output data for schur_product_max_all call
				schur_product_out_data schr_out_data;			
				schr_out_data.result_prev_layer_data = NULL; // THIS MUST be set during each loop
				schr_out_data.result_prev_layer_count = 0;
				
				schr_out_data.result_layer0_data = NULL; // THIS MUST be set during each loop
				// space for edges to layer0 must be in continious memory so set offset to number of other subparts saved
				schr_out_data.result_layer0_count = 0;
				
				schr_out_data.result_prev_layer_data = subparts_reserved_space ;
				schr_out_data.result_layer0_data = subparts_reserved_space + final_candidate.result_edges_prev_lay_size;
				
				uint best_max_schur_offset = final_candidate.max_schur_mem_offset;
				for (uint edge_i = edges_start_offset; edge_i < edges_end_offset; edge_i++) {				
					
					// each work-item processes one edge and data it is pointing to
					ocl_part_data_2 edge = layer_k_lib_edges[edge_i];
					
					// get type for this edge
					uint type = layer_k1_lib_parts[edge.node.offset].type;

					// calculate center position where part being processed over whole kernle is
					uint2 center_img_pos = convert_uint2(convert_float2((uint2)(final_candidate.x, final_candidate.y) + (uint2)(edge.x, edge.y)) * (float2)c);

					#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
					if (new_layer1_data->x ==  STOP_AT_X && new_layer1_data->y ==  STOP_AT_Y)
						schr_in_data.debug = 1337;
					else
						schr_in_data.debug = 0;
					#endif

					// set actual convolution_link_threshold (i.e. add_link_threshold) for this edge
					if (convolution_link_threshold > 0)
						schr_in_data.add_link_threshold = max_schur_memory_obj[best_max_schur_offset++] * convolution_link_threshold;
					
					// use schur_product_max_all with same settings as in candidates_do_matching to save all subparts to global memory
					schur_product_max_all(&edge, &layer_k_lib_apps[edge.app.offset], type, center_img_pos, &schr_in_data, &schr_out_data);									
				}	
				
				#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
				schr_in_data.debug = 0;
				#endif

				// now we still need to copy all subparts to center
				ocl_part_data_2 edge;
				edge.x = 0; edge.y = 0;
				
				// copy gaussion mask from global to private
				edge.distr_size = dummy_m_size.y * dummy_m_size.x;
				for (int i = 0; i < edge.distr_size; i++) {
					edge.distr[i] = gaussian_mask[i];
				}			
				
				uint type = layer_k1_img_parts[final_candidate.img_part_offset].m;
				
				// prepare input data for final schur_product_max_all call
				schr_in_data.use_response_type = SCHUR_PROD_USE_R_RESPONSE;

				// set actual convolution_link_threshold (i.e. add_link_threshold) for this edge
				if (convolution_link_threshold > 0)
					schr_in_data.add_link_threshold = max_schur_memory_obj[best_max_schur_offset++] * convolution_link_threshold;

				schur_product_max_all(&edge, NULL, type, (uint2)(final_candidate.x, final_candidate.y), &schr_in_data, &schr_out_data);
				
				// finaly we need to save how many subparts have we created
				new_layer1_data->edges_loc[OCL_TO_PREV_LAYER_EDGE].size = schr_out_data.result_prev_layer_count;
				new_layer1_data->edges_loc[OCL_TO_LAYER_0].size = schr_out_data.result_layer0_count;
			
				// update local offset
				save_local_offset++;
			}
			// update offset
			offset++;

			// and update to next candidate if it is valid one
			if (candidates_position_offset + offset < part_candidate_size) {
				final_candidate = part_candidates[candidates_position_offset + offset];
					
				// stop when next candidate is not on same position as this one
				if (current_pos_x != final_candidate.new_x || current_pos_y != final_candidate.new_y) {
					process_next_candidate = false;
				} 
			} else {
				process_next_candidate = false;
			}
		}			
	}
}

float2 schur_product_max_all(const __private ocl_part_data_2* edge, const __constant ocl_app_data* edge_app_data,
								const uint type, const uint2 center_img_pos,
								const __private schur_product_in_data* in_data, __private schur_product_out_data* out_data){


	// prepare values for loop
	float2 best_part = (float2)(-1, 0); 

	// variance can be changed based on parameter g_response_var_factor
	float variance = (in_data->g_response_var_factor > 0 ? in_data->g_response_var_factor : 1) * edge->gdistr.variance;

	float mean = edge->gdistr.mean;

	float gdistr_n_factor = variance <= 0 ? 0 : M_SQRT1_2_E * 1 / native_sqrt(M_PI_E * variance);
	uint result_best_id_offset = out_data != NULL ? out_data->result_prev_layer_count : 0 ;
	uint result_layer0_offset = out_data != NULL ? out_data->result_layer0_count : 0;	
	
	uint2 m_size = convert_uint2(sqrt((float2)edge->distr_size));
	uint2 m_size_half = m_size/(uint2)2;

	for (uint i = 0; i < m_size.y; i++) {
		for (uint j = 0; j < m_size.x; j++) {
			
			// calculate global (x,y) position for image parts for this work-item
			uint2 global_img_pos = center_img_pos + ( (uint2)(j, i) - m_size_half) ;
			
			uint pos = global_img_pos.y * in_data->img_size_layer_k1.x + global_img_pos.x;		

			// coordinates for all image parts for this position (x,y) start at this offset and size
			uint parts_offset = in_data->layer_k1_img_parts_coord[pos].offset; 
			uint parts_count = in_data->layer_k1_img_parts_coord[pos].size;

			uint n = i * m_size.x + j;
			
			// get matrix distribution value (factor) for this work-item within work-group size
			float m_factor = edge->distr[n];
			
			// in case we only check first part
			if (in_data->check_all_img_parts == false && parts_count > 0) {
				parts_count = min(parts_count, (uint)1);
			}			
				
			// now find each part at this specific location
			for (uint k = 0; k < parts_count; k++) {
				ocl_layer1_data part = in_data->layer_k1_img_parts[parts_offset + k]; // equivalent of nd in layer1_result::schur_product_max_all
					
				// verify type of part by matching with edge->app
				bool is_valid_type = false;
				
				float app_factor = 1;
				
				#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
				//if (in_data->debug == 1337)
					//printf("ocl schur prod: x: %d; y: %d; part.m: %d vs type: %d\n", part.x, part.y, part.m ,type);
				#endif
				
				// decide if is valid part type and calculate app factor for part 
				if (edge->app.size > 0 && edge_app_data != NULL) {					
					for (int l = 0; l < edge->app.size; l++) {
						if (part.m == edge_app_data[l].type) {
							is_valid_type = true;
							app_factor = edge_app_data[l].value;
							break;
						}
					}
				} else {
					is_valid_type = part.m == type ? true : false;
				}
				
				#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
				//if (in_data->debug == 1337)
					//printf("ocl schur prod: is_valid_type: %d\n", is_valid_type);
				#endif
				
				// if is valid type then compute responses based on selected method
				if (is_valid_type) {
					float response_times_m_factor = 0;
					// for part.response we can use three different values:
						// part.response.R
						// static value 1 
						// dist.pdf_val1(::log(part.response.R) - log(org_img_part.response.R) )
							// where dist is normal_distribution with mean and variance from edge

					switch (in_data->use_response_type) {
						// cases used with candidates_eliminate_forbidden
						case SCHUR_PROD_USE_VAL_RESPONSE: {	 // when using in candidates_eliminate_forbidden (equivalent to &layer1_result::vspf)
							
							response_times_m_factor = ocl_data_get_val(part);
							break;
						}					
						
						// cases used with candidates_do_matching
						case SCHUR_PROD_USE_IDENTITY_RESPONSE: { // when identity_g_response == true (equivalent to &layer1_result::idspf)
							response_times_m_factor = 1;
						}							
						case SCHUR_PROD_USE_R_RESPONSE: {	// when simple_g_response == true (equivalent to &layer1_result::rspf)
							response_times_m_factor = part.response.R;
							break;
						}
						case SCHUR_PROD_USE_SIMPLE_G_DISTR_RESPONSE: { // when ignore_g_distribution == true (equivalent to &layer1_result::sgspf)
							response_times_m_factor = part.response.G / part.response.RR;
							break;
						}
						case SCHUR_PROD_USE_G_DISTR_RESPONSE: { // when non of above is true (equivalent to g_response_spf gspf(ndval_log, ed->gdistr; ))
						// i.e. when simple_g_response != true && ignore_g_distribution != true && identity_g_response != true
						
							// implement normal distribution
							if (variance <= 0) {
								response_times_m_factor = 0;
							} else {
								float a = (log(part.response.R) - in_data->org_img_part_log_r_response) - mean;
								response_times_m_factor = exp(-a*a/2/variance);								
							}
							
							response_times_m_factor = response_times_m_factor * part.response.G / part.response.RR;
							break;
						}					
					}
					
					#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
					if (in_data->debug == 1337) {
						printf("ocl schur prod: x: %d; y: %d; part.m: %d vs type: %d\n", part.x, part.y, part.m ,type);
						printf("ocl schur prod: response_times_m_factor: %f; m_factor: %f; app_factor: %f\n", response_times_m_factor, m_factor, app_factor);
						printf("ocl schur prod: r: %f > thersh: %f == > %d\n", response_times_m_factor  * m_factor * app_factor, in_data->thresh,response_times_m_factor > in_data->thresh);
					}
					#endif
					
					response_times_m_factor = response_times_m_factor  * m_factor * app_factor;
					
					if (response_times_m_factor > in_data->thresh) {
						
						
						// when found correct one then save its positin and value if best one 
						if (response_times_m_factor > best_part.y) {
							best_part.x = parts_offset + k;
							best_part.y = response_times_m_factor;
						}
					
						// save index in case we also need list of all parts above some thresh
						// save only if value is above actual convolution_link_thershold (i.e. add_link_threshold = convolution_link_thershold * best schur_prod)
						#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
						if (in_data->debug == 1337) {
							printf("ocl schur prod: reesponse_times_m_factor > convolution_link_thershold : %f > %f \n", response_times_m_factor, in_data->add_link_threshold);
						}
						#endif

						if (out_data != NULL && response_times_m_factor > in_data->add_link_threshold) {
							if (out_data->result_prev_layer_data > 0) {
								out_data->result_prev_layer_data[result_best_id_offset].x = edge->x;
								out_data->result_prev_layer_data[result_best_id_offset].y = edge->y;
								out_data->result_prev_layer_data[result_best_id_offset].r = response_times_m_factor;
								out_data->result_prev_layer_data[result_best_id_offset].node.layer = in_data->layer_k1;
								out_data->result_prev_layer_data[result_best_id_offset].node.offset = parts_offset + k;
								result_best_id_offset++;
							}					
							out_data->result_prev_layer_count += 1;
								 
								 
							if (in_data->add_reconstruction_edges != 0) {
								// for part we are currently checking get each edge to layer0 and copy it for this new part
								
								// use seperate code for first layer 
								if (in_data->layer_k1 == 0) {				
									// in this case we only need to copy same edge as in result_prev_layer_data data
									// since layer k - 1 will NOT have any edges and edges in result_prev_layer_data actualy point
									// to layer0
									if (out_data->result_layer0_data > 0) {
										
										out_data->result_layer0_data[result_layer0_offset].x = 0;
										out_data->result_layer0_data[result_layer0_offset].y = 0;
										out_data->result_layer0_data[result_layer0_offset].r = 1;
										out_data->result_layer0_data[result_layer0_offset].type = 0;
										out_data->result_layer0_data[result_layer0_offset].node.layer = in_data->layer_k1;
										out_data->result_layer0_data[result_layer0_offset].node.offset = out_data->result_prev_layer_data[result_best_id_offset - 1].node.offset;
										
										result_layer0_offset++;
									}
									
									// output how many layer0 edges (should be/has been) written
									out_data->result_layer0_count += 1;
								} else {
									uint to_layer0_edges_offset = part.edges_loc[OCL_TO_LAYER_0].offset; 
									uint to_layer0_edges_count = part.edges_loc[OCL_TO_LAYER_0].size;
									
									uint edges_added_count = 0;
									
									if (out_data->result_layer0_data > 0) {
										for (uint l = 0; l < to_layer0_edges_count; l++) {
											ocl_edge_data_ip2 e = in_data->layer_k1_img_edges[to_layer0_edges_offset + l];
											
											bool add_lay0_edge = true;
											// verify that same part does not yet exists
											for (uint o = 0; o < result_layer0_offset; o++) {
												if (out_data->result_layer0_data[o].node.offset == e.node.offset) {
													add_lay0_edge = false;
													break;
												}
											}
											
											if (add_lay0_edge == true) {
												out_data->result_layer0_data[result_layer0_offset].x = 0;
												out_data->result_layer0_data[result_layer0_offset].y = 0;
												out_data->result_layer0_data[result_layer0_offset].r = 1;
												out_data->result_layer0_data[result_layer0_offset].type = e.type;
												out_data->result_layer0_data[result_layer0_offset].node.offset = e.node.offset;
												out_data->result_layer0_data[result_layer0_offset].node.layer = e.node.layer;
											
												result_layer0_offset++;
												edges_added_count++;
											}
										}									
									} else {
										// in case we do not have memory yet we can only use to_layer0_edges_count
										edges_added_count = to_layer0_edges_count;
									}
									
									// output how many layer0 edges (should be/has been) written
									out_data->result_layer0_count += edges_added_count;
								}
							}
						}
					}
				}
			}
		}
	}	
	
	return best_part;
}


inline
bool r_check_neighbor(const __global ocl_layer1_data* img_s_nodes,
					const __global ocl_layer1_data_coordinates* img_coord,
					const uint2 img_size,
					const int x, const int y, float RESPONSE) {
	bool is_valid;
	int neighbor_pos = y * img_size.x  + x;
	if (img_coord[neighbor_pos].size > 0) {
		is_valid = img_s_nodes[img_coord[neighbor_pos].offset].response.R < RESPONSE;
	} else {
		is_valid = true;
	}
	return is_valid;
}

inline
bool g_check_neighbor(const __global ocl_layer1_data* img_s_nodes,
					const __global ocl_layer1_data_coordinates* img_coord,
					const uint2 img_size,
					const int x, const int y, float RESPONSE) {
	bool is_valid;
	int neighbor_pos = y * img_size.x  + x;
	if (img_coord[neighbor_pos].size > 0) {
		is_valid = img_s_nodes[img_coord[neighbor_pos].offset].response.G < RESPONSE;
	} else {
		is_valid = true;
	}
	return is_valid;
}

__kernel
void inhibit_result(const __global ocl_layer1_data* img_s_nodes,
					const __global ocl_layer1_data_coordinates* img_coord,
					__global ocl_layer1_data_coordinates* img_coord_inhib,
					const uint2 img_size,
					__global int *error) {
	
	// CAN ONLY EXECUTE IN 1 dimension
	if (get_work_dim() != 1) {
		*error = 1;
		return; // error
	}
	
	
	// if using more work-items then there is parts to be processed then just skip invalid ones
	if (get_global_id(0) < img_size.x * img_size.y && *error == 0) {	

		int offset = img_coord[get_global_id(0)].offset;
		bool is_valid = false;
		
		// check if we have valid part on this position
		if (img_coord[get_global_id(0)].size > 0) {
			// get part from memory (all we actualy need is its position and G_RESPONSE (or R_RESPONSE if layer0)
			int x = img_s_nodes[offset].x;
			int y = img_s_nodes[offset].y;
			
			if (img_s_nodes[offset].z == 0) {
				float RESPONSE = img_s_nodes[offset].response.R;
				
				// now check for each part in local neighborhood if has better R_RESPONSE
				is_valid = true;
				is_valid = is_valid && r_check_neighbor(img_s_nodes, img_coord, img_size, x - 1, y - 1, RESPONSE);
				is_valid = is_valid && r_check_neighbor(img_s_nodes, img_coord, img_size, x + 0, y - 1, RESPONSE);
				is_valid = is_valid && r_check_neighbor(img_s_nodes, img_coord, img_size, x + 1, y - 1, RESPONSE);
				is_valid = is_valid && r_check_neighbor(img_s_nodes, img_coord, img_size, x - 1, y + 0, RESPONSE);
				is_valid = is_valid && r_check_neighbor(img_s_nodes, img_coord, img_size, x + 1, y + 0, RESPONSE);
				is_valid = is_valid && r_check_neighbor(img_s_nodes, img_coord, img_size, x - 1, y + 1, RESPONSE);
				is_valid = is_valid && r_check_neighbor(img_s_nodes, img_coord, img_size, x + 0, y + 1, RESPONSE);
				is_valid = is_valid && r_check_neighbor(img_s_nodes, img_coord, img_size, x + 1, y + 1, RESPONSE);
			} else {
				float RESPONSE = img_s_nodes[offset].response.G;
				
				// now check for each part in local neighborhood if has better G_RESPONSE
				is_valid = true;
				is_valid = is_valid && g_check_neighbor(img_s_nodes, img_coord, img_size, x - 1, y - 1, RESPONSE);
				is_valid = is_valid && g_check_neighbor(img_s_nodes, img_coord, img_size, x + 0, y - 1, RESPONSE);
				is_valid = is_valid && g_check_neighbor(img_s_nodes, img_coord, img_size, x + 1, y - 1, RESPONSE);
				is_valid = is_valid && g_check_neighbor(img_s_nodes, img_coord, img_size, x - 1, y + 0, RESPONSE);
				is_valid = is_valid && g_check_neighbor(img_s_nodes, img_coord, img_size, x + 1, y + 0, RESPONSE);
				is_valid = is_valid && g_check_neighbor(img_s_nodes, img_coord, img_size, x - 1, y + 1, RESPONSE);
				is_valid = is_valid && g_check_neighbor(img_s_nodes, img_coord, img_size, x + 0, y + 1, RESPONSE);
				is_valid = is_valid && g_check_neighbor(img_s_nodes, img_coord, img_size, x + 1, y + 1, RESPONSE);
			}
		}	
		
		// if this part is local max then save it to img_coord_inhib else set values to 0		
		img_coord_inhib[get_global_id(0)].size = is_valid ? 1 : 0;
		img_coord_inhib[get_global_id(0)].offset = is_valid ? offset : 0;
	}
}
