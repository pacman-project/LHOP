// include files
#include <utils.cl>

#if __DEVICE_LOCAL_MEM_TYPE__ > 0
	#define HAS_FAST_LOCAL_MEM 1
#endif

__kernel
void make_2x2sum(__global float* input, int2 image_size, int2 input_offset,
					__local float* local_data, int2 local_data_size,
					__global float* output, __global int* error) {
	
	// only two dimensions are allowed
	if (get_work_dim() != 2) {
		*error = 1;
		return;
	}
	
	#define VECTOR_SIZE 1
	int get_global_id_1 = get_global_id(1);
	int get_global_id_0 = get_global_id(0) * VECTOR_SIZE;

	
	if ((int)get_global_id_1 >= 0 && (int)get_global_id_1 < image_size.y - 1 &&
		(int)get_global_id_0 >= 0 && (int)get_global_id_0 < image_size.x - 1 && *error == 0) {


	// if device has fast local memory we can use it (CPU uses global memory for local mem so this would only slow it down
// (HAS_FAST_LOCAL_MEM should be set during compile time in HOST run-time)
#ifdef HAS_FAST_LOCAL_MEM
	
		int get_local_id_1 = get_local_id(1);
		int get_local_id_0 = get_local_id(0) * VECTOR_SIZE;
		
		int2 get_global_id_ = (int2)(get_global_id_0, get_global_id_1) + input_offset;
		int2 get_local_id_ = (int2)(get_local_id_0, get_local_id_1);
		
		// load local data with additionl one pixel on lower and right border
		load_local_data(input, image_size, local_data, local_data_size, get_global_id_, get_local_id_, VECTOR_SIZE, (int4)(1,1,0,0));
	
		int img_local_pos_offset_start = getImgPos( ((int)get_local_id_1), local_data_size.x , (int)get_local_id_0);

		float sum = 0;
		
		sum += local_data[img_local_pos_offset_start + 0];
		sum += local_data[img_local_pos_offset_start + 1];
		
		sum += local_data[img_local_pos_offset_start + max(0, local_data_size.x - 1) + 1];
		sum += local_data[img_local_pos_offset_start + max(0, local_data_size.x - 1) + 2];
	
#else
		int img_global_pos_offset_start = getImgPos( ((int)get_global_id_1  + input_offset.y), image_size.x, (int)get_global_id_0 + input_offset.x);
		
		float sum = 0;
		
		sum += input[img_global_pos_offset_start + 0];
		sum += input[img_global_pos_offset_start + 1];
		
		sum += input[img_global_pos_offset_start + max(0, image_size.x - 1) + 1];
		sum += input[img_global_pos_offset_start + max(0, image_size.x - 1) + 2];
#endif
		
		// save value back to global memory
		int pos = getImgPos( (int)get_global_id_1 + 1, image_size.x, (int)get_global_id_0 + 1);
		output[pos] = sum;
		
	} else {
#ifdef HAS_FAST_LOCAL_MEM
		// this barrier should be here since AMD implementation will not let threads finish unless all threads in same block get to barrier
		barrier(CLK_LOCAL_MEM_FENCE);
#endif
	}

	#undef VECTOR_SIZE
}

__attribute__((vec_type_hint(float4)))
__kernel
void get_max(__global float* input, int input_size, int offset,
			__local float* local_data, int local_data_size,
			__global float* result, __global int* error) {

	#define VECTOR_SIZE 4

	#define vector_float _vector_n(float,VECTOR_SIZE)	// this should produce float16, float8, float4 or float2 - based on value of VECTOR_SIZE
	#define vector_load _vector_n(vload,VECTOR_SIZE)

	if (get_work_dim() != 1) {
		*error = 1;
		return;
	}

	// local_data_size must be at least of size VECTOR_SIZE
	if (local_data_size < VECTOR_SIZE) {
		*error = 2;
		return;
	}

	if (*error == 0) {

		vector_float max_value = (vector_float)-MAXFLOAT;	
		
		__global float* input_offset = input + offset;

		int load_pos;
		for (load_pos = get_global_id(0); load_pos < input_size/VECTOR_SIZE; load_pos+=local_data_size) {
		//for (int i = 0; i < input_size/VECTOR_SIZE; i+=local_data_size) {
			//if ((int)get_global_id(0) + i >= 0 && (int)get_global_id(0) + i  < input_size) 
			{
				max_value = fmax(max_value, vector_load(load_pos, input_offset));
				//max_value = fmax(max_value, vector_load(get_global_id(0) + i, input_offset));
			}
		}
		
		// additionaly check for all values outside of input_size/VECTOR_SIZE (in case input_size is not dividable by VECTOR_SIZE)
		//for (int i = VECTOR_SIZE * (input_size/VECTOR_SIZE); i < input_size; i++) {
		for (int i = load_pos * VECTOR_SIZE; i < input_size; i++) {
			max_value.s0 = fmax(max_value.s0, input_offset[i]);
		}
		
		max_value.s0 = fmax(max_value.s0, max_value.s1);
		max_value.s0 = fmax(max_value.s0, max_value.s2);
		max_value.s0 = fmax(max_value.s0, max_value.s3);
#if VECTOR_SIZE > 4
		max_value.s0 = fmax(max_value.s0, max_value.s4);
		max_value.s0 = fmax(max_value.s0, max_value.s5);
		max_value.s0 = fmax(max_value.s0, max_value.s6);
		max_value.s0 = fmax(max_value.s0, max_value.s7);
	#if VECTOR_SIZE > 8	
		max_value.s0 = fmax(max_value.s0, max_value.s8);
		max_value.s0 = fmax(max_value.s0, max_value.s9);
		max_value.s0 = fmax(max_value.s0, max_value.sA);
		max_value.s0 = fmax(max_value.s0, max_value.sB);
		max_value.s0 = fmax(max_value.s0, max_value.sC);
		max_value.s0 = fmax(max_value.s0, max_value.sD);
		max_value.s0 = fmax(max_value.s0, max_value.sE);
	#endif
#endif
		
		local_data[get_local_id(0)] = max_value.s0;
		
		barrier(CLK_LOCAL_MEM_FENCE);
		
		if (get_local_id(0) == get_local_size(0) -1) {
			
		
			max_value = (vector_float)-MAXFLOAT;
			
			for (int i = 0; i < local_data_size/VECTOR_SIZE; i++) {
				max_value = fmax(max_value, vector_load(i, local_data));
			}
			
			// additionaly check for all values outside of local_data_size/VECTOR_SIZE (in case local_data_size is not dividable by VECTOR_SIZE)
			for (int i = VECTOR_SIZE * (local_data_size/VECTOR_SIZE); i < local_data_size; i++) {
				max_value.s0 = fmax(max_value.s0, local_data[i]);
			}
					
			max_value.s0 = fmax(max_value.s0, max_value.s1);
			max_value.s0 = fmax(max_value.s0, max_value.s2);
			max_value.s0 = fmax(max_value.s0, max_value.s3);
#if VECTOR_SIZE > 4		
			max_value.s0 = fmax(max_value.s0, max_value.s4);
			max_value.s0 = fmax(max_value.s0, max_value.s5);
			max_value.s0 = fmax(max_value.s0, max_value.s6);
			max_value.s0 = fmax(max_value.s0, max_value.s7);
	#if VECTOR_SIZE > 8	
			max_value.s0 = fmax(max_value.s0, max_value.s8);
			max_value.s0 = fmax(max_value.s0, max_value.s9);
			max_value.s0 = fmax(max_value.s0, max_value.sA);
			max_value.s0 = fmax(max_value.s0, max_value.sB);
			max_value.s0 = fmax(max_value.s0, max_value.sC);
			max_value.s0 = fmax(max_value.s0, max_value.sD);
			max_value.s0 = fmax(max_value.s0, max_value.sE);
	#endif		
#endif			
			*result = max_value.s0;
		}
	}
		
	#undef vector_float
	#undef vector_load
	#undef VECTOR_SIZE
}

__kernel
void prepare_results(__global float* filtered_images, __global float* input_max_imgs, __global float* input_max_imgs_index, __global float* input_sum2x2, int2 image_size, int2 new_image_size, // input data and size of it (all input data should be of same size)					
					__global float* input_image_mask, int2 input_image_mask_size,					
					__global float* maximum_img_max, __global float* maximum_sum2x2, 					
					__local float* local_data_max_imgs, __local float* local_data_sum2x2, int2 local_data_size,											
					float layer1_threshold, int layer1_3x3bound, int number_images, int border, float response_percent,
					__global ocl_layer1_data_coordinates* output, __global ocl_coord_sorting_pointer* output_keys,
					__global int* number_created, __global int* number_positions_used,
					__global int* error) {
	

	// only two dimensions are allowed
	if (get_work_dim() != 2) {
		*error = 1;
		return;
	}
	
	int check_neighboor_size = 3;
	int check_neighboor_size_half = check_neighboor_size/2;
	
	int get_global_id_1 = get_global_id(1);
	int get_global_id_0 = get_global_id(0);
	
	//int get_global_id_1 = check_neighboor_size_half + get_global_id(1);
	//int get_global_id_0 = check_neighboor_size_half + get_global_id(0);

	if ((int)get_global_id_1 >= 0 && (int)get_global_id_1 < image_size.y &&
		(int)get_global_id_0 >= 0 && (int)get_global_id_0 < image_size.x && *error == 0) {
					
		bool create_part = false;
		int create_parts_count = 0;
		
		if ((int)get_global_id_1 >= check_neighboor_size_half && (int)get_global_id_1 < image_size.y - check_neighboor_size_half &&
			(int)get_global_id_0 >= check_neighboor_size_half && (int)get_global_id_0 < image_size.x - check_neighboor_size_half) {

	// if device has fast local memory we can use it (CPU uses global memory for local mem so this would only slow it down
	// (HAS_FAST_LOCAL_MEM should be set during compile time in HOST run-time)
	#ifdef HAS_FAST_LOCAL_MEM

			int get_local_id_1 = check_neighboor_size_half + get_local_id(1);
			int get_local_id_0 = check_neighboor_size_half + get_local_id(0);
			
			int2 get_global_id_ = (int2)(get_global_id_0, get_global_id_1);
			int2 get_local_id_ = (int2)(get_local_id_0, get_local_id_1);
			
			// load local data for max_image and sum2x2 image with additional borders of one pixel size on each side
			load_local_data(input_max_imgs, image_size, local_data_max_imgs, local_data_size, get_global_id_, get_local_id_, 1, (int4)(1,1,1,1));
			load_local_data(input_sum2x2, image_size, local_data_sum2x2, local_data_size, get_global_id_, get_local_id_, 1, (int4)(1,1,1,1));
		
			// set position offset for local data
			int img_pos_offset_start = getImgPos( ((int)get_local_id_1 - check_neighboor_size_half), local_data_size.x , (int)get_local_id_0 - check_neighboor_size_half);
			int skipped_pos_x = max(0, local_data_size.x - check_neighboor_size);
			
			// load middle values from local data
			float mid = local_data_max_imgs[getImgPos( (int)get_local_id_1, local_data_size.x , (int)get_local_id_0)];
			float mid2 = local_data_sum2x2[getImgPos( (int)get_local_id_1, local_data_size.x , (int)get_local_id_0)];
			
	#else
			// no need to copy data only set position offset for global data
			int img_pos_offset_start = getImgPos( ((int)get_global_id_1 - check_neighboor_size_half), image_size.x, (int)get_global_id_0 - check_neighboor_size_half);
			int skipped_pos_x = max(0, image_size.x - check_neighboor_size);
			
			// load middle values from global data
			float mid = input_max_imgs[getImgPos( (int)get_global_id_1, image_size.x, (int)get_global_id_0)];
			float mid2 = input_sum2x2[getImgPos( (int)get_global_id_1, image_size.x, (int)get_global_id_0)];
	#endif

			float threshold = *maximum_img_max * layer1_threshold;
				
			float wfactor = input_image_mask_size.x == 0 ? 0 : (float)input_image_mask_size.x/image_size.x;
			float hfactor = input_image_mask_size.y == 0 ? 0 : (float)input_image_mask_size.y/image_size.y;
			
			bool is_valid_mask = input_image_mask_size.x > 0 && input_image_mask_size.y;
			bool is_pos_within_mask_bounds = get_global_id_0 < input_image_mask_size.x && get_global_id_1 < input_image_mask_size.y;
			
			int mask_pos = getImgPos( (int)round(get_global_id_1 * hfactor), input_image_mask_size.x, (int)round(get_global_id_0 * wfactor));
			
			
			// first count how many parts we will create for this position
			// do it outside of if clause for faster memory access
			int filter_pos = getImgPos((int)get_global_id_1, image_size.x,  (int)get_global_id_0);		
			
			for (int i = 0; i < number_images; i++) {			
				float value = filtered_images[filter_pos + (image_size.x * image_size.y) * i];
				// check only with mid values since only this value will go into part as R_RESPONSE
				create_parts_count += (value >= mid * response_percent ? 1 : 0);
			}
			
			if (is_valid_mask == false || (is_pos_within_mask_bounds && input_image_mask[mask_pos] > 1e-6)) {
				
				if (mid > threshold) {
					
					int sortingArrayPos = 0;
					// calculate how many of neighbors + mid in max_img are bigger then mid
					int gcount = 0;
					int img_pos_offset = max(0, img_pos_offset_start);
					
					for (int i = 0; i < check_neighboor_size; i++) {
						for (int j = 0; j < check_neighboor_size; j++) {

	#ifdef HAS_FAST_LOCAL_MEM
							gcount += local_data_max_imgs[img_pos_offset] >= mid ? 1 : 0;
	#else
							gcount += input_max_imgs[img_pos_offset] >= mid ? 1 : 0;
	#endif
							img_pos_offset++;
						}
						img_pos_offset += skipped_pos_x;
					}

					if (gcount <= layer1_3x3bound) {
						// create one node - reserve memory slot for this one
						atomic_add(number_created, create_parts_count);
						sortingArrayPos = atomic_inc(number_positions_used);
						create_part = true;
					} else {
						// calculate how many of neighbors + mid in sum2x2_img are bigger then mid
						int gcount2 = 0;
						img_pos_offset = max(0, img_pos_offset_start);
						
						for (int i = 0; i < check_neighboor_size; i++) {
							for (int j = 0; j < check_neighboor_size; j++) {

	#ifdef HAS_FAST_LOCAL_MEM
								gcount2 += local_data_sum2x2[img_pos_offset] >= mid2 ? 1 : 0;
	#else
								gcount2 += input_sum2x2[img_pos_offset] >= mid2 ? 1 : 0;
	#endif
								img_pos_offset++;							
							}
							img_pos_offset += skipped_pos_x;
						}	
						if (gcount2 <= layer1_3x3bound) {
							// create one node - reserve memory slot for this one
							atomic_add(number_created, create_parts_count);
							sortingArrayPos = atomic_inc(number_positions_used);
							create_part = true;
						} 
					}					
					
					if (create_part == true) {
						
						output_keys[sortingArrayPos].sort_value = mid;
						output_keys[sortingArrayPos].coord_ind = getImgPos((int)get_global_id_1, image_size.x, (int)get_global_id_0);
						output[sortingArrayPos].size = create_parts_count;
						output[sortingArrayPos].offset = 0;
					}
				}				
			}

		} else {
	#ifdef HAS_FAST_LOCAL_MEM
			// this barrier should be here since AMD implementation will not let threads finish unless all threads in same block get to barrier
			barrier(CLK_LOCAL_MEM_FENCE);
	#endif			
		}		
		
	}
	#undef VECTOR_SIZE
}

__kernel
void make_results(__global float* filtered_images, int2 image_size, int2 new_image_size, int number_images,
					__global ocl_layer1_data_coordinates* compact_coords, __global ocl_coord_sorting_pointer* compact_coords_keys, int compact_coords_size,
					__global ocl_layer1_data* result, const int border, const float response_percent, float power_correction, int normalization_part_index, 
					__global int* error) {
					
	// only one dimension is allowed
	if (get_work_dim() != 1) {
		*error = 1;
		return;
	}

	if (get_global_id(0) < compact_coords_size)  {
		float normalization_value = compact_coords_keys[normalization_part_index].sort_value;
		bool powercorr = (power_correction != 1.0);

		int old_pos = compact_coords_keys[get_global_id(0)].coord_ind;
		
		int x = fmod((float)old_pos , (float)image_size.x);
		int y = old_pos / image_size.x;

		// get max filtered value at this position that is used to filter out some parts
		//float mid = input_max_imgs[old_pos];
		float mid = compact_coords_keys[get_global_id(0)].sort_value;
		
		// if this image position has valid part then proceed
		if (compact_coords[get_global_id(0)].size > 0) { // this should always be true so there is no need to read this value
		
			// first get result position for this part position from pre-set offset values
			int result_offset = compact_coords[get_global_id(0)].offset;

			int save_part_offset = 0;
			// then make parts from each filtered image
			for (int i = 0; i < number_images; i++) {
				
				float r_resp = filtered_images[old_pos + (image_size.x * image_size.y) * i];
				
				if (r_resp >= mid * response_percent) {
					// create part
					ocl_layer1_data saved_part;
				
					// set indexes
					saved_part.x = x + border;
					saved_part.y = y + border;
					saved_part.z = 0;
					
					saved_part.m = i;
					
					r_resp = fmin(r_resp / normalization_value, 1.0f);
					if (powercorr) r_resp = pow(r_resp, power_correction);

					saved_part.response.R = r_resp;
					saved_part.response.RR = 1;
					saved_part.response.G = 1;
					saved_part.attr = 0;
					
					// set defualt attribute
					ocl_set_attr(&saved_part.attr, GRID_NODE_ATTR);
									
					for (int j = 0; j < OCL_EDGE_COUNT; j++) {
						saved_part.edges_loc[j].offset = 0;
						saved_part.edges_loc[j].size = 0;
					}
					
					// and save it to offset position
					result[result_offset + save_part_offset++] = saved_part;
				}
			}
			
		}
	}
	
}
