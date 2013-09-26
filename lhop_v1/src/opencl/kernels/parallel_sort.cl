// include files
#include <utils.cl>


/**
 * NVIDIA Bitonic Sort implementation
 */

inline void ComparatorPrivate(
    ocl_coord_sorting_pointer *keyA,
	ocl_layer1_data_coordinates *valA,
    ocl_coord_sorting_pointer *keyB,
	ocl_layer1_data_coordinates *valB,
    uint arrowDir
){
    if( (keyA->sort_value > keyB->sort_value) == arrowDir ){
        ocl_coord_sorting_pointer t; 
		t = *keyA; *keyA = *keyB; *keyB = t;
		ocl_layer1_data_coordinates v;
		v = *valA; *valA = *valB; *valB = v;
    }
}

inline void ComparatorLocal(
    __local ocl_coord_sorting_pointer *keyA,
	__local ocl_layer1_data_coordinates *valA,
    __local ocl_coord_sorting_pointer *keyB,
	__local ocl_layer1_data_coordinates *valB,
    uint arrowDir
){
    if( (keyA->sort_value > keyB->sort_value) == arrowDir ){
        ocl_coord_sorting_pointer t; 
		t = *keyA; *keyA = *keyB; *keyB = t;
		ocl_layer1_data_coordinates v;
		v = *valA; *valA = *valB; *valB = v;
    }
}
////////////////////////////////////////////////////////////////////////////////
// Bitonic sort kernel for large arrays (not fitting into local memory)
////////////////////////////////////////////////////////////////////////////////
//Bottom-level bitonic sort
__kernel void bitonicSortLocal(__global ocl_coord_sorting_pointer* buffer_keys,
								__global ocl_layer1_data_coordinates* buffer_vals,
								__local ocl_coord_sorting_pointer* local_buffer_keys,
								__local ocl_layer1_data_coordinates* local_buffer_vals,
								const uint local_size
){
//	if (get_global_id(0) == 0)
//		printf("bitonicSortLocal org: global size: %d\n", get_global_size(0));
	buffer_keys += get_group_id(0) * local_size + get_local_id(0);
	buffer_vals += get_group_id(0) * local_size + get_local_id(0);

	local_buffer_keys[get_local_id(0) +				   0] = buffer_keys[			   0];
	local_buffer_keys[get_local_id(0) + (local_size / 2)] = buffer_keys[(local_size / 2)];
	local_buffer_vals[get_local_id(0) +				   0] = buffer_vals[			   0];
	local_buffer_vals[get_local_id(0) + (local_size / 2)] = buffer_vals[(local_size / 2)];

    uint comparatorI = get_global_id(0) & ((local_size / 2) - 1);

    for(uint size = 2; size < local_size; size <<= 1){
        //Bitonic merge
        uint dir = (comparatorI & (size / 2)) != 0;
        for(uint stride = size / 2; stride > 0; stride >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            uint pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));
            ComparatorLocal(&local_buffer_keys[pos + 0], &local_buffer_vals[pos + 0],
							&local_buffer_keys[pos + stride], &local_buffer_vals[pos + stride], dir);
        }
    }

    //Odd / even arrays of LOCAL_SIZE_LIMIT elements
    //sorted in opposite directions
    {
        uint dir = (get_group_id(0) & 1);
        for(uint stride = local_size / 2; stride > 0; stride >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            uint pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));
            ComparatorLocal(&local_buffer_keys[pos + 0], &local_buffer_vals[pos + 0],
							&local_buffer_keys[pos + stride], &local_buffer_vals[pos + stride], dir);
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);
	buffer_keys[               0] = local_buffer_keys[get_local_id(0) +                0];
    buffer_keys[(local_size / 2)] = local_buffer_keys[get_local_id(0) + (local_size / 2)];
	buffer_vals[               0] = local_buffer_vals[get_local_id(0) +                0];
    buffer_vals[(local_size / 2)] = local_buffer_vals[get_local_id(0) + (local_size / 2)];
}

//Bitonic merge iteration for 'stride' >= LOCAL_SIZE_LIMIT
__kernel void bitonicMergeGlobal(
	__global ocl_coord_sorting_pointer* buffer_keys,
	__global ocl_layer1_data_coordinates* buffer_vals,
    int arrayLength,
    uint size,
    uint stride,
    uint sortDir
){
//	if (get_global_id(0) == 0)
//		printf("bitonicMergeGlobal org: size: %d, stride: %d, sort_dir %d, array_len: %d, global_size: %d\n", size, stride, sortDir, arrayLength, get_global_size(0));
    uint global_comparatorI = get_global_id(0);
    uint        comparatorI = global_comparatorI & ((uint)arrayLength / 2 - 1);

    //Bitonic merge
    uint dir = sortDir ^ ( (comparatorI & (size / 2)) != 0 );
    uint pos = 2 * global_comparatorI - (global_comparatorI & (stride - 1));

    ocl_coord_sorting_pointer keyA = buffer_keys[pos +      0];
    ocl_coord_sorting_pointer keyB = buffer_keys[pos + stride];
	
	ocl_layer1_data_coordinates valA = buffer_vals[pos +      0];
    ocl_layer1_data_coordinates valB = buffer_vals[pos + stride];

    ComparatorPrivate(&keyA, &valA, &keyB, &valB, dir);

	buffer_keys[pos +      0] = keyA;
	buffer_keys[pos + stride] = keyB;
    buffer_vals[pos +      0] = valA;
    buffer_vals[pos + stride] = valB;    
}

//Combined bitonic merge steps for
//'size' > LOCAL_SIZE_LIMIT and 'stride' = [1 .. LOCAL_SIZE_LIMIT / 2]
__kernel void bitonicMergeLocal(
	__global ocl_coord_sorting_pointer* buffer_keys,
	__global ocl_layer1_data_coordinates* buffer_vals,
	__local ocl_coord_sorting_pointer* local_buffer_keys,
	__local ocl_layer1_data_coordinates* local_buffer_vals,
	const uint local_size,
    uint arrayLength,
    uint stride,
    uint size,
    uint sortDir
){
//	if (get_global_id(0) == 0)
//		printf("bitonicMergeLocal org: size: %d, stride: %d\n", size, stride);
	buffer_keys += get_group_id(0) * local_size + get_local_id(0);
	buffer_vals += get_group_id(0) * local_size + get_local_id(0);

	local_buffer_keys[get_local_id(0) +				   0] = buffer_keys[			   0];
	local_buffer_keys[get_local_id(0) + (local_size / 2)] = buffer_keys[(local_size / 2)];
	local_buffer_vals[get_local_id(0) +				   0] = buffer_vals[			   0];
	local_buffer_vals[get_local_id(0) + (local_size / 2)] = buffer_vals[(local_size / 2)];

    //Bitonic merge
    uint comparatorI = get_global_id(0) & ((arrayLength / 2) - 1);
    uint         dir = sortDir ^ ( (comparatorI & (size / 2)) != 0 );
    for(; stride > 0; stride >>= 1){
        barrier(CLK_LOCAL_MEM_FENCE);
        uint pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));
		ComparatorLocal(&local_buffer_keys[pos + 0], &local_buffer_vals[pos + 0],
						&local_buffer_keys[pos + stride], &local_buffer_vals[pos + stride], dir);
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    buffer_keys[               0] = local_buffer_keys[get_local_id(0) +                0];
    buffer_keys[(local_size / 2)] = local_buffer_keys[get_local_id(0) + (local_size / 2)];
	buffer_vals[               0] = local_buffer_vals[get_local_id(0) +                0];
    buffer_vals[(local_size / 2)] = local_buffer_vals[get_local_id(0) + (local_size / 2)];
}


/**
 * Normal implementation of merge sort within each position of result.
 * Results in different position are executed in parallel but within each position merge sort is done in non-parallel.
 */

inline
bool compare_parts(__global ocl_layer1_data* d1,
					__global ocl_layer1_data* d2) {
	//if (d1->z == 0) 
	//	return d1->response.R > d2->response.R;
	//else
	//	return ocl_data_get_val(*d1) > ocl_data_get_val(*d2);
	if (ocl_data_get_val(*d1) == ocl_data_get_val(*d2))
		return d1->response.R > d2->response.R;
	else
		return ocl_data_get_val(*d1) > ocl_data_get_val(*d2);
}


__kernel 
void merge_sort_within_pos_result(__global ocl_layer1_data* img_s_nodes,
									const int img_s_nodes_size,
									__global ocl_layer1_data_coordinates* img_coord,
									const int2 img_size,									
									__global ocl_layer1_data* temp_data,
									__global int *error) {
	// CAN ONLY EXECUTE IN 1 dimension
	if (get_work_dim() != 1) {
		*error = 1;
		return; // error
	}
	
	if (get_global_id(0) < img_s_nodes_size && img_coord[get_global_id(0)].size > 1) {
	//if (get_global_id(0) < img_s_nodes_size) {

		//int coord_pos = getImgPos(img_s_nodes[get_global_id(0)].y, img_size.x, img_s_nodes[get_global_id(0)].x);
		
		//if (img_coord[coord_pos].size > 1) 
		{
			
			// count how many parts we have
			//int size = img_coord[coord_pos].size;
			int size = img_coord[get_global_id(0)].size;
			
			// get position for this part and its temp data (this can be outside of if clause due to possible better memory access)
			__global ocl_layer1_data* parts = img_s_nodes + img_coord[get_global_id(0)].offset;
			__global ocl_layer1_data* parts_temp_data = temp_data + img_coord[get_global_id(0)].offset;
				
			//__global ocl_layer1_data* parts = img_s_nodes + img_coord[coord_pos].offset;
			//__global ocl_layer1_data* parts_temp_data = temp_data + img_coord[coord_pos].offset;
			
			if (size > 1) {
				
				int sort_layer_count = 1;
				while (sort_layer_count < size) {
					
					for (int offset = 0; offset < size; offset+=2*sort_layer_count) {
						
						// merge two sorted arrays into one (save to tmp and then copy it back)
						int left = offset;
						int right = min(offset + sort_layer_count, size);
						
						int left_max_offset = min(left + sort_layer_count, size);
						int right_max_offset =  min(right + sort_layer_count, size);
						
						for (int i = 0; i < 2*sort_layer_count; i++) {
							if (left < left_max_offset && (right >= right_max_offset || compare_parts(&parts[left], &parts[right]))) {
								parts_temp_data[i] = parts[left];
								left++;
							} else if (right < right_max_offset) {
								parts_temp_data[i] = parts[right];
								right++;
							}
						}
						
						int max_i = offset + 2*sort_layer_count <= size ? 2*sort_layer_count : size - offset;
						for (int i = 0; i < max_i; i++) {
							parts[offset + i] = parts_temp_data[i];
						}		
					}
					sort_layer_count *= 2;
				}
			}
		}
	}
}

__kernel 
void merge_sort_within_pos_result_vlen(__global ocl_layer1_data* img_s_nodes,
									__global int* img_s_nodes_size,
									__global ocl_layer1_data_coordinates* img_coord,
									const int2 img_size,
									__global ocl_layer1_data* temp_data,
									__global int *error) {
	merge_sort_within_pos_result(img_s_nodes, *img_s_nodes_size, img_coord, img_size, temp_data, error);
}



/**
 * Parallel implementation fo merge sort (NOT IMPLEMENTED)
 *

__kernel 
void parallel_merge_sort_within_pos_result(__global ocl_layer1_data* img_s_nodes,
											const int img_s_nodes_size,
											__global ocl_layer1_data_coordinates* img_coord,
											const int2 img_size,									
											__global ocl_layer1_data* temp_data,
											__global int *error) {
	// CAN ONLY EXECUTE IN 1 dimension
	if (get_work_dim() != 1) {
		*error = 1;
		return; // error
	}
	
	if (get_global_id(0) < img_s_nodes_size) {

		//int coord_pos = getImgPos(img_s_nodes[get_global_id(0)].y, img_size.x, img_s_nodes[get_global_id(0)].x);
		
		//if (img_coord[coord_pos].size > 1) 
		{
			
			// count how many parts we have
			//int size = img_coord[coord_pos].size;
			int size = img_coord[get_global_id(0)].size;
			
			// get position for this part and its temp data (this can be outside of if clause due to possible better memory access)
			__global ocl_layer1_data* parts = img_s_nodes + img_coord[get_global_id(0)].offset;
			__global ocl_layer1_data* parts_temp_data = temp_data + img_coord[get_global_id(0)].offset;
				
			//__global ocl_layer1_data* parts = img_s_nodes + img_coord[coord_pos].offset;
			//__global ocl_layer1_data* parts_temp_data = temp_data + img_coord[coord_pos].offset;
			
			if (size > 10) {
				
				int sort_layer_count = 1;
				while (sort_layer_count < size) {
					
					for (int offset = 0; offset < size; offset+=2*sort_layer_count) {
						
						// merge two sorted arrays into one (save to tmp and then copy it back)
						int left = offset;
						int right = min(offset + sort_layer_count, size);
						
						int left_max_offset = min(left + sort_layer_count, size);
						int right_max_offset =  min(right + sort_layer_count, size);
						
						for (int i = 0; i < 2*sort_layer_count; i++) {
							if (left < left_max_offset && (right >= right_max_offset || compare_parts(&parts[left], &parts[right]))) {
								parts_temp_data[i] = parts[left];
								left++;
							} else if (right < right_max_offset) {
								parts_temp_data[i] = parts[right];
								right++;
							}
						}
						
						int max_i = offset + 2*sort_layer_count <= size ? 2*sort_layer_count : size - offset;
						for (int i = 0; i < max_i; i++) {
							//parts[offset + i].response.R += 2;
							parts[offset + i] = parts_temp_data[i];
						}		
					}
					sort_layer_count *= 2;
				}
			}
		}
	}
}
*/