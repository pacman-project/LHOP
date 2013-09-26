
#ifndef _UTILS_CL_H_
#define _UTILS_CL_H_

#ifndef CL_VERSION_1_1

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable

#endif 

#pragma OPENCL EXTENSION cl_amd_printf : enable

#define getGlobalID() get_global_id(0) * (get_global_size(1) * get_global_size(2)) + get_global_id(1) * get_global_size(2) + get_global_id(2)

// we cannot use mad24(y, size, x) for img pos since it would limit ous to image of size 2048 * 2048 
#define getImgPos(y, size, x) ((y) * (size)) + x 

#ifndef printf 
	#define printf(fm, ...) ; 
#endif


#define M_SQRT1_2_E 0.707107
#define M_PI_E 3.142857

#ifndef NULL
#define NULL 0
#endif


#define _vector_n_(name, n) name##n
#define _vector_n(name, n) _vector_n_(name,n)


/*
 * Const index keys for attributes. (they must match with values in host source)
 */

#define MULTIPLY_CONVOLUTION 0
#define LAYER1_TYPE_STRUCT 1
#define LAYER1_TYPE_APP 2
#define LAYER1_TYPE_DOG 3
#define LAYER1_TYPE_LOGGABOR 4

#define ATTR_VISITED 1024
#define ATTR_CHECKED 2048
#define NODE_DELETED_ATTR 0x8000


#define NODE_REDUCED_ATTR 1
#define HAS_NEXT_LAYER 4
#define FROM_PREV_LAYER 64
#define IN_NODE_ATTR 0x1000
#define BOUDARY_NODE_ATTR 0x2000
#define OUT_NODE_ATTR 0x4000

#define IMG_NODE_ATTR 2
#define GRID_NODE_ATTR 8

#define MAX_LAYER_NUMBER 16

#define ALLOW_UPDATE_ATTR 0x40000

#define TEXTURE_NODE 128

/**
 * Structures for layer1_data, part_data and their edges
 */
#define OCL_LYR_CENTER_BACK_EDGE 0
#define OCL_LYR_SRC_EDGE 1
#define OCL_LYR_FORBIDDEN_EDGE 2
#define OCL_TO_PREV_LAYER_EDGE 3
#define OCL_TO_SIMILAR_EDGE 4
#define OCL_TO_LAYER_0 5

// must eksplicitly define how many edge types we use
#define OCL_EDGE_COUNT 6

// opencl layer1_data (only information for part is here) (68 bytes)
typedef struct {
	int x, y, z;
	struct {
		float R, RR, G;
	} response;
	int m;

	uint attr;
	// start position and size of elements in global
	// opencl memory of connections (edge) for same layer
	struct {
		int offset, size;
	} edges_loc[OCL_EDGE_COUNT];
} ocl_layer1_data;

// quick replacement for function where val can be quickly replaced with any value from ocl_layer1_data
#define ocl_data_get_val(data) ((data).response.G)
#define ocl_data_get_vval(data) ((data).z == 0 ? (data).response.R : (data).response.G)

#define OCL_EDGE_DATA_EMPTY_TYPE 0
#define OCL_EDGE_DATA_IP2_TYPE 1

// opencl edge_data_ip2 (16 bytes)
typedef struct {
	int x, y;
	float r;
	int type; // type of edge connection (0 - empty, 1 - edge_data_ip2)
	// information for part_data that this is pointing to 
	struct {
		int offset;
		int layer;
	} node;
} ocl_edge_data_ip2;

// struction to hold pointer to parts for each image position (x,y) (2 bytes)
typedef struct { 
	int offset, size;
} ocl_layer1_data_coordinates;

// pointer for sorting of coordinates
typedef struct {
	float sort_value; // usually R_RESPONSE value
	int coord_ind;	// usually original coordinate index value (when used in crate_layer0.cl) or pointer to where candidates for this position are curently saved (when used in add_layer1.cl)
} ocl_coord_sorting_pointer;

// opencl part_data (72 bytes)
typedef struct {
	int cmx, cmy;
	int layer;
	int type;
	struct {
		float R, RR, G;
	} thresh;

	int bpcount; // number of "basic" parts
	uint attr;
	// start position and size of elements in global
	// opencl memory of connections (edge) for same layer
	struct {
		int offset, size;
	} edges_loc[OCL_EDGE_COUNT];
} ocl_part_data;

#define OCL_EDGE_DATA_EMPTY_TYPE 0
#define OCL_EDGE_DATA_PART_DATA_2_TYPE 1
#define OCL_EDGE_DATA_PART_DATA_2R_TYPE 2
#define OCL_EDGE_DATA_PART_DATA_2C_TYPE 3

// max size of distribution (WARNING: when adjusting this constant do not forget same one in cl_structures.h)
#define OCL_EDGE_DATA_2_DISTRIBUTION_MAX_SIZE 9*9 

// opencl part_data_2 for edges (136 bytes)
typedef struct {

	// information about values that are encoded 
	int cast_type; 

	int x, y;
	int type;
	int distr_size;
	float distr[OCL_EDGE_DATA_2_DISTRIBUTION_MAX_SIZE];
	struct {
		float mean, variance;
	} gdistr;

	struct {
		int offset, size;
	} app; // pointer to ocl_app_data

	// information for part_data that this is pointing to 
	struct {
		int offset;
		int layer;
	} node;
} ocl_part_data_2;

typedef struct {
	int type;
	float value;
} ocl_app_data;

typedef struct {
	uint lib_part_offset;
	uint img_part_offset;
	uint skip;
	float r_sum;
	float realization_ratio;	
	float g_prod;
	int x,y;
	int new_x,new_y;
	int result_edges_prev_lay_size;
	int result_edges_layer0_size;
	int part_type;
	int img_part_type;
	uint max_schur_mem_offset;
} layer_candidate;

#define ocl_layer_candidate_get_val(candidate) ((candidate).g_prod)

#define ocl_new_attr(attr, a) *attr = a
#define ocl_set_attr(attr, a) *attr |= a
#define ocl_clear_attr(attr, a) *attr &= ~a
#define ocl_is_attr_set(attr, a) ((attr & a) == a)
#define ocl_is_attr_set2(attr, a)  ((attr & a) != 0)


#define SCHUR_PROD_USE_VAL_RESPONSE 0
#define SCHUR_PROD_USE_IDENTITY_RESPONSE 1
#define SCHUR_PROD_USE_R_RESPONSE 2
#define SCHUR_PROD_USE_G_DISTR_RESPONSE 3
#define SCHUR_PROD_USE_SIMPLE_G_DISTR_RESPONSE 4


inline
void copy_to_local_data(__global float* image, int2 image_size, 
						__local float* local_data, int2 local_data_size,
						int2 get_global_id_, int2 get_local_id_,
						int2 copy_offset, int2 copy_size) {
	
	for (int j = 0; j < copy_size.y; j++) {
		for (int i = 0; i < copy_size.x; i++) {
			int src_offset = getImgPos( (int)get_global_id_.y + copy_offset.y + j, image_size.x, (int)get_global_id_.x + copy_offset.x + i );
			int dst_offset = getImgPos( (int)get_local_id_.y + copy_offset.y + j, local_data_size.x, (int)get_local_id_.x + copy_offset.x + i );
			
			local_data[dst_offset] = image[src_offset];
		}
	}
}

/**
 * Copies from global rectangular region to local data based on global and local id.
 * For one local work-item it loads only one memory chunk (of size vector_size) and
 * waits for all other work-items in same work gorup (barrier(CLK_LOCAL_MEM_FENCE);).
 * It can also load additional border data based on values in extra_border_load:
 *		extra_border_load.x = right border size
 *		extra_border_load.y = lower border size
 *		extra_border_load.z = left border size 
 *		extra_border_load.w = upper border size
 */ 
inline
void load_local_data(__global float* global_data, int2 global_data_size, 
						__local float* local_data, int2 local_data_size,
						int2 get_global_id_, int2 get_local_id_, 
						int vector_size, int4 extra_border_load) {
	
	//extra_border_load.x - right
	//extra_border_load.y - lower
	//extra_border_load.z - left 
	//extra_border_load.w - upper

	int2 copy_offset = (int2)(0,0);
	int2 copy_size = (int2)(0,0);

	copy_offset = (int2)(0,0);
	copy_size = (int2)(vector_size, 1);
	copy_to_local_data(global_data, global_data_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
					
	// if upper border
	if (extra_border_load.w > 0 && get_local_id(1) == 0) {
		copy_offset = (int2)(0, -extra_border_load.w);
		copy_size = (int2)(vector_size, extra_border_load.w);
		copy_to_local_data(global_data, global_data_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
		
		// if upper left corrner
		if (extra_border_load.z > 0 && get_local_id(0) == 0) {
			copy_offset = (int2)(- extra_border_load.z, -extra_border_load.w);
			copy_size = (int2)(extra_border_load.z, extra_border_load.w);
			copy_to_local_data(global_data, global_data_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
		}		
		
		// if upper right corrner
		if (extra_border_load.x > 0 && (get_local_id(0) == get_local_size(0) - 1 || get_global_id_.x + vector_size >= global_data_size.x - extra_border_load.x)) {
			copy_offset = (int2)(vector_size, -extra_border_load.w);
			copy_size = (int2)(extra_border_load.x, extra_border_load.w);
			copy_to_local_data(global_data, global_data_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
		}
	}

	// if lower border
	if (extra_border_load.y > 0 && (get_local_id(1) == get_local_size(1) - 1 || get_global_id_.y + vector_size >= global_data_size.y - extra_border_load.y)) {
		copy_offset = (int2)(0,1);
		copy_size = (int2)(vector_size, extra_border_load.y);
		copy_to_local_data(global_data, global_data_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);

		// if lower left corrner
		if (extra_border_load.z > 0 && get_local_id(0) == 0) {
			copy_offset = (int2)(- extra_border_load.z,1);
			copy_size = (int2)(extra_border_load.z, extra_border_load.y);
			copy_to_local_data(global_data, global_data_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
		}
		
		// if lower right corrner
		if (extra_border_load.x > 0 && (get_local_id(0) == get_local_size(0) - 1 || get_global_id_.x + vector_size >= global_data_size.x - extra_border_load.x)) {
			copy_offset = (int2)(vector_size,1);
			copy_size = (int2)(extra_border_load.x , extra_border_load.y);
			copy_to_local_data(global_data, global_data_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
		}
	}

	// if left border				
	if (extra_border_load.z > 0 && get_local_id(0) == 0) {
		copy_offset = (int2)(- extra_border_load.z, 0);
		copy_size = (int2)(extra_border_load.z, 1);
		copy_to_local_data(global_data, global_data_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
	}		
	
	// if right border 
	if (extra_border_load.x > 0 && (get_local_id(0) == get_local_size(0) - 1 || get_global_id_.x + vector_size >= global_data_size.x - extra_border_load.x)) {
		copy_offset = (int2)(vector_size, 0);
		copy_size = (int2)(extra_border_load.x, 1);
		copy_to_local_data(global_data, global_data_size, local_data, local_data_size, get_global_id_, get_local_id_, copy_offset, copy_size);
	}			
	
	barrier(CLK_LOCAL_MEM_FENCE);
	
}

__kernel
void set_memory_zero(__global int16* memory, int size, int unaligned_bytes) {
	for (int get_global_id_0 = get_global_id(0); get_global_id_0 < size; get_global_id_0 += get_local_size(0)) {
		memory[get_global_id_0] = (int16)0;
	}
	if (get_global_id(0) == get_global_size(0) - 1 && unaligned_bytes > 0) {
		__global int* unalinged_memory = (__global int*)&memory[size];
		// use constant size in for loop so compiler can unroll it if it has this capability
		for (int i = 0 ; i < sizeof(int16); i++) {
			if (i < unaligned_bytes) {
				unalinged_memory[i] = 0;
			}
		}
	}	
}


__kernel
void set_memory(__global int* memory, int value, int size) {
	if (get_global_id(0) < size) {
		memory[get_global_id(0)] = value;
	}
}

#endif
