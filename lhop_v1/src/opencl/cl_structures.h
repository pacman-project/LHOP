#pragma once

#ifndef _OPENCL_STRUCTURES_
#define _OPENCL_STRUCTURES_

#include <CL/opencl.h>

/**
 * Defines for layer1_creators types (struct, app, dog ...)
 */
#define LAYER1_TYPE_STRUCT 1
#define LAYER1_TYPE_APP 2
#define LAYER1_TYPE_DOG 3
#define LAYER1_TYPE_LOGGABOR 4

/**
 * All data structures here must be compliant with C99 standard and with OpenCL specifications.
 */

// edge (connection) types (MUST start with 0 since it is used for direct array access)
#define OCL_LYR_CENTER_BACK_EDGE 0
#define OCL_LYR_SRC_EDGE 1
#define OCL_LYR_FORBIDDEN_EDGE 2
#define OCL_TO_PREV_LAYER_EDGE 3
#define OCL_TO_SIMILAR_EDGE 4
#define OCL_TO_LAYER_0 5

// must eksplicitly define how many edge types we use
#define OCL_EDGE_COUNT 6

/**
 * Structure all layer1_data must hold three memory array for each layer:
 *
 * 1. One for each layer1_data i.e. all parts found in image for specific layer:
 *		- Each part must have its positions (x,y), layer (z), part type (m) and responses (R, RR, G)
 *		- Memory object must be order by parts based on their location i.e. parts with same
 *		  location (x,y) must be in consecutive memory positions.  
 *		- Each layer1_data may have one or more (or none) connections to other layer1_data in different 
 *		  layer (i.e. in different memory object). Each connection can have different type.
 *		- For each connection type in one layer1_data there is reported position (offset) where
 *		  data connections for this type start and count (size) of connections. Values are
 *		  referred to memory object for same layer.
 *	
 *	[p1(x1,y1), p2(x1,y1), p3(x1,y1), ... | p4(x2,y2), p5(x2,y2), p6(x2,y2) ... | .. ]
 *
 * 2. One memory object for quick access to parts based on thier location (x,y) for specific layer:
 *		- Memory array MUST be of size img_height * img_width and its values represent pointers (offset, size)
 *		  to array of parts (in layer1_data memory object) for specific location and layer
 *		- Values referre only to memory object of same layer.
 *		- If there are no parts for specific location (x,y) then its size is 0 and offset value is undefined.
 *
 *  [ (offset(x1,y1), size(x1,y1)), (offset(x2,y2), size(x2,y2)), ... ]
 *
 * 2. One memory object for quick access to parts based on thier location (x,y) for inhibited parts for specific layer:
 *		- Memory array MUST be of size img_height * img_width and its values represent pointers (offset, size)
 *		  to part (in layer1_data memory object) for specific location and layer.
 *		- Values referre only to memory object of same layer.
 *		- If there are no parts for specific location (x,y) then its size is 0 and offset value is undefined.
 *		- INFO: for inhibited parts there can be only one part per one location (x,y) therfore size can only be 0 or 1
 *
 *	[ (offset(x1,y1), size(x1,y1) = { 0 | 1 } ), (offset(x2,y2), size(x2,y2) = { 0 | 1 } ), ... ]
 */

// opencl layer1_data (only information for part is here) (68 bytes)
typedef struct {
	int x, y, z;
	struct {
		cl_float R, RR, G;
	} response;
	int m;

	cl_uint attr;

	// start position and size of elements in global
	// opencl memory of connections (edge) for same layer
	struct {
		int offset, size;
	} edges_loc[OCL_EDGE_COUNT];
} ocl_layer1_data;

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

// struction to hold pointer to parts for each image position (x,y) (8 bytes)
typedef struct { 
	int offset, size;
} ocl_layer1_data_coordinates;

// pointer for sorting of coordinates
typedef struct {
	float sort_value; // usually R_RESPONSE value
	int coord_ind;	// usually original coordinate index value (when used in crate_layer0.cl) or pointer to where candidates for this position are curently saved (when used in add_layer1.cl)
} ocl_coord_sorting_pointer;

/*
 * Structure for lib_part must hold two memory arrays for each layer:
 * 
 * 1. One for part_data i.e. all parts in library for specific layer:
 *		- Each part_data MUST be in same position as its type number.
 *		- Each part_data may have one or more (or none) connections to other parts_data in different 
 *		  layer (i.e. in different memory object). Each connection can have different type.
 *		- For each connection type in one part_data there is reported position (offset) where
 *		  data connections for this type start and count (size) of connections. Values are
 *		  referred to memory object for same layer
 *		
 *	 [p1, p1, p3 ....]
 *
 * 2. One for part_data2 i.e. all connections associated with parts in specific layer:
 *		- Values in memory object for specific layer must be sorted by part and by connection type
 *		- For each connection we also hold information of pointed part (i.e. part_data that this
 *		  connection is pointing to). We hold position (offset) of this part in memory object and
 *		  layer number for which this part belongs to (i.e. offset can referre to memory object of
 *		  different layer)
 *		- CAUTION: All values in ocl_part_data_2 can also be empty if there is no information for 
 *		  connection type - BUT node.offset and node.layer MUST always be set
 *
 *	 [e1(p1), e2(p1), .... | e3(p2), e4(p2) ... | .... ]
 */

// opencl part_data (72 bytes)
typedef struct {
	int cmx, cmy;
	int layer;
	int type;
	struct {
		float R, RR, G;
	} thresh;

	int bpcount; // number of "basic" parts
	cl_uint attr;
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

// max size of distribution (WARNING: when adjusting this constant do not forget same one in utils.cl)
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
	cl_uint lib_part_offset;
	cl_uint img_part_offset;
	cl_uint skip;
	float r_sum;
	float realization_ratio;	
	float g_prod;
	int x,y;
	int new_x,new_y;
	//int result_ids_size;
	int result_edges_prev_lay_size;
	int result_edges_layer0_size;
	int part_type;
	int img_part_type;
	cl_uint max_schur_mem_offset;
} ocl_layer_candidate;

inline void ocl_new_attr(cl_uint* attr, cl_uint a) { *attr = a; }
inline void ocl_set_attr(cl_uint* attr, cl_uint a) { *attr |= a; }
inline void ocl_clear_attr(cl_uint* attr, cl_uint a) { *attr &= ~a; }
inline bool ocl_is_attr_set(cl_uint attr, cl_uint a) { return (attr & a) == a; }
inline bool ocl_is_attr_set2(cl_uint attr, cl_uint a) { return (attr & a) != 0; }

#define SCHUR_PROD_USE_VAL_RESPONSE 0
#define SCHUR_PROD_USE_IDENTITY_RESPONSE 1
#define SCHUR_PROD_USE_R_RESPONSE 2
#define SCHUR_PROD_USE_G_DISTR_RESPONSE 3
#define SCHUR_PROD_USE_SIMPLE_G_DISTR_RESPONSE 4

#endif /* _PART_LIB_ */
