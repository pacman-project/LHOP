// Constant definitions used accross all modules
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

/*********************************************
 * = ADD NEW GLOBAL CONSTANST HERE
 *
 *
 * !!! MAKE SURE THERE ARE NOT ID CLASHES !!!
 *********************************************/

// Graph edge connection types
//    - Define new connections here by incrementing last one
//    - Make sure IDs are unique and consistent across all possible modules (coordinate with other developers)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum EdgeConnection {
	// part_lib connections
	TO_NEIGHBOOR = 0,

	//TO_REGION = 1,

	// Edges to all subparts (one layer down) other then center (in library; class part_lib).
	TO_LYR_SOURCE = 2, // connections to all parts in lower layer (one layer down) for each of the subparts other then center

	//TO_LYR_BACK = 3,

	// Edge to center (in library; class part_lib).
	TO_LYR_CENTER = 4, // inverse of lyrCenterBack (never used - maybe ?? :D)

	// Edges to all nodes (i.e. parts) in the next layer (one level up) which have same 
	// part for its center (in library; class part_lib).
	TO_LYR_CENTER_BACK = 5, // connection to all parts in upper layer (one level one) which have this part for its center

	// Edges to all subparts fro, the previous layer (one layer down) that were used in the 
	// inference process of this part (in class layer1_result)
	TO_PREV_LAYER = 6, // connection to all parts in lower layer (one layer down) that are used in this pard (in layer1_result for class layer1_data)
	TO_LYR_PREV = 7, // (in library; class part_lib).

	// Edges to forbidden subparts -- parts which are not allowed to be present at specific
	// positions -- (one layer down, in library; class part_lib). 
	TO_LYR_FORBIDDEN = 8,

	//TO_DIRECTION0 = 9,
	//TO_DIRECTION1 = 10,
	//TO_DIRECTION2 = 11,
	//TO_DIRECTION3 = 12,
	//TO_DIRECTION4 = 13, 
	//TO_DIRECTION5 = 14,
	//TO_DIRECTION6 = 15,
	//TO_DIRECTION7 = 16,
	//TO_DIRECTION8 = 17,
	//TO_DIRECTION9 = 18,
	//TO_DIRECTION10 = 19,
	//TO_DIRECTION11 = 20,
	//TO_DIRECTION12 = 21,
	//TO_DIRECTION13 = 22,
	//TO_DIRECTION14 = 23,
	//TO_DIRECTION15 = 24,

	// layer1_result connections
	TO_SIMILAR = 25,
	TO_NEXT_LAYER = 26, // used only in optimize_layer_2
	//TO_PREV_LAYER_I = 27,
	//TO_SEQUENCE_NEIGHBOR = 28,
	//TO_LYR_CHILD = 29,
	//TO_LYR_PARENT = 30,
	//TO_PROJECTION = 31,

	TO_NEXT_LAYER0 = 32, // used only in add_reconstruction_edges/add_reconstruction_edges_link
	//TO_NEXT_LAYER1 = 33,
	//TO_NEXT_LAYER2 = 34,
	//TO_NEXT_LAYER3 = 35,
	//TO_NEXT_LAYER4 = 36,
	//TO_NEXT_LAYER5 = 37,
	//TO_NEXT_LAYER6 = 38,

	TO_LAYER0 = 39, // edges to layer 0 also refered to as reconstruction edges
	TO_LYR_SRC_M = 40, // used as get only in get_similar_types/save_visview_centers but never set (maybe it can be removed)
	TO_HYPO_NODE = 41, // not used
	TO_LYR_CENTER1 = 42, // not uesed
	TO_LYR_SIMILAR = 43,

	// Edges to root parts ("similarity root node")
	TO_LYR_SIMROOT = 44,

	// Edges to variable-shape root parts
	TO_LYR_VS_ROOT = 45, // used only in add_laye1 but never set
	TO_LYR_VS_MEMBER = 46, // not used since VS_PART_ATTR is never set

	//TO_LYR_SIMVAL = 47,
	//TO_MISSING = 47,
	TO_MINMAX_LAYER0 = 48 // used only in save_visview (maybe it can be removed)
};

// Streamable types
//    - Define IDs of new types here by incrementing last one
//    - Make sure IDs are unique and consistent across all possible modules (coordinate with other developers)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
 * Streamable IDs used to identify specific class.
 * Each class must register its class name with this number (this must be done where the class is implemented i.e. in module)
 * streamable::register_class<StreambaleImplClass,StreamableId::<new_class_id>);
 */
enum StreamableId {
	NODE_DATA_T_INT = 0,
	NODE_DATA_T_CHAR = 1,
	NODE_DATA_T_DOUBLE = 2,
	EDGE_DATA = 3,
	PART_DATA = 4,
	PART_LIB = 5,
	IMG_NODE_DATA = 6,
	LAYER1_DATA = 7,
	NODE = 8,
	GRAPH = 9,
	IMG_GRAPH = 10,
	LAYER1_RESULT = 11,
	LAYER1_RESULT_STRUCT = 12,
	LAYER1_RESULT_APP = 13,
	LAYER1_RESULT_DOG = 14,
	PART_DATA_2 = 15,
	PART_DATA_2c = 16,
	//COMPRESSED_TREE = 17,
	//LIBRARY_TREE = 18,
	EDGE_DATA_T_DOUBLE = 19,
	//EQUIVALENCE_CLASS = 20,
	//CLASS_DISTRIBUTION = 21,
	//OBJECT_DISTRIBUTION = 22,
	//CLUSTER = 23,
	CIRCULAR_TOLERANCE_INT = 24,
	//CLUSTER_EDGE_DATA = 25,
	//CLUSTER_EDGE_DATA_R = 26,
	EDGE_DATA_NAME = 27,
	//RPART_DATA = 28,
	//STREAMABLE_T = 29,
	//PART_DATA_2R = 30,
	CPART_DATA = 31,
	PART_DATA_2A = 32,
	//OPART_DATA = 33,
	LAYER1_RESULT_LOGGABOR = 34,
	//SPART_DATA = 35,
	VS_PART_DATA = 36,
	PART_DATA_SIM = 37,
	LAYER1_RESULT_COLOR = 38,
	LAYERN_CREATOR = 39,
	MAPREDUCE_RESULT = 40,
	MAPREDUCE_ITEMS = 41,
	MAP_LEARNING_MAPREDUCE = 42,
	PART_LEARNING_MAPREDUCE = 43,
	LAYERN_CREATOR_MAPREDUCE = 44,
	K_BIN = 45,
	MAP_LEARNING = 46,
	PART_LEARNING = 47,
	ONLINE_GEO_LEARNING = 48,
	ONLINE_PATH_MAP = 49,
	PCA_LEARNING = 50,
	STREAMED_POINTER = 51
};

#endif /* _CONSTANTS_H_ */
