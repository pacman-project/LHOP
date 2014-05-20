
#include "core/legacy/part_lib.h"
#include "core/legacy/layer_1_result.h"

#include "core/legacy/inference/layer_n_creators.h"

#include "core/legacy/learning/layer_learning.h"
#include "core/legacy/learning/optimization.h"

#include "utils/graphs/graph.h"
#include "utils/serialization/streaming.h"
#include "utils/mapreduce.h"

#include "init.h"

void register_classes_core() {
	/***************************************************************************
	 * !!!! DO NOT MODIFY THIS FUNCTION UNLESS YOU'VE CHANGED STATIC_CORE  !!!!!
	 ***************************************************************************/
	
	// Graphs
	/////////////////////////////////////////////////////////////
	streamable::register_class<graph,StreamableId::GRAPH>();
	
	streamable::register_class<node,StreamableId::NODE>();
	
	streamable::register_class<node_data_t<int>,StreamableId::NODE_DATA_T_INT>();
	streamable::register_class<node_data_t<char>,StreamableId::NODE_DATA_T_CHAR>();
	streamable::register_class<node_data_t<double>,StreamableId::NODE_DATA_T_DOUBLE>();
	
	streamable::register_class<edge_data,StreamableId::EDGE_DATA>();
	streamable::register_class<edge_data_t<double>,StreamableId::EDGE_DATA_T_DOUBLE>();

	streamable::register_class<img_graph,StreamableId::IMG_GRAPH>();
	streamable::register_class<img_node_data,StreamableId::IMG_NODE_DATA>();
	
	// LHOP Library 
	/////////////////////////////////////////////////////////////
	streamable::register_class<part_data,StreamableId::PART_DATA>();
	streamable::register_class<part_lib,StreamableId::PART_LIB>();

	streamable::register_class<part_data_2,StreamableId::PART_DATA_2>();
	streamable::register_class<part_data_2a,StreamableId::PART_DATA_2A>();
	streamable::register_class<part_data_2c,StreamableId::PART_DATA_2c>();
	streamable::register_class<cpart_data,StreamableId::CPART_DATA>();
	streamable::register_class<vs_part_data,StreamableId::VS_PART_DATA>();
	streamable::register_class<part_data_sim,StreamableId::PART_DATA_SIM>();

	// LHOP layer files
	/////////////////////////////////////////////////////////////	
	streamable::register_class<layer1_result,StreamableId::LAYER1_RESULT>();
	streamable::register_class<layer1_data,StreamableId::LAYER1_DATA>();
	
	streamable::register_class<edge_data_name,StreamableId::EDGE_DATA_NAME>();
	
	// LHOP Inference
	/////////////////////////////////////////////////////////////
	streamable::register_class<layern_creator,StreamableId::LAYERN_CREATOR>();
	
	// LHOP Learning
	/////////////////////////////////////////////////////////////	
	streamable::register_class<map_learning,StreamableId::MAP_LEARNING>();
	streamable::register_class<part_learning,StreamableId::PART_LEARNING>();
	streamable::register_class<online_geo_learning,StreamableId::ONLINE_GEO_LEARNING>();
	streamable::register_class<online_path_map,StreamableId::ONLINE_PATH_MAP>();
	streamable::register_class<pca_learning,StreamableId::PCA_LEARNING>();

	// MISC
	/////////////////////////////////////////////////////////////	
	streamable::register_class<circular_tolerance<int>,StreamableId::CIRCULAR_TOLERANCE_INT>();
	streamable::register_class<K_bin,StreamableId::K_BIN>();
	streamable::register_class<streamed_pointer,StreamableId::STREAMED_POINTER>();
	
	// Map Reduce
	/////////////////////////////////////////////////////////////	
	streamable::register_class<mapreduce_result,StreamableId::MAPREDUCE_RESULT>();
	streamable::register_class<mapreduce_items,StreamableId::MAPREDUCE_ITEMS>();

	streamable::register_class<map_learning_mapreduce,StreamableId::MAP_LEARNING_MAPREDUCE>();
	streamable::register_class<part_learning_mapreduce,StreamableId::PART_LEARNING_MAPREDUCE>();
	streamable::register_class<layern_creator_mapreduce,StreamableId::LAYERN_CREATOR_MAPREDUCE>();
}