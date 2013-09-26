/// Initialization of atoms and serialization


#include "layer_1_result.h"
#include "part_lib.h"
#include "layer_n_creators.h"
#include "layer_learning.h"
#include "optimization.h"
#include "../utils/mapreduce.h"

// initialize custom class allocation

//DEFINE_CLASS_ALLOCATOR_2(node, 2048);
//DEFINE_CLASS_ALLOCATOR_2(edge_data_name, 2048);
//DEFINE_CLASS_ALLOCATOR_2(layer1_data, 2048);
//DEFINE_CLASS_ALLOCATOR(equivalence_class);


/// Initialize atoms (to have predetermined indices.
void init_atoms()
{
    atom("toNeighbor");
    atom("toRegion");	
    atom("lyrSrc"); // connections to all parts in lower layer (one layer down) for each of the subparts other then center
    atom("lyrBack"); 
    atom("lyrCenter"); // inverse of lyrCenterBack (never used - maybe ?? :D)
    atom("lyrCenterBack"); // connection to all parts in upper layer (one level one) which have this part for its center
    atom("toPrevLayer"); // connection to all parts in lower layer (one layer down) that are used in this pard (in layer1_result for class layer1_data)
    atom("lyrPrev");
    atom("lyrForbidden");
    atom("toDirection0");
    atom("toDirection1");
    atom("toDirection2");
    atom("toDirection3");
    atom("toDirection4");
    atom("toDirection5");
    atom("toDirection6");
    atom("toDirection7");
    atom("toDirection8");
    atom("toDirection9");
    atom("toDirection10");
    atom("toDirection11");
    atom("toDirection12");
    atom("toDirection13");
    atom("toDirection14");
    atom("toDirection15");
	// connections to all parts in same layer that are similar (in library for class part_data)
    atom("toSimilar");
    atom("toNextLayer");
    atom("toPrevLayerI");
    atom("toSequenceNeighbor");
    atom("layToChild");
    atom("layToParent");
    atom("toProjection");
    atom("toNextLayer0");
    atom("toNextLayer1");
    atom("toNextLayer2");
    atom("toNextLayer3");
    atom("toNextLayer4");
    atom("toNextLayer5");
    atom("toNextLayer6");
    atom("toLayer0");
    atom("lyrSrcM");
    atom("toHypoNode");
    atom("lyrCenter1");
    atom("lyrSimilar");
    atom("lyrSimRoot");
    atom("lyrVSRoot");
    atom("lyrVSMember");
    atom("lyrSimVal");
    atom("toMissing");
}


/// initialize classes for streaming.
void init_streaming()
{
    streamable::get_class_id_store().register_class(new node_data_t<int>());
    streamable::get_class_id_store().register_class(new node_data_t<char>());
    streamable::get_class_id_store().register_class(new node_data_t<double>());
    streamable::get_class_id_store().register_class(new edge_data());
    streamable::get_class_id_store().register_class(new part_data());
    streamable::get_class_id_store().register_class(new part_lib(0));
    streamable::get_class_id_store().register_class(new img_node_data());
    streamable::get_class_id_store().register_class(new layer1_data());
    streamable::get_class_id_store().register_class(node::new_node());
    streamable::get_class_id_store().register_class(new graph());
    streamable::get_class_id_store().register_class(new img_graph());
    streamable::get_class_id_store().register_class(new layer1_result());
    streamable::get_class_id_store().register_class(new layer1_result_struct());
    streamable::get_class_id_store().register_class(new layer1_result_app());
    streamable::get_class_id_store().register_class(new layer1_result_dog());
    streamable::get_class_id_store().register_class(new part_data_2());
    streamable::get_class_id_store().register_class(new part_data_2c());
    streamable::get_class_id_store().register_class(nullptr);    // new compressed_tree()
    streamable::get_class_id_store().register_class(nullptr);    // new library_tree()
    streamable::get_class_id_store().register_class(new edge_data_t<double>());
    streamable::get_class_id_store().register_class(nullptr);    // new equivalence_class()
    streamable::get_class_id_store().register_class(nullptr);    // new class_distribution()
    streamable::get_class_id_store().register_class(nullptr);    // new object_distribution()
    streamable::get_class_id_store().register_class(nullptr);    // new cluster()
    streamable::get_class_id_store().register_class(new circular_tolerance<int>());
    streamable::get_class_id_store().register_class(nullptr);    // new cluster_edge_data()
    streamable::get_class_id_store().register_class(nullptr);    // new cluster_edge_data_r()
    streamable::get_class_id_store().register_class(new edge_data_name());
    streamable::get_class_id_store().register_class(new rpart_data());
    streamable::get_class_id_store().register_class(new streamable_t<int>());
    streamable::get_class_id_store().register_class(new part_data_2r());
    streamable::get_class_id_store().register_class(new cpart_data());
    streamable::get_class_id_store().register_class(new part_data_2a());
    streamable::get_class_id_store().register_class(new opart_data());
	streamable::get_class_id_store().register_class(new layer1_result_loggabor());
    streamable::get_class_id_store().register_class(new spart_data());
    streamable::get_class_id_store().register_class(new vs_part_data());
    streamable::get_class_id_store().register_class(new part_data_sim());
    streamable::get_class_id_store().register_class(new layer1_result_color());
	streamable::get_class_id_store().register_class(new layern_creator());
    streamable::get_class_id_store().register_class(new mapreduce_result());
    streamable::get_class_id_store().register_class(new mapreduce_items());
    streamable::get_class_id_store().register_class(new map_learning_mapreduce());
    streamable::get_class_id_store().register_class(new part_learning_mapreduce());
	streamable::get_class_id_store().register_class(new layern_creator_mapreduce());
    streamable::get_class_id_store().register_class(new K_bin());
    streamable::get_class_id_store().register_class(new map_learning());
    streamable::get_class_id_store().register_class(new part_learning());
    streamable::get_class_id_store().register_class(new online_geo_learning());
    streamable::get_class_id_store().register_class(new online_path_map());
    streamable::get_class_id_store().register_class(new pca_learning());
	streamable::get_class_id_store().register_class(new streamed_pointer());
}

// attributes for nodes
//const unsigned NODE_REDUCED_ATTR = 1   
//const unsigned IMG_NODE_ATTR = 2
//const unsigned HAS_NEXT_LAYER = 4
//const unsigned GRID_NODE_ATTR = 8
//const unsigned ISLAND_ATTR = 16
//const unsigned EARTH_CENTER_ATTR = 32
//const unsigned FROM_PREV_LAYER = 64
//const unsigned TEXTURE_NODE = 128 
//const unsigned CLUSTER_PART_ATTR = 256;
//const unsigned ATTR_VISITED = 1024   -----+
//const unsigned ATTR_CHECKED = 2048        | =
//const unsigned ATTR_MATCHED = 1024   -----+
//const unsigned IN_NODE_ATTR = 0x1000
//const unsigned BOUDARY_NODE_ATTR = 0x2000
//const unsigned OUT_NODE_ATTR = 0x4000
//const unsigned NODE_DELETED_ATTR = 0x8000
//const unsigned HYPO_NODE_ATTR = 0x10000
//const unsigned CATEGORY_NODE_ATTR = 0x20000
//const unsigned ALLOW_UPDATE_ATTR = 0x40000
//const unsigned HYPO_NODE_CANDIDATE_ATTR = 0x80000
//const unsigned R_PART_ATTR = 0x100000
//const unsigned OBJ_PART_ATTR = 0x200000
//const unsigned HIDDEN_NODE_ATTR = 0x400000
//const unsigned ATTR_MARKED = 0x800000
//const unsigned ATTR_MARKED1 = 0x1000000
