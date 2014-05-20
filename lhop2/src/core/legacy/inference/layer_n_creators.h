
// layer_n_creators
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _LAYER_N_CREATORS_
#define _LAYER_N_CREATORS_

class layern_creator;

#include <ctime>
#include "utils/structures.h"
#include "utils/ocv.h"
#include "utils/utils.h"
#include "utils/config_dictionary.h"

#include "core/legacy/constants.h"

#include "core/legacy/layer_1_result.h"

// constants
///////////////////////////////////////////////////////////////////////////////

const unsigned ALLOW_UPDATE_ATTR = 0x40000;

const unsigned TEXTURE_NODE = 128;

const int invalid_edge_name = INT_MAX;

// layern_creator
///////////////////

class layern_creator : public streamable {
protected:
    struct new_edge_data_t {
        node* n;      // destination node
        int name;     // edge name
        int spindex;  // name of the subpart
        double r;     // response of the subpart

        new_edge_data_t(node* pn, int pname, int pindex, double pr) 
            : n(pn), name(pname), spindex(pindex), r(pr) { }
    };

    typedef vector<new_edge_data_t> new_edges_data_t;

    typedef dnpair (layer1_result::*product_s_mt)(int, int, int, int, int, matrix<double>&, int, int, double,
        schur_product_function*, double, double, bin_function<int,int,double>&, layer1_result::sp_result_t&);
    typedef dnpair (layer1_result::*product_ns_mt)(int, int, int, int, int, matrix<double>&, int, int, double,
        schur_product_function*, double, layer1_result::sp_result_t&);
    typedef dnpair (layer1_result::*product_ts_mt)(int, int, int, int, int, matrix<double>&, const map<int, double>&, int, double,
        schur_product_function*, double, layer1_result::sp_result_t&);

    part_lib* library;
	bool is_library_my;
	
public:

    // parameters
    double layer_contraction;               // next layer contraction factor
    bool manual_thresholds;
    double threshold_factor;
    double candidate_r_threshold;       
    double candidate_r_threshold_percent;
    double candidate_g_threshold;
    double candidate_g_threshold_percent;
    double r_response_threshold;            // add to next layer threshold (R_RESPONSE)
    double g_response_threshold;            // add to next layer threshold (G_RESPONSE)
    double s_response_threshold;            // threshold for S_RESPONSE
    double x_response_threshold;            // threshold for X_RESPONSE
    double g_response_threshold_percent;    // add to next layer threshold (G_RESPONSE), percent of the best
    bool identity_g_response;               // = false
    bool simple_g_response;                 // = false
    bool ignore_g_distribution;             // = false
    int g_response_operation;               // 0 = addition, 1 = multiplication; default = 1
    double g_response_var_factor;           // = 0.0 (no multiplication of g-response variance!)
    double type_thresh;                     // = 1.0 (legacy)
    double convolution_threshold;
    double convolution_link_threshold;
    double proj_max_forb_threshold;
    double forb_quot_threshold;
    double r_response_pow;
    double g_response_pow;
    double realization_ratio_threshold;     // = 0.5
    bool shape_check;                       // = false
    bool link_missing_support;              // = false
    double continuity_factor;               // = -1 (no continuity check if ~ <= 0)
    bool normalize_histogram;               // = true - experimental - to be removed!
    bool copy_prev_layer;                   // = false
    int texture_parts;                      // = 4
    int texture_radius;                     // = 4
    bool ignore_texture;                    // = false 
    double min_factor;                      // = 1.0
    double max_factor;                      // = 1.0
    int reconstruction_type;                // = 1    (1 = normal; 2 = 
    bool add_reconstruction_edges;          // = false
    bool add_activation_edges;              // = true
    double hypo_unrealized_threshold;       // = 0.4 (unrealized/all <= ~ --> add)
    double hypo_ratio_threshold;            // = 1E+6 (#of hypo nodes/#reconstruction nodes <= ~ --> add)
    int null_tolerance_limit;               // 3
    int positive_tolerance_threshold;       // 2
    int projection_radius;                  // 3
    double reconstruction_factor;           // 0.5
    int rec_null_tolerance_limit;           // 4
    int rec_tolerance_threshold;            // 1
    int tolerance_radius;                   // 4
    int min_part_distance;                  // -1
    bool depth_first_search;                // true (perform dfs on RP
    double hypo_start;                      // 0.5
    int hypo_nbhood;                        // 3
    double hypo_factor;                     // 1.2
    double hypo_val_threshold;              // 0.4
    int hypo_voters_threshold;              // 2
    bool add_edge_names;                    // false
    int new_positions_threshold;            // INT_MAX
    int same_position_threshold;            // INT_MAX
    set<int> allowed_parts;                 // empty

    int variation_dimension;                // how many "directions" we make
    int variation_factor;                   // 

    int inhibit_layer_response;             // default value is FALSE (no inhibition); G_RESPOSE, R_RESPONSE, RR_RESPONSE are possible values
    int inhibit_layer_max;                  // default value is INT_MAX
    double inhibit_layer_thresh;            // thresh for intersection of "boxes", default value is 0.4
    bool inhibit_layer_delete;              // delete inhibited nodes or not, default value is true

    bool strict_svm_checking;               // if svm file does not exist for certain object, discard it, default is false
    int scale_merge;                        // how many scales are merged, default value = 1 (no merge)

public:
	bool use_opencl;

    // info
    double min_cover_quot;                  // minimal cover quotient
    double quot_sum;
    int count;             

	// private constructor for serialization
	layern_creator() : library(nullptr) {}

    layern_creator(const char* cfgfile);
	layern_creator(const ConfigDictionary& cfg);
    ~layern_creator();

    part_lib* get_library();
    void set_library(part_lib* lib, bool is_my_lib = false);
    void set_library(part_lib* lib, const ConfigDictionary& cfg, bool is_my_lib = false);

    void add_layer(layer1_result* res, int layer, int stoplayer);
    void add_layer(vector<layer1_result*>& res, int layer, int stoplayer);
    void add_layer(layer1_result* res, const scmap_t& scmap, int layer, int stoplayer);
    void add_layer(vector<layer1_result*>& res, vector<scmap_t>& scmap, int layer, int stoplayer);
    bool category_layer_creator() { return reconstruction_type == 3; }

protected:

    void init();
	void cfg_init(const ConfigDictionary& cfg);
    double get_thresh(lib_data* d, int name, double defval);

    void add_layer1(layer1_result* res, int layer);    
	void add_layer1(layer1_result* res, const scmap_t& scmap, int layer);

    void add_layer3(layer1_result* res, int layer);

	bool add_layer7(layer1_result* res, int layer, const irectangle2& region);  
    bool add_layer7(layer1_result* res, const scmap_t& scmap, int layer, const irectangle2& region);  

	template<class I> itriple get_new_edges(vector<new_edge_data_t>& new_edges, int k, I begin, I end);
	void add_edges(node* newn, vector<new_edge_data_t>& new_edges /*, unsigned attr = 0*/);
    void add_edges(layer1_result* res, node* newn, vector<new_edge_data_t>& new_edges /*, unsigned attr = 0 */);
    void set_candidate_thresholds(layer1_result* res, int k);

    void mark_texture(layer1_result* res, int k);
    void inhibit_result(layer1_result* res, int k);

	void copy_layer(layer1_result* res, int k);
    void check_geometry(double& benergy, vector<double>& dvector, double& scdistance, 
        const path_map_t& lmap, path_map_t& imap, const ipoint2& center, bool calcbenergy);
    void get_new_support_edges(new_edges_data_t& edges, const vector<node*>& nodes);

	virtual streamable* make_instance() const { return new layern_creator(); }
	virtual void read_from_stream(istreamer& is);
	virtual void write_to_stream(ostreamer& os);
};



// I: iterator over a container of layer1_result::sp_result_pair_t = pair<node*, ipoint2>
// Return value: (#toPrev, #toLy0, #toHypo)
template<class I> itriple layern_creator::get_new_edges(vector<new_edge_data_t>& new_edges, int k, I begin, I end)
{
    typedef triple<node*, int, ipoint2> triple_t;

    /*static*/ int to_layer0 = EdgeConnection::TO_LAYER0;
    /*static*/ int to_prev_layer = EdgeConnection::TO_PREV_LAYER;
    /*static*/ int to_hyponode = EdgeConnection::TO_HYPO_NODE;

    itriple result(0, 0, 0);

    if (!add_activation_edges) 
        return result;

    if (!add_reconstruction_edges) {
        for (; begin != end; ++begin) new_edges.push_back(new_edge_data_t(begin->n, to_prev_layer, begin->index, begin->r));
        result.first = (int)new_edges.size();
    } else {
        for (I iter = begin; iter != end; ++iter) {
            node* n = iter->n;
            int p = iter->index;
            double r = iter->r;

            new_edges.push_back(new_edge_data_t(n, to_prev_layer, p, r));
            ++result.first;
            if (n->is_attr_set(HYPO_NODE_ATTR)) {
                new_edges.push_back(new_edge_data_t(n, to_hyponode, invalid_edge_name, 1.0));
                ++result.third;
            }
            if (k == 0) {
                new_edges.push_back(new_edge_data_t(n, to_layer0, invalid_edge_name, 1.0));
                ++result.second;
            }
        }
        if (k > 0) {
            set<node*> nbset;

            for (I iter = begin; iter != end; ++iter) {
                node* n = iter->n;
                foreach_neighbor(n, to_layer0, nbiter) nbset.insert(neighbor_node(nbiter)); 
            }
            for (set<node*>::iterator iter = nbset.begin(); iter != nbset.end(); ++iter) 
                new_edges.push_back(new_edge_data_t(*iter, to_layer0, invalid_edge_name, 1.0));
            result.second += (int)nbset.size();
            nbset.clear();
            for (I iter = begin; iter != end; ++iter) {
                node* n = iter->n;
                foreach_neighbor(n, to_hyponode, nbiter) nbset.insert(neighbor_node(nbiter)); 
            }
            for (set<node*>::iterator iter = nbset.begin(); iter != nbset.end(); ++iter) 
                new_edges.push_back(new_edge_data_t(*iter, to_hyponode, invalid_edge_name, 1.0));
            result.third += (int)nbset.size();
        }
    }
    return result;
}

#endif /* _LAYER_N_CREATORS_ */
