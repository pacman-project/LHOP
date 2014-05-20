// layer_learning learning algorithms:
//   * s_learning     learning of layers (statistics, g-responses, em-step, optimization)
//   * obj_learning   learning of object layer (models)
//   * t_learning     learning of thresholds

#pragma once

#ifndef _OBJECT_LEARNING_
#define _OBJECT_LEARNING_

#include "core/legacy/layer_1_result.h"
#include "core/legacy/inference/layer_n_creators.h"

#include "core/legacy/learning/layer_learning.h"

#include "utils/serialization/streamed_pointer.h"

#include <tuple>
#include <stdio.h>

// map_learning
/////////////////

// obj_learning (action=learn_objects)
///////////////////////////////////////////////////////////////////////////////

class obj_learning {
public:
    // typedefs
    //typedef pair<vector<ipoint2>, ipoint2> cluster_data_t; // ({(lyr. dif., type),...}, pos.)

    typedef vector<part_lib::cluster_data_t> obj_data_t;
    typedef pair<streamed_pointer, irectangle2> validation_data_t;
    typedef list<pair<double, obj_data_t> > objects_t;

    // parameters
    part_lib* library;
    part_lib* libraryD;
    int layer;
    double contraction;
    vector_parameter<double> r_response_threshold;
    vector_parameter<double> g_response_threshold;
    vector_parameter<double> rr_response_threshold;
    vector_parameter<int> max_objects;
    vector_parameter<int> max_add;
    vector_parameter<int> max_cluster_n;
    vector_parameter<int> min_cluster_n;
    vector_parameter<int> cluster_size;
    vector_parameter<double> cluster_member_threshold;
    vector_parameter<double> cover_threshold;
    vector_parameter<double> intersection_threshold;
    double cover_threshold0;
    int max_depth;
    double validation_threshold; // accept part if true/false >= validation_threshold
    int gaussian_dim;
    double gaussian_sigma;
    int reduce_radius;

    // validation
    list<validation_data_t> vset;
    layern_creator* creator;

    // results
    list<obj_data_t> objects;

    // methods
    obj_learning(const ConfigDictionary& cfg);
    ~obj_learning();

    void object_from_result(layer1_result* res, const scmap_t& scmap, const irectangle2& gtr);
    int add_to_library(part_lib* library, const list<obj_data_t>& objects, const string& name);
    int add_to_library(part_lib* library, const obj_data_t& object, const string& name);
    int add_to_library(const string& name);
	void add_validation_data(const streamed_pointer& ptr, const irectangle2& gtruth);

	void set_library(part_lib* plibrary);
protected:	
    void cfg_init(const ConfigDictionary& cfg);
    void object_from_result(objects_t& result, layer1_result* res, const scmap_t& scmap,
        int layer, const set<ipoint2>& covered, const irectangle2& gtr);

    bool validate_object(layer1_result* res, const irectangle2& gtr, const obj_data_t& object);
    bool validation_function_1(int tr, int fa);
    
	void reduce_object_data(obj_data_t& od);
};

// o_learning (action=learn_objects2)
///////////////////////////////////////////////////////////////////////////////

struct olv_object_part {
    typedef vector<node*> structure_t;  

    structure_t str;
    set<ipoint2> rec;

    olv_object_part() : str(), rec() { }

    void augment(node* n, const set<ipoint2>& s) 
    { 
        str.push_back(n); 
        rec.insert(s.begin(), s.end());
    }

    int type(int i) const { return ((layer1_data*)str[i]->data)->m; }
    ipoint2 pos(int i) const 
    { 
        layer1_data* nd = (layer1_data*)str[i]->data;
        return ipoint2(nd->x, nd->y);
    }
    int size() const { return (int)str.size(); }
    void print(ostream& os) const;
};

class o_learning {
protected:
    typedef map<iipair, ip2_vector> dupmap_t;
    typedef vector<pair<double, int> > mstat_t;
    typedef pair<streamed_pointer, list<irectangle2> > vset_item_t;
    typedef list<vset_item_t> vset_t;
    typedef vector<pair<double, int> > statmap_t;

    // parameters
    string catname;
    int max_models;
    int max_cluster_n;
    int min_cluster_n;
    int cluster_size;
    double cluster_member_threshold;
    double intersection_threshold;
    double hit_ratio_threshold;
    double hit_threshold;
    int redundancy_threshold;
    int type_bite_threshold;
    int gaussian_dim;
    double gaussian_sigma;
	ConfigDictionary inference_cfg;

    // algorithm fields, statistics, ...
    part_lib* library, * srclibrary;

    int srclayer;
    double contraction;
    
    vset_t vset; // validation set

    matrix<double> dist;

    dupmap_t dupmap;
    dupmap_t dupstat;
    statmap_t statmap;
    mstat_t mstat;
    vectorn<int>::type2 rhits;

public:
    o_learning(const ConfigDictionary& cfg);
    ~o_learning();
    
    void reset();

    // Adds to validation set (makes a copy and saves to disk)
    void add_to_validation_set(layer1_result* res, const list<irectangle2>& gtrs);
    void learn_models(layer1_result* res, const scmap_t& scmap);
    void update_duplet_statistics(layer1_result* res);
    part_lib* get_library() { return library; }
    void finalize();

protected:
    void add_to_library(part_lib* plb, list<olv_object_part>& models, layer1_result* res, const scmap_t& scmap);
    void augment_model(list<olv_object_part>& models, int& modelcount, olv_object_part model, const mstat_t& mstat,
        const map<node*, int>& indexmap, layer1_result* res, int maxmodels);
    void candidate_positions(list<pair<node*, set<ipoint2> > >& v, olv_object_part& model, 
        const map<node*, int>& indexmap, layer1_result* res, int m);
    bool check_duplet(olv_object_part& model, layer1_result* res, int m, const ipoint2& p);
    void get_mstat(mstat_t& mstat);
    void init_cfg(const ConfigDictionary& cfg);
    void make_duplets();
    void make_models(list<olv_object_part>& models, layer1_result* res);    

    void update_statmap(layer1_result* res, const list<irectangle2>& gtrs);
    void validate(vector<int>& tokeep, vectorn<int>::type2& rhits, part_lib* library, 
        int libsize0, int libsize1);
};


#endif /* _OBJECT_LEARNING_ */

