// layer_learning learning algorithms:
//   * s_learning     learning of layers (statistics, g-responses, em-step, optimization)
//   * obj_learning   learning of object layer (models)
//   * t_learning     learning of thresholds

#pragma once

#ifndef _LAYER_LEARNING_
#define _LAYER_LEARNING_

#include <map>
#include <set>

#include "core/legacy/layer_1_result.h"

#include "utils/structures.h"


////////////////////////
// em-step related only


// 2d-distribution
struct distribution2 {
    typedef map<ipoint2, double> map_t;

    map_t m;

    distribution2() : m() { }

    void update(const ipoint2& p, double d) 
    {   
        pair<map_t::iterator, bool> ibpair = m.insert(map_t::value_type(p, d));
        if (!ibpair.second) ibpair.first->second += d;
    }

    void update(const distribution2& d) 
    {
        for (map_t::const_iterator iter = d.m.begin(); iter != d.m.end(); ++iter) {
            update(iter->first, iter->second);
        }
    }

    pair<ipoint2, double> get_peak();
};

// histogram
struct histogram {
    typedef map<int, double> map_t;

    map_t m;

    histogram() : m() { }

    void update(int i, double d) 
    {
        pair<map_t::iterator, bool> ibpair = m.insert(map_t::value_type(i, d));
        if (!ibpair.second) ibpair.first->second += d;
    }

    void update(const histogram& h) 
    {
        for (map_t::const_iterator iter = h.m.begin(); iter != h.m.end(); ++iter) {
            pair<map_t::iterator, bool> ibpair = m.insert(map_t::value_type(iter->first, iter->second));
            if (!ibpair.second) ibpair.first->second += iter->second;
        }
    }

    friend ostream& operator<<(ostream& os, const histogram& h) 
    { 
        for (histogram::map_t::const_iterator iter = h.m.begin(); iter != h.m.end(); ++iter) 
            os << '[' << iter->first << ',' << iter->second << ']';
        return os; 
    }
};

// used by EM step only
void update_mean_map(map<int, map<int, pair<distribution2, histogram> > >& result, 
    layer1_result* res, int layer, part_lib* library, double rthresh, double gthresh);
void update_means(map<int, map<int, pair<distribution2, histogram> > >& result, int layer, part_lib* library);

////////////////////////////
// non-EM related

// used by object learning and layer optimization
void node_set_to_region_set(set<ipoint2>& result, const set<node*>& nset, const vector<set<ipoint2> >& rvector, int layer);


// used by get_maxima only
bool get_region_sym(matrix<double>& result, matrix<double>& m, int i0, int i1, int j0, int j1);
void max_normalize(matrix<double>& m, double factor = 1.0);

// used by layer optimization only
void select_nodes_by_type(vector<node*>& result, layer1_result* res, const vector<node*>& nodes, const set<int>& types);

// used by part and map learning during update
void get_region_map(map<node*, set<int> >& rmap, layer1_result* res, int layer, const vector<matrix<bool> >& regions, const vector<ipoint2>& rcenters);


// used by object learning only
void get_positions(set<ipoint2>& result, const set<node*>& s, int z);


#endif /* _LAYER_LEARNING_ */

