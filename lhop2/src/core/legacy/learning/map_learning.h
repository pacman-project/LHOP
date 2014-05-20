// layer_learning learning algorithms:
//   * s_learning     learning of layers (statistics, g-responses, em-step, optimization)
//   * obj_learning   learning of object layer (models)
//   * t_learning     learning of thresholds

#pragma once

#ifndef _MAP_LEARNING_
#define _MAP_LEARNING_

#include "core/legacy/layer_1_result.h"
#include "core/legacy/inference/layer_n_creators.h"

#include "core/legacy/learning/layer_learning.h"

#include "utils/serialization/streamed_pointer.h"


#include <tuple>
#include <stdio.h>

// map_learning
/////////////////

class part_learning;

class map_learning : public streamable {
protected:
    typedef matrix<double> stat_t;             // matrix of dimension stat_dim x stat_dim
    typedef map<ipoint2, stat_t> stat_map_t;         

    stat_map_t stat;    // statistics ("maps")
    int stat_dim;       // dimension of stat_t matrices

    vector<matrix<bool> > regions;  // regions of layer 1 masks -- initialized from library
    vector<ipoint2> rcenters;       // centers of regions

    // thresholds -- initialized from cfg
    int source_layer;           // index of source layer, mandatory

    double center_val_threshold;        // only nodes with val >= ~ are taken for centers;
                                        //   default val = 0.0 -- use ~_rel
    double center_val_threshold_rel;    // = 0.5; relative version of center_val_threshold, w/r to
    double nb_val_threshold_rel;        // = 0.6; map update only with neighbors of val: 
                                        //   cval*~ <= val <= cval*(2 - ~)
    double nbthresh_min, nbthresh_max;  // = 0.0, 0.0; i.e. use nb_val_threshold_rel
    double seq_min_intersection_percent;    // = 0.0
    double seq_max_intersection_percent;    // = 0.5     

    // finding maxima
    int max_max;    // = 4; max number of maxima
    double max_val_threshold;   // = 0.01; set all map values <= ~*max to 0 (min_update_count_percent)
    bool individual_max;        // = false; max (see map_val_threshold) is set for each map individually
    int max_nbhood_mask;        // = 5; size of "distribution"
    double max_sigma;           // = 0.0; if <= 0, then distribution is cut from the map
    int max_radius;             // = 2; min distance between two maxima and between each maximum and the center

public:
    typedef struct { ipoint2 pos; matrix<double> dist; } max_item_t;
    typedef vector<max_item_t> max_vector_t;     // vector of maxima (sorted by max values)
    typedef map<ipoint2, max_vector_t> max_map_t;

public:
    map_learning();
    map_learning(const ConfigDictionary& cfg);
    ~map_learning();

    void update(layer1_result* res);
    void reset(part_lib* library);
    void dispose();
    void cfg_init(const ConfigDictionary& cfg);

    void prepare_for_update(layer1_result* res);
    void display_statistics(const string& file_template, double blur_sigma = 0.75) const;
    void get_maxima(max_map_t& result) const;

    static void read_maxima_from_stream(max_map_t& maxima, vector<matrix<bool> >& regions, vector<ipoint2>& rcenters,
        const string& fname);

    friend class part_learning;
	
	void merge(map_learning* mlearner) {
		// all statistics is saved in stat variable - we can simply add sum all values at the same position
		
		for (auto iter = mlearner->stat.begin(); iter != mlearner->stat.end(); ++iter) {
			// for each value find matrix of this map_learning and sum it
			const ipoint2& iter_point = iter->first;
			stat_t& iter_stat = iter->second;

			stat_t& current_stat = stat[iter_point];
			if (current_stat.size() > 0)
				current_stat += iter_stat;
			else
				current_stat = iter_stat;
		}
	}

	// streamable implementations
	virtual streamable* make_instance() const { return new map_learning(); }
	virtual void read_from_stream(istreamer& is) {
		is.read(stat);
		is.read(stat_dim);
		is.read(regions);		
		is.read(rcenters);
		is.read(source_layer);
		is.read(center_val_threshold);
		is.read(center_val_threshold_rel);
		is.read(nb_val_threshold_rel);
		is.read(nbthresh_min);
		is.read(nbthresh_max);
		is.read(seq_min_intersection_percent);
		is.read(seq_max_intersection_percent);
		is.read(max_max);
		is.read(max_val_threshold);
		is.read(individual_max);
		is.read(max_nbhood_mask);
		is.read(max_sigma);
		is.read(max_radius);		
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write(stat);
		os.write(stat_dim);
		os.write(regions);
		os.write(rcenters);
		os.write(source_layer);
		os.write(center_val_threshold);
		os.write(center_val_threshold_rel);
		os.write(nb_val_threshold_rel);
		os.write(nbthresh_min);
		os.write(nbthresh_max);
		os.write(seq_min_intersection_percent);
		os.write(seq_max_intersection_percent);
		os.write(max_max);
		os.write(max_val_threshold);
		os.write(individual_max);
		os.write(max_nbhood_mask);
		os.write(max_sigma);
		os.write(max_radius);

	}
	virtual void copy_to(streamable* p, cloner& cl) 
    { 		
        map_learning* dest = (map_learning*)p;

		dest->stat = stat;
		dest->stat_dim = stat_dim;		
		dest->regions = regions;
		dest->rcenters = rcenters;
		dest->source_layer = source_layer;
		dest->center_val_threshold = center_val_threshold;
		dest->center_val_threshold_rel = center_val_threshold_rel;
		dest->nb_val_threshold_rel = nb_val_threshold_rel;
		dest->nbthresh_min = nbthresh_min;
		dest->nbthresh_max = nbthresh_max;
		dest->seq_min_intersection_percent = seq_min_intersection_percent;
		dest->seq_max_intersection_percent = seq_max_intersection_percent;
		dest->max_max = max_max;
		dest->max_val_threshold = max_val_threshold;
		dest->individual_max = individual_max;
		dest->max_nbhood_mask = max_nbhood_mask;
		dest->max_sigma = max_sigma;
		dest->max_radius = max_radius;
    }



protected:
    void get_maxima(max_vector_t& maxima, const matrix<double>& stat, double max_thresh) const;

};

#endif /* _MAP_LEARNING_ */

