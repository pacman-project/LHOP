// layer_learning learning algorithms:
//   * s_learning     learning of layers (statistics, g-responses, em-step, optimization)
//   * obj_learning   learning of object layer (models)
//   * t_learning     learning of thresholds

#pragma once

#ifndef _PART_LEARNING_
#define _PART_LEARNING_

#include "core/legacy/layer_1_result.h"
//#include "core/inference/layer_n_creators.h"

#include "core/legacy/learning/layer_learning.h"
#include "core/legacy/learning/map_learning.h"

#include "utils/serialization/streamed_pointer.h"

#include <tuple>
#include <stdio.h>

// part_learning
///////////////////////////////////////////////////////////////////////////////

class map_learning;

class online_path_map : public streamable {
protected:
    typedef pair<hdpoint_t, int> stat_t;
    typedef map<vector<int>, stat_t> map_t;

    map_t m;
public:
    online_path_map() : m() { }

    void update(const vector<int>& key, const hpoint_t& new_value) 
    {
        stat_t& existing_value = m[key];
        double n = ++existing_value.second;

        if (n == 1) {
            existing_value.first.p = (dpoint2)new_value.p;
            existing_value.first.h = new_value.h;
        } else {
            dpoint2 deltap = (dpoint2)new_value.p - existing_value.first.p;
            histogram_t deltah = new_value.h - existing_value.first.h;
        

            existing_value.first.p += deltap/n;
            deltah /= n;
            existing_value.first.h += deltah;
        }
    }

    void update(const path_map_t& pm) 
    {
        for (auto pmiter = pm.begin(); pmiter != pm.end(); ++pmiter) 
            update(pmiter->first, pmiter->second);
    }

    void get_mean(path_map_t& pm)
    {
        pm.clear();
        for (auto miter = m.begin(); miter != m.end(); ++miter) {
            hpoint_t& hp = pm[miter->first];

            hp.p = int_round(miter->second.first.p);
            hp.h = miter->second.first.h;
        }
    }

	void merge(const online_path_map& opm) 
	{
		for (map_t::const_iterator opmiter = opm.m.begin(); opmiter != opm.m.end(); ++opmiter) {
			const vector<int>& srckey = opmiter->first;
			const stat_t& srcps = opmiter->second;
	        stat_t& ps = m[srckey];

			if (ps.second == 0) {
				ps.first = srcps.first;
				ps.second = srcps.second;
			} else {
				int n = ps.second;
				int m = srcps.second;
				double fn = (double)n/(m + n);
				double fm = (double)m/(m + n);
				int maxi = ps.first.h.size();

				for (int i = 0; i < maxi; ++i) {
					ps.first.h[i] = fn*ps.first.h[i] + fm*srcps.first.h[i];
				}
				ps.first.p = ps.first.p*fn + srcps.first.p*fm;
				ps.second = m + n;
			}
		}
	}
	void print() {
		for (auto iter = m.begin(); iter != m.end(); ++iter) {
			printf("[");
			for (auto iter_1 = iter->first.begin(); iter_1 != iter->first.end(); ++iter_1) 
				printf("%d,",*iter_1);
			
			printf("]=");			
			printf("(%d,(%f,%f,[", iter->second.second, iter->second.first.p.x, iter->second.first.p.y);
			for (auto iter_1 = iter->second.first.h.begin();  iter_1 != iter->second.first.h.end(); ++iter_1) 
				printf("%f,",*iter_1);
			printf("])),");
		}
	}
	virtual streamable* make_instance() const { return new online_path_map(); }
	virtual void read_from_stream(istreamer& is) {
		int m_size;
		is.read(m_size);
		for (int i = 0; i < m_size; i++) {
			vector<int> key;
			stat_t value;		
			is.read(key); // vector<int>
			
			is.read(value.first.p); // dpoint2
			is.read(value.first.h); // histogram_t == vector<double>
			is.read(value.second); // int

			m[key] = value;
		}
	}

	virtual void write_to_stream(ostreamer& os) {
		os.write((int)m.size());
		for (auto iter = m.begin(); iter != m.end(); ++iter) {
			os.write(iter->first); // vector<int>
			os.write(iter->second.first.p); // dpoint2
			os.write(iter->second.first.h); // histogram_t == vector<double>
			os.write(iter->second.second); // int
		}
	}
};

class online_geo_learning : public streamable {
protected:
    vector<online_path_map> geo;
public:
    online_geo_learning() : geo() { geo.reserve(5); }

    void update(int i, const path_map_t& pm) 
    {
        if (i >= geo.size()) geo.resize(i + 1);
        geo[i].update(pm);
    }

    void get_mean(path_map_t& pm, int i)
    {
        geo[i].get_mean(pm);
    }

	void merge(const online_geo_learning& ogl) 
	{
		if (geo.size() < ogl.geo.size()) 
			geo.resize(ogl.geo.size());
		for (int i = 0; i < ogl.geo.size(); ++i)
			geo[i].merge(ogl.geo[i]);
	}
	void print() {
		for (auto iter = geo.begin(); iter != geo.end(); ++iter) 
			iter->print();
	}
	virtual streamable* make_instance() const { return new online_geo_learning(); }
	virtual void read_from_stream(istreamer& is) {
		int geo_size;
		is.read(geo_size);
		
		geo.resize(geo_size);
		for (int i = 0; i < geo_size; i++) {
			geo[i].read_from_stream(is);
		}
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write((int)geo.size());
		for (auto iter = geo.begin(); iter != geo.end(); ++iter)
			iter->write_to_stream(os);
	}
};

class part_learning : public streamable {
protected:
    enum { SC_SIMILARITY = 0, PCA_SIMILARITY = 1 };

public:
    struct candidates_item_t { 
        double val; ipoint2 coo; node* n; 
        candidates_item_t(double v = 0.0, const ipoint2& c = ipoint2::zero, node* pn = nullptr) : val(v), coo(c), n(pn) { }
        bool operator<(const candidates_item_t& ci) const { return val < ci.val; }
        bool operator>(const candidates_item_t& ci) const { return val > ci.val; }
    };

    struct real_part_t { 
        double val;
        vector<node*> nodes;

        real_part_t() : val(0.0), nodes() { }
    };

    typedef vector<double> histogram_t;

    struct similarity_item {
        int index;
        double val;
        vector<int> permutation;

        similarity_item() : index(0), val(0.0), permutation() { }
        similarity_item(int i, int v, const vector<int>& p) : index(i), val(v), permutation(p) { }
    };
	
    typedef vector<ipoint2> sym_part_t;     // symbollic representation of a part
                                            //   x == type, y == position in maxima 
                                            // at[0] is center
    
	typedef vector<path_map_t> geometry_t;
    typedef list<geometry_t> occurrences_t;

    struct stat_item_t { 
        double count;  // used by sort_stat
        online_geo_learning geol;
        pca_learning pcal;
        
        stat_item_t() : count(0.0), geol(), pcal() { }
    };

    typedef map<sym_part_t, stat_item_t> stat_t;

    typedef list<candidates_item_t> candidates_t;

protected:
    stat_t stat;
    map_learning::max_map_t maxima;

    vector<matrix<bool> > regions;  // regions of layer 1 masks -- initialized from map_learner
    vector<ipoint2> rcenters;       // centers of regions
    double layer_contraction;       // calculated in update, def. value is -1.0

    // settings and thresholds -- initialized from cfg file
    int source_layer;       // mandatory

    // finding maxima
    int max_max;    // = 4; max number of maxima
    double map_val_threshold;   // = 0.01; set all map values <= ~*max to 0 (min_update_count_percent)
    bool individual_max;        // = false; max (see map_val_threshold) is set for each map individually

    // part update
    double center_val_threshold;            // = -1.0
    double center_val_threshold_rel;        // = 0.5
    int min_seq_size;                       // = 2
    int max_seq_size;                       // = 3
    double seq_min_intersection_percent;    // = 0.0
    double seq_max_intersection_percent;    // = 0.5     
    int max_candidates;                     // = 5
    int max_stat_size;                      // = INT_MAX

    // similarity calculations
    int similarity_type;                    // SC_SIMILARITY or PCA_SIMILARITY
                                            // default value is determined whether there is a part with
                                            // vs_part_data on source_layer
    int max_matching_size;                  // = 30;  to invoke inhibition in 'geometry_distance'
public:
    part_learning();
    part_learning(const ConfigDictionary& cfg);
    ~part_learning();

    int get_source_layer() { return source_layer; }
    int get_stat_size() { return (int)stat.size(); }

    void init_cfg(const ConfigDictionary& cfg);
    void reset(const map_learning& mlearner);
    void reset(const string& fname);
    void reset();
    void dispose();
    void update(layer1_result* res, const map<ipoint2, histogram_t>& scmap);
    void display_maxima(const string& file_template, double blur_sigma = 0.75) const;

    void add_to_library(part_lib* library, int part_max_number, int sorting_type, int cluster_size, double dthresh, double scthresh);

	void save_stat();

	void merge(part_learning* plearning) {
		for (auto iter = plearning->stat.begin(); iter != plearning->stat.end(); ++iter) {
			const sym_part_t& iter_points = iter->first;
			stat_item_t& iter_stat = iter->second;
			
			stat_item_t& current_stat = stat[iter_points];
			current_stat.count += iter_stat.count;
			current_stat.geol.merge(iter_stat.geol);
			current_stat.pcal.merge(iter_stat.pcal);
		}
		layer_contraction = plearning->layer_contraction;
	}

	// streamable implementations
	virtual streamable* make_instance() const { return new part_learning(); }
	virtual void read_from_stream(istreamer& is) {
		is.read(regions);
		is.read(rcenters);
		is.read(layer_contraction);
		is.read(source_layer);
		is.read(max_max);
		is.read(map_val_threshold);
		is.read(individual_max);
		is.read(center_val_threshold);
		is.read(center_val_threshold_rel);
		is.read(min_seq_size);
		is.read(max_seq_size);
		is.read(seq_min_intersection_percent);
		is.read(seq_max_intersection_percent);
		is.read(max_candidates);
		is.read(max_stat_size);
		is.read(similarity_type);
		is.read(max_matching_size);

		// read maxima
		int maxima_size;
		is.read(maxima_size);
		for (int i = 0; i < maxima_size; i++) {
			ipoint2 key;
			map_learning::max_vector_t value;
			int value_size;

			is.read(key);
			is.read(value_size);

			maxima[key] = map_learning::max_vector_t(value_size);
			map_learning::max_vector_t& item_vector = maxima[key];
			for (int j = 0; j < value_size; j++) {
				map_learning::max_item_t& item = item_vector[j];
				is.read(item.pos);
				is.read(item.dist);
			}
		}

		// write stat
		int stat_count;
		is.read(stat_count);
		for (int i = 0; i < stat_count; i++) {
			sym_part_t key;
			stat_item_t item;
			is.read(key);

			stat[key] = item;

			stat_item_t &item_ref = stat[key];
			
			is.read(item_ref.count);
			item_ref.geol.read_from_stream(is);
			item_ref.pcal.read_from_stream(is);
		}
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write(regions);
		os.write(rcenters);
		os.write(layer_contraction);
		os.write(source_layer);
		os.write(max_max);
		os.write(map_val_threshold);
		os.write(individual_max);
		os.write(center_val_threshold);
		os.write(center_val_threshold_rel);
		os.write(min_seq_size);
		os.write(max_seq_size);
		os.write(seq_min_intersection_percent);
		os.write(seq_max_intersection_percent);
		os.write(max_candidates);
		os.write(max_stat_size);
		os.write(similarity_type);
		os.write(max_matching_size);

		// write maxima
		os.write((int)maxima.size());
		for (auto iter = maxima.begin(); iter != maxima.end(); ++iter) {
			os.write(iter->first);
			os.write((int)iter->second.size());
			
			for (auto iter_1 = iter->second.begin(); iter_1 != iter->second.end(); ++iter_1) {
				os.write(iter_1->pos);
				os.write(iter_1->dist);
			}
		}

		// write stat
		os.write((int)stat.size());
		for (auto iter = stat.begin(); iter != stat.end(); ++iter) {
			os.write(iter->first);
			
			os.write(iter->second.count);
			iter->second.geol.write_to_stream(os);
			iter->second.pcal.write_to_stream(os);
		}
	}

	virtual void copy_to(streamable* p, cloner& cl) 
    { 		
        part_learning* dest = (part_learning*)p;

		dest->stat = stat;
		dest->maxima = maxima;
		dest->regions = regions;
		dest->rcenters = rcenters;
		dest->layer_contraction = layer_contraction;
		dest->source_layer = source_layer;
		dest->max_max = max_max;
		dest->map_val_threshold = map_val_threshold;
		dest->individual_max = individual_max;
		dest->center_val_threshold = center_val_threshold;
		dest->center_val_threshold_rel = center_val_threshold_rel;
		dest->min_seq_size = min_seq_size;
		dest->max_seq_size = max_seq_size;
		dest->seq_min_intersection_percent = seq_min_intersection_percent;
		dest->seq_max_intersection_percent = seq_max_intersection_percent;
		dest->max_candidates = max_candidates;
		dest->max_stat_size = max_stat_size;
		dest->similarity_type = similarity_type;
		dest->max_matching_size = max_matching_size;

	}
	void print() {	
		printf("layer_contraction is %f\n",layer_contraction);
		printf("source_layer is %d\n",source_layer);
		printf("max_max is %d\n",max_max);
		printf("map_val_threshold is %f\n",map_val_threshold);
		printf("individual_max is %d\n",individual_max);
		printf("center_val_threshold is %f\n",center_val_threshold);
		printf("center_val_threshold_rel is %f\n",center_val_threshold_rel);
		printf("min_seq_size is %d\n",min_seq_size);
		printf("max_seq_size is %d\n",max_seq_size);
		printf("seq_min_intersection_percent is %f\n",seq_min_intersection_percent);
		printf("seq_max_intersection_percent is %f\n",seq_max_intersection_percent);
		printf("max_candidates is %d\n",max_candidates);
		printf("max_stat_size is %d\n",max_stat_size);
		printf("similarity_type is %d\n",similarity_type);
		printf("max_matching_size is %d\n",max_matching_size);
		
		printf("rcenters of size %d:\n\t", rcenters.size()); for (auto iter = rcenters.begin(); iter != rcenters.end(); ++iter) printf("(%d,%d), ", (*iter).x, (*iter).y);
		printf("regions of size %d:", regions.size()); for (auto iter = regions.begin(); iter != regions.end(); ++iter) { printf("\n\t"); iter->print(); }
		
		printf("stat of size %d:", stat.size()); 
		for (auto iter = stat.begin(); iter != stat.end(); ++iter) {
			printf("item key: "); for (auto iter_1 = iter->first.begin(); iter_1 != iter->first.end(); ++iter_1) printf("(%d,%d), ", (*iter_1).x, (*iter_1).y);
			printf("item value:\n\t");
			printf("\tcount = %f\n",iter->second.count);
			printf("\tgeol = "); iter->second.geol.print();
			printf("\tpcal = "); iter->second.pcal.print();
		}
    
		printf("maxima of size: %d\n", maxima.size());
		for (auto iter = maxima.begin(); iter != maxima.end(); ++iter) {
			printf("%d,%d = [",iter->first.x,iter->first.y);
			for (auto iter_1 = iter->second.begin(); iter_1 != iter->second.end(); ++iter_1) {
				printf("(%d,%d,[",iter_1->pos.x,iter_1->pos.y);
				for (auto iter_2 = iter_1->dist.begin();iter_2 != iter_1->dist.end(); ++iter_2) {
					printf("%f, ", *iter_2);
				}
				printf("]),");
			}
			printf("]\n");
		}
	}
protected:
    void get_best_candidates(candidates_t& bestc, layer1_result* res, 
        const map<node*, set<int> >& rmap, node* cn, const set<int>& fr);
    void extend_sequence(layer1_result* res, const scmap_t& scmap,
        const map<node*, set<int> >& rmap, node* cn, 
        const sym_part_t& seq, const real_part_t& rseq, const set<int>& fr);
    void update_sequence(const sym_part_t& seq, const real_part_t& rseq, const scmap_t& scmap, node* cn);
	
    void sort_stat(vector<pair<double, const sym_part_t*> >& pairs, int sorting_type);
    void merge_stat_sc(vector<vector<similarity_item> >& result, 
        map<sym_part_t, vector<int> >& perm, 
        map<sym_part_t, geometry_t>& stat2,
        int max_part_number, int sorting_type, int cluster_size, int space_size, double ethresh, double scthresh);
    void merge_stat_pca(vector<vector<similarity_item> >& result, 
        map<sym_part_t, vector<int> >& permap,
        map<sym_part_t, pca_data>& stat3,
        int space_size, int sorting_type, double norm);
    void get_average_stat(map<sym_part_t, geometry_t>& stat2);
    void get_pca_stat(map<sym_part_t, pca_data>& stat3);


};
#endif /* _PART_LEARNING_ */

