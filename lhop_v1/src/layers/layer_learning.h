// layer_learning learning algorithms:
//   * s_learning     learning of layers (statistics, g-responses, em-step, optimization)
//   * obj_learning   learning of object layer (models)
//   * t_learning     learning of thresholds

#pragma once

#ifndef _LAYER_LEARNING_
#define _LAYER_LEARNING_

#include "layer_1_result.h"
#include "layer_n_creators.h"
#include <tuple>
#include <stdio.h>

// s_learning and structures
///////////////////////////////////////////////////////////////////////////////

struct nb_data_q {     // --> s_learning ???
    double n;
    int x, y;

    nb_data_q(double vn = 0, int vx = 0, int vy = 0) : 
        n(vn), x(vx), y(vy) { }
    nb_data_q(const nb_data_q& d) : n(d.n), x(d.x), y(d.y) { }
    
    bool operator<(const nb_data_q& d) const { return n < d.n; }
    bool operator>(const nb_data_q& d) const { return n > d.n; }

    friend ostream& operator<<(ostream& os, const nb_data_q& d) 
    { 
        os << '{' << d.n << ',' << d.x << ',' << d.y << '}' << ' ';
        return os; 
    }
};

struct nb_data {
    double n;
    int x, y;
    matrix<double> nb;
    //set<node*> rset;

    nb_data(double vn = 0, int vx = 0, int vy = 0, int msize = 3) : 
        n(vn), x(vx), y(vy), nb(msize, msize, 0) /*, rset()*/ { }
    nb_data(const nb_data_q& d, int msize = 3) :
        n(d.n), x(d.x), y(d.y), nb(msize, msize, 0) /*, rset()*/ { }
    nb_data(const nb_data& d) : n(d.n), x(d.x), y(d.y), nb(d.nb) /*, rset(d.rset)*/ { }
    
    bool operator<(const nb_data& d) const { return n < d.n; }
    bool operator>(const nb_data& d) const { return n > d.n; }

    friend ostream& operator<<(ostream& os, const nb_data& d) 
    { 
        os << '{' << d.n << ',' << d.x << ',' << d.y << '}' << ' ';
        return os; 
    }
};

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
    map_learning(const config_dictionary& cfg);
    ~map_learning();

    void update(layer1_result* res);
    void reset(part_lib* library);
    void dispose();
    void cfg_init(const config_dictionary& cfg);

    void prepare_for_update(layer1_result* res);
    void display_statistics(const string& file_template, double blur_sigma = 0.75) const;
    void get_maxima(max_map_t& result) const;

    static void read_maxima_from_stream(max_map_t& maxima, vector<matrix<bool> >& regions, vector<ipoint2>& rcenters,
        const string& fname);
    void write_maxima_to_stream(const string& fname);

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

	// streamable implementations TODO: finalize
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

// part_learning
///////////////////////////////////////////////////////////////////////////////

class online_path_map : public streamable {
protected:
    //typedef map<vector<int>, hdpoint_t> dpath_map_t;
    typedef pair<hdpoint_t, int> stat_t;
    typedef map<vector<int>, stat_t> map_t;

    map_t m;
public:
    online_path_map() : m() { }

    void update(const vector<int>& key, const hpoint_t& hp) 
    {
        stat_t& ps = m[key];
        double n = ++ps.second;

        if (n == 1) {
            ps.first.p = (dpoint2)hp.p;
            ps.first.h = hp.h;
        } else {
            dpoint2 deltap = (dpoint2)hp.p - ps.first.p;
            histogram_t deltah = hp.h - ps.first.h;
        

            ps.first.p += deltap/n;
            deltah /= n;
            ps.first.h += deltah;
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

    //struct hpoint_t {    ---> part_data_2; needed in library!
    //    ipoint2 p; 
    //    histogram_t h; 

    //    hpoint_t() : p(), h() { }
    //};


    typedef vector<ipoint2> sym_part_t;     // symbollic representation of a part
                                            //   x == type, y == position in maxima 
                                            // at[0] is center
    //typedef map<vector<int>, hpoint_t> edge_map_t;  
    typedef vector<path_map_t> geometry_t;
    typedef list<geometry_t> occurrences_t;

    struct stat_item_t { 
        double count; 
        //occurrences_t occ; 
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
    part_learning(const config_dictionary& cfg);
    ~part_learning();

    int get_source_layer() { return source_layer; }
    int get_stat_size() { return (int)stat.size(); }

    void init_cfg(const config_dictionary& cfg);
    void reset(const map_learning& mlearner);
    void reset(const string& fname);
    void reset();
    void dispose();
    void update(layer1_result* res, const map<ipoint2, histogram_t>& scmap);
    void display_maxima(const string& file_template, double blur_sigma = 0.75) const;
    void add_to_library(part_lib* library, int part_max_number);
    void add_to_library(part_lib* library, int part_max_number, int sorting_type, int cluster_size, double dthresh, double scthresh);
    void similarity_matrix(matrix<double>& m, int max_part_number, double thresh);
    void export_stat(const string& fname, int max_part_number);
    void merge_stat_sc(vector<vector<similarity_item> >& result, int max_part_number, int cluster_size,
        int space_size, int sorting_type, double ethresh, double scthresh);
    //void merge_library(part_lib* library, int layer, int max_part_number, int space_size, int cluster_size,
    //    double ethresh, double scthresh);
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
		
		//map_learning::max_map_t maxima;
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
    void make_path_map(path_map_t& m, const scmap_t& scmap, node* cn, node* n);
    path_map_t make_path_map(const scmap_t& scmap, node* cn, node* n);
    void make_path_map_rec(list<pair<vector<int>, node*> >& l);
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
    void update_geometry(node* p, const geometry_t& g);
    void update_pca(node* p, const pca_data& pcd);

    // support functions for similarity calculations
    vector<int> geometry_distance(double& gmatching, vector<dpoint2>& dvector, double& scmatching, 
        const geometry_t& igeo, const geometry_t& jgeo, bool calc_benergy);

};

// obj_learning
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
    obj_learning(const config_dictionary& cfg);
    ~obj_learning();

    void object_from_result(layer1_result* res, const scmap_t& scmap, const irectangle2& gtr);
    int add_to_library(part_lib* library, const list<obj_data_t>& objects, const string& name);
    int add_to_library(part_lib* library, const obj_data_t& object, const string& name);
    int add_to_library(const string& name);
	void add_validation_data(const streamed_pointer& ptr, const irectangle2& gtruth);
    //void add_validation_data(layer1_result* res, const irectangle2& gtruth);
    void set_creator(layern_creator* creator);
	void set_library(part_lib* plibrary);
protected:	
    void cfg_init(const config_dictionary& cfg);
    void add_additional_lib_nodes(part_lib* library, obj_data_t& obj, int m);
    void object_from_result(objects_t& result, layer1_result* res, const scmap_t& scmap,
        int layer, const set<ipoint2>& covered, const irectangle2& gtr);
    bool validate_object(layer1_result* res, const irectangle2& gtr, 
        const vector<node*>& cluster_nodes, const obj_data_t& object);
    bool validate_object(layer1_result* res, const irectangle2& gtr, const obj_data_t& object);
    bool validation_function_1(int tr, int fa);
    //void normalize_object_data(obj_data_t& od);
    void reduce_object_data(obj_data_t& od);
};

// o_learning
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
	config_dictionary inference_cfg;

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
    o_learning(const config_dictionary& cfg);
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
    void init_cfg(const config_dictionary& cfg);
    void make_duplets();
    void make_models(list<olv_object_part>& models, layer1_result* res);    
    void print_models(const list<olv_object_part>& models);
    void update_statmap(layer1_result* res, const list<irectangle2>& gtrs);
    void validate(vector<int>& tokeep, vectorn<int>::type2& rhits, part_lib* library, 
        int libsize0, int libsize1);
};

// ell_learning
///////////////////////////////////////////////////////////////////////////////

class ell_learning {
public:
    typedef float real_t;
    typedef vector<real_t> real_vector_t;
    typedef list<real_vector_t> list_t;

protected:
    // (object, no#, subpart#) -> positive/negative statistics vector 
    typedef triple<int, int, ipoint2> key_t;
    typedef map<key_t, ell_learning_stat_item> map_t;
    typedef map<string, int> cat_map_t;
    typedef map<int, int> index_map_t;

    int category_layer; 
    double cover_threshold;
    string gtextension;
    part_lib* library;

    cat_map_t catmap;
    index_map_t posindexmap, negindexmap;
    map_t posstat, negstat;
public:
    ell_learning(const config_dictionary& cfg);
    ~ell_learning();

    void update(layer1_result* res, const list<pair<irectangle2, string> >& gtr);
    void update(layer1_result* res, const list<pair<irectangle2, int> >& gtr);
    part_lib* get_library() { return library; }
    void display_stat();
    void learn(list_t& positive, list_t& negative, int type);
    //void learn(const string& fname, int type);
    void learn_all(string out_dir);
protected:
    void learn(list_t& features, int& feature_count, int type, const map_t& stat);
    void cfg_init(const config_dictionary &cfg);
    void add_to_map(bool positive, int type, int index, const ipoint2& name, const ellipse& ell);
    void init_catmap();
    int get_new_sample_index(bool positive, int type);
};

// t_learning (threshold learning)
///////////////////////////////////////////////////////////////////////////////

class t_learning {
protected:
    // parameters
    int category_layer;           // mandatory
    bool use_groundthruth;        // = true
    bool false_is_negative;       // = true
    bool object_layer_only;       // = false
    double groundtruth_threshold; // = 0.7       "|intersection|/|union|"
    double tf_ratio;              // = 0.75
    bool remove_bad_parts;        // = true
	string gtextension;
    part_lib* library;

    set<int> toremove;
public:
    t_learning(const config_dictionary& cfg);
    ~t_learning();

    part_lib* get_library() { return library; }

    void learn(string pdir, const list<string>& positive, 
        string ndir, const list<string>& negative, string& sm, int threshold_name, int threshold_type);
    void begin();
    void end();

protected:
    void cfg_init(const config_dictionary& cfg);

    //void get_object_parts(set<int>& result, string& sm);
    double calc_obj_threshold(const vector<double>& pos, const vector<double>& neg);
    double calc_obj_neg_threshold(const vector<double>& pos, const vector<double>& neg);
    void update_tree_thresholds(layer1_result* res, node* n);
    void update_tree_thresholds(layer1_result* res, const vector<node*>& nodes);
};


struct clique_data {
    set<int> parts;
    int part0;
    cv::Mat mean;
    cv::Mat eigenvectors;
    cv::Mat eigenvalues;
    double sizefactor;
};

// pca_merging
///////////////////////////////////////////////////////////////////////////////

class pca_merging {
protected:
    typedef map<int, shape_learning_p> statmap_t;

    statmap_t stat;
    int layer;      // Layer 
    int iradius;    // Radius for inhibition
public:
    pca_merging(int l, int r);

    void update_stat(layer1_result* res, const response_filter& filter, double area_thresh);
    void update_library(part_lib* library, double mergenorm);
};


// global functions
///////////////////////////////////////////////////////////////////////////////

// g-response learning
////////////////////////

void update_g_distribution(map<pair<int, int>, online_distribution>& result, layer1_result* res, int layer, part_lib* library);

void set_g_distribution(part_lib* library, int layer, map<pair<int, int>, online_distribution>& result, double minvar, bool overwrite);

void learn_g_distributions(part_lib* library, int layer, const config_dictionary& srccfg, const string& nspace,
    const string& dir, const string& patt);

// em-step
////////////

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
    rmatrix cut_matrix(const ipoint2& p, int delta, double defval = 0.0);
    rmatrix to_matrix(double defval = 0.0);
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

void update_mean_map(map<int, map<int, pair<distribution2, histogram> > >& result, 
    layer1_result* res, int layer, part_lib* library, double rthresh, double gthresh);

void update_means(map<int, map<int, pair<distribution2, histogram> > >& result, int layer, part_lib* library);

// optimization
///////////////////////////////////////////////////////////////////////////////

// optmization_set
////////////////////

/// Used in the optimization process and it should be able to return layer1_result
/// using 'get' function. 
struct optimization_set {
    optimization_set() { }
    virtual ~optimization_set() { }
    virtual int size() const = 0;
    virtual layer1_result* get(int i) = 0;
};


// optimization_set_streamed
//////////////////////////////

/// Internally used by optimization process
struct optimization_set_streamed : public optimization_set {
protected:
    vector<streamed_pointer> container;
    layer1_result* current;
public:
    optimization_set_streamed() : container(), current(nullptr) { } 
    virtual ~optimization_set_streamed() { if (current != nullptr) delete current; }
    void add(const streamed_pointer& sp) { container.push_back(sp); }
    void clear() { container.clear(); }
    virtual int size() const { return (int)container.size(); }
    virtual layer1_result* get(int i) 
    { 
        //if (current != nullptr) delete current;
        //return (current = (layer1_result*)container[i].get()); 
		return ((layer1_result*)container[i].get()); 
    }
};



set<int> optimize_layer(part_lib* library, list<streamed_pointer>& workingset, int layer, const set<int>& init_set,
    double cover_thresh, double int_thresh, int bite_size);

set<int> optimize_layer_2(part_lib* library, list<streamed_pointer>& workingset, int layer, const set<int>& init_set,
    double cover_thresh, double int_thresh, int bite_size);

set<int> optimize_layer_locally(part_lib* library, list<streamed_pointer>& workingset, int layer, const set<int>& init_set,
    int n_parts, double gr_thresh, double s_thresh, int bite_size);

void add_similarity_edges(part_lib* library, int layer, int max_geo_size, double thresh);

#endif /* _LAYER_LEARNING_ */

