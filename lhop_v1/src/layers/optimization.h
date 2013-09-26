// optimization classes

#pragma once

#ifndef __OPTIMIZATION_H__
#define __OPTIMIZATION_H__

#include "../layers/layer_1_result.h"
#include "../layers/layer_n_creators.h"
#include "../layers/layer_learning.h"
#include "../utils/mapreduce.h"

// defines
///////////////////////////////////////////////////////////////////////////////

//#define _SAVEMAPS
#define PRINT_INFO(_info) cout << _info << endl;
#define PRINT_INFO_NUM(_info) cout << _info << ' ';
#define PRINT_DOT() cout << '.';
//#define PRINT_INFO(_info) 
//#define PRINT_DOT()


// classes
///////////////////////////////////////////////////////////////////////////////

//template<class V, class T> class local_optimization {
//public:
//    virtual V value(const T& state) = 0;
//    virtual void get_neighbor(const T& state, T& nbr) = 0;
//    virtual void move_to_neighbor(const T& state, const T& ) = 0;
//};

// layer1_optimization_new
///////////////////////////////////////////////////////////////////////////////

/// Internally used by optimization process
//struct optimization_set_sptr : public optimization_set {
//protected:
//    vector<layer1_result_ptr> container;
//    //layer1_result* current;
//public:
//    optimization_set_sptr() : container() { } 
//    virtual ~optimization_set_sptr() { }
//    void add(layer1_result_ptr sp) { container.push_back(sp); }
//    void clear() { container.clear(); }
//    virtual int size() const { return (int)container.size(); }
//    virtual layer1_result* get(int i) { return (layer1_result*)container[i]; }
//};
//

class layer1_optimization_new {
protected:
    part_lib* library;
    int layer;
    //optimization_set_streamed workingset;
	list<streamed_pointer> workingset;
    int parts_to_keep;
    set<int> final_parts;
    int loop;

    vector_parameter<double> cover_thresh;      // cfg: "covered_threshold", def. value: 0.85
    vector_parameter<double> int_thresh;        // cfg: "intersection_threshold", def. value: 0.2
    vector_parameter<int> bite_size;            // cfg: "bite_size", def. value: 10
    //vector_parameter<double> angle_cos_thresh;  // cfg: "angle_cos_threshold", def. value: 0.97
    //vector_parameter<double> a_thresh;          // cfg: "a_threshold", def. value 3
    //vector_parameter<double> e_thresh;          // cfg: "e_threshold", def. value 0.1
    vector_parameter<int> optimization;                           // cfg: "optimize", 0: no optimization, 1: global optimization (default)
                                                //          2: global optimization of type 2
    string library_image_file;
    bool show_labels;

public:
    layer1_optimization_new(const config_dictionary& cfg, int player);

    /// Optimization on or off?
    bool on() { return optimization != 0; }

    /// Returns the layer being optimized.
    int get_layer() { return layer; }

    /// Returns the library being optimized.
    part_lib* get_library() { return library; }

    /// Does nothing...
    void set_steps(int s) {  }

    /// Sets the library and the number of existing parts to keep.
    void set_library(part_lib* plibrary, int keep);

    /// Resets the structure (except of the 'library')
    /// and adds elements to 'test_set'
    /// See protected member for more details.
    template <class I> void add_to_test_set(I begin, I end);

    /// Perform optimization (calls optimize_layer function).
    /// The result is saved to final_parts.
    /// Returns 0.0 (returning double is for compatibility with the other version of 
    ///    layer1_optimization).
    double make_steps();

    /// Keeps only final_parts in the library and in the working set.
    void keep_final_parts();

    /// Print final_parts to cout.
    void print_final_parts();

    /// Saves graphic representation of the library (to library_image_file).
    void display_final_parts(const string& file, bool show_labels);
    void display_final_parts();

    /// Increments loop counter and parameters.
    void inc_loop();

protected:

    /// Fill parameters from cfg
    void init(const config_dictionary& cfg);

    /// Adds image to optimization_set
    void add_to_test_set(streamed_pointer tres);

    /// Clears the data set
    void reset();
};

template <class I> void layer1_optimization_new::add_to_test_set(I begin, I end)
{
    reset();
    for (; begin != end; ++begin) 
        add_to_test_set(*begin);
}

// cr_optimization
///////////////////////////////////////////////////////////////////////////////////

class optimization_data;

/* 
 *  Usage:
 *  - set the following (possibly in different order)
 *     ~ creators
 *     ~ learners
 *     ~ optimizers
 *     ~ 'start_layer' & 'end_layer'
 *  - set library (possibly delete everything above 'start_layer') !!!
 *  - add images using 'add_to_test_set'
 *  - call inter_layer_optimization
 *  - the result is the new library and elements of the 'test_set'
 *   
 */

/* Implementation of inter_layer_optimization
   - layer_preparation

*/

// cr_optimization_base
///////////////////////////////////////////////////////////////////////////////

class cr_optimization_base {

public:
    cr_optimization_base() : mapreduce_deployer(nullptr) { }

    // Should be called whenever the optimizer is ready
    // to accept new library and new test set
    virtual void reset() { }

    // Set library and optionally (if 'keep_parts' == false) 
    // delete all parts >= 'start_layer' in the library
    virtual void set_library(part_lib* plibrary, bool keep_parts = true) = 0;

    // Adds 'res' to 'init_test_set'
    // Makes a copy and deletes everything above 'start_layer'
    virtual void add_to_test_set(streamed_pointer res) = 0;
    template <class I> void add_range_to_test_set(I begin, I end);

    virtual part_lib* get_best_library() = 0;
    virtual void get_best_test_set(list<streamed_pointer>& container) = 0;
    virtual double execute() = 0;
    virtual int get_start_layer() = 0;
    //virtual void get_initial_parts(vector<int>&) = 0;

    virtual void print_best_parts() { }

	virtual void set_mapreduce(base_deployer_mapreudce* mapreduce) {
		this->mapreduce_deployer = mapreduce;
	}
protected:

	base_deployer_mapreudce* mapreduce_deployer;

    // Virtual methods for reaching "current objects".
    //virtual int get_current_layer() const = 0;
    //virtual part_lib* get_current_library() const = 0;
    //virtual list<layer1_result_ptr>& get_current_train_set() = 0;
    //virtual layer1_optimization* get_current_optimizer() const = 0;
    //virtual layern_creator* get_current_creator() const = 0;
    //virtual s_learning* get_current_learner() const = 0;
};

template <class I> void cr_optimization_base::add_range_to_test_set(I begin, I end)
{
    for (; begin != end; ++begin) 
        add_to_test_set(*begin);
}

// cr_set_optimization
///////////////////////////////////////////////////////////////////////////////

class cr_set_optimization : public cr_optimization_base {
protected:
    list<streamed_pointer> init_test_set;  
    list<streamed_pointer> current_train_set;
    list<streamed_pointer> best_test_set;
    list<cr_optimization_base*> optimizers;
    part_lib* init_library;
    part_lib* current_library;
    part_lib* best_library;
    int steps;

public:
    cr_set_optimization(optimization_data* oopt, int nsteps, const vector<cr_optimization_base*>& oset);

    virtual void reset();
    virtual void set_library(part_lib* plibrary, bool keep_parts = true);
    virtual void add_to_test_set(streamed_pointer res);

    virtual part_lib* get_best_library() { return best_library; }
    virtual void get_best_test_set(list<streamed_pointer>& container) 
    { 
        container.assign(best_test_set.begin(), best_test_set.end()); 
    }

    virtual int get_start_layer() { return optimizers.front()->get_start_layer(); }
    virtual void get_initial_parts(vector<int>& parts) { parts.clear(); }
    virtual double execute();

    // virtual void print_parts();
    void add_to_init_test_set(optimization_data* oopt);

	virtual void set_mapreduce(base_deployer_mapreudce* mapreduce) {
		for (auto iter = optimizers.begin(); iter != optimizers.end(); ++iter)
			(*iter)->set_mapreduce(mapreduce);
	}
};

// cr_layer_optimization
///////////////////////////////////////////////////////////////////////////////

class cr_layer_optimization : public cr_optimization_base {
protected:
    list<streamed_pointer> current_train_set;
    layer1_optimization_new* optimizer;
    layern_creator* creator;
    map_learning* mlearner;
    class part_learning* plearner;
    part_lib* current_library;

    // parameters of the optimization process
    int steps, loops;
    int uncovered_radius;
    double init_types_threshold;
	int em_steps;
    vector_parameter<int> part_max_number;
    vector_parameter<int> sorting_type; // 0 = occurrence sorting (default), 1 = area sorting
    int cluster_size;
    double merge_distance_threshold;
    double merge_sc_threshold;

    string maxima_images;
    string map_images;
    string library_image;
    string library_image_sc;
    string lib_export_name;
    bool show_labels;
    bool display_stat;
    bool video;

    // status of the optimization process
    int current_loop;

public:
    cr_layer_optimization(optimization_data* oopt, const config_dictionary& cfg, int layer);

    virtual void reset();
    virtual void set_library(part_lib* plibrary, bool keep_parts = true);
    virtual void add_to_test_set(streamed_pointer res);
    virtual part_lib* get_best_library() { return current_library; }
    virtual void get_best_test_set(list<streamed_pointer>& container);
    virtual int get_start_layer() { return optimizer->get_layer(); }
    //virtual void get_initial_parts(vector<int>& parts);
    virtual double execute();
    // virtual void print_parts();
	
	static void prepare_result(layer1_result* res, int start_layer);

protected:
    //void get_initial_parts(set<int>& parts, layer1_result* res);
    //void get_initial_parts(vector<int>& parts, double types_threshold);
    
    void part_learning();
    void g_distribution_learning();
    void optimization(int parts_to_keep);
    void EM_step(int parts_to_keep);
};


// optimization_data
///////////////////////////////////////////////////////////////////////////////////

class optimization_data {
public:
    typedef list<streamed_pointer>::iterator iterator_t;
protected:

    vector<layer1_optimization_new*> optimizers;
    vector<layern_creator*> creators;
    vector<map_learning*> mlearners;
    vector<part_learning*> plearners;
    part_lib* library;

    list<streamed_pointer> test_set;

public:    
	optimization_data(const config_dictionary& cfg, part_lib* lib);
	virtual ~optimization_data();

    void add_to_test_set(streamed_pointer res);

    layer1_optimization_new* get_optimizer(int layer) { return optimizers[layer]; }
    layern_creator* get_creator(int layer) { return creators[layer]; }
    map_learning* get_map_learner(int layer) { return mlearners[layer]; }
    part_learning* get_part_learner(int layer) { return plearners[layer]; }
    part_lib* get_library() { return library; }
    iterator_t get_test_set_begin() { return test_set.begin(); }
    iterator_t get_test_set_end() { return test_set.end(); }
    bool test_set_empty() { return test_set.empty(); }

protected:
    void initialize(const config_dictionary& cfg);

};

// public functions
///////////////////////////////////////////////////////////////////////////////

typedef map<int, pair<int, int> > validation_result_t;

void validate_parts(validation_result_t& result, layer1_result* res, const set<int>& parts,
    const list<irectangle2>& gtruth, int layer, double thresh);
    


class map_learning_mapreduce : public base_mapreduce {
	map_learning* mlearner;

	bool dispose_mlearner;
public:
	map_learning_mapreduce() : mlearner(nullptr), dispose_mlearner(1) {}
	virtual ~map_learning_mapreduce() {
		if (dispose_mlearner && mlearner != nullptr) delete mlearner;
	}
	////////////////////////////////////////////////
	// main functionality is written in map and reduce
	virtual streamable* map(streamable* item);
	virtual streamable* reduce(list<streamable*> &item_list);
	
	virtual string map_get_key(streamable* item) { return "merged_result"; }

	////////////////////////////////////////////////
	// support functions
	void set_map_learning(map_learning* m) { this->mlearner = m; dispose_mlearner = false; }

	// streamable implementations
	virtual streamable* make_instance() const { return new map_learning_mapreduce(); };
	virtual void read_from_stream(istreamer& is) {
		streamable* ptr;
		is.read(ptr); mlearner = (map_learning*)ptr;
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write(mlearner);
	}
};

class part_learning_mapreduce : public base_mapreduce {
	K_bin* bin;
	
	part_learning* plearner;
	layern_creator* creator;

	int current_loop;
	int layer;

	bool dispose_bin;
	bool dispose_plearner;
	bool dispose_creator;
public:
	part_learning_mapreduce() : plearner(nullptr), creator(nullptr), bin(nullptr), current_loop(0),layer(0), dispose_plearner(1), dispose_creator(1), dispose_bin(1) {}
	virtual ~part_learning_mapreduce() {
		if (dispose_bin && bin != nullptr) delete bin;
		if (dispose_plearner && plearner != nullptr) delete plearner;
		if (dispose_creator && creator != nullptr) delete creator;
	}
	
	////////////////////////////////////////////////
	// main functionality is written in map and reduce
	virtual streamable* map(streamable* item) ;
	virtual streamable* reduce(list<streamable*> &item_list);

	virtual string map_get_key(streamable* item) { return "merged_result"; }

	////////////////////////////////////////////////
	// support functions
	void set_bin(K_bin* b) { this->bin = b; dispose_bin = false; }
	void set_part_learning(part_learning* p) { this->plearner = p; dispose_plearner = false; }
	void set_creator(layern_creator* c) { this->creator = c; dispose_creator = false;}
	void set_current_loop(int c) { this->current_loop = c; }
	void set_layer(int l) { this->layer = l; }

	// streamable implementations
	virtual streamable* make_instance() const { return new part_learning_mapreduce(); }
	virtual void read_from_stream(istreamer& is) {
		streamable* ptr;

		is.read(ptr); plearner = (part_learning*)ptr;
		is.read(ptr); creator = (layern_creator*)ptr;
		is.read(ptr); bin = (K_bin*)ptr;
		is.read(current_loop);
		is.read(layer);
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write(plearner);
		os.write(creator);
		os.write(bin);
		os.write(current_loop);
		os.write(layer);
	}
};

class layern_creator_mapreduce : public base_mapreduce {
	layern_creator* creator;
	int layer;
	bool delete_layers_geq;
	bool add_reconstruction_edges_link;

	bool dispose_creator;
public:
	layern_creator_mapreduce() : creator(nullptr), dispose_creator(1), layer(0), delete_layers_geq(0), add_reconstruction_edges_link(0) {}
	virtual ~layern_creator_mapreduce() {
		if (dispose_creator && creator != nullptr) delete creator;
	}
	////////////////////////////////////////////////
	// main functionality is written in map and reduce
	virtual streamable* map(streamable* item);
	virtual streamable* reduce(list<streamable*> &item_list);

	virtual string map_get_key(streamable* item);

	// this mapreduce method does not have reduce
	virtual bool has_reduce() { return false; }
	// force mapreduce to wrap results with streamed_pointer* (to save on memory)
	virtual bool should_wrap_results() { return true; }

	////////////////////////////////////////////////
	// support functions
	void set_creator(layern_creator* c) { this->creator = c; dispose_creator = false; }
	void set_layer(int l) { this->layer = l; }
	void set_delete_layers_geq(bool d) { this->delete_layers_geq = d; }
	void set_add_reconstruction_edges_link(bool a) { this->add_reconstruction_edges_link = a; }

	// streamable implementations
	virtual streamable* make_instance() const { return new layern_creator_mapreduce(); }
	virtual void read_from_stream(istreamer& is) {
		streamable* ptr;
		is.read(ptr); creator = (layern_creator*)ptr;
		is.read(layer);
		is.read(delete_layers_geq);
		is.read(add_reconstruction_edges_link);
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write(creator);
		os.write(layer);
		os.write(delete_layers_geq);
		os.write(add_reconstruction_edges_link);
	}
};


#endif /* __OPTIMIZATION_H__ */
