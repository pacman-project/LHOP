// layer_learning learning algorithms:
//   * s_learning     learning of layers (statistics, g-responses, em-step, optimization)
//   * obj_learning   learning of object layer (models)
//   * t_learning     learning of thresholds

#pragma once

#ifndef _OPTIMIZATION_LEARNING_
#define _OPTIMIZATION_LEARNING_

#include "core/legacy/constants.h"
#include "core/legacy/layer_1_result.h"
#include "core/legacy/learning/layer_learning.h"

#include "utils/mapreduce.h"
#include "utils/serialization/streamed_pointer.h"


#include <tuple>
#include <stdio.h>

// optimization
///////////////////////////////////////////////////////////////////////////////


class layer1_optimization_new {
protected:
    part_lib* library;
    int layer;
	list<streamed_pointer> workingset;
    int parts_to_keep;
    set<int> final_parts;
    int loop;

    vector_parameter<double> cover_thresh;      // cfg: "covered_threshold", def. value: 0.85
    vector_parameter<double> int_thresh;        // cfg: "intersection_threshold", def. value: 0.2
    vector_parameter<int> bite_size;            // cfg: "bite_size", def. value: 10

    vector_parameter<int> optimization;                           // cfg: "optimize", 0: no optimization, 1: global optimization (default)
                                                //          2: global optimization of type 2
    string library_image_file;
    bool show_labels;

public:
    layer1_optimization_new(const ConfigDictionary& cfg, int player);

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
    void init(const ConfigDictionary& cfg);

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
		return ((layer1_result*)container[i].get()); 
    }
};



set<int> optimize_layer(part_lib* library, list<streamed_pointer>& workingset, int layer, const set<int>& init_set,
    double cover_thresh, double int_thresh, int bite_size);

set<int> optimize_layer_2(part_lib* library, list<streamed_pointer>& workingset, int layer, const set<int>& init_set,
    double cover_thresh, double int_thresh, int bite_size);

#endif /* _OPTIMIZATION_LEARNING_ */

