/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the HOPINTERFACE_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// HOPINTERFACE_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifndef __HOPINTERFACE
#define __HOPINTERFACE

#ifdef WIN32

#ifdef hop_EXPORTS
#define HOPINTERFACE_API __declspec(dllexport)
#else
#define HOPINTERFACE_API __declspec(dllimport)
#endif
#else
#define HOPINTERFACE_API 
#endif


#define HOP_FORMAT_BLOB     0
#define HOP_FORMAT_PROTOBUF 1
#define HOP_FORMAT_COMPRESS 16

#include <list>
#include "utils/mapreduce.h"
#include "utils/structures.h"
#include "utils/streaming.h"
#include "layers/layers.h"
#include "layers/hoc.h"

// classes
///////////////////////////////////////////////////////////////////////////////


struct HOPINTERFACE_API hop_exception {
    const char* str;

    hop_exception(const char* s) : str(s) { cout << s << endl; }
	virtual ~hop_exception() {
		if (str) delete str;
	}
};

class HOPINTERFACE_API hop_blob;

class HOPINTERFACE_API hop_struct {
protected:
    struct counted_ptr {
        void* ptr;
        int count;

        counted_ptr(void* p) : ptr(p), count(1) { }
    };

    counted_ptr* cptr;

public:
    hop_struct() : cptr(0) { }
    explicit hop_struct(void* data);
    //hop_struct(const hop_struct& hd) : cptr(0) { copy(hd); }
    virtual ~hop_struct() { }

    void* data() const { return cptr ? cptr->ptr : 0; }
    hop_struct& operator=(const hop_struct& hd) { copy(hd); return *this; }
    virtual hop_blob to_blob();
    virtual void to_file(const char* fname, int flags);
protected:

    void copy(const hop_struct& hd);
    virtual void dispose() = 0;
    virtual void write_to_stream(void* os); // pointer to ostreamer object
    virtual void read_from_stream(void* is); // pointer to istreamer object
};

// Arbitrary data.
// Note that memory is automatically deallocated using free!
class HOPINTERFACE_API hop_blob : public hop_struct {
protected:
    int size;

public:
    hop_blob(const hop_blob& b) : hop_struct(0), size(b.size) { copy(b); };

    // Data of 's' bytes
    hop_blob(void* data, int s);
    virtual ~hop_blob();

    int get_size() const { return size; }
    unsigned char* get_data_ptr() { return (unsigned char*)data(); }

	virtual void to_file(const char* fname, int flags) {
		if (data()) {
			ofstream out_f(fname, std::ofstream::binary);
			out_f.write((const char*)data(), get_size());
			out_f.close();
		}
	}

protected:
    virtual void dispose();
};

// Streamable implementation of hop_struct
class HOPINTERFACE_API hop_streamable : public hop_struct {
public:
	hop_streamable(const hop_streamable& stream_obj) : hop_struct(0) { copy(stream_obj); }
	explicit hop_streamable(void* stream_obj = 0) : hop_struct(stream_obj) {}
	
	hop_streamable(hop_blob& blob);
protected:
    virtual void dispose();
	virtual void write_to_stream(void* os);
};

// Interface to img class
class HOPINTERFACE_API hop_image : public hop_streamable {
public:
    hop_image(const hop_image& img) : hop_streamable(0) { copy(img); }

    // Construct object from img*
    hop_image(void* pimg);
    // Construct object from blob
    hop_image(hop_blob& blob);
    virtual ~hop_image();

    // Height of the image
    int get_width();

    // Width of the image
    int get_height();

    double* get_data_ptr();

    static hop_image from_bytes_rgb32(int width, int height, void* bytes);
protected:
    virtual void dispose();
    //virtual void write_to_stream(void* os);
};

// Interface to library of parts (class part_lib)
struct HOPINTERFACE_API hop_library_part {
    int layer;
    int type;
    // more information to come...

    hop_library_part(int l, int t) : layer(l), type(t) { }
};


class HOPINTERFACE_API hop_library : public hop_streamable {
public:
    static hop_library_part invalid_part;

public:
    hop_library(const hop_library& r) : hop_streamable(0) { copy(r); };

    // Construct object form part_lib*
    hop_library(void* plibrary);
    hop_library(hop_blob& blob);
    virtual ~hop_library();

    // Number of layers
    int get_layer_count();

    // Number of parts on layer 'l'. Returns -1 if layer l doean not exist.
    int get_part_count(int l);

    // Get part no. 'i' on layer 'l'. Returns invalid_part if 'l' or 'i' are 
    // out of range.
    hop_library_part get_part(int l, int i);

	void load_svm_models(char* filename);

	void print_layer_info(int layer);
protected:
    virtual void dispose();
    //virtual void write_to_stream(void* os);
};

struct HOPINTERFACE_API hop_node {
    int x, y;
    int layer;
    int type;
    double r_response;
	double g_response;
	double rr_response;
	double s_response;

    hop_node() : layer(-1), x(-1), y(-1), type(0), r_response(0.0), g_response(0.0), rr_response(0.0), s_response(0.0)  { }
    hop_node(int l, int px, int py, int t, double r, double g, double rr, double s) 
        : layer(l), x(px), y(py), type(t), r_response(r), g_response(g), rr_response(rr), s_response(s)  { }
    bool operator==(const hop_node& n) const 
        { return x == n.x && y == n.y && layer == n.layer && type == n.type; }
};

struct HOPINTERFACE_API hop_nodes {
    int size;
    hop_node* items;

    hop_nodes() : size(0), items(0) { }
};

struct HOPINTERFACE_API hop_indices {
    int size;
    int* items;

    hop_indices() : size(0), items(0) { }
};

class HOPINTERFACE_API hop_result : public hop_streamable {
public:
    static hop_node invalid_node;

public:
    hop_result(const hop_result& r) : hop_streamable(0) { copy(r); };

    // Constructs object from layer1_result*
    explicit hop_result(void* presult = 0);
    hop_result(hop_blob& blob);
    virtual ~hop_result();

	int get_original_image_width() const;
	int get_border() const;

    // Returns width of image at layer l
    int get_width(int l) const;

    // Returns height of image at layer l
    int get_height(int l) const;
    
    // Returns the number of layers
    int get_layer_count() const;

    // Returns node at particular position ('x', 'y', 'l')
    hop_node get_node_at(int l, int x, int y) const;

    // Returns the number of nodes at layer l
    int get_node_count(int l) const;

    // Returns the i-th node at layer l
    hop_node get_node(int l, int i) const;

    // Gets an array of "child nodes" of node n
    hop_nodes get_child_nodes(const hop_node& n) const;

    // Gets an array of indices of child nodes
    hop_indices get_child_nodes(int l, int i) const;

    // Saves the structure to a file -- compression is turned on for this structure!
    virtual void to_file(const char* fname, int flags);
protected:
    virtual void dispose();
    //virtual void write_to_stream(void* os);
};

class HOPINTERFACE_API hop_histogram_descriptor : public hoc::histogram_descriptor {
};


class HOPINTERFACE_API hop_histogram_generator {
	hoc::hoc_histogram_generator* cptr;

public:
    hop_histogram_generator() : cptr(0) { }
    virtual ~hop_histogram_generator() {
		delete cptr;
	}
	void initialize(const hop_library& lib, const char* params) {
		// use "hoc" namespace as root namespace
		config_dictionary cfg_root;	
		cfg_root.from_string(params);
		cfg_root.update_namespace_references();

		config_dictionary cfg;
		cfg.from_namespace_priority(cfg_root, 1, "hoc");

		cptr = new hoc::hoc_histogram_generator((part_lib*)lib.data());
		cptr->init_cfg(cfg);
	}

	void* data() { return cptr; }
};

class HOPINTERFACE_API hop_detection {
public:
	irectangle2 box;
	int category_id;
	vector<float> responses;
};

class HOPINTERFACE_API hop_streamed_pointer : public hop_streamable {
public:
	hop_streamed_pointer() : hop_streamable(nullptr) {}
	hop_streamed_pointer(const hop_streamed_pointer& sp) : hop_streamable(0) { copy(sp); }
	hop_streamed_pointer(void* presult) : hop_streamable(presult) { }	
	hop_streamed_pointer(hop_blob& blob) : hop_streamable(blob) { }

	// constructor that creates streamed_pointer from streamable data
	hop_streamed_pointer(const hop_streamable& s) : hop_streamable(new streamed_pointer((streamable*)s.data())) { }
	
	const char* get_filename() {
		streamed_pointer* ptr = (streamed_pointer*)data();
		string val = ptr->namep == nullptr ? nullptr : ptr->namep->name.c_str(); 

		char* val_c = new char[val.size()+1];
		memcpy(val_c, val.c_str(), val.size());
		val_c[val.size()] = 0;
		return val_c;
	}
	bool init_for_file(const char* name) {
		if (cptr == nullptr) {
			cptr = new counted_ptr(streamed_pointer::from_file(name));
			return true;
		} else {
			return false;
		}
	}
	virtual hop_streamable get() { 
		return hop_streamable(((streamed_pointer*)data())->get()); 
	}

	static hop_streamed_pointer from_file(const char* filename, bool is_file_disposable = true) {
		return hop_streamed_pointer(streamed_pointer::from_file(filename, is_file_disposable));
	}

};

class HOPINTERFACE_API hop_learning_results : public hop_streamable {
public:
	hop_learning_results(const hop_learning_results& r) : hop_streamable(0) { copy(r); };

	hop_learning_results(void* presult) : hop_streamable(presult) { }
	hop_learning_results(hop_blob& blob) : hop_streamable(blob) { }

	hop_library get_result() { return hop_library(((laylearning::lib_learning_results*)data())->lib); }

};

class HOPINTERFACE_API hop_learning_data_item : public hop_streamable {
public:
	hop_learning_data_item() : hop_streamable(new laylearning::lib_learning_controler::learning_data_item()) {}
	hop_learning_data_item(const hop_learning_data_item& r) : hop_streamable(0) { copy(r); };

	hop_learning_data_item(void* presult) : hop_streamable(presult) { }
	hop_learning_data_item(hop_blob& blob) : hop_streamable(blob) { }

	void add_input_image(const hop_streamable& obj) {
		laylearning::lib_learning_controler::learning_data_item* item = ((laylearning::lib_learning_controler::learning_data_item*)data());
		(*item)[laylearning::lib_learning_controler::INPUTDATA_IMAGE] = obj.data();
	}
	void add_groundtruth() {
		laylearning::lib_learning_controler::learning_data_item* item = ((laylearning::lib_learning_controler::learning_data_item*)data());
		(*item)[laylearning::lib_learning_controler::INPUTDATA_GROUNDTRUTH] = new list<irectangle2>();
	}
	void add_validation_pos_image(const hop_streamable& obj) {
		laylearning::lib_learning_controler::learning_data_item* item = ((laylearning::lib_learning_controler::learning_data_item*)data());
		(*item)[laylearning::lib_learning_controler::INPUTDATA_VALIDATION_POS_IMAGE] = obj.data();
	}
	void add_validation_neg_image(const hop_streamable& obj) {
		laylearning::lib_learning_controler::learning_data_item* item = ((laylearning::lib_learning_controler::learning_data_item*)data());
		(*item)[laylearning::lib_learning_controler::INPUTDATA_VALIDATION_NEG_IMAGE] = obj.data();
	}
};

class HOPINTERFACE_API hop_mapreduce_items : public hop_streamable {
public:
	hop_mapreduce_items(const hop_mapreduce_items& r) : hop_streamable(0) { copy(r); };

	hop_mapreduce_items(void* presult) : hop_streamable(presult) { }
	hop_mapreduce_items(hop_blob& blob) : hop_streamable(blob) { }

	virtual std::vector<hop_streamable> get_items() {
		list<streamable*>& items = ((mapreduce_items*)data())->items;

		std::vector<hop_streamable> hop_items;

		// reserve memory
		hop_items.reserve(items.size());

		// copy pointers to streamable
		for (auto iter = items.begin(); iter != items.end(); ++iter) {
			hop_items.push_back(hop_streamable(*iter));
		}

		return hop_items;
	}
};

class HOPINTERFACE_API hop_mapreduce_result : public hop_streamable {
public:
	hop_mapreduce_result(const hop_mapreduce_result& r) : hop_streamable(0) { copy(r); };

	hop_mapreduce_result(void* presult) : hop_streamable(presult) { }
	hop_mapreduce_result(hop_blob& blob) : hop_streamable(blob) { }

	void add(hop_streamable obj) { ((mapreduce_result*)data())->add((streamable*)obj.data()); }

	void set_result_wrapping(bool b) { ((mapreduce_result*)data())->wrap_results = b; }
	bool has_result_wrapping() { return ((mapreduce_result*)data())->wrap_results; }
};

class HOPINTERFACE_API hop_base_mapreduce : public hop_streamable {

public:
	hop_base_mapreduce(const hop_base_mapreduce& r) : hop_streamable(0) { copy(r); };

	hop_base_mapreduce(void* presult) : hop_streamable(presult) { }
	hop_base_mapreduce(hop_blob& blob) : hop_streamable(blob) { }

	virtual hop_streamable map(hop_streamable &item) {
		// get mapreduce implementation data
		base_mapreduce* mapreduce = (base_mapreduce*)this->data();

		// get streamable* item 
		streamable* s_item = (streamable*)item.data();

		// process item with map function
		streamable* result = mapreduce->map(s_item);

		// return result wrapped as streamable
		return hop_streamable((void*)result);
	}

	virtual hop_streamable reduce(std::vector<hop_streamable> &item_list) {
		// get mapreduce implementation data
		base_mapreduce* mapreduce = (base_mapreduce*)this->data();

		// compose list of hop_streamable into list of streamable*
		list<streamable*> s_item_list;

		for (auto iter = item_list.begin(); iter != item_list.end(); ++iter) {
			s_item_list.push_back((streamable*)iter->data());
		}

		// process item with reduce function
		streamable* result = mapreduce->reduce(s_item_list);
		
		// return result wrapped as streamable
		return hop_streamable((void*)result);
	}
	
	virtual const char* map_get_key(hop_streamable &item) { 
		string val = ((base_mapreduce*)data())->map_get_key((streamable*)(item.data())); 
		char* val_c = new char[val.size()+1];
		memcpy(val_c, val.c_str(), val.size());
		val_c[val.size()] = 0;
		return val_c;
	}

	virtual bool has_reduce() { return ((base_mapreduce*)data())->has_reduce(); }
	virtual bool should_wrap_results() { return ((base_mapreduce*)data())->should_wrap_results(); }

	virtual hop_mapreduce_result initilaize_results() { return hop_mapreduce_result(((base_mapreduce*)data())->initilaize_results()); }
};


class HOPINTERFACE_API hop_deployer_mapreudce : public base_deployer_mapreudce {
public:
	// This method will be called from c++ when request to submit job will be made.
	// Method should redirect call to Java overloaded submit method
	virtual mapreduce_result* submit(mapreduce_items* parallel_data, base_mapreduce* mapreduce_func) {

		hop_mapreduce_items items(parallel_data);
		hop_base_mapreduce func(mapreduce_func);

		hop_mapreduce_result result = submit(items, func);

		return (mapreduce_result*)result.data();
	}

	virtual hop_mapreduce_result submit(hop_mapreduce_items& parallel_data, hop_base_mapreduce& mapreduce_func) {
		//throw new_libhop_exception("Must be implemented by asscending class");

		// in case this is not overloaded by java call original base method (a lot of overhead due to wrappers but this should not be used often anyway)
		return hop_mapreduce_result(base_deployer_mapreudce::submit((mapreduce_items*)parallel_data.data(),(base_mapreduce*)mapreduce_func.data()));
	}
};

class HOPINTERFACE_API hop_learning_controler : public laylearning::lib_learning_controler {
public:

	hop_learning_controler()  {
	
	}
	virtual ~hop_learning_controler() {

	}

	void initialize(const char* params) {
        config_dictionary cfg_all;
		cfg_all.from_string(params);
		cfg_all.update_namespace_references();

		config_dictionary cfg_learn;
		cfg_learn.from_namespace_priority(cfg_all, 1, "learning");

		lib_learning_controler::initialize(cfg_learn);
	}
	void cleanup() {
		lib_learning_controler::cleanup();
	}

	void add_learning_data(const hop_learning_data_item& data_item) {
		learning_data_item* data_item_ptr = (learning_data_item*)data_item.data();
		lib_learning_controler::add_learning_data(*data_item_ptr);
	}
	  
	void set_starting_library(hop_library& lib) { 
		lib_learning_controler::set_starting_library((part_lib*)lib.data());
	}
	
	const char* get_action() { 
		return lib_learning_controler::get_action().c_str(); 
	}

	virtual hop_learning_results perform_learning() {
		printf("number of max openmp threads: %d\n", omp_get_max_threads( ));
		return hop_learning_results(lib_learning_controler::perform_learning());
	}

	virtual void set_mapreduce_implementation(hop_deployer_mapreudce* deployer) {
		lib_learning_controler::set_mapreduce_implementation(deployer);
	}

	void use_openmp_mapreduce_implementation() {
		lib_learning_controler::set_mapreduce_implementation(new openmp_deployer_mapreudce(), true);
	}
};

// functions
///////////////////////////////////////////////////////////////////////////////

// Inference from image (1st layer). 'result' is an array of results,
// each corresponding to one scale. Return value is size of the array.
HOPINTERFACE_API int hop_inference(hop_result*& result, hop_image& image, const char* params);
HOPINTERFACE_API int hop_inference(hop_result*& result, hop_image& image, const std::list<irectangle2> mask_regions, const char* params);

// Inference on 'result' with library 'lib' (higher layers)
HOPINTERFACE_API void hop_inference(hop_result& result, hop_library& lib, const char* params);


#ifdef __cplusplus    // If used by C++ code, 
extern "C" {          // we need to export the C interface
#endif

// Redirect stdout and stderr
HOPINTERFACE_API void hop_redirect_stdout(const char* filename);
HOPINTERFACE_API void hop_redirect_stderr(const char* filename);

// Enables unbufferd stdout to immediately write to stdout (instead of waiting for end of line)
HOPINTERFACE_API void hop_force_unbuffer_stdout();

// Reads image from file
HOPINTERFACE_API hop_image hop_read_image(const char* fname);

// Grayscale image from array of doubles with values between 0.0 and 1.0!
HOPINTERFACE_API hop_image hop_image_from_array(int width, int height, double* arr);

// Saves image to file
HOPINTERFACE_API void hop_save_image(hop_image& img, const char* fname);

// Reads library from file
HOPINTERFACE_API hop_library hop_read_library(const char* fname);

// Reads result from file
HOPINTERFACE_API hop_result hop_read_result(const char* fname);

// Inference on directory defined in cfgfile or params
HOPINTERFACE_API bool hop_1_inference(const char* cfgfile, const char* pattern, const char* params);

// Inference on directory defined in cfgfile or params
HOPINTERFACE_API bool hop_n_inference(const char* cfgfile, const char* pattern, const char* params);

// Leraning on directory defined in cfgfile or params
HOPINTERFACE_API bool hop_learning(const char* cfgfile, const char* pattern, const char* params);

// Leraning on directory defined in cfgfile or params
HOPINTERFACE_API bool hop_learning_tools(const char* cfgfile, const char* pattern, const char* params);

// Display result (depends on parameters)
HOPINTERFACE_API hop_image hop_display(hop_result& result, const char* params);

HOPINTERFACE_API hop_result hop_merge_scales(std::vector<hop_result>& result, const char* params);

// Create HoC (histogram of compositions) descriptors
HOPINTERFACE_API int hop_create_histograms(std::vector<hop_histogram_descriptor>*& result, hop_result& res_wraper, std::list<irectangle2> export_windows, hop_histogram_generator& hoc_generator);

// Obtain bounding box list of all detections at highest processed layer
//HOPINTERFACE_API int hop_get_detections(std::list<irectangle2>*& result, hop_result& res_wraper, hop_library& lib_wraper);
HOPINTERFACE_API int hop_get_detections(std::list<hop_detection>*& result, hop_result& res_wraper, hop_library& lib_wraper, const int only_category = -1);

// Learning (for distributed processing)
/*HOPINTERFACE_API hop_mlearner hop_generate_mlearner(hop_library& lib, const char* params);

HOPINTERFACE_API hop_mlearner hop_mlearner_update(hop_result& image_layers, hop_mlearner& mlearner);
HOPINTERFACE_API hop_mlearner hop_mlearner_merge(hop_mlearner& updated_learner, hop_mlearner& merged_learner);

HOPINTERFACE_API hop_plearner hop_generate_plearner(hop_mlearner& merged_mlearner, const char* params);

HOPINTERFACE_API hop_plearner hop_plearner_update(hop_result& image_layers, hop_plearner& plearner);
HOPINTERFACE_API hop_plearner hop_plearner_merge(hop_plearner& updated_learner, hop_plearner& merged_learner);

HOPINTERFACE_API hop_library hop_generate_library(hop_plearner& merged_plearner, const char* params);
*/
// Retrunes timestamp of build as string
HOPINTERFACE_API const char* hop_time_stamp();

// Retrunes timestamp of build as string
HOPINTERFACE_API void hop_enable_opencl(const char* c_kernel_paths, const char* c_used_devices);

#ifdef __cplusplus    // If used by C++ code, 
}
#endif


// This class is exported from the hop-interface.dll
//class HOPINTERFACE_API Chopinterface {
//public:
//	Chopinterface(void);
//	// TODO: add your methods here.
//}; 
//
//extern HOPINTERFACE_API int nhopinterface;
//
//HOPINTERFACE_API int fnhopinterface(void);

#endif

