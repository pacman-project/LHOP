/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// layer_1_creators

#pragma once
#ifndef _LAYER_1_CREATORS_
#define _LAYER_1_CREATORS_

#include "../utils/img.h"
#include "../graphs/img_graph.h"
#include "../utils/matrix.h"
#include "../utils/misc.h"
#include "../utils/utils.h"
#include "layer_1_result.h"
#include "layers.h"

// layer1_creator_...
///////////////////////

//const unsigned NODE_REDUCED_ATTR = 1;
//const unsigned HAS_NEXT_LAYER = 4;
//const unsigned FROM_PREV_LAYER = 64;
//const unsigned IN_NODE_ATTR = 0x1000;
//const unsigned BOUDARY_NODE_ATTR = 0x2000;
//const unsigned OUT_NODE_ATTR = 0x4000;

class layer1_creator {
protected:
    //vector<img*> images;
    vector<img*> filter_kernels;
    vector<iimg*> regions;

    // image processing
    int scale_mask_size; 
    double scale_sigma;
	int scale_limit;
    int max_scales;
    int init_size;
    double scale_factor;

    //img* max_image;
    //iimg* max_image_src;

    //layer1_result* result;

public:
	vector<string> use_opencl_devices;
	bool use_opencl;
	bool opencl_verify_result;

	bool make_dummy_result;

    int to_neighbor;              
    double layer1_region_threshold;    
    double layer1_threshold;            
    double power_correction;
    double normalization_percent;
    bool normalize;
    double response_percent;
    int layer1_3x3bound;                
    int layer1_neighb_radius;         
    int border;

	void init_filters();
	void clear_filters();
    layer1_creator();
    ~layer1_creator();

    //virtual layer1_result* create_result(img*, img* = nullptr);
	//int create_result(vector<layer1_result*>& results, bimg*, img*);
	//int create_result(vector<layer1_result*>& results, img*, const img*);
	//int create_result(vector<layer1_result*>& results, const img&, const img&);
	virtual int create_result_vector(vector<layer1_result*>& results, const img&, const img&);
	virtual int create_result_vector(vector<layer1_result*>& results, const img&, const vector<irectangle2>&);

    virtual void cfg_init(const config_dictionary& cfg);
protected:

	virtual layer1_result* create_result(img& im, img& maskim, int originalw);

	int ocl_create_result(layer1_result* result, const img* im, const img* maskimg);

    void image_processing(img* im);

    bool init_default_from_result(layer1_result* res);

    void init_result_with_default(layer1_result* res);

    /// result = new ...
    virtual layer1_result* init_result();

    /// creates vector of masks -> "masks"
    virtual void make_filter_kernels();

    /// applies filter kernels to img -> "images"
    virtual void make_images(vector<img*>& images, img& im, img& immask);

	/// OpenCL will need to know how many images will be created (to reserve memory)
	virtual int make_images_get_count();

    /// create "maximal image" -> max_image
    /// create a matrix with indices specifying which image gave maximal value -> max_img_src
    void make_max_image(img*& max_image, iimg*& max_image_src, const vector<img*>& images);

    /// applies "3 x 3 filter" on max_image -> layer1_result
    void make_result(layer1_result* result, img* max_image, iimg* max_image_src, const img& maskimg);

    /// "dilutes" layer1_result -> layer1_result_inhib
    void inhibit_result(layer1_result* result);

    /// finalizes the result; see comments before the implementation
    virtual void finalize(layer1_result* result, vector<img*>& images, img* max_image, iimg* max_image_src);

    virtual unsigned lib_type();

    void delete_filter_kernels();
    //void delete_images();

	CStopWatch clock_1, clock_2, clock_3, clock_4, clock_5;
#ifdef OPENCL	
	
	void ocl_set_offsets_recursively(const best_cmd_queue_info& cmd_queue_info, cl_mem coords, int coords_size, int additional_offsets, int work_group_size_used, 
									cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt);
	void ocl_make_images(const best_cmd_queue_info& cmd_queue_info, layer1_result* result, const img* im, cl_mem made_images, cl_mem max_image, cl_mem max_image_index, cl_mem kernel_err,
							int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt);
	void ocl_make_sum2x2(const best_cmd_queue_info& cmd_queue_info, const img* im, cl_mem max_image, cl_mem sum_2x2, cl_mem kernel_err,
							int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt);
	void ocl_get_maximus(const best_cmd_queue_info& cmd_queue_info, const img* im, cl_mem input, cl_mem output_max, cl_mem kernel_err,
							int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt);

	void ocl_prepare_results(const best_cmd_queue_info& cmd_queue_info, layer1_result* result, const img* im, const img* maskimg,
							 cl_mem made_images, cl_mem max_image, cl_mem max_image_index, cl_mem sum_2x2, cl_mem maximum_max_image, cl_mem maximum_sum_2x2,
							 cl_mem layer1_coords_compact_sort_keys, cl_mem layer1_coords_compact, cl_mem number_created_result, cl_mem number_positions_used_result, cl_mem kernel_err,
							 int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt);

	void ocl_sort_compact_coords(const best_cmd_queue_info& cmd_queue_info, const img* im, const cl_int coords_size,
							 cl_mem layer1_coords_compact, cl_mem layer1_coords_compact_sort_keys, cl_mem kernel_err,
							 int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt);
	
	void ocl_update_results_offsets(const best_cmd_queue_info& cmd_queue_info, const img* im, const cl_int coords_size, cl_mem layer1_coords, cl_mem kernel_err,
									int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt,
									cl_event* profile_start, cl_event* profile_stop);

	void ocl_make_results(const best_cmd_queue_info& cmd_queue_info, layer1_result* result, const img* im, int number_positions_used,
							cl_mem made_images, cl_mem layer1_coords_compact, cl_mem layer1_coords_compact_sort_keys, cl_mem layer1_result_data, cl_mem kernel_err,
							int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt);
	
	void ocl_sort_results(const best_cmd_queue_info& cmd_queue_info, layer1_result* result, const img* im, cl_mem layer1_result_data, cl_int compact_count, cl_mem layer1_coords_compact, cl_mem merge_sort_temp_data, cl_mem kernel_err,
							int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt);
	
	void ocl_inhibit_results(const best_cmd_queue_info& cmd_queue_info, layer1_result* result, const img* im, 
								cl_mem layer1_result_data, cl_mem layer1_coords, cl_int adjusted_layer1_coords_size, cl_mem layer1_coords_inhib, cl_mem kernel_err,
								int memory_copy_type, cl_uint wait_event_count, cl_event* wait_evnt, cl_event* evnt);
#endif

public:
    /// Functions to display intermediate results, get specific data...

    virtual void get_part_data(vector<part_data*>&, const config_dictionary&) {  }
    virtual void get_masks(vector<img*>&) { }
    part_lib* get_part_lib(const config_dictionary& cfg);
};

class layer1_creator_struct : public layer1_creator {
public:
    double gabor_lambda;                // = 6.0
    double gabor_gamma;                 // = 0.75
    double gabor_bw;                    // = 1.75
    int gabor_step;                     // = 30

public:
    layer1_creator_struct();
    ~layer1_creator_struct();

    void init_from_result(layer1_result_struct* res);
protected:
    virtual layer1_result* init_result();
    virtual void make_filter_kernels();
    virtual void make_images(vector<img*>& images, img& im, img& immask);
	virtual int make_images_get_count();
    virtual unsigned lib_type() { return STRUCT_LIB_ATTR; }

public:
    static void get_filter_kernels(vector<img*>&, double, double, double, int);
    virtual void get_filter_kernels(vector<img*>&);
    virtual void get_part_data(vector<part_data*>&, const config_dictionary&);

    virtual void cfg_init(const config_dictionary& cfg);
};

class layer1_creator_color : public layer1_creator {
public:
    int gabor_size;                // = 7
    int n_rotations;               // = 6
	double gauss_sigma;			   // = 5
	double vmd_kappa;			   // = 5.0
	double gray_weight;            // = 3.0
	bool hard_resize;			   // true
	double resize_sigma;		   // 1.0

protected:
	vector<double> masks_norm;
	img* gaussFilter, * gaussFilterT;
	img* gaussDxFilter, * gaussDxFilterT;

public:
    layer1_creator_color();
	virtual ~layer1_creator_color();

	virtual int create_result_vector(vector<layer1_result*>& results, const img&, const img&);
    virtual int create_result_vector(vector<layer1_result*>& results, const img&, const vector<irectangle2>&);	
    void init_from_result(layer1_result_color* res);
protected:
    virtual layer1_result* create_result(img& im, img& maskim, int originalw);
    virtual layer1_result* init_result();
    virtual void make_filter_kernels();
	virtual void make_images(vector<img*>& images, img& im) { throw; }
	virtual int make_images_get_count() { throw; }
    virtual unsigned lib_type() { return STRUCT_LIB_ATTR; }

	double Omega(int ch, int col);
	int Omega_size();
	void SO(vector<img*>& M, img& R, img& G, img& B);
	void edges(img*& imgMag, img*& imgDir, img*& imgMax, vector<img*>& M, img& gray);
	layer1_result* get_result(img* imgMag, img* imgDir, img* imgMax);

public:
	static void get_masks(vector<img*>&, int, int);
	virtual void get_masks(vector<img*>&);
    virtual void get_part_data(vector<part_data*>&, const config_dictionary&);

	virtual void cfg_init(const config_dictionary& cfg);
};

class layer1_creator_app : public layer1_creator {
public:
    double gabor_lambda;                // = 10.0
    double gabor_gamma;                 // = 0.8
    double gabor_bw;                    // = 8.0

    int gabor_step;                     // = 45
public:
    layer1_creator_app();
    ~layer1_creator_app();

    void init_from_result(layer1_result_app* res);
protected:
    virtual layer1_result* init_result();
    virtual void make_filter_kernels();
    virtual void make_images(vector<img*>& images, img& im, img& immask);
	virtual int make_images_get_count();
    virtual unsigned lib_type() { return APP_LIB_ATTR; }

public:
    static void get_filter_kernels(vector<img*>&, double, double, double, int);
    virtual void get_part_data(vector<part_data*>&, const config_dictionary&);
    virtual void cfg_init(const config_dictionary& cfg);
};

const unsigned DOG_LIB_ATTR = 2;

class layer1_creator_dog : public layer1_creator {
public:
    double sigma_inner;                 // = sqrt(2)
    double sigma_outer;                 // = 2
    double mask_size_factor;            // = 2.8     
public:
    layer1_creator_dog();
    ~layer1_creator_dog();

    void init_from_result(layer1_result_dog* res);
protected:
    virtual layer1_result* init_result();
    virtual void make_filter_kernels();
    virtual void make_images(vector<img*>& images, img& im, img& immask);
	virtual int make_images_get_count();
    virtual unsigned lib_type() { return DOG_LIB_ATTR; }

public:
    static void get_filter_kernels(vector<img*>&, double, double, double);
    virtual void get_part_data(vector<part_data*>&, const config_dictionary&);
    virtual void cfg_init(const config_dictionary& cfg);
};

const unsigned LOG_LIB_ATTR = 3;


class layer1_creator_loggabor : public layer1_creator {
public:
    layer1_creator_loggabor();
    ~layer1_creator_loggabor();

    void init_from_result(layer1_result_loggabor* res);
protected:
    virtual layer1_result* init_result();
    virtual void make_filter_kernels();
    virtual void make_images(vector<img*>& images, img& im, img& immask);
	virtual int make_images_get_count();
    virtual unsigned lib_type() { return LOG_LIB_ATTR; }

public:
    static void get_filter_kernels(vector<img*>&, double dummy = 0); // dummy is here for linux compatibility

	virtual void get_part_data(vector<part_data*>&, const config_dictionary&);
    virtual void cfg_init(const config_dictionary& cfg);
};


// global stuff
/////////////////


#endif /* _LAYER_1_ */
