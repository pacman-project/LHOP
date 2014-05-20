/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// layer_1_creators

#pragma once
#ifndef _LAYER_1_CREATORS_
#define _LAYER_1_CREATORS_

#include "utils/img.h"
#include "utils/graphs/img_graph.h"
#include "utils/matrix.h"
#include "utils/misc.h"
#include "utils/utils.h"
#include "utils/config_dictionary.h"

#include "utils/class_register.h"

#include "core/legacy/layer_1_result.h"

// layer1_creator_...
///////////////////////

class layer1_creator : public IRegistrableClass { // feature extraction (not inference)
protected:
    vector<img*> filter_kernels;
    vector<iimg*> regions;

    // image processing
    int scale_mask_size; 
    double scale_sigma;
	int scale_limit;
    int max_scales;
    int init_size;
    double scale_factor;


public:
	bool use_opencl;

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

    virtual int create_result_vector(vector<layer1_result*>& results, const img&, const img&);
	virtual int create_result_vector(vector<layer1_result*>& results, const img&, const vector<irectangle2>&);

    virtual void cfg_init(const ConfigDictionary& cfg);


protected:

	virtual layer1_result* create_result(img& im, img& maskim, int originalw);

    void image_processing(img* im); // unused

    bool init_default_from_result(layer1_result* res);

    void init_result_with_default(layer1_result* res);

    /// result = new ...
    virtual layer1_result* init_result();

    /// creates vector of masks -> "masks"
    virtual void make_filter_kernels();

    /// applies filter kernels to img -> "images"
    virtual void make_images(vector<img*>& images, img& im, img& immask);

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

	
    virtual void get_part_data(vector<part_data*>&, const ConfigDictionary&) {  }
    virtual void get_masks(vector<img*>&) { }

public:
    /// Functions to display intermediate results, get specific data...

    part_lib* get_part_lib(const ConfigDictionary& cfg);

	/// Factory interface that must be implemented for any class that provides feature extraction functionality
	class IFactory : public IRegistrableClassFactory {
	 public:
		 virtual layer1_creator* newInstance() const = 0;
	};
};


// global stuff
/////////////////


#endif /* _LAYER_1_ */
