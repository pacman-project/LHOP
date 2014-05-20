/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// feature_extraction

#pragma once
#ifndef _SHAPES_COLOR_FEATURE_EXTRACTION_
#define _SHAPES_COLOR_FEATURE_EXTRACTION_

#include "core/feature_extraction/abstract_feature_extraction.h"

/// @addtogroup module_shape Module - shapes
/// @{


// ColorEdgesFeatureExtraction
///////////////////////////////////////////////////////////////////////////////

class ColorEdgesFeatureExtraction : public AbstractFeatureExtraction  {
private:
    int gabor_size;                // = 7
    int n_rotations;               // = 6
	double gauss_sigma;			   // = 5
	double vmd_kappa;			   // = 5.0
	double gray_weight;            // = 3.0
	bool hard_resize;			   // true
	double resize_sigma;		   // 1.0

public:
	ColorEdgesFeatureExtraction(const IConfiguration& config, 
					const std::vector<IInputPreprocessing*>& preprocessing_steps, 
					const std::vector<IOutputPostprocessing*>postprocessing_steps) 
		: AbstractFeatureExtraction(config, preprocessing_steps,postprocessing_steps) {
	}


};

/// @}


#endif /* _SHAPES_COLOR_FEATURE_EXTRACTION_ */
