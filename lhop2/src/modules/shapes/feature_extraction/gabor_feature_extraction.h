/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// feature_extraction

#pragma once
#ifndef _SHAPES_GABOR_FEATURE_EXTRACTION_
#define _SHAPES_GABOR_FEATURE_EXTRACTION_

#include "core/feature_extraction/filterbank_feature_extraction.h"

static const char* const gabor_registration_names[] = { "gabor", "gabor.app", "gabor.dog", "gabor.loggabor"};

/// @addtogroup module_shape Module - shapes
/// @{

// General GaborFeatureExtraction (used by default)
///////////////////////////////////////////////////////////////////////////////

class GaborFeatureExtraction : public AbstractFilterbankFeatureExtraction  {
private:
	// used when generating filterbank
	double gabor_lambda;                // = 6.0
    double gabor_gamma;                 // = 0.75
    double gabor_bw;                    // = 1.75
    int gabor_step;                     // = 30

	// used when generating vocabulary part 
	double region_lambda;
	double region_gamma;
	double region_bw;

	double layer1_region_threshold;
public:
	GaborFeatureExtraction(const IConfiguration& config, 
					const std::vector<IInputPreprocessing*>& preprocessing_steps, 
					const std::vector<IOutputPostprocessing*>postprocessing_steps) 
		: AbstractFilterbankFeatureExtraction(config, preprocessing_steps, postprocessing_steps) {

		gabor_lambda = config.getDouble("gabor_lambda", 6.0);
		gabor_gamma = config.getDouble("gabor_gamma", 0.75);
		gabor_bw = config.getDouble("gabor_bw", 1.75);
		gabor_step = config.getInt("gabor_step", 30);

		region_lambda = config.getDouble("region_lambda", 7.0);
        region_gamma = config.getDouble("region_gamma", 1.0);
        region_bw = config.getDouble("region_bw", 2.75);

		layer1_region_threshold = config.getDouble("layer1_region_threshold", 0.1);
	}
	virtual std::vector<FeatureTypeDefinition> getFeatureTypeDefinitions();

	virtual std::vector<img*> generateFilterbank();
	virtual std::vector<img*> convolveImageWithFilterbank(const img* input_image, const std::vector<img*> filterbanks);

	class Factory : public IFactory {
	 public:
		 virtual AbstractFeatureExtraction* newInstance(const IConfiguration& config) const {			 
			 return new GaborFeatureExtraction(config, getPreprocessingSteps(config), getPostprocessingSteps(config));
		 }
		 virtual string assignedRegistrationName() const {  return gabor_registration_names[0]; }
		 virtual string getAssociatedInputObjectType() const { return IMAGE_INPUT_OBJECT; };
	};
};

// Different variations of gabor feature extraction (or gaussian):
///////////////////////////////////////////////////////////////////////////////

// AppGaborFeatureExtraction
///////////////////////////////////////////////////////////////////////////////
class AppGaborFeatureExtraction : public AbstractFilterbankFeatureExtraction  {
private:
	// used when generating filterbank
	double gabor_lambda;                // = 10.0
    double gabor_gamma;                 // = 0.8
    double gabor_bw;                    // = 8.0
    int gabor_step;                     // = 45

	// used when generating vocabulary part 
	double region_lambda;
	double region_gamma;
	double region_bw;

	double layer1_region_threshold;
public:
	AppGaborFeatureExtraction(const IConfiguration& config, 
					const std::vector<IInputPreprocessing*>& preprocessing_steps, 
					const std::vector<IOutputPostprocessing*>postprocessing_steps) 
		: AbstractFilterbankFeatureExtraction(config, preprocessing_steps,postprocessing_steps) {
		
		gabor_lambda = config.getDouble("gabor_lambda", 10.0);
		gabor_gamma = config.getDouble("gabor_gamma", 0.8);
		gabor_bw = config.getDouble("gabor_bw", 8.0);
		gabor_step = config.getInt("gabor_step", 45);

		config.getDouble("region_lambda", 7.0);
        config.getDouble("region_gamma", 1.0);
        config.getDouble("region_bw", 2.75);

		layer1_region_threshold = config.getDouble("layer1_region_threshold", 0.1);
	}
	virtual std::vector<FeatureTypeDefinition> getFeatureTypeDefinitions();

	virtual std::vector<img*> generateFilterbank();
	virtual std::vector<img*> generateFilterbank(double gabor_lambda, double gabor_gamma, double gabor_bw);
	virtual std::vector<img*> convolveImageWithFilterbank(const img* input_image, const std::vector<img*> filterbanks);

	class Factory : public IFactory {
	 public:
		 virtual AbstractFeatureExtraction* newInstance(const IConfiguration& config) const {
			 return new AppGaborFeatureExtraction(config, getPreprocessingSteps(config), getPostprocessingSteps(config));
		 }
		 virtual string assignedRegistrationName() const {  return gabor_registration_names[1]; }
		 virtual string getAssociatedInputObjectType() const { return IMAGE_INPUT_OBJECT; };
	};
};

// DogGaborFeatureExtraction
///////////////////////////////////////////////////////////////////////////////

class DogGaborFeatureExtraction : public AbstractFilterbankFeatureExtraction  {
public:
	// used when generating filterbank
    double sigma_inner;                 // = sqrt(2)
    double sigma_outer;                 // = 2
    double filter_kernel_size_factor;   // = 2.8     
	
	// used when generating vocabulary part 
	double region_sigma_inner;
	double region_sigma_outer;
	double region_mask_size_factor;

	double layer1_region_threshold;
public:
	DogGaborFeatureExtraction(const IConfiguration& config, 
					const std::vector<IInputPreprocessing*>& preprocessing_steps, 
					const std::vector<IOutputPostprocessing*>postprocessing_steps) 
		: AbstractFilterbankFeatureExtraction(config, preprocessing_steps,postprocessing_steps) {
		
		sigma_inner = config.getDouble("sigma_inner", sqrt(2.0));
		sigma_outer = config.getDouble("sigma_outer", 2);
		filter_kernel_size_factor = config.getDouble("mask_size_factor", 2.8);

		region_sigma_inner = config.getDouble("region_sigma_inner", 1.41);
        region_sigma_outer = config.getDouble("region_sigma_outer", 2.0);
        region_mask_size_factor = config.getDouble("region_mask_size_factor", 2.8);

		layer1_region_threshold = config.getDouble("layer1_region_threshold", 0.1);
	}

	virtual std::vector<FeatureTypeDefinition> getFeatureTypeDefinitions();

	virtual std::vector<img*> generateFilterbank();
	virtual std::vector<img*> generateFilterbank(double sigma_inner, double sigma_outer, double filter_kernel_size_factor);
	virtual std::vector<img*> convolveImageWithFilterbank(const img* input_image, const std::vector<img*> filterbanks);

	class Factory : public IFactory {
	 public:
		 virtual AbstractFeatureExtraction* newInstance(const IConfiguration& config) const {
			 return new DogGaborFeatureExtraction(config, getPreprocessingSteps(config), getPostprocessingSteps(config));
		 }
		 virtual string assignedRegistrationName() const {  return gabor_registration_names[2]; }
		 virtual string getAssociatedInputObjectType() const { return IMAGE_INPUT_OBJECT; };
	};
};

// LogGaborFeatureExtraction
///////////////////////////////////////////////////////////////////////////////

class LogGaborFeatureExtraction : public AbstractFilterbankFeatureExtraction  {
	double layer1_region_threshold;
public:
	LogGaborFeatureExtraction(const IConfiguration& config, 
					const std::vector<IInputPreprocessing*>& preprocessing_steps, 
					const std::vector<IOutputPostprocessing*>postprocessing_steps) 
		: AbstractFilterbankFeatureExtraction(config, preprocessing_steps,postprocessing_steps) {

		layer1_region_threshold = config.getDouble("layer1_region_threshold", 0.1);
	}
	virtual std::vector<FeatureTypeDefinition> getFeatureTypeDefinitions();

	virtual std::vector<img*> generateFilterbank();

	class Factory : public IFactory {
	 public:
		 virtual AbstractFeatureExtraction* newInstance(const IConfiguration& config) const {
			 return new LogGaborFeatureExtraction(config, getPreprocessingSteps(config), getPostprocessingSteps(config));
		 }
		 virtual string assignedRegistrationName() const {  return gabor_registration_names[3]; }
		 virtual string getAssociatedInputObjectType() const { return IMAGE_INPUT_OBJECT; };
	};
};


/// @}



#endif /* _SHAPES_GABOR_FEATURE_EXTRACTION_ */
