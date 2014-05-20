/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// feature_extraction

#include "gabor_feature_extraction.h"

#include "utils/class_register.h"

void register_gabor_feature_extraction_classes() {
	/************************************************************************* 
	 *           !!!! DO NOT MODIFY OR REMOVE THIS FUNCTION   !!!!!
	 ************************************************************************/

	ClassRegister::registerFactory<AbstractFeatureExtraction::IFactory, GaborFeatureExtraction::Factory>();	
	ClassRegister::registerFactory<AbstractFeatureExtraction::IFactory, AppGaborFeatureExtraction::Factory>();	
	ClassRegister::registerFactory<AbstractFeatureExtraction::IFactory, DogGaborFeatureExtraction::Factory>();	
	ClassRegister::registerFactory<AbstractFeatureExtraction::IFactory, LogGaborFeatureExtraction::Factory>();	
}

// GaborFeatureExtraction
///////////////////////////////////////////////////////////////////////////////

std::vector<FeatureTypeDefinition> GaborFeatureExtraction::getFeatureTypeDefinitions() {
	int number_filter_kernels = 180/gabor_step;
	
	std::vector<FeatureTypeDefinition> feature_types(number_filter_kernels);
	
	for (int angle = 0, i = 0; i < number_filter_kernels; angle += gabor_step, ++i) {
        
		// get original image filter from original settings and activation filter from region settings
		img* image_filter = img::gabor_mask(gabor_lambda, angle, 0.0, gabor_gamma, gabor_bw);
		img* activation_filter = img::gabor_mask(region_lambda, angle, 0.0, region_gamma, region_bw);

		// create feature type using image filter as reconstruction image and non-zero mask of activation filter as activation points
		feature_types[i].id = i;
		feature_types[i].activation_points = *activation_filter->get_bool_matrix(layer1_region_threshold);
		feature_types[i].reconstruction_image = *image_filter;

		delete image_filter;
		delete activation_filter;
    }

	return feature_types;
}

std::vector<img*> GaborFeatureExtraction::generateFilterbank() {
	
    int nfilter_kernels2 = 180/gabor_step;
	
	vector<img*> filter_kernels(2*nfilter_kernels2);
    
    for (int angle = 0, i = 0; i < nfilter_kernels2; angle += gabor_step, ++i) {
        filter_kernels[i] = img::gabor_mask(gabor_lambda, angle, 0.0, gabor_gamma, gabor_bw);
        filter_kernels[nfilter_kernels2 + i] = img::gabor_mask(gabor_lambda, angle, 90.0, gabor_gamma, gabor_bw);
    }

	return filter_kernels;
}
std::vector<img*> GaborFeatureExtraction::convolveImageWithFilterbank(const img* input_image, const std::vector<img*> filterbanks) {
	int nmasks2 = (int)filterbanks.size()/2;

    std::vector<img*> results(nmasks2);
	
	#pragma omp parallel for
    for (int i = 0; i < nmasks2; ++i) {        
		img* cosim = input_image->get_convolve(*filterbanks[i]);
        cosim->sqr();

		img* sinim = input_image->get_convolve(*filterbanks[nmasks2 + i]);
        sinim->sqr();
        *cosim += *sinim;
        cosim->sqrt();
        results[i] = cosim;

        delete sinim;
    }
	return results;
}


// AppGaborFeatureExtraction
///////////////////////////////////////////////////////////////////////////////

std::vector<FeatureTypeDefinition> AppGaborFeatureExtraction::getFeatureTypeDefinitions() {
	std::vector<img*> image_filters = generateFilterbank(gabor_lambda, gabor_gamma, gabor_bw);
	std::vector<img*> activation_filters = generateFilterbank(region_lambda, region_gamma, region_bw);
	
	std::vector<FeatureTypeDefinition> feature_types(image_filters.size());

	for (int i = 0; i < image_filters.size(); ++i) {
		// create feature type using image filter as reconstruction image and non-zero mask of activation filter as activation points
		feature_types[0].id = i;
		feature_types[0].activation_points = *activation_filters[i]->get_bool_matrix(layer1_region_threshold);
		feature_types[0].reconstruction_image = *image_filters[i];

		delete image_filters[i];
		delete activation_filters[i];
	}

	return feature_types;
}

std::vector<img*> AppGaborFeatureExtraction::generateFilterbank() {
	return generateFilterbank(gabor_lambda, gabor_gamma, gabor_bw);
}
std::vector<img*> AppGaborFeatureExtraction::generateFilterbank(double gabor_lambda, double gabor_gamma, double gabor_bw) {

	int nfilter_kernels = 360/gabor_step;

	vector<img*> filter_kernels(nfilter_kernels);
    for (int angle = 0, i = 0; i < nfilter_kernels; angle += gabor_step, ++i) {
        filter_kernels[i] = img::gabor_mask(gabor_lambda, angle, 90.0, gabor_gamma, gabor_bw);
    }

	return filter_kernels;
}
std::vector<img*> AppGaborFeatureExtraction::convolveImageWithFilterbank(const img* input_image, const std::vector<img*> filterbanks) {
	int nmasks2 = filterbanks.size()/2;
	
	std::vector<img*> results(filterbanks.size());
    for (int i = 0; i < nmasks2; ++i) {
        img* sinim = input_image->get_convolve(*filterbanks[i]);
        img* sinim2 = new img(*sinim);
        sinim->positive();
        sinim2->negative2();
        results[i] = sinim;
        results[nmasks2 + i] = sinim2;
    }

	return results;
}

// DogGaborFeatureExtraction
///////////////////////////////////////////////////////////////////////////////

std::vector<FeatureTypeDefinition> DogGaborFeatureExtraction::getFeatureTypeDefinitions() {
	std::vector<img*> image_filters = generateFilterbank(sigma_inner, sigma_outer, filter_kernel_size_factor);
	std::vector<img*> activation_filters = generateFilterbank(region_sigma_inner, region_sigma_outer, region_mask_size_factor);
	
	std::vector<FeatureTypeDefinition> feature_types(image_filters.size());

	for (int i = 0; i < image_filters.size(); ++i) {
		// create feature type using image filter as reconstruction image and non-zero mask of activation filter as activation points
		feature_types[0].id = i;
		feature_types[0].activation_points = *activation_filters[i]->get_bool_matrix(layer1_region_threshold);
		feature_types[0].reconstruction_image = *image_filters[i];

		delete image_filters[i];
		delete activation_filters[i];
	}

	return feature_types;
}

std::vector<img*> DogGaborFeatureExtraction::generateFilterbank() {
	return generateFilterbank(sigma_inner, sigma_outer, filter_kernel_size_factor);
}
std::vector<img*> DogGaborFeatureExtraction::generateFilterbank(double sigma_inner, double sigma_outer, double filter_kernel_size_factor) {
	int filter_kernel_w = 2*(int)(filter_kernel_size_factor*sigma_outer + 0.5) + 1;

    img* filter_kernel_in = img::gaussian_mask(filter_kernel_w, filter_kernel_w, sigma_inner);
    img* filter_kernel_out = img::gaussian_mask(filter_kernel_w, filter_kernel_w, sigma_outer);

    std::vector<img*> filter_kernels(2);

    *filter_kernel_in -= *filter_kernel_out;
    filter_kernels[0] = filter_kernel_in;
    filter_kernels[1] = new img(*filter_kernel_in);
    filter_kernels[1]->neg();

    delete filter_kernel_out;

	return filter_kernels;
}

std::vector<img*> DogGaborFeatureExtraction::convolveImageWithFilterbank(const img* input_image, const std::vector<img*> filterbanks) {

	std::vector<img*> results(2);

    results[0] = input_image->get_convolve(*filterbanks[0]);
    results[1] = new img(*results[0]);
    results[0]->positive();
    results[1]->negative2();

	return results;
}

// LogGaborFeatureExtraction
///////////////////////////////////////////////////////////////////////////////


std::vector<FeatureTypeDefinition> LogGaborFeatureExtraction::getFeatureTypeDefinitions() {
	// use the same filterbank as for image filtering/convolution
	std::vector<img*> filterbank = generateFilterbank();
	
	std::vector<FeatureTypeDefinition> feature_types(filterbank.size());
	for (int i = 0; i < filterbank.size(); ++i) {
		// create feature type using image filter as reconstruction image and its non-zero mask as activation points
		feature_types[i].id = i;
		feature_types[i].activation_points = *filterbank[i]->get_bool_matrix(layer1_region_threshold);
		feature_types[i].reconstruction_image = *filterbank[i];

		delete filterbank[i];
    }

	return feature_types;
}

std::vector<img*> LogGaborFeatureExtraction::generateFilterbank() {
	std::vector<img*> results;
	int nfilter_kernels = 6;

    results.resize(nfilter_kernels);
    for (int i = 0; i < 6; ++i) {
        results[i] = img::log_gabor_mask(5, i + 1);
    }

	return results;
}