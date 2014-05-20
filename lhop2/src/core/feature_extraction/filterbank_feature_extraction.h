/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// filterbank_feature_extraction

#pragma once
#ifndef _CORE_FILTERBANK_FEATURE_EXTRACTION_
#define _CORE_FILTERBANK_FEATURE_EXTRACTION_

#include "core/feature_extraction/abstract_feature_extraction.h"

/// @addtogroup core
/// @{
/// @addtogroup feature_extraction
/// @{

/**
 * Base class for extracting features from an image using a filterbank. The descended class should provide filters 
 * by implementing std::vector<img*> generateFilterbank(); Each returned filter corresponds to one feature type used 
 * to create first layer parts.
 *
 * The process of extractFeatureActivations() is implemented as follows:
 *  - casts input into AbstractImageInputObject (requires input in form of an image)
 *  - obtains/loads filterbank using abstract obtainFilterbank() method
 *  - performs convolution on the image with all filters in filterbank
 *  - creates resulting std::vector<FeatureActivationArray> by mapping result of one filter to one FeatureActivationArray
 * 
 * Descended class must/can override following methods:
 *  - std::vector<img*> generateFilterbank(): required to be implemented by the descendend class
 *  - std::vector<img*> convolveImageWithFilterbank(const img* input_image, const std::vector<img*> filterbanks): 
 *       current implementation performs filtering with each filter in filterbank but can be overridden to
 *       implement more efficient convolution (e.g. result from one filter can be simply modified values of another)
 */
class AbstractFilterbankFeatureExtraction : public AbstractFeatureExtraction {
private:
	std::vector<img*> loaded_filterbank;
public:
	AbstractFilterbankFeatureExtraction(const IConfiguration& config, 
					const std::vector<IInputPreprocessing*>& preprocessing_steps, 
					const std::vector<IOutputPostprocessing*>postprocessing_steps) 
		: AbstractFeatureExtraction(config, preprocessing_steps, postprocessing_steps) {
	}
	virtual ~AbstractFilterbankFeatureExtraction() {
		for (int i = 0; i < loaded_filterbank.size(); ++i)
			if (loaded_filterbank[i] != nullptr)
				delete loaded_filterbank[i];
	}
	virtual std::vector<img*> generateFilterbank() = 0;
	virtual std::vector<img*> convolveImageWithFilterbank(const img* input_image, const std::vector<img*> filterbanks) {

		std::vector<img*> result;

		for (auto filter_iter = filterbanks.begin(); filter_iter != filterbanks.end(); ++filter_iter) {
			result.push_back(input_image->get_convolve(**filter_iter));
		}   

		return result;
	}

	std::vector<FeatureActivationArray> extractFeatureActivations(const AbstractInputObject* input) {
		// we require input object in a form of an image
		const AbstractImageInputObject* image_input_object = dynamic_cast<const AbstractImageInputObject*>(input);

		// get filterbanks (and load them if they are not loaded yet)
		std::vector<img*> filterbank = loaded_filterbank.size() == 0 ? generateFilterbank() : loaded_filterbank;
		
		// cache filterbank in the memory
		if (loaded_filterbank.size() == 0)
			loaded_filterbank = filterbank;

		// load image 
		std::shared_ptr<img> image_input = image_input_object->getImage();

		// convolve image with each filterbank
		std::vector<img*> convolved_images = convolveImageWithFilterbank(image_input.get(), filterbank);

		// transform array of convolution results into array of FeatureActivationArrays
		std::vector<FeatureActivationArray> result(convolved_images.size());

		VocabularyLayer& vocabulary_initial_layer = getInitialVocabulary()->getLayer(0);

		for (int i = 0; i < convolved_images.size(); i++) {
			VocabularyPart& vocabulary_part = vocabulary_initial_layer.getPartByTypeId(i);

			result[i] = FeatureActivationArray(&vocabulary_part, image_input->width, image_input->height);

			//for (auto iter = convolved_images[i]->begin(); iter != convolved_images[i]->end(); ++iter) {
			//	std::cout << *iter << std::endl;
			//}

			// copy convolution results to feature candidates
			auto convolution_iter = convolved_images[i]->begin();			// image iterator
			auto result_features_iter = result[i].beginFeatureIterator();	// our feature iterator			

			int ii = 0; 
			while (convolution_iter != convolved_images[i]->end() && result_features_iter != result[i].endFeatureIterator()) {
				// just copy feature 
				result_features_iter->value = *convolution_iter;

				//if (std::abs(result_features_iter->value) > 0.000001) {
				//	std::cout << result[i].beginLocationAt(ii).getLocationX(result[i]) << " " << result[i].beginLocationAt(ii).getLocationY(result[i]) << " " << result_features_iter->value << std::endl;
				//}

				// move to next position
				++convolution_iter; 
				++result_features_iter;

				ii++;
			}

			
		}
		
		return result;
	}
};

/// @}
/// @}

#endif /* _CORE_FILTERBANK_FEATURE_EXTRACTION_ */
