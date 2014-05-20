/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// feature_extraction

#include "abstract_feature_extraction.h"

#include "utils/class_register.h"

void register_feature_extraction_classes() {
	/************************************************************************* 
	 *           !!!! DO NOT MODIFY OR REMOVE THIS FUNCTION   !!!!!
	 ************************************************************************/
}

/**
 * Simple pre-processing of the input object. Object is pre-processed with the IInputPreprocessing implementation in
 * AbstractFeatureExtraction::preprocessing_steps list in sequention order (order of steps in list is important). Each pre-processing
 * step can produce one or more new input objects. Result of one pre-processing is used as input in the next step.
 *
 * The input type must match with the input type required by the IInputPreprocessing implementation.
 *
 * Memory of intermediate results (between pre-process steps) is handled by this function but inpunt and the final results MUST
 * be handled by the caller.
 *
 * @param input_object: initial AbstractInputObject* object
 * @returns vector of AbstractInputObject* that have been fully preprocessed by all IInputPreprocessing in the AbstractFeatureExtraction::preprocessing_steps list
 */ 
std::vector<AbstractInputObject*> AbstractFeatureExtraction::performPreprocessing(const AbstractInputObject* input_object) {

	// start with initial input object
	std::vector<AbstractInputObject*> input_object_list;
	input_object_list.push_back((AbstractInputObject*)input_object);

	// for each pre-processing step transform preprocessed_objects input_object_list -> preprocessed_object_list
	for (auto step_iter = preprocessing_steps.begin(); step_iter != preprocessing_steps.end(); ++step_iter) {
		IInputPreprocessing* preprocess_step = (*step_iter);

		// do preprocessing
		std::vector<AbstractInputObject*> preprocessed_object_list = preprocess_step->doPreprocessing(input_object_list);

		// first we handle objects that need to be deleted from memory
		// we should delete any objects that are in input_object_list but are not any more in the preprocessed_object_list
		// unless this is initial input_object which should be handled by the caller 
		for (auto initial_iter = input_object_list.begin(); initial_iter != input_object_list.end(); ++initial_iter) {
			AbstractInputObject* initial_obj = *initial_iter;
			
			// skip if it is the input_object
			if (initial_obj == input_object)
				continue;

			// if we can find the same object in the preprocessed_object_list we continue to next object
			bool exists_in_preprocessed_list = false;
			for (auto preprocessed_list = preprocessed_object_list.begin(); preprocessed_list != preprocessed_object_list.end(); ++preprocessed_list) {
				AbstractInputObject* preprocessed_obj = *preprocessed_list;

				if (initial_obj == preprocessed_obj) {
					exists_in_preprocessed_list = true;
					break;
				}
			}

			// otherwise we need to delete this object
			if (exists_in_preprocessed_list == false) {
				delete initial_obj;
			}
		}

		// use preprocessed_object_list for the next step
		input_object_list = preprocessed_object_list;
	}

	// result will be stored in the input_object_list
	return input_object_list;
}

/**
 * Simple post-processing of the output objects. Object is post-processed with the IOutputPostprocessing implementation in
 * AbstractFeatureExtraction::postprocessing_steps list in sequention order (order of steps in list is important). Each post-processing
 * step can produce one or more new output objects. Result of one post-processing is used as input in the next step.
 *
 * The input type must match with the input type required by the IOutputPostprocessing implementation.
 *
 * Memory of intermediate results (between post-process steps) is handled by this function but final results MUST
 * be handled by the caller.
 *
 * @param output_object: initial AbstractInputObject* object
 * @returns vector of AbstractOutputObject* that have been fully postprocessed by all IOutputPostprocessing in the AbstractFeatureExtraction::postprocessing_steps list
 */ 
std::vector<AbstractOutputObject*> AbstractFeatureExtraction::performPostprocessing(const std::vector<AbstractOutputObject*>& output_object) {
	
	std::vector<AbstractOutputObject*> output_object_list = output_object;

	// for each post-processing step transform preprocessed_objects input_object_list -> preprocessed_object_list
	for (auto step_iter = postprocessing_steps.begin(); step_iter != postprocessing_steps.end(); ++step_iter) {
		IOutputPostprocessing* postprocess_step = (*step_iter);

		// do postprocessing
		std::vector<AbstractOutputObject*> postprocessed_object_list = postprocess_step->doPostprocessing(output_object_list);

		// first we handle objects that need to be deleted from memory
		// we should delete any objects that are in output_object_list but are not any more in the postprocessed_object_list
		for (auto initial_iter = output_object_list.begin(); initial_iter != output_object_list.end(); ++initial_iter) {
			AbstractOutputObject* initial_obj = *initial_iter;
			
			// if we can find the same object in the postprocessed_object_list we continue to next object
			bool exists_in_postprocessed_list = false;
			for (auto postprocessed_list = postprocessed_object_list.begin(); postprocessed_list != postprocessed_object_list.end(); ++postprocessed_list) {
				AbstractOutputObject* postprocessed_obj = *postprocessed_list;

				if (initial_obj == postprocessed_obj) {
					exists_in_postprocessed_list = true;
					break;
				}
			}

			// otherwise we need to delete this object
			if (exists_in_postprocessed_list == false) {
				delete initial_obj;
			}
		}

		// use postprocessed_object_list for the next step
		output_object_list = postprocessed_object_list;
	}

	// result will be stored in the output_object_list
	return output_object_list;
}

FeatureActivationArray AbstractFeatureExtraction::findMaxActivations(std::vector<FeatureActivationArray>& feature_activations, float* global_max_activation_return) {
	if (feature_activations.size() <= 0) 
		throw new_libhop_exception("Empty input received in the AbstractFeatureExtraction::findMaxActivations");
	
	// prepare output array of the same size
	FeatureActivationArray max_activations(nullptr, feature_activations[0].getWidth(), feature_activations[0].getHeight());

	float global_max_activation = -FLT_MAX;
	// iterate through all locations of activations (this goes over x and y) // TODO: check if this is thread safe and add OpenMP 
	for (auto location = max_activations.beginLocation(); location != max_activations.endLocation(); ++location) {
		float max_val = -FLT_MAX;
		// at each location get max value over all types of features
		for (auto iter = feature_activations.begin(); iter != feature_activations.end(); ++iter) { 
			max_val = std::max<float>(max_val, iter->getCandidateAt(location).value);
		}

		// save max value into max_activations for current location
		max_activations.getCandidateAt(location).value = max_val;

		// also get global max value
		global_max_activation = std::max<float>(global_max_activation, max_val);
	}

	if (global_max_activation_return != nullptr)
		*global_max_activation_return = global_max_activation;

	return max_activations;
}

/**
 * | 0 1 ...|                    | x   x   ...| 
 * | 2 3 ...| -> sum[0,1,2,3] -> | x [sum] ...|
 * | ...    |                    | ...        |
 */
FeatureActivationArray AbstractFeatureExtraction::findSum2x2MaxActivations(FeatureActivationArray& max_activations, float* global_max_activation_sum2x2_return) {
    
	// copy max data to max_sum2x2 data
	FeatureActivationArray max_feature_activations_sum2x2(max_activations);
	
	// create 4 iterators representing a 2x2 window
	// | 0 1 ...| 
	// | 2 3 ...|  
	// | ...    |
	// we will be moving this window over the array of max_activations
	std::vector<FeatureActivationResponse>::iterator max_activations_iter[9];

	max_activations_iter[0] = max_activations.beginFeatureIteratorAt(0, 0); 
    max_activations_iter[1] = max_activations.beginFeatureIteratorAt(1, 0); 
    max_activations_iter[2] = max_activations.beginFeatureIteratorAt(0, 1); 
    max_activations_iter[3] = max_activations.beginFeatureIteratorAt(1, 1);

	// create interator for saving result of sum 2x2
	// | 0  1  ...| 
	// | 2 [3] ...| -> write sum of 2x2 into lower right corner
	// | ...      |
	std::vector<FeatureActivationResponse>::iterator max_activations_sum2x2_iter = max_feature_activations_sum2x2.beginFeatureIteratorAt(1, 1);
	
	int width = max_activations.getWidth();
	int height = max_activations.getHeight();

    float max_sum = 0.0;
	for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
			// get sum of 2x2 window
			float sum = 0.0;
			for (int k = 0; k < 4; ++k) 
				sum += max_activations_iter[k]->value;

            // save to result position
			max_activations_sum2x2_iter->value = sum;

			// also get global max value
			max_sum = std::max<float>(max_sum, sum);

			// move the 2x2 window iterator
			for (int k = 0; k < 4; ++k)
				++max_activations_iter[k];

			// increment the result iterator
			++max_activations_sum2x2_iter;
		}
		// skip incrementing at the last loop since we will get out-of-bounds exceptions
		if (y < height -2) {
			// move the 2x2 windows iterator (by two locations to avoid edges) 
			for (int k = 0; k < 4; ++k)
				max_activations_iter[k]+=2;

			// increment the result iterator
			++max_activations_sum2x2_iter;
		}
	}

	if (global_max_activation_sum2x2_return != nullptr)
		*global_max_activation_sum2x2_return = max_sum;

    return max_feature_activations_sum2x2;       
}

AbstractFeatureExtraction::IntermediateData& AbstractFeatureExtraction::applyGlobalTresholding(IntermediateData& data, float layer1_threshold) {
	if (data.feature_activations.size() <= 0) 
		throw new_libhop_exception("Empty input received in the AbstractFeatureExtraction::applyGlobalTresholding");

	// calculate threshold by multiplying user threshold with the global max activation value
	double threshold = layer1_threshold * data.global_max_activations;
	
	FeatureActivationArray& max_activations = data.max_feature_activations;

	// iterate through all locations of max activations (this goes over x and y) // TODO: check if this is thread safe and add OpenMP 	
	// and shutdown activations (set to invalid) that are below threshold
	for (auto location = max_activations.beginLocation(); location != max_activations.endLocation(); ++location) {
		
		if (max_activations.getCandidateAt(location).value < threshold)			
			max_activations.setCandidateAtValid(location, false);		
	}

	// next we need to disable the corresponding activations of other feature types (not just max)
	// we can directly copy the is_valid memory block array to all the other feature types
	for (auto iter = data.feature_activations.begin(); iter != data.feature_activations.end(); ++iter)
		iter->activations_array.is_valid = max_activations.activations_array.is_valid;	

	return data;
}

AbstractFeatureExtraction::IntermediateData& AbstractFeatureExtraction::applyLocalInhibitionOfMax(IntermediateData& data, float layer1_3x3bound){
	if (data.feature_activations.size() <= 0) 
		throw new_libhop_exception("Empty input received in the AbstractFeatureExtraction::applyLocalInhibitionOfMax");
	
	FeatureActivationArray& max_activations = data.max_feature_activations;
	FeatureActivationArray& max_feature_activations_sum2x2 = data.max_feature_activations_sum2x2;

	// create 9 iterators representing a 3x3 region/window; a central one and one for each neighboor:
	// | 0  1  2 ...| 
	// | 3 [4] 5 ...|  
	// | 6  7  8 ...| 
	// | ...		|
	// we will be moving this window over the array of max_activations and max_activations_sum2x2
	std::vector<FeatureActivationResponse>::iterator max_activations_iter[9];
	std::vector<FeatureActivationResponse>::iterator max_activations_sum2x2_iter[9];

	// first row of window ([0 1 2])
	max_activations_iter[0] = max_activations.beginFeatureIteratorAt(0, 0); max_activations_sum2x2_iter[0] = max_feature_activations_sum2x2.beginFeatureIteratorAt(0, 0); 
    max_activations_iter[1] = max_activations.beginFeatureIteratorAt(1, 0); max_activations_sum2x2_iter[1] = max_feature_activations_sum2x2.beginFeatureIteratorAt(1, 0); 
    max_activations_iter[2] = max_activations.beginFeatureIteratorAt(2, 0); max_activations_sum2x2_iter[2] = max_feature_activations_sum2x2.beginFeatureIteratorAt(2, 0); 
	// second row of window ([3 4 5])
    max_activations_iter[3] = max_activations.beginFeatureIteratorAt(0, 1); max_activations_sum2x2_iter[3] = max_feature_activations_sum2x2.beginFeatureIteratorAt(0, 1);
    max_activations_iter[4] = max_activations.beginFeatureIteratorAt(1, 1); max_activations_sum2x2_iter[4] = max_feature_activations_sum2x2.beginFeatureIteratorAt(1, 1); // this is the middle point 
    max_activations_iter[5] = max_activations.beginFeatureIteratorAt(2, 1); max_activations_sum2x2_iter[5] = max_feature_activations_sum2x2.beginFeatureIteratorAt(2, 1);
	// last row of window ([6 7 8])
    max_activations_iter[6] = max_activations.beginFeatureIteratorAt(0, 2); max_activations_sum2x2_iter[6] = max_feature_activations_sum2x2.beginFeatureIteratorAt(0, 2);
    max_activations_iter[7] = max_activations.beginFeatureIteratorAt(1, 2); max_activations_sum2x2_iter[7] = max_feature_activations_sum2x2.beginFeatureIteratorAt(1, 2);
    max_activations_iter[8] = max_activations.beginFeatureIteratorAt(2, 2); max_activations_sum2x2_iter[8] = max_feature_activations_sum2x2.beginFeatureIteratorAt(2, 2);

	std::vector<int>::iterator max_is_valid_iter = max_activations.beginIsValidIteratorAt(1, 1); // is valid iterator for middle point only

	int width = max_activations.getWidth();
	int height = max_activations.getHeight();

	// iterate through all locations of max activations but skip edge locations
	for (int y = 1; y < height - 1; ++y) {
        for (int x = 1; x < width - 1; ++x) {
			// process only ones that are still active
			int is_valid = *max_is_valid_iter ;
			if (*max_is_valid_iter == true) {
			
				float mid = max_activations_iter[4]->value;
				int gcount = 0;

				// check how many neighbors are stronger than me				
				for (int k = 0; k < 9; ++k) {
					if (max_activations_iter[k]->value >= mid) 
						++gcount; 
				}
                if (gcount > layer1_3x3bound) {
					int gcount2 = 0;

					// if too many then also check how many neighbors are stronger then me in 2x2 sum
                    for (int k = 0; k < 9; ++k) {
						if (max_activations_sum2x2_iter[k]->value >= max_activations_sum2x2_iter[4]->value) 
							++gcount2; 
					}
                    if (gcount2 > layer1_3x3bound) {
						// if I am still not stronger then
						// shutdown activation as it gets inhibited by the neihbors
						*max_is_valid_iter = false;
                    }
                }
				
			}
			// move the window by one
			for (int k = 0; k < 9; ++k) { 
				++max_activations_iter[k]; 
				++max_activations_sum2x2_iter[k];
			}
			// incerement is_valid iterator
			++max_is_valid_iter;
		}

		// skip incrementing at the last loop since we will get out-of-bounds exceptions
		if (y < height -2) {
			// move the windows (by two locations to avoid edges)
			for (int k = 0; k < 9; ++k) { 
				max_activations_iter[k]+=2; 
				max_activations_sum2x2_iter[k]+=2;
			}
			// incerement is_valid iterator (also twice to avoid edges)
			max_is_valid_iter+=2;
		}
	}

	// next we need to disable the corresponding activations of other feature types (not just max)
	// we can directly copy the is_valid memory block array to all the other feature types
	for (auto iter = data.feature_activations.begin(); iter != data.feature_activations.end(); ++iter)
		iter->activations_array.is_valid = max_activations.activations_array.is_valid;

	return data;
}
AbstractFeatureExtraction::IntermediateData& AbstractFeatureExtraction::applyGlobalNormalization(IntermediateData& data) {
	if (data.feature_activations.size() <= 0) 
		throw new_libhop_exception("Empty input received in the AbstractFeatureExtraction::applyGlobalNormalization");

	// normalize with max activation value over all features and over all locations (replace with multiplication for efficency reasons)
	float global_normalization_value = 1/data.global_max_activations;

	// iterate over all feature types
	for (auto feature_type_iter = data.feature_activations.begin(); feature_type_iter != data.feature_activations.end(); ++feature_type_iter) {

		FeatureActivationArray& feature_type_activations = *feature_type_iter;

		// iterate over all location of this feature type
		for (auto location = feature_type_activations.beginLocation(); location != feature_type_activations.endLocation(); ++location) {
			// apply normalization only to valid loations (should be faster)
			if (feature_type_activations.isCandiateAtValid(location) == true) {				
				// activation_value[x,y] = activation_value[x,y] / global max activations value
				feature_type_activations.getCandidateAt(location).value *= global_normalization_value;
			}
		}
	}
	return data;
}

FeatureActivationArray& applyPowerCorrectionSingleFeature(FeatureActivationArray& feature_type_activations, float power_correction) {
	// iterate over all locations of this feature type
	for (auto location = feature_type_activations.beginLocation(); location != feature_type_activations.endLocation(); ++location) {
		// apply normalization only to valid loations (should be faster)
		if (feature_type_activations.isCandiateAtValid(location) == true) {				
			// activation_value[x,y] = activation_value[x,y] ^ power_correction 
			feature_type_activations.getCandidateAt(location).value = pow(feature_type_activations.getCandidateAt(location).value, power_correction);
		}
	}

	return feature_type_activations;
}
AbstractFeatureExtraction::IntermediateData& AbstractFeatureExtraction::applyPowerCorrection(IntermediateData& data, float power_correction){
	if (data.feature_activations.size() <= 0) 
		throw new_libhop_exception("Empty input received in the AbstractFeatureExtraction::applyGlobalNormalization");

	if (power_correction == 1)
		return data;

	// iterate over all feature types
	for (auto feature_type_iter = data.feature_activations.begin(); feature_type_iter != data.feature_activations.end(); ++feature_type_iter) {
		// and apply power correction to each feature type
		applyPowerCorrectionSingleFeature(*feature_type_iter, power_correction);
	}

	// also apply power correction to maximal feature (because it is used later on)
	// TODO: why does this not work ??
	//applyPowerCorrectionSingleFeature(data.max_feature_activations, power_correction);
	//applyPowerCorrectionSingleFeature(data.max_feature_activations_sum2x2, power_correction);

	//data.global_max_activations = pow(data.global_max_activations, power_correction);
	//data.global_max_activations_sum2x2 = pow(data.global_max_activations_sum2x2, power_correction);

	return data;
}

AbstractFeatureExtraction::IntermediateData& AbstractFeatureExtraction::applyLocationFiltering(IntermediateData& data, float response_percent){
	if (data.feature_activations.size() <= 0) 
		throw new_libhop_exception("Empty input received in the AbstractFeatureExtraction::applyGlobalNormalization");

	FeatureActivationArray& max_activations = data.max_feature_activations;
	
	// iterate over all feature types
	for (auto feature_type_iter = data.feature_activations.begin(); feature_type_iter != data.feature_activations.end(); ++feature_type_iter) {

		FeatureActivationArray& feature_type_activations = *feature_type_iter;

		// iterate over all locations of this feature type
		for (auto location = feature_type_activations.beginLocation(); location != feature_type_activations.endLocation(); ++location) {
			// apply normalization only to valid loations (should be faster)
			if (feature_type_activations.isCandiateAtValid(location) == true) {

				float max_activation_value = max_activations.getCandidateAt(location).value;
				float activation_value = feature_type_activations.getCandidateAt(location).value;
				
				// if max_activations[x,y] * response_percent > feature_type_activations[x,y] then shutdown activation
				if (max_activation_value * response_percent > activation_value)
					feature_type_activations.setCandidateAtValid(location, false);
			}
		}
	}
	return data;
}

#include "utils/uuid.h"

InferenceTree* AbstractFeatureExtraction::constructLayer(IntermediateData& data, int border) {

	int width = data.max_feature_activations.getWidth();
	int height = data.max_feature_activations.getHeight();

	cv::Size2i new_size(width + 2*border, height + 2*border);

	const int layer = 0;

	InferenceTree* result_tree = new InferenceTree();

	UUIDGenerator& inferred_part_uuid_gen = result_tree->getInferredPartUUIDGenerator();

	InferenceLayer* first_layer = new InferenceLayer(*result_tree, new_size, layer);

	result_tree->insertNewLayer(first_layer);

	FeatureActivationArray& max_activations = data.max_feature_activations;
	
	// iterate over all feature types to count how many parts we need to allocate
	/*int number_parts = 0;
	for (auto feature_type_iter = data.feature_activations.begin(); feature_type_iter != data.feature_activations.end(); ++feature_type_iter) {
		for (auto location = feature_type_iter->beginLocation(); location != feature_type_iter->endLocation(); ++location) {
			number_parts += feature_type_iter->isCandiateAtValid(location);
		}
	}*/

	for (auto feature_type_iter = data.feature_activations.begin(); feature_type_iter != data.feature_activations.end(); ++feature_type_iter) {
		
		VocabularyPart& corresponding_vocabulary_part = feature_type_iter->getTypeId();

		FeatureActivationArray& feature_type_activations = *feature_type_iter;

		auto location = feature_type_activations.beginLocation(); 		

		for (int y = 1; y < height - 1; ++y) {
			for (int x = 1; x < width - 1; ++x) {
				if (feature_type_activations.isCandiateAtValid(location)) {
					float activation_value = feature_type_activations.getCandidateAt(location).value;
					
					// get x,y location with added border
					cv::Point2i location(x + border, y + border);
					
					ResponsesArray default_responses;
					default_responses.set(ResponseType::R_RESPONSE_1, activation_value);
					default_responses.set(ResponseType::G_RESPONSE_1, 1.0);
					default_responses.set(ResponseType::RR_RESPONSE_1, 1.0);
					default_responses.set(ResponseType::S_RESPONSE_1, 0.0);

					UUIDType uuid = inferred_part_uuid_gen.generateUUID();

					first_layer->insertNewPart(new InferredPart(uuid, layer, location, corresponding_vocabulary_part.getUUID(), default_responses));
				}
				++location;
			}
			// skip border locations
			location+=2; 
		}
	}

	// TODO: save information about config/inital image into layer1_result (border and org image) !!!

	return result_tree;
}


std::shared_ptr<VocabularyTree> AbstractFeatureExtraction::getInitialVocabulary() {
	// return exising vocabulary if present
	if (vocabulary != nullptr) {
		return vocabulary;
	}

	vocabulary = std::shared_ptr<VocabularyTree>(new VocabularyTree());

	UUIDGenerator& vocabulary_part_uuid_gen = vocabulary->getVocabularyPartUUIDGenerator();

	// obtain definitions of all features from the actual implemention
	std::vector<FeatureTypeDefinition> feature_types = getFeatureTypeDefinitions();
		
	double contraction_factor = 0;

	VocabularyLayer* first_layer = new VocabularyLayer(*vocabulary,contraction_factor, 0);

	for (auto iter = feature_types.begin(); iter != feature_types.end(); ++iter) {
		FeatureTypeDefinition& feature_type = *iter;

		VocabularyPart::mask_t mask;
		VocabularyPart::region_t region;
			
		feature_type.reconstruction_image.to_point_map(mask, 0.0); // TODO: need to change this somehow
		feature_type.activation_points.to_point_set(region, false);

		int layer = 0;
		UUIDType uuid = vocabulary_part_uuid_gen.generateUUID();

		first_layer->insertNewPart(new VocabularyPart(uuid, layer, feature_type.id, mask, region));
	}

		
	vocabulary->insertNewLayer(first_layer);

	return vocabulary;
}