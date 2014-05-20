// abstract_layer_inference
///////////////////////////////////////////////////////////////////////////////

#include "abstract_inference.h"

void register_abstract_inference_classes() {
	/************************************************************************* 
	 *           !!!! DO NOT MODIFY OR REMOVE THIS FUNCTION   !!!!!
	 ************************************************************************/
	ClassRegister::registerFactory<AbstractLayerInference::IFactory, MultipleLayersInference::Factory>();
}


InferenceLayer* AbstractLayerInference::promoteCandidatesIntoLayer(InferredPartCandidateList* candidates_list, int layer, double layer_contraction) {
	// pre-compute (multiplication is faster then division)
	double layer_contraction_inverse = 1 / (double)layer_contraction;

	cv::Size2i new_layer_size = candidates_list->current_layer_inference.getSize();

	// perform layer contraction
	new_layer_size.width = cvRound(new_layer_size.width * layer_contraction_inverse);
	new_layer_size.height = cvRound(new_layer_size.height * layer_contraction_inverse);

	InferenceLayer* new_inference_layer = new InferenceLayer(candidates_list->current_layer_inference.getInferenceTree(), new_layer_size, layer);

	std::vector<IAttachableClass*> removed_candidates;

	auto& candidates = candidates_list->candidates;

	//	#pragma omp parallel for
	for (int i = 0; i < candidates.size(); ++i) {
		InferredPartCandidate& part_candidate = candidates[i];
		//for (auto iter = candidates_list->candidates.begin(); iter != candidates_list->candidates.begin(); ++iter) {
		//	InferredPartCandidate& part_candidate = *iter;
		if (part_candidate.invalid == false) {
			// promote candidate into new part

			cv::Point2i central_part = part_candidate.current_inferred_part->getLocation();
			cv::Point2i vocabulary_part_center_of_mass = part_candidate.candidate_vocabulary_part->getRelativeCenterOfMass();

			// location of new part is location of previous part + center of mass of vocabulary part			
			cv::Point2i new_part_location = central_part + vocabulary_part_center_of_mass;

			// perform layer contraction
			new_part_location.x = cvRound(new_part_location.x * layer_contraction_inverse);
			new_part_location.y = cvRound(new_part_location.x * layer_contraction_inverse);

			// construct and push new part into the layer
			// must used the same UUID since functionalities will point to this ID after they are moved from InferredPartCandidateList to InferrenceTree
			new_inference_layer->insertNewPart(new InferredPart(part_candidate.getUUID(), 
				layer, 
				new_part_location, 
				part_candidate.candidate_vocabulary_part->getUUID(),
				part_candidate.calculated_responses));
		} else {
			// mark it as removed
			removed_candidates.push_back(&part_candidate);
		}
	}

	// make sure we eliminate any reference pointers within any extension that is using it
	candidates_list->deleteAttachedReferences<InferenceTree>(removed_candidates);

	// finally copy all the extension attached to the InferredPartCandidateList during the construction into the new InferenceLayer*
	candidates_list->moveExtensionsTo(new_inference_layer->getInferenceTree());

	return new_inference_layer;
}

AbstractLayerInference* MultipleLayersInference::Factory::newInstance(const IConfiguration& config, const std::shared_ptr<VocabularyTree> vocabulary) const {

	int layer_index = config.getInt("layer_index", 0);
	int start_layer = -1;
	int end_layer = -1;

	if (layer_index > 0) {
		start_layer = config.getInt("start_layer", -1);
		end_layer = config.getInt("end_layer", -1);	

		// if both start_layer and end_layer are in config and at the same time layer_index is defined, then return error
		if (start_layer > 0 || end_layer > 0) {
			throw custom_libhop_exception(ConfigException, "Error: Ambiguous definition of layer. Found 'layer_index' and 'start_layer','end_layer' in config. Use either 'layer_index' or 'start_layer' and 'end_layer'.");
		} else {
			// if layer_index was define, then start_layer and end_layer should not be in config, so continue normaly
			start_layer = layer_index;
			end_layer = layer_index;	
		}
	} else {
		start_layer = config.getInt("start_layer");
		end_layer = config.getInt("end_layer");	
	}

	// validate start and end layer values
	if (start_layer <= 1) {			
		std::cout  << "Warning: Found start_layer that is <= 1 but create_layern cannot process 1. layer ... setting start_layer to 2" << std::endl;
		start_layer = 2;
	}
	// TODO: how should we handle this (maybe we can create NullInference that does nothing??)
	if (end_layer < start_layer) {
		// DO NOT throw exception so no error is reporter (required by FileImageInferenceJob in rosette/leoparts)
		std::cout << "Found start_layer > end_layer. Unable to continue processing" << std::endl;
		return nullptr;
	}

	std::vector<AbstractLayerInference*> layer_inferences_list;
	
	for (int i = start_layer; i <= end_layer; ++i) {
		std::string ly_namespace = string("ly") + i;
		std::string key = string("end_layer") + i;

		// get i-th layer configuration
		IConfiguration* layer_inference_config = config.getNamespace(ly_namespace.c_str());
		
		// find out what kind of inference type this layer should use (TODO: maybe we can check with the vocabulary if this will be valid type ???)
		std::string layer_inference_type = layer_inference_config->getString("type");

		// get factory for this inference layer (based on a type)
		AbstractLayerInference::IFactory* layer_inference_factory = ClassRegister::get().retrieveFactory<AbstractLayerInference::IFactory>(layer_inference_type);

		// create layer inference
		AbstractLayerInference* layer_inference = layer_inference_factory->newInstance(*layer_inference_config, vocabulary);

		// and push it into the list
		layer_inferences_list.push_back(layer_inference);
	}

	return new MultipleLayersInference(layer_inferences_list);
}


LayerOutputObject*  MultipleLayersInference::performInference(const AbstractLayerInputObject* input_object) {

	// save groundtruths so we do not have to load them each time
	std::vector<GroundtruthObject> input_groundtruths = input_object->getGroundtruths();

	const AbstractLayerInputObject* intermediate_input = input_object;

	LayerOutputObject* intermediate_output = nullptr;	

	// call individual layer inference processes and handle intput/output between them
	for (auto iter = layer_inferences.begin(); iter != layer_inferences.end(); ++iter) {
		AbstractLayerInference* layer_inference = *iter;

		// we need to delete any previous memory for output
		if (intermediate_output != nullptr) delete intermediate_output;

		// process with current layer
		intermediate_output = layer_inference->performInference(intermediate_input);

		// prepare as input for next layer
		intermediate_input = new MemoryLayerInputObject(intermediate_output->getLayerObject(), input_groundtruths, intermediate_output->getGroupMap());
	}	

	// also delete intermediate input object but only if its not the same as original input_object
	if (intermediate_input != input_object) {
		delete intermediate_input;
	}

	return intermediate_output;
}