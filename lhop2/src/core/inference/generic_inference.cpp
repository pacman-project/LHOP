// generic_layer_inference
///////////////////////////////////////////////////////////////////////////////

#include "generic_inference.h"

#include "core/inference/indexing_blocks/center_indexing.h"
#include "core/inference/matching_blocks/schur_product_matching.h"
#include "core/inference/selection_blocks/shape_check_selection.h"

void register_generic_inference_classes() {
	/************************************************************************* 
	 *           !!!! DO NOT MODIFY OR REMOVE THIS FUNCTION   !!!!!
	 ************************************************************************/
	ClassRegister::registerFactory<AbstractLayerInference::IFactory, GenericLayersInference::Factory>();
}

AbstractLayerInference* GenericLayersInference::Factory::newInstance(const IConfiguration& config, const std::shared_ptr<VocabularyTree> vocabulary) const {
	
	// TODO: should we have configuration for individual blocks in their factories (do they even have any factories)??

	///////////////////////////////////////////////////
	/// Indexing settings

	// settings for the candidate thresholding block
	double candidate_r_threshold = config.getDouble("candidate_r_threshold");
	double candidate_g_threshold = config.getDouble("candidate_g_threshold");
	
	std::vector<int> allowed_parts_vector = config.getIntArray("allowed_parts", std::vector<int>());
	std::set<int> allowed_parts(allowed_parts_vector.begin(), allowed_parts_vector.end());


	CandidateThresholdingBlock candidate_thresholding_block(candidate_r_threshold, candidate_g_threshold, allowed_parts);

	///////////////////////////////////////////////////
	/// Schur product matching settings

	// settings for the schur product block
	double convolution_threshold = config.getDouble("convolution_threshold", 0.1);

	IdentitySchurProductFunction schur_product_function;
	SchurProductBlock<IdentitySchurProductFunction> schur_product_block(schur_product_function, convolution_threshold);
	
	// settings for the response calculation block
	double r_response_pow = config.getDouble("r_response_pow", 1.0);
	double g_response_pow = config.getDouble("g_response_pow", 1.0);
	bool g_response_operation = config.getDouble("g_response_operation", true);

	RR_ResponseCalculationBlock rr_calculation_block;	
	G_ResponseCalculationBlock g_calculation_block(rr_calculation_block, g_response_pow, g_response_operation);
	R_ResponseCalculationBlock r_calculation_block(r_response_pow);
	
	ResponseCalculationBlock response_calculation_block(r_calculation_block, g_calculation_block);

	// settings for the candidate's responses validation block
	ResponseCandidateValidationBlock candidate_validation_block;
	
	///////////////////////////////////////////////////
	/// Selection per location settings

	/// constructing final building blocks
	IInferredCandidatesIndexingBlock* part_indexing = new IndexingBlock<CandidateThresholdingBlock>(candidate_thresholding_block);
	IInferredCandidatesFilteringBlock* part_matching = new SchurProductMatchingBlock<IdentitySchurProductFunction,
																	ResponseCalculationBlock,
																	ResponseCandidateValidationBlock>(schur_product_block, 
																										response_calculation_block, 
																										candidate_validation_block);
	IInferredCandidatesFilteringBlock* part_selection_per_location = new SelectionPerLocationBlock();

	int layer = config.getDouble("layer"); 
	
	// by default all configuration index layer with 1 as first layer so convert to 0-based index
	layer = layer - 1;

	// TODO: read layer contraction from config instead of using it from the vocabulary as vocabulary may not have correct values
	//int layer_contraction = vocabulary->getLayer(layer).getContrationFactor();
	double layer_contraction = config.getDouble("layer_contraction");

	return new GenericLayersInference(part_indexing, part_matching, part_selection_per_location, 
										layer, layer_contraction, vocabulary);
}


LayerOutputObject* GenericLayersInference::performInference(const AbstractLayerInputObject* input_object) {

	std::shared_ptr<InferenceTree> input_parse_tree = input_object->getLayerObject();
		
	InferenceLayer& current_inference_layer = input_parse_tree->getLayer(layer - 1);

	// generate possible candidates for new layer using indexing block
	InferredPartCandidateList* candidates = part_indexing->obtainInferenceCandidates(current_inference_layer, *vocabulary);

	// perform matching on each candidate using matching block
	candidates = part_matching->filterInferenceCandidates(candidates, *vocabulary); 

	// perform selection on remaining candidates using selection per location block
	candidates = part_selection_per_location->filterInferenceCandidates(candidates, *vocabulary); 

	// promote candidates into final layer parts
	InferenceLayer* new_inference_layer = promoteCandidatesIntoLayer(candidates, this->layer, this->layer_contraction);

	input_parse_tree->insertNewLayer(new_inference_layer); // TODO: what if layer with the same number already existed ??

	return new LayerOutputObject(input_parse_tree, input_object->getGroupMap());
}

