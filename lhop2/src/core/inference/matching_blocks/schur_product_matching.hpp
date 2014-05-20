// schur_product_matching
///////////////////////////////////////////////////////////////////////////////

#include "core/inference/matching_blocks/schur_product_matching.h"

#include "core/structures/subparts/subparts_extension.h"
#include "core/structures/shape_context/appearance_extension.h"

template <class T1, class T2, class T3>
InferredPartCandidateList* SchurProductMatchingBlock<T1,T2,T3>::filterInferenceCandidates(InferredPartCandidateList* inference_candidates, const VocabularyTree& vocabulary) {

	// get access to subparts of the vocabulary
	VocabularySubpartsAccess vocabulary_subparts_access = vocabulary.getAccess<VocabularySubpartsAccess>();
	
	// also set appearance access for the schur product block (it will use any appearance weights found during the schur product)
	VocabularyAppearanceAccess vocabulary_appearance_access = vocabulary.getAccess<VocabularyAppearanceAccess>();
	// set vocabulary_appearance_access pointer for the schur_product (but as it points to stack it will only be valid till the end of this call !!!)
	schur_product.setVocabularyAppearacnceAccess(vocabulary_appearance_access);

	// get pointer to list of inferred parts on current layer
	InferenceLayer& current_layer_inference = inference_candidates->current_layer_inference;

	// insert and get modifier for subparts of the inference tree in construction (i.e. InferenceSubpartsModifier modifier for the InferredPartCandidateList*)
	inference_candidates->insertExtension(new InferenceSubpartsExtension(*inference_candidates));

	// we will create subparts links so we need modifier for inference subparts
	InferenceSubpartsModifier inference_subparts_modifier = inference_candidates->getAccess<InferenceSubpartsModifier>();

	// get reference (!! has to be reference/pointer) to the generator for the inference subpart data
	// ERROR: !!!! we should bet UUID generator from the actual InferenceTree so we can use correct global IDs 
	UUIDGenerator& subpart_uuid_gen = inference_subparts_modifier.getInferrenceSubpartUUIDGenerator();

	auto& candidates = inference_candidates->candidates;
	
//	#pragma omp parallel for
	for (int i = 0; i < candidates.size(); ++i) {
		InferredPartCandidate& part_candidate = candidates[i];
	//for (auto iter = candidates.begin(); iter != candidates.end(); ++iter) {
	//	InferredPartCandidate& part_candidate = *iter;

		InferredPart& current_center_part = *(part_candidate.current_inferred_part);
		VocabularyPart& candidate_vocabulary_part = *(part_candidate.candidate_vocabulary_part);

		// save information about each subpart that will be realized
		std::vector<SchurProductResult*> schur_product_result_list;

		// find all subparts for this candidate vocabulary part through appropriate access class
		std::vector<VocabularySubpartData*> candidate_vocabulary_subparts = vocabulary_subparts_access.getSubparts(candidate_vocabulary_part);
		
		////////////////////////////////////////////////////////////////////////
		/// Matching section
		
		// for each vocabulary subpart perform schur matching around the current inferred part
		for (auto supbart_iter = candidate_vocabulary_subparts.begin(); supbart_iter != candidate_vocabulary_subparts.end(); ++supbart_iter) {
			VocabularySubpartData& vocabulary_subpart = **supbart_iter;

			// TODO: matching will happen even with the central node - how do we handle it? 
			// we do not need to match central node but we can use schur_product.match to find all the central parts
			// maybe we can have dummy distribution map for the central nodes and we do not have to handle anything from here

			// compute with the provided schur product implementation
			schur_product_result_list.push_back(schur_product.match(current_layer_inference, current_center_part, vocabulary_subpart));
		}
		
		////////////////////////////////////////////////////////////////////////
		/// Response calculation section

		// calculate response values from matching results of all subparts
		part_candidate.calculated_responses = response_calculation_block.updateResponse(part_candidate.calculated_responses, schur_product_result_list, part_candidate);

		////////////////////////////////////////////////////////////////////////
		/// Verification and save section

		// verify that the constructed part is valid (using thresholding of the calculated response values)
		if (candidate_validation_block.isCandidateValid(part_candidate) == true) {
			
			if (convolution_link_threshold < 1) {
				// if user-provided convolution_link_threshold is valid then
				// for each realized subpart filter-out links to subparts that do not have high enough schur product value	

				for (auto supbart_iter = schur_product_result_list.begin(); supbart_iter != schur_product_result_list.end(); ++supbart_iter) {
					SchurProductResult& subpart_schur_product = **supbart_iter;

					SchurProductPart* max_subpart = subpart_schur_product.getMaxSubpart();
				
					if (max_subpart != nullptr) {
						// as threshold use max schur product value multiplied with the user provided convolution_link_threshold value
						double relative_convolution_link_threshold = convolution_link_threshold * max_subpart->schur_product_value;

						// remove any element of subpart_schur_product that has is schur_product_value below relative_convolution_link_threshold
						std::remove_if(subpart_schur_product.valid_subparts.begin(), subpart_schur_product.valid_subparts.end(),
										// lambda function: returns true if schur_product_value is below relative_convolution_link_threshold
										[relative_convolution_link_threshold](SchurProductPart& schur_part) {
											return schur_part.schur_product_value <= relative_convolution_link_threshold;
										});
					}

				}
			}

			std::vector<InferenceSubpartData*> subpart_links;

			// save additional information needed for the promotion of the candidate into the actual part
			// use InferenceSubpartsModifier to insert the links
			for (auto supbart_iter = schur_product_result_list.begin(); supbart_iter != schur_product_result_list.end(); ++supbart_iter) {
				SchurProductResult& subpart_schur_product = **supbart_iter;

				int subpart_index = subpart_schur_product.vocabulary_subpart.index;
				
				// save all links returned by the schur product and associate them with the appropriate vocabulary subpart (i.e. using subpart index)
				for (auto iter_schur = subpart_schur_product.valid_subparts.begin(); iter_schur != subpart_schur_product.valid_subparts.end(); ++iter_schur) {
					subpart_links.push_back(new InferenceSubpartData(subpart_uuid_gen.generateUUID(),
																	*iter_schur->part,					// -> inferred part that survived the convolution thresholding
																	subpart_index,						// -> index of subpart relative to the vocabulary part we are inferring
																	iter_schur->schur_product_value));	// -> schur product value of the part
				}
				// insert subparts but just use candidate as reference because InferredPart will not be created until the end of inference
				// (we then need to make sure we use the same UUID of candidate also for new inferred part !!)
				inference_subparts_modifier.insertSubparts(part_candidate, subpart_links);
			}

			// TODO: when/where do we add reconstruction links (i.e. links to the first layer parts) ???
		} else {
			// if not valid then mark it as invalid
			part_candidate.invalid = true;
		}
	}

	return inference_candidates;
}

///////////////////////////////////////////////////////////////
/// Schur product implementations used by candidate matching

template <class T>
SchurProductResult* SchurProductBlock<T>::match(InferenceLayer& inferred_layer, InferredPart& inferred_center_part, VocabularySubpartData& candidate_vocabulary_part) {
	
	SchurProductResult* result = new SchurProductResult(candidate_vocabulary_part);

	int c = 1;
	
	// obtain distribution matrix (i.e. a matrix of gaussian map) produced during the learning stage of this subpart
    const cv::Mat& distribution_matrix = candidate_vocabulary_part.map_distribution;
	
	int dist_width = distribution_matrix.cols;
	int dist_height = distribution_matrix.rows;

	// get rectangular region around which we will do the matching
	// by moving from the current central part by an offset defined by the candidate_vocabulary_part
	cv::Point2i center_loc = (inferred_center_part.getLocation() + candidate_vocabulary_part.offset) * c;

    int x = center_loc.x - dist_width/2;
	int y = center_loc.y - dist_height/2;

	cv::Size2i layer_size = inferred_layer.getSize();

	// bounds check	
	if (x >= layer_size.width || x < 0 || 
		y >= layer_size.height || y < 0) 
		return result;

	int num_skiped_row_parts = layer_size.width - (int)dist_width;
	num_skiped_row_parts -= 1; // decrement by one, since ++region_iterator will handle that part

	// get appearance model for vocabulary subpart - we will ignore any types of inferred parts that are not in the appearance model
	auto subpart_appearance = vocabulary_appearance_access->getSubpartAppearance(candidate_vocabulary_part);

	// position the iterator (of parts in the current layer) at the start of the region
	InferenceLayer::Iterator region_iterator = inferred_layer.beginIteratorAt(x, y);
	
	// move both the distribution matrix index and the region over the layer parts index one by one
	for (int j = 0; j < dist_height; ++j, region_iterator+=num_skiped_row_parts) {
		for (int i = 0; i < dist_width; ++i, ++region_iterator) {
			// get all parts at this location
			std::vector<InferredPart*> all_parts_at_location = inferred_layer.getPartsAt(region_iterator);

			// perform calculation for each part
			for (auto part_iter = all_parts_at_location.begin(); part_iter != all_parts_at_location.end(); ++part_iter) {
				InferredPart& part = **part_iter;

				UUIDType corresponding_vocabulary_part_uuid = part.getCorrespondingVocabularyPartUUID();

				// ignore part if its type is not in the appearance model
				auto similar_apperance = subpart_appearance.find(corresponding_vocabulary_part_uuid);
				if (similar_apperance == subpart_appearance.end())
					continue;

                // calculate value as schur_product * distribution_matrix * app weight
				double value = schur_product_function.getValue(inferred_center_part, part) * distribution_matrix.at<float>(i,j) * similar_apperance->second.second;

				// ignore part if response is below convolution threshold
				if (value <= convolution_threshold)
					continue;

				// save the subpart and its schur product response if it passed all the criteria
				result->valid_subparts.push_back(SchurProductPart(&part, value));
            }
        }
    }

    return result;
}

///////////////////////////////////////////////////////////////
/// Types of schur product functions

double IdentitySchurProductFunction::getValue(const InferredPart& central_part, const InferredPart& part) const { 
	return 1.0; 
}

double R_ResponseSchurProductFunction::getValue(const InferredPart& central_part, const InferredPart& part) const {
	return part.getResponse(ResponseType::R_RESPONSE_1);
}

double V_ResponseSchurProductFunction::getValue(const InferredPart& central_part, const InferredPart& part) const {
	return part.getResponse(ResponseType::G_RESPONSE_1);
}

void G_ResponseSchurProductFunction::initialize(const InferredPart& current_center_part, const VocabularySubpartData& vocabulary_subpart){
	this->quotient = log(current_center_part.getResponse(ResponseType::R_RESPONSE_1)); 	// TODO: this is computed too many times (should cache the value somehow !!!)
	this->dist = vocabulary_subpart.gdistr;
	
	if (g_response_var_factor > 0.0) {
		this->dist.reset_variance(this->dist.get_variance() * g_response_var_factor);
	}
}

double G_ResponseSchurProductFunction::getValue(const InferredPart& central_part, const InferredPart& part) const {
	double r_response = part.getResponse(ResponseType::R_RESPONSE_1);
	double g_response = part.getResponse(ResponseType::G_RESPONSE_1);
	double rr_response = part.getResponse(ResponseType::RR_RESPONSE_1);
	return dist.pdf_val1(::log(r_response) - quotient) * g_response/rr_response; 
}

double Simple_G_ResponseSchurProductFunction::getValue(const InferredPart& part) const {
	return part.getResponse(ResponseType::G_RESPONSE_1)/part.getResponse(ResponseType::RR_RESPONSE_1);
}

///////////////////////////////////////////////////////////////
/// Response calculations used by candidate matching (TODO: maybe move into another file ??)


ResponsesArray& R_ResponseCalculationBlock::updateResponse(ResponsesArray& current_responses, const std::vector<SchurProductResult*>& schur_product_results, const InferredPartCandidate& part_candidate) {
	// calculate R_RESPONSE as sum over R_RESPONSE of all realized subparts (use max part from all schur product parts)

	ResponsesArray::ValueType r_sum = 0;

	for (auto iter_subpart = schur_product_results.begin(); iter_subpart != schur_product_results.end(); ++iter_subpart) {
		SchurProductResult& subpart_schur_product = (SchurProductResult&)**iter_subpart;
		
		SchurProductPart* max_subpart = subpart_schur_product.getMaxSubpart();

		if (max_subpart != nullptr)
			r_sum += max_subpart->part->getResponse(ResponseType::R_RESPONSE_1);
	}
	int subpart_count = schur_product_results.size();

    r_sum =  r_sum / (subpart_count + 1);
    r_sum = ::pow(r_sum, r_response_pow);

	current_responses.set(ResponseType::R_RESPONSE_1, r_sum);

	return current_responses;
}

ResponsesArray& RR_ResponseCalculationBlock::updateResponse(ResponsesArray& current_responses, const std::vector<SchurProductResult*>& schur_product_results, const InferredPartCandidate& part_candidate) {
	// count number of basic parts realized
	ResponsesArray::ValueType realized_basic_parts_count = 0;
	
	for (auto iter_subpart = schur_product_results.begin(); iter_subpart != schur_product_results.end(); ++iter_subpart) {
		SchurProductResult& subpart_schur_product = (SchurProductResult&)**iter_subpart;
		
		// count number of realized basic parts
		SchurProductPart* max_subpart = subpart_schur_product.getMaxSubpart();

		if (max_subpart != nullptr)
			realized_basic_parts_count += subpart_schur_product.vocabulary_subpart.subpart.getBasicPartCount() * max_subpart->part->getResponse(ResponseType::RR_RESPONSE_1);
	}

	realized_basic_parts_count = realized_basic_parts_count /  part_candidate.candidate_vocabulary_part->getBasicPartCount();

	current_responses.set(ResponseType::RR_RESPONSE_1, realized_basic_parts_count);

	return current_responses;
}

ResponsesArray& G_ResponseCalculationBlock::updateResponse(ResponsesArray& current_responses, const std::vector<SchurProductResult*>& schur_product_results, const InferredPartCandidate& part_candidate) {
	
	// use RR_ResponseCalculationBlock to calculate the RR_RESPONSE first
	current_responses = rr_response_block.updateResponse(current_responses, schur_product_results, part_candidate);

	// then calculate G_RESPONSE
	ResponsesArray::ValueType g_val = (g_response_operation == 0) ? 0.0 : 1.0; 

	for (auto iter_subpart = schur_product_results.begin(); iter_subpart != schur_product_results.end(); ++iter_subpart) {
		SchurProductResult& subpart_schur_product = (SchurProductResult&)**iter_subpart;
		
		SchurProductPart* max_subpart = subpart_schur_product.getMaxSubpart();

		if (max_subpart != nullptr) {
			if (g_response_operation == 0) 
				g_val += max_subpart->schur_product_value;
			else 
				g_val *= max_subpart->schur_product_value;
		}
	}
	
	// make sure RR_RESPONSE is caluclated first by using calling rr_response_block before
	g_val = g_val * current_responses.get(ResponseType::RR_RESPONSE_1);

	g_val = ::pow(g_val, g_response_pow);

	current_responses.set(ResponseType::G_RESPONSE_1, g_val);

	return current_responses;
}

///////////////////////////////////////////////////////////////
/// Candidate validation check based on calculated response values

bool R_ResponseCandidateValidationBlock::isCandidateValid(const InferredPartCandidate& part_candidate) {
	ResponsesArray::ValueType r_response = part_candidate.calculated_responses.get(ResponseType::R_RESPONSE_1);
	ResponsesArray::ValueType candidate_vocabulary_r_threshold = part_candidate.candidate_vocabulary_part->getResponseThreshold(ResponseType::R_RESPONSE_1);

	return r_response > candidate_vocabulary_r_threshold ? true : false;
}

bool RR_ResponseCandidateValidationBlock::isCandidateValid(const InferredPartCandidate& part_candidate) {
	ResponsesArray::ValueType rr_response = part_candidate.calculated_responses.get(ResponseType::RR_RESPONSE_1);
	ResponsesArray::ValueType candidate_vocabulary_rr_threshold = part_candidate.candidate_vocabulary_part->getResponseThreshold(ResponseType::RR_RESPONSE_1);

	return rr_response > candidate_vocabulary_rr_threshold ? true : false;
}

bool G_ResponseCandidateValidationBlock::isCandidateValid(const InferredPartCandidate& part_candidate) {
	ResponsesArray::ValueType g_response = part_candidate.calculated_responses.get(ResponseType::G_RESPONSE_1);
	ResponsesArray::ValueType candidate_vocabulary_g_threshold = part_candidate.candidate_vocabulary_part->getResponseThreshold(ResponseType::G_RESPONSE_1);

	return g_response > candidate_vocabulary_g_threshold ? true : false;
}
