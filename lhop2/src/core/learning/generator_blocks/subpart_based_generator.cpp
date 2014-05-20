
// subpart_based_generator
///////////////////////////////////////////////////////////////////////////////

#include "subpart_based_generator.h"

#include "core/structures/subparts/support_extension.h"

std::string VocabularySubpartDefinition::getUniqueIdentifier() const {
	std::stringstream ss;
	ss << vocabulary_subpart_uuid;
	// make sure to check if local_maxima is defined (central part does not have local maxima)
	if (local_maxima != nullptr)
		ss << local_maxima->offset;
	return ss.str();
}
VocabularyPartCandidateDefinition::VocabularyPartCandidateDefinition(UUIDType central_part_uuid, const VocabularySubpartsCombination& subparts_combination) 
	: central_part_uuid(central_part_uuid) {
	// copy to list of subpart definitions
	// ORDER IS IMPORTANT !! but we can simply add them since subparts_combination.subparts 
	// is already sorted by the subpart definition
	subparts_definitions.reserve(subparts_combination.subparts.size());
	for (auto iter = subparts_combination.subparts.begin(); iter != subparts_combination.subparts.end(); ++iter)
		subparts_definitions.push_back(iter->subpart_definition);

	// generate unique identifier based on ordered list of subpart definitions
	// use vocabulary_part_uuid and local maxima offset as identifier of each subpart
	std::stringstream ss;
	for (auto iter = subparts_definitions.begin(); iter != subparts_definitions.end(); ++iter) {
		ss << iter->getUniqueIdentifier();
	}
	unique_identifier = ss.str();
}

VocabularyPartCandidateList* SubpartBasedCandidatesGenerator::obtainVocabularyCandidates(VocabularyLayer& current_vocabulary_layer, const std::vector<AbstractLayerInputObject*> input_object) {

	// collect co-occurrence duplets statistics (map learning)
	// For each combination of two vocabulary parts (a pairs) we get a distribution
	// matrix that describes a frequency of spatial occurrences of one type of part (surrounding part)
	// around another type of part (central part).
	// the first method called by generic learning
	for (auto iter = input_object.begin(); iter != input_object.end(); ++iter) {
		const AbstractLayerInputObject* learning_object = *iter;

		// read learning object
		std::shared_ptr<InferenceTree> inference_tree = learning_object->getLayerObject();

		// update vocabulary part co-occurrence statistics from a specific layer of this learning sample
		cooccurrence_statistics->updateStatistics(inference_tree->getLayer(current_vocabulary_layer));
	}
	// after cooccurrence_statistics is updated we find local maximas for each pair of parts 
	// this is achived by calls to the CooccurrenceStatistics::PairwisePart::getLocalMaximas()
	// during the updateCandidatePartStatistics call

	// combine co-occurrence duplets to create new parts with multiple subparts
	for (auto iter = input_object.begin(); iter != input_object.end(); ++iter) {
		const AbstractLayerInputObject* learning_object = *iter;

		// read learning object
		std::shared_ptr<InferenceTree> inference_tree = learning_object->getLayerObject();

		// update new vocabulary part candidate statistics from a specific layer of this learning sample
		updateCandidatePartStatistics(inference_tree->getLayer(current_vocabulary_layer));
	}	

	// promote all generated VocabularyPartCandidateDefinition into VocabularyPartCandidate and return it
	return promoteVocabularyPartCandidateDefinitions(current_vocabulary_layer);
}

void SubpartBasedCandidatesGenerator::updateCandidatePartStatistics(const InferenceLayer& inference_layer) {
	
	InferenceSupportAccess layer1_support_parts_access = ((InferenceLayer&)inference_layer).getInferenceTree().getAccess<InferenceSupportAccess>(InferenceSupportExtension::Factory());

	// prepare for update for each individual statistics update
	for (auto iter = individual_statistics_collectors.begin(); iter != individual_statistics_collectors.end(); ++iter) {
		IPartCandidateStatisticsCollector* individual_statistics = *iter;

		// TODO: having virtual methods within so many nested loops will be performance bottleneck !!!!!
		individual_statistics->prepareUpdate(inference_layer);
	}

	// iterate over all location of the inference_tree 
	for (auto loc_iter = inference_layer.beginIterator(); loc_iter != inference_layer.endIterator(); ++loc_iter) {
		// get list of parts at that location
		auto location_parts = inference_layer.getPartsAt(loc_iter);

		// and iterate over all parts at that location
		for (auto part_iter = location_parts.begin(); part_iter != location_parts.end(); ++part_iter) {
			
			const InferredPart& center_part = **part_iter;

			// use part only if its vval >= center_part_threshold: TODO: change vval into something meaningfull
			if (center_part.getVVal() < center_part_threshold)
				continue;

			double subpart_support_threshold = 0.5 * center_part.getVVal(); // TODO: why is it needed ? should 0.5 be in config ?

			////////////////////////////////////////////////////////////////////////////////////////////////
			// Get list of subparts supported around this central part by doing projection one
			// layer lower for each possible subpart and checking if subpart does exist in that layer
			const std::vector<VocabularySubpartSupport> supported_subparts_list = getSupportedSubpartsSortedBySupportValue(inference_layer, center_part);

			// Now generate and collect statistics for all possible combinations of supported subparts using center_part as its center.
			// The main do/while loop processes each possible new candidate part (that consists from central part and multiple supported
			// subparts) and generate those new candidates by extending current candidate with any subpart supported by the current sample (inference_layer).
			// All new part candidates are collected in stack subparts_combinations_stack and main do/while loop pops first one out
			// at each start of the loop. After new candidates are generated we push them into the start of the stack and allow them to be processed 
			// in the proceeding loops.
			std::stack<VocabularySubpartsCombination> subparts_combinations_stack;

			// first insert empty combination as initial subparts combination
			subparts_combinations_stack.push(VocabularySubpartsCombination());
			do {
				// pop first one from "stack", this creates depth-first search which should consume less memory since depth 
				// is more limited then width 
				// max depth = worst case is number of max subparts
				// max width = worst case is size of supported_subparts_list (max number of vocabulary parts at the current layer * max number of local maximas per co-occurrence matrix)
				VocabularySubpartsCombination& subparts_combination = subparts_combinations_stack.top();
				
				// if there is optimal number of subparts (between min_subpart_limit and max_subpart_limit) we can collect part statistics for it
				if (subparts_combination.subparts.size() > 0 && // ignore initial combination with empty list
					subparts_combination.subparts.size() >= min_subpart_count && 
					subparts_combination.subparts.size() <= max_subpart_count) {
					// update statistics
					updateSubpartCombinationStatisitcs(subparts_combination, inference_layer, center_part);
				}

				// if we have not yet exceeded the limit of subparts then create new possible subparts
				// combinations by extending current part candidate with any of the supported subparts found around current central part of the sample
				if (subparts_combination.subparts.size() < max_subpart_count) {
					// add only specific number of new combinations (to restrict explosion of parts) or all if max number is 0 or less
					// (add only the best ones - they should already be sorted by the getSupportedSubpartsSortedBySupportValue(..))
					int remaining_space = max_new_subpart_hypothesis_added > 0 ? max_new_subpart_hypothesis_added : supported_subparts_list.size();

					double subpart_support_relative_threshold = -1;

					// go through all subpart hypothesis supported by this central part and create new combination
					for (int i = 0; i < supported_subparts_list.size(); ++i) {
						const VocabularySubpartSupport& new_subpart = supported_subparts_list[i];

						// first verify that this subpart has not been already selected into the current subpart combination
						if (subparts_combination.subparts.count(new_subpart) > 0)
							continue;

						// verify value of the supported part is above threshold
						// wrt layer immediately below
						if (new_subpart.supported_value <= subpart_support_threshold)
							continue;

						// verify that the intersection between found subpart and all the other subparts is between threshold
						// use intersection of layer 1 supporting parts (layer 1 supporting parts == all parts at layer 1 if we follow subparts down to 1. layer)
						std::set<InferredPart*> new_subpart_ly1_support = layer1_support_parts_access.getOnlyDifferentInitialLayerSupportParts(*new_subpart.supported_part);

						int intersection_value = InferenceSupportExtension::calculatePartsIntersection(new_subpart_ly1_support, subparts_combination.ly1_support);

						if (intersection_value < min_intersection_threhsold || intersection_value > max_intersection_threhsold) 
							continue;

						// also, add only 90% of best remaining subparts (TODO: same parameter was hard-coded into old code - is it usefully ? )
						// calculate subpart_support_relative_threshold first if not found yet
						// WARNING: order of checks is important !! make subpart_support_relative_threshold only after all the others 
						// checks have been done
						if (subpart_support_relative_threshold < 0)
							// since supported_subparts_list is sorted by highest supported_value we can initialize 
							// subpart_support_relative_threshold at first loop
							subpart_support_relative_threshold = 0.9 * new_subpart.supported_value; 
						else if (new_subpart.supported_value <= subpart_support_relative_threshold)
							continue;

						// check if we have added more then allowed subparts for each loop
						if (--remaining_space <= 0) // TODO: should it be --remaining_space or remaining_space-- ?? should it be < or <= ??
							break;

						// create new hypothesis combination using previous combination and new supported subpart hypothesis
						subparts_combinations_stack.push(VocabularySubpartsCombination(subparts_combination, new_subpart, new_subpart_ly1_support));
					}
				}

				// move to next element in stack
				subparts_combinations_stack.pop();
			} while (subparts_combinations_stack.empty() == false); // finish if no combinations remain for checking
		}
	}
}

std::vector<VocabularySubpartSupport> SubpartBasedCandidatesGenerator::getSupportedSubpartsSortedBySupportValue(const InferenceLayer& inference_layer, const InferredPart& center_part) {
	
	CooccurrenceStatistics::CentralPart& distribution_statistics = cooccurrence_statistics->getCentralPartStatistics(center_part.getCorrespondingVocabularyPartUUID());

	// create priority queue for all vocabulary subpart definitions that whose support around this central part is found
	std::priority_queue<VocabularySubpartSupport,
		std::vector<VocabularySubpartSupport>,
		VocabularySubpartSupport::ComparerBySupportedValue> supported_subparts_queue;

	// find subparts from co-occurrence distribution matrix (using its maxima)
	for (auto pairwise_stat_iter = distribution_statistics.sourunding_parts.begin(); pairwise_stat_iter != distribution_statistics.sourunding_parts.begin(); ++pairwise_stat_iter) {
		UUIDType surrounding_part_uuid = pairwise_stat_iter->first;
		CooccurrenceStatistics::PairwisePart& surrounding_part_statistics = pairwise_stat_iter->second;

		auto local_maximas = surrounding_part_statistics.getLocalMaximas();

		for (int i = 0; i < local_maximas.size(); ++i) {
			VocabularySubpartDefinition subpart_definition(surrounding_part_uuid, &local_maximas[i]);

			// check if this subpart is found around the center part in the sample (i.e. if a sample part supports this new subpart definition)
			// and save its supporting part (from inference layer) and its value
			VocabularySubpartSupport subpart_support = verifySubpartSupport(subpart_definition, inference_layer, center_part);

			// make sure subpart is supported with at least one valid part
			if (subpart_support.supported_part == nullptr)
				continue;

			// and then add it to the priority queue that will sort them by their supported part value
			supported_subparts_queue.push(subpart_support);							
		}
	}

	// copy from queue to vector so we can iterate over them
	std::vector<VocabularySubpartSupport> supported_subparts_vector;
	while (supported_subparts_queue.empty() == false) {
		supported_subparts_vector.push_back(supported_subparts_queue.top()); supported_subparts_queue.pop();
	}

	return supported_subparts_vector;
}

VocabularySubpartSupport SubpartBasedCandidatesGenerator::verifySubpartSupport(const VocabularySubpartDefinition& subpart_definition, const InferenceLayer& inference_layer, const InferredPart& inferred_center_part) {
	
	// obtain distribution matrix (i.e. a matrix of gaussian map) produced during the learning stage of this subpart
	const cv::Mat distribution_matrix = subpart_definition.local_maxima->distribution_matrix;

	int dist_width = distribution_matrix.cols;
	int dist_height = distribution_matrix.rows;

	// get rectangular region around which we will do the matching
	// by moving from the current central part by an offset defined by the candidate_vocabulary_part
	cv::Point2i center_loc = (inferred_center_part.getLocation() + subpart_definition.local_maxima->offset) ;

	int x = center_loc.x - dist_width/2;
	int y = center_loc.y - dist_height/2;

	cv::Size2i layer_size = inference_layer.getSize();

	// bounds check	
	if (x >= layer_size.width || x < 0 || 
		y >= layer_size.height || y < 0) 
		return VocabularySubpartSupport(subpart_definition);

	int num_skiped_row_parts = layer_size.width - (int)dist_width;
	num_skiped_row_parts -= 1; // decrement by one, since ++region_iterator will handle that part

	// position the iterator (of parts in the current layer) at the start of the region
	InferenceLayer::Iterator region_iterator = inference_layer.beginIteratorAt(x, y);

	double supported_value = 0;
	InferredPart* supported_part;	

	// move both the distribution matrix index and the region over the layer parts index one by one
	for (int j = 0; j < dist_height; ++j, region_iterator+=num_skiped_row_parts) {
		for (int i = 0; i < dist_width; ++i, ++region_iterator) {
			// get all parts at this location
			std::vector<InferredPart*> all_parts_at_location = inference_layer.getPartsAt(region_iterator);

			// perform calculation for each part
			for (auto part_iter = all_parts_at_location.begin(); part_iter != all_parts_at_location.end(); ++part_iter) {
				InferredPart& part = **part_iter;

				// ignore any part of different type then the surrounding hypothesis subpart
				if (subpart_definition.vocabulary_subpart_uuid != part.getCorrespondingVocabularyPartUUID())
					continue;

				// calculate value as G_RESPONSE * distribution_matrix
				double value = part.getResponse(G_RESPONSE_1) * distribution_matrix.at<float>(i,j);

				// ignore parts with response below 50% of center part threshold (TODO: move 0.5 into const or config variable)
				if (value <= center_part_threshold * 0.5) 
					continue;
			
				// save only part with max supported value
				if (supported_value < value) {
					supported_value = value;
					supported_part = &part;
				}
			}
		}
	}

	return VocabularySubpartSupport(subpart_definition, supported_part, supported_value);
}

void SubpartBasedCandidatesGenerator::updateSubpartCombinationStatisitcs(const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& center_part) {

	VocabularyPartCandidateDefinition part_candidate_definition(center_part.getCorrespondingVocabularyPartUUID(), subparts_combination);

	const std::string part_candidate_unique_identifier = part_candidate_definition.getUniqueIdentifier();

	// register this vocabulary part candidate definition so we can later promote it into vocabulary part candidate 
	// (if not already registered)
	auto existing_definition = all_part_candidates_definitions.find(part_candidate_unique_identifier);
	if (existing_definition == all_part_candidates_definitions.end())
		all_part_candidates_definitions[part_candidate_unique_identifier] = part_candidate_definition;

	// then delegate statistics to specific classes for statistics gathering
	for (auto iter = individual_statistics_collectors.begin(); iter != individual_statistics_collectors.end(); ++iter) {
		IPartCandidateStatisticsCollector* individual_statistics = *iter;
		
		// TODO: having virtual methods within so many nested loops will be performance bottleneck !!!!!
		individual_statistics->updatePartCandidateStatisitcs(part_candidate_definition, subparts_combination, inference_layer, center_part);
	}
}

VocabularyPartCandidateList* SubpartBasedCandidatesGenerator::promoteVocabularyPartCandidateDefinitions(VocabularyLayer& vocabulary_layer) {

	VocabularyPartCandidateList* vocabulary_part_candidate_list = new VocabularyPartCandidateList();

	// use generator of the VocabularyPart since VocabularyPartCandidate will be later 
	// promoted into the VocabularyPart and will have to pass on any attached extension
	UUIDGenerator& vocabulary_part_candidate_generator = vocabulary_layer.getVocabularyTree().getVocabularyPartUUIDGenerator();

	// based on set of collected subparts combinations create new candidate part
	for (auto iter = all_part_candidates_definitions.begin(); iter != all_part_candidates_definitions.end(); ++iter) {
		const std::string& part_candidate_identifyer = iter->first;
		VocabularyPartCandidateDefinition& part_candidate_def = iter->second;

		// create part candidate with its 
		VocabularyPartCandidate* vocabulary_part_candidate = new VocabularyPartCandidate(vocabulary_part_candidate_generator.generateUUID());

		// insert part candidate
		vocabulary_part_candidate_list->candidates.push_back(vocabulary_part_candidate);

		part_candidate_def.created_vocabulary_part_candidate = vocabulary_part_candidate;
	}

	// then perform any other work such as
		// SubpartPartCandidateStatisticsCollector:			create subparts 
		// ShapeContextPartCandidateStatisticsCollector:	promote shape context information into shape context extension 
		// PCAPartCandidateStatisticsCollector:				promote PCA information into PCA extension ()
		// <TODO>:											generate part similarity (based on shape context or PCA)
	for (auto iter = individual_statistics_collectors.begin(); iter != individual_statistics_collectors.end(); ++iter) {
		IPartCandidateStatisticsCollector* individual_statistics = *iter;

		individual_statistics->finalizePartCandidateStatisitcs(all_part_candidates_definitions, vocabulary_part_candidate_list, vocabulary_layer);
	}

	return vocabulary_part_candidate_list;
}


void SubpartPartCandidateStatisticsCollector::prepareUpdate(const InferenceLayer& inference_layer) {
	// nothing to do here
}

void SubpartPartCandidateStatisticsCollector::updatePartCandidateStatisitcs(const VocabularyPartCandidateDefinition& part_candidate, const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& center_part) {

	// multiply max support value from each subpart hypothesis
	float value = 1;
	for (auto iter = subparts_combination.subparts.begin(); iter != subparts_combination.subparts.end(); ++iter)	{
		value = value * iter->supported_value;
	}

	// simply count multiplied values
	part_candidate_counter[part_candidate.getUniqueIdentifier()] = value;
}

void SubpartPartCandidateStatisticsCollector::finalizePartCandidateStatisitcs(std::unordered_map<std::string, VocabularyPartCandidateDefinition>& all_part_candidates_definitions, VocabularyPartCandidateList* vocabulary_part_candidate_list, VocabularyLayer& vocabulary_layer) {

	// TODO: we can use collected statistics from subparts to also eliminate part candidate if needed

	// get subparts modifier as we need to add subparts
	VocabularySubpartsModifier vocabulary_subparts_modifier = vocabulary_part_candidate_list->getAccess<VocabularySubpartsModifier>();

	// TODO: ERROR: we should be using subparts from vocabulary_tree and not from vocabulary_part_candidate_list !!!!!!!
	UUIDGenerator& vocabulary_subpart_generator = vocabulary_subparts_modifier.getVocabularySubpartUUIDGenerator();

	for (auto iter = all_part_candidates_definitions.begin(); iter != all_part_candidates_definitions.end(); ++iter) {
		const std::string& part_candidate_identifyer = iter->first;
		VocabularyPartCandidateDefinition& part_candidate_def = iter->second;

		// TODO: should we check if part with the same subparts definition already exists ???? 
		// (old implementation check it and returned reference to existing part)

		// then add each subpart to its extension
		std::vector<VocabularySubpartData*> subparts_data;

		// add central subpart
		subparts_data.push_back(new VocabularySubpartData(vocabulary_subpart_generator.generateUUID(),
														vocabulary_layer.getPart(part_candidate_def.central_part_uuid),
														// use index 0, empty map distribution and (0,0) offset 
														0, cv::Mat(), cv::Point2i(0,0), normal_distribution1()));

		// add other subparts
		int i = 1; // 
		for (auto iter = part_candidate_def.subparts_definitions.begin(); iter != part_candidate_def.subparts_definitions.begin(); ++iter) {
			VocabularySubpartDefinition& subparts_definition = *iter;

			VocabularySubpartData* new_subpart_data = new VocabularySubpartData(vocabulary_subpart_generator.generateUUID(),
				vocabulary_layer.getPart(subparts_definition.vocabulary_subpart_uuid),
				i++, 
				subparts_definition.local_maxima->distribution_matrix, 
				subparts_definition.local_maxima->offset,
				normal_distribution1());

			subparts_data.push_back(new_subpart_data);

			// also save created subpart data to subparts_definition
			subparts_definition.created_subpart_data = new_subpart_data;

		}

		// assign subparts to vocabulary part candidate
		vocabulary_subparts_modifier.insertSubparts(*part_candidate_def.created_vocabulary_part_candidate, subparts_data);

		// TODO: how are we going to handle forbidden parts
		// maybe the concepts of forbidden parts have to be handled differently since its 
		// only job is to prevent inference of part that is only a subset of another part
		// i.e. when part a is subset of part b then part b has higher priority of inference then part a
		// in the old implementation this was handled as forbidden part which are defined here

	}
}

