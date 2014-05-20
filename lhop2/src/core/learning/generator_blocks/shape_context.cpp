
// shape_context
///////////////////////////////////////////////////////////////////////////////

#include "shape_context.h"

////////////////////////////////////////////////////////////////////////
////// Shape context

InferredPart* ShapeContextStatistics::followCentralPartToLy1(const InferredPart& root, const InferenceLayer& inference_layer) {
	// TODO: !!! this is to SLOW since it can be called to many times !! we need to cache all access classes as members
	InferenceSubpartsAccess suparts_access = inference_layer.getInferenceTree().getAccess<InferenceSubpartsAccess>();

	InferredPart* current_node = (InferredPart*)&root;
	while (current_node->getLayer() > 0) {
		std::vector<InferenceSubpartData*> subparts = suparts_access.getSubparts(*current_node);
		// find central subpart
		for (auto iter = subparts.begin(); iter != subparts.end(); ++iter){
			if ((*iter)->index == 0) {
				current_node = &(*iter)->subpart;
				break;
			}
		}
	}
	return current_node;
}

void ShapeContextStatistics::prepareUpdate(const InferenceLayer& inference_layer) {
	// get access to Belongie histogram extension
	//BelongieHistogramAccess belongie_histogram_access = inference_layer.getInferenceTree().getAccess<BelongieHistogramAccess>(BelongieHistogramExtension::Factory());

	// TODO: we should get all accesses here but we cannot have BelongieHistogramAccess as member since it does not have default constructor
	// if we create it as pointer then we still have performance problems as all calls will be virtual !!!
}

void ShapeContextStatistics::updatePartCandidateStatisitcs(const VocabularyPartCandidateDefinition& part_candidate, const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& center_part) {

	// TODO: !!! this is to SLOW since it can be called to many times !! we need to cache all access classes as members
	BelongieHistogramAccess belongie_histogram_access = ((InferenceLayer&)inference_layer).getInferenceTree().getAccess<BelongieHistogramAccess>(BelongieHistogramExtension::Factory());
	InferenceSupportAccess ly1_support_access = ((InferenceLayer&)inference_layer).getInferenceTree().getAccess<InferenceSupportAccess>(InferenceSubpartsExtension::Factory());
	
	cv::Point2i ly1_center_location = followCentralPartToLy1(center_part, inference_layer)->getLocation();

	// get shape context data for this candidate part based on its unique identifier
	PartShapeContextData& candidate_part_shape_context = part_candidates_shape_context[part_candidate.getUniqueIdentifier()];
	
	for (auto iter = subparts_combination.subparts.begin(); iter != subparts_combination.subparts.end(); ++iter) {
		
		SubpartShapeContextData& subpart_shape_context = candidate_part_shape_context.subparts[iter->subpart_definition.getUniqueIdentifier()];

		// get all possible ly1 parts (i.e. support leafs) from iter->supported_part as root 
		const std::vector<InferenceSupportData*>& ly1_support_leafs = ly1_support_access.getAllInitialLayerSupportParts(*iter->supported_part);
		
		// collect belongie histogram over each support leaf
		for (auto support_iter = ly1_support_leafs.cbegin(); support_iter != ly1_support_leafs.cend(); ++support_iter) {
			const InferenceSupportData& ly1_support = **support_iter;

			cv::Point2i ly1_support_location = ly1_support.support_subpart.getLocation();

			// get belongie histogram
			BelongieHistogramData* leaf_belongie_histogram = belongie_histogram_access.getHistogram(ly1_support.support_subpart, (InferenceTree&)inference_layer.getInferenceTree());

			subpart_shape_context.leafs[ly1_support.support_path].update(leaf_belongie_histogram->histogram, ly1_support_location - ly1_center_location); 
		}
	}
}

void ShapeContextStatistics::ShapeContextData::update(const cv::Mat belongie_histogram, const cv::Point2f part_center) {
	counter++;

	if (counter == 1) {
		this->avg_belongie_data.center_point = part_center;
		this->avg_belongie_data.histogram = belongie_histogram;		
	} else {
		cv::Point2f delta_center = part_center - this->avg_belongie_data.center_point;
		cv::Mat delta_histogram = belongie_histogram - this->avg_belongie_data.histogram;

		float counter_inverse = 1/(float)counter;
		this->avg_belongie_data.center_point += delta_center * counter_inverse;
		this->avg_belongie_data.histogram += delta_histogram * counter_inverse;
	}
}
void ShapeContextStatistics::ShapeContextData::merge(const ShapeContextData& second_shape_context) {
	if (counter == 0) {
		this->counter = second_shape_context.counter;
		this->avg_belongie_data.center_point = second_shape_context.avg_belongie_data.center_point;
		this->avg_belongie_data.histogram = second_shape_context.avg_belongie_data.histogram;
	} else {
		int n = this->counter;
		int m = second_shape_context.counter;
		float fn = (float)n/(m + n);
		float fm = (float)m/(m + n);
		
		this->avg_belongie_data.histogram = this->avg_belongie_data.histogram * fn + second_shape_context.avg_belongie_data.histogram * fm;
		this->avg_belongie_data.center_point = this->avg_belongie_data.center_point * fn + second_shape_context.avg_belongie_data.center_point * fm;
		this->counter = m + n;
	}
}

void ShapeContextStatistics::finalizePartCandidateStatisitcs(std::unordered_map<std::string, VocabularyPartCandidateDefinition>& all_part_candidates_definitions, VocabularyPartCandidateList* vocabulary_part_candidate_list, VocabularyLayer& vocabulary_layer) {

	ShapeContextModifier shape_context_modifier = vocabulary_part_candidate_list->getAccess<ShapeContextModifier>(ShapeContextExtension::Factory());

	// go through all part candidates that we have encountered/created during the update process
	for (auto iter = all_part_candidates_definitions.begin(); iter != all_part_candidates_definitions.end(); ++iter) {
		const std::string& part_candidate_identifyer = iter->first;
		VocabularyPartCandidateDefinition& part_candidate_def = iter->second;

		PartShapeContextData& candidate_part_shape_context = part_candidates_shape_context[part_candidate_identifyer];

		// for each part candidate go over all of its subparts and attach mean belongie histogram and center point to the vocabulary subpart data
		for (auto iter_subpart = part_candidate_def.subparts_definitions.begin(); iter_subpart != part_candidate_def.subparts_definitions.end(); ++iter_subpart) {
			const VocabularySubpartData& created_subpart_data = *(iter_subpart->created_subpart_data);

			SubpartShapeContextData& subpart_shape_context = candidate_part_shape_context.subparts[iter_subpart->getUniqueIdentifier()];

			// copy averaged belongie data (histogram and point) as geometry data
			// TODO: to many coping being done !!!!!
			std::map<PartSupportPath, BelongieHistogramData*> geometry_data;
			for (auto iter_leaf = subpart_shape_context.leafs.begin(); iter_leaf != subpart_shape_context.leafs.end(); ++iter_leaf)
				geometry_data[iter_leaf->first] = new BelongieHistogramData(iter_leaf->second.getAverageBelongieData());
			
			shape_context_modifier.insertGeometry(created_subpart_data, new ShapeContextGeometryData(geometry_data));
		}
	}
}

////////////////////////////////////////////////////////////////////////
////// Similarity based on shape context

void ShapeContextSimilarity::prepareUpdate(const InferenceLayer& inference_layer) {
	shape_context_statistics.prepareUpdate(inference_layer);
	
}

void ShapeContextSimilarity::updatePartCandidateStatisitcs(const VocabularyPartCandidateDefinition& part_candidate, const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& center_part) {
	shape_context_statistics.updatePartCandidateStatisitcs(part_candidate, subparts_combination, inference_layer, center_part);
	
	// TODO: do we need any specific statistics from the shape similarity ?
}

void ShapeContextSimilarity::finalizePartCandidateStatisitcs(std::unordered_map<std::string, VocabularyPartCandidateDefinition>& all_part_candidates_definitions, VocabularyPartCandidateList* vocabulary_part_candidate_list, VocabularyLayer& vocabulary_layer) {
	// need to first insert shape context information
	shape_context_statistics.finalizePartCandidateStatisitcs(all_part_candidates_definitions, vocabulary_part_candidate_list, vocabulary_layer);

	ShapeContextAccess shape_context_access = vocabulary_part_candidate_list->getAccess<ShapeContextAccess>();

	// based on part_learning::merge_stat_sc() 
	// calculate geometry distance between each combination of parts
	// and create cluster of similar parts in iterative manner (greedy)
}


