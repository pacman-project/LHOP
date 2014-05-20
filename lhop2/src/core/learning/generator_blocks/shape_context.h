
// shape_context
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_LEARNING_SHAPE_CONTEXT_
#define _CORE_LEARNING_SHAPE_CONTEXT_

#include <queue>
#include <stack>
#include <functional>
#include <unordered_set>
#include <unordered_map>

#include "core/learning/abstract_learning.h"

#include "core/learning/generator_blocks/subpart_based_generator.h"

#include "core/structures/subparts/subparts_extension.h"
#include "core/structures/subparts/support_extension.h"
#include "core/structures/shape_context/shape_context_extension.h"


/// @addtogroup core
/// @{
/// @addtogroup learning
/// @{
/// @addtogroup learning_building_blocks
/// @{
/// @addtogroup learning_building_blocks_shape_context
/// @{

////////////////////////////////////////////////////////////////////////
////// Classes for shape context 
/**
 * TODO: write detail description
 */
class ShapeContextStatistics : public IPartCandidateStatisticsCollector {
	//BelongieHistogramAccess* belongie_histogram_access;
	class ShapeContextData {
		BelongieHistogramData avg_belongie_data;
		int counter;
	public:
		ShapeContextData() : counter(0), avg_belongie_data(cv::Mat(), cv::Point()) {}
		void update(const cv::Mat belongie_histogram, const cv::Point2f part_center);
		void merge(const ShapeContextData& second_shape_context);
		BelongieHistogramData& getAverageBelongieData() { return avg_belongie_data; }
	};	
	struct SubpartShapeContextData {
		// maps unique_id(path between root and ly1 support leaf) -> ShapeContextData
		std::map<PartSupportPath, ShapeContextData> leafs;
	};	
	struct PartShapeContextData {
		// maps unique_id(VocabularySubpartDefinition) -> SubpartShapeContextData
		std::unordered_map<std::string, SubpartShapeContextData> subparts;
	};
	// maps unique_id(VocabularyPartCandidateDefinition) -> PartShapeContextData
	std::unordered_map<std::string,PartShapeContextData> part_candidates_shape_context;
public:
	virtual void prepareUpdate(const InferenceLayer& inference_layer);

	virtual void updatePartCandidateStatisitcs(const VocabularyPartCandidateDefinition& part_candidate, const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& center_part);

	virtual void finalizePartCandidateStatisitcs(std::unordered_map<std::string, VocabularyPartCandidateDefinition>& all_part_candidates_definitions, VocabularyPartCandidateList* vocabulary_part_candidate_list, VocabularyLayer& vocabulary_layer);

private:
	InferredPart* followCentralPartToLy1(const InferredPart& root, const InferenceLayer& inference_layer);
};


/**
 * Class for computing similarity between parts based on shape context statistics. 
 * It relies upon ShapeContextStatistics and first relays all calls to this class.
 * In finalizePartCandidateStatisitcs it computes geometry distance between all
 * combination of parts and creates cluster of similar parts.
 */
class ShapeContextSimilarity : public IPartCandidateStatisticsCollector {
	ShapeContextStatistics shape_context_statistics;
public:

	virtual void prepareUpdate(const InferenceLayer& inference_layer);

	virtual void updatePartCandidateStatisitcs(const VocabularyPartCandidateDefinition& part_candidate, const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& center_part);

	virtual void finalizePartCandidateStatisitcs(std::unordered_map<std::string, VocabularyPartCandidateDefinition>& all_part_candidates_definitions, VocabularyPartCandidateList* vocabulary_part_candidate_list, VocabularyLayer& vocabulary_layer);

};

////////////////////////////////////////////////////////////////////////
////// Classes for PCA

/**
 * TODO: write detail description
 */
/*class PCAStatistics : public IPartCandidateStatisticsCollector {
public:
	virtual void prepareUpdate(const InferenceLayer& inference_layer) {
		// nothing to do here
	}	

	virtual void updatePartCandidateStatisitcs(const VocabularyPartCandidateDefinition& part_candidate, const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& center_part) {
		// TODO: implement
	}

	virtual void finalizePartCandidateStatisitcs(std::unordered_map<std::string, VocabularyPartCandidateDefinition>& all_part_candidates_definitions, VocabularyPartCandidateList* vocabulary_part_candidate_list, VocabularyLayer& vocabulary_layer) {
		// TODO: implement
	}
};

class PCASimilarity : public IPartCandidateStatisticsCollector {
	PCAStatistics pca_statistics;
public:
	virtual void prepareUpdate(const InferenceLayer& inference_layer) {
		pca_statistics.prepareUpdate(inference_layer);
	}	

	virtual void updatePartCandidateStatisitcs(const VocabularyPartCandidateDefinition& part_candidate, const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& center_part) {
		pca_statistics.updatePartCandidateStatisitcs(part_candidate, subparts_combination, inference_layer, center_part);
		// TODO: implement
	}

	virtual void finalizePartCandidateStatisitcs(std::unordered_map<std::string, VocabularyPartCandidateDefinition>& all_part_candidates_definitions, VocabularyPartCandidateList* vocabulary_part_candidate_list, VocabularyLayer& vocabulary_layer) {
		pca_statistics.finalizePartCandidateStatisitcs(all_part_candidates_definitions, vocabulary_part_candidate_list, vocabulary_layer);
		// TODO: implement
	}
};
*/


/// @}
/// @}
/// @}
/// @}


#endif /* _CORE_LEARNING_SHAPE_CONTEXT_ */
