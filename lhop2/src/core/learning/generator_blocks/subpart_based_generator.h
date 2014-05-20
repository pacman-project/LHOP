
// subpart_based_generator
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_LEARNING_COOCCURRENCE_BASED_GEN_
#define _CORE_LEARNING_COOCCURRENCE_BASED_GEN_

#include <queue>
#include <stack>
#include <functional>
#include <unordered_set>
#include <unordered_map>

#include "core/learning/abstract_learning.h"

#include "core/learning/generator_blocks/cooccurrence_statistics.h"

#include "core/structures/subparts/subparts_extension.h"


/// @addtogroup core
/// @{
/// @addtogroup learning
/// @{
/// @addtogroup learning_building_blocks
/// @{
/// @addtogroup learning_building_blocks_cooccurrence
/// @{

struct VocabularySubpartSupport;
struct VocabularySubpartsCombination;

/**
 * Simple definition class for one vocabulary subpart composed from:
 *  - uuid of the surrounding vocabulary part (vocabulary_part_uuid)
 *  - pointer to the local maxima that was used to create this subpart definition
 *
 * Has overloaded less operator (to enable sorting and equal operations). Two objects
 * are equal if vocabulary_part_uuid and local maxima offset are identical.
 */
struct VocabularySubpartDefinition {
	UUIDType vocabulary_subpart_uuid;
	CooccurrenceStatistics::LocalMaximaDistribution* local_maxima;

	// used only after the promoteVocabularyPartCandidateDefinitions(..)
	// so other IPartCandidateStatisticsCollector can find corresponding
	// data during the finalizePartCandidateStatisitcs(..)
	VocabularySubpartData* created_subpart_data;

	VocabularySubpartDefinition() : vocabulary_subpart_uuid(), local_maxima(nullptr) {}
	VocabularySubpartDefinition(const UUIDType& vocabulary_subpart_uuid, CooccurrenceStatistics::LocalMaximaDistribution* local_maxima)
		: vocabulary_subpart_uuid(vocabulary_subpart_uuid), local_maxima(local_maxima) {
	}

	bool operator<(const VocabularySubpartDefinition& obj) const {
		if (vocabulary_subpart_uuid != obj.vocabulary_subpart_uuid)
			return vocabulary_subpart_uuid < obj.vocabulary_subpart_uuid; // sort by vocabulary part UUID first

		const cv::Point2i a_loc = local_maxima->offset;
		const cv::Point2i b_loc = obj.local_maxima->offset;

		// then sort by offset location
		if (a_loc.x != b_loc.x)
			return a_loc.x < b_loc.x; // first by x location
		else
			return a_loc.y < b_loc.y; // then by y location
	}

	std::string getUniqueIdentifier() const;

};


/**
 * An definition for specific VocabularyPartCandidate instance. This class
 * uses a vector of VocabularySubpartDefinition to uniquely identify each
 * subparts combination. For each VocabularySubpartsCombination one
 * corresponding VocabularyPartCandidateDefinition is created during update process.
 */
struct VocabularyPartCandidateDefinition {
	std::string unique_identifier;
	UUIDType central_part_uuid;

	std::vector<VocabularySubpartDefinition> subparts_definitions;

	// used only after the promoteVocabularyPartCandidateDefinitions(..)
	// so other IPartCandidateStatisticsCollector can find corresponding
	// part candidate during the finalizePartCandidateStatisitcs(..)
	VocabularyPartCandidate* created_vocabulary_part_candidate;

	VocabularyPartCandidateDefinition() {}
	VocabularyPartCandidateDefinition(UUIDType central_part_uuid, const VocabularySubpartsCombination& subparts_combination);

	/**
	 * Unique identifier for specific combination of subparts. Combination is
	 * defined by a set of subparts definition converted into string that relies on
	 * the sored order of subpart definitions to correctly identify combination
	 * of subparts. A tuple (vocabulary part UUID and local maxima offset) is
	 * used as identifier of each subpart definition.
	 */
	std::string getUniqueIdentifier() const { return unique_identifier; }
};


/**
 * TODO: write detail description
 */
class IPartCandidateStatisticsCollector {
public:
	virtual void prepareUpdate(const InferenceLayer& inference_layer) = 0;

	virtual void updatePartCandidateStatisitcs(const VocabularyPartCandidateDefinition& part_candidate, const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& center_part) = 0;

	virtual void finalizePartCandidateStatisitcs(std::unordered_map<std::string, VocabularyPartCandidateDefinition>& all_part_candidates_definitions, VocabularyPartCandidateList* vocabulary_part_candidate_list, VocabularyLayer& vocabulary_layer) = 0;
	
};

/**
 * TODO: write detail description
 */
class SubpartBasedCandidatesGenerator : public IVocabularyCandidatesGeneratorBlock {
private:
	// parameters:
	const double center_part_threshold;
	const int min_subpart_count;
	const int max_subpart_count;
	const int max_new_subpart_hypothesis_added;

	const double min_intersection_threhsold;
	const double max_intersection_threhsold;

	// building blocks:
	CooccurrenceStatistics* cooccurrence_statistics;
	std::vector<IPartCandidateStatisticsCollector*> individual_statistics_collectors;

	// list of all observed new vocabulary part definitions
	std::unordered_map<std::string, VocabularyPartCandidateDefinition> all_part_candidates_definitions;
public:
	SubpartBasedCandidatesGenerator(CooccurrenceStatistics* cooccurrence_statistics, 
										std::vector<IPartCandidateStatisticsCollector*> individual_statistics_collectors,
										double center_part_threshold,
										int min_subpart_count, int max_subpart_count,
										int max_new_subpart_hypothesis_added,
										double min_intersection_threhsold,
										double max_intersection_threhsold) 
		: cooccurrence_statistics(cooccurrence_statistics), individual_statistics_collectors(individual_statistics_collectors),
		  center_part_threshold(center_part_threshold), 
		  min_subpart_count(min_subpart_count), max_subpart_count(max_subpart_count),
		  max_new_subpart_hypothesis_added(max_new_subpart_hypothesis_added),
		  min_intersection_threhsold(min_intersection_threhsold),
		  max_intersection_threhsold(max_intersection_threhsold) {
	}

	VocabularyPartCandidateList* obtainVocabularyCandidates(VocabularyLayer& current_vocabulary_layer, const std::vector<AbstractLayerInputObject*> input_object);

protected:

	virtual void updateCandidatePartStatistics(const InferenceLayer& inference_tree);

	/**
	 * Retrieve a sorted list of all the supported subparts (definitions) around inferred part center_part using the following procedure:
	 *  - use CooccurrenceStatistics to retrieve a list of all possible subparts around specific vocabulary part (of the same type as center_part)
	 *  - match each possible subpart that statistics says it is frequent around the center_part to verify that the sample (inference_layer) supports subpart
	 *  - all possible matched/supported subparts are sorted according to its matched/supported value and returned as vector of VocabularySubpartSupport
	 */
	virtual std::vector<VocabularySubpartSupport> getSupportedSubpartsSortedBySupportValue(const InferenceLayer& inference_layer, const InferredPart& center_part);
	
	/**
	 * Verification method to check if specific subpart (definition( is supported by the sample (inference_layer) around specific center_part.
	 * Matching is done as simple schur product. We obtain offset and distribution matrix from CooccurrenceStatistics for this particular subpart
	 * (relative to particular central vocabulary subpart) and apply it to area around center_part to find the matching inferred part in sample.
	 * We return only best part found (one with max dist(i,j) * part.G_RESPONSE) and return it as VocabularySubpartSupport with supported_part
	 * and supported_value.
	 */
	virtual VocabularySubpartSupport verifySubpartSupport(const VocabularySubpartDefinition& subpart_definition, const InferenceLayer& inference_layer, const InferredPart& center_part);

	/**
	 * Updates any statistics for the specific combination of subparts (combination of subpart definitions) that is supported by the sample at the center_part.
	 * An unique ID is assigned to the subparts combination and then statistics gathering is delegated to any implementation that wishes to collect any kind 
	 * of statistics that can be found within this center_part for this specific subparts combination.
	 */
	virtual void updateSubpartCombinationStatisitcs(const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& inferred_center_part);

	/**
	 * Promote all subparts combinations that have good enough statistics into vocabulary part candidates. 
	 * Input parameter VocabularyTree is required only to get UUID for future VocabularyPart (i.e. for VocabularyPartCandidate)
	 */
	virtual VocabularyPartCandidateList* promoteVocabularyPartCandidateDefinitions(VocabularyLayer& vocabulary_layer);
};


/**
 * Simple class that holds support for one subpart definition (subpart found around specific center part) composed from:
 *  - subpart definition (VocabularySubpartDefinition)
 *  - pointer to the inferred part in the sample that supports this subpart definition
 *  - value of the support from the inferred part that supports this subpart definition
 * 
 * Has two possible comparer classes:	 
 *  - ComparerBySupportedValue: uses support value as compare value
 *  - ComparerBySubpartDefinition: uses comparer of the VocabularySubpartDefinition
 */			
struct VocabularySubpartSupport {	
	VocabularySubpartDefinition subpart_definition;

	InferredPart* supported_part;
	double supported_value;

	VocabularySubpartSupport() : subpart_definition(), supported_part(nullptr), supported_value(-1) {}
	VocabularySubpartSupport(const VocabularySubpartDefinition& subpart_definition) : subpart_definition(subpart_definition), supported_part(nullptr), supported_value(-1) {}
	VocabularySubpartSupport(const VocabularySubpartDefinition& subpart_definition, InferredPart* supported_part, double supported_value) 
		: subpart_definition(subpart_definition), supported_part(supported_part), supported_value(supported_value) {
	}

	struct ComparerBySupportedValue {
		bool operator()(const VocabularySubpartSupport& a, const VocabularySubpartSupport& b ) { return a.supported_value < b.supported_value; }
	};
	struct ComparerBySubpartDefinition {
		bool operator()(const VocabularySubpartSupport& a, const VocabularySubpartSupport& b ) const { return a.subpart_definition < b.subpart_definition; }
	};
};

/**
 * Simple class with a set of supported subparts and their definitions. Each instance of this class is basis for
 * one new vocabulary part candidate. Theirs subpart definitions represent candidates for subparts of the new vocabulary part.
 */
struct VocabularySubpartsCombination {
	// use ComparerBySubpartDefinition comparer to identify two identical subpart hypothesis 
	// WARNING: it is important to maintain the correct order of subpart in the std::set
	// since the order of subparts will be used to generate unique ID in the VocabularyPartCandidateDefinition
	std::set<VocabularySubpartSupport, VocabularySubpartSupport::ComparerBySubpartDefinition> subparts;

	// cache ly1_support to avoid excessive computation of 
	std::set<InferredPart*> ly1_support;

	VocabularySubpartsCombination() {}
	VocabularySubpartsCombination(const VocabularySubpartsCombination& current_combination, const VocabularySubpartSupport& added_subpart, const std::set<InferredPart*>& added_subpart_ly1_support) {
		subparts = current_combination.subparts; // make a copy
		subparts.insert(added_subpart);

		ly1_support = current_combination.ly1_support;
		ly1_support.insert(added_subpart_ly1_support.begin(),added_subpart_ly1_support.end());
	}
};

/**
 * Method updatePartCandidateStatisitcs:
 * For each part_candidate (identified by the same UUID) we simply multiply
 * values from its support part (supported_value) and collected them in part_candidate_counter
 * by summing them together.
 *
 * Method finalizePartCandidateStatisitcs:
 * Insert subpart data into VocabularySubpartsFunctionality for each subpart definition.
 */
class SubpartPartCandidateStatisticsCollector : public IPartCandidateStatisticsCollector {
	std::unordered_map<std::string, float> part_candidate_counter;
public:
	virtual void prepareUpdate(const InferenceLayer& inference_layer);

	virtual void updatePartCandidateStatisitcs(const VocabularyPartCandidateDefinition& part_candidate, const VocabularySubpartsCombination& subparts_combination, const InferenceLayer& inference_layer, const InferredPart& center_part);

	virtual void finalizePartCandidateStatisitcs(std::unordered_map<std::string, VocabularyPartCandidateDefinition>& all_part_candidates_definitions, VocabularyPartCandidateList* vocabulary_part_candidate_list, VocabularyLayer& vocabulary_layer);
};

/// @}
/// @}
/// @}
/// @}


#endif /* _CORE_LEARNING_COOCCURRENCE_BASED_GEN_ */
