// indexing_block/center_indexing
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_CENTER_INDEXING_
#define _CORE_CENTER_INDEXING_

#include "core/inference/abstract_inference.h"

/// @addtogroup core
/// @{
/// @addtogroup inference
/// @{
/// @addtogroup inferrence_building_blocks 
/// @{
/// @addtogroup indexing
/// @{

// forward definitions
class ICandidateThresholdingBlock;
class CandidateThresholdingBlock;
class NullCandidateThresholdingBlock;


//////////////////////////////////////////////////////
/// Indexing and initial candidate thresholding implementations
/// @{

/**
 * Interface for any candidate thresholding block used by the IndexingBlock. 
 */ 
class ICandidateThresholdingBlock {
public:
	virtual bool isCandidateCenterPartValid(const InferredPart& center_part) = 0;
	virtual bool isCandidateIndexingPartValid(const VocabularyPart& indexing_part) = 0;
};

/*
 * Default candidate thresholding:
 *  - isCandidateCenterPartValid: allows central parts with R_RESPONSE and G_RESPONSE over certain threshold
 *  - isCandidateIndexingPartValid: allows indexing parts of specific type
 */
class CandidateThresholdingBlock : public ICandidateThresholdingBlock {
	double candidate_r_threshold;
	double candidate_g_threshold;
	set<int> allowed_parts;
public:
	CandidateThresholdingBlock(double candidate_r_threshold,
								double candidate_g_threshold,
								set<int> allowed_parts)
		: candidate_r_threshold(candidate_r_threshold), candidate_g_threshold(candidate_g_threshold), allowed_parts(allowed_parts) {
	}

	bool isCandidateCenterPartValid(const InferredPart& center_part);
	bool isCandidateIndexingPartValid(const VocabularyPart& indexing_part);
};

/**
 * Empty thresholding implementation of the ICandidateThresholdingBlock that
 * does not do any thresholding and always returns true.
 */
class NullCandidateThresholdingBlock : public ICandidateThresholdingBlock {
public:
	bool isCandidateCenterPartValid(const InferredPart& center_part) { return true; }
	bool isCandidateIndexingPartValid(const VocabularyPart& indexing_part)  { return true; }
};

/// @}

/**
 * Indexing block that creates initial list of candidates by finding corresponding vocabualry parts
 * and creating one candidate for each vocabulary part that is indexed as using current part for its center.
 * 
 * Can also use candidate thresholding by supplying any ICandidateThresholdingBlock implementation as template
 * parameter. By default it uses CandidateThresholdingBlock class.
 */
template <class TCandidateThresholdingBlock = CandidateThresholdingBlock>
class IndexingBlock : public IInferredCandidatesIndexingBlock {
	TCandidateThresholdingBlock candidate_thresholding;
public:
	IndexingBlock(TCandidateThresholdingBlock& candidate_thresholding) : candidate_thresholding(candidate_thresholding) {
	}

	virtual InferredPartCandidateList* obtainInferenceCandidates(InferenceLayer& current_layer_inference, const VocabularyTree& vocabulary);
};

// include template definitions
#include "core/inference/indexing_blocks/center_indexing.hpp"

/// @}
/// @}
/// @}
/// @}

#endif /* _CORE_CENTER_INDEXING_ */
