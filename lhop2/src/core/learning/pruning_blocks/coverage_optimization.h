
// coverage_optimization
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_LEARNING_COVERAGE_OPTIMIZATION_
#define _CORE_LEARNING_COVERAGE_OPTIMIZATION_

#include "core/learning/abstract_learning.h"

/// @addtogroup core
/// @{
/// @addtogroup learning
/// @{
/// @addtogroup learning_building_blocks
/// @{
/// @addtogroup learning_building_blocks_coverage
/// @{

class CoverageOptimizationMCMC : public IVocabularyCandidatesPruningBlock {
public:
	virtual VocabularyPartCandidateList* pruneVocabularyCandidates(VocabularyPartCandidateList* vocabulary_candidates, const std::vector<AbstractLayerInputObject*> input_object);
};


/// @}
/// @}
/// @}
/// @}


#endif /* _CORE_LEARNING_COVERAGE_OPTIMIZATION_ */
