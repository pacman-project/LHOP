// shape_check_selection
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_SHAPE_CHECK_SELECTION_
#define _CORE_SHAPE_CHECK_SELECTION_

#include "core/inference/abstract_inference.h"

/// @addtogroup core
/// @{
/// @addtogroup inference
/// @{
/// @addtogroup inferrence_building_blocks
/// @{
/// @addtogroup shape_check
/// @{

////////////////////////////////////////////////////////////////
// Candidate selection per one location implementations

class SelectionPerLocationBlock : public IInferredCandidatesFilteringBlock {
public:
	virtual InferredPartCandidateList* filterInferenceCandidates(InferredPartCandidateList* inference_candidates, const VocabularyTree& vocabulary) ;
};

/// @}
/// @}
/// @}
/// @}

#endif /* _CORE_SHAPE_CHECK_SELECTION_ */
