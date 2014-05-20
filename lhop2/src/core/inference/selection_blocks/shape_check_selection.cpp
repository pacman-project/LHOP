// shape_check_selection
///////////////////////////////////////////////////////////////////////////////

#include "shape_check_selection.h"


InferredPartCandidateList* SelectionPerLocationBlock::filterInferenceCandidates(InferredPartCandidateList* inference_candidates, const VocabularyTree& vocabulary) {

	// first generate "scmap" (previously computed by the get_sc_map)
		// scmap a list of histograms
		// for each first layer location that contains part (feature) we create a histogram from its surrounding parts

	// chech shape of each candidate
	{
		// retrive subparts of this candidate 
		// we require InferenceSubpartsAccess to retrive this information and thefore has to be present in the inference_candidates
		
		// from subparts we are only intereseted in the max subpart found for each different subpart index

		// convert information of subparts into map between TreePath and its shape histogram (path_map_t)

		// we perform the checking process now
		{
			// get geometry information of the vocabulary part from shape context extension (currently geometry_extension.h)
		}


		// from the process we create two new responses S_RESPONSE and X_RESPONSE

	}
	return inference_candidates;
}