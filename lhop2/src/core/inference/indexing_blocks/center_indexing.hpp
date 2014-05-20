// indexing_block/center_indexing
///////////////////////////////////////////////////////////////////////////////

#include "core/inference/indexing_blocks/center_indexing.h"

#include "core/structures/indexing/indexing_extension.h"

bool CandidateThresholdingBlock::isCandidateCenterPartValid(const InferredPart& center_part) {
	// verify that the responses of this part are not below candidate_r(g)_threshold (see above)
	// otherwise this part is not considered to be a valid center
	return center_part.getResponse(ResponseType::R_RESPONSE_1) >= candidate_r_threshold 
		&& center_part.getResponse(ResponseType::G_RESPONSE_1) >= candidate_g_threshold ? true : false;			
}

bool CandidateThresholdingBlock::isCandidateIndexingPartValid(const VocabularyPart& indexing_part) {
	// make sure indexing part is in the allowed part list
	return allowed_parts.empty() || allowed_parts.find(indexing_part.getTypeId()) != allowed_parts.end() ? true : false;
}

template <class T>
InferredPartCandidateList* IndexingBlock<T>::obtainInferenceCandidates(InferenceLayer& current_layer_inference, const VocabularyTree& vocabulary) {
	InferredPartCandidateList* candidates_list = new InferredPartCandidateList(current_layer_inference);
	
	// get reference (has to be pointer/reference) to the generator for inferred part uuids
	UUIDGenerator& inferred_part_uuid_gen = current_layer_inference.getInferenceTree().getInferredPartUUIDGenerator();

	// get access to indexing extension
	VocabularyIndexingAccess indexing_access = vocabulary.getAccess<VocabularyIndexingAccess>();

	auto& candidates = candidates_list->candidates;

	// go over all inferred part of the current layer
	for (auto loc = current_layer_inference.beginIterator(); loc != current_layer_inference.endIterator(); ++loc) {
		std::vector<InferredPart*>& parts_at = current_layer_inference.getPartsAt(loc);

		for (auto iter = parts_at.begin(); iter != parts_at.end(); ++iter) {
			InferredPart& center_candidate_part = **iter; // *iter is InferredPart*

			// verify that candidate center is valid (i.e. is above certain threshold)
			if (candidate_thresholding.isCandidateCenterPartValid(center_candidate_part) == false)
				continue;

			// get corresponding vocabulary part
			const VocabularyPart& vocabulary_part = center_candidate_part.getCorrespondingVocabularyPart(vocabulary);

			// find all central parts indexed to the vocabulary part (use VocabularyIndexingAccess to find this information)
			std::vector<VocabularyPart*> indexed_parts = indexing_access.getIndexedCentralParts(vocabulary_part);

			// for each indexed vocabulary part create one new candidate
			for (auto indexed_iter = indexed_parts.begin(); indexed_iter != indexed_parts.end(); ++indexed_iter) {
				VocabularyPart& indexed_part = **indexed_iter;	// *indexed_iter is VocabularyPart*		

				// make sure indexing part is valid 
				if (candidate_thresholding.isCandidateIndexingPartValid(indexed_part) == false)
					continue;

				// save InferredPartCandidate with original candidate and indexed candidate into the list of candidates
				candidates.push_back(InferredPartCandidate(inferred_part_uuid_gen.generateUUID(), &center_candidate_part, &indexed_part));
			}
		}
	}

	return candidates_list;
}
