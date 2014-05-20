// generic_layer_inference
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_GENERIC_LAYER_INFERENCE_
#define _CORE_GENERIC_LAYER_INFERENCE_

#include "core/inference/abstract_inference.h"

/// @addtogroup core
/// @{
/// @addtogroup inference
/// @{

/**
 * Generic layer building (inference) is constructed from three main building blocks:
 *  - indexing - available implementations: IndexingBlock
 *  - part matching - available implementations: SchurProductMatchingBlock (which is also constructed from a few other building blocks)
 *  - part selection per location - available implementations: /
 *
 * Main method is performInference(...) which takes abstract layer input object (existing InferenceTree) and returns 
 * layer output object with one added layer. Layer is added based on the provided vocabulary of parts (VocabularyTree) 
 * which must be supplied through the constructor.
 *
 * The process works by first obtaining the list of candidates through part indexing block (IIndexingBlock). All the 
 * candidates are then filtered with the matching block and finally filtered with the part selection per location block.  
 *
 * Final promotion of candidate parts into the layer parts is done with the (protected) promoteCandidatesIntoLayer(...) method which 
 * will create InferredPart from CandidatePart using its location, calculated responses and any additionally created 
 * extension attached to the parts. The individual building blocks can attach extension directly to the InferredPartCandidate and
 * the promoteCandidatesIntoLayer(...) method will transfer that information directly to the new InferredPart.
 */
class GenericLayersInference : public AbstractLayerInference {
private:
	// building blocks
	IInferredCandidatesIndexingBlock* part_indexing;
	IInferredCandidatesFilteringBlock* part_matching;
	IInferredCandidatesFilteringBlock* part_selection_per_location;

	const int layer;
	const double layer_contraction;
	
	const std::shared_ptr<VocabularyTree> vocabulary;
public:

	GenericLayersInference(IInferredCandidatesIndexingBlock* part_indexing, 
							IInferredCandidatesFilteringBlock* part_matching, 
							IInferredCandidatesFilteringBlock* part_selection_per_location,
							int layer, double layer_contraction, 
							std::shared_ptr<VocabularyTree> vocabulary)
		: part_indexing(part_indexing), part_matching(part_matching), part_selection_per_location(part_selection_per_location),
		  layer(layer), layer_contraction(layer_contraction), vocabulary(vocabulary) {
	}

	virtual ~GenericLayersInference() {
		if (part_indexing != nullptr) delete part_indexing;
		if (part_matching != nullptr) delete part_matching;
		if (part_selection_per_location != nullptr) delete part_selection_per_location;
	}

	virtual LayerOutputObject* performInference(const AbstractLayerInputObject* input_object);

	class Factory : public IFactory {
	 public:
		 /**
		  * Main factory for the GenericLayersInference objects. It reads all the settings from 
		  * the provided configuration object.
		  */
		 virtual AbstractLayerInference* newInstance(const IConfiguration& config, const std::shared_ptr<VocabularyTree> vocabulary) const;

		 /**
		  * Registration name of this inference. It is used to identify this inference when requested 
		  * by the user in the IConfiguration.
		  */
		 virtual std::string assignedRegistrationName() const {
			return "generic";
		}
	};
};

//////////////////////////////////////////////////
// Misc building blocks

/// @addtogroup building_blocks
/// @{


class EliminateForbiddenBlock : public IInferredCandidatesFilteringBlock {
public:
	virtual InferredPartCandidateList* filterInferenceCandidates(InferredPartCandidateList* inference_candidates, const VocabularyTree& vocabulary) {}
};

/// @}

/// @}
/// @}


#endif /* _CORE_GENERIC_LAYER_INFERENCE_ */
