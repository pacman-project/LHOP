
// abstract_layer_inference
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_ABSTRACT_LAYER_INFERENCE_
#define _CORE_ABSTRACT_LAYER_INFERENCE_

#include <functional>

#include "utils/class_register.h"
#include "utils/configuration.h"

#include "core/input_output/input_object.h"
#include "core/input_output/output_object.h"

#include "core/structures/vocabulary.h"
#include "core/structures/parse_tree.h"

#include "core/structures/responses_array.h"

/// @addtogroup core
/// @{
/// @addtogroup inference
/// @{


// forward definitions
struct InferredPartCandidateList;
struct InferredPartCandidate;

//////////////////////////////////////////////////////
/// Abstract inference and building block interfaces 

class AbstractLayerInference : public IRegistrableClass {
public:
	virtual LayerOutputObject* performInference(const AbstractLayerInputObject* input_object) = 0;

	class IFactory : public IRegistrableClassFactory {
	 public:
		virtual AbstractLayerInference* newInstance(const IConfiguration& config, const std::shared_ptr<VocabularyTree> vocabulary) const = 0;		
	};
protected:
	
	/**
	 * Creates new layer for the inference tree and promotes all candidates from the candidates_list->candidates 
	 * that are still alive (i.e. have invalid == false) into layer parts.
	 * 
	 * InferenceLayer* is created based on size of the previous layer (candidates_list->current_layer_inference)
	 * but is contracted by the layer_contraction factor (size is divided by the layer_contraction factor).
	 * 
	 * Part is created at the location of its center part and offseted by the center of mass of the corresponding 
	 * vocabulary part. Location is additionally contracted (divided) by the layer_contraction factor.
	 * All the functionalities attached to the candidate are moved to the InferenceTree (using InferredPartCandidateList.moveExtensionsTo()).
	 * Any references to  the deleted (not promoted) candidates are first removed (using InferredPartCandidateList.deleteAttachedReferences<>()) 
	 * before they are transferred to the InferenceTree.
	 */
	InferenceLayer* promoteCandidatesIntoLayer(InferredPartCandidateList* candidates, int layer, double layer_contraction);
};

/// @addtogroup inferrence_building_blocks
/// @{

class IInferredCandidatesIndexingBlock {
public:
	virtual InferredPartCandidateList* obtainInferenceCandidates(InferenceLayer& current_layer_inference, const VocabularyTree& vocabulary) = 0;
};

class IInferredCandidatesFilteringBlock {
public:
	virtual InferredPartCandidateList* filterInferenceCandidates(InferredPartCandidateList* inference_candidates, const VocabularyTree& vocabulary) = 0;

};

/// @}

//////////////////////////////////////////////////////
// Multiple Layer Inference implementation

class MultipleLayersInference : public AbstractLayerInference {
private:
	std::vector<AbstractLayerInference*> layer_inferences;
public:
	MultipleLayersInference(std::vector<AbstractLayerInference*> layer_inferences) : layer_inferences(layer_inferences) {
	}

	virtual LayerOutputObject* performInference(const AbstractLayerInputObject* input_object);

	class Factory : public IFactory {
	 public:
		virtual AbstractLayerInference* newInstance(const IConfiguration& config, const std::shared_ptr<VocabularyTree> vocabulary) const;

		virtual std::string assignedRegistrationName() const {
			return "inference.multiple-layer";
		}
	};
};

//////////////////////////////////////////////////////
/// Structures used during the inference (for candidates)

struct InferredPartCandidateList : public ExtensionHolder {
	// extends from the ExtensionHolder so we can attach any extension for parts still in construction
	// (i.e. any extension using InferredPart will now attach itself to InferredPartCandidate until they are promoted into the InferredPart)

	InferenceLayer& current_layer_inference;
	std::vector<InferredPartCandidate> candidates;

	// default constructor
	InferredPartCandidateList(InferenceLayer& current_layer_inference) : current_layer_inference(current_layer_inference) {}
};

struct InferredPartCandidate : public IAttachableClass {
	InferredPart* current_inferred_part;
	VocabularyPart* candidate_vocabulary_part;

	ResponsesArray calculated_responses;

	bool invalid;

	// default constructor
	InferredPartCandidate(UUIDType uuid, InferredPart* current_inferred_part, VocabularyPart* candidate_vocabulary_part) :
		IAttachableClass(uuid), current_inferred_part(current_inferred_part), candidate_vocabulary_part(candidate_vocabulary_part), invalid(false) {
	}
};


/// @}
/// @}

#endif /* _CORE_ABSTRACT_LAYER_INFERENCE_ */
