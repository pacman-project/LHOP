
// abstract_learning
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_ABSTRACT_LAYER_LEARNING_
#define _CORE_ABSTRACT_LAYER_LEARNING_

#include <vector>
#include <memory>

#include "utils/class_register.h"
#include "utils/configuration.h"

#include "core/input_output/input_object.h"
#include "core/input_output/output_object.h"

#include "core/structures/vocabulary.h"
#include "core/structures/parse_tree.h"

#include "core/structures/responses_array.h"

struct VocabularyPartCandidateList;
struct VocabularyPartCandidate;

/// @addtogroup core
/// @{
/// @addtogroup learning
/// @{

/**
 * Any class providing learning of vocabulary layer parts must extend from this class.
 * The class should "connect" to the vocabulary at specific layer (vocabulary and which 
 * layer should be defined in the constructor) and should implement updateVocabulary(..) 
 * method that updates statistics used to generate vocabulary parts. 
 *
 * The class must perform generation of parts in-place into the existing vocabulary which
 * is also its output.
 * 
 */
class AbstractLayerLearning : public IRegistrableClass {
public:
	/**
	 * Main update should update statistics used to build the vocabulary parts. The output
	 * of this method is the same as 
	 *
	 * Current version still supports batch-mode (using more the one input objects) 
	 * but should be replace with online-mode (updating only one input object at a time)
	 * 
	 * TODO: how to bring different types of input signals (not just groundtruth) of the each input object into the system ??
	 */
	virtual void updateVocabulary(const std::vector<AbstractLayerInputObject*> input_object) = 0;

	/**
	 * Retrieves output vocabulary.
	 */
	virtual VocabularyOutputObject* getOutputVocabulary() = 0;

	class IFactory : public IRegistrableClassFactory {
	public:
		virtual AbstractLayerLearning* newInstance(const IConfiguration& config, std::shared_ptr<VocabularyTree> vocabulary) const = 0;
	};
};

/// @addtogroup learning_building_blocks
/// @{

class IVocabularyCandidatesGeneratorBlock {
public:
	virtual VocabularyPartCandidateList* obtainVocabularyCandidates(VocabularyLayer& current_vocabulary_layer, const std::vector<AbstractLayerInputObject*> input_object) = 0;
};

class IVocabularyCandidatesPruningBlock {
public:
	virtual VocabularyPartCandidateList* pruneVocabularyCandidates(VocabularyPartCandidateList* vocabulary_candidates, const std::vector<AbstractLayerInputObject*> input_object) = 0;

};

/// @}


//////////////////////////////////////////////////////
// Multiple Layer Inference implementation

class MultipleLayersBatchLearning : public AbstractLayerLearning {
private:
	std::vector<AbstractLayerLearning*> layer_learning;
public:
	MultipleLayersBatchLearning(std::vector<AbstractLayerLearning*> layer_learning) : layer_learning(layer_learning) {
	}

	virtual void updateVocabulary(const std::vector<AbstractLayerInputObject*> input_object);

	/**
	 * Retrieves output vocabulary from the last layer
	 */
	virtual VocabularyOutputObject* getOutputVocabulary() {
		return layer_learning.back()->getOutputVocabulary();
	}

	class Factory : public IFactory {
	public:
		virtual AbstractLayerLearning* newInstance(const IConfiguration& config, std::shared_ptr<VocabularyTree> vocabulary) const;

		virtual std::string assignedRegistrationName() const {
			return "learning.multiple-layer";
		}
	};
};


//////////////////////////////////////////////////////
/// Structures used during the learning (for candidates)



struct VocabularyPartCandidateList : public ExtensionHolder {
	// extends from the ExtensionHolder so we can attach any extension for parts still in construction
	// (i.e. any extension using VocabularyPart will now attach itself to VocabularyPartCandidate until they are promoted into the VocabularyPart)
	std::vector<VocabularyPartCandidate*> candidates;

	// default constructor
	VocabularyPartCandidateList() {}
};

struct VocabularyPartCandidate : public IAttachableClass {
	bool invalid;

	// default constructor
	VocabularyPartCandidate(UUIDType uuid) :
		IAttachableClass(uuid), invalid(false) {
	}
};



/// @}
/// @}


#endif /* _CORE_ABSTRACT_LAYER_LEARNING_ */
