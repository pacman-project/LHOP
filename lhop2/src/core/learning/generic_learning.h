
// generic_learning
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_GENERIC_LAYER_LEARNING_
#define _CORE_GENERIC_LAYER_LEARNING_

#include "core/learning/abstract_learning.h"

/// @addtogroup core
/// @{
/// @addtogroup learning
/// @{

/**
 * Batch-mode learning for each generic layer
 */
class GenericLayerBatchLearning : public AbstractLayerLearning {
protected:
	const int layer;
	std::shared_ptr<VocabularyTree> vocabulary;

	IVocabularyCandidatesGeneratorBlock* part_generator;
	IVocabularyCandidatesPruningBlock* part_pruning;
public:
	GenericLayerBatchLearning(std::shared_ptr<VocabularyTree> vocabulary, int layer, 
							IVocabularyCandidatesGeneratorBlock* part_generator,
							IVocabularyCandidatesPruningBlock* part_pruning) 
		: vocabulary(vocabulary), layer(layer), part_generator(part_generator), part_pruning(part_pruning) {
	}

	void updateVocabulary(const std::vector<AbstractLayerInputObject*> input_object);

	virtual VocabularyOutputObject* getOutputVocabulary() {
		return new VocabularyOutputObject(vocabulary);
	}

	class Factory : public IFactory {
	public:
		virtual AbstractLayerLearning* newInstance(const IConfiguration& config, std::shared_ptr<VocabularyTree> vocabulary) const;

		virtual std::string assignedRegistrationName() const {
			return "learning.generic-layer";
		}
	};

protected:
	virtual void promoteCandidatesIntoVocabulary(VocabularyPartCandidateList* vocabulary_candidates);
};



/// @}
/// @}

#endif /* _CORE_ABSTRACT_LAYER_LEARNING_ */
