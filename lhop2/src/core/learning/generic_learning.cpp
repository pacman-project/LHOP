
// generic_learning
///////////////////////////////////////////////////////////////////////////////

#include "generic_learning.h"


void register_generic_learning_classes() {
	/************************************************************************* 
	 *           !!!! DO NOT MODIFY OR REMOVE THIS FUNCTION   !!!!!
	 ************************************************************************/
	ClassRegister::registerFactory<AbstractLayerLearning::IFactory, GenericLayerBatchLearning::Factory>();
}

AbstractLayerLearning* GenericLayerBatchLearning::Factory::newInstance(const IConfiguration& config, std::shared_ptr<VocabularyTree> vocabulary) const {
	return nullptr;
}

void GenericLayerBatchLearning::updateVocabulary(const std::vector<AbstractLayerInputObject*> input_object) {
	// current version does batch learning only (no online learning)
	// so generate parts using following steps:

	// 1. generate parts using co-occurrence statistics
	VocabularyPartCandidateList* vocabulary_candidates = part_generator->obtainVocabularyCandidates(vocabulary->getLayer(layer - 1), input_object);

	// 2. optimize number of parts (using greedy coverage algorithm)
	vocabulary_candidates = part_pruning->pruneVocabularyCandidates(vocabulary_candidates, input_object);

	// 3. perform EM-step optimization 
	// ..

	// finally promote candidate parts into real vocabulary parts
	promoteCandidatesIntoVocabulary(vocabulary_candidates);
}

void GenericLayerBatchLearning::promoteCandidatesIntoVocabulary(VocabularyPartCandidateList* vocabulary_candidates) {
	// TODO:
	// create new vocabulary layer if it does not exists yet

	// promote all vocabulary candidates into real vocabulary parts
/*
       foreach partcand in list
       gen vocab part
       copy uuid and state variables and extensions
*/
}
