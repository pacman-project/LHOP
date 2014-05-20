
// abstract_learning
///////////////////////////////////////////////////////////////////////////////

#include "abstract_learning.h"


void register_abstract_learning_classes() {
	/************************************************************************* 
	 *           !!!! DO NOT MODIFY OR REMOVE THIS FUNCTION   !!!!!
	 ************************************************************************/
	ClassRegister::registerFactory<AbstractLayerLearning::IFactory, MultipleLayersBatchLearning::Factory>();
}

AbstractLayerLearning* MultipleLayersBatchLearning::Factory::newInstance(const IConfiguration& config, std::shared_ptr<VocabularyTree> vocabulary) const {
	// TODO: generate multiple layer learning from configuration
	return nullptr;
}