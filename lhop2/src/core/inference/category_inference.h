
//  category_layer_inference
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _CORE_CATEGORY_LAYER_INFERENCE_
#define _CORE_CATEGORY_LAYER_INFERENCE_

#include "core/inference/abstract_inference.h"

/// @addtogroup core
/// @{
/// @addtogroup inference
/// @{


class CategoryLayerInference : public AbstractLayerInference {

public:

	virtual LayerOutputObject* performInference(const AbstractLayerInputObject* input_object);

	class Factory : public IFactory {
	 public:
		virtual AbstractLayerInference* newInstance(const IConfiguration& config, const std::shared_ptr<VocabularyTree> vocabulary) const {
			return new CategoryLayerInference();
		}

		/**
		  * Registration name of this inference. It is used to identify this inference when requested 
		  * by the user in the IConfiguration.
		  */
		virtual std::string assignedRegistrationName() const {
			return "categorization.category";
		}
	};
};

class ObjectLayerInference : public AbstractLayerInference {

public:

	virtual LayerOutputObject* performInference(const AbstractLayerInputObject* input_object) ;

	class Factory : public IFactory {
	 public:
		virtual AbstractLayerInference* newInstance(const IConfiguration& config, const std::shared_ptr<VocabularyTree> vocabulary) const {
			return new ObjectLayerInference();
		}
		/**
		  * Registration name of this inference. It is used to identify this inference when requested 
		  * by the user in the IConfiguration.
		  */
		virtual std::string assignedRegistrationName() const {
			return "categorization.object";
		}
	};
};


/// @}
/// @}


#endif /* _CORE_CATEGORY_LAYER_INFERENCE_ */
