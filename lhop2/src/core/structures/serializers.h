/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// serializers

#pragma once
#ifndef _CORE_SERIALIZERS_
#define _CORE_SERIALIZERS_

#include "utils/serialization/serialization.h"

#include "core/structures/responses_array.h"
#include "core/structures/vocabulary.h"
#include "core/structures/parse_tree.h"

#include "core/structures/indexing/indexing_extension.h"
#include "core/structures/subparts/subparts_extension.h"

////////////////////////////////////////////////////////////////////////////////////////////////
//// ExtensionHolder

/// @addtogroup core
/// @{
/// @addtogroup main_structures
/// @{
/// @addtogroup serialization
/// @{


/**
 * Abstract class for ExtensionHolder serializer. Has special premission (is a friend class)  to access
 * private variables of the ExtensionHolder and exposes them as protected getters.
 */ 
class AbstractExtensionHolderSerializer : public AbstractSerializer {
public:
	AbstractExtensionHolderSerializer(std::shared_ptr<SerializedObjectStorage> storage) : AbstractSerializer(storage) {
	}
protected:
	// expose methods needed for the serialization
	std::set<IExtension*>& getFunctionalities(ExtensionHolder& holder) const { return holder.extensions; }
	const std::set<IExtension*>& getFunctionalities(const ExtensionHolder& holder) const { return holder.extensions; }
};


class DefaultExtensionHolderSerializer : public AbstractExtensionHolderSerializer {
public:
	DefaultExtensionHolderSerializer(std::shared_ptr<SerializedObjectStorage> storage) : AbstractExtensionHolderSerializer(storage) {
	}

	virtual void serialize(const ISerializableObject& object, ostreamer& output_stream);
	virtual ISerializableObject* deserialize(istreamer& input_stream, ISerializableObject* object = nullptr);

	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new DefaultExtensionHolderSerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "DefaultExtensionHolderSerializer"; }
	};
};


////////////////////////////////////////////////////////////////////////////////////////////////
//// ResponsesArray

/**
 * Abstract class for ResponsesArray serializer. Has special premission (is a friend class) to access
 * private variables of the ResponsesArray and exposes them as protected getters.
 */ 
class AbstractResponsesArraySerializer : public AbstractSerializer {
public:
	AbstractResponsesArraySerializer(std::shared_ptr<SerializedObjectStorage> storage) : AbstractSerializer(storage) {
	}
protected:
	// expose methods needed for the serialization
	std::vector<ResponsesArray::ValueType>& getValues(ResponsesArray& responses) const { return responses.values; }
	const std::vector<ResponsesArray::ValueType>& getValues(const ResponsesArray& responses) const { return responses.values; }
};

class DefaultResponsesArraySerializer : public AbstractResponsesArraySerializer {
public:
	DefaultResponsesArraySerializer(std::shared_ptr<SerializedObjectStorage> storage) : AbstractResponsesArraySerializer(storage) {
	}

	virtual void serialize(const ISerializableObject& object, ostreamer& output_stream);
	virtual ISerializableObject* deserialize(istreamer& input_stream, ISerializableObject* object = nullptr);

	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new DefaultResponsesArraySerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "DefaultResponseArraySerializer"; }
	};
};

////////////////////////////////////////////////////////////////////////////////////////////////
//// InferenceTree

class DefaultInferenceTreeSerializer : public DefaultExtensionHolderSerializer {
	DefaultUUIDGeneratorSerializer uuid_generator_serializer;
	DefaultResponsesArraySerializer response_array_serializer;
public:
	DefaultInferenceTreeSerializer(std::shared_ptr<SerializedObjectStorage> storage) 
		: DefaultExtensionHolderSerializer(storage), response_array_serializer(storage), uuid_generator_serializer(storage) {
	}

	virtual void serialize(const InferenceTree& inference_tree, ostreamer& output_stream);
	void serializeInferenceLayer(const InferenceLayer& inference_layer, ostreamer& output_stream);
	void serializeInferredPart(const InferredPart& inference_layer, ostreamer& output_stream);
	
	virtual InferenceTree* deserialize(istreamer& input_stream, InferenceTree* object = nullptr);
	InferenceLayer* deserializeInferenceLayer(istreamer& input_stream, InferenceTree& inference_tree);
	InferredPart* deserializeInferredPart(istreamer& input_stream);

	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new DefaultInferenceTreeSerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "DefaultInferenceTreeSerializer"; }
	};
};


////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularyTree

// TODO


////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularyIndexingFunctionality

class DefaultVocabularyIndexingSerializer : public AbstractSerializer {
public:
	DefaultVocabularyIndexingSerializer(std::shared_ptr<SerializedObjectStorage> storage) : AbstractSerializer(storage) {
	}

	virtual void serialize(const ISerializableObject& object, ostreamer& output_stream);
	virtual ISerializableObject* deserialize(istreamer& input_stream, ISerializableObject* object = nullptr);
	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new DefaultVocabularyIndexingSerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "DefaultVocabularyIndexingAccessSerializer"; }
	};
};

////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularySubpartFunctionality

// TODO


////////////////////////////////////////////////////////////////////////////////////////////////
//// InferenceSubpartFunctionality

class DefaultInferenceSubpartsExtSerializer : public AbstractSerializer {
	DefaultUUIDGeneratorSerializer uuid_generator_serializer;
public:
	DefaultInferenceSubpartsExtSerializer(std::shared_ptr<SerializedObjectStorage> storage) 
		: AbstractSerializer(storage), uuid_generator_serializer(storage)  {
	}
	virtual void serialize(const ISerializableObject& object, ostreamer& output_stream);
	virtual ISerializableObject* deserialize(istreamer& input_stream, ISerializableObject* object = nullptr) ;

	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new DefaultInferenceSubpartsExtSerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "DefaultInferenceSubpartsExtSerializer"; }
	};
};

/// @}
/// @}
/// @}


#endif /* _CORE_SERIALIZERS_ */

