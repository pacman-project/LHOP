/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// serializers

#pragma once
#ifndef _CORE_SERIALIZERS_TEXT_
#define _CORE_SERIALIZERS_TEXT_

#include "utils/serialization/serialization.h"

#include "core/structures/responses_array.h"
#include "core/structures/vocabulary.h"
#include "core/structures/parse_tree.h"

#include "core/structures/indexing/indexing_extension.h"
#include "core/structures/subparts/subparts_extension.h"

#include "core/structures/serializers.h"

////////////////////////////////////////////////////////////////////////////////////////////////
//// ExtensionHolder

/// @addtogroup core
/// @{
/// @addtogroup main_structures
/// @{
/// @addtogroup serialization
/// @{

class TextExtensionHolderSerializer : public DefaultExtensionHolderSerializer {
public:
	TextExtensionHolderSerializer(std::shared_ptr<SerializedObjectStorage> storage) : DefaultExtensionHolderSerializer(storage) {
	}

	virtual void serialize(const ISerializableObject& object, ostreamer& output_stream);
	virtual ISerializableObject* deserialize(istreamer& input_stream, ISerializableObject* object = nullptr);

	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new TextExtensionHolderSerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "TextFunctionalityHolderSerializer"; }
	};
};


////////////////////////////////////////////////////////////////////////////////////////////////
//// ResponsesArray

class TextResponsesArraySerializer : public DefaultResponsesArraySerializer {
public:
	TextResponsesArraySerializer(std::shared_ptr<SerializedObjectStorage> storage) : DefaultResponsesArraySerializer(storage) {
	}

	virtual void serialize(const ISerializableObject& object, ostreamer& output_stream);
	virtual ISerializableObject* deserialize(istreamer& input_stream, ISerializableObject* object = nullptr);

	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new TextResponsesArraySerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "TextResponsesArraySerializer"; }
	};
};

////////////////////////////////////////////////////////////////////////////////////////////////
//// InferenceTree

class TextInferenceTreeSerializer : public TextExtensionHolderSerializer {
	TextResponsesArraySerializer response_array_serializer;
public:
	TextInferenceTreeSerializer(std::shared_ptr<SerializedObjectStorage> storage) : TextExtensionHolderSerializer(storage), response_array_serializer(storage) {
	}

	virtual void serialize(const InferenceTree& inference_tree, ostreamer& output_stream);
	void serializeInferenceLayer(const InferenceLayer& inference_layer, ostreamer& output_stream);
	void serializeInferredPart(const InferredPart& inference_layer, ostreamer& output_stream);
	
	virtual InferenceTree* deserialize(istreamer& input_stream, InferenceTree* object = nullptr);
	InferenceLayer* deserializeInferenceLayer(istreamer& input_stream, InferenceTree& inference_tree);
	InferredPart* deserializeInferredPart(istreamer& input_stream);

	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new TextInferenceTreeSerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "TextInferenceTreeSerializer"; }
	};
};


////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularyTree

// TODO


////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularyIndexingFunctionality


class TextVocabularyIndexingSerializer : public DefaultVocabularyIndexingSerializer {
public:
	TextVocabularyIndexingSerializer(std::shared_ptr<SerializedObjectStorage> storage) : DefaultVocabularyIndexingSerializer(storage) {
	}

	virtual void serialize(const ISerializableObject& object, ostreamer& output_stream);
	virtual ISerializableObject* deserialize(istreamer& input_stream, ISerializableObject* object = nullptr);
	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new TextVocabularyIndexingSerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "TextVocabularyIndexingAccessSerializer"; }
	};
};

////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularySubpartFunctionality

// TODO


////////////////////////////////////////////////////////////////////////////////////////////////
//// InferenceSubpartFunctionality


class TextInferenceSubpartsFunctionalitySerializer : public DefaultInferenceSubpartsExtSerializer {
public:
	TextInferenceSubpartsFunctionalitySerializer(std::shared_ptr<SerializedObjectStorage> storage) : DefaultInferenceSubpartsExtSerializer(storage) {
	}
	virtual void serialize(const ISerializableObject& object, ostreamer& output_stream);
	virtual ISerializableObject* deserialize(istreamer& input_stream, ISerializableObject* object = nullptr) ;

	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new TextInferenceSubpartsFunctionalitySerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "TextInferenceSubpartsFunctionalitySerializer"; }
	};
};


/// @}
/// @}
/// @}

#endif /* _CORE_SERIALIZERS_TEXT_ */

