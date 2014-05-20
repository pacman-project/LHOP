/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// serializers_text

#include "serializers_text.h"

////////////////////////////////////////////////////////////////////////////////////////////////
//// ExtensionHolder

void TextExtensionHolderSerializer::serialize(const ISerializableObject& object, ostreamer& output_stream) {
	const ExtensionHolder& holder = (const ExtensionHolder&)object;

	output_stream << " ===== Functionalities ===== " << std::endl;

	auto func_list = getFunctionalities(holder);

	// serialize each extension
	for (auto iter = func_list.begin(); iter != func_list.end(); ++iter) {
		IExtension* extension = *iter;

		serializeAbstractObject(*extension, output_stream);
	}
}

ISerializableObject* TextExtensionHolderSerializer::deserialize(istreamer& input_stream, ISerializableObject* object) {		
	ExtensionHolder* holder =  object != nullptr ? (ExtensionHolder*)object : new ExtensionHolder();
	
	// skip

	return holder;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//// ResponsesArray

void TextResponsesArraySerializer::serialize(const ISerializableObject& object, ostreamer& output_stream) {
	const ResponsesArray& response_array = (const ResponsesArray&)object;
	
	auto response_values = getValues(response_array);

	output_stream << "Response Array: ";
	for (int i = 0; i < response_values.size(); ++i) {
		output_stream <<  response_values[i] << ", ";
	}	
}
ISerializableObject* TextResponsesArraySerializer::deserialize(istreamer& input_stream, ISerializableObject* object) {
	ResponsesArray* response_array = object != nullptr ? (ResponsesArray*)object : new ResponsesArray();
	
	// skip

	return response_array;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//// InferenceTree

 void TextInferenceTreeSerializer::serialize(const InferenceTree& inference_tree, ostreamer& output_stream) {

	output_stream << " ===== Inference Tree ===== " << std::endl;

	int number_of_layer = inference_tree.getNumberOfLayers();

	// write number of layers
	output_stream << "Number of layers: " << number_of_layer << std::endl;

	// serialize each layer
	for (int i = 0; i < number_of_layer; ++i) {
		
		output_stream << "Layer: " << i << std::endl;

		serializeInferenceLayer(inference_tree.getLayer(i), output_stream);
	}

	// serialize with parent after all parts are saved
	TextExtensionHolderSerializer::serialize(inference_tree, output_stream);
}

void TextInferenceTreeSerializer::serializeInferenceLayer(const InferenceLayer& inference_layer, ostreamer& output_stream) {

	output_stream << " ===== Inference Layer ===== " << std::endl;

	// write size of layer
	output_stream << "Size of layer: " << inference_layer.getSize().width << " " << inference_layer.getSize().height << std::endl;

	// write number of parts 
	int number_all_parts = 0;
	for (auto location = inference_layer.beginIterator(); location != inference_layer.endIterator(); ++location) {
		const std::vector<InferredPart*>& location_parts = inference_layer.getPartsAt(location);
		for (auto iter = location_parts.begin(); iter != location_parts.end(); ++iter) {
			++number_all_parts;
		}
	}

	output_stream << "Number of layer parts: " << number_all_parts << std::endl;

	for (auto location = inference_layer.beginIterator(); location != inference_layer.endIterator(); ++location) {
		const std::vector<InferredPart*>& location_parts = inference_layer.getPartsAt(location);
		for (auto iter = location_parts.begin(); iter != location_parts.end(); ++iter) {
			// serialize each part
			serializeInferredPart(**iter, output_stream);
		}
	}
}

void TextInferenceTreeSerializer::serializeInferredPart(const InferredPart& inferred_part, ostreamer& output_stream) {
	// serialize individual propreties
	output_stream << "\tUUID " << inferred_part.getUUID()
					<< ", location " << inferred_part.getLocation() 
					<< ", layer " << inferred_part.getLayer() 
					<< ", corresponding vocabulary part " << inferred_part.getCorrespondingVocabularyPartUUID() << " ";

	// serialize response array
	response_array_serializer.serialize(inferred_part.getResponseArray(), output_stream);

	output_stream << std::endl;
}

	
 InferenceTree* TextInferenceTreeSerializer::deserialize(istreamer& input_stream, InferenceTree* object) {
	// create inference tree object first
	InferenceTree* inference_tree = object != nullptr ? (InferenceTree*)object : new InferenceTree();
	
	// skip

	return inference_tree;
}

InferenceLayer* TextInferenceTreeSerializer::deserializeInferenceLayer(istreamer& input_stream, InferenceTree& inference_tree) {
	return nullptr;
}

InferredPart* TextInferenceTreeSerializer::deserializeInferredPart(istreamer& input_stream) {
	return nullptr;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularyTree

// TODO


////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularyIndexingFunctionality

void TextVocabularyIndexingSerializer::serialize(const ISerializableObject& object, ostreamer& output_stream) {
	const VocabularyIndexingExtension& func = (const VocabularyIndexingExtension&)object;
	
	output_stream << " ===== Vocabulary Indexing ===== " << std::endl;

	output_stream << "Number of indexes: " << func.indexed_center_links.size() << std::endl;

	int i = 0;
	for (auto iter = func.indexed_center_links.begin(); iter != func.indexed_center_links.end(); ++iter, ++i) {
		output_stream << "Index: " << i;
		output_stream << ", Key " << iter->first << ", Values:";
		for (auto iter_part = iter->second.begin(); iter_part != iter->second.begin(); ++iter_part) {
			output_stream << " " << (*iter_part)->getUUID();
		}
		output_stream << std::endl;
	}
}
ISerializableObject* TextVocabularyIndexingSerializer::deserialize(istreamer& input_stream, ISerializableObject* object){
	// TODO: write deserialization
	// TODO: how are we going to get VocabularyPart from UUIDs (solution: use storage->retriveRegisteredObject(..))
	return nullptr;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularySubpartFunctionality

// TODO


////////////////////////////////////////////////////////////////////////////////////////////////
//// InferenceSubpartFunctionality


void TextInferenceSubpartsFunctionalitySerializer::serialize(const ISerializableObject& object, ostreamer& output_stream){
	const InferenceSubpartsExtension& func = (const InferenceSubpartsExtension&)object;

	output_stream << " ===== Inference Subparts ===== " << std::endl;

	output_stream << "Number of indexes: " << func.subparts_links.size() << std::endl;

	int i = 0;
	for (auto iter = func.subparts_links.begin(); iter != func.subparts_links.end(); ++iter, ++i) {
		output_stream << "Index: " << i;
		output_stream << ", Key " << iter->first;
		output_stream << ", Values:";
		for (auto iter_part = iter->second.begin(); iter_part != iter->second.begin(); ++iter_part) {

			output_stream << "\tUUID " << (*iter_part)->getUUID() << ", r " << (*iter_part)->r << ", index " << (*iter_part)->index;
			
			output_stream << " Points to InferredPart UUID " << (*iter_part)->subpart.getUUID() << ", ";
		}
		output_stream << std::endl;
	}
				
}
ISerializableObject* TextInferenceSubpartsFunctionalitySerializer::deserialize(istreamer& input_stream, ISerializableObject* object) {
	return nullptr;
}
