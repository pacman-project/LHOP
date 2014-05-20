/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// serializers

#include "serializers.h"

////////////////////////////////////////////////////////////////////////////////////////////////
//// ExtensionHolder

void DefaultExtensionHolderSerializer::serialize(const ISerializableObject& object, ostreamer& output_stream) {
	const ExtensionHolder& holder = (const ExtensionHolder&)object;
	
	// register holder to the storage of serialized object (so that other serializers can find pointer to it)
	storage->registerSerializedObject((void*)&holder, output_stream);

	// get private variable of the holder (through protected getter that is allowed to access this variable)
	auto& func_list = getFunctionalities(holder);

	// write number of functionalities
	output_stream.write((int)func_list.size());

	// serialize each extension
	for (auto iter = func_list.begin(); iter != func_list.end(); ++iter) {
		IExtension* extension = *iter;

		serializeAbstractObject(*extension, output_stream);
	}
}

ISerializableObject* DefaultExtensionHolderSerializer::deserialize(istreamer& input_stream, ISerializableObject* object) {		
	ExtensionHolder* holder =  object != nullptr ? (ExtensionHolder*)object : new ExtensionHolder();

	// register holder to the storage of unserialized object (so that other serializers can find pointer to it)
	storage->registerUnserializedObject((void*)&holder, input_stream);

	// read number of functionalities
	int number_functionalities;
	input_stream.read(number_functionalities);

	// deserialize each extension and insert it into holder
	while (number_functionalities-- > 0) {
		holder->insertExtension((IExtension*)deserializeAbstractObject(input_stream));
	}

	return holder;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//// ResponsesArray

void DefaultResponsesArraySerializer::serialize(const ISerializableObject& object, ostreamer& output_stream) {
	const ResponsesArray& response_array = (const ResponsesArray&)object;

	auto& response_values = getValues(response_array);

	output_stream.write((int)response_values.size());

	for (auto iter = response_values.begin(); iter != response_values.end(); ++iter) {
		output_stream.write(*iter);
	}
}
ISerializableObject* DefaultResponsesArraySerializer::deserialize(istreamer& input_stream, ISerializableObject* object) {
	ResponsesArray* response_array = object != nullptr ? (ResponsesArray*)object : new ResponsesArray();

	int number_of_values;
	input_stream.read(number_of_values);

	auto& response_values = getValues(*response_array);

	response_values.resize(number_of_values);

	for (int i = 0; i < number_of_values; ++i)
		input_stream.read(response_values[i]);

	return response_array;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//// InferenceTree

 void DefaultInferenceTreeSerializer::serialize(const InferenceTree& inference_tree, ostreamer& output_stream) {

	int number_of_layer = inference_tree.getNumberOfLayers();

	// serialize UUID generators
	uuid_generator_serializer.serialize(inference_tree.getInferredPartUUIDGenerator(), output_stream);

	// write number of layers
	output_stream.write(number_of_layer);

	// serialize each layer
	for (int i = 0; i < number_of_layer; ++i) {
		serializeInferenceLayer(inference_tree.getLayer(i), output_stream);
	}

	// serialize with parent after all parts are saved
	DefaultExtensionHolderSerializer::serialize(inference_tree, output_stream);
}

void DefaultInferenceTreeSerializer::serializeInferenceLayer(const InferenceLayer& inference_layer, ostreamer& output_stream) {
	// write size of layer
	output_stream.write(inference_layer.getSize().width);
	output_stream.write(inference_layer.getSize().height);
	output_stream.write(inference_layer.getLayerIndex());

	// write number of parts 
	int number_all_parts = 0;
	for (auto location = inference_layer.beginIterator(); location != inference_layer.endIterator(); ++location) {
		const std::vector<InferredPart*>& location_parts = inference_layer.getPartsAt(location);
		for (auto iter = location_parts.begin(); iter != location_parts.end(); ++iter) {
			++number_all_parts;
		}
	}
	output_stream.write(number_all_parts);

	for (auto location = inference_layer.beginIterator(); location != inference_layer.endIterator(); ++location) {
		const std::vector<InferredPart*>& location_parts = inference_layer.getPartsAt(location);
		for (auto iter = location_parts.begin(); iter != location_parts.end(); ++iter) {
			// serialize each part
			serializeInferredPart(**iter, output_stream);
		}
	}
}

void DefaultInferenceTreeSerializer::serializeInferredPart(const InferredPart& inferred_part, ostreamer& output_stream) {
	// serialize individual propreties
	output_stream.write(inferred_part.getUUID());
	output_stream.write(inferred_part.getLocation().x);
	output_stream.write(inferred_part.getLocation().y);
	output_stream.write(inferred_part.getLayer());
	output_stream.write(inferred_part.getCorrespondingVocabularyPartUUID());

	// serialize response array
	response_array_serializer.serialize(inferred_part.getResponseArray(), output_stream);

	// register part to the storage of serialized object (so that other serializers can find pointer to it)
	storage->registerSerializedObject((void*)&inferred_part, output_stream);
}

	
 InferenceTree* DefaultInferenceTreeSerializer::deserialize(istreamer& input_stream, InferenceTree* object) {
	// create inference tree object first
	InferenceTree* inference_tree = object != nullptr ? (InferenceTree*)object : new InferenceTree();

	// unserialize UUID generators
	uuid_generator_serializer.deserialize(input_stream, &inference_tree->getInferredPartUUIDGenerator());
	
	int number_of_layers;
	input_stream.read(number_of_layers);

	// deserialize individual layers
	while (number_of_layers-- > 0) {
		inference_tree->insertNewLayer(deserializeInferenceLayer(input_stream, *inference_tree));
	}

	// deserialize parent
	DefaultExtensionHolderSerializer::deserialize(input_stream, inference_tree);

	return inference_tree;
}

InferenceLayer* DefaultInferenceTreeSerializer::deserializeInferenceLayer(istreamer& input_stream, InferenceTree& inference_tree) {
	// read size of layer
	cv::Size2i size;
	input_stream.read(size.width);
	input_stream.read(size.height);
	
	int layer_index;
	input_stream.read(layer_index);

	InferenceLayer* inference_layer = new InferenceLayer(inference_tree, size, layer_index);

	// read number of parts
	int number_all_parts = 0;
	input_stream.read(number_all_parts);

	// deserialize individual parts
	while (number_all_parts-- > 0) {
		inference_layer->insertNewPart(deserializeInferredPart(input_stream));
	}

	return inference_layer;
}

InferredPart* DefaultInferenceTreeSerializer::deserializeInferredPart(istreamer& input_stream) {
	UUIDType uuid;
	cv::Point2i location;
	int layer;
	UUIDType corresponding_vocabulary_part_uuid;
	ResponsesArray response_array;

	input_stream.read(uuid);
	input_stream.read(location.x);
	input_stream.read(location.y);
	input_stream.read(layer);
	input_stream.read(corresponding_vocabulary_part_uuid);
		
	response_array_serializer.deserialize(input_stream, &response_array);

	InferredPart* part = new InferredPart(uuid, layer, location, corresponding_vocabulary_part_uuid, response_array);

	 // register part to the storage of unserialized object (so that other serializers can find pointer to it)
	storage->registerUnserializedObject((void*)&part, input_stream);

	return part;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularyTree

// TODO


////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularyIndexingFunctionality


void DefaultVocabularyIndexingSerializer::serialize(const ISerializableObject& object, ostreamer& output_stream){
	const VocabularyIndexingExtension& func = (const VocabularyIndexingExtension&)object;
	
	// register holder to the storage of serialized object 
	storage->registerSerializedObject((void*)&func.holder, output_stream);

	output_stream.write((int)func.indexed_center_links.size());

	for (auto iter = func.indexed_center_links.begin(); iter != func.indexed_center_links.end(); ++iter) {
		output_stream.write(iter->first);
		output_stream.write((int)iter->second.size());
		for (auto iter_part = iter->second.begin(); iter_part != iter->second.begin(); ++iter_part) {
			output_stream.write((*iter_part)->getUUID());
		}
	}
}
ISerializableObject* DefaultVocabularyIndexingSerializer::deserialize(istreamer& input_stream, ISerializableObject* object){
	// TODO: write deserialization
	// TODO: how are we going to get VocabularyPart from UUIDs (solution: use storage->retriveRegisteredObject(..))
	return nullptr;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//// VocabularySubpartFunctionality

// TODO


////////////////////////////////////////////////////////////////////////////////////////////////
//// InferenceSubpartFunctionality


void DefaultInferenceSubpartsExtSerializer::serialize(const ISerializableObject& object, ostreamer& output_stream){
	const InferenceSubpartsExtension& func = (const InferenceSubpartsExtension&)object;

	// register holder to the storage of serialized object 
	storage->registerSerializedObject((void*)&func.holder, output_stream);

	// serialize UUID generators
	uuid_generator_serializer.serialize(func.uuid_generators_inference_subpart_data, output_stream);
	
	output_stream.write((int)func.subparts_links.size());

	for (auto iter = func.subparts_links.begin(); iter != func.subparts_links.end(); ++iter) {
		output_stream.write(iter->first);
		output_stream.write((int)iter->second.size());
		for (auto iter_part = iter->second.begin(); iter_part != iter->second.begin(); ++iter_part) {
			output_stream.write((*iter_part)->getUUID());
			output_stream.write((*iter_part)->r);
			output_stream.write((*iter_part)->index);
			
			// subpart (InferredPart) points to object within InferenceTree and is therfore the responsibility of 
			// InferenceTree serializer to serialize it - we only register it here to obtain an id that we can save the reference
			
			// register part to the storage of serialized object (so that other serializers can find pointer to it)
			storage->registerSerializedObject((void*)&(*iter_part)->subpart, output_stream);
		}
	}
				
}
ISerializableObject* DefaultInferenceSubpartsExtSerializer::deserialize(istreamer& input_stream, ISerializableObject* object) {
	
	// first retrive reference to holder (as it is required by the constructor of the InferenceSubpartsFunctionality)
	ExtensionHolder* holder = (ExtensionHolder*)storage->retriveUnserializedObject(input_stream);

	InferenceSubpartsExtension* func = object != nullptr ? (InferenceSubpartsExtension*)object : new InferenceSubpartsExtension(*holder);

	// deserialize UUID generators
	uuid_generator_serializer.deserialize(input_stream, &func->uuid_generators_inference_subpart_data);

	int subpart_links_size;
	input_stream.read(subpart_links_size);

	while (subpart_links_size-- > 0) {
		int reference_key;
		input_stream.read(reference_key);

		int number_of_subparts;
		input_stream.read(number_of_subparts);
		while (number_of_subparts-- > 0) {
			UUIDType uuid;
			double r;
			int index;
			int subpart_reference_id;

			input_stream.read(uuid);
			input_stream.read(r);
			input_stream.read(index);
		
			// retrive already unserialized object (must have already been unserailized and registerd by its own serilizer called before this one TODO: how do we control order of the serializers)
			InferredPart* reference_subpart = (InferredPart*)storage->retriveUnserializedObject(input_stream);
			
			func->subparts_links[(size_t)reference_key].push_back(new InferenceSubpartData(uuid, *reference_subpart, index, r));
		}
	}
	
	return func;
}
