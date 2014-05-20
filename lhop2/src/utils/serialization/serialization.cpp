// support for storing/restoring classes to/from stream

#include "serialization.h"


void AbstractSerializer::serializeAbstractObject(const ISerializableObject& object, ostreamer& output_stream) {
	AbstractSerializer::IFactory* serializer_factory = object.getSerializer();
	
	if (serializer_factory != nullptr) {

		// save name of the serializer
		output_stream.write(serializer_factory->assignedRegistrationName());

		// get actual serializer
		AbstractSerializer* serializer = serializer_factory->newInstance(storage);

		// use serializer for this specific class to handle serialization
		serializer->serialize(object, output_stream);

		delete serializer;
		delete serializer_factory;
	}
}
ISerializableObject* AbstractSerializer::deserializeAbstractObject(istreamer& input_stream) {
	// get name of the serializer
	std::string name;
	input_stream.read(name);

	// create serializer from its factory (obtain factory first)
	AbstractSerializer::IFactory* serializer_factory = ClassRegister::get().retrieveFactory<AbstractSerializer::IFactory>(name);

	AbstractSerializer* serializer = serializer_factory->newInstance(storage);

	// unserialize the object
	ISerializableObject* object = serializer->deserialize(input_stream);

	delete serializer_factory;
	delete serializer;

	// TODO: object should get registered with the shared serialization storage

	return object;
}


int SerializedObjectStorage::registerSerializedObject(void* ptr, ostreamer& output_stream) {
	// check if pointer already registered
	auto existing_object = serialized_objects.find(ptr);

	int id;
	
	if (existing_object != serialized_objects.end()) {
		id = existing_object->second;
	} else {
		// generate new ID for this pointer
		id = counter++;

		serialized_objects[ptr] = id;
	}
	// save ID into the stream
	output_stream.write(id);

	return id;
}
void* SerializedObjectStorage::retriveSerializedObject(int id) {
	// TODO: skip
	return nullptr;
}

int SerializedObjectStorage::registerUnserializedObject(void* ptr, istreamer& input_stream) {
	// read id from stream
	int id;
	input_stream.read(id);

	// find if any pointer already registered to this ID
	auto existing_object = unserialized_objects.find(id);

	if (existing_object == unserialized_objects.end()) {
		unserialized_objects[id] = ptr;
	} else {
		// TODO: how do we handle if someone already registered
	}

	return id;
}
void* SerializedObjectStorage::retriveUnserializedObject(istreamer& input_stream) {
	// read id from stream 
	int id;
	input_stream.read(id);

	// retrieve pointer from map
	auto existing_object = unserialized_objects.find(id);

	if (existing_object == unserialized_objects.end()) {
		throw new_libhop_exception("Cannot retrieve reference to an object that has not been unserialized first.");
	}

	return existing_object->second;
}
