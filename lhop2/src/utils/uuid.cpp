
#include "uuid.h"


////////////////////////////////////////////////////////////////////////////////////////////////
//// Serializers for utility classes

AbstractSerializer::IFactory* UUIDGenerator::getSerializer() const  {
	return new DefaultUUIDGeneratorSerializer::Factory();
}

void DefaultUUIDGeneratorSerializer::serialize(const ISerializableObject& object, ostreamer& output_stream) {
	UUIDGenerator& uuid_generator = (UUIDGenerator&)object;
	
	// we only need to write current value
	output_stream.write(getCurrentValue(uuid_generator));
}

ISerializableObject* DefaultUUIDGeneratorSerializer::deserialize(istreamer& input_stream, ISerializableObject* object) {
	// we only need to read current value
	UUIDType current_value;
	input_stream.read(current_value);
	
	UUIDGenerator* uuid_generator = (UUIDGenerator*)object;

	if (object == nullptr) {
		object = new UUIDGenerator(current_value);
	} else {
		// get pointer to current value of the uuid_generator
		UUIDType& current_value_ptr = getCurrentValue(*uuid_generator);

		// and assign value read from stream to that pointer
		current_value_ptr = current_value;
	}
	return object;
}