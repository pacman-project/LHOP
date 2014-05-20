// support for storing/restoring classes to/from stream

#pragma once
#ifndef _SERIALIZATION_H_
#define _SERIALIZATION_H_

#include <string>
#include <memory>

#include "utils/class_register.h"
#include "utils/serialization/streaming.h"

class ISerializableObject;
class AbstractSerializer;
class SerializedObjectStorage;

class AbstractSerializer : public IRegistrableClass {
protected:
	std::shared_ptr<SerializedObjectStorage> storage;
public:
	AbstractSerializer(std::shared_ptr<SerializedObjectStorage> storage) : storage(storage) {}

	// must be implemented by the descendent
	virtual void serialize(const ISerializableObject& object, ostreamer& output_stream) = 0;
	virtual ISerializableObject* deserialize(istreamer& input_stream, ISerializableObject* object = nullptr) = 0;

	class IFactory : public IRegistrableClassFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const = 0;
		virtual std::string assignedRegistrationName() const = 0;
	};

protected:
	/**
	 * Serializes any object extending from ISerializableObject for which you do not known the
	 * exact type during the compile time of the serialization and therfore its serializer class 
     * has to be retrived from the object itself. The serializer class id will also be written
	 * into the stream and will be used by the deserializeAbstractObject(..) to perform correct 
     * de-serialization.
	 */
	void serializeAbstractObject(const ISerializableObject& object, ostreamer& output_stream);
	/*
     * Performs de-serialization of an abstract object whoses class type was not known during the 
	 * compile time. The serializer class should be saved during the serialization and will be used
     * during the deserialize process.
	 */
	ISerializableObject* deserializeAbstractObject(istreamer& input_stream) ;
};
 

class ISerializableObject {
public:
	virtual AbstractSerializer::IFactory* getSerializer() const = 0;
};

class SerializedObjectStorage {
	int counter;
	
	std::map<void*, int> serialized_objects; // used only by the registerSerializedObject() and retriveSerializedObject()
	std::map<int, void*> unserialized_objects;// used only by the registerUnserializedObject() and retriveUnserializedObject()
public:
	SerializedObjectStorage() : counter(0) {
	}

	/**
     * Assignes ID to the pointer, saves it into its map and
	 * writes it into the stream.
     */
	int registerSerializedObject(void* ptr, ostreamer& output_stream);
	/**
     * Is it needed ?
	 */
	void* retriveSerializedObject(int id);

    /**
     * Registers unserilized pointer to an object by 
     * reading its ID from stream and returning it if needed.
     */
	int registerUnserializedObject(void* ptr, istreamer& input_stream);
	/**
	 * Retrieves pointer to already unserialized object
	 * based on its id.
	 */
	void* retriveUnserializedObject(istreamer& input_stream);
};


#endif /* _SERIALIZATION_H_ */
