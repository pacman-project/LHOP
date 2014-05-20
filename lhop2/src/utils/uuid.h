
#pragma once
#ifndef _UTILS_UUID_
#define _UTILS_UUID_

#include <cstdint>

#include <stdlib.h>
#include <time.h>

#include "utils/serialization/serialization.h"

class AbstractUUIDGeneratorSerializer;
class DefaultUUIDGeneratorSerializer;

// use new C++ 64 bit integer (unsigned) TODO: if replace with something if some compiler do not support it yet
typedef uint64_t UUIDType;

/**
 * Sequential UUID generator using the initial value.
 * Obtain new UUID by invoking generateUUID() which will return current assigned 
 * value and increment the internal counter;
 * 
 * Use static newInstance<T>() function to construct a generator for specific 
 * class T.
 * 
 */ 
class UUIDGenerator : public ISerializableObject {
	friend class AbstractUUIDGeneratorSerializer;
	UUIDType current_value;
public:
	UUIDGenerator(UUIDType current_value) : current_value(current_value) {
	}

	UUIDType generateUUID() {
		return current_value++;
	}

	/**
	 * Create UUID generator for specific class T by assigning
	 * first 16 bits to class ID, next 8 bit assigned randomly
	 * and next 40 bits allocated for sequential generation.
	 *
	 * We obtaint ID for class T from its hash code using type(T).hash_code() and
	 * directly cast it into 16 bits (we discard all the other bits).
	 *
	 * TODO: hopefully the 16 bit hash + 8 random bits will not clash, but
	 * maybe we should be more inteligent in the future !!!
	 */
	template <class T>
	static UUIDGenerator newInstance() {
		size_t class_hash = typeid(T).hash_code();

		// TODO: maybe we can use every second (or third) bit from class_hash
		short class_id = (short)class_hash;

		// generate randomly 8 bits
		srand (time(nullptr));
		char random_bits = rand() % 256;

		// combine class_id + random bits into starting_uuid_value
		UUIDType starting_uuid_value = 0;
		
		starting_uuid_value |= class_id;
		
		starting_uuid_value <<= sizeof(random_bits) * 8;
		starting_uuid_value |= random_bits;

		// calculate the size of UUIDType and shift it so that class_id + random_bits will be tle left most bits
		starting_uuid_value <<= (sizeof(starting_uuid_value) - (sizeof(class_id) + sizeof(random_bits)) ) *8;

		// do not return pointer since UUIDGenerator is a small class with only one value
		// and copy of pointer or of the whole class should require the same amount of time
		return UUIDGenerator(starting_uuid_value);
	}

	virtual AbstractSerializer::IFactory* getSerializer() const;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//// Serializer classes for UUIDGenerator

#include "utils/serialization/serialization.h"

/**
 * Abstract serializer for UUIDGenerator that exposes UUIDGenerator::current_value through 
 * protected method getCurrentValue(..)
 */ 
class AbstractUUIDGeneratorSerializer : public AbstractSerializer {
public:
	AbstractUUIDGeneratorSerializer(std::shared_ptr<SerializedObjectStorage> storage) : AbstractSerializer(storage) {
	}
protected:
	UUIDType& getCurrentValue(UUIDGenerator& uuid_generator) { return uuid_generator.current_value; }
};

/**
 * Default serializer for the UUIDGenerator class.
 */
class DefaultUUIDGeneratorSerializer : public AbstractUUIDGeneratorSerializer {
public:
	DefaultUUIDGeneratorSerializer(std::shared_ptr<SerializedObjectStorage> storage) : AbstractUUIDGeneratorSerializer(storage) {
	}
	virtual void serialize(const ISerializableObject& object, ostreamer& output_stream);
	virtual ISerializableObject* deserialize(istreamer& input_stream, ISerializableObject* object = nullptr) ;

	class Factory : public IFactory {
	public:
		virtual AbstractSerializer* newInstance(std::shared_ptr<SerializedObjectStorage> storage) const { return new DefaultUUIDGeneratorSerializer(storage); }
		virtual std::string assignedRegistrationName() const { return "UUIDGenerator"; }
	};
};


#endif


