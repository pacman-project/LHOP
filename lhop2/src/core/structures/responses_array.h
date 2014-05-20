// part_responses

#pragma once
#ifndef _CORE_RESPONSES_ARRAY_
#define _CORE_RESPONSES_ARRAY_

#include <vector>

#include "utils/serialization/serialization.h"


/// @addtogroup core
/// @{
/// @addtogroup main_structures
/// @{

class AbstractResponsesArraySerializer;

enum ResponseType { R_RESPONSE_1, G_RESPONSE_1, RR_RESPONSE_1, S_RESPONSE_1, X_RESPONSE_1 };

/**
 * Response array containing a list of responses (based on ResponseTypes). It is used by the
 * InferredPart as array of responses calculated during the inference process and by the 
 * VocabularyPart as array of thresholds.
 */ 
class ResponsesArray : public ISerializableObject {
	// grant access to abstract serializer that will expose private methods as protected getters
	friend class AbstractResponsesArraySerializer;
public:
	// by default use double as type (change to float if do not need that accuracy)
	typedef double ValueType;
	static const ResponsesArray DEFAULT_RESPONSES;

private:
	std::vector<ValueType> values;
public:
	ResponsesArray() : values(5,-1) {
	}
	ValueType get(ResponseType response_type) const {
		return values[response_type];
	}

	void set(ResponseType response_type, ValueType value) {
		values[response_type] = value;
	}

	virtual AbstractSerializer::IFactory* getSerializer() const;
};

/// @}
/// @}

#endif /* _CORE_RESPONSES_ARRAY_ */

