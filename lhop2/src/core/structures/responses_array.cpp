// part_responses

#include "responses_array.h"

#include "core/structures/serializers.h"

const ResponsesArray ResponsesArray::DEFAULT_RESPONSES = ResponsesArray();


AbstractSerializer::IFactory* ResponsesArray::getSerializer() const {
	return new DefaultResponsesArraySerializer::Factory();
}
