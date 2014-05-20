// attached_extension

#include "attached_extension.h"

#include "core/structures/serializers.h"

AbstractSerializer::IFactory* ExtensionHolder::getSerializer() const {
	return new DefaultExtensionHolderSerializer::Factory();
}


