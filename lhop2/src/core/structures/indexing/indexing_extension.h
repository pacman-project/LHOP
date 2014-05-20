/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// indexing_extension

#pragma once
#ifndef _CORE_INDEXING_EXTENSION_
#define _CORE_INDEXING_EXTENSION_

#include <vector>
#include <unordered_map>

#include "core/structures/parse_tree.h"
#include "core/structures/vocabulary.h"
#include "core/structures/attached_extension.h"

/// @addtogroup core
/// @{
/// @addtogroup main_structures
/// @{

/// @addtogroup ext
/// @{
/// @addtogroup ext_vocabulary_tree
/// @{

////////////////////////////////////////////////////////////////////////////////
//// Vocabulary parts indexing and its access

class IVocabularyIndexingAccess;
class IVocabularyIndexingModifier;
class VocabularyIndexingAccess;
class VocabularyIndexingModifier;
class VocabularyIndexingExtension;

class IVocabularyIndexingAccess : public IExtensionAccess {
public:
	virtual std::vector<VocabularyPart*> getIndexedCentralParts(const VocabularyPart& part) = 0;
};

class IVocabularyIndexingModifier : public IExtensionModifier {
public:
	virtual void insertIndexedCentralParts(const IAttachableClass& part, const std::vector<VocabularyPart*>& indexes_for_insertion) = 0;
	virtual void deleteIndexedCentralParts(const IAttachableClass& part, const std::vector<VocabularyPart*>& indexes_for_deletion) = 0;
};

/**
 * Functionality provides an reverse index information about which vocabulary parts are
 * used as an central subpart for another vocabulary part. It is attached to the
 * VocabularyTree, specifically VocabularyPart. It does not store any other information
 * except for the pointer to the VocabularyPart containing reference part as it center.
 * Access is provide through  VocabularyIndexingAccess class (read-only access) and 
 * VocabularyIndexingModifier class (write-only access).
 */
class VocabularyIndexingExtension : public SimpleAbstractExtension<VocabularyIndexingExtension,
																	VocabularyIndexingAccess,
																	VocabularyIndexingModifier,
																	VocabularyPart> {
public:
	ExtensionHolder& holder;

	// map: UUID of vocabulary part to list of vocabulary parts (the list are all parts indexed as users of UUID as its central part
	std::unordered_map<UUIDType, std::vector<VocabularyPart*>> indexed_center_links;


	VocabularyIndexingExtension(ExtensionHolder& holder) : holder(holder) {
	}


	virtual void moveTo(IExtension& dst_extension);
	
	virtual AbstractSerializer::IFactory* getSerializer() const;
};

class VocabularyIndexingAccess : public IVocabularyIndexingAccess {
	VocabularyIndexingExtension& ext;
public:
	VocabularyIndexingAccess(VocabularyIndexingExtension& ext) : ext(ext) {

	}
	virtual std::vector<VocabularyPart*> getIndexedCentralParts(const VocabularyPart& part) {
		return ext.indexed_center_links[part.getUUID()];
	}
};


class VocabularyIndexingModifier : public IVocabularyIndexingModifier {
	VocabularyIndexingExtension& ext;
public:
	VocabularyIndexingModifier(VocabularyIndexingExtension& ext) : ext(ext) {

	}
	/**
	 * Input argument must contain only pointers to InferredPart.
	 * Method deletes any links that used the parts in the argument list as
	 * theirs reference point. 
	 */ 
	virtual void deleteAllReferences(const std::vector<IAttachableClass*>& reference_to);

	/**
	 * Creates an index link between input part and provided list of subparts. The link represent an reversed index
	 * i.e. all parts in the indexes_for_insertion use part as its central part.
	 */
	virtual void insertIndexedCentralParts(const IAttachableClass& part, const std::vector<VocabularyPart*>& indexes_for_insertion);

	/**
	 * Removes provided links associated with the input part.
	 */
	virtual void deleteIndexedCentralParts(const IAttachableClass& part, const std::vector<VocabularyPart*>& indexes_for_deletion);


	
};

/// @}
/// @}
/// @}
/// @}

#endif /* _CORE_INDEXING_EXTENSION_ */

