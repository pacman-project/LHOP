/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// support_extension

#pragma once
#ifndef _CORE_SUPPORT_EXTENSION_
#define _CORE_SUPPORT_EXTENSION_

#include <typeinfo>
#include <typeindex>
#include <vector>
#include <unordered_map>

#include "core/structures/parse_tree.h"
#include "core/structures/vocabulary.h"
#include "core/structures/attached_extension.h"
#include "core/structures/subparts/subparts_extension.h"

/// @addtogroup core
/// @{
/// @addtogroup main_structures
/// @{

/// @addtogroup ext
/// @{

////////////////////////////////////////////////////////////////////////////////
//// Inference subparts extension

/// @addtogroup ext_support
/// @{

// forward definitions
class IInferenceSupportAccess;
class InferenceSupportAccess;
class InferenceSupportModifier;
class InferenceSupportExtension;
class InferenceSupportData;

////////////////////////////////////////////////////////////////////////////////
//// Implementation for support of InferenceTree 


class IInferenceSupportAccess : public IExtensionAccess {
public:
	virtual const std::vector<InferenceSupportData*>& getAllInitialLayerSupportParts(const IAttachableClass& part) = 0;
	virtual std::set<InferredPart*> getOnlyDifferentInitialLayerSupportParts(const IAttachableClass& part) = 0;
};


// PartSupportPath == a list of InferredPart UUIDs that define 
// parts from root to its leafs at the first layer (i.e. support parts)
// each UUIDType is one node at that path
typedef std::vector<UUIDType> PartSupportPath;

/**
 * 
 */ 
struct InferenceSupportData : public IAttachableClass {
	// first layer subpart reference (i.e. a supporting element of a part)
	InferredPart& support_subpart;
	
	// defines path between root node (part) and its leaf on the first layer
	PartSupportPath support_path;

	InferenceSupportData(UUIDType uuid, InferredPart& support_subpart, const PartSupportPath& support_path) 
		: IAttachableClass(uuid), support_subpart(support_subpart), support_path(support_path) {

	}
};


/**
 * 
 */
class InferenceSupportExtension : public SimpleAbstractExtension<InferenceSupportExtension,
																	InferenceSupportAccess,
																	InferenceSupportModifier,
																	InferredPart> {
public:
	// UUID generators
	UUIDGenerator uuid_generators_inference_support_data; // for class InferenceSupportData 

	ExtensionHolder& holder;

	// we map InferredPart* -> std::vector<InferenceSupportData*>
	// where we use InferredPart->getUUID() as key identifier
	std::unordered_map<UUIDType, std::vector<InferenceSupportData*>> support_subparts_links;

	// default constructor
	InferenceSupportExtension(ExtensionHolder& holder) 
		: holder(holder), uuid_generators_inference_support_data(UUIDGenerator::newInstance<InferenceSupportData>()) {
	}
	
	virtual void moveTo(IExtension& dst_extension); // moves all data into dst_extension and clear this one

	class Factory : public IFactory {
	public:
		virtual IExtension* newInstance(ExtensionHolder& holder) const { return new InferenceSupportExtension(holder); }
	};

	/**
	 * Returns intersection factor between two sets of supporting parts.
	 * Intersection factor is defined as count(intersection_parts(support_a,support_b)) / min(size(support_a), size(support_b))
	 */
	static float calculatePartsIntersection(const std::set<InferredPart*> support_a, const std::set<InferredPart*> support_b);

};

////////////////////////////////////////////////
/// Access class for retrieving support of InferenceTree

/**
 * Class provides access to the subparts stored within the InferenceTree object.
 * Access is read-only from the getSubparts(..)
 */ 
class InferenceSupportAccess : public IInferenceSupportAccess {	
	InferenceSupportExtension& ext;
	InferenceSubpartsAccess subparts_access;
public:
	InferenceSupportAccess(InferenceSupportExtension& ext) : ext(ext), subparts_access(ext.holder.getAccess<InferenceSubpartsAccess>()) {
	}

	/**
	 * Retrieves list of all supporting parts (subparts) including their paths to the root 
	 * that are associated with the specific input part. 
	 */
	virtual const std::vector<InferenceSupportData*>& getAllInitialLayerSupportParts(const IAttachableClass& part);		

	/**
	 * Get only a set of different InferredParts that act as supporting parts on first layer 
	 * (we do not know how many time each supporting part is used in the parse tree).
	 * 
	 * TODO: currently comparison is done by pointer value and not by actual value of the part (uuid, or location)
	 *       so be careful when using it !!!
	 */
	virtual std::set<InferredPart*> getOnlyDifferentInitialLayerSupportParts(const IAttachableClass& part);


};


////////////////////////////////////////////////
/// Modifier class for support of InferenceTree

/**
 * Class is only needed for deletion of its reference.
 */
class InferenceSupportModifier : public IExtensionModifier {
	InferenceSupportExtension& ext;
public:
	InferenceSupportModifier(InferenceSupportExtension& ext) : ext(ext) {
	}
	
	/**
	 * Input argument must contain only pointers to InferredPart.
	 * Method deletes any links that used the parts in the argument list as
	 * theirs reference point. 
	 */ 
	virtual void deleteAllReferences(const std::vector<IAttachableClass*>& reference_to);
};

/// @}


////////////////////////////////////////////////////////////////////////////////
//// Implementation for initial layer support of VocabularyPart


/// @addtogroup ext_vocabulary_tree
/// @{


/// @}

/// @}
/// @}


#endif /* _CORE_SUPPORT_EXTENSION_ */

