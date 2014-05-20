/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// subparts_extension

#pragma once
#ifndef _CORE_SUBPARTS_EXTENSION_
#define _CORE_SUBPARTS_EXTENSION_

#include <typeinfo>
#include <typeindex>
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

////////////////////////////////////////////////////////////////////////////////
//// Inference subparts extension

/// @addtogroup ext_inference_tree
/// @{

// forward definitions
class IInferenceSubpartsAccess;
class IInferenceSubpartsModifier;
class InferenceSubpartsAccess;
class InferenceSubpartsModifier;
class InferenceSubpartsExtension;
class InferenceSubpartData;

class IInferenceSubpartsAccess : public IExtensionAccess {
public:
	virtual std::vector<InferenceSubpartData*> getSubparts(const IAttachableClass& part) = 0;		
};

class IInferenceSubpartsModifier : public IExtensionModifier {
public:
	virtual void insertSubparts(const IAttachableClass& part, const std::vector<InferenceSubpartData*>& subparts_for_insertion) = 0;
	virtual void deleteSubparts(const IAttachableClass& part, const std::vector<InferenceSubpartData*>& subparts_for_deletion) = 0;
	
	// get reference to the UUID generator for InferrenceSubpartData class
	virtual UUIDGenerator& getInferrenceSubpartUUIDGenerator() = 0;
};


class IInferenceFirstLayerSubpartsAccess : public IExtensionAccess {
public:
	virtual std::vector<InferenceSubpartData*> getSubparts(const InferredPart& part) = 0;
};

/**
 * Main structure saved in the InferenceSubpartsExtension. It points to another inferred part
 * and represents a link to subpart with defined index (corresponding to vocabulary subpart index)
 * and response r (== response of matching (SchurProductBlock) which is usually schur product with the distribution map)
 */ 
struct InferenceSubpartData : public IAttachableClass {
	InferredPart& subpart;
	int index;
	double r;

	InferenceSubpartData(UUIDType uuid, InferredPart& subpart, int index, double r) : IAttachableClass(uuid), subpart(subpart), index(index), r(r) {
	}
};

////////////////////////////////////////////////////////////////////////////////
//// Inference subparts implementation

/**
 * Extension for associating one InferredPart as an subpart of another InferredPart.
 * It is attached to the InferenceTree, specifically it uses InferredPart as its reference.
 * The data stored within this connection is InferenceSubpartData. Access is provide through
 * InferenceSubpartsAccess class (read-only access) and InferenceSubpartsModifier class (write-only access).
 */
class InferenceSubpartsExtension : public SimpleAbstractExtension<InferenceSubpartsExtension,
																InferenceSubpartsAccess,
																InferenceSubpartsModifier,
																InferredPart> {
public:
	// UUID generators
	UUIDGenerator uuid_generators_inference_subpart_data; // for class InferenceSubpartData 

	ExtensionHolder& holder;

	// we map InferredPart* -> std::vector<InferenceSubpartData*>
	// where we use InferredPart->getUUID() as key identifier
	std::unordered_map<UUIDType, std::vector<InferenceSubpartData*>> subparts_links;

	// default constructor
	InferenceSubpartsExtension(ExtensionHolder& holder) 
		: holder(holder), uuid_generators_inference_subpart_data(UUIDGenerator::newInstance<InferenceSubpartData>()) {
	}
	
	// moves all data into dst_extension and clear this one
	virtual void moveTo(IExtension& dst_extension);
	
	virtual AbstractSerializer::IFactory* getSerializer() const;
	
	class Factory : public IFactory {
	public:
		virtual IExtension* newInstance(ExtensionHolder& holder) const { return new InferenceSubpartsExtension(holder); }
	};
};

////////////////////////////////////////////////
/// Access class for subparts of InferenceTree

/**
 * Class provides access to the subparts stored within the InferenceTree object.
 * Access is read-only from the getSubparts(..)
 */ 
class InferenceSubpartsAccess : public IInferenceSubpartsAccess {	
	InferenceSubpartsExtension& ext;
public:
	InferenceSubpartsAccess(InferenceSubpartsExtension& ext) : ext(ext) {
	}

	/**
	 * Retrieves list of links/subparts associated with the specific input part.
	 */
	virtual std::vector<InferenceSubpartData*> getSubparts(const IAttachableClass& part) {
		return ext.subparts_links[part.getUUID()];
	}
};

////////////////////////////////////////////////
/// Modifier class for subparts of InferenceTree

/**
 * Class provided modification access to the subparts stored within the InferenceTree object.
 * Access is write-only from the methods: insertSubparts(..) and deleteSubparts(..)
 */
class InferenceSubpartsModifier : public IInferenceSubpartsModifier {
	InferenceSubpartsExtension& ext;
public:
	InferenceSubpartsModifier(InferenceSubpartsExtension& ext) : ext(ext) {
	}
	
	/**
	 * Input argument must contain only pointers to InferredPart.
	 * Method deletes any links that used the parts in the argument list as
	 * theirs reference point. 
	 */ 
	virtual void deleteAllReferences(const std::vector<IAttachableClass*>& reference_to);

	/**
	 * Creates links between input part and provided subpart. Each subpart must point to another part and must
	 * contain the necessary data. 
	 */
	virtual void insertSubparts(const IAttachableClass& part, const std::vector<InferenceSubpartData*>& subparts_for_insertion);
	/**
	 * Removes provided links associated with the input part.
	 */
	virtual void deleteSubparts(const IAttachableClass& part, const std::vector<InferenceSubpartData*>& subparts_for_deletion);

	virtual UUIDGenerator& getInferrenceSubpartUUIDGenerator();
};

/// @}


////////////////////////////////////////////////////////////////////////////////
//// Vocabulary subparts and its access

/// @addtogroup ext_vocabulary_tree
/// @{


// forward definitions
class IVocabularySubpartsAccess;
class IVocabularySubpartsModifier;
class VocabularySubpartsAccess;
class VocabularySubpartsModifier;
class VocabularySubpartsExtension;
class VocabularySubpartData;


class IVocabularySubpartsAccess : public IExtensionAccess {
public:
	virtual std::vector<VocabularySubpartData*> getSubparts(const IAttachableClass& part) = 0;
};

class IVocabularySubpartsModifier : public IExtensionModifier {
public:
	virtual void insertSubparts(const IAttachableClass& part, const std::vector<VocabularySubpartData*>& subparts_for_insertion) = 0;
	virtual void deleteSubparts(const IAttachableClass& part, const std::vector<VocabularySubpartData*>& subparts_for_deletion) = 0;
	
	// get reference to the UUID generator for VocabularySubpartData class
	virtual UUIDGenerator& getVocabularySubpartUUIDGenerator() = 0;
};

/**
 * Main structure saved by the VocabularySubpartsFunctionality. It points to another vocabulary part
 * and represents a link to subpart with defined index, distribution map, normal g distribution and 
 * location point (x,y) which defines offset from the main central subpart.
 * 
 * TODO: we use all structure information only about non-central subparts
 */ 
struct VocabularySubpartData : public IAttachableClass {
	const VocabularyPart& subpart;
	int index;
	cv::Mat map_distribution;
	normal_distribution1 gdistr; // TODO: check which values are required
	
	cv::Point2i offset;

	VocabularySubpartData(UUIDType uuid, const VocabularyPart& subpart, int index, cv::Mat map_distribution, cv::Point2i offset, normal_distribution1 gdistr)
		: IAttachableClass(uuid), subpart(subpart), index(index), map_distribution(map_distribution), offset(offset), gdistr(gdistr)  { 
	}
};

/**
 * Functionality for associating one VocabularyPart as an subpart of another VocabularyPart.
 * It is usually attached to the VocabularyTree, specifically it uses VocabularyPart as its reference.
 * The data stored within this connection is VocabularySubpartData. Access is provide through
 * VocabularySubpartsAccess class (read-only access) and VocabularySubpartsModifier class (write-only access).
 */
class VocabularySubpartsExtension : public SimpleAbstractExtension<VocabularySubpartsExtension,
																	VocabularySubpartsAccess,
																	VocabularySubpartsModifier,
																	VocabularyPart>  {
public:
	// UUID generators
	UUIDGenerator uuid_generators_vocabulary_subpart_data; // for class VocabularySubpartData 

	ExtensionHolder& holder;

	// we map VocabularyPart* -> std::vector<VocabularySubpartData*>
	// where we use VocabularyPart->getUUID() as key identifier
	std::unordered_map<UUIDType, std::vector<VocabularySubpartData*>> subparts_links;

	// default constructor
	VocabularySubpartsExtension(ExtensionHolder& holder) 
		: holder(holder), uuid_generators_vocabulary_subpart_data(UUIDGenerator::newInstance<VocabularySubpartData>()) { 
	}
	
	// moves all data into dst_extension and clear this one
	virtual void moveTo(IExtension& dst_extension); 

	virtual AbstractSerializer::IFactory* getSerializer() const {
		return nullptr; // TODO: return serializer 
	}

	class Factory : public IFactory {
	public:
		virtual IExtension* newInstance(ExtensionHolder& holder) const { return new VocabularySubpartsExtension(holder); }
	};
};

////////////////////////////////////////////////
/// Access class for subparts of InferenceTree

/**
 * Class provides access to the subparts stored within the VocabularyTree object.
 * Access is read-only from the getSubparts(..)
 */ 
class VocabularySubpartsAccess : public IVocabularySubpartsAccess {
	VocabularySubpartsExtension& ext;
public:
	VocabularySubpartsAccess(VocabularySubpartsExtension& ext) : ext(ext) {
	}

	virtual std::vector<VocabularySubpartData*> getSubparts(const IAttachableClass& part) {
		return ext.subparts_links[part.getUUID()];
	}
};

////////////////////////////////////////////////
/// Modifier class for subparts of InferenceTree

/**
 * Class provided modification access to the subparts stored within the VocabularyTree object.
 * Access is write-only from the methods: insertSubparts(..) and deleteSubparts(..)
 */
class VocabularySubpartsModifier : public IVocabularySubpartsModifier {
	VocabularySubpartsExtension& ext;
public:
	VocabularySubpartsModifier(VocabularySubpartsExtension& ext) : ext(ext) {
	}
	
	virtual void deleteAllReferences(const std::vector<IAttachableClass*>& reference_to);

	/**
	 * Creates links between input part and provided subpart. Each subpart must point to another part and must
	 * contain the necessary data. 
	 */
	virtual void insertSubparts(const IAttachableClass& part, const std::vector<VocabularySubpartData*>& subparts_for_insertion);
	/**
	 * Removes provided links associated with the input part.
	 */
	virtual void deleteSubparts(const IAttachableClass& part, const std::vector<VocabularySubpartData*>& subparts_for_deletion);

	/** 
	 * Get UUID generator for VocabularySubpartData
	 */
	virtual UUIDGenerator& getVocabularySubpartUUIDGenerator();
};

/// @}

/// @}
/// @}


#endif /* _CORE_SUBPARTS_EXTENSION_ */

