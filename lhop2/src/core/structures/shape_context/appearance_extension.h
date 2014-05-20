/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// appearance_extension

#pragma once
#ifndef _CORE_APPEARANCE_EXTENSION_
#define _CORE_APPEARANCE_EXTENSION_

#include <typeinfo>
#include <typeindex>

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
/// @addtogroup ext_vocabulary_tree
/// @{

////////////////////////////////////////////////////////////////////////////////
//// Vocabulary subparts appearance and its access

class IVocabularyAppearanceAccess : public IExtensionAccess {
public:
	virtual map<UUIDType, pair<vector<UUIDType>, double>> getSubpartAppearance(const VocabularySubpartData& subpart_data) = 0;
};

class IVocabularyAppearanceModifier : public IExtensionModifier {
public:
	// TODO:
};

// forward definitions
class VocabularyAppearanceAccess;
class VocabularyAppearanceModifier;
class VocabularyAppearanceExtension;


class VocabularyAppearanceExtension : public SimpleAbstractExtension<VocabularyAppearanceExtension,
																		VocabularyAppearanceAccess,
																		VocabularyAppearanceModifier,
																		VocabularySubpartData> {
public:
	ExtensionHolder& holder;
	
	// TODO: implement appearance extension

	VocabularyAppearanceExtension(ExtensionHolder& holder) : holder(holder) {
	}

	virtual void moveTo(IExtension& dst_extension);

	virtual AbstractSerializer::IFactory* getSerializer() const;
};


////////////////////////////////////////////////
/// Access class for appearance of the VocabularyTree

class VocabularyAppearanceAccess : public IVocabularyAppearanceAccess {
	VocabularyAppearanceExtension& ext;
public:
	VocabularyAppearanceAccess(VocabularyAppearanceExtension& ext) : ext(ext) {
	}

	virtual map<UUIDType, pair<vector<UUIDType>, double>> getSubpartAppearance(const VocabularySubpartData& subpart_data) {
		// TODO
		map<UUIDType, pair<vector<UUIDType>, double>> idenity;
		idenity[subpart_data.subpart.getUUID()] = std::pair<vector<UUIDType>, double>(vector<UUIDType>(), 1);
		return idenity;
	}
};

////////////////////////////////////////////////
/// Modifier class for appearance of the VocabularyTree

class VocabularyAppearanceModifier : public IVocabularyAppearanceModifier {
	VocabularyAppearanceExtension& ext;
public:
	VocabularyAppearanceModifier(VocabularyAppearanceExtension& ext) : ext(ext) {
	}
	
	virtual void deleteAllReferences(const std::vector<IAttachableClass*>& reference_to) {
		// TODO
	}

};


/// @}
/// @}
/// @}
/// @}

#endif /* _CORE_APPEARANCE_EXTENSION_ */

