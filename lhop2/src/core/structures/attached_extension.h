// attached_extension

#pragma once
#ifndef _CORE_ATTCHED_EXTENSION_
#define _CORE_ATTCHED_EXTENSION_

#include <map>
#include <vector>
#include <set>
#include <typeinfo>
#include <typeindex>

#include "utils/uuid.h"
#include "utils/serialization/serialization.h"

/// @addtogroup core
/// @{
/// @addtogroup main_structures
/// @{
/// @addtogroup ext
/// @{

/**
 * Extensions are a way to augment the vocabulary/inference trees (or any other objects that allows it)
 * with the additional information. It is used to provide e.g. subpart links
 * to each part, to implement part indexing in the vocabulary, to add
 * subpart geometry, to add appearance models or to extract layer 1 support parts
 *
 * The main extension is implemented from the IExtension interface
 * and any object that enables the use of extensions must extend the ExtensionHolder class.
 *
 * Individual extension must expose an access layer to anybody using it through the
 * interface IExtensionAccess. The extended IExtension class must return appropriate
 * implementation from the getAccess() method.
 *
 * Example usage:
 *  VocabularyTree vocabulary(...);
 *  VocabularyPart* vocabulary_part = ...;
 *
 *  // get access to the vocabulary indexing
 *  VocabularyIndexingAccess indexing_access = vocabulary.getAccess<VocabularyIndexingAccess>(); 
 *
 *  // obtained indexed parts for the provided vocabulary part
 *  std::vector<VocabularyPart*> indexed_parts = indexing_access.getIndexedCentralParts(vocabulary_part);
 */

// forward definitions
class DefaultExtensionHolderSerializer;
class ExtensionHolder;

class ExtensionNotPresentException : public libhop_exception {
public:
	ExtensionNotPresentException(const string& file, const size_t line, const string& msg = string()): libhop_exception(file, line, msg, "Extension not present") {}
	ExtensionNotPresentException(const string& file, const size_t line, const std::exception& from_ex): libhop_exception(file, line, from_ex, "Extension not present") {}
};


/**
 * Interface for any class that will have extension attached to itself.
 * Each attachable class must have an unique identifier (uuid) that can be used as to/from pointer reference.
 */ 
class IAttachableClass {
	UUIDType uuid;
public:
	IAttachableClass(UUIDType uuid) : uuid(uuid) {}

	UUIDType getUUID() const { return uuid; }
};

/**
 * Interface for class that provides access to the extension.
 * Make this class as simple as possible with only acting as an wrapper to the 
 * actual extension.
 * NOTICE: all methods of this class can be used in performance critical code
 * so avoid any expensive calls (virtual calls) and minimize any overhead.
 */
class IExtensionAccess {
public:
	virtual ~IExtensionAccess() {} // force it into polymorphic type
};

/**
 * Interface for class that is able to modify IExtension.
 * IExtensionAccess should only provide simple access methods
 * while IExtensionModifier should provide access to modifications,
 * such as deletion, insertion etc.
 * By default it provides method for deletion of references and is called
 * when ExtensionHolder->deleteAttachedReferences() is being invoked.
 */ 
class IExtensionModifier : public IExtensionAccess {
public:
	/**
	 * Method must implement deletion of its elements when the reference is deleted from outside of the extension.
	 * This will be called only if IExtension registered (through IExtension::getAttachableClassesTypeInfo()) 
     * that it uses specific IAttachableClass class as reference. The method can safely cast elements of reference_to 
     * to that class (or must determine by itself correct class if it registered to multiple IAttachableClass classes).
	 */ 
	virtual void deleteAllReferences(const std::vector<IAttachableClass*>& reference_to) = 0;
};

/**
 * Interface for any extension class. It provides an access class by returning 
 * IExtensionAccess implementation in the getAccess() call and a modifier class by
 * returning IExtensionModifier implementation in the getModifier() call.
 */
class IExtension : public ISerializableObject {
public:
	virtual IExtensionAccess* getAccess() = 0;
	virtual IExtensionModifier* getModifier() = 0;

	// IExtension must provide typeid of its access and modifier implementation class and
	// typeids of all IAttachableClass classes that it can attach to.
	virtual std::type_index getAccessTypeInfo() = 0;
	virtual std::type_index getModifierTypeInfo() = 0;
	virtual std::vector<std::type_index> getAttachableClassesTypeInfo() = 0;

	virtual void moveTo(IExtension& dst_extension) = 0;

	class IFactory {
	public:
		virtual IExtension* newInstance(ExtensionHolder& holder) const = 0;
	};
};

/**
 * Holder of all extension attached to the desired class.
 * Each extension inserted is registered by typeid(*extension->getAccess()) name.
 * I.e. it uses the name of its access class.
 */ 
class ExtensionHolder : public ISerializableObject {
	// grant access to abstract serializer that will expose private methods as protected getters
	friend class AbstractExtensionHolderSerializer;
private:
	std::set<IExtension*> extensions;

	// maps: typeid(TFunctionalityAccess).hash_code --> IExtension*
	std::map<size_t, IExtensionAccess*> extension_accesses;
	
	// maps: type(IAttachableClass).hash_code -> array of IExtension*
	std::map<size_t, std::vector<IExtension*>> extension_references;
public:	

	virtual ~ExtensionHolder() {
		// TODO: delete all extensions and their memory
	}

	/**
	 * Inserts new extension into the holder. Will return true if inserted 
	 * or false if one with the same access/modifier classes already exists
	 */ 
	bool insertExtension(IExtension* ext) {
		
		size_t access_type_id = ext->getAccessTypeInfo().hash_code();
		size_t modifier_type_id = ext->getModifierTypeInfo().hash_code();

		// check if access or modifier class already exists
		if (extension_accesses.find(access_type_id) != extension_accesses.end() || 
			extension_accesses.find(modifier_type_id) != extension_accesses.end()) {
			// TODO: what to do if only one type exist but the other does not ?? (this should be an indicator of an error ?)
			return false;
		}

		// save extension
		extensions.insert(ext);

		// register extension access and modifier based on its hashcode of typeid
		extension_accesses[access_type_id] = ext->getAccess();
		extension_accesses[modifier_type_id] = ext->getModifier();

		// register extension to the list of possible attachable classes based on its hashcode of typeid
		std::vector<std::type_index> attachable_classes_types = ext->getAttachableClassesTypeInfo();

		for (auto iter = attachable_classes_types.begin(); iter != attachable_classes_types.end(); ++iter) {
			extension_references[iter->hash_code()].push_back(ext);
		}

		return true;
	}

	/**
	 * Get access to the extension associated with provided access class.
	 * Extension must already be present in the class otherwise this method will throw an exception.
	 */ 
	template <class TExtensionAccess>
	TExtensionAccess getAccess() const {
		// get hash code of the extension access
		size_t hash_code = typeid(TExtensionAccess).hash_code();

		// find extension instance in the internal list
		auto extension_iter = extension_accesses.find(hash_code);
		if (extension_iter == extension_accesses.end()) {
			throw custom_libhop_exception(ExtensionNotPresentException,"Cannot find requested extension access - extension is not present");
		}

		// and return a copy of a value
		return *dynamic_cast<TExtensionAccess*>(extension_iter->second);
	}
	/**
	 * Get access to the extension associated with provided access class.
	 * If extension is not present in the class this method will use provided factory to create it.
	 */ 
	template <class TExtensionAccess>
	TExtensionAccess getAccess(const IExtension::IFactory& factory) {
		try {
			return getAccess<TExtensionAccess>();
		} catch (const ExtensionNotPresentException& ex) {
			// create extension using provided factory and insert it
			insertExtension(factory.newInstance(*this));

			// try getting access again
			return getAccess<TExtensionAccess>();
		}		
	}

	/**
     * Get list of IExtensionModifier* of all extensions that can be attached to the TAttachableClass.
	 * When IExtension is inserted it must also provide a list of all TAttachableClass classes that it will
	 * attach to. This function find the inverse of that i.e. it finds all the extensions that use specific 
     * TAttachableClass class and returns theirs modifiers
	 */
	template <class TAttachableClass>
	std::vector<IExtensionModifier*> getAttachableReferences() {
		std::vector<IExtensionModifier*> result;

		// get hash code of the attachable access
		size_t hash_code = typeid(TAttachableClass).hash_code();

		auto extension_ref_iter = extension_references.find(hash_code);
		if (extension_ref_iter != extension_references.end()) {

			for (auto iter = extension_ref_iter->second.begin(); iter != extension_ref_iter->second.end(); ++iter) {
				result.push_back((*iter)->getModifier());
			}
		}
		return result;
	}
	
	/**
     * Method for deletion of references that point to specific IAttachableClass class.
     * Since IAttachableClass does not know about any pointer to it, it must call deleteAttachedReferences()
     * that will handle deletion of any references that point to the IAttachableClass being deleted.
     * I.e. function calls deleteAllReferences(objects_for_deletion) on any modifier returned by getAttachableReferences()
     */
	template <class TAttachableClass> 
	void deleteAttachedReferences(const std::vector<IAttachableClass*>& objects_for_deletion) {
		// first find if all functionalities that can be attached to this part (get theirs modifier interfaces)
		std::vector<IExtensionModifier*> attached_extension_modifiers = getAttachableReferences<TAttachableClass>();

		// make sure any extension using this part as a reference will delete its reference
		for (auto iter = attached_extension_modifiers.begin(); iter != attached_extension_modifiers.end(); ++iter) {
			(*iter)->deleteAllReferences(objects_for_deletion);
		}
	}

	/**
	 * Moves all extensions to a new ExtensionHolder provided in as the argument.
	 * Existing IExtension* pointers become invalid as they they are merged into the
	 * functionalities of the dst.
	 * Internal lists are cleared and all other functions will not return any valid objects.
	 */ 
	void moveExtensionsTo(ExtensionHolder& dst) {
		
		for (auto iter = extensions.begin(); iter != extensions.end(); ++iter) {
			IExtension* extension = *iter;

			// merge or add new extension
			auto existing_extension = dst.extensions.find(extension);
			if (existing_extension == dst.extensions.end()) {
				dst.insertExtension(extension);
			} else {
				// just merge it with the existing extension
				extension->moveTo(**existing_extension); // existing_extension is IExtension* and moveTo() uses IExtension&
				
				// we can delete current extension
				delete extension;
			}
		}

		// don't do anything with the extension references as dst already contains correct pointers 
		// (if the same extension was already present) or was inserted with the insertExtension() call
		// just clear the lists
		extensions.clear();
		extension_accesses.clear();
		extension_references.clear();
	}

	virtual AbstractSerializer::IFactory* getSerializer() const;
};


/**
 * Convenience class as template for simple extension that uses 
 * specific classes as defined by template parameter:
 *  - TAccessClass (IExtensionAccess): constructor MUST take one parameter
 *  - TModifierClass (IExtensionModifier)
 *  - TAttachableClass (IAttachableClass)
 * By default it does not use any serializer (returns nullptr as serializer factory)
 */
template <class TExtension, class TAccessClass, class TModifierClass, class TAttachableClass> 
class SimpleAbstractExtension : public IExtension {
public:
	virtual IExtensionAccess* getAccess() { return new TAccessClass((TExtension&)*this); }
	virtual IExtensionModifier* getModifier() { return new TModifierClass((TExtension&)*this); }

	// IExtension must provide typeid of its access and modifier implementation class and
	// typeids of all IAttachableClass classes that it can attach to.
	virtual std::type_index getAccessTypeInfo() { return type_index(typeid(TAccessClass)); }
	virtual std::type_index getModifierTypeInfo() { return type_index(typeid(TModifierClass)); }
	virtual std::vector<std::type_index> getAttachableClassesTypeInfo() {
  		std::vector<std::type_index> result;
		result.push_back(type_index(typeid(TAttachableClass)));
		return result;
	};
	virtual AbstractSerializer::IFactory* getSerializer() const { return nullptr;}
};

/// @}
/// @}
/// @}

#endif /* _CORE_ATTCHED_EXTENSION_ */

