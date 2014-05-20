#pragma once
#ifndef _CORE_OUTPUT_OBJECT__
#define _CORE_OUTPUT_OBJECT__

#include "utils/img.h"

#include "core/structures/parse_tree.h"
#include "core/structures/vocabulary.h"

#include "core/input_output/groupable_object.h"

/// @addtogroup core
/// @{
/// @addtogroup input_output
/// @{

/**
 * Base abstract class for any output object returned by inference, learning, toolset or any other modules
 */ 
class AbstractOutputObject : public GroupableObject {
public:
	virtual ~AbstractOutputObject() {} // add virtual destructor to force polymorphisem (otherwise we cannot use dynamic_cast<..>())
};


/**
 * Base output object for saving images. 
 * (Class may currently not be used anywhere)
 *
 * Implementation must provide method saveImage(const img*).
 */ 
class ImageOutputObject : public AbstractOutputObject {
protected:
	std::shared_ptr<img> image;
public:
	ImageOutputObject(std::shared_ptr<img> image, std::multimap<GroupableObject::Type, GroupableObject::Member> group) : image(image) {
		this->group = group;
	}
	ImageOutputObject(const ImageOutputObject& obj) : image(obj.image)  {
		this->group = obj.group;
	}

	
	std::shared_ptr<img> getImage() {
		return image;
	}
};


/**
 * Base output object for saving layer1_result (TODO: change name) object
 * Class is used by feature extraction and inference
 *
 * Implementation must provide method saveLayerObject(const layer1_result*).
 */
class LayerOutputObject : public AbstractOutputObject {
protected:
	std::shared_ptr<InferenceTree> ly_object;
public:
	LayerOutputObject(std::shared_ptr<InferenceTree> ly_object) : ly_object(ly_object) {
	}
	LayerOutputObject(std::shared_ptr<InferenceTree> ly_object, std::multimap<GroupableObject::Type, GroupableObject::Member> group) : ly_object(ly_object) {
		this->group = group;
	}

	LayerOutputObject(const LayerOutputObject& obj) : ly_object(obj.ly_object)  {
		this->group = obj.group;
	}

	std::shared_ptr<InferenceTree> getLayerObject() {
		return ly_object;
	}
};


/**
 * Base output object for saving vocabulary object (former part_lib).
 * Class is used mostly by learning and by other utility toolsets that display information about vocabulary.
 *
 * Implementation must provide method saveVocabularyObject(const part_lib*).
 */ 
class VocabularyOutputObject : public AbstractOutputObject {
protected:
	std::shared_ptr<VocabularyTree> vocabulary;
public:
	VocabularyOutputObject(std::shared_ptr<VocabularyTree> vocabulary) : vocabulary(vocabulary) {
	}
	VocabularyOutputObject(std::shared_ptr<VocabularyTree> vocabulary, std::multimap<GroupableObject::Type, GroupableObject::Member> group) : vocabulary(vocabulary) {
		this->group = group;
	}
	VocabularyOutputObject(const VocabularyOutputObject& obj) : vocabulary(obj.vocabulary)  {
		this->group = obj.group;
	}

	std::shared_ptr<VocabularyTree> getVocabularyObject() {
		return vocabulary;
	}
};


/// @}
/// @}

#endif /* _CORE_OUTPUT_OBJECT__ */
