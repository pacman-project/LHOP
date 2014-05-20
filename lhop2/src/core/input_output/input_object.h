#pragma once
#ifndef _CORE_INPUT_OBJECT__
#define _CORE_INPUT_OBJECT__

#include <memory>

#include "utils/img.h"
#include "core/legacy/layer_1_result.h"
#include "core/legacy/part_lib.h"

#include "core/structures/parse_tree.h"
#include "core/structures/vocabulary.h"

#include "core/input_output/groupable_object.h"

/// @addtogroup core
/// @{
/// @addtogroup input_output
/// @{

/**
 * Names used for the registration of each implementation that extends any AbstractInputObject.
 * e.g. class FileImageInputObject can register its factory with IMAGE_INPUT_OBJECT and e.g.
 * FilterbankFeatureExtraction can post a request for its input in IMAGE_INPUT_OBJECT.
 */
#define INPUT_OBJECT_WITH_GROUNDTRUTH_NAME "InputObjectWithGroundtruth"
#define IMAGE_INPUT_OBJECT "ImageInputObject"
#define LAYER_INPUT_OBJECT "LayerInputObject"
#define VOCABULARY_INPUT_OBJECT "VocabularyInputObject"

/**
 * Base abstract class for any input object required by inference, learning, toolset or any other modules
 * We define AbstractInputObject as groupable i.e. each input object can belonged to one or more groups. Using groups
 * allows us to split object into multiple input objects and to merge the results from split results.
 */ 
class AbstractInputObject : public GroupableObject {
public:
	virtual ~AbstractInputObject() {} // add virtual destructor to force polymorphisem (otherwise we cannot use dynamic_cast<..>())
};

// TODO: move somewhere into utilities
struct GroundtruthObject {
	irectangle2 region;
	std::string label;
};

/*
 * Base input object for retrieving groundtruth.
 * Class is extended by at least AbstractImageInputObject and AbstractLayerInputObject which have groundtruth associated with the object
 */ 
class InputObjectWithGroundtruth : public AbstractInputObject {
public:
	virtual std::vector<GroundtruthObject> getGroundtruths(const std::string& only_label = "") const = 0;
};

/**
 * Base input object for retrieving images. 
 * Class is used by feature extraction of the shape module (could also be used by texture module or color module or any other module using images as input).
 *
 * Implementation must provide method getImage() returing img* object.
 */ 
class AbstractImageInputObject : public virtual InputObjectWithGroundtruth {
public:
	virtual std::shared_ptr<img> getImage() const = 0;
};

/**
 * Simple image input object that holds image and groundtruths in memory as img class and as vector of GroundtruthObject
 */
class MemoryImageInputObject : public AbstractImageInputObject {
	const std::shared_ptr<img> image;
	std::vector<GroundtruthObject> groundtruths;
public:
	MemoryImageInputObject(std::shared_ptr<img> image, const std::vector<GroundtruthObject>& groundtruths, std::multimap<GroupableObject::Type, GroupableObject::Member> group) 
		:  image(image), groundtruths(groundtruths) { 
		this->group = group;
	}

	virtual std::shared_ptr<img> getImage() const { return image; }
	virtual std::vector<GroundtruthObject> getGroundtruths(const std::string& only_label = "") const { 
		if (only_label.empty()) {
			return groundtruths;
		} else {
			std::vector<GroundtruthObject> lable_groundtruths;
			
			for (auto iter = groundtruths.end(); iter != groundtruths.end(); ++iter)
				if (iter->label.compare(only_label) == 0)
					lable_groundtruths.push_back(*iter);

			return lable_groundtruths;
		}
	} 
};

/**
 * Base input object for retrieving InferenceTree input object
 * Class is used by the inference and learning.
 *
 * Implementation must provide method getLayerObject() returing smart pointer to the InferenceTree object.
 */
class AbstractLayerInputObject : public virtual InputObjectWithGroundtruth {
public:
	virtual std::shared_ptr<InferenceTree> getLayerObject() const = 0;
};


/**
 * Simple layer input object that holds inference tree and groundtruths in memory as img class and as vector of GroundtruthObject
 */
class MemoryLayerInputObject : public AbstractLayerInputObject {
	const std::shared_ptr<InferenceTree> ly_obj;
	std::vector<GroundtruthObject> groundtruths;
public:
	MemoryLayerInputObject(std::shared_ptr<InferenceTree> ly_obj, const std::vector<GroundtruthObject>& groundtruths, std::multimap<GroupableObject::Type, GroupableObject::Member> group) 
		:  ly_obj(ly_obj), groundtruths(groundtruths) { 
		this->group = group;
	}

	virtual std::shared_ptr<InferenceTree> getLayerObject() const { return ly_obj; }
	virtual std::vector<GroundtruthObject> getGroundtruths(const std::string& only_label = "") const { 
		if (only_label.empty()) {
			return groundtruths;
		} else {
			std::vector<GroundtruthObject> lable_groundtruths;
			
			for (auto iter = groundtruths.end(); iter != groundtruths.end(); ++iter)
				if (iter->label.compare(only_label) == 0)
					lable_groundtruths.push_back(*iter);

			return lable_groundtruths;
		}
	} 
};

/**
 * Base input object for retrieving vocabulary tree
 * Class is used mostly by learning and by other utility toolsets that display information about vocabulary.
 *
 * Implementation must provide method getVocabularyObject() returing smart pointer to the VocabularyTree object.
 */ 
class AbstractVocabularyInputObject : public AbstractInputObject {
public:
	virtual std::shared_ptr<VocabularyTree> getVocabularyObject() const = 0;
};

/// @}
/// @}

#endif /* _CORE_INPUT_OBJECT__ */
