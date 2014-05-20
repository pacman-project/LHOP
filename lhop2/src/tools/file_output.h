#pragma once
#ifndef _TOOLSET_FILE_OUTPUT__
#define _TOOLSET_FILE_OUTPUT__

#include <sstream>

#include "core/input_output/output_object.h"

#include "utils/class_register.h"

/// @addtogroup toolset Toolset
/// @{

/// @addtogroup file_input_output File input/output structures
/// @{


class IFileOutput : public IRegistrableClass  {
public:
	virtual void save(const std::string& basename) = 0;
};

/**
 * Abstract IFileOutput class for any groupable object (i.e. for any object that extends GroupableObject)
 */ 
class FileGroupableOutputObject : public IFileOutput {
protected:
	/**
	 * Create additional extension naming convention from groups of the object.
	 * Used group name and group member id to create string from template "_[group name]_[group member id]".
	 * Template is added for each group including for multiple same type groups.
	 * E.g. having groups: (scale,1), (scale_0), (tile,4) will produce extension: _scale_1_scale_0_tile_4
	 */ 
	std::string getGroupsExtensions(GroupableObject* group_object) {
		
		std::multimap<GroupableObject::Type, GroupableObject::Member> group = group_object->getGroupMap();

		std::stringstream ss;
		for (auto iter = group.begin(); iter != group.end(); ++iter) {
			ss << "_" << iter->first.name << "_" << iter->second.group_member_id;
		}

		return ss.str();
	}
};

/**
 * Class for saving images to a file. 
 */ 
class FileImageOutputObject : public FileGroupableOutputObject {
private:
	ImageOutputObject* image_output;
public:
	FileImageOutputObject(ImageOutputObject* image_output) : image_output(image_output) { 
		if (image_output == nullptr)
			throw new_libhop_exception("Cannon use null pointer in constructor of FileImageOutputObject");
	}

	virtual void save(const std::string& basename) {
		// set .png as default output
		std::string extension = getGroupsExtensions(image_output) + ".png";

		// replace whole extention with grouping name and with .png extentions
		string output_filename = change_extension(basename,extension);

		// get image
		std::shared_ptr<img> image = image_output->getImage();
		
		// save image
		if (image != nullptr)
			image->save(output_filename);
		//else TODO: should we throw an exception or just silently ignore ???
	}
};


#include <sstream>
#include <fstream>
#include "core/structures/serializers_text.h"

/**
 * Class for saving layer1_result (TODO: rename) to the file
 */
class FileLayerOutputObject : public FileGroupableOutputObject {
private:
	LayerOutputObject* layer_output;
public:
	FileLayerOutputObject(LayerOutputObject* layer_output) : layer_output(layer_output) {
		if (layer_output == nullptr)
			throw new_libhop_exception("Cannon use null pointer in constructor of FileLayerOutputObject");
	}

	virtual void save(const std::string& basename) {
		
		// get layer
		std::shared_ptr<InferenceTree> layer = layer_output->getLayerObject();

		// set .ly1 as default output
		std::string extension = getGroupsExtensions(layer_output) + ".ly" + layer->getNumberOfLayers();

		// replace whole extention with grouping name and with .png extentions
		std::string output_filename = change_extension(basename,extension);
				
		ofstreamer out_file;

		out_file.open(change_extension(basename,getGroupsExtensions(layer_output) + "_text_ly" + layer->getNumberOfLayers()));
		TextInferenceTreeSerializer(std::shared_ptr<SerializedObjectStorage>(new SerializedObjectStorage())).serialize(*layer, out_file);
		out_file.close();

		ofstreamer out_file_1;
		out_file_1.open(output_filename);
		DefaultInferenceTreeSerializer(std::shared_ptr<SerializedObjectStorage>(new SerializedObjectStorage())).serialize(*layer, out_file_1);
		out_file_1.close();

		ifstreamer in_file;
		in_file.open(output_filename);
		InferenceTree* deserialized_parse_tree = DefaultInferenceTreeSerializer(std::shared_ptr<SerializedObjectStorage>(new SerializedObjectStorage())).deserialize(in_file);
		in_file.close();

		ofstreamer out_file_2;
		out_file_2.open(change_extension(basename,getGroupsExtensions(layer_output) + "_text_unserialized_ly" + layer->getNumberOfLayers()));
		TextInferenceTreeSerializer(std::shared_ptr<SerializedObjectStorage>(new SerializedObjectStorage())).serialize(*deserialized_parse_tree, out_file_2);
		out_file_2.close();
		
		// save layer
		// TODO: implement saving
		//if (layer != nullptr)
		//	layer->save(output_filename);
		//else TODO: should we throw an exception or just silently ignore ???
	}
};


/**
 * Class for saving vocabulary (former part_lib) to the file
 */ 
class FileVocabularyOutputObject : public FileGroupableOutputObject {
private:
	VocabularyOutputObject* vocabulary_output;
public:
	FileVocabularyOutputObject(VocabularyOutputObject* vocabulary_output) : vocabulary_output(vocabulary_output) {
		if (vocabulary_output == nullptr)
			throw new_libhop_exception("Cannon use null pointer in constructor of FileVocabularyOutputObject");
	}

	virtual void save(const std::string& basename) {
		// set .plb as default output
		std::string extension = getGroupsExtensions(vocabulary_output) + ".plb";

		// replace whole extention with grouping name and with .png extentions
		string output_filename = change_extension(basename,extension);

		// get vocabulary 
		std::shared_ptr<VocabularyTree> vocabulary = vocabulary_output->getVocabularyObject();

		// save vocabulary
		// TODO: implement saving
		// if (vocabulary != nullptr)
		//	vocabulary->save(output_filename);
		//else TODO: should we throw an exception or just silently ignore ???
	}
};

/// @}
/// @}


#endif /* _TOOLSET_FILE_OUTPUT__ */
