#pragma once
#ifndef _TOOLSET_FILE_INPUT__
#define _TOOLSET_FILE_INPUT__

#include "core/input_output/input_object.h"

#include "utils/class_register.h"

#include "core/structures/serializers.h"


/// @addtogroup toolset Toolset
/// @{

/// @addtogroup file_input_output File input/output structures
/// @{


/**
 * Base class for any input object that can be loaded from file
 */ 
class IFileInput : public IRegistrableClass  {
public:
	virtual ~IFileInput() {} // add virtual destructor to force polymorphisem (otherwise we cannot use dynamic_cast<..>())

	/// Factory interface that must be implemented for any class that provides file loading functionality
	class IFactory : public IRegistrableClassFactory {
	 public:
		 virtual IFileInput* newInstance(const std::string filename) const = 0;
	};
};


/**
 * Class for loading groundtruth from a file.
 *
 * Two constructors are available:
 *  - FileInputObjectWithGroundtruth(const std::string& groundtruth_filename): will read from groundtruth_filename
 *  - FileInputObjectWithGroundtruth(const std::string& filename, const std::string& extension): will changed extentsion of filename with new extension and read from constructed filename
 */ 
class FileInputObjectWithGroundtruth : public virtual InputObjectWithGroundtruth, public IFileInput {
private:
	std::string groundtruth_filename;
public:
	FileInputObjectWithGroundtruth(const std::string& groundtruth_filename) : groundtruth_filename(groundtruth_filename) {}
	FileInputObjectWithGroundtruth(const std::string& filename, const std::string& extension) : groundtruth_filename(change_extension(filename, extension)) {}

	/**
	 * Reads groundtruth file.
	 */ 
	virtual std::vector<GroundtruthObject> getGroundtruths(const std::string& only_label = "") const;

	class Factory : public IFactory {
	 public:
		 virtual IFileInput* newInstance(const std::string filename) const { 
			 return new FileInputObjectWithGroundtruth(filename); 			 
		 }
		 virtual string assignedRegistrationName() const { return INPUT_OBJECT_WITH_GROUNDTRUTH_NAME; }
	};
};
/**
 * Class for loading images from a file. 
 * Can also read associated groundtruth file with .groundtruth extension.
 */ 
class FileImageInputObject : public FileInputObjectWithGroundtruth, public AbstractImageInputObject {
private:
	std::string image_filename;
public:
	FileImageInputObject(const std::string& image_filename, const std::string& gt_ext) 
		: image_filename(image_filename), FileInputObjectWithGroundtruth(image_filename, gt_ext) {
	}

	/**
	 * Reads image from filename using img::read_colors
	 */
	virtual std::shared_ptr<img> getImage() const {
		img* image = new img();
		img::read_colors(image_filename, *image);
		
		image->normalize(0,1);

		return std::shared_ptr<img>(image);
	}

	class Factory : public IFactory {
	 public:
		 virtual IFileInput* newInstance(const std::string filename) const { 
			 return newInstance(filename, ".groundtruth");
		 }
		 virtual IFileInput* newInstance(const std::string filename, const std::string& gt_ext) const { 
			 return new FileImageInputObject(filename, gt_ext); 
		 }
		 virtual string assignedRegistrationName() const { return IMAGE_INPUT_OBJECT; }
	};
};


/**
 * Class for loading layer (layer1_result) object from a file. 
 * Can also read associated groundtruth file with .groundtruth extension.
 */
class FileLayerInputObject : public FileInputObjectWithGroundtruth, public AbstractLayerInputObject {
private:
	std::string layer_filename;
public:
	FileLayerInputObject(const std::string& layer_filename, const std::string& gt_ext) 
		: layer_filename(layer_filename), FileInputObjectWithGroundtruth(layer_filename, gt_ext) {
	}

	virtual std::shared_ptr<InferenceTree> getLayerObject() const {
		// legacy: TODO need to change
		/*layer1_result* res;
		try {
			res = (layer1_result*)streamable::read(layer_filename);
		} catch (...) {
			res = nullptr;
		}*/
		
		ifstreamer in_file;
		in_file.open(layer_filename);

		InferenceTree* layer_deserialized = DefaultInferenceTreeSerializer(std::shared_ptr<SerializedObjectStorage>(new SerializedObjectStorage())).deserialize(in_file);

		return std::shared_ptr<InferenceTree>(layer_deserialized);
	}

	class Factory : public IFactory {
	 public:
		 virtual IFileInput* newInstance(const std::string filename) const { 
			 return newInstance(filename, ".groundtruth");
		 }
		 virtual IFileInput* newInstance(const std::string filename, const std::string& gt_ext) const { 
			 return new FileLayerInputObject(filename, gt_ext); 
		 }
		 virtual string assignedRegistrationName() const { return LAYER_INPUT_OBJECT; }
	};
};


/**
 * Class for loading vocabulary (part_lib) object from a file. 
 */ 
class FileVocabularyInputObject : public IFileInput, public AbstractVocabularyInputObject {
private:
	std::string vocabulary_filename;
public:
	FileVocabularyInputObject(const std::string& vocabulary_filename) : vocabulary_filename(vocabulary_filename){}

	virtual std::shared_ptr<VocabularyTree> getVocabularyObject() const {
		// legacy: TODO need to change
		part_lib* vocabulary;
		try {
			vocabulary = (part_lib*)streamable::read(vocabulary_filename);
		} catch (...) {
			vocabulary = nullptr;
		}
		return std::shared_ptr<VocabularyTree>(new VocabularyTreeFromPartLib(vocabulary));
	}

	class Factory : public IFactory {
	 public:
		 virtual IFileInput* newInstance(const std::string filename) const {  return new FileVocabularyInputObject(filename); }
		 virtual string assignedRegistrationName() const { return VOCABULARY_INPUT_OBJECT; }
	};
};

/// @}

/// @}

#endif /* _TOOLSET_FILE_INPUT__ */
