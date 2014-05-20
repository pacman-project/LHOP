#pragma once
#ifndef _CORE_OUTPUT_POSTPROCESSING__
#define _CORE_OUTPUT_POSTPROCESSING__

#include <vector>

#include "core/input_output/output_object.h"

/// @addtogroup core
/// @{
/// @addtogroup input_output
/// @{
/// @addtogroup postprocessing
/// @{


/**
 * Main inteface for any post-processing.
 *
 * Classes must implement method vector<AbstractInputObject*> doPostprocessing(vector<AbstractInputObject*>) that takes a list of input objects
 * and returns a new list of preprocessed input objects. Classes MUST NOT delete any original input AbstractInputObject* objects.
 */ 
class IOutputPostprocessing {
public:
	virtual std::vector<AbstractOutputObject*> doPostprocessing(const std::vector<AbstractOutputObject*>& output_list) const = 0 ;
};


/**
 * Abstract base class (interface) for any postprocessing on image itself.
 */
class ImagePostprocessing : public IOutputPostprocessing {
	
};

/**
 * Abstract base class for any postprocessing on layer object.
 */
class LayerPostprocessing : public IOutputPostprocessing {
protected:
	/**
	 * Casting to LayerOutputObject* needed for any layer postprocessing.
	 */
	LayerOutputObject* castToLayerObject(AbstractOutputObject* input_object) const {
		LayerOutputObject* output_layer = dynamic_cast<LayerOutputObject*>(input_object);

		if (output_layer == nullptr)
			throw new_libhop_exception("LayerPostprocessing can only preprocess LayerOutputObject* objects");

		return output_layer;
	}
};


/**
 * Combine all layer files belonging to the same group usig the layer1_result->add_results()
 */ 
class CombineLayerPartsPostprocessing : public LayerPostprocessing {
private:
	std::string group_type;
public:
	CombineLayerPartsPostprocessing(const std::string& group_type = "tile") : group_type(group_type) {}

	virtual std::vector<AbstractOutputObject*> doPostprocessing(const std::vector<AbstractOutputObject*>& output_list) const;
};

/**
 * Merge all layer files belonging to the same group usig the layer1_result->merge()
 */ 
class MergeScalesPostprocessing : public LayerPostprocessing {
private:
	std::string group_type;
public:
	MergeScalesPostprocessing(const std::string& group_type = "scale") : group_type(group_type) {}

	virtual std::vector<AbstractOutputObject*> doPostprocessing(const std::vector<AbstractOutputObject*>& output_list) const;
};


/**
 * Abstract base class for any postprocessing on vocabulary object.
 */
class VocabularyPostprocessing : public IOutputPostprocessing {
};

/// @}
/// @}
/// @{

#endif /* _CORE_OUTPUT_POSTPROCESSING__ */