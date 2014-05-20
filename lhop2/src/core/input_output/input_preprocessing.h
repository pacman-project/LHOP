#pragma once
#ifndef _CORE_INPUT_PREPROCESSING__
#define _CORE_INPUT_PREPROCESSING__

#include <vector>

#include "core/input_output/input_object.h"

/// @addtogroup core
/// @{
/// @addtogroup input_output
/// @{
/// @addtogroup preprocessing
/// @{

/**
 * Main interface for any preprocessing.
 *
 * Classes must implement method vector<AbstractInputObject*> doPreprocessing(vector<AbstractInputObject*>) that takes a list of input objects
 * and returns a new list of preprocessed input objects. Classes MUST NOT delete any original input AbstractInputObject* objects.
 */ 
class IInputPreprocessing {
public:
	virtual std::vector<AbstractInputObject*> doPreprocessing(const std::vector<AbstractInputObject*> input_list) const = 0 ;
};

/**
 * Abstract base class (interface) for any preprocessing on groundtruths.
 */
class GroundruthPreprocessing : public IInputPreprocessing {
protected:
	/**
	 * Casting to InputObjectWithGroundtruth* needed for any groundtruh preprocessing.
	 */
	InputObjectWithGroundtruth* castToImageObject(AbstractInputObject* input_object) const {
		InputObjectWithGroundtruth* input_grountruth = dynamic_cast<AbstractImageInputObject*>(input_object);

		if (input_grountruth == nullptr)
			throw new_libhop_exception("GroundruthPreprocessing can only preprocess InputObjectWithGroundtruth* objects");

		return input_grountruth;
	}
};
/**
 * Abstract base class (interface) for any preprocessing on image itself.
 */
class ImagePreprocessing : public IInputPreprocessing {
protected:
	/**
	 * Casting to AbstractImageInputObject* needed for any image preprocessing.
	 */
	AbstractImageInputObject* castToImageObject(AbstractInputObject* input_object) const {
		AbstractImageInputObject* input_image = dynamic_cast<AbstractImageInputObject*>(input_object);

		if (input_image == nullptr)
			throw new_libhop_exception("ImagePreprocessing can only preprocess AbstractImageInputObject* objects");

		return input_image;
	}
};

/**
 * Performs image resize on each input. Resize value is defined as:
 *	if resize_value > 0:
 *   - then resize_value defines max side in pixels
 *  if resize_value < 0:
 *   - then resize_value defines factor (multiplied by 100) of scaling (i.e. -120 => scale factor of 1.2)
 *
 * Requires AbstractImageInputObject in input list.
 */ 
class ResizeImagePreprocessing : public ImagePreprocessing {
	int resize_value;
public:
	ResizeImagePreprocessing(int resize_value) : resize_value(resize_value) {}

	virtual std::vector<AbstractInputObject*> doPreprocessing(const std::vector<AbstractInputObject*> input_list) const;
};

/**
 * Performs bluring on each input image.
 *
 * Requires AbstractImageInputObject in input list
 */
class BlurImagePreprocessing : public ImagePreprocessing {
	int scale_mask_size;
	double scale_sigma;
public:
	BlurImagePreprocessing(int scale_mask_size, double scale_sigma) : scale_mask_size(scale_mask_size), scale_sigma(scale_sigma) {}

	virtual std::vector<AbstractInputObject*> doPreprocessing(const std::vector<AbstractInputObject*> input_list) const;
};

/**
 * Performs scaling of image using scale_factor until scale_limit in pixels is reached (min side of image must not be below this limit)
 * or until max number scales have been produced. As soon as one or the other argument is reached the scaling stops.
 *
 * Requires AbstractImageInputObject in input list.
 */ 
class ScaleImagePreprocessing : public ImagePreprocessing {
private:
	int scale_limit;
    int max_scales;

    double scale_factor;

public:
	ScaleImagePreprocessing(double scale_factor, int scale_limit, int max_scales) : scale_factor(scale_factor), scale_limit(scale_limit), max_scales(max_scales) {}

	virtual std::vector<AbstractInputObject*> doPreprocessing(const std::vector<AbstractInputObject*> input_list) const;
};

/**
 * Splits colored image in seperate channels and returns new input for each image channel (original image can be discarded).
 */
class SeperateImageColorsPreprocessing : public ImagePreprocessing {
public:
	// each channel will be assigned to one of those groups (for recombining)
	enum { GROUP_CHANNEL_R, GROUP_CHANNEL_G, GROUP_CHANNEL_B };

	virtual std::vector<AbstractInputObject*> doPreprocessing(const std::vector<AbstractInputObject*> input_list) const;
};

/**
 * TODO: not implemented but it should apply bit-mask to an input image (TODO: where do we get bitmask ?)
 */
class ApplyImageMaskPreprocessing : public ImagePreprocessing {
public:
	virtual std::vector<AbstractInputObject*> doPreprocessing(const std::vector<AbstractInputObject*> input_list) const;
};

/**
 * Horizontaly flipes the image and their coresponding groundtruths.
 */
class FlipImagePreprocessing : public ImagePreprocessing {
public:
	virtual std::vector<AbstractInputObject*> doPreprocessing(const std::vector<AbstractInputObject*> input_list) const;
};

/**
 * Based on groundtruth regions attached to the image it extracts and returns only specific regions of the image (original image can be discarded).
 * Can also extract only based on specific groundtruth label (if set to empty string "" it will take all regions regardless of its label)
 */
class ExtractImageGroundtruthPreprocessing : public ImagePreprocessing  {
private:
	std::string only_label;
public:
	ExtractImageGroundtruthPreprocessing(const std::string& only_label = "") : only_label(only_label) {}

	virtual std::vector<AbstractInputObject*> doPreprocessing(const std::vector<AbstractInputObject*> input_list) const;
};

/**
 * Splits image into tiles for easier processing (in case of big images).
 * 
 * Spliting is controled by paramters n, overlep and rthresh, where n defines n x n tiles,
 * overlap defines in number of overlaping pixels in neighboring tiles and rthresh defines
 * minimal groundtruth-vs-tile intersection/union for gt to be included into specific tile 
 * (uses PASCAL way of calcuating intersection/union).
 */
class SplitImagePreprocessing : public ImagePreprocessing {
private:
	int n; int overlap;
	double rthresh;
public:
	SplitImagePreprocessing(int n, int overlap, double rthresh) : n(n), overlap(overlap), rthresh(rthresh) {}

	virtual std::vector<AbstractInputObject*> doPreprocessing(const std::vector<AbstractInputObject*> input_list) const;
};

/**
 * Abstract base class for any preprocessing on layer object.
 */
class LayerPreprocessing : public IInputPreprocessing {
};

/**
 * Abstract base class for any preprocessing on vocabulary object.
 */
class VocabularyPreprocessing : public IInputPreprocessing {
};

/// @}
/// @}
/// @}

#endif /* _CORE_INPUT_PREPROCESSING__ */
