/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// feature_extraction

#pragma once
#ifndef _CORE_FEATURE_EXTRACTION_
#define _CORE_FEATURE_EXTRACTION_

#include "utils/img.h"
#include "utils/graphs/img_graph.h"
#include "utils/matrix.h"
#include "utils/misc.h"
#include "utils/utils.h"
//#include "utils/config_dictionary.h"
#include "utils/configuration.h"

#include "utils/class_register.h"

#include "core/structures/parse_tree.h"
#include "core/structures/vocabulary.h"

#include "core/input_output/input_object.h"
#include "core/input_output/output_object.h"

#include "core/input_output/input_preprocessing.h"
#include "core/input_output/output_postprocessing.h"

/// @addtogroup core
/// @{
/// @addtogroup feature_extraction
/// @{

// Base classes for feature extraction 
//////////////////////////////////////////////

/**
 * A struct for feature type definition. Each feature type is defined by:
 *   - ID (that that starts with 0) 
 *        A unique number identifying this feature type 
 *   - feature type activation points (...):
 *        A list of offset points that are used to compute response of the activation at each location (x,y)
 *   - reconstruction image (img):
 *        An image representation of this feature type. This will be used for visualization of each activation/part.
 */ 
struct FeatureTypeDefinition {
	int id;
	matrix<bool> activation_points;
	img reconstruction_image;
};

/**
 * A representations of one activation/candidate for a final first layer part.
 * Must be as simple as possible due to memory requirements.
 */
struct FeatureActivationResponse {
	float value;
};

/**
 * Class holding an 2D array of FeatureActivationResponses (which currently holds only one float)
 * extracted from an input object for specific feature type.
 *
 * Any AbstractFeatureExtraction class MUST implement extractFeatureActivations and return one or more 
 * FeatureActivationArrays that correspond to response values of this specific feature type.
 * The type_id corresponds connects to the feature type used when extracting this activation array.
 *
 * This is in effect a wrapper (with iterators) over a simple 2D array of floats (i.e. 2D array of activations) 
 * plus definition of feature type
 */ 
class FeatureActivationArray {
	friend class AbstractFeatureExtraction;
protected:
	// corresponding feature type ID
	VocabularyPart* vocabulary_feature_part; // TODO: this should be VocabularyPart& but cannot use it due to default constructor for std::vector and due to struct IntermediateData
	
	int array_x_size;
	int array_y_size;

	// simple struct for 2D array of activations
	struct FeatureCandidatesArrays {
		std::vector<int> is_valid; // vector<bool> is memory optimized (only 1 bit per element)
		std::vector<FeatureActivationResponse> list; // 2D array containing list of features
		

		FeatureCandidatesArrays() {}
		// move assignment to avoid excessive copying
		FeatureCandidatesArrays(FeatureCandidatesArrays&& obj) : is_valid(std::move(obj.is_valid)), list(std::move(obj.list)) {}
		// copy assignment has to be explicitly declared when declaring also move operator
		FeatureCandidatesArrays(const FeatureCandidatesArrays& obj) : is_valid(obj.is_valid), list(obj.list) {}
	} activations_array;

public:	
	/**
	 * Simple location allows for efficient iteration over all elements of the feature type.
	 */
	class Location {
		friend class FeatureActivationArray;
		int index;
	public:
		Location() : index(0) {}
		Location(const FeatureActivationArray& feeature_type, int i) : index(i) {}
		Location(const FeatureActivationArray& feeature_type, int x, int y) : index(y * feeature_type.array_x_size + x) {}

		int getLocationX(const FeatureActivationArray& feeature_type) { return index % feeature_type.array_x_size; }
		int getLocationY(const FeatureActivationArray& feeature_type) { return index / feeature_type.array_x_size; }

		// overload for operation used in iteration
		inline bool operator==(const Location& b) { return index == b.index; }
		inline bool operator!=(const Location& b) { return index != b.index; }
		inline Location& operator++() { index++; return *this; }
		inline Location& operator--() { index--; return *this; }

		inline Location& operator+=(int i) { index+=i; return *this; }
		inline Location& operator-=(int i) { index-=i; return *this; }
	};

private:
	// simple end marker for iterators
	Location location_end;
public:	
	FeatureActivationArray() 
		: vocabulary_feature_part(nullptr), array_x_size(0), array_y_size(0), location_end(*this,0) { // default used only for std::vector
	} 
	FeatureActivationArray(VocabularyPart* vocabulary_feature_part, int array_x_size, int array_y_size) 
		: vocabulary_feature_part(vocabulary_feature_part), array_x_size(array_x_size), array_y_size(array_y_size), 
		  location_end(*this, array_x_size*array_y_size) {
		
		FeatureActivationResponse empty_feature; empty_feature.value = 0;
		
		activations_array.list = std::vector<FeatureActivationResponse>(array_x_size * array_y_size, empty_feature); // using memset might be more efficeint
		activations_array.is_valid = std::vector<int>(array_x_size * array_y_size, true);
	}

	// move assignment to avoid excessive copying
	FeatureActivationArray(FeatureActivationArray&& obj) 
		: vocabulary_feature_part(vocabulary_feature_part), array_x_size(obj.array_x_size), array_y_size(obj.array_y_size),
		location_end(obj.location_end), activations_array(std::move(obj.activations_array)) {			
	}
	// copy assignment has to be explicitly declared when declaring also move operator
	FeatureActivationArray(const FeatureActivationArray& obj) : vocabulary_feature_part(vocabulary_feature_part), array_x_size(obj.array_x_size), array_y_size(obj.array_y_size), location_end(obj.location_end), activations_array(obj.activations_array) {			
	}

	FeatureActivationResponse& getCandidateAt(const Location& loc) { return activations_array.list[loc.index]; }
	bool isCandiateAtValid(const Location& loc) { return activations_array.is_valid[loc.index]; }
	void setCandidateAtValid(const Location& loc, bool valid) { activations_array.is_valid[loc.index] = valid; }
	
	Location beginLocation() const { return Location(*this, 0); }
	const Location& endLocation() const { return location_end; }

	Location beginLocationAt(int i) const { return Location(*this, i); }
	Location beginLocationAt(int x, int y) const { return beginLocationAt(y * array_x_size + x); }


	std::vector<FeatureActivationResponse>::iterator beginFeatureIterator() { return activations_array.list.begin(); }
	std::vector<FeatureActivationResponse>::iterator endFeatureIterator() { return activations_array.list.end(); }

	std::vector<FeatureActivationResponse>::iterator beginFeatureIteratorAt(int i) { return beginFeatureIterator() + i; }
	std::vector<FeatureActivationResponse>::iterator beginFeatureIteratorAt(int x, int y) { return beginFeatureIteratorAt(y * array_x_size + x); }
	
	std::vector<int>::iterator beginIsValidIterator() { return activations_array.is_valid.begin(); }
	std::vector<int>::iterator endIsValidIterator() { return activations_array.is_valid.end(); }

	std::vector<int>::iterator beginIsValidIteratorAt(int i) { return beginIsValidIterator() + i; }
	std::vector<int>::iterator beginIsValidIteratorAt(int x, int y) { return beginIsValidIteratorAt(y * array_x_size + x); }

	int getWidth() const { return array_x_size; }
	int getHeight() const { return array_y_size; }

	VocabularyPart& getTypeId() const { return *vocabulary_feature_part; }

	// move assignment to avoid excessive copying (of std::vectors in feature_activations)
	FeatureActivationArray& operator=(FeatureActivationArray&& obj) {
		vocabulary_feature_part = obj.vocabulary_feature_part; array_x_size = obj.array_x_size; array_y_size = obj.array_y_size; location_end = obj.location_end;
		activations_array.is_valid = std::move(obj.activations_array.is_valid);
		activations_array.list = std::move(obj.activations_array.list);
		return *this;
	}
};

/**
 * Base class for any feature extraction. It contains 8 implemented methods for processing feature list of a FeatureActivationArray 
 * plus one main abstract method (extractFeatureActivations) that return initial list of FeatureActivationArray. Each FeatureActivationArray
 * corresponds to activations from one feature type and must be created in the implemented extractFeatureActivations().
 *
 * Current descended implementation:
 *  - AbstractFilterbankFeatureExtraction (still abstract)
 *      - GaborFeatureExtraction
 *      - AppGaborFeatureExtraction
 *      - DogGaborFeatureExtraction
 *      - LogGaborFeatureExtraction
 *      - ColorEdgesFeatureExtraction
 */ 
class AbstractFeatureExtraction : public IRegistrableClass { 
private:
	float layer1_threshold;
	float layer1_3x3bound;
	bool global_normalize;
	float power_correction;
	float response_percent;
	int border_size;

	struct IntermediateData {
		std::vector<FeatureActivationArray> feature_activations;

		// can be calculated by findMaxActivations
		FeatureActivationArray max_feature_activations;		
		float global_max_activations;

		// sum 2x2 of max activations
		// can be calculated by computeSum2x2MaxActivations
		FeatureActivationArray max_feature_activations_sum2x2;
		float global_max_activations_sum2x2;
	};

	std::vector<IInputPreprocessing*> preprocessing_steps;
	std::vector<IOutputPostprocessing*> postprocessing_steps;

	std::shared_ptr<VocabularyTree> vocabulary;
public:
    AbstractFeatureExtraction(const IConfiguration& config, 
					const std::vector<IInputPreprocessing*>& preprocessing_steps, 
					const std::vector<IOutputPostprocessing*>postprocessing_steps,
					std::shared_ptr<VocabularyTree> vocabulary = nullptr) 
		: preprocessing_steps(preprocessing_steps), postprocessing_steps(postprocessing_steps), vocabulary(vocabulary) {

		// load settings
		layer1_threshold = config.getDouble("layer1_threshold", 0.1);
		layer1_3x3bound = config.getDouble("layer1_3x3bound", 3);
		global_normalize = config.getBool("normalize", true);
		power_correction = config.getDouble("power_correction", 0.7);
		response_percent = config.getDouble("response_percent", 0.8);
		border_size =  config.getInt("border_size");
	}

	virtual ~AbstractFeatureExtraction() {
		for (auto iter = preprocessing_steps.begin(); iter != preprocessing_steps.end(); ++iter) delete *iter;
		for (auto iter = postprocessing_steps.begin(); iter != postprocessing_steps.end(); ++iter) delete *iter;
	}

	/**
	 * Entry point for feature extraction. Takes AbstractInputObject and creates LayerOutputObject with the initial layer.
	 * 
	 * Performs extraction in following steps:
	 *  1. Calls virtual method extractFeatureActivations() to retrieve list of part activations(candidates) for each type of feature 
	 *  2. Performs possible thresholding, inhibition, global or local normalization, etc. on activated parts
	 *  3. Promotes remaining activations into initial parts on first layer of layer1_result and returns them in a wrapper (i.e. LayerOutputObject)
	 */
	LayerOutputObject* performExtraction(const AbstractInputObject* input) {
		IntermediateData data;

		// Retrive list of part activations/candidates from abstract method extractFeatureActivations() (MUST be implemented by the descended class)
		// The goal of this method is to transform features of the input domain (shape, texture, color, motion etc.) 
		// into the LHOP part activations/candidates where for one location we can have multiple types of features
		data.feature_activations = extractFeatureActivations(input);
		
		clock_t start,end;
		
		// first for each location find the strongest feature type 
		// i.e. perform max operation over all different features at one location
		data.max_feature_activations = std::move(findMaxActivations(data.feature_activations, &data.global_max_activations));

		// also perform maximum operation over sum of 2x2
		data.max_feature_activations_sum2x2 = findSum2x2MaxActivations(data.max_feature_activations, &data.global_max_activations_sum2x2);
		
		// threshold the activations i.e. shutdown candidates/activations that are below global threshold
		applyGlobalTresholding(data, layer1_threshold);

		// perform local inhibition using max activations and sum2x2 of max activations i.e. shutdown candidates/activations whos neighbors are stronger		
		applyLocalInhibitionOfMax(data, layer1_3x3bound);

		// perform global normalization of each activation (in future this might be replaced with local normalization)
		if (normalize)
			applyGlobalNormalization(data);

		// perform power correction of each activation
		applyPowerCorrection(data, power_correction);

		// recalculate max features as theirs values changed due to the applied power correction
		data.max_feature_activations = std::move(findMaxActivations(data.feature_activations, &data.global_max_activations));
		data.max_feature_activations_sum2x2 = findSum2x2MaxActivations(data.max_feature_activations, &data.global_max_activations_sum2x2);

		// perform filtering of individual activations whose value is below threshold relative to maximal activation at the same location
		applyLocationFiltering(data, response_percent);

		// finally create layer1_result and promote remaining activations/candidates into first layer parts 
		InferenceTree* result = constructLayer(data, border_size);
		
		return new LayerOutputObject(std::shared_ptr<InferenceTree>(result));
	}
	
	/**
	 * Obtains pointer to initial vocabulary either by using any provided vocabulary from the constructor 
	 * or by generating new vocabulary based on the getFeatureTypeDefinitions() implemented by the 
	 * descended class.
	 */
	virtual std::shared_ptr<VocabularyTree> getInitialVocabulary();

	/**
	 * Performs all pre-processing steps on the input object.
	 */
	std::vector<AbstractInputObject*> performPreprocessing(const AbstractInputObject* input_object);

	/**
	 * Performs all post-processing steps on the list of output objects.
	 * Each output object is a pair of its input object and output data.
	 */
	std::vector<AbstractOutputObject*> performPostprocessing(const std::vector<AbstractOutputObject*>& output_object);

	/// Factory interface that must be implemented for any class that provides feature extraction extension
	class IFactory : public IRegistrableClassFactory {
	 public:
		virtual AbstractFeatureExtraction* newInstance(const IConfiguration& config) const = 0;		
		virtual string getAssociatedInputObjectType() const = 0;

		std::vector<IInputPreprocessing*> getPreprocessingSteps(const IConfiguration& config) const {
			std::vector<IInputPreprocessing*> steps;

			// resize
			int init_size = config.getInt("init_size", -100);

			if (init_size != -100)
				steps.push_back(new ResizeImagePreprocessing(init_size));

			// blur			
			float scale_sigma = config.getDouble("scale_sigma", 0.0);
			int scale_mask_size = (int)(5.0 * scale_sigma);

			if (scale_sigma > 0) 
				steps.push_back(new BlurImagePreprocessing(scale_mask_size, scale_sigma));

			// scaling
			int scale_limit = config.getInt("scale_limit", 200);
			int max_scales = config.getInt("max_scales", 99);
			float scale_factor = config.getDouble("scale_factor", 1/::pow(2.0, 1.0/3.0));
			
			if (scale_factor != 1)
				steps.push_back(new ScaleImagePreprocessing(scale_factor, scale_limit, max_scales));

			return steps;
		}
		std::vector<IOutputPostprocessing*> getPostprocessingSteps(const IConfiguration& config) const {
			return std::vector<IOutputPostprocessing*>();
		}
	};
protected:
	////////////////////////////////////////////////////////////////////////////////////
	////////////// Abstract methods required to be implemented //////////////////////////

	/**
     * Abstract methods that must be implemented by the descended class.
     *
     * Computes activations to all feature type used in the process of performExtraction().
     * Each element of the resulting vector will be one feature type that can have multiple activations 
     * spread over different locations.
     */
	virtual std::vector<FeatureActivationArray> extractFeatureActivations(const AbstractInputObject* input) = 0;	
	
	/**
     * Abstract methods that must be implemented by the descended class.
     *
     * Returns a list of definitions for all feature types that the implementation extracts.
     */
	virtual std::vector<FeatureTypeDefinition> getFeatureTypeDefinitions() = 0;	

	//////////////////////////////////////////////////////////////////////
	////////////// Processing/filtering methods //////////////////////////

	/**
	 * Searches for max activation on each location (x,y) over all types of features.
	 *
	 * Returns max activations as a new feature type (with type id -1) where value of a candidate/feature is max value at the same location
	 * over all feature types of the input array. Additionally returns global max activation value if global_max_activation_return ptr is not null.
	 */ 
	FeatureActivationArray findMaxActivations(std::vector<FeatureActivationArray>& feature_candidate, float* global_max_activation_return = nullptr);
	
	/**
	 * Method computes sum of 2x2 windows from max activations (assuming max activations is 2-dimensional array).
	 * Takes FeatureActivationArray of max activations as input (max_activations) and return the sum2x2 of max_activations. Additionally returs 
     * global max activation value of sum2x2 if global_max_activation_sum2x2_return ptr is not null.
	 * 
     * NOTICE: Returned array has undefined values at the edges (on all sides).
     *
	 * max_activations -> sum2x2 ->  max_activations_sum2x2
	 * | 0 1 ...|                    | x   x   ...| 
	 * | 2 3 ...| -> sum[0,1,2,3] -> | x [sum] ...|
	 * | ...    |                    | ...        |
     * 
	 */
	FeatureActivationArray findSum2x2MaxActivations(FeatureActivationArray& max_activations, float* global_max_activation_sum2x2_return);

	/**
	 * Uses data.max_feature_activations to shutdown all location whose max activation is below a threshold (same values as threshold are passed through).
	 * Threshold is calculated as layer1_threshold * global_max_activation_value.
	 *
	 * For efficiency it disables only max activations and then copies is_valid array to other feature types. In output all feature types will have identical
	 * is_valid array therefore any previous individual modification of feature type's is_valid array will be overridden.
	 */
	IntermediateData& applyGlobalTresholding(IntermediateData& data, float layer1_threshold);
	/** 
	 * Performs local inhibition (neighbor of 3x3) of activations based on data.max_feature_activations and data.max_feature_activations_sum2x2.
	 * Activation is not inhibited if number of stronger neighbors in 3x3 area of either max activations or either of max 2x2 sum is below threshold (layer1_3x3bound).
	 *
	 * For efficiency it disables only max activations and then copies is_valid array to other feature types. In output all feature types will have identical
	 * is_valid array therefore any previous individual modification of feature type's is_valid array will be overridden.
	 */
	IntermediateData& applyLocalInhibitionOfMax(IntermediateData& data, float layer1_3x3bound);

	/**
	 * Performs global normalization by dividing each activation value (over all location and all feature types) with global max activation value (data.global_max_activations):
	 *		feature[x,y] = feature[x,y] / max(over all features types)
	 *
     * Method will skip computation on location that are invalid (whose isCandiateAtValid() will return false).
     * Operation is done in-place and results is the same object as input.
     */ 
	IntermediateData& applyGlobalNormalization(IntermediateData& data);

    /**
	 * Performs power correction to each activation value (over all location and all feature types):
     *		feature[x,y] = feature[x,y] ^ power_correction_value
	 *
     * Method will skip computation on location that are invalid (whose isCandiateAtValid() will return false).
     * Operation is done in-place and results is the same object as input.
     */ 
	IntermediateData& applyPowerCorrection(IntermediateData& data, float power_correction);

	/**
	 * Shutdowns all activations if their activation value is below max value of activation at the same location (x,y) multiplied by response_percent factor:
     *		if feature[x,y] < max[all features at x,y] * response_percent then shutdown activation
     *
     * Method will skip computation on location that are invalid (whose isCandiateAtValid() will return false).
     * Operation is done in-place and results is the same object as input.
     */ 
	IntermediateData& applyLocationFiltering(IntermediateData& data, float response_percent);

	/** 
	 * Construct resulting InferenceTree* with the same width/height of the input image + added border.
     * Each remaining activation is used to initialize R_RESPONSE of an added part.
	 */
	InferenceTree* constructLayer(IntermediateData& data, int border);
};

/// @}
/// @}

#endif /* _CORE_FEATURE_EXTRACTION_ */
