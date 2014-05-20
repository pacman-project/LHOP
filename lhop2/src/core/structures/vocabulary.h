/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// vocabulary

#pragma once
#ifndef _CORE_VOCABULARY_
#define _CORE_VOCABULARY_

#include <vector>
#include <unordered_map>

//#include <opencv2/opencv.hpp>

#include "core/legacy/part_lib.h"

#include "core/structures/responses_array.h"
#include "core/structures/attached_extension.h"

/// @addtogroup core
/// @{
/// @addtogroup main_structures
/// @{
/// @addtogroup vocabulary_tree
/// @{

class VocabularyTree;
class VocabularyLayer;
class VocabularyPart;

class VocabularyTree : public ExtensionHolder {
private:
	// UUID generators
	UUIDGenerator uuid_generators_vocabulary_part; // for class InferredPart 

	//static double similarity_threshold = 1E-6;

	std::vector<VocabularyLayer*> layers_list;
public:
	VocabularyTree() : layers_list(0), uuid_generators_vocabulary_part(UUIDGenerator::newInstance<VocabularyPart>()) {
		std::cout << " calling constructor of VocabularyTree with layers list size: " << layers_list.size() << std::endl;
	}
	// TODO: write destructor but add it to the .cpp

    /// Returns the number of layers in the library.
    int getNumberOfLayers() { return layers_list.size(); }
    int getLastLayerIndex() { return layers_list.size() - 1; }

	// get reference to the specific layer of the vocabulary
	VocabularyLayer& getLayer(int layer) { return *layers_list[layer]; }
	const VocabularyLayer& getLayer(int layer) const { return *layers_list[layer]; }

	void insertNewLayer(VocabularyLayer* layer) { layers_list.push_back(layer); }

	UUIDGenerator& getVocabularyPartUUIDGenerator() { return uuid_generators_vocabulary_part; }
};

class VocabularyLayer {
private:
	VocabularyTree& vocabulary_tree;

	const double contraction_factor;
	const int layer_index; // TODO: how do we get this layer index info ??

	std::unordered_map<UUIDType,VocabularyPart*> parts;
public:
	VocabularyLayer(VocabularyTree& vocabulary_tree, double contraction_factor, int layer_index) : vocabulary_tree(vocabulary_tree), contraction_factor(contraction_factor), layer_index(layer_index) {
	}

	// TODO: write destructor but add it to the .cpp

	double getContrationFactor() { return contraction_factor; }
    int getNumberOfParts() { return parts.size(); } 
	int getLayerIndex() const { return layer_index;} 	

	VocabularyTree& getVocabularyTree() { return vocabulary_tree; }
	const VocabularyTree& getVocabularyTree() const { return vocabulary_tree; }

	VocabularyPart& getPart(UUIDType uuid);
	const VocabularyPart& getPart(UUIDType uuid) const;

	VocabularyPart& getPartByTypeId(int type_id);

	void insertNewPart(VocabularyPart* part);
};

class VocabularyPart : public IAttachableClass {
public:
	typedef map<ipoint2, double> mask_t;
    typedef set<ipoint2> region_t;
private:
	const int bpcount;  // number of "basic" parts; non-streamable, filled by get_basic_part_number()

	const int layer;
	const int type_id;
	ResponsesArray thresholds;

	// TODO: move this into the new extension
	// this should be renamed into reconstruction points/image
	// 
	// mask == original filter used to construct the part
    const mask_t mask;  // Non-serializable, except for layer 0, filled recursively by get_mask

	// this should be renamed into activation points i.e. points of the image where this part fired from (was activated from) 
	// (for first layer parts this can be area of the filter with locations of the non-zero values)
	// region == mask of the region filter i.e. locations of the region filter values with the positive factor
    const region_t region; // Non-serializable, except for layer 0, filled recursively by get_region

	const cv::Point2i relative_center_of_mass; // center of mass relative to the central node
public:
	VocabularyPart(UUIDType uuid, int layer, int type_id, mask_t mask, region_t region, const ResponsesArray& thresholds_ = ResponsesArray::DEFAULT_RESPONSES, cv::Point2i relative_center_of_mass = cv::Point2i(0))
		: IAttachableClass(uuid), layer(layer), type_id(type_id), mask(mask), region(region), thresholds(thresholds_), relative_center_of_mass(relative_center_of_mass), bpcount(0) {
	}
	int getTypeId() const { return type_id; }
	cv::Point2i getRelativeCenterOfMass() const { return relative_center_of_mass; }

	ResponsesArray::ValueType getResponseThreshold(ResponseType response_type) const { return thresholds.get(response_type); }
	void setResponseThreshold(ResponseType response_type, ResponsesArray::ValueType value) { thresholds.set(response_type,value); }


	int getBasicPartCount() const {
		// TODO: implement getBasicPartCount() in the VocabularyPart
		return bpcount;
	}
	
	
};

///////////////////////////////////////////////////////////////////////////
////// Scaffolding class for loading VocabularyTree from the old part_lib

#include "core/legacy/part_lib.h"

class VocabularyTreeFromPartLib : public VocabularyTree {
	part_lib* lib;
public:
	VocabularyTreeFromPartLib(part_lib* lib) : lib(lib), VocabularyTree() {
		initializeFromPartLib(lib);
	}

private:
	void initializeFromPartLib(part_lib* lib);
};

/// @}
/// @}
/// @}

#endif /* _CORE_VOCABULARY_ */

