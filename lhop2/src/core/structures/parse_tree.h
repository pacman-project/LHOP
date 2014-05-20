/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// parse_tree

#pragma once
#ifndef _CORE_PARSE_TREE_
#define _CORE_PARSE_TREE_

//#include <opencv2/opencv.hpp>

#include "core/structures/responses_array.h"
#include "core/structures/attached_extension.h"
#include "core/structures/vocabulary.h"

#include "utils/uuid.h"
#include "utils/serialization/serialization.h"



/// @addtogroup core
/// @{
/// @addtogroup main_structures
/// @{
/// @addtogroup inference_tree Inference tree
/// @{


// forward definitions
class InferenceTree;
class InferenceLayer;
class InferredPart;

/**
 * An InferenceTree class (aka parse tree or layer1_result) contains a list of parts inferred from
 * an input object (e.g. an image or a video sequence) divided into separate layers.
 * 
 * For each layer we hold an InferenceLayer class that contains information about parts.
 */
class InferenceTree : public ExtensionHolder {
	// UUID generators
	UUIDGenerator uuid_generators_inferred_part; // for class InferredPart 

	std::vector<InferenceLayer*> layers_list;
public:
	InferenceTree() : layers_list(0), uuid_generators_inferred_part(UUIDGenerator::newInstance<InferredPart>()) {		
	}
	virtual ~InferenceTree();

	InferenceLayer& getLayer(int layer) { return *layers_list[layer]; }
	const InferenceLayer& getLayer(int layer) const { return *layers_list[layer]; }

	// Get layer that was inferred from specific vocabulary layer.
	InferenceLayer& getLayer(const VocabularyLayer& vocabulary_layer);
	// Get layer that was inferred from specific vocabulary layer.
	const InferenceLayer& getLayer(const VocabularyLayer& vocabulary_layer) const;

	void insertNewLayer(InferenceLayer* layer_parts) {
		layers_list.push_back(layer_parts);
	}

	int getNumberOfLayers() const { return layers_list.size(); }

	/**
	 * TODO: implement if needed (should work the same as layer1_result::merge())
	 * what is the difference between InferenceTree::addExistingTree() ??
	 */
	void mergeFrom(InferenceTree& inference_tree) {
		// TODO: implement merging
		// this should substitute layer1_result::merge()
		throw new_libhop_exception("Merging of InferenceTree objects not yet implemented in the InferenceTree");
	}
	/**
	 * TODO: implement if needed (should work the same as layer1_result::add_result())
	 * what is the difference between InferenceTree::mergeFrom() ??
	 */
	void addExistingTree(std::vector<InferenceTree*> exising_inference_trees) {
		// TODO: implement merging
		// this should substitute layer1_result::add_result() 
		throw new_libhop_exception("Merging of InferenceTree objects not yet implemented in the InferenceTree");
	}	

	UUIDGenerator& getInferredPartUUIDGenerator() {
		return uuid_generators_inferred_part;
	}

	const UUIDGenerator& getInferredPartUUIDGenerator() const {
		return uuid_generators_inferred_part;
	}

	virtual AbstractSerializer::IFactory* getSerializer() const;
};

/**
 * An InferenceLayer class contains a list of parts assigned to only specific layer.
 * We contain parts in 2D array and hold multiple parts per each location.
 */
class InferenceLayer {
private:
	InferenceTree& inference_tree;

	const int layer_index;

	const cv::Size2i size;

	std::vector<std::vector<InferredPart*>> part_array; // TODO: change into unordered map for memory efficiency

public:
	class Iterator {
		friend class InferenceLayer;
		int index;
	public:
		Iterator() : index(0) {}
		Iterator(const InferenceLayer& layer, int i) : index(i) {}
		Iterator(const InferenceLayer& layer, int x, int y) : index(y * layer.size.width + x) {}
		Iterator(const InferenceLayer& layer, const cv::Point2i& loc) : index(loc.y * layer.size.width + loc.x) {}

		int getLocationX(const InferenceLayer& layer) { return index % layer.size.width; }
		int getLocationY(const InferenceLayer& layer) { return index / layer.size.width; }

		// overload for operation used in iteration
		inline bool operator==(const Iterator& b) { return index == b.index; }
		inline bool operator!=(const Iterator& b) { return index != b.index; }
		inline Iterator& operator++() { index++; return *this; }
		inline Iterator& operator--() { index--; return *this; }

		inline Iterator& operator+=(int i) { index+=i; return *this; }
		inline Iterator& operator-=(int i) { index-=i; return *this; }

	};
private:
	Iterator iterator_end;
public:
	InferenceLayer(InferenceTree& inference_tree, cv::Size2i size, int layer_index) 
		: inference_tree(inference_tree), size(size), layer_index(layer_index), part_array(size.width * size.height), iterator_end(*this, size.width * size.height) {
	}

	virtual ~InferenceLayer();

	InferenceTree& getInferenceTree() { return inference_tree; }
	const InferenceTree& getInferenceTree() const { return inference_tree; }

	cv::Size2i getSize() const { return size; }	
	int getLayerIndex() const { return layer_index;} 

	void insertNewPart(InferredPart* part);

	std::vector<InferredPart*>& getPartsAt(Iterator& iter);
	const std::vector<InferredPart*>& getPartsAt(Iterator& iter) const;

	std::vector<InferredPart*> getPartsInRegion(const cv::Rect region) const;

	void deleteParts(const std::vector<InferredPart*>& parts_for_deletion);

	Iterator beginIterator() const { return Iterator(); }
	Iterator endIterator() const { return iterator_end; }

	Iterator beginIteratorAt(int x, int y) const { return Iterator(*this, x, y); }
};


/**
 * An InferredPart class is a representation of one inferred part from its input object.
 * It is defined by its layer number, 2D array location (x,y), and its type ID. It also contains
 * a list of responses defined in PartResponses class.
 *
 * Additionally, multiple IExtension classes provide additional extension that can be accessed
 * through theirs access interfaces:
 *  - InferenceSubpartsAccess
 *  - InferenceFirstLayerSubpartsAccess
 */ 
class InferredPart : public IAttachableClass {
private:
	const cv::Point2i location; // x,y
	const int layer; // is this needed here since we know layer from the InferenceLayer ?	
	ResponsesArray responses; // R_RESPONSE, G_RESPONSE, RR_RESPONSE ...

	const UUIDType corresponding_vocabulary_part_uuid;
public:
	InferredPart(UUIDType uuid, int layer, cv::Point2i location, UUIDType corresponding_vocabulary_part_uuid, const ResponsesArray& responses = ResponsesArray::DEFAULT_RESPONSES)
		: IAttachableClass(uuid), layer(layer), location(location), corresponding_vocabulary_part_uuid(corresponding_vocabulary_part_uuid), responses(responses) {
		// TODO: generate uuid based on internal counter ??
	}
	int getLayer() const { return layer; }
	cv::Point2i getLocation() const { return location; }
	
	const ResponsesArray& getResponseArray() const { return responses; }

	ResponsesArray::ValueType getResponse(ResponseType response_type) const { return responses.get(response_type); }
	void setResponse(ResponseType response_type, double value) { responses.set(response_type, value); }

	// vval() comes from the old code and is defined as G_RESPONSE_1 at first layer and R_RESPONSE_1 at other layers
	ResponsesArray::ValueType getVVal() const {
		return layer != 0 ? getResponse(R_RESPONSE_1) : getResponse(G_RESPONSE_1);
	}

	UUIDType getCorrespondingVocabularyPartUUID() const {
		return corresponding_vocabulary_part_uuid;
	}
	const VocabularyPart& getCorrespondingVocabularyPart(const VocabularyTree& vocabulary) {
		return vocabulary.getLayer(layer).getPart(corresponding_vocabulary_part_uuid);
	}

	
};


///////////////////////////////////////////////////////////////////////////
////// Scaffolding class for loading ParseTree from the old layer1_result

#include "core/legacy/layer_1_result.h"

class ParseTreeFromLayer1Result : public InferenceTree {
	layer1_result* res;
public:
	ParseTreeFromLayer1Result(layer1_result* res) : res(res) {}
};

/// @}
/// @}
/// @}


#endif /* _CORE_PARSE_TREE_ */

