/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// parse_tree


#include "parse_tree.h"

#include "core/structures/serializers.h"

InferenceTree::~InferenceTree() {
	for (auto iter = layers_list.begin(); iter != layers_list.end(); ++iter) delete *iter;
}

InferenceLayer& InferenceTree::getLayer(const VocabularyLayer& vocabulary_layer) { 
	return getLayer(vocabulary_layer.getLayerIndex()); 
}
const InferenceLayer& InferenceTree::getLayer(const VocabularyLayer& vocabulary_layer) const { 
	return getLayer(vocabulary_layer.getLayerIndex()); 
}

InferenceLayer::~InferenceLayer() {
	for (auto iter_loc = part_array.begin(); iter_loc != part_array.end(); ++iter_loc) {
		for (auto iter_part = iter_loc->begin(); iter_part != iter_loc->end(); ++iter_part) {
			delete *iter_part;
		}
	}
}

void InferenceLayer::insertNewPart(InferredPart* part) {
	cv::Point2i loc = part->getLocation();
	int i = loc.y * this->size.width + loc.x;
	part_array[i].push_back(part);
}

std::vector<InferredPart*>& InferenceLayer::getPartsAt(InferenceLayer::Iterator& iter) {
	return part_array[iter.index];
}

const std::vector<InferredPart*>& InferenceLayer::getPartsAt(InferenceLayer::Iterator& iter) const {
	return part_array[iter.index];

}

std::vector<InferredPart*> InferenceLayer::getPartsInRegion(const cv::Rect region) const {
	std::vector<InferredPart*> part_list;

	int num_skiped_row_parts = this->size.width - (int)region.width;
	num_skiped_row_parts -= 1; // decrement by one, since ++region_iterator will handle that part

	InferenceLayer::Iterator region_iterator = this->beginIteratorAt(region.x, region.y);

	for (int j = 0; j < region.height; ++j, region_iterator+=num_skiped_row_parts) {
		for (int i = 0; i < region.width; ++i, ++region_iterator) {
			// get all parts at this location 
			std::vector<InferredPart*> all_parts_at_location = this->getPartsAt(region_iterator);

			// and copy pointers to returned array (part_list = [part_list + all_parts_at_location])
			part_list.insert(part_list.end(), all_parts_at_location.begin(), all_parts_at_location.end());
			
		}
	}

	return part_list;
}

void InferenceLayer::deleteParts(const std::vector<InferredPart*>& parts_for_deletion) {
	// delete any references that are attached to this part
	const std::vector<IAttachableClass*> parts_for_deletion_as_references(parts_for_deletion.begin(), parts_for_deletion.end());

	inference_tree.deleteAttachedReferences<InferredPart>(parts_for_deletion_as_references);

	// delete parts from the layer
	for (auto iter = parts_for_deletion.begin(); iter != parts_for_deletion.end(); ++iter) {
		InferredPart* part = *iter;
			
		// find storage for the location of current part
		auto loc_iter = InferenceLayer::Iterator(*this, part->getLocation());

		std::vector<InferredPart*>& parts_at = getPartsAt(loc_iter);
		
		// erase part from the storage
		std::remove_if(parts_at.begin(), parts_at.end(), 
						// lambda function: current_part == part
						[part](InferredPart* current_part) {
							return current_part == part;
						});
	}

}

AbstractSerializer::IFactory* InferenceTree::getSerializer() const {
	return new DefaultInferenceTreeSerializer::Factory();
}
