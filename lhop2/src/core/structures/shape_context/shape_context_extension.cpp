/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// shape_context_extension

#include "shape_context_extension.h"


////////////////////////////////////////////////////////////////////////////////
//// Belongie histogram caching as extension

void BelongieHistogramExtension::moveTo(IExtension& dst_extension) {
	BelongieHistogramExtension& dst_subpart_extension = dynamic_cast<BelongieHistogramExtension&>(dst_extension);

	// copy subparts_links to destination
	dst_subpart_extension.histogram_map.insert(histogram_map.begin(), histogram_map.end());

	// and just clear current array (do not delete any BelongieHistogramData pointers as they are still alive in the dst_extension) 
	histogram_map.clear();
}

void BelongieHistogramModifier::deleteAllReferences(const std::vector<IAttachableClass*>& reference_to) {
	// delete subpart 		

	// type should be InferredPart otherwise we cannot handle it
	for (auto iter = reference_to.begin(); iter != reference_to.end(); ++iter) {
		InferredPart& part = *(InferredPart*)*iter;

		// find this part in the internal map 
		auto found_part_iter = ext.histogram_map.find(part.getLocation());

		if (found_part_iter != ext.histogram_map.end()) {

			/// first delete part support data
			// delete any references to the subparts (if any class was using InferenceSupportData as reference)
			// NO NEED as it is not attachable !!
			// ext.holder.deleteAttachedReferences<BelongieHistogramData>(std::vector<IAttachableClass*>(found_part_iter->second.begin(), found_part_iter->second.end()));

			// then delete BelongieHistogramData pointer
			delete found_part_iter->second;

			// and remove part from internal list
			ext.histogram_map.erase(found_part_iter);
		}
	}
}

BelongieHistogramData* BelongieHistogramAccess::getHistogram(const IAttachableClass& part, InferenceTree& inference_tree) {
	// check if location is cached in histogram_map 

	// get_sc_map() and get_histogram() within layer_1_result.cpp

	// part has to be on 1. layer

	// get location of part
	
	// get 1. layer of parts from inference

	// get part within the region of your location

	// go through all locations that have at least one part
	// and create histogram from those parts for each subregion
		// you can use K_bin as histogram

	// cache histogram to histogram_map

	return nullptr;
}



////////////////////////////////////////////////////////////////////////////////
//// Shape context information attached to the VocabularyTree

void ShapeContextExtension::moveTo(IExtension& dst_extension) {
	ShapeContextExtension& dst_subpart_extension = dynamic_cast<ShapeContextExtension&>(dst_extension);

	// copy subparts_links to destination
	dst_subpart_extension.geometry_data.insert(geometry_data.begin(), geometry_data.end());

	// and just clear current array (do not delete any ShapeContextGeometryData pointers as they are still alive in the dst_extension) 
	geometry_data.clear();
}

AbstractSerializer::IFactory* ShapeContextExtension::getSerializer() const {
	return nullptr; // TODO !!!
}

void ShapeContextModifier::deleteAllReferences(const std::vector<IAttachableClass*>& reference_to) {
	// delete subpart 		

	// type should be InferredPart otherwise we cannot handle it
	for (auto iter = reference_to.begin(); iter != reference_to.end(); ++iter) {
		VocabularySubpartData& subpart = *(VocabularySubpartData*)*iter;

		// find this part in the internal map 
		auto found_part_iter = ext.geometry_data.find(subpart.getUUID());

		if (found_part_iter != ext.geometry_data.end()) {

			/// first delete part support data
			// delete any references to the subparts (if any class was using InferenceSupportData as reference)
			// NOT NEEDED as it is not attachable !!
			// ext.holder.deleteAttachedReferences<ShapeContextGeometryData>(std::vector<IAttachableClass*>(found_part_iter->second.begin(), found_part_iter->second.end()));

			// then delete ShapeContextGeometryData pointer
			delete found_part_iter->second;

			// and remove part from internal list
			ext.geometry_data.erase(found_part_iter);
		}
	}
}


void* ShapeContextExtension::calculateGeometryDistance(const ShapeContextGeometryData& geometry_a, const ShapeContextGeometryData& geometry_b, bool calc_benergy, int max_matching_size) {
	// TODO: need to implement as geometry_distance(..) and part_geometry_matching(..) from old code in part_learning.cpp!!!

}