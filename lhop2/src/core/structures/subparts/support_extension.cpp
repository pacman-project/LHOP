/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// support_extension

#include "support_extension.h"

#include "core/structures/serializers.h"

////////////////////////////////////////////////////////////////////////////////
//// Initial layer (first layer) support extension

void InferenceSupportExtension::moveTo(IExtension& dst_extension) {
	InferenceSupportExtension& dst_subpart_extension = dynamic_cast<InferenceSupportExtension&>(dst_extension);

	// copy subparts_links to destination
	dst_subpart_extension.support_subparts_links.insert(support_subparts_links.begin(), support_subparts_links.end());

	// and just clear current array (do not delete any InferenceSubpartData pointers as they are still alive in the dst_extension) 
	support_subparts_links.clear();
}

void InferenceSupportModifier::deleteAllReferences(const std::vector<IAttachableClass*>& reference_to) {
	// delete subpart 		

	// type should be InferredPart otherwise we cannot handle it
	for (auto iter = reference_to.begin(); iter != reference_to.end(); ++iter) {
		InferredPart& part = *(InferredPart*)*iter;
	
		// find this part in the internal map 
		auto found_part_iter = ext.support_subparts_links.find(part.getUUID());

		if (found_part_iter != ext.support_subparts_links.end()) {
		
			/// first delete part support data
			// delete any references to the subparts (if any class was using InferenceSupportData as reference)
			ext.holder.deleteAttachedReferences<InferenceSupportData>(std::vector<IAttachableClass*>(found_part_iter->second.begin(), found_part_iter->second.end()));
			
			// then delete all InferenceSupportData pointer
			for (auto iter = found_part_iter->second.begin(); iter != found_part_iter->second.end(); ++iter) {
				delete *iter;
			}

			// and remove part from internal list
			ext.support_subparts_links.erase(found_part_iter);
		}
	}
}

const std::vector<InferenceSupportData*>& InferenceSupportAccess::getAllInitialLayerSupportParts(const IAttachableClass& part) {
	
	std::vector<InferenceSupportData*>& existing_support_parts = ext.support_subparts_links[part.getUUID()];

	if (existing_support_parts.empty()) {
		InferredPart& inferred_part = (InferredPart&)part;

		// fill support parts from previous layer parts
		std::vector<InferenceSubpartData*> subpart_list = subparts_access.getSubparts(inferred_part);

		if (inferred_part.getLayer() <= 0) {
			// find initial layer support for each subpart and add them to the support list
			for (auto iter = subpart_list.begin(); iter != subpart_list.begin(); ++iter) {
				InferenceSubpartData* subpart_data = *iter;
				// recursively get support for subpart 
				std::vector<InferenceSupportData*> subpart_support_list = this->getAllInitialLayerSupportParts(*subpart_data);
			
				// and merge it into the existing_support_parts list as new InferenceSupportData
				for (auto support_iter = subpart_support_list.begin(); support_iter != subpart_support_list.end(); ++support_iter) {
					InferenceSupportData* support_subpart = *support_iter;

					// add UUID of current 
					std::vector<UUIDType> support_subpart_leaf_to_root_path = support_subpart->support_path;
					support_subpart_leaf_to_root_path.push_back(subpart_data->getUUID());

					// add to list of supporting subparts
					// supporting subpart will also bes saved into the extension.support_subparts_links 
					// for future accesses since existing_support_parts is reference to cache
					existing_support_parts.push_back(new InferenceSupportData(ext.uuid_generators_inference_support_data.generateUUID(),
																				support_subpart->support_subpart, 
																				support_subpart_leaf_to_root_path));
				}
			}
		} else {
			// this is first layer part - add directly to as support part
			existing_support_parts.push_back(new InferenceSupportData(ext.uuid_generators_inference_support_data.generateUUID(),
																		inferred_part, 
																		std::vector<UUIDType>()));
		}
	}

	return existing_support_parts;
}

std::set<InferredPart*> InferenceSupportAccess::getOnlyDifferentInitialLayerSupportParts(const IAttachableClass& part) {

	// TODO: implement more efficiently !!! by caching all points instead of using getAllInitialLayerSupportParts !!!

	const std::vector<InferenceSupportData*>& all_supporing_parts = getAllInitialLayerSupportParts(part);

	// compare parts by their UUIDs
	struct InferredPartComparer {
		bool operator()(const InferredPart* a, const InferredPart* b) { return a->getUUID() < b->getUUID(); }
	};

	std::set<InferredPart*> only_different_supporting_parts;

	// just extract all different layer 1 parts TODO: is comparison by pointer OK ? 
	for (auto iter = all_supporing_parts.begin(); iter != all_supporing_parts.end(); ++iter) {
		only_different_supporting_parts.insert(&(*iter)->support_subpart);
	}

	return only_different_supporting_parts;
}

// we need ocv.h only for std::less<cv::Point2i> !!
#include "utils/ocv.h"

float InferenceSupportExtension::calculatePartsIntersection(const std::set<InferredPart*> support_a, const std::set<InferredPart*> support_b) {

	// TODO: current implementation is not the same as old one !!!!!
	//		 we should also use vocabulary region data points and calculate support based on those points !!!
	//       -> replace getLocation() with getLocation() + getCorrespondingVocabularyPart().getRegions()  !!!

	std::set<cv::Point2i> support_a_points, support_b_points;

	for (auto iter = support_a.cbegin(); iter != support_a.cend(); ++iter) {
		support_a_points.insert((*iter)->getLocation());
	}

	for (auto iter = support_b.cbegin(); iter != support_b.cend(); ++iter) {
		support_b_points.insert((*iter)->getLocation());
	}

	auto inter_a = support_a_points.cbegin();
	auto inter_b = support_b_points.cbegin();

	float result = 0;

	std::less<cv::Point2i> cvPointLessOperator;
	while (inter_a != support_a_points.end() && inter_b != support_b_points.end()) {
		if (cvPointLessOperator(*inter_a, *inter_b)) ++inter_a; 
		else if (cvPointLessOperator(*inter_b, *inter_a)) ++inter_b;
		else {
			++result;
			++inter_a; ++inter_b;
		}
	}

	return result / std::min<int>(support_a_points.size(),support_b_points.size());
}

////////////////////////////////////////////////////////////////////////////////
//// Implementation for support of VocabularyPart
