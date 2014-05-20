/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// indexing_extension

#include "indexing_extension.h"

#include "core/structures/serializers.h"

void VocabularyIndexingExtension::moveTo(IExtension& dst_extension){
	VocabularyIndexingExtension& dst_indexing_extension = dynamic_cast<VocabularyIndexingExtension&>(dst_extension);

	// copy indexed_center_links to destination for each part 
	for (auto iter = indexed_center_links.begin(); iter != indexed_center_links.begin(); ++iter) {
		auto existing_links = dst_indexing_extension.indexed_center_links[iter->first];
		existing_links.insert(existing_links.begin(), iter->second.begin(), iter->second.end());
	}

	// and just clear current array (do not delete any InferenceSubpartData pointers as they are still alive in the dst_extension) 
	indexed_center_links.clear();
}

AbstractSerializer::IFactory* VocabularyIndexingExtension::getSerializer() const {
	// TODO: how do we handle multiple serializers ??
	return new DefaultVocabularyIndexingSerializer::Factory();
}


void VocabularyIndexingModifier::deleteAllReferences(const std::vector<IAttachableClass*>& reference_to) {
	// delete subpart 		

	// type should be InferredPart otherwise we cannot handle it
	for (auto iter = reference_to.begin(); iter != reference_to.end(); ++iter) {
		IAttachableClass& part = (IAttachableClass&)**iter;
	
		// find this part in the internal map 
		auto found_part_iter = ext.indexed_center_links.find(part.getUUID());

		if (found_part_iter != ext.indexed_center_links.end()) {
		
			/// first delete its subparts
			insertIndexedCentralParts(part, found_part_iter->second);

			// then remove part from internal list
			ext.indexed_center_links.erase(found_part_iter);
		}
	}
}

void VocabularyIndexingModifier::insertIndexedCentralParts(const IAttachableClass& part, const std::vector<VocabularyPart*>& indexes_for_insertion){
	std::vector<VocabularyPart*>& existing_subparts_of_part = ext.indexed_center_links[part.getUUID()];
	existing_subparts_of_part.insert(existing_subparts_of_part.end(), indexes_for_insertion.begin(), indexes_for_insertion.end());
}

void VocabularyIndexingModifier::deleteIndexedCentralParts(const IAttachableClass& part, const std::vector<VocabularyPart*>& indexes_for_deletion){

	//////////////////////////////////////////////////
	// Delete references pointing to this subparts
	std::vector<IAttachableClass*> abstract_indexes_for_deletion(indexes_for_deletion.begin(), indexes_for_deletion.end());

	// delete any references to the subparts (if any class was using InferenceSubpartData as reference)
	ext.holder.deleteAttachedReferences<VocabularyPart>(abstract_indexes_for_deletion);

	//////////////////////////////////////////////////
	// Delete subparts		
	std::vector<VocabularyPart*>&  existing_indexes_of_part = ext.indexed_center_links[part.getUUID()];

	// remove each subpart of subparts_for_insertion from existing_subparts_of_part list
	std::remove_if(existing_indexes_of_part.begin(), existing_indexes_of_part.end(), 
					
					// lambda functions: returns true if existing_subpart is contained in the subparts_for_deletion list
					[&indexes_for_deletion](VocabularyPart* existing_index) {
						
						return std::any_of(indexes_for_deletion.begin(), indexes_for_deletion.end(),
											
											// lambda function: equal operation
											[&existing_index](VocabularyPart* deletion_index){
												return existing_index == deletion_index;
											});
						});

}
