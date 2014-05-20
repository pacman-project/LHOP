/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// subparts_extension

#include "subparts_extension.h"

#include "core/structures/serializers.h"

////////////////////////////////////////////////////////////////////////////////
//// Inference subparts extension

void InferenceSubpartsExtension::moveTo(IExtension& dst_extension) {
	InferenceSubpartsExtension& dst_subpart_extension = dynamic_cast<InferenceSubpartsExtension&>(dst_extension);

	// copy subparts_links to destination
	dst_subpart_extension.subparts_links.insert(subparts_links.begin(), subparts_links.end());

	// and just clear current array (do not delete any InferenceSubpartData pointers as they are still alive in the dst_extension) 
	subparts_links.clear();
}
void InferenceSubpartsModifier::deleteAllReferences(const std::vector<IAttachableClass*>& reference_to) {
	// delete subpart 		

	// type should be InferredPart otherwise we cannot handle it
	for (auto iter = reference_to.begin(); iter != reference_to.end(); ++iter) {
		InferredPart& part = *(InferredPart*)*iter;
	
		// find this part in the internal map 
		auto found_part_iter = ext.subparts_links.find(part.getUUID());

		if (found_part_iter != ext.subparts_links.end()) {
		
			/// first delete its subparts
			deleteSubparts(part, found_part_iter->second);

			// then remove part from internal list
			ext.subparts_links.erase(found_part_iter);
		}
	}
}

void InferenceSubpartsModifier::insertSubparts(const IAttachableClass& part, const std::vector<InferenceSubpartData*>& subparts_for_insertion) {
	std::vector<InferenceSubpartData*>& existing_subparts_of_part = ext.subparts_links[part.getUUID()];
	existing_subparts_of_part.insert(existing_subparts_of_part.end(), subparts_for_insertion.begin(), subparts_for_insertion.end());
}
void InferenceSubpartsModifier::deleteSubparts(const IAttachableClass& part, const std::vector<InferenceSubpartData*>& subparts_for_deletion) {

	//////////////////////////////////////////////////
	// Delete references pointing to this subparts
	std::vector<IAttachableClass*> abstract_subparts_for_deletion(subparts_for_deletion.begin(), subparts_for_deletion.end());

	// delete any references to the subparts (if any class was using InferenceSubpartData as reference)
	ext.holder.deleteAttachedReferences<InferenceSubpartData>(abstract_subparts_for_deletion);

	//////////////////////////////////////////////////
	// Delete subparts		
	std::vector<InferenceSubpartData*>&  existing_subparts_of_part = ext.subparts_links[part.getUUID()];

	// remove each subpart of subparts_for_insertion from existing_subparts_of_part list
	std::remove_if(existing_subparts_of_part.begin(), existing_subparts_of_part.end(), 
					
					// lambda functions: returns true if existing_subpart is contained in the subparts_for_deletion list
					[&subparts_for_deletion](InferenceSubpartData* existing_subpart) {
						
						return std::any_of(subparts_for_deletion.begin(), subparts_for_deletion.end(),
											
											// lambda function: equal operation
											[&existing_subpart](InferenceSubpartData* deletion_subpart){
												return existing_subpart == deletion_subpart;
											});
						});

	// TOOD: should we delete memory for all subparts_for_deletion ??
}

AbstractSerializer::IFactory* InferenceSubpartsExtension::getSerializer() const {
	return new DefaultInferenceSubpartsExtSerializer::Factory();
}

UUIDGenerator& InferenceSubpartsModifier::getInferrenceSubpartUUIDGenerator() {
	return ext.uuid_generators_inference_subpart_data;
}

////////////////////////////////////////////////////////////////////////////////
//// Vocabulary subparts extension


void VocabularySubpartsExtension::moveTo(IExtension& dst_extension) {
	VocabularySubpartsExtension& dst_subpart_extension = dynamic_cast<VocabularySubpartsExtension&>(dst_extension);

	// copy subparts_links to destination
	dst_subpart_extension.subparts_links.insert(subparts_links.begin(), subparts_links.end());

	// and just clear current array (do not delete any InferenceSubpartData pointers as they are still alive in the dst_extension) 
	subparts_links.clear();
}


void VocabularySubpartsModifier::deleteAllReferences(const std::vector<IAttachableClass*>& reference_to) {
	// delete subpart 		

	// type should be VocabularyPart otherwise we cannot handle it
	for (auto iter = reference_to.begin(); iter != reference_to.end(); ++iter) {
		VocabularyPart& part = *(VocabularyPart*)*iter;
	
		// find this part in the internal map 
		auto found_part_iter = ext.subparts_links.find(part.getUUID());

		if (found_part_iter != ext.subparts_links.end()) {
		
			/// first delete its subparts
			deleteSubparts(part, found_part_iter->second);

			// then remove part from internal list
			ext.subparts_links.erase(found_part_iter);
		}
	}
}

void VocabularySubpartsModifier::insertSubparts(const IAttachableClass& part, const std::vector<VocabularySubpartData*>& subparts_for_insertion) {
	std::vector<VocabularySubpartData*>& existing_subparts_of_part = ext.subparts_links[part.getUUID()];
	existing_subparts_of_part.insert(existing_subparts_of_part.end(), subparts_for_insertion.begin(), subparts_for_insertion.end());
}
void VocabularySubpartsModifier::deleteSubparts(const IAttachableClass& part, const std::vector<VocabularySubpartData*>& subparts_for_deletion) {

	//////////////////////////////////////////////////
	// Delete references pointing to this subparts
	std::vector<IAttachableClass*> abstract_subparts_for_deletion(subparts_for_deletion.begin(), subparts_for_deletion.end());

	// delete any references to the subparts (if any class was using InferenceSubpartData as reference)
	ext.holder.deleteAttachedReferences<VocabularySubpartData>(abstract_subparts_for_deletion);

	//////////////////////////////////////////////////
	// Delete subparts		
	std::vector<VocabularySubpartData*>&  existing_subparts_of_part = ext.subparts_links[part.getUUID()];

	// remove each subpart of subparts_for_insertion from existing_subparts_of_part list
	std::remove_if(existing_subparts_of_part.begin(), existing_subparts_of_part.end(), 
					
					// lambda functions: returns true if existing_subpart is contained in the subparts_for_deletion list
					[&subparts_for_deletion](VocabularySubpartData* existing_subpart) {
						
						return std::any_of(subparts_for_deletion.begin(), subparts_for_deletion.end(),
											
											// lambda function: equal operation
											[&existing_subpart](VocabularySubpartData* deletion_subpart){
												return existing_subpart == deletion_subpart;
											});
						});

}


UUIDGenerator& VocabularySubpartsModifier::getVocabularySubpartUUIDGenerator() {
	return ext.uuid_generators_vocabulary_subpart_data;
}