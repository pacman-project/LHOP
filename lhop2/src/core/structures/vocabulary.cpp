/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// vocabulary


#include "vocabulary.h"

VocabularyPart& VocabularyLayer::getPart(UUIDType uuid) {
	return *parts[uuid];
}

const VocabularyPart& VocabularyLayer::getPart(UUIDType uuid) const {
	auto iter = parts.find(uuid);
	if (iter != parts.end())
		return *iter->second;
	else
		throw new_libhop_exception("Access to vocabulary part with unknown id: uuid/type_id does not exist at this layer");
}

VocabularyPart& VocabularyLayer::getPartByTypeId(int type_id) {
	for (auto iter = parts.begin(); iter != parts.end(); ++iter) {
		if (iter->second->getTypeId() == type_id)
			return *iter->second;
	}
	throw new_libhop_exception("Part with specific type ID not found on this layer.");
}

void VocabularyLayer::insertNewPart(VocabularyPart* part) {
	parts[part->getUUID()] = part;
}

////////////////////////////////////////////////////////////////////////
//// VocabularyTreeFromPartLib

#include "core/legacy/constants.h"
#include "core/structures/indexing/indexing_extension.h"
#include "core/structures/subparts/subparts_extension.h"
#include "core/structures/shape_context/appearance_extension.h"



void VocabularyTreeFromPartLib::initializeFromPartLib(part_lib* lib) {
	
	std::unordered_map<node*, VocabularyPart*> new_parts_mapping;

	UUIDGenerator& vocabulary_part_uuid_gen = getVocabularyPartUUIDGenerator();

	// first insert raw parts for all layers
	// and associated original node* with new VocabularyPart*
	for (int layer = 0; layer < lib->layer_count; ++layer) {

		VocabularyLayer* vocabulary_layer = new VocabularyLayer(*this, lib->contractions[layer], layer);

		std::vector<node*>& layer_parts = lib->parts[layer];

		for (auto iter = layer_parts.begin(); iter != layer_parts.end(); ++iter) {
			node* p = *iter;

			// TODO: how do we handle other types
			// currently we handle only part_data
			part_data* pd = dynamic_cast<part_data*>(p->data);
			
			UUIDType uuid = vocabulary_part_uuid_gen.generateUUID();

			ResponsesArray thresholds;
			thresholds.set(ResponseType::R_RESPONSE_1, pd->td.get_thresh(R_THRESH));
			thresholds.set(ResponseType::G_RESPONSE_1, pd->td.get_thresh(G_THRESH));
			thresholds.set(ResponseType::S_RESPONSE_1, pd->td.get_thresh(S_THRESH));
			thresholds.set(ResponseType::RR_RESPONSE_1, pd->td.get_thresh(RR_THRESH));

			VocabularyPart* new_vocabulary_part = new VocabularyPart(uuid, layer, pd->type, pd->mask, pd->region, thresholds, cv::Point2i(pd->cmx, pd->cmy));
			vocabulary_layer->insertNewPart(new_vocabulary_part);

			// map existing node of part with new vocabulary part (so that we can later easily add links/edges between parts)
			new_parts_mapping[p] = new_vocabulary_part;
		}

		this->insertNewLayer(vocabulary_layer);
	}
	
	this->insertExtension(new VocabularyIndexingExtension(*this));
	this->insertExtension(new VocabularySubpartsExtension(*this));
	this->insertExtension(new VocabularyAppearanceExtension(*this));

	std::unordered_map<edge_data*, VocabularySubpartData*> new_subparts_mapping;
	
	VocabularyIndexingModifier vocabulary_indexing_modifier = this->getAccess<VocabularyIndexingModifier>();
	VocabularySubpartsModifier vocabulary_subparts_modifier = this->getAccess<VocabularySubpartsModifier>();

	UUIDGenerator& vocabulary_subpart_uuid_gen = vocabulary_subparts_modifier.getVocabularySubpartUUIDGenerator();

	// then add all connections for all layers
	for (int layer = 0; layer < lib->layer_count; ++layer) {

		VocabularyLayer& vocabulary_layer = this->getLayer(layer);

		std::vector<node*>& layer_parts = lib->parts[layer];

		for (auto iter = layer_parts.begin(); iter != layer_parts.end(); ++iter) {
			node* p = *iter;

			part_data* pd = dynamic_cast<part_data*>(p->data);
			
			// find corresponding new vocabulary part
			VocabularyPart* vocabulary_part = new_parts_mapping[p];

			std::vector<VocabularyPart*> indexed_central_parts;
			std::vector<VocabularySubpartData*> subparts;

			// insert different types of connections
			for (auto neighbor_iter = p->neighbors.begin(); neighbor_iter != p->neighbors.end(); ++neighbor_iter) {
				int neighboor_type = neighbor_iter->first;
				
				node* neighboor_node = neighbor_iter->second.first;
				edge_data* neighboor_data = neighbor_iter->second.second;

				// find corresponding neighbor new vocabulary part
				VocabularyPart* neighboor_vocabulary_part = new_parts_mapping[neighboor_node];

				switch (neighboor_type) {
					case EdgeConnection::TO_NEIGHBOOR: {
						break;
					}
					case EdgeConnection::TO_LYR_PREV: {
						// TO_LYR_PREV == TO_LYR_SOURCE + TO_LYR_CENTER 
						// use VocabularySubpartsModifier						
						break;
					}					
					case EdgeConnection::TO_LYR_SOURCE: {
						// Edges to all subparts (one layer down) other then center (in library; class part_lib). connections to all parts in lower layer (one layer down) for each of the subparts other then center	
						part_data_2* part_data = (part_data_2*)neighboor_data;

						// TODO: save part_data->geo to geometry extension
						// TODO: save part_data->app to appearance extension
						
						UUIDType uuid = vocabulary_subpart_uuid_gen.generateUUID();

						VocabularySubpartData* subpart_data = new VocabularySubpartData(uuid, *neighboor_vocabulary_part, part_data->index, part_data->distr.toCvMat(), cv::Point2i(part_data->x, part_data->y), part_data->gdistr);
						subparts.push_back(subpart_data);

						new_subparts_mapping[neighboor_data] = subpart_data;
						break;
					}					
					case EdgeConnection::TO_LYR_CENTER: {
						// Edge to center (in library; class part_lib). inverse of lyrCenterBack (never used - maybe ?? :D)
						// use VocabularySubpartsModifier
						part_data_2a* part_data = (part_data_2a*)neighboor_data;

						// TODO: save part_data->geo to geometry extension
						// TODO: save part_data->app to appearance extension
						
						UUIDType uuid = vocabulary_subpart_uuid_gen.generateUUID();

						VocabularySubpartData* subpart_data = new VocabularySubpartData(uuid, *neighboor_vocabulary_part, 0, cv::Mat(), cv::Point2i(), normal_distribution1());
						subparts.push_back(subpart_data);

						new_subparts_mapping[neighboor_data] = subpart_data;
						break;
					}
					case EdgeConnection::TO_LYR_CENTER_BACK: {
						// Edges to all nodes (i.e. parts) in the next layer (one level up) which have same  part for its center (in library; class part_lib). Connectios to all parts in upper layer (one level one) which have this part for its center
						// use VocabularyIndexingModifier						
						indexed_central_parts.push_back(neighboor_vocabulary_part);
						break;
					}

					/////////////////////////////////////////////////////////////////////////////////////////					
					case EdgeConnection::TO_LYR_FORBIDDEN: {
						// Edges to forbidden subparts -- parts which are not allowed to be present at specific positions -- (one layer down, in library; class part_lib). 
						// TODO: write extension for forbidden subparts
						break;
					}
					case EdgeConnection::TO_LYR_SIMILAR: {
						// TODO: write extension for similarity between parts
						break;
					}
					
					case EdgeConnection::TO_LYR_SIMROOT: {
						// TODO: write extension for similarity between parts
						break;
					}
					/////////////////////////////////////////////////////////////////////////////////////////
					default: {
					}
				}
			}

			vocabulary_indexing_modifier.insertIndexedCentralParts(*vocabulary_part, indexed_central_parts);
			vocabulary_subparts_modifier.insertSubparts(*vocabulary_part, subparts);
		}
	}


}