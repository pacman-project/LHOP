
#include "output_postprocessing.h"

/**
 * Class enables grouping of the GroupableObject to allow for merging of the results.
 * Execute in the following procedure:
 *
 * ObjectGrouping group(initial_objects);
 * while (group.hasNext(group_type_str) {
 *   std::map<int,GroupableObject*> group_members = group.getNext();
 *   /// merge group members into results: group_members -> results
 *   group.postMergeResults(results)
 * }
 * group.getFinalResults();
 */ 
class ObjectGrouping {
private:
	std::string group_type;

	std::vector<GroupableObject*> remaining_output_list;

	std::vector<GroupableObject*> final_results;

	// partial resuls
	std::map<int,GroupableObject*> current_group_members;
	int highest_group_id;
public:
	ObjectGrouping(std::vector<AbstractOutputObject*> initial_list, const std::string& group_type) : group_type(group_type) {
		for (auto iter = initial_list.begin(); iter != initial_list.end(); ++iter) {
			remaining_output_list.push_back(castToGroupableObject(*iter));
		}
	}

	/**
	 * Casting to GroupableObject* needed for any grouping.
	 */
	static GroupableObject* castToGroupableObject(AbstractOutputObject* input_object){
		GroupableObject* output_groupable = dynamic_cast<GroupableObject*>(input_object);

		if (output_groupable == nullptr)
			throw new_libhop_exception("GroupablePostprocessing can only preprocess GroupableObject* objects");

		return output_groupable;
	}

	bool hasNext() {
		return remaining_output_list.size() > 0 ? true : false;
	}

	std::map<int,GroupableObject*> getNext() {
		// find all possible group ids
		highest_group_id = INT_MIN;

		std::map<int,std::map<int,GroupableObject*>> group_members;

		// find all OutputObjects that belong to the group_type
		for (auto iter = remaining_output_list.begin(); iter != remaining_output_list.end(); ++iter) {
			GroupableObject* output_layer = *iter;

			std::multimap<GroupableObject::Type, GroupableObject::Member> groups = output_layer->getGroupMap();

			// find all groups ids of group_type
			bool is_group_type_member = false;
			for (auto group_iter = groups.begin(); group_iter != groups.end(); ++group_iter) {
				GroupableObject::Type input_group_type = group_iter->first;
				GroupableObject::Member input_group_member = group_iter->second;

				if (group_type.compare(input_group_type.name) == 0) {
					group_members[input_group_member.group_id].insert(std::pair<int,GroupableObject*>(input_group_member.group_member_id,output_layer));

					highest_group_id = max<int>(highest_group_id, input_group_member.group_id);
					is_group_type_member = true;
				}
			}
			if (is_group_type_member == false)
				final_results.push_back((AbstractOutputObject*)output_layer); // pass objects that are not of interest into results
		}
		
		// clear remaining_output_list as it will be input for next loop
		remaining_output_list.clear();

		// 
		for (auto group_iter = group_members.begin(); group_iter != group_members.end(); ++group_iter) {
			std::map<int,GroupableObject*>& loop_group_members = group_iter->second;

			// return only group with highest id and pass all the others to remaining list
			if (group_iter->first != highest_group_id) {
				// copy all others to the remaining output list 
				for (auto iter = loop_group_members.begin(); iter != loop_group_members.end(); ++iter) {
					remaining_output_list.push_back(iter->second);
				}
			}
		}

		// save group members with the highest ID
		current_group_members = group_members[highest_group_id];

		// return group members next for merging (i.e. the ones with highest group ID)
		return current_group_members;
	}

	void postMergeResults(std::vector<GroupableObject*> merged_result) {		
		for (auto iter = merged_result.begin(); iter != merged_result.end(); ++iter) {
			// remove grouping from the resuls so they will not be processed any more
			(*iter)->removeGroupId(group_type, highest_group_id);
			
			// put into the list for next iteration
			remaining_output_list.push_back(*iter);
		}
	}

	std::vector<GroupableObject*> getFinalResult() {
		return final_results;
	}
};

std::vector<AbstractOutputObject*> CombineLayerPartsPostprocessing::doPostprocessing(const std::vector<AbstractOutputObject*>& output_list) const {

	// use ObjectGrouping to group output_list and merge only groups returned by the ObjectGrouping
	ObjectGrouping group(output_list, group_type);

	while (group.hasNext()) {
		// get next group ready for merging
		std::map<int,GroupableObject*> group_members = group.getNext();

		// perform merging
		std::shared_ptr<InferenceTree> merged_results(nullptr);
		std::vector<InferenceTree*> tmp_result(1);
				
		for (auto iter = group_members.begin(); iter != group_members.end(); ++iter) {
			// do merging
			LayerOutputObject* output_layer = castToLayerObject((AbstractOutputObject*)iter->second);

			std::shared_ptr<InferenceTree> layer_obj = output_layer->getLayerObject();

			if (merged_results == nullptr) {
				merged_results = layer_obj;
			} else  {				
				tmp_result[0] = layer_obj.get();
				
				merged_results->addExistingTree(tmp_result);
				//merged_results->add_results(tmp_result, 0);
			}
		}	
		
		// create new result and push it to the group in case the same object also will have to be processed twice
		std::vector<GroupableObject*> group_results;
		group_results.push_back(new LayerOutputObject(merged_results, output_list.front()->getGroupMap()));

		group.postMergeResults(group_results);
	}

	// get final results and copy them into the output list
	std::vector<GroupableObject*> grouped_results = group.getFinalResult();;

	std::vector<AbstractOutputObject*> results;
	for (auto iter = grouped_results.begin(); iter != grouped_results.end(); ++iter) {
		results.push_back((AbstractOutputObject*) *iter);
	}
	
	return results;
}

std::vector<AbstractOutputObject*> MergeScalesPostprocessing::doPostprocessing(const std::vector<AbstractOutputObject*>& output_list) const {

	// use ObjectGrouping to group output_list and merge only groups returned by the ObjectGrouping
	ObjectGrouping group(output_list, group_type);

	while (group.hasNext()) {
		// get next group ready for merging
		std::map<int,GroupableObject*> group_members = group.getNext();

		// perform merging
		std::shared_ptr<InferenceTree> merged_results(nullptr);
		vector<layer1_result*> tmp_result(1);
				
		for (auto iter = group_members.begin(); iter != group_members.end(); ++iter) {
			// do merging
			LayerOutputObject* output_layer = castToLayerObject((AbstractOutputObject*)iter->second);

			std::shared_ptr<InferenceTree> layer_obj = output_layer->getLayerObject();

			if (merged_results == nullptr) {
				merged_results = layer_obj;
			} else  {
				merged_results->mergeFrom(*layer_obj);
				//merged_results->merge(layer_obj.get(), layer_obj->border);
			}
		}	
		
		// create new result and push it to the group in case the same object also will have to be processed twice
		std::vector<GroupableObject*> group_results;
		group_results.push_back(new LayerOutputObject(merged_results, output_list.front()->getGroupMap()));

		group.postMergeResults(group_results);
	}

	// get final results and copy them into the output list
	std::vector<GroupableObject*> grouped_results = group.getFinalResult();;

	std::vector<AbstractOutputObject*> results;
	for (auto iter = grouped_results.begin(); iter != grouped_results.end(); ++iter) {
		results.push_back((AbstractOutputObject*) *iter);
	}
	
	return results;
}