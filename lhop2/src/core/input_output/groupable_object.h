#pragma once
#ifndef _CORE_GRPUPABLE_OBJECT__
#define _CORE_GRPUPABLE_OBJECT__

/// @addtogroup core
/// @{
/// @addtogroup input_output
/// @{


/**
 * Abstract representation of any object that can be grouped in one or more groups. 
 * Object's group is identified with group_id and object itself is identified by group_member_id within a group 
 * (variable are defined in struct GroupableObject::Member)
 *
 * This class can be used by pre/post-processing when one input object gets split into multiple input objects
 * thus producing multiple output object from one input. Use GroupableObject to identify splitted objects and to
 * find correct groups of objects for appropriate merging where needed.
 */ 
class GroupableObject {
public:
	struct Member {
		int group_id;
		int group_member_id;
	};

	struct Type {
		std::string name;

		bool operator< (const Type& obj) const {
			return name < obj.name;
		}
	};
protected:
	std::multimap<Type, Member> group;
public:
	std::multimap<Type, Member> getGroupMap() const {
		return group;
	}
	void addGroup(Type type, Member member) {
		group.insert(std::pair<Type, Member>(type,member));
	}
	void removeGroupId(std::string type, int group_id) {
		for (auto iter = group.begin(); iter != group.end(); ++iter) {
			if (type.compare(iter->first.name) == 0 && iter->second.group_id == group_id) {
				group.erase(iter);
				break;
			}
		}
	}
};

/// @}
/// @}

#endif /* _CORE_GRPUPABLE_OBJECT__ */
