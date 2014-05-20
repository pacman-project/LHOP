// main_toolset.h

#pragma once
#ifndef _TOOLS_LHOP_MAIN_TOOLSET_H_
#define _TOOLS_LHOP_MAIN_TOOLSET_H_

#include <string>
#include <list>

#include "utils/toolset.h"
#include "utils/class_register.h"
#include "utils/structures.h"

using namespace std;

/// @addtogroup toolset Toolset
/// @{


/**
 * Toolset used to extract first layer features.
 *  - call tag (assigned functionality name): "feature_extraction" 
 */ 
class FeatureExtractionToolset : public IToolset {
	const ClassRegister& pluggable_manager;
public:
	FeatureExtractionToolset(const ClassRegister& pluggable_manager) : pluggable_manager(pluggable_manager) {
	}

	virtual bool areArgumentsValid(int argc, char* argv[]);
	virtual void main(int argc, char* argv[]);	

	virtual string getShortDescription();
	virtual string getLongDescription();
	virtual string getUsageDescription();

	class Factory : public IFactory {

		virtual IToolset* newInstance(const ClassRegister& pluggables_manager) const { 
			return new FeatureExtractionToolset(pluggables_manager);
		};
		virtual string assignedRegistrationName() const { return "feature_extraction"; };

	};
};

class InferenceToolset : public IToolset {
	const ClassRegister& pluggable_manager;
public:
	InferenceToolset(const ClassRegister& pluggable_manager) : pluggable_manager(pluggable_manager) {
	}


	virtual bool areArgumentsValid(int argc, char* argv[]);
	virtual void main(int argc, char* argv[]);

	virtual string getShortDescription();
	virtual string getLongDescription();
	virtual string getUsageDescription();

	class Factory : public IFactory {

		virtual IToolset* newInstance(const ClassRegister& pluggables_manager) const { 
			return new InferenceToolset(pluggables_manager);
		};
		virtual string assignedRegistrationName() const { return "inference"; }
	};
};

class LearningToolset : public IToolset {
	const ClassRegister& pluggable_manager;
public:
	LearningToolset(const ClassRegister& pluggable_manager) : pluggable_manager(pluggable_manager) {
	}


	virtual bool areArgumentsValid(int argc, char* argv[]);
	virtual void main(int argc, char* argv[]);	

	virtual string getShortDescription();
	virtual string getLongDescription();
	virtual string getUsageDescription();

	class Factory : public IFactory {

		virtual IToolset* newInstance(const ClassRegister& pluggables_manager) const { 
			return new LearningToolset(pluggables_manager);
		};
		virtual string assignedRegistrationName() const { return "learning"; }
	};
};

// legacy toolsets

class LegacyFeatureExtractionToolset : public IToolset {
	const ClassRegister& pluggable_manager;
public:
	LegacyFeatureExtractionToolset (const ClassRegister& pluggable_manager) : pluggable_manager(pluggable_manager) {
	}

	virtual bool areArgumentsValid(int argc, char* argv[]);
	virtual void main(int argc, char* argv[]);	

	virtual string getShortDescription();
	virtual string getLongDescription();
	virtual string getUsageDescription();

	class Factory : public IFactory {

		virtual IToolset* newInstance(const ClassRegister& pluggables_manager) const { 
			return new LegacyFeatureExtractionToolset(pluggables_manager);
		};
		virtual string assignedRegistrationName() const { return "legacy_feature_extraction"; };

	};
};

class LegacyInferenceToolset : public IToolset {
	const ClassRegister& pluggable_manager;
public:
	LegacyInferenceToolset(const ClassRegister& pluggable_manager) : pluggable_manager(pluggable_manager) {
	}


	virtual bool areArgumentsValid(int argc, char* argv[]);
	virtual void main(int argc, char* argv[]);

	virtual string getShortDescription();
	virtual string getLongDescription();
	virtual string getUsageDescription();

	class Factory : public IFactory {

		virtual IToolset* newInstance(const ClassRegister& pluggables_manager) const { 
			return new LegacyInferenceToolset(pluggables_manager);
		};
		virtual string assignedRegistrationName() const { return "legacy_inference"; }
	};
};

/// groundtruth read/write functions
///////////////////////////////////////////////////

/// fname must have an extension; it is automatically changed to .grundtruth
void read_groundtruth(list<irectangle2>& rectangles, const string& fname, 
					  const string& cat_name_only = "", const string& gt_extension = ".groundtruth");

void read_groundtruth(list<std::pair<irectangle2,string> >& rectangles, const string& fname,
					  const string& cat_name_only = "", const string& gt_extension = ".groundtruth",
                      bool display_warnings = true);

void save_groundtruth(const list<std::pair<irectangle2,string> >& rectangles, const string& outdir, const string& srcname, 
    const string& prefix, const string& suffix, bool changeext = true);

/// @}
#endif /* _TOOLS_LHOP_MAIN_TOOLSET_H_*/