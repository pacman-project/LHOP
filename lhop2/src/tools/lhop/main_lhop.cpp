// main_lhop.cpp : Defines the main entry point for lhop
//

#include "utils/utils.h"
#include "utils/toolset.h"
#include "utils/class_register.h"

#include "tools/main_toolset.h"

#include <stdio.h>
#include <string.h>

using namespace std;

/**
 * Registration method that MUST never be called but MUST BE PRESENT !!!
 * (ClassRegister::registerFactory uses some template voodoo to register instance in 
 */
void register_main_toolset_classes() {
	/************************************************************************* 
	 *           !!!! DO NOT MODIFY OR REMOVE THIS FUNCTION   !!!!!
	 ************************************************************************/

	ClassRegister::registerFactory<IToolset::IFactory, FeatureExtractionToolset::Factory>();
	ClassRegister::registerFactory<IToolset::IFactory, LearningToolset::Factory>();
	ClassRegister::registerFactory<IToolset::IFactory, InferenceToolset::Factory>();

	ClassRegister::registerFactory<IToolset::IFactory, LegacyFeatureExtractionToolset::Factory>();
	ClassRegister::registerFactory<IToolset::IFactory, LegacyInferenceToolset::Factory>();
}

void printToolsetList() {

}

// Run main toolset defined by the first program argument (second in argv.
// Return 1 on success or 0 failure
bool runToolset(int argc, char* argv[]) {

	const ClassRegister& plugables_manager = ClassRegister::get();

	const string& toolset_command = argv[1];

	// retreive implementation factory for specific user supplied toolset command
	IToolset::IFactory* toolset_factory = plugables_manager.retrieveFactory<IToolset::IFactory>(toolset_command, false);

	if (toolset_factory == nullptr) {
		// display error
		cout << "Error: invalid command '" << toolset_command << "'" << endl;
		// and return failure
		return false;
	}

	// then create new toolset using standard factory constructors
	IToolset* toolset = toolset_factory->newInstance(plugables_manager);
	
	bool display_help = false;

	// check for help argument
	for (int i = 2; i < argc; i++) {
		if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
			cout << "found help request" << endl;
			display_help = true;
			break;
		}
	}
	// also check if arguments are valid (let implementing toolset decide if arguments are valid)
	if (display_help == false)
		display_help = toolset->areArgumentsValid(argc, argv) ? false : true;

	if (display_help) {
		// display help (usage + long description)
		cout << "Usage of '" << toolset_command << "' command:" << endl
			<< "\t\t" << toolset->getUsageDescription() << endl << endl
			<< "Command description:" << endl << endl
			<< "\t\t" << toolset->getLongDescription() << endl << endl;
	} else {	
		// run toolset
		toolset->main(argc, argv);
	}

	// clean memory
	if (toolset != nullptr)	delete toolset;
	if (toolset_factory != nullptr) delete toolset_factory;

	// return success
	return true;
}

int main(int argc, char* argv[])
{
	// display compile timestamp
	cout << "lhop (" __DATE__ " / " __TIME__ ")" << endl << endl;
	try {
		// initiate random (maybe should add into some initi script ??)
		random_seed((int)time(nullptr));

		// main swicth	
		switch (argc) {
			// print toolset list if no toolset command supplied 
			case 1: printToolsetList(); break;
			default: 
				// run main function of requested toolset (or print help if invalid of arguments supplied)
				if (runToolset(argc, argv) == 0)
					printToolsetList();;
		}
	
    } catch (const libhop_exception& e) {
		cout << e.what() << endl;
    } catch (const exception& e) {
		cout << "General exception '" << typeid(e).name() << "' with message: '" << e.what() << "'" << endl;
	} 	
	return 0;
}

