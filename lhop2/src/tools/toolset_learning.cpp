// learning.cpp : Defines the entry point feature extraction functionality
//

#include "tools/main_toolset.h"

#include "core/legacy/layers.h"
#include "core/legacy/learning/optimization.h"

using namespace laylearning;

void optimization_load_input(lib_learning_controler &learning_controler, const char* pattern) {

	// get configuration
	const ConfigDictionary& cfg =  learning_controler.get_cfg();

	// read library
    string libname;
    cfg.getValue(libname, "library", true);

    part_lib* library =  part_lib::read(libname);

	if (library == nullptr)
		throw custom_libhop_exception(ConfigException, string("Library '" + libname + "' not found"));

	learning_controler.set_starting_library(library);

	// read input files
	string dir = cfg.getValueString("src_dir", "");
	list<string> files;
	layer1_result* res;
	int count = 0;
	int maxcount = cfg.getValueInt("file_limit", INT_MAX);
	int start_layer = cfg.getValueInt("start_layer", 1);
	int end_layer = cfg.getValueInt("end_layer", 2);

	end_dir(dir);
	if (!cfg.getValueBool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.getValueString("pattern", "")))
		list_directory(files, dir + pattern);
	cout << "from_file = " << (cfg.getValueBool("from_file", false) ? "true" : "false") << endl;
	cout << "pattern = " << pattern << endl;
	cout << "src_dir = " << dir << endl;
	cout << "start_layer = " << start_layer << endl;
	cout << "end_layer = " << end_layer << endl;

	vector<streamed_pointer> input_data;
	vector<map<int,void*>> extra_data;

	for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
		if (++count > maxcount) break;

		const string filename = dir + *file;
		char* filename_c = new char[filename.size() + 1];
		filename.copy(filename_c, filename.size()); filename_c[filename.size()] = 0;

		map<int,void*> data_item;
		data_item[lib_learning_controler::INPUTDATA_IMAGE] = streamed_pointer::from_file(filename, false); // create streamed_pointer to layer1_result without loading the file
		data_item[lib_learning_controler::INPUTDATA_FILENAME] = filename_c;

		learning_controler.add_learning_data(data_item);		
	}
}

void optimization_save_results(lib_learning_results* results, lib_learning_controler &learning_controler) {

	const ConfigDictionary& cfg = learning_controler.get_cfg();

	lib_optimization_results* opt_results = static_cast<lib_optimization_results*>(results);

	int final_layer = opt_results->final_layer_number;

	// Save results
    string lib_export_name = cfg.getValueString("lib_export_name", "lib");
    string res_export_name = cfg.getValueString("res_export_name", "res");
	string name = cfg.getValueString("out_library", lib_export_name + (final_layer) + ".plb");
    int i = 0;

    PRINT_INFO("");
    PRINT_INFO("SAVING");

	opt_results->lib->save(name);

    if (cfg.getValueBool("save_test_set", false)) {

		list<streamed_pointer>& sresult = opt_results->processed_layers;

        for (list<streamed_pointer>::iterator iter = sresult.begin(); iter != sresult.end(); ++iter) {
            layer1_result* res = (layer1_result*)iter->get();

            name = res_export_name + (i++) + string(".ly") + (final_layer );
            PRINT_INFO(name);
            res->delete_edges_complement(EdgeConnection::TO_PREV_LAYER);
            res->save(name);
            delete res;
        }
    }
}


void random_selection(list<string>& files, int n, const string& file)
{
    if (n < 0 || n >= (int)files.size()) return;

    vector<int> perm;

    random_permutation(perm, (int)files.size());
    
    vector<int> todrop(perm.begin() + n, perm.end());
    ofstream os(file.c_str());
    list<string>::iterator iter;
    int count = 0;

    sort(todrop.begin(), todrop.end());
    iter = files.begin();
    while (iter != files.end()) {
        if (!binary_search(todrop.begin(), todrop.end(), count++)) ++iter;
        else {
            if (!os.fail()) os << *iter << endl;
            iter = files.erase(iter);
        }
    }
}
/**
 * learn_objects must receive list of following input data:
 *  - (INPUTDATA_VALIDATION_POS_IMAGE, INPUTDATA_GROUNDTRUTH, INPUTDATA_FILENAME) positive validation file
 *  - (INPUTDATA_VALIDATION_NEG_IMAGE, INPUTDATA_FILENAME) negative validation file
 *  - (INPUTDATA_IMAGE, INPUTDATA_GROUNDTRUTH, INPUTDATA_FILENAME) main training images
 * All images are layer1_results in streamable_pointer* form
 */
void learn_objects_load_input(lib_learning_controler &learning_controler, const char* patt) {
	
	const ConfigDictionary& cfg =  learning_controler.get_cfg();

	string pattern(patt);

	string dir = cfg.getValueString("src_dir", "");
	string vdir = cfg.getValueString("validation_src_dir", "");
	string vfilepos = cfg.getValueString("positive_validation_files", "");
	string vfileneg = cfg.getValueString("negative_validation_files", "");
	string libimgfile = cfg.getValueString("library_image_file", "");
	int rndchoice = cfg.getValueInt("random_select", -1);
	string catname = cfg.getValueString("category_name", "");
	int maxcount = cfg.getValueInt("file_limit", INT_MAX);

	// read library
    string libname;
    cfg.getValue(libname, "library", true);

    part_lib* library = part_lib::read(libname);

	if (library == nullptr)
		throw custom_libhop_exception(ConfigException, string("Library '" + libname + "' not found"));

	learning_controler.set_starting_library(library);


	list<string> files;
	list<string> vfilespos;
	list<string> vfilesneg;
	layer1_result* res = nullptr;

	end_dir(dir);
	end_dir(vdir);
	if (!cfg.getValueBool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.getValueString("pattern", "")))
		list_directory(files, dir + pattern);
	if (rndchoice > 0) random_selection(files, rndchoice, cfg.getValueString("dropped_files", ""));
	if (!vfilepos.empty() || !vfileneg.empty()) {
		file_list_from_file(vfilespos, vfilepos, vdir, cfg.getValueString("positive_pattern", ""));
		file_list_from_file(vfilesneg, vfileneg, vdir, cfg.getValueString("negative_pattern", ""));		
	}
	cout << "from_file = " << (cfg.getValueBool("from_file", false) ? "true" : "false") << endl;
	cout << "pattern = " << pattern << endl;
	cout << "src_dir = " << dir << endl;

	// read validation files
	for (list<string>::iterator vfile = vfilespos.begin(); vfile != vfilespos.end(); ++vfile) {

		const string filename = vdir + *vfile;
		char* filename_c = new char[filename.size() + 1];
		filename.copy(filename_c, filename.size()); filename_c[filename.size()] = 0;

		list<irectangle2>* gtrs = new list<irectangle2>();
		read_groundtruth(*gtrs, filename, catname);
		
		map<int,void*> data_item;
		data_item[lib_learning_controler::INPUTDATA_VALIDATION_POS_IMAGE] = streamed_pointer::from_file(filename, false); // create streamed_pointer to layer1_result without loading the file
		data_item[lib_learning_controler::INPUTDATA_FILENAME] = filename_c;
		data_item[lib_learning_controler::INPUTDATA_GROUNDTRUTH] = gtrs;

		learning_controler.add_learning_data(data_item);
	}
	for (list<string>::iterator vfile = vfilesneg.begin(); vfile != vfilesneg.end(); ++vfile) {		

		const string filename = vdir + *vfile;
		char* filename_c = new char[filename.size() + 1];
		filename.copy(filename_c, filename.size()); filename_c[filename.size()] = 0;

		map<int,void*> data_item;
		data_item[lib_learning_controler::INPUTDATA_VALIDATION_NEG_IMAGE] = streamed_pointer::from_file(filename, false); // create streamed_pointer to layer1_result without loading the file
		data_item[lib_learning_controler::INPUTDATA_FILENAME] = filename_c;

		learning_controler.add_learning_data(data_item);
	}

	// read main files
	int count = 0;

	for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
		if (++count > maxcount) break;

		const string filename = dir + *file;
		char* filename_c = new char[filename.size() + 1];
		filename.copy(filename_c, filename.size()); filename_c[filename.size()] = 0;

		list<irectangle2>* gtrs = new list<irectangle2>();
		read_groundtruth(*gtrs, dir + *file, catname);

		map<int,void*> data_item;
		data_item[lib_learning_controler::INPUTDATA_IMAGE] = streamed_pointer::from_file(filename, false); // create streamed_pointer to layer1_result without loading the file
		data_item[lib_learning_controler::INPUTDATA_FILENAME] = filename_c; 
		data_item[lib_learning_controler::INPUTDATA_GROUNDTRUTH] = gtrs;

		learning_controler.add_learning_data(data_item);
	}
}

void learn_objects_save_results(lib_learning_results* results, lib_learning_controler &learning_controler) {
	
	const ConfigDictionary& cfg = learning_controler.get_cfg();
	
	string libname = cfg.getValueString("out_library", "olib.plb");
	
	results->lib->save(libname);

	int layer;
	cfg.getValue(layer, "layer", true); // in case obj_learning::init_cfg will change layer definition then this has to be changed too
	string libimgfile = cfg.getValueString("library_image_file", "");
	
	if(!libimgfile.empty()) {
		cout << "Saving lib image: layer " << layer << endl;
		for(int il=1; il<=layer+3; il++) {
			string str;
			str += (string)"-" + il + ".png";
			cout << "Saving lib image (layer " << il << ") to " << change_extension(libimgfile, str).c_str() << endl;
			results->lib->save_all(change_extension(libimgfile, str).c_str(), il, 0, -1,
				cfg.getValueBool("show_labels", false));
		}
	}

}


/**
 * learn_objects must recieve list of following input data:
 *  - (INPUTDATA_VALIDATION_POS_IMAGE, INPUTDATA_GROUNDTRUTH, INPUTDATA_FILENAME) validation files
 *  - (INPUTDATA_IMAGE, INPUTDATA_GROUNDTRUTH, INPUTDATA_FILENAME) main training images
 * All images are layer1_results in streamable_pointer* form
 */
void learn_objects2_load_input(lib_learning_controler &learning_controler, const char* patt) {
	
	string pattern(patt);

	const ConfigDictionary& cfg = learning_controler.get_cfg();

	string dir = cfg.getValueString("src_dir", "");
	string vdir = cfg.getValueString("validation_src_dir", "");
	string vfiles = cfg.getValueString("validation_files", "");
	
	string catname;
	list<string> files;

	cfg.getValue(catname, "category_name", true);

	end_dir(dir);
	if (!cfg.getValueBool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.getValueString("pattern", "")))
		list_directory(files, dir + pattern);
	for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
		
		const string filename = dir + *file;
		char* filename_c = new char[filename.size() + 1];
		filename.copy(filename_c, filename.size()); filename_c[filename.size()] = 0;
		
		list<irectangle2>* gtrs = new list<irectangle2>();
		read_groundtruth(*gtrs, vdir + *file, catname);

		map<int,void*> data_item;
		data_item[lib_learning_controler::INPUTDATA_IMAGE] = streamed_pointer::from_file(filename, false);
		data_item[lib_learning_controler::INPUTDATA_FILENAME] = filename_c;
		
		if (gtrs->size() <= 0) {
			// groundtruth file is empty - we need to calculate bbox based on min/max location of shape nodes (we need to load layer1_result)
			layer1_result* res = (layer1_result*)(((streamed_pointer*)data_item[lib_learning_controler::INPUTDATA_IMAGE])->get());

			if (res != nullptr)
				gtrs->push_back(get_nodes_bounding_rectangle(res->shape_nodes[0].begin(), res->shape_nodes[0].end()));

			delete res;
		}	

		data_item[lib_learning_controler::INPUTDATA_GROUNDTRUTH] = gtrs;

		learning_controler.add_learning_data(data_item);
	}

	files.clear();
	file_list_from_file(files, vfiles, vdir, cfg.getValueString("validation_pattern", ""));
	for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
		const string filename = vdir + *file;
		char* filename_c = new char[filename.size() + 1];
		filename.copy(filename_c, filename.size()); filename_c[filename.size()] = 0;

		list<irectangle2>* gtrs = new list<irectangle2>();
		read_groundtruth(*gtrs, vdir + *file, catname);

		map<int,void*> data_item;
		data_item[lib_learning_controler::INPUTDATA_VALIDATION_POS_IMAGE] = streamed_pointer::from_file(filename, false);
		data_item[lib_learning_controler::INPUTDATA_FILENAME] = filename_c;
		
		if (gtrs->size() <= 0) {
			// groundtruth file is empty - we need to calculate bbox based on min/max location of shape nodes (we need to load layer1_result)
			layer1_result* res = (layer1_result*)(((streamed_pointer*)data_item[lib_learning_controler::INPUTDATA_VALIDATION_POS_IMAGE])->get());

			if (res != nullptr)
				gtrs->push_back(irectangle2(2*res->x_size(0), 2*res->x_size(0), 3*res->x_size(0), 3*res->y_size(0)));

			delete res;
		}	

		data_item[lib_learning_controler::INPUTDATA_GROUNDTRUTH] = gtrs;

		learning_controler.add_learning_data(data_item);
	}
}
void learn_objects2_save_results(lib_learning_results* results, lib_learning_controler &learning_controler) {
	const ConfigDictionary& cfg = learning_controler.get_cfg();

	string libname = cfg.getValueString("out_library", "olib.plb");

	results->lib->save(libname);
}


string LearningToolset::getShortDescription() {
	return "";
}
string LearningToolset::getLongDescription() {
	return "";
}
string LearningToolset::getUsageDescription() {
	return "";
}

bool LearningToolset::areArgumentsValid(int argc, char* argv[]) {
	// TODO: do more intelligent verification
	return argc > 2  ? true : false;
}


#include "core/learning/abstract_learning.h"

#include "tools/file_input.h"
#include "tools/file_output.h"

#include "tools/main_toolset.h"

#include "utils/utils.h"
#include "utils/configuration.h"

// new learning implementation
void new_learning(int argc, char* argv[]) {
	// parse input arguments
	const char* cfgfile = argv[2];						// path to configuration file
	const char* pattern = argc > 3 ? argv[3] : "";		// optional pattern
	const char* params = argc > 4 ? argv[4] : "";		// optional string of additional configuration values

	DictionaryConfiguration learning_config;

	// load from string first but do not update namespace references (i.e. leave .from_namespace = ... values as is)
	learning_config.appendConfig(StringDictionaryConfiguration(params)); 
	// then also load from file (append will NOT override current values) 
	// and update namespace references at the same time (i.e. replace .from_namespace = ... with the values from the actual namespace)
	learning_config.appendConfig(FileDictionaryConfiguration(cfgfile));  

	// get inference.ly1 namespace
	IConfiguration* config_ptr = learning_config.getNamespace("learning");

	if (config_ptr == nullptr)
		throw new_libhop_exception("Missing learning namespace in configuration");

	IConfiguration& config = *config_ptr;

	std::string srcdir = fixPathEndSeparator(config.getString("src_dir"));
	std::string out_library = config.getString("out_library");

	std::list<std::string> files;

	if (!config.getBool("from_file") || !list_from_file(files, pattern, srcdir))
		list_directory(files, srcdir + "/" + string(pattern));

	std::string library_name = config.getString("library");

	// load vocabulary from file
	FileVocabularyInputObject vocabulary(library_name);

	// create multilayer learning that will be able to handle each individual layer 
	// (irrespectively of the type of learning - configuration will define which type will be used)
	std::string multiple_layer_learning_type = MultipleLayersBatchLearning::Factory().assignedRegistrationName(); // just find the registration name of the MultipleLayersBatchLearning

	// find factory for specific learning type
	AbstractLayerLearning::IFactory* factory = ClassRegister::get().retrieveFactory<AbstractLayerLearning::IFactory>(multiple_layer_learning_type);

	// create inference object (will load all the necessary settings)
	AbstractLayerLearning* layers_learning = factory->newInstance(config, vocabulary.getVocabularyObject());

	// create a list of all training files
	std::vector<AbstractLayerInputObject*> input_learning_objects;

	FileLayerInputObject::Factory file_input_object_factory;

	for (auto fiter = files.begin(); fiter != files.end(); ++fiter) {
		// create file input object from a filename using a file object factory
		input_learning_objects.push_back((AbstractLayerInputObject*)file_input_object_factory.newInstance(srcdir + *fiter));
	}

	// perform learning by updating existing vocabulary
	layers_learning->updateVocabulary(input_learning_objects);

	// save vocabulary
	// TODO: implement saving output vocabulary
	FileVocabularyOutputObject file_vocabulary_output(layers_learning->getOutputVocabulary());

	file_vocabulary_output.save(out_library);
	

	// delete factories and inference objects
	delete factory;
	delete layers_learning;
	delete config_ptr;
}

void LearningToolset::main(int argc, char* argv[]) {
	
	// parse input arguments
	const char* cfg_file = argv[2];
	const char* patt = argc > 3 ? argv[3] : "";
	const char* params = argc > 4 ? argv[4] : "";

	new_learning(argc, argv);
	
	///// old learning

	lib_learning_controler learning_controler;

	// all config values
    ConfigDictionary cfg_all(cfg_file);
	cfg_all.fromString(params);
	cfg_all.updateNamespaceReferences();

	ConfigDictionary cfg_learn;
	cfg_learn.fromNamespacePriority(cfg_all, 1, "learning");

	// initilaize controler 
	learning_controler.initialize(cfg_learn);

	base_deployer_mapreudce* default_mr = new openmp_deployer_mapreudce();

	learning_controler.set_mapreduce_implementation(default_mr);

	// load input data (layer1_result files for training, validation etc) based on action type
	const string action = learning_controler.get_action();
	if (action == "optimization") optimization_load_input(learning_controler, patt);			
    else if (action == "learn_objects") learn_objects_load_input(learning_controler, patt);
	else if (action == "learn_objects2") learn_objects2_load_input(learning_controler, patt);

	// perform library learning
	lib_learning_results* results = learning_controler.perform_learning();

	// save library 
	if (action == "optimization") optimization_save_results(results, learning_controler);
    else if (action == "learn_objects") learn_objects_save_results(results, learning_controler);
	else if (action == "learn_objects2") learn_objects2_save_results(results, learning_controler);

	// do cleanup
	learning_controler.cleanup();

	delete results;
	//delete default_mr;

}

