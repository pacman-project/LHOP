/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */

#ifdef WIN32
#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
#include <windows.h>
#endif

#include "hop.h"

#include <stdio.h>

#include <ctime>
#include <tuple>
#include "layers/layers.h"
#include "layers/optimization.h"

using namespace laylearning;

void optimization_load_input(lib_learning_controler &learning_controler, const char* pattern) {

	// get configuration
	const config_dictionary& cfg =  learning_controler.get_cfg();

	// read library
    string libname;

	part_lib* library = nullptr;

    cfg.get_value(libname, "library", true);
    read_library(libname, library);

	if (library == nullptr)
		throw custom_libhop_exception(config_exception, string("Library '" + libname + "' not found"));

	learning_controler.set_starting_library(library);

	// read input files
	string dir = cfg.get_value_string("src_dir", "");
	list<string> files;
	layer1_result* res;
	int count = 0;
	int maxcount = cfg.get_value_int("file_limit", INT_MAX);
	int start_layer = cfg.get_value_int("start_layer", 1);
	int end_layer = cfg.get_value_int("end_layer", 2);

	end_dir(dir);
	if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
		list_directory(files, dir + pattern);
	cout << "from_file = " << (cfg.get_value_bool("from_file", false) ? "true" : "false") << endl;
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

	const config_dictionary& cfg = learning_controler.get_cfg();

	lib_optimization_results* opt_results = static_cast<lib_optimization_results*>(results);

	int final_layer = opt_results->final_layer_number;

	// Save results
    string lib_export_name = cfg.get_value_string("lib_export_name", "lib");
    string res_export_name = cfg.get_value_string("res_export_name", "res");
	string name = cfg.get_value_string("out_library", lib_export_name + (final_layer) + ".plb");
    int i = 0;

    PRINT_INFO("");
    PRINT_INFO("SAVING");

	opt_results->lib->save(name);

    if (cfg.get_value_bool("save_test_set", false)) {

		list<streamed_pointer>& sresult = opt_results->processed_layers;

        for (list<streamed_pointer>::iterator iter = sresult.begin(); iter != sresult.end(); ++iter) {
            layer1_result* res = (layer1_result*)iter->get();

            name = res_export_name + (i++) + string(".ly") + (final_layer );
            PRINT_INFO(name);
            res->delete_edges_complement(atom("toPrevLayer").get_index());
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
 * learn_objects must recieve list of following input data:
 *  - (INPUTDATA_VALIDATION_POS_IMAGE, INPUTDATA_GROUNDTRUTH, INPUTDATA_FILENAME) positive validation file
 *  - (INPUTDATA_VALIDATION_NEG_IMAGE, INPUTDATA_FILENAME) negative validation file
 *  - (INPUTDATA_IMAGE, INPUTDATA_GROUNDTRUTH, INPUTDATA_FILENAME) main training images
 * All images are layer1_results in streamable_pointer* form
 */
void learn_objects_load_input(lib_learning_controler &learning_controler, const char* patt) {
	
	const config_dictionary& cfg =  learning_controler.get_cfg();

	string pattern(patt);

	string dir = cfg.get_value_string("src_dir", "");
	string vdir = cfg.get_value_string("validation_src_dir", "");
	string vfilepos = cfg.get_value_string("positive_validation_files", "");
	string vfileneg = cfg.get_value_string("negative_validation_files", "");
	string libimgfile = cfg.get_value_string("library_image_file", "");
	int rndchoice = cfg.get_value_int("random_select", -1);
	string catname = cfg.get_value_string("category_name", "");
	int maxcount = cfg.get_value_int("file_limit", INT_MAX);

	// read library
    string libname;

	part_lib* library = nullptr;

    cfg.get_value(libname, "library", true);
    read_library(libname, library);

	if (library == nullptr)
		throw custom_libhop_exception(config_exception, string("Library '" + libname + "' not found"));

	learning_controler.set_starting_library(library);


	list<string> files;
	list<string> vfilespos;
	list<string> vfilesneg;
	layer1_result* res = nullptr;

	end_dir(dir);
	end_dir(vdir);
	if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
		list_directory(files, dir + pattern);
	if (rndchoice > 0) random_selection(files, rndchoice, cfg.get_value_string("dropped_files", ""));
	if (!vfilepos.empty() || !vfileneg.empty()) {
		file_list_from_file(vfilespos, vfilepos, vdir, cfg.get_value_string("positive_pattern", ""));
		file_list_from_file(vfilesneg, vfileneg, vdir, cfg.get_value_string("negative_pattern", ""));		
	}
	cout << "from_file = " << (cfg.get_value_bool("from_file", false) ? "true" : "false") << endl;
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
	
	const config_dictionary& cfg = learning_controler.get_cfg();
	
	string libname = cfg.get_value_string("out_library", "olib.plb");
	
	results->lib->save(libname);

	int layer;
	cfg.get_value(layer, "layer", true); // in case obj_learning::init_cfg will change layer definition then this has to be changed too
	string libimgfile = cfg.get_value_string("library_image_file", "");
	
	if(!libimgfile.empty()) {
		cout << "Saving lib image: layer " << layer << endl;
		for(int il=1; il<=layer+3; il++) {
			string str;
			str += (string)"-" + il + ".png";
			cout << "Saving lib image (layer " << il << ") to " << change_extension(libimgfile, str).c_str() << endl;
			results->lib->save_all(change_extension(libimgfile, str).c_str(), il, 0, -1,
				cfg.get_value_bool("show_labels", false));
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

	const config_dictionary& cfg = learning_controler.get_cfg();

	string dir = cfg.get_value_string("src_dir", "");
	string vdir = cfg.get_value_string("validation_src_dir", "");
	string vfiles = cfg.get_value_string("validation_files", "");
	
	string catname;
	list<string> files;

	cfg.get_value(catname, "category_name", true);

	end_dir(dir);
	if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
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
				gtrs->push_back(node_set_bounding_rectangle(res->shape_nodes[0].begin(), res->shape_nodes[0].end()));

			delete res;
		}	

		data_item[lib_learning_controler::INPUTDATA_GROUNDTRUTH] = gtrs;

		learning_controler.add_learning_data(data_item);
	}

	files.clear();
	file_list_from_file(files, vfiles, vdir, cfg.get_value_string("validation_pattern", ""));
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
	const config_dictionary& cfg = learning_controler.get_cfg();

	string libname = cfg.get_value_string("out_library", "olib.plb");

	results->lib->save(libname);
}

HOPINTERFACE_API bool hop_learning(const char* cfg_file, const char* patt, const char* params)
{
	bool ok = true;
    try {
		lib_learning_controler learning_controler;

		// all config values
        config_dictionary cfg_all(cfg_file);
		cfg_all.from_string(params);
		cfg_all.update_namespace_references();

		config_dictionary cfg_learn;
		cfg_learn.from_namespace_priority(cfg_all, 1, "learning");

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

    } catch (const libhop_exception& e) {
		cout << e.what() << endl;
		ok = false;
    } catch (const exception& e) {
		cout << "General exception '" << typeid(e).name() << "' with message: '" << e.what() << "'" << endl;
		ok = false;
	}
	return ok;
}
