// layers.cpp : Defines the exported functions for the DLL application.
//

#include "layers.h"

#include "core/legacy/layer_1_result.h"

#include "core/legacy/inference/layer_1_creators.h"
#include "core/legacy/inference/layer_n_creators.h"

#include "core/legacy/learning/optimization.h"
#include "core/legacy/learning/layer_learning.h"

#include "utils/class_register.h"

namespace lay1create {

// export function for layer 1 (and its supporting functions)
///////////////////////////////////////////////////////////////////////////////


void create_layer1(vector<layer1_result*>& result, const ConfigDictionary& cfg, const vector<img>& srcimg, const void* mask, const bool using_regions)
{
    try {
        string type;
		string outdir = "";
		string libname = "";

        result.clear();
        cfg.getValue(type, "type", true);
		cfg.getValue(libname, "part_lib_name");
		cfg.getValue(outdir, "lib_out_dir");
		end_dir(outdir);	

		vector<layer1_result*> tmp_result(1);

		// create factory for layer1_creator
		layer1_creator::IFactory* factory = ClassRegister::get().retrieveFactory<layer1_creator::IFactory>(type);

		for (vector<img>::const_iterator img_it = srcimg.begin(); img_it != srcimg.end(); img_it++) {
			vector<layer1_result*> channel_result;

			// create appropriate class instance
			layer1_creator* cr = factory->newInstance();

			// if we have valid instance we need to init config and filter, and then create actual result 
			if (cr != nullptr) {
				cr->cfg_init(cfg);

				cr->init_filters();
				
				if (using_regions)
					cr->create_result_vector(channel_result, *img_it, *(vector<irectangle2>*)mask);
				else
					cr->create_result_vector(channel_result, *img_it, *(img*)mask);

				cr->clear_filters();

				// save part library is output path for it is valid				
				if (libname.size() > 0) {
					part_lib* library = cr->get_part_lib(cfg);
					library->save(outdir + libname);
					delete library;
				}

				delete cr;
			}
			if (result.empty()) {
				result = channel_result;
			} else {
				// merge with main results for each scale 
				for (int i = 0; i < result.size(); i++) {
					tmp_result[0] = channel_result[i];
					result[i]->add_results(tmp_result, 0);
				}
			}
		}
		
    } catch (const libhop_exception& e) {
		cout << e.what() << endl;
    } catch (const exception& e) {
		cout << "General exception '" << typeid(e).name() << "' with message: '" << e.what() << "'" << endl;
	} 
}

void create_layer1(vector<layer1_result*>& result, const ConfigDictionary& cfg, const vector<img>& srcimg, const img& mask_image) {
	create_layer1(result, cfg, srcimg, &mask_image, false);
}
void create_layer1(vector<layer1_result*>& result, const ConfigDictionary& cfg, const vector<img>& srcimg, const vector<irectangle2>& mask_regions) {
	create_layer1(result, cfg, srcimg, &mask_regions, true);
}

}

namespace layncreate {

// export function for layer n (and its supporting functions)
///////////////////////////////////////////////////////////////////////////////

clock_t create_layern(vector<layer1_result*>& result, const ConfigDictionary& cfg, part_lib* library)
{
    clock_t startt, endt, totalt = 0;

    try {

		int layer_index = cfg.getValueInt("layer_index", 0);
		int start_layer = -1, end_layer = -1;

		cfg.getValue(start_layer, "start_layer", layer_index > 0 ? false : true);
		cfg.getValue(end_layer, "end_layer", layer_index > 0 ? false : true);	

		// either ('layer_index') or ('start_layer' and 'end_layer') can be defined but not both
		if (layer_index > 0) {
			// if both start_layer and end_layer are in config and at the same time layer_index is defined, then return error
			if (start_layer > 0 || end_layer > 0) {
				throw new_libhop_exception("Error: Ambiguous definition of layer. Found 'layer_index' and 'start_layer','end_layer' in config. Use either 'layer_index' or 'start_layer' and 'end_layer'.");
			} else {
				// if layer_index was define, then start_layer and end_layer should not be in config, so continue normaly
				start_layer = layer_index;
				end_layer = layer_index;	
			}
		}

		// validate start and end layer values
		if (start_layer <= 1) {			
			cout  << "Warning: Found start_layer that is <= 1 but create_layern cannot process 1. layer ... setting start_layer to 2" << endl;
			start_layer = 2;
		}
		if (end_layer < start_layer) {
			throw new_libhop_exception("Found start_layer > end_layer. Unable to continue processing");
		}

		// verify that library is valid
		if (library == nullptr)
			throw new_libhop_exception("Invalid input library (library is nullptr)");

		vector<layern_creator*> layern_creator_list(end_layer+1, nullptr);

		// create creators for each layer
		for (int i = start_layer; i <= end_layer; ++i) {
			string ly_namespace = string("ly") + i;

			ConfigDictionary cfgi;
			cfgi.fromNamespacePriority(cfg, 1, ly_namespace.c_str());

			// construct creator and set proper library
			layern_creator_list[i] = new layern_creator(cfgi);
			layern_creator_list[i]->set_library(library);
		}		

		// process for each layer and each result 
		for (int i = start_layer; i <= end_layer; ++i) {			
			startt = clock();
			for (vector<layer1_result*>::iterator riter = result.begin(); riter != result.end(); ++riter) {
				layern_creator_list[i]->add_layer(*riter, i, 0);
			}
			endt = clock();
			totalt += endt - startt;

			if (!result.empty())
				cout << '(' << result[0]->cover_quotient(i - 2) << ')' << '[' << (double)(endt - startt)/CLOCKS_PER_SEC << " sec]";
		}

		cout << " (time = " << (double)totalt/CLOCKS_PER_SEC << ")" << endl;


		for (int i = start_layer; i <= end_layer; ++i) {
			layern_creator_list[i]->set_library(nullptr);
			delete layern_creator_list[i];
		}

    } catch (const libhop_exception& e) {
		cout << e.what() << endl;
		totalt = -1;
    } catch (const exception& e) {
		cout << "General exception '" << typeid(e).name() << "' with message: '" << e.what() << "'" << endl;
		totalt = -1;
	}
    return totalt;
}

}

namespace laylearning {

void lib_learning_controler::initialize(const ConfigDictionary& cfg_learn) {

	action = cfg_learn.getValueString("action", "optimization");

	// then get only values needed for specific action
	cfg.fromNamespacePriority(cfg_learn, 1, action.c_str());
}
     
void lib_learning_controler::cleanup(){
	for (auto item_iter = learning_data.begin(); item_iter != learning_data.end(); item_iter++) {
		for (auto data_iter = (*item_iter).begin(); data_iter != (*item_iter).end(); data_iter++) {
			if (data_iter->second != nullptr)
				delete data_iter->second;
		}
	}
	if (start_lib != nullptr)
		delete start_lib;

	if (dispose_deployer && mapreduce_deployer != nullptr)
		delete mapreduce_deployer;
}

void lib_learning_controler::add_learning_data(const learning_data_item& data) {
	learning_data.push_back(data);
}
      

// Optimizes layers from start_layer to min_layer
lib_learning_results* lib_learning_controler::perform_optimization()
{

	string type = cfg.getValueString("optimization_type", "default"); // this is not used any more
	clock_t start, end;

	start = clock();
		
	cfg.fromFile(cfg.getValueString("import_configuration", "").c_str());

	optimization_data odata(cfg, this->start_lib);

	if (this->learning_data.empty())
		throw custom_libhop_exception(ConfigException, string("Can not perform optimization on an empty set"));
		
	// load input data
	for (auto iter = learning_data.begin(); iter != learning_data.end(); iter++) {
		streamed_pointer* ptr = static_cast<streamed_pointer*>(iter->at(INPUTDATA_IMAGE));

		if (ptr == nullptr) continue;

		const char* filename = (const char*)((*iter)[INPUTDATA_FILENAME]);
		cout << "Processing " << (filename != nullptr ? filename : "<missing filename>") << endl;

		odata.add_to_test_set(*ptr);
	}

	int overall_steps = cfg.getValueInt("overall_steps", 2);
	vector<cr_optimization_base*> lopts;
	int final_layer;

	vector<int> optseq;
	if (cfg.isDefined("start_layer")) {
		int start_layer = cfg.getValueInt("start_layer", 1);
		int end_layer = cfg.getValueInt("end_layer", 2);

		for (int l = start_layer; l <= end_layer; ++l) {
			optseq.push_back(l);
		}
	} else {
		cfg.getValue(optseq, "optimization_sequence", true);		
	}

	for (int i = 0; i < (int)optseq.size(); ++i) {
		string ly_namespace = string("ly") + optseq[i];
		ConfigDictionary ly_cfg;
		ly_cfg.fromNamespacePriority(cfg, 1, ly_namespace.c_str());
		lopts.push_back(new cr_layer_optimization(&odata, ly_cfg, optseq[i]));
	}
	final_layer = optseq.back();

	cr_set_optimization cropt(&odata, overall_steps, lopts);

	// Fill the optimizer 
	cropt.set_library(this->start_lib);
	cropt.add_to_init_test_set(&odata);

	cropt.set_mapreduce(this->mapreduce_deployer);

	// Execute 
	cropt.execute();
    
	// Get the results	
	list<streamed_pointer> ll;
	lib_optimization_results* results = new lib_optimization_results(cropt.get_best_library(), ll, final_layer);

	list<streamed_pointer>& ll1 = results->processed_layers;
	cropt.get_best_test_set(ll1);

	end = clock();

	cout << "Processing time: " << (double)(end - start)/CLOCKS_PER_SEC << endl;
	
	return results;
}

lib_learning_results* lib_learning_controler::learn_object_parts()
{
	obj_learning learner(cfg);

	learner.set_library(this->start_lib);

	string libimgfile = cfg.getValueString("library_image_file", "");
	string vcreatorcfg = cfg.getValueString("validation_creator", "");
	string catname = cfg.getValueString("category_name", "");

	layer1_result* res = nullptr;
	K_bin bin(12, 2, 4, 7);

	// load validation files
	for (auto iter = learning_data.begin(); iter != learning_data.end(); ++iter) {
		streamed_pointer* ptr = (streamed_pointer*)((*iter)[INPUTDATA_VALIDATION_POS_IMAGE]);
		
		if (ptr == nullptr) continue;

		const char* filename = (const char*)((*iter)[INPUTDATA_FILENAME]);
		cout << "Loading positive validation file " << (filename != nullptr ? filename : "<missing filename>") << endl;

		list<irectangle2>* gtrs = (list<irectangle2>*)(iter->at(INPUTDATA_GROUNDTRUTH));
		learner.add_validation_data(*ptr, gtrs->front()); 
	}
	for (auto iter = learning_data.begin(); iter != learning_data.end(); ++iter) {
		streamed_pointer* ptr = (streamed_pointer*)((*iter)[INPUTDATA_VALIDATION_NEG_IMAGE]);

		if (ptr == nullptr) continue;

		const char* filename = (const char*)((*iter)[INPUTDATA_FILENAME]);
		cout << "Loading negative validation file " << (filename != nullptr ? filename : "<missing filename>") << endl;

		learner.add_validation_data(*ptr, irectangle2()); 
	}


	// learn objects
	int count = 0;

	for (auto iter = learning_data.begin(); iter != learning_data.end(); ++iter) {
		streamed_pointer* ptr = static_cast<streamed_pointer*>((*iter)[INPUTDATA_IMAGE]);
		
		if (ptr == nullptr) continue;

		res = (layer1_result*)ptr->get();

		const char* filename = (const char*)((*iter)[INPUTDATA_FILENAME]);

		if (res == nullptr) {
			cout << "Skipping " << (filename != nullptr ? filename : ptr->get_name_only()) << " as file is missing or invalid" << endl;
			continue;
		} else {
			cout << "Processing  " << (filename != nullptr ? filename : "<missing filename>") << endl;
		}

		if (res) {
			list<irectangle2>& gtrs = *static_cast<list<irectangle2>*>(iter->at(INPUTDATA_GROUNDTRUTH));
			clock_t start, end;
			scmap_t scmap;

			get_sc_map(scmap, res, bin, true);
			start = clock();
			learner.object_from_result(res, scmap, gtrs.empty() ? irectangle2() : gtrs.front());
			end = clock();
			cout << "Processing time: " << (double)(end - start)/CLOCKS_PER_SEC << endl;

			delete res;
		}
	}
	int howmany = learner.add_to_library(catname);

	cout << howmany << " object parts added." << endl;

	// copy address to created library
	part_lib* lib = learner.library;

	// remove lib reference from learner (since learner gets deleted after end of this call, but library does not yet gets saved)
	learner.library = nullptr;

	return new lib_learning_results(lib);
}

lib_learning_results* lib_learning_controler::learn_object_parts2()
{
	string pattern = ""; // TODO : FIX THIS 

	cout << "Doing object part learning (learn_objects2)" << endl << endl;
	
	o_learning olearner(cfg);

	string catname;
	list<string> files;

	cfg.getValue(catname, "category_name", true);

	// Load train files 1st time

	cout << "Duplet Statistics" << endl;

	for (auto iter = learning_data.begin(); iter != learning_data.end(); ++iter) {
		streamed_pointer* ptr = static_cast<streamed_pointer*>((*iter)[INPUTDATA_IMAGE]);
		
		if (ptr == nullptr) continue;

		layer1_result* res = (layer1_result*)ptr->get();

		const char* filename = (const char*)((*iter)[INPUTDATA_FILENAME]);

		if (res == nullptr) {
			cout << "File " << (filename != nullptr ? filename : ptr->get_name_only()) << " can not be opened!" << endl;
		} else {
			cout << "Processing file " << (filename != nullptr ? filename : ptr->get_name_only());

			olearner.update_duplet_statistics(res);

			delete res;
			cout << " done." << endl;
		}
	}

	// Validation files

	cout << "Validation sets" << endl;

	// "Self"-validation
	if (cfg.getValueBool("self_validation", false)) {
		for (auto iter = learning_data.begin(); iter != learning_data.end(); ++iter) {
			streamed_pointer* ptr = static_cast<streamed_pointer*>((*iter)[INPUTDATA_IMAGE]);
		
			if (ptr == nullptr) continue;

			layer1_result* res = (layer1_result*)ptr->get();

			const char* filename = (const char*)((*iter)[INPUTDATA_FILENAME]);

			if (res == nullptr) {
				cout << "File " << (filename != nullptr ? filename : ptr->get_name_only()) << " can not be opened!" << endl;
			} else {
				cout << "Processing file " << (filename != nullptr ? filename : ptr->get_name_only());

				list<irectangle2>* gtrs = static_cast<list<irectangle2>*>((*iter)[INPUTDATA_GROUNDTRUTH]);

				olearner.add_to_validation_set(res, *gtrs);

				delete res;
				cout << " done." << endl;
			}
		}

	}

	// load validation files with groundtruth
	files.clear();
	for (auto iter = learning_data.begin(); iter != learning_data.end(); ++iter) {
		streamed_pointer* ptr = static_cast<streamed_pointer*>((*iter)[INPUTDATA_VALIDATION_POS_IMAGE]);
		
		if (ptr == nullptr) continue;

		layer1_result* res = (layer1_result*)ptr->get();

		const char* filename = (const char*)((*iter)[INPUTDATA_FILENAME]);
		if (res == nullptr) {
			cout << "File " << (filename != nullptr ? filename : ptr->get_name_only()) << " can not be opened!" << endl;
		} else {
			cout << "Processing file " << (filename != nullptr ? filename : ptr->get_name_only());
			list<irectangle2>* gtrs = static_cast<list<irectangle2>*>((*iter)[INPUTDATA_GROUNDTRUTH]);

			olearner.add_to_validation_set(res, *gtrs);

			delete res;
			cout << " done." << endl;
		}
	}
    
	// Make models

	cout << "Make models" << endl;

	K_bin bin(12, 2, 4, 7);
    
	files.clear();
	for (auto iter = learning_data.begin(); iter != learning_data.end(); ++iter) {
		streamed_pointer* ptr = static_cast<streamed_pointer*>((*iter)[INPUTDATA_IMAGE]);
		
		if (ptr == nullptr) continue;

		layer1_result* res = (layer1_result*)ptr->get();

		const char* filename = (const char*)((*iter)[INPUTDATA_FILENAME]);
		if (res == nullptr) {
			cout << "File " << (filename != nullptr ? filename : ptr->get_name_only()) << " can not be opened!" << endl;
		} else {
			scmap_t scmap;

			cout << "Processing file " << (filename != nullptr ? filename : ptr->get_name_only());

			get_sc_map(scmap, res, bin, true);
			olearner.learn_models(res, scmap);

			delete res;
			cout << " done." << endl;
		}
		//break;
	}
	olearner.finalize();
	
	return new lib_learning_results(olearner.get_library());
}

}