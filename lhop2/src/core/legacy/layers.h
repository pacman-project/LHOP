// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the LAYERS_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// LAYERS_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
//#ifdef LAYERS_EXPORTS
//#define LAYERS_API __declspec(dllexport)
//#else
//#define LAYERS_API __declspec(dllimport)
//#endif

#pragma once
#ifndef _LAYERS_
#define _LAYERS_

#include <vector>
#include <string>
#include <ctime>

#include "utils/matrix.h"
#include "utils/mapreduce.h"
#include "utils/serialization/streamed_pointer.h"
#include "utils/config_dictionary.h"

#include "core/legacy/layer_1_result.h"
#include "core/legacy/inference/layer_n_creators.h"

// exported functions
///////////////////////

namespace lay1create {
	/// Create layer 1 from image
	/// \param result a list of processed images (scales)
	/// \param cfg a configuration dictionary
	/// \param srcimg image to process
	/// \param maskimg mask image (may be empty)
	void create_layer1(vector<layer1_result*>& result, const ConfigDictionary& cfg, const vector<img>& srcimg, const img& mask_image);
	void create_layer1(vector<layer1_result*>& result, const ConfigDictionary& cfg, const vector<img>& srcimg, const vector<irectangle2>& mask_regions);
}

namespace layncreate {
	clock_t create_layern(vector<layer1_result*>& result, const ConfigDictionary& cfg, part_lib* library);
}

namespace laylearning {

// main control class for library learning 
////////////////////////////////////////////////
class lib_learning_results {
public:
	part_lib* lib;

	lib_learning_results(part_lib* lib) : lib(lib) {}
	virtual ~lib_learning_results() {}
};


class lib_learning_controler {
public:
	typedef map<int,void*> learning_data_item;

	enum { INPUTDATA_IMAGE = 0, INPUTDATA_GROUNDTRUTH = 1, INPUTDATA_VALIDATION_POS_IMAGE = 2, INPUTDATA_VALIDATION_NEG_IMAGE = 3, INPUTDATA_FILENAME = 4 };

private:
	string action;
	ConfigDictionary cfg;

	vector<learning_data_item> learning_data;

	part_lib* start_lib;

	base_deployer_mapreudce* mapreduce_deployer;
	bool dispose_deployer;
public:
	lib_learning_controler() : start_lib(nullptr), mapreduce_deployer(nullptr), dispose_deployer(1)  {}

	void initialize(const ConfigDictionary& cfg_learn);     
	void cleanup();

	virtual void add_learning_data(const learning_data_item& data) ;
	  
	virtual void set_starting_library(part_lib* lib) { this->start_lib = lib; }

	ConfigDictionary& get_cfg() { return cfg; }
	string get_action() { return this->action; }

	lib_learning_results* perform_learning() {
		// do learning based on action defined in config 
		if (action == "optimization") return this->perform_optimization();
        else if (action == "learn_objects") return this->learn_object_parts();
        else if (action == "learn_objects2") return this->learn_object_parts2();
	}	

	virtual void set_mapreduce_implementation(base_deployer_mapreudce* deployer, bool dispose_on_cleanup = 0) {
		this->mapreduce_deployer = deployer;
	}	

private:
	/**
	 * Main learning method learns from 2. to n-th layer using optimization process
	 */ 
	lib_learning_results* perform_optimization();
	/**
	 * Old code for object layer learning (n+1-th layer i.e. last layer)
	 */
	lib_learning_results* learn_object_parts();
	/**
	 * Latest code for object layer learning (n+1-th layer i.e. last layer)
	 */
	lib_learning_results* learn_object_parts2();
};

// results classes for different type of learning 
////////////////////////////////////////////////


class lib_optimization_results : public lib_learning_results {
public:
	list<streamed_pointer> processed_layers;
	int final_layer_number;

	lib_optimization_results(part_lib* lib, list<streamed_pointer> processed_layers, int layer) : lib_learning_results(lib), processed_layers(processed_layers), final_layer_number(layer) {}
	virtual ~lib_optimization_results() {
	}
};

}

#endif