// inference.cpp : Defines the entry point feature extraction functionality
//



#include <stdio.h>
#include <string>
#include <cctype>

#include "core/inference/abstract_inference.h"

#include "tools/file_input.h"
#include "tools/file_output.h"

#include "tools/main_toolset.h"

#include "utils/utils.h"
#include "utils/configuration.h"


string InferenceToolset::getShortDescription() {
	return "Infere vocabulary part from existing layer object with N-layers";
}
string InferenceToolset::getLongDescription() {
	return "TODO";
}
string InferenceToolset::getUsageDescription() {
	return "TODO";
}

bool InferenceToolset::areArgumentsValid(int argc, char* argv[]) {
	// TODO: do more inteligent verification
	return argc > 2  ? true : false;
}

void inferre_layer_from_files(const IConfiguration& config, const char* pattern) {
	
	std::string srcdir = fixPathEndSeparator(config.getString("src_dir"));
	std::string outdir = fixPathEndSeparator(config.getString("out_dir"));

	std::list<std::string> files;

	if (!config.getBool("from_file") || !list_from_file(files, pattern, srcdir))
		list_directory(files, srcdir + "/" + string(pattern));
	
	std::string library_name = config.getString("part_lib_name");

	// load vocabualry from file
	FileVocabularyInputObject vocabulary(library_name);

	// create multilayer inference that will be able to handle each individual layer 
	// (irrespectevly of the type of inference - configuration will define which type will be used)
	// TODO: we can later use specific multi-layer inferences (e.g. for specific modalities such as textures, 3D, color etc.) if requied by the library
	//		 but for now just simply use default MultipleLayersInference

	std::string multiple_layer_inference_type = MultipleLayersInference::Factory().assignedRegistrationName(); // just find the registration name of the MultipleLayersInference

	// find factory for specific inference type
	AbstractLayerInference::IFactory* factory = ClassRegister::get().retrieveFactory<AbstractLayerInference::IFactory>(multiple_layer_inference_type);

	// create inference object (will load all the neccessary settings)
	AbstractLayerInference* layers_inference = factory->newInstance(config, vocabulary.getVocabularyObject());

	FileLayerInputObject::Factory file_input_object_factory;

	// Process files
	for (auto fiter = files.begin(); fiter != files.end(); ++fiter) {
		
		// create file input object from filename
		FileLayerInputObject* file_input_object = (FileLayerInputObject*)file_input_object_factory.newInstance(srcdir + *fiter);

		// call inference on the input object
		LayerOutputObject* output_tree = layers_inference->performInference(file_input_object);
		
		// create file output object from result
		FileLayerOutputObject file_output_tree(output_tree);

		// save to output dir
		file_output_tree.save(outdir + *fiter);

		// delete input and output
		delete file_input_object;
		delete output_tree;
	}

	// delete factories and inference objects
	delete factory;
	delete layers_inference;
}

/// Inference of *multiple layers* on files described by pattern -- using config cfg
void InferenceToolset::main(int argc, char* argv[]) {
	
	// parse input arguments
	const char* cfgfile = argv[2];						// path to configuration file
	const char* fpatt = argc > 3 ? argv[3] : "";		// optional pattern
	const char* params = argc > 4 ? argv[4] : "";		// optional string of additional configuration values

	DictionaryConfiguration config;

	// load from string first but do not update namespace references (i.e. leave .from_namespace = ... values as is)
	config.appendConfig(StringDictionaryConfiguration(params)); 
	// then also load from file (append will NOT override current values) 
	// and update namespace references at the same time (i.e. replace .from_namespace = ... with the values from the actual namespace)
	config.appendConfig(FileDictionaryConfiguration(cfgfile));  
	
	// get inference.ly1 namespace
	IConfiguration* inference_config = config.getNamespace("inference");
	
	if (inference_config == nullptr)
		throw new_libhop_exception("Missing inference namespace in configuration");

	inferre_layer_from_files(*inference_config, fpatt);

	delete inference_config;
}