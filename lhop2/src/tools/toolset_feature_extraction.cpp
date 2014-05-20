// feature_extraction.cpp : Defines the entry point feature extraction functionality
//

#include <stdio.h>
#include <string>
#include <cctype>

#include "core/feature_extraction/abstract_feature_extraction.h"

#include "tools/main_toolset.h"

#include "tools/file_input.h"
#include "tools/file_output.h"

#include "utils/misc.h"
#include "utils/configuration.h"

void extract_features_from_files(const IConfiguration& config, const char* pattern) {

	// load input (source) and output dir
	string srcdir = fixPathEndSeparator(config.getString("src_dir"));
	string outdir = fixPathEndSeparator(config.getString("out_dir"));
	string outprefix = config.getString("out_prefix");

	list<string> files;

	// get the list of input files from reading text file that contians input filename in each line (if from_file is set to true)
	if (!config.getBool("from_file", false) || !list_from_file(files, pattern, srcdir))
		list_directory(files, srcdir + "/" + string(pattern));

    cout << "Creating layer 1 from '" << srcdir << "'" << endl;
    cout << "  to '" << outdir << "'" << endl;
    cout << files.size() << " files found (from file = " << (config.getBool("from_file", false) ? "true" : "false") <<
        ", pattern = " << pattern << ")" << endl;

	string feature_extraction_type = config.getString("type");

	// legacy support
	/////////////////////////////////////////////
	if (feature_extraction_type.compare("struct") == 0) 
		feature_extraction_type = "gabor";
	else if (feature_extraction_type.compare("app") == 0) 
		feature_extraction_type = "gabor.app";
	else if (feature_extraction_type.compare("dog") == 0) 
		feature_extraction_type = "gabor.dog";
	else if (feature_extraction_type.compare("loggabor") == 0) 
		feature_extraction_type = "gabor.loggabor";
	/////////////////////////////////////////////

	// get factory for specific implementation of feature extraction (shape, motion, 3d, color, texture ...)
	AbstractFeatureExtraction::IFactory* factory = ClassRegister::get().retrieveFactory<AbstractFeatureExtraction::IFactory>(feature_extraction_type);

	if (factory == nullptr)
		throw new_libhop_exception("Unable to find factory object for the requested feature extraction type (" + feature_extraction_type + ")");

	// find which input type is associated with this implementation e.g. it can be ImageInput, 3DInput, MotionInput ... (factory should provide the necessary information)
	string input_object_type = factory->getAssociatedInputObjectType();

	// get File factory for associated input object i.e. retrive factory that will be able to create object which have IFileInput interface and are of correct input type for associated feature extraction
	// [IFileInput + feature extraction type (shape, motion, 3d ..)] should return factories such as FileImageInput, File3DInput, FileMotionInput ...
	IFileInput::IFactory* file_input_factory = ClassRegister::get().retrieveFactory<IFileInput::IFactory>(input_object_type);

	if (file_input_factory == nullptr)
		throw new_libhop_exception("Unable to find FileInput factory object for input type (" + input_object_type + ") associated with the requested feature extraction type (" + feature_extraction_type + ")");

	
	// create class for doing actual feature extraction
	// use factory based on feature extraction type (shape, motion, 3d, color, texture ...)
	AbstractFeatureExtraction* feature_extraction = factory->newInstance(config);

	// process each file
	for (auto fiter = files.begin(); fiter != files.end(); ++fiter) {

		// create file input object from filename
		// the file_input_factory must create correct input type for feature extraction type (FileImageInput, File3DInput or FileMotionInput ..)
		// but we can directly cast it into AbstractInputObject* as feature_extraction only knows about AbstractInputObject
		AbstractInputObject* input_object = dynamic_cast<AbstractInputObject*>(file_input_factory->newInstance(srcdir + *fiter));
		
		// first handle pre-processing (i.e. spliting images into RGB, spliting into tiles, resizeing, bluring, scaling ...)
		// we take one input object created from file and we can produce multiple new input objects
		// preprocessing steps are controled by the AbstractFeatureExtraction implementation (its Factory class should have constructed proper preprocessing objects from the configuration)
		std::vector<AbstractInputObject*> preprocessed_input_objects =  feature_extraction->performPreprocessing(dynamic_cast<AbstractInputObject*>(input_object)); 

		std::vector<AbstractOutputObject*> output_objects;
		// do actual feature extraction for each pre-processed input object
		for (auto input_iter = preprocessed_input_objects.begin(); input_iter != preprocessed_input_objects.end(); ++input_iter) {
			// perform extraction
			LayerOutputObject* output_layer = feature_extraction->performExtraction(*input_iter);

			// save output as <input,output> tuple
			output_objects.push_back(output_layer);

			// delete input (but only if it is not the same as initial input_object created before pre-processing)
			if (*input_iter != input_object)
				delete *input_iter;
		}

		// perform any post-processing (e.g. combining back from tiles, combining from RGB ...)
		output_objects =  feature_extraction->performPostprocessing(output_objects);

		// save the results into proper output folder
		for (auto output_iter = output_objects.begin(); output_iter != output_objects.end(); ++output_iter) {
			// do dynamic cast just in case improper post-processing might have been used
			LayerOutputObject* output_layer = dynamic_cast<LayerOutputObject*>(*output_iter); 

			// create file output object from result
			FileLayerOutputObject file_output_layer(output_layer);

			// save to output dir
			file_output_layer.save(outdir + outprefix + *fiter);

			// delete output
			delete output_layer;
		}

		// clear memory
		delete input_object;
	}

	delete feature_extraction;
}

string FeatureExtractionToolset::getShortDescription() {
	return "";
}
string FeatureExtractionToolset::getLongDescription() {
	return "";
}
string FeatureExtractionToolset::getUsageDescription() {
	return "";
}

bool FeatureExtractionToolset::areArgumentsValid(int argc, char* argv[]) {
	// TODO: do more inteligent verification
	return argc > 2  ? true : false;
}

void FeatureExtractionToolset::main(int argc, char* argv[]) {

	// parse input arguments
	const char* cfgfile = argv[2];
	const char* patt = argc > 3 ? argv[3] : "";
	const char* params = argc > 4 ? argv[4] : "";

	DictionaryConfiguration config;

	// load from string first but do not update namespace references (i.e. leave .from_namespace = ... values as is)
	config.appendConfig(StringDictionaryConfiguration(params)); 
	// then also load from file (append will NOT override current values) 
	// and update namespace references at the same time (i.e. replace .from_namespace = ... with the values from the actual namespace)
	config.appendConfig(FileDictionaryConfiguration(cfgfile));  
	
	// get inference.ly1 namespace
	IConfiguration* inference_ly1_config = config.getNamespace("inference.ly1");
	
	if (inference_ly1_config == nullptr)
		throw new_libhop_exception("Missing inference.ly1 namespace in configuration");

	extract_features_from_files(*inference_ly1_config, patt);
	
	delete inference_ly1_config;
}