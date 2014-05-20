
#include "file_input.h"

void register_file_input_classes() {
	/************************************************************************* 
	 *           !!!! DO NOT MODIFY OR REMOVE THIS FUNCTION   !!!!!
	 ************************************************************************/

	ClassRegister::registerFactory<IFileInput::IFactory, FileInputObjectWithGroundtruth::Factory>();
	ClassRegister::registerFactory<IFileInput::IFactory, FileImageInputObject::Factory>();
	ClassRegister::registerFactory<IFileInput::IFactory, FileLayerInputObject::Factory>();
	ClassRegister::registerFactory<IFileInput::IFactory, FileVocabularyInputObject::Factory>();

}

//// FileInputObjectWithGroundtruth methods
//////////////////////////////////////////////////

std::vector<GroundtruthObject> FileInputObjectWithGroundtruth::getGroundtruths(const std::string& only_label) const {
	std::vector<GroundtruthObject> result;

	
    std::ifstream is(groundtruth_filename.c_str());

	bool warning_displayed = false;
    if (is.is_open()) {
        while (true) {
			string line;
			// read and parse whole line as stream of string
			getline(is,line);
			stringstream sline(line);

			rectangle2<double> rect;
			string label;

			// first read rectangle values
			sline >> rect;
            if (is.fail()) break;
			
			// then read category name and trim it
			sline >> label;
			label = trimString(label);
			
			if (only_label.length() > 0 && label != only_label) {
				// skip all ground truth values that do not belog to category cat_name_only
				continue;
			}
			result.push_back(GroundtruthObject());

			result.back().region = irectangle2((int)rect.ll.x, (int)rect.ll.y, (int)rect.ur.x, (int)rect.ur.y);
			result.back().label = label;
        }
	} else {
		cout << "\nNo groundtruth file found" << endl;
	}

	return result;
}

