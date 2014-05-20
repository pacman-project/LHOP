// dictionary for configuration propreties

#include "utils/convert.h"
#include "utils/misc.h"
#include "configuration.h"

#include <cstring>

using namespace std;


std::vector<std::string> splitString(const string& str, char c) {
    string::size_type start = 0, size = str.size(), pos;
    
    std::vector<std::string> result;
    while (start < size && (pos = str.find_first_of(c, start)) != string::npos) {
        result.push_back(str.substr(start, pos - start));
        start = pos + 1;
    }
    if (start < size) 
		result.push_back(str.substr(start, size - start));

	return result;
}

// ConfigDictionary
//////////////////////

bool DictionaryConfiguration::hasKey(const std::string& key) const {
	return (dictionary.find(key) != dictionary.end()) ? true : false;
}

/// Getter/setters for integer
////////////////////////////////////
int DictionaryConfiguration::getInt(const std::string& key) const {
    std::string value_str = getString(key);

    try {
        return Convert::ToInteger<int>(value_str.c_str(), (unsigned)value_str.size());
    } catch (InvalidConversionException) {
        throw custom_libhop_exception(ConfigException_, "Invalid value for key '" + key + " (expecting int)'");
    }
}

void DictionaryConfiguration::setInt(const std::string& key, int val) {
    ostringstream result;
    result << val;
    setString(key, result.str());
}

std::vector<int> DictionaryConfiguration::getIntArray(const std::string& key) const {
    std::string value_str = getString(key);

	try {
		std::vector<int> result;
        std::vector<std::string> strings = splitString(value_str, ',');
        
        for (size_t i = 0; i < strings.size(); ++i)
            result.push_back(Convert::ToInteger<int>(strings[i].c_str(), (unsigned)strings[i].size()));
		
		return result;
    } catch (InvalidConversionException) {
        throw custom_libhop_exception(ConfigException_, "Invalid value for key '" + key + "' (expecting int)");
    }    
}

/// Getter/setters for double
////////////////////////////////////


double DictionaryConfiguration::getDouble(const std::string& key) const {
    std::string value_str = getString(key);

    try {
        return Convert::ToDouble(value_str.c_str(), (unsigned)value_str.size());;
    } catch (InvalidConversionException) {
        throw custom_libhop_exception(ConfigException_, "Invalid value for key '" + key + "' (expecting double)");
    }
}

void DictionaryConfiguration::setDouble(const std::string& key, double val) {
    ostringstream result;
    result << val;
    setString(key, result.str());
}

std::vector<double> DictionaryConfiguration::getDoubleArray(const std::string& key) const {
    std::string value_str = getString(key);
     
    try {
		std::vector<double> result;
        std::vector<std::string> strings = splitString(value_str, ',');

        for (size_t i = 0; i < strings.size(); ++i)
            result.push_back(Convert::ToDouble(strings[i].c_str(), (unsigned)strings[i].size()));

        return result;
    } catch (InvalidConversionException) {
        throw custom_libhop_exception(ConfigException_, "Invalid value for key '" + key + " (expecting double)'");
    }    
}


/// Getter/setters for boolean
////////////////////////////////////


bool DictionaryConfiguration::getBool(const std::string& key) const {
	std::string value_str = getString(key);

	if (value_str == "true" || value_str == "false")
		return value_str == "true";
	else
		throw custom_libhop_exception(ConfigException_, "Invalid value for key '" + key + " (expecting bool as 'true' or 'false')'");
}


void DictionaryConfiguration::setBool(const std::string& key, bool val) {
    setString(key, val ? "true" : "false");
}

std::vector<bool> DictionaryConfiguration::getBoolArray(const std::string& key) const {
	std::string value_str = getString(key);
    
	std::vector<bool> result;
    std::vector<std::string> strings = splitString(value_str, ',');

    for (size_t i = 0; i < strings.size(); ++i) {		
		if (strings[i] == "true" || strings[i] == "false")
			result.push_back(strings[i] == "true");
		else
			throw custom_libhop_exception(ConfigException_, "Invalid value for key '" + key + " (expecting bool as 'true' or 'false')'");
    }

    return result;    
}

/// Getter/setters for string
////////////////////////////////////

std::string DictionaryConfiguration::getString(const std::string& key) const {
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end())
        throw custom_libhop_exception(ConfigException_, "Key '" + key + "' not specified");

    return iter->second;
}

void DictionaryConfiguration::setString(const std::string& key, const std::string& val) {
	map<string,string>::iterator iter;

	if ((iter = dictionary.find(key)) == dictionary.end()) 
		dictionary.insert(pair<string, string>(key, val));
	else
		iter->second = val;
}

std::vector<std::string> DictionaryConfiguration::getStringArray(const std::string& key) const {
	std::string value_str = getString(key);

	std::vector<std::string> result = splitString(value_str, ',');
    for (size_t i = 0; i < result.size(); ++i) 
		result[i] = trimString(result[i]);

	return result; 
}

/// Namespace manipulation methods
///////////////////////////////////////

bool isNamespaceKey(const string& key, const string& ns, string& nskey)
{
    size_t pos = key.find('.',0);

	if (ns.length() == 0 && pos == string::npos) {
		// in case we requested no namespace (ns is empty) and key has no namespace
		nskey = key;
		return true;
	} else if (ns.length() != 0 && pos != string::npos && key.find(ns,0) == 0) {
		// in case key has namespace and it matches requested namespace ns
		nskey = key.substr(ns.length() + 1);
        return true;
	} else {
        return false;
    }
}

IConfiguration* DictionaryConfiguration::getNamespace(const std::string& namespace_name) const  {
	// split full namespace name into individual namespaces
	std::vector<std::string> namespace_list = splitString(namespace_name, '.'); 

	std::map<std::string, std::string> namespace_dictionary = dictionary; 
	// find key/values for each individual namespace 
	for (auto ns_iter = namespace_list.begin(); ns_iter != namespace_list.end(); ++ns_iter) {

		std::map<std::string, std::string> new_namespace_dictionary; 

		std::string& current_namespace = *ns_iter;
		
		// find all keys that need to be copied 
		// but use namespace_dictionary from previous iteration (in first loop it must be initialized to this->dictionary)
		for (auto key_value_iter = namespace_dictionary.begin(); key_value_iter != namespace_dictionary.end(); ++key_value_iter) {
			std::string namespace_key;

			// add key/value if it has no namespace or if it has current_namespace
			if (isNamespaceKey(key_value_iter->first, "", namespace_key) ||  isNamespaceKey(key_value_iter->first, current_namespace, namespace_key)) {
				// force update if key already exists since higher namespaces have higher priority
				new_namespace_dictionary[namespace_key] = key_value_iter->second;
			}
		}

		// replace resulting dictionary with new dictionary
		namespace_dictionary = new_namespace_dictionary;
	}

	return new DictionaryConfiguration(namespace_dictionary);
}

void DictionaryConfiguration::updateNamespaceReferences() {
	const char* NS_REFERENCE_KEY = "from_namespace";
	size_t reference_key_size = strlen(NS_REFERENCE_KEY);
	// find all keys that are references to other namespaces
	map<string,string> cfg_for_add;
	set<string> cfg_for_delete;
	for (auto it = dictionary.begin(); it != dictionary.end(); it++) {
		size_t pos = it->first.rfind(NS_REFERENCE_KEY);
		if (pos != string::npos && it->first.length() - pos == reference_key_size) {
			// found reference key; 

			// get namespace of current key and namspace this value is pointing to
			const string key_nsname = it->first.substr(0,pos);
			string ref_nsname = it->second;

			// delete this key (so there can be no problems with recursive inclusion)
			//cfg_for_delete.insert(it->first);

			// if target namespace contains :: then we need to look into new file for target namespace values
			bool delete_target = false;
			DictionaryConfiguration* target_cfg = this;
			
			size_t seperator_pos = ref_nsname.find("::");
			if (seperator_pos != string::npos) {
				const string target_file = ref_nsname.substr(0,seperator_pos);
				const string target_nsname = ref_nsname.substr(seperator_pos + 2);

				if (target_file.length() > 0) {
					target_cfg = new FileDictionaryConfiguration(target_file);
					delete_target = true; 
				}
				ref_nsname = target_nsname;
			}			

			// now include all namespace from value of this key
			for (auto sditer = dictionary.begin(); sditer != dictionary.end(); ++sditer) {
				string ref_nskey;

				if (isNamespaceKey(sditer->first, ref_nsname, ref_nskey) && cfg_for_delete.find(ref_nskey) == cfg_for_delete.end()) {
					const string new_key = key_nsname + ref_nskey;
					// now include only if key is not defined (std::map will only insert if key does not exists)
					cfg_for_add.insert(std::pair<string,string>(new_key, sditer->second));
				}
			}
			if (delete_target == true)
				delete target_cfg;
		}
	}
	// erase existing references
	for (set<string>::iterator it = cfg_for_delete.begin(); it != cfg_for_delete.end(); it++) {
		dictionary.erase(*it);
	}
	// add new referenced namespace values	
	for (auto it = cfg_for_add.begin(); it != cfg_for_add.end(); it++) {
		dictionary[it->first] = it->second;
	}
}

/// FileDictionaryConfiguration
///////////////////////////////////////////////////

void FileDictionaryConfiguration::parseFile(const std::string& filename) {
    ifstream is;    

    is.open(filename, ios::in);
    while (is) {
		std::string line;
		string::size_type pos;

        getline(is, line);
        if ((pos = line.find("#")) != string::npos) {
            line.erase(line.begin() + pos, line.end());
        }
        if ((pos = line.find("=")) == string::npos) continue;
		std::string left = trimString(line.substr(0, pos));
        std::string right = trimString(line.substr(pos + 1, line.size() - pos - 1));
        
		dictionary[left] = right;
    }
    is.close();
}


void FileDictionaryConfiguration::saveToFile(const string& filename) {
	ofstream ofs(filename.c_str());
	
	for(auto iter = dictionary.begin(); iter != dictionary.end(); ++iter) {
		ofs << iter->first << " = " << iter->second << endl;
	}

	ofs.close();
}


/// StringDictionaryConfiguration
///////////////////////////////////////////////////

void StringDictionaryConfiguration::parseString(const std::string& config_string, const char delimiter)
{
    size_t pos;
        
    std::string s = trimString(config_string);
    while (!s.empty()) {
		std::string line;
        if ((pos = s.find(delimiter)) == string::npos) {
            line = s;
            s = "";
        } else {
            line = s.substr(0, pos);
            s.erase(0, pos + 1);
        }
		// use # ase comment if delimiter is newline
        if (delimiter == '\n' && (pos = line.find("#")) != string::npos) {
            line.erase(line.begin() + pos, line.end());
        }
        if ((pos = line.find("=")) == string::npos) continue;
        std::string left = trimString(line.substr(0, pos));
        std::string right = trimString(line.substr(pos + 1, line.size() - pos - 1));

        dictionary[left] = right;
    }
}


std::string StringDictionaryConfiguration::toString() {
	ostringstream oss;
	
	for(auto iter = dictionary.begin(); iter != dictionary.end(); ++iter) {
		oss << iter->first << " = " << iter->second << "; ";
	}

	return oss.str();
}
