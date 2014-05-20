// dictionary for configuration propreties

#include "utils/convert.h"
#include "utils/misc.h"
#include "config_dictionary.h"

#include <cstring>

using namespace std;


void split_string(const string& str, char c, vector<string>& v)
{
    string::size_type start = 0, size = str.size(), pos;
    
    v.clear();
    while (start < size && (pos = str.find_first_of(c, start)) != string::npos) {
        v.push_back(str.substr(start, pos - start));
        start = pos + 1;
    }
    if (start < size) v.push_back(str.substr(start, size - start));
}

// ConfigDictionary
//////////////////////

ConfigDictionary::ConfigDictionary() : dictionary() { }

ConfigDictionary::ConfigDictionary(const char* fname) : dictionary() 
{
    parseFile(fname);
}

ConfigDictionary::ConfigDictionary(string fname) : dictionary() 
{
    parseFile(fname.c_str());
}

ConfigDictionary::ConfigDictionary(const ConfigDictionary& cd) :
    dictionary(cd.dictionary.begin(), cd.dictionary.end())
{ }

void ConfigDictionary::fromString(const char* str)
{
    parseString(str);
}

void ConfigDictionary::fromFile(const char* fname)
{
    parseFile(fname);
}

bool isNamespaceKey(string& nskey, const string& key, const string& ns)
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

void ConfigDictionary::fromNamespace(const ConfigDictionary& cfg, const string& nsname)
{
    for (DictionaryMap::const_iterator sditer = cfg.dictionary.begin(); sditer != cfg.dictionary.end(); ++sditer) {
        string nskey;

        if (isNamespaceKey(nskey, sditer->first, nsname)) 
            insert(nskey, sditer->second);            
    }
}

void ConfigDictionary::fromNamespacePriority(const ConfigDictionary& cfg, int ns_size, ...) {
	
	va_list namespace_priority_list;
	va_start(namespace_priority_list, ns_size); 

	ConfigDictionary src_cfg = cfg;	
	for (int i = 0; i < ns_size; i++) {
		ConfigDictionary dst_cfg;
		// first get from root namespace of this one 
		dst_cfg.fromNamespace(src_cfg, "");
		// then get name of this namespace
		char* nsname = va_arg(namespace_priority_list,char*);
		// print warning if not namespace found 
		if (strlen(nsname) > 0 && src_cfg.doesNamespaceExists(nsname) == false)
			cout << "Warning: No namespace '" << nsname << "' found in config." << endl;
		// and copy from values from this namespace
		dst_cfg.fromNamespace(src_cfg, nsname);
		// finaly save this new cfg
		src_cfg = dst_cfg;
		
	}
	va_end(namespace_priority_list);

	dictionary.insert(src_cfg.dictionary.begin(),src_cfg.dictionary.end());
}

void ConfigDictionary::updateNamespaceReferences() {
	const char* NS_REFERENCE_KEY = "from_namespace";
	size_t reference_key_size = strlen(NS_REFERENCE_KEY);
	// find all keys that are references to other namespaces
	ConfigDictionary cfg_for_add;
	set<string> cfg_for_delete;
	for (DictionaryMap::iterator it = dictionary.begin(); it != dictionary.end(); it++) {
		size_t pos = it->first.rfind(NS_REFERENCE_KEY);
		if (pos != string::npos && it->first.length() - pos == reference_key_size) {
			// found reference key; 

			// get namespace of current key and namspace this value is pointing to
			const string key_nsname = it->first.substr(0,pos);
			string ref_nsname = it->second;

			// delete this key (so there can be no problems with recursive inclusion)
			cfg_for_delete.insert(it->first);


			// if target namespace contains :: then we need to look into new file for target namespace values
			bool delete_target = false;
			ConfigDictionary* target_cfg = this;
			
			size_t seperator_pos = ref_nsname.find("::");
			if (seperator_pos != string::npos) {
				const string target_file = ref_nsname.substr(0,seperator_pos);
				const string target_nsname = ref_nsname.substr(seperator_pos + 2);

				if (target_file.length() > 0) {
					target_cfg = new ConfigDictionary(target_file);
					delete_target = true; 
				}
				ref_nsname = target_nsname;
			}			

			// now include all namespace from value of this key
			for (DictionaryMap::const_iterator sditer = dictionary.begin(); sditer != dictionary.end(); ++sditer) {
				string ref_nskey;

				if (isNamespaceKey(ref_nskey, sditer->first, ref_nsname) && cfg_for_delete.find(ref_nskey) == cfg_for_delete.end()) {
					const string new_key = key_nsname + ref_nskey;
					// now include only if key is not defined
					if (isDefined(new_key) == false)
						cfg_for_add.insert(new_key, sditer->second);
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
	for (DictionaryMap::iterator it = cfg_for_add.dictionary.begin(); it != cfg_for_add.dictionary.end(); it++) {
		insert(it->first, it->second);
	}
}
bool ConfigDictionary::doesNamespaceExists(const string& nsname) {
	for (DictionaryMap::iterator it = dictionary.begin(); it != dictionary.end(); it++) {
		string nskey;
		if (isNamespaceKey(nskey, it->first, nsname))
			return true;
	}
	return false;
}

bool ConfigDictionary::toFile(const string& filename)
{
	try
	{
		ofstream ofs(filename.c_str());
		for(DictionaryMap::iterator i = dictionary.begin(); i != dictionary.end(); ++i)
		{
			ofs << i->first << " = " << i->second << endl;
		}
		ofs.close();
		return true;
	}
	catch(...)
	{
		return false;
	}
}
/*
void ConfigDictionary::replaceKeys(const char* fname)
{
    dictionary.clear();
    fromFile(fname);
}*/

bool ConfigDictionary::getValue(int& i, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(ConfigException, "Key '" + key + "' not specified");
        return false;
    }
    try {
        i = Convert::ToInteger<int>(iter->second.c_str(), (unsigned)iter->second.size());
        return true;
    } catch (InvalidConversionException) {
        if (thr) throw custom_libhop_exception(ConfigException, "Invalid value for key '" + key + "'");
        return false;
    }
}

bool ConfigDictionary::getValue(vector<int>& v, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(ConfigException, "Key '" + key + "' not specified");
        return false;
    }
    try {
        vector<string> strings;

        v.clear();
        split_string(iter->second, ',', strings);
        for (size_t i = 0; i < strings.size(); ++i) {
            int res = Convert::ToInteger<int>(strings[i].c_str(), (unsigned)strings[i].size());
            v.push_back(res);
        }
        return true;
    } catch (InvalidConversionException) {
        if (thr) throw custom_libhop_exception(ConfigException, "Invalid value for key '" + key + "'");
        return false;
    }    
}

bool ConfigDictionary::getValue(vector<string>& v, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(ConfigException, "Key '" + key + "' not specified");
        return false;
    }
    v.clear();
    split_string(iter->second, ',', v);
    for (size_t i = 0; i < v.size(); ++i) v[i] = trimString(v[i]);
    return true;
}

void ConfigDictionary::setValue(const string& key, const string& val)
{
	map<string,string>::iterator iter;

	if ((iter = dictionary.find(key)) == dictionary.end()) 
		dictionary.insert(pair<string, string>(key, val));
	else
		iter->second = val;
}

void ConfigDictionary::setValue(const string& key, int val)
{
    ostringstream result;
    result << val;
    setValue(key, result.str());
}

void ConfigDictionary::setValue(const string& key, double val)
{
    ostringstream result;
    result << val;
    setValue(key, result.str());
}

bool ConfigDictionary::getValue(vector<double>& v, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(ConfigException, "Key '" + key + "' not specified");
        return false;
    }
    try {
        vector<string> strings;

        v.clear();
        split_string(iter->second, ',', strings);
        for (size_t i = 0; i < strings.size(); ++i) {
            double res = Convert::ToDouble(strings[i].c_str(), (unsigned)strings[i].size());
            v.push_back(res);
        }
        return true;
    } catch (InvalidConversionException) {
        if (thr) throw custom_libhop_exception(ConfigException, "Invalid value for key '" + key + "'");
        return false;
    }    
}

bool ConfigDictionary::getValue(vector<bool>& v, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(ConfigException, "Key '" + key + "' not specified");
        return false;
    }
    try {
        vector<string> strings;

        v.clear();
        split_string(iter->second, ',', strings);
        for (size_t i = 0; i < strings.size(); ++i) {
            v.push_back(strings[i] == "true");
        }
        return true;
    } catch (InvalidConversionException) {
        if (thr) throw custom_libhop_exception(ConfigException, "Invalid value for key '" + key + "'");
        return false;
    }    
}

bool ConfigDictionary::getValue(double& d, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(ConfigException, "Key '" + key + "' not specified");
        return false;
    }
    try {
        d = Convert::ToDouble(iter->second.c_str(), (unsigned)iter->second.size());
        return true;
    } catch (InvalidConversionException) {
        if (thr) throw custom_libhop_exception(ConfigException, "Invalid value for key '" + key + "'");
        return false;
    }
}

bool ConfigDictionary::getValue(string& s, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(ConfigException, "Key '" + key + "' not specified");
        return false;
    }
    s = iter->second;
    return true;
}

bool ConfigDictionary::getValue(bool& b, const string& key, bool thr) const
{
	string s;
	bool result = getValue(s, key, thr);

	if (!result) return false;
	b = (s == "true");
    return true;
}


double ConfigDictionary::getValueDouble(const string& key, double def, bool thr) const
{
    double d;

    if (getValue(d, key, thr)) return d; else return def;
}

int ConfigDictionary::getValueInt(const string& key, int def, bool thr) const
{
    int i;

    if (getValue(i, key, thr)) return i; else return def;
}

string ConfigDictionary::getValueString(const string& key, const string& def, bool thr) const
{
    string s;

    if (getValue(s, key, thr)) return s; else return def;
}

bool ConfigDictionary::getValueBool(const string& key, bool def, bool thr) const
{
    string s;

    if (getValue(s, key, thr)) return s == "true"; else return def;
}


bool ConfigDictionary::isDefined(const string& key) const
{
    return (dictionary.find(key) != dictionary.end());
}


void ConfigDictionary::parseString(string s, const char delimiter)
{
    size_t pos;
    string line, left, right;
        
    s = trimString(s);
    while (!s.empty()) {
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
        left = line.substr(0, pos);
        right = line.substr(pos + 1, line.size() - pos - 1);
        left = trimString(left); 
		right = trimString(right);

        insert(left, right);
    }
}

void ConfigDictionary::parseFile(const char* fname)
{
    ifstream is;
    string s, left, right;
    string::size_type pos;

    is.open(fname, ios::in);
    while (is) {
        getline(is, s);
        if ((pos = s.find("#")) != string::npos) {
            s.erase(s.begin() + pos, s.end());
        }
        if ((pos = s.find("=")) == string::npos) continue;
        left = s.substr(0, pos);
        right = s.substr(pos + 1, s.size() - pos - 1);
        left = trimString(left); 
		right = trimString(right);
        
        insert(left, right);
    }
    is.close();
}

void ConfigDictionary::insert(const string& key, const string& val)
{
    pair<DictionaryMap::iterator, bool> ib =
        dictionary.insert(pair<string, string>(key, val));

    if (!ib.second) ib.first->second = val;
}
