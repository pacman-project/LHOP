// config_dictionary.h

#pragma once
#ifndef _config_dictionary_H_
#define _config_dictionary_H_

#include <string>
#include <map>

#include "utils/serialization/streaming.h"
#include "utils/exceptions.h"

using namespace std;

// ConfigException
/////////////////////

class ConfigException : public libhop_exception {
public:
	ConfigException(const string& file, const size_t line, const string& msg = "Unknown ConfigException"): libhop_exception(file, line, msg, "ConfigException") {}
};


#include <stdarg.h>

// ConfigDictionary 
//////////////////////
class ConfigDictionary {	
protected:
    typedef map<string, string> DictionaryMap;

    DictionaryMap dictionary;
public:
    ConfigDictionary();
    ConfigDictionary(const char*);
    ConfigDictionary(string);
    ConfigDictionary(const ConfigDictionary&);

    //void replaceKeys(const char*);
    void fromFile(const char*);
	bool toFile(const string& filename);
    void fromString(const char*);
    void fromNamespace(const ConfigDictionary& cfg, const string& nsname);
	void fromNamespacePriority(const ConfigDictionary& cfg, int ns_size, ...);
	void updateNamespaceReferences();
	bool doesNamespaceExists(const string& nsname);
    bool getValue(double&, const string&, bool = false) const;
	bool getValue(bool&, const string&, bool = false) const;
    bool getValue(int&, const string&, bool = false) const;
    bool getValue(string&, const string&, bool = false) const;
    bool getValue(vector<int>&, const string&, bool = false) const;
    bool getValue(vector<double>&, const string&, bool = false) const;
    bool getValue(vector<bool>&, const string&, bool = false) const;
	bool getValue(vector<string>& v, const string& key, bool = false) const;
    template <typename T> void getValue(vector_parameter<T>&, const string&, const T&) const;
    double getValueDouble(const string&, double, bool = false) const;
    int getValueInt(const string&, int, bool = false) const;
    string getValueString(const string&, const string&, bool = false) const;
    bool getValueBool(const string&, bool, bool = false) const;
    bool isDefined(const string&) const;
	void setValue(const string& key, const string& val);
    void setValue(const string& key, int val);
    void setValue(const string& key, double val);
	
protected:
    void parseFile(const char*);
    void parseString(string s, const char delimiter = ';');
    void insert(const string& key, const string& val);
};

template <typename T> void ConfigDictionary::getValue(vector_parameter<T>& par, const string& key, 
        const T& def) const
{
    vector<T> v;

    getValue(v, key, false);
    if (v.empty()) v.push_back(def);
    par.set_val(v);

}
#endif /* _config_dictionary_H_*/