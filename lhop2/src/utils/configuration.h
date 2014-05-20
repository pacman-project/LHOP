// config_dictionary.h

#pragma once
#ifndef _UTILS_CONFIGURATION_H_
#define _UTILS_CONFIGURATION_H_

#include <string>
#include <map>

#include "utils/serialization/streaming.h"
#include "utils/exceptions.h"

using namespace std;

// ConfigException
/////////////////////

class ConfigException_ : public libhop_exception {
public:
	ConfigException_(const string& file, const size_t line, const string& msg = "Unknown ConfigException"): libhop_exception(file, line, msg, "ConfigException") {}
};


class IConfiguration {
public:
	////////////////////////////////////////////
	/// Abstract methods required by interface	
	virtual bool hasKey(const std::string& key) const = 0;

	virtual double getDouble(const std::string& key) const = 0;
	virtual bool getBool(const std::string& key) const = 0;
	virtual int getInt(const std::string& key) const = 0;
	virtual std::string getString(const std::string& key) const = 0;

	virtual std::vector<double> getDoubleArray(const std::string& key) const = 0;
	virtual std::vector<bool> getBoolArray(const std::string& key) const = 0;
	virtual std::vector<int> getIntArray(const std::string& key) const = 0;
	virtual std::vector<std::string> getStringArray(const std::string& key) const = 0;

	virtual void setDouble(const std::string& key, double val) = 0;
	virtual void setBool(const std::string& key, bool val) = 0;
	virtual void setInt(const std::string& key, int val) = 0;
	virtual void setString(const std::string& key, const std::string& val) = 0;

	virtual IConfiguration* getNamespace(const std::string& namespace_name) const = 0;

	////////////////////////////////////////////
	/// Additional overloaded getters with default values	
	virtual double getDouble(const std::string& key, double default_value) const {
		return hasKey(key) ? getDouble(key) : default_value;
	}
	virtual bool getBool(const std::string& key, bool default_value) const {
		return hasKey(key) ? getBool(key) : default_value;
	}
	virtual int getInt(const std::string& key, int default_value) const {
		return hasKey(key) ? getInt(key) : default_value;
	}
	virtual std::string getString(const std::string& key, std::string default_value) const {
		return hasKey(key) ? getString(key) : default_value;
	}
	virtual std::vector<double> getDoubleArray(const std::string& key, std::vector<double> default_value) const {
		return hasKey(key) ? getDoubleArray(key) : default_value;
	}
	virtual std::vector<bool> getBoolArray(const std::string& key, std::vector<bool> default_value) const {
		return hasKey(key) ? getBoolArray(key) : default_value;
	}
	virtual std::vector<int> getIntArray(const std::string& key, std::vector<int> default_value) const {
		return hasKey(key) ? getIntArray(key) : default_value;
	}
	virtual std::vector<std::string> getStringArray(const std::string& key, std::vector<std::string> default_value) const {
		return hasKey(key) ? getStringArray(key) : default_value;
	}
};

#include <stdarg.h>

class DictionaryConfiguration : public IConfiguration {
protected:
    std::map<std::string, std::string> dictionary;
public:
	DictionaryConfiguration() {}
	DictionaryConfiguration(const std::map<std::string, std::string>& key_values) : dictionary(key_values) {
		updateNamespaceReferences();
	}
	DictionaryConfiguration(const DictionaryConfiguration& config) : dictionary(config.dictionary) {
	}
	
	/**
	 * Updates current dictionary with key/values from config but only if key does not exist in current dictionary.
	 * Will also update namespace refernces if set to true (by default is true).
	 */ 
	void appendConfig(const DictionaryConfiguration& config, bool update_namespace_references = true) {
		dictionary.insert(config.dictionary.begin(), config.dictionary.end());

		if (update_namespace_references)
			updateNamespaceReferences();
	}

	/////////////////////////////////////////////////////////////////////
	/// General setters/getters required by IConfiguration interface	 
	virtual bool hasKey(const std::string& key) const;

	virtual double getDouble(const std::string& key) const;
	virtual bool getBool(const std::string& key) const;
	virtual int getInt(const std::string& key) const;
	virtual std::string getString(const std::string& key) const;

	virtual std::vector<double> getDoubleArray(const std::string& key) const;
	virtual std::vector<bool> getBoolArray(const std::string& key) const;
	virtual std::vector<int> getIntArray(const std::string& key) const;
	virtual std::vector<std::string> getStringArray(const std::string& key) const;

	virtual void setDouble(const std::string& key, double val);
	virtual void setBool(const std::string& key, bool val);
	virtual void setInt(const std::string& key, int val);
	virtual void setString(const std::string& key, const std::string& val);

	/////////////////////////////////////////////////////////////////////
	/// Namespace methods required by IConfiguration interface	
	virtual IConfiguration* getNamespace(const std::string& namespace_name) const ;

protected:
	void updateNamespaceReferences();
};


/**
 * FileDictionaryConfiguration
 *
 * Loads dictionary configurations from a file. One line in file corresponds to one pair of key/value
 * in following format: "key = value # optional coment". Char '#' represent a comment and everything 
 * from a '#' to the end of line is ignored.
 */
class FileDictionaryConfiguration : public DictionaryConfiguration {
public:
	FileDictionaryConfiguration(const std::string& filename) {
		parseFile(filename);
		updateNamespaceReferences();
	}
	FileDictionaryConfiguration(const DictionaryConfiguration& config) : DictionaryConfiguration(config) {}

	void saveToFile(const string& filename);
private:	
	void parseFile(const std::string& filename);
};


/**
 * StringDictionaryConfiguration
 *
 * Loads dictionary configurations from a string with ';' as default delimiter between values.
 * Example of config_string: 'start_layer = 1; end_layer = 4; ly1.inference.init_size = 200;'
 */
class StringDictionaryConfiguration : public DictionaryConfiguration {
public:
	StringDictionaryConfiguration(const std::string& config_string, const char delimiter = ';') {
		parseString(config_string, delimiter);
		updateNamespaceReferences();
	}
	StringDictionaryConfiguration(const DictionaryConfiguration& config) : DictionaryConfiguration(config) {}

	std::string toString();
private:
    void parseString(const std::string& config_string, const char delimiter = ';');
};

#endif /* _UTILS_CONFIGURATION_H_*/