// miscellanous classes and functions

#include <string.h>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <map>
#include <sstream>
#if defined WIN32 | defined WIN64
#include <windows.h>
#endif
#include "platform.h"
#include "convert.h"
#include "misc.h"



using namespace std;

// functions
///////////////////////////////////////////////////////////////////////////////

//write float to matlab matrix
void data2mat(char* filename, float* ptr, const int dimx, const int dimy)
{
	ofstream ofs;
	ofs.open(filename);
	for(int y = 0; y < dimy; ++y)
	{
		for(int x = 0; x < dimx; ++x)
		{
			if(x != 0) ofs << ',';
			ofs << *ptr;
			ptr++;
		}
		ofs << endl;
	}
	ofs.close();
}

void data2mat(char* filename, double* ptr, const int dimx, const int dimy)
{
	ofstream ofs;
	ofs.open(filename);
	for(int y = 0; y < dimy; ++y)
	{
		for(int x = 0; x < dimx; ++x)
		{
			if(x != 0) ofs << ',';
			ofs << *ptr;
			ptr++;
		}
		ofs << endl;
	}
	ofs.close();
}

void trim(string& str)
{
    string::size_type pos = 0, len = str.length();

    while (pos < len && str[pos] <= ' ') ++pos;
    str.erase(0, pos);
    pos = str.length();
    while (pos > 0 && str[pos - 1] <= ' ') --pos;
    str.erase(pos, str.length() - pos);
}

    //string::size_type pos = str.find_last_not_of(' ');

    //if (pos != string::npos) {
    //    str.erase(pos + 1);
    //    pos = str.find_first_not_of(' ');
    //    if (pos != string::npos) str.erase(0, pos);
    //}
    //else str.erase(str.begin(), str.end());


void end_dir(string& dir_name)
{
    if (dir_name.size() > 0 && dir_name[dir_name.size() - 1] != FILE_SEPARATOR) dir_name += FILE_SEPARATOR;
}

string change_extension(const string& fname, const string& ext, const string& extsep)
{
    string result;
    string::size_type dot_pos;

    if ((dot_pos = fname.rfind(extsep)) == string::npos) result = fname + ext;
    else result = fname.substr(0, dot_pos) + ext;
    return result;
}

string get_extension(const string& fname, const string& extsep)
{
    string::size_type dot_pos;

    if ((dot_pos = fname.rfind(extsep)) == string::npos) return "";
    else return fname.substr(dot_pos);
}


string change_extension(const string& fname, const string& ext)
{
    return change_extension(fname, ext, ".");
}

int scale_number(const string& fname)
{
	string s = change_extension(fname, "");
	if (s.empty()) return -1; else return (int)(s.back() - '0');
}

int hash_function(const string& s)
{
    int result = 0;

    for (string::const_iterator iter = s.begin(); iter != s.end(); ++iter) {
        char c = *iter;
        result += (int)(c - ' ');
    }
    return result;
}

string get_prefix(string str, const string& sep, int n /* = 1*/)
{
    string::size_type sep_pos = str.size(), off = 0;

    while (n-- > 0) {
        if ((sep_pos = str.find(sep, off)) == string::npos) return string();
        else off = sep_pos + 1;
    }
    return str.substr(0, sep_pos);
}

string operator+(const string& s, int i)
{
    ostringstream result;
    result << s << i;
    return result.str();
}

string fill_left(int d, int width, char fillchar /* = '0' */)
{
    stringstream s;

    s.width(width);
    s.fill(fillchar);
    s << d;
    return s.str();
}

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

void HSLtoRGB(const double h, const double s, const double l, double* result)
{
    double q = (l < 0.5) ? (l*(1 + s)) : (l + s - l*s);
    double p = 2*l - q;
    double t[3] = { h + 1.0/3.0, h, h - 1.0/3.0 };
    double res;
    int C;
    
    for (C = 0; C < 3; ++C) {
        if (t[C] < 0.0) t[C] += 1.0;
        if (t[C] > 1.0) t[C] -= 1.0;
    }
    for (C = 0; C < 3; ++C) {
        if (t[C] < 1.0/6.0) res = p + ((q - p)*6.0*t[C]);
        else if (t[C] < 0.5) res = q;
        else if (t[C] < 2.0/3.0) res = p + ((q - p)*6.0*(2.0/3.0 - t[C]));
        else res = p;
        result[C] = res;
    }
}

double round(double d, int places)
{
    int p10 = 1;

    for (int i = 0; i < places; ++i) p10 *= 10;
    return (d < 0.0) ? ceil(d * p10 - 0.5)/p10 : floor(d * p10 + 0.5)/p10;
} 

set<int> set_range(int len)
{
    set<int> result;

    for (int i = 0; i < len; ++i) result.insert(i);
    return result;
}

set<int> set_range(int min, int max)
{
    set<int> result;

    while (min <= max) result.insert(min++); 
    return result;
}

vector<int> vector_range(int min, int max)
{
    vector<int> result;

    while (min <= max) result.push_back(min++); 
    return result;
}

// method definitions
///////////////////////////////////////////////////////////////////////////////

// memory allocation
//////////////////////

// mem_allocator
//////////////////

mem_allocator::mem_allocator(size_t ssize, size_t scount /* = 1024 */, const char* n) :
    name(n),
    mem(0), 
    slot_size(max(ssize + sizeof(char*), sizeof(info_block))), 
    slot_count(scount), 
    mem_count(1)
{ 
    mem = (char**)malloc(sizeof(char*) * mem_count);
    init_mem_range(mem[0]);
    free_slot = mem[0] + slot_size;
}

mem_allocator::~mem_allocator()
{
    for (size_t i = 0; i < mem_count; ++i)  free(mem[i]);
    free(mem);
}

void mem_allocator::init_mem_range(char*& memp)
{
    memp = (char*)malloc(slot_size * slot_count);
    char* p = memp + slot_size;
    info_block* info;
    info_block0* info0 = (info_block0*)memp;

    info0->size = slot_size - sizeof(char*);
    info0->bit_size = ilog2_x86(info0->size);
    for (size_t i = 2; i < slot_count; ++i) {
        info = (info_block*)p;
        info->mem0 = memp;
        info->mem0_next = memp;
        info->next = i;
        p += slot_size;
    }
    info = (info_block*)p;
    info->mem0 = memp;
    info->mem0_next = memp;
    info->next = UINT_MAX;
}

void mem_allocator::expand_mem()
{
    char** newmem = (char**)malloc(sizeof(char*) * (mem_count + 1));

    memcpy(newmem, mem, sizeof(char*) * mem_count);
    free(mem);
    mem = newmem;
    init_mem_range(mem[mem_count]);
    free_slot = mem[mem_count] + slot_size;
    ++mem_count;
    //cout << "expanded: size " << slot_size << " with " << slot_count << " slots in " << mem_count << " mems" << endl;
}

// mem_allocator_set
//////////////////////

const char* pow_strings[] = {"1", "2", "4", "8", "16", "32", "64", "128", "256",
    "512", "1024", "2048"};
const int table_sizes[] = {0, 0, 0, 64, 512, 1024, 512, 512, 128, 64};



mem_allocator_set::mem_allocator_set(size_t pcount /* = 10 */) :
    count(max((unsigned int)pcount, 10U))
{
    allocators = new mem_allocator*[count];
    memset(allocators, 0, sizeof(mem_allocator*) * count);
    for (size_t i = 4; i < count; ++i) {
        allocators[i] = new mem_allocator(1 << i, table_sizes[i], pow_strings[i]);
    }

}

mem_allocator_set::~mem_allocator_set()
{
    for (size_t i = 0; i < count; ++i) 
        if (allocators[i] != nullptr) delete allocators[i];
    delete allocators;
}

// allocator_base (static member definition)
//////////////////////////////////////////////

mem_allocator_set allocator_base::allocators;

// int_map
////////////

void int_map::reset(int domain, int codomain) 
{ 
    f.resize(domain);
    finv.resize(codomain);
    make_identity();
}


void int_map::reset(int domain_range, const vector<int>& v)
{
    f.resize(domain_range);
    fill(f.begin(), f.end(), -1);
    finv.resize(v.size());
    for (int i = 0; i < (int)v.size(); ++i) {
        if (v[i] >= 0 && v[i] < domain_range) {
            f[v[i]] = i;
            finv[i] = v[i];
        }
    }
}

// config_dictionary
//////////////////////

config_dictionary::config_dictionary() : dictionary() { }

config_dictionary::config_dictionary(const char* fname) : dictionary() 
{
    parse_file(fname);
}

config_dictionary::config_dictionary(string fname) : dictionary() 
{
    parse_file(fname.c_str());
}

config_dictionary::config_dictionary(const config_dictionary& cd) :
    dictionary(cd.dictionary.begin(), cd.dictionary.end())
{ }

void config_dictionary::from_string(const char* str)
{
    parse_string(str);
}

void config_dictionary::from_string_nl(const char* str)
{
    parse_string_nl(str);
}

void config_dictionary::from_file(const char* fname)
{
    parse_file(fname);
}

bool is_namespace_key(string& nskey, const string& key, const string& ns)
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

void config_dictionary::from_namespace(const config_dictionary& cfg, const string& nsname)
{
    for (dictionary_t::const_iterator sditer = cfg.dictionary.begin(); sditer != cfg.dictionary.end(); ++sditer) {
        string nskey;

        if (is_namespace_key(nskey, sditer->first, nsname)) 
            insert(nskey, sditer->second);            
    }
}

void config_dictionary::from_namespace_priority(const config_dictionary& cfg, int ns_size, ...) {
	/*// include keys by reverse order (first namespace has bigger priority, last one has lower)
	for (vector<string>::const_reverse_iterator r_iter = nsname_priority.rbegin(); r_iter != nsname_priority.rend(); r_iter++) {
		from_namespace(cfg, *r_iter);
	}
	// finally no namespace indicates biggest priority so it should override all others
	from_namespace(cfg, "");*/
	va_list namespace_priority_list;
	va_start(namespace_priority_list, ns_size); 

	config_dictionary src_cfg = cfg;	
	for (int i = 0; i < ns_size; i++) {
		config_dictionary dst_cfg;
		// first get from root namespace of this one 
		dst_cfg.from_namespace(src_cfg, "");
		// then get name of this namespace
		char* nsname = va_arg(namespace_priority_list,char*);
		// print warning if not namespace found 
		if (strlen(nsname) > 0 && src_cfg.does_namespace_exist(nsname) == false)
			cout << "Warning: No namespace '" << nsname << "' found in config." << endl;
		// and copy from values from this namespace
		dst_cfg.from_namespace(src_cfg, nsname);
		// finaly save this new cfg
		src_cfg = dst_cfg;
		
	}
	va_end(namespace_priority_list);

	dictionary.insert(src_cfg.dictionary.begin(),src_cfg.dictionary.end());
}

void config_dictionary::update_namespace_references() {
	const char* NS_REFERENCE_KEY = "from_namespace";
	int reference_key_size = strlen(NS_REFERENCE_KEY);
	// find all keys that are references to other namespaces
	config_dictionary cfg_for_add;
	set<string> cfg_for_delete;
	for (dictionary_t::iterator it = dictionary.begin(); it != dictionary.end(); it++) {
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
			config_dictionary* target_cfg = this;
			
			size_t seperator_pos = ref_nsname.find("::");
			if (seperator_pos != string::npos) {
				const string target_file = ref_nsname.substr(0,seperator_pos);
				const string target_nsname = ref_nsname.substr(seperator_pos + 2);

				if (target_file.length() > 0) {
					target_cfg = new config_dictionary(target_file);
					delete_target = true; 
				}
				ref_nsname = target_nsname;
			}			

			// now include all namespace from value of this key
			for (dictionary_t::const_iterator sditer = dictionary.begin(); sditer != dictionary.end(); ++sditer) {
				string ref_nskey;

				if (is_namespace_key(ref_nskey, sditer->first, ref_nsname) && cfg_for_delete.find(ref_nskey) == cfg_for_delete.end()) {
					const string new_key = key_nsname + ref_nskey;
					// now include only if key is not defined
					if (is_defined(new_key) == false)
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
	for (dictionary_t::iterator it = cfg_for_add.dictionary.begin(); it != cfg_for_add.dictionary.end(); it++) {
		insert(it->first, it->second);
	}
}
bool config_dictionary::does_namespace_exist(const string& nsname) {
	for (dictionary_t::iterator it = dictionary.begin(); it != dictionary.end(); it++) {
		string nskey;
		if (is_namespace_key(nskey, it->first, nsname))
			return true;
	}
	return false;
}

bool config_dictionary::to_file(const string& filename)
{
	try
	{
		ofstream ofs(filename.c_str());
		for(dictionary_t::iterator i = dictionary.begin(); i != dictionary.end(); ++i)
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

void config_dictionary::replace_keys(const char* fname)
{
    dictionary.clear();
    parse_file(fname);
}

bool config_dictionary::get_value(int& i, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(config_exception, "Key '" + key + "' not specified");
        return false;
    }
    try {
        i = Convert::ToInteger<int>(iter->second.c_str(), (unsigned)iter->second.size());
        return true;
    } catch (InvalidConversionException) {
        if (thr) throw custom_libhop_exception(config_exception, "Invalid value for key '" + key + "'");
        return false;
    }
}

bool config_dictionary::get_value(vector<int>& v, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(config_exception, "Key '" + key + "' not specified");
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
        if (thr) throw custom_libhop_exception(config_exception, "Invalid value for key '" + key + "'");
        return false;
    }    
}

bool config_dictionary::get_value(vector<string>& v, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(config_exception, "Key '" + key + "' not specified");
        return false;
    }
    v.clear();
    split_string(iter->second, ',', v);
    for (size_t i = 0; i < v.size(); ++i) trim(v[i]);
    return true;
}

void config_dictionary::set_value(const string& key, const string& val)
{
	map<string,string>::iterator iter;

	if ((iter = dictionary.find(key)) == dictionary.end()) 
		dictionary.insert(pair<string, string>(key, val));
	else
		iter->second = val;
}

void config_dictionary::set_value(const string& key, int val)
{
    ostringstream result;
    result << val;
    set_value(key, result.str());
}

void config_dictionary::set_value(const string& key, double val)
{
    ostringstream result;
    result << val;
    set_value(key, result.str());
}

bool config_dictionary::get_value(vector<double>& v, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(config_exception, "Key '" + key + "' not specified");
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
        if (thr) throw custom_libhop_exception(config_exception, "Invalid value for key '" + key + "'");
        return false;
    }    
}

bool config_dictionary::get_value(vector<bool>& v, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(config_exception, "Key '" + key + "' not specified");
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
        if (thr) throw custom_libhop_exception(config_exception, "Invalid value for key '" + key + "'");
        return false;
    }    
}

bool config_dictionary::get_value(double& d, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(config_exception, "Key '" + key + "' not specified");
        return false;
    }
    try {
        d = Convert::ToDouble(iter->second.c_str(), (unsigned)iter->second.size());
        return true;
    } catch (InvalidConversionException) {
        if (thr) throw custom_libhop_exception(config_exception, "Invalid value for key '" + key + "'");
        return false;
    }
}

bool config_dictionary::get_value(string& s, const string& key, bool thr) const
{
    map<string,string>::const_iterator iter;

    if ((iter = dictionary.find(key)) == dictionary.end()) {
        if (thr) throw custom_libhop_exception(config_exception, "Key '" + key + "' not specified");
        return false;
    }
    s = iter->second;
    return true;
}

bool config_dictionary::get_value(bool& b, const string& key, bool thr) const
{
	string s;
	bool result = get_value(s, key, thr);

	if (!result) return false;
	b = (s == "true");
    return true;
}


double config_dictionary::get_value_double(const string& key, double def, bool thr) const
{
    double d;

    if (get_value(d, key, thr)) return d; else return def;
}

int config_dictionary::get_value_int(const string& key, int def, bool thr) const
{
    int i;

    if (get_value(i, key, thr)) return i; else return def;
}

string config_dictionary::get_value_string(const string& key, const string& def, bool thr) const
{
    string s;

    if (get_value(s, key, thr)) return s; else return def;
}

bool config_dictionary::get_value_bool(const string& key, bool def, bool thr) const
{
    string s;

    if (get_value(s, key, thr)) return s == "true"; else return def;
}


bool config_dictionary::is_defined(const string& key) const
{
    return (dictionary.find(key) != dictionary.end());
}

void config_dictionary::parse_string(string s)
{
    size_t pos;
    string line, left, right;
        
    trim(s);
    while (!s.empty()) {
        if ((pos = s.find(';')) == string::npos) {
            line = s;
            s = "";
        } else {
            line = s.substr(0, pos);
            s.erase(0, pos + 1);
        }
        if ((pos = line.find("=")) == string::npos) continue;
        left = line.substr(0, pos);
        right = line.substr(pos + 1, line.size() - pos - 1);
        trim(left); trim(right);

        insert(left, right);
    }
}

void config_dictionary::parse_string_nl(string s)
{
    size_t pos;
    string line, left, right;
        
    trim(s);
    while (!s.empty()) {
        if ((pos = s.find('\n')) == string::npos) {
            line = s;
            s = "";
        } else {
            line = s.substr(0, pos);
            s.erase(0, pos + 1);
        }
        if ((pos = line.find("#")) != string::npos) {
            line.erase(line.begin() + pos, line.end());
        }
        if ((pos = line.find("=")) == string::npos) continue;
        left = line.substr(0, pos);
        right = line.substr(pos + 1, line.size() - pos - 1);
        trim(left); trim(right);

        insert(left, right);
    }
}

void config_dictionary::parse_file(const char* fname)
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
        trim(left); trim(right);
        
        insert(left, right);
    }
    is.close();
}

void config_dictionary::insert(const string& key, const string& val)
{
    pair<dictionary_t::iterator, bool> ib =
        dictionary.insert(pair<string, string>(key, val));

    if (!ib.second) ib.first->second = val;
}

// matrix operations
///////////////////////////////////////////////////////////////////////////////

// Eigenvalues and eigenvectors of a symmetric 2 x 2 matrix A
// A = [ a  b ]
//     [ b  d ]
void eigensystem(double& l1, double& l2, dpoint2& v1, dpoint2& v2, double a, double b, double d)
{
    if (b == 0) {
        l1 = a; l2 = b;
        v1.set(1.0, 0.0);
        v2.set(0.0, 1.0);
    }

    double x, D, v;
    
    x = a - d;
    D = ::sqrt(4 * b * b + x * x);
    l1 = (a + d - D)/2;
    l2 = (a + d + D)/2;
    
    v = (a - d - D)/2.0/b;
    v1.set(v, 1.0);
    v1.div(::sqrt(v1.norm2()));

    v = (a - d + D)/2.0/b;
    v2.set(v, 1.0);
    v2.div(::sqrt(v2.norm2()));

}


//unicode string conversion
/*
std::string wstrtostr(const std::wstring &wstr)
{
    // Convert a Unicode string to an ASCII string
    std::string strTo;
    char *szTo = new char[wstr.length() + 1];
    szTo[wstr.size()] = '\0';
    WideCharToMultiByte(CP_ACP, 0, wstr.c_str(), -1, szTo, (int)wstr.length(), nullptr, nullptr);
    strTo = szTo;
    delete[] szTo;
    return strTo;
}

std::wstring strtowstr(const std::string &str)
{
    // Convert an ASCII string to a Unicode String
    std::wstring wstrTo;
    wchar_t *wszTo = new wchar_t[str.length() + 1];
    wszTo[str.size()] = L'\0';
    MultiByteToWideChar(CP_ACP, 0, str.c_str(), -1, wszTo, (int)str.length());
    wstrTo = wszTo;
    delete[] wszTo;
    return wstrTo;
}

void replace_substr(string& orig, const string& substr_old, const string& substr_new)
{
	int i = 0;
	while(i != string::npos || i <= orig.size())
	{
		i = orig.find(substr_old, i);
		if(i == string::npos) return;
		orig.replace(i, substr_old.size(), substr_new);
		i+= substr_new.size();
	}
}

//escape all perl regex special char
void escape_for_perlregex(string& se)
{
	replace_substr(se, "\\", "\\\\");	//MUST BE 1st
	replace_substr(se, ".", "\\.");
	replace_substr(se, "[", "\\[");
	replace_substr(se, "{", "\\{");
	replace_substr(se, "(", "\\(");
	replace_substr(se, ")", "\\)");
	replace_substr(se, "*", "\\*");
	replace_substr(se, "+", "\\+");
	replace_substr(se, "?", "\\?");
	replace_substr(se, "|", "\\|");
	replace_substr(se, "^", "\\^");
	replace_substr(se, "$", "\\$");	
}
*/

