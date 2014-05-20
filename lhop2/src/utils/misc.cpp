// miscellanous classes and functions

#include <fstream>
#include <vector>

#if defined WIN32 | defined WIN64
#include <windows.h>
#endif

#include "utils/platform.h"
#include "misc.h"



using namespace std;

// functions
///////////////////////////////////////////////////////////////////////////////
std::string trimString(const std::string& str)
{
	std::string result = str;

    string::size_type pos = 0, len = result.length();

    while (pos < len && result[pos] <= ' ') ++pos;
    result.erase(0, pos);
    pos = result.length();
    while (pos > 0 && result[pos - 1] <= ' ') --pos;
    result.erase(pos, result.length() - pos);

	return result;
}

string fixPathEndSeparator(const string& dir_name) {
    if (dir_name.size() > 0 && dir_name[dir_name.size() - 1] != FILE_SEPARATOR) 
		return dir_name + FILE_SEPARATOR;
	else
		return dir_name;
}

// TODO: replace with fixEndDirSeparator
void end_dir(string& dir_name) {
	dir_name = fixPathEndSeparator(dir_name);
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

string operator+(const string& s, int i)
{
    ostringstream result;
    result << s << i;
    return result.str();
}

double round(double d, int places)
{
    int p10 = 1;

    for (int i = 0; i < places; ++i) p10 *= 10;
    return (d < 0.0) ? ceil(d * p10 - 0.5)/p10 : floor(d * p10 + 0.5)/p10;
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
