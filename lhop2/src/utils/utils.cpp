/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// miscellanous classes and functions

#include "utils/platform.h"
#include "utils.h"
#include "utils/misc.h"
#include "utils/matching/matching.h"
#include "utils/graphs/graph_utils.h"

#include <string.h>
#include <cstdlib> 
#if defined WIN32 | defined WIN64
#include <windows.h>
#include <Psapi.h>
#else
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <fnmatch.h>
#endif


#include <iostream>
#include <opencv2/opencv.hpp>


using namespace std;

// region
///////////////////////////////////////////////////////////////////////////////

void region2::add(const region2& r, const ipoint2& d, double factor)
{
    for (const_iter_t iter = r.points.begin(); iter != r.points.end(); ++iter) 
        points.insert((*iter + d) * factor);
}

void region2::add(const region2& r, const ipoint2& d)
{
    for (const_iter_t iter = r.points.begin(); iter != r.points.end(); ++iter) 
        points.insert(*iter + d);
}

void region2::assign(const matrix<double>& m, double thresh)
{
    points.clear();
    for_each_xy(m, i, j) {
        if (m(i, j) >= thresh) points.insert(ipoint2(i, j));
    }
}


// K_bin and shape context histogram
///////////////////////////////////////////////////////////////////////////////

// nangles: number of angles
// r: sorted vector of radii boundaries
K_bin::K_bin(int nangles, const vector<int>& rv)
{
    init_matrix(nangles, rv);
}

K_bin::K_bin()
{
    vector<int> rv;

    rv.push_back(20);
    init_matrix(6, rv);
}

K_bin::K_bin(int nangles, int r1)
{
    vector<int> rv;

    rv.push_back(r1);
    init_matrix(nangles, rv);
}

K_bin::K_bin(int nangles, int r1, int r2)
{
    vector<int> rv;

    rv.push_back(r1);
    rv.push_back(r2);
    init_matrix(nangles, rv);
}

K_bin::K_bin(int nangles, int r1, int r2, int r3)
{
    vector<int> rv;

    rv.push_back(r1);
    rv.push_back(r2);
    rv.push_back(r3);
    init_matrix(nangles, rv);
}

void K_bin::reset(int nangles, const vector<int>& rv)
{
    init_matrix(nangles, rv);
}

// nangles: number of angles
// r: sorted vector of radii boundaries
//
//   radii:  r[size - 1]  r[0]    0     r[0]   r[size - 1]
//              |.:........:......:......:..........:.|
//
//
void K_bin::init_matrix(int nangles, const vector<int>& rv)
{
    double angle = 2.0*M_PI/nangles;
    double adelta = angle/2.0;
    int msize = 2*rv.back() + 3;

    m.resize(msize, msize);
    nbins = (int)rv.size()*nangles;
    
    int mc = m.width/2;

    for (int i = 0; i < msize; ++i) {
        for (int j = 0; j < msize; ++j) {
            int x = i - mc, y = j - mc;

            if (x == 0 && y == 0) {
                m(i, j) = 0;
            } else {
                double a = atan2((double)y, (double)x) + adelta;

                if (a < 0) a += 2.0*M_PI;
                
                int aindex = (int)(a/angle);

                if (aindex >= nangles) aindex = nangles - 1; // just in case

                double r = sqrt((double)(x*x + y*y));
                int rindex = (int)rv.size();

                for (int k = 0; k < (int)rv.size(); ++k) {
                    if (r < rv[k]) { 
                        rindex = k; 
                        break; 
                    }
                }
                m(i, j) = rindex*nangles + aindex;
            }
        }
    }
}

// Return bin or -1
int K_bin::get_bin(const ipoint2& pt, const ipoint2& center /* = ipoint2::zero */) const
{
    ipoint2 p = pt + ipoint2((int)m.width/2 - center.x, (int)m.height/2 - center.y);

    if (p.x >= 0 && p.x < (int)m.width && p.y >= 0 && p.y < (int)m.height) {
        int b = m(p.x, p.y);

        if (b < nbins) return b;
    }
    return -1;
}

void K_bin::get_histogram(vector<int>& h, const ip2_vector& pts, const ipoint2& center /* = ipoint2::zero */) const
{
    ipoint2 pp((int)m.width/2 - center.x, (int)m.height/2 - center.y);

    h.resize(nbins, 0);
    for (ip2_vector::const_iterator iter = pts.begin(); iter != pts.end(); ++iter) {
        ipoint2 p = *iter + pp;

        if (p.x >= 0 && p.x < (int)m.width && p.y >= 0 && p.y < (int)m.height) {
            int bin = m(p.x, p.y);
            
            if (bin < nbins) ++h[bin];
        }
    }
}

void K_bin::print() const
{
    m.print_mathematica(cout);
}


void normalize_histogram(vector<double>& hn, const vector<int>& h)
{
    hn.resize(h.size(), 0.0);

    int count = 0;

    for (int i = 0; i < h.size(); ++i) 
        count += h[i];
    if (count == 0) 
        return;
    for (int i = 0; i < h.size(); ++i)
        hn[i] = (double)h[i]/count;

}
double sc_histogram_distance(const vector<double>& h1, const vector<double>& h2)
{
    size_t maxk = min<int>(h1.size(), h2.size());
    double result = 0.0;

    for (size_t k = 0; k < maxk; ++k) {
        double den = h1[k] + h2[k];
        
        if (den != 0.0) {
            double x = h1[k] - h2[k];

            result += x*x/den;
        }
    }
    return result/2.0;
}


// M_bin
///////////////////////////////////////////////////////////////////////////////

M_bin::M_bin(const string& fname) : m(), nbins(), centers(), areas(), max_distance2(), user_value(nullptr)
{
    matrix<int> fm = read_int_matrix(fname);

    if (!fm.empty()) init_matrix(fm);
}

// Get bin at (i, j), note that (0, 0) is center of "matrix".
int M_bin::get_bin(int i, int j) const
{
    int x = (int)m.width/2 + i, y = (int)m.height/2 + j;

    if (x < 0 || x >= (int)m.width || y < 0 || y >= (int)m.height)
        return -1;
    return m(x, y);
}

void M_bin::init_matrix(const matrix<int>& pm)
{
    m = pm;
    nbins = 0;
    for_each_element (m, i) {
        int b = m[i];

        if (b >= 0 && nbins <= b) nbins = b + 1;
    }
    centers.resize(nbins);
	areas.resize(nbins);
	max_distance2.resize(nbins);

    ipoint2 c((int)m.width/2, (int)m.height/2);

    for (int i = 0; i < nbins; ++i) {
        int count = 0;
        ipoint2& p = centers[i];
        
		int max_distance_px = 0;

        p = ipoint2::zero;
        for_each_xy_int (m, x, y) {
            int b = m(x, y);

            if (b == i) { p.add(x - c.x, y - c.y); ++count; }
        }

        if (count > 0) p /= count;

		float dx = 0, dy = 0;

		for_each_xy_int (m, x, y) {
            int b = m(x, y);
            if (b == i) {
				max_distance_px = max<int>(max_distance_px, p.distance2(ipoint2(x,y) - c));
			}
        }
		max_distance2[i] = max_distance_px;
		areas[i] = count;
    }
	
}

M_bin* M_bin::get_resized(int w, int h) const{
	printf("resizing bin from %d,%d to %d,%d\n",width(), height(), w,h);
	matrix<int> new_m(w,h);

	float width_scale = m.width / (float)new_m.width;
	float height_scale = m.height / (float)new_m.height;

	for (int y = 0; y < (int)(new_m).height; ++y) {
		for (int x = 0; x < (int)(new_m).width; ++x) {
			int old_x = max<int>((int)min<int>((int)ceil(x * width_scale - 0.5), (int)m.width-1), (int)0);
			int old_y = max<int>((int)min<int>((int)ceil(y * height_scale - 0.5), (int)m.height-1), (int)0);
			new_m(x,y) = m(old_x, old_y);
		}
	}
	M_bin* new_bin = new M_bin(new_m);
	new_bin->nbins = nbins;
	return new_bin;
}
// functions
///////////////////////////////////////////////////////////////////////////////

vector<ipoint2> inhibit_point_set(const vector<ipoint2>& pts, int radius)
{
    if (pts.empty()) return vector<ipoint2>();
    if (radius <= 0) return pts;

    irectangle2 box;
    vector<ipoint2> result;

    for (vector<ipoint2>::const_iterator piter = pts.begin(); piter != pts.end(); ++piter) 
        box.eat(*piter);

    matrix<bool> m(box.x_dim() + 1, box.y_dim() + 1);

    m.fill(true);
    for (vector<ipoint2>::const_iterator piter = pts.begin(); piter != pts.end(); ++piter) {
        ipoint2 p = *piter - box.ll;

        if (m(p.x, p.y)) {
            result.push_back(*piter);
            m.set_region_circ(p.x, p.y, radius, false);
        }
    }
    return result;
}

vector<pair<int, ipoint2> > inhibit_point_set(const vector<pair<int, ipoint2> >& pts, int radius)
{
    typedef pair<int, ipoint2> item_t;
    typedef vector<item_t> result_t;
    typedef map<int, vector<ipoint2> > map_t;

    if (pts.empty()) return result_t();
    if (radius <= 0) return pts;

    map_t ptsm;
    result_t result;

    result.reserve(pts.size());
    for (result_t::const_iterator piter = pts.begin(); piter != pts.end(); ++piter) {
        ptsm[piter->first].push_back(piter->second);
    }
    for (map_t::iterator miter = ptsm.begin(); miter != ptsm.end(); ++miter) {
        vector<ipoint2> ipts = inhibit_point_set(miter->second, radius);

        for (vector<ipoint2>::iterator viter = ipts.begin(); viter != ipts.end(); ++viter) 
            result.push_back(item_t(miter->first, *viter));
    }
    return result;
}
void gaussian_mask(int dimx, int dimy, double sigma, matrix<double>& mask)
{
    mask.resize(dimx, dimy);

    int centerx = dimx/2;
    int centery = dimy/2;

    for_each_xy(mask, x, y) { 
        double d = (x - centerx)*(x - centerx) + (y - centery)*(y - centery);
        mask(x, y) = ::exp(-d/2.0/sigma/sigma);
    }
}
// Sum of dist(v1[i], v2[perm[i]]), i = 0, ..., v1.size()
// Returns -1 if perm[i] is not a valid index for v2 for some i
double distance_cost(const vector<int>& perm, const ip2_vector& v1, const ip2_vector& v2)
{
    double result = 0.0;

    if (perm.size() != v1.size()) return -1.0;
    for (int i = 0; i < (int)v1.size(); ++i) {
        int p = perm[i];

        if (p < 0 || p >= v2.size()) return -1.0;
        result += sqrt((double)v1[i].distance2(v2[p]));
    }
    return result;
}

// Calculates assignment (permutation) of points v1 to points v2 such that the sum of 
// distances is minimal
// Solves by examining all possibilities, so do not use it if v1 and v2 are large!
// Returns -1.0 if |v1| != |v2|
// Implement Hungarian method ;)
double min_distance_matching(vector<int>& perm, const ip2_vector& v1, const ip2_vector& v2)
{
    if (v1.size() != v2.size()) return -1.0;

    vector<int> tmpp = vector_range(0, (int)v1.size() - 1);
    double bestval = numeric_limits<double>::max();

    do {
        double val = distance_cost(tmpp, v1, v2);

        if (val < bestval) {
            bestval = val;
            perm = tmpp;
        }
    } while (next_permutation(tmpp.begin(), tmpp.end()));
    return bestval;
}

vector<int> identity_permutation(int n)
{
    vector<int> result(n, 0);

    for (int i = 0; i < n; ++i) 
        result[i] = i;
    return result;
}

void random_seed(int seed)
{
	srand((unsigned)seed);
}

void random_permutation(vector<int>& perm, int n)
{
    int* arr = new int[n];
    int i, p, t;

    for (i = 0; i < n; ++i) arr[i] = i; 
    for (i = n - 1; i > 0; --i) {
        p = rand()%(i + 1);
        t = arr[p];
        arr[p] = arr[i];
        arr[i] = t;
    }
    array_to_vector(perm, arr, n);
    delete[] arr;
}

int random_discrete(const vector<double>& d)
{
    if (d.empty())
        return -1;

    double sum = 0;

    for (size_t i = 0; i < d.size(); ++i) sum += d[i];

    double rnd = random_real()*sum;
    int result = 0;

    sum = d[0];
    while (sum < rnd) sum += d[++result];
    return result;
}

// "Normalize" vector (translate & scale).
pair<double, dpoint2> translate_and_scale(vector<dpoint2>& v)
{
    if (v.empty()) return pair<double, dpoint2>();

    dpoint2 avg = dpoint2::center(v.begin(), v.end());
    double scale = 0.0;

    for (vector<dpoint2>::iterator iter = v.begin(); iter != v.end(); ++iter) {
        scale += iter->distance2(avg);
    }
    scale = ::sqrt(scale/(double)v.size());
    for (vector<dpoint2>::iterator iter = v.begin(); iter != v.end(); ++iter) {
        *iter -= avg;
        iter->div(scale);
    }
    return pair<double, dpoint2>(scale, avg);
}

pair<double, dpoint2> translate_and_scale(vector<pair<int, ipoint2> >& pts)
{
    vector<ipoint2> v = extract_second<ipoint2>(pts.begin(), pts.end());
    auto result = translate_and_scale(v, 100.0);

    replace_second(pts, v);
    return result;
}

pair<double, dpoint2> translate_and_scale(vector<ipoint2>& v, double factor /* = 100 */)
{
    vector<dpoint2> pts = cast_vector<dpoint2, ipoint2>(v);
    auto result = translate_and_scale(pts);

    v = cast_vector<ipoint2, dpoint2>(pts, factor);
    return result;
}

// Rotate vector 'v' for angle 'a' about origin.
void rotate(vector<dpoint2>& v, double a)
{
    double cosa = cos(a);
    double sina = sin(a);

    for (vector<dpoint2>::iterator iter = v.begin(); iter != v.end(); ++iter) {
        iter->set(cosa*iter->x - sina*iter->y, sina*iter->x + cosa*iter->y);
    }
}

// Optimal rotation of v -> ref.
double optimal_rotation(const vector<dpoint2>& v, const vector<dpoint2>& ref)
{
    if (v.empty() || v.size() != ref.size()) 
        return 0.0;

    double numer = 0.0, denom = 0.0;
    
    for (vector<dpoint2>::const_iterator viter = v.begin(), refiter = ref.begin(); viter != v.end(); ++viter, ++refiter) {
        numer += viter->x*refiter->y - viter->y*refiter->x;
        denom += viter->x*refiter->x + viter->y*refiter->y;
    }
    return atan2(numer, denom);
}

bool delete_file(const string& path)
{
    return remove(path.c_str()) == 0;
}

unsigned long get_process_id()
{
#if defined WIN32 | defined WIN64
    return GetCurrentProcessId();
#else
	return getpid();
#endif
}

unsigned long get_working_set_size()
{
#if defined WIN32 | defined WIN64
    PROCESS_MEMORY_COUNTERS lMemInfo; 

    GetProcessMemoryInfo(GetCurrentProcess(), &lMemInfo, sizeof(lMemInfo));
    return (unsigned long)lMemInfo.WorkingSetSize;
#else
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    //http://www.gnu.org/s/libc/manual/html_node/Resource-Usage.html
    return (unsigned long) (usage.ru_idrss + usage.ru_ixrss + usage.ru_isrss) * 1024;
#endif
}

void to_unix_path(string& s)
{
    for_each(s.begin(), s.end(), [](char& c) { if (c == '\\') c = '/'; });
}

string wilcard_to_regex(const string& s)
{
    string result;

    for (string::const_iterator iter = s.begin(); iter != s.end(); ++iter) {
        if (*iter == '*') result += ".*";
#if defined WIN32 | defined WIN64
        else if (*iter == '\\') result += "/";
#endif
        else if (*iter == '?') result += '.';
        else if (*iter == '.') result += "\\.";
        else result += *iter;
    }
    return result;
}

#if defined WIN32 | defined WIN64

void list_directory_no_delete(list<string>& result, const string& pattern, const string& prefix)
{
    WIN32_FIND_DATA w32fd;
    HANDLE hFind = INVALID_HANDLE_VALUE;

    hFind = FindFirstFile(pattern.c_str(), &w32fd);
    if (hFind == INVALID_HANDLE_VALUE) 
        return;
    do {
        if (!(w32fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) 
            result.push_back(prefix + w32fd.cFileName);
    } while (FindNextFile(hFind, &w32fd) != 0);
    FindClose(hFind);
}

void list_directory(list<string>& names, const string& rpath)
{
    names.clear();
    list_directory_no_delete(names, rpath, "");
}

const char* get_current_path() {
	int buffer_length = GetCurrentDirectory(0,nullptr);

#ifdef _UNICODE
	LPTSTR buffer = (LPTSTR)malloc(sizeof(wchar_t) * buffer_length);
#else
	LPTSTR buffer = (LPTSTR)malloc(sizeof(char) * buffer_length);
#endif
	
	GetCurrentDirectory(buffer_length, buffer);
#ifdef _UNICODE
	// THIS MAY NOT HAVE BEEN TESTED !!
	int wBuffer_len = WideCharToMultiByte(CP_ACP, 0, buffer, -1, nullptr, 0, nullptr, nullptr);

	LPTSTR wBuffer = (LPTSTR)malloc(sizeof(char) * wBuffer_len);
	WideCharToMultiByte(CP_ACP, 0, buffer, -1, wBuffer, wBuffer_len, nullptr, nullptr);

	return (char*)wBuffer;
#else
	return (char*)buffer;
#endif
}

#else

void list_directory(list<string>& names, const string& rpath)
{
    DIR *dp;
    struct dirent *ep;

	size_t sep_pos = rpath.find_last_of("/\\");
	
	string folder, filename;
	if (sep_pos < 0) {
		folder = "./";
		filename = rpath;
	} else {
		folder = rpath.substr(0,sep_pos);
		filename = rpath.substr(sep_pos+1);
	}
    const char* pattern = filename.c_str();

    dp = opendir(folder.c_str());

    if (dp != nullptr)
	{
        while (ep = readdir (dp)) {
			if (fnmatch(pattern, ep->d_name, 0) == 0)
				names.push_back(string(ep->d_name));
		}
        (void) closedir (dp);
    } else
        perror ("Couldn't open the directory");
}

const char* get_current_path() {
	return (char*)get_current_dir_name();
}


#endif

bool create_directory(const string& folder) {	
	bool ok = false;
#if defined WIN32 | defined WIN64		
	ok = CreateDirectory(folder.c_str(), nullptr) > 0 ? true : false; // returns 0 on error
#else
	ok = mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0 ? true : false; // return 0 on success
#endif
	return ok;
}
bool create_directories(const string& folder) {	
	
	bool ok = false;
	// create dir
	if (folder.length() > 0 ) {
		// if parent dir does not exists then call create directory on parent
		ok = create_directory(folder);

		if (ok == false) {
			const string parent = get_parent_of_path(folder);		
			if (parent.length() > 0 && create_directories(parent))
				ok = create_directory(folder);
		}
	}
	return ok;
}

const string get_filename_of_path(const string& path) {
	size_t sep_pos = path.find_last_of("/\\");
	return (sep_pos == string::npos) ? path : path.substr(sep_pos+1);
}

const string get_parent_of_path(const string& path) {
	size_t sep_pos = path.find_last_of("/\\");
	return (sep_pos == string::npos) ? "" : path.substr(0,sep_pos);
}

bool file_list_from_file(list<string>& l, const string& fname, const string& dir, const string& ext /* = "" */)
{
    if (ext.empty())
        list_from_file(l, fname, "");
    else {
        list<string> pl;

        list_from_file(pl, fname, "");
        for (list<string>::iterator iter = pl.begin(); iter != pl.end(); ++iter) {
            char buf[1024];
            list<string> fnl;

            sprintf(buf, ext.c_str(), iter->c_str());
            list_directory(fnl, dir + buf);
            l.insert(l.end(), fnl.begin(), fnl.end());
        }
    }
    return !l.empty();
}

bool file_lists_from_list(list<list<string> >& l, const list<string>& pl, const string& dir, const string& ext)
{
    for (list<string>::const_iterator iter = pl.begin(); iter != pl.end(); ++iter) {
        char buf[1024];
        list<string> fnl;

        sprintf(buf, ext.c_str(), iter->c_str());
        list_directory(fnl, dir + buf);
        l.push_back(fnl);
    }
    return !l.empty();
}

bool file_lists_from_file(list<list<string> >& l, const string& fname, const string& dir, const string& ext)
{
    list<string> pl;

    list_from_file(pl, fname, "");
    return file_lists_from_list(l, pl, dir, ext);
}
matrix<int> read_int_matrix(const string& s)
{
    int mwidth = 0;
    vector<int> v;

    ifstream is(s.c_str());
    
    bool fail = is.fail();
    int wcount = 0;

    while (!fail && !is.eof()) {
        char c;
        bool newline = false;

        do {
            c = is.get();
            if (c == '\n') newline = true;
        } while (c == ' ' || c == '\t' || c == '\n');
        if (!is.eof()) {

            if (!newline) ++wcount;
            else {
                if (mwidth == 0) { 
                    if (wcount != 0) mwidth = wcount; 
                } else if (mwidth != wcount) fail = true; 
                wcount = 1;
            }

            if (!fail) {
                int elt;

                is.putback(c);
                is >> elt;

                fail = is.fail();
                if (!fail) v.push_back(elt);
            }
        }
    }
    is.close();

    if (fail || mwidth == 0 || ((int)v.size()) % mwidth != 0) 
        return matrix<int>();

    matrix<int> m(mwidth, (int)v.size()/mwidth);
    
    m.assign(v.begin(), v.end());
    return m;
}


matrix<float> read_float_matrix(const string& s)
{
    int mwidth = 0;
    vector<float> v;

    ifstream is(s.c_str());
    
    bool fail = is.fail();
    int wcount = 0;

    while (!fail && !is.eof()) {
        char c;
        bool newline = false;

        do {
            c = is.get();
            if (c == '\n') newline = true;
        } while (c == ' ' || c == '\t' || c == '\n');
        if (!is.eof()) {

            if (!newline) ++wcount;
            else {
                if (mwidth == 0) { 
                    if (wcount != 0) mwidth = wcount; 
                } else if (mwidth != wcount) fail = true; 
                wcount = 1;
            }

            if (!fail) {
                float elt;

                is.putback(c);
                is >> elt;

                fail = is.fail();
                if (!fail) v.push_back(elt);
            }
        }
    }
    is.close();

    if (fail || mwidth == 0 || ((int)v.size()) % mwidth != 0) 
        return matrix<float>();

    matrix<float> m(mwidth, (int)v.size()/mwidth);
    
    m.assign(v.begin(), v.end());
    return m;
}


// Returns a vector of integers (a 'permutation') such that 
// pts1[perm[i]] <-> pts2[i] is optimal assignment.
vector<int> point_matching(const vector<dpoint2>& pts1, const vector<dpoint2>& pts2)
{
	if (pts1.empty() || pts2.empty()) return vector<int>();

    matrix<double> dist(pts2.size(), pts1.size());
    matrix<int> match(pts2.size(), pts1.size());
    vector<int> result(pts2.size());

    for (int i = 0; i < (int)pts1.size(); ++i) 
        for (int j = 0; j < (int)pts2.size(); ++j) 
            dist(j, i) = (double)pts1[i].distance2(pts2[j]);
    best_matching(match.ptr(0, 0), dist.height, dist.width, dist.ptr(0, 0));
    for (int i = 0; i < (int)pts1.size(); ++i) {
        for (int j = 0; j < (int)pts2.size(); ++j) 
            if (match(j, i) == 1) {
                result[j] = i;
            }
    }
    return result;
}

// Mean shift algorithm
// 'path': result; front() == init
// 'pv': vector of points
// 'init': initial value
// 'radius': neighborhood of current 'x'
// 'c': normal kernel parameter in e^{c t^2}
// 'eps', 'maxsteps': iteration limits
// returns last 'eps' i.e.: |x_n - x_{n-1}|
double mean_shift(vector<dpoint2>& path, const vector<dpoint2>& pv, const dpoint2& init, 
    double radius, double c, double eps, int maxsteps)
{
    double radius2 = radius*radius;
    double eps2 = eps*eps;
    double delta2 = numeric_limits<double>::infinity();
    dpoint2 x = init;

    path.clear();
    path.push_back(x);
    for (int s = 0; s < maxsteps; ++s) {
        dpoint2 num = dpoint2::zero;
        double denom = 0.0;

        for (vector<dpoint2>::const_iterator iter = pv.begin(); iter != pv.end(); ++iter) {
            const dpoint2& p = *iter;
            double dxn = (p - x).norm2();

            if (dxn <= radius2) {
                double k = ::exp(c*dxn);
                num += p*k;
                denom += k;
            }
        }
        num.x /= denom;
        num.y /= denom;
        delta2 = x.distance2(num);
        x = num;
        path.push_back(x);
        if (delta2 < eps2)
            break;
        //cout << x << endl;
    }
    return sqrt(delta2);
}
