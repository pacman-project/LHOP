// miscellaneous solution specific classes and functions

#pragma once

#ifndef __UTILS_H__
#define __UTILS_H__

#define _USE_MATH_DEFINES  // necessary for <math.h> to define M_PI,.... constants, etc. 
#include <math.h>
#include <climits>
#include <stdio.h>
#include <stdlib.h>
#include "utils/misc.h"
#include "utils/matrix.h"
#include "utils/structures.h"
#include <algorithm>

#ifdef max
#undef max
#endif

using namespace std;


// defines
///////////////////////////////////////////////////////////////////////////////

#define START_TIMER(t) clock_t t = clock()
#define STOP_AND_SHOW(t, msg) cout << "Evaluation time of \"" << msg << "\" is "\
    << (double)(clock() - t)/CLOCKS_PER_SEC << " sec." << endl


// constants
///////////////////////////////////////////////////////////////////////////////

const double angle_degree = M_PI/180.0;
const double angle_degree1 = 180.0/M_PI; 

// declarations
///////////////////////////////////////////////////////////////////////////////

struct wpoint;

// functions
///////////////////////////////////////////////////////////////////////////////

bool delete_file(const string& path); 

void list_directory(list<string>& names, const string& path);

bool create_directories(const string& folder);

const char* get_current_path();
const string get_parent_of_path(const string& path);
const string get_filename_of_path(const string& path);

bool file_list_from_file(list<string>& l, const string& fname, const string& dir, const string& ext = "");

bool file_lists_from_list(list<list<string> >& l, const list<string>& pl, const string& dir, const string& ext);

bool file_lists_from_file(list<list<string> >& l, const string& fname, const string& dir, const string& ext);

unsigned long get_process_id();

unsigned long get_working_set_size();

vector<ipoint2> inhibit_point_set(const vector<ipoint2>& pts, int radius);
vector<pair<int, ipoint2> > inhibit_point_set(const vector<pair<int, ipoint2> >& pts, int radius);

void gaussian_mask(int dimx, int dimy, double sigma, matrix<double>& mask);

matrix<int> read_int_matrix(const string& s);
matrix<float> read_float_matrix(const string& s);

double distance_cost(const vector<int>& perm, const ip2_vector& v1, const ip2_vector& v2);

double min_distance_matching(vector<int>& perm, const ip2_vector& v1, const ip2_vector& v2);

vector<int> identity_permutation(int n);


// Random

void random_seed(int seed);

void random_permutation(vector<int>& perm, int n);

int random_discrete(const vector<double>& d);

// "Procrustes"

pair<double, dpoint2> translate_and_scale(vector<dpoint2>& v);
pair<double, dpoint2> translate_and_scale(vector<pair<int, ipoint2> >& pts);
pair<double, dpoint2> translate_and_scale(vector<ipoint2>& v, double factor = 100);

void rotate(vector<dpoint2>& v, double a);

double optimal_rotation(const vector<dpoint2>& v, const vector<dpoint2>& ref);

int filled_region_size(const vector<ipoint2>& pts);


// Point matching (Hungarian method)

vector<int> point_matching(const vector<dpoint2>& pts1, const vector<dpoint2>& pts2);

// Misc

double mean_shift(vector<dpoint2>& path, const vector<dpoint2>& pv, const dpoint2& init, 
    double radius, double c, double eps, int maxsteps);


// templates

template<class T> int intersection_size(const set<T>& s1, const set<T>& s2);


// classes 
///////////////////////////////////////////////////////////////////////////////

// region
///////////

struct region2 {
public:
    typedef set<ipoint2> container_t;
    typedef container_t::iterator iter_t;
    typedef container_t::const_iterator const_iter_t;

protected:
    container_t points;

public:
    region2() : points() { }
    region2(const region2& r) : points(r.points) { }

    void insert(const ipoint2& p) { points.insert(p); }
    int intersection_size(const region2& r) const { return ::intersection_size(points, r.points); }

    void add(const region2& r) { points.insert(r.points.begin(), r.points.end()); }
    void add(const region2& r, const ipoint2& delta);
    void add(const region2& r, const ipoint2& delta, double factor);

    void assign(const matrix<double>& m, double thresh); 
};

// wpoint
///////////

struct wpoint {
    int i, j;
    double w;

    wpoint(int iv = 0, int jv = 0, double wv = 0.0) : i(iv), j(jv), w(wv) { }

    bool operator<(const wpoint& q) const 
    {
        if (i < q.i) return true; 
        else if (i > q.i) return false; 
        else return j < q.j;
    }

    bool operator>(const wpoint& q) const 
    {
        if (i > q.i) return true; 
        else if (i < q.i) return false; 
        else return j > q.j;
    }

    bool operator==(const wpoint& q) const 
    {
        return i == q.i && j == q.j;
    }

    wpoint operator-(const iipair& p) 
    {
        wpoint result = *this;

        result.i -= p.first;
        result.j -= p.second;
        return result;
    }


};

// grid_point
///////////////

class grid_point : public point2<int> {
public:
    static matrix<double> default_filter;

    int type;
    matrix<double>* filter;

    grid_point() : point2<int>(), type(0), filter(&default_filter) { }
    grid_point(int x, int y, int t) : point2<int>(x, y), type(t), filter(&default_filter) { }
    grid_point(int x, int y, int t,  matrix<double>* f) : point2<int>(x, y), type(t), filter(f) { }
    grid_point(const grid_point& p) : point2<int>(p), type(p.type), filter(p.filter) { }
};


class normal_distribution1 {
protected:
    double mean;
    double variance;
    double nfactor;

protected:
    void reset_nfactor() { nfactor = variance <= 0.0 ? 0.0 : 1/sqrt(2*M_PI*variance); }
public:
    normal_distribution1(double m = 0.0, double s2 = 1.0) : mean(m), variance(s2) { reset_nfactor(); }
    normal_distribution1(const pair<double, double>& p) : mean(p.first), variance(p.second) 
        { reset_nfactor(); }

    double get_mean() { return mean; }
    double get_variance() { return variance; }
    void reset_mean(double m) { mean = m; }
    void reset_variance(double v) { variance = v; reset_nfactor(); }
    void reset(double m, double v) { mean = m; variance = v; reset_nfactor(); }

    // Value of the probability density function
    double pdf_val(double x) const
    { 
        if (variance <= 0) return 0.0;
        double a = x - mean;
        return nfactor*exp(-a*a/2/variance);
    }

    // Value of the probability density function without the "normalization" factor, 
    // i.e. value is 1.0 for x = mean.
    double pdf_val1(double x) const
    {
        if (variance <= 0) return 0.0;
        double a = x - mean;
        return exp(-a*a/2/variance);
    }
};

class von_mises_distribution {
protected:
    double mu;
    double kappa;
    double nfactor;

protected:
    void reset_nfactor() 
    { 
        // normalizing constant is 2 pi BesselI(0, kappa);
        // here we use BesselI(0, x) = e^x/sqrt(2 pi kappa)
        nfactor = kappa <= 0 ? (2*M_PI) : (sqrt(2*M_PI)*exp(kappa)/kappa); 
    }

public:
    von_mises_distribution() : mu(0), kappa(0) { reset_nfactor(); }
    von_mises_distribution(double m, double k) : mu(m), kappa(k) { reset_nfactor(); }

    von_mises_distribution(const pair<double, double>& p) : mu(p.first), kappa(p.second) 
        { reset_nfactor(); }

    double get_mu() { return mu; }
    double get_kappa() { return kappa; }
    void reset_mu(double m) { mu = m; }
    void reset_kappa(double k) 
    { 
        kappa = k; 
        reset_nfactor();
    }

    void reset(double m, double k) { mu = m; kappa = k; reset_nfactor(); }

    // Value of the probability density function
    double pdf_val(double x) const
    { 
        return exp(kappa*cos(x - mu))/nfactor;
    }

    // Value of the probability density function st. the value is 1 for x = mu, 
    double pdf_val1(double x) const
    {
        return exp(kappa*cos(x - mu) - kappa);
    }
};


template<class T> struct online_mean {
protected:
    int n;
    T mean;
public:
    online_mean() : n(0), mean() { }
    void new_data(const T& x) 
    {
        T delta = x - mean;
        ++n;
        mean += delta/n;
    }
    int sample_count() const { return n; }
    T get_mean() const { return mean; }
};

typedef online_mean<dpoint2> dpoint2_mean;

// K_bin
///////////////////////////////////////////////////////////////////////////////

struct K_bin : public streamable {
protected:
    matrix<int> m;
    int nbins;
public:
    K_bin();
    K_bin(int nangles, const vector<int>& rv);
    K_bin(int nangles, int r1);
    K_bin(int nangles, int r1, int r2);
    K_bin(int nangles, int r1, int r2, int r3);

    void get_histogram(vector<int>& h, const ip2_vector& pts, const ipoint2& center = ipoint2::zero) const;
    int get_bin(const ipoint2& p, const ipoint2& center = ipoint2::zero) const;
    int bin_count() const { return nbins; }
    void print() const;
    int width() const { return (int)m.width; }
    int height() const { return (int)m.width; }
    matrix<int> get_matrix() const { return m; }
    void reset(int nangles, const vector<int>& rv);

	virtual streamable* make_instance() const { return new K_bin(); }
	virtual void read_from_stream(istreamer& is) {
		is.read(nbins);
		is.read(m);
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write(nbins);
		os.write(m);
	}

protected:
    void init_matrix(int nangles, const vector<int>& rv);
};

void normalize_histogram(vector<double>& hn, const vector<int>& h);

template<class T> void normalize_histogram(vector<T>& hn)
{
    T sum = 0;

    for (auto iter = hn.begin(); iter != hn.end(); ++iter)
        sum += *iter;
    if (sum != 0)
        for (auto iter = hn.begin(); iter != hn.end(); ++iter)
            *iter /= sum;
}

double sc_histogram_distance(const vector<double>& h1, const vector<double>& h2);

// M_bin
///////////////////////////////////////////////////////////////////////////////

// "Matrix" bin

class M_bin_user_object {
public:
    virtual ~M_bin_user_object() {}
};

struct M_bin {
protected:
    matrix<int> m;
    int nbins;
    vector<ipoint2> centers;  // Centers of bins; note that (0, 0) is center of matrix 'm'
	vector<int> areas;  // Centers of bins; note that (0, 0) is center of matrix 'm'
	vector<float> max_distance2;  // distance to farthest element belonging to bin	
	M_bin_user_object* user_value;
public:
	M_bin(const matrix<int>& pm) : user_value(nullptr) { init_matrix(pm); }
    M_bin(const K_bin& kb) : user_value(nullptr) { init_matrix(kb.get_matrix()); }
    M_bin(const string& fname);
    
	virtual ~M_bin() { if (user_value != nullptr) delete user_value; }

    void init(const matrix<int>& pm) { init_matrix(pm); }
    void init(const K_bin& kb) { init_matrix(kb.get_matrix()); }
    int bin_count() const { return nbins; }
    int get_bin(int i, int j) const;
    ipoint2 bin_center(int b) const { return (b >= 0 && b < nbins) ? centers[b] : ipoint2(INT_MAX, INT_MAX); }
	int bin_area(int b) const { return (b >= 0 && b < nbins) ? areas[b] : 0; }
	float bin_max_distance2(int b) const { return (b >= 0 && b < nbins) ? max_distance2[b] : 0; }
    int width() const { return (int)m.width; }
    int height() const { return (int)m.height; }
	M_bin* get_resized(int w, int h) const;
	
	M_bin_user_object* get_user_value() {
		return user_value;
	}
	void set_user_value(M_bin_user_object* user_value) {
		this->user_value = user_value;
	}
protected:
    void init_matrix(const matrix<int>& pm);
};


// global templates
///////////////////////////////////////////////////////////////////////////////

template<class T, class S> void cast_vector(vector<T>& dest, const vector<S>& src)
{
    dest.resize(src.size());

    typename vector<S>::const_iterator siter = src.begin();
    typename vector<T>::iterator diter = dest.begin();

    while (siter != src.end()) *(diter++) = (T)(*(siter++));
}

template<class T, class S> void cast_vector(vector<T>& dest, const vector<S>& src, const S& factor)
{
    dest.resize(src.size());

    typename vector<S>::const_iterator siter = src.begin();
    typename vector<T>::iterator diter = dest.begin();

    while (siter != src.end()) *(diter++) = (T)(factor*(*(siter++)));
}

template<class T, class S> vector<T> cast_vector(const vector<S>& src)
{
    vector<T> dest(src.size());
    typename vector<S>::const_iterator siter = src.begin();
    typename vector<T>::iterator diter = dest.begin();

    while (siter != src.end()) *(diter++) = (T)(*(siter++));
    return dest;
}

template<class T, class S, class F> vector<T> cast_vector(const vector<S>& src, const F& factor)
{
    vector<T> dest(src.size());
    typename vector<S>::const_iterator siter = src.begin();
    typename vector<T>::iterator diter = dest.begin();

    while (siter != src.end()) *(diter++) = (T)((*(siter++))*factor);
    return dest;
}

template<class T> vector<T> take(const vector<T>& v, const vector<int>& ind)
{
    vector<T> result(ind.size());

    for (size_t i = 0; i < ind.size(); ++i) {
        result[i] = v[ind[i]];
    }
    return result;
}

template<class T, class I> vector<T> extract_second(I begin, I end)
{
    vector<T> result;

    while (begin != end) {
        result.push_back(begin->second);
        ++begin;
    }
    return result;
}

template<class T, class I> vector<T> extract_first(I begin, I end)
{
    vector<T> result;

    while (begin != end) {
        result.push_back(begin->first);
        ++begin;
    }
    return result;
}

template<class T, class I, class F> vector<T> extract(I begin, I end, F fun)
{
    vector<T> result;

    while (begin != end) {
        result.push_back(fun(*begin));
        ++begin;
    }
    return result;
}

template<class A, class B> void replace_second(vector<pair<A, B> >& v1, 
    const vector<B>& v2)
{
    size_t maxi = min(v1.size(), v2.size());

    for (size_t i = 0; i < maxi; ++i) 
        v1[i].second = v2[i];
    if (v1.size() > maxi) 
        v1.resize(maxi);
}

template<class T, class S> void divide_vector(vector<T>& v, const S& d)
{
    for (typename vector<T>::iterator iter = v.begin(); iter != v.end(); ++iter) 
        (*iter) /= (T)d;
}

template<class T, class S> void mult_vector(vector<T>& v, const S& m)
{
    for (typename vector<T>::iterator iter = v.begin(); iter != v.end(); ++iter) 
        (*iter) *= m;
}

template<class T> T median_value(vector<T> val)
{
    sort(val.begin(), val.end());
    return val[val.size()/2];
}

template<class T> T median_value(vector<T> val, double perc)
{
    if (val.empty()) return T();
    sort(val.begin(), val.end());
    size_t pos = (size_t)(perc*(val.size() - 1));
    return val[pos];
}

template<class Key, class T> void get_key_set(set<Key>& kset, const map<Key, T>& m)
{
    kset.clear();
    for (typename map<Key, T>::const_iterator miter = m.begin(); miter != m.end(); ++miter) {
        kset.insert(miter->first);
    }
}

template<class T> bool list_from_file(list<T>& l, const string& fname, const string& prefix)
{
    ifstream is(fname.c_str());

    l.clear();
    while (!is.fail()) {
        T s;
        is >> s;
        if (!is.fail()) l.push_back(s);
    }
    return !l.empty();
}

template<class T, class D> double hierarchical_clustering(vector<vector<T> >& clist, const vector<T>& list, D distance)
{
    clist.clear();
    for (typename vector<T>::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
        clist.push_back(vector<T>());
        clist.back().push_back(*iter);
    }
    return hierarchical_clustering_iter(clist, distance);
}

template<class T, class D> double hierarchical_clustering_iter(vector<vector<T> >& clist, D distance)
{
    int n = (int)clist.size();
    double min = numeric_limits<double>::max();

    if (n <= 1) return min;

    int mini = -1, minj = -1;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double d = distance(clist[i], clist[j]);
            
            if (d < min) { min = d; mini = i; minj = j; }
        }
    }
    clist.push_back(vector<T>());

    vector<T>& last = clist.back();

    last.insert(last.end(), clist[mini].begin(), clist[mini].end());
    last.insert(last.end(), clist[minj].begin(), clist[minj].end());
    if (mini > minj) {
        clist.erase(clist.begin() + mini);
        clist.erase(clist.begin() + minj);
    } else {
        clist.erase(clist.begin() + minj);
        clist.erase(clist.begin() + mini);
    }

    return min;
}

inline int distance2(int x1, int y1, int x2, int y2)
{
    int xd = x2 - x1;
    int yd = y2 - y1;

    return xd*xd + yd*yd;
}

template <typename T> struct max_function {
    T operator()(const T& a, const T& b) const { return max<T>(a, b); }
};

template<class T> void array_to_vector(vector<T>& v, const T arr[], int size)
{
    v.resize(size);
    typename vector<T>::iterator iter = v.begin();

    for (int i = 0; i < size; ++i) 
        *(iter++) = arr[i];
}

template<class T> int intersection_size(const set<T>& s1, const set<T>& s2)
{
    typename set<T>::const_iterator i1 = s1.begin(), i2 = s2.begin();
    int result = 0;

    while (i1 != s1.end() && i2 != s2.end()) {
        if (*i1 < *i2) ++i1; 
        else if (*i2 < *i1) ++i2;
        else {
            ++result;
            ++i1; ++i2;
        }
    }
    return result;
}

// not used
template<class T> void set_intersection(set<T>& intersection, const set<T>& s1, const set<T>& s2)
{
    typename set<T>::const_iterator i1 = s1.begin(), i2 = s2.begin();

    intersection.clear();
    while (i1 != s1.end() && i2 != s2.end()) {
        if (*i1 < *i2) ++i1; 
        else if (*i2 < *i1) ++i2;
        else {
            intersection.insert(*i1);
            ++i1; ++i2;
        }
    }
}

template<class T> void set_difference(set<T>& d, const set<T>& s1, const set<T>& s2)
{
    typename set<T>::const_iterator i1 = s1.begin(), i2 = s2.begin();

    while (i1 != s1.end() && i2 != s2.end()) {
        if (*i1 < *i2) { d.insert(*i1); ++i1; }
        else if (*i2 < *i1) ++i2;
        else { ++i1; ++i2; }
    }
    while (i1 != s1.end()) {
        d.insert(*i1);
        ++i1;
    }
}

template<class T, class I> void listable_to_vector(vector<T>& v, I begin, I end)
{
    v.clear();
    for (;begin != end; ++begin) {
        v.push_back(*begin);
    }
}

template<class T, class I> T average(I begin, I end)
{
    int n = 0;
    T result(0);

    for (; begin != end; ++begin, ++n) result += *begin;
    return (n == 0) ? T(0) : result/n;
}

// Square of the diameter of a set of point<T> objects.
template<class T, class I> T diameter2(I begin, I end)
{
    T result = 0;

    for (I i = begin; i != end; ++i) {
        for (I j = begin; j != end; ++j) {
            T d = i->distance2(*j);

            if (d > result) result = d;
        }
    }
    return result;
}

template<class T> void convolve_matrix(matrix<double>& result, matrix<T>& im, const matrix<double>& dm)
{
    result.resize(im.width, im.height, 0.0);

    if (dm.width != dm.height) return;

    int ddim = (int)dm.width;
    int ddim2 = (int)dm.width/2;
    int dd = ddim*ddim;

    T** sptr = new T*[dd];
    int i, j, k;

    k = 0;
    for (j = 0; j < ddim; ++j) 
        for (i = 0; i < ddim; ++i) 
            sptr[k++] = &im(i, j);

    double* dptr = &result(ddim2, ddim2);
    int mwidth = (int)im.width - ddim2;
    int mheight = (int)im.height - ddim2;
    double conv;

    for (j = ddim2; j < mheight; ++j) {
        for (i = ddim2; i < mwidth; ++i) {
            conv = 0.0;
            for (k = 0; k < dd; ++k) conv += dm[k] * (*sptr[k]);
            *dptr = conv;
            for (k = 0; k < dd; ++k) ++(sptr[k]);
            ++dptr;
        }
         for (k = 0; k < dd; ++k) sptr[k] += 2*ddim2;
         dptr += 2*ddim2;
    }
    delete sptr;
}

template<class T> void ipoint2_set_normalization_parameters(ipoint2& avg, double& scale, T begin, T end)
{
    avg.set(0, 0);
    scale = 0.0;

    if (begin == end)
        return;

    int count = 0;

    for (T iter = begin; iter != end; ++iter) {
        ipoint2 p = *iter;
        
        avg += p;
        ++count;
    }
    avg /= count;
    for (T iter = begin; iter != end; ++iter) {
        ipoint2 p = *iter - avg;
        
        scale += p.x*p.x + p.y*p.y;
    }
    scale = sqrt(scale/count);
}

template<class T> void normalize_ipoint2_set(T begin, T end, int sfactor)
{
    if (begin == end)
        return;

    double scale = 0.0;
    ipoint2 avg(0, 0);
    int count = 0;

    for (T iter = begin; iter != end; ++iter) {
        ipoint2 p = *iter;
        
        avg += p;
        ++count;
    }
    avg /= count;
    for (T iter = begin; iter != end; ++iter) {
        ipoint2 p = *iter - avg;
        
        scale += p.x*p.x + p.y*p.y;
    }
    scale = sqrt(scale/count);
    for (T iter = begin; iter != end; ++iter) {
        ipoint2& p = *iter; 
        
        p -= avg;
        p.x = (int)(sfactor*p.x/scale);
        p.y = (int)(sfactor*p.y/scale);
    }

}

template<class T> matrix<double> convolve_matrix(matrix<T>& m, const matrix<double>& k)
{
    matrix<double> result;

    convolve_matrix(result, m, k);
    return result;
}
 
// Auxiliary function for ordering (see below)
template <class T> struct ordering_less : binary_function<T, T, bool> {
    bool operator() (const T& x, const T& y) const { return *x.first < *y.first; }
};

// Auxiliary function for ordering_pred (see below)
template <class T, class Pred> struct ordering_less_pred {
    Pred p;
    ordering_less_pred(Pred pp) : p(pp) { }
    bool operator() (const T& x, const T& y) const { return p(*x.first, *y.first); }
};

// Returns ordering of vector range sorted in non-descending order with STL sort function.
template<class T> vector<int> ordering(typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end)
{
    typedef pair<const T*, int> pair_t;
    typedef vector<pair_t> vector_t;

    vector_t pv;
    int i = 0;

    for (typename vector<T>::const_iterator iter = begin; iter != end; ++iter) {
        pv.push_back(pair_t(&(*iter), i++));
    }
    sort(pv.begin(), pv.end(), ordering_less<pair_t>());
    
    vector<int> result;

    result.reserve(pv.size());
    for (typename vector_t::iterator viter = pv.begin(); viter != pv.end(); ++viter) 
        result.push_back(viter->second);
    return result;
}

// Returns ordering of vector range sorted in non-descending order with STL sort function.
template<class T, class Pred> 
vector<int> ordering(typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end, Pred p)
{
    typedef pair<const T*, int> pair_t;
    typedef vector<pair_t> vector_t;

    vector_t pv;
    int i = 0;

    for (typename vector<T>::const_iterator iter = begin; iter != end; ++iter) {
        pv.push_back(pair_t(&(*iter), i++));
    }
    sort(pv.begin(), pv.end(), ordering_less_pred<pair_t, Pred>(p));
    
    vector<int> result;

    result.reserve(pv.size());
    for (typename vector_t::iterator viter = pv.begin(); viter != pv.end(); ++viter) 
        result.push_back(viter->second);
    return result;
}


// Permutes range of vector starting at 'begin' according to permutation 'p'.
template<class T> void permute_range(typename vector<T>::iterator begin, const vector<int>& p)
{
    vector<int> mem(p.size(), 0);

    for (int i = 0; i < (int)p.size(); ++i) {
        int pos = i;
        int ppos = p[pos];

        mem[pos] = 1;
        while (mem[ppos] == 0) {
            iter_swap(begin + pos, begin + ppos);
            mem[ppos] = 1;
            pos = ppos;
            ppos = p[ppos];
        }
    }
}


// Returns a HOP_REAL random number from interval [0, 1).
inline double random_real() 
{
    return (double)rand() / ((double)RAND_MAX + 1);
}

inline double random_real(double min, double max) 
{
    return random_real()*(max - min) + min;
}

// Returns a random integer from interval [min, max).
inline int random_int(int min, int max)
{
    return (int)(random_real()*(max - min) + min);
}

inline double normal_random()
{
    double U = (double)(rand() + 1)/((double)RAND_MAX + 2);
    double V = (double)(rand() + 1)/((double)RAND_MAX + 2);

    return sqrt(-2.0*log(U))*cos(2*M_PI*V);
}

inline double normal_random(double mu, double sigma)
{
    return mu + sigma*normal_random();
}

// Utility function for performance timing in microseconds in windows only

#if defined WIN32 | defined WIN64

#include <windows.h>

typedef struct {
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
} stopWatch;

class CStopWatch {

private:
	double sum_time;
	stopWatch timer;
	LARGE_INTEGER frequency;
	double LIToSecs( LARGE_INTEGER & L)  {
		return ((double)L.QuadPart /(double)frequency.QuadPart) ;
	}
public:
	CStopWatch() {
		sum_time = 0;
		timer.start.QuadPart=0;
		timer.stop.QuadPart=0; 
		QueryPerformanceFrequency( &frequency ) ;
	}

	void startTimer() { QueryPerformanceCounter(&timer.start) ; }
	void stopTimer( ) { QueryPerformanceCounter(&timer.stop) ; }

	void sumElapsedTime() {
		sum_time += getElapsedTime();
	}
	double getElapsedTime()  {
		LARGE_INTEGER time;
		time.QuadPart = timer.stop.QuadPart - timer.start.QuadPart;
		return LIToSecs( time) ;
	}
	double getSumTime()  {
		return sum_time;
	}
	void resetSumTime() {
		sum_time = 0;
	}

};
#else

#include<ctime>
#include<time.h>

#include <stdlib.h>
#include <sys/time.h>

class CStopWatch {

private:
	double sum_time;
	struct timeval start;
	struct timeval stop;
	struct timezone tz;
	clock_t start1,stop1;
public:
	CStopWatch() {
		sum_time = 0;
		start.tv_usec = 0;
		stop.tv_usec = 0;
	}
	void startTimer() {
		gettimeofday(&start, &tz);
		start1 = clock();
	}
	void stopTimer( ) {
		gettimeofday(&stop, &tz);
		stop1 = clock();
	}
	double getElapsedTime()  {
		//return (double)(stop.tv_usec-start.tv_usec)*1e-6;
		return (double)(stop1-start1)/CLOCKS_PER_SEC;
	}
	void sumElapsedTime() {
		sum_time += getElapsedTime();
	}
	double getSumTime()  {
		return sum_time;
	}

	void resetSumTime() {
		sum_time = 0;
	}
};
#endif


#endif /* __UTILS_H__ */




