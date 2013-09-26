// misc.h

#pragma once
#ifndef _MISC_H_
#define _MISC_H_

#define NOMINMAX
//#define PI 3.14159265

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <limits.h>

#if defined WIN32 | defined WIN64
#include <windows.h>
#endif

#include "structures.h"
#include "streaming.h"
#include "exceptions.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// defines & typedefs
///////////////////////////////////////////////////////////////////////////////

#define foreach_point_in_disk(_x, _y, _radius) \
    for (int _x = -(_radius); _x <= (_radius); ++_x) \
        for (int _y = -(_radius); _y <= (_radius); ++_y) \
            if (_x*_x + _y*_y <= (_radius)*(_radius))

template <typename T>
struct vectorn {
    typedef vector<vector<vector<T> > > type3;
    typedef vector<vector<T> > type2;
    typedef vector<T> type1;
};

#define foreach_vector3(_v, _i, _j, _k) \
    for (int _i = 0; _i < _v.size(); ++_i) \
        for (int _j = 0; _j < _v[_i].size(); ++_j) \
            for (int _k = 0; _k < _v[_i][_j].size(); ++_k) 

#define foreach23_vector3(_v, _i, _j, _k) \
    for (int _j = 0; _j < _v[_i].size(); ++_j) \
        for (int _k = 0; _k < _v[_i][_j].size(); ++_k) 

template<typename T> T& at(typename vectorn<T>::type3& v3, int i, int j, int k, const T& defval)
{

    if (i >= v3.size()) v3.resize(i + 1, typename vectorn<T>::type2());
    if (j >= (v3[i]).size()) (v3[i]).resize(j + 1, typename vectorn<T>::type1());
    if (k >= (v3[i][j]).size()) (v3[i][j]).resize(k + 1, defval);
    return v3[i][j][k];
}

template<typename T> T& at(typename vectorn<T>::type2& v2, int i, int j, const T& defval)
{
    if (i >= v2.size()) v2.resize(i + 1, typename vectorn<T>::type1());
    if (j >= (v2[i]).size()) (v2[i]).resize(j + 1, defval);
    return v2[i][j];
}

// streamable_t<T>
///////////////////////////////////////////////////////////////////////////////

template<class T> class streamable_t : public streamable {
public:
    T data;

    streamable_t() { }
    streamable_t(const T& d) : data(d) { }
    streamable_t(const streamable_t& s) : data(s.data) { }

    virtual void copy_to(streamable* p, cloner& cl) { streamable::copy_to(p, cl); ((streamable_t*)p)->data = data; }
    virtual streamable* make_instance() const { return new streamable_t(); }
    virtual void read_from_stream(istreamer& is) { streamable::read_from_stream(is); is.read(data); }
    virtual void write_to_stream(ostreamer& os) { streamable::write_to_stream(os); os.write(data); }
};

typedef streamable_t<int> streamable_int;

// predicates
///////////////

class bfs_predicate {
public:
    bfs_predicate() { }
    virtual bool value(int type, int dist, int nbcount) const { return true; }
};

template<class A, class B, class R> struct bin_function {
    virtual R operator()(const A&, const B&) const = 0;
};

template<class T, class U> struct pair1_greater {
    bool operator()(const pair<T, U>& p1, const pair<T, U>& p2) const 
    {
        return p1.first > p2.first;
    }
};

template<class T> struct geq_predicate {
    T val;

    geq_predicate(const T& v) : val(v) { }
    bool operator()(const T& t) { return t >= val; }
};

template<class T> struct leq_predicate {
    T val;

    leq_predicate(const T& v) : val(v) { }
    bool operator()(const T& t) { return t <= val; }
};


// 

static inline int ilog2_x86(unsigned int x) 
{
    int retval;
 
#ifdef WIN64
      retval = (int)log((double)x);
#elif defined WIN32
	__asm { bsr eax, x
              mov retval, eax
      }
#else
      asm ( "bsr %%eax, %0 ; \n mov %0, %%eax ; \n" : "=r"(retval) 
           : "r"(x)
           : "%eax" );
#endif
	return retval;
}

// mem_allocator
///////////////////////////////////////////////////////////////////////////////

#define CLASS_ALLOCATION()\
    private:\
        static mem_allocator _allocator;\
    public:\
        void* operator new(size_t count) { return _allocator.alloc_mem(); }\
        void operator delete(void* object) { _allocator.free_mem(object); }
        
#define DEFINE_CLASS_ALLOCATOR(_class) mem_allocator _class::_allocator(sizeof(_class), 1024, #_class);
#define DEFINE_CLASS_ALLOCATOR_2(_class, _size) mem_allocator _class::_allocator(sizeof(_class), _size, #_class);

// D E B U G

//static clock_t T_alloc_mem = 0;
//static clock_t T_free_mem = 0;
//static clock_t T_alloc_mem2 = 0;
//static clock_t T_free_mem2 = 0;
//static clock_t T_start, T_end;




//     +------- free_slot
//     |      +------- position of the data
//     V      V
//    -+------+-----------+------+-----------------+
//     | mem0 | mem0_next | next |......           |
//    -+------+-----------+------+-----------------+
//
//     \_______ info_block ______/
//
//     \_______________ slot_size _________________/

class mem_allocator {
protected:
    struct info_block { char* mem0; char* mem0_next; size_t next; };
    struct info_block0 { size_t size; size_t bit_size; };

    const char* name;
    char** mem;
    char* free_slot;   
    size_t slot_size; 
    size_t slot_count;
    size_t mem_count;

public:
    mem_allocator(size_t ssize, size_t scount, const char* n);

	~mem_allocator();

	void* alloc_mem()
	{
		void* r;
		r = alloc_mem_r();
		return r;
	}

    void* alloc_mem_r() 
    {
        //cout << " allocating " << name << endl;		
        if (free_slot == nullptr) return nullptr; // throw something?
        info_block* info = (info_block*)free_slot;
        
        if (info->next == UINT_MAX) {
            expand_mem(); 
            return alloc_mem_r();
        } else {
            void* result = &(info->mem0_next);
            free_slot = info->mem0_next + slot_size * info->next;
            return result;
        }
    }

    void free_mem(void* p)
    {
        info_block* info = (info_block*)((char*)p - sizeof(info->mem0));
        info_block* free_info = (info_block*)free_slot;

        info->mem0_next = free_info->mem0;
        info->next = (free_slot - free_info->mem0)/slot_size;
        free_slot = (char*)info;
	}

    static void free_mem(void* p, mem_allocator** allocators) 
    {
        info_block* info = (info_block*)((char*)p - sizeof(info->mem0));
        info_block0* info0 = (info_block0*)info->mem0;
        mem_allocator* alloc = allocators[info0->bit_size];
        info_block* free_info = (info_block*)alloc->free_slot;

        info->mem0_next = free_info->mem0;
        info->next = (alloc->free_slot - free_info->mem0)/alloc->slot_size;
        alloc->free_slot = (char*)info;
    }

protected:
    void init_mem_range(char*& memp);  
    void expand_mem();

};

// mem_allocator_set
//////////////////////

class mem_allocator_set {
protected:
    mem_allocator** allocators;
    size_t count;

public:
    mem_allocator_set(size_t pcount = 10);
    ~mem_allocator_set();

    void* alloc_mem(size_t c) 
    {
        int i = ilog2_x86(c) + 1;
        if (i >= count) throw new_libhop_exception("Increase the number of mem_allocators.");
        return allocators[i]->alloc_mem();
    }

    void free_mem(void* p) 
    {
        mem_allocator::free_mem(p, allocators);
    }

};

// allocator_base
///////////////////

class allocator_base {
private:
    static mem_allocator_set allocators;

public:
    void* operator new(size_t count) { return allocators.alloc_mem(count); }
    void operator delete(void* object) { allocators.free_mem(object); }
};


// triple
/////////////

template<class T, class U, class V> 
struct triple {
    T first;
    U second;
    V third;

    triple() : first(T()), second(U()), third(V()) { }
    triple(const T& f, const U& s, const V& t) : 
        first(f), second(s), third(t) { }
    
    bool operator==(const triple& t) const { 
        return (first == t.first && second == t.second && third == t.third);
    }
    
    bool operator<(const triple& t) const {
        return (first < t.first || 
            !(t.first < first) && (second < t.second || !(t.second < second) && third < t.third));
    }
    
    bool operator>(const triple& t) const {
        return (t < *this);
    }
    
    bool operator!=(const triple& t) const {
        return !(*this == t);
    }

    friend ostream& operator<<(ostream& os, const triple& t)
    {
        os << t.first << ' ' << t.second << ' ' << t.third;
        return os;
    }

    static bool first_less(const triple& t1, const triple& t2) { return t1.first < t2.first; }
    static bool second_less(const triple& t1, const triple& t2) { return t1.second < t2.second; }
    static bool third_less(const triple& t1, const triple& t2) { return t1.third < t2.third; }

};

typedef triple<int,int,int> itriple;


// int_map
///////////////////////////////////////////////////////////////////////////////

class int_map {
protected:
    vector<int> f, finv;

public:
    int_map(int domain = 0, int codomain = 0) : f(domain), finv(codomain) { } 

    void reset(int domain, int codomain);
    void reset(int domain_range, const vector<int>& v);

    int operator()(int i) { return f[i]; }
    int inv(int i) { return finv[i]; }
    int domain() { return (int)f.size(); }
    int codomain() { return (int)finv.size(); }


protected:

    void make_identity() 
    { 
        int md = (int)f.size(), mcd = (int)finv.size();
        for (int i = 0; i < md; ++i) { f[i] = i % mcd; finv[i % mcd] = i; }
    }

};

// iipair & support
/////////////////////

typedef pair<int, int> iipair;

inline iipair operator-(const iipair& p, const iipair& q)
{
    return iipair(p.first - q.first, p.second - q.second);
}

inline iipair operator+(const iipair& p, const iipair& q)
{
    return iipair(p.first + q.first, p.second + q.second);
}

inline iipair& operator+=(iipair& p, const iipair& q)
{

    p.first += q.first;
    p.second += q.second;
    return p;
}

inline iipair& operator-=(iipair& p, const iipair& q)
{

    p.first -= q.first;
    p.second -= q.second;
    return p;
}

struct iipair_equal : public binary_function<iipair, iipair, bool> {
    int dist;
  
    iipair_equal(int d = 0) : dist(d) { }

    bool operator()(const iipair& p, const iipair& q) const
        { return abs((double) (p.first - q.first)) <= dist && abs((double) (p.second - q.second)) <= dist; }
};

// rectangle
//////////////

struct rectangle {
    int x0, y0;
    int x1, y1;

    rectangle(int vx0 = 0.0, int vy0 = 0.0, int vx1 = 0.0, int vy1 = 0.0) 
        : x0(vx0), y0(vy0), x1(vx1), y1(vy1) 
    { }

    bool inside(int x, int y) 
    { 
        return x0 <= x && x < x1 && y0 <= y && y < y1; 
    }
};

inline istream& operator>>(istream& is, rectangle& r) 
{
    is >> r.x0 >> r.y0;
    is >> r.x1 >> r.y1;
    return is;
}


// vector_parameter
/////////////////////

template<typename T> struct vector_parameter {
protected:
    T val;
    size_t val_i;
    vector<T> val_v;
public:
    vector_parameter() { val_v.resize(1); val = val_v[0]; val_i = 0; }
    //vector_parameter(const vector<T>& v) : val_v(v) { if 

    void set_val(const vector<T>& v) 
    {
        val_v.assign(v.begin(), v.end());
        if (val_v.empty()) 
            throw new_libhop_exception("Empty vector");
        val = v[0];
        val_i = 0;
    }

    void set_val(const T& v)
    {
        vector<T> vv;

        vv.push_back(v);
        set_val(vv);
    }

    T get_val() const { return val; }
    operator T() const { return val; }
    T operator[](int i) const 
    { 
        if (i < 0) return val_v[0]; 
        else if (i < (int)val_v.size()) return val_v[i];
        else return val_v.back();
    }

    void inc() { if (++val_i < val_v.size()) val = val_v[val_i]; }
    
};


// config_exception
/////////////////////

class config_exception : public libhop_exception {
public:
	config_exception(const string& file, const size_t line, const string& msg = "Unknown config_exception"): libhop_exception(file, line, msg, "config_exception") {}
};

// config_dictionary 
//////////////////////
#include <stdarg.h>

class config_dictionary {	
protected:
    typedef map<string, string> dictionary_t;

    dictionary_t dictionary;
public:
    config_dictionary();
    config_dictionary(const char*);
    config_dictionary(string);
    config_dictionary(const config_dictionary&);

    void replace_keys(const char*);
    void from_file(const char*);
	bool to_file(const string& filename);
    void from_string(const char*);
    void from_string_nl(const char*);
    void from_namespace(const config_dictionary& cfg, const string& nsname);
	void from_namespace_priority(const config_dictionary& cfg, int ns_size, ...);
	void update_namespace_references();
	bool does_namespace_exist(const string& nsname);
    bool get_value(double&, const string&, bool = false) const;
	bool get_value(bool&, const string&, bool = false) const;
    bool get_value(int&, const string&, bool = false) const;
    bool get_value(string&, const string&, bool = false) const;
    bool get_value(vector<int>&, const string&, bool = false) const;
    bool get_value(vector<double>&, const string&, bool = false) const;
    bool get_value(vector<bool>&, const string&, bool = false) const;
	bool get_value(vector<string>& v, const string& key, bool = false) const;
    template <typename T> void get_value(vector_parameter<T>&, const string&, const T&) const;
    double get_value_double(const string&, double, bool = false) const;
    int get_value_int(const string&, int, bool = false) const;
    string get_value_string(const string&, const string&, bool = false) const;
    bool get_value_bool(const string&, bool, bool = false) const;
    bool is_defined(const string&) const;
	void set_value(const string& key, const string& val);
    void set_value(const string& key, int val);
    void set_value(const string& key, double val);
	
protected:
    void parse_file(const char*);
    void parse_string(string);
    void parse_string_nl(string);
    void insert(const string& key, const string& val);
};

template <typename T> void config_dictionary::get_value(vector_parameter<T>& par, const string& key, 
        const T& def) const
{
    vector<T> v;

    get_value(v, key, false);
    if (v.empty()) v.push_back(def);
    par.set_val(v);
}

// tolerance objects
//////////////////////

template <class T> class tolerance_object : public streamable {
public:
    virtual double value(const T& x, const T& y) = 0;
    virtual matrix<double>* get_matrix() { return nullptr; }

    virtual streamable* make_instance() const = 0;
    virtual void read_from_stream(istreamer& is) { streamable::read_from_stream(is); }
    virtual void write_to_stream(ostreamer& os) { streamable::write_to_stream(os); }
};

template <class T> class circular_tolerance : public tolerance_object<T> {
protected:
    T radius2;
    matrix<double>* m;
public:
    circular_tolerance() : m(nullptr) { }
    circular_tolerance(T r) : radius2(r*r), m(nullptr) { }
    ~circular_tolerance() { if (m != nullptr) delete m; }

    void set_radius(T r) 
    { 
        if (r*r == radius2) return;
        radius2 = r*r; 
        if (m != nullptr) delete m; 
        m = nullptr;
    }

    virtual double value(const T& x, const T& y) { return x*x + y*y > radius2 ? 0.0 : 1.0; }
    virtual matrix<double>* get_matrix() 
    { 
        if (m == nullptr) {
            m = new matrix<double>();
            m->shape_circle((int)(sqrt((double)radius2)), 1.0, 0.0);
        }
        return m; 
    }

    virtual streamable* make_instance() const { return new circular_tolerance<T>(); };
    virtual void read_from_stream(istreamer& is) { streamable::read_from_stream(is); is.read(radius2); }
    virtual void write_to_stream(ostreamer& os) { streamable::write_to_stream(os); os.write(radius2); }
};

// average_learning
///////////////////////////////////////////////////////////////////////////////

// Oper is class with defined operators add_to(S, S) ("S += S"), div_int(S, int) ("S/int")
template<class T, class S, class Oper> class average_learning {
protected:
    typedef pair<S, int> stat_t;
    typedef map<T, stat_t> map_t;

    map_t map_;
    Oper op;
public:
    average_learning() : map_(), op() { }
    
    void update(const T& key, const S& p) 
    {
        stat_t& ps = map_[key];

        op.add_to(ps.first, p);
        ps.second++;
    }

    void update(const std::map<T, S>& m) 
    {
        for (typename std::map<T, S>::const_iterator miter = m.begin(); miter != m.end(); ++miter) 
            update(miter->first, miter->second);
    }


    void get_average_map(std::map<T, S>& result) const
    {
        typedef std::map<T, S> result_t;
	    typedef typename result_t::value_type result_value_t;

        result.clear();
        for (typename map_t::const_iterator iter = map_.begin(); iter != map_.end(); ++iter) 
            result.insert(result_value_t(iter->first, op.div_int(iter->second.first, iter->second.second)));
    }
};

// openmp lock
////////////////

#ifdef _OPENMP
struct omp_lock {
private:
    omp_lock_t l;
public:
    omp_lock() { omp_init_lock(&l); }
    ~omp_lock() { omp_destroy_lock(&l); }
    void lock() { omp_set_lock(&l); }
    void unlock() { omp_unset_lock(&l); }
};
#else
struct omp_lock {
    void lock() { }
    void unlock() { }
};
#endif

struct scope_lock {
private: 
    omp_lock& l;
public:
    explicit scope_lock(omp_lock& pl) : l(pl) { l.lock(); }
    ~scope_lock() { l.unlock(); }
};

// functions
///////////////////////////////////////////////////////////////////////////////

double round(double d, int places);

set<int> set_range(int len);

set<int> set_range(int min, int max);

vector<int> vector_range(int min, int max);

// inline & template functions
///////////////////////////////////////////////////////////////////////////////

inline int int_round(double d)
{
    return (int)(d > 0.0 ? (d + 0.5) : (d - 0.5));
}

inline ipoint2 int_round(const dpoint2& p)
{
    return ipoint2(int_round(p.x), int_round(p.y));
}

inline int mod(int a, int b) 
{
    if (a >= 0) return a % b; else return b*(1 - a/b) + a;
}

//inline int sqr(int x) { return x*x; }

  
template<class I> iipair center_of_mass(I begin, I end)
{
    iipair result(0, 0);
    int n = 0;

    while (begin != end) {
        result += *begin;
        ++n;
        ++begin;
    }
    if (n == 0) return result;
    return iipair(int_round((double)result.first/n), int_round((double)result.second/n));
}

// Returns true if the elements of one range are equal to the elements of
// the second range where equality is given by a binary predicate.
// If the predicate is not equality, the function works properly if
// for each two elements the sets of elements that are considered the same by
// the predicate are disjoint.
template<class I, class J, class P> 
    bool predicate_equal(I begin1, I end1, J begin2, J end2, P equal_pred)
{
    if (end1 - begin1 != end2 - begin2) return false;

    typename std::list< typename J::value_type > Jlist(begin2, end2);
    typename std::list< typename J::value_type >::iterator j;
    bool exists;

    for (I i = begin1; i != end1; ++i) {
        exists = false;
        j = Jlist.begin();
        while (!exists && j != Jlist.end()) {
            if (equal_pred(*i, *j)) { exists = true; Jlist.erase(j++); } else ++j;
        }
        if (!exists) return false;
    }
    return true;
}

// Returns true if the elements of one range are equal to the elements of
// the second range.
template<class I, class J> 
    bool ranges_equal(I begin1, I end1, J begin2, J end2)
{
    if (end1 - begin1 != end2 - begin2) return false;

    typename std::list< typename J::value_type > Jlist(begin2, end2);
    typename std::list< typename J::value_type >::iterator j;
    bool exists;

    for (I i = begin1; i != end1; ++i) {
        exists = false;
        j = Jlist.begin();
        while (!exists && j != Jlist.end()) {
            if (*i == *j) { exists = true; Jlist.erase(j++); } else ++j;
        }
        if (!exists) return false;
    }
    return true;
}

// Returns true if the elements of one range are subset of those of the second range.
template<class I, class J> 
    bool is_subset(I begin1, I end1, J begin2, J end2, bool strict = false)
{
    typename std::list< typename J::value_type > Jlist(begin2, end2);
    typename std::list< typename J::value_type >::iterator j;
    bool exists;

    for (I i = begin1; i != end1; ++i) {
        exists = false;
        j = Jlist.begin();
        while (!exists && j != Jlist.end()) {
            if (*i == *j) { exists = true; Jlist.erase(j++); } else ++j;
        }
        if (!exists) return false;
    }
    return (!strict || !Jlist.empty());
}

// diff = range2 \ range1
// returns true if range2 is subset of range1.   
template<class I, class J, class K> 
bool range_difference(I begin1, I end1, J begin2, J end2, K& diff)
{
    diff.assign(begin2, end2);
    typename std::list< typename J::value_type >::iterator j;
    bool exists, subset = true;

    for (I i = begin1; i != end1; ++i) {
        exists = false;
        j = diff.begin();
        while (!exists && j != diff.end()) {
            if (*i == *j) { exists = true; diff.erase(j++); } else ++j;
        }
        if (!exists && subset) subset = false;
    }
    return subset;
}


// Transforms p such that it represents the next permutation.
// Returns true if the next permutation exists and false otherwise.
template<class T> bool next_permutation(vector<T>& p)
{
    int j, k, l;
    int n = (int)p.size() - 1;
    int count = 0;
    T tmp;

    j = n - 1;
    while (j >= 0 && p[j] >= p[j + 1]) --j;
    if (j < 0) { // no next permutation
        return false;
    } else {
        l = n;
        while (p[j] >= p[l]) --l;
        tmp = p[j]; p[j] = p[l]; p[l] = tmp;
        k = j + 1; l = n;
        while (k < l) {
            tmp = p[k]; p[k] = p[l]; p[l] = tmp;
            ++k;
            --l;
        }
    }
    return true;
}

template<class T> void resize_vector(vector<T>& v, int size)
{
    vector<T> result(size);
    int vsize = (int)v.size();

    if (vsize > 0) {
        double f = (double)vsize/size;
        
        for (int i = 0; i < size; ++i) {
            int pos = (int)(i*f);
            if (pos < vsize) result[i] = v[pos]; else result[i] = v.back();
        }
    }
    v.swap(result);
}

template<class T> vector<T> get_resized_vector(const vector<T>& v, int size)
{
    vector<T> result = v;
    resize_vector(result, size);
    return result;
}

template<class T> void permute(vector<T>& v, const vector<int>& perm)
{
    if (perm.size() != v.size()) return;

    vector<T> pv(v.size());

    for (int i = 0; i < (int)perm.size(); ++i) {
        pv[i] = v[perm[i]];
    }
    v.swap(pv);
}

template<class T> void permute_and_resize(vector<T>& v, const vector<int>& perm)
{
    vector<T> pv(perm.size());

    for (int i = 0; i < (int)perm.size(); ++i) {
        pv[i] = v[perm[i]];
    }
    v.swap(pv);
}

template<class T> bool size_less(const vector<T>& v1, const vector<T>& v2) 
{ 
    return v1.size() < v2.size();
}


// functions
///////////////////////////////////////////////////////////////////////////////

void data2mat(char* filename, float* ptr, const int dimx, const int dimy);
void data2mat(char* filename, double* ptr, const int dimx, const int dimy);

void trim(string& str);

void end_dir(string& dir_name);

string change_extension(const string& fname, const string& ext, const string& extsep);

int scale_number(const string& fname);

string change_extension(const string& fname, const string& ext);

string get_extension(const string& fname, const string& extsep);

int hash_function(const string& s);

string get_prefix(string str, const string& sep, int n = 1);

string operator+(const string& s, int i);

string fill_left(int d, int width, char fillchar = '0');

void HSLtoRGB(const double h, const double s, const double l, double* result);

template <class In, class Out, class Pred> Out copy_if(In first, In last, Out res, Pred p)
{
	while(first != last)
	{
		if(p(*first)) *res++ = *first;
		++first;
	}
	return res;
}

// ad-hoc matrix functions
// inverse of a 3x3 matrix<T>
template <class T> void matrix3inv(const matrix<T>& A, matrix<T>& invA)
{
	invA = matrix<T>(3,3);
	T idet = +A(0,0)*(A(1,1)*A(2,2)-A(1,2)*A(2,1))
    -A(1,0)*(A(0,1)*A(2,2)-A(2,1)*A(0,2))
    +A(2,0)*(A(0,1)*A(1,2)-A(1,1)*A(0,2));
	idet = 1/idet;
	invA(0,0) =  (A(1,1)*A(2,2)-A(1,2)*A(2,1))*idet;
	invA(1,0) = -(A(1,0)*A(2,2)-A(2,0)*A(1,2))*idet;
	invA(2,0) =  (A(1,0)*A(2,1)-A(2,0)*A(1,1))*idet;
	invA(0,1) = -(A(0,1)*A(2,2)-A(2,1)*A(0,2))*idet;
	invA(1,1) =  (A(0,0)*A(2,2)-A(2,0)*A(0,2))*idet;
	invA(2,1) = -(A(0,0)*A(2,1)-A(0,1)*A(2,0))*idet;
	invA(0,2) =  (A(0,1)*A(1,2)-A(0,2)*A(1,1))*idet;
	invA(1,2) = -(A(0,0)*A(1,2)-A(0,2)*A(1,0))*idet;
	invA(2,2) =  (A(0,0)*A(1,1)-A(0,1)*A(1,0))*idet;
}

template <class Cont> void to_matlab_matrix(const Cont cont, const int dimx, ostream& os)
{
	int i = 0;
	for(typename Cont::const_iterator iter = cont.begin();
		iter != cont.end(); ++iter)
	{
		++i;
		os << *iter;
		if(i % dimx == 0)
		{
			os << endl;
			i = 0;
		}
		else
		{
			os << ' ';
		}
	}
}

// Eigenvalues and eigenvectors of a symmetric 2 x 2 matrix A
// A = [ a  b ]
//     [ b  d ]
void eigensystem(double& l1, double& l2, dpoint2& v1, dpoint2& v2, double a, double b, double d);


/*
template <class T> void from_matlab_matrix(const string& filename, vector<T>& v)
{
	ifstream ifs(filename.c_str());
	v.clear();
	while(!ifs.eof())
	{
		string s;
		getline(ifs, s);
		replace(s.begin(), s.end(), ',', ' ');
		istringstream iss(s);
		copy(istream_iterator<T>(iss), istream_iterator<T>(), back_inserter(v));
	}
}
*/
//unicode to normal strings in windows
std::string wstrtostr(const std::wstring &wstr);
std::wstring strtowstr(const std::string &str);
//void replace_substr(string& orig, const string& substr_old, const string& substr_new);
//void escape_for_perlregex(string& se);


#endif /* _MISC_H_ */
