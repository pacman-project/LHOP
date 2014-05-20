/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// part_lib

#pragma once
#ifndef _PART_LIB_
#define _PART_LIB_

#include "utils/img.h"
#include "utils/graphs/img_graph.h"
#include "utils/matrix.h"
#include "utils/misc.h"
#include "utils/utils.h"
#include "utils/ocv.h"
#include "utils/hopmath.h"

#include <memory>

#if defined WIN32 | defined WIN64
#include <hash_map>
using namespace stdext;
#else
#include <ext/hash_map>
using namespace __gnu_cxx;
#endif

// forward declaration 
////////////////////////

class part_lib;

// constants
//////////////

const unsigned STRUCT_LIB_ATTR = 0;
const unsigned APP_LIB_ATTR = 1;


// typedefs
/////////////

typedef pair<double, double> ddpair;

// part_lib & data for nodes ... 
//////////////////////////////////

const unsigned CLUSTER_PART_ATTR = 0x100;
const unsigned HYPO_NODE_ATTR = 0x10000; // never set
const unsigned CATEGORY_NODE_ATTR = 0x20000;
const unsigned HYPO_NODE_CANDIDATE_ATTR = 0x80000;
const unsigned R_PART_ATTR = 0x100000;
const unsigned OBJ_PART_ATTR = 0x200000;
const unsigned VS_PART_ATTR = 0x400000; // never set

// thresh_data
////////////////

// list of thresholds
const int R_THRESH = 0;    // threshold for R_RESPONSE 
const int G_THRESH = 2;    // threshold for G_RESPONSE (layer1_data.val)
const int S_THRESH = 3;
const int RR_THRESH = 5;   // threshold for RR_RESPONSE (realization ratio)

struct thresh_data {
protected:
    typedef hash_map<int, double> container_t;

    container_t tmap;
public:
    thresh_data() : tmap() { }

    double get_thresh(int name, double defval = 0.0, double factor = 1.0) 
    { 
        container_t::iterator iter = tmap.find(name);
        if (iter == tmap.end()) return defval; else return factor*(iter->second);
    }

    void set_thresh(int name, double val) 
    {
        pair<container_t::iterator, bool> ibpair = tmap.insert(container_t::value_type(name, val));
        if (!ibpair.second) ibpair.first->second = val;
    }

    // Set min value of the existing thresh and val (or sets thresh to val if 
    // thresh does not exist)
    void update_thresh(int name, double val) 
    {
        pair<container_t::iterator, bool> ibpair = tmap.insert(container_t::value_type(name, val));
        if (!ibpair.second) ibpair.first->second = std::min<double>(val, ibpair.first->second);
    }

    bool empty() { return tmap.empty(); }

    void write_to_stream(ostreamer& os);
    void read_from_stream(istreamer& is);
};


// edge_path
//////////////

struct edge_path_item {
    int index; 
    int type; 

    edge_path_item() : index(0), type(0) { }
    edge_path_item(int i, int t) : index(i), type(t) { }
};
// this thing uses only edge_path_data_t which is in layer_1_result (i.e. layer structure not library structure)
struct edge_path {
    typedef vector<edge_path_item> path_t;

    path_t path;

    edge_path() : path() { }
    edge_path(const edge_path& p) : path(p.path) { }

    //const path_t operator()() { return path; }
    void push_back(int i, int t) { path.push_back(edge_path_item(i, t)); }
    void clear() { path.clear(); }
    int size() const { return (int)path.size(); }
    vector<int> edges() const
    { 
        vector<int> result;

        result.reserve(path.size());
        for (path_t::const_iterator piter = path.begin(); piter != path.end(); ++piter) 
            result.push_back(piter->index);
        return result;
    }
    void write_to_stream(ostreamer& os) 
    {
        os.write((int)path.size());
        os.write((char*)&path.at(0), sizeof(edge_path_item)*path.size());
    }
    void read_from_stream(istreamer& is)
    {
        int size;

        is.read(size);
        path.resize(size);
        is.read((char*)&path.at(0), sizeof(edge_path_item)*size);
    }
};

// lib_data
/////////////

struct lib_data : public node_data {
protected:
    int bpcount;  // number of "basic" parts; non-streamable, filled by get_basic_part_number()
public:
    int layer;
    int type;
    thresh_data td;

    lib_data(int l = 0, int t = -1) : bpcount(-1), layer(l), type(t), td() { }
    lib_data(const lib_data& ld) : bpcount(ld.bpcount), layer(ld.layer), type(ld.type), td(ld.td) { }

    double get_thresh(int name, double defval, double factor = 1.0) { return td.get_thresh(name, defval, factor); }

    virtual void copy_to(streamable* p, cloner& cl)
    {
        node_data::copy_to(p, cl);
        ((lib_data*)p)->layer = layer;
        ((lib_data*)p)->type = type;
        ((lib_data*)p)->td = td;
    }

    virtual streamable* make_instance() const { return new lib_data(); }
    virtual void read_from_stream(istreamer& is) 
    { 
        node_data::read_from_stream(is); 
        is.read(type); is.read(layer); 
        td.read_from_stream(is);
    }
    virtual void write_to_stream(ostreamer& os) 
    { 
        node_data::write_to_stream(os); 
        os.write(type); os.write(layer); 
        td.write_to_stream(os);
    }
    virtual img get_image(node* n, unsigned libtype, part_lib* library) { return img(); }

protected:
    int get_basic_part_number(node* n); 

    friend class part_lib;
};

// part_data
//////////////

typedef vector<double> histogram_t;

struct hpoint_t {
    ipoint2 p;
    histogram_t h;
};

struct hdpoint_t {
    dpoint2 p;
    histogram_t h;
};


typedef map<vector<int>, hpoint_t> path_map_t;

struct path_map_stat {
protected:
    typedef map<vector<int>, int> count_map_t;

    path_map_t pm;
    count_map_t cm;
public:
    path_map_stat() : pm(), cm() { }

    void add(const path_map_t& m); 
    path_map_t get() const;
    void get(path_map_t& m) const;
};

struct part_data : public lib_data {
protected:
    typedef map<ipoint2, double> mask_t;
    typedef set<ipoint2> region_t;

    static vector<double> contractions;
public:
	// this should be renamed into reconstruction points/image
	// 
	// mask == original filter used to construct the part
    mask_t mask;  // Non-serializable, except for layer 0, filled recursively by get_mask

	// this should be renamed into activation points i.e. points of the image where this part fired from (was activated from) 
	// (for first layer parts this can be area of the filter with locations of the non-zero values)
	// this is used only int the ly1 support to calculate intersection between two parts
	// region == mask of the region filter i.e. locations of the region filter values with the positive factor
    region_t region; // Non-serializable, except for layer 0, filled recursively by get_region

    int cmx, cmy; // center of mass relative to the central node

    part_data() : mask(), region(), cmx(0), cmy(0), lib_data(0, -1) {}
    part_data(const img& m, const matrix<bool>& r, int l, int t = -1); 
    part_data(const part_data& pd);

    virtual void copy_to(streamable* p, cloner& cl)
    {
        lib_data::copy_to(p, cl);
        ((part_data*)p)->cmx = cmx;
        ((part_data*)p)->cmy = cmy;
        ((part_data*)p)->mask = mask;
        ((part_data*)p)->region = region;
    }
    virtual streamable* make_instance() const { return new part_data(); }
    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);

    virtual img get_image(node* n, unsigned libtype, part_lib* library); // vizualzation purpuse only

    ipoint2 get_mask(matrix<double>& m, node* n, part_lib* library); // used by vizualzation and similarity
    irectangle2 get_mask_box(node* n, part_lib* library); // used by vizualzation and get_mask (therefor also used by similarity)

	ipoint2 get_region_matrix(matrix<bool>& m, node* n);
    void get_region_set(set<ipoint2>& s, node* n);
    void reset_mask();
    void reset_region();

    static void set_contractions(const vector<double>& cv);
    static void set_contractions(double carr[], int size);
	static double total_contraction(int layer);

protected:
    void make_mask(node* n, unsigned libtype, part_lib* library);
    virtual void merge_mask_maps(map<ipoint2, double>& result, const map<ipoint2, double>& src, 
        const ipoint2& delta, const double& factor, unsigned libtype);
    void make_region(node* n);
};

// vs_part_data
/////////////////

struct pca_data {
    cv::Mat mean;
    cv::Mat eigenvectors;
    cv::Mat eigenvalues;
    double sizefactor;

    pca_data() : mean(), eigenvectors(), eigenvalues(), sizefactor(1.0) { }
    pca_data(const cv::Mat& m, const cv::Mat& eve, const cv::Mat& eva, double sf) 
        : mean(m), eigenvectors(eve), eigenvalues(eva), sizefactor(sf) { }
    pca_data(const pca_data& d) 
        : mean(d.mean), eigenvectors(d.eigenvectors), eigenvalues(d.eigenvalues), sizefactor(d.sizefactor) { }

    pca_data& operator=(const pca_data& d) 
    { 
        mean = d.mean; eigenvectors = d.eigenvectors; eigenvalues = d.eigenvalues; 
        sizefactor = d.sizefactor;
        return *this; 
    }

    void write_to_stream(ostreamer& os) const;
    void read_from_stream(istreamer& is);
};

struct svm_data {
    shared_ptr<cv::SVM> svm;
    cv::Mat mean, sigma;

    svm_data() : svm(nullptr), mean(), sigma() { }
    svm_data(cv::SVM* svmp) : svm(svmp), mean(), sigma() { }
    svm_data(cv::SVM* svmp, const cv::Mat& m, const cv::Mat& s) : svm(svmp), mean(m), sigma(s) { }

    void write_to_stream(ostreamer& os) const;
    void read_from_stream(istreamer& is);
};

class vs_part_data : public part_data, public cv_storable_data {
public:
    pca_data pcad;
    svm_data svmd;  // shape (unused)
    svm_data svmt;  // thresholds
    svm_data svmh;  // "hoc"
public:
    vs_part_data() : pcad(), svmd(), svmt(), svmh() { }
    vs_part_data(const part_data& pd) : part_data(pd), pcad(), svmd(), svmt(), svmh() { }

    void set_data(const pca_data& d) { pcad = d; }
    string name() const;

    virtual void copy_to(streamable* p, cloner& cl)
    {
        part_data::copy_to(p, cl);
        ((vs_part_data*)p)->pcad = pcad;
        ((vs_part_data*)p)->svmd = svmd;
        ((vs_part_data*)p)->svmt = svmt;
        ((vs_part_data*)p)->svmh = svmh;
    }

    virtual streamable* make_instance() const { return new vs_part_data(); }
    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);
    virtual void write_to_storage(CvFileStorage* fs);
    virtual void read_from_storage(CvFileStorage* fs);
};

class pca_learning  : public streamable {
protected:
    typedef vector<ipoint2> point_stat_t;
    typedef list<point_stat_t> samples_t;

    samples_t samples;
public:
    pca_learning() : samples() { }
    
    void update(const vector<ipoint2>& pts) 
    {
        samples.push_back(pts);
    }

	void merge(const pca_learning& pcal)
	{
		for (samples_t::const_iterator iter = pcal.samples.begin(); iter != pcal.samples.end(); ++iter) {
			samples.push_back(*iter);
		}
	}
	void print() {
	}
    pca_data get_pca_data(int maxcomp = 5) const;

	virtual streamable* make_instance() const { return new pca_learning(); }
	virtual void read_from_stream(istreamer& is) {
		is.read(samples);
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write(samples);
	}
};

// cpart_data ("class part data")
///////////////////////////////////////////////////////////////////////////////

class cpart_data : public part_data {
public:
    string name;

    cpart_data(const string& n = "") : part_data(), name(n) { }

    virtual void copy_to(streamable* p, cloner& cl) { part_data::copy_to(p, cl); ((cpart_data*)p)->name = name; }
    virtual streamable* make_instance() const { return new cpart_data(); }
    virtual void read_from_stream(istreamer& is) 
    { 
        part_data::read_from_stream(is); 
        is.read(name); 
    }

    virtual void write_to_stream(ostreamer& os) 
    { 
        part_data::write_to_stream(os); 
        os.write(name); 
    }

};

// part_data_2a
///////////////////

class part_data_2a : public edge_data {
public:
    typedef map<int, pair<vector<int>,  double> > app_map_t;	// type -> (newindex, 1.0)
    typedef map<int, path_map_t> geo_map_t;

    app_map_t app; // appearance
    geo_map_t geo; // Geometric "average" of the subpart 
    
    part_data_2a() : app(), geo() { }
    part_data_2a(int m) : app(), geo() { app.insert(app_map_t::value_type(m, app_map_t::mapped_type(vector<int>(), 1.0))); }
    part_data_2a(const app_map_t& amap) : app(amap), geo() { }
    part_data_2a(const part_data_2a& pd) : app(pd.app), geo(pd.geo) { }

    void set_geo(const path_map_t& pm)
    {
        if (app.empty()) {
            cout << "Trying to update geometry of edge with empty 'app'." << endl;
            throw;
        }
        geo[app.begin()->first] = pm;
    }

    bool operator==(const part_data_2a& d) { return app == d.app; }
    
    virtual void copy_to(streamable* p, cloner& cl)
    {
        edge_data::copy_to(p, cl);
        ((part_data_2a*)p)->app = app;
        ((part_data_2a*)p)->geo = geo;
    }

    virtual streamable* make_instance() const { return new part_data_2a(); }
    virtual void read_from_stream(istreamer& is); 
    virtual void write_to_stream(ostreamer& os); 
};

// part_data_2
////////////////

class part_data_2 : public part_data_2a {
public:
    int index;      // index of edge -- used in edge_path
    int x, y;       // position relative to center node (center node has position (0, 0))
    matrix<double> distr;   // distribution
    ddpair gdistr;  // mean and variance of the distribution for G_RESPONSE
    thresh_data td; // thresholds (if set)

    part_data_2() : part_data_2a(), index(-1), x(0), y(0), distr(), gdistr(0.0, 0.0), td() { }
    part_data_2(int vx, int vy, const matrix<double>& d, int m, int i) 
        : part_data_2a(m), index(i), x(vx), y(vy), distr(d), gdistr(0.0, 0.0), td() { }
    part_data_2(const iipair& p, const matrix<double>& d, int m, int i) 
        : part_data_2a(m), index(i), x(p.first), y(p.second), distr(d), gdistr(0.0, 0.0), td() { }
    part_data_2(const part_data_2& pd) : part_data_2a(pd), index(pd.index), x(pd.x), y(pd.y), distr(pd.distr), gdistr(0.0, 0.0), td() { }

    virtual void copy_to(streamable* p, cloner& cl)
    {
        edge_data::copy_to(p, cl);
        ((part_data_2*)p)->index = index;
        ((part_data_2*)p)->x = x;
        ((part_data_2*)p)->y = y;
        ((part_data_2*)p)->distr = distr;
        ((part_data_2*)p)->gdistr = gdistr;
        ((part_data_2*)p)->app = app;
        ((part_data_2*)p)->td = td;
        ((part_data_2*)p)->geo = geo;
    }

    virtual streamable* make_instance() const { return new part_data_2(); }
    
    virtual void read_from_stream(istreamer& is) 
	{ 
        part_data_2a::read_from_stream(is); 
        is.read(index);
        is.read(x); is.read(y); 
        is.read(distr); 
        is.read(gdistr);
		if (is.get_version() <= 3.3) {	// "Bug" in previous versions; it is serialized already in part_data_2a
			map<int, double> app1;

			is.read(app1);	
		}
        td.read_from_stream(is);
    }

    virtual void write_to_stream(ostreamer& os) 
	{ 
        part_data_2a::write_to_stream(os); 
        os.write(index);
        os.write(x); os.write(y); 
        os.write(distr); 
        os.write(gdistr); 
        td.write_to_stream(os);
    }

};


// part_data_2c
/////////////////

class part_data_2c : public edge_data {
public:
    int x, y;
	
    part_data_2c(int vx = 0, int vy = 0) : x(vx), y(vy) { }
    part_data_2c(const iipair& p) : x(p.first), y(p.second) { }
    part_data_2c(const ipoint2& p) : x(p.x), y(p.y) { }
    part_data_2c(const part_data_2c& pd) : x(pd.x), y(pd.y) { }

    bool operator==(const part_data_2c& d) { return x == d.x && y == d.y; }
    
    virtual void copy_to(streamable* p, cloner& cl)
    {
        edge_data::copy_to(p, cl);
        ((part_data_2c*)p)->x = x;
        ((part_data_2c*)p)->y = y;
    }

    virtual streamable* make_instance() const { return new part_data_2c(); }
    virtual void read_from_stream(istreamer& is) 
        { edge_data::read_from_stream(is); is.read(x); is.read(y); }
    virtual void write_to_stream(ostreamer& os) 
        { edge_data::write_to_stream(os); os.write(x); os.write(y); }
};

// part_data_sim
//////////////////

class part_data_sim : public edge_data {
public:
    map<vector<int>, int> perm;
    double val;

    part_data_sim(double v = 0.0) : perm(), val(v) { }
    part_data_sim(const map<vector<int>, int>& p, double v = 0.0) : perm(p), val(v) { }

    virtual void copy_to(streamable* p, cloner& cl)
    {
        edge_data::copy_to(p, cl); 
        ((part_data_sim*)p)->perm = perm;
        ((part_data_sim*)p)->val = val;
    }

    virtual streamable* make_instance() const { return new part_data_sim(); }
    virtual void read_from_stream(istreamer& is) 
        { edge_data::read_from_stream(is); is.read(perm); is.read(val); }
    virtual void write_to_stream(ostreamer& os) 
        { edge_data::write_to_stream(os); os.write(perm); os.write(val); }
};

// part_data_str
//////////////////

struct part_str {
    int ldelta, type;
    int x, y;
    part_data_2* data; // optional

    part_str(int ld, int t, int xc, int yc, part_data_2* d) : 
        ldelta(ld), type(t), x(xc), y(yc), data(d) { }
    part_str(const part_str& p) : 
        ldelta(p.ldelta), type(p.type), x(p.x), y(p.y), data(p.data) { }

    bool operator==(const part_str& p) const { 
        return (ldelta == p.ldelta && type == p.type && x == p.x && y == p.y); 
    }
};

// part_lib
/////////////

class part_lib : public graph {
protected:
    bool updated;
    unsigned attr;

public:
    typedef vector<node*>::iterator niterator;
    typedef pair<irectangle2, set<int> > rs_pair_t;

    // Data for description of object parts
    struct cluster_data_t { 
        vector<ipoint2> types;      // {(lyr. diff., type), ...}
        ipoint2 pos;                // position
        map<int, path_map_t> geo;   // geo 

        bool operator<(const cluster_data_t& cd) const { return types < cd.types; }
    };


    vector<vector<node*> > parts;
    double contractions[MAX_LAYER_NUMBER];
    streamable* layer_data[MAX_LAYER_NUMBER];
    int layer_count;
    static double similarity_threshold;


    part_lib(unsigned libtype = 0);
    virtual ~part_lib();

    unsigned get_attr() { return attr; }
   
	
	/// Returns the number of parts
    int layer_size(int index = 0) { if (index < layer_count) return (int)parts[index].size(); else return 0; }

    /// Returns the number of layers in the library.
    int max_layer() { return layer_count; }

    int max_layer_index() { return layer_count - 1; }


	bool is_subpart(const vector<part_str>& left, const vector<part_str>& right, list<part_str>& diff);
	
	int part_exists(const vector<node*>& lyr, const vector<iipair>& str, const vector<iipair>& coo);
    
    int find_equivalent_object_part(int lyr, const vector<cluster_data_t>& objdata, int masksize);
    int find_equivalent_object_part(int lyr, const vector<pair<vector<ipoint2>, ipoint2> >& objdata, int masksize);
    
    int add_part(int lyr, part_data* pd, const vector<iipair>& str, const vector<iipair>& coo,
        const vector<matrix<double>*>& distr, double eqthresh = -1.0, double contraction = 1.0, int tol = 0, int eqtol = 0);
    int add_object_part(int lyr, part_data* pd, const vector<cluster_data_t>& objdata,
        const vector<matrix<double>*>& distr, double contraction = 1.0);
    int add_category(int lyr, const vector<int>& members, const string& name);
    string get_category(node* p);
    
    void delete_parts(int lyr, vector<int> v);
    void delete_parts(int lyr);
    void delete_parts_geq(int lyr);

	void delete_unused_parts();

    void keep_parts(int lyr, const vector<int>& v);
    void keep_clusters(int lyr, const set<int>& v);
    
    int make_similarity_clusters(int lyr, double thresh);
    void add_similarity_edges(int layer, double contraction);
    void get_similar_types(int tlayer, int t, int clayer, int ct, int x, int y, double typethresh, set<int>& appmap);
    double similarity(node* part1, node* part2, double rfactor, int depth, double M(double[], int, int));
    static double defaultMfunction(double* arr, int n, int lyr);
	
    pair<node*, part_data_2*> get_neighbor_pair(node* p, int index);
    
    double contraction(int l1, int l2);
    
    bool is_vs_layer(int layer);
protected:

    int add_object_part(int lyr, part_data* pd, const vector<pair<vector<ipoint2>, ipoint2> >& objdata,
        const vector<matrix<double>*>& distr, double contraction = 1.0);
    void extend_subparts(int lyr, node* n, niterator begin, niterator end);
    void get_structure(node* n, vector<part_str>& str);
    void get_structure(node* n, int name, vector<part_str>& str);
public:
    void get_regions(int lyr, vector<matrix<bool> >& regions, vector<ipoint2>& region_centers);
    void get_regions(int lyr, vector<set<ipoint2> >& regions);

	void init_basic_part_numbers(int layer);

	int get_basic_part_number(node* p); 

    vector<int> get_root_parts_map(int layer);
	
	void update_var_fields();

    void edge_paths(vector<edge_path>& paths, node* p);
	
    double get_thresh(int name, int lyr, int type, double defval);

	///////////////////////////////////////
	//// Visualization stuff
    img lib_image(int lyr, int min, int max, bool show_labels, bool make_border, 
        bool mark_center, bool one_row);

    void save_all(const char* fname, int lyr, int min = 0, int max = -1, bool show_labels = true, bool make_border = false, 
        bool mark_center = false, bool one_row = false);
    void save_all_sc(const char* fname, int lyr, int min = 0, int max = -1, bool show_labels = true, bool one_row = false,
        bool save_variations = true);
    void save_sc(const string& outdir, bool show_labels, bool sc, bool save_variations = true);
    void save_all_sc_mma(const char* fname, int lyr, int min = 0, int max = -1);
    void save_all(const char* fname, int lyr, const vector<int>& l, bool show_labels, bool mark_center = false);

    void save_part(const char* fname, int lyr, int part, double threshold = 0.0);

    void save_visview_filters(const string& dir);
    void save_visview_centers(const string& dir, int lyr);
    void save_visview_parts(const string& dir, int lyr);
    void save_visview_contraction(const string& dir, int lyr);
    void save_visview_centers(const string& dir);
    void save_visview_parts(const string& dir);
    void save_visview_contractions(const string& dir);

    void drop_images();
	

	void display_layer_info(int layer);

	///////////////////////////////////////
	//// Streamable methods	(for saving/loading into memory/file)
	static part_lib* read(const string& name);
    virtual void save(const string& file, int zlib_compression_level = Z_NO_COMPRESSION);

	virtual void copy_to(streamable* p, cloner& cl);
    virtual streamable* make_instance() const { return new part_lib(0); }
    
	virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);

    void read_cv_storable_data(CvFileStorage* fs);
    void write_cv_storable_data(CvFileStorage* fs);
};

// global stuff
/////////////////

// geometry/shape check/shape context functions
void add_to_path_map(path_map_t& pm, const path_map_t& pm2);

path_map_t reduce_path_map(const path_map_t& pm, int radius);

path_map_t synchronize_path_map(const path_map_t& pm, const path_map_t& srcpm);

void path_map_union(path_map_t& pm, const path_map_t& pm1, int ext);

void get_library_geo(path_map_t& pm, node* pn);
void get_library_geo(path_map_t& geo, node* pn, const map<int, int>& imap);

vector<ipoint2> get_library_geo(node* pn);

vector<pair<int, ipoint2> > get_library_geo_pieces(node* pn, int radius);

void get_dissimilar_cluster(vector<node*>& cluster, const vector<node*>& parts, int maxsize);

vector<ipoint2> get_path_map_points(const path_map_t& pm);

void inhibit_path_map(path_map_t& pm, int maxsize);

bool is_sim_root(node* p);

vector<dpoint2> random_sample_from_part(node* p);

vector<pair<dpoint2, int> >  random_sample_part(part_lib* library, node* n);


cv::Mat back_projection(const cv::Mat& data, const pca_data& pcd);
vector<dpoint2> sample_from_pca(const pca_data& pd);


#endif /* _PART_LIB_ */
