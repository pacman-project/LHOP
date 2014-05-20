/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// layer_1_result

#pragma once
#ifndef _LAYER_1_RESULT_
#define _LAYER_1_RESULT_

#include <algorithm>
#include <utility>

#include "utils/img.h"
#include "utils/graphs/img_graph.h"
#include "utils/matrix.h"
#include "utils/misc.h"
#include "utils/utils.h"
#include "utils/config_dictionary.h"

#include "core/legacy/constants.h"

#include "core/legacy/part_lib.h"

// apparently windows uses _isnan and gnu uses isnan
#ifndef WIN32
#define _isnan ::isnan
#endif

// constants
//////////////

const unsigned NODE_REDUCED_ATTR = 1;
const unsigned HAS_NEXT_LAYER = 4;
const unsigned FROM_PREV_LAYER = 64;
const unsigned IN_NODE_ATTR = 0x1000;
const unsigned BOUDARY_NODE_ATTR = 0x2000;
const unsigned OUT_NODE_ATTR = 0x4000;
const unsigned HIDDEN_NODE_ATTR = 0x400000;
const unsigned ATTR_MARKED = 0x800000;
const unsigned ATTR_MARKED1 = 0x1000000;


void set_debug_output(bool val);
bool get_debug_output();

// typedefs
/////////////

typedef pair<double, double> ddpair;

typedef map<ipoint2, vector<double> > scmap_t;

// layer1_data
////////////////

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!! Important: change 'response_from_string' and 'response_type'   !!!
// !!! functions when adding or removing responses                    !!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
const int R_RESPONSE = 0;         // "response r"
const int G_RESPONSE = 2;         // response (g-response)
const int S_RESPONSE = 3;         // s-response
const int RR_RESPONSE = 5;        // Realization ratio #of covered basic parts/#of all basic parts
const int X_RESPONSE = 6;         // "test" response

int response_from_string(string str);

double reverse_s_response(double s, double maxs = 5.0);

struct response_map {
protected:
    static const double empty_val; // = NaN

    int size;
    vector<double> v;
public:
    response_map() : size(0), v(8, empty_val) { }
    
    bool response_defined(int name) { return name < (int)v.size() && !_isnan((double)v[name]); }

    void set_response(int name, double value) 
    { 
        if (name >= (int)v.size()) v.resize(name + 1, empty_val);
        if (_isnan(v[name])) ++size;
        v[name] = value;
    }

    double get_response(int name, double defval) const
    { 
        if (_isnan(v[name])) return defval; else return v[name]; 
    }

    double get_response(int name) const
    {
        if (!_isnan(v[name])) return v[name];
        else 
            throw new_libhop_exception(string("Response ") + name + string(" not defined."));
        
    }

    double operator()(int name) const { return get_response(name); }

    void write_to_stream(ostreamer& os);
    void read_from_stream(istreamer& is);
};

struct response_sort_f {
    int rname;

    response_sort_f(int r) : rname(r) { }
    bool operator()(node* n, node* m) const;
};

struct neg_response_sort_f {
    int rname;

    neg_response_sort_f(int r) : rname(r) { }
    bool operator()(node* n, node* m) const;
};


class layer1_data : public img_node_data {
public:
    response_map r;
    int m;

    layer1_data() {}
    layer1_data(double rval, int pm) : r(), m(pm) { r.set_response(R_RESPONSE, rval); }
    layer1_data(const response_map& pr, int pm) : r(pr), m(pm) { }
    layer1_data(const layer1_data& d) : r(d.r), m(d.m), img_node_data(d)  { }

	double val() const { return r(G_RESPONSE); }
    double vval() const { return z == 0 ? r(R_RESPONSE) : r(G_RESPONSE); }

    virtual void gc_clone(const node_data* d, const map<node*,node*>& cmap)
    {
        img_node_data::gc_clone(d, cmap);
        m = ((layer1_data*)d)->m;
        r = ((layer1_data*)d)->r;
    }

    virtual void copy_to(streamable* p, cloner& cl)
    {
        img_node_data::copy_to(p, cl);
        ((layer1_data*)p)->m = m;
        ((layer1_data*)p)->r = r;
    }

    virtual streamable* make_instance() const { return new layer1_data(); }
    virtual void read_from_stream(istreamer& is)
    {
        img_node_data::read_from_stream(is);

        if (is.get_version() >= 2.4) {
            is.read(m);
            r.read_from_stream(is);
        } else {
            double val;

            is.read(val);
            is.read(m);
            if (is.get_version() >= 2.1) r.read_from_stream(is);
            if (z == 0) {
                r.set_response(R_RESPONSE, val);
                r.set_response(G_RESPONSE, 1.0);
                r.set_response(RR_RESPONSE, 1.0);
                r.set_response(S_RESPONSE, 0.0);
            }
        } 
    }

    virtual void write_to_stream(ostreamer& os)
    {
        img_node_data::write_to_stream(os);
        os.write(m);
        r.write_to_stream(os);
    }

    virtual void print(ostream& os) const 
    { 
        os << m << " (" << val() << ')'; 
    }

    /// compares by value
    static bool less1(layer1_data& a, layer1_data& b)
    { 
        return a.val() < b.val();
    }

    /// compares nodes with associated layer1_data
    static bool less1n(node* a, node* b)
    { 
		if (((layer1_data*)a->data)->val() == ((layer1_data*)b->data)->val()) {
			return ((layer1_data*)a->data)->r(R_RESPONSE) + 10E-6 < ((layer1_data*)b->data)->r(R_RESPONSE);
		} else {
			return ((layer1_data*)a->data)->val() + 10E-6 < ((layer1_data*)b->data)->val();
		}
    }

    static bool greater1n(node* a, node* b)
    { 
		if (((layer1_data*)a->data)->val() == ((layer1_data*)b->data)->val()) {
			return ((layer1_data*)a->data)->r(R_RESPONSE) > ((layer1_data*)b->data)->r(R_RESPONSE);
		} else {
			return ((layer1_data*)a->data)->val() > ((layer1_data*)b->data)->val();
		}
    }

    static void sort_by_response(vector<node*>& v, int rname)
    {
        sort(v.begin(), v.end(), response_sort_f(rname));
    }
    
    static void sort_by_neg_response(vector<node*>& v, int rname)
    {
        sort(v.begin(), v.end(), neg_response_sort_f(rname));
    }

};

inline bool response_sort_f::operator()(node* n, node* m) const 
{
    layer1_data* nd = (layer1_data*)n->data;
    layer1_data* md = (layer1_data*)m->data;
    
    return nd->r(rname) > md->r(rname);
}

inline bool neg_response_sort_f::operator()(node* n, node* m) const 
{
    layer1_data* nd = (layer1_data*)n->data;
    layer1_data* md = (layer1_data*)m->data;
    
    return nd->r(rname) < md->r(rname);
}


// edge_data_
//////////////////

class edge_data_name : public edge_data {
public:
    int index;
	double r;

    edge_data_name(const edge_data_name& ed) : edge_data(), index(ed.index), r(ed.r) { }
    edge_data_name(int i, double pr) : index(i), r(pr) { }
	edge_data_name()  { } 

    virtual void clone(const edge_data_name* d) 
    { 
        index = ((edge_data_name*)d)->index;  
        r = ((edge_data_name*)d)->r;
    }
    
    virtual void copy_to(streamable* p, cloner& cl) 
    {
        edge_data::copy_to(p, cl);
        ((edge_data_name*)p)->index = index;
        ((edge_data_name*)p)->r = r;
    }

    virtual streamable* make_instance() const { return new edge_data_name(); }
    virtual void read_from_stream(istreamer& is) 
    { 
        edge_data::read_from_stream(is);
        is.read(index);
        if (is.get_version() >= 2.5) is.read(r); else r = 1.0;
    }

    virtual void write_to_stream(ostreamer& os) 
    { 
        edge_data::write_to_stream(os); 
        os.write(index);
        os.write(r);
    }

};


typedef edge_data_tw<edge_path> edge_path_data_t;


// layer_info
///////////////

struct layer_info {
    int covered;

    layer_info(int c = 0) : covered(c) { }
    void read_from_stream(istreamer& is) { is.read(covered); }
    void write_to_stream(ostreamer& os) { os.write(covered); }
};  

// layer1_result
//////////////////

const int VVE_SIMPLE_FILENAMES = 1;
const int VVE_TRUE_COORDINATES = 2;
const int VVE_BOUNDING_BOXES = 4;
const int VVE_MATLAB_SCALES = 8;

typedef pair<double, node*> dnpair;

struct sp_values {
    double g;
    double energy;
    node* n;

    sp_values() : g(0.0), energy(0.0), n(nullptr) { }
};

struct unary_double_power : public unary_double_function {
    double power;
    unary_double_power(double p = 1.0) : power(p) { }
    virtual double operator()(const double& x) const { return (power == 1.0) ? x : ::pow(x, power); }
};

typedef unary_function_t<layer1_data*, double> schur_product_function;

struct schur_product_params {
    double thresh;          // convolution thresh
    double dthresh;         // geometry checking thresh (distance)
    double scthresh;        // geometry checking thresh (shape context)

	// not used any more ?
    ipoint2 center;         // center (layer1) of 

    part_data_2* pdata;     // pointer to edge data in library

    schur_product_function* f;  

    schur_product_params(double t, double dt, double sct, const ipoint2& c, part_data_2* pd, schur_product_function* spf) 
        : thresh(t), dthresh(dt), scthresh(sct), center(c), pdata(pd), f(spf) { }
};

struct r_response_spf : public schur_product_function {
    virtual result_t operator()(const argument_t& nd) const { return nd->r(R_RESPONSE); }
};

struct v_response_spf : public schur_product_function {
    virtual result_t operator()(const argument_t& nd) const { return nd->val(); }
};

struct g_response_spf : public schur_product_function {
    double quotient;
    normal_distribution1 dist;

    g_response_spf(double q, const ddpair& d) : quotient(q), dist(d) { }
    virtual result_t operator()(const argument_t& nd) const 
    { 
        return dist.pdf_val1(::log(nd->r(R_RESPONSE)) - quotient)*nd->r(G_RESPONSE)/nd->r(RR_RESPONSE); 
    }
};

struct simple_g_response_spf : public schur_product_function {

    virtual result_t operator()(const argument_t& nd) const { 
        return nd->r(G_RESPONSE)/nd->r(RR_RESPONSE);
    }
};


struct identity_spf : public schur_product_function {
    virtual result_t operator()(const argument_t& nd) const { return 1.0; }
};

struct response_filter {
protected:
    typedef pair<double, bool> map_val_t; // bool: true: OK if val >= resp (e.g. R, RR, G), 
                                          // false: OK if val <= resp (e.g. S)
    typedef map<int, map_val_t> map_t;
    map_t rmap;
    
public:
    response_filter() : rmap() { }
    response_filter(int r, double v, bool less) : rmap() { rmap.insert(map_t::value_type(r, map_val_t(v, less))); }
    response_filter(int r1, double v1, bool less1, int r2, double v2, bool less2) : rmap() {
        rmap.insert(map_t::value_type(r1, map_val_t(v1, less1))); 
        rmap.insert(map_t::value_type(r2, map_val_t(v2, less2)));
    }

    void add_response(int r, double v, bool less) { rmap.insert(map_t::value_type(r, map_val_t(v, less))); }

    bool check(layer1_data* nd) const { 
        for (map_t::const_iterator iter = rmap.begin(); iter != rmap.end(); ++iter) {
            if (iter->second.second && nd->r(iter->first) < iter->second.first) return false;
            if (!iter->second.second && nd->r(iter->first) > iter->second.first) return false;
        }
        return true;
    }

    bool check(layer1_data* nd, part_lib* library) const { 
        if (library == 0) return check(nd);
        else {
            for (map_t::const_iterator iter = rmap.begin(); iter != rmap.end(); ++iter) {
                if (iter->second.second && nd->r(iter->first) < library->get_thresh(iter->first, nd->z, nd->m, iter->second.first))
                    return false;
                if (!iter->second.second && nd->r(iter->first) > library->get_thresh(iter->first, nd->z, nd->m, iter->second.first)) 
                    return false;
            }
            return true;
        }
    }
};


class layer1_result : public img_graph  {
public:
    // Schur product result type
    struct sp_result_data_t {
        node* n;
        int index;    // name (index of subpart)
        double r;     // response * map value

        sp_result_data_t(node* pn, int pindex, double pr) : n(pn), index(pindex), r(pr) { }
    };

    typedef vector<sp_result_data_t> sp_result_t;

    typedef img& (img::*combine_t)(const img& src, int desti, int destj, HOP_REAL factor);
    
    struct lpr_matrix_type {
        int m;
        lpr_matrix_type* next;

        lpr_matrix_type() : m(-1), next(nullptr) { }
        lpr_matrix_type(int pm, lpr_matrix_type* pnext) : m(pm), next(pnext) { }
        ~lpr_matrix_type() { if (next != nullptr) delete next; }

        bool find_m(int pm) 
        { 
            lpr_matrix_type* p = this;
            do {
                if (p->m == pm) return true;
                p = p->next;
            } while (p != nullptr);
            return false;
        }
        void add_m(int pm) 
        { 
            if (m == pm) return;
            else if (m < 0) m = pm; 
            else if (next != nullptr) next->add_m(pm);
            else next = new lpr_matrix_type(pm, nullptr);
        }
    };

    struct box_data_t {
        irectangle2 box;
        response_map r;
        double val;
        int size;
        int m;
        int nm; // neighbor's type

        box_data_t(const irectangle2& pbox, const response_map& pr, double pval, int psize, int pm, int pnm) 
            : box(pbox), r(pr), val(pval), size(psize), m(pm), nm(pnm) { }

        double getValue() const { return val; }

        struct box_data_f : public unary_function<box_data_t, irectangle2> {
            irectangle2& operator()(box_data_t& bd) const { return bd.box; }
        };

		bool operator< (const box_data_t& bd) const { return val < bd.val; }
		bool operator> (const box_data_t& bd) const { return val > bd.val; }
    };

    // Static functions -- "schur_product_function"s
    static r_response_spf rspf;
    static v_response_spf vspf;
    static simple_g_response_spf sgspf;
    static identity_spf idspf;


    int to_neighbor;              // = EdgeConnection::TO_NEIGHBOOR
    double layer1_region_threshold;    
    double layer1_threshold;            
    int layer1_3x3bound;                
    int layer1_neighb_radius;         
    int border;
    int original_width;

    unsigned attr;                // additional info: NOT streamable 

	// Parts sorted by layers. For each layer only one part for one location (x,y)
	// is present in vector, other parts in same location can be accessed by layer1_data->next.
    vector<vector<node*> > shape_nodes;
    vector<vector<node*> > shape_nodes_inhib;
    vector<layer_info> info;


    layer1_result() : img_graph(), shape_nodes(), shape_nodes_inhib(), info(), attr(0) { }
	virtual ~layer1_result();

	void get_layer_nodes(set<node*>& result, int z, const vector<int>& parts, double thresh, double thresh2, double thresh3, double thresh4, dpoint2* within_bounds = nullptr);
	void get_layer_nodes(vector<node*>& result, int z);
    void get_layer_nodes(vector<node*>& result, int z, double bestthresh);
	
	void delete_layers_geq(int z, bool set_has_next_layer_attr = true);
	
	void delete_nodes(const set<node*>& ns);
    void delete_nodes(unsigned attr);

    void inhibit(int z);
    void update_and_inhibit(int k);

    void add_results(const vector<layer1_result*>& results, int z);
    void add_results(layer1_result* res1, layer1_result* res2, int z);
    
	// used only by optimization
	void add_reconstruction_edges(int z); 
    void add_reconstruction_edges_link(int z); 
    void add_reconstruction_edges_fwd(int z);

    int max_layer_index();
    int layer_count() { return (int)shape_nodes.size(); }

	void get_parts(set<int>& parts, int z); // unused by libhop but used by vishop
	int get_nnodes(int z) { if(z < (int)shape_nodes.size() && z >= 0) return (int)shape_nodes[z].size(); return 0; }
    void synchronize_with_library(part_lib* library);
	
	void merge(layer1_result* res, int border = 0, double mfactor = 1.0);
    double get_q_quantile_val(int layer, double invq);

	irectangle2 get_box_with_cached_link(node* n);
	irectangle2 get_box(node* n);
	
    void get_boxes(list<box_data_t>& boxes, part_lib* library, int z, int response, bool use_lib_responses, bool revPCA, const response_filter& rfilter, double factor, int border0);
	void get_boxes(list<box_data_t>& boxes, int z, const set<int>& types, double inhibit, int resp_type = G_RESPONSE );
    
	void check_with_groundtruth(vector<double>& result, const list<irectangle2>& gtr, int z, const set<int>& types, double inhibit);
	
    double convolve_all(int x, int y, matrix<double>& m, int mt, double thresh, int z = 0);
    pair<node*, double> convolve_max_all(int x, int y, matrix<double>& m, int mt, double thresh, int z = 0);

	node* find_node_in_region(const point2<int>& p, int z, int deltax, int deltay, int m, unsigned attr, double valthresh = 0.0);
	
	double schur_product_max(int x, int y, const matrix<double>& m, int type, int z, schur_product_function* f, double thresh, vector<node*>& res);
	dnpair schur_product_max_all_new(sp_result_t& res, int xo, int yo, int z, double c, schur_product_params& spd);
    dnpair schur_product_max_all_ol(int xo, int yo, int dx, int dy, int index, matrix<double>& m, const map<int, pair<vector<int>, double> >& appmap, int z, schur_product_function* f, double thresh, sp_result_t& res);
    dnpair schur_product_max_all(int x, int y, matrix<double>& m, int type, int z, schur_product_function* f, double thresh);

	void get_nodes(vector<node*>& result, int z, const set<int>& mset);
	void get_contractions(vector<double>& contractions);

    void inc_covered(int z) { if (z < (int)info.size()) info[z].covered++; }
    double cover_quotient(int z) 
    { 
        if (z < (int)info.size() && z < (int)shape_nodes.size()) 
            return (double)info[z].covered/(double)shape_nodes[z].size();
        else
            return 0.0;
    }

    void set_attribute(unsigned a) { attr |= a; }
    void clear_attribute(unsigned a) { attr &= ~a; }
    bool is_attribute_set(unsigned a) { return (attr & a) != 0; }

	void type_statistics(vector<int>& stat, const vector<node*>& nodes, bool next); // never used

	///////////////////////////////////////
	/// Comperators
	bool data_less(img_node_data* d1, img_node_data* d2) 
    { 
		if (abs(((layer1_data*)d1)->val() - ((layer1_data*)d2)->val()) < 1e-6) {
            return ((layer1_data*)d1)->m > ((layer1_data*)d2)->m;
		} else {
			return ((layer1_data*)d1)->val() + 1E-6 < ((layer1_data*)d2)->val();
		}
    }

    bool data_greater(img_node_data* d1, img_node_data* d2) 
    { 
		if (abs(((layer1_data*)d1)->val() - ((layer1_data*)d2)->val()) < 1e-6) {
            return ((layer1_data*)d1)->m < ((layer1_data*)d2)->m;
		} else {
			return ((layer1_data*)d1)->val() > ((layer1_data*)d2)->val();
		}
    }
    bool data_equal(img_node_data* d1, img_node_data* d2) 
        { return ((layer1_data*)d1)->m == ((layer1_data*)d2)->m; }

	///////////////////////////////////////
	/// Visualization and export functions	
    img* get_image(set<node*>& nodes, int x_size, int y_size, bool paintuncov, const color& uncovered, const color& defcol,  const vector<int>& parts, const vector<color>& colors, bool bigpoints = false);
    img* get_image(int z, int zt, bool paintuncov, const color& uncovered, const color& defcol, const vector<int>& parts, const vector<color>& colors, bool bigpoints = false);
	
	img* get_image_inhib(int z, int zt, bool paintuncov, const color& uncovered, const color& defcol, const vector<int>& parts, const vector<color>& colors);
	
    img* get_image_boxed(part_lib* library, int z, const vector<int>& parts, const vector<color>& c, const color& defcol, bool inhibit);
    img* get_part_image(int z, int zt, node* p);
	
	img* get_image_reconstructed(part_lib* library, int z, int zt, const vector<int>& parts, 
        bool colored = false, bool drawbox = false, double factorpow = 1.0, bool all = false);
    img* get_image_reconstructed(part_lib* library, int z, const vector<int>& parts);
	
	// required by get_image_reconstructed to draw image
	// (overrided only by layer1_result_app and replaced by &img::combine_sum)
	combine_t get_combine_function() { return &img::combine_max; }
	
    void save_visview_filters(const string& dir);
	
    void save_visview(const string& dir, const string& name, const string& fname, part_lib* library, int z, int zmax, const map<node*,int>& prevlayermap, int save_mode);
    void save_visview(const string& dir, const string& nameint, const string& fname, part_lib* library, int maxlyr, int save_mode);
	
	void write_edge_info(ostream& os);
	
	// never used
	void debug_print_node_value_less(const double val, int layer);
		
	///////////////////////////////////////
	//// Streamable methods	(for saving/loading into memory/file)
    virtual void copy_to(streamable* p, cloner& cl);
    virtual streamable* make_instance() const { return new layer1_result(); }

    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);

protected:
    virtual void copy_to(graph* dest, const map<node*, node*>& cmap);
    void fill_shape_node_vectors();

};

// global stuff
///////////////////////////////////////////////////////////////////////////////

void read_layer1_result(layer1_result*& res, const string& fname);

layer1_result* read_layer1_result(const string& fname);

void save_layer1_result(layer1_result* res, const string& fname);

void save_node_set_mma(const string& fname, const set<node*>& nodes);

void keep_best_boxes(list<layer1_result::box_data_t>& boxes, int k);

void inhibit_layer(layer1_result* res, int z, int response, int maxn, double thresh);

// geometry/shape check/shape context functions
void link_path(layer1_result* res, node* n, int edge_name, int link_name);

void get_sc_map(scmap_t& scmap, layer1_result* res, const K_bin& bin, bool normalize);

node* follow_center_link(layer1_result* res, node* n);

int part_geometry_matching(double& benergy, vector<dpoint2>& dvector, double& scdistance, const path_map_t& im, const path_map_t& jm, bool calcbenergy);

int part_geometry_matching(double& benergy, vector<dpoint2>& dvector, double& scdistance, node* p, node* q, bool calcbenergy);

void part_geometry_matching_sym(double& benergy, vector<dpoint2>& dvector, double& scdistance, const path_map_t& im, const path_map_t& jm, bool calcbenergy);

ipoint2 get_path_map(path_map_t& pmap, layer1_result* res, const scmap_t& scmap, node* n, bool link);

void get_node_geo_p(vector<pair<int, ipoint2> >& ptsm, layer1_result* res, node* n);

void get_node_geo(vector<pair<int, ipoint2> >& ptsm, layer1_result* res, node* n);

vector<ipoint2> robust_PCA(layer1_result* res, node* n, node* p);

node* get_closest_node(layer1_result* res, int layer, int type, const dpoint2& p, int minr, int maxr, unsigned attr);

irectangle2 bounding_rectangle_of_projection(layer1_result* res, const vector<node*>& nodes);

void sample_tree(set<node*>& result, const set<node*>& ns);

void sample_tree(set<node*>& result, node* n);

void cluster_detections_ms(map<node*, node*>& clmap, layer1_result* res, 
    const list<node*>& nodes, double sigma, int steps);

void cluster_detections_ms(map<node*, vector<node*> >& clusters, layer1_result* res, 
    const list<node*>& nodes, double sigma, int steps);

void cluster_detections_ms(map<node*, vector<node*> >& clusters, layer1_result* res, int layer, const response_filter& rsf,
    double sigma, int steps);

double layer0_cover_ratio(layer1_result* res, int l);

double check_svmt(const svm_data& svmd, int nchildren, node* n, bool dfvalue);

// templates
///////////////////////////////////////////////////////////////////////////////

inline int node_type(node* n)
{
    return ((layer1_data*)n->data)->m;
}


inline int node_layer(node* n)
{
    return ((layer1_data*)n->data)->z;
}

inline ipoint2 node_coordinates(node* n) 
{
    layer1_data* nd = (layer1_data*)n->data;

    return ipoint2(nd->x, nd->y);
}

#endif /* _LAYER_1_RESULT_ */

