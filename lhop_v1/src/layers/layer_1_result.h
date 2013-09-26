/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// layer_1_result

#pragma once
#ifndef _LAYER_1_RESULT_
#define _LAYER_1_RESULT_

//#ifdef LAYERS_EXPORTS
//#define LAYERS_API __declspec(dllexport)
//#else
//#define LAYERS_API __declspec(dllimport)
//#endif

// #include <ctime>
#include <algorithm>
#include <utility>

#include "../utils/img.h"
#include "../graphs/img_graph.h"
#include "../utils/matrix.h"
#include "../utils/misc.h"
#include "../utils/utils.h"

#include "part_lib.h"

// apparently windows uses _isnan and gnu uses isnan
#ifndef WIN32
#define _isnan ::isnan
#endif

// constants
//////////////

// !! NOTE : DO NOT forget to add them to opencl source when added new attribute key (values must in both definitions must match)
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
int response_type(int r);
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
#ifdef OPENCL
	ocl_layer1_data* ocl_struct;
#else
	void* ocl_struct;
#endif

    layer1_data() { }
    layer1_data(double rval, int pm) : r(), m(pm), ocl_struct(nullptr) { r.set_response(R_RESPONSE, rval); }
    layer1_data(const response_map& pr, int pm) : r(pr), m(pm), ocl_struct(nullptr) { }
    layer1_data(const layer1_data& d) : r(d.r), m(d.m), img_node_data(d), ocl_struct(nullptr)  { }

#ifdef OPENCL
	layer1_data(ocl_layer1_data* ocl_struct) {
		
		this->ocl_struct = ocl_struct;

		this->m = ocl_struct->m;

		this->r.set_response(R_RESPONSE, ocl_struct->response.R);
		this->r.set_response(RR_RESPONSE, ocl_struct->response.RR);
		this->r.set_response(G_RESPONSE, ocl_struct->response.G);
	}
#endif
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
        //if (((layer1_data*)a->data)->z == 0) return ((layer1_data*)a->data)->r(R_RESPONSE) < ((layer1_data*)b->data)->r(R_RESPONSE);
        //else return ((layer1_data*)a->data)->val() < ((layer1_data*)b->data)->val();
		if (((layer1_data*)a->data)->val() == ((layer1_data*)b->data)->val()) {
			return ((layer1_data*)a->data)->r(R_RESPONSE) + 10E-6 < ((layer1_data*)b->data)->r(R_RESPONSE);
		} else {
			return ((layer1_data*)a->data)->val() + 10E-6 < ((layer1_data*)b->data)->val();
		}
    }

    static bool greater1n(node* a, node* b)
    { 
        //if (((layer1_data*)a->data)->z == 0) return ((layer1_data*)a->data)->r(R_RESPONSE) > ((layer1_data*)b->data)->r(R_RESPONSE);
        //else return ((layer1_data*)a->data)->val() > ((layer1_data*)b->data)->val();
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

#ifdef OPENCL
	ocl_edge_data_ip2* ocl_struct;
#else
	void* ocl_struct;
#endif
    edge_data_name(const edge_data_name& ed) : edge_data(), index(ed.index), r(ed.r) { }
    edge_data_name(int i, double pr) : index(i), r(pr), ocl_struct(nullptr) { }
	edge_data_name() : ocl_struct(nullptr) { } 

#ifdef OPENCL
	edge_data_name(ocl_edge_data_ip2* ocl_struct) {
        index = ocl_struct->x;
		r = ocl_struct->r;
	}
#endif
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

    //bool check(int r, double v, bool neval = true) const {
    //    map_t::const_iterator iter = rmap.find(r);
    //    if (iter == rmap.end()) return neval; 
    //    else if (less) return v >= iter->second;
    //    else return v <= iter->second;
    //}

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

        double get_value() const { return val; }

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


    int to_neighbor;              // = atom("toNeighbor").index  
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

	CStopWatch clock_make_data_from_ocl;
#ifdef OPENCL
	// opencl memory objects
	vector<pair<ocl_layer1_data*, int> > ocl_shape_nodes;
	vector<pair<ocl_edge_data_ip2*, int> > ocl_edges;

	// Position (offset) of parts in memory object based on implicit
	// location. Size of each array for each layer must be of size (img_height * img_width) where
	// img_height and img_width can contract for each layer diffrently.
	vector<pair<ocl_layer1_data_coordinates*, int> > ocl_shape_nodes_coord;
	// number of non zero positions per layer
	vector<int> ocl_shape_nodes_coord_non_zero_count;
	// List of inhibited nodes for each location (only positions of parts in memory object for same layer).
	// There can be only ONE part in each position so size should only be 0 or 1.
	vector<pair<ocl_layer1_data_coordinates*, int> > ocl_shape_nodes_inhib_coord;
	// number of non zero positions per layer
	vector<int> ocl_shape_nodes_inhib_coord_non_zero_count;
#endif

    layer1_result() : img_graph(), shape_nodes(), shape_nodes_inhib(), info(), attr(0) { }
	virtual ~layer1_result();

    img* get_image(set<node*>& nodes, int x_size, int y_size,
        bool paintuncov, const color& uncovered, const color& defcol, 
        const vector<int>& parts, const vector<color>& colors, bool bigpoints = false);
    img* get_image(int z, int zt, bool paintuncov, const color& uncovered, const color& defcol, 
        const vector<int>& parts, const vector<color>& colors, bool bigpoints = false);
	img* get_image_new(int z, const vector<int>& parts, const color& defcol, const vector<color>& colors, bool bigpoints = false);
    img* get_image_inhib(int z, int zt, bool paintuncov, const color& uncovered, const color& defcol, 
        const vector<int>& parts, const vector<color>& colors);
    img* get_image_inhib_simple(int z);
    img* get_image_inhib_simple_01(int z);
    img* get_image_boxed(int z, const vector<int>& parts, const vector<color>& c, const color& defcol, bool inhibit);
    img* get_part_image(int z, int zt, node* p);
    img* get_part_reconstructed(int z, int zt, node* p);
    img* get_image_reconstructed(int z, int zt, const vector<int>& parts, 
        bool colored = false, bool drawbox = false, double factorpow = 1.0, bool all = false);
    img* get_image_reconstructed(part_lib* library, int z, const vector<int>& parts);
    void get_reconstruction_nodes(set<node*>& result, int z, const vector<int>& parts);
    virtual combine_t get_combine_function() { return &img::combine_max; }
    template<class Man> int get_reconstruction_nodes(Man& ex, int z, const vector<int>& parts,
        double thresh, double thresh2, double thresh3, int index);
	template<class Man> void get_hypo_node_count(Man& man,
		node* n);
    void part_reconstruction(map<int, img>& result, int z);
    void get_layer_nodes(set<node*>& result, int z, const vector<int>& parts, 
        double thresh, double thresh2, double thresh3, double thresh4, dpoint2* within_bounds = nullptr);
	void get_layer_nodes(map<int,set<node*> >& nodes, int z, const vector<int>& parts, 
        double thresh, double thresh2, double thresh3, double thresh4, dpoint2* within_bounds = nullptr);
	void get_layer_nodes(vector<node*>& result, int z);
    void get_layer_nodes(vector<node*>& result, int z, double bestthresh);
    void distance_clustering(vector<vector<ipoint2> >& clusters, int layer, int cluster_n);

    void save_visview_filters(const string& dir);

    void dilute_layer(int z, const vector<int>& parts);
    void delete_layers_geq(int z, bool set_has_next_layer_attr = true);
    void remove_grid_nodes(int x, int y, int z, unsigned attr, double valthresh);
    void filter_nodes(int z, const irectangle2& rect);
    virtual void delete_nodes(const set<node*>& ns);
    virtual void delete_nodes(unsigned attr);

    void inhibit(int z);
    void update_and_inhibit(int k);
    virtual void add_results(const vector<layer1_result*>& results, int z);
    virtual void add_results(layer1_result* res1, layer1_result* res2, int z);
    void add_reconstruction_edges(int z); 
    void add_reconstruction_edges_link(int z); 
    void add_reconstruction_edges_back(int z); 
    void add_reconstruction_edges_fwd(int z);
    void add_reconstruction_edges_fwd_link(int z);
    void add_reconstruction_edges_leq_fwd(int z);

    void adjust_rectangle(irectangle2& rect, int srclayer, int destlayer);
    double size_factor(int srclayer, int destlayer);
    void get_box(irectangle2& box, int z);
    ipoint2 get_center_of_mass(int z);

    int max_layer_index();
    int layer_count() { return (int)shape_nodes.size(); }
    int get_max_part_index(int z);
    void get_parts(set<int>& parts, int z);
	int get_nnodes(int z) { if(z < (int)shape_nodes.size() && z >= 0) return (int)shape_nodes[z].size(); return 0; }
    void synchronize_with_library(part_lib* library);
	void shift_part_id(layer1_result* res, const int layer, const int offset);
    void merge(layer1_result* res, int border = 0, double mfactor = 1.0);
    double get_q_quantile_val(int layer, double invq);

	virtual bool data_less(img_node_data* d1, img_node_data* d2) 
    { 
		if (abs(((layer1_data*)d1)->val() - ((layer1_data*)d2)->val()) < 1e-6) {
			//return ((layer1_data*)d1)->r(R_RESPONSE) + 1E-6 < ((layer1_data*)d2)->r(R_RESPONSE);
            return ((layer1_data*)d1)->m > ((layer1_data*)d2)->m;
		} else {
			return ((layer1_data*)d1)->val() + 1E-6 < ((layer1_data*)d2)->val();
		}
    }

    virtual bool data_greater(img_node_data* d1, img_node_data* d2) 
    { 
		if (abs(((layer1_data*)d1)->val() - ((layer1_data*)d2)->val()) < 1e-6) {
			//return ((layer1_data*)d1)->r(R_RESPONSE) > ((layer1_data*)d2)->r(R_RESPONSE);
            return ((layer1_data*)d1)->m < ((layer1_data*)d2)->m;
		} else {
			return ((layer1_data*)d1)->val() > ((layer1_data*)d2)->val();
		}
    }
    virtual bool data_equal(img_node_data* d1, img_node_data* d2) 
        { return ((layer1_data*)d1)->m == ((layer1_data*)d2)->m; }

    virtual void get_masks(vector<img*>&) { }
    virtual void get_part_data(vector<part_data*>&, config_dictionary& cfg) { }
    virtual void get_real_masks(vector<img*>& masks) { get_masks(masks); }
    virtual void get_regions(vector<matrix<bool>*>&, config_dictionary& cfg) { }
    virtual int node_region(node*) { return 0; }

    //void check_geometry(double& benergy, vector<dpoint2>& dvector, double& scdistance, 
    //    const scmap_t& scmap, node* n, const path_map_t& gmap, const ipoint2& center);

	irectangle2 get_box_with_cached_link(node* n);
	irectangle2 get_box(node* n);
    irectangle2 get_box_0(node* n);
    void get_boxes(list<box_data_t>& boxes, part_lib* library, int z, int response, bool use_lib_responses, bool revPCA,
        const response_filter& rfilter, double factor, int border0);
	void get_predicted_boxes(list<box_data_t>& boxes, int z, int response, const response_filter& rfilter, 
        const vector<irectangle2>& pred_boxes, double factor, int border0);
    void get_boxes(list<box_data_t>& boxes, int z, const set<int>& types, double inhibit, int resp_type = G_RESPONSE );
    //void get_boxes(list<pair<irectangle2, node*> >& boxes, int z, int response, int maxn, double factor);
    // ==> inhibit_layer function
    void get_boxes_inhibited(list<box_data_t>& boxes, int z, const set<int>& types, 
        double inhibit, double thresh);

    void check_with_groundtruth(vector<double>& result, const list<irectangle2>& gtr, 
        int z, const set<int>& types, double inhibit);
    vector<node*> filter_with_groundtruth(vector<node*>& nodes, const list<irectangle2>& gtr, double thresh);
	vector<node*> filter_with_groundtruth(vector<node*>& nodes, const list<irectangle2>& gtr, double threshp, double threshn);
    void count_hits(map<int, pair<int, int> >& stat, const list<irectangle2>& gtrs, 
        const list<box_data_t>& boxes, double thresh);
    void count_hits(vector<bool>& hits, int& misses, vector<pair<double, bool> >& confidence,
        const list<irectangle2>& gtrs, const list<box_data_t>& boxes, double thresh);
    void count_hits(vector<bool>& hits, int& misses, vector<pair<double, bool> >& confidence,
        const list<pair<irectangle2, int> >& gtrs, const list<box_data_t>& boxes, double thresh);
    void count_hits_inhibited(vector<bool>& hits, int& misses, vector<pair<double, bool> >& confidence,
        const list<irectangle2>& gtrs, int z, const set<int>& types, double thresh);
    void get_edge_info(set<itriple>& info);

    ddpair top_response_distribution(int layer, int response);
    ddpair response_distribution(int layer, int response);

    void layer_projection_matrix(matrix<lpr_matrix_type>& result, const vector<node*>& nodes, int edgename, int r);

    double convolve(int x, int y, matrix<double>& m, int mt, int z = 0);
    pair<node*, double> convolve_max(int x, int y, matrix<double>& m, int mt, double thresh, int z = 0);
    double convolve_all(int x, int y, matrix<double>& m, int mt, double thresh, int z = 0);
    pair<node*, double> convolve_max_all(int x, int y, matrix<double>& m, int mt, double thresh, int z = 0);

    static int get_closest_node_distance2(vector<node*>& nodes, const point2<int>& p);
    node* find_node_in_region(const point2<int>& p, int z, int deltax, int deltay, 
        int m, unsigned attr, double valthresh = 0.0);
    node* find_node_in_region(node* n, int deltax, int deltay, int m, unsigned attr, double valthresh = 0.0);
    void find_nodes_in_region(vector<node*>& result, int x, int y, int z, int dx, int dy,
        int m, unsigned attr, double valthresh = 0.0);
    node* m_node_at(int x, int y, int z, int m);
    void sorted_nodes_at(vector<node*>& v, int x, int y, int z, int response);

    double schur_product_max(int x, int y, const matrix<double>& m, int type, int z, 
        schur_product_function* f, double thresh, vector<node*>& res);
    dnpair schur_product_max(int x, int y, const matrix<double>& m, int type, int z, 
        schur_product_function* f, double thresh);
    double schur_product_max(int x, int y, const matrix<double>& m, int type, int z, 
        schur_product_function* f, double thresh, double typethresh, bin_function<int, int, double>& typecompare, 
        vector<node*>& res);
    dnpair schur_product_max_all(int xo, int yo, int dx, int dy, int index, matrix<double>& m, int type, int z, double c,
        schur_product_function* f, double thresh, double typethresh, bin_function<int, int, double>& typecompare, sp_result_t& res);
    dnpair schur_product_max_all_best(int xo, int yo, int dx, int dy, int index, matrix<double>& m, int type, int z, double c,
        schur_product_function* f, double thresh, double typethresh, bin_function<int, int, double>& typecompare, sp_result_t& res);
    dnpair schur_product_max_all(int xo, int yo, int dx, int dy, int index, matrix<double>& m, const map<int, double>& appmap, int z, double c,
        schur_product_function* f, double thresh, sp_result_t& res);
    dnpair schur_product_max_all_new(sp_result_t& res, int xo, int yo, int z, double c, schur_product_params& spd);
    dnpair schur_product_max_all_ol(int xo, int yo, int dx, int dy, int index, matrix<double>& m, const map<int, pair<vector<int>, double> >& appmap, int z,
        schur_product_function* f, double thresh, sp_result_t& res);
    dnpair schur_product_max_all_best(int xo, int yo, int dx, int dy, int index, matrix<double>& m, const map<int, double>& appmap, int z, double c,
        schur_product_function* f, double thresh, sp_result_t& res);
    dnpair schur_product_max_all(int xo, int yo, int dx, int dy, int index, matrix<double>& m, int type, int z, double c,  
        schur_product_function* f, double thresh, sp_result_t& res);
    double schur_product_max_all_dir(int xo, int yo, int dx, int dy, int index, double facmin, double facmax, 
        matrix<double>& m, int type, int z, schur_product_function* f, double thresh, sp_result_t& res, double& factor);
    double schur_product_max_all_dir(int xo, int yo, int dx, int dy, int index, double facmin, double facmax, 
        matrix<double>& m, const set<int>& typeset, int z, schur_product_function* f, double thresh, sp_result_t& res, double& factor);
    dnpair schur_product_max_all_best(int xo, int yo, int dx, int dy, int index, matrix<double>& m, int type, int z, double c,
        schur_product_function* f, double thresh, sp_result_t& res);
    dnpair schur_product_max_all(int x, int y, matrix<double>& m, int type, int z, schur_product_function* f, double thresh);
    bool box_empty(int z, const irectangle2& rect);
    node* find_node_in_rect(int z, const irectangle2& rect, const map<int, double>& appmap);
    int find_nodes_in_rect(sp_result_t& result, int z, const irectangle2& rect, const map<int, double>& appmap, int max = INT_MAX);
    void get_nodes(vector<node*>& result, int z, const set<int>& mset);
    void get_nodes(vector<node*>& result, int z0, int z1, const set<int>& mset);
    void get_contractions(vector<double>& contractions);

    void inc_covered(int z) { if (z < (int)info.size()) info[z].covered++; }
    double cover_quotient(int z) 
    { 
        if (z < (int)info.size() && z < (int)shape_nodes.size()) 
            return (double)info[z].covered/(double)shape_nodes[z].size();
        else
            return 0.0;
    }

    virtual void copy_to(streamable* p, cloner& cl);
    virtual streamable* make_instance() const { return new layer1_result(); }
    virtual void write_vgr_label(ostream& os, node* n, int count);
    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);

    void save_visview(const string& dir, const string& name, const string& fname, part_lib* library, int z, int zmax, const map<node*,int>& prevlayermap, int save_mode);
    void save_visview(const string& dir, const string& nameint, const string& fname, part_lib* library, int maxlyr, int save_mode);
    void write_clu(ostream& os, const vector<node*> nodes);
    void write_mma2(ostream& os, const vector<node*>& nodes, int edgename, bool alltypes = false);
    void write_mma2(ostream& os);
    void write_mma3(ostream& os, int m);
    void write_mma4(ostream& os, int slayer, int dlayer);
    void write_matlab(ostream& os, int m);
    void type_statistics(vector<int>& stat, const vector<node*>& nodes, bool next);
    void type_statistics(vector<int>& stat, int layer, bool next);
    void get_covering_statistics(const vector<node*>& nodes, const vector<int>& parts, vector<int>& stat);
    void write_edge_info(ostream& os);
    void write_node_info(ostream& os, const vector<node*>& nodes);
    static void write_edge_info(ostream& os, const set<itriple>& info);
    template<class I> static void select_best_val(list<node*>& result, I begin, I end, int count);

    void set_attribute(unsigned a) { attr |= a; }
    void clear_attribute(unsigned a) { attr &= ~a; }
    bool is_attribute_set(unsigned a) { return (attr & a) != 0; }

	/* void printLyDiff(layer1_result* r2);   REMOVED */
	void debug_print_node_value_less(const double val, int layer);

#ifdef OPENCL
	/**
	  * Opencl functions
	  */

	// This function copies all data from layer1_data into memory object suitable for OpenCL processing.
	void ocl_make_data(int layer, bool overrride_existing_ocl_data = false);
	void make_data_from_ocl(int layer, bool add_edge_names, bool overrride_existing_data = false);

	float ocl_cover_quotient(int z) {
		if (z < (int)info.size() && z < (int)ocl_shape_nodes_coord_non_zero_count.size()) 
            return (float)info[z].covered/(float)ocl_shape_nodes_coord_non_zero_count[z];
        else
            return 0.0;
	}
#endif
protected:
    virtual void copy_to(graph* dest, const map<node*, node*>& cmap);
    void fill_shape_node_vectors();

};

template<class Man> int layer1_result::get_reconstruction_nodes(Man& man, int z, 
    const vector<int>& parts, double thresh, double thresh2, double thresh3, int index)
{
    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);

    man.reset();
    if (shape_nodes[z].empty()) return 0;

    vector<int> sparts(parts.begin(), parts.end());
    vector<node*>& s_nodes = shape_nodes[z];
    bool first_only = thresh < 0.0;
    int type;

    sort(sparts.begin(), sparts.end());
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        if (first_only) {
            type = ((layer1_data*)n->data)->m;
            if (binary_search(sparts.begin(), sparts.end(), type)) man.add_src_node(n); 
        } else {
            do {
                layer1_data* nd = (layer1_data*)n->data;

                if (nd->r(R_RESPONSE) < thresh || nd->r.get_response(G_RESPONSE, 1.0) < thresh2 ||
                        nd->r.get_response(RR_RESPONSE, 1.0) < thresh3) break;
                if (binary_search(sparts.begin(), sparts.end(), nd->m)) man.add_src_node(n); 
                n = nd->next;
            } while (n != nullptr);
        } 
    }
    int result = man.size();

    man.keep_index(index);
    recurse2(man, man.src, atom("toPrevLayer").get_index());
    return result;
}

template<class Man> void layer1_result::get_hypo_node_count(Man& man,
    node* n)
{
    man.reset();
	man.add_src_node(n);
    man.keep_index(-1);
    recurse2(man, man.src, atom("toPrevLayer").get_index());
}

//template<class Man> int layer1_result::get_reconstruction_nodes(Man& man, int z, 
//    const set<node*>& nodes, double thresh, int index /* = -1 */)
//{
//    if (z < 0) z = (int)shape_nodes.size() - 1;
//    else z = min(z, (int)shape_nodes.size() - 1);
//
//    man.reset();
//    if (nodes.empty()) return 0;
//
//    for (set<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
//        node* n = *iter;
//        if (first_only) {
//            man.add_src_node(n); 
//        } else {
//            do {
//                layer1_data* nd = (layer1_data*)n->data;
//                if (nd->val < thresh) break;
//				man.add_src_node(n); 
//                n = nd->next;
//            } while (n != nullptr);
//        } 
//    }
//    int result = man.size();
//    man.keep_index(index);
//    recurse2(man, man.src, atom("toPrevLayer").get_index());
//    return result;
//}

template<class I> void layer1_result::select_best_val(list<node*>& result, I begin, I end, int count)
{
    typedef pair<double, node*> vector_pair_t;
    
    vector<vector_pair_t> v;

    while (begin != end) {
        node* n = *begin;
        layer1_data* nd = (layer1_data*)n->data;

        v.push_back(vector_pair_t(nd->val, n));
        ++begin;
    }
    sort(v.begin(), v.end(), greater<vector_pair_t>());

    int n = 0;

    for (vector<vector_pair_t>::iterator iter = v.begin(); iter != v.end() && n < count; ++iter, ++n) {
        result.push_back(iter->second);
    }
}


// layer1_result_struct
/////////////////////////

class layer1_result_struct : public layer1_result {
public:
    double gabor_lambda;                
    double gabor_gamma;                 
    double gabor_bw;                    
    int gabor_step;    

    layer1_result_struct() : layer1_result() { }

    virtual void get_masks(vector<img*>&); 
    virtual void get_real_masks(vector<img*>&);
    virtual void get_part_data(vector<part_data*>&, config_dictionary&);
    virtual void get_regions(vector<matrix<bool>*>&, config_dictionary&);

    virtual void copy_to(streamable* p, cloner& cl);
    virtual streamable* make_instance() const { return new layer1_result_struct(); }
    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);

protected:
    virtual void copy_to(graph* dest, const map<node*, node*>& cmap);
};

// layer1_result_color
////////////////////////

class layer1_result_color : public layer1_result {
public:
    int gabor_size;                
    int n_rotations;     

    layer1_result_color() : layer1_result() { }

    virtual void get_masks(vector<img*>&); 
    virtual void get_real_masks(vector<img*>&);
    virtual void get_part_data(vector<part_data*>&, config_dictionary&);
    virtual void get_regions(vector<matrix<bool>*>&, config_dictionary&);

    virtual void copy_to(streamable* p, cloner& cl);
    virtual streamable* make_instance() const { return new layer1_result_color(); }
    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);

protected:
    virtual void copy_to(graph* dest, const map<node*, node*>& cmap);
};

// layer1_result_app
//////////////////////

class layer1_result_app : public layer1_result {
public:
    double gabor_lambda;                
    double gabor_gamma;                 
    double gabor_bw;                    
    int gabor_step;    

    layer1_result_app() : layer1_result() { }

    virtual void get_masks(vector<img*>&); 

    virtual combine_t get_combine_function() { return &img::combine_sum; }
    virtual void copy_to(streamable* p, cloner& cl);
    virtual streamable* make_instance() const { return new layer1_result_app(); }
    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);
};

class layer1_result_dog : public layer1_result {
public:
    double sigma_inner;                 
    double sigma_outer;                 
    double mask_size_factor;            

    layer1_result_dog() : layer1_result() { }

    virtual void get_masks(vector<img*>&);

    virtual void copy_to(streamable* p, cloner& cl);
    virtual streamable* make_instance() const { return new layer1_result_dog(); }
    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);
};


class layer1_result_loggabor : public layer1_result {
public:
    layer1_result_loggabor() : layer1_result() { }

    virtual void get_masks(vector<img*>&); 

    virtual streamable* make_instance() const { return new layer1_result_loggabor(); }
};

// ell_learning
///////////////////////////////////////////////////////////////////////////////

struct ell_learning_stat_item {
    //double cosa, sina;
    dpoint2 angle;
    double a;
    double e;

    //ell_learning_stat_item(double pcosa = 0.0, double psina = 0.0, double pa = 0.0, double pe = 0.0) 
    //    : cosa(pcosa), sina(psina), a(pa), e(pe) { }
    ell_learning_stat_item(const dpoint2& pangle = dpoint2::zero, double pa = 0.0, double pe = 0.0) 
        : angle(pangle), a(pa), e(pe) { }

    ell_learning_stat_item(const ellipse& ell) 
    {
        //cosa = cos(ell.angle);
        //sina = sin(ell.angle);
        angle.x = cos(2 * ell.angle);
        angle.y = sin(2 * ell.angle);
        a = 3 * ell.a;
        e = ell.eccentricity();
    }

    static int size() { return 4; }

    template<class T> void push_back_to(vector<T>& v) const
    {
        //v.push_back((T)cosa);
        //v.push_back((T)sina);
        v.push_back((T)angle.x);
        v.push_back((T)angle.y);
        v.push_back((T)a);
        v.push_back((T)e);
    }

};

// global stuff
///////////////////////////////////////////////////////////////////////////////

void read_layer1_result(layer1_result*& res, const string& fname);

layer1_result* read_layer1_result(const string& fname);

void read_layer1_result_visview(layer1_result*& res, const string& fname);

void save_layer1_result(layer1_result* res, const string& fname);

void save_node_set_mma(const string& fname, const set<node*>& nodes);

void inhibit_boxes(list<layer1_result::box_data_t>& boxes, double t);

void keep_best_boxes(list<layer1_result::box_data_t>& boxes, int k);

void check_hits(vector<const layer1_result::box_data_t*>& hitboxes, vector<const layer1_result::box_data_t*>& missboxes,
    vector<bool>& hits, int& misses, const list<pair<irectangle2, int> >& gtrs, 
    const list<layer1_result::box_data_t>& boxes, double thresh);

void inhibit_layer(layer1_result* res, int z, int response, int maxn, double thresh);

void link_path(layer1_result* res, node* n, int edge_name, int link_name);

void get_sc_map(scmap_t& scmap, layer1_result* res, const K_bin& bin, bool normalize);

node* follow_center_link(layer1_result* res, node* n);

int part_geometry_matching(double& benergy, vector<dpoint2>& dvector, double& scdistance, 
    const path_map_t& im, const path_map_t& jm, bool calcbenergy);

int part_geometry_matching(double& benergy, vector<dpoint2>& dvector, double& scdistance, 
    node* p, node* q, bool calcbenergy);

void part_geometry_matching_sym(double& benergy, vector<dpoint2>& dvector, double& scdistance, 
    const path_map_t& im, const path_map_t& jm, bool calcbenergy);

ipoint2 get_path_map(path_map_t& pmap, layer1_result* res, const scmap_t& scmap, node* n, bool link);

void get_node_geo_p(vector<pair<int, ipoint2> >& ptsm, layer1_result* res, node* n);

void get_node_geo(vector<pair<int, ipoint2> >& ptsm, layer1_result* res, node* n);

vector<ipoint2> robust_PCA(layer1_result* res, node* n, node* p);

//void get_boxes(list<layer1_result::box_data_t>& boxes, layer1_result* res, const map<node*, vector<node*> >& clustering, 
//    int cattype, int response, double factor, int border0);

void filter_nodes(list<node*>& nodes, layer1_result* res, int layer, const response_filter& rsf);

vector<node*> filter_nodes(const vector<node*>& nv, const response_filter& rsf);

node* get_closest_node(layer1_result* res, int layer, int type, const dpoint2& p, int minr, int maxr,
    unsigned attr);

irectangle2 bounding_rectangle_of_projection(layer1_result* res, const vector<node*>& nodes);

void sample_tree(set<node*>& result, const set<node*>& ns);

void sample_tree(set<node*>& result, node* n);

void fill_hoc_feature_vector(svm2::vector_t& v, layer1_result* res, int layer, int layer_size, const K_bin& bin, const ipoint2& bcenter,
    const irectangle2& rect);

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

// Ellipse from container of nodes. We assume that nodes are from the
// same layer.
template<class T> void fit_ellipse_to_nodes(ellipse& ell, T begin, T end)
{
    set<ipoint2> pts;

    for (T i = begin; i != end; ++i) {
        node* n = *i;
        img_node_data* nd = (img_node_data*)n->data;

        pts.insert(ipoint2(nd->x, nd->y));
    }
    fit_ellipse(ell, pts);
}
/*
template<class T> void nodes_in_rectangle(set<node*>& nset, const irectangle2& rect, T begin, T end)
{
    nset.clear();

    if (rect.invalid()) return;
    for (T i = begin; i != end; ++i) {
        node* n = *i;
        img_node_data* nd = (img_node_data*)n->data;

        if (rect.inside(nd->x, nd->y)) nset.insert(ipoint2(nd->x, nd->y));
    }
}*/

template<class T> irectangle2 bounding_rectangle_of_nodes(T begin, T end)
{
    irectangle2 result;

    for (T i = begin; i != end; ++i) {
        node* n = *i;
        img_node_data* nd = (img_node_data*)n->data;

        result.eat(nd->x, nd->y);
    }
    return result;
}


template<class T> void node_set_to_point_set(set<ipoint2>& result, T begin, T end)
{
    for (T i = begin; i != end; ++i) {
        node* n = *i;
        img_node_data* nd = (img_node_data*)n->data;

        result.insert(ipoint2(nd->x, nd->y));
    }
}


#endif /* _LAYER_1_RESULT_ */

