// layer_n_creators
///////////////////////////////////////////////////////////////////////////////

#include <queue>
#include "layer_1.h"
#include "layer_n_creators.h"

#ifdef OPENCL
#include "opencl/clContextManager.h"
#include "opencl/cl_utils.h"
#endif

#include <limits>
#include "graphs/graph_utils.h"

// local typedefs & declarations
///////////////////////////////////////////////////////////////////////////////

class layer1_data_h : public layer1_data {
protected:
    set<itriple> voters;
public:
    layer1_data_h(double pval, int pm) : layer1_data(pval, pm), voters() { }
    layer1_data_h(double pval, int pm, int vc, int vx, int vy) : layer1_data(pval, pm), voters() 
        { add_voter(vc, vx, vy); }

    void add_voter(int vc, int vx, int vy) { voters.insert(itriple(vc, vx, vy)); }
    int voters_number() { return (int)voters.size(); }
};

class h_remove_predicate {
protected:
    double valthr;
    int vnthr;
public:
    h_remove_predicate(double pvalthr = 0.0, int pvnthr = 0) : valthr(pvalthr), vnthr(pvnthr) { }

    bool operator()(node* n) const {
        if (!n->is_attr_set(HYPO_NODE_CANDIDATE_ATTR)) return false;

        layer1_data_h* nd = (layer1_data_h*)n->data;
        bool result = nd->val() < valthr || nd->voters_number() < vnthr;
        
        n->data = new layer1_data(*nd);
        n->clear_attr(HYPO_NODE_CANDIDATE_ATTR);
        n->set_attr(HYPO_NODE_ATTR);
        delete nd;
        return result;
    }
};

// layern_creator
///////////////////////////////////////////////////////////////////////////////

// layern_creator support structures
//////////////////////////////////////

// proj_data
//////////////

/// Structure used to store necessary data for new nodes in the process of inference (add_layerX functions)

struct proj_data_base {
    response_map r;     // new response map (values for R{G,RR}_RESPONSE
    int type;           // type
    int x, y;           // position delta (w/r to the positon of the central node)
    layer1_result::sp_result_t* d; // child nodes, i.e. nodes neighbors of the new node
                                   // with toPrevLayer edges. This is a vector of pairs (n, p) 
                                   // where n is child node and p is the name of this node (edge_data_name).
                                   // The name is the coordinate of the subpart (in library).

    proj_data_base(const response_map& vr, int vtype, int vx, int vy, layer1_result::sp_result_t* vd) : 
        r(vr), type(vtype), x(vx), y(vy), d(vd) { }
    ~proj_data_base() {  }

    double val() { return r(G_RESPONSE); }
    bool operator<(const proj_data_base& pd) const { return r(G_RESPONSE) < pd.r(G_RESPONSE); }
    bool operator>(const proj_data_base& pd) const { return r(G_RESPONSE) > pd.r(G_RESPONSE); }
};


struct proj_data_1 : public proj_data_base {
    path_map_t ipmap;

    proj_data_1(const response_map& vr, int vtype, int vx, int vy, layer1_result::sp_result_t* vd, const path_map_t& vipmap) : 
        proj_data_base(vr, vtype, vx, vy, vd), ipmap(vipmap) { }
};

struct proj_data_7 : public proj_data_base {
    vector<node*> msupp;  // "missing" support

    proj_data_7(const response_map& vr, int vtype, int vx, int vy, layer1_result::sp_result_t* vd, const vector<node*>& ms) : 
        proj_data_base(vr, vtype, vx, vy, vd), msupp(ms) { }
};

class exact_comparer : public bin_function<int, int, double> {
public:
    virtual double operator()(const int& t1, const int& t2) const { return (t1 == t2) ? 1.0 : 0.0; }
};

class similarity_comparer : public bin_function<int, int, double> {
    part_lib* library;
    int layer; 
public:
    similarity_comparer(part_lib* l, int ly) : library(l), layer(ly) { }
    virtual double operator()(const int& t1, const int& t2) const 
    { 
        return library->part_similarity(layer, t1, t2); 
    }
};

// layern_creator definitions
///////////////////////////////

//vector<layern_creator*> layern_creator::creators(MAX_LAYER_NUMBER, nullptr);


bool ly1r_type_equal(img_node_data* d1, img_node_data* d2) 
{
    return ((layer1_data*)d1)->m == ((layer1_data*)d2)->m;
}


layern_creator::layern_creator(const char* cfgfile):
	use_opencl(false), opencl_verify_result(false)
{
    config_dictionary cfg(cfgfile);
    init();
    cfg_init(cfg);
}

layern_creator::layern_creator(const config_dictionary& cfg):
	use_opencl(false), opencl_verify_result(false)
{
    init();
    cfg_init(cfg);
}

layern_creator::~layern_creator()
{
	if (library != nullptr && is_library_my == true)  {
		delete library;
	}
}

part_lib* layern_creator::get_library() { 
	return library; 
}

void layern_creator::set_library(part_lib* lib, bool is_my_lib) { 
	// destroy old lib first
	if (is_library_my == true && lib != nullptr && library != lib)
		delete library;

	library = lib; 
	is_library_my = is_my_lib;
#ifdef OPENCL
	if (use_opencl && lib != nullptr && lib->has_opencl_data() == false)
		lib->ocl_make_data();
#endif
}

void layern_creator::set_library(part_lib* lib, const config_dictionary& cfg, bool is_my_lib) { 
	// destroy old lib first
	if (is_library_my == true && lib != nullptr && library != lib)
		delete library;

	library = lib; 
	is_library_my = is_my_lib;

    if (library != nullptr) {
        for (list<node*>::iterator niter = library->nodes.begin(); niter != library->nodes.end(); ++niter) {
            node* p = *niter;
            lib_data* pd = (lib_data*)p->data;

            if (p->is_attr_set(OBJ_PART_ATTR)) {
                string catname = library->get_category(p);

                if (cfg.is_defined("s_response_threshold." + catname)) {
                    pd->td.set_thresh(S_RESPONSE, cfg.get_value_double("s_response_threshold." + catname, 2000.0));
                }
                if (cfg.is_defined("g_response_threshold." + catname)) {
                    pd->td.set_thresh(G_RESPONSE, cfg.get_value_double("g_response_threshold." + catname, 0.0));
                }
                if (cfg.is_defined("realization_ratio_threshold." + catname)) {
                    pd->td.set_thresh(RR_RESPONSE, cfg.get_value_double("realization_ratio_threshold." + catname, 0.0));
                }

            }
        }
    }

#ifdef OPENCL
	if (use_opencl && lib != nullptr && lib->has_opencl_data() == false)
		lib->ocl_make_data();
#endif
}

void layern_creator::mark_texture(layer1_result *res, int k)
{
    vector<node*>& s_nodes = res->shape_nodes[k];
    vector<node*>::iterator iter, niter;
    vector<node*> neighbors;
    
    for (iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        set<int> types;

        neighbors.clear();
        res->get_neighbors_circular(neighbors, n, texture_radius, texture_radius, 0, false);
        for (niter = neighbors.begin(); niter != neighbors.end(); ++niter) {
            layer1_data* nd = (layer1_data*)(*niter)->data;
            types.insert(nd->m);
        }
        if ((int)types.size() >= texture_parts) {
            do {
                n->set_attr(TEXTURE_NODE);
                n = ((layer1_data*)n->data)->next;
            } while (n != nullptr);
        }
    }

}

// Copy uncovered previous layer, i.e. with index k-1 to the layer with index k.
// New nodes get FROM_PREV_LAYER attribute.
// USE ONLY ON THE LAST LAYER
void layern_creator::copy_layer(layer1_result* res, int k)
{
    int k1 = k - 1;

    if (k1 < 0) return;
    
    vector<node*>& s_nodes1 = res->shape_nodes[k1];
    int to_prev_layer = atom("toPrevLayer").get_index();
    int layersize = (int)library->parts[k].size();

    for (vector<node*>::iterator iter = s_nodes1.begin(); iter != s_nodes1.end(); ++iter) {
        node* n = *iter;

        if (n->is_attr_set(HAS_NEXT_LAYER)) continue;
        
        layer1_data* nd = (layer1_data*)n->data;
        int newx = int_round(nd->x/layer_contraction); // !floor!
        int newy = int_round(nd->y/layer_contraction); // !floor!
            
        if (res->node_at(newx, newy, k) == nullptr) {
            layer1_data* newd = new layer1_data(nd->r, nd->m + layersize);
            node* newn = res->add_grid_node(newd, newx, newy, k);

            newn->set_attr(FROM_PREV_LAYER);
            n->set_attr(HAS_NEXT_LAYER);
            res->add_edge(newn, n, to_prev_layer);
        }
    }
}


void layern_creator::inhibit_result(layer1_result* res, int k)
{
    res->update_and_inhibit(k);

    if (g_response_threshold_percent <= 0.0 || g_response_threshold_percent > 1.0) 
        return;

    // Remove excessive hypotheses 
    vector<node*>& nodes = res->shape_nodes[k];
    set<node*> toremove;

    for (vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        node* m = n;
        layer1_data* nd = (layer1_data*)n->data;
        double thresh = nd->val()*g_response_threshold_percent;
        bool tail = false;

        do {
            nd = (layer1_data*)n->data;
            if (nd->val() < thresh) {
                toremove.insert(n);
                if (!tail) {
                    ((layer1_data*)m->data)->next = nullptr;
                    tail = true;
                }
            }
            m = n;
            n = nd->next;
        } while (n != nullptr);
    }
    res->delete_isolated_nodes(toremove);
}

void layern_creator::inhibit_result2(layer1_result* res, int k)
{
    int w = res->x_size(k);
    int h = res->y_size(k);
    layer1_data* d;
    int w1 = w - 1, h1 = h - 1;

    vector<node*>& s_nodes = res->shape_nodes[k];
    vector<node*>& s_nodes_inhib = res->shape_nodes_inhib[k];

    list<node*>::iterator niter;
    vector<node*>::iterator iter;

    // fill s_nodes
    s_nodes.clear();
    for (niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
        d = (layer1_data*)(*niter)->data;
        if (d->z == k && d->x > 0 && d->x < w1 && d->y > 0 && d->y < h1
                && (*niter)->is_attr_set(IMG_NODE_ATTR)) { 
            s_nodes.push_back(*niter);
            // cout << '(' << d->x << ',' << d->y << ')' << ' ' << d->m << ' ';
            // cout << (*niter)->degree(atom("toPrevLayer").get_index()) << endl;
        }
    }
    sort(s_nodes.begin(), s_nodes.end(), layer1_data::less1n); 

    // fill s_nodes_inhib (== s_nodes) and 
    // reduce the number of "toPrevLayer" edges
    int edgename = atom("toPrevLayer").get_index();

    s_nodes_inhib.clear();
    for (iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        multimap<int, edge_pair> newneighbors;

        n->set_attr(NODE_REDUCED_ATTR);
        s_nodes_inhib.push_back(n);
        forall_neighbors(n, niter) {
            if (neighbor_index(niter) != edgename) newneighbors.insert(*niter);
            else {
                node* nn = neighbor_node(niter);

                if (!nn->is_attr_set(NODE_REDUCED_ATTR)) newneighbors.insert(*niter);
            }
        }
        n->neighbors = newneighbors;
    }
}

void layern_creator::init()
{
	is_library_my = false;
    library = nullptr;
    min_cover_quot = 1.0;
    quot_sum = 0.0;
    count = 0;
}

void layern_creator::cfg_init(const config_dictionary& cfg)
{
    layer_contraction = cfg.get_value_double("layer_contraction", 2.0);
    manual_thresholds = cfg.get_value_bool("manual_thresholds", false);
    threshold_factor = cfg.get_value_double("threshold_factor", 1.0);
    if (!cfg.is_defined("candidate_threshold") && cfg.is_defined("proj_prod_threshold"))
        candidate_r_threshold = cfg.get_value_double("proj_prod_threshold", 0.2);
    else
        candidate_r_threshold = cfg.get_value_double("candidate_r_threshold", 0.2);
    candidate_g_threshold = cfg.get_value_double("candidate_g_threshold", 0.2);
    if (!cfg.is_defined("convolution_threshold") && cfg.is_defined("proj_max_threshold"))
        convolution_threshold = cfg.get_value_double("proj_max_threshold", 0.1);
    else
        convolution_threshold = cfg.get_value_double("convolution_threshold", 0.1);
    convolution_link_threshold = cfg.get_value_double("convolution_link_threshold", 0.0);
    if (!cfg.is_defined("r_response_threshold") && cfg.is_defined("proj_add_threshold"))
        r_response_threshold = cfg.get_value_double("proj_add_threshold", 0.3);
    else    
        r_response_threshold = cfg.get_value_double("r_response_threshold", 0.3);
    candidate_r_threshold_percent = cfg.get_value_double("candidate_r_threshold_percent", numeric_limits<double>::quiet_NaN());
    candidate_g_threshold_percent = cfg.get_value_double("candidate_g_threshold_percent", numeric_limits<double>::quiet_NaN());
    r_response_pow = cfg.get_value_double("r_response_pow", 1.0);
    g_response_pow = cfg.get_value_double("g_response_pow", 1.0);
    g_response_threshold = cfg.get_value_double("g_response_threshold", 0.0);
    g_response_threshold_percent = cfg.get_value_double("g_response_threshold_percent", 0.0);
    identity_g_response = cfg.get_value_bool("identity_g_response", false);
    simple_g_response = cfg.get_value_bool("simple_g_response", false);
    ignore_g_distribution = cfg.get_value_bool("ignore_g_distribution", false);
    g_response_operation = cfg.get_value_int("g_response_operation", 1);
    g_response_var_factor = cfg.get_value_double("g_response_var_factor", 0.0);
    proj_max_forb_threshold = cfg.get_value_double("proj_max_forb_threshold", 0.1);
    forb_quot_threshold = cfg.get_value_double("forb_quot_threshold", 0.2);
    copy_prev_layer = cfg.get_value_bool("copy_prev_layer", false);
    texture_parts = cfg.get_value_int("texture_parts", 4);
    ignore_texture = cfg.get_value_bool("ignore_texture", false);
    texture_radius = cfg.get_value_int("texture_radius", 4);
    type_thresh = cfg.get_value_double("type_threshold", 1.0);
    min_factor = cfg.get_value_double("min_factor", 1.0);
    max_factor = cfg.get_value_double("max_factor", 1.0);
    reconstruction_type = cfg.get_value_int("reconstruction_type", 1);
    add_reconstruction_edges = cfg.get_value_bool("add_reconstruction_edges", false);
    add_activation_edges = cfg.get_value_bool("add_activation_edges", true);
    hypo_unrealized_threshold = cfg.get_value_double("hypo_unrealized_threshold", 0.4); 
    hypo_ratio_threshold = cfg.get_value_double("hypo_ratio_threshold", 1E+6);
    realization_ratio_threshold = cfg.get_value_double("realization_ratio_threshold", 0.5);
    shape_check = cfg.get_value_bool("shape_check", false);
    link_missing_support = cfg.get_value_bool("link_missing_support", false);
    continuity_factor = cfg.get_value_double("continuity_factor", -1.0);
    normalize_histogram = cfg.get_value_bool("normalize_histogram", true);
    null_tolerance_limit = cfg.get_value_int("null_tolerance_limit", 3);
    positive_tolerance_threshold = cfg.get_value_int("positive_tolerance_threshold", 2);
    projection_radius = cfg.get_value_int("projection_radius", 3);
    reconstruction_factor = cfg.get_value_double("reconstruction_factor", 0.5);
    rec_null_tolerance_limit = cfg.get_value_int("reconstruction_null_tolerance_limit", 4);
    rec_tolerance_threshold = cfg.get_value_int("reconstruction_tolerance_threshold", 1);
    tolerance_radius = cfg.get_value_int("tolerance_radius", -1);
    min_part_distance = cfg.get_value_int("min_part_distance", -1);

    string ilr = cfg.get_value_string("inhibit_layer_response", "FALSE");

    if (ilr == "G_RESPONSE") inhibit_layer_response = G_RESPONSE;
    else if (ilr == "RR_RESPONSE") inhibit_layer_response = RR_RESPONSE;
    else if (ilr == "R_RESPONSE") inhibit_layer_response = R_RESPONSE;
    else inhibit_layer_response = -1;

    inhibit_layer_max = cfg.get_value_int("inhibit_layer_max", INT_MAX);
    inhibit_layer_thresh = cfg.get_value_double("inhibit_layer_thresh", 0.4);
    inhibit_layer_delete = cfg.get_value_bool("inhibit_layer_delete", true);

    hypo_start = cfg.get_value_double("hypo_start", 0.5);
    hypo_nbhood = cfg.get_value_int("hypo_nbhood", 3);
    hypo_factor = cfg.get_value_double("hypo_factor", 1.2);
    hypo_val_threshold = cfg.get_value_double("hypo_val_threshold", 0.4);
    hypo_voters_threshold = cfg.get_value_int("hypo_voters_threshold", 2);
    
    add_edge_names = cfg.get_value_bool("add_edge_names", false);
    new_positions_threshold = cfg.get_value_int("new_positions_threshold", INT_MAX);
    same_position_threshold = cfg.get_value_int("same_position_threshold", INT_MAX);
    variation_dimension = cfg.get_value_int("variation_dimension", 1);
    variation_factor = cfg.get_value_int("variation_factor", 1);
    scale_merge = cfg.get_value_int("scale_merge", 1);

    s_response_threshold = cfg.get_value_double("s_response_threshold", 2000.0);
    x_response_threshold = cfg.get_value_double("x_response_threshold", 0.0);

    vector<int> allow;

    cfg.get_value(allow, "allowed_parts");
    allowed_parts.insert(allow.begin(), allow.end());
	
	use_opencl = cfg.get_value_bool("use_opencl", false);
	opencl_verify_result = cfg.get_value_bool("opencl_verify_result", false);

	cfg.get_value(use_opencl_devices, "use_opencl_devices");

#ifdef OPENCL

	if (use_opencl == true) {
		if (OpenCL::context_manager->is_initialized == false) {
			cout << "Warning: Found 'use_opencl = true' in config but OpenCL Context IS NOT initialized (maybe missing 'opencl_enable = true' in config)." << endl << "Continuing without OpenCL acceleration." << endl;
			use_opencl = false;
		} else {

			// convert library to opencl data on first use
			if (library != nullptr && library->has_opencl_data() == false) {
				library->ocl_make_data();
			}
		}
	}
#else
	if (use_opencl) {
		cout << "Warning: Found 'use_opencl = true' in config but binary not build with OpenCL support. Will use original implementation !!" << endl;
		use_opencl = false;
	}
#endif
}

void check_layer(layer1_result* res, int layer, const string& print)
{
    int k = layer - 1;

    vector<node*>& s_nodes = res->shape_nodes[k];

    cout << " total of " << s_nodes.size() << " nodes on layer (index) " << k << endl;
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        if (n->is_attr_set(HYPO_NODE_ATTR)) {
            cout << print << " HYPO_NODE_FOUND, data of type " << typeid(*(n->data)).name() << endl;
        } //else {
        //    cout << "NOT A HYPO NODE, data of type " << typeid(*(n->data)).name() << endl;
        //}
    }
}

void layern_creator::add_layer(layer1_result* res, const scmap_t& scmap, int layer, int stoplayer)
{
	// verify that library is build
	if (library == nullptr) {
		cout << "Error: Unable to call layern_creator::add_layer without library !!" << endl;
		throw std::exception();
	}
	
	// use this flag in case opencl is not pressent
	if (use_opencl == true) {
#ifdef OPENCL
		// make opencl data for one layer before since we will need that data for opencl calls
		if (/*res->ocl_shape_nodes == nullptr || */res->ocl_shape_nodes.size() <= (layer - 1) -1 || 
			res->ocl_shape_nodes[(layer - 1) -1].second < res->shape_nodes[(layer - 1) -1].size())
			res->ocl_make_data((layer - 1) -1);
		switch (reconstruction_type) {
			case 2: 
			case 3: 
			case 4: 
			case 5: 				
			case 6:				
			case 7:
			case 8:
			case 9: //add_layer7(res, layer, irectangle2()); break;   // call with an invalid rectangle
			case 10:
				cout << "unsupported opnecl reconstruction_type !!! " << endl;
				throw new std::exception();
			default: ocl_add_layer1(res, layer);
		}
#else
		cout << "Binary not build with OpenCL support !!" << endl;
		throw new std::exception();
#endif			
	} else {
		switch (reconstruction_type) {
			case 2: add_layer2(res, layer); break;
			case 3: add_layer3(res, layer); break;
			case 4: add_layer4(res, layer); break;
			case 5: 
				//check_layer(res, layer - 1, "PRE");
				add_hypothetical_nodes(res, layer - 1);
				//check_layer(res, layer - 1, "POST");
				add_layer5(res, layer); 
				break;
			case 6:
				add_hypothetical_nodes(res, layer - 1);
				break;
			case 7: add_layer5(res, layer); break;
			case 8: add_layer6(res, layer); break;
			case 9: 
                try {
                    add_layer7(res, scmap, layer, irectangle2());   // call with an invalid rectangle
                } catch (libhop_exception& ex) {
                    cout << "lihop_exception '" << ex.get_message() << "'." << endl;
                } catch (exception& ex) {
                    cout << "exception '" << ex.what() << "'." << endl;
                } catch (...) {
                    cout << "Unknown exception." << endl;
                }
                break;
			case 10: add_layer8(res, layer, stoplayer); break;
			default: 
                add_layer1(res, scmap, layer);
                //if (!shape_check) add_layer1(res, layer);
                //else {
                //    K_bin bin(12, 2, 4, 7); // layern_creator field...
                //    scmap_t scmap;

                //    get_sc_map(scmap, res, bin, normalize_histogram);
                //}

		}		
	}

/*    int to_prev = atom("toPrevLayer");
    for (layer1_result::iter_t aiter = res->nodes.begin(); aiter != res->nodes.end(); ++aiter) {
        node* n = *aiter; 
        layer1_data* nd = (layer1_data*)n->data;

        foreach_neighbor(n, to_prev, naiter) {
            edge_data_name* ed = dynamic_cast<edge_data_name*>(neighbor_edge_data(naiter));
            layer1_data* nnd = (layer1_data*)neighbor_node_data(naiter);

            if (ed == nullptr) cout << "From " << nd->x << ", " << nd->y << " ---> " << nnd->x << ", " << nnd->y << endl;
        }
    }
*/
}

void layern_creator::add_layer(layer1_result* res, int layer, int stoplayer)
{
    K_bin bin(12, 2, 4, 7); 
    scmap_t scmap;

    if (shape_check) 
        get_sc_map(scmap, res, bin, normalize_histogram);
    add_layer(res, scmap, layer, stoplayer);
}

// Adds layer to each member (scale) of res; scales are merged according to scale_merge parameter.
// res: a vector of layer1_results on different scales; res[0] is supposed to be the largest scale
// layer: layer to add
// Note: merged members at the end of res are disposed!
void layern_creator::add_layer(vector<layer1_result*>& res, vector<scmap_t>& scmap, int layer, int stoplayer)
{
    //cout << "Reslen: " << res.size() << endl;
    if (scale_merge <= 1) {
        for (int i = 0; i < (int)res.size(); ++i) 
            add_layer(res[i], scmap[i], layer, 0);
    } else {
        for (int i = 0; i <= (int)res.size() - scale_merge; ++i) {
            for (int s = 1; s < scale_merge; ++s)
                res[i]->merge(res[i + s], 100);
            add_layer(res[i], scmap[i], layer, 0); //, stoplayer);
            //cout << "*** " << layer  << endl;
        }
        for (int i = (int)res.size() - scale_merge + 1; i < (int)res.size(); ++i)
            delete res[i];
        res.resize(res.size() - scale_merge + 1);
    }
}

void layern_creator::add_layer(vector<layer1_result*>& res, int layer, int stoplayer)
{
    vector<scmap_t> scmap(res.size());
    K_bin bin(12, 2, 4, 7); 

    if (shape_check) {
        for (int i = 0; i < (int)res.size(); ++i) 
            get_sc_map(scmap[i], res[i], bin, normalize_histogram);
    }
    add_layer(res, scmap, layer, stoplayer);
}


void layern_creator::set_candidate_thresholds(layer1_result* res, int k)
{
    vector<node*>& s_nodes = res->shape_nodes[k];

    if (s_nodes.empty()) 
        return;

	if (!_isnan(candidate_r_threshold_percent)) {
		double max_r;

        if (k == 0)
            max_r = ((layer1_data*)s_nodes.front()->data)->r(R_RESPONSE);
        else {
            max_r = 0.0;
		    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); iter++) {
			    double resp = ((layer1_data*)(*iter)->data)->r(R_RESPONSE);
			    if (resp > max_r) max_r = resp;
		    }		
        }
		candidate_r_threshold = max_r * candidate_r_threshold_percent;
    }
    if (!_isnan(candidate_g_threshold_percent)) {
		double max_r = 0;

        if (k == 0) max_r = 1.0;
        else max_r = ((layer1_data*)s_nodes.front()->data)->r(G_RESPONSE);
		candidate_g_threshold = max_r * candidate_g_threshold_percent;
    }
}

double layern_creator::get_thresh(lib_data* d, int name, double defval)
{
    return manual_thresholds ? defval : d->get_thresh(name, defval, threshold_factor);
}

void filter_links(layer1_result::sp_result_t& conn, const layer1_result::sp_result_t& newconn, double thresh)
{
    if (thresh <= 0.0) 
        conn.insert(conn.end(), newconn.begin(), newconn.end());
    else {
        for (layer1_result::sp_result_t::const_iterator iter = newconn.begin(); iter != newconn.end(); ++iter) 
            if (iter->r >= thresh) conn.push_back(*iter);
    }
}

bool sp_result_compare(const layer1_result::sp_result_data_t& spr1, const layer1_result::sp_result_data_t& spr2)
{
    return spr1.r > spr2.r;
}

// Check if geometry of the part ('imap') matches that of the "model" ('lmap')
// returns 'benergy' (bending TPS energy), 'dvector' (distances between points),
// and 'scdistance' (average of shape context distances).
// Note: Bending energy is calculated only if 'calcbenergy' is true
void layern_creator::check_geometry(double& benergy, vector<double>& dvector, double& scdistance, 
    const path_map_t& lmap, path_map_t& imap, const ipoint2& center,
    bool calcbenergy)
{  
    for (path_map_t::iterator pmiter = imap.begin(); pmiter != imap.end(); ++pmiter)
        pmiter->second.p -= center;

    vector<dpoint2> dispv;
    int unmatched = part_geometry_matching(benergy, dispv, scdistance, lmap, imap, calcbenergy);
    //double f = (double)lmap.size()/(lmap.size() - unmatched); // divide by zero?

    for (vector<dpoint2>::iterator dpiter = dispv.begin(); dpiter != dispv.end(); ++dpiter) 
        dvector.push_back(sqrt(dpiter->norm2()));
    //benergy *= f; 
    scdistance /= (double)lmap.size();
}

// Used internally in add_layer{1,7}
struct inference_data {
    int index;
    node* n;

    inference_data(int pi, node* pn) : index(pi), n(pn) { }
};

void get_library_geo(path_map_t& pm, node* pn, const vector<inference_data>& id)
{
    map<int, int> imap;

    for (vector<inference_data>::const_iterator iter = id.begin(); iter != id.end(); ++iter) {
        imap.insert(pair<int, int>(iter->index, node_type(iter->n)));
    }
    get_library_geo(pm, pn, imap);
}

void get_node_geo(path_map_t& pm, layer1_result* res, const scmap_t& scmap, const vector<inference_data>& idata)
{
    pm.clear();
    for (vector<inference_data>::const_iterator iter = idata.begin(); iter != idata.end(); ++iter) {
        path_map_t tpm;

        get_path_map(tpm, res, scmap, iter->n, false);
        path_map_union(pm, tpm, iter->index);
    }
}

// Extract (labeled) point set 'ptsm' from inference data; 'ptsm' is translated and scaled
pair<double, dpoint2> get_node_geo(vector<pair<int, ipoint2> >& ptsm, layer1_result* res, const vector<inference_data>& idata)
{
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");

    for (vector<inference_data>::const_iterator iter = idata.begin(); iter != idata.end(); ++iter) {
        node* m = iter->n;
        layer1_data* md = (layer1_data*)m->data;
        set<node*> nset;
        set<ipoint2> pset;

        //res->recurse_and_link(m, toprev, to0, nset);
        m->add_to_neighbor_set(to0, nset);

        node_set_to_point_set(pset, nset.begin(), nset.end());

        for (set<ipoint2>::iterator pseti = pset.begin(); pseti != pset.end(); ++pseti) {
            ptsm.push_back(pair<int, ipoint2>(iter->index, *pseti));
        }
    }
    ptsm = inhibit_point_set(ptsm, 5);
    return translate_and_scale(ptsm);
}


void print_geo(ostream& os, const path_map_t& pm)
{
    os << '{';
    for (path_map_t::const_iterator iter = pm.begin(); iter != pm.end(); ++iter) {
        if (iter != pm.begin()) os << ',';
        os << '{' << iter->second.p.x << ',' << iter->second.p.y << '}';
    }
    os << '}';
    os << endl;
}


bool check_continuity(const path_map_t& pm, double factor)
{
    vector<ipoint2> pts = get_path_map_points(pm);

    return continuous_point_set(pts, factor);
}

pair<double, double> check_svm_geometry_p(const vector<pair<int, ipoint2> >& ptsm, const svm_data& pcad, node* p)
{
    if (pcad.svm == nullptr)
        return pair<double, double>(0.0, 0.0);

    lib_data* pd = (lib_data*)p->data;
    vector<pair<int, ipoint2> > ipts = get_library_geo_pieces(p, pd->layer);
    vector<ipoint2> libpts = extract_second<ipoint2>(ipts.begin(), ipts.end());

    translate_and_scale(libpts);

    vector<int> match = piecewise_point_matching_p(ptsm, ipts);
    double avgdelta = 0.0;
    int count = 0;

    // Calculate average delta
    for (int i = 0; i < (int)ipts.size(); ++i) {
        int j = match[i];

        if (j >= 0 && j < (int)ptsm.size()) {
            avgdelta += sqrt((double)(ptsm[j].second.distance2(ipts[i].second)));
            ++count;
        }
    }
    if (count > 0) avgdelta /= count;

    // Create vector, count matched "pieces"
    set<int> matchedpcs, allpcs;
    vector<ipoint2> iptsm(ipts.size());

    for (int i = 0; i < (int)iptsm.size(); ++i) {
        int j = match[i];

        if (j >= 0 && j < (int)ptsm.size()) {
            iptsm[i] = ptsm[j].second;
            matchedpcs.insert(ipts[i].first);
        } else { 
            iptsm[i] = ipts[i].second + ipoint2(avgdelta, avgdelta);
        }
        allpcs.insert(ipts[i].first);
    }
    translate_and_scale(iptsm);


    // Test/train data vector: (x1, y1, squared size1, x2, y2, squared size2, ...)
    int nsamples = (int)libpts.size();
    int nfeatures = 3*nsamples;
    cv::Mat test(1, nfeatures, CV_32FC1, cv::Scalar_<float>(0.0));
    float* ptr = test.ptr<float>(0);

    for (int i = 0; i < nsamples; ++i) {
        int x = iptsm[i].x - libpts[i].x;
        int y = iptsm[i].y - libpts[i].y;

        *(ptr++) = (float)x;
        *(ptr++) = (float)y;
        *(ptr++) = (float)(x*x + y*y);
    }

    float result = pcad.svm->predict(test, true);
    
    //cout << "Pred: " << result << " ";

    return pair<double, double>(result, (double)matchedpcs.size()/allpcs.size());
}

/*double check_svmt(const svm_data& svmd, double gresp, double rrresp, double sresp,
    int nchildren, const vector<inference_data>& children, bool dfvalue)
{
    if (svmd.svm == nullptr) 
        return 1;

    int nfeatures = 3 + nchildren;
    cv::Mat test(1, nfeatures, CV_32FC1, cv::Scalar_<float>(0.0));

    test.at<float>(0, 0) = gresp;
    test.at<float>(0, 1) = rrresp;
    test.at<float>(0, 2) = sresp;
    for (int i = 0; i < (int)children.size(); ++i) {
        layer1_data* nd = (layer1_data*)(children[i].n->data);

        test.at<float>(0, 3 + children[i].index) = nd->r(G_RESPONSE);
    }

    return svmd.svm->predict(test, dfvalue);
}*/

// Return value:
// .first: distance to projection
// .second: % of matched pieces
// 'ptsmiss' are points of the model which do not match with points of ptsm.
// We assume that 'ptsm' are translated and scaled points (e.g. obtained using get_node_geo)
pair<double, double> check_pca_geometry_p(vector<dpoint2>& ptsmiss, const vector<pair<int, ipoint2> >& ptsm, const pca_data& pcad, node* p)
{
    static map<int, vector<pair<int, ipoint2> > > libptsmap;

    lib_data* pd = (lib_data*)p->data;
    vector<pair<int, ipoint2> >& ipts = libptsmap[pd->type];

    if (ipts.empty()) {
        ipts = get_library_geo_pieces(p, pd->layer);
        translate_and_scale(ipts);
    }

    vector<int> match = piecewise_point_matching_p(ptsm, ipts);

    //double avgdelta = 0.0;
    //int count = 0;

    // Calculate average delta
    //for (int i = 0; i < (int)ipts.size(); ++i) {
    //    int j = match[i];

    //    if (j >= 0 && j < (int)ptsm.size()) {
    //        avgdelta += sqrt((double)(ptsm[j].second.distance2(ipts[i].second)));
    //        ++count;
    //    }
    //}
    //if (count > 0) avgdelta /= count;

    // Create vector, count matched "pieces"
    
    set<int> matchedpcs, allpcs;
    vector<dpoint2> dptsm(ipts.size());

    ptsmiss.clear();
    for (int i = 0; i < (int)dptsm.size(); ++i) {
        int j = match[i];

        if (j >= 0 && j < (int)ptsm.size()) {
            dptsm[i] = (dpoint2)ptsm[j].second;
//            sptsm.push_back(((dpoint2)ipts[i].second)/100.0);
            matchedpcs.insert(ipts[i].first);
        } else {
            dptsm[i] = (dpoint2)ipts[i].second; // + dpoint2(random_real(-1, 1), random_real(-1, 1));
            ptsmiss.push_back(dptsm[i]);
        }
        allpcs.insert(ipts[i].first);
    }
    translate_and_scale(dptsm);

    //pair<double, dpoint2> st = translate_and_scale(sptsm);

    //for_each(dptsm.begin(), dptsm.end(), [st](dpoint2& p) { p = p*st.first + st.second; });

    //for (int i = 0; i < (int)ipts.size(); ++i) {
    //    int j = match[i];

    //    if (j < 0 || j >= (int)ptsm.size())
    //        dptsm[i] = ((dpoint2)ipts[i].second)/100.0;
    //}

    //translate_and_scale(dptsm);

    // debug
    //img iml, imi;
    //vector<ipoint2> iptsconv = extract_second<ipoint2>(ipts.begin(), ipts.end());

    //point_img(iml, iptsconv, 10, true);
    //iml.save("c:\\work\\tmp\\libimg.png");
    //point_img(imi, cast_vector<ipoint2, dpoint2>(dptsm), 10, true);
    //imi.save("c:\\work\\tmp\\infimg.png");
    // debug
    //throw;

    cv::Mat data = flatten(dptsm);

    if (data.cols != pcad.mean.cols) {
        cout << "PCA data size mismatch. Wrong library?" << endl;
        throw;
    }

    cv::Mat backp = back_projection(data, pcad);
    //double norm = cv::norm(backp, data, cv::NORM_L2);
    //double prob = erf_value(data, pcad);

    //return pair<double, double>(norm + 1 - prob , (double)matchedpcs.size()/allpcs.size());
    return pair<double, double>(subspace_distance(data, pcad.mean, pcad.eigenvectors, pcad.eigenvalues, 1.0),
        (double)matchedpcs.size()/allpcs.size());
}

pair<double, double> check_pca_geometry_p(const vector<pair<int, ipoint2> >& ptsm, const pca_data& pcad, node* p)
{
    vector<dpoint2> ptsmiss;

    return check_pca_geometry_p(ptsmiss, ptsm, pcad, p);
}

double check_pca_geometry(const path_map_t& pm, const pca_data& pcad, part_data_sim* ed)
{
    //const dpoint2 invalid_dp2(10000.0, 10000.0);

    //if (ed == nullptr || ed->perm.empty()) {
    //    cout << "Invalid part similarity data!" << endl;
    //    throw;
    //}

    //vector<dpoint2> pts(ed->perm.size(), invalid_dp2);

    //for (map<vector<int>, int>::const_iterator iter = ed->perm.begin(); iter != ed->perm.end(); ++iter) {
    //    path_map_t::const_iterator pmiter = pm.find(iter->first);

    //    if (pmiter != pm.end()) pts[iter->second] = (dpoint2)(pmiter->second.p);
    //}

    //// Fill missing points
    //vector<dpoint2> pts0 = pts;

    //for (int i = 0; i < (int)pts.size(); ++i) {
    //    if (pts[i] == invalid_dp2) {
    //        dpoint2 prev = invalid_dp2, next = invalid_dp2;

    //        for (int j = i - 1; j >= 0; --j) 
    //            if (pts0[j] != invalid_dp2) {
    //                prev = pts0[j];
    //                break;
    //            }

    //        for (int j = i + 1; j < (int)pts0.size(); ++j) 
    //            if (pts0[j] != invalid_dp2) {
    //                next = pts0[j];
    //                break;
    //            }
    //        if (prev == invalid_dp2) pts[i] = next;
    //        else if (next == invalid_dp2) pts[i] = prev;
    //        else pts[i] = (next + prev)/2.0;
    //    }
    //}

    //translate_and_scale(pts);

    //cv::Mat data = flatten(pts);
    //cv::Mat backp = back_projection(data, pcad);

    //return cv::norm(backp, data, cv::NORM_L2);    

    vector<dpoint2> pts = cast_vector<dpoint2, ipoint2>(get_path_map_points(pm));
    vector<dpoint2> v = partition(pcad.mean);

    if (pts.size() < v.size()) return numeric_limits<double>::max();

    //continuous_point_set

    pts = get_resized_vector(pts, (int)v.size());
    translate_and_scale(pts);

    permute(pts, point_matching(pts, v));

    cv::Mat data = flatten(pts);
    //cv::Mat backp = back_projection(data, pcad);
    //double prob = erf_value(data, pcad);
    //double norm = cv::norm(backp, data, cv::NORM_L2);

    //return norm + 1 - prob;
    return subspace_distance(data, pcad.mean, pcad.eigenvectors, pcad.eigenvalues, 1.0);
}

double check_geometry2(const vector<pair<int, ipoint2> >& ptsm, node* p)
{
    vector<pair<int, ipoint2> > ipts = get_library_geo_pieces(p, 5);

    translate_and_scale(ipts);

    vector<int> match = piecewise_point_matching_p(ipts, ptsm);
    vector<ipoint2> ptsm2 = extract_second<ipoint2>(ptsm.begin(), ptsm.end());

    vector<dpoint2> iptsd = cast_vector<dpoint2, ipoint2>(extract_second<ipoint2>(ipts.begin(), ipts.end()));
    vector<dpoint2> ptsm2d = cast_vector<dpoint2, ipoint2>(ptsm2);

    translate_and_scale(iptsd);
    translate_and_scale(ptsm2d);

    double result = 0.0;
    for (int i = 0; i < (int)ptsm2d.size(); ++i)
        result += ptsm2d[i].distance2(iptsd[match[i]]);

    return sqrt(result);
}

void layern_creator::add_layer1(layer1_result* res, int layer)
{
    const scmap_t scmap; // Empty map; no geometry checking!

    add_layer1(res, scmap, layer);
}

// Debug checking:
// - add toLayer0 edges;
// - check if all toLayer0 paths end in the 1st layer
void check_consistency(layer1_result* res, int z)
{
    for (auto niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;

        if (nd->z == z) {
            set<node*> nset;

            res->recurse_and_link(n, atom("toPrevLayer"), atom("toLayer0"), nset);
            for (auto siter = nset.begin(); siter != nset.end(); ++siter) {
                if (node_layer(*siter) != 0) {
                    layer1_data* nd = (layer1_data*)(*siter)->data;
                    cerr << "toLayer0 edge does not reach layer 0" << endl;
                    cerr << "m = " << nd->m << endl;
                    cerr << "z = " << nd->z << endl;
                    cerr << "x = " << nd->x << endl;
                    cerr << "y = " << nd->y << endl;
                    cerr << "val = " << nd->vval() << endl;
                    throw;
                }
            }
        }
    }
}

void layern_creator::add_layer1(layer1_result* res, const scmap_t& scmap, int layer)
{
    typedef double product_function_t(int, int, matrix<double>&, int, int, 
        double, double, bin_function<int,int,double>&, vector<node*>&);
    typedef vector<pair<int, proj_data_1*> > parallel_candidates_t;

    if (layer < 1) return;

    int k = layer - 1;
    int k1 = k - 1;

    int x_size = int_round(res->x_size(k1)/layer_contraction); // !floor!
    int y_size = int_round(res->y_size(k1)/layer_contraction); // !floor!
    
	// Sets appropriate thresholds for candidate_r_threshold and candidate_g_threshold.
    // They can be set by user to a value used for all nodes or, in case when 
    // candidate_r_threshold_percent (candidate_g_threshold_percent) are set, as a percent of the 
    // best value at the current position (position of n).
    set_candidate_thresholds(res, k1);

    // Mark "texture nodes" on "prevoius layer"; texture nodes do not index into the inferred layer 
    if (ignore_texture) mark_texture(res, k1);

    if (!res->grid(k)) res->new_grid(x_size, y_size, k);
    while ((int)res->shape_nodes.size() <= k) res->shape_nodes.push_back(vector<node*>());
    while ((int)res->shape_nodes_inhib.size() <= k) res->shape_nodes_inhib.push_back(vector<node*>());
    while ((int)res->info.size() <= k) res->info.push_back(layer_info());
    if (!res->grid(k1)) res->init_grid(k1);

	// Edges to all nodes (i.e. parts) in the next layer (one level up) which have same 
    // part for its center (in library; class part_lib).
    int to_center_back = atom("lyrCenterBack"); 

	// Edges to all subparts (one layer down) other then center (in library; class part_lib).
    int to_part = atom("lyrSrc");

    // Edge to center (in library; class part_lib).
    int to_center = atom("lyrCenter"); 

	// Edges to forbidden subparts -- parts which are not allowed to be present at specific
    // positions -- (one layer down, in library; class part_lib). 
    int to_forbidden = atom("lyrForbidden");

    // Edges to all subparts fro, the previous layer (one layer down) that were used in the 
    // inference process of this part (in class layer1_result)
    int to_prev_layer = atom("toPrevLayer");

	// Edges to root parts ("similarity root node")
    int to_sim_root = atom("lyrSimRoot");

    // Edges to layer 0 (reconstruction edges)
    int to_layer_0 = atom("toLayer0");

    // Edges to variable-shape root parts
    int to_vs_root = atom("lyrVSRoot");
    matrix<double> distr;

	// Parts from the library for layers layer and layer - 1 
    vector<node*>& parts = library->parts[k];
    vector<node*>& parts1 = library->parts[k1];

	// Nodes from image from layer layer - 1
    vector<node*> s_nodes1;
    
    res->get_layer_nodes(s_nodes1, k1);

    if (s_nodes1.empty()) return;

    // Add reconstruction edges in case of shape checking (since this is not thread-safe, it must be
    // done in advance)
    if (shape_check) {
        for (auto niter = s_nodes1.begin(); niter != s_nodes1.end(); ++niter) {
            link_path(res, *niter, to_prev_layer, to_layer_0);
        }
    }
    library->init_basic_part_numbers(k);

    vector<list<proj_data_1> > candidates(s_nodes1.size());

	// For each node found in layer - 1 of image do ...
    // Note that we need to process all nodes from s_nodes as well as the 
    // nodes found at the same positions (those linked using img_graph.next pointers).
    #pragma omp parallel for
    for (int ni = 0; ni < s_nodes1.size(); ++ni) {
		// Part data of the current node n
        node* n = s_nodes1[ni];
        layer1_data* nd = (layer1_data*)n->data;
        node* nc = follow_center_link(res, n);
        layer1_data* ncd = (layer1_data*)nc->data;
        schur_product_params spp(convolution_threshold, 1, 1, ipoint2(ncd->x, ncd->y), nullptr, nullptr);
        part_data_2 centerpd(0, 0, matrix<double>(), 0, 0); // dummy part_data for central part

        gaussian_mask(5, 5, 2.0, centerpd.distr);

		// If the responses of this part are below candidate_r(g)_threshold (see above)
        // then this part is not considered to be a valid center
        if (nd->r(R_RESPONSE) < candidate_r_threshold || nd->r(G_RESPONSE) < candidate_g_threshold ||
                n->is_attr_set(HAS_NEXT_LAYER) || (ignore_texture && n->is_attr_set(TEXTURE_NODE))) {
            continue;
        }

        double ndval_log = log(nd->r(R_RESPONSE)); 
        int x_size_c = int_round(nd->x/layer_contraction); // !floor!
        int y_size_c = int_round(nd->y/layer_contraction); // !floor!
        
		// New part coordinates must be within this bounds to be used..
        if (x_size_c <= 0 || x_size_c + 1 >= x_size || y_size_c <= 0 || y_size_c + 1 >= y_size) {
            continue;
        }

		// Let cpart be a part in the library with the same type as the n (== nd->m).
        node* cpart = parts1[nd->m];
        part_data* cpartd = (part_data*)cpart->data;
        
// used for debuging with opencl comparison
//#define STOP_AT_X 238
//#define STOP_AT_Y 144
//#define STOP_AT_TYPE 5

		// Go through all parts which have cpart for its center (on the next layer) 
        foreach_neighbor(cpart, to_center_back, iter) {
			node* pc = neighbor_node(iter);        // neighbor of cpart
            part_data* pcd = (part_data*)pc->data; // data of pc
            int bpn = library->get_basic_part_number(pc); // basic part number
            double bpr = 0; // count how many basic parts are realized
            //int forb_part_count = 0;
            double max_forb_sum = 0.0;
			int max_forb_count = 0;	
            vector<inference_data> idata;


            #if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
			int new_x = int_round((nd->x + pcd->cmx)/layer_contraction); // !floor!
			int new_y = int_round((nd->y + pcd->cmy)/layer_contraction); // !floor!
			
			if (new_x == STOP_AT_X && new_y == STOP_AT_Y && pcd->type == STOP_AT_TYPE)
					cout << "";
			#endif

            // Check whether this part is allowed to be added; allowed_parts is a parameter.
            // Empty value means no restriction.
            if (!allowed_parts.empty() && allowed_parts.find(pcd->type) == allowed_parts.end())
                continue;

			// Counts how many forbidden subparts (of pc) are found.
            foreach_neighbor(pc, to_forbidden, iter2) {
                node* p = iter2->second.first;
                part_data* pd = (part_data*)p->data;
                part_data_2* ed = (part_data_2*)iter2->second.second;
                dnpair max;
                
				// Gets max value from edge matrix multiplyed by calcualed by layer1_result::rspf .
                // function
                max = res->schur_product_max_all(nd->x + ed->x, nd->y + ed->y, ed->distr, pd->type, 
                    k1, &layer1_result::vspf, convolution_threshold /* unused !! */);
                if (max.first > proj_max_forb_threshold) { ++max_forb_count; max_forb_sum += max.first; }
            }
			// If there are too many forbbiden parts found then skip this part 
            // (try to match another part from lib).
            if (max_forb_count == 0 || max_forb_sum/max_forb_count < forb_quot_threshold) {

				// This part is not forbidden, continue with matching.

                // conn is a vector of matched subparts (nodes) from the previous layer.
                // More precisely, a vector of pairs (n, p) where n is a node and 
                // p indicates to which subpart node n belongs. This is done
                // by a point p (ipoint2) representing the coordinates of the 
                // subpart as specified in the library.
                layer1_result::sp_result_t* conn = new layer1_result::sp_result_t();
                double r_sum = 0.0;             // == R_RESPONSE
                double realization_ratio = 0.0; // == RR_RESPONSE
                double g_val = (g_response_operation == 0) ? 0.0 : 1.0;  // G_RESPONSE
                int realized_count = 0;
                int part_count = 0;
                response_map r;

				#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
				int new_x = int_round((nd->x + pcd->cmx)/layer_contraction); // !floor!
				int new_y = int_round((nd->y + pcd->cmy)/layer_contraction); // !floor!
				if (new_x == STOP_AT_X && new_y == STOP_AT_Y && pcd->type == STOP_AT_TYPE)
					cout << "";
				#endif

                /*
                set_debug_output(false);
                if (k == 3) {
                    set_debug_output(false);
                    foreach_neighbor (pc, to_sim_root, sriter) {
                        lib_data* rpd = (lib_data*)neighbor_node_data(sriter);
                        if (rpd->type == 134) {
                            set_debug_output(true);
                            cout << "DEBUG ON (" << pcd->type << ')' << endl;
                        }
                    }
                }*/

                // For each subpart of the currently selected part do the matching...
                foreach_neighbor(pc, to_part, iter2) {
                    node* p = neighbor_node(iter2); 
                    part_data* pd = (part_data*)p->data;
                    part_data_2* ed = (part_data_2*)neighbor_edge_data(iter2); 
                    g_response_spf gspf(ndval_log, ed->gdistr);
                    layer1_result::sp_result_t newconn;
                    dnpair max;

                    spp.pdata = ed;

					// Do not consider G_RESPONSE, used in the learning phase.
					if (identity_g_response) {
                        spp.f = &layer1_result::idspf;
                        max = res->schur_product_max_all_new(newconn, nd->x, nd->y, k1, 1.0, spp);
					} else if (simple_g_response) {
                        spp.f = &layer1_result::rspf;
                        max = res->schur_product_max_all_new(newconn, nd->x, nd->y, k1, 1.0, spp);
                    } else if (ignore_g_distribution) {
                        spp.f = &layer1_result::sgspf;
                        max = res->schur_product_max_all_new(newconn, nd->x, nd->y, k1, 1.0, spp);
                    } else {
                        if (g_response_var_factor > 0.0) 
                            gspf.dist.reset_variance(gspf.dist.get_variance() * g_response_var_factor);
                        spp.f = &gspf;
                        max = res->schur_product_max_all_new(newconn, nd->x, nd->y, k1, 1.0, spp);
                    }

					// Check if (max) value found by schur_product_max_all is above threshold
                    if (max.first >= convolution_threshold) { 
                        layer1_data* nnd = (layer1_data*)max.second->data;

						#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
						if (new_x == STOP_AT_X && new_y == STOP_AT_Y && pcd->type == STOP_AT_TYPE)
							printf("org: selected node: x: %d; y: %d; m: %d; r_resp: %f; current r_sum: %f; g_val: %f\n", nnd->x, nnd->y, nnd->m, nnd->r(R_RESPONSE), r_sum);
						#endif

						// Update r_sum and g_sum (used in the calculation of R_RESPONSE and G_RESPONSE)
                        // and bpr (used in the calculation of RR_RESPONSE)
                        if (convolution_link_threshold >= 1.0) 
                            conn->push_back(layer1_result::sp_result_data_t(max.second, ed->index, max.first));
                        else 
                            filter_links(*conn, newconn, convolution_link_threshold * max.first);
                        bpr += library->get_basic_part_number(p)*nnd->r(RR_RESPONSE);
                        ++realized_count; 
                        r_sum += nnd->r(R_RESPONSE); 
                        if (g_response_operation == 0) g_val += max.first;
                        else g_val *= max.first;

                        if (shape_check && !scmap.empty()) {
                            idata.push_back(inference_data(ed->index, max.second));
                        }
                    }
                    ++part_count;
                }
				
				#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
				if (new_x == STOP_AT_X && new_y == STOP_AT_Y && pcd->type == STOP_AT_TYPE)
					printf("org: first bpr: %f; bpn: %d; r_sum: %f; g_val: %f\n", bpr, bpn, r_sum, g_val);
				#endif

                // Calculate RR_RESPONSE
                bpr += library->get_basic_part_number(cpart)*nd->r(RR_RESPONSE);
				realization_ratio = bpr/bpn;

                // Calculate R_RESPONSE (from r_sum) and G_RESPONSE (from g_val)
                r_sum += nd->r(R_RESPONSE);
                r_sum /= (part_count + 1);
                r_sum = ::pow(r_sum, r_response_pow);
                g_val = ::pow(g_val * realization_ratio, g_response_pow);

				#if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
				if (new_x == STOP_AT_X && new_y == STOP_AT_Y && pcd->type == STOP_AT_TYPE)
					printf("org:  bpr: %f; bpn: %d; r_sum: %f; g_val: %f\n", bpr, bpn, r_sum, g_val);
				#endif

                //if (layer > 4){
                //    delete conn;
                //    goto endfor;
                //}

				// If R(G)_RESPONSEs are good (w/r to thresholds)
				// then add candidate to the list of possible new nodes.
                // They are inserted to the inference graph below.
                if ((float)r_sum >= (float)get_thresh(pcd, R_THRESH, r_response_threshold) && 
                        (float)realization_ratio >= (float)get_thresh(pcd, RR_THRESH, realization_ratio_threshold) &&
                        (float)g_val >= (float)get_thresh(pcd, G_THRESH, g_response_threshold)) {
                    layer1_result::sp_result_t newconn;
                    dnpair max;
                    part_data_2a* pced = (part_data_2a*)pc->get_edge_data(to_center);
					
					// Use this function to add connections to the center part also.
                    // (we use the schur_product_functions to do this with a "dummy"
                    // map matrix.
                    centerpd.app.clear();
					centerpd.app.insert(part_data_2::app_map_t::value_type(nd->m, part_data_2::app_map_t::mapped_type(vector<int>(), 1.0)));
                    centerpd.geo = pced->geo;
                    spp.f = &layer1_result::rspf;
                    spp.pdata = &centerpd;
					max = res->schur_product_max_all_new(newconn, nd->x, nd->y, k1, 1.0, spp);

                    if (max.first != 0.0) {
                        //throw new_libhop_exception("Center part not added!");
                        path_map_t ipmap;

                        if (convolution_link_threshold >= 1.0) 
                            conn->push_back(layer1_result::sp_result_data_t(max.second, 0, max.first));
                        else
                            filter_links(*conn, newconn, convolution_link_threshold * max.first);

                        if (shape_check && !scmap.empty()) {
                            idata.push_back(inference_data(0, max.second));
                            get_node_geo(ipmap, res, scmap, idata);
                        }

				        // Push this candidate to the list of candidates
                        // In fact, it is a priority_queue, but for no obvious reason...
                        // See proj_data class for more about info stored in the queue.
			            r.set_response(R_RESPONSE, r_sum);
                        r.set_response(RR_RESPONSE, realization_ratio);
                        r.set_response(G_RESPONSE, g_val);
                        //r.set_response(E_RESPONSE, benergy);
                        candidates[ni].push_back(proj_data_1(r, pcd->type, pcd->cmx, pcd->cmy, conn, ipmap));
                        
                        if (get_debug_output()) cout << "ACCEPTED" << endl;
                    } else
                        delete conn; 

                } else 
                    // Delete conn (a vector of connections) since this candidate does not match the
                    // criteria.
                    delete conn; 
            }
        }
    }

    // Make a vector of candidates to be inserted at the same position
    map<ipoint2, parallel_candidates_t> cooccurrent;   // new position -> (candidate index, proj_data*)

    for (int ci = 0; ci < (int)candidates.size(); ++ci) {
        layer1_data* nd = (layer1_data*)(s_nodes1[ci]->data);

        for (auto citer = candidates[ci].begin(); citer != candidates[ci].end(); ++citer) {
            int new_x = int_round((nd->x + citer->x)/layer_contraction); // !floor!
            int new_y = int_round((nd->y + citer->y)/layer_contraction); // !floor!

            cooccurrent[ipoint2(new_x, new_y)].push_back(pair<int, proj_data_1*>(ci, &(*citer)));
        }
    }

    cout << "|";

    vector<parallel_candidates_t> par_candidates;

    for (auto coociter = cooccurrent.begin(); coociter != cooccurrent.end(); ++coociter) 
        par_candidates.push_back(coociter->second);

    // Result of the following block (shape checking): 
    //  - vector of new 'layer1_result's and 
    //  - vector of new edges

    vector<list<new_edges_data_t> > new_edges(par_candidates.size());
    vector<map<int, pair<layer1_data*, new_edges_data_t*> > > new_data(par_candidates.size());

    #pragma omp parallel for
    for (int pci = 0; pci < (int)par_candidates.size(); ++pci) {

        for (auto citer = par_candidates[pci].begin(); citer != par_candidates[pci].end(); ++citer) {
            int ci = citer->first;
            node* n = s_nodes1[ci];
            layer1_data* nd = (layer1_data*)n->data;
            node* nc = follow_center_link(res, n);
            layer1_data* ncd = (layer1_data*)nc->data;
            const proj_data_1& best = *citer->second;

            new_edges[pci].push_back(vector<new_edge_data_t>());

            itriple new_edges_info = get_new_edges(new_edges[pci].back(), k1, best.d->begin(), best.d->end());  // toPrev, toLy0, toHypo

            // if parameter recostruction_edges is true then "toLayer0" edges are added 
            // from the new node. 
            if (!add_reconstruction_edges || new_edges_info.second == 0 || 
                    (double)new_edges_info.third/new_edges_info.second <= hypo_ratio_threshold) {

                // Calculate the new position of the node.
                int new_x = int_round((nd->x + best.x)/layer_contraction); // !floor!
                int new_y = int_round((nd->y + best.y)/layer_contraction); // !floor!
                node* libp = parts[best.type];
                map<node*, part_data_sim*> roots;

                foreach_neighbor (libp, to_vs_root, sriter) {
                    roots.insert(pair<node*, part_data_sim*>(neighbor_node(sriter), nullptr));
                }
                if (roots.empty()) {
                    foreach_neighbor (libp, to_sim_root, sriter) {
                        roots.insert(pair<node*, part_data_sim*>(neighbor_node(sriter), 
                            dynamic_cast<part_data_sim*>(neighbor_edge_data(sriter))));
                    }
                }
                // For each of the roots of the type
                for (map<node*, part_data_sim*>::iterator rtiter = roots.begin(); rtiter != roots.end(); ++rtiter) {
                //foreach_neighbor (libp, to_sim_root, sriter) {
                    //edge_data_t<double>* rped = (edge_data_t<double>*)neighbor_edge_data(sriter);
                    //lib_data* rpd = (lib_data*)neighbor_node_data(sriter);
                    lib_data* rpd = (lib_data*)rtiter->first->data;
                    vs_part_data* vspd = dynamic_cast<vs_part_data*>(rtiter->first->data);
                    double benergy = 0.0, scdistance = 0.0;
                    bool insert = true;

                    if (shape_check) {
                        if (vspd != nullptr) {
                            scdistance = check_pca_geometry(best.ipmap, vspd->pcad, rtiter->second);
                            //cout << "(n=" << norm << ")";
                            insert = scdistance <= s_response_threshold;
                        } else if (!scmap.empty()) {
                            //insert = true; //rpd->type == best.type;
                            vector<double> dv;
                            path_map_t lpmap, ipmap2;

                            get_library_geo(lpmap, libp);
                            lpmap = reduce_path_map(lpmap, 3); // 3!!!??
                            ipmap2 = synchronize_path_map(best.ipmap, lpmap);

                            //ofstream osl("c:\\work\\lpmap.m");
                            //ofstream osi("c:\\work\\ipmap.m");
                            //print_geo(osl, lpmap);
                            //print_geo(osi, ipmap2);
                            //osl.close();
                            //osi.close();

                            check_geometry(benergy, dv, scdistance, lpmap, ipmap2, ipoint2(ncd->x, ncd->y),
                                x_response_threshold > 0.0); //, *conn);
                            insert = scdistance <= s_response_threshold && benergy <= x_response_threshold;
                        }
                    }
                    if (insert) {
                        layer1_data* ly1d = new layer1_data(best.r, rpd->type); 

                        ly1d->x = new_x;
                        ly1d->y = new_y;
                        ly1d->z = k;
                        ly1d->r.set_response(S_RESPONSE, scdistance);
                        ly1d->r.set_response(X_RESPONSE, benergy);

                        // check if the same type (rpd->type) exists at the current position and keep 
                        // the one with the higher value (G_RESPONSE)

                        bool exists = new_data[pci].find(ly1d->m) != new_data[pci].end();
                        pair<layer1_data*, new_edges_data_t*>& depair = new_data[pci][ly1d->m];

                        if (!exists) {
                            depair.first = ly1d;
                            depair.second = &new_edges[pci].back();
                        } else if (depair.first->vval() >= ly1d->vval()) {
                            delete ly1d;
                        } else {
                            delete depair.first;
                            depair.first = ly1d;
                            depair.second = &new_edges[pci].back();
                        }

                        if (!n->is_attr_set(HAS_NEXT_LAYER)) {  // ?
                            res->inc_covered(k1);
                            n->set_attr(HAS_NEXT_LAYER);
                        }

                        
                        // Find node with the type of the new node (best.type) at this position.
                        /*pair<node*, node**> fpair = res->find_at(ly1d, new_x, new_y, k);
                        node* newn = fpair.first;

				        #if defined STOP_AT_X && defined STOP_AT_Y && defined STOP_AT_TYPE
				        if (new_x == STOP_AT_X && new_y == STOP_AT_Y && best.type == STOP_AT_TYPE)
					        cout << "";
				        #endif

                        if (newn == nullptr) {
                            // There is no node of this type, insert
                            newn = res->add_grid_node(ly1d, new_x, new_y, k);
        					
                            if (!n->is_attr_set(HAS_NEXT_LAYER)) {
                                res->inc_covered(k1);
                                n->set_attr(HAS_NEXT_LAYER);
                            }

                        } else {
                            if (!res->data_greater(ly1d, (layer1_data*)newn->data)) 
                                // R_RESPONSE is smaller, do not insert!
                                newn = nullptr; 
                            else {
                                // There is a node of the same type, remove the old one
                                res->remove_grid_node(fpair);
                                ((layer1_data*)newn->data)->r = ly1d->r;
                                newn->clear_neighbors();
                                res->add_grid_node(newn);
                            } 
                            delete ly1d;

                        }
                        if (newn != nullptr) {
                            // We add edges from the new node to all "constructing" nodes.
                            add_edges(res, newn, new_edges);
                        }*/
                    }
                }
            }
            delete best.d;
        }

    }

    // We insert the candidates into the graph; at each position we insert nodes 
    // sorted by best value (G_RESPONSE)

    for (int ndi = 0; ndi < (int)new_data.size(); ++ndi) {
        vector<pair<layer1_data*, new_edges_data_t*> > v;

        v.reserve(new_data[ndi].size());
        for (auto nditer = new_data[ndi].begin(); nditer != new_data[ndi].end(); ++nditer) 
            v.push_back(nditer->second);
        sort(v.begin(), v.end(), [](const pair<layer1_data*, new_edges_data_t*>& a, const pair<layer1_data*, new_edges_data_t*>& b) ->
            bool { 
                double aval = a.first->vval(), bval = b.first->vval();
                return abs(aval - bval) < 1E-6 ? a.first->m > b.first->m : aval < bval;
            });
        // take last 'same_position_threshold' elements

        if ((int)v.size() > same_position_threshold) {
            for (int i = 0; i < (int)v.size() - same_position_threshold; ++i)
                delete v[i].first;
            v.erase(v.begin(), v.begin() + (v.size() - same_position_threshold));
        }

        for (auto viter = v.begin(); viter != v.end(); ++viter) {
            node* newn = res->add_grid_node(viter->first, viter->first->x, viter->first->y, k);

            /*if (viter != v.begin() && viter->first->vval() < (viter - 1)->first->vval()) {
                cerr << viter->first->vval() << " < " << (viter - 1)->first->vval() << endl;
                throw;
            }*/

            int error = 0;
            for (auto diter = viter->second->begin(); diter != viter->second->end(); ++diter) {
                if (diter->n == nullptr) { error = 1; }
                if (!newn->is_attr_set(IMG_NODE_ATTR)) { 
                    error = 2; cout << viter->first->vval() << '/' << (viter - 1)->first->vval() << ' ';
                    cout << viter->first->m << '/' << (viter - 1)->first->m << ' ';
                }
                if (node_layer(diter->n) != k1) { error = 3; }
            }
            if (error > 0) {
                cout << "Error " << error;
                throw;
            }
            if (viter->second->empty()) cerr << "Empty set of edges." << endl;
            add_edges(res, newn, *viter->second);
        }
    }
	
    // Just for statistics...
    double quot = res->cover_quotient(k1);

    count++;
    quot_sum += quot;
    if (quot < min_cover_quot) min_cover_quot = quot;

    if (copy_prev_layer) copy_layer(res, k);

    inhibit_result(res, k);
    res->delete_edges(atom("toLayer0"));
    //check_consistency(res, k); 
    //res->delete_edges(atom("toLayer0"));

}

void layern_creator::add_layer2(layer1_result* res, int layer)
{
    throw new_libhop_exception("layern_creator::add_layer2 is not supported.");
}

// Use this function to add category layer only!
void layern_creator::add_layer3(layer1_result* res, int layer)
{
    int k = layer - 1;
    int k1 = k - 1;

    if (k <= 0) return;

    int x_size = res->x_size(k1);
    int y_size = res->y_size(k1);
    
    res->new_grid(x_size, y_size, k);
    while ((int)res->shape_nodes.size() <= k) res->shape_nodes.push_back(vector<node*>());
    while ((int)res->shape_nodes_inhib.size() <= k) res->shape_nodes_inhib.push_back(vector<node*>());
    while ((int)res->info.size() <= k) res->info.push_back(layer_info());

    vector<node*>& s_nodes1 = res->shape_nodes[k1];
    vector<node*>& s_nodes = res->shape_nodes[k];
    vector<node*>& s_nodes_inhib = res->shape_nodes_inhib[k];
    vector<node*>& parts1 = library->parts[k1];

    int to_prev_layer = atom("toPrevLayer");
    int to_center = atom("lyrCenterBack");
        
    for (vector<node*>::iterator niter = s_nodes1.begin(); niter != s_nodes1.end(); ++niter) {
        node* n = *niter;

        do {
            layer1_data* nd = (layer1_data*)n->data;
            node* cpart = parts1[nd->m];
            part_data* cpartd = (part_data*)cpart->data;

            foreach_neighbor(cpart, to_center, iter) {
			    node* pc = neighbor_node(iter);        // neighbor of cpart
                part_data* pcd = (part_data*)pc->data; // data of pc

                if (!pc->is_attr_set(CATEGORY_NODE_ATTR)) 
                    continue;

                layer1_data* ly1d = new layer1_data(nd->r, pcd->type);
                pair<node*, node**> fpair = res->find_at(ly1d, nd->x, nd->y, k);
                node* newn = fpair.first;

                //if (newn == nullptr) {
                    newn = res->add_grid_node(ly1d, nd->x, nd->y, k);
				    if (!n->is_attr_set(HAS_NEXT_LAYER)) {
                        res->inc_covered(k1);
                        n->set_attr(HAS_NEXT_LAYER);
                    }
                //} else {
                    //if (!res->data_greater(ly1d, (layer1_data*)newn->data)) 
                    //    newn = nullptr; 
                    //else {
                        //res->remove_grid_node(fpair);
                        //((layer1_data*)newn->data)->r = ly1d->r;
                        //newn->clear_neighbors();
                        //res->add_grid_node(newn);
                    //} 
                    //delete ly1d;
                //}
                if (newn != nullptr) {
                    res->add_edge(newn, n, to_prev_layer);
                }
            }
            n = nd->next;
        } while (n != nullptr);
    }
    inhibit_result(res, k);
}

void layern_creator::add_layer4(layer1_result* res, int layer)
{
    throw new_libhop_exception("layern_creator::add_layer4 is not supported");
}

// Version which does not check a similarity between parts
// - res = graph
// - layer = name of the layer (i.e. index + 1)
// Remarks: 
void layern_creator::add_layer5(layer1_result* res, int layer)
{
    int k = layer - 1;
    int k1 = k - 1;

    if (k <= 0) return;
    //if (res->grid(k)) return;

    int x_size = int_round(res->x_size(k1)/layer_contraction); // !floor!
    int y_size = int_round(res->y_size(k1)/layer_contraction); // !floor!
    int max_x = res->x_size(k1);
    int max_y = res->y_size(k1);

    
    res->new_grid(x_size, y_size, k);
    while ((int)res->shape_nodes.size() <= k) res->shape_nodes.push_back(vector<node*>());
    while ((int)res->shape_nodes_inhib.size() <= k) res->shape_nodes_inhib.push_back(vector<node*>());
    while ((int)res->info.size() <= k) res->info.push_back(layer_info());

    vector<node*>& s_nodes1 = res->shape_nodes[k1];

    //cout << "------ k = " << k << "    " << s_nodes1.size() << "-------" << endl;

    if (s_nodes1.empty()) return;

    vector<node*>& parts1 = library->parts[k1];
    int to_center = atom("lyrCenterBack").get_index();
    int to_part = atom("lyrSrc").get_index();
    int to_prev_layer = atom("toPrevLayer").get_index();
    int to_forbidden = atom("lyrForbidden").get_index();
    matrix<double> distr;
    vector<node*>::iterator niter = s_nodes1.begin();
    node* n = *niter;
    int COUNTER = 0;
    int x_size_c, y_size_c;

    gaussian_mask(5, 5, 2.0, distr);
    while (n != nullptr) {
        layer1_data* nd = (layer1_data*)n->data;
        node* cpart = parts1[nd->m];
        part_data* cpartd = (part_data*)cpart->data;
        priority_queue<proj_data_base> candidates;

        if (n->is_attr_set(HAS_NEXT_LAYER) || 
            (ignore_texture && n->is_attr_set(TEXTURE_NODE))) goto increment;

        x_size_c = int_round(nd->x/layer_contraction); // !floor!
        y_size_c = int_round(nd->y/layer_contraction); // !floor!

        if (x_size_c <= 0 || x_size_c + 1 >= x_size || 
            y_size_c <= 0 || y_size_c + 1 >= y_size) goto increment;

        foreach_neighbor(cpart, to_center, iter) {
            node* pc = iter->second.first;
            part_data* pcd = (part_data*)pc->data;
            //int forb_part_count = 0;
            double max_forb_sum = 0.0;
            int max_forb_count = 0;
            dnpair max;

            foreach_neighbor(pc, to_forbidden, iter2) {
                node* p = iter2->second.first;
                part_data* pd = (part_data*)p->data;
                part_data_2* ed = (part_data_2*)iter2->second.second;
                
                max = res->schur_product_max(nd->x + ed->x, nd->y + ed->y, ed->distr, pd->type, 
                    k1, R_RESPONSE, candidate_r_threshold);
                if (max.first > proj_max_forb_threshold) { ++max_forb_count; max_forb_sum += max.first; }
                //++forb_part_count;
            }

            if (max_forb_count == 0 || max_forb_sum/max_forb_count < forb_quot_threshold) {
                layer1_result::sp_result_t* conn = new layer1_result::sp_result_t();
                double max_sum = 0.0;
                double realization_ratio = 0.0;
                int max_thresh_count = 0;
                int part_count = 0;
                response_map r;

                foreach_neighbor(pc, to_part, iter2) {
                    node* p = iter2->second.first;
                    part_data* pd = (part_data*)p->data;
                    part_data_2* ed = (part_data_2*)iter2->second.second;
                    dnpair max;

                    max = res->schur_product_max_all(nd->x, nd->y, ed->x, ed->y, ed->index, ed->distr, pd->type, 
                        k1, 1.0, R_RESPONSE, candidate_r_threshold, *conn);
                    if (max.first >= convolution_threshold) {
                        ++max_thresh_count; max_sum += max.first; 
                    }
                    ++part_count;
                }
                realization_ratio = (part_count == 0) ? 1.0 : (double)max_thresh_count/((double)part_count);
                max_sum += nd->val();
                max_sum /= (part_count + 1);
                max_sum = ::pow(max_sum, r_response_pow);

                if (max_sum > get_thresh(pcd, R_THRESH, r_response_threshold) && 
                        realization_ratio >= get_thresh(pcd, RR_THRESH, realization_ratio_threshold)) {
                    res->schur_product_max_all(nd->x, nd->y, 0, 0, 0, distr, nd->m, k1, 1.0, R_RESPONSE, 
                        candidate_r_threshold, *conn);
                    r.set_response(R_RESPONSE, max_sum);
                    r.set_response(RR_RESPONSE, realization_ratio);
                    candidates.push(proj_data_base(r, pcd->type, pcd->cmx, pcd->cmy, conn));
                } else 
                    delete conn;
            }
        }

        // connect with the best
        while (!candidates.empty()) {
            const proj_data_base& best = candidates.top();
            vector<new_edge_data_t> new_edges;
            itriple new_edges_info = get_new_edges(new_edges, k1, best.d->begin(), best.d->end());  // toPrev, toLy0, toHypo

            if (k == 5 && new_edges_info.second != 0) cout << (double)new_edges_info.third/new_edges_info.second << ' ';

            if (!add_reconstruction_edges || new_edges_info.second == 0 ||
                    (double)new_edges_info.third/new_edges_info.second <= hypo_ratio_threshold) {
                int new_x = int_round((nd->x + best.x)/layer_contraction); // !floor!
                int new_y = int_round((nd->y + best.y)/layer_contraction); // !floor!

                if (new_x >= 0 && new_x < x_size && new_y >= 0 && new_y < y_size) {
                    layer1_data* ly1d = new layer1_data(best.r, best.type);
                    pair<node*, node**> fpair = res->find_at(ly1d, new_x, new_y, k);
                    node* newn = fpair.first;

                    if (newn == nullptr) {
                        newn = res->add_grid_node(ly1d, new_x, new_y, k);
                    } else {
                        if (!res->data_greater(ly1d, (layer1_data*)newn->data)) newn = nullptr;
                        else {
                            res->remove_grid_node(fpair);
                            ((layer1_data*)newn->data)->r = ly1d->r;
                            newn->clear_neighbors();
                            res->add_grid_node(newn);
                        } 
                        delete ly1d;
                    }
                    if (newn != nullptr) 
                        add_edges(res, newn, new_edges);

                    if (!n->is_attr_set(HAS_NEXT_LAYER)) {
                        res->inc_covered(k1);
                        n->set_attr(HAS_NEXT_LAYER);
                    }
                }
            }
            delete best.d;
            candidates.pop();
        }
    increment:
        n = nd->next;
        if (n == nullptr && ++niter != s_nodes1.end()) n = *niter;
    }

    double quot = res->cover_quotient(k1);

    count++;
    quot_sum += quot;
    if (quot < min_cover_quot) min_cover_quot = quot;

    if (copy_prev_layer) copy_layer(res, k);

    inhibit_result(res, k);
}

void layern_creator::add_layer6(layer1_result* res, int layer)
{
    throw new_libhop_exception("layern_creator::add_layer6 is not supported.");
}

void get_shape_xml(vector<float>& result, layer1_result* res, spart_data* pd, layer1_result::sp_result_t& conn)
{
    typedef map<vector<int>, ipoint2> vmap_t;

    int to_src = atom("lyrSrc");
    int to_prev = atom("toPrevLayer");
    int to_0 = atom("toLayer0");
    set<ipoint2> pts;

    for (layer1_result::sp_result_t::iterator iter = conn.begin(); iter != conn.end(); ++iter) {
        set<node*> nset;
        set<ipoint2> pset;

        link_path(res, iter->n, to_prev, to_0);
        iter->n->get_neighbor_set(to_0, nset);
        node_set_to_point_set(pset, nset.begin(), nset.end());
        pts.insert(pset.begin(), pset.end());
    }
    ipoint2 mean;
    double scale;

    ipoint2_set_normalization_parameters(mean, scale, pts.begin(), pts.end());
    
    const shape_model& sm = pd->get_model();
    int sum = 0;
    int sumcount = 0;
    vmap_t vmap;

    for (layer1_result::sp_result_t::iterator iter = conn.begin(); iter != conn.end(); ++iter) {
        int count = iter->n->count_neighbors(to_0);
        
        foreach_neighbor(iter->n, to_0, niter) {
            layer1_data* nd = (layer1_data*)neighbor_node_data(niter);
            edge_path_data_t* ed = (edge_path_data_t*)neighbor_edge_data(niter);
            vector<int> v = ed->data.edges();

            v.push_back(iter->index);
            vmap.insert(vmap_t::value_type(v, 
                ipoint2((int)(100*(nd->x - mean.x)/scale), (int)(100*(nd->y - mean.y)/scale))));
            //int d = sm.distance2(v, ipoint2((int)(100*(nd->x - mean.x)/scale), (int)(100*(nd->y - mean.y)/scale)));

            //if (random_int(0, count) <= 0.25*count) {
            //    if (d != INT_MAX) sum += d; else { sum += 5000; }  //// !!!!!????!!!
            //    ++sumcount; 
            //}

        }
    }

    const shape_model::map_t& smap = pd->get_model().get_map();
    int ds = (int)smap.size();
    int index = 0;

    result.clear();
    result.resize(3*ds, 0);

    for (shape_model::map_t::const_iterator diter = smap.begin(); diter != smap.end(); ++diter) {
        vmap_t::iterator siter = vmap.find(diter->first);

        if (siter != vmap.end()) {
            result[index] = 1;
            result[index + ds] = siter->second.x;
            result[index + 2*ds] = siter->second.y;
        }
        ++index;
    }

}

double check_shape(layer1_result* res, spart_data* pd, layer1_result::sp_result_t& conn, bool strict)
{
    int to_src = atom("lyrSrc");
    int to_prev = atom("toPrevLayer");
    int to_0 = atom("toLayer0");
    set<ipoint2> pts;

    for (layer1_result::sp_result_t::iterator iter = conn.begin(); iter != conn.end(); ++iter) {
        set<node*> nset;
        set<ipoint2> pset;

        link_path(res, iter->n, to_prev, to_0);
        iter->n->get_neighbor_set(to_0, nset);

                    //n->get_neighbor_set(to0, supp);
                    //node_set_to_point_set(psupp, supp.begin(), supp.end());

                    //double q = (double)psupp.size()/p0.size();

                    //if (q >= cover_threshold) {
                    //    ipoint2 mean;
                    //    double scale;

                    //    ipoint2_set_normalization_parameters(mean, scale, psupp.begin(), psupp.end());
                    //    foreach_neighbor(n, to0, nniter) {
                    //        edge_path_data_t* ned = (edge_path_data_t*)neighbor_edge_data(nniter);
                    //        layer1_data* nnd = (layer1_data*)neighbor_node_data(nniter);
                    //        ipoint2 p((int)(100*(nnd->x - mean.x)/scale), (int)(100*(nnd->y - mean.y)/scale));

                    //        makermap[nd->m].update(ned->data, p);
                    //    }
                    //}




        node_set_to_point_set(pset, nset.begin(), nset.end());
        pts.insert(pset.begin(), pset.end());
    }
    ipoint2 mean;
    double scale;

    ipoint2_set_normalization_parameters(mean, scale, pts.begin(), pts.end());

    //select random 10% points and calculate distance
    const shape_model& sm = pd->get_model();
    int sum = 0;
    int sumcount = 0;

    for (layer1_result::sp_result_t::iterator iter = conn.begin(); iter != conn.end(); ++iter) {
        int count = iter->n->count_neighbors(to_0);
        
        foreach_neighbor(iter->n, to_0, niter) {
            layer1_data* nd = (layer1_data*)neighbor_node_data(niter);
            edge_path_data_t* ed = (edge_path_data_t*)neighbor_edge_data(niter);
            vector<int> v = ed->data.edges();

            v.push_back(iter->index);

            int d = sm.distance2(v, ipoint2((int)(100*(nd->x - mean.x)/scale), (int)(100*(nd->y - mean.y)/scale)));

            if (d > sum) sum = d;
            //if (random_int(0, count) <= 0.25*count) {
            //    if (d != INT_MAX) sum += d; else { sum += 5000; }  //// !!!!!????!!!
            //    ++sumcount; 
            //}
        }
    }

    //return (double)sum/sumcount;
    return sum;
    
    ///////
    /*sm.get_matrix().save_mathematica("matrix.m");
    ofstream os("pts.m");
    os << '{';
    for (set<ipoint2>::iterator i = pts.begin(); i != pts.end(); ++i) {
        if (i != pts.begin()) os << ',';
        os << '{' << i->x << ',' << i->y << '}';
    }
    os << '}';
    os.close(); 

    int dummy;
    cout << (double)hits/sm.region_count();
    cout << "; Waiting to enter some character: ";
    cin >> dummy;*/
    //////

    //return (double)hits/sm.region_count(); 
}

void ell_svm_vector(vector<float>& v, layer1_result* res, node* p, layer1_result::sp_result_t& conn)
{
    typedef map<int, set<node*> > map_t;

    int to_src = atom("lyrSrc");
    int to_prev = atom("toPrevLayer");
    int to_0 = atom("toLayer0");
    int feature_count = 4; //4*(p->count_neighbors(to_src) + 1);
    map_t spmap;
    set<node*> objnodes;

    for (layer1_result::sp_result_t::iterator iter = conn.begin(); iter != conn.end(); ++iter) {
        map_t::iterator miter = spmap.find(iter->index);
        set<node*> nset;

        res->recurse_and_link(iter->n, to_prev, to_0, nset);
        //iter->n->get_neighbor_set(to_0, nset);
        objnodes.insert(nset.begin(), nset.end());
        if (miter == spmap.end()) spmap.insert(map_t::value_type(iter->index, nset));
        else miter->second.insert(nset.begin(), nset.end());
    }

    v.clear();
    v.reserve(feature_count + 1);
    ellipse ell;

    fit_ellipse_to_nodes(ell, objnodes.begin(), objnodes.end());
    ell_learning_stat_item(ell).push_back_to(v);
    foreach_neighbor(p, to_src, iter) {
        part_data* spd = (part_data*)neighbor_node_data(iter);
        part_data_2* sped = (part_data_2*)neighbor_edge_data(iter);

        map_t::iterator miter = spmap.find(sped->index);
        if (miter == spmap.end()) ell_learning_stat_item().push_back_to(v);
        else {
            fit_ellipse_to_nodes(ell, miter->second.begin(), miter->second.end());
            ell_learning_stat_item(ell).push_back_to(v);
        }

    }

}

float dist(const vector<float>& v1, const vector<float>& v2)
{
    float result = 0.0;

    for (int i = 0; i < (int)v1.size() && i < (int)v2.size(); ++i) {
        float x = v1[i] - v2[i];

        result += x*x;
    }
    return sqrt(result);
}

float dist(const list<vector<float> >& l, const vector<float>& v)
{
    typedef list<vector<float> > list_t;

    float mind = numeric_limits<float>::infinity();

    for (list_t::const_iterator iter = l.begin(); iter != l.end(); ++iter) {
        float d = dist(*iter, v);
        if (d < mind) mind = d;
    }
    return mind;
}

float predict(const list<vector<float> >& pos, const list<vector<float> >& neg, const vector<float>& v)
{
    float pdist = dist(pos, v);
    float ndist = dist(neg, v);
    if (pdist <= ndist) return 0.0; else return 1.0;
}

void get_support_nodes(vector<node*>& supp, layer1_result* res, const vector<ipoint2>& misspos, int delta)
{
    if (!res->grid(0)) 
        res->init_grid(0);

    for (auto piter = misspos.begin(); piter != misspos.end(); ++piter) {
        //node* n = res->get_closest_node(piter->x, piter->y, 0, delta, delta);
        vector<node*> result = res->get_grid_nodes_circular(piter->x, piter->y, 0, delta, delta);

        supp.insert(supp.end(), result.begin(), result.end());
    }
}

void layern_creator::get_new_support_edges(new_edges_data_t& edges, const vector<node*>& nodes) 
{
    const int ename = atom("toPrevLayer");

    for (auto niter = nodes.begin(); niter != nodes.end(); ++niter) {
        edges.push_back(new_edge_data_t(*niter, ename, invalid_edge_name, 1.0));
    }
}

bool layern_creator::add_layer7(layer1_result* res, int layer, const irectangle2& region)
{
    scmap_t scmap;

    return add_layer7(res, scmap, layer, region);
}

int debug_number = 0;

// Version which handles "object layer" 
// If 'region' is not invalid, it considers only nodes in the specified rectangle.
// Return value: true if some node was added, false otherwise.
bool layern_creator::add_layer7(layer1_result* res, const scmap_t& scmap, int layer, const irectangle2& region)
{
    typedef vector<pair<int, proj_data_7*> > parallel_candidates_t;

    int to_center_back = atom("lyrCenterBack");
    int to_part = atom("lyrSrc");
    int to_prev_layer = atom("toPrevLayer");
    int to_layer_0 = atom("toLayer0");
    int to_prev_layerI = atom("toPrevLayerI");

    int k = layer - 1;
    int k1 = k - 1;

    if (k <= 0) return false;

    vector<node*>& parts1 = library->parts[k1];
    vector<node*>& parts = library->parts[k];

    if (res->max_layer_index() < k1 || res->shape_nodes[k1].empty() || 
        parts1.empty() || parts.empty()) return false;

    while ((int)res->shape_nodes.size() <= k) res->shape_nodes.push_back(vector<node*>());
    while ((int)res->shape_nodes_inhib.size() <= k) res->shape_nodes_inhib.push_back(vector<node*>());
    while ((int)res->info.size() <= k) res->info.push_back(layer_info());

    set_candidate_thresholds(res, k1);

    // Set the contraction; either use the value from cfg file or, in case of
    // "rlayer" sets the predetermined correct values.
    double contraction = layer_contraction;
    
    int x_size = int_round(res->x_size(k1)/contraction); // !floor!
    int y_size = int_round(res->y_size(k1)/contraction); // !floor!

    vector<ipoint2> directions;
    int dirnum = variation_dimension/2;

    for (int i = -dirnum; i <= dirnum; ++i)
        for (int j = -dirnum; j <= dirnum; ++j) 
            directions.push_back(ipoint2(variation_factor*i, variation_factor*j));


    if (!res->grid(k)) res->new_grid(x_size, y_size, k);
    if (!res->grid(k1)) res->init_grid(k1);

    
    //library->get_similar_types(k1, nd->m, type_thresh, types);
    //titer = types.begin(); // type list has always size >= 1!
    //ntype = *titer;

    //node* cpart = parts1[ntype];

    vector<node*> s_nodes1;

    res->get_layer_nodes(s_nodes1, k1);

    if (s_nodes1.empty()) return false;

    // Add reconstruction edges in case of shape checking (since this is not thread-safe, it must be
    // done in advance)
    if (shape_check) {
        for (auto niter = s_nodes1.begin(); niter != s_nodes1.end(); ++niter) {
            //link_path(res, *niter, to_prev_layer, to_layer_0);
            res->recurse_and_link(*niter, to_prev_layer, to_layer_0);
        }
    }
    library->init_basic_part_numbers(k);

    //////node* n = *niter;
    //////layer1_data* nd = (layer1_data*)n->data;
    //////int ntype;
    //////int x_size_c = int_round(nd->x/contraction); // !floor!
    //////int y_size_c = int_round(nd->y/contraction); // !floor!

    bool result = false;
    vector<list<proj_data_7> > candidates(s_nodes1.size());

    cout << "|";

    #pragma omp parallel for
    for (int ni = 0; ni < (int)s_nodes1.size(); ++ni) {
        node* n = s_nodes1[ni];
        layer1_data* nd = (layer1_data*)n->data;
        int ntype = nd->m;

        if (nd->r(R_RESPONSE) < candidate_r_threshold || nd->r(G_RESPONSE) < candidate_g_threshold ||
                (!region.invalid() && !region.inside(nd->x, nd->y))) 
            continue;

        node* cpart = parts1[ntype];
        part_data* cpartd = (part_data*)cpart->data;
        node* nc = follow_center_link(res, n);
        layer1_data* ncd = (layer1_data*)nc->data;

        foreach_neighbor (cpart, to_center_back, cpn_iter) {
            node* pc = neighbor_node(cpn_iter);

            if (!allowed_parts.empty() && allowed_parts.find(((lib_data*)pc->data)->type) == allowed_parts.end())
                continue;

            part_data_2c* pced = (part_data_2c*)neighbor_edge_data(cpn_iter);
            part_data* pcd = (part_data*)pc->data;
            int bpn = library->get_basic_part_number(pc); // basic part number

            if (min_part_distance >= 0 && res->find_node_in_region(ipoint2(nd->x + pcd->cmx - pced->x, nd->y + pcd->cmy - pced->y), k, 
                min_part_distance, min_part_distance, pcd->type, 0) != nullptr) continue;

            int part_count = pc->count_neighbors(to_part);
            int part_size = part_count;
            int failnum = (int)((1.0 - get_thresh(pcd, RR_THRESH, realization_ratio_threshold))*part_count);
                    
            for (double factor = min_factor; factor < max_factor + 10E-6; factor += 0.2) { // make 0.2 parameter...
                for (int d = 0; d < (int)directions.size(); ++d) {
                    double bpr = 0; // count how many basic parts are realized
                    double r_sum = 0.0;
                    double g_val = (g_response_operation == 0) ? 0.0 : 1.0;
                    int part_count = 1.0;
                    double realization_ratio = 0.0;
                    //int realized_count = 0;
                    //int fail = 0;
                    layer1_result::sp_result_t* conn = new layer1_result::sp_result_t();
                    const ipoint2& dir = directions[d];
                    response_map r;
                    vector<inference_data> idata;
                    
                    foreach_neighbor(pc, to_part, iter2) {
                        node* p = neighbor_node(iter2);
                        lib_data* pd = (lib_data*)p->data;
                        part_data_2* ed = (part_data_2*)neighbor_edge_data(iter2);
                        layer1_result::sp_result_t newconn;
                        double c = library->contraction(pd->layer + 2, k1 + 1);
                        dnpair max;

                        max = res->schur_product_max_all_ol(
                            (int)(c*(nd->x + factor*(ed->x - pced->x + dir.x))), 
                            (int)(c*(nd->y + factor*(ed->y - pced->y + dir.y))),
                            ed->x, ed->y, ed->index, 
                            ed->distr, ed->app, pd->layer, &layer1_result::sgspf, convolution_threshold, newconn);
                        if (max.first >= convolution_threshold) {
                            layer1_data* nnd = (layer1_data*)max.second->data;
                            if (convolution_link_threshold >= 1.0) 
                                conn->push_back(layer1_result::sp_result_data_t(max.second, ed->index, max.first));
                            else
                                filter_links(*conn, newconn, convolution_link_threshold * max.first);
                            r_sum += nnd->r(R_RESPONSE);
                            bpr += library->get_basic_part_number(p)*nnd->r(RR_RESPONSE);
                            if (g_response_operation == 0) g_val += max.first; 
                            else g_val *= max.first;

                            if (shape_check && !scmap.empty()) {
                                idata.push_back(inference_data(ed->index, max.second));
                            }
                        }
                        ++part_count;
                    }
                    // Calculate RR_RESPONSE
					realization_ratio = bpr/bpn;

                    //realization_ratio = (part_count == 0) ? 1.0 : (double)realized_count/((double)part_count);
                    r_sum += nd->r(R_RESPONSE);
                    r_sum /= part_count;
                    r_sum = ::pow(r_sum, r_response_pow);
                    g_val = ::pow(g_val * realization_ratio, g_response_pow);

                    if (r_sum > get_thresh(pcd, R_THRESH, r_response_threshold) && 
                            realization_ratio >= get_thresh(pcd, RR_THRESH, realization_ratio_threshold) &&
                            g_val >= get_thresh(pcd, G_THRESH, g_response_threshold)) {
                        int type = pcd->type;
                        vs_part_data* vspcd = dynamic_cast<vs_part_data*>(pc->data);
                        bool ok = true;
                        double benergy = 0.0;
                        double scdistance = 0.0;
                        double maxdistance = 0.0;
                        vector<node*> missnodes;

                        if (shape_check && !scmap.empty()) {
                            if (vspcd != nullptr && !vspcd->pcad.mean.empty()) {
                            //if (vspcd != nullptr && vspcd->svmd.svm != nullptr) {
                                vector<pair<int, ipoint2> > ipts;
                                pair<double, dpoint2> tsdata = get_node_geo(ipts, res, idata);
                                vector<dpoint2> misspos;
                                        
                                ddpair result = check_pca_geometry_p(misspos, ipts, vspcd->pcad, pc);
                                //ddpair result = check_svm_geometry_p(ipts, vspcd->svmd, pc);
                                scdistance = result.first;
                                ok = scdistance <= get_thresh(pcd, S_THRESH, s_response_threshold);
                                if (ok && link_missing_support && !misspos.empty()) {
                                    vector<ipoint2> missposi;

                                    for (auto mpiter = misspos.begin(); mpiter != misspos.end(); ++mpiter) {
                                        missposi.push_back(ipoint2(mpiter->x*tsdata.first/100.0 + tsdata.second.x,
                                            mpiter->y*tsdata.first/100.0 + tsdata.second.y));
                                    }
                                    get_support_nodes(missnodes, res, missposi, max(int_round(tsdata.first * 0.05), 3));
                                }
                            } else {
                                vector<pair<int, ipoint2> > ipts;
                                pair<double, dpoint2> tsdata = get_node_geo(ipts, res, idata);

                                scdistance = check_geometry2(ipts, pc);
                                ok = scdistance <= get_thresh(pcd, S_THRESH, s_response_threshold);
                            }
                        }
                        /*
                        if (vspcd != nullptr && vspcd->svmt.svm != nullptr) {
                            double f = check_svmt(vspcd->svmt, g_val, realization_ratio,
                                scdistance, part_size, idata, true);
                            //ok = f < 0.5;
                            ok = f > -0.5;
                            benergy = f + 0.51;
                        }*/

                        if (ok) {
                            r.set_response(R_RESPONSE, r_sum);
                            r.set_response(RR_RESPONSE, realization_ratio);
                            r.set_response(G_RESPONSE, g_val); 
                            r.set_response(S_RESPONSE, scdistance); 
                            r.set_response(X_RESPONSE, benergy);
                            candidates[ni].push_back(proj_data_7(r, type, pcd->cmx - pced->x, pcd->cmy - pced->y, conn, missnodes));
                        } else 
                            delete conn;
                    } else
                        delete conn;
                }
                    
            }
            
        } // end of "foreach 'to_center_back' neighbor"
    }

    cout << "|";

    // Prepare for "fast insertion"
    map<ipoint2, parallel_candidates_t> cooccurrent;   // new position -> (candidate index, proj_data*)
    vector<double> values;
    double value_threshold;

    if (new_positions_threshold < INT_MAX)
        values.reserve(candidates.size()*2);
    for (int ci = 0; ci < (int)candidates.size(); ++ci) {
        layer1_data* nd = (layer1_data*)(s_nodes1[ci]->data);

        for (auto citer = candidates[ci].begin(); citer != candidates[ci].end(); ++citer) {
            int new_x = int_round((nd->x + citer->x)/contraction); // !floor!
            int new_y = int_round((nd->y + citer->y)/contraction); // !floor!

            cooccurrent[ipoint2(new_x, new_y)].push_back(pair<int, proj_data_7*>(ci, &(*citer)));
            if (new_positions_threshold < INT_MAX) 
                values.push_back(citer->val());
        }
    }

    if (values.size() < new_positions_threshold) value_threshold = 0.0;
    else {
        sort(values.begin(), values.end(), greater<double>());
        value_threshold = values[new_positions_threshold - 1] - 1E-6;
    }

    vector<parallel_candidates_t> par_candidates;

    for (auto coociter = cooccurrent.begin(); coociter != cooccurrent.end(); ++coociter) 
        par_candidates.push_back(coociter->second);

    // Result of the following block (shape checking): 
    //  - vector of new 'layer1_result's and 
    //  - vector of new edges

    vector<list<new_edges_data_t> > new_edges(par_candidates.size());
    vector<map<int, pair<layer1_data*, new_edges_data_t*> > > new_data(par_candidates.size());
    
    for (int pci = 0; pci < (int)par_candidates.size(); ++pci) {

        for (auto citer = par_candidates[pci].begin(); citer != par_candidates[pci].end(); ++citer) {
            int ci = citer->first;
            node* n = s_nodes1[ci];
            layer1_data* nd = (layer1_data*)n->data;
            const proj_data_7& best = *citer->second;

            new_edges[pci].push_back(vector<new_edge_data_t>());

            itriple new_edges_info = get_new_edges(new_edges[pci].back(), k1, best.d->begin(), best.d->end());  

            get_new_support_edges(new_edges[pci].back(), best.msupp);
            if (citer->second->val() > value_threshold && (!add_reconstruction_edges || new_edges_info.second == 0 || 
                    (double)new_edges_info.third/new_edges_info.second <= hypo_ratio_threshold)) {
                int new_x = int_round((nd->x + best.x)/contraction); // !floor!
                int new_y = int_round((nd->y + best.y)/contraction); // !floor!

                if (new_x >= 0 && new_x < x_size && new_y >= 0 && new_y < y_size) {
                    layer1_data* ly1d = new layer1_data(best.r, best.type);

                    ly1d->x = new_x;
                    ly1d->y = new_y;
                    ly1d->z = k;

                    bool exists = new_data[pci].find(ly1d->m) != new_data[pci].end();
                    pair<layer1_data*, new_edges_data_t*>& depair = new_data[pci][ly1d->m];

                    if (!exists) {
                        depair.first = ly1d;
                        depair.second = &new_edges[pci].back();
                    } else if (depair.first->vval() >= ly1d->vval()) {
                        delete ly1d;
                    } else {
                        delete depair.first;
                        depair.first = ly1d;
                        depair.second = &new_edges[pci].back();
                    }

                    if (!n->is_attr_set(HAS_NEXT_LAYER)) {  // ?
                        res->inc_covered(k1);
                        n->set_attr(HAS_NEXT_LAYER);
                    }

                }
            }
            delete best.d;
        }
    }

    // We insert the candidates into the graph; at each position we insert nodes 
    // sorted by best value (G_RESPONSE)

    for (int ndi = 0; ndi < (int)new_data.size(); ++ndi) {
        vector<pair<layer1_data*, new_edges_data_t*> > v;

        v.reserve(new_data[ndi].size());
        for (auto nditer = new_data[ndi].begin(); nditer != new_data[ndi].end(); ++nditer) 
            v.push_back(nditer->second);
        sort(v.begin(), v.end(), [](const pair<layer1_data*, new_edges_data_t*>& a, const pair<layer1_data*, new_edges_data_t*>& b) 
            { return a.first->vval() < b.first->vval(); });
        for (auto viter = v.begin(); viter != v.end(); ++viter) {
            node* newn = res->add_grid_node(viter->first, viter->first->x, viter->first->y, k);

            add_edges(res, newn, *viter->second);
        }
    }
        

    double quot = res->cover_quotient(k1);

    count++;
    quot_sum += quot;
    if (quot < min_cover_quot) min_cover_quot = quot;

    if (copy_prev_layer) copy_layer(res, k);

    if (inhibit_layer_response >= 0) {
        inhibit_layer(res, k, inhibit_layer_response, inhibit_layer_max, inhibit_layer_thresh);
        if (inhibit_layer_delete) 
            res->delete_nodes(HIDDEN_NODE_ATTR);
    }

    //cout << "Best dist: " << bestshapedist << " of type: " << besttype << endl;

    inhibit_result(res, k);
    //res->delete_edges_leq(to_layer_0, k - 1);
    res->delete_edges(to_layer_0);
    return result;
}

// Add node(s) to layer 'n->z + 1' with 'n' as a center
void layern_creator::add_layer8(layer1_result* res, node* n, int endk)
{
	throw new_libhop_exception("Error: add_layer8 no longer supported");
	/*
    layer1_data* nd = (layer1_data*)n->data;
    int k1 = nd->z;
    int k = k1 + 1;

    if (endk < k) return;

    node* cpart = library->parts[k1][nd->m];
    int to_center_back = atom("lyrCenterBack");
    int to_part = atom("lyrSrc");
    if (!res->grid(k)) res->new_grid(res->x_size(k1), res->y_size(k1), k);
    int back_count = cpart->count_neighbors(to_center_back);

    foreach_neighbor (cpart, to_center_back, iter) {
        node* pc = neighbor_node(iter);
        part_data_2c* pced = (part_data_2c*)neighbor_edge_data(iter);
        int newx = -1, newy = -1;
        int newm = -1;
        double newval;
        response_map newr;
        layer1_result::sp_result_t conn;
        
        if (pc->is_attr_set(R_PART_ATTR)) {
            rpart_data* pcd = (rpart_data*)pc->data;

            newx = nd->x - pced->x;
            newy = nd->y - pced->y;

            if (res->node_at(newx, newy, endk, 1, 1) != nullptr) continue;
            if (res->find_node_in_region(ipoint2(newx, newy), k, 1, 1, pcd->type, 0) != nullptr) continue;

            if (back_count == 1 && cpart->is_attr_set(R_PART_ATTR)) {
                // Do nothing...
                newm = pcd->type;
                newval = nd->val(); 
            } else {
                int part_count = pc->count_neighbors(to_part);
                int failnum = max<int>(0, part_count - (int)(get_thresh(pcd, RR_THRESH, realization_ratio_threshold) * pcd->minsp));
                int fail = 0;

                foreach_neighbor(pc, to_part, iter2) {
                    node* p = neighbor_node(iter2);
                    part_data* pd = (part_data*)p->data;
                    part_data_2r* ed = (part_data_2r*)neighbor_edge_data(iter2);
                    map<int, double> appmap;
                    //int foundn;
                    node* foundn;

                    library->get_similar_types(pd->layer, pd->type, pcd->layer, pcd->type, ed->x, ed->y, 1.0, appmap);
                    appmap.insert(pair<int, double>(pd->type, 1.0));

                    foundn = res->find_node_in_rect(pd->layer, 
                        ed->rect + ipoint2(nd->x + ed->x - pced->x, nd->y + ed->y - pced->y),
                        appmap);
                    //foundn = res->find_nodes_in_rect(conn, pd->layer, 
                    //    ed->rect + ipoint2(nd->x + ed->x - pced->x, nd->y + ed->y - pced->y), 
                    //    typeset  set the parameter 'max' to 1 if we need performance! );

                    if (foundn == nullptr) {
                        ++fail;
                        if (fail > failnum) break;
                    }
                }
                if (fail <= failnum) { 
                    newm = pcd->type;
                    newval = nd->val(); 
                }
            }
            
        } else if (pc->is_attr_set(OBJ_PART_ATTR)) {
            part_data* pcd = (part_data*)pc->data;

            newx = nd->x + pcd->cmx - pced->x;
            newy = nd->y + pcd->cmy - pced->y;

            if (res->node_at(newx, newy, endk, 1, 1) != nullptr) continue;
            if (res->find_node_in_region(ipoint2(newx, newy), k, 1, 1, pcd->type, 0) != nullptr) continue;

            int part_count = pc->count_neighbors(to_part);
            int failnum = max<int>(0, (int)((1.0 - get_thresh(pcd, RR_THRESH, realization_ratio_threshold))*part_count));
                    
            for (double factor = min_factor; factor < max_factor + 10E-6; factor += 0.2) { // make 0.2 parameter...
                double realization_ratio = 0.0;
                int max_thresh_count = 0;
                int fail = 0;

                newval = 0.0;
                foreach_neighbor(pc, to_part, iter2) {
                    node* p = neighbor_node(iter2);
                    lib_data* pd = (lib_data*)p->data;
                    part_data_2* ed = (part_data_2*)neighbor_edge_data(iter2);
                    map<int, double> appmap;
                    dnpair max;

                    library->get_similar_types(pd->layer, pd->type, pcd->layer, pcd->type, ed->x, ed->y, type_thresh, appmap);
                    appmap.insert(pair<int, double>(pd->type, 1.0));

                    max = res->schur_product_max_all(nd->x, nd->y, 
                        (int)(factor*(ed->x - pced->x)), (int)(factor*(ed->y - pced->y)), ed->index,
                        ed->distr, appmap, pd->layer, 1.0, R_RESPONSE, candidate_r_threshold, conn);

                    if (max.first < convolution_threshold) ++fail;
                    else { ++max_thresh_count; newval += max.first; }
                    if (fail > failnum) break;
                }
                realization_ratio = (part_count == 0) ? 1.0 : (double)max_thresh_count/((double)part_count);
                newval += nd->val();
                newval /= (part_count + 1);
                newval *= ::pow(newval, r_response_pow);

                if (newval > get_thresh(pcd, R_THRESH, r_response_threshold) && 
                        realization_ratio >= get_thresh(pcd, RR_THRESH, realization_ratio_threshold)) {
                    newm = pcd->type;
                    newr.set_response(R_RESPONSE, newval);
                    newr.set_response(RR_RESPONSE, realization_ratio);
                    break;
                }
            }
        }

        if (newm >= 0 && newx >= 0 && newx < res->x_size(k) && newy >= 0 && newy < res->y_size(k)) {
            // Add node
            vector<new_edge_data_t> new_edges;
        
            get_new_edges(new_edges, k1, conn.begin(), conn.end());  

            layer1_data* ly1d = new layer1_data(newr, newm);
            pair<node*, node**> fpair = res->find_at(ly1d, newx, newy, k);
            node* newn = fpair.first;

            if (newn == nullptr) {
                newn = res->add_grid_node(ly1d, newx, newy, k);
            } else {
                if (!res->data_greater(ly1d, (layer1_data*)newn->data)) newn = nullptr;
                else {
                    res->remove_grid_node(fpair);
                    ((layer1_data*)newn->data)->r = ly1d->r;
                    newn->clear_neighbors();
                    res->add_grid_node(newn);
                } 
                delete ly1d;
            }
            if (newn != nullptr) {
                add_edges(res, newn, new_edges);
                res->add_edge(newn, n, atom("toPrevLayerI"));

                if (!n->is_attr_set(HAS_NEXT_LAYER)) {
                    res->inc_covered(k1);
                    n->set_attr(HAS_NEXT_LAYER);
                }
                creators[k + 1]->add_layer8(res, newn, endk);
            }
        }

    }
*/
}

void layern_creator::add_layer8(layer1_result* res, int layer, int stoplayer)
{
	throw new_libhop_exception("Error: add_layer8 no longer supported");
	/*
    int k = layer - 1;
    int k1 = k - 1;
    int ell = stoplayer - 1;

    if (res->max_layer_index() < k1 || res->shape_nodes[k1].empty() || 
        library->max_layer_index() < ell || ell < k) return;

    while ((int)res->shape_nodes.size() <= ell) res->shape_nodes.push_back(vector<node*>());
    while ((int)res->shape_nodes_inhib.size() <= ell) res->shape_nodes_inhib.push_back(vector<node*>());
    while ((int)res->info.size() <= ell) res->info.push_back(layer_info());

    vector<node*>& s_nodes1 = res->shape_nodes[k1];
    
    for (vector<node*>::iterator niter = s_nodes1.begin(); niter != s_nodes1.end(); ++niter) {
        node* n = *niter;

        while (n != nullptr) {
            add_layer8(res, n, ell);
            n = ((layer1_data*)n->data)->next;
        }

    }
    for (int i = k; i <= ell; ++i)
        inhibit_result(res, i);
		*/
}

// Adds hypothetical nodes
// - res = graph
// - layer = name of the layer (i.e. index + 1) to add the nodes
// Remarks:
void layern_creator::add_hypothetical_nodes(layer1_result* res, int layer)
{
    int k1 = layer - 1;

    if (k1 < 0 || k1 >= (int)res->shape_nodes.size()) return;
    if (!res->grid(k1)) res->init_grid(k1);

    vector<node*>& parts1 = library->parts[k1];
    vector<node*>& s_nodes1 = res->shape_nodes[k1];
    int to_center = atom("lyrCenterBack").get_index();
    int to_part = atom("lyrSrc").get_index();
    int to_prev_layer = atom("toPrevLayer").get_index();
    int to_forbidden = atom("lyrForbidden").get_index();
    int maxx = res->x_size(k1);
    int maxy = res->y_size(k1);
    set<ipoint2> hypo_positions;
    list<ipoint2> new_hypo_positions;

    for (vector<node*>::iterator niter = s_nodes1.begin(); niter != s_nodes1.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;
        node* cpart = parts1[nd->m];
        part_data* cpartd = (part_data*)cpart->data;
        priority_queue<proj_data_base> candidates;
        multimap<int, edge_pair>::iterator citer, citerend;
        
        foreach_neighbor(cpart, to_center, iter) {
            node* pc = iter->second.first;
            part_data* pcd = (part_data*)pc->data;
            double max_sum = 0.0;
            double realization_ratio = 0.0;
            int max_thresh_count = 0;
            int part_count = 0;
            list<edge_pair*> unrealized_parts;

            foreach_neighbor(pc, to_part, iter2) {
                node* p = iter2->second.first;
                part_data* pd = (part_data*)p->data;
                part_data_2* ed = (part_data_2*)iter2->second.second;
                dnpair max;

                max = res->schur_product_max_all(nd->x + ed->x, nd->y + ed->y, ed->distr, pd->type, 
                    k1, R_RESPONSE, candidate_r_threshold);
                if (max.first < convolution_threshold) 
                    unrealized_parts.push_back(&(iter2->second));
                ++part_count;
            }

            double ratio = (double)unrealized_parts.size()/(part_count + 1);

            if (!unrealized_parts.empty() && (ratio <= hypo_unrealized_threshold)) {
                for (list<edge_pair*>::iterator epi = unrealized_parts.begin(); epi != unrealized_parts.end(); ++epi) {
                    edge_pair* ep = *epi;
                    int m = ((part_data*)ep->first->data)->type;
                    part_data_2* ed = (part_data_2*)ep->second;
                    int newx = nd->x + ed->x;
                    int newy = nd->y + ed->y;

                    if (newx < 0 || newx >= maxx || newy < 0 || newy >= maxy) 
                        continue;

                    layer1_data_h* ly1d = new layer1_data_h(hypo_start - 1.0/(part_count + 1.0), m, pcd->type, ed->x, ed->y);     
                    pair<node*, node**> fpair = res->find_at(ly1d, newx, newy, k1);
                    node* newn = fpair.first;
                    vector<node*> nbhood;

                    if (newn != nullptr && !newn->is_attr_set(HYPO_NODE_CANDIDATE_ATTR)) {        // though this would be strange
                        delete ly1d;
                        continue; 
                    }

                    res->find_nodes_in_region(nbhood, newx, newy, k1, hypo_nbhood, hypo_nbhood, m, 0);    

                    for (vector<node*>::iterator nbiter = nbhood.begin(); nbiter != nbhood.end(); ++nbiter) {
                        if ((*nbiter)->is_attr_set(HYPO_NODE_CANDIDATE_ATTR)) {
                            layer1_data_h* hnd = (layer1_data_h*)(*nbiter)->data;   
                            hnd->r.set_response(R_RESPONSE, min(hnd->val() * hypo_factor, 1.0));
                            hnd->add_voter(pcd->type, ed->x, ed->y);
                        }
                    }
                    if (newn != nullptr) 
                        delete ly1d;
                    else {
                        if (res->node_at(newx, newy, k1) == nullptr) 
                            new_hypo_positions.push_back(ipoint2(newx, newy));      
                        ly1d->r.set_response(R_RESPONSE, min(ly1d->val() * pow(hypo_factor, (double)nbhood.size()), 1.0));   
                        newn = res->add_grid_node(ly1d, newx, newy, k1);
                        newn->set_attr(HYPO_NODE_CANDIDATE_ATTR);
                        hypo_positions.insert(ipoint2(newx, newy));
                    }
                }
            }
        }
    }
    
    h_remove_predicate r_predicate(hypo_val_threshold, hypo_voters_threshold); 

    for (set<ipoint2>::iterator iter = hypo_positions.begin(); iter != hypo_positions.end(); ++iter) {
        const ipoint2& p = *iter;

        res->sort_at(p.x, p.y, k1);
        res->img_graph::remove_grid_nodes(p.x, p.y, k1, r_predicate);
        
    }
    int newhp = 0;
    for (list<ipoint2>::iterator iter = new_hypo_positions.begin(); iter != new_hypo_positions.end(); ++iter) {
        ipoint2& p = *iter;
        node* n = res->node_at(p.x, p.y, k1);
        
        if (n != nullptr) { s_nodes1.push_back(n); ++newhp; }
    }
    cout << '(' << newhp << " hyponodes added) ";
    res->delete_isolated_nodes(NODE_DELETED_ATTR);
}

void layern_creator::add_edges(node* newn, vector<new_edge_data_t>& new_edges /*, unsigned attr = 0*/)
{
    for (vector<new_edge_data_t>::iterator iter = new_edges.begin(); iter != new_edges.end(); ++iter) {
        if (!add_edge_names || iter->spindex == invalid_edge_name) {
            graph::add_edge(newn, iter->n, iter->name);
            //iter->n->set_attr(attr);
        } else {
			graph::add_edge_2(newn, iter->n, new edge_data_name(iter->spindex, iter->r), iter->name);
            //iter->n->set_attr(attr);
        }
    }
}

void layern_creator::add_edges(layer1_result* res, node* newn, vector<new_edge_data_t>& new_edges /*, unsigned attr = 0*/)
{
    for (vector<new_edge_data_t>::iterator iter = new_edges.begin(); iter != new_edges.end(); ++iter) {
        if (!add_edge_names || iter->spindex == invalid_edge_name) {
            res->add_edge(newn, iter->n, iter->name);
            //iter->n->set_attr(attr);
        } else {
            res->add_edge_2(newn, iter->n, new edge_data_name(iter->spindex, iter->r), iter->name);
            //iter->n->set_attr(attr);
        }
    }
}


void layern_creator::read_from_stream(istreamer& is) {
	streamable* ptr;
	is.read(ptr); library = (part_lib*)ptr;
	is_library_my = true;
		
	is.read(layer_contraction);
	is.read(manual_thresholds);
	is.read(threshold_factor);
	is.read(candidate_r_threshold);
	is.read(candidate_r_threshold_percent);
	is.read(candidate_g_threshold);
	is.read(candidate_g_threshold_percent);
	is.read(r_response_threshold);
	is.read(g_response_threshold);
	is.read(s_response_threshold);
	is.read(x_response_threshold);
	is.read(g_response_threshold_percent);
	is.read(identity_g_response);
	is.read(simple_g_response);
	is.read(ignore_g_distribution);
	is.read(g_response_operation);
	is.read(g_response_var_factor);
	is.read(type_thresh);
	is.read(convolution_threshold);
	is.read(convolution_link_threshold);
	is.read(proj_max_forb_threshold);
	is.read(forb_quot_threshold);
	is.read(r_response_pow);
	is.read(g_response_pow);
	is.read(realization_ratio_threshold);
	is.read(shape_check);
	is.read(link_missing_support);
	is.read(continuity_factor);
	is.read(normalize_histogram);
	is.read(copy_prev_layer);
	is.read(texture_parts);
	is.read(texture_radius);
	is.read(ignore_texture);
	is.read(min_factor);
	is.read(max_factor);
	is.read(reconstruction_type);
	is.read(add_reconstruction_edges);
	is.read(add_activation_edges);
	is.read(hypo_unrealized_threshold);
	is.read(hypo_ratio_threshold);
	is.read(null_tolerance_limit);
	is.read(positive_tolerance_threshold);
	is.read(projection_radius);
	is.read(reconstruction_factor);
	is.read(rec_null_tolerance_limit);
	is.read(rec_tolerance_threshold);
	is.read(tolerance_radius);
	is.read(min_part_distance);
	is.read(depth_first_search);
	is.read(hypo_start);
	is.read(hypo_nbhood);
	is.read(hypo_factor);
	is.read(hypo_val_threshold);
	is.read(hypo_voters_threshold);
	is.read(add_edge_names);
	is.read(new_positions_threshold);
	is.read(same_position_threshold);
	is.read(allowed_parts);
	is.read(variation_dimension);
	is.read(variation_factor);
	is.read(inhibit_layer_response);
	is.read(inhibit_layer_max);
	is.read(inhibit_layer_thresh);
	is.read(inhibit_layer_delete);
	is.read(strict_svm_checking);
	is.read(scale_merge);

	is.read(min_cover_quot);
	is.read(quot_sum);
	is.read(count);

	use_opencl = false;
	opencl_verify_result = true;

}
void layern_creator::write_to_stream(ostreamer& os) {
	os.write(library);

	os.write(layer_contraction);
	os.write(manual_thresholds);
	os.write(threshold_factor);
	os.write(candidate_r_threshold);
	os.write(candidate_r_threshold_percent);
	os.write(candidate_g_threshold);
	os.write(candidate_g_threshold_percent);
	os.write(r_response_threshold);
	os.write(g_response_threshold);
	os.write(s_response_threshold);
	os.write(x_response_threshold);
	os.write(g_response_threshold_percent);
	os.write(identity_g_response);
	os.write(simple_g_response);
	os.write(ignore_g_distribution);
	os.write(g_response_operation);
	os.write(g_response_var_factor);
	os.write(type_thresh);
	os.write(convolution_threshold);
	os.write(convolution_link_threshold);
	os.write(proj_max_forb_threshold);
	os.write(forb_quot_threshold);
	os.write(r_response_pow);
	os.write(g_response_pow);
	os.write(realization_ratio_threshold);
	os.write(shape_check);
	os.write(link_missing_support);
	os.write(continuity_factor);
	os.write(normalize_histogram);
	os.write(copy_prev_layer);
	os.write(texture_parts);
	os.write(texture_radius);
	os.write(ignore_texture);
	os.write(min_factor);
	os.write(max_factor);
	os.write(reconstruction_type);
	os.write(add_reconstruction_edges);
	os.write(add_activation_edges);
	os.write(hypo_unrealized_threshold);
	os.write(hypo_ratio_threshold);
	os.write(null_tolerance_limit);
	os.write(positive_tolerance_threshold);
	os.write(projection_radius);
	os.write(reconstruction_factor);
	os.write(rec_null_tolerance_limit);
	os.write(rec_tolerance_threshold);
	os.write(tolerance_radius);
	os.write(min_part_distance);
	os.write(depth_first_search);
	os.write(hypo_start);
	os.write(hypo_nbhood);
	os.write(hypo_factor);
	os.write(hypo_val_threshold);
	os.write(hypo_voters_threshold);
	os.write(add_edge_names);
	os.write(new_positions_threshold);
	os.write(same_position_threshold);
	os.write(allowed_parts);
	os.write(variation_dimension);
	os.write(variation_factor);
	os.write(inhibit_layer_response);
	os.write(inhibit_layer_max);
	os.write(inhibit_layer_thresh);
	os.write(inhibit_layer_delete);
	os.write(strict_svm_checking);
	os.write(scale_merge);
	os.write(min_cover_quot);
	os.write(quot_sum);
	os.write(count);       
}