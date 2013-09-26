// test.cpp : Defines the entry point for the console application.
//


#ifndef _WIN32_WINNT		// Allow use of features specific to Windows XP or later.                   
#define _WIN32_WINNT 0x0501	// Change this to the appropriate value to target other versions of Windows.
#endif						

#include <stdio.h>
//#include <tchar.h>

#if defined WIN32 | defined WIN64
#include <windows.h>
#include <psapi.h>
#include <crtdbg.h>
#endif

// TODO: reference additional headers your program requires here

#include "layers/layer_1.h"




// #include <iostream>
// #include <ctime>
// #include "../img/img.h"
#include "layers/optimization.h"
#include "layers/initialization.h"


#include <ctime>

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <sys/stat.h>
#include <hash_map>
#include <queue>
#include "utils/ocv.h"


//#define _SAVEMAPS
#define PRINT_INFO(_info) cerr << _info << endl;
#define PRINT_DOT() cerr << '.';
//#define PRINT_INFO(_info) 
//#define PRINT_DOT()
//#define LOG_SP(_info) log_stream << _info << endl; log_stream.flush();
#define LOG_SP(_info) cerr << _info << endl; 
//#define LOG_SP(_info)


using namespace std;


// constants
///////////////////////////////////////////////////////////////////////////////

// layer1_data_ms

//class layer1_data_ms : public layer1_data {
//    CLASS_ALLOCATION();
//public:
//    set<int> mset;
//
//    layer1_data_ms() : layer1_data(), mset() { }
//    layer1_data_ms(double pval, int pm, const set<int>& pset) : layer1_data(pval, pm), mset(pset) { }
//    layer1_data_ms(const layer1_data_ms& d) : layer1_data(d), mset(d.mset) { }
//};
//
//DEFINE_CLASS_ALLOCATOR_2(layer1_data_ms, 1024);



// learning
///////////////////////////////////////////////////////////////////////////////

class stat_learning {
public:
    typedef matrix<int> nbhood_t;

protected:
    //typedef set<int> region_t;

    matrix<nbhood_t> nbhoods;
    //vector<region_t*> c_regions;
    //vector<region_t*> d_regions;

public:
    stat_learning(int cs, int ds, int nbsize);

    void reset();
    void reset(int cs, int ds, int nbsize);
    //void reset_regions(const vector<region_t*>& cr, const vector<region_t*>& cr

    int c_size() const { return (int)nbhoods.width; }   // "center" dimension
    int d_size() const { return (int)nbhoods.height; }  // the "other" dimension
    int nb_size() const { return (nbhoods.size() == 0) ? 0 : (int)nbhoods[0].width; }
    nbhood_t& at(int c, int d) { return nbhoods(c, d); }
    virtual string name_at(int c, int d) { return string("nbhood_") + c + string("-") + d; }

    virtual int update(int c, int d, int x, int y);
};

class library_path_walker : public path_walker {
public:
    typedef vector<pair<int, ipoint2> > library_path_t;
    typedef list<library_path_t> result_t;

    int llimit;

    library_path_walker(int ll) : path_walker(), llimit(ll) { }

    bool check_neighbor(container_t& neighbors, const container_t::value_type& n, edge_data* ed, node* nn)
    {
        if (((part_data*)(nn->data))->layer < llimit) return false;
        else return path_walker::check_neighbor(neighbors, n, ed, nn);
    }

    void get_result(result_t& tresult)
    {
        tresult.clear();
        for (list<path_t>::iterator iter = result.begin(); iter != result.end(); ++iter) {
            path_t& p = *iter;
            
            tresult.push_back(library_path_t());
            for (path_t::iterator piter = p.begin(); piter != p.end(); ++piter) {
                part_data* pd = (part_data*)piter->first->data;
                part_data_2* ed = (part_data_2*)piter->second;

                if (ed == nullptr) tresult.back().push_back(pair<int, ipoint2>(pd->type, ipoint2::zero));
                else tresult.back().push_back(pair<int, ipoint2>(pd->type, ipoint2(ed->x, ed->y)));
            }
            
        }
    }
};

class layer1_result_walker : public path_walker {
public:
    typedef vector<pair<int, ipoint2> > layer1_result_path_t;
    typedef multimap<node*, layer1_result_path_t> result_t;

protected:
    int llimit;

public:

    layer1_result_walker(int ll) : path_walker(), llimit(ll) { }

    bool check_neighbor(container_t& neighbors, const container_t::value_type& n, edge_data* ed, node* nn)
    {
        if (((part_data*)(nn->data))->layer < llimit) return false;
        else return path_walker::check_neighbor(neighbors, n, ed, nn);
    }

    void insert_leaf(const container_t::value_type& n) 
    { 
        if (((layer1_data*)n.back().first->data)->z == llimit)
            result.push_back(n); 
    }

    void get_result(result_t& tresult)
    {
        tresult.clear();
        for (list<path_t>::iterator iter = result.begin(); iter != result.end(); ++iter) {
            path_t& p = *iter;
            result_t::iterator iiter;
            
            iiter = tresult.insert(result_t::value_type(p.back().first, layer1_result_path_t()));
            for (path_t::iterator piter = p.begin(); piter != p.end(); ++piter) {
                layer1_data* pd = (layer1_data*)piter->first->data;
                edge_data_tw<ipoint2>* ed = (edge_data_tw<ipoint2>*)piter->second;

                if (ed == nullptr) iiter->second.push_back(pair<int, ipoint2>(pd->m, ipoint2::zero));
                else iiter->second.push_back(pair<int, ipoint2>(pd->m, ed->data));
            }
            
        }
    }
};

class part_stat_learning : public stat_learning {
public:
    // path = ((node, edge name (i.e. coordinates)),...) = (...(v_i, e_i),...)
    //        where edge e_i = v_i-1 -> v_i and e_0 = (0, 0)
    typedef vector<pair<int, ipoint2> > path_t;  
    typedef list<path_t> path_collection_t; 
    typedef map<path_t, int> path_map_t;
    typedef vector<int> path_vector_t;

protected:
    // center layer, "other" layer, projection layer; player must be <= dlayer, clayer
    int clayer, dlayer, player;  
    path_collection_t cpaths, dpaths;
    path_map_t cpmap, dpmap;

    int nbsize;
    double center_val_threshold;
public:
    part_stat_learning(const config_dictionary& cfg);

    virtual string name_at(int c, int d);
    int update(layer1_result* res);
protected:
    void prepare_for_update(layer1_result* res);
    void init(const config_dictionary& cfg);
    void init_from_library(part_lib* newlib);
};

class stat_learning_exception : public std::exception {
    stat_learning_exception() : std::exception() { }
};


// part_stat_learning
///////////////////////

void part_stat_learning::init(const config_dictionary& cfg)
{
    // init layers
    clayer = cfg.get_value_int("c_layer", 0);
    dlayer = cfg.get_value_int("d_layer", 0);
    player = cfg.get_value_int("p_layer", 0);
    nbsize = cfg.get_value_int("nb_size", 15);

    // init library
    string libname;
    part_lib* library;

    cfg.get_value(libname, "library", true);
    read_library(libname, library);

    init_from_library(library);

    if (library != nullptr) delete library;
}

// Parse all possible paths from s(ource)layer to p(rojection)layer 
// wrt the given library.
void get_library_paths(library_path_walker::result_t& result, part_lib* library, int slayer, int player)
{
    if (player < 0 || slayer < player || slayer > library->max_layer_index()) return;

    library_path_walker walker(player);
    library_path_walker::container_t start;
    vector<int> edge_names(2);

    edge_names[0] = atom("lyrSrc");
    edge_names[1] = atom("lyrCenter");
    path_walker::add_to_container(start, library->parts[slayer].begin(), library->parts[slayer].end());
    library->recurse2(walker, start, edge_names);
    walker.get_result(result);

}

void part_stat_learning::init_from_library(part_lib* library)
{
    if (library == nullptr) return;

    // make all possible paths from clayer -> player and from dlayer -> player
    // library_path_walker::result_t cpaths, dpaths;

    get_library_paths(cpaths, library, clayer, player);
    get_library_paths(dpaths, library, dlayer, player);

    reset((int)cpaths.size(), (int)dpaths.size(), nbsize);

    int count;
    path_collection_t::iterator iter;

    for (iter = cpaths.begin(), count = 0; iter != cpaths.end(); ++iter, ++count) 
        cpmap.insert(path_map_t::value_type(*iter, count));
    for (iter = dpaths.begin(), count = 0; iter != dpaths.end(); ++iter, ++count) 
        dpmap.insert(path_map_t::value_type(*iter, count));

}

// returns a |dlayer| x |dlayer| matrix of 
void make_region_matrix(matrix<region2*>& matrix, layer1_result* res, int dlayer, int slayer)
{
    vector<img*> masks;
    vector<region2> regions;
    vector<ipoint2> centers;
    int name = atom("toPrevLayer").get_index();

    res->get_masks(masks);
    regions.resize(masks.size());
    centers.resize(masks.size());
    
    for (size_t i = 0; i < masks.size(); ++i) {
        img* m = masks[i];

        regions[i].assign(*m, 0.1);
        centers[i].x = m->width/2;
        centers[i].y = m->height/2;
        delete m;
    }

    vector<node*>& nodes = res->shape_nodes[slayer];
    double factor = (double)res->x_size(dlayer)/res->x_size(0);
    double factor2 = (double)res->x_size(dlayer)/res->x_size(slayer);

    matrix.resize(res->x_size(dlayer), res->y_size(dlayer), nullptr);
    for (vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;
        set<node*> proj;
        region2*& region = matrix((int)(nd->x*factor), (int)(nd->y*factor));
        
        if (region != nullptr) continue;
        
        region = new region2();
        res->recurse_from_node(n, name, proj);
        for (set<node*>::iterator piter = proj.begin(); piter != proj.end(); ++piter) {
            node* p = *piter;
            layer1_data* pd = (layer1_data*)p->data;

            if (pd->z == 0) 
                region->add(regions[pd->m], ipoint2(pd->x, pd->y) - centers[pd->m],
                    factor);
        }
    }

}

part_stat_learning::part_stat_learning(const config_dictionary& cfg) :
    stat_learning(0, 0, 0)
{ 
    init(cfg);
}

void part_stat_learning::prepare_for_update(layer1_result* res)
{
    res->set_attr(res->shape_nodes[player].begin(), res->shape_nodes[player].end(), 
        ALLOW_UPDATE_ATTR);
    res->connect_neighbors_circular_2(res->shape_nodes[player], nb_size()/2, nb_size()/2, 0, 
        ALLOW_UPDATE_ATTR, atom("toNeighbor"));
    res->set_attribute(ALLOW_UPDATE_ATTR);
}

string part_stat_learning::name_at(int c, int d)
{
    stringstream str;

    str << ::at(cpaths, c) << " --- " << ::at(dpaths, d) << '\0';
    return string(str.str());
}

int part_stat_learning::update(layer1_result* res)
{
    if (clayer < 0 || clayer > res->max_layer_index() || dlayer < 0 || dlayer > res->max_layer_index() 
        || player > clayer || player > dlayer) 
        return -1;

    if (!res->is_attribute_set(ALLOW_UPDATE_ATTR))
        prepare_for_update(res);

    layer1_result_walker walker(player);
    layer1_result_walker::container_t start;
    layer1_result_walker::result_t cmap, dmap;
    vector<node*>& cnodes = res->shape_nodes[clayer];
    vector<node*>& dnodes = res->shape_nodes[dlayer];
    vector<node*>& pnodes = res->shape_nodes[player];
    int name = atom("toNeighbor");

    // get paths from cnodes
    walker.add_to_container(start, cnodes.begin(), cnodes.end());
    res->recurse2(walker, start, atom("toPrevLayer"));
    walker.get_result(cmap);

    // get paths from dnodes
    start.clear();
    walker.reset();
    walker.add_to_container(start, dnodes.begin(), dnodes.end());
    res->recurse2(walker, start, atom("toPrevLayer"));
    walker.get_result(dmap);

    cerr << "cpaths" << endl;
    cerr << cmap.size() << endl;
    cerr << cmap.begin()->second << endl;

    cerr << "dpaths" << endl;
    cerr << dmap.size() << endl;
    cerr << dmap.begin()->second << endl;

    // for all "center nodes" cn
    for (vector<node*>::iterator iter = pnodes.begin(); iter != pnodes.end(); ++iter) {
        node* cn = *iter;
        layer1_data* cd = (layer1_data*)cn->data;
        ipoint2 cnpos(cd->x, cd->y);
        layer1_result_walker::result_t::iterator cmiter;

        // for all paths ending in cn
        for (cmiter = cmap.find(cn); cmiter != cmap.end() && cmiter->first == cn; ++cmiter) {
            path_t& path = cmiter->second;
            path_map_t::iterator cpmiter = cpmap.find(path);

            if (cpmiter == cpmap.end()) continue;

            int ci = cpmiter->second;

            // for all nodes n in the neighborhood of cn
            foreach_neighbor(cn, name, niter) {
                node* n = neighbor_node(niter);
                layer1_data* nd = (layer1_data*)n->data;
                ipoint2 posdelta(nd->x - cnpos.x, nd->y - cnpos.y);
                layer1_result_walker::result_t::iterator dmiter;
                
                // for all paths ending in n
                for (dmiter = dmap.find(n); dmiter != dmap.end() && dmiter->first == n; ++dmiter) {
                    path_t& npath = dmiter->second;
                    path_map_t::iterator dpmiter = dpmap.find(npath);

                    if (dpmiter == dpmap.end()) continue;

                    int di = dpmiter->second;

                    // update!
                    stat_learning::update(ci, di, posdelta.x, posdelta.y);
                }

            }

        }

    }
    
    return 1;
}

//int part_stat_learning::update(layer1_result* res)
//{
//    if (clayer < 0 || clayer > res->max_layer_index() || dlayer < 0 || dlayer > res->max_layer_index()) 
//        return -1;
//
//    if (!res->is_attribute_set(ALLOW_UPDATE_ATTR))
//        prepare_for_update(res);
//
//    double factor = (double)res->x_size(clayer)/res->x_size(dlayer);
//    vector<node*>& cnodes = res->shape_nodes[clayer];
//    int name = atom("toNeighbor").get_index();
//    matrix<region2*> rmatrix;
//    int nbsize2 = sqr(nb_size());
//
//    make_region_matrix(rmatrix, res, clayer, dlayer);
//    for (vector<node*>::iterator iter = cnodes.begin(); iter != cnodes.end(); ++iter) {
//        node* cn = *iter;
//        layer1_data* cd = (layer1_data*)cn->data;
//
//        if (!cn->is_attr_set(ALLOW_UPDATE_ATTR) || cd->val < center_val_threshold)
//            continue;
//
//        ipoint2 cnpos(cd->x, cd->y);
//        region2* cr = rmatrix(cnpos.x, cnpos.y);
//
//        foreach_neighbor(cn, name, niter) {
//            node* n = neighbor_node(niter);
//            layer1_data* nd = (layer1_data*)neighbor_node_data(niter);
//            ipoint2 npos((int)(factor*nd->x), (int)(factor*nd->y));
//
//            if (cnpos.distance2(npos) > nbsize2) continue;
//
//            if (cr->intersection_size(*rmatrix(npos.x, npos.y)) < 4)
//                stat_learning::update(cd->m, nd->m, npos.x - cnpos.x, npos.y - cnpos.y);
//        }
//        
//    }
//
//    for_each_iter(rmatrix, matrix<region2*>::iterator, iter) {
//        if (*iter != nullptr) delete *iter;
//    }
//
//    return 1;
//}

// stat_learning
//////////////////

stat_learning::stat_learning(int cs, int ds, int nbsize) : 
    nbhoods(cs, ds, nbhood_t(nbsize, nbsize, 0))
{
}

void stat_learning::reset()
{
    for_each_iter(nbhoods, matrix<nbhood_t>::iterator, iter) iter->fill(0);
}

void stat_learning::reset(int cs, int ds, int nbsize)
{
    nbhoods.resize(cs, ds, nbhood_t(nbsize, nbsize, 0));
}

int stat_learning::update(int c, int d, int x, int y)
{
    if (c >= 0 && c < c_size() && d >= 0 && d < d_size()) {
        int nbs = nb_size();
        int nbs2 = nbs/2;
        int nx = nbs2 + x, ny = nbs2 + y;

        if (nx >= 0 && nx < nbs && ny >= 0 && ny < nbs)
            return ++(nbhoods(c, d).at(nx, ny));
    }
    return -1;
}

    //typedef set<layer2_sequence>::iterator set_iterator;

    //part_lib* library;
    //int source_layer_index;                 // index (i.e. layer name - 1) of the layer which 
    //                                        // serves as source
    //int nb_size;                            // = 17
    //int to_neighbor;                        // = atom("toNeighbor").get_index()
    //int to_region;                          // = atom("toRegion").get_index()
    //double region_intersection_threshold;   // = 0.3 (if  int <= thresh the object survives)
    //double region2_intersection_threshold;  // = 0.3
    //int max_nbhood;                         // = 5
    //int max_nbhood_mask;                    // = 5
    //int max_border;                         // = max_nbhood_mask/2
    //double max_threshold_percent;           // = 0.5
    //int max_max;                            // = -1
    //double seq_threshold;                   // = 0.5
    //double edge_threshold;                  // = 0.5
    //int seq_min_intersection;               // = 0
    //int seq_max_intersection;               // = 1000    
    //double center_val_threshold;            // = 0.2
    //double nb_val_threshold_rel;            // = 0.6
    //double nbthresh_min, nbthresh_max;      // = 0.0, 0.0 i.e. use nb_val_threshold_rel
    //int blur_dim;                           // = 5
    //double blur_sigma;                      // = 0.75
    //double sequence_count_threshold;        // = 0.5
    //img* blur_mask;                         // = gaussian_mask(blur_dim, blur_dim, blur_sigma)
    //int min_part_length;                    // = 2
    //int max_part_length;                    // = 3
    //double part_count_factor;               // = 1.0
    //double part_seq_size_factor;            // = 1.0
    //int part_max_number;                    // = 20
    //int part_eq_tolerance;                  // = 0
    //double layer_contraction;               // = 1.0
    //int center_island_radius;               // = -1
    //double region_eq_threshold;             // = -1.0
    //int region_eq_tolerance;                // = 0
    //double island_threshold_percent;        // = max_threshold_percent
    //double min_update_count_percent;        // = 0.1
    //bool individual_max;                    // false
    //double distribution_sigma;              // < 0.0 take the whole distribution from the map (default); 
    //                                        // > 0.0 sigma of the Gaussian mask - dimension of the mask is 
    //                                        //       max_nbhood_mask x max_nbhood_mask
    //int rnb_radius;                         // = 15
    //int rnb_min_radius;                     // = 2

    //string nbhoods_image_template;          // = "" ( == do not save)
    //string nbhoods_matrix_template;         // = "" ( == do not save)
    //string maxima_image_template;           // = "" ( == do not save)

    //vector<matrix<bool>*> regions;          // regions of layer 1 masks (pointers to library original)
    //vector<iipair> region_centers;          // centers of regions
    //vector<matrix<bool>*> sregions;         // regions "of source_layer_index" (pointers to library originals)
    //vector<iipair> sregion_centers;         // centers of sregions
    //matrix<matrix<int> > nbhoods;            // tmp result 1, counts occurrences
    //matrix<vector<nb_data> > nb_maxima;      // tmp result 2, maxima
    //vector<set<int> > nb_rsets;              // tmp result 3, sets of sregions
    //int_map type_map;                       // type_map(t) -> index of the type
    //                                        // type_map.inv(i) -> type associated with index i

//    bool debug;
//    unsigned max_seq_updates;               
//    set<layer2_sequence> sequences;
//
//    s_learning(layer1_result*);
//    s_learning(config_dictionary& cfg);
//    ~s_learning();
//
//    void init(layer1_result*);
//    void cfg_init(config_dictionary& cfg);
//    void destroy();
//    void reset(part_lib* newlib, const vector<int>& initial_parts);
//    
//    //void choose_parts(const vector<int>& tstat);
//    //template<class I> void choose_parts(I begin, I end);
//    void make_region_sets(matrix<set<int> >& ir, layer1_result* l);
//    void prepare_for_update(layer1_result* l);
//    void update_nbhoods(layer1_result* l);
//    template<class I> void update_nbhoods(I begin, I end);
//    void find_maxima(matrix<double>& im, vector<nb_data>& max, int type);
//    void find_maxima();
//    void filter_sequences();
//    void apply_to_result(layer1_result* res);
//    template<class I> void apply_to_result(I begin, I end);
//    void add_to_library(const vector<int>& indices, int maxadd);
//    void add_to_library2(list<pair<int, grid_points> >& gpl, vector<int>& indices, int maxadd);
//
//    //void write_init(const char* fname);
//    //void write_init(ostreamer& os);
//    void write_nbhoods(const string& file);
//    void write_nbhoods(ostreamer& os);
//    //void read_init(istreamer& is);
//    void read_nbhoods(const string& fname);
//    void read_nbhoods(istreamer& is);
//
//    void display_nbhoods(const string& file_template, const string& matrix_template);
//    void display_maxima(const string& file_template);
//
//    void add_edges(layer1_result* res);
//    void add_rnb_edges(layer1_result* res, int dest_layer_index, int iradius, int oradius);
//
//    //matrix< get_matrix(int i, int j);
//    //void update
//
//protected:
//    void get_region_set_from_node(layer1_result* l, node* n, set<int>& result);
//



// hierarchy --- structure for TEST purposes only
///////////////////////////////////////////////////

class hierarchy : graph {
public:
    vector<layer1_result*> images;
    vector<node*> layer_nodes[MAX_LAYER_NUMBER];

    hierarchy(const vector<layer1_result*> pimages);

    void write_mma(ostream& os);
    void write_net(ostream& os);
protected:
    int build_hierarchy(int layer);
};


hierarchy::hierarchy(vector<layer1_result*> pimages) : 
    images(pimages)
{
    int layer = 1;

    while (build_hierarchy(layer) > 0) layer++;
}
    
// build the hierarchy from the specified layer 
// returns the number of added edges
int hierarchy::build_hierarchy(int layer)
{
    int to_prev_layer = atom("toPrevLayer").get_index();
    int to_proj = atom("toProjection").get_index();
    int result = 0;

    for (int i = 0; i < (int)images.size(); ++i) {
        layer1_result* res = images[i];
    
        if (layer >= (int)res->shape_nodes.size()) continue;

        vector<node*>& nodes = res->shape_nodes[layer];
        vector<node*>& lnodes = layer_nodes[layer];

        for (vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
            node* n = *iter;

            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if (nd->m >= (int)lnodes.size()) lnodes.resize(nd->m + 1);
                if (lnodes[nd->m] == nullptr) lnodes[nd->m] = add_node();

                node* ln = lnodes[nd->m];
                set<node*> proj;
                
                res->recurse_from_node(n, to_prev_layer, proj);
                for (set<node*>::iterator niter = proj.begin(); niter != proj.end(); ++niter) {
                    add_edge(ln, *niter, to_proj);
                    ++result;
                }

                n = nd->next;
            } 
        }
        
    }
    return result;
}

void hierarchy::write_mma(ostream& os)
{
    int to_proj = atom("toProjection").get_index();
    map<node*, int> node_enum;
    int ecount = 0;

    os << '{';
    for (int i = 1; i < MAX_LAYER_NUMBER && !layer_nodes[i].empty(); ++i) {
        vector<node*>& nodes = layer_nodes[i];

        os << endl;
        for (vector<node*>::iterator niter = nodes.begin(); niter != nodes.end(); ++niter) {
            node* n = *niter;

            if (n == nullptr) continue;

            pair<map<node*, int>::iterator, bool> insertr;
            int ni;
            
            insertr = node_enum.insert(pair<node*, int>(n, (int)node_enum.size()));
            ni = insertr.first->second;
            foreach_neighbor(n, to_proj, nbiter) {
                insertr = node_enum.insert(pair<node*, int>(neighbor_node(nbiter), (int)node_enum.size()));
                if (ecount++ > 0) os << ',';
                os << '{' << ni << ',' << insertr.first->second << '}';
            }
        }
    }
    os << '}' << endl;
}

void hierarchy::write_net(ostream& os)
{
    int to_proj = atom("toProjection").get_index();
    map<node*, int> node_enum;
    list<iipair> edges;
    int ecount = 0;

    for (int i = 1; i < MAX_LAYER_NUMBER && !layer_nodes[i].empty(); ++i) {
        vector<node*>& nodes = layer_nodes[i];

        for (vector<node*>::iterator niter = nodes.begin(); niter != nodes.end(); ++niter) {
            node* n = *niter;

            if (n == nullptr) continue;

            pair<map<node*, int>::iterator, bool> insertr;
            int ni;
            
            insertr = node_enum.insert(pair<node*, int>(n, (int)node_enum.size()));
            ni = insertr.first->second;
            foreach_neighbor(n, to_proj, nbiter) {
                insertr = node_enum.insert(pair<node*, int>(neighbor_node(nbiter), (int)node_enum.size()));
                edges.push_back(iipair(ni, insertr.first->second));
            }
        }
    }

    // Write vertices
    int vnum = (int)node_enum.size();

    os << "*Vertices " << vnum << endl;
    for (int i = 1; i <= vnum; ++i) {
        os << i << ' ' << '\"' << i << '\"' << endl;
    }

    // Write edges
    os << "*Edges " << edges.size() << endl;
    for (list<iipair>::iterator i = edges.begin(); i != edges.end(); ++i) {
        iipair& e = *i;
        os << e.first + 1 << ' ' << e.second + 1 << ' ' << '1' << endl;
    }

}

void export_hierarchy(const char* fnames, const char* type)
{
    list<string> files;
    layer1_result* res;
    vector<layer1_result*> results;

    list_directory(files, fnames);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        read_layer1_result(res, *file);
        results.push_back(res);
        cerr << "Loading " << *file << endl;
    }

    hierarchy h(results);
    ofstream os;

    cerr << "Writing hierarchy in ";

    if (type[0] == 'm' || type[0] == 'M') {
        os.open("hierarchy.m");
        cerr << "Mathematica ...";
        h.write_mma(os);
    } else {
        os.open("hierarchy.net");
        cerr << "Pajek ...";
        h.write_net(os);
    }

    cerr << endl;
    os.close();

    for (vector<layer1_result*>::iterator i = results.begin(); i != results.end(); ++i)
        delete *i;
}

// optimization_module ----> to optimization.h/cpp
////////////////////////

// optimization_data --> optimization_data 

class optimization_module_base {
protected:
    int steps;
    optimization_data* data;
public:
    optimization_module_base(int nsteps, optimization_data* pdata);
    virtual int execute() = 0;
};

optimization_module_base::optimization_module_base(int nsteps, optimization_data* pdata)
{
    data = pdata;
    steps = nsteps;
}




// implementations
///////////////////////////////////////////////////////////////////////////////



// layer1_optimization usage test
///////////////////////////////////

//void optimize_layer(const config_dictionary& cfg, const string& pattern)
//{
//    // init "optimizer" structure
//    part_lib* library;
//    int layer;
//
//    read_library(cfg.get_value_string("library", ""), library);
//    layer = cfg.get_value_int("layer", 1);
//    
//    layer1_optimization optimizer(cfg, library, layer);
//
//    // load files
//    string dir = cfg.get_value_string("src_dir", "");
//    list<string> files;
//    layer1_result* res;
//    list<layer1_result*> res_list;
//
//    end_dir(dir);
//    list_directory(files, dir + pattern);
//    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
//        cerr << "Loading " << *file << endl;
//
//        read_layer1_result(res, dir + *file);
//        if (res) {
//            optimizer.add_to_test_set(res);
//            res_list.push_back(res);
//        }
//    }
//    // perform optimization
//    /*
//    int step = 0;
//    char fname[1000];
//    int answer;
//
//    while (true) {
//        cout << "0 = stop, x > 0 = make x steps: "; cin > > answer;
//
//        if (answer == 0) break;
//        while (answer-- > 0) {
//            sprintf(fname, "c:\\work\\step%d.bmp", step++);
//            optimizer.make_step();
//            optimizer.print_parts();
//            optimizer.save_parts(fname);
//        }
//    } 
//    */
//    int step = 0;
//    int libno = 0;
//    char fname[1000];
//    int answer;
//
//    while (true) {
//        cout << "-1 = save new library, 0 = stop, x > 0 = make x steps: "; cin > > answer;
//
//
//        if (answer == 0) break;
//        if (answer == -1) {
//            vector<int> parts;
//            part_lib* new_library;
//
//            read_library(cfg.get_value_string("library", ""), new_library);
//            optimizer.get_parts(parts);
//            new_library->keep_parts(layer, parts);
//            sprintf(fname, "library%d.plb", libno);
//            new_library->save(fname);
//            delete new_library;
//        } else {
//            while (answer-- > 0) optimizer.make_step_1();
//
//            sprintf(fname, "step_random%d.bmp", step++);
//            optimizer.print_parts();
//            optimizer.save_parts(fname);
//        }
//    } 
//
//    // delete layer1_results
//    for (list<layer1_result*>::iterator iter = res_list.begin(); iter != res_list.end(); ++iter) {
//        delete *iter;
//    }
//    delete library;
//}


/*
void optimize_layer2(const config_dictionary& cfg, const string& pattern)
{
    typedef map<string, int> t_map;

    cr_optimization optimizer;

    string dir = cfg.get_value_string("src_dir", "");
    int rndchoice = cfg.get_value_int("random_select", -1);

    list<string> files;
    layer1_result* res;
    //list<layer1_result*> res_list;
    t_map typemap;

    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, pattern, dir))
        list_directory(files, dir + pattern);
    if (rndchoice > 0) 
        random_selection(files, rndchoice, cfg.get_value_string("dropped_files", ""));
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cerr << "Loading " << *file << endl;
        read_layer1_result(res, dir + *file);
        if (res) {
            string pref = get_prefix(*file, "_");
            pair<t_map::iterator, bool> ipair = typemap.insert(t_map::value_type(pref, (int)typemap.size()));

            res->type_index = ipair.first->second;
            optimizer.add_to_test_set(layer1_result_ptr(res));
            //res_list.push_back(res);
        }
    }
    optimizer.inter_layer_optimization(1, 10);
}
*/



/*void optimization_test1(const config_dictionary& cfg, optimization_data& oopt)
{
    vector<cr_optimization_base*> lopts;
    int start_layer = 1, end_layer = 3;

    for (int l = start_layer; l <= end_layer; ++l) {
        lopts.push_back(new cr_layer_optimization(&oopt, 100, l));
    }

    cr_set_optimization cropt(&oopt, 2, lopts);

    cropt.set_library(oopt.get_library());
    cropt.add_range_to_test_set(oopt.get_test_set_begin(), oopt.get_test_set_end());

    cropt.execute();   
    
    part_lib* lresult = cropt.get_best_library();
    list<layer1_result_ptr> sresult;
    string name;
    int i = 0;
    string lib_export_name = cfg.get_value_string("lib_export_name", "lib");
    string res_export_name = cfg.get_value_string("res_export_name", "res");

    cropt.get_best_test_set(sresult);

    PRINT_INFO("");
    PRINT_INFO("SAVING");
    name = lib_export_name + (end_layer + 1) + ".plb";
    lresult->save(name);
    for (list<layer1_result_ptr>::iterator iter = sresult.begin(); iter != sresult.end(); ++iter) {
        name = res_export_name + (i++) + string(".ly") + (end_layer + 1);
        PRINT_INFO(name);

        (*iter)->delete_edges_complement(atom("toPrevLayer").get_index());
        (*iter)->save(name);
    }

}*/


void test_overall_optimization(const config_dictionary& cfg, const string& pattern)
{
    string dir = cfg.get_value_string("src_dir", "");
    list<string> files;
    layer1_result* res;
    optimization_data oopt(cfg);

    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, pattern, dir))
        list_directory(files, dir + pattern);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cerr << "Processing " << *file << endl;
        read_layer1_result(res, dir + *file);
        if (res) 
            oopt.add_to_test_set(streamed_pointer(res));
    }
//    optimization_test1(cfg, oopt);
}

void save_confidence_vector(vector<pair<double, bool> >& confidence, const string& fname)
{
    typedef pair<double, bool> confidence_t;

    ofstream os(fname.c_str());

    if (!os.is_open()) return;

    vector<confidence_t>::iterator iter;
    bool first;

    sort(confidence.begin(), confidence.end(), greater<confidence_t>());

    // Write fp vector.
    first = true;
    os << "fp={"; 
    for (iter = confidence.begin(); iter != confidence.end(); ++iter) {
        if (first) first = false; else os << ',';
        if (iter->second) os << "0.0"; else os << iter->first;
    }
    os << '}' << ';' << endl;

    // Write tp vector.
    first = true;
    os << "tp={"; 
    for (iter = confidence.begin(); iter != confidence.end(); ++iter) {
        if (first) first = false; else os << ',';
        if (iter->second) os << iter->first; else os << "0.0"; 
    }
    os << '}' << ';' << endl;

    os.close();
}

void test_multicheck(const config_dictionary& cfg, const char* fname)
{
    typedef pair<double, bool> confidence_t;

    string dir = cfg.get_value_string("src_dir", "");
    list<string> patterns;
    layer1_result* res;
    part_lib* library;
    list<irectangle2> gtruths;
    string extension;
    string libname;
    int misses;
    int layer;
    double thresh;
    int total_misses = 0;
    int total_objects = 0;
    int objects_detected = 0;
    int object_layer, object_index;
    set<int> parts;
    vector<confidence_t> confidence;
 
    cfg.get_value(thresh, "threshold", true);
    cfg.get_value(object_layer, "object_layer", true);
    cfg.get_value(object_index, "object_index", true);
    cfg.get_value(extension, "extension", true);
    cfg.get_value(libname, "library", true);
    layer = object_layer - 1;
    read_library(libname, library);
    library->get_parent_parts(parts, object_layer, object_index);
    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(patterns, fname, dir))
        patterns.push_back(fname);
    for (list<string>::iterator pattern = patterns.begin(); pattern != patterns.end(); ++pattern) {
        list<string> files;
        vector<bool> hits;
        int maxmisses = 0;
        
        list_directory(files, dir + *pattern + string("_*") + extension);
        for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
            read_layer1_result(res, dir + *file);
            read_groundtruth(gtruths, dir + *file);
            if (res != nullptr) {
                cerr << '.';
                if (!gtruths.empty()) {
                    res->count_hits_inhibited(hits, misses, confidence, gtruths, layer, parts, thresh);
                    if (misses > maxmisses) maxmisses = misses;
                }
                delete res;
            }
        }
        cout << *pattern << "  hits: " << hits << "  misses: " << maxmisses << endl;
        total_objects += (int)hits.size();
        objects_detected += (int)count(hits.begin(), hits.end(), true);
        total_misses += maxmisses;
    }
    cout << "Objects detected: " << objects_detected << "  Total objects: " << total_objects << endl;
    cout << "Total misses: " << total_misses << endl;
    
    save_confidence_vector(confidence, cfg.get_value_string("out_file", "tp-fp.m"));
}

void test_check(const config_dictionary& cfg, const char* fname)
{
    typedef pair<double, bool> confidence_t;

    string dir = cfg.get_value_string("src_dir", "");
    bool save = cfg.get_value_bool("save_boxes", false);
    bool saveimg = cfg.get_value_bool("save_image", false);
    list<string> patterns;
    string filepattern;
    int misses;
    int layer;
    double thresh;
    int total_misses = 0;
    int total_objects = 0;
    int objects_detected = 0;
    vector<confidence_t> confidence;
    vector<int> vec;
    set<int> parts;
	map<string, int> catmap;
	string libname;
	part_lib* library;
    response_filter rfilter;
    
    cfg.get_value(vec, "parts");
    cfg.get_value(thresh, "threshold", true);
    cfg.get_value(layer, "layer", true);
    cfg.get_value(filepattern, "pattern", true);
	cfg.get_value(libname, "library", true);

    rfilter.add_response(R_RESPONSE, cfg.get_value_double("r_response_threshold", 0.0));
    rfilter.add_response(G_RESPONSE, cfg.get_value_double("g_response_threshold", 0.0));
    rfilter.add_response(RR_RESPONSE, cfg.get_value_double("rr_response_threshold", 0.0));

	read_library(libname, library);
	if (library == nullptr) {
		cerr << "Can not open library (library is needed to extract category names)" << endl;
		return;
	} else {
		if (layer < 0 || layer > library->max_layer_index()) {
			cerr << "Layer #" << layer << " does not exist." << endl;
			return;
		}
		for (vector<node*>::iterator iter = library->parts[layer].begin(); iter != library->parts[layer].end(); ++iter) {
			cpart_data* pd = dynamic_cast<cpart_data*>((*iter)->data);

			if (pd) 
				catmap.insert(pair<string, int>(pd->name, pd->type));
		}
        delete library;
	}
	if (catmap.empty()) {
		cerr << "No category names can be extracted from layer " << layer << " in the library " << libname << "!" << endl;
		return;
	}
	
    parts.insert(vec.begin(), vec.end());
    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(patterns, fname, dir))
        patterns.push_back(fname);
    for (list<string>::iterator pattern = patterns.begin(); pattern != patterns.end(); ++pattern) {
        list<string> files;

        list<pair<irectangle2, string> > gtruths;
        vector<const layer1_result::box_data_t*> hitp;
        vector<bool> hits;
        int misses = 0;
        list<layer1_result::box_data_t> boxes;
        int x_size0 = 0;
        double factor;
        img* im = nullptr;
        int maxwidth = 0;
        int origwidth = 0;
        int border = 0;
        char buf[1024];

        
        sprintf(buf, filepattern.c_str(), pattern->c_str());
        list_directory(files, dir + buf);
        for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
		    layer1_result* res;

			read_layer1_result(res, dir + *file);
            if (res != nullptr) {
                if (x_size0 > 0) factor = (double)x_size0/res->x_size(0);
                else {
                    read_groundtruth(gtruths, dir + *file);
                    factor = 1.0;
                    x_size0 = res->x_size(0);
                }
                res->get_boxes(boxes, layer, RR_RESPONSE, rfilter, factor);

                if (res->x_size(0) > maxwidth) {
                    maxwidth = res->x_size(0);
                    origwidth = res->original_width;
                    border = res->border;
                }
                if (saveimg && (im == nullptr || im->width < res->x_size(0))) {
                    if (im != nullptr) delete im;
                    im = res->get_image_reconstructed(0, 0, vector<int>(), true);
                }
                cout << '.';
                delete res;
            }
        }

        if (!gtruths.empty()) {
			list<pair<irectangle2, int> > igtruths;

			for (list<pair<irectangle2, string> >::iterator giter = gtruths.begin(); giter != gtruths.end(); ++giter) {
				map<string, int>::iterator miter = catmap.find(giter->second);
				if (miter != catmap.end()) 
					igtruths.push_back(pair<irectangle2, int>(giter->first, miter->second));
			}
            inhibit_boxes(boxes, thresh);
            hitp = check_hits(hits, misses, igtruths, boxes, thresh);
        }

        if (save) {
            ofstream os((dir + *pattern + "_boxes.txt").c_str());

            for (list<layer1_result::box_data_t>::iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
                irectangle2 box = biter->box;

                if (origwidth > 0) {
                    box -= ipoint2(border, border);
                    box.resize(origwidth/maxwidth);
                }
                os << biter->box;
                os << ' ' << ((find(hitp.begin(), hitp.end(), &(*biter)) == hitp.end()) ? 0 : 1);
                os << endl;
            }
            os.close();
        }
        if (im != nullptr) {
            for (list<layer1_result::box_data_t>::iterator biter = boxes.begin(); biter != boxes.end(); ++biter)
                if (find(hitp.begin(), hitp.end(), &(*biter)) == hitp.end())
                    im->draw_box(biter->box, COL_TO_REAL(255, 0, 0), 3);
                else 
                    im->draw_box(biter->box, COL_TO_REAL(0, 255, 0), 3);

            im->save(dir + *pattern + ".png");
            delete im;
        }
        cout << *pattern << "  hits: " << hits << "  misses: " << misses << endl;
        total_objects += (int)hits.size();
        objects_detected += (int)count(hits.begin(), hits.end(), true);
        total_misses += misses;
    }
    cout << "Objects detected: " << objects_detected << "  Total objects: " << total_objects << endl;
    cout << "Total misses: " << total_misses << endl;
    
    //save_confidence_vector(confidence, cfg.get_value_string("out_file", "tp-fp.m"));
}

void part_hits_statistics(const config_dictionary& cfg, const char* fname)
{
    typedef map<int, pair<int, int> > stat_t;

    string dir = cfg.get_value_string("src_dir", "");
    string gtruthext = cfg.get_value_string("groundtruth_extension", ".groundtruth");
    list<string> files;
    layer1_result* res;
    list<irectangle2> gtruths;
    stat_t stat;
    int layer;
    double thresh;
 
    cfg.get_value(thresh, "threshold", true);
    cfg.get_value(layer, "layer", true);
    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, dir))
        list_directory(files, dir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        read_layer1_result(res, dir + *file);
        if (res != nullptr) {
            cerr << "Processing " << *file << " ...";
            read_groundtruth(gtruths, change_extension(dir + *file, gtruthext));
            if (!gtruths.empty()) {
                list<layer1_result::box_data_t> boxes;

                res->get_boxes(boxes, layer, set<int>(), 0.0);
                res->count_hits(stat, gtruths, boxes, thresh);
            }
            delete res;
            cerr << endl;
        }
    }
    ofstream f(cfg.get_value_string("out_file", "hit-statistics.txt").c_str());

    for (stat_t::iterator iter = stat.begin(); iter != stat.end(); ++iter) {
        f << iter->first << ' ' << iter->second.first << ' ' << iter->second.second << endl;
    }
    f.close();
}

void merge_scales(const config_dictionary& cfg, const char* fname)
{
    string dir = cfg.get_value_string("src_dir", "");
    list<string> patterns;
    string extension;
    string suffix = cfg.get_value_string("suffix", "_merged");
 
    cfg.get_value(extension, "extension", true);
    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(patterns, fname, dir))
        patterns.push_back(fname);
	int border = 0;
	cfg.get_value(border, "border_size");
    for (list<string>::iterator pattern = patterns.begin(); pattern != patterns.end(); ++pattern) {
        list<string> files;
        layer1_result* result = nullptr;
        
        list_directory(files, dir + *pattern + string("_?") + extension);
        for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
            layer1_result* res;

            read_layer1_result(res, dir + *file);
            if (res != nullptr) {
                if (result == nullptr) {
                    result = res;
                    cerr << "Adding to " << *file << endl;
                } else { 
                    cerr << "    " << *file;
                    result->merge(res, border);
                    cerr << endl;
                    delete res;
                }
            }
        }
        if (result != nullptr) {
            save_layer1_result(result, dir + *pattern + suffix + extension);
            delete result;
        }
    }
}



void test_3(const config_dictionary& cfg, const char* fname)
{
/*    layer1_result* res;

    read_layer1_result(res, fname);
    if (res != nullptr) {
        res->add_reverse_edges(atom("toPrevLayer"), atom("toPrevLayerRev"));
        res->save("c:\\temp\\rev.lyx");
        delete res;
    }
*/   
    part_lib* library;

    read_library(fname, library);
    if (library != nullptr) {
        int slayer = cfg.get_value_int("s_layer", 2);
        int player = cfg.get_value_int("p_layer", 1);
        int part = cfg.get_value_int("part", 0);
        library_path_walker walker(player);
        library_path_walker::container_t start;
        library_path_walker::result_t result;

        for (vector<node*>::iterator piter = library->parts[slayer].begin();
                piter != library->parts[slayer].end(); ++piter) {
            path_walker::add_to_container(start, *piter);    
        }

        library->recurse2(walker, start, atom("lyrPrev"));
        walker.get_result(result);

        for (library_path_walker::result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
            cerr << *iter << endl;;
        }

        delete library;
    }
}


void slearning(const config_dictionary& cfg, const char* fname)
{
    //stat_learning slearner(10, 10, 21);
    //
    //for (int i = 0; i < 100000; ++i) {
    //    int c = random_int(0, 10), d = random_int(0, 10);
    //    int x = random_int(-10, 11), y = random_int(-10, 11);
    //    slearner.update(c, d, x, y);
    //}

    int csize = cfg.get_value_int("c_size", 10);
    int dsize = cfg.get_value_int("d_size", 10);
    int nbsize = cfg.get_value_int("nb_size", 23);
    int clayer = cfg.get_value_int("c_layer", 1);
    int dlayer = cfg.get_value_int("d_layer", 0);
    int player = cfg.get_value_int("p_layer", 0);
    string dir = cfg.get_value_string("src_dir", "");

    //part_stat_learning slearner(clayer, dlayer, player, nbsize);
    part_stat_learning slearner(cfg);

    layer1_result* res;

    end_dir(dir);
    read_layer1_result(res, dir + fname);

    if (res == nullptr) return;

    clock_t start = clock();
    slearner.update(res);
    clock_t end = clock();

    cerr << "Update time " << ((double)(end - start))/CLOCKS_PER_SEC << " sec." << endl;

    for (int x = 0; x < min(10, slearner.c_size()); ++x) {
        for (int y = 0; y < min(10, slearner.d_size()); ++y) {
            stat_learning::nbhood_t& nb = slearner.at(x, y);
            img im(nb.width, nb.height, 0.0);

            for_each_xy_int(nb, i, j) {
                im(i, j) = nb(i, j);
            }
            cerr << '[' << x << ',' << y << ']' << ' ';
            try { 
                cerr << slearner.name_at(x, y);
            } catch (...) 
            {
            }
            cerr << endl;
            string name = string("c:\\temp\\test") + x + "_" + y + ".png";
            im.save_jet_colormap(name.c_str(), -1000);
        }
    }

    delete res;
}


bool pred(node* n)
{
    return ((layer1_data*)n->data)->z <= 1;
}


class mem_test_class {
public:
    int x;
    double y;
    size_t z;

    mem_test_class(int px, double py) : x(px), y(py) { }


};

class mem_test_class_2 : public mem_test_class, public allocator_base {
public:
    double z;

    mem_test_class_2(double pz) : mem_test_class(12, 45.6), z(pz) { }
};

void test1_old()
{
#if defined WIN32 | defined WIN64
    stdext::hash_multimap<int, double> aaa;


    mem_test_class* cl1 = new mem_test_class(1, 100.0);
    mem_test_class_2* cl2 = new mem_test_class_2(11.0);
    //mem_test_class_2 cl2;

    cerr << "mem_test_class " << "  " << sizeof(*cl1) << endl;
    cerr << "mem_test_class_2 " << "  " << sizeof(*cl2) << endl;

    mem_allocator alloc(sizeof(double), 50, "test");
    double* arr[100];

    for (int i = 0; i < 96; ++i) {
        arr[i] = (double*)alloc.alloc_mem();
        *(arr[i]) = i + 0.1;
    }
    alloc.free_mem(arr[7]);

    delete cl1;
#endif
}   


struct rectangle_hierarchy {
protected:
    struct rh_layer {
        typedef pair<irectangle2, irectangle2*> list_item_t;
        list<list_item_t> rlist;

        rh_layer() : rlist() { }
    };

    vector<rh_layer> layers;
public:
    rectangle_hierarchy(int leafdim, int lcount);
};

rectangle_hierarchy::rectangle_hierarchy(int leafdim, int lcount) : 
    layers()
{
        
}


template<class I> ipoint2 tpoint_center(I begin, I end)
{
    ipoint2 result;
    int count = 0;

    while (begin != end) {
        result += ipoint2(begin->first, begin->second);
        ++count;
        ++begin;
    }
    return result/count;
}

template<class T, class U, class V, class Pred> void partition2(vector<vector<triple<T, U, V> > >& prt, vector<triple<T, U, V> > srcpts, Pred f)
{
    if (srcpts.size() == 1) prt.push_back(srcpts);
    else {
        sort(srcpts.begin(), srcpts.end(), f);

        size_t n = srcpts.size();
        size_t n2 = n/2;
        typename vector<triple<T, U, V> >::const_iterator iter = srcpts.begin();
        size_t i = 0;

        prt.push_back(vector<triple<T, U, V> >());
        for (; i < n2; ++i, ++iter) prt.back().push_back(*iter);
        prt.push_back(vector<triple<T, U, V> >());
        for (; i < n; ++i, ++iter) prt.back().push_back(*iter);
    }
}

template<class T, class U, class V> void partition_tpoints(vector<vector<triple<T, U, V> > >& prt, 
    const vector<triple<T, U, V> >& srcpts)
{
    prt.clear();

    if (srcpts.size() <= 4) 
        prt.push_back(srcpts);
    else if (srcpts.size() <= 8)
        partition2(prt, srcpts, triple<T, U, V>::first_less);
    else {
        vector<vector<triple<T, U, V> > > tmp;

        partition2(tmp, srcpts, triple<T, U, V>::first_less);
        partition2(prt, tmp[0], triple<T, U, V>::second_less);
        partition2(prt, tmp[1], triple<T, U, V>::second_less);
    }
}

void learn_layer(part_lib* library, list<layer1_result*>& rlist, int layer)
{
    typedef triple<int, int, int> vtriple;

	set<vector<vtriple> > rparts;
    int srclayer = layer - 1;

    for (list<layer1_result*>::iterator iter = rlist.begin(); iter != rlist.end(); ++iter) {
        layer1_result* res = *iter;

        if (res->max_layer_index() < srclayer) continue;

        if (!res->grid(layer)) res->init_grid(srclayer);
        vector<node*>& s_nodes = res->shape_nodes[srclayer];
        vector<vtriple> pts;

        // make triples (x, y, type)
        for (vector<node*>::iterator niter = s_nodes.begin(); niter != s_nodes.end(); ++niter) {
            layer1_data* nd = (layer1_data*)(*niter)->data;
                
            pts.push_back(vtriple(nd->x, nd->y, nd->m));
        }

        if (layer >= 2) {
            sort(pts.begin(), pts.end());
            rparts.insert(pts);
        } else { // layer == 1
            vector<vector<vtriple> > prt;

            partition_tpoints(prt, pts);
            for (vector<vector<vtriple> >::iterator piter = prt.begin(); piter != prt.end(); ++piter) {
                sort(piter->begin(), piter->end());
                rparts.insert(*piter);
            }
        }

    }

    // vvv   D E B U G   vvv
    //for (set<vector<itriple> >::iterator siter = rparts.begin(); siter != rparts.end(); ++siter) {
    //    vector<itriple>& v = *siter;

    //    for (vector<itriple>::iterator viter = v.begin(); viter != v.end(); ++viter) {
    //        cerr << '(' << viter->first << ',' << viter->second << ';' << viter->third << ')' << ' ';
    //    }
    //    cerr << endl;
    //}
    // ^^^   D E B U G   ^^^

    // Add parts to the library
    while (library->layer_count <= layer) ++library->layer_count;

    vector<node*>& parts = library->parts[layer];
    vector<node*>& srcparts = library->parts[layer - 1];

    for (set<vector<vtriple> >::const_iterator siter = rparts.begin(); siter != rparts.end(); ++siter) {
        const vector<vtriple>& v = *siter;
        ipoint2 ccoo(v.front().first, v.front().second);
        ipoint2 center = tpoint_center(v.begin(), v.end());
        node* p = library->add_node(new rpart_data(layer, (int)parts.size(), center.x - ccoo.x, center.y - ccoo.y, 0));

        parts.push_back(p);

        // The first one is center ...
        library->add_edge(p, srcparts[v.front().third], atom("lyrCenter"), atom("lyrCenterBack")); 
        library->add_edge(p, srcparts[v.front().third], atom("lyrPrev"));
        for (vector<itriple>::const_iterator viter = v.begin() + 1; viter != v.end(); ++viter) {
            library->add_edge_2(p, srcparts[viter->third], 
                new part_data_2c(viter->first - ccoo.x, viter->second - ccoo.y), atom("lyrSrc"));
            library->add_edge(p, srcparts[viter->third], atom("lyrPrev"));
        }
        
    }

}

// Add layer to layer1_result for learning purposes.
void add_layer(layer1_result* res, part_lib* library, int layer)
{
    if (library->max_layer_index() < layer) return;

    while ((int)res->shape_nodes.size() <= layer) res->shape_nodes.push_back(vector<node*>());
    while ((int)res->shape_nodes_inhib.size() <= layer) res->shape_nodes_inhib.push_back(vector<node*>());

    int srclayer = layer - 1;
    vector<node*>& s_nodes = res->shape_nodes[srclayer];
    vector<node*>& srcparts = library->parts[srclayer];
    int to_cback = atom("lyrCenterBack");
    int to_center = atom("lyrCenter");
    int to_src = atom("lyrSrc");
    int to_prev = atom("toPrevLayer");

    res->new_grid(res->x_size(srclayer), res->y_size(srclayer), layer);
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;
        node* p = srcparts[nd->m];

        foreach_neighbor(p, to_cback, citer) {
            node* c = neighbor_node(citer);
            rpart_data* cd = (rpart_data*)c->data;
            set<node*> conn;

            int x0 = nd->x + cd->cx; 
            int y0 = nd->y + cd->cy;

            //if (x0 % 2 != 0 || y0 % 2 != 0) continue; // ***

            // check all "subparts"
            bool ok = true;
            int count = 0;

            conn.insert(n);
            foreach_neighbor(c, to_src, siter) {
                node* s = neighbor_node(siter);
                rpart_data* sd = (rpart_data*)neighbor_node_data(siter);
                part_data_2c* ed = (part_data_2c*)neighbor_edge_data(siter);
                node* m = res->node_at_safe(nd->x + ed->x, nd->y + ed->y, srclayer);
                layer1_data* md;

                if (!m || (md = (layer1_data*)m->data)->m != sd->type) { ok = false; break; }
                conn.insert(m);
                ++count;
            }
            if (ok) {
                node* newn = res->add_grid_node(new layer1_data(0.25 + 0.25 * count, cd->type), 
                    x0, y0, layer);

                for (set<node*>::iterator iter = conn.begin(); iter != conn.end(); ++iter) 
                    res->add_edge(newn, *iter, to_prev);
            } 

        }
        
    }
    
    for (list<node*>::iterator iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;

        if (nd->z == layer && n->is_attr_set(IMG_NODE_ATTR)) {
            res->shape_nodes[layer].push_back(n);
            res->shape_nodes_inhib[layer].push_back(n);
        }
    }

}

template<class T> void add_layer(T begin, T end, part_lib* library, int layer)
{
    for (; begin != end; ++begin) add_layer(*begin, library, layer);
}


// member of a rlayer, i.e. rpart is vector<ipoint2>
// and is a *sorted* vector of active "subrectangles"
void find_rparts(vector<vector<ipoint2> >& result, const list<vector<ipoint2> >& vlist)
{
	set<vector<ipoint2> > rparts;

	for (list<vector<ipoint2> >::const_iterator iter = vlist.begin(); iter != vlist.end(); ++iter) {
		const vector<ipoint2>& v = *iter;
		irectangle2 box;

		for (vector<ipoint2>::const_iterator viter = v.begin(); viter != v.end(); ++viter) 
			box.eat(*viter);

		ipoint2 minp((box.ll.x % 2 == 0) ? box.ll.x : box.ll.x - 1, 
			(box.ll.y % 2 == 0) ? box.ll.y : box.ll.y - 1);
		imatrix m(box.x_dim() + 2, box.y_dim() + 2, 0);

		for (vector<ipoint2>::const_iterator viter = v.begin(); viter != v.end(); ++viter) {
			ipoint2 p = *viter - minp;

			m(p.x, p.y) = 1;
		}

		for (int i = 0; i < (int)m.width - 1; i += 2) {
			for (int j = 0; j < (int)m.height - 1; j += 2) {
				vector<ipoint2> rpart;
				
				for (int k = 0; k < 2; ++k) 
					for (int l = 0; l < 2; ++l) 
						if (m(i + k, j + l) > 0) rpart.push_back(ipoint2(k, l));
				rparts.insert(rpart);
			}
		}

	}
	result.assign(rparts.begin(), rparts.end());

}



// 'src' must be sorted to allow for binary_search!
// 'rpart' is always sorted!
void apply_rpart(vector<ipoint2>& result, const vector<ipoint2>& rpart, const vector<ipoint2>& src)
{
	
}

void apply_rparts(vector<ipoint2>& result, const vector<vector<ipoint2> >& rparts, const vector<ipoint2>& src0)
{
	vector<ipoint2> src(src0.begin(), src0.end());

	sort(src.begin(), src.end());
	//for (vector<ipoint2>::iterator
}

void save_rlist(const char* fname, list<layer1_result*>& rlist, int layer)
{
    ofstream os(fname);

    os << '{';
    for (list<layer1_result*>::iterator iter = rlist.begin(); iter != rlist.end(); ++iter) {
        layer1_result* r = *iter;

        if (iter != rlist.begin()) os << ',' << endl;
        os << '{';
        for (vector<node*>::iterator niter = r->shape_nodes[layer].begin(); niter != r->shape_nodes[layer].end(); ++niter) {
            layer1_data* nd = (layer1_data*)(*niter)->data;
            if (niter != r->shape_nodes[layer].begin()) os << ',';
            os << '{' << nd->x << ',' << nd->y << ',' << nd->m << '}';
        }
        os << '}';
    }
    os << '}' << endl;

    os.close();
}

bool connected_2c(node* n1, node* n2, part_data_2c* ed, int name)
{
    foreach_neighbor(n1, name, niter) {
        if (neighbor_node(niter) == n2 && 
            *((part_data_2c*)neighbor_edge_data(niter)) == *ed) return true;
    }
    return false;
}

void combine_libraries(part_lib* newlib, part_lib* rlib, part_lib* library, vector<int>& partition, 
    vector<int>& o2n, vector<ipoint2>& clist)
{
    int srclayer = newlib->max_layer_index();
    int layer = srclayer + 1;
    int objlayer = srclayer + 1;
    int tosrc = atom("lyrSrc");
    int tocenter = atom("lyrCenter");
    int tocenter1 = atom("lyrCenter1");
    int tocenterback = atom("lyrCenterBack");
    int toprev = atom("lyrPrev");
    int tosrcM = atom("lyrSrcM");
    
    // add layer srclayer + 1 -- which is special & connect it to all 
    ++newlib->layer_count;

    vector<node*>* srcrparts = &rlib->parts[0];
    vector<node*>* srcparts = &newlib->parts[srclayer];
    vector<node*>* destparts = &newlib->parts[layer];

    newlib->layer_data[layer] = new streamable_int(0);
    for (vector<node*>::iterator siter = srcrparts->begin(); siter != srcrparts->end(); ++siter) {
        node* rp = *siter;
        rpart_data* rpd = (rpart_data*)rp->data;
        node* newn = newlib->add_node(new rpart_data(layer, (int)destparts->size() /* == rpd->type */, 0, 0, 0), R_PART_ATTR); 
        
        destparts->push_back(newn);
        for (vector<node*>::iterator iter = srcparts->begin(); iter != srcparts->end(); ++iter) {
            node* p = *iter;
            lib_data* pd = (lib_data*)p->data;
            
            if (partition[pd->type] == rpd->type) {
                newlib->add_edge(newn, p, tocenter, tocenterback);
                newlib->add_edge(newn, p, toprev);
            }
        }
    }

    // Add subsequent layers
    for (int srcrlayer = 1; srcrlayer < rlib->layer_count; ++srcrlayer) {
        ++newlib->layer_count;
        srcrparts = &rlib->parts[srcrlayer];
        srcparts = &newlib->parts[layer];
        destparts = &newlib->parts[++layer];

        newlib->layer_data[layer] = new streamable_int(srcrlayer);
        for (vector<node*>::iterator siter = srcrparts->begin(); siter != srcrparts->end(); ++siter) {
            node* rp = *siter;
            rpart_data* rpd = (rpart_data*)rp->data;
            node* newn = newlib->add_node(new rpart_data(layer, (int)destparts->size(), rpd->cx, rpd->cy, 0), R_PART_ATTR);
            int ctype = ((rpart_data*)rp->get_neighbor(tocenter)->data)->type;

            destparts->push_back(newn);
            newlib->add_edge(newn, (*srcparts)[ctype], tocenter, tocenterback);
            newlib->add_edge(newn, (*srcparts)[ctype], toprev);
            foreach_neighbor(rp, tosrc, niter) {
                int ctype = ((rpart_data*)neighbor_node_data(niter))->type;

                newlib->add_edge_2(newn, (*srcparts)[ctype], 
                    new part_data_2c(*((part_data_2c*)neighbor_edge_data(niter))), tosrc);
                newlib->add_edge(newn, (*srcparts)[ctype], toprev);
            }
        }
    }

    // And finally, we have to add the top
    ++newlib->layer_count;
    srcrparts = &library->parts[objlayer];
    srcparts = &newlib->parts[layer];
    destparts = &newlib->parts[++layer];
    vector<node*>* origsrcparts = &newlib->parts[srclayer];

    // Copy part data and add edges: lyrCenterBack form layer - 1 to layer 
    //                    lyrCenter and lyrSrc to srclayer
    for (vector<node*>::iterator iter = srcrparts->begin(); iter != srcrparts->end(); ++iter) {
        node* p = *iter;
        part_data* pd = (part_data*)p->data;
        part_data* newd = new part_data(*pd);
        node* newn = newlib->add_node(newd);
        node* c = p->get_neighbor(tocenter);
        part_data* cd = (part_data*)c->data;
        vector<node*> path;

        newd->layer = layer;
        destparts->push_back(newn);
        newlib->add_edge((*srcparts)[o2n[newd->type]], newn, tocenterback, tocenter1);
        newlib->add_edge_2(newn, (*origsrcparts)[cd->type], new part_data_2c(clist[pd->type]), tocenter);

        // path from (*srcparts)[o2n[newd->type]] to (*origsrcparts)[cd->type]
        // 
        if (newlib->find_path(path, (*srcparts)[o2n[newd->type]], toprev, (*origsrcparts)[cd->type])) {
            size_t psize = path.size();

            if (psize > 2) {
                node* src = path[psize - 3];
                node* dest = path[psize - 2];

                if (src->is_neighbor(dest, tosrc)) {  // equiv to !is_neighbor(dest, tocenter)
                    foreach_neighbor(p, tosrcM, niter) {
                        lib_data* nd = (lib_data*)neighbor_node_data(niter);
                        part_data_2c* ed = (part_data_2c*)neighbor_edge_data(niter);

                        if (!connected_2c(src, dest, ed, tosrcM))
                            newlib->add_edge_2(src, dest, new part_data_2c(*ed), tosrcM);
                    }
                }

            }
        }

        foreach_neighbor(p, tosrc, niter) {
            lib_data* nd = (lib_data*)neighbor_node_data(niter);
            part_data_2* ed = (part_data_2*)neighbor_edge_data(niter);
            newlib->add_edge_2(newn, (*origsrcparts)[nd->type], new part_data_2(*ed), tosrc);
        }
        foreach_neighbor(p, tosrcM, niter) {
            lib_data* nd = (lib_data*)neighbor_node_data(niter);
            part_data_2c* ed = (part_data_2c*)neighbor_edge_data(niter);
            newlib->add_edge_2(newn, (*origsrcparts)[nd->type], new part_data_2c(*ed), tosrcM);
        }
        foreach_neighbor(p, toprev, niter) {
            lib_data* nd = (lib_data*)neighbor_node_data(niter);
            part_data_2c* ed = (part_data_2c*)neighbor_edge_data(niter);
            newlib->add_edge_2(newn, (*origsrcparts)[nd->type], new part_data_2c(*ed), toprev);
        }
     }


    

    
}

int simsum = 0;
int allparts = 0;

// Get partition of layer O - 1 with respect to the object layer O.
int object_layer_partition(vector<int>& partition, part_lib* library, int O)
{
    vector<node*>& parts = library->parts[O];
    vector<node*>& partsm1 = library->parts[O - 1];
    int srcname = atom("lyrSrc");

    partition.resize(partsm1.size());
    for (size_t i = 0; i < partition.size(); ++i) 
        partition[i] = i;

    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* pc = *iter;
        lib_data* pcd = (lib_data*)pc->data;

        foreach_neighbor(pc, srcname, niter) {
            node* p = neighbor_node(niter);
            lib_data* pd = (lib_data*)p->data;
            part_data_2c* ed = (part_data_2c*)neighbor_edge_data(niter);
            set<int> typeset;
            int t = INT_MAX;

            library->get_similar_types(pd->layer, pd->type, pcd->layer, pcd->type, ed->x, ed->y, 0.1, typeset);        
            typeset.insert(pd->type);
            
            simsum += typeset.size();
            ++allparts;
            for (set<int>::iterator iter = ++typeset.begin(); iter != typeset.end(); ++iter) {
                int nt = *iter;
                if (partition[nt] != nt && partition[nt] < t) t = partition[nt];
            }
            if (t == INT_MAX) 
                t = partition[*(typeset.begin())];
            for (set<int>::iterator iter = ++typeset.begin(); iter != typeset.end(); ++iter)
                partition[*iter] = t;
        }
    }

    map<int, int> old2new;
    map<int, int>::iterator miter;
    int ccount = 0;

    for (size_t i = 0; i < partition.size(); ++i) {
        miter = old2new.find(partition[i]);
        if (miter != old2new.end()) partition[i] = miter->second;
        else {
            old2new.insert(pair<int, int>(partition[i], ccount));
            partition[i] = ccount++;
        }
    }
    return ccount;

}

void test1_merge()
{
    part_lib* library;
    read_library("C:\\work\\data\\multiclass\\LabelMe\\test\\olibv-cat22.plb", library);
    //read_library("C:\\work\\data\\multiclass\\mug-ferrari\\olibv-original.plb", library);
    vector<int> partition;
    int c;

    c = object_layer_partition(partition, library, 5);

    cerr << "Sum = " << simsum << " avg = " << (double)simsum/allparts << endl;
    for (int i = 0; i < c; ++i) {
        cerr << "Cluster " << i << ": ";
        for (size_t j = 0; j < partition.size(); ++j) 
            if (partition[j] == i) cerr << j << ' ';
        cerr << endl;
    }

    delete library;
}


struct magic_function : public unary_double_function {
    virtual double operator()(double a) const { return a*a; }
};


double test_function(double a, const unary_double_function& uf)
{
    return uf(a) + 1;
}


struct function_magique : public unary_double_function {
    virtual double operator()(const double& d) const { return ::sqrt(d); }
};

struct test_class {
    virtual int get() const { return 1; }
};

struct test_class_2 : public test_class {
    int x;
    test_class_2(int xx) : x(xx) { }
    virtual int get() const { return x; }
};

void test_print(const vector<test_class>& cl)
{
//    cerr << "Result of get = " << cl[0]->get() << endl;
}

// Cluster points until maximal diameter^2 of cluster does not exceed maxdiam2.
void diameter_clustering(vector<vector<ipoint2> >& clusters, int maxdiam2)
{
    cerr << clusters.size() << endl;
    while ((int)clusters.size() > 1) {
        vector<vector<ipoint2> > memclusters = clusters;
        int diam2 = 0;

        cerr << ".";
        hierarchical_clustering_iter(clusters, ipoint2_set_distance);
        for (int c = 0; c < (int)clusters.size(); ++c) {
            int d = diameter2<int>(clusters[c].begin(), clusters[c].end());
            if (d > diam2) diam2 = d;
        }
        if (diam2 > maxdiam2){
            clusters = memclusters;
            break;
        }
    }
}





/*

foreach (n in layer 5 and in validation res')
    stat(n.m) += |rec(n) isct gtr|/|gtr|
sort stat(m)
foreach (model image res)
    model = {}
    models = {}
    augment_model(models, res, model, stat(m), dupmap, maxmodels, validation set, ...)

********

augment_model(models, model, stat(m), maxmodels, validation set)
    match model to res
    forall m (sorted by stat)
        match m according to duplets to res
        if match is valid and |rec(match) isct rec(model)|/rec(match) < thresh then
            add model + m to models
            augment_modes(models, model + m, stat(m), maxmodels, validation set, ...)



*/


/*void test1()
{
    layer1_result* res, * res2;
    clock_t start, end;

    read_layer1_result(res, "D:\\work\\data\\multiclass\\mug-poeticon\\layerx\\final\\test_Video13-002_0.ly7");
    start = clock();
    res->add_reconstruction_edges_fwd(4);
    end = clock();
    cout << "Time for leq_fwd: " << (double)(end - start)/CLOCKS_PER_SEC << "s" << endl;

    read_layer1_result(res2, "D:\\work\\data\\multiclass\\mug-poeticon\\layerx\\final\\test_Video13-002_0.ly7");
    start = clock();

    int to_prev = atom("toPrevLayer");
    int to_0 = atom("toLayer0");

    for (list<node*>::iterator iter = res2->nodes.begin(); iter != res2->nodes.end(); ++iter) {
        if (node_layer(*iter) == 4)
            res2->recurse_and_link(*iter, to_prev, to_0);
    }
    end = clock();
    cout << "Time for r&l: " << (double)(end - start)/CLOCKS_PER_SEC << "s" << endl;
    int total = 0, diff = 0;

    for (list<node*>::iterator iter = res->nodes.begin(), iter2 = res2->nodes.begin(); iter != res->nodes.end() && iter2 != res2->nodes.end(); ++iter, ++iter2) {
        if (node_layer(*iter) != 4) 
            continue;

        set<node*> nset, nset2;
        set<ipoint2> pset, pset2;

        (*iter)->get_neighbor_set(to_0, nset);
        (*iter2)->get_neighbor_set(to_0, nset2);
        node_set_to_point_set(pset, nset.begin(), nset.end());
        node_set_to_point_set(pset2, nset2.begin(), nset2.end());
        if (pset != pset2) {
            cerr << "Difference detected!" << endl;
            cerr << "Sizes: " << pset.size() << ", " << pset2.size() << endl;
            diff++;
        }
        total++;
    }
    cerr << "Total: " << total << "  difference: " << diff << endl;

    delete res2;
    delete res;

}*/


// Return the number of points (shape positions) where 'res1' is different
// from 'res' on layer 'layer'.
int layer_0_difference(layer1_result* res, layer1_result* res1, double tolerance)
{
    if (!res->grid(0)) res->init_grid(0);

    vector<node*>& s_nodes = res1->shape_nodes[0];
    int result = 0;

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;
        bool found = false;

        for (int dx = -2; dx < 3 && !found; ++dx) 
            for (int dy = -2; dy < 3 && !found; ++dy) {
                node* n0 = res->node_at_safe(nd->x + dx, nd->y + dy, 0);        

                while (n0 != nullptr) {
                    layer1_data* n0d = (layer1_data*)n0->data;

                    if (abs(n0d->r(R_RESPONSE) - nd->r(R_RESPONSE)) > tolerance)
                        break;
                    if (n0d->m == nd->m) {
                        found = true;
                        break;
                    }
                    n0 = n0d->next;
                }
            }
        if (!found) ++result;
    }
    return result;
}


void test1()
{
    layer1_result* res, * res1;
    clock_t start, end;
    int* leak = new int;

    *leak = 100;
    read_layer1_result(res, "D:\\work\\data\\multiclass\\poeticon-video13\\layer1\\Video13-003_0.ly1");
    read_layer1_result(res1, "D:\\work\\data\\multiclass\\poeticon-video13\\layer1\\Video13-004_0.ly1");
    start = clock();
    //inhibit_layer(res, 6, RR_RESPONSE, 100, 0.4);
    cout << "Result: " << layer_0_difference(res, res1, 0.5) << "/" << res1->shape_nodes[0].size() << endl;
    end = clock();
    cout << "Time for operation: " << (double)(end - start)/CLOCKS_PER_SEC << "s" << endl;

    delete res1;
    read_layer1_result(res1, "D:\\work\\data\\multiclass\\poeticon-video13\\layer1\\Video13-012_0.ly1");

    start = clock();
    //inhibit_layer(res, 6, RR_RESPONSE, 100, 0.4);
    cout << "Result: " << layer_0_difference(res, res1, 0.5) << "/" << res1->shape_nodes[0].size() << endl;
    end = clock();
    cout << "Time for operation: " << (double)(end - start)/CLOCKS_PER_SEC << "s" << endl;

    delete res;
    delete res1;

}

/*
void test1()
{
    layer1_result* res, * res2;
    clock_t start, end;

    read_layer1_result(res, "D:\\work\\data\\multiclass\\mug-poeticon\\layerx\\final\\test_Video13-002_0.ly7");
    start = clock();
    res->add_reconstruction_edges_fwd(4);
    end = clock();
    cout << "Time for leq_fwd: " << (double)(end - start)/CLOCKS_PER_SEC << "s" << endl;

    read_layer1_result(res2, "D:\\work\\data\\multiclass\\mug-poeticon\\layerx\\final\\test_Video13-002_0.ly7");
    start = clock();

    int to_prev = atom("toPrevLayer");
    int to_0 = atom("toLayer0");

    for (list<node*>::iterator iter = res2->nodes.begin(); iter != res2->nodes.end(); ++iter) {
        if (node_layer(*iter) == 4)
            res2->recurse_and_link(*iter, to_prev, to_0);
    }
    end = clock();
    cout << "Time for r&l: " << (double)(end - start)/CLOCKS_PER_SEC << "s" << endl;
    int total = 0, diff = 0;

    for (list<node*>::iterator iter = res->nodes.begin(), iter2 = res2->nodes.begin(); iter != res->nodes.end() && iter2 != res2->nodes.end(); ++iter, ++iter2) {
        if (node_layer(*iter) != 4) 
            continue;

        set<node*> nset, nset2;
        set<ipoint2> pset, pset2;

        (*iter)->get_neighbor_set(to_0, nset);
        (*iter2)->get_neighbor_set(to_0, nset2);
        node_set_to_point_set(pset, nset.begin(), nset.end());
        node_set_to_point_set(pset2, nset2.begin(), nset2.end());
        if (pset != pset2) {
            cerr << "Difference detected!" << endl;
            cerr << "Sizes: " << pset.size() << ", " << pset2.size() << endl;
            diff++;
        }
        total++;
    }
    cerr << "Total: " << total << "  difference: " << diff << endl;

    delete res2;
    delete res;
*/
    //vectorn<int>::type3 v;

    //at(v, 10, 2, 1, 0) = 100;

    //cerr << at(v, 1, 1, 5, 0) << endl;
    //cerr << at(v, 10, 2, 1, 0) << endl;

    //foreach_vector3(v, i, j, k) {
    //    cerr << i << ',' << j << ',' << k << ':' << v[i][j][k] << endl;
    //}


    //list<string> files;
    //string dir = "D:\\work\\data\\multiclass\\mug-ferrari\\layerx\\";
    //list<streamed_pointer> pl;


    //list_from_file(files, "C:\\work\\data\\multiclass\\mug-ferrari\\train_files_5_debug.txt", dir);
    //for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
    //    layer1_result* res;

    //    read_layer1_result(res, dir + *file);
    //    if (res != nullptr) {
    //        pl.push_back(streamed_pointer(res));

    //        delete res;
    //    }
    //}
    //
    //streamed_pointer sp = pl.back();

    //pl.clear();


    ////////////////////////////////////////////////////////////////////////////////////////////

    //list<string> testlist;
    //vector<int> v;

    //testlist.push_back("work"); // 0 
    //testlist.push_back("data"); // 1
    //testlist.push_back("multiclass"); // 2
    //testlist.push_back("mug");  // 3 
    //testlist.push_back("ferrari");  // 4
    //testlist.push_back("layerx");  // 5

    //v.push_back(2);
    ////v.push_back(0);
    ////v.push_back(4);
    ////v.push_back(3);
    ////v.push_back(7);

    //keep_items(testlist, v);

    //for (list<string>::iterator i = testlist.begin(); i != testlist.end(); ++i) {
    //    cerr << *i << endl;
    //}

    //return;

    //string dir = "D:\\work\\data\\multiclass\\mug-ferrari\\layerx\\";
    //list<string> files;
    //config_dictionary cfg("C:\\work\\data\\multiclass\\mug-ferrari\\new_o_learning.cfg");

    //o_learning olearner(cfg);

    //cerr << "Duplet Statistics" << endl;

    //list_from_file(files, "C:\\work\\data\\multiclass\\mug-ferrari\\train_files_5_debug.txt", dir);
    //for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
    //    layer1_result* res;

    //    read_layer1_result(res, dir + *file);
    //    if (res == nullptr) {
    //        cerr << "File " << *file << " can not be opened!" << endl;
    //    } else {
    //        cerr << "Processing file " << *file;

    //        olearner.update_duplet_statistics(res);

    //        delete res;
    //        cerr << " done." << endl;
    //    }
    //}

    //cerr << "Validation sets" << endl;

    //files.clear();
    //list_from_file(files, "C:\\work\\data\\multiclass\\mug-ferrari\\tl_train_files_5_debug.txt", dir);
    //for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
    //    layer1_result* res;

    //    read_layer1_result(res, dir + *file);
    //    if (res == nullptr) {
    //        cerr << "File " << *file << " can not be opened!" << endl;
    //    } else {
    //        cerr << "Processing file " << *file;
    //        list<irectangle2> gtrs;

    //        read_groundtruth(gtrs, dir + *file, "mug", ".groundtruth");
    //        olearner.add_to_validation_set(res, gtrs);

    //        delete res;
    //        cerr << " done." << endl;
    //    }
    //}

    //cerr << "Make models" << endl;

    //files.clear();
    //list_from_file(files, "C:\\work\\data\\multiclass\\mug-ferrari\\train_files_5_debug.txt", dir);
    //for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
    //    layer1_result* res;

    //    read_layer1_result(res, dir + *file);
    //    if (res == nullptr) {
    //        cerr << "File " << *file << " can not be opened!" << endl;
    //    } else {
    //        cerr << "Processing file " << *file;

    //        olearner.learn_models(res);

    //        delete res;
    //        cerr << " done." << endl;
    //    }
    //    //break;
    //}
    //olearner.finalize();

    //olearner.get_library()->save("C:\\work\\data\\multiclass\\mug-ferrari\\libx.plb");

    //
    //cout << "{";
    //for (tupmap_t::iterator iter = dupmap.begin(); iter != dupmap.end(); ++iter) {
    //    if (iter != dupmap.begin()) cout << ","; 
    //    cout << "{";

    //    for (ip2_vector::iterator viter = iter->second.begin(); viter != iter->second.end(); ++viter) {
    //        if (viter != iter->second.begin()) cout << ",";
    //        cout << "{" << viter->x << "," << viter->y << "}";
    //    }
    //    cout << "}" << endl;

    //    int ss;
    //    cin > > ss;
    //}
    //cout << "}" << endl;
//}

void test1_clusters()
{
    vector<vector<ipoint2> > clusters;
    
    int pts[][2] = {{-9, -9}, {-9, 
        2}, {-8, -1}, {-6, -9}, {-4, -7}, {-4, -1}, {-3, -10}, {0, 
        4}, {2, -7}, {2, 5}, {3, -8}, {3, 4}, {3, 8}, {3, 9}, {6, 7}, {6, 
        10}, {7, -8}, {8, -7}, {10, -5}};

    for (int i = 0; i < 19; ++i) {
        clusters.push_back(vector<ipoint2>());
        clusters.back().push_back(ipoint2(pts[i][0], pts[i][1]));
    }
    while ((int)clusters.size() > 1) {
        vector<vector<ipoint2> > memclusters = clusters;
        int diam2 = 0;

        hierarchical_clustering_iter(clusters, ipoint2_set_distance);
        for (int c = 0; c < (int)clusters.size(); ++c) {
            int d = diameter2<int, vector<ipoint2>::iterator>(clusters[c].begin(), clusters[c].end());
            if (d > diam2) diam2 = d;
        }
        if (diam2 > 9*9){
            clusters = memclusters;
            break;
        }
    }

    cout << "{";
    for (int c = 0; c < clusters.size(); ++c) {
        if (c > 0) cout << ","; 
        cout << "{";
        for (int i = 0; i < clusters[c].size(); ++i) {
            if (i > 0) cout << ",";
            cout << "{" << clusters[c][i].x << "," << clusters[c][i].y << "}";
        }
        cout << "}";
    }
    cout << "}" << endl;




    //layer1_result* res0, * res1;

    //read_layer1_result(res0, "D:\\work\\data\\multiclass\\mug-poeticon\\layerx\\test_Picture07_0.ly7");
    //read_layer1_result(res1, "D:\\work\\data\\multiclass\\mug-poeticon\\layerx\\test_Picture07_1.ly7");
    //res0->merge(res1, 100);
    //save_layer1_result(res0, "c:\\work\\merged.ly7");

    //delete res0;
    //delete res1;
}


//void test1()
//{
//    list<string> fnames;
//    string srcdir = "G:\\recognition\\edge\\teledyne\\layer3\\";
//    int layer_index = 1;
//    int slices = 4;
//    int levels = 2;
//
//    fnames.push_back("43949781_257_513_513_769_0.ly3");
//    fnames.push_back("43949781_257_513_513_769_1.ly3");
//    fnames.push_back("43949781_257_513_513_769_2.ly3");
//
//

    //cerr << "TEST1" << endl;

    //vector<test_class_2> vec;

    //test_print(vec);
    
    /*const char* fname = "C:\\science\\data\\lowlevel\\simpsons\\bart1.gif";

    ifstream f(fname, ios::in | ios::binary);
    struct stat results;

    stat(fname, &results);
    cerr << "File size: " << results.st_size << endl;   
    
    char* buffer = (char*)malloc(sizeof(char) * results.st_size);
    if (!f.read(buffer, results.st_size)) cerr << "An error occurred!" << endl;
    else cerr << "Successfully read the above value of bytes!" << endl;
    Magick::Blob blb;

    blb.updateNoCopy(buffer, sizeof(char) * results.st_size);


    Magick::Image im(blb);
    im.write("c:\\bart.png");

    free(buffer);*/

    
//}

void test1_some_old_version()
{
#if defined WIN32 | defined WIN64
    cerr << "TEST1" << endl;

    response_map rmap;
    int hit = 0;

    rmap.set_response(0, 1.10);
    rmap.set_response(2, 1.12);
    rmap.set_response(3, 1.13);
    rmap.set_response(7, 1.14);
    rmap.set_response(9, 1.15);

    clock_t start = clock();
    for (int i = 0; i < 10; ++i) {
        int x = random_int(0, 9);
        if (rmap.get_response(x, -1.0) == 1.10) ++hit;
    }
    clock_t end = clock();

    cerr << "hit = " << hit << endl;
    cerr << "Processing time " << ((double)(end - start))/CLOCKS_PER_SEC << " sec." << endl;


    online_distribution normal;

    normal.new_data(1.0);
    normal.new_data(2.0);
    normal.new_data(3.0);
    normal.new_data(4.0);
    normal.new_data(5.0);
    normal.new_data(4.0);
    normal.new_data(3.0);
    normal.new_data(2.0);
    normal.new_data(1.0);

    cerr << "Mean: " << normal.get_mean() << endl;
    cerr << "Variance: " << normal.get_variance() << endl;

    //cerr << "Magic function: " << test_function(2.0, magic_function()) << endl;

    unary_double_power f;

    hit = 0;
    cerr << endl << endl << "Power performance" << endl;
    start = clock();

    for (int i = 0; i < 10; ++i) {
        double x = random_real();
        if (x < 3.0) x = 1.0;
        //double y = (x == 1.0) ? x : ::pow(2.0, x);
        f.power = x;
        double y = f(2.0);
        //cerr << y << ' ';
        if (y < 1.5) ++hit;
    }
    end = clock();

    cerr << "hit = " << hit << endl;
    cerr << "Processing time " << ((double)(end - start))/CLOCKS_PER_SEC << " sec." << endl;

    function_magique mf;

    unary_double_function* df = &mf;

    cerr << "Magique?: " << (*df)(2.0) <<  endl;

    MEMORYSTATUSEX statex;

    statex.dwLength = sizeof (statex);

    GlobalMemoryStatusEx(&statex);
    cerr << (1U << 31) << endl;
    cerr << statex.ullTotalPhys/1024 << " Kb of total pysical memory" << endl;
    cerr << statex.ullAvailPhys/1024 << " Kb of physical memory free" << endl;
    cerr << statex.ullTotalPageFile/1024 << " Kb of total page file memory" << endl;
    cerr << statex.ullAvailPageFile/1024 << " Kb of available page file memory" << endl;
    cerr << "Memory load = " << statex.dwMemoryLoad << endl;
    
    PROCESS_MEMORY_COUNTERS pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));

    cerr << "This process: WorkingSetSize: " << pmc.WorkingSetSize/1024 << endl;

    list<layer1_result_ptr> rlist;

    streamable_pointer::set_max_mem_objects(150);

    matrix<double> m(10, 30, 1.0);

    for (int i = 0; i < 3; ++i) {
        layer1_result* res;
        read_layer1_result(res, "G:\\recognition\\edge\\simpsons\\layer2\\homer20simpson_0.ly2");
        res->add_reconstruction_edges(1);
        rlist.push_back(layer1_result_ptr(res));
        //if (i % 100 == 0) {
        //    int x;
        //    cin > > x;
        //    if (x > 1) break;
        //}
        //layer1_result* r2 = rlist.back();
    }
    //for (list<layer1_result_ptr>::iterator rliter = rlist.begin(); rliter != rlist.end(); ++rliter) {
    //    layer1_result* r2 = *rliter;
    //    r2->save("C:\\work\\x.ly2");
    //}

#endif

}




void test1_xxx()
{
    layer1_result* res;

    read_layer1_result(res, "tmpimg_0.bmp");
    if (res == nullptr)
        cerr << "There has been an error reading lyX file!" << endl;
    /*Magick::Image img;
    
    try {
        img = Magick::Image("c:\\work\\data\\polar\\cawm_0.ly1");
    } catch (Magick::Exception& e) {
        cerr << e.what() << endl;
    }
    cerr << "Image valid? " << img.isValid() << endl;*/

    
}

void test1_2009()
{


    /*

    list<irectangle2> rl;

    irectangle2(0, 0, 6, 6).split(rl, 3, 3);
    cout << '{';
    for (list<irectangle2>::iterator iter = rl.begin(); iter != rl.end(); ++iter) {
        cout << ',';
        iter->mma_write(cout);
    }
    cout << '}';
    */

    //layer1_result* res;

    //read_layer1_result(res, "D:\\work\\data\\multiclass\\mug-ferrari\\layerx\\mugs_learn_image_71_1.ly5");

    //delete res;


    
/*    part_lib* lib;

    read_library("c:\\work\\data\\multiclass\\mug-ferrari\\olibv-original.plb", lib);
    for (int i = 0; i < lib->layer_count; ++i) {
        cerr << "Layer " << i << " size = " << lib->parts[i].size() << endl;
    }
    for (vector<node*>::iterator iter = lib->parts[4].begin(); iter != lib->parts[4].end(); ++iter) {
        node* p = *iter;
        set<node*> s;
        p->get_neighbor_set(atom("lyrCenterBack"), s);
        cerr << "lyrCenterBack: " << s.size() << endl;
        s.clear();
        p->get_neighbor_set(atom("lyrSrc"), s);
        cerr << "lyrSrc: " << s.size() << endl;
    }

    delete lib;
    */

    part_lib* library;

    //read_library("C:\\work\\data\\multiclass\\LabelMe\\test\\olibv-cat22.plb", library);
    read_library("C:\\work\\data\\multiclass\\mug-ferrari\\olibv-original.plb", library);

    if (library == nullptr) return;

    vector<node*>& partsm1 = library->parts[4];
	vector<node*>& parts = library->parts[5];
    int backname = atom("lyrCenterBack");
    int srcname = atom("lyrSrc");
    int centername = atom("lyrCenter");
    //int rect0dim = 4;
	list<layer1_result*> rlist;
    int count = 0;

    int partition_classes = (int)partsm1.size();
    vector<int> partition(partsm1.size(), 0);

    // make a simple partition by "mod"
    for (size_t i = 0; i < partsm1.size(); ++i) partition[i] = i % partition_classes;


    //  int maxxsize = 0, maxysize = 0;
    //  for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
    //      node* p = *iter;
    //  	irectangle2 box;
    //// Get the bounding box of the subparts
		    //box.eat(0, 0);
    //      foreach_neighbor(p, srcname, siter) {
    //          part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
    //          
	//	box.eat(ed->x/rect0dim, ed->y/rect0dim);
    //      }
    //      if (box.x_dim() > maxxsize) maxxsize = box.x_dim();
    //      if (box.y_dim() > maxysize) maxysize = box.y_dim();
    //  }
    //  cerr << "Max x size = " << maxxsize << "  max y size = " << maxysize << endl;

      //  return;

	// For each part of layer X make a corresponding layer1_result
    for (vector<node*>::iterator iter = parts.begin() + 8; iter != parts.end() && count < 8; ++iter, ++count) {
		// break; // REMOVE!!!!!!!!!!!!!!!!!!!!!!!!

        node* p = *iter;
        lib_data* pd = (lib_data*)p->data; //%%%
		layer1_result* r = new layer1_result();
		irectangle2 box;
		node* n;

		// Get the bounding box of the subparts
		box.eat(0, 0);
        foreach_neighbor(p, srcname, siter) {
            part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
            
			box.eat(ed->x, ed->y);
        }

		// normalize the size of box s.t. (0, 0) is on even 
		//if (box.ll.x < 0) --box.ll.x;
		//if (box.ll.y < 0) --box.ll.y;
		//if (box.ll.x % 2 != 0) --box.ll.x;
		//if (box.ll.y % 2 != 0) --box.ll.y;

		// add nodes to r; shape_nodes = shape_nodes_inhib are filled
		// s.t. the first element is the "center"
        node* c = p->get_neighbor(centername);
        lib_data* cd = (lib_data*)c->data;
        //set<int> typeset;

        r->new_grid(box.x_dim() + 1, box.y_dim() + 1, 0);
        r->shape_nodes.push_back(vector<node*>());
        r->shape_nodes_inhib.push_back(vector<node*>());
        //typeset.insert(partition[cd->type]);
        n = r->add_grid_node(new layer1_data(1.0, partition[cd->type]), -box.ll.x, -box.ll.y, 0);
		r->shape_nodes[0].push_back(n);
		r->shape_nodes_inhib[0].push_back(n);
        foreach_neighbor(p, srcname, siter) {
            part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
            lib_data* pnd = (lib_data*)neighbor_node_data(siter);
            
            //typeset.clear();
            //library->get_similar_types(pnd->layer, pnd->type, pd->layer, pd->type, ed->x, ed->y, 1.0, typeset);
            //typeset.insert(partition[pnd->type]);
            n = r->add_grid_node(new layer1_data(1.0, partition[pnd->type]), 
				ed->x - box.ll.x,
                ed->y - box.ll.y, 
				0);
			r->shape_nodes[0].push_back(n);
			r->shape_nodes_inhib[0].push_back(n);
        }

		rlist.push_back(r);
    }

    // D E B U G
    // save rlist to mathematica

    save_rlist("c:\\temp\\init.m", rlist, 0);
    // D E B U G
	
	// "Faked" test input
	// vvvvvvvvvvvvvvvvvvvvv
	if (false) {
        layer1_result* r = new layer1_result();
		node* n;
        
        r->new_grid(8, 8, 0);
        r->shape_nodes.push_back(vector<node*>());
        r->shape_nodes_inhib.push_back(vector<node*>());

		n = r->add_grid_node(new layer1_data(1.0, 0), 4, 4, 0);
		r->shape_nodes[0].push_back(n);
		r->shape_nodes_inhib[0].push_back(n);

		n = r->add_grid_node(new layer1_data(1.0, 0), 5, 5, 0);
		r->shape_nodes[0].push_back(n);
		r->shape_nodes_inhib[0].push_back(n);

		n = r->add_grid_node(new layer1_data(1.0, 0), 1, 6, 0);
		r->shape_nodes[0].push_back(n);
		r->shape_nodes_inhib[0].push_back(n);

		n = r->add_grid_node(new layer1_data(1.0, 0), 1, 2, 0);
		r->shape_nodes[0].push_back(n);
		r->shape_nodes_inhib[0].push_back(n);

		n = r->add_grid_node(new layer1_data(1.0, 0), 5, 2, 0);
		r->shape_nodes[0].push_back(n);
		r->shape_nodes_inhib[0].push_back(n);

		rlist.push_back(r);
	}
	// ^^^^^^^^^^^^^^^^^^^^^

    part_lib* rlib = new part_lib(0);

    // add layer 0 to the library ---------------------------------------------
    rlib->layer_count = 0;
    for (int i = 0; i < partition_classes; ++i) 
        rlib->parts[0].push_back(rlib->add_node(new rpart_data(0, i, 0, 0, 0)));

    rlib->save_all("c:\\temp\\0.png", 1, 0, -1, false, false);

    // add layer 1 to the library ---------------------------------------------
    learn_layer(rlib, rlist, 1);

//    rlib->save_all("c:\\temp\\1.png", 2, 0, 50, false, true);        

    // add layer 1 to the results
    for (list<layer1_result*>::iterator riter = rlist.begin(); riter != rlist.end(); ++riter) 
        add_layer(*riter, rlib, 1);
        
    save_rlist("c:\\temp\\1.m", rlist, 1);

    // add layer 2 to the library ---------------------------------------------
    learn_layer(rlib, rlist, 2);

  //  rlib->save_all("c:\\temp\\2.png", 3, 0, 50, false, true);        

    // add layer 1 to the results
    for (list<layer1_result*>::iterator riter = rlist.begin(); riter != rlist.end(); ++riter) 
        add_layer(*riter, rlib, 2);
        
    save_rlist("c:\\temp\\2.m", rlist, 2);

    // cleanup ----------------------------------------------------------------

	for (list<layer1_result*>::iterator riter = rlist.begin(); riter != rlist.end(); ++riter) {
		delete *riter;
	}
    delete rlib;
    delete library;
}

void test1_semiold()
{
    //_CrtSetBreakAlloc(3690);



    cerr << "log2 of 1244 = " << ilog2_x86(8) + 1 << endl;
    int j = 0;
    clock_t start, stop;

    start = clock();
    for (unsigned i = 0; i < 10000; ++i) {
        j += ilog2_x86(i);
    }
    stop = clock();
    cerr << j;
    cerr << " Processing time ilog2_x86: " << (double)(stop - start)/CLOCKS_PER_SEC << endl;

}

// Choses a minimal (approximately) subset of the set 's'
// which covers the set 's' according to the neighborhood
// given by the edges with name 'ename'.
//void reduce_node_set(set<node*>& result, const set<node*>& s, int ename)
//{
//    set<node*> uncovered(s);
//
//    result.clear();
//    while (!uncovered.empty()) {
//        node* bestn = nullptr;
//        int bestnbsize = -1;
//
//        for (set<node*>::iterator iter = uncovered.begin(); iter != uncovered.end(); ++iter) {
//            node* n = *iter;
//            int nbsize = 0;
//
//            foreach_neighbor(n, ename, niter) {
//                if (uncovered.find(neighbor_node(niter)) != uncovered.end()) ++nbsize;
//            }
//            if (nbsize > bestnbsize) { bestnbsize = nbsize; bestn = n; }
//        }   
//
//        result.insert(bestn);
//        foreach_neighbor(bestn, ename, niter) {
//            uncovered.erase(neighbor_node(niter));
//        }
//        uncovered.erase(bestn);
//    }
//}


//template<class T> bool operator<(const vector<T>& a, const vector<T>& b)
//{
//    vector<T>::const_iterator itera = a.begin(), iterb = b.begin();
//
//    for (; itera != a.end() && iterb != b.end(); ++itera, ++iterb) {
//        if (*itera < *iterb) return true; 
//        if (*itera > *iterb) return false;
//    }
//    return iterb != b.end();
//}

irectangle2* find_rectangle(list<irectangle2>& rlist, int x, int y)
{
    for (list<irectangle2>::iterator iter = rlist.begin(); iter != rlist.end(); ++iter) {
        if (iter->inside_uropen(x, y)) return &(*iter);
    }
    return nullptr;
}

// r_part_data
////////////////

struct r_part_data : public lib_data {
public:
    vector<irectangle2> rlist;

    r_part_data() : rlist(), lib_data() { }
    r_part_data(int l, int t, const vector<irectangle2>& rl) : rlist(rl), lib_data(l, t) { }

    virtual void copy_to(streamable* p, cloner& cl) 
    {
        lib_data::copy_to(p, cl);
        ((r_part_data*)p)->rlist = rlist;
    }

    virtual streamable* make_instance() const { return new r_part_data(); }    
    virtual void write_to_stream(ostreamer& os)
    {
        lib_data::write_to_stream(os);
        os.write(rlist);
    }
    virtual void read_from_stream(istreamer& is)
    {
        lib_data::read_from_stream(is);
        is.read(rlist);
    }

};

// ---> to part_lib
void add_r_part_to_library(part_lib* library, int lyr, int ctype, const vector<irectangle2>& rset, const vector<node*>& nset)
{
    if (lyr < 1) return;
    while (lyr + 1 >= library->layer_count) ++library->layer_count;

    int prevname = atom("lyrPrev");
    int srcname = atom("lyrSrc");
    int srcmname = atom("lyrSrcM");
    vector<node*>& parts = library->parts[lyr];
    vector<node*>& partsp1 = library->parts[lyr + 1];
    vector<node*>& partsm1 = library->parts[lyr - 1];
    node* c = partsm1[ctype];
    node* p = library->add_node(new r_part_data(lyr, (int)parts.size(), rset), R_PART_ATTR);

    // Add r_part node
    parts.push_back(p);
    library->add_edge(p, c, atom("lyrCenter"), atom("lyrCenterBack"));
    library->add_edge_2(p, c, new part_data_2c(0, 0), prevname);

    //set<node*> ns;
    //c->get_neighbor_set(atom("lyrCenterBack"), ns); cerr << ns.size() << endl;
    
    // Add nodes do lyr + 1
    for (vector<node*>::const_iterator iter = nset.begin(); iter != nset.end(); ++iter) {
        node* oldp = *iter;
        node* newp = library->add_node(new part_data(*((part_data*)oldp->data)));
        part_data* newpd = (part_data*)newp->data;

        // Set type copy the neighbors
        newpd->type = (int)partsp1.size();
        newpd->layer = lyr + 1;
        partsp1.push_back(newp);
        //library->add_edge(newp, p, atom("lyrCenter"), atom("lyrCenterBack"));
        library->add_edge(p, newp, atom("lyrCenterBack"), atom("lyrCenter1"));
        library->add_edge(newp, c, atom("lyrCenter"));

        foreach_neighbor(oldp, srcname, niter) {
            lib_data* nd = (lib_data*)neighbor_node_data(niter);
            part_data_2* ed = (part_data_2*)neighbor_edge_data(niter);
            
            library->add_edge_2(newp, partsm1[nd->type], new part_data_2(*ed), srcname);
        }
        foreach_neighbor(oldp, srcmname, niter) {
            lib_data* nd = (lib_data*)neighbor_node_data(niter);
            part_data_2c* ed = (part_data_2c*)neighbor_edge_data(niter);
            
            library->add_edge_2(newp, partsm1[nd->type], new part_data_2c(*ed), srcmname);
        }
        foreach_neighbor(oldp, prevname, niter) {
            lib_data* nd = (lib_data*)neighbor_node_data(niter);
            part_data_2c* ed = (part_data_2c*)neighbor_node(niter);

            library->add_edge_2(newp, partsm1[nd->type], new part_data_2c(*ed), prevname);
        }
        

    }

}

struct type_rectangle {
    set<int> tset;
    irectangle2 rect;

    type_rectangle(int t, const irectangle2& r) : tset(), rect(r) { tset.insert(t); }
    type_rectangle(const set<int>& ts, const irectangle2& r) : tset(ts), rect(r) { }
    type_rectangle(const type_rectangle& tr) : tset(tr.tset), rect(tr.rect) { }

    bool operator<(const type_rectangle& r2) const
    {
        if (rect.ll < r2.rect.ll) return true;
        else if (rect.ll > r2.rect.ll) return false;
        else return rect.ur < r2.rect.ur;
    }

    friend ostream& operator<<(ostream& os, const type_rectangle& r) 
    {
        os << '[' << *r.tset.begin() << ',' << r.rect << ']';
        return os;
    }

};

struct tr_vector {
    typedef vector<type_rectangle> trvector_t;
    typedef trvector_t::iterator trviter_t;
    typedef trvector_t::const_iterator trvciter_t;
    typedef pair<int, ipoint2> info_t;  // original type & original center

    set<info_t> info;
    trvector_t vec;

    tr_vector(int t, const ipoint2& p) : info(), vec() { info.insert(info_t(t, p)); }

    void display_info() const;
    void insert_into(tr_vector& v2) const;
    void grow(int g);
    bool subvector(const tr_vector& v2) const;
    bool subvector(const tr_vector& v2, const map<ipoint2, int>& cooc, int maxcooc, int coocthresh) const;
};

typedef list<tr_vector> trvlist_t;

// Insert this[i] into v2[j] if this[i] is a subset of v2[j]
void tr_vector::insert_into(tr_vector& v2) const
{
    for (trvciter_t iter1 = vec.begin(); iter1 != vec.end(); ++iter1) {
        for (trviter_t iter2 = v2.vec.begin(); iter2 != v2.vec.end(); ++iter2) {
            if (iter2->rect.inside(iter1->rect)) 
                iter2->tset.insert(iter1->tset.begin(), iter1->tset.end());
        }
    }
    v2.info.insert(info.begin(), info.end());
}

void tr_vector::display_info() const
{
    cerr << "# of items: " << info.size() << ": ";
    for (set<info_t>::const_iterator iter = info.begin(); iter != info.end(); ++iter) {
        cerr << *iter;
    }
}

void tr_vector::grow(int g)
{
    for (trviter_t iter = vec.begin(); iter != vec.end(); ++iter) 
        iter->rect.grow(g);
}

// Tests whether "v1 subset v2". 
bool tr_vector::subvector(const tr_vector& v2) const
{
    for (trvciter_t iter1 = vec.begin(); iter1 != vec.end(); ++iter1) {
        bool inside = false;

        for (trvciter_t iter2 = v2.vec.begin(); iter2 != v2.vec.end(); ++iter2) {
            if (iter2->rect.inside(iter1->rect)) { inside = true; break; }
        }
        if (!inside) return false;
    }
    return true;
}

// Returns:
//                                \sum cooc(i, j)
//                   (i, j) \in s1 \times s2 and i <= j
//
// cooc(i, i) is defined as 'maxcooc'.
int co_occurrence_sum(const set<int>& s1, const set<int>& s2, const map<ipoint2, int>& cooc, int maxcooc)
{
    int result = 0;

    for (set<int>::const_iterator iter1 = s1.begin(); iter1 != s1.end(); ++iter1) {
        int i = *iter1;
        int maxi = 0;

        for (set<int>::const_iterator iter2 = s2.begin(); iter2 != s2.end(); ++iter2) {
            int j = *iter2;
            
            if (i == j) {
                maxi = maxcooc;
                break;
            }
            if (i < j) {
                map<ipoint2, int>::const_iterator miter = cooc.find(ipoint2(i, j));

                if (miter != cooc.end() && miter->second > maxi) maxi = miter->second;
            }
        }
        result += maxi;
    }
    return result;
}

// Tests whether "v1 subset v2" and checks similarity of types (not yet...)
bool tr_vector::subvector(const tr_vector& v2, const map<ipoint2, int>& cooc, int maxcooc, int coocthresh) const
{
    for (trvciter_t iter1 = vec.begin(); iter1 != vec.end(); ++iter1) {
        bool inside = false;

        for (trvciter_t iter2 = v2.vec.begin(); iter2 != v2.vec.end(); ++iter2) {
            if (iter2->rect.inside(iter1->rect) && co_occurrence_sum(iter1->tset, iter2->tset, cooc, maxcooc) >= coocthresh) {
                inside = true; 
                break; 
            }
        }
        if (!inside) return false;
    }
    return true;
}

bool type_vector_compare(const tr_vector& v1, const tr_vector& v2)
{
    return v1.vec.size() < v2.vec.size();
}

bool type_rectangle_compare(const type_rectangle& r1, const type_rectangle& r2)
{
    return r1.tset.size() < r1.tset.size();
}

part_lib* insert_c2f_layer(part_lib* srclib, int layer, trvlist_t& trvlist)
{
    part_lib* newlib = (part_lib*)srclib->get_copy();

    newlib->insert_empty_layer(layer);

    int srclayer = layer - 1;
    int c2flayer = layer;
    int objlayer = layer + 1;
    vector<node*>& c2fparts = newlib->parts[c2flayer];
    vector<node*>& objparts = newlib->parts[objlayer];
    vector<node*>& srcparts = newlib->parts[srclayer];
    vector<node*>& origparts = srclib->parts[layer];

    // Add "objlayer" first (without) "backedges"
    //for (vector<node*>::iterator iter = origparts.begin(); iter != origparts.end(); ++iter) {
    //    node* origp = *iter;
    //    part_data* origpd = (part_data*)origp->data;
    //    part_data* newpd = new part_data(*origpd);
    //    node* newp = newlib->add_node(newpd, origp->attr);
    //    
    //    newpd->layer = objlayer;
    //    objparts.push_back(newp);
    //    forall_neighbors(origp, niter) {
    //        int name = neighbor_index(niter);
    //        node* n = neighbor_node(niter);
    //        lib_data* nd = (lib_data*)neighbor_node_data(niter);
    //        edge_data* ed = neighbor_edge_data(niter);

    //        if (name == atom("lyrCenter")) 
    //            newlib->add_edge(newp, srcparts[nd->type], name);
    //        else if (name == atom("lyrSrc")) 
    //            newlib->add_edge_2(newp, srcparts[nd->type], new part_data_2(*((part_data_2*)ed)), name);
    //        else if (name == atom("lyrPrev") || name == atom("lyrSrcM"))
    //            newlib->add_edge_2(newp, srcparts[nd->type], new part_data_2c(*((part_data_2c*)ed)), name);
    //    }
    //}

    // Delete "lyrCenterBack" edges from "srclayer"
    newlib->delete_edges(srcparts.begin(), srcparts.end(), atom("lyrCenterBack"));

    // Add "c2flayer"
    for (trvlist_t::iterator iter = trvlist.begin(); iter != trvlist.end(); ++iter) {
        tr_vector& tr = *iter;
        rpart_data* pd = new rpart_data(c2flayer, (int)c2fparts.size(), 0, 0, 0);
        node* p = newlib->add_node(pd, R_PART_ATTR);

        c2fparts.push_back(p);

        // Add "down" edges 
        for (tr_vector::trviter_t riter = tr.vec.begin(); riter != tr.vec.end(); ++riter) {
            type_rectangle& r = *riter;
            ipoint2 rcenterpt = r.rect.center();
            bool first = true;

            for (set<int>::iterator titer = r.tset.begin(); titer != r.tset.end(); ++titer) {
                if (!first) newlib->add_edge_2(p, srcparts[*titer], new part_data_2c(rcenterpt), atom("lyrSrcM"));
                else {
                    newlib->add_edge_2(p, srcparts[*titer], 
                        new part_data_2r(rcenterpt, (r.rect - rcenterpt).grow(2)), atom("lyrSrc"));
                    first = false;
                }
                newlib->add_edge_2(srcparts[*titer], p, new part_data_2c(rcenterpt), atom("lyrCenterBack"));
            }
            
        }

        // Add "up" edges
        for (set<tr_vector::info_t>::iterator iter = tr.info.begin(); iter != tr.info.end(); ++iter) 
            newlib->add_edge_2(p, objparts[iter->first], new part_data_2c(0, 0), atom("lyrCenterBack"));

    }

    return newlib;    
}

// The map 'result' "counts" how many times i and j, i < j, represented as a pair (i, j) co-occur 
// in the (sub)part "sets" of object parts in the library.
int co_occurrence_map(map<ipoint2, int>& result, part_lib* library, int layer)
{
    vector<node*>& parts = library->parts[layer];
    int srcname = atom("lyrSrc");
    int srcMname = atom("lyrSrcM");
    int max = 0;
    int count = 0;

    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* p = *iter;
        lib_data* pd = (lib_data*)p->data;

        foreach_neighbor(p, srcname, niter) {
            lib_data* pnd = (lib_data*)neighbor_node_data(niter);
            part_data_2* ed = (part_data_2*)neighbor_edge_data(niter);
            set<int> typeset;
            
            /*pt -= box.ll;*/
            library->get_similar_types(pnd->layer, pnd->type, pd->layer, pd->type, ed->x, ed->y, 1.0, typeset);
            typeset.insert(pnd->type);
            for (set<int>::iterator siter1 = typeset.begin(); siter1 != typeset.end(); ++siter1) {
                int i = *siter1;

                for (set<int>::iterator siter2 = typeset.begin(); siter2 != typeset.end(); ++siter2) {
                    int j = *siter2;

                    if (i < j) {
                        map<ipoint2, int>::iterator iter = result.find(ipoint2(i, j));

                        if (iter == result.end()) {
                            result.insert(pair<ipoint2, int>(ipoint2(i, j), 1));
                            if (max == 0) max = 1;
                        } else {
                            ++iter->second;
                            if (iter->second > max) max = iter->second;
                        }
                    }
                }
            }
        }
        if (count++ % 10 == 0) cerr << '.';
    }
    cerr << endl;
    return max;
}

// Insert a layer between regular layer and object layer 'layer'
// in the 'library'. All necessary parameters are given through 
// config_dictionary 'cfg'.
part_lib* insert_c2f_layer(const config_dictionary& cfg, part_lib* library, int layer)
{
    cerr << "Doing insert_c2f_layer" << endl;

    int grow = cfg.get_value_int("grow", 8);
    double subthresh = cfg.get_value_double("subset_threshold", 0.3);
    bool display = cfg.get_value_bool("display_merged_parts", false);

    if (library == nullptr || library->max_layer_index() < layer) return nullptr;

    vector<node*>& parts = library->parts[layer];
    vector<node*>& partsm1 = library->parts[layer - 1];
    int srcname = atom("lyrSrc");
    int centername = atom("lyrCenter");
    trvlist_t origlist;

    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* p = *iter;
        lib_data* pd = (lib_data*)p->data;
        //node* c = p->get_neighbor(centername);
        //lib_data* cd = (lib_data*)c->data;
        //irectangle2 box;

        //     foreach_neighbor(p, srcname, siter) {
        //         part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
        //         
        //          box.eat(ed->x, ed->y);
        //     }

        origlist.push_back(tr_vector(pd->type, ipoint2::zero));
        //origlist.back().vec.push_back(type_rectangle(cd->type, irectangle2(-box.ll /* ! */, -box.ll /* ! */)));
        foreach_neighbor(p, srcname, siter) {
            part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
            lib_data* pnd = (lib_data*)neighbor_node_data(siter);
            ipoint2 pt(ed->x, ed->y);
            set<int> typeset;
            
            //pt -= box.ll; /* ! */
            library->get_similar_types(pnd->layer, pnd->type, pd->layer, pd->type, ed->x, ed->y, 1.0, typeset);
            typeset.insert(pnd->type);
            origlist.back().vec.push_back(type_rectangle(typeset, irectangle2(pt, pt)));
        }
        sort(origlist.back().vec.begin(), origlist.back().vec.end());
    }
    
    trvlist_t tmplist(origlist);
    trvlist_t reslist;
    map<ipoint2, int> cooc;
    int maxcooc;

    cerr << "Making cooc map ";
    maxcooc = co_occurrence_map(cooc, library, layer);

    cerr << "Merging parts ";
    tmplist.sort(type_vector_compare);
    while (!tmplist.empty()) {
        tr_vector tmpv = tmplist.front();

        tmpv.grow(grow);
        for (trvlist_t::iterator iter = tmplist.begin(); iter != tmplist.end();) {
            if (!iter->subvector(tmpv, cooc, maxcooc, (int)(maxcooc * subthresh))) ++iter;
            else {
                iter->insert_into(tmpv);
                iter = tmplist.erase(iter);
            }
        }
        reslist.push_back(tmpv);
    }

    // Display some statistics:
    int minc = INT_MAX, maxc = 0, sum = 0;
    for (trvlist_t::iterator iter = reslist.begin(); iter != reslist.end(); ++iter)  {
        int i = iter->info.size();
        if (i > maxc) maxc = i;
        if (i < minc) minc = i;
        sum += i;
    }
    cerr << endl;
    cerr << "No. of all parts (original): " << origlist.size() << endl;
    cerr << "No. of coarse (merged) parts: " << reslist.size() << endl;
    cerr << "Sizes of coarse parts: max: " << maxc << " min: " << minc << " avg: " << (double)sum/reslist.size() << endl;

    // Display ...
    if (display) {
        for (trvlist_t::iterator iter = reslist.begin(); iter != reslist.end(); ++iter) {
            iter->display_info();
            cerr << endl;
        }
    }

    part_lib* newlib = insert_c2f_layer(library, layer, reslist);

    return newlib;
}

void coarse2fineT1(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing 'coarse2fineT1'" << endl;
    //part_lib* library;

    //read_library(fname, library);

    //if (library == nullptr) return;

    int layer = cfg.get_value_int("layer", 5);
    int matrixdim = cfg.get_value_int("matrix_dim", 200);
    int dmax = cfg.get_value_int("max_distance", 3);
    string dir = cfg.get_value_string("src_dir", "");
    int misses_count = 0;
    int gtruth_count = 0;
    int boxes_count = 0;

    list<string> files;
    //matrix<int> result(matrixdim, matrixdim, 0);

    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, dir))
        list_directory(files, dir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;
        list<irectangle2> gtruths;
        list<layer1_result::box_data_t> boxes;
        vector<pair<double, bool> > confidence; // not used
        vector<bool> hitvector;
        int hits, misses;
        double power;

        cerr << "Reading '" << *file << "'";

        read_layer1_result(res, dir + *file);
        read_groundtruth(gtruths, dir + *file);

        if (res == nullptr) continue;

        cerr << " done. ";

        res->get_boxes(boxes, layer, set<int>(), 0.8);
        res->count_hits(hitvector, misses, confidence, gtruths, boxes, 0.5);
        hits = count(hitvector.begin(), hitvector.end(), true);

        cerr << "Hits: " << hits << "; misses: " << misses << endl;
        boxes_count += (int)boxes.size();
        misses_count += misses;
        gtruth_count += (int)gtruths.size();

        //if (res->max_layer_index() >= layer) {
        //    vector<node*>& s_nodes = res->shape_nodes[layer];
        //    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        //        node* n = *iter;
        //        layer1_data* nd = (layer1_data*)n->data;
        //        vector<node*> nnodes;
        //        
        //        res->get_neighbors_rectangular(nnodes, n, dmax, dmax, 0);
        //        for (vector<node*>::iterator niter = nnodes.begin(); niter != nnodes.end(); ++niter) {
        //            node* nn = *niter;
        //            layer1_data* nnd = (layer1_data*)nn->data;
        //            if (nd->m < (int)result.width && nnd->m < (int)result.height) 
        //                ++result(nd->m, nnd->m);
        //        }
        //    }
        //}
        delete res;
    }

    cerr << "Final result: # of objects: " << gtruth_count << "; # of misses: " << misses_count << 
        "; # of boxes: " << boxes_count << endl;
    cerr << "power (beta): " << (boxes_count == 0 ? 0.0 : (double)(boxes_count - misses_count)/boxes_count) << endl;
    
    //result.save_mathematica(cfg.get_value_string("out_name", "matrix.m").c_str());

    
    //string outname = cfg.get_value_string("out_library", "newlib.plb");
    //int layers = cfg.get_value_int("c2flayers", 1);
    //part_lib* result = nullptr;

    // - Process each test image with layer 'layer' in "rectangle mode" 
    //   (set proj_prod_threshold, proj_add_threshold to 0 in case 'layer' is object layer)
    // - Create matrix N_ij where n_ij = # of occurrences of part i and j where 
    //   distance between part i and j are <= dmax; dmax should be small!
    


    //delete library;
}

void push_back_rpart(list<pair<irectangle2, set<int> > >& result, const irectangle2& r, const set<int>& s)
{
    typedef pair<irectangle2, set<int> > rpart_t;
    
    int rcost = r.area()*(ilog2_x86((unsigned)s.size()) + 1);

    for (list<rpart_t>::iterator iter = result.begin(); iter != result.end(); ++iter) {
        int icost = iter->first.area()*(ilog2_x86((unsigned)iter->second.size()) + 1);
        irectangle2 b = iter->first.bounding_rectangle(r);
        set<int> bset(iter->second);
        bset.insert(s.begin(), s.end());
        int bcost = b.area()*(ilog2_x86((unsigned)bset.size()) + 1);

        if (bcost < rcost + icost) {
            iter->first = b;
            iter->second = bset;
            return;
        }
    }
    result.push_back(rpart_t(r, s));
}

int rectangle_cost(int xdim, int ydim, unsigned ssize)
{
    return xdim * ydim * (ilog2_x86(ssize) + 1);
}

int part_cost(part_lib* library, int lyr, int i, int objgrow)
{
    if (lyr < 0 || i < 0 || library->max_layer_index() < lyr || i >= library->parts[lyr].size())
        return 0;

    node* pc = library->parts[lyr][i];
    int lyrsrc = atom("lyrSrc");
    int result = 0;

    if (pc->is_attr_set(OBJ_PART_ATTR)) {
        part_data* pcd = (part_data*)pc->data;

        foreach_neighbor (pc, lyrsrc, niter) {
            node* p = neighbor_node(niter);
            lib_data* pd = (lib_data*)neighbor_node_data(niter);
            part_data_2* ed = (part_data_2*)neighbor_edge_data(niter);
            set<int> typeset;
            
            library->get_similar_types(pd->layer, pd->type, pcd->layer, pcd->type, ed->x, ed->y, 1.0, typeset);
            typeset.insert(pd->type);

            result += rectangle_cost((int)ed->distr.width - 1 + 2*objgrow, (int)ed->distr.height - 1 + 2*objgrow, typeset.size());
        }
    } else if (pc->is_attr_set(R_PART_ATTR)) {
        rpart_data* pcd = (rpart_data*)pc->data;
        
        foreach_neighbor (pc, lyrsrc, niter) {
            node* p = neighbor_node(niter);
            lib_data* pd = (lib_data*)neighbor_node_data(niter);
            part_data_2r* ed = (part_data_2r*)neighbor_edge_data(niter);
            set<int> typeset;
            
            library->get_similar_types(pd->layer, pd->type, pcd->layer, pcd->type, ed->x, ed->y, 1.0, typeset);
            typeset.insert(pd->type);
                            
            result += rectangle_cost(ed->rect.x_dim(), ed->rect.y_dim(), typeset.size());
        }
    }
    return result;
}

// Merge 'parts' from library 'lib' at layer 'lyr'. The result is returned in partlist;
// 'costs' is a triple (mincost, maxcost, newcost) and 'minsp' is minimal number of subparts in
// parts in 'parts'.
 void merge_parts(list<pair<irectangle2, set<int> > >& partlist, itriple& costs, int& minsp, part_lib* lib, int lyr, const set<int>& parts, int grow)
{
    typedef pair<irectangle2, set<int> > new_part_t;

    //list<new_part_t> partlist;
    int lyrsrc = atom("lyrSrc");
    
    minsp = INT_MAX;
    partlist.clear();

    // Merge 'parts' into new part.
    for (set<int>::const_iterator siter = parts.begin(); siter != parts.end(); ++siter) {
        node* pc = lib->parts[lyr][*siter];

        if (pc->is_attr_set(OBJ_PART_ATTR)) {
            part_data* pcd = (part_data*)pc->data;
            int ncount = 0;

            foreach_neighbor (pc, lyrsrc, niter) {
                node* p = neighbor_node(niter);
                lib_data* pd = (lib_data*)neighbor_node_data(niter);
                part_data_2* ed = (part_data_2*)neighbor_edge_data(niter);
                set<int> typeset;
                irectangle2 rect(-(int)ed->distr.width/2, -(int)ed->distr.height/2, 
                    (int)ed->distr.width/2, (int)ed->distr.height/2);
                
                rect.grow(grow);
                lib->get_similar_types(pd->layer, pd->type, pcd->layer, pcd->type, ed->x, ed->y, 1.0, typeset);
                typeset.insert(pd->type);
                rect += ipoint2(ed->x, ed->y);
                                
                push_back_rpart(partlist, rect, typeset);
                ++ncount;
            }
            if (ncount < minsp) minsp = ncount;
        } else if (pc->is_attr_set(R_PART_ATTR)) {
            rpart_data* pcd = (rpart_data*)pc->data;
            
            foreach_neighbor (pc, lyrsrc, niter) {
                node* p = neighbor_node(niter);
                lib_data* pd = (lib_data*)neighbor_node_data(niter);
                part_data_2r* ed = (part_data_2r*)neighbor_edge_data(niter);
                set<int> typeset;
                
                lib->get_similar_types(pd->layer, pd->type, pcd->layer, pcd->type, ed->x, ed->y, 1.0, typeset);
                typeset.insert(pd->type);
                                
                push_back_rpart(partlist, ed->rect + ipoint2(ed->x, ed->y), typeset);
            }
            if (pcd->minsp < minsp) minsp = pcd->minsp;
        }
    }

    // Check the cost of the new part
    int maxcost = 0, mincost = INT_MAX;
    int newcost = 0;

    for (set<int>::const_iterator siter = parts.begin(); siter != parts.end(); ++siter) {
        int cost = part_cost(lib, lyr, *siter, grow);

        maxcost = max<int>(maxcost, cost);
        mincost = min<int>(mincost, cost);
    }
    for (list<new_part_t>::iterator iter = partlist.begin(); iter != partlist.end(); ++iter) {
        newcost += rectangle_cost(iter->first.x_dim(), iter->first.y_dim(), iter->second.size());
    }
   
    //cout << mincost << "," << newcost << endl;
    costs = itriple(mincost, maxcost, newcost);
}

// Merge 'parts' on layer 'lyr' from 'lib' and insert them into the 'newlib' on layer 'newlyr'
bool merge_parts(part_lib* newlib, part_lib* lib, int newlyr, int lyr, const set<int>& parts, int grow, double factor)
{
    typedef pair<irectangle2, set<int> > new_part_t;

    list<new_part_t> partlist;
    itriple costs;
    int minsp;

    merge_parts(partlist, costs, minsp, lib, lyr, parts, grow); 

    //  mincost = costs.first, maxcost = costs.second, newcost = costs.third
    // 'mincost < factor*newcost'
    if (costs.first < factor*costs.third) 
        return false;

    // Insert part into the 'newlib'.
    node* newp = lib->add_node(new rpart_data(newlyr, newlib->parts[newlyr].size(), 0, 0, minsp), R_PART_ATTR);

    newlib->parts[newlyr].push_back(newp);

    // Down edges
    for (list<new_part_t>::iterator iter = partlist.begin(); iter != partlist.end(); ++iter) {
        int j = 0;
        ipoint2 p = iter->first.center();
        irectangle2 r = iter->first - p;

        for (set<int>::iterator siter = iter->second.begin(); siter != iter->second.end(); ++siter) {
            node* n = newlib->parts[newlyr - 1][*siter];

            if (j == 0) newlib->add_edge_2(newp, n, new part_data_2r(p, r), atom("lyrSrc"));
            else newlib->add_edge_2(newp, n, new part_data_2c(p), atom("lyrSrcM"));
            newlib->add_edge_2(n, newp, new part_data_2c(p), atom("lyrCenterBack"));
            ++j;
        }

    }

    // Up edges
    for (set<int>::const_iterator siter = parts.begin(); siter != parts.end(); ++siter) {
        node* n = newlib->parts[newlyr + 1][*siter];

        newlib->add_edge_2(newp, n, new part_data_2c(0, 0), atom("lyrCenterBack"));
    }
    
    return true;
}

// Merge parts from 'lib' (at layer 'lyr') to 'newlib' (at layer 'newlyr').
// Merges all parts with mincost/newcost >= factor. Merging is performed greedy starting 
// from pair with maximal mincost/newcost ratio.
int merge_parts(part_lib*& newlib, part_lib* lib, int newlyr, int lyr, int grow, double factor)
{
    typedef pair<irectangle2, set<int> > new_part_t;
    typedef pair<double, vector<int> > item_t;


    // Make all pairs
    int layersize = (int)lib->parts[lyr].size();
    list<new_part_t> partlist;
    itriple costs;
    int minsp;
    list<item_t> pairs;

    for (int i = 0; i < layersize; ++i) {
        for (int j = i + 1; j < layersize; ++j) {
            set<int> parts;
            parts.insert(i);
            parts.insert(j);
            
            merge_parts(partlist, costs, minsp, lib, lyr, parts, grow);
            if (costs.first >= factor*costs.third) {
                pairs.push_back(item_t((double)costs.first/costs.third, vector<int>(parts.begin(), parts.end())));
            }
        }
    }

    // Sort pairs and add them to the library!
    pairs.sort(greater<item_t>());

    newlib = (part_lib*)lib->get_copy();
    newlib->insert_empty_layer(newlyr);
    newlib->delete_edges(newlib->parts[newlyr - 1].begin(), newlib->parts[newlyr - 1].end(), atom("lyrCenterBack"));
    vector<bool> coverv(layersize, false);
    int result = 0;

    // Make all pairs.
    for (list<item_t>::iterator iter = pairs.begin(); iter != pairs.end(); ++iter) {
        set<int> toMerge(iter->second.begin(), iter->second.end());
        bool disjoint = true;

        for (set<int>::iterator i = toMerge.begin(); disjoint && i != toMerge.end(); ++i) {
            if (coverv[*i]) disjoint = false;
        }
        if (disjoint && merge_parts(newlib, lib, newlyr, lyr, toMerge, grow, 0.0)) {
            ++result;
            for (set<int>::iterator i = toMerge.begin(); i != toMerge.end(); ++i) 
                coverv[*i] = true;
        }
    }

    // Add the remaining "singleton" parts.
    for (int i = 0; i < layersize; ++i) {
        if (coverv[i]) continue;

        set<int> toMerge;
        toMerge.insert(i);

        merge_parts(newlib, lib, newlyr, lyr, toMerge, grow, 0.0);
        coverv[i] = true;
    }
    
    return result;
}


void merge2partsT(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing merge2partsT - train version" << endl;

    part_lib* library;

    read_library(fname, library);

    if (library == nullptr) return;

    int layer = cfg.get_value_int("layer", 5);
    string outname = cfg.get_value_string("out_library", "newlib.plb");
    int maxpart = cfg.get_value_int("max_part", INT_MAX);
    int grow = cfg.get_value_int("grow", 1);
    double factor = cfg.get_value_double("cost_factor", 0.7);

    part_lib* result = (part_lib*)library->get_copy();
    int layersize = (int)library->parts[layer].size();

    result->insert_empty_layer(layer);
    result->delete_edges(result->parts[layer - 1].begin(), result->parts[layer - 1].end(), atom("lyrCenterBack"));
    maxpart = min<int>(layersize, maxpart);

    // Make all pairs.
    for (int i = 0; i < maxpart; ++i) {
        for (int j = i + 1; j < maxpart; ++j) {
            set<int> toMerge;
            toMerge.insert(i); 
            toMerge.insert(j);

            if (merge_parts(result, library, layer, layer, toMerge, grow, factor)) 
                cerr << "  " << i << ',' << j;
        }
    }

    // Add the remaining "singleton" parts.
    for (int i = maxpart; i < layersize; ++i) {
        set<int> toMerge;
        toMerge.insert(i);

        merge_parts(result, library, layer, layer, toMerge, grow, factor);
    }
    result->save(outname);

    // Display some statistics...
    //int libcost = 0;

    //for (int i = 0; i < layersize; ++i) 
    //    
    //cerr << "Cost of " << part1 << ": " << part_cost(library, layer, part1) << 
    //    "; cost of " << part2 << ": " << part_cost(library, layer, part2) << endl;
    //cerr << "Cost of merged part: " << part_cost(result, layer, 0) << endl;

    delete library;
}

void merge2partsM(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing merge2partsM - manual list of parts" << endl;

    part_lib* library;

    read_library(fname, library);

    if (library == nullptr) return;

    int layer = cfg.get_value_int("layer", 5);
    string outname = cfg.get_value_string("out_library", "newlib.plb");
    int grow = cfg.get_value_int("grow", 1);
    vector<int> pairs;

    cfg.get_value(pairs, "pairs");

    part_lib* result = (part_lib*)library->get_copy();
    int layersize = (int)library->parts[layer].size();

    result->insert_empty_layer(layer);
    result->delete_edges(result->parts[layer - 1].begin(), result->parts[layer - 1].end(), atom("lyrCenterBack"));
    vector<bool> coverv(layersize, false);

    // Make all pairs.
    for (size_t i = 0; i + 1 < pairs.size(); ++i) {
        set<int> toMerge;
        toMerge.insert(pairs[i]); 
        toMerge.insert(pairs[i + 1]);

        if (merge_parts(result, library, layer, layer, toMerge, grow, 0.0)) 
            coverv[pairs[i]] = coverv[pairs[i + 1]] = true;
    }

    // Add the remaining "singleton" parts.
    for (size_t i = 0; i < layersize; ++i) {
        if (coverv[i]) continue;

        set<int> toMerge;
        toMerge.insert(i);

        merge_parts(result, library, layer, layer, toMerge, grow, 0.0);
        coverv[i] = true;
    }
    result->save(outname);

    delete library;
}

void merge2parts(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing merge2parts" << endl;

    part_lib* library;
    part_lib* newlibrary;

    read_library(fname, library);

    if (library == nullptr) return;

    int layer = cfg.get_value_int("src_layer", 5);
    int newlayer = cfg.get_value_int("out_layer", layer);
    string outname = cfg.get_value_string("out_library", "newlib.plb");
    int grow = cfg.get_value_int("grow", 1);
    double factor = cfg.get_value_double("cost_factor", 0.7);
    int partsadded;

    partsadded = merge_parts(newlibrary, library, newlayer, layer, grow, factor);

    cerr << partsadded << " merged parts added." << endl;

    newlibrary->save(outname);
    delete newlibrary;
    delete library;
}

void cross_layer_svm(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing cross_layer_svm" << endl;

    int layer = cfg.get_value_int("src_layer", 6);
    int svmlayer = cfg.get_value_int("src_layer_s", 1);
    string srcdir = cfg.get_value_string("src_dir", "");
    string svmsrcdir = cfg.get_value_string("src_dir_s", "");
    string svmext = cfg.get_value_string("svm_extension", ".ly2");
    string svmfile = cfg.get_value_string("svm_file", "");
    string outdir = cfg.get_value_string("out_dir", "");
    int maxcount = cfg.get_value_int("file_limit", INT_MAX);
    vector<int> svm_angles;
    vector<int> svm_radii;
    int svm_type_count;

    cfg.get_value(svm_angles, "svm_angles", true);
    cfg.get_value(svm_radii, "svm_radii", true);
    cfg.get_value(svm_type_count, "svm_type_count", true);

    list<string> files;
    int count = 0;
    svm_predictor predictor(svmfile);

    if (!predictor.valid()) {
        cerr << "Can not load the svm file '" << svmfile << "'!" << endl;
        return;
    }

    irectangle2 box = predictor.get_average_box();

    box.grow(2);
    box -= box.center();

    end_dir(outdir);
    end_dir(srcdir);
    end_dir(svmsrcdir);

    vector<int> amap, rmap;
    int acount, rcount;
    vector<float> featureMap;

    get_feature_maps(amap, acount, rmap, rcount, svm_angles, svm_radii);

    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, srcdir))
        list_directory(files, srcdir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        if (++count > maxcount) break;

        layer1_result* res;
        layer1_result* svmres;

        cerr << "Processing " << *file;
        read_layer1_result(res, srcdir + *file);
        read_layer1_result(svmres, change_extension(srcdir + *file, svmext));
        if (res == nullptr || svmres == nullptr) { 
            cerr << " error reading file(s)!";
        } else {
            vector<node*>& s_nodes = res->shape_nodes[layer];
            double size_ratio = (double)res->x_size(svmlayer)/res->x_size(layer);
            ofstream os(change_extension(outdir + *file, ".m").c_str());

            if (os.fail()) {
                cerr << "Can not open output file '" << (outdir + *file) << "'!" << endl;
            } else {
                bool writecomma = false;

                os << '{';
                for (vector<node*>::iterator niter = s_nodes.begin(); niter != s_nodes.end(); ++niter) {
                    node* n = *niter;

                    while (n != nullptr) {
                        layer1_data* nd = (layer1_data*)n->data;
                        ipoint2 p(nd->x, nd->y);
                        vector<node*> result;
                        irectangle2 rbox;
                        float prediction;

                        p *= size_ratio;
                        rbox = box*size_ratio;

                        svmres->get_nodes_in_rect(result, svmlayer, rbox + p);
                        get_feature_map_radial(featureMap, result, p, amap, acount, rmap, rcount, svm_type_count);
                        prediction = predictor.predict(featureMap);

                        if ((int)prediction == 0) {
                            if (writecomma) os << ','; else writecomma = true;
                            os << '{' << nd->x << ',' << nd->y << ',' << nd->m << '}';
                        }

                        n = nd->next;
                    }
                }
                os << '}' << endl;
                os.close();
            }

            cerr << " done";
        }
        if (res != nullptr) delete res;
        if (svmres != nullptr) delete svmres;
        cerr << endl;
    }
    
}

void svm_test(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing svm_test" << endl;

    int layer = cfg.get_value_int("src_layer", 1);
    //string srcdir = cfg.get_value_string("src_dir", "");
    string srcdir = "D:\\work\\data\\multiclass\\svmtest\\";

    map<int, double> trainMap;
    list<string> pfiles;
    list<string> nfiles;
    part_lib* library;

    pfiles.push_back("train_mugs-eth_image1_0.ly5");
    pfiles.push_back("train_mugs-eth_image2_0.ly5");
    pfiles.push_back("train_mugs-eth_image3_0.ly5");
    pfiles.push_back("train_mugs-eth_image4_0.ly5");
    //pfiles.push_back("train_mugs-eth_image5_0.ly5");

    nfiles.push_back("train_pliers-eth_image1_0.ly5");
    nfiles.push_back("train_pliers-eth_image2_0.ly5");
    nfiles.push_back("train_pliers-eth_image3_0.ly5");
    nfiles.push_back("train_pliers-eth_image4_0.ly5");
    //nfiles.push_back("train_pliers-eth_image5_0.ly5");

    read_library(srcdir + "lib5.plb", library);

    vector<int> radii;
    vector<int> angles;
    int type_count = 30; /// <----- !!!!!!!
    int feature_count = ((int)radii.size() + 1)*((int)angles.size() + 1)*type_count;

    svm svmTrain(feature_count);

    cerr << "Loading positive examples...";
    for (list<string>::iterator iter = pfiles.begin(); iter != pfiles.end(); ++iter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *iter);
        if (res == nullptr) continue;

        map<int, double> featureMap;
        vector<node*>& nodes = res->shape_nodes[layer];
        ipoint2 center = node_set_center(nodes.begin(), nodes.end());
        
        get_feature_map_radial(featureMap, nodes, center, angles, radii, type_count); 
        svmTrain.add_positive_sample(featureMap);

        delete res;
    }

    cerr << endl << "Loading negative examples...";
    for (list<string>::iterator iter = nfiles.begin(); iter != nfiles.end(); ++iter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *iter);
        if (res == nullptr) continue;

        map<int, double> featureMap;
        vector<node*>& nodes = res->shape_nodes[layer];
        ipoint2 center = node_set_center(nodes.begin(), nodes.end());
        
        get_feature_map_radial(featureMap, nodes, center, angles, radii, type_count); 
        svmTrain.add_negative_sample(featureMap);

        delete res;
    }

    cerr << endl << "Training SVM...";
    svmTrain.train("c:\\work\\svm.xml", irectangle2());

    // Do test
    ///////////////////////////////////////////////////////////////////////////

    list<string> tfiles;

    tfiles.push_back("train_mugs-eth_image4_0.ly5");
    tfiles.push_back("train_mugs-eth_image5_0.ly5");
    tfiles.push_back("train_pliers-eth_image5_0.ly5");

    svm_predictor svmPredict("c:\\work\\svm.xml");

    cerr << "Classifying examples..." << endl;
    for (list<string>::iterator iter = tfiles.begin(); iter != tfiles.end(); ++iter) {
        layer1_result* res;
        double pres;

        read_layer1_result(res, srcdir + *iter);
        if (res == nullptr) continue;

        map<int, double> featureMap;
        vector<node*>& nodes = res->shape_nodes[layer];
        ipoint2 center = node_set_center(nodes.begin(), nodes.end());
        
        get_feature_map_radial(featureMap, nodes, center, angles, radii, type_count); 
        pres = svmPredict.predict(featureMap);

        cerr << *iter << ": " << pres << endl;

        delete res;
    }

}


// Calculate the "power" of the merged parts on 'layer'.
// "Algorithm": for each detection n of type t on 'layer' which is 
//   "covered" with some detection of type u on 'layer' + 1 increase mmap(t, u) by 1.
//   mmap(t, -1) counts the number of all occurrences of t.
void powercalc(map<iipair, int>& mmap, layer1_result* res, int layer, double factor)
{
    if (res == nullptr || layer < 0 || res->max_layer_index() < layer + 1) return;

    res->init_grid(layer + 1);
    
    vector<node*>& s_nodes = res->shape_nodes[layer];
    vector<node*>& s_nodes1 = res->shape_nodes[layer + 1];
    int toprev = atom("toPrevLayer");
    int toprevI = atom("toPrevLayerI");

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        
        while (n != nullptr) {
            layer1_data* nd = (layer1_data*)n->data;

            map<iipair, int>::iterator miter = mmap.find(iipair(nd->m, -1));

            if (miter != mmap.end()) ++miter->second;
            else mmap.insert(pair<iipair, int>(iipair(nd->m, -1), 1));

            node* supn = res->node_at(nd->x, nd->y, nd->z + 1);
            while (supn != nullptr) {
                layer1_data* supnd = (layer1_data*)supn->data;

                if (supn->is_neighbor(n, toprevI)) {
                    map<iipair, int>::iterator miter = mmap.find(iipair(nd->m, supnd->m));

                    if (miter != mmap.end()) ++miter->second;
                    else mmap.insert(pair<iipair, int>(iipair(nd->m, supnd->m), 1));
                }
                supn = supnd->next;
            }
            n = nd->next;
        }
    }
}

struct powercalc_info {
    //powercalc_info() : cost(0), structure(), classes() { }
    //powercalc_info(int c, const vector<iipair>& s, const set<int>& cs) : cost(c), structure(s), classes(cs) { }
    int cost;
    vector<iipair> structure;    // structure w.r. to the sup-layer; (index, power)
    set<int> classes;            // a set of classes containing the part
};


void powercalc(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing powercalc" << endl;

    string dir = cfg.get_value_string("src_dir", "");
    int layer = cfg.get_value_int("layer", 5);
    int count = 0;
    int maxcount = cfg.get_value_int("file_limit", INT_MAX);



    // Calculate "powers" of parts on a set of processed images
    list<string> files;
    map<iipair, int> power;    // maps ("part index", "covered index") to # of occurrences

    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, dir))
        list_directory(files, dir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        if (++count > maxcount) break;

        layer1_result* res;

        cerr << "Processing " << *file;
        read_layer1_result(res, dir + *file);
        if (res != nullptr) {
            powercalc(power, res, layer, 1.0);
            delete res;
        }
        cerr << " done" << endl;
    }

    // Calculate "costs"
    string libname = cfg.get_value_string("library", "");

    part_lib* library;
    map<int, powercalc_info> info;
    int backname = atom("lyrCenterBack");
    
    read_library(libname, library);
    if (library != nullptr) {
        vector<node*>& parts = library->parts[layer];

        for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
            node* p = *iter;
            lib_data* pd = (lib_data*)p->data;
            powercalc_info pinfo;

            // Calc cost
            pinfo.cost = part_cost(library, pd->layer, pd->type, 1);

            // Calc sup-layer structure
            foreach_neighbor (p, backname, niter) {
                lib_data* pnd = (lib_data*)neighbor_node_data(niter);
                int nc = part_cost(library, pnd->layer, pnd->type, 1);

                pinfo.structure.push_back(iipair(pnd->type, nc));
            }

            // Calc class structure
            set<node*> cnodes;

            library->recurse_from_node(p, backname, cnodes);
            for (set<node*>::iterator siter = cnodes.begin(); siter != cnodes.end(); ++siter)
                pinfo.classes.insert(((lib_data*)(*siter)->data)->type);
            info.insert(pair<int, powercalc_info>(pd->type, pinfo));
        }
    }

    // Output

    cerr << "Result:" << endl;
    for (map<iipair, int>::iterator iter = power.begin(); iter != power.end(); ++iter) {
        cout << iter->first.first << ',' << iter->first.second << ',' << iter->second;
        
        if (iter->first.second < 0) {
            map<int, powercalc_info>::iterator iiter = info.find(iter->first.first);

            if (iiter != info.end()) {
                cout << ',' << iiter->second.cost;

                cout << ',' << iiter->second.structure.size();
                for (size_t i = 0; i < iiter->second.structure.size(); ++i) 
                    cout << ',' << iiter->second.structure[i].first << ',' << iiter->second.structure[i].second;
                cout << ',' << iiter->second.classes.size();
                for (set<int>::iterator i = iiter->second.classes.begin(); i != iiter->second.classes.end(); ++i)
                    cout << ',' << *i;
            }
        }
        cout << endl;
    }
}

void coarse2fine(const config_dictionary& cfg, const char* fname)
{
    part_lib* library;

    read_library(fname, library);

    if (library == nullptr) return;

    int layer = cfg.get_value_int("layer", 5);
    string outname = cfg.get_value_string("out_library", "newlib.plb");
    int layers = cfg.get_value_int("c2flayers", 1);
    part_lib* result = nullptr;

    for (int l = 0; l < layers; ++l) {
        string c2fcfgkey = string("c2fcfg") + (l + 1);
        string c2fcfg;
        
        cfg.get_value(c2fcfg, c2fcfgkey, true);
        
        config_dictionary cf(c2fcfg);
        part_lib* newlib = insert_c2f_layer(cf, library, layer);
        
        if (result == nullptr) result = newlib;
        else {
            part_lib* tmplib = result->c2f_merge(newlib, layer + l - 1, layer);

            delete result;
            result = tmplib;
        }
    }

    if (result != nullptr) {
        result->save(outname);
        delete result;
    }
    delete library;    
}

void insert_rlayers(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing linsert" << " ******************OBSOLETE************** " << endl;

    part_lib* library;
    int layer = cfg.get_value_int("layer", 5);
    int rect0dim = cfg.get_value_int("rect_dim", 4);
    string outname = cfg.get_value_string("out_library", "newlib.plb");
    //int partition_classes = cfg.get_value_int("partition_classes", 2);

    read_library(fname, library);
    if (library == nullptr || library->max_layer_index() < layer) return;

    vector<node*>& parts = library->parts[layer];
    vector<node*>& partsm1 = library->parts[layer - 1];
    int srcname = atom("lyrSrc");
    int centername = atom("lyrCenter");
	list<layer1_result*> rlist;
    vector<ipoint2> clist;
    int partition_classes = partsm1.size();
    vector<int> partition(partsm1.size());


    // make a simple partition by "mod" operation
    for (size_t i = 0; i < partsm1.size(); ++i) 
        partition[i] = i;
    //partition_classes = object_layer_partition(partition, library, layer);

    // Make "layer1_result" structure from each object part composition
    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {

        node* p = *iter;
		layer1_result* r = new layer1_result();
		irectangle2 box;
		node* n;

		// Get the bounding box of the subparts
		box.eat(0, 0);
        foreach_neighbor(p, srcname, siter) {
            part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
            
			box.eat(ed->x/*/rect0dim + ((ed->x < 0) ? -1 : 0)*/, 
                ed->y/*/rect0dim + ((ed->y < 0) ? -1 : 0)*/);
        }

		// normalize the size of box s.t. (0, 0) is on even 
		//if (box.ll.x < 0) --box.ll.x;
		//if (box.ll.y < 0) --box.ll.y;
		//if (box.ll.x % 2 != 0) --box.ll.x;
		//if (box.ll.y % 2 != 0) --box.ll.y;

		// add nodes to r; shape_nodes = shape_nodes_inhib are filled
		// s.t. the first element is the "center"
        node* c = p->get_neighbor(centername);
        lib_data* cd = (lib_data*)c->data;

        r->new_grid(box.x_dim() + 1, box.y_dim() + 1, 0);
        r->shape_nodes.push_back(vector<node*>());
        r->shape_nodes_inhib.push_back(vector<node*>());
        n = r->add_grid_node(new layer1_data(1.0, partition[cd->type]), -box.ll.x, -box.ll.y, 0);
        clist.push_back(-box.ll);
		r->shape_nodes[0].push_back(n);
		r->shape_nodes_inhib[0].push_back(n);
        foreach_neighbor(p, srcname, siter) {
            part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
            lib_data* pnd = (lib_data*)neighbor_node_data(siter);
            int x = ed->x/*/rect0dim + ((ed->x < 0) ? -1 : 0)*/ - box.ll.x;
            int y = ed->y/*/rect0dim + ((ed->y < 0) ? -1 : 0)*/ - box.ll.y;
            
            if (r->node_at(x, y, 0) == nullptr) {
                n = r->add_grid_node(new layer1_data(1.0, partition[pnd->type]), x, y, 0);
			    r->shape_nodes[0].push_back(n);
			    r->shape_nodes_inhib[0].push_back(n);
            }
        }

		rlist.push_back(r);
    }
    

    cerr << "Starting with " << partition_classes << " parts" << endl;
    // Make new "rectangle" library 
    part_lib* rlib = new part_lib(0);

    // add layer 0 to the library ---------------------------------------------
    rlib->layer_count = 0;
    for (int i = 0; i < partition_classes; ++i) 
        rlib->parts[0].push_back(rlib->add_node(new rpart_data(0, i, 0, 0, 0)));

    // add layer l = 1 to the library -----------------------------------------
    learn_layer(rlib, rlist, 1);
    add_layer(rlist.begin(), rlist.end(), rlib, 1);


    // add layer l = 2 to the library -----------------------------------------
    learn_layer(rlib, rlist, 2);
    add_layer(rlist.begin(), rlist.end(), rlib, 2);


    // o2n[i] is the index of the ith part in rlist
    vector<int> o2n;
    vector<ipoint2>::iterator cliter = clist.begin();

    for (list<layer1_result*>::iterator riter = rlist.begin(); riter != rlist.end(); ++riter, ++cliter) {
        layer1_data* nd = (layer1_data*)((*riter)->shape_nodes[2].front()->data);
        o2n.push_back(nd->m);
        *cliter -= ipoint2(nd->x, nd->y);
    }

    // "Transplant" the rectangle library to the new, final, library
    part_lib* newlib = (part_lib*)library->get_copy();
    newlib->delete_parts_geq(layer);
    newlib->layer_count = layer;
    combine_libraries(newlib, rlib, library, partition, o2n, clist);

    // Save images -- optionally
    if (cfg.get_value_bool("save_layers", false)) {
        for (int l = 0; l < rlib->layer_count; ++l) {
            cerr << "Saving layer " << l << " with " << rlib->parts[l].size() << " parts." << endl;
            string name = string("c:\\temp\\") + l + string(".png");
            rlib->save_all(name.c_str(), l + 1, 
                cfg.get_value_int("saving_start", 0),
                cfg.get_value_int("saving_end", -1), false, true);
        }
    }

    // Save
    newlib->save(outname);

    // Cleaning
	for (list<layer1_result*>::iterator riter = rlist.begin(); riter != rlist.end(); ++riter) {
		delete *riter;
	}
    delete library;
    delete newlib;
    delete rlib;
}

void layer_statistics(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing layer_statistics" << endl;

    string dir = cfg.get_value_string("src_dir", "");
    list<string> files;
    layer1_result* res;
 
    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, dir))
        list_directory(files, dir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        read_layer1_result(res, dir + *file);
        
        cerr << "Statistics for: " << *file << endl;
        for (int i = 0; i <= res->max_layer_index(); ++i) {
            int count = 0;
            int minat = INT_MAX;
            int maxat = 0;
            int avgat = 0;
            double avgval = 0;
            double maxval = 0.0;
            double minval = 1.0;

            for (vector<node*>::iterator iter = res->shape_nodes[i].begin(); iter != res->shape_nodes[i].end(); ++iter) {
                node* n = *iter;
                layer1_data* nd = (layer1_data*)n->data;
                vector<node*> nodesat;

                res->nodes_at(nodesat, n);

                //for (vector<node*>::iterator atiter = nodesat.begin(); atiter != nodesat.end(); ++atiter) {
                //    layer1_data* atnd = (layer1_data*)(*atiter)->data;
                    double v = nd->vval();

                    if (v < minval) minval = v;
                    if (v > maxval) maxval = v;
                    avgval += v;
                //}

                int atsize = (int)nodesat.size();
                if (atsize > maxat) maxat = atsize; 
                if (atsize < minat) minat = atsize;
                count += atsize;
            }
            if (res->shape_nodes[i].size() > 0) {
                avgat = count/(int)res->shape_nodes[i].size();
                avgval /= res->shape_nodes[i].size();
            }
            cerr << "Layer " << i << ": " << count << " parts" << endl;
            cerr << "  minat: " << minat << ", maxat: " << maxat << ", avgat: " << avgat;
            cerr << "  minval: " << minval << ", maxval: " << maxval << ", avgval " << avgval << endl;

        }

        if (res != nullptr) {
            delete res;
        }
    }

}

struct nodes_at_sort_f {
    int nbname;

    nodes_at_sort_f() : nbname(atom("toLayer0")) { }
    bool operator()(node* n, node* m) 
    {
        layer1_data* nd = (layer1_data*)n->data;
        layer1_data* md = (layer1_data*)m->data;

        return (0.9*nd->r(R_RESPONSE) + 0.1*nd->r(G_RESPONSE)) > (0.9*md->r(R_RESPONSE) + 0.1*md->r(G_RESPONSE));

        //int ncount = n->count_neighbors(nbname);
        //int mcount = m->count_neighbors(nbname);
        //
        //if (ncount > mcount) return true;
        //else if (ncount < mcount) return false;
        //else {
        //    layer1_data* nd = (layer1_data*)n->data;
        //    layer1_data* md = (layer1_data*)m->data;

        //    return (nd->r(R_RESPONSE) + nd->r(G_RESPONSE))/2.0 > (md->r(R_RESPONSE) + md->r(G_RESPONSE))/2.0;
        //}
        //layer1_data* nd = (layer1_data*)n->data;
        //layer1_data* md = (layer1_data*)m->data;

        //return (nd->r(R_RESPONSE) + nd->r(G_RESPONSE))/2.0 > (md->r(R_RESPONSE) + md->r(G_RESPONSE))/2.0;
        //return (nd->r(R_RESPONSE))/1.0 > (md->r(R_RESPONSE))/1.0;
    }
};

template<class T, class U> pair<T, U> sort(const pair<T, U>& p)
{
    if (p.first <= p.second) return p; else return pair<T, U>(p.second, p.first);
}

void select_nodes_by_type_x(vector<node*>& result, layer1_result* res, const vector<node*>& nodes, 
    const set<int>& types)
{
    for (vector<node*>::const_iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        
        do {
            layer1_data* nd = (layer1_data*)n->data;

            if (types.find(nd->m) != types.end()) {
                result.push_back(n);
                break;
            }
            n = nd->next;
        } while (n != nullptr);
    }
}

void best_response_distribution(const config_dictionary& cfg, const char* fname)
{
    streamable_pointer::set_max_mem_objects(30);

    string srcdir = cfg.get_value_string("src_dir", "");
    string libname;
    int layer;
    list<string> files;
    part_lib* library;

    cfg.get_value(layer, "layer", true);
    cfg.get_value(libname, "library", true);
    read_library(libname, library);

    if (library == nullptr) {
        cerr << "Library does not exist." << endl;
        return;
    }
    if (library->max_layer_index() > layer || layer < 0) {
        cerr << "Layer " << layer << " does not exist in the library!" << endl;
        delete library;
        return;
    }

    list<layer1_result_ptr> workingset;

    end_dir(srcdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, srcdir))
        list_directory(files, srcdir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        cerr << "Processing: " << *file;
        read_layer1_result(res, srcdir + *file);
        if (res == nullptr) 
            cerr << " error reading file" << endl;
        else {
            res->add_reconstruction_edges(layer);
            res->init_grid(layer);
            workingset.push_back(layer1_result_ptr(res));
            cerr << " done" << endl;
        }
    }

    nodes_at_sort_f sortf;
    double min_uncovered = numeric_limits<double>::max();
    double covered_thresh = cfg.get_value_double("covered_threshold", 0.8);
    double start_factor = cfg.get_value_double("start_factor", 0.2);
    double factor_power = cfg.get_value_double("factor_power", 2.0);
    double factor = start_factor;

    set<int> current_set;
    vector<int> selected;

    do {
        vector<int> statistics(library->layer_size(layer), 0);

        cerr << endl;
        cerr << "Gathering statistics" << endl << "--------------------" << endl;

        for (list<layer1_result_ptr>::iterator wsiter = workingset.begin(); wsiter != workingset.end(); ++wsiter) {
            layer1_result* res = *wsiter;
            vector<node*>& s_nodes = res->shape_nodes_inhib[layer];

            // Calc covering ratio
            vector<node*> visited;
            set<node*> covered;
            set<node*> maxcovered;

            select_nodes_by_type_x(visited, res, s_nodes, current_set);
            res->recurse(visited, atom("toLayer0"), covered);
            res->recurse(s_nodes, atom("toLayer0"), maxcovered);

            double covered_ratio = (double)covered.size()/(double)maxcovered.size();

            cerr << "CR = " << covered_ratio << endl;

            if (covered_ratio >= covered_thresh) continue;

            // Update statistics
            for (vector<node*>::iterator sniter = s_nodes.begin(); sniter != s_nodes.end(); ++sniter) {
                node* n = *sniter;
                vector<node*> nvec;
                bool found = false;

                res->nodes_at(nvec, n);

                for (vector<node*>::iterator nviter = nvec.begin(); nviter != nvec.end(); ++nviter) {
                    layer1_data* nnd = (layer1_data*)(*nviter)->data;
                    if (current_set.find(nnd->m) != current_set.end()) {
                        found = true;
                        break;
                    }
                }
                if (found) continue;

                sort(nvec.begin(), nvec.end(), sortf);
                
                node* bestn = nvec.front();
                layer1_data* bestnd = (layer1_data*)bestn->data;

                ++statistics[bestnd->m];
            }
        }

        vector<int>::iterator riter;
        vector<int> tmpstat = statistics;

        cerr << endl;
        cerr << "Statistics:" << endl << "-----------" << endl;

        selected.clear();
        while (*(riter = max_element(tmpstat.begin(), tmpstat.end())) > 0) {
            int pos = riter - tmpstat.begin();

            cerr << pos << ": #" << *riter << endl;

            tmpstat[pos] = 0;
            selected.push_back(pos);
        }

        cerr << endl;
        cerr << "Selection:" << endl << "-----------" << endl;

        for (int i = 0; i < (int)selected.size() && statistics[selected[i]] >= factor*statistics[selected[0]]; ++i) {
            current_set.insert(selected[i]);  

            cerr << selected[i] << ' ';
        }

        cerr << endl;

        factor = pow(factor, factor_power);

    } while (!selected.empty());

    vector<int> v(current_set.begin(), current_set.end());
    if (library != nullptr) {
        library->save_all("outlib.png", layer + 1, v, cfg.get_value_bool("show_labels", false));
        delete library;
    }

    if (cfg.is_defined("out_library")) {
        library->keep_parts(layer, v);
        library->save(cfg.get_value_string("out_library", ""));
    }
}

/*
void best_response_distribution(const config_dictionary& cfg, const char* fname)
{
    typedef vector<int> result_t;
    typedef map<iipair, int> result2_t;

    string srcdir = cfg.get_value_string("src_dir", "");
    string libname;
    int layer;
    list<string> files;
    part_lib* library;

    cfg.get_value(layer, "layer", true);
    cfg.get_value(libname, "library", true);
    read_library(libname, library);

    if (library == nullptr) {
        cerr << "Library does not exist." << endl;
        return;
    }
    if (library->max_layer_index() > layer || layer < 0) {
        cerr << "Layer " << layer << " does not exist in the library!" << endl;
        delete library;
        return;
    }

    nodes_at_sort_f sortf;
    result_t result(library->layer_size(layer), 0);
    result2_t result2;

    end_dir(srcdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, srcdir))
        list_directory(files, srcdir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        cerr << "Processing: " << *file;
        read_layer1_result(res, srcdir + *file);
        if (res == nullptr) 
            cerr << " error reading file" << endl;
        else {
            cerr << endl;
            //ddpair dtop = res->top_response_distribution(layer, R_RESPONSE);
            //ddpair dall = res->response_distribution(layer, R_RESPONSE);
            //cerr << "    \"top\" r-distribution: m = " << dtop.first << ", var = " << dtop.second << "  " << endl;
            //cerr << "    r-distribution: m = " << dall.first << ", var = " << dall.second << "  " << endl;
            //dtop = res->top_response_distribution(layer, G_RESPONSE);
            //dall = res->response_distribution(layer, G_RESPONSE);
            //cerr << "    \"top\" g-distribution: m = " << dtop.first << ", var = " << dtop.second << "  " << endl;
            //cerr << "    g-distribution: m = " << dall.first << ", var = " << dall.second << "  " << endl;
            res->add_reconstruction_edges(layer);
            res->init_grid(layer);

            for (int i = 0; i < res->shape_nodes_inhib[layer].size(); ++i) {
                vector<node*> v;
                node* n = res->shape_nodes_inhib[layer][i];
                layer1_data* nd = (layer1_data*)n->data;

                //res->sorted_nodes_at(v, nd->x, nd->y, nd->z, R_RESPONSE);
                //cerr << "Best r-response: " << ((layer1_data*)(v[0]->data))->r(R_RESPONSE) << endl;
                //res->sorted_nodes_at(v, nd->x, nd->y, nd->z, G_RESPONSE);
                //cerr << "Best g-response: " << ((layer1_data*)(v[0]->data))->r(G_RESPONSE) << endl;


                vector<node*> nvec;

                res->nodes_at(nvec, nd->x, nd->y, nd->z);
                sort(nvec.begin(), nvec.end(), sortf);
                
                node* bestn = nvec.front();
                layer1_data* bestnd = (layer1_data*)bestn->data;
                set<node*> bestset;

                ++result[bestnd->m];

                bestn->get_neighbor_set(atom("toLayer0"), bestset);
                for (vector<node*>::iterator niter = ++nvec.begin(); niter != nvec.end(); ++niter) {
                    node* n = *niter;
                    layer1_data* nd = (layer1_data*)n->data;
                    set<node*> nset;

                    n->get_neighbor_set(atom("toLayer0"), nset);
                    if (is_subset(nset.begin(), nset.end(), bestset.begin(), bestset.end())) {
                        double rdif = bestnd->r(R_RESPONSE) - nd->r(R_RESPONSE);
                        double gdif = bestnd->r(G_RESPONSE) - nd->r(G_RESPONSE);

                        //cerr << "m = " << nd->m << " issubset of bestm = " << bestnd->m << 
                        //    "  r-dif = " << rdif << "  g-dif = " << gdif << endl;
                        if (rdif < 0.1 && gdif < 0.1) {
                            result2_t::iterator r2iter = 
                                result2.insert(result2_t::value_type(sort(iipair(nd->m, bestnd->m)), 0)).first;
                            ++r2iter->second;
                            ++result[nd->m];
                        }
                    }
                }
                
                //for (vector<node*>::iterator iter = nvec.begin(); iter != nvec.end(); ++iter) {
                //    layer1_data* nnd = (layer1_data*)(*iter)->data;

                //    cerr << "m = " << nnd->m << ", covering: " << (*iter)->count_neighbors(atom("toLayer0")) << 
                //        ", r-response: " << nnd->r(R_RESPONSE) << ", g-response: " <<
                //        nnd->r(G_RESPONSE) << endl;
                //}
                //cerr << endl;
                //int dummy;
                //cin > > dummy;
            }

            delete res;
            cerr << " done" << endl;
        }
    }

    vector<int> final;
    result_t::iterator riter;

    cerr << "Result:" << endl << "---------------" << endl;
    while (*(riter = max_element(result.begin(), result.end())) > 0) {
        int pos = riter - result.begin();

        cerr << pos << ": #" << *riter << endl;

        for (result2_t::iterator r2iter = result2.begin(); r2iter != result2.end(); ++r2iter) {
            if (r2iter->first.first == pos) { 
                result[r2iter->first.second] -= r2iter->second;
                r2iter->second = 0;
            } else if (r2iter->first.second == pos) {
                result[r2iter->first.first] -= r2iter->second;
                r2iter->second = 0;
            }
        }
        result[pos] = 0;
        final.push_back(pos);
    }

    //vector<iipair> result3;
    //for (result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
    //    result3.push_back(iipair(iter->second, iter->first));
    //}
    //sort(result3.begin(), result3.end(), greater<iipair>());

    //cerr << endl;
    //cerr << "Result:" << endl;
    //cerr << "--------------------------------" << endl;
    //for (result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
    //    cerr << iter->first << ": " << iter->second << endl;
    //}

    //vector<int> order;

    //for (vector<iipair>::iterator iter = result2.begin(); iter != result2.end(); ++iter) {
    //    cerr << iter->second << ": " << iter->first << endl;
    //    order.push_back(iter->second);
    //}

    if (library != nullptr) {
        library->save_all("outlib.png", layer + 1, final, cfg.get_value_bool("show_labels", false));
        delete library;
    }
}*/


void compare_layer1_results(const config_dictionary& cfg, const char* fname)
{
    string file1, file2;
    layer1_result* res1, * res2;

    cfg.get_value(file1, "file1", true);
    cfg.get_value(file2, "file2", true);
    
    cerr << "comparing for \"" << file1 << " <= " << file2 << "\"..." << endl;
    
    read_layer1_result(res1, file1);
    if (res1 == nullptr) { 
        cerr << "Error loading: " << file1 << endl; 
        return;
    }
    read_layer1_result(res2, file2);
    if (res2 == nullptr) {
        cerr << "Error loading: " << file2 << endl;
        delete res1;
        return;
    }

    int maxlayer = std::min<int>(res1->max_layer_index(), res2->max_layer_index());

    for (int i = 0; i <= maxlayer; ++i) {
        cerr << "Checking layer " << i << endl;

        if (!res1->grid(i)) res1->init_grid(i);
        if (!res2->grid(i)) res2->init_grid(i);
        for (list<node*>::iterator iter = res1->nodes.begin(); iter != res1->nodes.end(); ++iter) {
            node* n = *iter;
            layer1_data* nd = (layer1_data*)n->data;
            set<node*> set1, set2;

            if (!n->is_attr_set(IMG_NODE_ATTR) || nd->z != i) continue;
            res1->nodes_at(set1, nd->x, nd->y, nd->z);
            res2->nodes_at(set2, nd->x, nd->y, nd->z);

            for (set<node*>::iterator s1iter = set1.begin(); s1iter != set1.end(); ++s1iter) {
                node* n1 = *s1iter;
                layer1_data* n1d = (layer1_data*)n1->data;
                bool found  = false;

                for (set<node*>::iterator s2iter = set2.begin(); !found && s2iter != set2.end(); ++s2iter) {
                    node* n2 = *s2iter;
                    layer1_data* n2d = (layer1_data*)n2->data;

                    found = n2d->m == n1d->m && n1d->val() <= n2d->val();
                }
                if (!found) 
                    cerr << "Type " << n1d->m << " at (" << n1d->x << ',' << n1d->y << ',' << n1d->z << ')' << 
                        " does not exist in file2 at the same position!" << endl;
            }
        }
    }
    
    delete res1;
    delete res2;
}


void library_layer_statistics(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing layer_statistics" << endl;

    part_lib* library;

    read_library(fname, library);
    if (library == nullptr) {
        cerr << "Can not open library '" << fname << "'" << endl;
        return;
    }

    for (int i = 0; i <= library->max_layer_index(); ++i) {
        vector<node*>& parts = library->parts[i];

        cout << "Layer: " << i << endl;
        for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
            node* p = *iter;
            lib_data* pd = (lib_data*)p->data;

            cout << "  Part " << pd->type << ": R_THRESH = " << pd->get_thresh(R_THRESH, -1.0) <<
                ",  G_THRESH = " << pd->get_thresh(G_THRESH, -1.0) << ", RR_THRESH = " << pd->get_thresh(RR_THRESH, -1.0) << endl;
        }
        cout << endl;
    }

    delete library;
}

void back_statistics(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing library_statistics" << endl;

    part_lib* library;
    int backname = atom("lyrCenterBack");
    int objectlayer = INT_MAX - 1;

    read_library(fname, library);
    if (library == nullptr) {
        cerr << "Can not open library '" << fname << "'" << endl;
        return;
    }

	// Basic library statistics
	//
	cout << "Contractions vector: ";
	for (int i = 0; i < library->layer_count + 1; ++i) {
		if (i > 0) cout << ", ";
		cout << library->contractions[i];
	}
	cout << endl;

    // Back statistics 
    //
    for (int i = 0; i < library->layer_count; ++i) {
        vector<node*>& parts = library->parts[i];
        int max = 0, min = INT_MAX;
        double avg = 0.0;
        
        for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
            node* n = *iter;
            int count = 0;

            foreach_neighbor(n, backname, niter) {
                ++count;
            }
            if (count > max) max = count;
            if (count < min) min = count;
            avg += count;
        }
        avg /= parts.size();
        if (parts.size() > 0 && parts[0]->is_attr_set(OBJ_PART_ATTR)) 
            objectlayer = i;

        cout << "Layer " << i << ": size = " << parts.size() << "; min up = " << min;
        cout << "; max up = " << max << "; avg up = " << avg << endl;

    }

    // Category layer statistics
    //
    if (objectlayer + 1 <= library->max_layer_index()) {
        vector<node*>& parts = library->parts[objectlayer + 1];
        int prevname = atom("lyrPrev");
        int count = 0;

        cout << "Layer " << objectlayer + 1 << " is class layer: number of classes = " << parts.size() << endl;
        for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
            node* p = *iter;
            //set<vector<ipoint2> > posset;

            //foreach_neighbor(p, prevname, oiter) {
            //    vector<ipoint2> pos;
            //    irectangle2 box;

            //    box = library->subpart_positions(pos, neighbor_node(oiter));
            //    for (vector<ipoint2>::iterator i = pos.begin(); i != pos.end(); ++i) *i -= box.ll;
            //    sort(pos.begin(), pos.end());
            //    posset.insert(pos);
            //}
            cout << "  Class #" << ++count << " (" << ((cpart_data*)(p->data))->name << "): number of objects = " << p->count_neighbors(prevname) /*posset.size()*/ << endl;
        }
    }

    delete library;
}

// Does the following:
//  - Delete all parts on 'layer' having 1 subpart (e.g. only central part)
void clean_library(const config_dictionary& cfg, const char* fname)
{
    int thresh = cfg.get_value_int("threshold", 2);

    cerr << "Doing cleaning of the library" << endl;

    part_lib* library;

    read_library(fname, library);
    if (library == nullptr) {
        cerr << "Can not open library '" << fname << "'" << endl;
        return;
    }

    int layer = cfg.get_value_int("library", 5);
    vector<node*>& parts = library->parts[layer];
    vector<int> todelete;
    
    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* p = *iter;
        lib_data* pd = (lib_data*)p->data;
        set<node*> nb;

        p->get_neighbor_set(atom("lyrSrc"), nb);
        if (nb.size() <= thresh)
            todelete.push_back(pd->type);
    }

    library->delete_parts(layer, todelete);
    library->save(cfg.get_value_string("out_library", "newlib.plb"));
    cerr << todelete.size() << " parts on layer " << layer << " deleted." << endl;

    delete library;
}

void split_layer(const config_dictionary& cfg, const char* fname)
{
    typedef vector<irectangle2*> rectangle_set_t;

    part_lib* library;
    int layer = cfg.get_value_int("layer", 5);
    int splitn = cfg.get_value_int("split_n", 3);
    string outname = cfg.get_value_string("out_library", "newlib.plb");

    read_library(fname, library);
    if (library == nullptr || library->max_layer_index() < layer) return;

    part_lib* newlib = (part_lib*)library->get_copy();
    vector<node*>& parts = library->parts[layer];
    vector<node*>& partsm1 = library->parts[layer - 1];
    int backname = atom("lyrCenterBack");
    int srcname = atom("lyrSrc");
    int maxmap = 0; // D E B U G
    int maxup = 0; // D E B U G
    int avguporig = 0;
    int avgupnew = 0;
    double maxratio = 0.0;

    newlib->delete_parts_geq(layer);

    for (vector<node*>::iterator iter = partsm1.begin(); iter != partsm1.end(); ++iter) {
        node* cn = *iter;
        lib_data* cnd = (lib_data*)cn->data;
        list<irectangle2> rlist;
        irectangle2 box;
        int pcount = 0; // D E B U G
 
        // Find the bounding rectangle of all parts with center 'cn'.
        foreach_neighbor(cn, backname, citer) {
            node* p = neighbor_node(citer);

            foreach_neighbor(p, srcname, siter) {
                part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
                
                box.eat(ed->x, ed->y);
            }
            pcount++; // D E B U G
        }
        
        // Split the rectangle
        box.ur.x += 1;
        box.ur.y += 1;
        box.split(rlist, splitn, splitn);

        // D E B U G
        avguporig += pcount;
        int maxcontents = 0;
        if (pcount > maxup) maxup = pcount;
        //cerr << "|^|: " << pcount << ";  Bounding rectangle: " << box;
        
        typedef map<vector<irectangle2*>, vector<node*> > map_t;

        map_t set2nodes;

        // Make the 'set2nodes' map; i.e. giving the correspondence between sets S of 
        // rectangles and sets of parts with center 'cn' and subparts in *all* rectangles 
        // from S.
        foreach_neighbor(cn, backname, citer) {
            node* p = neighbor_node(citer);
            set<irectangle2*> rset;

            foreach_neighbor(p, srcname, siter) {
                part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
                irectangle2* rectp = find_rectangle(rlist, ed->x, ed->y);
                
                if (rectp != nullptr) rset.insert(rectp); // Though not being in any rectangle would be strange!
            }
    
            vector<irectangle2*> rsetv(rset.begin(), rset.end());

            map_t::iterator mapiter = set2nodes.insert(map_t::value_type(rsetv, vector<node*>())).first;
            mapiter->second.push_back(p);
            
            if ((int)mapiter->second.size() > maxcontents) maxcontents = mapiter->second.size();
        }

        // D E B U G
        //cerr << "   : " << set2nodes.size() << "  cont: " << maxcontents << endl;
        if ((int)set2nodes.size() > maxmap) maxmap = set2nodes.size();
        if ((double)maxmap/maxup > maxratio) maxratio = (double)maxmap/maxup;
        avgupnew += (int)set2nodes.size();
        // END D E B U G

        // Resize the rectangles a bit (to allow for some variance)
        for (list<irectangle2>::iterator riter = rlist.begin(); riter != rlist.end(); ++riter) {
            riter->grow(5);
        }

        // Add new parts to library
        for (map_t::iterator mapiter = set2nodes.begin(); mapiter != set2nodes.end(); ++mapiter) {
            const map_t::key_type& rpset = mapiter->first;
            map_t::mapped_type& nset = mapiter->second;
            vector<irectangle2> rset;

            for (map_t::key_type::const_iterator riter = rpset.begin(); riter != rpset.end(); ++riter) {
                rset.push_back(**riter);
            }
            add_r_part_to_library(newlib, layer, cnd->type, rset, nset);
        }

    }

    // D E B U G
    cerr << "Max original up: " << maxup << "  max merged up: " << maxmap << 
        "  max ratio (up/merged): " << maxratio << endl;
    cerr << "Avg originalup: " << (double)avguporig/partsm1.size() << "  avg newup: " << 
        (double)avgupnew/partsm1.size() << endl;

    newlib->save(outname);

    delete newlib;
    delete library;
}

void test2(const string& fname)
{
    vector<itriple> test;
    vector<vector<itriple> > result;

    test.push_back(itriple(1, 2, 100));
    //test.push_back(itriple(4, 5, 101));
    //test.push_back(itriple(-4, 2, 102));
    //test.push_back(itriple(13, 15, 103));
    //test.push_back(itriple(11, -5, 104));

    partition_tpoints(result, test);

    cerr << result << endl;

    return;

    part_lib* library;

    read_library(fname, library);
    if (library == nullptr) return;

    ofstream os("c:\\temp\\test2output.m");
    vector<node*>& parts = library->parts[4];
    int namecb = atom("lyrCenterBack");

    os << '{';
    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* p4 = *iter;
        part_data* p4d = (part_data*)p4->data;
        bool first = true;
        
        if (iter != parts.begin()) os << ',';
        os << '{';
        foreach_neighbor(p4, namecb, niter) {
            if (first) first = false; else os << ',';

            if (neighbor_node(niter)->is_attr_set(R_PART_ATTR)) {
                r_part_data* nd = (r_part_data*)neighbor_node_data(niter);
                
                os << '{' << nd->type << ',' << nd->layer << ',' << '{';
                for (vector<irectangle2>::iterator riter = nd->rlist.begin(); riter != nd->rlist.end(); ++riter) {
                    if (riter != nd->rlist.begin()) os << ',';
                    riter->mma_write(os);
                }
                os << '}' << '}';
            } else {
                part_data* nd = (part_data*)neighbor_node_data(niter);
                os << '{' << nd->type << ',' << nd->layer << '}' << endl;
            }
        }
        os << '}' << endl;
    }
    os << '}';
    os.close();
    delete library;
}


void test2_semiold(const string& fname)
{
    layer1_result* res;

    read_layer1_result(res, fname);
    if (res != nullptr) {
        cout << "Read successfully \'" << fname << '\'' << endl;

        int slayer = 4;
        int ename = atom("toPrevLayer");

        vector<node*>& s_nodes = res->shape_nodes[slayer];

        for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
            node* n = *iter;

            foreach_neighbor(n, ename, niter) {
                layer1_data* nd = (layer1_data*)neighbor_node_data(niter);
                edge_data_ip2* ed = (edge_data_ip2*)neighbor_edge_data(niter);

                if (ed != nullptr) {
                  //  ed->data
                }
            }
        }


        delete res;
    }
}

void test2_x(const string& fname)
{
    layer1_result* res;

    read_layer1_result(res, fname);
    if (res != nullptr) {
        cout << "Read successfully \'" << fname << '\'' << endl;

        int nbname = atom("toNeighbor");
        vector<node*>& s_nodes = res->shape_nodes[0];
        set<node*> result;
        set<node*> initset(s_nodes.begin(), s_nodes.end());

        if (!res->grid(0)) res->init_grid(0);
        res->connect_neighbors_circular(s_nodes, 5, 5, 0, nbname);
        //reduce_node_set(result, initset, atom("toNeighbor"));

        img imresult(res->x_size(0), res->x_size(0));

        for (set<node*>::iterator iter = result.begin(); iter != result.end(); ++iter) {
            layer1_data* nd = (layer1_data*)(*iter)->data;
            imresult.draw_big_point(nd->x, nd->y, 1.0);
        }

        imresult.save("c:\\temp\\test.png");
    
        //for (int i = 0; i < 10000; ++i) {
        //    layer1_result* resx = (layer1_result*)res->get_copy_s();
        //    resx->save(string("c:\\temp\\x") + i + string(".lyx"));
        //    delete resx;
        //}
        delete res;
    }
}

//void test2(const string& fname)
//{
//    layer1_result* res;
//
//    read_layer1_result(res, fname);
//    if (res != nullptr) {
//        cout << "Read successfully \'" << fname << '\'' << endl;
//    
//        int count = 0;
//
//        for (vector<node*>::iterator iter = res->shape_nodes[5].begin(); iter != res->shape_nodes[5].end(); ++iter) {
//            node* n = *iter;
//
//            do {
//                layer1_data* nd = (layer1_data*)n->data;
//                set<node*> nset;
//                int all = 0, layer0 = 0, hypo = 0;
//
//                res->recurse_from_node(n, atom("toPrevLayer").get_index(), nset);
//                for (set<node*>::iterator i = nset.begin(); i != nset.end(); ++i) {
//                    node* nn = *i;
//                    layer1_data* nnd = (layer1_data*)nn->data;
//                    
//                    ++all;
//                    if (nn->is_attr_set(HYPO_NODE_ATTR))
//                        hypo += 1;
//                    if (nnd->z == 0) 
//                        ++layer0;
//                }
//                cout << "Node #" << count << "  type = " << nd->m << "  value = " << nd->val << endl;
//                cout << "  all = " << all << "  layer0 = " << layer0 << "  hypo = " << hypo;
//                cout << "  hypo/layer0 = " << (double)hypo/layer0 << endl;
//                ++count;
//                n = nd->next;
//            } while (n != nullptr);
//        }
//
//    }
//}


/*void test1()
{
    layer1_result* res;

    read_layer1_result(res, "D:\\work\\data\\optimization\\ferrari\\layerxf\\fer_mugs_muki_2.ly6");
    if (res == nullptr) return;
    for (layer1_result::iter_t iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
        node* n = *iter;
        if (n->is_attr_set2(HYPO_NODE_ATTR)) {
            layer1_data* nd = (layer1_data*)n->data;
            cerr << nd->z << ' ' << nd->m << endl;
        }
    }
    delete res;
}*/

void __test1()
{
    cerr << "TEST1" << endl;

/*    list<string> files;
    list<layer1_result*> reslist;
    clock_t start, end;

    list_directory(files, "C:\\work\\data\\learning\\pos0\\layer2\\*.ly2");

    list<string>::iterator iter = files.begin();
    int i = 0;
    while (i < 1) {
        layer1_result* res;

        read_layer1_result(res, "C:\\work\\data\\learning\\pos0\\layer3\\learn_circle_1_0.ly3");
        //reslist.push_back(res);
        if (!res->grid(1)) res->init_grid(1);

        cerr << "Start copy" << endl;
        start = clock();
        //for (int j = 0; j < 40; ++j) {

            //res->delete_nodes(j, 1);
        //}
        layer1_result* g = (layer1_result*)res->get_induced_copy(pred);
        end = clock();
        cerr << "End copy" << endl;
        cerr << "Copy made in " << (double)(end - start)/CLOCKS_PER_SEC << " sec" << endl;
        g->save("c:\\work\\copy.ly3");

        delete res;
        delete g;

        ++iter; ++i;
    }
    for (list<layer1_result*>::iterator j = reslist.begin(); j != reslist.end(); ++j) {
        delete *j;
    }

    cerr << endl;
    cerr << -1 << "      " << (-1 << 16) << endl;
    cerr << 1 << "      " << (1 << 16) << endl;


    // ----------- PERMUTATION -----------
    vector<int> perm;

    for (int i = 1; i <= 4; ++i) perm.push_back(i);
    do {
        cerr << perm << endl;
    } while (next_permutation(perm));
*/
    // ----------- CLUSTERING -----------

    vector<ipoint2> vec;
    vec.push_back(ipoint2(1, 2));
    vec.push_back(ipoint2(1, 5));
    vec.push_back(ipoint2(16, 19));
    vec.push_back(ipoint2(2, 3));
    vec.push_back(ipoint2(17, 18));
    vec.push_back(ipoint2(1, 3));
    vec.push_back(ipoint2(20, 7));

    vector<vector<ipoint2> > result;
    double d;

    d = hierarchical_clustering(result, vec, ipoint2_set_distance);
    while (result.size() > 1) {
        cout << "Min distance: " << d << endl;
        cout << result << endl;
        d = hierarchical_clustering_iter(result, ipoint2_set_distance);
    }

    // -------- FERRARI GROUNDTRUTH -----------
    ifstream is("E:\\work\\data\\ferrari\\Mugs\\jazzburger_mugs.groundtruth");
    img im(string("E:\\work\\data\\ferrari\\Mugs\\jazzburger.jpg"));
    double x1, y1, x2, y2;

    if (is.is_open()) {
        is >> x1 >> y1 >> x2 >> y2;
        im.draw_box(irectangle2((int)x1, (int)y1, (int)x2, (int)y2), 1.0);
        im.save_normalized("test.bmp");
    }

    // ------ STRING PREFIX -------
    string teststr("fer_mug03_1.ly1");
    cout << get_prefix(teststr, "_", 2) << endl;
}

//cr_optimization1.cfg mm_mug*_1.ly1

void print_statistics(layer1_result* res, int layer)
{
    int name = atom("toNextLayer0").get_index() + layer;
    set<node*> covering_nodes;
    set<node*> uncovered;
    set<node*> tmp;
    set<ipoint2> coo_set;

    res->partition_by_neighbors(covering_nodes, uncovered, 
            res->shape_nodes[0].begin(), res->shape_nodes[0].end(), name);

    cerr << "layer #" << layer << endl;
    cerr << "-----------------------------" << endl;
    cerr << "# of uncovered: " << uncovered.size() << " of " << res->shape_nodes[0].size();
    cerr << " = " << (double)(uncovered.size())/res->shape_nodes[0].size() << "%" << endl;
    cerr << "# of covered nodes: " << covering_nodes.size() << endl;

    for (set<node*>::iterator iter = uncovered.begin(); iter != uncovered.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;
        coo_set.insert(ipoint2(nd->x, nd->y));
    }

    cerr << "uncovered points: " << coo_set.size() << endl;


    cerr << "how many uncovered nodes are covered by nodes of a particular layer?:" << endl;
    for (int i = 1; i <= layer; ++i) {
        name = atom("toNextLayer0").get_index() + i;
        tmp.clear();
        res->get_neighbors(tmp, uncovered.begin(), uncovered.end(), name);
        cerr << "  by layer " << i << ": " << tmp.size() << endl;
    }
}

struct check_type {
    set<int> types;

    bool operator()(node* n) const
    {
        layer1_data* nd = (layer1_data*)n->data;
        return types.find(nd->m) != types.end();
    }
};


void print_statistics(layer1_result* res, int layer, int p, int q)
{
    int name = atom("toNextLayer0").get_index() + layer;
    set<node*> tmpp, tmpq;
    check_type predicate;

    predicate.types.insert(p);
    res->select_by_neighbors(tmpp, res->shape_nodes[0].begin(), res->shape_nodes[0].end(),
        name, predicate);
    cerr << "#nbs of " << p << " = " << tmpp.size() << endl;
    predicate.types.clear();
    predicate.types.insert(q);
    res->select_by_neighbors(tmpq, res->shape_nodes[0].begin(), res->shape_nodes[0].end(),
        name, predicate);
    cerr << "#nbs of " << q << " = " << tmpq.size() << endl;
    cerr << "#nbs of " << q << " and " << p << " = " << intersection_size(tmpp, tmpq) << endl;

}

// test version 
// to be "cleaned up" of the following
//   - add_reconstruction_edges
//
void change_covered_responses(layer1_result* res, int layer, double thresh, double newvalf) 
{
    int name = atom("toNextLayer0").get_index() + layer;
    int to_0 = atom("toLayer0").get_index();
    
    for (int l = 1; l <= layer; ++l)
        res->add_reconstruction_edges(l);

    set<node*> covered, uncovered;
    vector<node*>& s_nodes = res->shape_nodes[layer];

    // get covered nodes
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        do {
            n->add_to_neighbor_set(to_0, covered);
            n = ((layer1_data*)n->data)->next;
        } while (n != nullptr);
    }
    
    // reduce responses on prevoius layers
    for (int l = layer - 1; l > 0; --l) {
        vector<node*>& s_nodes = res->shape_nodes[l];

        for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
            node* n = *iter;
            layer1_data* nd;
            set<node*> neighbors;

            do {
                nd = (layer1_data*)n->data;
                n->get_neighbor_set(to_0, neighbors);
                if (intersection_size(neighbors, covered) > (int)(thresh * neighbors.size()))
                    nd->r.set_response(R_RESPONSE, nd->val() * newvalf);
                n = nd->next;
            } while (n != nullptr);
        }
    }
    
}

void X_X_test2(const string& fname)
{
    layer1_result* res;
    list<string> files;

    list_directory(files, fname);

    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cerr << "Processing " << *file << "..." << endl;
        read_layer1_result(res, *file);
        if (res) {

            layer1_result* res2;
            {
                START_TIMER(t);
                res2 = (layer1_result*)res->get_copy_s();
                STOP_AND_SHOW(t, "copy through memstream");
                img* im = res2->get_image_reconstructed(0, 0, vector<int>());
                im->save("bla.png");                  
                delete im;
                delete res2;
            }
            {
                START_TIMER(t);
                res2 = (layer1_result*)res->get_induced_copy(graph::copy_predicate_t());
                STOP_AND_SHOW(t, "copy using graph copy function");
                delete res2;
            }
            //cerr << "Size: " << os.pcount();
            delete res;
        }
    }
}

/*void test2(string fname)
{
    config_dictionary cfg;

    cfg.from_string(fname.c_str());
    cerr << "tkey = " << cfg.get_value_string("tkey", "undefined") << endl;
}*/

/*void test2(const string& fname)
{
    layer1_result* res;
    list<string> files;
    list<layer1_result_ptr> results;
    int count = 0;

    list_directory(files, fname);
    layer1_result_ptr::set_max_mem_objects(50);
    cerr << "Max streamable pointers in memory: " << layer1_result_ptr::get_max_mem_objects() << endl;
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cerr << "Processing " << *file << "..." << endl;
        read_layer1_result(res, *file);
        results.push_back(layer1_result_ptr(res));
        cerr << "Memory consumed: " << get_working_set_size() << endl; 
    }
    layer1_result_ptr some_ptr(results.back());
    for (list<layer1_result_ptr>::iterator iter = results.begin(); iter != results.end(); ++iter) {
        layer1_result* res = (layer1_result*)iter->get_ptr();
        img* im = res->get_image_reconstructed(0, 0, vector<int>());
        string s = string("c:\\temp\\bla") + (count++) + ".png";

        im->save(s);
        delete im;
    }
    for (list<layer1_result_ptr>::iterator iter = results.begin(); iter != results.end(); ++iter) {
        cerr << "Count: " << iter->get_count() << ", in memory: " << iter->in_memory() << endl;
    }
}*/

/*void test2(const string& fname)
{
    layer1_result* res;
    list<string> files;
    int layer = 1;

    list_directory(files, fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        read_layer1_result(res, *file);
        if (res != nullptr) {
            cerr << "Processing " << *file << "..." << endl;

            set<node*> covered, uncovered, layer0, result;
            set<node*> all(res->shape_nodes[layer].begin(), res->shape_nodes[layer].end());

            res->get_neighbors(covered, res->shape_nodes[layer + 1].begin(), 
                res->shape_nodes[layer + 1].end(), atom("toPrevLayer").get_index());
            set_difference(uncovered, all, covered);
            res->add_reconstruction_edges(layer);
            res->get_neighbors(layer0, uncovered.begin(), uncovered.end(), atom("toLayer0").get_index());
            res->get_neighbors(result, layer0.begin(), layer0.end(), atom("toNextLayer0").get_index() + layer);

            for (set<node*>::iterator iter = result.begin(); iter != result.end(); ++iter) {
                layer1_data* nd = (layer1_data*)(*iter)->data;

                cout << nd->x << ' ' << nd->y << ' ' << nd->z << ' ' << nd->m << endl;
            }
            cout << "0 0 0 0" << endl;
            for (set<node*>::iterator iter = layer0.begin(); iter != layer0.end(); ++iter) {
                layer1_data* nd = (layer1_data*)(*iter)->data;

                cout << nd->x << ' ' << nd->y << ' ' << nd->z << ' ' << nd->m << endl;
            }

            delete res;
        }
    }
}*/

void x_test2(const string& fname)
{
    layer1_result* res;
    list<string> files;

    list_directory(files, fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        read_layer1_result(res, *file);
        cerr << "Processing " << *file << "..." << endl;
        if (res != nullptr) {
            change_covered_responses(res, 4, 0.80, 0.2);
            res->save("reduced_" + *file);
            delete res;
        }
    }
}

void ___test2(const string& fname)
{
    img tmp = img::string_to_img(fname);
    
    tmp.save("c:\\work\\str.bmp");
    tmp = tmp.get_resized(-200, -300);
    tmp.save("c:\\work\\strr.bmp");
}

//void __test2(const string& fname)
//{
//    layer1_result* res;
//    part_lib* library;
//    config_dictionary cfg("C:\\work\\data\\optimization\\cr_optimize4.cfg");
//
//    read_layer1_result(res, fname);
//    read_library("C:\\work\\data\\optimization\\lib4.plb", library);
//
//    layer1_optimization optimizer(cfg, 4);
//
//    if (res != nullptr) {
//        print_statistics(res, 4);
//
//        //optimizer.add_to_test_set(res);
//        
//        for (int i = 0; i < 1; ++i) {
//            optimizer.make_step_1();
//            optimizer.print_final_parts();
//            PRINT_INFO("Cover value: " << optimizer.get_final_value());
//        }
//        cerr << endl;
//        optimizer.keep_final_parts();
//        res->save("c:\\work\\data\\optimization\\x.ly5");
//        print_statistics(res, 4, 0, 1);
//        delete res;
//    }
//}

void test2_() 
{
    layer1_result* res;
    
    read_layer1_result(res, "C:\\work\\data\\optimization\\res20.ly5");
    if (res != nullptr) {
        for (int layer = 1; layer <= res->max_layer_index(); ++layer) {
            print_statistics(res, layer);
        }

        delete res;
    }
}

void save_catlayer(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing save_catlayer" << endl;

    part_lib* library;

    read_library(fname, library);
    if (library == nullptr) return;

    vector<node*>& parts = library->parts[library->max_layer_index()];

    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* p = *iter;
        lib_data* pd = (lib_data*)p->data;

        cout << pd->type;
        forall_neighbors(p, niter) {
            lib_data* pnd = (lib_data*)neighbor_node_data(niter);

            cout << ',' << pnd->type;
        }
        cout << endl;
    }

    delete library;

}

void change_labels(const config_dictionary& cfg, const char* fname)
{
    cerr << "Doing change_labels" << endl;  
    
    layer1_result* res;
    list<string> files;
    string srcdir = cfg.get_value_string("src_dir", "");
    string destdir = cfg.get_value_string("out_dir", "");
    int layer = cfg.get_value_int("layer", -1);
    vector<int> typemap;


    cfg.get_value(typemap, "type_map", true);
    end_dir(srcdir);
    end_dir(destdir);
    list_directory(files, srcdir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        string name = srcdir + *file;

        cerr << "Processing " << *file;
        read_layer1_result(res, name);
        if (res == nullptr) cerr << " loading failed!";
        else {
            int l = (layer < 0 || layer > res->max_layer_index()) ? res->max_layer_index() : layer;
            vector<node*>& s_nodes = res->shape_nodes[l];

            for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
                node* n = *iter;

                while (n != nullptr) {
                    layer1_data* nd = (layer1_data*)n->data;

                    if (nd->m >= 0 && nd->m < typemap.size()) nd->m = typemap[nd->m];
                    n = nd->next;
                }
            }
            res->save(destdir + *file);
            cerr << " done.";
            delete res;
        }
        cerr << endl;
    }

}

void test_compression(const config_dictionary& cfg, const char* fname) {
	
	string dir = cfg.get_value_string("test_dir","");
	string file = cfg.get_value_string("test_file_4","test1.ly3");


	long start_org = 0, end_org = 0, start_comprss = 0, end_comprss = 0;

	start_org = clock();
	
		cout << "reading original file" << endl;	
		layer1_result* res = (layer1_result*)streamable::read(dir + file);	

		cout << "saving as original file" << endl;
		res->save(dir + file + "._org_resaved");
	
	end_org = clock();

	cout << "done in " << (double)(end_org - start_org) /CLOCKS_PER_SEC << endl;
	cout << "compressing original file" << endl;
	res->save(dir + file + ".z", -1);

	start_comprss = clock();
	
		cout << "reading compressed file" << endl;
		layer1_result* res1 = (layer1_result*)streamable::read(dir + file + ".z");

		cout << "saving as original file" << endl;
		res1->save(dir + file + "._new");
	
	end_comprss = clock();
	
	cout << "done in " << (double)(end_comprss - start_comprss) /CLOCKS_PER_SEC << endl;


	cout << "done" << endl;

}

void test_rect_stream(const config_dictionary& cfg, const char* fname) {
	
	string file;
	cfg.get_value(file,"file",true);
	
	list<irectangle2> rectangles;
	read_groundtruth(rectangles, file);

	
	ifstream is(file.c_str());

	rectangle2<double> rect;

	

	is >> rect;
	

	cout << (int)rect.ll.x << " " << (int)rect.ll.y << " " << (int)rect.ur.x << " " <<(int)rect.ur.y << endl;
	cout << "is ok: " << is.fail() << endl;

	

	string ff;
	is >> ff;
	cout << "new cat: " << ff << endl;
	cout << "is ok: " << is.fail() << endl;
}

void test_suplea_learning(const config_dictionary& cfg, const char* patt)
{
    cerr << "Doing test_suplea_learning" << endl;  
    
    layer1_result* res;
    list<string> files;
    string srcdir = cfg.get_value_string("src_dir", "");
    vector<int> typemap;
    s_learning sl(cfg);

    end_dir(srcdir);
    list_directory(files, srcdir + patt);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        string name = srcdir + *file;

        cerr << "Processing " << *file;
        read_layer1_result(res, name);
        if (res == nullptr) cerr << " loading failed!";
        else {
	        sl.prepare_for_update(res);
	        sl.update_nbhoods(res);
            delete res;
        }
        cerr << endl;
    }
    sl.find_maxima();
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        string name = srcdir + *file;

        cerr << "Processing " << *file;
        read_layer1_result(res, name);
        if (res == nullptr) cerr << " loading failed!";
        else {
            sl.apply_to_result(res); 
            delete res;
        }
        cerr << endl;
    }
    sl.filter_sequences();
    sl.add_to_library(vector<int>(), sl.part_max_number);
    learn_g_distributions(sl, cfg, "g_distribution", srcdir, patt);

    sl.library->save_all("c:\\work\\testlib.png", sl.source_layer_index + 2, 0, -1, true);
    sl.library->save("c:\\work\\testlib.plb");

}

// Duplet learning

class d_learning {
protected:
    typedef map<ipoint2, imatrix> stat_map_t;

    stat_map_t statmap;
    int rfr2;     // (receptive field radius)^2 = (rfdim/2)^2

    // parameters
    int rfdim;    // dimension of the receptive field (rfdim x rfdim); default = 5

public:
    d_learning(const config_dictionary& cfg);
	bool update(int i, const ipoint2& ip, int j, const ipoint2& jp);

protected:
    void set_rfdim(int d);
    void init_cfg(const config_dictionary& cfg);
};

d_learning::d_learning(const config_dictionary& cfg)
{
    init_cfg(cfg);
}

void d_learning::set_rfdim(int d) 
{
    if (d <= 0) d = 5;
    if (d % 2 == 0) d += 1;
    rfdim = d;
    rfr2 = d*d/4;
}

void d_learning::init_cfg(const config_dictionary& cfg)
{
    set_rfdim(cfg.get_value_int("rf_dimension", 5));
}

// Update statistics: type i at ip vs. type j at jp.
// Return value: true if stat was updated and false otherwise (distance too big for rf).
bool d_learning::update(int i, const ipoint2& ip, int j, const ipoint2& jp)
{
    if (ip.distance2(jp) > rfr2) return false;

    stat_map_t::iterator iter = statmap.find(ipoint2(i, j));

    if (iter == statmap.end())
        iter = statmap.insert(stat_map_t::value_type(ipoint2(i, j), imatrix(rfdim, rfdim, 0))).first;

    int x = jp.x - ip.x + rfdim/2;
    int y = jp.y - jp.x + rfdim/2;

    iter->second.at(x, y) += 1;
    return true;
}



void test0(const char* cfg_file, const char* pattern, const char* params)
{
    try { 
        config_dictionary cfg(cfg_file);
        cfg.from_string(params);

        string mode = cfg.get_value_string("mode", "onelayer");

        if (mode == "onelayer") cerr << "mode = onelayer not supported" << endl;
        //else if (mode == "cr") optimize_layer2(cfg, pattern);
        else if (mode == "last") cerr << "mode = last ---> slearning" << endl;
        else if (mode == "overall") test_overall_optimization(cfg, pattern);
        else if (mode == "validate") cerr << "mode = validate ---> slearning" << endl;
        else if (mode == "check") test_check(cfg, pattern);
        else if (mode == "multicheck") test_multicheck(cfg, pattern);
        else if (mode == "merge") merge_scales(cfg, pattern);
        else if (mode == "statistics") part_hits_statistics(cfg, pattern);
        else if (mode == "slearning") slearning(cfg, pattern);
        else if (mode == "test") test_3(cfg, pattern);
        else if (mode == "splitlayer") split_layer(cfg, pattern);
        else if (mode == "insert_rlayers") insert_rlayers(cfg, pattern);
        else if (mode == "back_statistics") back_statistics(cfg, pattern);
        else if (mode == "clean_library") clean_library(cfg, pattern);
        else if (mode == "layer_statistics") layer_statistics(cfg, pattern);
        else if (mode == "coarse2fine") coarse2fine(cfg, pattern);
        else if (mode == "coarse2fineT1") coarse2fineT1(cfg, pattern);
        else if (mode == "save_catlayer") save_catlayer(cfg, pattern);
        else if (mode == "change_labels") change_labels(cfg, pattern);
        else if (mode == "merge2partsM") merge2partsM(cfg, pattern); // manual 
        else if (mode == "merge2partsT") merge2partsT(cfg, pattern); // train
        else if (mode == "merge2parts") merge2parts(cfg, pattern);
        else if (mode == "powercalc") powercalc(cfg, pattern);
        else if (mode == "svm_test") svm_test(cfg, pattern);
        else if (mode == "cross_layer_svm") cross_layer_svm(cfg, pattern);
        else if (mode == "library_layer_statistics") library_layer_statistics(cfg, pattern);
        else if (mode == "compare_layer1_results") compare_layer1_results(cfg, pattern);
        else if (mode == "best_response_distribution") best_response_distribution(cfg, pattern);
		else if (mode == "compressed_streaming") test_compression(cfg, pattern);
		else if (mode == "test_rect_stream") test_rect_stream(cfg, pattern);
        else if (mode == "suplea_learning") test_suplea_learning(cfg, pattern);

    } catch (const libhop_exception& e) {
		cerr << e.what() << endl;
    } catch (const exception& e) {
		cerr << "General exception '" << typeid(e).name() << "' with message: '" << e.what() << "'" << endl;
	} 
}

int main(int argc, char* argv[])
{
#if defined WIN32 | defined WIN64
    _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
    _CrtSetBreakAlloc(338);
#endif
    init_atoms();
    init_streaming();
    
    //int* boo = new int;

    cerr << "test (" __DATE__ " / " __TIME__ ")" << endl;
    srand((unsigned)time(0)); 
    //srand(2009);

    cerr << argc << " arguments given" << endl;

    switch (argc) {
        case 1 : test1(); break;
        case 2 : test2(argv[1]); break;
        /*case 2 : export_hierarchy(argv[1], "m"); break;*/
        case 3 : test0(argv[1], argv[2], ""); break;
        case 4 : test0(argv[1], argv[2], argv[3]); break;
        //case 4 : test
        /*case 3 : export_hierarchy(argv[1], argv[2]); break;*/
        default: 
            cout << "test files ['m' | 'n']" << endl;
            cout << "  (m = Mathematica, n = Pajek)" << endl;
    }
    return 0;
}

