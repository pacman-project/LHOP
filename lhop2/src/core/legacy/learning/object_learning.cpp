
#include <opencv2/opencv.hpp>
#include <queue>
#include "object_learning.h"
#include "utils/graphs/graph_utils.h"

// global functions	
///////////////////////////////////////////////////////////////////////////////



// obj_learning
///////////////////////////////////////////////////////////////////////////////


obj_learning::obj_learning(const ConfigDictionary& cfg) :
    library(nullptr),
    libraryD(nullptr),
    creator(nullptr)
{
    cfg_init(cfg);
}

obj_learning::~obj_learning()
{
    if (library != nullptr) delete library;
    if (libraryD != nullptr) delete libraryD;
    if (creator != nullptr) delete creator;
}

void obj_learning::set_library(part_lib* plibrary)
{
    library = plibrary;
    if (library) {
        libraryD = (part_lib*)library->get_copy_s();
        libraryD->delete_parts_geq(layer + 1);
	}
}

void obj_learning::cfg_init(const ConfigDictionary& cfg)
{
    string libname;
    
    cfg.getValue(layer, "layer", true);

    cfg.getValue(r_response_threshold, "r_response_threshold", 0.0);
    cfg.getValue(g_response_threshold, "g_response_threshold", 0.0);
    cfg.getValue(rr_response_threshold, "rr_response_threshold", 0.0);
    contraction = cfg.getValueDouble("contraction", 1.0);
    cfg.getValue(max_objects, "max_objects", 4);
    cfg.getValue(max_add, "max_add", 0);
    if (max_add[0] <= 0) max_objects.set_val(max_add[0]);
    cfg.getValue(max_cluster_n, "max_cluster_n", 6);
    cfg.getValue(min_cluster_n, "min_cluster_n", 4);
    cfg.getValue(cluster_size, "cluster_size", 10);
    gaussian_dim = cfg.getValueInt("gaussian_dim", 5);
    gaussian_sigma = cfg.getValueDouble("gaussian_sigma", 2.0);
    cfg.getValue(cluster_member_threshold, "cluster_member_threshold", 0.5);
    cfg.getValue(cover_threshold, "cover_threshold", 0.5);
    cover_threshold0 = cfg.getValueDouble("layer0_cover_ratio_threshold", 0.0);
    validation_threshold = cfg.getValueDouble("validation_threshold", 1.5);
    cfg.getValue(intersection_threshold, "intersection_threshold", 0.2);
    max_depth = cfg.getValueInt("max_depth", 0);
    reduce_radius = cfg.getValueInt("reduce_radius", 3);
}

// recurse2 Man
//  typedef container_t
//  node* get_node(container_t::const_iterator) - returns the node part of the data
//  bool check_neighbor(const container_t&, node* n, node* nn) - returns true if the 
//     neighbor nn is "good" and adds it to the container
//  void insert_leaf(node* n) - insert n to the result
//

struct ofr_recurse {
    typedef set<node*> container_t;
    typedef container_t::const_iterator const_iterator_t;
    typedef container_t::iterator iterator_t;

    container_t result;
    layer1_result* res;

    ofr_recurse(layer1_result* r) : res(r) 
    { 
        if (!res->grid(0)) res->init_grid(0); 
    }

    void reset() { result.clear(); }

    node* get_node(const container_t::value_type& val) { return val; }

    bool check_neighbor(container_t& neighbors, const container_t::value_type& n, edge_data* ed, node* nn) 
    { 
        neighbors.insert(nn); 
        return true; 
    }

    void insert_leaf(const container_t::value_type& n) 
    { 
        layer1_data* nd = (layer1_data*)n->data;

        if (nd->z == 0)
            result.insert(res->node_at(nd->x, nd->y, 0));
    }
};


struct list_item_t {
    int size;
    double response;
    node* n;

    list_item_t(int s, double r, node* pn) : size(s), response(r), n(pn) { }
    bool operator<(const list_item_t& i2) const { return response < i2.response; }
    bool operator>(const list_item_t& i2) const { return response > i2.response; }
};

void obj_learning::object_from_result(objects_t& result, layer1_result* res, const scmap_t& scmap,
    int layer, const set<ipoint2>& covered, const irectangle2& gtr)
{
    typedef set<node*> node_set_t;
    typedef set<ipoint2> point_set_t;
    typedef map<node*, point_set_t> map_t;
    typedef pair<ipoint2, set<int> > tabu_item_t;

    int depth = this->layer - layer;

    if (layer < 0 || layer >= (int)res->shape_nodes.size() || depth > max_depth) return;

    if (!res->grid(layer)) res->init_grid(layer);

    vector<node*>& s_nodes = res->shape_nodes[layer];
    int ename = EdgeConnection::TO_PREV_LAYER;
    //int nbname = EdgeConnection::TO_NEIGHBOOR;
    map_t node2set;
    list<tabu_item_t> tabu;
    list<list_item_t> size2node;
    ofr_recurse rmanager(res);
    ofr_recurse::container_t start;
    point_set_t proj0;
    point_set_t pset;
    int projsize = 0;

    if (s_nodes.empty()) return;

    vector<point_set_t> ly0sets;

    library->get_regions(1, ly0sets);

    // Initialize node->set map, projection
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        while (n != nullptr) {
            layer1_data* nd = (layer1_data*)n->data;

            if (nd->r(R_RESPONSE) >= r_response_threshold[depth] && nd->r(G_RESPONSE) >= g_response_threshold[depth] &&
                    nd->r(RR_RESPONSE) >= rr_response_threshold[depth]) {
                pset.clear();
                start.clear();
                start.insert(n);
                rmanager.reset();
                res->recurse2(rmanager, start, ename);

                node_set_to_region_set(pset, rmanager.result, ly0sets, 0);
                //get_positions(pset, rmanager.result, 0);
                int intsize = intersection_size(pset, covered);

                if ((double)intsize/pset.size() < intersection_threshold[depth])
                    node2set.insert(map_t::value_type(n, pset));
            }

            n = nd->next;
        }
    }

    if (node2set.empty()) return;

    pset.clear();
    start.clear();
    start.insert(s_nodes.begin(), s_nodes.end());
    rmanager.reset();
    res->recurse2(rmanager, start, ename);
    //get_positions(pset, rmanager.result, 0);
    node_set_to_region_set(pset, rmanager.result, ly0sets, 0);
    projsize = (int)pset.size() - intersection_size(covered, pset);
    get_positions(proj0, rmanager.result, 0);

    pset.clear();
    node_set_to_region_set(pset, set<node*>(res->shape_nodes[0].begin(), res->shape_nodes[0].end()), ly0sets, 0);

    // D E B U G  -  D I S P L A Y
    cout << "projection ratio: " << (double)projsize/pset.size() << endl;

    // Initialize size2node list & sort it in descending order
    for (map_t::iterator miter = node2set.begin(); miter != node2set.end(); ++miter) {
        const layer1_data* nd = (layer1_data*)(miter->first)->data;

        size2node.push_back(list_item_t((int)miter->second.size(), nd->r(G_RESPONSE) , miter->first));
    }
    size2node.sort(greater<list_item_t>());

    int objectcount = 0;
    int topsize = size2node.front().size;

    cout << "==============================" << endl;  // D E B U G  -  D I S P L A Y
    cout << "depth: " << depth <<  "; 'model reconstruction nodes'/'layer 0 nodes' ="; // D E B U G  -  D I S P L A Y

    while (objectcount < max_objects[depth]) {
        obj_data_t object;
        node_set_t nset;
        point_set_t rset;
        list<list_item_t>::iterator liter = size2node.begin();

        // Find best cover using the greedy method ignoring the first item in tabu list
        while (liter != size2node.end() && 
            rset.size() < cover_threshold[depth] * projsize && (int)nset.size() <= max_cluster_n[depth]) {
            node* topn = liter->n;
            layer1_data* topnd = (layer1_data*)topn->data;

            if (tabu.empty() || ((tabu.front().first.distance2(topnd->x, topnd->y) > 3*3) &&  // !!!!!??????
                    tabu.front().second.find(topnd->m) == tabu.front().second.end())) {
                point_set_t& topset = node2set.find(topn)->second;

                if (intersection_size(topset, rset) <= (int)(intersection_threshold[depth] * topset.size())) {
                    nset.insert(topn);
                    rset.insert(topset.begin(), topset.end());
                }
            }
            ++liter;
        }

        //double modelvalue = (double)rset.size()/projsize;
        double modelvalue = (double)intersection_size(rset, proj0)/proj0.size();

        // D E B U G  -  D I S P L A Y
        cout << ' ' << modelvalue;

        if (!tabu.empty()) tabu.pop_front();

        // Create 'object' structure from the nodes from the 'nset' and adding 
        // "equivalent" nodes
        for (node_set_t::iterator iter = nset.begin(); iter != nset.end(); ++iter) {
            part_lib::cluster_data_t cd;
            path_map_t pm, pml;
            node* n = *iter;
            layer1_data* nd = (layer1_data*)n->data;
            point_set_t& mnset = node2set.find(n)->second;
            vector<node*> nbhood;
            set<int> nbm;

            cd.types.push_back(ipoint2(-1 - depth, nd->m));
            
            // Update geometry
            get_path_map(pm, res, scmap, n, true);
            add_to_path_map(cd.geo[nd->m], pm);

            nbm.insert(nd->m);
            res->get_neighbors_circular(nbhood, n, 4, 4, 0, true);
            for (vector<node*>::iterator i = nbhood.begin(); i != nbhood.end() && (int)cd.types.size() < cluster_size; ++i) {
                node* n2 = *i;
                layer1_data* nd2 = (layer1_data*)n2->data;
                
                if (nbm.find(nd2->m) == nbm.end()) {
                    map_t::iterator fiter = node2set.find(n2);

                    if (fiter != node2set.end()) {
                        point_set_t& mnset2 = fiter->second;

                        if (intersection_size(mnset, mnset2) >= cluster_member_threshold * mnset.size()) {
                            cd.types.push_back(ipoint2(-1 - depth, nd2->m));

                            // Update geometry
                            get_path_map(pm, res, scmap, n2, true);
                            add_to_path_map(cd.geo[nd2->m], pm);

                            nbm.insert(nd2->m);
                        }
                    }
                }
            }
            cd.pos = ipoint2(nd->x, nd->y);
            
            object.push_back(cd);
            tabu.push_back(tabu_item_t(ipoint2(nd->x, nd->y), nbm));
}
        
        // Validate object
        bool valid = false;

        if (!vset.empty() && creator != nullptr) {
            valid = validate_object(res, gtr, object);
            if (valid) result.push_back(objects_t::value_type(modelvalue, object));
        } else {
            valid = ((int)object.size() >= min_cluster_n[depth]) && (modelvalue >= cover_threshold[depth]);
            if (valid) result.push_back(objects_t::value_type(modelvalue, object));
            else {
                // Recursive step
                objects_t newobjects;
                
                object_from_result(newobjects, res, scmap, layer - 1, rset, gtr);

                objects_t::iterator noiter = newobjects.begin();
                int count = 0;

                while (noiter != newobjects.end() && count < max_add[depth]) {
                    result.push_back(objects_t::value_type(modelvalue, object));
                    result.back().second.insert(result.back().second.end(), noiter->second.begin(), noiter->second.end());
                }
            }
        }

        ++objectcount;
    }
    cout << endl;
}

void obj_learning::object_from_result(layer1_result* res, const scmap_t& scmap, const irectangle2& gtr)
{
    vector<double> resc;

    res->get_contractions(resc);
    for (int i = 0; i < (int)resc.size() - 1; ++i) 
        if (round(resc[i], 1) != round(library->contractions[i + 2], 1)) 
            throw new_libhop_exception("Contraction exception");
    if (round(contraction, 1) != round(resc.back(), 1)) 
        throw new_libhop_exception("Contraction exception");

    double cr0 = layer0_cover_ratio(res, layer);

    if (cr0 < cover_threshold0) {
        cout << "Layer " << layer << "covers " << cr0 << " < " << cover_threshold0 << " of layer 0. Skipping." << endl;
        return;
    }
    res->delete_edges(EdgeConnection::TO_LAYER0);

    objects_t newobjects;
    
    object_from_result(newobjects, res, scmap, layer, set<ipoint2>(), gtr);
    
    newobjects.sort(greater<objects_t::value_type>());
    objects_t::iterator noiter = newobjects.begin();
    int i = 0; 

    while (noiter != newobjects.end() && i < max_add[0]) {
        cout << '(' << noiter->first << ')';
        objects.push_back(noiter->second);
        ++noiter;
        ++i;
    }
    cout << ' ' << i << " new objects added" << endl;

}

bool obj_learning::validate_object(layer1_result* res, const irectangle2& gtr, const obj_data_t& object)
{
    int tr = 0, fa = 0;
    part_lib* newlib = (part_lib*)libraryD->get_copy_s();
    int counter = 0;

    cout << " V";
    add_to_library(newlib, object, "temporary");
    creator->set_library(newlib);
    if (!gtr.invalid()) 
        vset.push_front(validation_data_t(
            streamed_pointer((layer1_result*)res->get_copy_s()), 
            gtr)
        );
    for (list<validation_data_t>::iterator iter = vset.begin(); iter != vset.end(); ++iter) {
		layer1_result* rnew = (layer1_result*)iter->first.get();
        vector<double> v;
        list<irectangle2> gtr;

        creator->add_layer(rnew, layer + 2, 0);
        
        counter++;
        gtr.push_back(iter->second);
        rnew->check_with_groundtruth(v, gtr, layer + 1, set<int>(), 0.0);
        
        if (gtr.front().invalid() && !v.empty()) { cout << '!'; }  // D E B U G
        if (!v.empty()) 
            if (*min_element(v.begin(), v.end()) > 0.6) {
                ++tr; 
                cout << 'T'; // D E B U G
            } else {
                ++fa;
                cout << 'F'; // D E B U G
            }
        else cout << '.';
        delete rnew;
    }
    if (!gtr.invalid()) 
        vset.pop_front();
    creator->set_library(nullptr);
    delete newlib;
    cout << ' ' << tr << '/' << fa << endl;
    return validation_function_1(tr, fa);
}

void obj_learning::add_validation_data(const streamed_pointer& ptr, const irectangle2& gtruth)
{
    vset.push_back(validation_data_t(ptr, gtruth));
}

bool obj_learning::validation_function_1(int tr, int fa)
{
    return (fa == 0) ? tr > 0 : (double)tr/fa > validation_threshold;
}

// Reduces geometry (path maps) with class parameter 'reduce_radius'
void obj_learning::reduce_object_data(obj_data_t& od)
{
    for (obj_data_t::iterator iter = od.begin(); iter != od.end(); ++iter) {
        for (map<int, path_map_t>::iterator pmiter = iter->geo.begin(); pmiter != iter->geo.end(); ++pmiter) {
            pmiter->second = reduce_path_map(pmiter->second, reduce_radius);
        }
    }
}
// Adds the parts to the library
// Returns the number of parts added
int obj_learning::add_to_library(part_lib* lib, const list<obj_data_t>& objects, const string& name)
{
    vector<int> newparts;

    for (list<obj_data_t>::const_iterator iter = objects.begin(); iter != objects.end(); ++iter) {
        if (iter->empty()) continue;

        matrix<double> mask;

        gaussian_mask(gaussian_dim, gaussian_dim, gaussian_sigma, mask);

        vector<matrix<double>*> maskv(iter->size(), &mask);

        if (lib->find_equivalent_object_part(layer + 1, *iter, mask.width) < 0) {
            obj_data_t od = *iter;
            //normalize_object_data(od);
            reduce_object_data(od);

            int newm = lib->add_object_part(layer + 2, new part_data(), od, maskv, contraction);

            if (newm >= 0) newparts.push_back(newm);
        } 
    }
    lib->add_category(layer + 3, newparts, name);
    return (int)newparts.size();
}

int obj_learning::add_to_library(part_lib* lib, const obj_data_t& object, const string& name)
{
    list<obj_data_t> objects;

    objects.push_back(object);
    return add_to_library(lib, objects, name);
}

int obj_learning::add_to_library(const string& name)
{
    return add_to_library(library, objects, name);
}

// o_learning
///////////////////////////////////////////////////////////////////////////////

void find_maxima(ip2_vector& maxima, const ip2_vector& v)
{
    typedef pair<double, ipoint2> queue_item_t;
    typedef priority_queue<queue_item_t> queue_t;

    const int border = 5;
    const int radius = 5;
    const int max_maxima = 7;
    const double max_percent = 0.05;

    static img* blurmask = img::gaussian_mask(border, border, 0.75);

    maxima.clear();
    if (v.empty()) return;

    irectangle2 box = irectangle2::bounding_rectangle(v.begin(), v.end());
    matrix<double> m(box.x_dim() + 2*border + 1, box.y_dim() + 2*border + 1, 0.0);

    for (ip2_vector::const_iterator iter = v.begin(); iter != v.end(); ++iter) {
        m(iter->x - box.ll.x + border, iter->y - box.ll.y + border) += 1;
    }
    m = convolve_matrix(m, *blurmask);

    double max = m.maximum();
    queue_t queue;

    m.change_values_leq(max_percent * max, 0.0);
    for_each_xy_int (m, i, j) {
        double lmax = m.maximum_c((int)i, (int)j, radius, radius);

        if (lmax > 10E-6 && m(i, j) == lmax) 
            queue.push(queue_item_t(lmax, ipoint2(i, j)));
    }

    matrix<int> bm(m.width, m.height, 0);

    bm.set_region_c(border - box.ll.x, border - box.ll.y, radius, radius, 1);
    while (!queue.empty() && (int)maxima.size() < max_maxima) {
        ipoint2 p = queue.top().second;

        if (bm(p.x, p.y) == 0) {
            bm.set_region_c(p.x, p.y, radius, radius, 1);
            maxima.push_back(ipoint2(p.x + box.ll.x - border, p.y + box.ll.y - border));
        }
        queue.pop();
    }
    
}


void olv_object_part::print(std::ostream &os) const
{
    os << '{';
    for (int i = 0; i < size(); ++i) {
        if (i != 0) os << ',';
        os << '{' << type(i) << ',' << '{' << pos(i).x << ',' << pos(i).y << '}' << '}';
    }
    os << '}';
}

inline double div0(double a, int b) { return b == 0 ? a : a/b; }

struct mstat_compare {
    vector<pair<double, int> >* pstat;

    mstat_compare(vector<pair<double, int> >* ps) : pstat(ps) { }
    bool operator()(node* n1, node* n2) const 
    { 
        int nt1 = node_type(n1), nt2 = node_type(n2);
        return div0(pstat->at(nt1).first, pstat->at(nt1).second) > div0(pstat->at(nt2).first, pstat->at(nt2).second); 
    }
};

o_learning::o_learning(const ConfigDictionary& cfg) 
    : library(nullptr), srclibrary(nullptr), dupstat(), dupmap(), statmap()
{
    init_cfg(cfg);
}

o_learning::~o_learning()
{
    if (library) delete library;
    if (srclibrary) delete srclibrary;
}

void o_learning::reset()
{
    dupstat.clear();
    dupmap.clear();
    statmap.clear();
    vset.clear();
    rhits.clear();
}

void o_learning::init_cfg(const ConfigDictionary& cfg)
{
    string libname;    
    cfg.getValue(libname, "library", true);
	
	library = part_lib::read(libname);

    srclibrary = (part_lib*)library->get_copy_s();

    cfg.getValue(srclayer, "src_layer", true);
    cfg.getValue(contraction, "contraction", true);    
    catname = cfg.getValueString("category_name", "");
    max_cluster_n = cfg.getValueInt("max_cluster_n", 6);
    min_cluster_n = cfg.getValueInt("min_cluster_n", 4);
    cluster_size = cfg.getValueInt("cluster_size", 10);
    cluster_member_threshold = cfg.getValueDouble("cluster_member_threshold", 0.7);
    intersection_threshold = cfg.getValueDouble("intersection_threshold", 0.7);
    hit_ratio_threshold = cfg.getValueDouble("hit_ratio_threshold", 0.8);
    hit_threshold = cfg.getValueDouble("hit_threshold", 0.1);
    redundancy_threshold = cfg.getValueInt("redundancy_threshold", 10);
    type_bite_threshold = cfg.getValueInt("type_bite_threshold", 5);
    max_models = cfg.getValueInt("max_models", 500);
    gaussian_dim = cfg.getValueInt("gaussian_dim", 5);
    gaussian_sigma = cfg.getValueDouble("gaussian_sigma", 2.0);

    double gdim = cfg.getValueInt("gaussian_dim", 7);
    double gsigma = cfg.getValueDouble("gaussian_sigma", 2.0);
    gaussian_mask(gdim, gdim, gsigma, dist);

	inference_cfg.fromNamespacePriority(cfg, 1, "validation");
}

void o_learning::update_duplet_statistics(layer1_result* res)
{
    for (list<node*>::iterator niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;

        if (nd->z != srclayer) continue;
        for (list<node*>::iterator miter = res->nodes.begin(); miter != res->nodes.end(); ++miter) {
            node* m = *miter;
            layer1_data* md = (layer1_data*)m->data;

            if (md->z != srclayer) continue;
            
            dupmap_t::iterator fiter = dupstat.find(iipair(nd->m, md->m));

            if (fiter == dupstat.end()) 
                fiter = dupstat.insert(dupmap_t::value_type(iipair(nd->m, md->m), ip2_vector())).first;
            fiter->second.push_back(ipoint2(md->x - nd->x, md->y - nd->y));
        }

    }
}

void o_learning::make_duplets()
{
    cout << "(dupstat size: " << dupstat.size() << ')';

    dupmap.clear();
    for (dupmap_t::const_iterator iter = dupstat.begin(); iter != dupstat.end(); ++iter) {
        ip2_vector v = iter->second;
        ip2_vector mv;

        find_maxima(mv, v);
        if (!mv.empty()) 
            dupmap.insert(dupmap_t::value_type(iter->first, mv));
    }
}

void o_learning::update_statmap(layer1_result* res, const list<irectangle2>& gtrs)
{
    for (list<node*>::iterator iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;

        if (nd->z != srclayer) 
            continue;

        set<node*> nset;
        set<ipoint2> pset;
        double stat = 0.0;

        n->get_neighbor_set(EdgeConnection::TO_LAYER0, nset);
        
        irectangle2 box = get_nodes_bounding_rectangle(nset.begin(), nset.end());

        for (list<irectangle2>::const_iterator riter = gtrs.begin(); riter != gtrs.end(); ++riter) {
            double r = (double)riter->intersection(box).area()/(double)box.area();

            if (r > stat) stat = r;
        }
        if (nd->m >= (int)statmap.size()) statmap.resize(nd->m + 1, statmap_t::value_type(0.0, 0));
        statmap[nd->m].first += stat;
        statmap[nd->m].second += 1;
    }
}

void o_learning::add_to_validation_set(layer1_result* res, const list<irectangle2>& gtrs)
{
    cout << "<R";
    cout << "><U";
    update_statmap(res, gtrs);
    cout << "><S";
    vset.push_back(vset_item_t(streamed_pointer(res), gtrs));
    cout << '>';
}

void o_learning::get_mstat(mstat_t& mstat)
{
    mstat.clear();

    if (statmap.size() < library->layer_size(srclayer))
        statmap.resize(library->layer_size(srclayer), statmap_t::value_type(0.0, 0));
    for (int i = 0; i < (int)statmap.size(); ++i) {
        mstat.push_back(mstat_t::value_type(div0(statmap[i].first, statmap[i].second), i));
    }
    sort(mstat.begin(), mstat.end(), greater<mstat_t::value_type>());
}


int maxvsize = 0;

// Removes all nodes in the list of nodes and their reconstruction sets
// which intersect (their reconstruction sets) 
void compress_node_list(list<pair<node*, set<ipoint2> > >& result, double thresh)
{
    typedef pair<node*, set<ipoint2> > result_item_t;


    list<result_item_t>::iterator iter = result.begin(); 
    
    while (iter != result.end()) {
        bool keep = true;

        for (list<result_item_t>::iterator liter = result.begin(); liter != iter; ++liter) {
            if ((double)intersection_size(iter->second, liter->second)/liter->second.size() >= thresh) {
                keep = false;
                break;
            }
        }
        if (keep) ++iter;
        else iter = result.erase(iter);
    }
}

void o_learning::candidate_positions(list<pair<node*, set<ipoint2> > >& v, olv_object_part& model, 
    const map<node*, int>& indexmap, layer1_result* res, int m)
{
    typedef pair<node*, set<ipoint2> > result_item_t;

    v.clear();
    for (int i = 0; i < model.size(); ++i) {
        dupmap_t::iterator fiter = dupmap.find(iipair(model.type(i), m));
        ipoint2 pos = model.pos(i);

        // if there is a duplet from siter-th part in the model
        if (fiter != dupmap.end()) {
            vector<node*> r;

            // match with this duplet's positions
            for (ip2_vector::iterator dupiter = fiter->second.begin(); dupiter != fiter->second.end(); ++dupiter) {
                vector<node*> candidates;

                res->schur_product_max(pos.x + dupiter->x, pos.y + dupiter->y,
                    dist, m, srclayer, &layer1_result::idspf, 0.1, candidates);
                if (!candidates.empty())
                    r.push_back(candidates.front());
            }
            for (vector<node*>::iterator riter = r.begin(); riter != r.end(); ++riter) {
                node* rn = *riter;
                set<node*> nset;
                set<ipoint2> pset;
                map<node*, int>::const_iterator fiter = indexmap.find(rn);
                map<node*, int>::const_iterator mfiter = indexmap.find(model.str.back());

                // Skip if the node was checked before...
                if (fiter != indexmap.end() && mfiter != indexmap.end() && fiter->second <= mfiter->second)
                    continue;

                (*riter)->get_neighbor_set(EdgeConnection::TO_LAYER0, nset);
                node_set_to_point_set(pset, nset.begin(), nset.end());
                if (intersection_size(model.rec, pset)/(double)pset.size() < 0.4 /* ��� */) {
                    v.push_back(result_item_t(rn, pset)); 
                }
            }
        }
    }

    compress_node_list(v, 0.6 /* ����� */);

    if (v.size() > maxvsize) maxvsize = v.size();
}

double max_intersection_percent(const list<set<ipoint2> >& pslist, const set<ipoint2>& ps)
{
    if (ps.empty()) return 0.0;

    double result = 0.0;

    for (list<set<ipoint2> >::const_iterator iter = pslist.begin(); iter != pslist.end(); ++iter) {
        double q = (double)intersection_size(*iter, ps)/ps.size();

        if (q > result) result = q;
    }
    return result;
}

void o_learning::augment_model(list<olv_object_part>& models, int& modelcount, olv_object_part model, const mstat_t& mstat, 
    const map<node*, int>& indexmap, layer1_result* res, int maxmodels)
{
    typedef list<pair<node*, set<ipoint2> > > position_list_t;

    if (model.size() >= max_cluster_n) 
        return;

    list<set<ipoint2> > recmem;
    int usedm = 0;
    int tailcount = 0;

    for (mstat_t::const_iterator miter = mstat.begin(); miter != mstat.end(); ++miter) {
        int m = miter->second;
        position_list_t pos;
        bool augmented = false;

        candidate_positions(pos, model, indexmap, res, m);

        for (position_list_t::iterator positer = pos.begin(); positer != pos.end(); ++positer) { 
            if (max_intersection_percent(recmem, positer->second) < intersection_threshold) {
                olv_object_part memmodel = model;

                model.augment(positer->first, positer->second);
                models.push_back(model);

                if (model.size() >= min_cluster_n) {
                    ++modelcount;
                    ++tailcount;
                }
                if (tailcount > maxmodels) {
                    return;
                }

                int memmc = modelcount;

                augment_model(models, modelcount, model, mstat, indexmap, res, maxmodels/4 + 2);

                tailcount += modelcount - memmc;

                if (tailcount > maxmodels) {
                    return;
                }
                model = memmodel;
                recmem.push_back(positer->second);

                augmented = true;
            }
        }
        if (augmented)
            ++usedm;

        // We consider only a few -- best -- types for adding parts
        if (usedm > 3 /* ����� */)  
            break;
    }
}

void o_learning::make_models(list<olv_object_part>& models, layer1_result* res)
{
    vector<node*> lnodes;
    map<node*, int> indexmap;

    cout << " (map update";
    if (dupmap.empty() && !dupstat.empty()) {
        make_duplets();
        dupstat.clear();
    }
    if (mstat.empty()) 
        get_mstat(mstat);
    cout << ")";

    // Prepare res and node set (vector)
    res->add_reconstruction_edges_fwd(srclayer);
    res->get_layer_nodes(lnodes, srclayer);

    //debug_check1(lnodes, statmap);
    sort(lnodes.begin(), lnodes.end(), mstat_compare(&statmap));

    // Make node index map
    for (vector<node*>::iterator niter = lnodes.begin(); niter != lnodes.end(); ++niter) {
        int index = (int)(niter - lnodes.begin());
        indexmap.insert(pair<node*, int>(*niter, index));
    }
    
    // Start making models 
    int used = 0;
    int modelcount = 0;
    list<set<ipoint2> > recmem;

    for (vector<node*>::iterator niter = lnodes.begin(); niter != lnodes.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;

        set<node*> nset;
        set<ipoint2> pset;

        n->get_neighbor_set(EdgeConnection::TO_LAYER0, nset);
        node_set_to_point_set(pset, nset.begin(), nset.end());

        if (max_intersection_percent(recmem, pset) < 0.7 /* ����� */) {
            olv_object_part model;

            model.augment(n, pset);
            augment_model(models, modelcount, model, mstat, indexmap, res, max_models/4 /* ������ see 3 below */);
            recmem.push_back(pset);
            ++used;
        }

        // We consider only a few -- best -- types for initial model part
        if (used > 3 /* ������ */)
            break;
    }

    for (list<olv_object_part>::iterator iter = models.begin(); iter != models.end() ; ) {
        if (iter->size() < min_cluster_n || iter->size() > max_cluster_n) 
            iter = models.erase(iter);
        else
            ++iter;
    }
}

// result: vector items are pairs ((-1, type), path_map_t)
// res: layer1_result with toLayer0 edges added
// n: node
// size: check nodes in circle with diameter size/2 
void get_overlapping_parts(vector<pair<ipoint2, path_map_t> >& result, layer1_result* res, const scmap_t& scmap,
    node* n, int size, int thresh)
{
    typedef pair<ipoint2, path_map_t> result_item_t;

    int name = EdgeConnection::TO_LAYER0;
    vector<node*> nbhood;
    set<node*> nset;
    set<ipoint2> psetn;
    map<int, path_map_stat> resmap;

    n->get_neighbor_set(name, nset);
    node_set_to_point_set(psetn, nset.begin(), nset.end());
    result.clear();

    res->get_neighbors_circular(nbhood, n, size/2, size/2, 0, true);
    nbhood.insert(nbhood.begin(), n);
    for (vector<node*>::iterator iter = nbhood.begin(); iter != nbhood.end(); ++iter) {
        node* nn = *iter;
        set<ipoint2> psetnn;
        
        nset.clear();
        nn->get_neighbor_set(name, nset);
        node_set_to_point_set(psetnn, nset.begin(), nset.end());
        if (intersection_size(psetn, psetnn) >= psetn.size() * thresh) {
            int m = node_type(nn);
            path_map_t pm;
            
            get_path_map(pm, res, scmap, nn, true);
            resmap[m].add(pm);
        }
    }
    for (map<int, path_map_stat>::iterator rmiter = resmap.begin(); rmiter != resmap.end(); ++rmiter) {
        if (rmiter->first == node_type(n)) 
            result.insert(result.begin(), result_item_t(ipoint2(-1, rmiter->first), rmiter->second.get()));
        else
            result.push_back(result_item_t(ipoint2(-1, rmiter->first), rmiter->second.get()));
    }
}

// Note that list of models is changed; models that are not added to the library are removed
void o_learning::add_to_library(part_lib* plb, list<olv_object_part>& models, layer1_result* res,
    const scmap_t& scmap)
{
    matrix<double> mask;
    list<olv_object_part>::iterator miter = models.begin();

    gaussian_mask(gaussian_dim, gaussian_dim, gaussian_sigma, mask);
    while (miter != models.end()) {
        const vector<node*>& str = miter->str;
        vector<part_lib::cluster_data_t> objdata;

        for (vector<node*>::const_iterator siter = str.begin(); siter != str.end(); ++siter) {
            node* n = *siter;
            layer1_data* nd = (layer1_data*)n->data;
            vector<pair<ipoint2, path_map_t> > simnodes;

            get_overlapping_parts(simnodes, res, scmap, n, 3, cluster_member_threshold);  // 3 !!!!!!

            if ((int)simnodes.size() > cluster_size) 
                simnodes.resize(cluster_size);
            if (simnodes.empty()) {
                cout << "simnodes empty" << endl;
                throw exception();
            }

            part_lib::cluster_data_t cluster;

            for (vector<pair<ipoint2, path_map_t> >::iterator sniter = simnodes.begin(); sniter != simnodes.end(); ++sniter) {
                cluster.types.push_back(sniter->first);
                cluster.geo.insert(pair<int, path_map_t>(sniter->first.y, sniter->second));
            }
            cluster.pos.set(nd->x, nd->y);

            objdata.push_back(cluster);
        }
        vector<matrix<double>*> maskv(objdata.size(), &mask);
        if (plb->find_equivalent_object_part(srclayer + 1, objdata, mask.width) < 0) {
            plb->add_object_part(srclayer + 2, new part_data(), objdata, maskv, contraction);
            ++miter;
        } else 
            miter = models.erase(miter);
        
    }
}

template<class T> T& at(vector<T>& v, int i, const T& defval)
{
    if (i >= (int)v.size()) v.resize(i + 1, defval);
    return v.at(i);
}


// hits[m] is number of hits of part m
// misses[m] is number of misses of part m
// rhits[m][i][j] is a number of hits of part m to j-th rectangle in i-th image
void count_hits(vector<int>& hits, vector<int>& misses, vectorn<int>::type3& rhits, layer1_result* res, 
    int resi, int layer, const list<irectangle2>& reclist, double thresh, int maxtypes)
{
    list<layer1_result::box_data_t> boxes;
    response_filter rfilter;

    res->get_boxes(boxes, nullptr, layer, RR_RESPONSE, false, false, rfilter, 1.0, 100);
    keep_best_boxes(boxes, maxtypes);
    for (list<layer1_result::box_data_t>::iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
        irectangle2& r = biter->box;
        int hit = -1;
        int rcount = 0;

        for (list<irectangle2>::const_iterator iter = reclist.begin(); iter != reclist.end(); ++iter) {
            if ((double)r.intersection(*iter).area()/r.union_area(*iter) >= thresh) hit = rcount;
            ++rcount;
        }

        if (hit < 0) ++at(misses, biter->m, 0); 
        else {
            ++at(hits, biter->m, 0); 
            ++at(rhits, biter->m, resi, hit, 0);
        }
    }
}

// tokeep: parts to keep
// library: library to perform validation on
// rhits[i][j] counts how many times the j-th rectangle of i-th validation image was hit.
// libsize0: old size of library 
// libsize1: new size of library
void o_learning::validate(vector<int>& tokeep, vectorn<int>::type2& rhits, part_lib* library, 
    int libsize0, int libsize1)
{
    // some more changes?
    cout << "realization_ratio_threshold: " << inference_cfg.getValueDouble("realization_ratio_threshold", 0.0) << endl;

    layern_creator creator(inference_cfg);
    vector<int> hits, misses;
    int boxcount = 0;
    int vrescount = 0;
    vectorn<int>::type3 rhitsl;

    creator.set_library(library);
    creator.allowed_parts = set_range(libsize0, libsize1 - 1);
    tokeep.clear();
    for (vset_t::iterator vsiter = vset.begin(); vsiter != vset.end(); ++vsiter) {
        bool skip = true;

        for (int j = 0; j < (int)vsiter->second.size(); ++j) {
            if (at(rhits, vrescount, j, 0) < redundancy_threshold) {
                skip = false;
                break;
            }
        }
        if (!skip) {
            layer1_result* res = (layer1_result*)vsiter->first.get();
            vector<int> rh;
            
            creator.add_layer(res, srclayer + 2, 0); 
            if (res->max_layer_index() >= srclayer + 1) {
                cout << "[" << res->shape_nodes[srclayer + 1].size() << "]";
                count_hits(hits, misses, rhitsl, res, vrescount, srclayer + 1, vsiter->second, 0.6, 
                    type_bite_threshold); // !!!!!
            }
            boxcount += (int)vsiter->second.size();

            delete res;
        }
        ++vrescount;

        cout << '.';
    }
    if (hits.size() > misses.size()) misses.resize(hits.size(), 0);
    else hits.resize(misses.size(), 0);

    // count "new hits"
    vector<int> newhits(hits.size(), 0);

    foreach_vector3(rhitsl, i, j, k) {
        int lc = rhitsl[i][j][k];
        
        if (lc > 0 && at(rhits, j, k, 0) < redundancy_threshold) {
            ++newhits[i];
            //else { cout << "red " << i << " "; }
        }
    }

    cout << endl;
    cout << "hits: " << hits << endl;
    cout << "misses: " << misses << endl;
    cout << "newhits: " << newhits << endl;
    cout << "boxcount: " << boxcount << endl;
    cout << "hit-thresh: " << hit_threshold*boxcount << endl;

    for (int i = 0; i < (int)hits.size(); ++i) {
        if (newhits[i] > 0 && 
                (double)hits[i]/(hits[i] + misses[i]) >= hit_ratio_threshold && 
                hits[i] > hit_threshold*boxcount) {
            foreach23_vector3(rhitsl, i, rhlj, rhlk) {
                ++at(rhits, rhlj, rhlk, 0);
            }
            tokeep.push_back(i - libsize0);
        }
    }
    cout << "tokeep: " << tokeep << endl;
}

template<class T> void delete_items(list<T>& l, vector<int> indices)
{
    int li = 0, vi = 0;
    typename list<T>::iterator liter = l.begin();

    sort(indices.begin(), indices.end());
    for (int i = 0; i < (int)indices.size(); ++i) {
        while (liter != l.end() && li < indices[i]) { 
            li++; ++liter; 
        }
        if (liter != l.end() && li == indices[i]) {
            liter = l.erase(liter);
            ++li;
        }
    }
}

template<class T> void keep_items(list<T>& l, vector<int> indices)
{
    int li = 0, vi = 0;
    typename list<T>::iterator liter = l.begin();

    sort(indices.begin(), indices.end());
    for (int i = 0; i < (int)indices.size(); ++i) {
        while (liter != l.end() && li < indices[i]) { 
            liter = l.erase(liter);
            li++;  
        }
        if (liter != l.end() && li == indices[i]) {
            ++liter;
            ++li;
        }
    }
    while (liter != l.end()) 
        liter = l.erase(liter);
}


void o_learning::learn_models(layer1_result* res, const scmap_t& scmap)
{
    list<olv_object_part> models;

    // Make models
    cout << "<MM";
    make_models(models, res);
    cout << '(' << models.size() << ")>";

    // Add models to library
    part_lib* tmplib = (part_lib*)srclibrary->get_copy_s();
    int libsize0 = tmplib->layer_size(srclayer + 1);

    // Note that res is needed to extract "equivalent" parts
    cout << "<ATL";
    res->delete_edges(EdgeConnection::TO_LAYER0); // link_path <-> add_reconstruction_edges_fwd clash
    add_to_library(tmplib, models, res, scmap); 

    int libsize1 = tmplib->layer_size(srclayer + 1);

    cout << "(" << libsize1 - libsize0 << ")>";

    // Nothing to do if no parts produced
    if (libsize1 == libsize0)
        return;
    
    // Validate library; infer + count
    vector<int> tokeep;
    
    cout << "<V";
    validate(tokeep, rhits, tmplib, libsize0, libsize1);
    cout << "(" << tokeep.size() << ")";
    cout << "><ATL2";
    cout << "M=" << models.size();

    keep_items(models, tokeep);
    //cout << "after keeping of " << tokeep.size() << " parts " << models.size() << endl;
    add_to_library(library, models, res, scmap);
    cout << ">(";

    cout << library->parts[srclayer + 1].size() << " libsize)";
}

void o_learning::finalize()
{
    vector<int> newparts;

    for (int i = srclibrary->layer_size(srclayer + 1); i < library->layer_size(srclayer + 1); ++i)
        newparts.push_back(i);

    library->add_category(srclayer + 3, newparts, catname);
}
