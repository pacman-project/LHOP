
#include <cv.h>
#include <queue>
#include "layer_learning.h"
#include "../graphs/graph_utils.h"

// global functions
///////////////////////////////////////////////////////////////////////////////

bool get_region_sym(matrix<double>& result, matrix<double>& m, int i0, int i1, int j0, int j1)
{
    int xcenter = (i0 + i1)/2, ycenter = (j0 + j1)/2;
    if (xcenter < 0 || xcenter >= (int)m.width || ycenter < 0 || ycenter >= (int)m.height)
        return false;

    result.resize(i1 - i0, j1 - j0);

    int i, j, truei, truej, di, dj;

    for (i = i0, di = 0; i < i1; ++i, ++di) {
        truei = (i < 0 || i >= (int)m.width) ? i1 - di - 1 : i;
        for (j = j0, dj = 0; j < j1; ++j, ++dj) {
            truej = (j < 0 || j >= (int)m.height) ? j1 - dj - 1 : j;
            result(di, dj) = m(truei, truej);
        }
    }
    return true;
}

void max_normalize(matrix<double>& m, double factor = 1.0)
{
    double max = m.maximum();

    if (max == 0.0) return;
    for_each_element (m, i) {
        m[i] /= max*factor;
    }
}

void convolve_int_matrix(matrix<int>& result, matrix<int>& im, const img& dm)
{
    result = im;

    if (dm.width != dm.height) return;

    int ddim = (int)dm.width;
    int ddim2 = (int)dm.width/2;
    int dd = ddim*ddim;

    int** sptr = new int*[dd];
    int i, j, k;

    k = 0;
    for (j = 0; j < ddim; ++j) 
        for (i = 0; i < ddim; ++i) 
            sptr[k++] = &im(i, j);

    int* dptr = &result(ddim2, ddim2);
    int mwidth = (int)im.width - ddim2;
    int mheight = (int)im.height - ddim2;
    HOP_REAL conv;

    for (j = ddim2; j < mheight; ++j) {
        for (i = ddim2; i < mwidth; ++i) {
            conv = 0.0;
            for (k = 0; k < dd; ++k) conv += dm[k] * (*sptr[k]);
            *dptr = (int)conv;
            for (k = 0; k < dd; ++k) ++(sptr[k]);
            ++dptr;
        }
         for (k = 0; k < dd; ++k) sptr[k] += 2*ddim2;
         dptr += 2*ddim2;
    }
    delete sptr;
}

void get_reconstruction_map(map<node*, set<node*> >& result, layer1_result* res, int layer)
{
    typedef map<node*, set<node*> > result_t;

    int prevname = atom("toPrevLayer");

    result.clear();
    for (list<node*>::iterator iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
        node* n = *iter;
        result_t::iterator riter = result.insert(result_t::value_type(n, result_t::mapped_type())).first;

        res->recurse_from_node(n, prevname, riter->second);
    }
}

void get_positions(set<ipoint2>& result, const set<node*>& s, int z)
{
    result.clear();
    for (set<node*>::const_iterator siter = s.begin(); siter != s.end(); ++siter) {
        const img_node_data* nd = (img_node_data*)(*siter)->data;

        if (nd->z == z) result.insert(ipoint2(nd->x, nd->y));        
    }
}

void node_set_to_region_set(set<ipoint2>& result, const set<node*>& nset, const vector<set<ipoint2> >& rvector, int layer)
{
    for (set<node*>::const_iterator iter = nset.begin(); iter != nset.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;
        const set<ipoint2>& rs = rvector[nd->m];

        if (nd->z == layer) {
            for (set<ipoint2>::const_iterator riter = rs.begin(); riter != rs.end(); ++riter) {
                result.insert(ipoint2(nd->x + riter->x, nd->y + riter->y));
            }
        }
    }
}

// map_learning
///////////////////////////////////////////////////////////////////////////////

// auxilliary functions

void get_region_map(map<node*, set<int> >& rmap, layer1_result* res, int layer,
    const vector<matrix<bool> >& regions, const vector<ipoint2>& rcenters)
{
    typedef map<node*, set<int> > result_t;

    int to_prev = atom("toPrevLayer");

    rmap.clear();
    for (list<node*>::iterator iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
        node* n = *iter;
        set<node*> proj;
        set<int>& rset = rmap[n];

        res->recurse_from_node(n, to_prev, proj);
        for (set<node*>::iterator siter = proj.begin(); siter != proj.end(); ++siter) {
            layer1_data* pd = (layer1_data*)(*siter)->data;

            if (pd->z == 0) {
                const ipoint2& rc = rcenters[pd->m];

                get_region_set(pd->x, pd->y, res->x_size(0), res->y_size(0), 
                    regions[pd->m], rc.x, rc.y, rset);
            }
        }
    }
}

void points_from_edge_map(ip2_vector& v, const path_map_t& m)
{
    v.clear();
    v.reserve(m.size());
    for (path_map_t::const_iterator iter = m.begin(); iter != m.end(); ++iter) {
        v.push_back(iter->second.p);
    }
}


void centers_from_geometry(ip2_vector& v, const part_learning::geometry_t& g)
{
    ip2_vector mv;

    for (part_learning::geometry_t::const_iterator iter = g.begin(); iter != g.end(); ++iter) {
        points_from_edge_map(mv, *iter);
        v.push_back(ipoint2::center(mv.begin(), mv.end()));
    }
}



map_learning::map_learning() :
    stat(),
    regions(),
    rcenters()
{
}

map_learning::map_learning(const config_dictionary& cfg) :
    stat(),
    regions(),
    rcenters()
{
    cfg_init(cfg);
}

map_learning::~map_learning()
{
}

void map_learning::cfg_init(const config_dictionary& cfg)
{
    if (cfg.is_defined("source_layer_index")) {
        source_layer = cfg.get_value_int("source_layer_index", 0) - 1;  // backward compatibility!
    } else {
        cfg.get_value(source_layer, "source_layer", true);
    }
    cfg.get_value(stat_dim, "nb_size", true);

    // map update
    center_val_threshold = cfg.get_value_double("center_val_threshold", 0.0);
    center_val_threshold_rel = cfg.get_value_double("center_val_threshold_rel", 0.5);
    nb_val_threshold_rel = cfg.get_value_double("nb_val_threshold_rel", 0.6);
    nbthresh_min = cfg.get_value_double("nbthresh_min", 0.0);
    nbthresh_max = cfg.get_value_double("nbthresh_max", 0.0);
    seq_min_intersection_percent = cfg.get_value_double("seq_min_intersection_percent", 0.0);
    seq_max_intersection_percent = cfg.get_value_double("seq_max_intersection_percent", 0.5);

    // finding maxima
    max_max = cfg.get_value_int("max_max", 4);
    max_val_threshold = cfg.get_value_double("max_val_threshold", 0.01);
    if (!cfg.is_defined("max_val_threshold"))
        max_val_threshold = cfg.get_value_double("min_update_count_percent", 0.01);
    individual_max = cfg.get_value_bool("individual_max", false);
    max_sigma = cfg.get_value_double("max_sigma", 0.0);
    max_nbhood_mask = cfg.get_value_int("max_nbhood_mask", 5);
    max_radius = cfg.get_value_int("max_radius", 2); 
}

void map_learning::reset(part_lib* library)
{
    library->get_regions(1, regions, rcenters);
    stat.clear();
}

void map_learning::dispose()
{
    stat.clear();
}

// Prepares res for update. Sets ALLOW_UPDATE_ATTR to all nodes.
void map_learning::prepare_for_update(layer1_result* res)
{
    if (res == nullptr) 
        return;

    res->set_attr(res->shape_nodes[source_layer].begin(),
        res->shape_nodes[source_layer].end(), ALLOW_UPDATE_ATTR);
    //res->connect_neighbors_circular(res->shape_nodes[source_layer], stat_dim/2, stat_dim/2, 
    //    ALLOW_UPDATE_ATTR, atom("toNeighbor"));
}

void map_learning::update(layer1_result* res)
{
    typedef map<node*, set<int> > rmap_t;

    if (res == nullptr)
        return;

    if (regions.empty()) 
        throw new_libhop_exception("map_learning::update: Map learner not initialized with library");

    vector<node*>& s_nodes = res->shape_nodes[source_layer];
    rmap_t rmap;

    get_region_map(rmap, res, source_layer, regions, rcenters);
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* cn = *iter;

        do {
            layer1_data* cd = (layer1_data*)cn->data;
            ipoint2 cp(cd->x, cd->y);
        
            if (cn->is_attr_set(ALLOW_UPDATE_ATTR) && cd->vval() >= center_val_threshold) {
                double nbt_min, nbt_max;

                if (nbthresh_min >= nbthresh_max) {
                    nbt_min = cd->vval() * nb_val_threshold_rel;
                    nbt_max = cd->vval() * (2 - nb_val_threshold_rel);
                } else {
                    nbt_min = nbthresh_min;
                    nbt_max = nbthresh_max;
                }
                set<int>& cset = rmap[cn];
                int xoff = stat_dim/2 - cd->x;
                int yoff = stat_dim/2 - cd->y;
                vector<node*> nbhood;

                res->get_nodes_circular(nbhood, cn, stat_dim/2, stat_dim/2, ALLOW_UPDATE_ATTR);

                // Update map with neighbors
                for (vector<node*>::iterator niter = nbhood.begin(); niter != nbhood.end(); ++niter) {
                    node* nb = *niter;
                    layer1_data* nbd = (layer1_data*)nb->data;
                    double nbval = nbd->vval();
                    set<int>& nbset = rmap[nb];
                    int is = intersection_size(cset, nbset);
                    int mins = min<int>((int)nbset.size(), (int)cset.size());

                    if ((is >= seq_min_intersection_percent * mins && is <= seq_max_intersection_percent * mins) &&
                            (nbt_min <= nbval && nbval <= nbt_max) &&
                            cp.distance2(nbd->x, nbd->y) > sqr(max_radius + 0.5)) {
                        matrix<double>& im = stat[ipoint2(cd->m, nbd->m)];
                        // Uncomment if update is distance dependent
                        //int x = nbd->x - cd->x, y = nbd->y - cd->y;
                        //int dist2 = x*x + y*y;

                        if (im.empty()) im.resize(stat_dim, stat_dim, 0.0);
                        im(xoff + nbd->x, yoff + nbd->y) += nbval; // * ::pow(dist2, -0.125)
                    }
                }

            } // end of "if (center->is_attr_set(ALLOW_UPDATE_ATTR) && cdata->val >= center_val_threshold)"

            cn = cd->next;
            break; // Only the best!
        } while (cn != nullptr);
    }
        
}

// Save statistics to files
// file_template should contain two %d format characters
void map_learning::display_statistics(const string& file_template, double blur_sigma /* = 0.75 */) const
{
    char fname[1000];
    int blur_dim = 4*blur_sigma + 1;
    img* blur_mask = img::gaussian_mask(blur_dim, blur_dim, blur_sigma);

    for (stat_map_t::const_iterator iter = stat.begin(); iter != stat.end(); ++iter) {
        matrix<double> tmp(iter->second);
        img im;

        tmp.add_border(blur_dim/2, blur_dim/2, 0.0);
        convolve_matrix(im, tmp, *blur_mask);

        sprintf(fname, file_template.c_str(), iter->first.x, iter->first.y);
        im.save_jet_colormap(fname, -1000);
    }
    delete blur_mask;
}

// Get vector of maxima
// stat: matrix of statistics
// max_maxima: max number of maxima
// max_thresh: each max must be >= max_thresh; 
void map_learning::get_maxima(max_vector_t& maxima, const matrix<double>& stat, double max_thresh) const
{
    typedef pair<double, ipoint2> queue_item_t;
    typedef priority_queue<queue_item_t> queue_t;

    const int border = 5;

    static img blurmask = img::gaussian_mask2(border, border, 0.75);

    maxima.clear();
    matrix<double> m(stat);
    matrix<double> dist;
    queue_t queue;
    
    m.add_border(border, border, 0.0);
    m = convolve_matrix(m, blurmask);
    m.change_values_leq(max_thresh, 0.0);

    if (max_sigma > 0.0)
        dist = img::gaussian_mask2(max_nbhood_mask, max_nbhood_mask, max_sigma);

    for_each_xy_int (m, i, j) {
        double lmax = m.maximum_c((int)i, (int)j, max_radius, max_radius);

        if (lmax > 1E-6 && m(i, j) == lmax) 
            queue.push(queue_item_t(lmax, ipoint2(i, j)));
    }

    matrix<int> bm(m.width, m.height, 0);

    //cout << queue.size();
    //if (!queue.empty()) 
    //    cout << "; " << queue.top();
    //cout << endl;
    bm.set_region_c(m.width/2, m.height/2, max_radius, max_radius, 1);
    while (!queue.empty() && (int)maxima.size() < max_max) {
        ipoint2 p = queue.top().second;

        //cout << "--" << p << endl;
        if (bm(p.x, p.y) == 0) {
            //cout << "***" << endl;
            bm.set_region_c(p.x, p.y, max_radius, max_radius, 1);
            maxima.push_back(map_learning::max_item_t());
            maxima.back().pos = ipoint2(p.x - m.width/2, p.y - m.height/2);
            if (!dist.empty())
                maxima.back().dist = dist;
            else {
                int radius2 = max_nbhood_mask/2;

                get_region_sym(maxima.back().dist, m, p.x - radius2, p.x + radius2 + 1,
                    p.y - radius2, p.y + radius2 + 1);
                max_normalize(maxima.back().dist);
            }
            //cout << ipoint2(x, y) << ' ';
        }
        queue.pop();
    }
}

void map_learning::get_maxima(max_map_t& result) const
{
    double max = 0.0;

    result.clear();
    if (!individual_max) {
        for (stat_map_t::const_iterator iter = stat.begin(); iter != stat.end(); ++iter) {
            const stat_t& m = iter->second;
            int r = m.maximum();

            if (r > max) max = r;
        }
        max *= max_val_threshold;
    }
    for (stat_map_t::const_iterator iter = stat.begin(); iter != stat.end(); ++iter) {
        const stat_t& m = iter->second;
        
        if (individual_max)
            max = max_val_threshold * m.maximum();
        //cout << "coordinates: " << iter->first << endl;
        get_maxima(result[iter->first], m, max);
    }
}

void map_learning::read_maxima_from_stream(max_map_t& maxima, vector<matrix<bool> >& regions, vector<ipoint2>& rcenters,
        const string& fname)
{
    ifstreamer is;
    
    is.open(fname);

    int maxsize;

    is.read(maxsize);
    for (int m = 0; m < maxsize; ++m) {
        ipoint2 key;
        max_vector_t& mapped = maxima[key];
        int vsize;

        is.read(key);
        is.read(vsize);
        mapped.reserve(vsize);
        for (int v = 0; v < vsize; ++v) {
            max_item_t item;

            is.read(item.pos);
            is.read(item.dist);
            mapped.push_back(item);
        }
    }
    is.read(regions);
    is.read(rcenters);

    is.close();
}

void map_learning::write_maxima_to_stream(const string& fname)
{
    ofstreamer os;
    max_map_t maxima;
    
    get_maxima(maxima);

    os.open(fname);

    os.write((int)maxima.size());
    for (max_map_t::iterator iter = maxima.begin(); iter != maxima.end(); ++iter) {
        os.write(iter->first);
        os.write((int)iter->second.size());
        for (max_vector_t::iterator viter = iter->second.begin(); viter != iter->second.end(); ++viter) {
            os.write(viter->pos);
            os.write(viter->dist);
        }
    }
    os.write(regions);
    os.write(rcenters);

    os.close();
}


// part_learning 
///////////////////////////////////////////////////////////////////////////////

struct position_compare_less {
    bool operator()(const ipoint2& p1, const ipoint2& p2)
    {
        if (abs(p1.x - p2.x) < 2*abs(p1.y - p2.y)) return p1.y < p2.y; else return p1.x < p2.x;
    }
};

part_learning::part_learning() :
    stat(), maxima(), layer_contraction(-1.0)
{ 
}

part_learning::part_learning(const config_dictionary& cfg) :
    stat(), maxima(), layer_contraction(-1.0)
{ 
    init_cfg(cfg);
}

part_learning::~part_learning()
{
}

void part_learning::init_cfg(const config_dictionary& cfg)
{
    if (cfg.is_defined("source_layer_index")) {
        source_layer = cfg.get_value_int("source_layer_index", 0) - 1;  // backward compatibility!
    } else {
        cfg.get_value(source_layer, "source_layer", true);
    }

    // part update
    center_val_threshold = cfg.get_value_double("center_val_threshold", -1.0);
    center_val_threshold_rel = cfg.get_value_double("center_val_threshold_rel", 0.2);
    if (cfg.is_defined("seq_max_intersection_percent2"))
        seq_min_intersection_percent = cfg.get_value_double("seq_min_intersection_percent2", 0.0);
    else
        seq_min_intersection_percent = cfg.get_value_double("seq_min_intersection_percent", 0.0);
    if (cfg.is_defined("seq_max_intersection_percent2"))
        seq_max_intersection_percent = cfg.get_value_double("seq_max_intersection_percent2", 0.5);
    else
        seq_max_intersection_percent = cfg.get_value_double("seq_max_intersection_percent", 0.5);
    max_candidates = cfg.get_value_int("max_candidates", 5);
    min_seq_size = cfg.get_value_int("min_part_length", 2);
    max_seq_size = cfg.get_value_int("max_part_length", 3);
    max_stat_size = cfg.get_value_int("max_stat_size", INT_MAX);

    // similarity calculation
    similarity_type = cfg.get_value_int("similarity_type", -1);

    // similarity calculation inhibition 
    max_matching_size = cfg.get_value_int("max_matching_size", 30);
}

void part_learning::reset(const map_learning& mlearner)
{
    mlearner.get_maxima(maxima);
    regions = mlearner.regions;
    rcenters = mlearner.rcenters;
}

void part_learning::reset(const string& fname)
{
    map_learning::read_maxima_from_stream(maxima, regions, rcenters, fname);    
}

void part_learning::reset()
{
    stat.clear();
}

void part_learning::dispose()
{
    stat.clear();
    maxima.clear();
    regions.clear();
    rcenters.clear();
}

// Update statistics. 'res' is inference graph and 'sc' is "shape context".
void part_learning::update(layer1_result* res, const scmap_t& scmap)
{
    typedef map<node*, set<int> > rmap_t;

    if (res == nullptr) 
        return;
    if (source_layer < 0 || source_layer > res->max_layer_index())
        throw new_libhop_exception("part_learning::update: Invalid source_layer."); 
    if (source_layer == 0) layer_contraction = 1.0; 
    else {
        if (layer_contraction < 0)
            layer_contraction = round((double)res->x_size(source_layer - 1)/res->x_size(source_layer), 1);
        else if (layer_contraction != round((double)res->x_size(source_layer - 1)/res->x_size(source_layer), 1))
            throw new_libhop_exception("part_learning::update: Contraction exception.");
    }

    bool rel_cvt = false; // "quantile" center_val_threshold from center_val_threshold or not...

    if (center_val_threshold < 0.0) {
        center_val_threshold = res->get_q_quantile_val(source_layer, 1 - center_val_threshold_rel);
        rel_cvt = true;
    }

    cout << " (cvt " << center_val_threshold << ")" << flush;

    vector<node*>& s_nodes = res->shape_nodes[source_layer];
    rmap_t rmap;

    get_region_map(rmap, res, source_layer, regions, rcenters);
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        while (n != nullptr && !n->is_attr_set(HAS_NEXT_LAYER)) {
            layer1_data* nd = (layer1_data*)n->data;

            if (nd->vval() >= center_val_threshold) {
                sym_part_t seq;
                real_part_t rseq;
                set<int> fr = rmap[n];

                // 1st element is center
                seq.push_back(ipoint2(nd->m, -1));
                rseq.val += nd->vval();
                rseq.nodes.push_back(n);
                extend_sequence(res, scmap, rmap, n, seq, rseq, fr);
            }   
            //break; // just take the best one (maybe 90% of the best one)
            n = nd->next;
        }
    }

    if (rel_cvt) 
        center_val_threshold = -1.0;
}

// Get list of candidates: triples {double v, ipoint2 p, node n} 
//   (response value, symbollic coord., corresp. node)
// 'res' is inference graph
// 'rmap' is reconstruction map for 'res'
// 'cn' is center node
// 'fr' is support of current sequence
void part_learning::get_best_candidates(candidates_t& bestc, layer1_result* res, 
    const map<node*, set<int> >& rmap, node* cn, const set<int>& fr)
{
    typedef map_learning::max_map_t map_t;

    layer1_data* cd = (layer1_data*)cn->data;
    double thresh = 0.5*cd->vval();

    for (map_t::iterator miter = maxima.lower_bound(ipoint2(cd->m, 0)); miter != maxima.end() && miter->first.x == cd->m; ++miter) {
        int vindex = 0;

        for (map_t::mapped_type::iterator viter = miter->second.begin(); viter != miter->second.end(); ++viter) {
            pair<node*, double> conv;

            conv = res->convolve_max_all(cd->x + viter->pos.x, cd->y + viter->pos.y, viter->dist, miter->first.y, 
                center_val_threshold*0.5, source_layer);
            if (conv.first != nullptr && conv.second > thresh) {
                const set<int>& br = rmap.find(conv.first)->second;
                int is = intersection_size(br, fr);

                if (is >= seq_min_intersection_percent*br.size() && is <= seq_max_intersection_percent*br.size()) {
                    bestc.push_back(candidates_item_t(conv.second, ipoint2(miter->first.y, vindex), conv.first));
                    bestc.sort(greater<candidates_item_t>());
                    if (bestc.size() > max_candidates)
                        bestc.resize(max_candidates);
                }
                
            }
            ++vindex;
        }
    }
}

// Extends part sequence ('seq') as much as possible
// 'res': inference graph
// 'scmap': "shape context" map
// 'rmap': reconstruction map (from node to its support on layer 0)
// 'cn': center node
// 'seq': current sequence
// 'rseq' is "realization" sequence (value of the realization and vector of nodes 
//    corresponding to each subpart)
// 'fr': current part support
void part_learning::extend_sequence(layer1_result* res, const scmap_t& scmap, 
    const map<node*, set<int> >& rmap, 
    node* cn, const sym_part_t& seq, const real_part_t& rseq, const set<int>& fr)
{
    candidates_t bestc;

    get_best_candidates(bestc, res, rmap, cn, fr);

    if (bestc.empty() || seq.size() > max_seq_size) {
        if (seq.size() >= min_seq_size && seq.size() <= max_seq_size) {
            update_sequence(seq, rseq, scmap, cn);
        }
    } else {
        double thresh = 0.9*bestc.front().val;

        for (candidates_t::iterator citer = bestc.begin(); citer != bestc.end() && citer->val >= thresh; ++citer) {
            set<int> fr1 = fr;
            sym_part_t seq1 = seq;
            real_part_t rseq1 = rseq;
            const set<int>& nr = rmap.find(citer->n)->second;

            fr1.insert(nr.begin(), nr.end());
            seq1.push_back(citer->coo);
            rseq1.val *= citer->val;
            rseq1.nodes.push_back(citer->n);
            extend_sequence(res, scmap, rmap, cn, seq1, rseq1, fr1);
        }
    }
}

void part_learning::display_maxima(const string& file_template, double blur_sigma /* = 0.75 */) const
{
    typedef map_learning::max_map_t map_t;

    char fname[1000];
    int blur_dim = 4*blur_sigma + 1;
    img* blur_mask = img::gaussian_mask(blur_dim, blur_dim, blur_sigma);
    irectangle2 box;

    for (map_t::const_iterator iter = maxima.begin(); iter != maxima.end(); ++iter) {
        for (map_t::mapped_type::const_iterator viter = iter->second.begin(); viter != iter->second.end(); ++viter)
            box.eat(viter->pos);
    }

    int width = 2*max(box.ur.x, -box.ll.x) + 11;
    int height = 2*max(box.ur.y, -box.ll.y) + 11;

    for (map_t::const_iterator iter = maxima.begin(); iter != maxima.end(); ++iter) {
        matrix<double> tmp(width, height, 0.0);
        img im;

        for (map_t::mapped_type::const_iterator viter = iter->second.begin(); viter != iter->second.end(); ++viter) {
            tmp(viter->pos.x + width/2, viter->pos.y + height/2) = 1;
        }

        tmp.add_border(blur_dim/2, blur_dim/2, 0.0);
        convolve_matrix(im, tmp, *blur_mask);

        sprintf(fname, file_template.c_str(), iter->first.x, iter->first.y);
        im.save_jet_colormap(fname, -1000);
    }
    delete blur_mask;
}

// Extends map 'm' with rules path |-> hpoint_t = {p, h} where p is a coordinate
// of the endpoint of the path, h is a histogram at p
// 'scmap' is a shape context map
// 'cn' is center node
// 'n' is starting node
void part_learning::make_path_map(path_map_t& m, const scmap_t& scmap, node* cn, node* n)
{
    typedef pair<vector<int>, node*> list_item_t;
    typedef list<list_item_t> list_t;

    list_t l;
    ipoint2 cpt = node_coordinates(cn);

    l.push_back(list_item_t(vector<int>(), n));
    make_path_map_rec(l);
    for (list_t::iterator liter = l.begin(); liter != l.end(); ++liter) {
        vector<int> p = liter->first;
        layer1_data* nd = (layer1_data*)liter->second->data;
        ipoint2 pt(nd->x, nd->y);

        pair<path_map_t::iterator, bool> ibpair = m.insert(path_map_t::value_type(p, hpoint_t()));

        if (ibpair.second) { // only insert the first path !!??
            ibpair.first->second.p = pt - cpt;

            scmap_t::const_iterator sciter = scmap.find(pt);

            if (sciter == scmap.end()) 
                throw new_libhop_exception("Shape context map is not defined at specific point.");
            ibpair.first->second.h = sciter->second;
        }
    }

}

path_map_t part_learning::make_path_map(const scmap_t& scmap, node* cn, node* n)
{
    path_map_t pm;

    make_path_map(pm, scmap, cn, n);
    return pm;
}

// Support function for make_path_map; l contains a list of incomplete paths and their end nodes.
void part_learning::make_path_map_rec(list<pair<vector<int>, node*> >& l)
{
    typedef pair<vector<int>, node*> list_item_t;
    typedef list<list_item_t> list_t;

    int name = atom("toPrevLayer");

    list_t newl;

    for (list_t::iterator liter = l.begin(); liter != l.end(); ++liter) {
        node* n = liter->second;
        
        foreach_neighbor (n, name, niter) {
            node* nn = neighbor_node(niter);
            edge_data_name* nned = (edge_data_name*)neighbor_edge_data(niter);
            vector<int> v = liter->first;

            if (nned == nullptr) 
                throw new_libhop_exception("Edge data not found; use add_edge_names = true");

            v.push_back(nned->index);
            newl.push_back(list_item_t(v, nn));
        }
    }
    if (!newl.empty()) {
        l.clear();
        l.splice(l.begin(), newl);
        make_path_map_rec(l);
    }

}

// Find last node in "toPrevLayer path" with index 0 (indicating central node).
node* get_center_node(node* cn)
{
    int name = atom("toPrevLayer");
    node* n = cn;

    if (cn == nullptr) return nullptr;
    while (true) {
        node* n0 = nullptr;

        foreach_neighbor (n, name, niter) {
            edge_data_name* ed = (edge_data_name*)neighbor_edge_data(niter);

            if (ed == nullptr) throw new_libhop_exception("Edge data not found; use add_edge_names = true");
            if (ed->index == 0) {
                n0 = neighbor_node(niter);
                break;
            }
        }
        if (n0 == nullptr) return n;
        n = n0;
    }
}



// Update statistics 'part_learning::stat'
// 'seq' is sequence to add to the stat; vector of symbollic coordinates (i, j) where i is part and j is 
//    maxima position, the first enty is data for center where j = -1
// 'rseq' is a sequence of 'seq' realization (value of the realization and vector of nodes 
//    corresponding to each subpart)
// 'scmap' is shape context map for the updated image
// 'cn' is center node
void part_learning::update_sequence(const sym_part_t& seq, const real_part_t& rseq, const scmap_t& scmap, node* cn)
{
    vector<ipoint2> geoseq; // Geometric coordinates of seq - to sort seq and rseq

    for (int i = 1; i < (int)seq.size(); ++i) {
        map_learning::max_item_t& mitem = maxima[ipoint2(seq[0].x, seq[i].x)][seq[i].y];
        geoseq.push_back(mitem.pos);
    }

    vector<int> ord = ordering<ipoint2>(geoseq.begin(), geoseq.end());
    sym_part_t sseq = seq;
    real_part_t srseq = rseq;

    permute_range<ipoint2>(sseq.begin() + 1, ord);
    permute_range<node*>(srseq.nodes.begin() + 1, ord);

    stat_t::mapped_type& s = stat[sseq];

    s.count += srseq.val; ///srseq.nodes.size(); 

    //if (s.occ.size() > max_stat_size)
    //    return;

    //s.occ.push_back(geometry_t());

    //geometry_t& emp = s.occ.back();
    node* cn0 = get_center_node(cn);
    vector<ipoint2> pts;

    if (cn0 == nullptr || node_layer(cn0) != 0) 
        throw new_libhop_exception("Path to central node on layer 0 does not exist.");
    //emp.resize(srseq.nodes.size());
    for (int i = 0; i < (int)srseq.nodes.size(); ++i) {
        node* n = srseq.nodes[i];
        path_map_t pm;
        
        make_path_map(pm, scmap, cn0, n);
        s.geol.update(i, pm);
        if (similarity_type == PCA_SIMILARITY) {
            vector<ipoint2> gpts = get_path_map_points(pm);
            pts.insert(pts.end(), gpts.begin(), gpts.end());
        }
    }
    if (!pts.empty())
        s.pcal.update(pts);
}

double cx_hull_area(const vector<ipoint2>& pts)
{
    vector<cv::Point2f> cvpts, outpts;

    for (auto piter = pts.begin(); piter != pts.end(); ++piter) {
        cvpts.push_back(cv::Point2f(piter->x, piter->y));
    }
    cv::convexHull(cvpts, outpts);
    return cv::contourArea(outpts);
}


// Returns a sorted vector of pairs where pair.first is count in stat_item_t
// and pair.second is pointer to key values of stat.
void part_learning::sort_stat(vector<pair<double, const sym_part_t*> >& pairs, int sorting_type)
{
    typedef pair<double, const sym_part_t*> vector_item_t;

    pairs.clear();
    pairs.reserve(stat.size());

    if (sorting_type == 1) {
        cout << "Sorting type = 1" << endl;
        map<sym_part_t, double> areamap;

        for (stat_t::iterator siter = stat.begin(); siter != stat.end(); ++siter) {
            int gsize = (int)siter->first.size();
            geometry_t geo2;
            ip2_vector v;
            vector<dpoint2> dv;

            geo2.resize(gsize);
            for (int g = 0; g < gsize; ++g)
                siter->second.geol.get_mean(geo2[g], g);
            centers_from_geometry(v, geo2);
            pairs.push_back(vector_item_t(cx_hull_area(v), &(siter->first)));
        }
    } else {
        cout << "Sorting type = 0" << endl;
        for (stat_t::iterator iter = stat.begin(); iter != stat.end(); ++iter) {
            pairs.push_back(vector_item_t(iter->second.count, &(iter->first)));
        }
    }
    sort(pairs.begin(), pairs.end(), greater<vector_item_t>());
}

// Set geometry data 'g' to library point 'p';
void part_learning::update_geometry(node* p, const geometry_t& g)
{
    int name = atom("lyrSrc");
    edge_pair cpair = p->get_neighbor_pair(atom("lyrCenter"));
    
    if (cpair.first != nullptr) {
        part_data_2a* ced = (part_data_2a*)cpair.second;
        ced->set_geo(g[0]);
    }
    foreach_neighbor (p, name, piter) {
        part_data_2* ped = (part_data_2*)neighbor_edge_data(piter);
        ped->set_geo(g[ped->index]);
    }
}

// Set pca data to library node. It throws an exception if data of p is not 
// of type vs_part_data
void part_learning::update_pca(node* p, const pca_data& pcd)
{
    vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);

    if (pd == nullptr) {
        cout << "Trying to set pca data to non-\"vs_part_data\" node." << endl;
        throw;
    }
    pd->set_data(pcd);
}

// Auxiliary function in add_to_library - adds lyrSimRoot edges according to part similarity 
// i.e. according to the result returned by merge_stat.
// 'layer' library layer
// 'mresult' result returned by merge_stat
// 'indexmap' maps old part indices (before inserting to library) to new part indices (library indices)
// (indices can be different due to skipping...)
void update_simroot_edges(part_lib* library, int layer, int oldlibsize, 
    const vector<vector<part_learning::similarity_item> >& mresult, const map<int, int>& indexmap)
{
    int rootname = atom("lyrSimRoot");
    int simname = atom("lyrSimilar");
    vector<node*>& parts = library->parts[layer];

    for (int i = 0; i < (int)mresult.size(); ++i) {
        const vector<part_learning::similarity_item>& mi = mresult[i];
        map<int, int>::const_iterator miter = indexmap.find(mi[0].index);

        if (miter != indexmap.end()) {
            node* p0 = parts[miter->second];

            for (int j = 0; j < (int)mi.size(); ++j) {
                int pindex = indexmap.find(mi[j].index)->second;

                if (pindex >= oldlibsize) {
                    node* p = parts[pindex];

                    //if (mi[j].permutation.empty()) 
                    //    library->add_edge_2(p, p0, new edge_data_t<double>(mi[j].val), nullptr, rootname, simname);
                    //else 
                    
                    library->add_edge_2(p, p0, new part_data_sim(mi[j].val), nullptr, rootname, simname);
                }
            }
        }
    }
}

// Use only for backward compatibility!!!
void part_learning::add_to_library(part_lib* library, int part_max_number)
{
    add_to_library(library, part_max_number, 0, 1, 0.0, 0.0);
}

void part_learning::add_to_library(part_lib* library, int part_max_number, int sorting_type, int cluster_size, double dthresh, double scthresh)
{
    typedef pair<double, const sym_part_t*> vector_item_t;
    typedef map<sym_part_t, geometry_t> stat2_t;
    typedef map<sym_part_t, pca_data> stat3_t;

    matrix<double> dummy;
    vector<vector_item_t> parts;
    stat2_t stat2; // for SC_SIMILARITY
    stat3_t stat3; // for PCA_SIMILARITY
    vector<vector<similarity_item> > mresult;  // previously: pair<int, double>
    map<sym_part_t, vector<int> > permap;
    set<int> indset;
    int libsize = (int)library->parts[source_layer + 1].size();
    int count = 0;
    
    sort_stat(parts, sorting_type);

    if (similarity_type == -1) 
        similarity_type = library->is_vs_layer(source_layer) ? PCA_SIMILARITY : SC_SIMILARITY;

    cout << "Similarity type: " << similarity_type << endl;

    if (similarity_type == SC_SIMILARITY)
        merge_stat_sc(mresult, permap, stat2, part_max_number, sorting_type, cluster_size, part_max_number, dthresh, scthresh); 
    else if (similarity_type == PCA_SIMILARITY)
        merge_stat_pca(mresult, permap, stat3, part_max_number, sorting_type, dthresh);
    else throw new_libhop_exception("Unknown similarity type.");

    for (int i = 0; i < (int)mresult.size(); ++i) {
        for (int j = 0; j < (int)mresult[i].size(); ++j) {
            indset.insert(mresult[i][j].index);
        }
    }

    map<int, int> indexm;
	
    for (set<int>::iterator isiter = indset.begin(); isiter != indset.end(); ++isiter) { 
        //vector<vector_item_t>::iterator viter = parts.begin(); viter != parts.end() && count < part_max_number; ++viter) {
        const sym_part_t* spp = parts[*isiter].second;
        sym_part_t sppp = *spp;

        permute_range<ipoint2>(sppp.begin(), permap[*spp]);

        vector<iipair> str;
        vector<iipair> coo;
        vector<matrix<double>*> distr;
        int ctype = sppp.front().x;
        int pn;

        str.push_back(iipair(-1, ctype));
        coo.push_back(iipair(0, 0));
        distr.push_back(&dummy);
        for (sym_part_t::const_iterator piter = sppp.begin() + 1; piter != sppp.end(); ++piter) {
            int type = piter->x;
            int maxi = piter->y;
            map_learning::max_item_t& mitem = maxima[ipoint2(ctype, type)][maxi];

            str.push_back(iipair(-1, type));
            coo.push_back(iipair(mitem.pos.x, mitem.pos.y));
            distr.push_back(&mitem.dist);
        }

        part_data* pd = (similarity_type == PCA_SIMILARITY) ? new vs_part_data() : new part_data();

        pn = library->add_part(source_layer + 2, pd, str, coo, distr, -1, 
            layer_contraction, 0, 0); /// 1.0 = contraction!!!!
        if (pn < 0) {
            indexm.insert(pair<int, int>(*isiter, -1 - pn));
            delete pd;
        } else {
            indexm.insert(pair<int, int>(*isiter, pn));
            ++count;
            if (similarity_type == PCA_SIMILARITY) {
                update_pca(library->parts[source_layer + 1][pn], stat3[*spp]);
            } else { 
                geometry_t geo = stat2[*spp];

                permute_range<path_map_t>(geo.begin(), permap[*spp]);
                update_geometry(library->parts[source_layer + 1][pn], geo);
            }
        }
        cout << ' ' << pn;
    }
	cout << endl;
    update_simroot_edges(library, source_layer + 1, libsize, mresult, indexm);
}

void part_learning::similarity_matrix(matrix<double>& m, int max_part_number, double thresh)
{
    typedef pair<double, const sym_part_t*> vector_item_t;

    vector<vector_item_t> parts;

    cerr << "Similarity_matrix not implemented (debug function)." << endl;
    throw;

    /*
    //cout << "presort" << endl;
    sort_stat(parts);
    //cout << "postsort" << endl;
    max_part_number = min<int>((int)parts.size(), max_part_number);
    //cout << max_part_number << endl;
    m.resize(max_part_number, max_part_number, 0);
    for (int i = 0; i < max_part_number; ++i) {
        stat_t::iterator iter = stat.find(*parts[i].second);

        if (iter == stat.end()) throw new_libhop_exception("!!!!!!!!!!!!!");
        occurrences_t& iocc = iter->second.occ;

        //cout << "i = " << i << " ";

        // compare i(ter)-th part to j(ter)-th part; for each occurrence of i-th part 
        // check distance to each occurrence of j-th part; NNNN if all distances are <= thresh, then
        // say that these two parts are equivalent.
        for (int j = i; j < max_part_number; ++j) {
            stat_t::iterator jter = stat.find(*parts[j].second);
            if (jter == stat.end()) throw new_libhop_exception("????????????????");
            occurrences_t& jocc = jter->second.occ;
            double maxdist = 0.0;

            //cout << "j = " << j << " ";

            for (occurrences_t::iterator oiter = iocc.begin(); oiter != iocc.end(); ++oiter) {
                geometry_t& igeo = *oiter;

                for (occurrences_t::iterator ojter = jocc.begin(); ojter != jocc.end(); ++ojter) {
                    geometry_t& jgeo = *ojter;
                    double gdist, scdist;
                    vector<dpoint2> dvec;

                    geometry_distance(gdist, dvec, scdist, igeo, jgeo, false);

                    //if (dist.first > maxdist)  // !!!!!?????
                    //    maxdist = dist.first;  // !!!!!?????
                    if (scdist > maxdist)  // !!!!!?????
                        maxdist = scdist;  // !!!!!?????

                }
            }
            m(i, j) = m(j, i) = maxdist;
            //if (maxdist <= thresh) {
            //    cout << i << " ~ " << j << " ";
            //}
        }
    }
    */
}



// Returns distance between two geometies.
// We assume that geo[0] is geometry of central part; 
// igeo[1:] and jgeo[1:] and returns the one with the smallest distance:
// 'benergy' is TPS bending energy
// 'scmatching' is average shape context matching
// return value is permutation on 'igeo' used to match subparts of 'igeo' to 'jgeo'
vector<int> part_learning::geometry_distance(double& benergy, vector<dpoint2>& dvector, double& scmatching, 
    const geometry_t& igeo, const geometry_t& jgeo, bool calc_benergy)
{
    ip2_vector iv, jv;
    vector<int> perm;
    double em, scm;
    
    centers_from_geometry(iv, igeo);
    centers_from_geometry(jv, jgeo);
    
    iv.erase(iv.begin());  // assume that center (always first elt.) maps to center
    jv.erase(jv.begin());

    min_distance_matching(perm, iv, jv);

    path_map_t ipm, jpm;

    path_map_union(ipm, igeo[0], 0);
    path_map_union(jpm, jgeo[0], 0);

    for (int i = 0; i < (int)perm.size(); ++i) {
        path_map_union(ipm, igeo[i + 1], perm[i] + 1);
        path_map_union(jpm, jgeo[i + 1], i + 1);
    }

    // ipm inhibition -- to reduce the number of points 
    // while calculating bending energy and scenergy
    inhibit_path_map(ipm, max_matching_size);

    part_geometry_matching(benergy, dvector, scmatching, ipm, jpm, calc_benergy);
    scmatching /= ipm.size();

    // Expand perm to return value (prepend 0) and add 1 to rest
    perm += 1;
    perm.insert(perm.begin(), 0);
    return perm;
}

// Export "stat" to fname
// Experimental!
void part_learning::export_stat(const string& fname, int max_part_number)
{
    cout << "export_stat not implemented (experimental - debug function)." << endl;
    
    throw;

    /*
    typedef pair<double, const sym_part_t*> vector_item_t;
    typedef map<vector<int>, ip2_vector> dist_map_t;
    typedef vector<dist_map_t> dist_vector_t;

    vector<vector_item_t> parts;
    ofstream os;

    os.open(fname.c_str());

    sort_stat(parts);
    max_part_number = min<int>((int)parts.size(), max_part_number);
    
    os << '{';

    for (int i = 0; i < max_part_number; ++i) {
        stat_t::iterator iter = stat.find(*parts[i].second);
        dist_vector_t geo(iter->first.size());
        occurrences_t& occ = iter->second.occ;
        int gsize = (int)geo.size();

        for (occurrences_t::iterator oiter = occ.begin(); oiter != occ.end(); ++oiter) {
            for (int g = 0; g < gsize; ++g) {
                for (path_map_t::iterator miter = oiter->at(g).begin(); miter != oiter->at(g).end(); ++miter) {
                    geo[g][miter->first].push_back(miter->second.p);
                }
            }
        }
        if (i != 0) os << ',' << '\n';
        os << '{';
        for (int g = 0; g < gsize; ++g) {
            if (g != 0) os << ',';
            os << '{';
            for (dist_map_t::iterator miter = geo[g].begin(); miter != geo[g].end(); ++miter) {
                if (miter != geo[g].begin()) os << ',';
                os << '{';
                for (int j = 0; j < (int)miter->second.size(); ++j) {
                    if (j != 0) os << ',';
                    os << '{' << miter->second[j].x << ',' << miter->second[j].y << '}';
                }
                os << '}';
            }
            os << '}';
        }
        os << '}';
    }
    os << '}';

    os.close();
    */
}

struct avl_operator_hpoint {
    void add_to(hpoint_t& p, const hpoint_t& q) const
    { 
        p.p += q.p; 
        if (p.h.empty()) p.h = q.h; else p.h += q.h; 
    }

    hpoint_t div_int(const hpoint_t& p, int n) const
    {
        hpoint_t result;

        result.p = p.p; result.p /= n;
        result.h = p.h; result.h /= (double)n;
        return result;
    }
};

// Take average geometry over all occurrences for each part.
void part_learning::get_average_stat(map<sym_part_t, geometry_t>& stat2)
{
    stat2.clear();

    for (stat_t::iterator siter = stat.begin(); siter != stat.end(); ++siter) {
        int gsize = (int)siter->first.size();
        geometry_t& geo2 = stat2[siter->first];

        geo2.resize(gsize);
        for (int g = 0; g < gsize; ++g)
            siter->second.geol.get_mean(geo2[g], g);
    }
/*  
    stat2.clear();
    for (stat_t::iterator siter = stat.begin(); siter != stat.end(); ++siter) {
        occurrences_t& occ = siter->second.occ;
        int gsize = (int)siter->first.size();
        geometry_t& geo2 = stat2[siter->first];

        geo2.resize(gsize);
        for (int g = 0; g < gsize; ++g) {
            average_learning<path_map_t::key_type, path_map_t::mapped_type, avl_operator_hpoint> avl;

            for (occurrences_t::iterator oiter = occ.begin(); oiter != occ.end(); ++oiter) 
                avl.update(oiter->at(g));
            avl.get_average_map(geo2[g]);

            //for (occurrences_t::iterator oiter = occ.begin(); oiter != occ.end(); ++oiter) {
            //    path_map_t& part = oiter->at(g);
            //    
            //    for (path_map_t::iterator emiter = part.begin(); emiter != part.end(); ++emiter) {
            //        hpoint_t& p2 = part2[emiter->first];

            //        p2.p += emiter->second.p;
            //        if (p2.h.empty()) p2.h = emiter->second.h;
            //        else p2.h += emiter->second.h;
            //    }
            //    ++ocount;
            //}
            //for (path_map_t::iterator em2iter = part2.begin(); em2iter != part2.end(); ++em2iter) {
            //    em2iter->second.p /= ocount;
            //    em2iter->second.h /= (double)ocount;
            //}
        }
    }
*/
}

// Makes clusters ('result') according to similariry 
// 'perm' contains requested permutation for each part.
// 'stat2' contains average geometry of parts; note: subparts are permuted!
// 'max_part_number' max number of clusters
// 'space_size' max number of parts taken from 'part_learning::stat'
// 'ethresh', 'scthresh' energy and shape context thresholds for similarity
// Todo: Change algorithm: 
//   select part which is similar to most parts, ...
void part_learning::merge_stat_sc(vector<vector<similarity_item> >& result, 
    map<sym_part_t, vector<int> >& permap,
    map<sym_part_t, geometry_t>& stat2,
    int max_part_number, int sorting_type, int cluster_size, int space_size, double ethresh, double scthresh)
{
    typedef map<sym_part_t, geometry_t> stat2_t; // "average" version of stat_t
    typedef pair<double, const sym_part_t*> vector_item_t;
    typedef vector<geometry_t*> merged_item_t;
    typedef vector<merged_item_t> merge_t;
    typedef vector<similarity_item> result_item_t;
    typedef vector<result_item_t> result_t;
    typedef map<sym_part_t, vector<int> > permap_t;

    vector<vector_item_t> parts;

    sort_stat(parts, sorting_type);
    get_average_stat(stat2);

    if (space_size < (int)parts.size())
        parts.resize(space_size);

    // Merge stat, take max_part_number merged parts

    //vector<bool> used(parts.size(), false);
    vector<vector<int> > permutations(parts.size(), vector<int>());
    int unused = (int)permutations.size();

    result.clear();
    while ((int)result.size() < max_part_number && unused > 0) {
        //cout << "U="<<unused;

        // Find the first unused part and inset it to 'last' set
        for (int i = 0; i < (int)permutations.size(); ++i) {
            if (permutations[i].empty() /*!used[i]*/) {
                result.push_back(result_item_t());
                result.back().push_back(similarity_item(i, 0.0, vector<int>()));
                //used[i] = true; 
                permutations[i] = identity_permutation((int)parts[i].second->size());
                --unused; 
                break;
            }
        }

        result_item_t& last = result.back();

        // Go through all parts, add them in a greedy way to the result.back().
        for (int i = 0; i < (int)parts.size(); ++i) {
            if (last[0].index == i) // we do not compare "seed" part (1st part) to itself;
                continue;

            stat2_t::iterator iter = stat2.find(*parts[i].second);
            bool insert = true;
            double dist2first = 0.0;
            vector<int> perm;

            // check if i-th part is close to all parts in 'last'.
            for (int j = 0; j < (int)last.size() && insert; ++j) {
                stat2_t::iterator jter = stat2.find(*parts[last[j].index].second);
                double benergy, scdist;
                vector<dpoint2> dvector;
                
                if (jter->second.size() != iter->second.size()) insert = false;
                else {
                    // match iter->second to jter->second
                    vector<int> tperm = perm;

                    perm = geometry_distance(benergy, dvector, scdist, iter->second, jter->second, ethresh < 1000); 
                    insert = (tperm.empty() || perm == tperm) && (benergy <= ethresh) && (scdist <= scthresh);
                    if (j == 0) 
                        dist2first = scdist;
                }
            }
            if (insert && (permutations[i].empty() || permutations[i] == perm)) {
                last.push_back(similarity_item(i, dist2first, vector<int>()));
                if (permutations[i].empty()/*!used[i]*/) { 
                    //used[i] = true; 
                    permutations[i] = perm;
                    --unused; 
                    //permute_range<path_map_t>(iter->second.begin(), perm);
                }
            }
        }
    }
    cout << "MEND";

    // Fill 'permap'
    for (int i = 0; i < (int)permutations.size(); ++i) {
        permap.insert(permap_t::value_type(*parts[i].second, permutations[i]));
    }

    // Adjust cluster sizes
    if (cluster_size < 1) cluster_size = 1;
    for (result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
        result_item_t& cluster = *iter;

        if ((int)cluster.size() > cluster_size)
            cluster.resize(cluster_size);
    }

    //for (int i = 0; i < (int)result.size(); ++i) {
    //    for (int j = 0; j < result[i].size(); ++j) {
    //        cout << result[i][j].index << ' ';
    //    }
    //    cout << endl;
    //}
}

/*class clique_data {
public:
    set<int> parts; // Members of the "clique" (including part0)
    int part0;      // Initial part (name of the "clique")
    cv::Mat mean;
    cv::Mat eigenvectors;

    clique_data& operator=(const clique_data& cd) 
    { 
        parts = cd.parts; part0 = cd.part0; mean = cd.mean; eigenvectors = cd.eigenvectors;
        return *this;
    }

    clique_data(const clique_data& cd) : parts(cd.parts), part0(cd.part0), 
        mean(cd.mean), eigenvectors(cd.eigenvectors) { }
    clique_data() : parts(), part0(0), mean(), eigenvectors() { cout << 'C'; }
    ~clique_data() { }
};*/

// Makes clusters ('result') according to pca similarity 
// Note that, due to "compatibility" with merge_stat_sc, some parameters are just
// filled with default values
// 'perm' filled with identities
// 'stat3' contains pca varaibility of parts
// 'space_size' max number of parts taken from 'part_learning::stat'
// 'norm' measure for similarity in terms of original - back-projection norm
void part_learning::merge_stat_pca(vector<vector<similarity_item> >& result, 
    map<sym_part_t, vector<int> >& permap,
    map<sym_part_t, pca_data>& stat3,
    int space_size, int sorting_type, double norm)
{
    typedef map<sym_part_t, pca_data> stat3_t; 
    typedef pair<double, const sym_part_t*> vector_item_t;
    typedef vector<geometry_t*> merged_item_t;
    typedef vector<merged_item_t> merge_t;
    typedef vector<similarity_item> result_item_t;
    typedef vector<result_item_t> result_t;
    typedef map<sym_part_t, vector<int> > permap_t;
    typedef vector<int> perm_t;

    vector<vector_item_t> parts;

    sort_stat(parts, 0);
    get_pca_stat(stat3);

    if (space_size < (int)parts.size())
        parts.resize(space_size);
    
    vector<vector<dpoint2> > geovec(parts.size(), vector<dpoint2>());
    vector<clique_data> cliques(parts.size());

    // 0. Fill 'cliques'
    for (int i = 0; i < (int)parts.size(); ++i) {
        const sym_part_t* sp = parts[i].second;

        cliques[i].parts.insert(i);
        cliques[i].part0 = i;
        cliques[i].mean = stat3[*sp].mean;
        cliques[i].eigenvectors = stat3[*sp].eigenvectors;
    }    

    // 1. Get part geometry (from mean) for each part
    for (vector<clique_data>::const_iterator cliter = cliques.begin(); cliter != cliques.end(); ++cliter) 
        geovec[cliter->part0] = partition(cliter->mean);

    // 2. Calculate norm of part geometry to its "back"-projection from clique space -- for each part and each clique
    //    matrix(p, c)
    matrix<double> cliquedist((int)parts.size(), (int)cliques.size());
    matrix<vector<int> > cliqueperm((int)parts.size(), (int)cliques.size());

    for (int pi = 0; pi < (int)parts.size(); ++pi) {
        vector<dpoint2> gv = geovec[pi];
        cout << '.';

        for (int ci = 0; ci < (int)cliques.size(); ++ci) {
            vector<dpoint2> v = partition(cliques[ci].mean);
            vector<dpoint2> gvresized = get_resized_vector(gv, (int)v.size());

            translate_and_scale(gvresized);

            vector<int> perm = point_matching(gvresized, v);

            //write_matching(string("c:\\work\\tmp\\match") + pi + string("-") + ci + string(".m"), gvresized, v);
            
            permute(gvresized, perm);
            
            cv::Mat data = flatten(gvresized);
            cv::Mat coeffs;
            cv::Mat bproj;
            cv::Mat result;
            
            gemm(data - cliques[ci].mean, cliques[ci].eigenvectors, 1, cv::Mat(), 0, coeffs, cv::GEMM_2_T);
            gemm(coeffs, cliques[ci].eigenvectors, 1, cliques[ci].mean, 1, result, 0);
            cliquedist(pi, ci) = cv::norm(result, data, cv::NORM_L2);
            cliqueperm(pi, ci) = perm;
        }
    }
    cout << endl;

    // Remove redundant cliques based on difference on back-projection
    vector<int> toremove;

    while (true) {

        // 3. Clique D is absorbed in clique C if D can be expressed with basis of C with small error
        vector<vector<int> > nabsorbed(cliques.size(), vector<int>());   // nabsorbed[i] counts cliques absorbed in clique i

        for (int ci = 0; ci < (int)cliques.size(); ++ci) {
            if (cliques[ci].parts.empty())
                continue;

            for (int di = 0; di < (int)cliques.size(); ++di) {
                const set<int>& D = cliques[di].parts;
                double maxdiff = 0.0;

                if (D.empty())
                    continue;

                for (set<int>::const_iterator iter = D.begin(); iter != D.end(); ++iter) {
                    int e = *iter;
                    double diff = cliquedist(e, ci);

                    if (diff > maxdiff) maxdiff = diff;
                }
                //cout << di << " is absorbed in " << ci << " with " << maxdiff << endl;
                if (ci != di && maxdiff <= norm) 
                    nabsorbed[ci].push_back(di);
                
            }
        }

        // 4. Find clique which absorbes most cliques
        vector<vector<int> >::iterator miter = max_element(nabsorbed.begin(), nabsorbed.end(), size_less<int>);

        // 4a. If there is no clique which absorbes some other clique, quit
        if (miter == nabsorbed.end() || miter->size() == 0)
            break;

        // 5. Keep such clique and remove absorbed cliques
        int keepi = (int)(miter - nabsorbed.begin());
        vector<int>& torm = *miter;

        for (vector<int>::iterator riter = torm.begin(); riter != torm.end(); ++riter) {
            cliques[keepi].parts.insert(cliques[*riter].parts.begin(), cliques[*riter].parts.end());
            toremove.push_back(*riter);
            cliques[*riter].parts.clear();
        }
    }
    sort(toremove.begin(), toremove.end(), greater<int>());
    for (vector<int>::iterator riter = toremove.begin(); riter != toremove.end(); ++riter) 
        cliques.erase(cliques.begin() + *riter);
    
    // Fill result
    result.clear();
    result.resize(cliques.size());

    for (int i = 0; i < (int)cliques.size(); ++i) {
        clique_data& clique = cliques[i];

        result[i].push_back(similarity_item(clique.part0, 0.0, vector_range(0, (int)geovec[clique.part0].size() - 1)));
        for (set<int>::iterator citer = clique.parts.begin(); citer != clique.parts.end(); ++citer) {
            if (*citer != clique.part0) 
                result[i].push_back(similarity_item(*citer, cliquedist(*citer, clique.part0), cliqueperm(*citer, clique.part0)));
        }
    }

    // Fill 'permap' (trivial)
    for (int i = 0; i < (int)parts.size(); ++i) {
        const sym_part_t* sp = parts[i].second;

        permap.insert(permap_t::value_type(*sp, vector_range(0, sp->size() - 1)));
    }

    // Adjust cluster sizes
    // ?
}

void part_learning::merge_stat_sc(vector<vector<similarity_item> >& result, 
    int max_part_number, int cluster_size, int space_size, int sorting_type, double ethresh, double scthresh)
{
    map<sym_part_t, geometry_t> stat2;
    map<sym_part_t, vector<int> > perm;
    
    merge_stat_sc(result, perm, stat2, max_part_number, cluster_size, space_size, sorting_type, ethresh, scthresh);
}

void part_learning::get_pca_stat(map<sym_part_t, pca_data>& stat3)
{
    typedef map<sym_part_t, pca_data> stat3_t;

    stat3.clear();
    for (stat_t::iterator siter = stat.begin(); siter != stat.end(); ++siter) {
        int gsize = (int)siter->first.size();
        pca_data& geo3 = stat3[siter->first];

        geo3 = siter->second.pcal.get_pca_data();
    }

    //for (stat_t::iterator siter = stat.begin(); siter != stat.end(); ++siter) {
    //    occurrences_t& occ = siter->second.occ;
    //    int gsize = (int)siter->first.size();
    //    pca_data& geo3 = stat3[siter->first];
    //    pca_learning pcal;

    //    for (occurrences_t::iterator oiter = occ.begin(); oiter != occ.end(); ++oiter) {
    //        vector<ipoint2> pts;

    //        for (int g = 0; g < gsize; ++g) {
    //            vector<ipoint2> gpts = get_path_map_points(oiter->at(g));

    //            pts.insert(pts.end(), gpts.begin(), gpts.end());
    //        }
    //        pcal.update(pts);
    //    }
    //    geo3 = pcal.get_pca_data();
    //}
}

//void set_part_geometry(node* p, const vector<int>& perm, const geometry_t& geo)
//{
//    int name = atom("lyrSrc");
//
//    foreach_naighbor(
//}

/*
// Merges stat (using merge_stat) and adds info about merging to 'library'.
void part_learning::merge_library(part_lib* library, int layer, int max_part_number, int cluster_size, 
        int space_size, double ethresh, double scthresh)
{
    typedef map<sym_part_t, geometry_t> stat2_t;
    typedef pair<double, const sym_part_t*> vector_item_t;
    typedef pair<int, double> merge_pair_t;
    typedef vector<merge_pair_t> merge_item_t;

    if (layer < 0 || layer > library->max_layer_index())
        throw new_libhop_exception("part_learning::merge_library: invalid layer.");

    stat2_t stat2;
    vector<vector_item_t> parts;
    vector<merge_item_t> merge;

    sort_stat(parts);
    merge_stat(merge, stat2, max_part_number, cluster_size, space_size, ethresh, scthresh);
    
    // Do "merging"
    // - add geometry data
    // - add lyrSimilar and lyrSimRoot edges
    //      lyrSimilar, edge data is similarity weight, type: edge_data_t<double>
    //      lyrSimBack, no edge data
    vector<node*>& libparts = library->parts[layer];

    if (parts.size() < libparts.size())
        throw new_libhop_exception("part_learning::merge_library: library layer and stat size mismatch.");

    for (int m = 0; m < (int)merge.size(); ++m) {
        merge_item_t& cluster = merge[m];


        for (int i = 0; i < (int)cluster.size(); ++i) {

        }
    }

}
*/

void part_learning::save_stat()
{
    for (auto siter = stat.begin(); siter != stat.end(); ++siter) {
        auto sypart = siter->first;
        double count = siter->second.count;

        for (auto spiter = sypart.begin(); spiter != sypart.end(); ++spiter) {
            if (spiter != sypart.begin()) cout << ',';
            cout << spiter->x << ',' << spiter->y;
        }
        cout << ',' << count << endl;
    }
}

// obj_learning
///////////////////////////////////////////////////////////////////////////////


obj_learning::obj_learning(const config_dictionary& cfg) :
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

void obj_learning::cfg_init(const config_dictionary& cfg)
{
    string libname;
    
    //cfg.get_value(libname, "library", true);
    cfg.get_value(layer, "layer", true);

    //read_library(libname, library);
    //if (library) {
    //    libraryD = (part_lib*)library->get_copy_s();
    //    libraryD->delete_parts_geq(layer + 1);
    //}

    cfg.get_value(r_response_threshold, "r_response_threshold", 0.0);
    cfg.get_value(g_response_threshold, "g_response_threshold", 0.0);
    cfg.get_value(rr_response_threshold, "rr_response_threshold", 0.0);
    contraction = cfg.get_value_double("contraction", 1.0);
    cfg.get_value(max_objects, "max_objects", 4);
    cfg.get_value(max_add, "max_add", 0);
    if (max_add[0] <= 0) max_objects.set_val(max_add[0]);
    cfg.get_value(max_cluster_n, "max_cluster_n", 6);
    cfg.get_value(min_cluster_n, "min_cluster_n", 4);
    cfg.get_value(cluster_size, "cluster_size", 10);
    gaussian_dim = cfg.get_value_int("gaussian_dim", 5);
    gaussian_sigma = cfg.get_value_double("gaussian_sigma", 2.0);
    cfg.get_value(cluster_member_threshold, "cluster_member_threshold", 0.5);
    cfg.get_value(cover_threshold, "cover_threshold", 0.5);
    cover_threshold0 = cfg.get_value_double("layer0_cover_ratio_threshold", 0.0);
    validation_threshold = cfg.get_value_double("validation_threshold", 1.5);
    cfg.get_value(intersection_threshold, "intersection_threshold", 0.2);
    max_depth = cfg.get_value_int("max_depth", 0);
    reduce_radius = cfg.get_value_int("reduce_radius", 3);
}

/*void obj_learning::object_from_result(layer1_result* res)
{
    typedef pair<int, node*> queue_t;

    if (layer < 0 || layer >= (int)res->shape_nodes.size()) return;

    if (!res->grid(layer)) res->init_grid(layer);

    vector<node*>& s_nodes = res->shape_nodes[layer];
    vector<vector<ipoint2> > clusters;
    int to_prev_layer = atom("toPrevLayer").get_index();
    obj_data_t object;

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;

        clusters.push_back(vector<ipoint2>());
        clusters.back().push_back(ipoint2(nd->x, nd->y));
    }
    while ((int)clusters.size() > cluster_n) {
        hierarchical_clustering_iter(clusters, ipoint2_set_distance);
    }
    for (vector<vector<ipoint2> >::iterator citer = clusters.begin(); citer != clusters.end(); ++citer) {
        vector<ipoint2>& cluster = *citer;
        //node* bestn = nullptr;
        //int bests = 0;
        priority_queue<queue_t> queue;

        for (vector<ipoint2>::iterator piter = cluster.begin(); piter != cluster.end(); ++piter) {
            ipoint2& p = *piter;
            irectangle2 box;
            set<node*> proj;
            node* n = res->node_at(p.x, p.y, layer);

            res->recurse_from_node(n, to_prev_layer, proj);
            box = node_set_bounding_rectangle(proj.begin(), proj.end());
            if (!box.invalid())
                queue.push(queue_t(box.area(), n));
            //if (bestn == nullptr || (int)proj.size() > bests) {
            //    bests = (int)proj.size();
            //    bestn = n;
            //}
        }
        //
        int c = 0;
        cluster_data_t cd;
        int top_val = -1;

        while ((int)cd.first.size() < cluster_size && !queue.empty()) {
            if (top_val < 0) top_val = queue.top().first;
            else if (queue.top().first < cluster_member_threshold * (double)top_val) 
                break;

            layer1_data* nd = (layer1_data*)queue.top().second->data;

            if (c++ == 0) { 
                cd.second.x = nd->x; 
                cd.second.y = nd->y; 
            }
            if (find(cd.first.begin(), cd.first.end(), nd->m) == cd.first.end())
                cd.first.push_back(nd->m);
            queue.pop();
        }
        object.push_back(cd);
    }
    objects.push_back(object);
}*/


// Choses a minimal (approximately) subset of the set 's'
// which covers the set 's' according to the neighborhood
// given by the edges with name 'ename'.
void reduce_node_set(set<node*>& result, const set<node*>& s, int ename)
{
    set<node*> uncovered(s);

    result.clear();
    while (!uncovered.empty()) {
        node* bestn = nullptr;
        int bestnbsize = -1;

        for (set<node*>::iterator iter = uncovered.begin(); iter != uncovered.end(); ++iter) {
            node* n = *iter;
            int nbsize = 0;

            foreach_neighbor(n, ename, niter) {
                if (uncovered.find(neighbor_node(niter)) != uncovered.end()) ++nbsize;
            }
            if (nbsize > bestnbsize) { bestnbsize = nbsize; bestn = n; }
        }   

        result.insert(bestn);
        foreach_neighbor(bestn, ename, niter) {
            uncovered.erase(neighbor_node(niter));
        }
        uncovered.erase(bestn);
    }
}

irectangle2 node_set_extent(const set<node*>& nset)
{
    irectangle2 result;
    
    for (set<node*>::const_iterator iter = nset.begin(); iter != nset.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;

        result.eat(nd->x, nd->y);
    }
    return result;
}

struct recurse0_check {
    bool operator()(node* n) { return n->is_attr_set(IMG_NODE_ATTR); }
};

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

//int debugNO = 1;

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

    // D E B U G 
    //string dfn = string("c:\\temp\\obj_learning_") + debugNO + string(".m");
    //debugNO++;
    //save_node_set_mma(dfn, covered);
    // D E B U G

    if (!res->grid(layer)) res->init_grid(layer);

    vector<node*>& s_nodes = res->shape_nodes[layer];
    int ename = atom("toPrevLayer");
    //int nbname = atom("toNeighbor");
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
        //set<node*> redset;
        //reduce_node_set(redset, miter->second, nbname);
        //irectangle2 ext = node_set_extent(miter->second);

        //cout << ' ' << (int)sqrt((double)ext.size2());
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
                /* liter->first >= cover_threshold * topsize*/ 
                /*(int)nset.size() < min_cluster_n !!! */
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
            //library->get_path_map(pml, library->parts[nd->z][nd->m]);
            get_path_map(pm, res, scmap, n, true);
            //augment_path_map(pm, pml, true);
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
                            //library->get_path_map(pml, library->parts[nd2->z][nd2->m]);
                            get_path_map(pm, res, scmap, n2, true);
                            //augment_path_map(pm, pml, true);
                            add_to_path_map(cd.geo[nd2->m], pm);

                            nbm.insert(nd2->m);
                        }
                    }
                }
            }
            cd.pos = ipoint2(nd->x, nd->y);
            
            object.push_back(cd);
            tabu.push_back(tabu_item_t(ipoint2(nd->x, nd->y), nbm));

            // D E B U G  -  D I S P L A Y
            //cout << '{';
            //for (vector<ipoint2>::iterator diter = cd.first.begin(); diter != cd.first.end(); ++diter) {
            //    if (diter != cd.first.begin()) cout << ',';
            //    cout << '(' << diter->x << ',' << diter->y << ')';
            //}
            //cout << "},(" << cd.second.x << ',' << cd.second.y << ')' << endl;
            // E N D   D E B U G  -  D I S P L A Y
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
                
                // D E B U G 
                //string dfn = string("c:\\temp\\obj_learning_") + debugNO + string("_") + depth + string(".m");
                //debugNO++;
                //save_node_set_mma(dfn, rset);
                // D E B U G

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
    res->delete_edges(atom("toLayer0"));

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


bool obj_learning::validate_object(layer1_result* res, const irectangle2& gtr, 
    const vector<node*>& cluster_nodes, const obj_data_t& object)
{
    if (!vset.empty() && creator != nullptr) 
        return validate_object(res, gtr, object);
    else {
        set<node*> reconstruction;
        int layer0_size = (int)res->shape_nodes[0].size();

        res->recurse(cluster_nodes, atom("toPrevLayer").get_index(), reconstruction);
        cout << "cover ratio: " << (double)reconstruction.size()/layer0_size << endl;
        if ((double)reconstruction.size()/layer0_size >= cover_threshold) 
            return true;
        else 
            return false;
    } 
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
        // D E B U G 
        //rnew->save(string("c:\\temp\\val") + (counter++) + string("_0.ly6"));
        counter++;
        gtr.push_back(iter->second);
        rnew->check_with_groundtruth(v, gtr, layer + 1, set<int>(), 0.0);
        // D E B U G
        //cout << "(" << v << ")";
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

/*void obj_learning::add_validation_data(layer1_result* res, const irectangle2& gtruth)
{
    vset.push_back(validation_data_t(streamed_pointer(res), gtruth));
}*/

void obj_learning::set_creator(layern_creator* newc)
{
    if (creator != nullptr) delete creator;
    creator = newc;
}

bool obj_learning::validation_function_1(int tr, int fa)
{
    return (fa == 0) ? tr > 0 : (double)tr/fa > validation_threshold;
}

// Normalizes geometry data of part_lib::cluster_data_t, should only be used before adding to library
// geometry data of subparts are sums of equivalent parts, therefore they need to be normalized (average)
// before adding to library!
/*void obj_learning::normalize_object_data(obj_data_t& od)
{
    for (obj_data_t::iterator iter = od.begin(); iter != od.end(); ++iter) {
        for (map<int, path_map_t>::iterator pmiter = iter->geo.begin(); pmiter != iter->geo.end(); ++pmiter) 
            divide_path_map(pmiter->second, (int)iter->types.size());
    }
}*/

// Reduces geometry (path maps) with class parameter 'reduce_radius'
void obj_learning::reduce_object_data(obj_data_t& od)
{
    for (obj_data_t::iterator iter = od.begin(); iter != od.end(); ++iter) {
        for (map<int, path_map_t>::iterator pmiter = iter->geo.begin(); pmiter != iter->geo.end(); ++pmiter) {
            pmiter->second = reduce_path_map(pmiter->second, reduce_radius);
        }
    }
}

void obj_learning::add_additional_lib_nodes(part_lib* library, obj_data_t& obj, int m)
{
    vector<node*>& parts1 = library->parts[layer];
    node* cn = library->parts[layer + 1][m];
    int name = atom("lyrSrcM").get_index();
    int to_part = atom("lyrSrc").get_index();
    int to_center_back = atom("lyrCenterBack").get_index();
    int to_center = atom("lyrCenter").get_index();
    vector<ipoint2>& centers = obj.front().types;

    // to center
    for (vector<ipoint2>::iterator i = centers.begin(); i != centers.end(); ++i) {
        library->add_edge_unique(cn, parts1[i->y], to_center, to_center_back);
    }

    // from center to "other" parts 
    for (obj_data_t::iterator oiter = ++obj.begin(); oiter != obj.end(); ++oiter) {
        vector<ipoint2>& v = oiter->types;
        part_data_2* ed = (part_data_2*)library->get_edge_data(cn, 
            library->parts[layer + 1 + v.front().x][v.front().y], to_part);

        for (vector<ipoint2>::iterator i = v.begin(); i != v.end(); ++i) {
            library->add_edge_2(cn, 
                library->parts[layer + 1 + i->x][i->y], new part_data_2c(ed->x, ed->y), name);
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

        //for (int i = 0; i < (int)origobject.size(); ++i) {
        //    obj_data_t object(origobject);

        //    object.erase(object.begin() + i);
        //    object.insert(object.begin(), origobject[i]);

        //    vector<iipair> str;
        //    vector<iipair> coo; 
        //    vector<matrix<double>*> distr;
        //    ipoint2 coo0; 
        //    matrix<double> mask;
        //    int newm;

        //    gaussian_mask(gaussian_dim, gaussian_dim, gaussian_sigma, mask);
        //    coo0 = object.front().second;
        //    coo.push_back(iipair(0, 0));
        //    str.push_back(object.front().first.front().to_pair());
        //    distr.push_back(&mask);
        //    for (obj_data_t::iterator oiter = ++object.begin(); oiter != object.end(); ++oiter) {
        //        coo.push_back((oiter->second - coo0).to_pair());
        //        str.push_back(oiter->first.front().to_pair());
        //        distr.push_back(&mask);
        //    }
        //    newm = lib->add_part(layer + 2, new part_data(), str, coo, distr, -1, contraction, 0, 0);
        //    if (newm >= 0) {
        //        newparts.push_back(newm);
        //        add_additional_lib_nodes(lib, object, newm);
        //    }
        //}
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

o_learning::o_learning(const config_dictionary& cfg) 
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

void o_learning::init_cfg(const config_dictionary& cfg)
{
    string libname;
    
    cfg.get_value(libname, "library", true);
    read_library(libname, library);
    srclibrary = (part_lib*)library->get_copy_s();

    cfg.get_value(srclayer, "src_layer", true);
    cfg.get_value(contraction, "contraction", true);    
    catname = cfg.get_value_string("category_name", "");
    max_cluster_n = cfg.get_value_int("max_cluster_n", 6);
    min_cluster_n = cfg.get_value_int("min_cluster_n", 4);
    cluster_size = cfg.get_value_int("cluster_size", 10);
    cluster_member_threshold = cfg.get_value_double("cluster_member_threshold", 0.7);
    intersection_threshold = cfg.get_value_double("intersection_threshold", 0.7);
    hit_ratio_threshold = cfg.get_value_double("hit_ratio_threshold", 0.8);
    hit_threshold = cfg.get_value_double("hit_threshold", 0.1);
    redundancy_threshold = cfg.get_value_int("redundancy_threshold", 10);
    type_bite_threshold = cfg.get_value_int("type_bite_threshold", 5);
    max_models = cfg.get_value_int("max_models", 500);
    gaussian_dim = cfg.get_value_int("gaussian_dim", 5);
    gaussian_sigma = cfg.get_value_double("gaussian_sigma", 2.0);

    double gdim = cfg.get_value_int("gaussian_dim", 7);
    double gsigma = cfg.get_value_double("gaussian_sigma", 2.0);
    gaussian_mask(gdim, gdim, gsigma, dist);

	inference_cfg.from_namespace_priority(cfg, 1, "validation");
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
        //vector<ip2_vector> clusters;
        ip2_vector v = iter->second;
        ip2_vector mv;

        //cout << '.';
        //sort(v.begin(), v.end());
        //unique(v.begin(), v.end());

        //for (ip2_vector::iterator viter = v.begin(); viter != v.end(); ++viter) {
        //    clusters.push_back(ip2_vector());
        //    clusters.back().push_back(*viter);
        //}
        //diameter_clustering(clusters, maxdiam2);

        //v.clear();
        //for (vector<ip2_vector>::iterator citer = clusters.begin(); citer != clusters.end(); ++citer) {
        //    ipoint2 avg = average<ipoint2>(citer->begin(), citer->end());

        //    v.push_back(avg);
        //}
        find_maxima(mv, v);
        if (!mv.empty()) 
            dupmap.insert(dupmap_t::value_type(iter->first, mv));
    }
}

void o_learning::print_models(const list<olv_object_part>& models)
{
    for (list<olv_object_part>::const_iterator i = models.begin(); i != models.end(); ++i) {
        i->print(cout);
        cout << endl;
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

        n->get_neighbor_set(atom("toLayer0"), nset);
        
        irectangle2 box = node_set_bounding_rectangle(nset.begin(), nset.end());

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
    //res->add_reconstruction_edges_leq_fwd(srclayer);
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

//double o_learning::check_duplet(olv_object_part& model, layer1_result* res, int m, const ipoint2& p)
//{
////dnpair layer1_result::schur_product_max_all(int xo, int yo, int dx, int dy, matrix<double>& m, int type, int z, double c,
////        schur_product_function* f, double thresh, sp_result_t& res)
//    if (!res->grid(srclayer)) 
//        res->init_grid(srclayer);
//    for (olv_object_part::structure_t::iterator siter = model.str.begin(); siter != model.str.end(); ++siter) {
//        vector<node*> r;
//
//        res->schur_product_max(
//    }
//
//}

int maxvsize = 0;

/* Throw out */

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

                (*riter)->get_neighbor_set(atom("toLayer0"), nset);
                node_set_to_point_set(pset, nset.begin(), nset.end());
                if (intersection_size(model.rec, pset)/(double)pset.size() < 0.4 /* £££ */) {
                    v.push_back(result_item_t(rn, pset)); 
                }
            }
        }
    }

    compress_node_list(v, 0.6 /* £££££ */);

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

    // cout << "maxmodels " << maxmodels << endl;
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
        if (usedm > 3 /* £££££ */)  
            break;
    }
}

void debug_check1(vector<node*>& lnodes, vector<pair<double, int> >& statmap)
{
    for (vector<node*>::iterator i = lnodes.begin(); i != lnodes.end(); ++i) {
        if (node_type(*i) >= statmap.size()) {
            cout << "statmap size/node type mismatch: ";
            cout << "node type: " << node_type(*i) << "; statmap.size(): " << statmap.size() << endl;
            throw exception();
        }
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

        n->get_neighbor_set(atom("toLayer0"), nset);
        node_set_to_point_set(pset, nset.begin(), nset.end());

        if (max_intersection_percent(recmem, pset) < 0.7 /* £££££ */) {
            olv_object_part model;

            model.augment(n, pset);
            augment_model(models, modelcount, model, mstat, indexmap, res, max_models/4 /* ££££££ see 3 below */);
            recmem.push_back(pset);
            ++used;
        }

        // We consider only a few -- best -- types for initial model part
        if (used > 3 /* ££££££ */)
            break;
    }

    for (list<olv_object_part>::iterator iter = models.begin(); iter != models.end() ; ) {
        if (iter->size() < min_cluster_n || iter->size() > max_cluster_n) 
            iter = models.erase(iter);
        else
            ++iter;
    }
    //cout << "  used: " << used << endl;
    //cout << "  maxv: " << maxvsize << endl;
    //print_models(models); 
}

// result: vector items are pairs ((-1, type), path_map_t)
// res: layer1_result with toLayer0 edges added
// n: node
// size: check nodes in circle with diameter size/2 
void get_overlapping_parts(vector<pair<ipoint2, path_map_t> >& result, layer1_result* res, const scmap_t& scmap,
    node* n, int size, int thresh)
{
    typedef pair<ipoint2, path_map_t> result_item_t;

    int name = atom("toLayer0");
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
    // vector<ipoint2> is a set of equivalent subparts and
    // ipoint2 is a coordinate of the subpart *relative to the center*; 
    // see part_lib::add_object_part
    //typedef pair<vector<ipoint2>, ipoint2> obj_item_t;
    //    struct cluster_data_t { 
    //    vector<ipoint2> types;      // {(lyr. diff., type), ...}
    //    ipoint2 pos;                // position
    //    map<int, path_map_t> geo;   // geo 

    //    bool operator<(const cluster_data_t& cd) const { return types < cd.types; }
    //};


    matrix<double> mask;
    list<olv_object_part>::iterator miter = models.begin();

    gaussian_mask(gaussian_dim, gaussian_dim, gaussian_sigma, mask);
    while (miter != models.end()) {
        const vector<node*>& str = miter->str;
        vector<part_lib::cluster_data_t> objdata;

        //cout << '.' << str.size();
        
        for (vector<node*>::const_iterator siter = str.begin(); siter != str.end(); ++siter) {
            node* n = *siter;
            layer1_data* nd = (layer1_data*)n->data;
            vector<pair<ipoint2, path_map_t> > simnodes;

            //cout << '[';
            get_overlapping_parts(simnodes, res, scmap, n, 3, cluster_member_threshold);  // 3 !!!!!!
            //cout << ']';
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
            //rhits[biter->m
            //if (rhits.size() != reclist.size()) 
            //    rhits.resize(reclist.size(), 0);
            //++rhits[hit];
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
    //inference_cfg.set_value("realization_ratio_threshold", "0.6");
    //inference_cfg.set_value("g_response_threshold", "0.001");
    //inference_cfg.set_value("r_response_threshold", "0.001");
    //inference_cfg.set_value("candidate_r_threshold", "0.001");
    //inference_cfg.set_value("candidate_g_threshold", "0.001");
    //inference_cfg.set_value("convolution_threshold", "0.001");

    // some more changes?
    cout << "realization_ratio_threshold: " << inference_cfg.get_value_double("realization_ratio_threshold", 0.0) << endl;

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

//void o_learning::filter_equivalent_models(list<olv_object_part>& models)
//{
//    list<
//}

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
    res->delete_edges(atom("toLayer0")); // link_path <-> add_reconstruction_edges_fwd clash
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


// ell_learning
///////////////////////////////////////////////////////////////////////////////

ell_learning::ell_learning(const config_dictionary &cfg) :
    library(nullptr)
{
    cfg_init(cfg);
}

ell_learning::~ell_learning()
{
    if (library != nullptr) delete library;
}

void ell_learning::cfg_init(const config_dictionary &cfg)
{
    string libname;

    cfg.get_value(category_layer, "category_layer", true);
    //cfg.get_value(category_name, "category_name", true);
    cover_threshold = cfg.get_value_double("cover_threshold", 0.7);
    cfg.get_value(libname, "library", true);
	gtextension = cfg.get_value_string("groundtruth_extension", ".groundtruth");

    read_library(libname, library);
    init_catmap();
}

void ell_learning::init_catmap()
{
    catmap.clear();
    for (vector<node*>::iterator iter = library->parts[category_layer].begin(); 
		    iter != library->parts[category_layer].end(); ++iter) {
		node* n = *iter;
		cpart_data* nd = (cpart_data*)n->data;
		
        catmap.insert(cat_map_t::value_type(nd->name, nd->type));
	}
}


//void ell_learning::learn(string dir, const list<string>& files, const string& catname)    
void ell_learning::update(layer1_result* res, const list<pair<irectangle2, int> >& gtr)
{
    typedef map<int, set<node*> > psmap_t;
    typedef list<pair<irectangle2, int> > gtr_t;

    if (res == nullptr && res->max_layer_index() < category_layer) return;

    int to_prev = atom("toPrevLayer");
    int to_0 = atom("toLayer0");
    vector<node*>& c_nodes = res->shape_nodes[category_layer];
    //irectangle2 objrect;

    cout << "<R";
    res->add_reconstruction_edges_fwd_link(category_layer - 1);
    //res->add_reconstruction_edges_leq_fwd(category_layer - 2);
    cout << ">";
    for (vector<node*>::iterator citer = c_nodes.begin(); citer != c_nodes.end(); ++citer) {
        node* cn = *citer;
        layer1_data* cnd = (layer1_data*)cn->data;

        //if (cnd->m != catindex)
        //    continue;

        foreach_neighbor(cn, to_prev, oiter) {
            // Object
            node* on = neighbor_node(oiter);
            layer1_data* ond = (layer1_data*)on->data;
            set<node*> objnodes;
            ellipse ell;
            psmap_t psmap;
            
            // Update statistics for subparts
            foreach_neighbor(on, to_prev, iter) {
                node* sn = neighbor_node(iter);
                layer1_data* snd = (layer1_data*)sn->data;
                edge_data_name* nned = (edge_data_name*)neighbor_edge_data(iter);
                set<node*> nset;
                psmap_t::iterator miter;
                
                if (nned == nullptr) {
                    cout << "No edge name, exiting.";
                    throw exception();
                }
                sn->get_neighbor_set(to_0, nset);
                objnodes.insert(nset.begin(), nset.end());
                miter = psmap.find(nned->index);
                if (miter == psmap.end()) psmap.insert(psmap_t::value_type(nned->index, nset));
                else miter->second.insert(nset.begin(), nset.end());
            }

            bool isfalse = true;
            ell_learning_stat_item stat;
            irectangle2 recrect = bounding_rectangle_of_nodes(objnodes.begin(), objnodes.end());
            
            for (gtr_t::const_iterator gtiter = gtr.begin(); gtiter != gtr.end(); ++gtiter) {
                irectangle2 intsec = gtiter->first.intersection(recrect);
                double factor = (double)(intsec.area())/(double)gtiter->first.area();

                if (factor >= cover_threshold && gtiter->second == cnd->m) {
                    isfalse = false;
                    break;
                }
            }

            int index = get_new_sample_index(!isfalse, ond->m);

            fit_ellipse_to_nodes(ell, objnodes.begin(), objnodes.end());
            add_to_map(!isfalse, ond->m, index, ipoint2(INT_MIN, INT_MIN), ell);

            for (psmap_t::iterator miter = psmap.begin(); miter != psmap.end(); ++miter) {
                fit_ellipse_to_nodes(ell, miter->second.begin(), miter->second.end());
                add_to_map(!isfalse, ond->m, index, miter->first, ell);
            }
        }
    }

}

void ell_learning::update(layer1_result* res, const list<pair<irectangle2, string> >& gtr)
{
    typedef list<pair<irectangle2, string> > gtr_t;

    list<pair<irectangle2, int> > gtr2;

    for (gtr_t::const_iterator gtiter = gtr.begin(); gtiter != gtr.end(); ++gtiter) {
        cat_map_t::iterator iter = catmap.find(gtiter->second);

        if (iter != catmap.end()) gtr2.push_back(pair<irectangle2, int>(gtiter->first, iter->second));
        else cout << "Category '" << gtiter->second << "' not found in library; skipping update." << endl;
    }
    update(res, gtr2);
}

void ell_learning::display_stat()
{

}

void ell_learning::learn(list_t& features, int& feature_count, int type, const map_t& stat)
{
    int to_src = atom("lyrSrc");
    map_t::const_iterator iter;
    ell_learning_stat_item empty;

    features.clear();
    node* p = library->parts[category_layer - 1][type];
    part_data* pd = (part_data*)p->data;
    int index = 0;
    feature_count = 4; //3*(p->count_neighbors(to_src) + 1);
    
    while ((iter = stat.find(map_t::key_type(type, index, ipoint2(INT_MIN, INT_MIN)))) != stat.end()) {
        features.push_back(real_vector_t());
        features.back().reserve(feature_count + 1);

        iter->second.push_back_to(features.back());
        foreach_neighbor(p, to_src, iter) {
            part_data* spd = (part_data*)neighbor_node_data(iter);
            part_data_2* sped = (part_data_2*)neighbor_edge_data(iter);

            map_t::const_iterator fiter = stat.find(map_t::key_type(type, index, ipoint2(sped->x, sped->y)));
            if (fiter == stat.end()) empty.push_back_to(features.back());
            else fiter->second.push_back_to(features.back());
        }
        ++index;
    }
}

void ell_learning::learn(list_t& positive, list_t& negative, int type)
{
    int feature_count;
    list_t::iterator liter;

    learn(positive, feature_count, type, posstat);
    learn(negative, feature_count, type, negstat);

    if (positive.empty() && negative.empty())
        return;
}

//void ell_learning::learn(const string& fname, int type)
//{
//    list_t pos, neg;
//    svm2 asvm;
//
//    learn(pos, neg, type);
//    if (pos.empty() || neg.empty()) {
//        //if (pos.empty()) cout << "Type " << type << " has no positive examples." << endl;
//        //if (neg.empty()) cout << "Type " << type << " has no negative examples." << endl;
//        return;
//    }
//    
//    cout << "Type " << type << " saved." << endl;
//}

#include <ml.h>

void ell_learning::learn_all(string out_file)
{
    cout << "ell_learning not supported" << endl;
    throw exception();

    vector<node*>& parts = library->parts[category_layer - 1];

    //end_dir(out_dir);

	//CvFileStorage* xml_out_storage = cvOpenFileStorage(out_file.c_str(), nullptr, CV_STORAGE_WRITE);

	//// copy existing svm models to new storage
	//if (library->svm_models_storage != nullptr) {
	//	cv::FileNode root_node(library->svm_models_storage, cvGetRootFileNode(library->svm_models_storage, 0));
	//	
	//	for (cv::FileNodeIterator it = root_node.begin(); it != root_node.end(); it++) {
	//		cvWriteFileNode(xml_out_storage, (*it).name().c_str(), (*it).node, 0);
	//	}
	//}

    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* p = *iter;
        int type = ((part_data*)p->data)->type;
        opart_data* opd = dynamic_cast<opart_data*>(p->data);
        list_t pos, neg;
        string fname;
        unsigned vtype; // verification type

        learn(pos, neg, type);

        cout << "type: " << type << ", |pos|, |neg|: " << pos.size() << ", " << neg.size() << endl;

        if (!pos.empty() && !neg.empty()) vtype = opart_data::SVM;
        else if (pos.empty() && neg.empty()) vtype = opart_data::UNKNOWN;
        else if (pos.empty()) vtype = opart_data::NEGATIVE;
        else vtype = opart_data::POSITIVE;

        if (vtype == opart_data::SVM) {
            svm2 asvm;

            fname = string("type_") + type; //+ string(".xml");
            //asvm.train(fname, pos, neg, xml_out_storage);
        }

        if (opd != nullptr) {
            opd->set_name(fname);
            opd->set_vtype(vtype);
        } else {
            part_data* pd = (part_data*)p->data;

            opd = new opart_data(*pd, fname, vtype);
            p->data = opd;
            delete pd;
        }

        //if (type == 152) {
        //    cout << "Mean: " << amean << endl;
        //    cout << "Positive: ";
        //    for (list_t::iterator liter = pos.begin(); liter != pos.end(); ++liter) {
        //        cout << *liter << "  " << endl;
        //    }
        //    cout << endl;
        //    cout << "Negative: ";
        //    for (list_t::iterator liter = neg.begin(); liter != neg.end(); ++liter) {
        //        cout << *liter << "  " << endl;
        //    }
        //    cout << endl;
        //}

        /// !!!!
        //if (pos.empty() || neg.empty()) {
        //    opd->td.set_thresh(RR_THRESH, 1.0);
        //    opd->td.set_thresh(G_THRESH, 1.0);
        //}
    }
 
	// release existing SVM storage 
	//if (library->svm_models_storage != nullptr) {
	//	cvReleaseFileStorage(&(library->svm_models_storage));
	//}

	//library->svm_models_storage = xml_out_storage;	
}

int ell_learning::get_new_sample_index(bool positive, int type)
{
    index_map_t::iterator iter;
    int result = 0;

    if (positive) {
        if ((iter = posindexmap.find(type)) == posindexmap.end()) posindexmap.insert(iipair(type, 0));
        else result = ++(iter->second);
    } else {
        if ((iter = negindexmap.find(type)) == negindexmap.end()) negindexmap.insert(iipair(type, 0));
        else result = ++(iter->second);
    }
    return result;
}

void ell_learning::add_to_map(bool positive, int type, int index, const ipoint2& name, const ellipse& ell)
{
    if (positive) {
        posstat.insert(map_t::value_type(key_t(type, index, name), ell_learning_stat_item(ell)));
    } else {
        negstat.insert(map_t::value_type(key_t(type, index, name), ell_learning_stat_item(ell)));
    }
}

// t_learning
///////////////////////////////////////////////////////////////////////////////

void get_subpart_values(/*map<pair<int, ipoint2>, double>& result, */ const vector<node*>& nodes)
{
    typedef map<pair<int, ipoint2>, double> result_t;
    typedef vector<node*> src_t;

    int to_prev = atom("toPrevLayer");

    for (src_t::const_iterator niter = nodes.begin(); niter != nodes.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;
        set<int> spset;

        cout << "type: " << nd->m << " (" << nd->r(RR_RESPONSE) << "), len:  " << n->count_neighbors(to_prev) << "; ";
        foreach_neighbor(n, to_prev, niter) {
            layer1_data* nnd = (layer1_data*)neighbor_node_data(niter);
            edge_data_name* nned = (edge_data_name*)neighbor_edge_data(niter);

            if (spset.find(nned->index) == spset.end()) {
                cout << nned->index << " (" << nnd->r(RR_RESPONSE) << ") ";
                spset.insert(nned->index);
            }
            //if (nned != nullptr) {
            //    pair<result_t::iterator, bool> ibpair;
                
                //ibpair = result.insert(result_t::value_type(result_t::key_type(nd->m, nned->data), nnd->r(RR_RESPONSE)));
                //if (!ibpair.second) 
                //    if (ibpair.first->second 

            //}
        }
        cout << endl;

    }
}

void get_values(map<iipair, vector<double> >& val, const vector<node*>& nodes)
{
    typedef map<iipair, vector<double> > result_t;
    typedef vector<node*> src_t;

    for (src_t::const_iterator niter = nodes.begin(); niter != nodes.end(); ++niter) {
        layer1_data* nd = (layer1_data*)(*niter)->data;
        pair<result_t::iterator, bool> ibpair;
        
        // -------- LEARNING --------
        val[iipair(nd->m, G_THRESH)].push_back(nd->r.get_response(G_RESPONSE));
        val[iipair(nd->m, RR_THRESH)].push_back(nd->r.get_response(RR_RESPONSE));
        val[iipair(nd->m, S_THRESH)].push_back(nd->r.get_response(S_RESPONSE, 0.0));
        // --------------------------

    }
}

void collect_subpart_statistics(map<int, map<iipair, int> >& result, const vector<node*>& nodes)
{
    typedef map<iipair, int> result_value_t;
    typedef map<int, result_value_t> result_t;   // type |-> ((depth, edgename.index) |-> count)
    typedef vector<node*> src_t;

    int ename = atom("toPrevLayer");
    
    for (src_t::const_iterator niter = nodes.begin(); niter != nodes.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;
        set<iipair> spset;

        foreach_neighbor(n, ename, nbiter) {
            edge_data_name* nned = (edge_data_name*)neighbor_edge_data(nbiter);
            layer1_data* nnd = (layer1_data*)neighbor_node_data(nbiter);

            if (nned != nullptr) {
                spset.insert(iipair(nnd->z - nd->z, nned->index));
                //result_t::iterator riter = result.insert(result_t::value_type(nd->m, result_value_t())).first;
                //result_value_t::iterator rviter = riter->second.insert(result_value_t::value_type(itriple(nnd->z - nd->z, nned->data.x, nned->data.y), 0)).first;

                //++rviter->second;
            }
        }
        for (set<iipair>::iterator siter = spset.begin(); siter != spset.end(); ++siter) {
            result_t::iterator riter = result.insert(result_t::value_type(nd->m, result_value_t())).first;
            result_value_t::iterator rviter = riter->second.insert(result_value_t::value_type(*siter, 0)).first;

            ++rviter->second;
        }
    }

}

t_learning::t_learning(const config_dictionary &cfg) :
    library(nullptr)
{
    cfg_init(cfg);
}

t_learning::~t_learning()
{
    if (library != nullptr) delete library;
}

void t_learning::cfg_init(const config_dictionary &cfg)
{
    string libname;

    cfg.get_value(category_layer, "category_layer", true);
    use_groundthruth = cfg.get_value_bool("use_groundthruth", true);
    false_is_negative = cfg.get_value_bool("false_is_negative", true);
    groundtruth_threshold = cfg.get_value_double("groundtruth_threshold", 0.7);
    object_layer_only = cfg.get_value_bool("object_layer_only", false);
    cfg.get_value(libname, "library", true);
    tf_ratio = cfg.get_value_double("tf_ratio", 0.75);

	gtextension = cfg.get_value_string("groundtruth_extension",".groundtruth");
    remove_bad_parts = cfg.get_value_bool("remove_bad_parts", true); 

    read_library(libname, library);
}

//void t_learning::get_object_parts(set<int>& result, string& sm)
//{
//    int nbname = atom("lyrPrev");
//
//	// find category node with name sm
//	node* n = nullptr;
//	for (vector<node*>::iterator iter = library->parts[category_layer].begin(); 
//		iter != library->parts[category_layer].end(); ++iter) {
//		n = *iter;
//		cpart_data* nd = (cpart_data*)n->data;
//		if (nd->name == sm) break;
//	}
//	if (n != nullptr) {
//		foreach_neighbor(n, nbname, iter) {
//			lib_data* pd = (lib_data*)neighbor_node_data(iter);
//			result.insert(pd->type);
//		}
//	}
//}

// Returns min threshold t s.t. |pos > t|/|pos + neg > t| >= tf_ratio
// and -1.0 if the highest possible value does not satisfy this condition
double t_learning::calc_obj_threshold(const vector<double>& pos, const vector<double>& neg)
{
    const double step = 0.01;

    if (neg.empty() && pos.empty()) return 1.0;
    if (neg.empty()) return *min_element(pos.begin(), pos.end());
    if (pos.empty()) return *max_element(neg.begin(), neg.end()) + step;

    double maxval = std::max(*max_element(pos.begin(), pos.end()), *max_element(neg.begin(), neg.end()));
    double result;
    geq_predicate<double> pred(maxval);
    int posc, negc;

    posc = count_if(pos.begin(), pos.end(), pred);
    negc = count_if(neg.begin(), neg.end(), pred);
    if (negc > 0 && (double)posc/(posc + negc) < tf_ratio) return -1.0; 
    result = pred.val;
    do {
        pred.val -= step;
        posc = count_if(pos.begin(), pos.end(), pred);
        negc = count_if(neg.begin(), neg.end(), pred);
        if (negc > 0 && (double)posc/(posc + negc) < tf_ratio) break;
        result = pred.val;
    } while (result > 0.0);
    if (result > 0) return result;
    else return *min_element(pos.begin(), pos.end());
}

// Returns min threshold t s.t. |pos > t|/|pos + neg > t| >= tf_ratio
// and -1.0 if the highest possible value does not satisfy this condition
// Version for "negative" thresholds; i.e "<= thresh" is allowed
double t_learning::calc_obj_neg_threshold(const vector<double>& pos, const vector<double>& neg)
{
    const double step = 0.01;
    const double maxvalue = 5.0;

    if (neg.empty() && pos.empty()) return 0.0;
    if (neg.empty()) return *max_element(pos.begin(), pos.end());
    if (pos.empty()) return *min_element(neg.begin(), neg.end()) - step;

    double minval = std::min<double>(*min_element(pos.begin(), pos.end()), *min_element(neg.begin(), neg.end()));
    double result;
    leq_predicate<double> pred(minval);
    int posc, negc;

    posc = count_if(pos.begin(), pos.end(), pred);
    negc = count_if(neg.begin(), neg.end(), pred);
    if (negc > 0 && (double)posc/(posc + negc) < tf_ratio) return -1.0; 
    result = pred.val;
    do {
        pred.val += step;
        posc = count_if(pos.begin(), pos.end(), pred);
        negc = count_if(neg.begin(), neg.end(), pred);
        if (negc > 0 && (double)posc/(posc + negc) < tf_ratio) break;
        result = pred.val;
    } while (result < maxvalue);
    if (result < maxvalue) return result;
    else return *max_element(pos.begin(), pos.end());
}

//double t_learning::calc_obj_threshold(const vector<node*>& pos, const vector<node*>& neg)
//{
//    vector<double> posv, negv;
//
//    for (vector<node*>::const_iterator iter = pos.begin(); iter != pos.end(); ++iter) {
//        layer1_data* nd = (layer1_data*)(*iter)->data;
//        posv.push_back(nd->val);
//    }
//    for (vector<node*>::const_iterator iter = neg.begin(); iter != neg.end(); ++iter) {
//        layer1_data* nd = (layer1_data*)(*iter)->data;
//        negv.push_back(nd->val);
//    }
//    return calc_obj_threshold(pos, neg);
//}

void t_learning::update_tree_thresholds(layer1_result* res, node* n)
{
    layer1_data* nd = (layer1_data*)n->data;
    double thr = library->get_thresh(G_THRESH, nd->z, nd->m, -1.0);

    if (thr < 0 || nd->val() < thr) return;

    int ename = atom("toPrevLayer");
    set<node*> snodes;

    res->subgraph_from_node(n, ename, snodes);
    for (set<node*>::iterator iter = snodes.begin(); iter != snodes.end(); ++iter) {
        node* sn = *iter;
        layer1_data* snd = (layer1_data*)sn->data;

        if (sn == n) {
            library->update_thresh(RR_THRESH, snd->z, snd->m, snd->r.get_response(RR_RESPONSE, 0.99));
        } else {
            library->update_thresh(R_THRESH, snd->z, snd->m, snd->r.get_response(R_RESPONSE, 0.99));
            library->update_thresh(G_THRESH, snd->z, snd->m, snd->r.get_response(G_RESPONSE, 0.99));
            library->update_thresh(RR_THRESH, snd->z, snd->m, snd->r.get_response(RR_RESPONSE, 0.99));
        }
    }
}

void t_learning::update_tree_thresholds(layer1_result* res, const vector<node*>& nodes)
{
    for (vector<node*>::const_iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        update_tree_thresholds(res, *iter);
    }
}


void t_learning::begin()
{
    toremove.clear();
}

void t_learning::end()
{
    if (remove_bad_parts) 
        library->delete_parts(category_layer - 1, vector<int>(toremove.begin(), toremove.end()));
    else 
        cout << "Warning: \"bad\" parts were not removed!" << endl;
}

// 'threshold_name' is one of R_THRESH, G_THRESH, RR_THRESH, S_THRESH, ... 
//    NOTE: *_THRESH should match *_RESPONSE values!
// 'threshold_type' is -1 if threshold is disciminative in negative way, i.e. if lower value is 
//    better (S) and +1 if higher value is better (R, G, RR)
void t_learning::learn(string pdir, const list<string>& positive, 
    string ndir, const list<string>& negative, string& sm, int threshold_name, int threshold_type)
{
    typedef map<iipair, vector<double> > map_t; // (object, thresh_name) -> responses
    typedef map<int, double> thresh_t;
    typedef map<int, map<itriple, int> > substat_t; // see collect_subpart_statistics

    end_dir(pdir);
    end_dir(ndir);
    map_t pvals, nvals;
    set<int> types;
    substat_t substat;

    // Get object parts
    library->get_object_parts(types, category_layer, sm);

    // Pass 1 - learn object layer thresholds
    // ------------------------------------------------------------------------

    cout << endl << "Pass 1 (learning object layer thresholds)" << endl << endl;

    // Read & process positive examples
    cout << "Positive examples" << endl;

    for (list<string>::const_iterator file = positive.begin(); file != positive.end(); ++file) {
        layer1_result* res;

        cout << "Processing " << *file;

        read_layer1_result(res, pdir + *file);
        if (res == nullptr) cout << " error reading file" << endl; 
        else {
            list<irectangle2> rectangles;
            vector<node*> pnodes;

            res->get_nodes(pnodes, category_layer - 1, types);
			if (pnodes.size() > 1000) {
				if (threshold_type > 0) layer1_data::sort_by_response(pnodes, threshold_name);		
                else layer1_data::sort_by_neg_response(pnodes, threshold_name);
				pnodes.resize(1000);
			}
            if (use_groundthruth) 
                read_groundtruth(rectangles, pdir + *file, sm);
            vector<node*> negatives = res->filter_with_groundtruth(pnodes, rectangles, groundtruth_threshold);
			if (false_is_negative) {
                //cout << "Negative:" << endl;
                //cout << "---------" << endl;
                get_values(nvals, negatives);
                //get_subpart_values(negatives);
			}
            //cout << "Positive:" << endl;
            //cout << "---------" << endl;
            get_values(pvals, pnodes);
            //get_subpart_values(pnodes);
            //collect_subpart_statistics(substat, pnodes);

            cout << " done" << endl;
            delete res;
        }
    }

    // Print "substat"
    //for (substat_t::iterator ssiter = substat.begin(); ssiter != substat.end(); ++ssiter) {
    //    cout << ssiter->first << ": ";
    //    for (substat_t::mapped_type::iterator miter = ssiter->second.begin(); miter != ssiter->second.end(); ++miter) {
    //        cout << '(' << miter->first.first << ',' << miter->first.second << ',' << miter->first.third << "): ";
    //        cout << miter->second << endl;
    //    }
    //}

    // Read & process negative examples
    cout << "Negative examples" << endl;

    for (list<string>::const_iterator file = negative.begin(); file != negative.end(); ++file) {
        layer1_result* res;

        cout << "Processing " << *file;

        read_layer1_result(res, ndir + *file);
        if (res == nullptr) cout << " error reading file" << endl; 
        else {
            vector<node*> nnodes;

            res->get_nodes(nnodes, category_layer - 1, types);
            get_values(nvals, nnodes);

            cout << " done" << endl;
            delete res;
        }
    }

    // Calc object thresholds and update (overwrite) the values in the library
    list<map_t::key_type> keys;

    for (map_t::iterator piter = pvals.begin(); piter != pvals.end(); ++piter) 
        keys.push_back(piter->first);
    for (map_t::iterator niter = nvals.begin(); niter != nvals.end(); ++niter)
        keys.push_back(niter->first);

    for (list<map_t::key_type>::iterator kiter = keys.begin(); kiter != keys.end(); ++kiter) {

        if (kiter->second != threshold_name) 
            continue;

        int m = kiter->first;
        map_t::iterator niter = nvals.find(iipair(m, kiter->second));
        map_t::iterator piter = pvals.find(iipair(m, kiter->second));
        double thr;

        if (threshold_type > 0)
            thr = calc_obj_threshold(
                (piter != pvals.end()) ? piter->second : vector<double>(),
                (niter != nvals.end()) ? niter->second : vector<double>()
            );
        else 
            thr = calc_obj_neg_threshold(
                (piter != pvals.end()) ? piter->second : vector<double>(),
                (niter != nvals.end()) ? niter->second : vector<double>()
            );

        //cout << piter->second << endl;
        //if (niter != nvals.end()) cout << niter->second << endl;
        //cout << thr << endl;
        //cout << endl;

        if (thr == -1.0) {
            thr = threshold_type > 0 ? 1.0 : 0.0;
            toremove.insert(m);
        }
        cout << "m = " << m << ": " << thr << endl;
        library->set_thresh(kiter->second, category_layer - 1, m, thr);
    }

    // Pass 2 - learn (update) part thresholds
    // ------------------------------------------------------------------------

    set<int> alltypes(types);

    types.clear();
    set_difference(types, alltypes, toremove);

    if (!object_layer_only) {
        cout << endl << "Pass 2 (learning part thresholds)" << endl << endl;

        // Read & process positive examples
        for (list<string>::const_iterator file = positive.begin(); file != positive.end(); ++file) {
            layer1_result* res;

            cout << "Processing " << *file;

            read_layer1_result(res, pdir + *file);
            if (res == nullptr) cout << " error reading file" << endl; 
            else {
                list<irectangle2> rectangles;
                vector<node*> pnodes;

                //if (use_groundthruth) read_groundtruth(rectangles, pdir + *file);
                res->get_nodes(pnodes, category_layer - 1, types);
                update_tree_thresholds(res, pnodes);
                
                //    For each node n in pnodes above a sufficient theshold
                //       For each child node cn in subtree of n
                //          update thresh of (cn->z, cn->m) to cn->value

                cout << " done" << endl;
                delete res;
            }
        }
    }

    // Display parts and their thresholds on object layer
    // ------------------------------------------------------------------------

    cout << endl;
    library->display_thresholds(library->parts[category_layer - 1]);

}


// Keeps only nodes with area(projection_bb ^ gtr)/area(projection_bb u gtr) >= thresh for some
// rectangle in gtr. Returns the nodes for which this is not true.
void filter_with_groundtruths(vector<pair<node*, int> >& positive, vector<node*>& negative,
    vector<node*>& nodes, layer1_result* res, const list<irectangle2>& gtr, double thresh)
{
    vector<node*> result;

    if (gtr.empty()) {
        negative = nodes;
        return;
    }

    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter) {
        irectangle2 box = res->get_box(*iter);

        if (box.invalid()) {
            negative.push_back(*iter);
        } else {
            bool inside = false;
            int i = 0;

            for (list<irectangle2>::const_iterator biter = gtr.begin(); !inside && biter != gtr.end(); ++biter, ++i) {
                double r = (double)box.intersection(*biter).area()/box.union_area(*biter);
                inside = r >= thresh;
            }
            if (!inside) negative.push_back(*iter);
            else positive.push_back(pair<node*, int>(*iter, i));
        }
    }
}

vector<pair<node*, int> > filter_nodes(const vector<pair<node*, int> >& nv, const response_filter& rsf)
{
	vector<pair<node*, int> > result;

    for (auto niter = nv.begin(); niter != nv.end(); ++niter) {
        node* n = niter->first;
        layer1_data* nd = (layer1_data*)n->data;

        if (rsf.check(nd)) 
            result.push_back(*niter);
    }
	return result;
}

void update_negative_value(vector<double>& v, int i, double val)
{
    if (i < v.size()) {
        if (v[i] < val) v[i] = val;
    } else {
        v.resize(i, -1.0);
        v[i] = val;
    }
    
}

vector<irectangle2> transform_boxes(const vector<irectangle2>& r, double factor, int border)
{
    vector<irectangle2> result;

    for (auto riter = r.begin(); riter != r.end(); ++riter) {
        result.push_back(*riter);
        result.back().ur *= factor;
        result.back().ll *= factor;
        result.back().ur += border;
        result.back().ll += border;
    }
    return result;
}

template<class T> void sort_by_response(T begin, T end, int response, bool ascending)
{
    if (ascending) {
        sort(begin, end, [response](node* n, node* m) -> bool { 
                layer1_data* nd = (layer1_data*)n->data, * md = (layer1_data*)m->data;
                return nd->r(response) < md->r(response);
            }
        );
    } else {
        sort(begin, end, [response](node* n, node* m) -> bool { 
                layer1_data* nd = (layer1_data*)n->data, * md = (layer1_data*)m->data;
                return nd->r(response) > md->r(response);
            }
        );
    }
}

template<class T> void trim(vector<T>& v, size_t newmaxsize)
{
    v.resize(min(v.size(), newmaxsize));
}


template<class T> double vector_diff(const vector<point2<T> >& v1, const vector<point2<T> >& v2)
{
    size_t maxi = std::min<size_t>(v1.size(), v2.size());

    if (maxi == 0) return numeric_limits<double>::infinity();

    double val = 0;

    for (size_t i = 0; i < maxi; ++i) {
        val += sqrt((double)(v1[i].distance2(v2[i])));
    }
    return val/maxi;
}

// pca_merging
///////////////////////////////////////////////////////////////////////////////

pca_merging::pca_merging(int l, int r) : stat()
{
    layer = l;
    iradius = r;
}

void pca_merging::update_stat(layer1_result* res, const response_filter& filter, double area_thresh)
{
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");

    if (res == nullptr || layer < 0 || layer > res->max_layer_index()) return;

    vector<node*>& s_nodes = res->shape_nodes[layer];
    vector<node*>& s_nodes0 = res->shape_nodes[0];
    set<ipoint2> allpts0;

    node_set_to_point_set(allpts0, s_nodes0.begin(), s_nodes0.end());

    irectangle2 box = irectangle2::bounding_rectangle(allpts0.begin(), allpts0.end());

    for (vector<node*>::iterator siter = s_nodes.begin(); siter != s_nodes.end(); ++siter) {
        node* n = *siter;

        while (n != nullptr) {
            layer1_data* nd = (layer1_data*)n->data;

            if (filter.check(nd)) {
                //set<node*> ly0n;
                //set<ipoint2> ly0p;

                //res->recurse_and_link(n, toprev, to0, ly0n);
                //node_set_to_point_set(ly0p, ly0n.begin(), ly0n.end());

                vector<pair<int, ipoint2> > pts; // = inhibit_point_set(vector<ipoint2>(ly0p.begin(), ly0p.end()), iradius);

                get_node_geo_p(pts, res, n);
                pts = inhibit_point_set(pts, nd->z);
                translate_and_scale(pts);

                vector<ipoint2> ptsvec = extract_second<ipoint2>(pts.begin(), pts.end());
                
                if (filled_region_size(ptsvec) > area_thresh*box.area()) {
                    stat[nd->m].update(pts);
                    cout << '<' << nd->m << '>';
                }
            }
            n = nd->next;
        }


    }
}

// We assume that part field in each item of 'cliques' vector has, initially, exactly one 
// element in it.
void reduce_pca_set(vector<clique_data>& cliques, double norm)
{
    set<int> partset;

    for (vector<clique_data>::const_iterator cliter = cliques.begin(); cliter != cliques.end(); ++cliter) 
        partset.insert(cliter->parts.begin(), cliter->parts.end());

    if (partset.empty()) 
        return;

    int maxpart = *max_element(partset.begin(), partset.end());

    vector<vector<dpoint2> > geovec(maxpart + 1, vector<dpoint2>());

    // 1. Get part geometry (from mean) for each part
    for (vector<clique_data>::const_iterator cliter = cliques.begin(); cliter != cliques.end(); ++cliter) 
        geovec[*cliter->parts.begin()] = partition(cliter->mean);

    // 2. Calculate norm of part geometry to its "back"-projection from clique space -- for each part and each clique
    //    matrix(p, c)
    matrix<double> cliquedist(maxpart + 1, (int)cliques.size());

    for (int pi = 0; pi <= maxpart; ++pi) {
        vector<dpoint2> gv = geovec[pi];
        cout << '.';

        if (gv.empty()) 
            continue;

        for (int ci = 0; ci < (int)cliques.size(); ++ci) {
            vector<dpoint2> v = partition(cliques[ci].mean);
            vector<dpoint2> gvresized = get_resized_vector(gv, (int)v.size());

            translate_and_scale(gvresized);

            vector<int> perm = point_matching(gvresized, v);

            //write_matching(string("c:\\work\\tmp\\match") + pi + string("-") + ci + string(".m"), gvresized, v);
            
            permute(gvresized, perm);
            
            cv::Mat data = flatten(gvresized);
            cv::Mat coeffs;
            cv::Mat bproj;
            cv::Mat result;
            
            gemm(data - cliques[ci].mean, cliques[ci].eigenvectors, 1, cv::Mat(), 0, coeffs, cv::GEMM_2_T);
            gemm(coeffs, cliques[ci].eigenvectors, 1, cliques[ci].mean, 1, result, 0);
            cliquedist(pi, ci) = cv::norm(result, data, cv::NORM_L2);
        }
    }
    cout << endl;

    // Remove redundant cliques based on difference on back-projection
    vector<int> toremove;

    while (true) {

        // 3. Clique D is absorbed in clique C if D can be expressed with basis of C with small error
        vector<vector<int> > nabsorbed(cliques.size(), vector<int>());   // nabsorbed[i] counts cliques absorbed in clique i

        for (int ci = 0; ci < (int)cliques.size(); ++ci) {
            if (cliques[ci].parts.empty())
                continue;

            for (int di = 0; di < (int)cliques.size(); ++di) {
                const set<int>& D = cliques[di].parts;
                double maxdiff = 0.0;

                if (D.empty())
                    continue;

                for (set<int>::const_iterator iter = D.begin(); iter != D.end(); ++iter) {
                    int e = *iter;
                    double diff = cliquedist(e, ci);

                    if (diff > maxdiff) maxdiff = diff;
                }
                //cout << di << " is absorbed in " << ci << " with " << maxdiff << endl;
                if (ci != di && maxdiff <= norm) 
                    nabsorbed[ci].push_back(di);
                
            }
        }

        cout << "na: " << nabsorbed.size() << endl;

        // 4. Find clique which absorbes most cliques
        vector<vector<int> >::iterator miter = max_element(nabsorbed.begin(), nabsorbed.end(), size_less<int>);

        // 4a. If there is no clique which absorbes some other clique, quit
        if (miter == nabsorbed.end() || miter->size() == 0)
            break;

        // 5. Keep such clique and remove absorbed cliques
        int keepi = (int)(miter - nabsorbed.begin());
        vector<int>& torm = *miter;

        for (vector<int>::iterator riter = torm.begin(); riter != torm.end(); ++riter) {
            cliques[keepi].parts.insert(cliques[*riter].parts.begin(), cliques[*riter].parts.end());
            toremove.push_back(*riter);
            cliques[*riter].parts.clear();
        }
    }
    sort(toremove.begin(), toremove.end(), greater<int>());
    for (vector<int>::iterator riter = toremove.begin(); riter != toremove.end(); ++riter) {
        cliques.erase(cliques.begin() + *riter);
    }
}

void reduce_pca_set(part_lib* library, int layer, double norm)
{
    if (library == nullptr || layer < 0 || layer > library->max_layer_index()) 
        return;

    vector<node*>& parts = library->parts[layer];
    vector<clique_data> cliques(parts.size());

    for (int i = 0; i < (int)parts.size(); ++i) {
        node* p = parts[i];
        vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);

        if (pd != nullptr) {
            cliques[i].parts.insert(pd->type);
            cliques[i].mean = pd->pcad.mean;
            cliques[i].eigenvectors = pd->pcad.eigenvectors;
        }
    }
    reduce_pca_set(cliques, 1.0);

    for (int i = 0; i < (int)cliques.size(); ++i) {
        for (set<int>::iterator ci = cliques[i].parts.begin(); ci != cliques[i].parts.end(); ++ci) {
            cout << *ci << ' ';
        }
        cout << endl;
    }

}

vector<pair<int, ipoint2> > get_scaled_point_set(const vector<pair<int, ipoint2> >& ipts, double factor)
{
    vector<pair<int, ipoint2> > result;

    result.reserve(ipts.size());
    for (int i = 0; i < (int)ipts.size(); ++i) {
        const ipoint2& p = ipts[i].second;
        ipoint2 sp(int_round(p.x*factor), int_round(p.y*factor));

        result.push_back(pair<int, ipoint2>(ipts[i].first, sp));
    }
    return result;
}

void pca_merging::update_library(part_lib* library, double mergenorm)
{
    int rootname = atom("lyrSimRoot");
    int simname = atom("lyrSimilar");

    if (library == nullptr || layer < 0 || layer > library->max_layer_index()) 
        return;

    vector<node*>& parts = library->parts[layer];

    // Clear all lyrSimRoot and lyrSimilar edges
    for (vector<node*>::iterator piter = parts.begin(); piter != parts.end(); ++piter) {
        node* p = *piter;

        p->delete_edges(rootname);
        p->delete_edges(simname);
    }

    // Update clique data from stat
    cout << "Updating clique data from stat." << endl;

    vector<clique_data> cliques;

    cliques.reserve(stat.size());
    for (statmap_t::iterator siter = stat.begin(); siter != stat.end(); ++siter) {
        cout << '.';

        int type = siter->first;
        vector<pair<int, ipoint2> > ipts = get_library_geo_pieces(parts[type], layer);
        double factor = translate_and_scale(ipts).first;

        pca_data pcad = siter->second.get_pca_data(ipts, factor, 5); // max 5 components  

        if (pcad.mean.cols > 0) {
            cliques.push_back(clique_data());
            cliques.back().parts.insert(type);
            cliques.back().part0 = type;
            cliques.back().mean = pcad.mean;
            cliques.back().eigenvectors = pcad.eigenvectors;

            vs_part_data* pvsd = dynamic_cast<vs_part_data*>(parts[type]->data);

            if (pvsd == nullptr) {
                part_data* pd = (part_data*)parts[type]->data;

                pvsd = new vs_part_data(*pd);
                delete pd;
                parts[type]->data = pvsd;
            } 
            pvsd->pcad = pcad;
        }

        svm_data svmd = siter->second.get_svm_data(ipts);

        if (svmd.svm != nullptr) {
            vs_part_data* pvsd = dynamic_cast<vs_part_data*>(parts[type]->data);

            if (pvsd == nullptr) {
                part_data* pd = (part_data*)parts[type]->data;

                pvsd = new vs_part_data(*pd);
                delete pd;
                parts[type]->data = pvsd;
            }
            pvsd->svmd = svmd;
        } else {    // 

        }

    }

    // Add trivial pca models to other parts
    for (int i = 0; i < (int)parts.size(); ++i) {
        vs_part_data* pvsd = dynamic_cast<vs_part_data*>(parts[i]->data);

        if (pvsd == nullptr) {
            vector<pair<int, ipoint2> > ipts = get_library_geo_pieces(parts[i], layer);
            part_data* pd = (part_data*)parts[i]->data;
            pvsd = new vs_part_data(*pd);
            shape_learning_p pcal;
            double sizefactor;
            
            sizefactor = translate_and_scale(ipts).first;

            vector<ipoint2> pts = extract_second<ipoint2>(ipts.begin(), ipts.end());
            vector<dpoint2> dpts = cast_vector<dpoint2, ipoint2>(pts);

            translate_and_scale(dpts);

            vector<dpoint2> dpts2 = dpts;

            mult_vector(dpts2, 1.2);
            dpts2 = dpts2 - dpts;
            
            pvsd->pcad.mean = flatten(dpts);
            pvsd->pcad.eigenvectors = flatten(dpts2);
            pvsd->pcad.eigenvalues = cv::Mat(1, 1, CV_64F, cv::Scalar(0.15));
            pvsd->pcad.sizefactor = sizefactor;

            delete pd;
            parts[i]->data = pvsd;
        }
    }

    // Merge cliques

    cout << endl;
    //cout << "Merging cliques." << endl;

    //reduce_pca_set(cliques, mergenorm);

    // Print result
    for (int i = 0; i < (int)cliques.size(); ++i) {
        cout << cliques[i].part0 << ": ";
        for (set<int>::iterator ci = cliques[i].parts.begin(); ci != cliques[i].parts.end(); ++ci) {
            cout << *ci << ' ';
        }
        cout << endl;
    }

    // Update library
    for (int i = 0; i < (int)cliques.size(); ++i) {
        node* p0 = parts[cliques[i].part0];

        for (set<int>::iterator ci = cliques[i].parts.begin(); ci != cliques[i].parts.end(); ++ci) {
            node* p = parts[*ci];

            library->add_edge_2(p, p0, new part_data_sim(1.0), nullptr, rootname, simname);
        }
    }

}

// distribution2
///////////////////////////////////////////////////////////////////////////////

pair<ipoint2, double> distribution2::get_peak() 
{
    double bestval = (numeric_limits<double>::min)();
    ipoint2 bestp;

    for (map_t::iterator miter = m.begin(); miter != m.end(); ++miter) {
        if (miter->second > bestval) { 
            bestval = miter->second; 
            bestp = miter->first; 
        }
    }
    return pair<ipoint2, double>(bestp, bestval);
}

rmatrix distribution2::cut_matrix(const ipoint2& p, int delta, double defval /* = 0.0 */)
{
    rmatrix result(2*delta + 1, 2*delta + 1, defval);
    int xmin = p.x - delta, xmax = p.x + delta;
    int ymin = p.y - delta, ymax = p.y + delta;
    
    for (map_t::iterator miter = m.begin(); miter != m.end(); ++miter) {
        const ipoint2& mp = miter->first;

        if (mp.x >= xmin && mp.x <= xmax && mp.y >= ymin && mp.y <= ymax)
            result(mp.x - xmin, mp.y - ymin) = miter->second;
    }
    return result;
}

rmatrix distribution2::to_matrix(double defval /* = 0.0 */)
{
    irectangle2 box;

    for (map_t::iterator miter = m.begin(); miter != m.end(); ++miter) 
        box.eat(miter->first);
    if (box.invalid()) return rmatrix();

    rmatrix result(box.x_dim() + 1, box.y_dim() + 1, defval);

    for (map_t::iterator miter = m.begin(); miter != m.end(); ++miter)
        result(miter->first.x - box.ll.x, miter->first.y - box.ll.y) = miter->second;
    return result;
}

// global exported functions
///////////////////////////////////////////////////////////////////////////////


// Result type:
//   int (part type) -> (ipoint2 (subpart id) -> online_distribution)
//   == (int, ipoint2) -> online_distribution
// Returns true if all s_nodes 
void update_g_distribution(map<pair<int, int>, online_distribution>& result, layer1_result* res, int layer, part_lib* library)
{
    typedef map<int, double> prmap_t;
    typedef map<pair<int, int>, online_distribution> result_t;

    vector<node*>& s_nodes = res->shape_nodes[layer];
    int srcname = atom("lyrSrc");
    int centername = atom("lyrCenter");
    int prevname = atom("toPrevLayer");

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;
        node* p = library->parts[nd->z][nd->m];
        prmap_t prmap; // part -> best response

        //  Fill prmap (part-best-response map)
        foreach_neighbor(n, prevname, niter) {
            layer1_data* nnd = (layer1_data*)neighbor_node_data(niter);
            edge_data_name* nned = (edge_data_name*)neighbor_edge_data(niter);

            if (nned == nullptr) 
                continue; //throw exception("Edge data not found; use add_edge_names = true");

            pair<prmap_t::iterator, bool> ibpair = prmap.insert(prmap_t::value_type(nned->index, nnd->r(R_RESPONSE))); // G_RESPONSE

            if (!ibpair.second && ibpair.first->second < nnd->r(R_RESPONSE)) ibpair.first->second = nnd->r(R_RESPONSE); // R_RESPONSE
        }

        // Check if center & all subparts are present
        if (p->get_neighbor(centername) != nullptr && prmap.find(0) != prmap.end()) {
            bool ok = true;

            // Check if all subparts are present
            foreach_neighbor(p, srcname, piter) {
                part_data_2* pnd = (part_data_2*)neighbor_edge_data(piter);
                if (prmap.find(pnd->index) == prmap.end()) {
                    ok = false;
                    break;
                }
            }
            if (ok) {
                double cresp = prmap.find(0)->second;

                for (prmap_t::iterator priter = prmap.begin(); priter != prmap.end(); ++priter) {
                    if (priter->first != 0) {
                        pair<result_t::iterator, bool> ibpair = 
                            result.insert(result_t::value_type(result_t::key_type(nd->m, priter->first), online_distribution()));

                        ibpair.first->second.new_data(log(priter->second) - log(cresp));
                    }
                }
            }
        }
    }

}

// Set g-distribution to parts in the library; it considers them set if variance is not 0.0.
void set_g_distribution(part_lib* library, int layer, map<pair<int, int>, online_distribution>& result,
    double minvar, bool overwrite)
{
    typedef map<pair<int, int>, online_distribution> result_t;

    int srcname = atom("lyrSrc");

    for (result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
        node* p1 = library->parts[layer][iter->first.first];
        pair<node*, part_data_2*> p2p = library->get_neighbor_pair(p1, iter->first.second);

        if (p2p.second->gdistr.second == 0.0 || overwrite) {
            p2p.second->gdistr.first = iter->second.get_mean();
            p2p.second->gdistr.second = max(iter->second.get_variance(), minvar);
        }
        //cout << '(' << p2p.second->gdistr.first << ',' << p2p.second->gdistr.second << ')' << endl;
    }
}

// srccfg is config dictionary where namespace is used to make configuration
// for the inference process; source_layer is layer index of the part_learning
void learn_g_distributions(part_lib* library, int source_layer, const config_dictionary& srccfg, const string& nspace,
    const string& dir, const string& patt)
{
    typedef map<pair<int, int>, online_distribution> dist_map_t;

    config_dictionary icfg; // inference configuration
    
    icfg.from_namespace(srccfg, nspace); // fill from srccfg.namespace 
    icfg.from_string( // set/override some keys
        "realization_ratio_threshold = 0.9;" // !?
        "add_edge_names = true;"
        "layer_contraction = 1.0;"
        "identity_g_response = true"
    );

    layern_creator creator(icfg);
    list<string> files;
    dist_map_t distmap;
    int layer = source_layer + 1;

    // Gather statistics
    creator.set_library(library);
    list_directory(files, dir + patt);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;
    
        cout << "Processing: " << *file;
        read_layer1_result(res, dir + *file);
        if (res == nullptr) 
            cout << " error reading file" << endl;
        else {
            cout << " <I";

            creator.add_layer(res, layer + 1, 0);

            cout << "><G";

            update_g_distribution(distmap, res, layer, library);

            cout << '>';

            delete res;
            cout << " done" << endl;
        }
    }

    // Set distributions
    set_g_distribution(library, layer, distmap, icfg.get_value_double("min_variance", 0.1), false);

    // Prevent library from being disposed by creator!
    creator.set_library(nullptr);
}

// 
void get_layer0_overlaps(list<pair<double, node*> >& ovector, layer1_result* res, const map<node*, set<node*> >& rmap, 
    node* n, int radius)
{
    typedef pair<double, node*> item_t;
    typedef map<int, pair<node*, double> > map_t;

    layer1_data* nd = (layer1_data*)n->data;

    if (!res->grid(nd->z)) res->init_grid(nd->z);

    const set<node*>& nrec = rmap.find(n)->second;
    map_t result;

    int mini = std::max<int>(0, nd->x - radius), maxi = std::min<int>(nd->x + radius + 1, res->x_size());
    int minj = std::max<int>(0, nd->y - radius), maxj = std::min<int>(nd->y + radius + 1, res->y_size());

    ovector.clear();
    for (int i = mini; i < maxi; ++i) {
        for (int j = minj; j < maxj; ++j) {
            node* nn = res->node_at(i, j, nd->z);

            while (nn != nullptr) {
                layer1_data* nnd = (layer1_data*)nn->data;
                const set<node*>& nnrec = rmap.find(nn)->second;

                double f = (double)intersection_size(nrec, nnrec)/nrec.size();

                if (f > 0.0) {
                    pair<map_t::iterator, bool> ibpair = result.insert(map_t::value_type(nnd->m, pair<node*, double>(nn, f)));

                    if (!ibpair.second && ibpair.first->second.second < f) {
                        ibpair.first->second.first = nn;
                        ibpair.first->second.second = f;
                    }
                }
                nn = nnd->next;
            }
        }
    }
    for (map_t::iterator iter = result.begin(); iter != result.end(); ++iter) 
        ovector.push_back(item_t(iter->second.second, iter->second.first));
    ovector.sort(greater<item_t>());
}

void update_mean_map(map<int, map<int, pair<distribution2, histogram> > >& result, 
    layer1_result* res, int layer, part_lib* library, double rthresh, double gthresh)
{
    typedef pair<distribution2, histogram> data_t;
    typedef map<int, data_t> subpart_map_t;
    typedef map<int, subpart_map_t> result_t;
    typedef result_t::mapped_type result_mt;
    typedef map<node*, set<node*> > reconstruction_map_t;
    
    vector<node*>& s_nodes = res->shape_nodes[layer];
    reconstruction_map_t rmap;
    int prevname = atom("toPrevLayer");
    int srcname = atom("lyrSrc");
    int centername = atom("lyrCenter");
    double factor = (double)res->x_size(layer)/res->x_size(layer - 1);

    if (factor < 0.9 || factor > 1.1) {
        cout << "Warning: Contraction factor is not 1.0. No maps will be updated." << endl;
        return;
    }

    //cout << 'R';
    get_reconstruction_map(rmap, res, layer - 1);
    //cout << 'S';
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        do {
            layer1_data* nd = (layer1_data*)n->data;
            node* part = library->parts[nd->z][nd->m];
            part_data* partd = dynamic_cast<part_data*>(part->data);
            int cmx = 0, cmy = 0;
            
            if (partd != nullptr) { cmx = partd->cmx; cmy = partd->cmy; }
            if (nd->r.get_response(R_RESPONSE) >= rthresh || nd->r.get_response(G_RESPONSE) >= gthresh) {
                subpart_map_t ppmap;

                foreach_neighbor(n, prevname, niter) {
                    node* nn = neighbor_node(niter);
                    layer1_data* nnd = (layer1_data*)neighbor_node_data(niter);
                    edge_data_name* nned = (edge_data_name*)neighbor_edge_data(niter);

                    if (nned == nullptr) 
                        continue; //throw exception("Edge data not found; use add_edge_names = true");
                    subpart_map_t::iterator iiter = 
                        ppmap.insert(subpart_map_t::value_type(nned->index, data_t(distribution2(), histogram()))).first;

                    // update distribution of subpart
                    iiter->second.first.update(ipoint2((int)(nnd->x - (nd->x - cmx)), (int)(nnd->y - (nd->y - cmy))), 
                        nned->r);

                    // update appearance of subpart
                    iiter->second.second.update(nnd->m, nned->r);

                    // iterate over all parts which have sufficient overlap with nn
                    list<pair<double, node*> > overlaps;
                    
                    get_layer0_overlaps(overlaps, res, rmap, nn, 2); // THRESHOLD!!!!!
                    for (list<pair<double, node*> >::iterator oiter = overlaps.begin(); oiter != overlaps.end(); ++oiter) {
                        if (oiter->first < 0.6) break;
                        else {
                            layer1_data* nnnd = (layer1_data*)oiter->second->data;
                            iiter->second.second.update(nnnd->m, oiter->first * nnnd->val());
                        }
                    }
                    
                }

                node* p = library->parts[nd->z][nd->m];

                if (p->get_neighbor(centername) != nullptr && ppmap.find(0) != ppmap.end()) {
                    bool ok = true;

                    // Check if all subparts are present
                    foreach_neighbor(p, srcname, piter) {
                        part_data_2* pnd = (part_data_2*)neighbor_edge_data(piter);
                        if (ppmap.find(pnd->index) == ppmap.end()) {
                            ok = false;
                            break;
                        }
                    }

                    if (ok) {
                        //dpoint2 centerp = ppmap.find(ipoint2::zero)->second.get_mean(); // Hopefully (0, 0)
                        result_t::iterator riter =
                            result.insert(result_t::value_type(nd->m, result_mt())).first;

                        for (subpart_map_t::iterator ppiter = ppmap.begin(); ppiter != ppmap.end(); ++ppiter) {
                            result_mt::iterator miter = 
                                riter->second.insert(result_mt::value_type(ppiter->first, result_mt::mapped_type())).first;
                            miter->second.first.update(ppiter->second.first);
                            miter->second.second.update(ppiter->second.second);
                        }
                    }
                }
            }
            n = nd->next;
            // break;
        } while (n != nullptr);
    }
}

void set_appearance(map<int, pair<vector<int>, double> >& app, int appindex, const histogram& h, double thresh)
{
    histogram::map_t::const_iterator hiter;
    double maxval = 0.0;

    for (hiter = h.m.begin(); hiter != h.m.end(); ++hiter) {
        if (hiter->second > maxval) maxval = hiter->second;
    }

    if (maxval == 0.0) {
        cout << "Warning: histogram has no entries. Skipping update of appearance" << endl;
        return;
    }

    app.clear();
    for (hiter = h.m.begin(); hiter != h.m.end(); ++hiter) {
        double r = hiter->second/maxval;

        if (r > thresh) 
            app.insert(pair<int, pair<vector<int>, double> >(hiter->first, pair<vector<int>, double>(vector<int>(), r)));
    }
}

void update_means(map<int, map<int, pair<distribution2, histogram> > >& result, int layer, part_lib* library)
{
    typedef map<int, map<int, pair<distribution2, histogram> > > result_t;

    int srcname = atom("lyrSrc");
    int centername = atom("lyrCenter");

    for (result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
        node* p = library->parts[layer][iter->first];

        foreach_neighbor(p, srcname, piter) {
            part_data_2* pnd = (part_data_2*)neighbor_edge_data(piter);
            map<int, pair<distribution2, histogram> >::iterator fiter = 
                iter->second.find(pnd->index);

            if (fiter != iter->second.end()) {
                // update new position (map missing!!!)
                ipoint2 mean = fiter->second.first.get_peak().first;

                pnd->x = int_round(mean.x);
                pnd->y = int_round(mean.y);

                // update appearance
				set_appearance(pnd->app, pnd->index, fiter->second.second, 0.2);
            }
            
        }
    }
}

void select_nodes_by_type(vector<node*>& result, layer1_result* res, const vector<node*>& nodes, 
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

struct nodes_at_sort_f {
    int nbname;

    nodes_at_sort_f() : nbname(atom("toLayer0")) { }
    bool operator()(node* n, node* m) const
    {
        layer1_data* nd = (layer1_data*)n->data;
        layer1_data* md = (layer1_data*)m->data;

        return (0.9*nd->r(G_RESPONSE) + 0.1*nd->r(R_RESPONSE)) > (0.9*md->r(G_RESPONSE) + 0.1*md->r(R_RESPONSE));

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

void centralize_set(set<ipoint2>& ps)
{
    set<ipoint2> tmp;
    ipoint2 center(0, 0);

    if (ps.size() < 2) return;
    for (set<ipoint2>::iterator iter = ps.begin(); iter != ps.end(); ++iter) 
        center += *iter;
    center.set(int_round((double)center.x/ps.size()), int_round((double)center.y/ps.size()));
    for (set<ipoint2>::iterator iter = ps.begin(); iter != ps.end(); ++iter) 
        tmp.insert(*iter - center);
    ps = tmp;
}

// functions
///////////////////////////////////////////////////////////////////////////////

struct functions {
    static double BesselI0_approx(double x);
};

double functions::BesselI0_approx(double x)
{
    if (x <= 0.0) return 1;
    else return exp(x)/sqrt(2*M_PI*x);
}


// von_mises_online_distribution
///////////////////////////////////////////////////////////////////////////////

struct von_mises_online_distribution {
protected:
    double csum, ssum;
    int n;
public:
    von_mises_online_distribution() : csum(0.0), ssum(0.0), n(0) { }

    void new_data(double a) { csum += cos(a); ssum += sin(a); ++n; }
    int samples() { return n; }
    void reset();
    pair<double, double> get_parameters() const;

protected:
    static double f6(double x);
};

// Approximate inverse of I1(x)/I0(x)
// Taken from: 
//   Annette J. Dobson, Approximations for the von Mises Concentration Statistic,
//   Journal of the Royal Statistical Society. Series C (Applied Statistics), Vol. 27, No. 3(1978), pp. 345-347
double von_mises_online_distribution::f6(double x)
{
    if (x >= 0.65) {
        return (9 - 8*x + 3*x*x)/8/(1 - x);
    } else {
        double x2 = x*x;
        double x3 = x2*x;
        return 2*x + x3 + 5*x2*x3/6;
    } 
}

void von_mises_online_distribution::reset() 
{
    csum = ssum = 0.0;
    n = 0;
}

// Returns estimated mu and kappa.
pair<double, double> von_mises_online_distribution::get_parameters() const
{
    if (n == 0) return pair<double, double>(0, 0);

    double c = csum/n;
    double s = ssum/n;
    double R = sqrt(c*c + s*s);

    return pair<double, double>(atan2(s, c), f6(R));
}

// ellipse_learner
///////////////////////////////////////////////////////////////////////////////

class ellipse_comparer;

class ellipse_learner {
protected:
    typedef set<ipoint2> pset_t;

    struct online_ell_parameters {
        von_mises_online_distribution ell_angle;
        online_distribution ell_a;
        online_distribution ell_e;

        online_ell_parameters() : ell_angle(), ell_a(), ell_e() { }
    };

    vector<pset_t> ly0sets; // "filter image"
    vector<online_ell_parameters> ell_stat_learn; // statistics

public:
    ellipse_learner(part_lib* library, int layer);
    
    void new_data(layer1_result* res, node* nn);
    void new_data(set<ipoint2>& prset, layer1_result* res, node* nn);
    void print();

    friend class ellipse_comparer;
};

void ellipse_learner::print()
{
    for (vector<online_ell_parameters>::iterator iter = ell_stat_learn.begin(); iter != ell_stat_learn.end(); ++iter) {
        cout << "angle: (" << iter->ell_angle.get_parameters().first << ", " << iter->ell_angle.get_parameters().second << "); ";
        cout << "a: (" << iter->ell_a.get_mean() << ", " << iter->ell_a.get_variance() << "); ";
        cout << "e: (" << iter->ell_e.get_mean() << ", " << iter->ell_e.get_variance() << endl;
    }
}

ellipse_learner::ellipse_learner(part_lib *library, int layer) :
    ell_stat_learn(library->layer_size(layer), online_ell_parameters())
{
    library->get_regions(1, ly0sets);
}

void ellipse_learner::new_data(layer1_result* res, node* nn)
{
    set<ipoint2> prset;

    new_data(prset, res, nn);
}

void ellipse_learner::new_data(set<ipoint2>& prset, layer1_result* res, node* nn)
{
    int tolayer0 = atom("toLayer0");
    layer1_data* nnd = (layer1_data*)nn->data;
    set<node*> nrset;
    double cx, cy, a, b, angle; // elipse parameters

    prset.clear();
    res->recurse_from_node(nn, tolayer0, nrset);
    node_set_to_region_set(prset, nrset, ly0sets, 0);

    set<ipoint2> prsetc = prset;

    centralize_set(prsetc);
    fit_ellipse(cx, cy, a, b, angle, prsetc);


    // ############
    //img im(prset, 0);
    //im.save("c:\\temp\\prset.bmp");
    //cout << "m = " << nnd->m << endl;
    //cout << "cx = " << cx << ", cy = " << cy << ", a = " << a << ", b = " << b << ", angle = " << angle << endl;
    //int dummy;
    //cin >> dummy;
    // ############

    if (a > 0.0) {
        double e = sqrt(a*a - b*b)/a;

        ell_stat_learn[nnd->m].ell_angle.new_data(2*angle);
        ell_stat_learn[nnd->m].ell_a.new_data(a);
        ell_stat_learn[nnd->m].ell_e.new_data(e);
    }

}

// ellipse_comparer
///////////////////////////////////////////////////////////////////////////////

class ellipse_comparer {
protected:
    //struct ell_parameters {
    //    double ell_angle, ell_a, ell_e;
    //};

    struct ell_parameters_distribution {
        von_mises_distribution ell_angle;
        normal_distribution1 ell_a;
        normal_distribution1 ell_e;

        ell_parameters_distribution() : ell_angle(), ell_a(), ell_e() { }
    };

    vector<ell_parameters_distribution> ell_statistics;

    double angle_thresh;
    double a_thresh;
    double e_thresh;
public:
    ellipse_comparer(const ellipse_learner& learner, double anglet, double at, double et);

    bool compare(int i, int j);

protected:
    void init(const ellipse_learner& learner);
};

ellipse_comparer::ellipse_comparer(const ellipse_learner& learner, double anglet, double at, double et) :
    ell_statistics(), angle_thresh(anglet), a_thresh(at), e_thresh(et) 
{
    init(learner);
}

void ellipse_comparer::init(const ellipse_learner& learner)
{
    ell_statistics.reserve(learner.ell_stat_learn.size());

    for (int i = 0; i < (int)learner.ell_stat_learn.size(); ++i) {
        const ellipse_learner::online_ell_parameters& item = learner.ell_stat_learn[i];

        von_mises_distribution angle_dist(item.ell_angle.get_parameters());
        normal_distribution1 a_dist(item.ell_a.get_mean(), item.ell_a.get_variance());
        normal_distribution1 e_dist(item.ell_e.get_mean(), item.ell_e.get_variance());

        ell_statistics.push_back(ell_parameters_distribution());
        ell_statistics.back().ell_angle = angle_dist;
        ell_statistics.back().ell_a = a_dist;
        ell_statistics.back().ell_e = e_dist;

        //cout << "i: " << i << endl;
        //cout << "    angle: (" << angle_dist.get_mu() << ',' << angle_dist.get_kappa() << ')' << endl;
        //cout << "        a: (" << a_dist.get_mean() << ',' << a_dist.get_variance() << ')' << endl;
        //cout << "        e: (" << e_dist.get_mean() << ',' << e_dist.get_variance() << ')' << endl;
    }
}

bool ellipse_comparer::compare(int i, int j)
{
    ell_parameters_distribution& dist0 = ell_statistics[i];
    ell_parameters_distribution& dist1 = ell_statistics[j];

    return cos(dist0.ell_angle.get_mu() - dist1.ell_angle.get_mu()) > angle_thresh &&
        abs(dist0.ell_a.get_mean() - dist1.ell_a.get_mean()) < a_thresh &&
        abs(dist0.ell_e.get_mean() - dist1.ell_e.get_mean()) < e_thresh;
}

class optimize_layer_stat : public streamable {
public:
	vector<int> statistics;
	optimize_layer_stat(int num_parts) : statistics(vector<int>(num_parts)) {}

	// streamable implementations
	virtual streamable* make_instance() const { return new optimize_layer_stat(0); }
	virtual void read_from_stream(istreamer& is) { is.read(statistics); }
	virtual void write_to_stream(ostreamer& os) { os.write(statistics);	}
};

class optimize_layer_mapreduce : public base_mapreduce {
	set<int> current_set; vector<set<ipoint2>> ly0sets;
	int num_parts; int layer; double cover_thresh; double int_thresh;
public:
	optimize_layer_mapreduce() { }
	optimize_layer_mapreduce(const set<int>& current_set, const vector<set<ipoint2>>& ly0sets, int num_parts, int layer, double cover_thresh, double int_thresh) :
		current_set(current_set), ly0sets(ly0sets), num_parts(num_parts), layer(layer), cover_thresh(cover_thresh), int_thresh(int_thresh) {}
	virtual ~optimize_layer_mapreduce() { }

	////////////////////////////////////////////////
	// main functionality is written in map and reduce
	// this is part of optimize_layer() that gathres statistics for each layer1_result
	virtual streamable* map(streamable* item) {
		// we know item should be streamed_pointer*
		layer1_result* res = (layer1_result*)((streamed_pointer*)item)->get();
		
		optimize_layer_stat* result = new optimize_layer_stat(num_parts);		
		
		int tolayer0 = atom("toLayer0");

		vector<node*>& s_nodes = res->shape_nodes[layer];

		// Calc covering ratio
		vector<node*> visited;
		set<node*> coveredn;
		set<node*> maxcoveredn;
		set<ipoint2> covered;
		set<ipoint2> maxcovered;

		select_nodes_by_type(visited, res, s_nodes, current_set);
		res->recurse(visited, tolayer0, coveredn);
		//get_positions(covered, coveredn, 0);
		node_set_to_region_set(covered, coveredn, ly0sets, 0);
		res->recurse(s_nodes, tolayer0, maxcoveredn);
		node_set_to_region_set(maxcovered, maxcoveredn, ly0sets, 0);

		double cover_ratio = (double)covered.size()/(double)maxcovered.size();

		//cout << "CR = " << cover_ratio << endl;

		if (cover_ratio >= cover_thresh) return result;

		// Update statistics
		for (vector<node*>::iterator sniter = s_nodes.begin(); sniter != s_nodes.end(); ++sniter) {
			node* n = *sniter;
			vector<node*> nvec;
			bool found = false;

			res->nodes_at(nvec, n);
			for (vector<node*>::iterator nviter = nvec.begin(); nviter != nvec.end(); ++nviter) {
				node* nn = *nviter;
				layer1_data* nnd = (layer1_data*)nn->data;
				if (current_set.find(nnd->m) == current_set.end()) {
					set<ipoint2> prset; 
					set<node*> nrset;

					res->recurse_from_node(nn, tolayer0, nrset);
					node_set_to_region_set(prset, nrset, ly0sets, 0);
					//elearner.new_data(prset, res, nn);
					int is = intersection_size(prset, covered);

					if (is < prset.size()*int_thresh)
						++result->statistics[nnd->m];
				}
			}

		}
		
		delete res;

		return result;
	}
	virtual streamable* reduce(list<streamable*> &item_list) {

		optimize_layer_stat* result = new optimize_layer_stat(num_parts);		

		for (auto iter = item_list.begin(); iter != item_list.end(); ++iter) {
			optimize_layer_stat* partial_result = (optimize_layer_stat*)*iter;
			for (int i = 0; i < num_parts; i++)
				result->statistics[i] += partial_result->statistics[i];

		}

		return result;
	}

	virtual string map_get_key(streamable* item) { return "merged_result"; }

	virtual bool has_reduce() { return true; }
	virtual bool should_wrap_results() { return false; }

	// streamable implementations
	virtual streamable* make_instance() const { return new optimize_layer_mapreduce(); }
	virtual void read_from_stream(istreamer& is) {
		is.read(current_set);
		is.read(ly0sets);
		is.read(num_parts);
		is.read(layer);
		is.read(cover_thresh);
		is.read(int_thresh);
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write(current_set);
		os.write(ly0sets);
		os.write(num_parts);
		os.write(layer);
		os.write(cover_thresh);
		os.write(int_thresh);
	}
};

/// Optimizes library on specific layer using workingset; algorithm parameters are 
/// * cover_thresh: threshold for covering ratio ( = |covered|/|maxcovered|)
///        if image is covered >= cover_threshold, skip it
/// * int_thresh: update part statistics if intersection of this part with already 
///        covered part of image >= int_thresh * size of part
/// * bite_size: how many parts are processed in one greedy step
/// * angle_thresh, a_thresh, e_thresh: thresholds for similarity of ellipses
///       (angle, major semi-axis, eccentricity)
/// Other
/// * init_set: initial parts (to be considered fixed)
/// Returns optimized set of parts (includes init_set).
set<int> optimize_layer(part_lib* library, list<streamed_pointer>& workingset, int layer, const set<int>& init_set,
    double cover_thresh, double int_thresh, int bite_size)
{
    typedef set<ipoint2> pset_t;

    if (library->layer_size(layer) == 0)
        return set<int>();

    nodes_at_sort_f sortf;
    set<int> current_set(init_set);
    vector<node*>& parts = library->parts[layer];
    vector<node*> bite;
    int tolayer0 = atom("toLayer0");
    vector<pset_t> ly0sets;

    library->get_regions(1, ly0sets);

	base_deployer_mapreudce* mapreduce_deployer = new openmp_deployer_mapreudce();

    do {        
		vector<int> statistics;
		int num_parts = library->layer_size(layer);		
        //ellipse_learner elearner(library, layer);

        //cout << endl;
        //cout << "Gathering statistics" << endl << "--------------------" << endl;
        cout << '.';		
		
		//////////////////////////////////////////////////////////////////////////////
		// CALL TO optimize_layer_mapreduce::map and optimize_layer_mapreduce::reduce
		// (check base_deployer_mapreudce for default calling implementation)
		{	
			mapreduce_items* mr_process_list = new mapreduce_items();

			for (auto iter = workingset.begin(); iter != workingset.end(); ++iter) {
				mr_process_list->items.push_back(&*iter); 
			}

			optimize_layer_mapreduce* collect_stat_mapreduce = new optimize_layer_mapreduce(current_set, ly0sets, num_parts, layer, cover_thresh, int_thresh);

			mapreduce_result* result = mapreduce_deployer->submit(mr_process_list, collect_stat_mapreduce);

			list<streamable*> stremable_results = result->get_result();

			if (stremable_results.size() == 1) {
				// results OK	
				optimize_layer_stat* collect_stat_result = (optimize_layer_stat*)stremable_results.front();
			
				statistics = collect_stat_result->statistics;
			
				delete collect_stat_result;
			} else {
				// results NOT OK				
				throw new_libhop_exception("optimize_layer() - mapreduce_deployer returned invalid results for optimize_layer_mapreduce (not all results from training set were returned)");
			}

			//for (auto iter = mr_process_list->items.begin(); iter != mr_process_list->items.end(); ++iter)
			//	delete *iter;

			delete mr_process_list;
			delete collect_stat_mapreduce;
		}
        //elearner.print();
        // #######################################
        //vector<img> images;
        //for (int index = 0; index < rset.size(); ++index) {
        //    images.push_back(img(rset[index], 0));
        //}
        //img imresult = img::concat(images, 0.0);
        //imresult.save("c:\\temp\\avg.png");
        // #######################################

        //ellipse_comparer ecomparer(elearner, angle_thresh, a_thresh, e_thresh);

        // Selection process
        vector<int> statordering = ordering<int>(statistics.begin(), statistics.end(), greater<int>());
        vector<node*> candidateparts;

        cout << '|' << statistics[statordering[0]];
        for (int i = 0; i < (int)statordering.size() && (int)candidateparts.size() < bite_size*8 &&   // !!!!8
                statistics[statordering[i]] > 0; ++i) {
            candidateparts.push_back(parts[statordering[i]]);
        }
        get_dissimilar_cluster(bite, candidateparts, bite_size);
        for (vector<node*>::iterator biter = bite.begin(); biter != bite.end(); ++biter) {
            lib_data* pd = (lib_data*)(*biter)->data;

            current_set.insert(pd->type);

        }

    } while (!bite.empty());

	delete mapreduce_deployer;

    // ###################
    cout << endl;
    for (set<int>::iterator iter = current_set.begin(); iter != current_set.end(); ++iter) 
        cout << *iter << ' ';
    cout << endl;
    // ###################

    return current_set;
}

/// Optimizes library on specific layer using workingset; algorithm parameters are 
/// * cover_thresh: threshold for covering ratio ( = |covered|/|maxcovered|)
///        if image is covered >= cover_threshold, skip it
/// * int_thresh: update part statistics if maximal intersection of this part with some 
///        already chosen part <= int_thresh*(size of this part)
/// * bite_size: how many parts are processed in one greedy step
/// * angle_thresh, a_thresh, e_thresh: thresholds for similarity of ellipses
///       (angle, major semi-axis, eccentricity)
/// Other
/// * init_set: initial parts (to be considered fixed)
/// Returns optimized set of parts (includes init_set).
set<int> optimize_layer_2(part_lib* library, list<streamed_pointer>& workingset, int layer, const set<int>& init_set,
    double cover_thresh, double int_thresh, int bite_size)
{
    cout << "optimize_layer_2" << endl;
    //cout << "Int_thresh = " << int_thresh << endl;

    typedef set<ipoint2> pset_t;

    if (library->layer_size(layer) == 0)
        return set<int>();

    nodes_at_sort_f sortf;
    set<int> current_set(init_set);
    vector<node*>& parts = library->parts[layer];
    vector<node*> bite;
    int tolayer0 = atom("toLayer0");
    int tonext = atom("toNextLayer");
    //int tonext = atom("toNextLayer0").get_index() + layer;
    vector<pset_t> ly0sets;

    library->get_regions(1, ly0sets);
    do {
        vector<double> statistics(library->layer_size(layer), 0.0);
        // ellipse_learner elearner(library, layer);

        //cout << endl;
        //cout << "Gathering statistics" << endl << "--------------------" << endl;
        cout << '.';
        for (auto iter = workingset.begin(); iter != workingset.end(); ++iter) {
			layer1_result* res = (layer1_result*)iter->get();
            vector<node*>& s_nodes = res->shape_nodes[layer];

            // Calc covering ratio
            vector<node*> visited;
            //set<node*> coveredn;
            //set<node*> maxcoveredn;
            //set<ipoint2> covered;
            //set<ipoint2> maxcovered;

            select_nodes_by_type(visited, res, s_nodes, current_set);
            //res->recurse(visited, tolayer0, coveredn);
            //get_positions(covered, coveredn, 0);
            //node_set_to_region_set(covered, coveredn, ly0sets, 0);
            //res->recurse(s_nodes, tolayer0, maxcoveredn);
            //node_set_to_region_set(maxcovered, maxcoveredn, ly0sets, 0);
            //{
                //set<ipoint2> pts;
                //get_positions(pts, maxcoveredn, 0);
                //cout << "#" << i << ": " << pts.size() << " quot: " << (double)pts.size()/res->shape_nodes[0].size();
                //cout << " cov: " << (double)covered.size()/(double)maxcovered.size() << endl;
            //}

            // Link visited nodes "backwards"
            for (vector<node*>::iterator viter = visited.begin(); viter != visited.end(); ++viter) {
                node* vn = *viter;
                set<node*> vnsupp;
                
                res->recurse_from_node(vn, tolayer0, vnsupp);
                for (set<node*>::iterator siter = vnsupp.begin(); siter != vnsupp.end(); ++siter) {
                    res->add_edge(*siter, vn, tonext);  
                }
            }

            //double cover_ratio = (double)covered.size()/(double)maxcovered.size();

            //cout << "CR = " << covered_ratio << endl;

            //if (cover_ratio >= cover_thresh) continue;

            // Update statistics
            for (vector<node*>::iterator sniter = s_nodes.begin(); sniter != s_nodes.end(); ++sniter) {
                node* n = *sniter;
                vector<node*> nvec;

                res->nodes_at(nvec, n);
                for (vector<node*>::iterator nviter = nvec.begin(); nviter != nvec.end(); ++nviter) {
                    node* nn = *nviter;
                    layer1_data* nnd = (layer1_data*)nn->data;

                    if (current_set.find(nnd->m) == current_set.end()) {
                        set<node*> nnsupp;
                        set<ipoint2> nnsuppr; 
                        set<node*> bnodes;
                        bool update = true;

                        res->recurse_from_node(nn, tolayer0, nnsupp);
                        node_set_to_region_set(nnsuppr, nnsupp, ly0sets, 0);
                        res->get_neighbors(bnodes, nnsupp.begin(), nnsupp.end(), tonext);

                        for (set<node*>::iterator siter = bnodes.begin(); siter != bnodes.end(); ++siter) {
                            node* bn = *siter;

                            if (bn == nullptr || bn == nn) 
                                continue;

                            set<node*> bnsupp;
                            set<ipoint2> bnsuppr;

                            res->recurse_from_node(bn, tolayer0, bnsupp);
                            node_set_to_region_set(bnsuppr, bnsupp, ly0sets, 0);
                            //cout << "/=" << (double)intersection_size(bnsuppr, nnsuppr)/nnsuppr.size() << " ";

                            if (intersection_size(bnsuppr, nnsuppr) >= nnsuppr.size()*int_thresh) {
                                update = false;
                                break;
                            }
                        }

                        if (update) {
                            //elearner.new_data(res, nn);
                            statistics[nnd->m] += nnd->vval();
                        }

                    }
                }

            }

			delete res;
        }
        
        //elearner.print();
        // #######################################
        //vector<img> images;
        //for (int index = 0; index < rset.size(); ++index) {
        //    images.push_back(img(rset[index], 0));
        //}
        //img imresult = img::concat(images, 0.0);
        //imresult.save("c:\\temp\\avg.png");
        // #######################################

        //ellipse_comparer ecomparer(elearner, angle_thresh, a_thresh, e_thresh);

        // Selection process
        vector<int> statordering = ordering<double>(statistics.begin(), statistics.end(), greater<double>());
        vector<node*> candidateparts;

        for (int i = 0; i < (int)statordering.size() && (int)candidateparts.size() < bite_size*8 &&   // !!!!8
                statistics[statordering[i]] > 0; ++i) {
            candidateparts.push_back(parts[statordering[i]]);
        }
        get_dissimilar_cluster(bite, candidateparts, bite_size);
        for (vector<node*>::iterator biter = bite.begin(); biter != bite.end(); ++biter) {
            lib_data* pd = (lib_data*)(*biter)->data;

            current_set.insert(pd->type);

        }

    } while (!bite.empty());

    // ###################
    cout << endl;
    for (set<int>::iterator iter = current_set.begin(); iter != current_set.end(); ++iter) 
        cout << *iter << ' ';
    cout << endl;
    // ###################

    return current_set;
}

/// Aghorithm for local optimization
/// * init_set: fixed parts
/// * n_parts: sample of parts, -1 = all
/// * gr_thresh: g-response threshold for update
/// * s_thresh: threshold for statistics
/// * bite_size: greedy step size 
set<int> optimize_layer_locally(part_lib* library, list<streamed_pointer>& workingset, int layer, const set<int>& init_set,
    int n_parts, double gr_thresh, double s_thresh, int bite_size)
{
    set<int> current_set(init_set);
    int lysize = library->layer_size(layer);
    vector<node*>& parts = library->parts[layer];
    vector<double> statistics;
    //vector<double>::iterator siter;
    double maxstatvalue = 0;
    double maxstat = 0;
    //int pass = 0; /// debug!!!!!

    do {
        //ellipse_learner elearner(library, layer);

        statistics.assign(lysize, 0.0);
        for (auto iter = workingset.begin(); iter != workingset.end(); ++iter) {
            layer1_result* res = (layer1_result*)iter->get();
            vector<node*>& s_nodes = res->shape_nodes[layer];
            int maxi = (n_parts < 0) ? (int)s_nodes.size() : min((int)s_nodes.size(), n_parts);

            res->init_grid(layer);
            for (int i = 0; i < maxi; ++i) {
                node* n = s_nodes[i];
                layer1_data* nd = (layer1_data*)n->data;
                vector<node*> nvec;
                double foundval = 0.0;

                for (int x = nd->x - 1; x < nd->x + 2; ++x) 
                    for (int y = nd->y - 1; y < nd->y + 2; ++y) 
                        res->nodes_at(nvec, x, y, layer);
                for (vector<node*>::iterator nviter = nvec.begin(); nviter != nvec.end(); ++nviter) {
                    node* nn = *nviter;
                    layer1_data* nnd = (layer1_data*)nn->data;
                    if (current_set.find(nnd->m) != current_set.end()) {
                        if (nnd->val() > foundval) foundval = nnd->val();
                    }
                }
                if (foundval < gr_thresh) {
                    if (nd->val() > gr_thresh) {
                        statistics[nd->m] += nd->val();
                    }
                }
                
            }

        }

        cout << endl;
        // ########
        //string fname = string("c:\\work\\stat_") + (++pass) + string(".m");
        //ofstream os(fname.c_str());
        //os << '{';
        //for (int i = 0; i < statistics.size(); ++i) {
        //    if (i > 0) os << ',';
        //    os << statistics[i];
        //}
        //os << '}';
        //os.close();
        // ########

        vector<node*> bite;
        vector<int> statordering = ordering<double>(statistics.begin(), statistics.end(), greater<double>());
        vector<node*> candidateparts;

        maxstatvalue = statistics[statordering[0]];
        if (maxstatvalue > maxstat) maxstat = maxstatvalue;
        for (int i = 0; i < (int)statordering.size() && (int)candidateparts.size() < bite_size*8 &&  // !!!!8
                statistics[statordering[i]] > maxstat * s_thresh; ++i) { 
            candidateparts.push_back(parts[statordering[i]]);
        }
        get_dissimilar_cluster(bite, candidateparts, bite_size);
        for (vector<node*>::iterator biter = bite.begin(); biter != bite.end(); ++biter) {
            lib_data* pd = (lib_data*)(*biter)->data;

            current_set.insert(pd->type);

        }

        //cout << endl << " max stat. element: " << *siter << " (thresh: " << maxstat * s_thresh << ')' << endl;
    } while (maxstatvalue > maxstat * s_thresh);

    // ###################
    cout << endl;
    for (set<int>::iterator iter = current_set.begin(); iter != current_set.end(); ++iter) 
        cout << *iter << ' ';
    cout << endl;
    // ###################
    
    return current_set;
}

// Adds "lyrSimVal" edges between parts of the same layer.
// Geometry (path map size) is reduced to 'max_geo_size'.
// Edges are added only if similarity is > thresh.
void add_similarity_edges(part_lib* library, int layer, int max_geo_size, double thresh)
{
    int ename = atom("lyrSimVal");

    if (layer < 0 || layer > library->max_layer_index())
        return;

    for (int pi = 0; pi < library->layer_size(layer); ++pi) {
        node* p = library->parts[layer][pi];
        vector<hpoint_t> phpts;
        vector<ipoint2> ppts;
        path_map_t pm;

        get_library_geo(pm, p);
        inhibit_path_map(pm, max_geo_size);

        for (auto pmiter = pm.begin(); pmiter != pm.end(); ++pmiter) 
            phpts.push_back(pmiter->second);
        ppts = extract<ipoint2>(phpts.begin(), phpts.end(), [](const hpoint_t& hp) { return hp.p; });
        for (int qi = 0; qi < library->layer_size(layer); ++qi) {
            node* q = library->parts[layer][qi];
            vector<hpoint_t> qhpts;
            vector<ipoint2> qpts;
            path_map_t qm;

            get_library_geo(qm, q);
            for (auto qmiter = qm.begin(); qmiter != qm.end(); ++qmiter) 
                qhpts.push_back(qmiter->second);
            qpts = extract<ipoint2>(qhpts.begin(), qhpts.end(), [](const hpoint_t& hp) { return hp.p; });

            vector<int> perm = point_matching(ppts, qpts);
            double scdistance = 0.0;
            for (int i = 0; i < (int)perm.size(); ++i) {
                scdistance += sc_histogram_distance(phpts[perm[i]].h, qhpts[i].h);
            }
            scdistance /= perm.size();
            library->add_edge_2_unique(p, q, new edge_data_t<double>(scdistance), ename);
        }
    }
}



