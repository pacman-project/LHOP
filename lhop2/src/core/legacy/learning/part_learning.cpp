
#include <opencv2/opencv.hpp>
#include <queue>

#include "part_learning.h"

#include "utils/graphs/graph_utils.h"

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


struct position_compare_less {
    bool operator()(const ipoint2& p1, const ipoint2& p2)
    {
        if (abs(p1.x - p2.x) < 2*abs(p1.y - p2.y)) return p1.y < p2.y; else return p1.x < p2.x;
    }
};

// part_learning 
///////////////////////////////////////////////////////////////////////////////


part_learning::part_learning() :
    stat(), maxima(), layer_contraction(-1.0)
{ 
}

part_learning::part_learning(const ConfigDictionary& cfg) :
    stat(), maxima(), layer_contraction(-1.0)
{ 
    init_cfg(cfg);
}

part_learning::~part_learning()
{
}

void part_learning::init_cfg(const ConfigDictionary& cfg)
{
    if (cfg.isDefined("source_layer_index")) {
        source_layer = cfg.getValueInt("source_layer_index", 0) - 1;  // backward compatibility!
    } else {
        cfg.getValue(source_layer, "source_layer", true);
    }

    // part update
    center_val_threshold = cfg.getValueDouble("center_val_threshold", -1.0);
    center_val_threshold_rel = cfg.getValueDouble("center_val_threshold_rel", 0.2);
    if (cfg.isDefined("seq_max_intersection_percent2"))
        seq_min_intersection_percent = cfg.getValueDouble("seq_min_intersection_percent2", 0.0);
    else
        seq_min_intersection_percent = cfg.getValueDouble("seq_min_intersection_percent", 0.0);
    if (cfg.isDefined("seq_max_intersection_percent2"))
        seq_max_intersection_percent = cfg.getValueDouble("seq_max_intersection_percent2", 0.5);
    else
        seq_max_intersection_percent = cfg.getValueDouble("seq_max_intersection_percent", 0.5);
    max_candidates = cfg.getValueInt("max_candidates", 5);
    min_seq_size = cfg.getValueInt("min_part_length", 2);
    max_seq_size = cfg.getValueInt("max_part_length", 3);
    max_stat_size = cfg.getValueInt("max_stat_size", INT_MAX);

    // similarity calculation
    similarity_type = cfg.getValueInt("similarity_type", -1);

    // similarity calculation inhibition 
    max_matching_size = cfg.getValueInt("max_matching_size", 30);
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
            // just take the best one (maybe 90% of the best one)
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

	// max_seq_size == max_part_length 
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

// Support function for make_path_map; l contains a list of incomplete paths and their end nodes.
void make_path_map_rec(list<pair<vector<int>, node*> >& l)
{
	typedef pair<vector<int>, node*> list_item_t;
	typedef list<list_item_t> list_t;

	int name = EdgeConnection::TO_PREV_LAYER;

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
// Extends map 'm' with rules path |-> hpoint_t = {p, h} where p is a coordinate
// of the endpoint of the path, h is a histogram at p
// 'scmap' is a shape context map
// 'cn' is center node
// 'n' is starting node
void make_path_map(path_map_t& m, const scmap_t& scmap, node* cn, node* n)
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

// Find last node in "toPrevLayer path" with index 0 (indicating central node).
// (i.e. finds central subpart on layer 1 of node cn)
node* get_center_node(node* cn)
{
    int name = EdgeConnection::TO_PREV_LAYER;
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
    sym_part_t sorted_seq = seq;
    real_part_t sorted_rseq = rseq;

    permute_range<ipoint2>(sorted_seq.begin() + 1, ord);
    permute_range<node*>(sorted_rseq.nodes.begin() + 1, ord);

    stat_t::mapped_type& s = stat[sorted_seq];

	// this is where part (with its subparts) is first time saved into statistics
    s.count += sorted_rseq.val; // used later to sort each new part

	// following lines are used only to generate geometry (shape context) and PCA information
	node* cn0 = get_center_node(cn);
    vector<ipoint2> pts;

    if (cn0 == nullptr || node_layer(cn0) != 0) 
        throw new_libhop_exception("Path to central node on layer 0 does not exist.");
    
    for (int i = 0; i < (int)sorted_rseq.nodes.size(); ++i) {
        node* n = sorted_rseq.nodes[i];
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
void update_geometry(node* p, const part_learning::geometry_t& g)
{
    int name = EdgeConnection::TO_LYR_SOURCE;
    edge_pair cpair = p->get_neighbor_pair(EdgeConnection::TO_LYR_CENTER);
    
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
void update_pca(node* p, const pca_data& pcd)
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
    int rootname = EdgeConnection::TO_LYR_SIMROOT;
    int simname = EdgeConnection::TO_LYR_SIMILAR;
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

                    library->add_edge_2(p, p0, new part_data_sim(mi[j].val), nullptr, rootname, simname);
                }
            }
        }
    }
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

// Returns distance between two geometies.
// We assume that geo[0] is geometry of central part; 
// igeo[1:] and jgeo[1:] and returns the one with the smallest distance:
// 'benergy' is TPS bending energy
// 'scmatching' is average shape context matching
// return value is permutation on 'igeo' used to match subparts of 'igeo' to 'jgeo'
vector<int> geometry_distance(double& benergy, vector<dpoint2>& dvector, double& scmatching, 
    const part_learning::geometry_t& igeo, const part_learning::geometry_t& jgeo, bool calc_benergy, int max_matching_size)
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

                    perm = geometry_distance(benergy, dvector, scdist, iter->second, jter->second, ethresh < 1000, max_matching_size); 
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
}


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

	struct clique_data {
		set<int> parts;
		int part0;
		cv::Mat mean;
		cv::Mat eigenvectors;
		cv::Mat eigenvalues;
		double sizefactor;
	};

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

void part_learning::get_pca_stat(map<sym_part_t, pca_data>& stat3)
{
    typedef map<sym_part_t, pca_data> stat3_t;

    stat3.clear();
    for (stat_t::iterator siter = stat.begin(); siter != stat.end(); ++siter) {
        int gsize = (int)siter->first.size();
        pca_data& geo3 = stat3[siter->first];

        geo3 = siter->second.pcal.get_pca_data();
    }
}

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

