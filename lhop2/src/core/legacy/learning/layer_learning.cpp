
#include <opencv2/opencv.hpp>
#include <queue>
#include "layer_learning.h"
#include "utils/graphs/graph_utils.h"

// global functions
///////////////////////////////////////////////////////////////////////////////


// make set of integers 
//  (x, y) - center
//  width, height - dimensions of the picture
//  r - region
//  (cx, cy) - center of the region
//  result - note that it is NOT deleted, elements are added to this set!!
void get_region_set(int x, int y, int width, int height, const matrix<bool>& r, int cx, int cy, 
                    set<int>& result)
{
    int rwidth = (int)r.width, rheight = (int)r.height;
    int x_start, y_start, x_end, y_end;
    int ix_start, iy_start;

    x_start = max(x - cx, 0); 
    x_end = min(x + rwidth - cx, width);
    y_start = max(y - cy, 0); 
    y_end = min(y + rheight - cy, height);
    ix_start = x_start - (x - cx);
    iy_start = y_start - (y - cy);

    //result.clear();
    for (int i = ix_start, x = x_start; x < x_end; ++i, ++x) {
        for (int j = iy_start, y = y_start; y < y_end; ++j, ++y) {
            if (r(i, j)) result.insert(x + y*width);
        }
    }
}


/**
 * For each node in res (at specified layer) get their support points, i.e.
 * get theirs layer 1 locations (as one dimensional location = y * width + x)
 */
void get_region_map(map<node*, set<int> >& rmap, layer1_result* res, int layer,
    const vector<matrix<bool> >& regions, const vector<ipoint2>& rcenters)
{
    typedef map<node*, set<int> > result_t;

    int to_prev = EdgeConnection::TO_PREV_LAYER;

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

void max_normalize(matrix<double>& m, double factor)
{
    double max = m.maximum();

    if (max == 0.0) return;
    for_each_element (m, i) {
        m[i] /= max*factor;
    }
}

void get_reconstruction_map(map<node*, set<node*> >& result, layer1_result* res, int layer)
{
    typedef map<node*, set<node*> > result_t;

    int prevname = EdgeConnection::TO_PREV_LAYER;

    result.clear();
    for (list<node*>::iterator iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
        node* n = *iter;
        result_t::iterator riter = result.insert(result_t::value_type(n, result_t::mapped_type())).first;

        res->recurse_from_node(n, prevname, riter->second);
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
    int prevname = EdgeConnection::TO_PREV_LAYER;
    int srcname = EdgeConnection::TO_LYR_SOURCE;
    int centername = EdgeConnection::TO_LYR_CENTER;
    double factor = (double)res->x_size(layer)/res->x_size(layer - 1);

    if (factor < 0.9 || factor > 1.1) {
        cout << "Warning: Contraction factor is not 1.0. No maps will be updated." << endl;
        return;
    }

    get_reconstruction_map(rmap, res, layer - 1);
    
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
                        continue; 
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

    int srcname = EdgeConnection::TO_LYR_SOURCE;
    int centername = EdgeConnection::TO_LYR_CENTER;

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

template<class T, class U> pair<T, U> sort(const pair<T, U>& p)
{
    if (p.first <= p.second) return p; else return pair<T, U>(p.second, p.first);
}
// functions
///////////////////////////////////////////////////////////////////////////////

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

