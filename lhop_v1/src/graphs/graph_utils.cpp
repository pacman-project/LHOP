
#include "graph_utils.h"
#include "utils/utils.h"
#include <algorithm>


// Returns a minimum spanning tree graph based on Euclidean distances between points 
// in 'pts' vector.
// Note: 
// - node data of node 'n' is of type 'node_data_t<int>' and contains index of a point
//   corresponding to 'n'
graph* point_graph(const vector<ipoint2>& pts)
{
    if (pts.empty()) return new graph();

    graph* g = new graph();
    vector<node*> nodes(pts.size(), nullptr);
    vector<int> dist(pts.size(), INT_MAX/2);    // dist[i] is distance from i to tree
    vector<int> closest(pts.size(), 0);         // closest[i] is the closest point to i in current tree

    while (true) {
        int mind = INT_MAX;
        int mini = -1;

        for (int i = 0; i < (int)dist.size(); ++i) {
            if (dist[i] < mind) {
                mind = dist[i];
                mini = i;
            }
        }
        if (mind == INT_MAX) 
            break;
        nodes[mini] = g->add_node(new node_data_t<int>(mini));
        if (mini != closest[mini]) 
            g->add_edge(nodes[mini], nodes[closest[mini]], 0, 0);
        for (int i = 0; i < (int)pts.size(); ++i) {
            if (dist[i] < INT_MAX) {
                int d = pts[mini].distance2(pts[i]);

                if (d < dist[i]) { dist[i] = d; closest[i] = mini; }
            }
        }
        dist[mini] = INT_MAX;
    }   
    return g;
}

img point_graph_image(graph* g, const vector<ipoint2>& pts, int border /* = 10 */)
{
    irectangle2 box = irectangle2::bounding_rectangle(pts.begin(), pts.end());

    if (box.invalid()) return img(2*border, 2*border, 0.0);

    img result(2*border + box.x_dim() + 1, 2*border + box.y_dim() + 1, 0.0);
    int xdelta = border - box.ll.x;
    int ydelta = border - box.ll.y;

    for (graph::iter_t iter = g->nodes.begin(); iter != g->nodes.end(); ++iter) {
        node* n = *iter;
        node_data_t<int>* nd = (node_data_t<int>*)n->data;

        forall_neighbors(n, niter) {
            node_data_t<int>* nnd = (node_data_t<int>*)neighbor_node_data(niter);

            result.draw_line(pts[nd->data].x + xdelta, pts[nd->data].y + ydelta,
                pts[nnd->data].x + xdelta, pts[nnd->data].y + ydelta, 1.0);
        }
    }
    return result;
}

// Longest path from 'start' node in a directed acyclic graph.
vector<node*> dag_longest_path(graph* g, node* start)
{
    vector<node*> result;
    vector<node*> path;

    g->clear_attr(ATTR_VISITED);
    path.reserve(g->nodes.size());
    path.push_back(start);
    while (!path.empty()) {
        node* n = path.back();
        node* un = nullptr;

        n->set_attr(ATTR_VISITED);

        // Find unvisited neighbor
        forall_neighbors (n, iter) {
            node* nn = neighbor_node(iter);

            if (!nn->is_attr_set(ATTR_VISITED)) {
                un = nn;
                break;
            }
        }
        if (un != nullptr) path.push_back(un);
        else {
            if (result.size() < path.size()) result = path;
            path.pop_back();
        } 
    }
    g->clear_attr(ATTR_VISITED);
    return result;
}

// Use to find longest path if 'g' is tree.
vector<node*> tree_longest_path(graph* g)
{
    vector<node*> result;

    if (!g->nodes.empty()) {
        result = dag_longest_path(g, g->nodes.front());
        result = dag_longest_path(g, result.back());
    }
    return result;
}

// Use to find longest path if 'g' is a forest.
vector<node*> forest_longest_path(graph* g)
{
    vector<node*> result;
    set<node*> availnodes(g->nodes.begin(), g->nodes.end());

    while (!availnodes.empty()) {
        vector<node*> tmpresult;
        
        tmpresult = dag_longest_path(g, *availnodes.begin());
        tmpresult = dag_longest_path(g, tmpresult.back());

        availnodes.erase(availnodes.begin());
        for (auto iter = tmpresult.begin(); iter != tmpresult.end(); ++iter) 
            availnodes.erase(*iter);
        if (tmpresult.size() > result.size()) 
            swap(tmpresult, result);
    }
    return result;
}


// Closest path in dfs tree from 'start' to 'nset'.
// Returns empty path if start is in nset.
vector<node*> dfs_closest_path(graph* g, const set<node*>& nset, node* start)
{
    if (nset.find(start) != nset.end())
        return vector<node*>();

    vector<node*> result;
    vector<node*> path;

    g->clear_attr(ATTR_VISITED);
    path.reserve(g->nodes.size());
    path.push_back(start);
    while (!path.empty()) {
        node* n = path.back();
        node* un = nullptr;

        n->set_attr(ATTR_VISITED);

        // Find unvisited neighbor
        forall_neighbors (n, iter) {
            node* nn = neighbor_node(iter);

            if (nset.find(nn) != nset.end()) {
                if (result.empty() || path.size() + 1 < result.size()) {
                    result = path;
                    result.push_back(nn);
                }
            } else if (!nn->is_attr_set(ATTR_VISITED)) {
                un = nn;
                break;
            }
        }
        if (un != nullptr) path.push_back(un);
        else path.pop_back();
    }
    g->clear_attr(ATTR_VISITED);
    return result;
}

// Auxillary function for 'point_matching_p'
void match_paths(vector<int>& result, double& val, 
    graph* g1, graph* g2, const vector<node*> lpath1, const vector<node*> lpath2,
    const vector<ipoint2>& pts1, const vector<ipoint2>& pts2)
{
    double f = (double)(lpath2.size() - 1)/(lpath1.size() - 1);
    set<node*> matched; // Matched nodes of g2.

    val = 0.0; // Sum of distances^2
    result.resize(pts2.size());
    
    // Match lpath1 to lpath2
    for (int i = 0; i < (int)lpath1.size(); ++i) {
        int j = int_round(f*i);
        node* n1 = lpath1[i];
        node* n2 = lpath2[j];
        int i1 = ((node_data_t<int>*)n1->data)->data;
        int i2 = ((node_data_t<int>*)n2->data)->data;

        val += pts1[i1].distance2(pts2[i2]);
        result[i2] = i1;
        matched.insert(n2);
    }

    // Match rest
    set<node*> unmatched;

    set_difference(unmatched, set<node*>(g2->nodes.begin(), g2->nodes.end()), matched);
    while (!unmatched.empty()) {
        node* n = *unmatched.begin();
        vector<node*> p = dfs_closest_path(g2, matched, n);

        for (int pi = 0; pi + 1 < (int)p.size(); ++pi) {
            node* pn = p[pi];
            int i2 = ((node_data_t<int>*)pn->data)->data;
            int i1 = ((node_data_t<int>*)p.back()->data)->data;
            
            result[i2] = result[i1];
            val += pts1[result[i1]].distance2(pts2[i2]);
            unmatched.erase(pn);
            matched.insert(pn);
        }
    }
}

void save_matching_dbg(const string& fname, const vector<int>& perm, const vector<ipoint2>& pts1, 
    const vector<ipoint2>& pts2)
{
    ofstream os(fname.c_str());
    int maxi = max((int)pts1.size(), (int)pts2.size());

    for (int i = 0; i < maxi; ++i) {
        if (i < perm.size()) os << (perm[i] + 1) << ',';
        else os << 0 << ',';
        if (i < pts1.size()) os << pts1[i].x << ',' << pts1[i].y << ',';
        else os << 0 << ',' << 0 << ',';
        if (i < pts2.size()) os << pts2[i].x << ',' << pts2[i].y;
        else os << 0 << ',' << 0;            
        os << endl;
    }
    os.close();
}

// Match points 'pts1' to points 'pts2' using matching of longest paths in the 
// minimum spanning trees; retruns mapping 'm' such that 
// pts1[m[i]] <-> pts2[i]. Note that 'm' is not necessarily a permutation!
vector<int> point_matching_p(const vector<ipoint2>& pts1, const vector<ipoint2>& pts2)
{
    graph* g1 = point_graph(pts1);
    vector<node*> lpath1 = tree_longest_path(g1);
    graph* g2 = point_graph(pts2);
    vector<node*> lpath2 = tree_longest_path(g2);

    vector<int> result;
    double val;
    vector<int> resultinv;
    double valinv;

    match_paths(result, val, g1, g2, lpath1, lpath2, pts1, pts2);
    //save_matching_dbg("c:\\work\\match.csv", result, cast_vector<dpoint2, ipoint2>(pts1), 
    //    cast_vector<dpoint2, ipoint2>(pts2));
    reverse(lpath1.begin(), lpath1.end());
    match_paths(resultinv, valinv, g1, g2, lpath1, lpath2, pts1, pts2);
    //save_matching_dbg("c:\\work\\matchinv.csv", resultinv, cast_vector<dpoint2, ipoint2>(pts1), 
    //    cast_vector<dpoint2, ipoint2>(pts2));

    delete g2;
    delete g1;

    if (valinv < val) return resultinv;
    else return result;
}

vector<int> point_matching_p(const vector<dpoint2>& pts1, const vector<dpoint2>& pts2, double factor /* = 100 */)
{
    return point_matching_p(cast_vector<ipoint2, dpoint2, double>(pts1, factor),
        cast_vector<ipoint2, dpoint2, double>(pts2, factor));
}

// Match points 'pts1' to points 'pts2' using graph matching "point_matching_p" for each of the pieces;
// "piece" is a subset of points with the same .first, piece numbers must start
vector<int> piecewise_point_matching_p(const vector<pair<int, ipoint2> >& pts1, const vector<pair<int, ipoint2> >& pts2)
{
    typedef pair<int, ipoint2> item_t;
    typedef vector<item_t> vector_t;
    typedef vector<int> result_t;
    typedef map<int, vector<int> > map_t;

    if (pts2.empty()) return result_t();

    result_t result(pts2.size(), -1);
    
    map_t pieces1, pieces2;

    for (int i = 0; i < (int)pts1.size(); ++i) 
        pieces1[pts1[i].first].push_back(i);
    for (int i = 0; i < (int)pts2.size(); ++i) 
        pieces2[pts2[i].first].push_back(i);
    for (map_t::iterator piter = pieces2.begin(); piter != pieces2.end(); ++piter) {
        vector<int>& piece2 = piter->second;
        vector<int>& piece1 = pieces1[piter->first];

        if (!piece1.empty()) {
            vector<ipoint2> p1;
            vector<ipoint2> p2;

            p1.reserve(piece1.size());
            p2.reserve(piece2.size());
            for (int i = 0; i < (int)piece1.size(); ++i) 
                p1.push_back(pts1[piece1[i]].second);
            for (int i = 0; i < (int)piece2.size(); ++i) 
                p2.push_back(pts2[piece2[i]].second);

            vector<int> match = point_matching_p(p1, p2);

            for (int i = 0; i < (int)match.size(); ++i)
                result[piece2[i]] = piece1[match[i]];
        }

    }
    return result;
}

vector<double> curve_gaps(graph* g, const vector<ipoint2>& pts)
{
    set<node*> visited;
    vector<double> result;

    for (graph::iter_t niter = g->nodes.begin(); niter != g->nodes.end(); ++niter) {
        node* n = *niter;
        node_data_t<int>* nd = (node_data_t<int>*)n->data;

        forall_neighbors (n, iter) {
            node* nn = neighbor_node(iter);
            node_data_t<int>* nnd = (node_data_t<int>*)nn->data;
            
            if (visited.find(nn) == visited.end()) 
                result.push_back(sqrt((double)(pts[nnd->data].distance2(pts[nd->data]))));
        }
        visited.insert(n);
    }
    return result;
}

bool continuous_point_set(const vector<ipoint2>& ipts, double factor /* = 3.0 */)
{
    graph* g = point_graph(ipts);
    vector<double> gaps = curve_gaps(g, ipts);
    double m = median_value(gaps, 0.8);

    delete g;
    return (*max_element(gaps.begin(), gaps.end())) < factor*m;
}

double curve_length(const vector<ipoint2>& ipts)
{
    graph* g = point_graph(ipts);
    vector<double> gaps = curve_gaps(g, ipts);
    double result = 0.0;

    for (vector<double>::iterator iter = gaps.begin(); iter != gaps.end(); ++iter) {
        result += *iter;
    }

    delete g;
    return result;
}


