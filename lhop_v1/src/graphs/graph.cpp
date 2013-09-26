// graph library includes definitions found in
//  - graph.h
//  - img_graph.h
////////////////////////////////////////////////

#define _USE_MATH_DEFINES  // necessary for <math.h> to define M_PI,.... constants, etc. 
#include <math.h>
#include <queue>
#include "../utils/atom.h"
#include "../utils/streaming.h"
#include "../utils/structures.h"
#include "../utils/utils.h"
#include "graph.h"


#include "../layers/part_lib.h"

// node
///////////////////////////////////////////////////////////////////////////////

int node::count_neighbors(int n)
{
    int count = 0;

    for (iterator iter = get_neighbors_begin(n); iter != neighbors.end() && iter->first == n; ++iter)
        ++count;
    return count;
}

node* node::get_neighbor(int n)
{
    iterator iter = get_neighbors_begin(n);
    if (iter == neighbors.end()) return nullptr;
    else return iter->second.first;
}

edge_data* node::get_edge_data(int n)
{
    iterator iter = get_neighbors_begin(n);
    if (iter == neighbors.end()) return nullptr;
    else return iter->second.second;
}

edge_data* node::get_edge_data(node* n, int name)
{
    foreach_neighbor(this, name, iter) {
        if (neighbor_node(iter) == n) return neighbor_edge_data(iter);
    }
    return nullptr;
}


edge_pair node::get_neighbor_pair(int n)
{
    iterator iter = get_neighbors_begin(n);
    if (iter == neighbors.end()) return edge_pair(nullptr, nullptr);
    else return iter->second;
}

void node::get_neighbor_set(int n, set<node*>& nset)
{
    nset.clear();
    iterator iter;
    for (iter = get_neighbors_begin(n); iter != neighbors.end() && iter->first == n; ++iter)
        nset.insert(iter->second.first);
}

void node::get_neighbors(set<node*>& n)
{
    n.clear();
    forall_neighbors (this, iter) {
        n.insert(neighbor_node(iter));
    }
}

void node::add_to_neighbor_set(int n, set<node*>& nset)
{
    iterator iter;
    for (iter = get_neighbors_begin(n); iter != neighbors.end() && iter->first == n; ++iter)
        nset.insert(iter->second.first);
}

void node::get_neighbor_pair_set(int n, set<edge_pair>& nset)
{
    nset.clear();
    iterator iter;
    for (iter = get_neighbors_begin(n); iter != neighbors.end() && iter->first == n; ++iter)
        nset.insert(iter->second);
}

bool node::is_neighbor(node* n, int name)
{
    foreach_neighbor(this, name, iter) {
        if (neighbor_node(iter) == n) return true;
    }
    return false;
}

int node::degree(int n)
{
    int result = 0;

    for (iterator iter = get_neighbors_begin(n); iter != neighbors.end() && iter->first == n; ++iter)
        ++result;
    return result;
}

void node::get_neighbors_iter(int n, node::iterator& begin, node::iterator& end)
{
    end = begin = neighbors.find(n);
    while (end != neighbors.end() && end->first == n) ++end;
}

void node::delete_edges(const set<node*>& ns)
{
    multimap<int, edge_pair> newneighbors;

    for (iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        if (ns.find(neighbor_node(iter)) == ns.end()) newneighbors.insert(*iter);
        else if (neighbor_edge_data(iter) != nullptr) delete neighbor_edge_data(iter);
    } 
    neighbors = newneighbors;
}

void node::delete_edges(int index, node* n)
{
    multimap<int, edge_pair> newneighbors;

    for (iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        if (neighbor_node(iter) != n || neighbor_index(iter) != index) newneighbors.insert(*iter);
        else if (neighbor_edge_data(iter) != nullptr) delete neighbor_edge_data(iter);
    } 
    neighbors = newneighbors;
}

void node::gc_clone_to(node* n, graph* g, const map<node*, node*>& cmap) const
{
    map<node*,node*>::const_iterator f;

    n->attr = attr;
    for (const_iterator i = neighbors.begin(); i != neighbors.end(); ++i) {
        node* nb = i->second.first;
        edge_data* ed = i->second.second;

        if ((f = cmap.find(nb)) != cmap.end()) {
            if (ed == nullptr) g->add_edge(n, f->second, i->first);
            else g->add_edge_2(n, f->second, ed->get_gc_clone(cmap), i->first);
        }
    }
    if (data) n->data = data->get_gc_clone(cmap);
}

void node::delete_edges(node* n)
{
    multimap<int, edge_pair> newneighbors;

    for (iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
       if (neighbor_node(iter) != n) newneighbors.insert(*iter);
       else if (neighbor_edge_data(iter) != nullptr) delete neighbor_edge_data(iter);
    }
    neighbors = newneighbors;
}

void node::delete_edges(unsigned attr)
{
    multimap<int, edge_pair> newneighbors;

    for (iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        if (!neighbor_node(iter)->is_attr_set(attr)) newneighbors.insert(*iter);
        else if (neighbor_edge_data(iter) != nullptr) delete neighbor_edge_data(iter);
    }
    neighbors = newneighbors;
}

void node::delete_edges(int index)
{
    multimap<int, edge_pair> newneighbors;

    for (iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        if (neighbor_index(iter) != index) newneighbors.insert(*iter);
        else if (neighbor_edge_data(iter) != nullptr) delete neighbor_edge_data(iter);
    }
    neighbors = newneighbors;
}

void node::delete_edges_complement(int index)
{
    multimap<int, edge_pair> newneighbors;

    for (iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        if (neighbor_index(iter) == index) newneighbors.insert(*iter);
        else if (neighbor_edge_data(iter) != nullptr) delete neighbor_edge_data(iter);
    }
    neighbors = newneighbors;
}

void node::rename_edges(int from, int to)
{
    multimap<int, edge_pair> newneighbors;

    for (iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        if (neighbor_index(iter) != from) newneighbors.insert(*iter);
        else newneighbors.insert(container::value_type(to, iter->second));
    }
    neighbors.swap(newneighbors);
}

void node::copy_to(streamable* p, cloner& cl)
{
    int key;
    streamable* n, * ed;

    for (iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        key = iter->first;
        n = cl.get_copy(iter->second.first);
        ed = cl.get_copy(iter->second.second);
        ((node*)p)->neighbors.insert(node_pair(key, edge_pair((node*)n, (edge_data*)ed)));
    }
    ((node*)p)->data = (node_data*)cl.get_copy(data);
    ((node*)p)->attr = attr;
}

void node::read_from_stream(istreamer& is)
{
    streamable::read_from_stream(is);

    unsigned nnum;
    int key;
    streamable* n,* ed;

    is.read(nnum);
    for (unsigned i = 0; i < nnum; ++i) {
        is.read(key);
        is.read(n);
        is.read(ed);
        neighbors.insert(node_pair(key, edge_pair((node*)n, (edge_data*)ed)));
    }
    is.read((streamable*&)data);
    is.read(attr);
}

void node::write_to_stream(ostreamer& os)
{
    streamable::write_to_stream(os);
    
    unsigned nnum = (unsigned)neighbors.size();

    os.write(nnum);
    for (iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
        os.write((int)iter->first);
        os.write(iter->second.first);
        os.write(iter->second.second);
    }
    os.write(data);
    os.write(attr);
}


// graph
///////////////////////////////////////////////////////////////////////////////

graph::~graph()
{
    if (!is_subgraph) delete_nodes();
}

void graph::get_edges(list<edge_data*>& edges, node* n1, node* n2, int name)
{
    edges.clear();
    foreach_neighbor(n1, name, iter) {
        if (iter->second.first == n2) edges.push_back(iter->second.second);
    }
}

edge_data* graph::get_edge_data(node* n1, node* n2, int name)
{
    foreach_neighbor(n1, name, iter) {
        if (neighbor_node(iter) == n2) return neighbor_edge_data(iter);
    }
    return nullptr;
}

edge_data** graph::get_edge_data_p(node* n1, node* n2, int name) 
{
    foreach_neighbor(n1, name, iter) {
        if (neighbor_node(iter) == n2) return &neighbor_edge_data(iter);
    }
    return nullptr;
}

void graph::add_reverse_edges(node* n, int name, int rname)
{
    foreach_neighbor(n, name, iter) {
        add_edge(neighbor_node(iter), n, rname);
    }
}

void graph::add_reverse_edges(int name, int rname)
{
    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        add_reverse_edges(*iter, name, rname);
    }
}

void graph::get_node_enumeration(const vector<node*>& nodes, map<node*, int>& enumeration)
{
    vector<node*>::const_iterator iter;
    int count = 0;

    for (iter = nodes.begin(); iter != nodes.end(); ++iter) {
        enumeration.insert(pair<node*,int>(*iter, count++));
    }
}

void graph::get_all_edges(set<int>& result)
{
    result.clear();
    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;

        forall_neighbors(n, iter) {
            result.insert(neighbor_index(iter));
        }
    }
}

void graph::write_vgr(ostream& os, const vector<node*>& nodes, int edgename)
{
    vector<node*>::const_iterator iter;
    double alpha = 2*M_PI/(double)nodes.size();
    map<node*, int> enumeration;
    int count;

    count = 0;
    for (iter = nodes.begin(); iter != nodes.end(); ++iter) {
        enumeration.insert(pair<node*,int>(*iter, count++));
    }
    count = 0;
    for (iter = nodes.begin(); iter != nodes.end(); ++iter) {
        double x = cos(alpha*count);
        double y = sin(alpha*count);

        os << ++count << ' ' << x << ' ' << y;
        foreach_neighbor(*iter, edgename, i) {
            node* m = i->second.first;
            map<node*, int>::iterator pos = enumeration.find(m);

            if (pos != enumeration.end()) os << ' ' << pos->second + 1;
        }
        os << endl;
    }
}

void graph::write_mma(ostream& os)
{
    map<node*, int> enumeration;
    int index = 0;
    bool first = true;

    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        enumeration.insert(pair<node*, int>(*iter, index++));
    }
    os << '{';
    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        int nindex = enumeration.find(n)->second;
        
        forall_neighbors(n, niter) {
            int nnindex = enumeration.find(neighbor_node(niter))->second;
            int eindex = neighbor_index(niter);

            if (first) first = false; else os << ',' << endl;
            os << '{' << nindex << ',' << nnindex << ',' << '\"' << atom::get_name(eindex) << '\"' << '}';
        }

    }
    os << '}';
}

graph* graph::get_induced_copy(const copy_predicate_t& pred)
{
    graph* result = (graph*)make_instance();
    map<node*, node*> cmap;

    for (list<node*>::iterator i = nodes.begin(); i != nodes.end(); ++i) {
        node* n = *i;

        if (pred(n)) cmap.insert(pair<node*, node*>(n, result->add_node()));
    }
    for (map<node*, node*>::iterator i = cmap.begin(); i != cmap.end(); ++i) {
        if (i->second != nullptr) i->first->gc_clone_to(i->second, result, cmap);
    }
    copy_to(result, cmap);
    return result;
}

void graph::copy_to(graph* dest, const map<node*, node*>& cmap)
{
}

void graph::delete_nodes()
{
	for (iter_t it = nodes.begin(); it != nodes.end(); it++) {
		delete (*it);
	}
	nodes.clear();
}

void graph::delete_nodes(unsigned attr)
{
    list<node*>::iterator iter;

    for (iter = nodes.begin(); iter != nodes.end(); ++iter) 
        (*iter)->delete_edges(attr);

    iter = nodes.begin();
    while (iter != nodes.end()) {
        if (!(*iter)->is_attr_set(attr)) ++iter; 
        else {
            delete *iter;
            iter = nodes.erase(iter);
        } 
    }
}

// Delete all nodes with the attribute attr.
// Remark: We assume that these nodes are (semi)isolated, i.e. there are no 
//   edges pointing to these nodes!
void graph::delete_isolated_nodes(unsigned attr)
{
    list<node*>::iterator iter = nodes.begin();

    while (iter != nodes.end()) {
        if (!(*iter)->is_attr_set(attr)) ++iter; 
        else {
            delete *iter;
            iter = nodes.erase(iter);
        } 
    }
}

// We assume ns is a set of (semi)isolated nodes, i.e. there are no 
//   edges pointing to these nodes!
void graph::delete_isolated_nodes(const set<node*>& ns)
{
    list<node*>::iterator iter = nodes.begin();

    while (iter != nodes.end()) {
        if (ns.find(*iter) == ns.end()) ++iter;
        else {
            delete *iter;
            iter = nodes.erase(iter);
        }
    }
}

void graph::delete_nodes(const set<node*>& ns) 
{
    list<node*>::iterator iter = nodes.begin();

    while (iter != nodes.end()) {
        if (ns.find(*iter) == ns.end()) ++iter;
        else iter = nodes.erase(iter);
    }
    for (iter = nodes.begin(); iter != nodes.end(); ++iter) {
        (*iter)->delete_edges(ns);
    }
    for (set<node*>::const_iterator siter = ns.begin(); siter != ns.end(); ++siter) {
        delete *siter;
    }
}

void graph::delete_node(node* n)
{
    list<node*>::iterator iter = nodes.begin();

    while (iter != nodes.end()) {
        if (*iter == n) iter = nodes.erase(iter);
        else {
            (*iter)->delete_edges(n);
            ++iter;
        }
    }
    delete n;
}

void graph::delete_edges(int index)
{
    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter)
        (*iter)->delete_edges(index);
}


void graph::delete_edges_complement(int index)
{
    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter)
        (*iter)->delete_edges_complement(index);
}

void graph::rename_edges(int from, int to)
{
    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) 
        (*iter)->rename_edges(from, to);
}

// Returns a bfs tree
// root - root node of tree 
// maxdist - depth of tree
// predicate - should return true on allowed edges
// tree - result; node data are of type pair<node*,int> where pair.first is a pointer to the original node and
//                pair.second is edge name pointing to this node
void graph::bfs_subtree(graph& tree, node* root, int maxdist, const bfs_predicate& predicate)
{
    typedef pair<int,node*> inpair;
    typedef pair<node*,int> nd_pair;
    typedef node_data_tx<nd_pair> tree_data;

    queue<inpair> q;
    list<node*> visited;

    node* tn = tree.add_node(new tree_data(nd_pair(root, -1)));
    q.push(inpair(0, tn));
    root->set_attr(ATTR_VISITED);
    visited.push_back(root);

    while (!q.empty()) {
        inpair& p = q.front();
        node* n = ((tree_data*)p.second->data)->data.first;

        if (p.first < maxdist) {
            int nbcount = 0;
            forall_neighbors(n, iter) {
                if (predicate.value(iter->first, p.first, nbcount)) { //if (allowededges.empty() || allowededges.find(iter->first) != allowededges.end()) {
                    node* m = iter->second.first;

                    if (!m->is_attr_set(ATTR_VISITED)) {
                        tn = tree.add_node(new tree_data(nd_pair(m, iter->first)));
                        tree.add_edge(p.second, tn, 1);
                        q.push(inpair(p.first + 1, tn));
                        m->set_attr(ATTR_VISITED);
                        visited.push_back(m);
                    }
                    ++nbcount;
                }
            }
            
        }
        q.pop();
    }

    for (list<node*>::iterator iter = visited.begin(); iter != visited.end(); ++iter) {
        (*iter)->clear_attr(ATTR_VISITED);
    }
}

struct subset_predicate {
    int maximum;

    subset_predicate(int m) : maximum(m) { }
    bool operator()(const list<node*>& s, node* n) const { return (int)s.size() < maximum; }
};

// Returns a list of bfs subtrees 
// root - root node of trees 
// maxdist - maximal depth of trees
// mindist - minimal depth of trees
// predicate - should return true on allowed edges (i.e. it can check edge name)
// maxdegree - maximal degree of a subtree
// tree - result; node data are of type pair<node*,int> where pair.first is a pointer to the original node and
//                pair.second is edge name pointing to this node
void graph::bfs_subtrees(list<graph>& trees, node* root, int maxdist, int mindist, int maxdegree,
                         const bfs_predicate& predicate)
{
    typedef set<graph*> gpset;
    typedef triple<int,node*,gpset> qtriple;
    typedef pair<node*,int> nd_pair;
    typedef node_data_tx<nd_pair> tree_data;

    list<qtriple> q;
    list<node*> visited;
    list<qtriple>::iterator qiter;

    trees.clear();
    trees.push_back(graph(true));
    trees.back().add_node(root);
    q.push_back(qtriple(0, root, gpset()));
    q.front().third.insert(&trees.back());
    root->set_attr(ATTR_VISITED);
    visited.push_back(root);
    subset_predicate ssp(maxdegree);

    while (!q.empty()) {
        qtriple& p = q.front();
        node* n = p.second;
        map<node*,gpset*> gpsmap;

        if (p.first < maxdist) {
            subset_creator<node*> ssc;
            int qsize = (int)q.size();

            forall_neighbors(n, iter) {
                if (predicate.value(iter->first, p.first, 0)) {
                    node* m = iter->second.first;

                    if (!m->is_attr_set(ATTR_VISITED)) {
                        m->set_attr(ATTR_VISITED);
                        visited.push_back(m);
                        ssc.add(m, ssp);
                        q.push_back(qtriple(p.first + 1, m, gpset()));
                        gpsmap.insert(pair<node*,gpset*>(m, &(q.back().third)));
                    }
                }
            }
            for (gpset::iterator giter = p.third.begin(); giter != p.third.end(); ++giter) {
                graph* gp = *giter;
                gpset newpset;

                for (subset_creator<node*>::iter_t ssiter = ssc.subsets.begin(); 
                        ssiter != ssc.subsets.end(); ++ssiter) {
                    list<node*>& sst = *ssiter;
                    
                    if (sst.size() == 0) continue;
                    trees.push_back(graph(*gp));

                    graph& newtree = trees.back();

                    for (list<node*>::iterator niter = sst.begin(); niter != sst.end(); ++niter) {
                        node* nx = *niter;
                        map<node*,gpset*>::iterator gpsiter;

                        newtree.add_node(nx);
                        if ((gpsiter = gpsmap.find(nx)) != gpsmap.end()) 
                            gpsiter->second->insert(&newtree);
                    }
                    newpset.insert(&newtree);
                }
                qiter = q.begin();
                for (int i = 1; i < qsize; ++i) {
                    ++qiter;

                    gpset& qgp = qiter->third;

                    if (qgp.find(gp) != qgp.end()) 
                        qgp.insert(newpset.begin(), newpset.end());
                }
                
            }

        }
        if (n == root) --ssp.maximum;
        q.pop_front();
    }

    for (list<node*>::iterator iter = visited.begin(); iter != visited.end(); ++iter) {
        (*iter)->clear_attr(ATTR_VISITED);
    }
}


/*bool graph::complete_from_point(node* n, list<node*>& result)
{
    list<node*> visited;

    n->set_attr(ATTR_CHECKED);
    result.clear();
    visited.push_back(n);
    result.push_back(n);
    complete_from_point_r(n, result, visited);
}

bool graph::complete_from_point(node* n, list<node*>& result, list<node*>& visited)
{
}

void graph::find_complete_graphs(list<list<node*> >& cglist, int minn) const
{
    cglist.clear();

    for (citer_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;

        forall_neighbors(n, niter) {
            node* nn = *niter;

            
        }
    }
}

*/

void graph::copy_to(streamable* p, cloner& cl)
{
    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        ((graph*)p)->nodes.push_back((node*)cl.get_copy(*iter));
    }
}

void graph::read_from_stream(istreamer &is)
{
    streamable::read_from_stream(is);

    unsigned nnodes;
    streamable* n;
    
    is.read(nnodes);
    for (unsigned i = 0; i < nnodes; ++i) {
        is.read(n);
        nodes.push_back((node*)n);
    }
}

void graph::write_to_stream(ostreamer& os)
{
    streamable::write_to_stream(os);

    unsigned nnodes = (unsigned)nodes.size();

    os.write(nnodes);
    for (list<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter)
        os.write(*iter);
}

// Recursively visits all nodes starting from node 'n', following edges with 'edge_name'.
// End nodes are returned in 'result' and all visited nodes are linked to 'result' with 'link_name'.
// If there is an edge named 'link_name' then this node is considered to be already "visited".
void graph::recurse_and_link(node* n, int edge_name, int link_name, set<node*>& result)
{
    if (n->has_neighbor(link_name)) 
        n->add_to_neighbor_set(link_name, result);
    else {
        foreach_neighbor(n, edge_name, iter) {
            set<node*> nset;

            recurse_and_link(neighbor_node(iter), edge_name, link_name, nset);
            result.insert(nset.begin(), nset.end());
        }
        if (result.empty()) result.insert(n);
        else add_edges(n, result.begin(), result.end(), 0, link_name);
    }
}


void graph::recurse_and_link(int layer, int edge_name, int link_name)
{
    for (auto niter = nodes.begin(); niter != nodes.end(); ++niter) {
        node* n = *niter;

        if (z_coo(n) == layer) recurse_and_link(n, edge_name, link_name);
    }
}


// See more general version of the function above.
void graph::recurse_and_link(node* n, int edge_name, int link_name)
{
    set<node*> result;

    recurse_and_link(n, edge_name, link_name, result);
}


// Finds a 'path' between 'srcn' and 'destn' using edges with name 'edge_name'.
// Use it with care ;)
bool graph::find_path(vector<node*>& path, node* srcn, int edge_name, node* destn) const
{
    list<vector<node*> > paths;

    paths.push_back(vector<node*>());
    paths.back().push_back(srcn);

    if (srcn == destn) {
        path = paths.back(); 
        return true;
    }
    while (!paths.empty()) {
        size_t lsize = paths.size();
 
        for (size_t i = 0; i < lsize; ++i) {
            vector<node*>& p = paths.front();
            node* n = p.back();

            foreach_neighbor(n, edge_name, niter) {
                node* nn = neighbor_node(niter);

                if (nn == destn) {
                    path = p;
                    path.push_back(nn);
                    return true;
                }
                paths.push_back(p);
                paths.back().push_back(nn);
            }
            paths.pop_front();
        }
    }
    return false;
}

void graph::set_attr(unsigned attr)
{
    set_attr(nodes.begin(), nodes.end(), attr);
}

void graph::clear_attr(unsigned attr)
{
    clear_attr(nodes.begin(), nodes.end(), attr);
}


/*bool graph::is_subgraph_of(const graph& g) const
{
    if (nodes.size() > g.nodes.size()) return false;

    citer_t iter1, iter2 = nodes.begin();
    vector<int> nodemap(nodes.size(), -1);

    while (iter != nodes.end()) {
        node* n = *iter;
        
        if (match_n2n(
    }
    
    
}*/

// global functions
///////////////////////////////////////////////////////////////////////////////

// Bron_Kerbosch algorithm for finding maximal cliques, with pivoting; 
// see http://arxiv.org/PS_cache/arxiv/pdf/1006/1006.5440v1.pdf, Figure 2.
// (See also http://en.wikipedia.org/wiki/Bron-Kerbosch_algorithm.)
void Bron_Kerbosch(list<set<node*> >& cliques, const set<node*>& P, const set<node*>& R, const set<node*>& X)
{
    if (P.empty() && X.empty()) {
        cliques.push_back(R);
        return;
    }

    // Choose pivot
    set<node*> S;
    node* u = nullptr;  // choose u to maximize |P \cap N(u)|
    set<node*> N;
    int maxint = -1;   

    S.insert(P.begin(), P.end());
    S.insert(X.begin(), X.end());
    for (set<node*>::iterator siter = S.begin(); siter != S.end(); ++siter) {
        (*siter)->get_neighbors(N);
        
        int isize = intersection_size(P, N);

        if (isize > maxint) { maxint = isize; u = *siter; }
    }

    S.clear();
    u->get_neighbors(N);
    set_difference(S, P, N);

    set<node*> PP = P;
    set<node*> XX = X;
    for (set<node*>::iterator siter = S.begin(); siter != S.end(); ++siter) {
        set<node*> newP, newR, newX;
        node* v = *siter;

        v->get_neighbors(N);
        set_intersection(newP, PP, N);
        newR.insert(R.begin(), R.end()); newR.insert(v);
        set_intersection(newX, XX, N);
        Bron_Kerbosch(cliques, newP, newR, newX);
        PP.erase(v);
        XX.insert(v);
    }
}

// Returns a list of sets of nodes representing maximal cliques.
// Uses Bron-Kerbosch algorithm; see function Bron_Kerbosch.
void get_maximal_cliques(list<set<node*> >& cliques, graph* g)
{
    cliques.clear();
    Bron_Kerbosch(cliques, set<node*>(g->nodes.begin(), g->nodes.end()), set<node*>(), set<node*>());
}

