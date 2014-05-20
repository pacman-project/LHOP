// graph library includes definitions found in
//  - graph.h
//  - graph_edges.h
//  - graph_nodes.h
//  - img_graph.h
////////////////////////////////////////////////


//#define _USE_MATH_DEFINES  // necessary for <math.h> to define M_PI,.... constants, etc. 
//#include <math.h>

#include "graph_nodes.h"
#include "graph_edges.h"
#include "graph.h"

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

