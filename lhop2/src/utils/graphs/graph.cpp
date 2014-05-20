// graph library includes definitions found in
//  - graph.h
//  - graph_edges.h
//  - graph_nodes.h
//  - img_graph.h
////////////////////////////////////////////////


#define _USE_MATH_DEFINES  // necessary for <math.h> to define M_PI,.... constants, etc. 
#include <math.h>
#include <queue>

#include "graph.h"

#include "utils/utils.h"
#include "utils/serialization/streaming.h"

// graph
///////////////////////////////////////////////////////////////////////////////


graph::~graph()
{
    if (!is_subgraph) delete_nodes();
}

edge_data* graph::get_edge_data(node* n1, node* n2, int name)
{
    foreach_neighbor(n1, name, iter) {
        if (neighbor_node(iter) == n2) return neighbor_edge_data(iter);
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

// See more general version of the function above.
void graph::recurse_and_link(node* n, int edge_name, int link_name)
{
    set<node*> result;

    recurse_and_link(n, edge_name, link_name, result);
}

void graph::set_attr(unsigned attr)
{
    set_attr(nodes.begin(), nodes.end(), attr);
}

void graph::clear_attr(unsigned attr)
{
    clear_attr(nodes.begin(), nodes.end(), attr);
}


// global functions
///////////////////////////////////////////////////////////////////////////////

