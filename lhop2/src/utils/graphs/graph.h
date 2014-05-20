// graph class
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "utils/serialization/streaming.h"
#include "utils/misc.h"
#include "utils/smemory.h"

#include "graph_nodes.h"
#include "graph_edges.h"

using namespace std;

// helper methods to allow node/neighbor traversal
//////////////////////////////////////////////////

#define node_neighbor_iterator(_iter) multimap<int, node_pair>::iterator _iter;
#define forall_neighbors(_nodep, _iter)\
    for (multimap<int, edge_pair>::iterator _iter = (_nodep)->neighbors.begin();\
    _iter != (_nodep)->neighbors.end();\
    ++_iter)
#define foreach_neighbor(_nodep, _index, _iter)\
    for (multimap<int, edge_pair>::iterator _iter = (_nodep)->get_neighbors_begin(_index);\
    (_iter != (_nodep)->neighbors.end()) && (_iter->first == _index);\
    ++_iter)
#define while_neighbor_has_index(_nodep, _index, _iter)\
    for (multimap<int, edge_pair>::iterator _iter = (_nodep)->get_neighbors_begin(_index);\
    (_iter != (_nodep)->neighbors.end()) && (_iter->first == _index); )
#define neighbor_node(_iter) ((_iter)->second.first)
#define neighbor_node_data(_iter) ((_iter)->second.first->data)
#define neighbor_index(_iter) ((_iter)->first)
#define neighbor_edge_data(_iter) ((_iter)->second.second)


// graph
//////////

class graph : public streamable, public custom_alloc {
public:
    omp_lock lock;

    struct copy_predicate_t {
        virtual bool operator()(node*) const { return true; }
    };

    typedef list<node*>::iterator iter_t;
    typedef list<node*>::const_iterator citer_t;

    list<node*> nodes;
    bool is_subgraph;

    graph(bool subgr = false) : nodes(), is_subgraph(subgr) { }
    graph(const graph& g) : nodes(g.nodes), is_subgraph(true) { }
    virtual ~graph();

    void add_node(node* n) { nodes.push_back(n); }
    node* add_node(node_data* nd, unsigned attr = 0) 
        { add_node(new node(nd, attr)); return nodes.back(); }
    node* add_node() { add_node(new node()); return nodes.back(); }

    bool empty() const { return nodes.empty(); }

    virtual void delete_nodes();
    virtual void delete_nodes(unsigned attr);
    virtual void delete_isolated_nodes(unsigned attr);
    virtual void delete_isolated_nodes(const set<node*>& ns);
    template<class I> void delete_nodes(I begin, I end);
    virtual void delete_nodes(const set<node*>& ns);
    virtual void delete_node(node* n);
    template<class C> static void remove_nodes(C& container, unsigned attr);
    void delete_edges(int index);
    template<class I> void delete_edges(I begin, I end, int name);
    void delete_edges_complement(int index);
    void rename_edges(int from, int to);

    static void add_edge(node* n1, node* n2, int name);
    bool add_edge_unique(node* n1, node* n2, int name);
    void add_edge(node* n1, node* n2, int name, int back_name);
    void add_edge_unique(node* n1, node* n2, int name, int back_name);
    bool add_edge_2_unique(node* n1, node* n2, edge_data* d, int name);
    static void add_edge_2(node* n1, node* n2, edge_data* d, int name);
    void add_edge_2(node* n1, node* n2, edge_data* d, edge_data* bd, int name, int back_name);
    template<class I> void add_edges_to_neighbors_of(node* n, I begin, I end, int name);
    template<class I> void add_edges(node* n, I begin, I end, unsigned attr, int name);

    edge_data* get_edge_data(node* n1, node* n2, int name);
    void add_reverse_edges(node* n, int name, int rname);
    void add_reverse_edges(int name, int rname);
    bool connected(node* n1, node* n2, int name) 
        { return n1->is_neighbor(n2, name); }

    template<class Man> void recurse2(Man& man, const typename Man::container_t& list, int edge_name);
    template<class Man> void recurse2(Man& man, const typename Man::container_t& list, const vector<int>& edge_names);
    template<class Cont> void recurse(const Cont& list, int edge_name, set<node*>& result);
    template<class Cont, class Cond> void recurse(const Cont& list, int edge_name, 
        Cond& check, set<node*>& result);
    template<class Cont, class Cond> void recurse0(const Cont& list, int edge_name, 
        Cond& check, set<node*>& result);
    template<class Cont, class Cond> void recurse05(const Cont& list, int edge_name, 
        const Cond& check, set<node*>& result);
    template<class Cont, class Cond> void recurse3(const Cont& list, const Cond& sel_neighbors, 
        set<node*>& result);
    void recurse_from_node(node* src, int edge_name, set<node*>& result);
    template<class Cond> void recurse_from_node(node* src, int edge_name, Cond check, 
        set<node*>& result);
    template<class Cond> void recurse0_from_node(node* src, int edge_name, Cond check, 
        set<node*>& result);
    void recurse_and_link(node* n, int edge_name, int link_name, set<node*>& result);
    void recurse_and_link(node* n, int edge_name, int link_name);
    
    template<class Filter> void recurse_and_link(node* n, int edge_name, int link_name, set<node*>& result);
    template<class Cont> void recurse_and_link(const Cont& list, int edge_name, int link_name, set<node*>& result);

    void find_complete_graphs(list<list<node*> >& cglist, int minn) const;
    template<class I> static void get_neighbors(set<node*>& result, I begin, I end, int edge_name);
    template<class I> void partition_by_neighbors(set<node*>& children, set<node*>& childless, 
        I begin, I end, int edge_name);
    template<class I> void get_childless(set<node*>& result, I begin, I end, int edge_name);
    template<class I, class P> void select_by_neighbors(set<node*>& result, I begin, I end, 
        int edge_name, const P& pred);
    template<class I, class F> void apply(I begin, I end, const F& f);
    template<class I> void set_attr(I begin, I end, unsigned attr);
    void set_attr(unsigned attr);
    template<class I> void clear_attr(I begin, I end, unsigned attr);
    void clear_attr(unsigned attr);
    template<class Cont> void subgraph(const Cont& list, int edge_name, set<node*>& result);
    void subgraph_from_node(node* n, int edge_name, set<node*>& result);

    virtual void copy_to(streamable* p, cloner& cl);
    virtual streamable* make_instance() const { return new graph(); }
    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);

protected:
    bool complete_from_point_r(node* n, list<node*>& result);
    bool complete_from_point(node* n, list<node*>& result);

    virtual void copy_to(graph* dest, const map<node*, node*>& cmap);

};

// path_walker 
////////////////

// To be used as "Man" in graph::recurse2 to recursively get all 
// paths form a set of nodes along edges with a certain name.

class path_walker {
public:
    typedef list<path_t> container_t; 
    typedef container_t::iterator iterator_t;
    typedef container_t::const_iterator const_iterator_t;

    list<path_t> result;

    path_walker() : result() { }

    void reset() { result.clear(); }

    node* get_node(const container_t::value_type& val) { return val.back().first; }

    bool check_neighbor(container_t& neighbors, const container_t::value_type& n, edge_data* ed, node* nn)
    {
        neighbors.push_back(n);
        neighbors.back().push_back(edge_pair(nn, ed));
        return true;
    }

    void insert_leaf(const container_t::value_type& n) { result.push_back(n); }

    template<class T> static void add_to_container(container_t& l, T begin, T end)
    {
        while (begin != end) { 
            l.push_back(path_t()); 
            l.back().push_back(edge_pair(*begin, nullptr));
            ++begin;
        }
    }

    static void add_to_container(container_t& l, node* n)
    {
        l.push_back(path_t());
        l.back().push_back(edge_pair(n, nullptr));
    }

};


// inline definitions
///////////////////////////////////////////////////////////////////////////////

// graph
//////////

inline void graph::add_edge(node *n1, node *n2, int name)
{
    n1->neighbors.insert(node_pair(name, edge_pair(n2, nullptr)));
}

inline bool graph::add_edge_unique(node *n1, node *n2, int name)
{
    if (n1->is_neighbor(n2, name)) return false; 
    add_edge(n1, n2, name);
    return true;
}

inline bool graph::add_edge_2_unique(node *n1, node *n2, edge_data* d, int name)
{
    if (n1->is_neighbor(n2, name)) return false; 
    add_edge_2(n1, n2, d, name);
    return true;
}


inline void graph::add_edge(node *n1, node *n2, int name, int back_name)
{
    n1->neighbors.insert(node_pair(name, edge_pair(n2, nullptr)));
    n2->neighbors.insert(node_pair(back_name, edge_pair(n1, nullptr)));
}

inline void graph::add_edge_unique(node *n1, node *n2, int name, int back_name)
{
    add_edge_unique(n1, n2, name);
    add_edge_unique(n2, n1, back_name);
}

inline void graph::add_edge_2(node* n1, node* n2, edge_data* d, int name)
{
    n1->neighbors.insert(node_pair(name, edge_pair(n2, d)));
}

inline void graph::add_edge_2(node* n1, node* n2, edge_data* d, edge_data* bd, 
                              int name, int back_name)
{
    n1->neighbors.insert(node_pair(name, edge_pair(n2, d)));
    n2->neighbors.insert(node_pair(back_name, edge_pair(n1, bd)));
}

template<class I> void graph::add_edges_to_neighbors_of(node* n, I begin, I end, int name)
{
    set<node*> nbset;

    for (; begin != end; ++begin) {
        node* sn = *begin;

        foreach_neighbor(sn, name, iter) {
            node* nn = neighbor_node(iter);
            if (nbset.insert(nn).second) add_edge(n, nn, name);
        }
    }
}

template<class I> void graph::add_edges(node* n, I begin, I end, unsigned attr, int name)
{
    set<node*> nbset;

    for (; begin != end; ++begin) {
        node* sn = *begin;

        if (sn->is_attr_set(attr)) add_edge(n, sn, name);        
    }
}

template<class I> void graph::delete_nodes(I begin, I end)
{
    set<node*> ns(begin, end);
    delete_nodes(ns);    
}

template<class C> void graph::remove_nodes(C& container, unsigned attr)
{
    typename C::iterator iter = container.begin();

    while (iter != container.end()) {
        if ((*iter)->is_attr_set(attr)) iter = container.erase(iter);
        else ++iter;
    }
}

template<class I> void graph::delete_edges(I begin, I end, int name)
{
    for (; begin != end; ++begin) {
        (*begin)->delete_edges(name);
    }
}

template<class Cont> void graph::recurse(const Cont& list, 
    int edge_name, set<node*>& result)
{
    set<node*> neighbors;
    
    for (typename Cont::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
        node* n = *iter;
        int count = 0;

        foreach_neighbor(n, edge_name, i) {
            neighbors.insert(i->second.first);
            ++count;
        }
        if (count == 0) result.insert(n);
    }
    if (neighbors.size() > 0)
        recurse(neighbors, edge_name, result);
}

template<class Cont, class Cond> void graph::recurse0(const Cont& list, 
    int edge_name, Cond& check, set<node*>& result)
{
    set<node*> neighbors;
    
    for (typename Cont::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
        node* n = *iter;
        int count = 0;

        foreach_neighbor(n, edge_name, i) {
            neighbors.insert(i->second.first);
            ++count;
        }
        if (count == 0 && check(n)) result.insert(n);
    }
    if (neighbors.size() > 0)
        recurse0(neighbors, edge_name, check, result);
}

template<class Cont, class Cond> void graph::recurse05(const Cont& list, 
    int edge_name, const Cond& check, set<node*>& result)
{
    set<node*> neighbors;
    
    for (typename Cont::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
        node* n = *iter;
        int count = 0;

        foreach_neighbor(n, edge_name, i) {
            node* nn = i->second.first;

            if (check(n, nn)) {
                neighbors.insert(nn);
                ++count;
            }
        }
        if (count == 0 && check(n)) result.insert(n);
    }
    if (neighbors.size() > 0)
        recurse05(neighbors, edge_name, check, result);
}

template<class Cont, class Cond> void graph::recurse(const Cont& list, 
    int edge_name, Cond& check, set<node*>& result)
{
    set<node*> neighbors;
    
    for (typename Cont::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
        node* n = *iter;
        int count = 0;

        foreach_neighbor(n, edge_name, i) {
            node* nn = i->second.first;

            if (check(n, nn)) {
                neighbors.insert(nn);
                ++count;
            }
        }
        if (count == 0) result.insert(n);
    }
    if (neighbors.size() > 0)
        recurse(neighbors, edge_name, check, result);
}

// recurse2 Man
//  typedef container_t
//  node* get_node(container_t::const_iterator) - returns the node part of the data
//  bool check_neighbor(const container_t&, node* n, node* nn) - returns true if the 
//     neighbor nn is "good" and adds it to the container
//  void insert_leaf(node* n) - insert n to the result
//
template<class Man> void graph::recurse2(Man& man, const typename Man::container_t& list, 
    int edge_name)
{
    typename Man::container_t neighbors;
    
    for (typename Man::const_iterator_t iter = list.begin(); iter != list.end(); ++iter) {
        const typename Man::container_t::value_type& item = *iter;
        node* n = man.get_node(item);
        int count = 0;

        foreach_neighbor(n, edge_name, i) {
            if (man.check_neighbor(neighbors, item, neighbor_edge_data(i), neighbor_node(i))) 
                ++count;
        }
        if (count == 0) man.insert_leaf(item);
    }
    if (neighbors.size() > 0)
        recurse2(man, neighbors, edge_name);
}

// recurse2 Man
//  typedef container_t
//  node* get_node(container_t::const_iterator) - returns the node part of the data
//  bool check_neighbor(const container_t&, node* n, node* nn) - returns true if the 
//     neighbor nn is "good" and adds it to the container
//  void insert_leaf(node* n) - insert n to the result
//
template<class Man> void graph::recurse2(Man& man, const typename Man::container_t& list, 
    const vector<int>& edge_names)
{
    typename Man::container_t neighbors;
    
    for (typename Man::const_iterator_t iter = list.begin(); iter != list.end(); ++iter) {
        const typename Man::container_t::value_type& item = *iter;
        node* n = man.get_node(item);
        int count = 0;

        for (vector<int>::const_iterator niter = edge_names.begin(); niter != edge_names.end(); ++niter) {
            foreach_neighbor(n, *niter, i) {
                if (man.check_neighbor(neighbors, item, neighbor_edge_data(i), neighbor_node(i))) 
                    ++count;
            }
        }
        if (count == 0) man.insert_leaf(item);
    }
    if (neighbors.size() > 0)
        recurse2(man, neighbors, edge_names);
}

template<class Cont, class Cond> void graph::recurse3(const Cont& list, 
    const Cond& sel_neighbors, set<node*>& result)
{
    set<node*> neighbors;
    
    for (typename Cont::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
        node* n = *iter;
        int count = sel_neighbors(neighbors, n);

        if (count == 0) result.insert(n);
    }
    if (neighbors.size() > 0)
        recurse3(neighbors, sel_neighbors, result);
}

// Recursively visits all nodes starting from nodes in 'list', following edges with 'edge_name'.
// End nodes are returned in 'result' and all visited nodes are linked to 'result' with 'link_name'.
// If there is an edge named 'link_name' then this node is considered to be already "visited".
template<class Cont> void graph::recurse_and_link(const Cont& list, 
    int edge_name, int link_name, set<node*>& result)
{
    set<node*> neighbors;
    
    for (typename Cont::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
        node* n = *iter;
        int count = 0;

        foreach_neighbor(n, edge_name, i) {
            neighbors.insert(i->second.first);
            ++count;
        }
        if (count == 0) result.insert(n);
    }
    if (neighbors.size() > 0)
        recurse(neighbors, edge_name, result);
}

// Recursively visits all nodes starting from node 'n', following edges with 'edge_name'.
// End nodes are returned in 'result' and all visited nodes are linked to 'result' with 'link_name'.
// If there is an edge named 'link_name' then this node is considered to be already "visited".
template<class Filter> 
void graph::recurse_and_link(node* n, int edge_name, int link_name, set<node*>& result)
{
    if (n->has_neighbor(link_name)) 
        n->add_to_neighbor_set(link_name, result);
    else {
		Filter result_filter;
		//set<node*> result_filter;
        foreach_neighbor(n, edge_name, iter) {
            set<node*> nset;

            recurse_and_link(neighbor_node(iter), edge_name, link_name, nset);			
            result_filter.insert(nset.begin(), nset.end());
        }
		if (result_filter.empty()) result.insert(n);
		else {
			add_edges(n, result_filter.begin(), result_filter.end(), 0, link_name);
			result.insert(result_filter.begin(), result_filter.end());
		}
    }
}

template<class Cont> void graph::subgraph(const Cont& list, int edge_name, set<node*>& result)
{
    set<node*> neighbors;
    
    result.insert(list.begin(), list.end());
    for (typename Cont::const_iterator iter = list.begin(); iter != list.end(); ++iter) {
        node* n = *iter;

        foreach_neighbor(n, edge_name, i) {
            neighbors.insert(neighbor_node(i));
        }
    }
    if (neighbors.size() > 0) 
        subgraph(neighbors, edge_name, result);
}

inline void graph::subgraph_from_node(node* n, int edge_name, set<node*>& result)
{
    vector<node*> singleton;
    singleton.push_back(n);
    subgraph(singleton, edge_name, result);
}

inline void graph::recurse_from_node(node* src, int edge_name, set<node*>& result)
{
    vector<node*> singleton;
    singleton.push_back(src);
    recurse(singleton, edge_name, result);
}

template<class Cond> void graph::recurse_from_node(node* src, int edge_name, Cond check, set<node*>& result)
{
    vector<node*> singleton;
    singleton.push_back(src);
    recurse(singleton, edge_name, check, result);
}

template<class Cond> void graph::recurse0_from_node(node* src, int edge_name, Cond check, set<node*>& result)
{
    vector<node*> singleton;
    singleton.push_back(src);
    recurse0(singleton, edge_name, check, result);
}

template<class I> void graph::get_neighbors(set<node*>& result, I begin, I end, int edge_name)
{
    while (begin != end) {
        node* n = *begin;

        foreach_neighbor(n, edge_name, iter) {
            result.insert(neighbor_node(iter));
        }
        ++begin;
    }
}

template<class I> void graph::partition_by_neighbors(set<node*>& children, set<node*>& childless, 
        I begin, I end, int edge_name)
{
    while (begin != end) {
        node* n = *begin;
        int count = 0;

        foreach_neighbor(n, edge_name, iter) {
            children.insert(neighbor_node(iter));
            ++count;
        }
        if (count == 0) childless.insert(n);
        ++begin;
    }    
}

template<class I> void graph::get_childless(set<node*>& result, I begin, I end, int edge_name)
{
    while (begin != end) {
        node* n = *begin;
        if (n->degree(edge_name) == 0) result.insert(n);
    }
}

template<class I, class P> void graph::select_by_neighbors(set<node*>& result, I begin, I end, 
    int edge_name, const P& pred)
{
    while (begin != end) {
        node* n = *begin;

        foreach_neighbor(n, edge_name, iter) {
            if (pred(neighbor_node(iter))) result.insert(n);
        }
        ++begin;
    }
}

template<class I, class F> void graph::apply(I begin, I end, const F& f)
{
    while (begin != end) {
        f(*begin);
        ++begin;
    }
}

template<class I> void graph::set_attr(I begin, I end, unsigned attr)
{
    for (; begin != end; ++begin) (*begin)->set_attr(attr);
}

template<class I> void graph::clear_attr(I begin, I end, unsigned attr)
{
    for (; begin != end; ++begin) (*begin)->clear_attr(attr);
}


#endif /* _GRAPH_H_ */

