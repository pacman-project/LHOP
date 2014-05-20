// graph nodes classes
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _GRAPH_NODES_H_
#define _GRAPH_NODES_H_

////#include "utils/atom.h"
#include "utils/serialization/streaming.h"
#include "utils/smemory.h"

using namespace std;

class graph;

class node;

class edge_data;

const unsigned NODE_DELETED_ATTR = 0x8000;

class node_data : public streamable, public custom_alloc {
public:
    node_data() { }

    virtual void clone(const node_data* d) { }
    node_data* get_clone() 
    { 
        node_data* res = new node_data(); 
        res->clone(this); 
        return res;
    }
    virtual void gc_clone(const node_data* d, const map<node*,node*>& cmap) { }
    node_data* get_gc_clone(const map<node*,node*>& cmap)
    {
        node_data* res = (node_data*)make_instance();
        res->gc_clone(this, cmap);
        return res;
    }

    virtual void copy_to(streamable* p, cloner& cl) { /* streamable::copy_to(p, cl); */ }
    virtual streamable* make_instance() const { return new node_data(); }
    virtual void read_from_stream(istreamer& is) { streamable::read_from_stream(is); }
    virtual void write_to_stream(ostreamer& os) { streamable::write_to_stream(os); }

    virtual void print(ostream& os) const { os << "-- Empty --"; } 

    friend ostream& operator<<(ostream& os, const node_data& d) { d.print(os); return os; }
};

template<class T> class node_data_t : public node_data {
public:
    T data;
    
    node_data_t() { }
    node_data_t(const T& d) : data(d) { }

    virtual void clone(const node_data* d) { data = ((node_data_t*)d)->data; }
    virtual void gc_clone(const node_data* d, const map<node*, node*>& cmap) { data = ((node_data_t*)d)->data; }

    virtual void copy_to(streamable* p, cloner& cl) 
    { 
        node_data::copy_to(p, cl);
        ((node_data_t*)p)->data = data;
    }

    virtual streamable* make_instance() const { return new node_data_t<T>(); }
    virtual void read_from_stream(istreamer& is)  
    { 
        node_data::read_from_stream(is); 
        is.read(data); 
    }
    
    virtual void write_to_stream(ostreamer& os) 
    { 
        node_data::write_to_stream(os); 
        os.write(data);
    }

    virtual void print(ostream& os) const
    {
        os << data << endl;
    }

};


template<class T> class node_data_tx : public node_data {
public:
    T data;
    
    node_data_tx() { }
    node_data_tx(const T& d) : data(d) { }

    virtual void clone(const node_data* d) { data = ((node_data_tx*)d)->data; }

    virtual streamable* make_instance() const 
        { throw streaming_exception(__FILE__, __LINE__,streaming_exception::NONSTREAMABLE_CLASS, "node_data_tx"); }
};

// node
/////////

const unsigned ATTR_VISITED = 1024;
const unsigned ATTR_CHECKED = 2048;

typedef pair<node*, edge_data*> edge_pair;
typedef pair<int, edge_pair> node_pair;
typedef vector<edge_pair> path_t; 

class node : public streamable, public custom_alloc {
public:
    typedef multimap<int, edge_pair> container;
    typedef multimap<int, edge_pair>::iterator iterator;
    typedef multimap<int, edge_pair>::const_iterator const_iterator;

    multimap<int, edge_pair> neighbors;
    node_data* data;
    unsigned attr;

protected:
    node(const node& n, unsigned vattr = 0) : neighbors(n.neighbors), data(n.data), attr(vattr) { }
    node(node_data* d, unsigned vattr = 0) : neighbors(), data(d), attr(vattr) { }

public:
    node(unsigned vattr = 0) : neighbors(), data(nullptr), attr(vattr) { }

    // Use only to make instance of a node outside a graph!
    //static node* new_node(unsigned vattr = 0) { return new node(vattr); }

    virtual ~node() 
    { 
        if (data != nullptr) delete data; 
        clear_neighbors();
    }

    virtual void gc_clone_to(node* n, graph* g, const map<node*, node*>& cmap) const;

    iterator get_neighbors_begin(int n) { return neighbors.find(n); }
    bool neighbors_iter_end(const iterator& iter, int n) { return iter == neighbors.end() || iter->first != n; }
    
	void get_neighbors_iter(int n, iterator& begin, iterator& end);

    int count_neighbors(int n);
    node* get_neighbor(int n);
    bool has_neighbor(int n) { return get_neighbors_begin(n) != neighbors.end(); }

	edge_data* get_edge_data(int n);
    edge_data* get_edge_data(node* n, int name);
    edge_pair get_neighbor_pair(int n);

    bool is_neighbor(node* n, int name);
    void get_neighbor_set(int n, set<node*>& nset);
    void add_to_neighbor_set(int n, set<node*>& nset);
	
    void get_neighbor_pair_set(int n, set<edge_pair>& nset);

    int degree(int n);

    void clear_neighbors() { 
		// posible ERROR !!!!
		//for (iterator iter = neighbors.begin(); iter != neighbors.end(); ++iter)
        //    if (iter->second.second != nullptr) delete iter->second.second;
        neighbors.clear(); 
	}

    void delete_edges(const set<node*>& ns);
    void delete_edges(node* n);
    void delete_edges(unsigned attr);
    void delete_edges(int index);
    void delete_edges(int index, node* n);
    void delete_edges_complement(int index);
    
    void rename_edges(int from, int to);

    void new_attr(unsigned a) { attr = a; }
    void set_attr(unsigned a) { attr |= a; }
    void clear_attr(unsigned a) { attr &= ~a; }
    bool is_attr_set(unsigned a) { return (attr & a) == a; }
    bool is_attr_set2(unsigned a) { return (attr & a) != 0; }

    virtual void copy_to(streamable* p, cloner& cl);
    virtual streamable* make_instance() const { return new node(); }
    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);

    friend ostream& operator<<(ostream& os, const node& n)
    {
        os << *(n.data);
        return os;
    }

    friend class graph;

};


#endif /* _GRAPH_NODES_H_ */

