// img graph class

#pragma once
#ifndef _IMG_GRAPH_H_
#define _IMG_GRAPH_H_

#define _USE_MATH_DEFINES  // necessary for <math.h> to define M_PI,.... constants, etc. 
#include <math.h>
#include <string.h>
#include "graph.h"
#include "../utils/matrix.h"

// constants
///////////////////////////////////////////////////////////////////////////////

const unsigned IMG_NODE_ATTR = 2;    // First in the grid
const unsigned GRID_NODE_ATTR = 8;   // All in the grid

const unsigned MAX_LAYER_NUMBER = 16;

// img_node_data
//////////////////

struct img_node_data : public node_data {
    int x, y, z;
    node* next;

    img_node_data(int vx = 0, int vy = 0, int vz = 0) : x(vx), y(vy), z(vz), next(nullptr) { }
    img_node_data(const img_node_data& nd) : x(nd.x), y(nd.y), z(nd.z), next(nd.next) { } 

    virtual void clone(const node_data* d) 
    { 
        img_node_data* ind = (img_node_data*)d;
        x = ind->x; y = ind->y; z = ind->z;
    }

    virtual void gc_clone(const node_data* d, const map<node*,node*>& cmap);

    int distance2(const img_node_data& d) const
    {
        int dx = x - d.x, dy = y - d.y;

        return dx*dx + dy*dy;
    }

    virtual void copy_to(streamable* p, cloner& cl) 
    {
        node_data::copy_to(p, cl);
        ((img_node_data*)p)->x = x;
        ((img_node_data*)p)->y = y;
        ((img_node_data*)p)->z = z;
        ((img_node_data*)p)->next = (node*)cl.get_copy(next);
    }

    virtual streamable* make_instance() const { return new img_node_data(); }
    virtual void read_from_stream(istreamer& is) 
    { 
        node_data::read_from_stream(is);    
        is.read(x); is.read(y); is.read(z);
        is.read((streamable*&)next);
    }

    virtual void write_to_stream(ostreamer& os)
    { 
        node_data::write_to_stream(os);
        os.write(x); os.write(y); os.write(z);
        os.write((streamable*)next);
    }

    virtual void print(ostream& os) const
    {
        os << '(' << x << ',' << x << ')';
    }

};

// img_node_data_t
////////////////////

template<class T> class img_node_data_t : public img_node_data {
public:
    T data;

    img_node_data_t() : img_node_data(0, 0, 0) { }
    img_node_data_t(const T& d) : img_node_data(0, 0, 0), data(d) { }


    virtual void copy_to(streamable* p, cloner& cl)
    {
        img_node_data::copy_to(p, cl);
        ((img_node_data_t*)p)->data = data;
    }

    virtual streamable* make_instance() const { return new img_node_data_t<T>(); }
    virtual void read_from_stream(istreamer& is)  
    { 
        img_node_data::read_from_stream(is); 
        is.read(data); 
    }
    
    virtual void write_to_stream(ostreamer& os) 
    { 
        img_node_data::write_to_stream(os); 
        os.write(data);
    }
    
    virtual void print(ostream& os) const
    {
        img_node_data::print(os);
        os << ' ' << data;
    }

};

// img_graph_layer
////////////////////

struct img_graph_layer {
    int x_size, y_size;
    node** grid;

    img_graph_layer() : x_size(0), y_size(0), grid(0) { }
    img_graph_layer(int x, int y) : x_size(0), y_size(0), grid(0) { init_grid(x, y); }
    ~img_graph_layer() { if (grid) delete[] grid; }
    
    node*& r_node_at(int i, int j) { return grid[x_size*j + i]; }
    node* operator()(int i, int j) { return grid[x_size*j + i]; }
    node* node_at(int i, int j) { return (grid) ? grid[x_size*j + i] : 0; }
    node* node_at_safe(int i, int j) 
        { return (!grid || i < 0 || i >= x_size || j < 0 || j >= y_size) ? 0 : grid[x_size*j + i]; }

    void init_grid(int x, int y) 
    { 
        if (grid && x_size*y_size == x*y) {
            x_size = x;
            y_size = y;
        } else {
            if (grid) delete[] grid;
            x_size = x;
            y_size = y;
            grid = new node*[x_size*y_size];
        }
        memset(grid, 0, x_size*y_size*sizeof(node*));
    }

    void init_grid() 
    {
        if (grid || x_size == 0 || y_size == 0) return;
        grid = new node*[x_size*y_size];
        memset(grid, 0, x_size*y_size*sizeof(node*));
    }


};

struct layer_predicate {
    int lyr;

    layer_predicate(int l = 0) : lyr(l) { }

    bool operator()(const node*, const node* n) const { return ((img_node_data*)(n->data))->z >= lyr; }
};

// img_graph
//////////////

class img_graph : public graph {
public:
    img_graph_layer layers[MAX_LAYER_NUMBER];

    img_graph();
    img_graph(int x_size, int y_size, int z = 0);
    img_graph(istreamer&);

protected:
    node*& r_node_at(int i, int j, int k = 0) { return layers[k].r_node_at(i, j); }
    node*& r_next_node(node* n) { return ((img_node_data*)n->data)->next; }

public:

    node* operator()(int i, int j, int k = 0) { return layers[k](i, j); }
    node* node_at(int i, int j, int k = 0) { return layers[k].node_at(i, j); }
    node* node_at(int i, int j, int k, int dx, int dy);
    node* node_at_safe(int i, int j, int k = 0) 
        { return layers[k].node_at_safe(i, j); }
    img_node_data* node_data(node* n) { return (img_node_data*)n->data; }
    void nodes_at(set<node*>& result, int i, int j, int k);
    void nodes_at(vector<node*>& result, int i, int j, int k);
    void nodes_at(vector<node*>& result, node* n);

    node** grid(int z = 0) { return layers[z].grid; }
    int x_size(int z = 0) { return layers[z].x_size; }
    int y_size(int z = 0) { return layers[z].y_size; }
    bool empty(int z = 0) { return x_size(z) == 0 || y_size(z) == 0; }
    void init_grid();
    void init_grid(int z);
    void new_grid(int x, int y, int z); 

    void change_data(node* n, img_node_data* newd, bool dispose);

    virtual void delete_nodes();
    virtual void delete_nodes(const set<node*>& ns);
    virtual void delete_node(node* n);

    void delete_edges_leq(int index, int l);
    node* add_grid_node(img_node_data* nd, int i, int j, int k = 0);
    void add_grid_node(node* nd);
    node* make_grid_node(node* n);
    void remove_grid_node(const pair<node*, node**>& dp);
    template<class P> void remove_grid_nodes(int i, int j, int k, P pred);
    virtual bool data_less(img_node_data*, img_node_data*) { return false; }
    virtual bool data_greater(img_node_data*, img_node_data*) { return false; }
    virtual bool data_equal(img_node_data*, img_node_data*) { return false; }
    pair<node*, node**> find_at(img_node_data* nd, int i, int j, int k = 0);
    void sort_at(int i, int j, int k = 0);
    void sort_at(node* n);
    static bool img_node_greater(img_graph* g, node* n1, node* n2);
    void all_grid_nodes(vector<node*>& result, const vector<node*> src);
    int count_at(int i, int j, int k = 0);

    void connect_neighbors(node* n, int tox, int toy, int attr, int name);
    void get_neighbors_circular(vector<node*>& result, node* n, int tox, int toy, unsigned attr, bool all);
    void get_neighbors_rectangular(vector<node*>& result, node* n, int tox, int toy, unsigned attr);
    void add_neighbors_circular(set<node*>& result, node* n, int tox, int toy, unsigned attr);
    void connect_neighbors_circular(node* n, int tox, int toy, unsigned attr, int name);
    void connect_neighbors_circular(vector<node*>& l, int tox, int toy, unsigned attr, int name);
    void connect_neighbors_circular_1(node* n, int iradius, int oradius, int attr, 
        const vector<int>& namemap, bool one);
    void connect_neighbors_circular_1(vector<node*>& l, int iradius, int oradius, int attr, 
        const vector<int>& namemap, bool one);
    void connect_neighbors_circular_2(node* n, int tox, int toy, int zdelta, unsigned attr, int name);
    void connect_neighbors_circular_2(vector<node*>& l, int tox, int toy, int zdelta, unsigned attr, int name);
    void connect_region(node* n, matrix<bool>& r, int name);
    void get_region(node* n, matrix<bool>& r, set<node**>& result);
    int get_neighbor_set_extents2(node* n, int name);
    double get_neighbor_set_extents(node* n, int name) { return sqrt((double)get_neighbor_set_extents2(n, name)); }
    int get_nodes_in_rect(vector<node*>& result, int z, const irectangle2& rect);
    void get_nodes_circular(vector<node*>& result, node* n, int tox, int toy, unsigned attr);
    node* get_closest_grid_node(int x, int y, int z, int dx, int dy);
    vector<node*> get_grid_nodes_circular(int x, int y, int z, int dx, int dy);

    virtual void write_vgr_label(ostream& os, node* n, int count);
    void write_vgr(ostream& os, const vector<node*>& nodes, int edgename);
    void write_mma2(ostream& os, const vector<node*>& nodes, int edgename);

    virtual void copy_to(streamable* p, cloner& cl);
    virtual streamable* make_instance() const { return new img_graph(); }
    virtual void read_from_stream(istreamer& is);
    virtual void write_to_stream(ostreamer& os);

protected:
    void copy_to(graph* dest, const map<node*, node*>& cmap);

};

class img_node_comparer {
protected:
    img_graph* g;
public:
    img_node_comparer(img_graph* pg) : g(pg) { }

    bool operator()(node* n1, node* n2) 
    { 
        return g->data_greater((img_node_data*)n1->data, (img_node_data*)n2->data); 
    }
};

// earth
//////////

class earth : public img_graph {
public:
    earth(const matrix<double>& m, int cir = -1);

    node* make_island(int sx, int sy, double thresh, int min_radius = 0);
    node* get_center(int sx, int sy);
    bool on_island(int sx, int sy);
    void print(ostream& os);
};

const unsigned ISLAND_ATTR = 16;
const unsigned EARTH_CENTER_ATTR = 32;


// inline & template definitions
///////////////////////////////////////////////////////////////////////////////

inline int x_coo(node* n) { return ((img_node_data*)n->data)->x; }
inline int y_coo(node* n) { return ((img_node_data*)n->data)->y; }
inline int z_coo(node* n) { return ((img_node_data*)n->data)->z; }
template<class T> T t_img_data(node* n) { return ((img_node_data_t<T>*)n->data)->data; }
inline img_node_data* cast_img_data(node* n) { return (img_node_data*)n->data; }

template<class Container> int img_node_set_distance2(const Container& c1, const Container& c2)
{
    if (c1.empty() || c1.empty()) return -1;

    int result = cast_img_data(*c1.begin())->distance2(*cast_img_data(*c2.begin()));

    for (typename Container::const_iterator iter1 = c1.begin(); iter1 != c1.end(); ++iter1) {
        img_node_data* d1 = (img_node_data*)(*iter1)->data;

        for (typename Container::const_iterator iter2 = c2.begin(); iter2 != c2.end(); ++iter2) {
            int dist2 = d1->distance2(*((img_node_data*)(*iter2)->data));
            
            if (dist2 < result) result = dist2;
        }
    }
    return result;
}

template<class P> void img_graph::remove_grid_nodes(int i, int j, int k, P predicate)
{
    node** np = &r_node_at(i, j, k);
    node* n = *np;

    while (n != nullptr) {
        if (predicate(n)) {
            n->set_attr(NODE_DELETED_ATTR);
            *np = ((img_node_data*)n->data)->next;
        } else {
            np = &((img_node_data*)n->data)->next;
        }
        n = *np;
    }
    if ((n = node_at(i, j, k)) != nullptr) n->set_attr(IMG_NODE_ATTR);
}


#endif
