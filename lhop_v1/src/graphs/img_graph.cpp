// graph library includes definitions found in
//  - img_graph.h
////////////////////////////////////////////////

#include <algorithm>
#include <queue>
#include "../utils/atom.h"
#include "../utils/streaming.h"
#include "graph.h"
#include "img_graph.h"
#include "../utils/structures.h"
#include "../utils/utils.h"

// img_node_data 
///////////////////////////////////////////////////////////////////////////////


void img_node_data::gc_clone(const node_data* d, const map<node*,node*>& cmap)
{ 
    img_node_data* srcd = (img_node_data*)d;
    x = srcd->x;
    y = srcd->y;
    z = srcd->z;

    if (srcd->next != nullptr) {
        node* nx = srcd->next;
        map<node*,node*>::const_iterator f;

        while (nx != nullptr && (f = cmap.find(nx)) == cmap.end()) {
            nx = ((img_node_data*)nx->data)->next;
        } 
        next = (nx == nullptr || f == cmap.end()) ? nullptr : f->second;
    }
}


// img_graph
///////////////////////////////////////////////////////////////////////////////

img_graph::img_graph() : graph() { }

img_graph::img_graph(int x, int y, int z) : graph()
{
    layers[z].init_grid(x, y);
}

void img_graph::connect_neighbors(node* n, int tox, int toy, int attr, int name)
{
    if (!n->is_attr_set(IMG_NODE_ATTR)) return;

    img_node_data* nd = (img_node_data*)n->data;
    int k = nd->z;

    if (!grid(k)) init_grid(k);
    
    int maxx = min(nd->x + tox, x_size(k) - 1);
    int maxy = min(nd->y + toy, y_size(k) - 1);
    int jstart = max(nd->y - toy, 0);
    node* m;

    for (int i = max(nd->x - tox, 0); i <= maxx; ++i) 
        for (int j = jstart; j <= maxy; ++j) {
            m = node_at(i, j, k);
            if (m != nullptr && m->is_attr_set(attr)) add_edge(n, m, name);
        }
}

void img_graph::get_neighbors_circular(vector<node*>& result, node* n, int tox, int toy, unsigned attr, bool all)
{
    if (!n->is_attr_set(GRID_NODE_ATTR)) return;

    img_node_data* nd = (img_node_data*)n->data;
    int k = nd->z;

    if (!grid(k)) init_grid(k);

    int centerx = nd->x, centery = nd->y;
    int maxx = min(centerx + tox, x_size(k) - 1);
    int maxy = min(centery + toy, y_size(k) - 1);
    int jstart = max(centery - toy, 0);
    int aa = tox*tox, bb = toy*toy;
    int aabb = aa*bb;
    node* m;

    for (int i = max(centerx - tox, 0); i <= maxx; ++i) 
        for (int j = jstart; j <= maxy; ++j) {
            if ((i - centerx)*(i - centerx)*bb + (j - centery)*(j - centery)*aa <= aabb) {
                m = node_at(i, j, k);

                while (m != nullptr) {
                    if (m->is_attr_set(attr)) result.push_back(m);
                    if (!all) break;
                    m = ((img_node_data*)m->data)->next;
                }
            }
        }
}

node* img_graph::node_at(int i, int j, int k, int dx, int dy)
{
    int maxi = min<int>(i + dx + 1, x_size(k));
    int maxj = min<int>(j + dy + 1, y_size(k));
    int mini = max<int>(0, i - 1);
    int minj = max<int>(0, j - 1);

    for (int ii = mini; ii < maxi; ++ii) {
        for (int jj = minj; jj < maxj; ++jj) {
            node* n = node_at(ii, jj, k);
            if (n != nullptr) return n;
        }
    }
    return nullptr;
}

void img_graph::nodes_at(set<node*>& result, int i, int j, int k)
{
    node* n = node_at_safe(i, j, k);

    while (n != nullptr) {
        result.insert(n);
        n = ((img_node_data*)n->data)->next;
    }
}

void img_graph::nodes_at(vector<node*>& result, int i, int j, int k)
{
    node* n = node_at_safe(i, j, k);

    while (n != nullptr) {
        result.push_back(n);
        n = ((img_node_data*)n->data)->next;
    }
}

void img_graph::nodes_at(vector<node*>& result, node* n)
{
    while (n != nullptr) {
        result.push_back(n);
        n = ((img_node_data*)n->data)->next;
    }
}

void img_graph::get_neighbors_rectangular(vector<node*>& result, node* n, int tox, int toy, unsigned attr)
{
    if (!n->is_attr_set(GRID_NODE_ATTR)) return;

    img_node_data* nd = (img_node_data*)n->data;
    int k = nd->z;

    if (!grid(k)) init_grid(k);

    int centerx = nd->x, centery = nd->y;
    int maxx = min<int>(centerx + tox, x_size(k) - 1);
    int maxy = min<int>(centery + toy, y_size(k) - 1);
    int jstart = max<int>(centery - toy, 0);
    node* m;

    for (int i = max<int>(centerx - tox, 0); i <= maxx; ++i) 
        for (int j = jstart; j <= maxy; ++j) {
            m = node_at(i, j, k);

            while (m != nullptr) {
                if (m->is_attr_set(attr)) result.push_back(m);
                m = ((img_node_data*)m->data)->next;
            }
        }
        
}

void img_graph::add_neighbors_circular(set<node*>& result, node* n, int tox, int toy, unsigned attr)
{
    if (!n->is_attr_set(IMG_NODE_ATTR)) return;

    img_node_data* nd = (img_node_data*)n->data;
    int k = nd->z;

    if (!grid(k)) init_grid(k);

    int centerx = nd->x, centery = nd->y;
    int maxx = min(centerx + tox, x_size(k) - 1);
    int maxy = min(centery + toy, y_size(k) - 1);
    int jstart = max(centery - toy, 0);
    int aa = tox*tox, bb = toy*toy;
    int aabb = aa*bb;
    node* m;

    for (int i = max(centerx - tox, 0); i <= maxx; ++i) 
        for (int j = jstart; j <= maxy; ++j) {
            if ((i - centerx)*(i - centerx)*bb + (j - centery)*(j - centery)*aa <= aabb) {
                m = node_at(i, j, k);
                if (m != nullptr && m->is_attr_set(attr)) 
                    result.insert(m);
            }
        }
}

void img_graph::connect_neighbors_circular(node* n, int tox, int toy, unsigned attr, int name)
{
    if (!n->is_attr_set(IMG_NODE_ATTR)) return;

    img_node_data* nd = (img_node_data*)n->data;
    int k = nd->z;

    if (!grid(k)) init_grid(k);

    int centerx = nd->x, centery = nd->y;
    int maxx = min(centerx + tox, x_size(k) - 1);
    int maxy = min(centery + toy, y_size(k) - 1);
    int jstart = max(centery - toy, 0);
    int aa = tox*tox, bb = toy*toy;
    int aabb = aa*bb;
    node* m;

    for (int i = max(centerx - tox, 0); i <= maxx; ++i) 
        for (int j = jstart; j <= maxy; ++j) {
            if ((i - centerx)*(i - centerx)*bb + (j - centery)*(j - centery)*aa <= aabb) {
                m = node_at(i, j, k);
                if (m != nullptr && m->is_attr_set(attr)) add_edge(n, m, name);
            }
        }
}

// The same as "connect_neighbors_circular" only that it connects n to nodes 
// from another layer (ie. nd->z + zdelta)
//
void img_graph::connect_neighbors_circular_2(node* n, int tox, int toy, int zdelta, unsigned attr, int name)
{
    if (!n->is_attr_set(IMG_NODE_ATTR)) return;

    img_node_data* nd = (img_node_data*)n->data;
    int k = max(nd->z + zdelta, 0);
    double xfactor = (double)x_size(k)/x_size(nd->z);
    double yfactor = (double)y_size(k)/y_size(nd->z);

    if (!grid(k)) init_grid(k);

    tox = (int)(xfactor * tox);
    toy = (int)(yfactor * toy);

    int centerx = (int)(xfactor * nd->x), centery = (int)(yfactor * nd->y);
    int maxx = min(centerx + tox, x_size(k) - 1);
    int maxy = min(centery + toy, y_size(k) - 1);
    int jstart = max(centery - toy, 0);
    int aa = tox*tox, bb = toy*toy;
    int aabb = aa*bb;
    node* m;

    for (int i = max(centerx - tox, 0); i <= maxx; ++i) 
        for (int j = jstart; j <= maxy; ++j) {
            if ((i - centerx)*(i - centerx)*bb + (j - centery)*(j - centery)*aa <= aabb) {
                m = node_at(i, j, k);
                if (m != nullptr && m->is_attr_set(attr)) add_edge(n, m, name);
            }
        }
}

void img_graph::connect_neighbors_circular(vector<node*>& l, int tox, int toy, unsigned attr, int name)
{
    vector<node*>::iterator iter;

    for (iter = l.begin(); iter != l.end(); ++iter) {
        connect_neighbors_circular(*iter, tox, toy, attr, name);
    }
}

void img_graph::connect_neighbors_circular_2(vector<node*>& l, int tox, int toy, int zdelta, unsigned attr, int name)
{
    vector<node*>::iterator iter;

    for (iter = l.begin(); iter != l.end(); ++iter) {
        connect_neighbors_circular_2(*iter, tox, toy, zdelta, attr, name);
    }
}

// n - center node
// radius - radius to cover 
// attr - connect only with the neighbors with this attribute set
// namemap - vector of size 360; if the angle to some neighbor is alpha, 
//     then the edge gets name namemap[alpha]; if namemap[alpha] < 0, no edge is added
// one - true if only one neighbor of each name is permitted
void img_graph::connect_neighbors_circular_1(node* n, int iradius, int oradius, int attr, const vector<int>& namemap, bool one)
{
    if (!n->is_attr_set(IMG_NODE_ATTR)) return;

    img_node_data* nd = (img_node_data*)n->data;
    int k = nd->z;

    if (!grid(k)) init_grid(k);
    
    int centerx = nd->x;
    int centery = nd->y;
    double degree1 = 180.0/M_PI;
    set<int> forbiddennames;
    node* m;

    for (int rad = iradius; rad <= oradius; ++rad) {
        int maxx = min(nd->x + rad, x_size(k) - 1);
        int maxy = min(nd->y + rad, y_size(k) - 1);
        int rr = rad*rad;
        int jstart = max(nd->y - rad, 0);

        for (int i = max(nd->x - rad, 0); i <= maxx; ++i) {
            int di = i - centerx;

            for (int j = jstart; j <= maxy; ++j) {
                int dj = j - centery;
    
                if (di*di + dj*dj <= rr) {
                    m = node_at(i, j, k);
                    if (m != nullptr && m != n && m->is_attr_set2(attr)) {
                        double alpha = (dj < 0) ? 
                            (2*M_PI + atan2((double)dj, (double)di)) :
                            atan2((double)dj, (double)di);
                        int name = namemap[(int)(alpha*degree1)];
                        
                        if (name >= 0) {
                            if (!one) add_edge(n, m, name);
                            else if (forbiddennames.find(name) == forbiddennames.end()) {
                                add_edge(n, m, name);
                                forbiddennames.insert(name);
                            }
                        }
                    
                    }
                }
            }
        }
    }
}

// calls connect_neighbors_1 for each node in the list l
void img_graph::connect_neighbors_circular_1(vector<node*>& l, int iradius, int oradius, int attr, 
    const vector<int>& namemap, bool one)
{
    vector<node*>::iterator iter;

    for (iter = l.begin(); iter != l.end(); ++iter) {
        connect_neighbors_circular_1(*iter, iradius, oradius, attr, namemap, one);
    }
}

void img_graph::get_nodes_circular(vector<node*>& result, node* n, int tox, int toy, unsigned attr)
{
    //if (!n->is_attr_set(IMG_NODE_ATTR)) return;
    img_node_data* nd = (img_node_data*)n->data;
    int k = nd->z;

    if (!grid(k)) init_grid(k);

    int centerx = nd->x, centery = nd->y;
    int maxx = min(centerx + tox, x_size(k) - 1);
    int maxy = min(centery + toy, y_size(k) - 1);
    int jstart = max(centery - toy, 0);
    int aa = tox*tox, bb = toy*toy;
    int aabb = aa*bb;
    node* m;

    for (int i = max(centerx - tox, 0); i <= maxx; ++i) 
        for (int j = jstart; j <= maxy; ++j) {
            if ((i - centerx)*(i - centerx)*bb + (j - centery)*(j - centery)*aa <= aabb) {
                m = node_at(i, j, k);
                while (m != nullptr) {
                    if (m->is_attr_set(attr)) result.push_back(m);
                    m = ((img_node_data*)m->data)->next;
                }
            }
        }
}

/// initialize grid nodes from node list
void img_graph::init_grid()
{
    for (unsigned k = 0; k < MAX_LAYER_NUMBER; ++k) {
        layers[k].init_grid();
        if (!grid(k)) break;
    }
    
    img_node_data* d;
    for (list<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        if ((*iter)->is_attr_set(IMG_NODE_ATTR)) {
            d = (img_node_data*)((*iter)->data);
            r_node_at(d->x, d->y, d->z) = *iter;
        }
    }

}

void img_graph::init_grid(int z)
{
    layers[z].init_grid();
    if (!grid(z)) return;
    
    img_node_data* d;
    for (list<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        if ((*iter)->is_attr_set(IMG_NODE_ATTR)) {
            d = (img_node_data*)((*iter)->data);
            if (d->z == z) r_node_at(d->x, d->y, d->z) = *iter;
        }
    }

}

void img_graph::new_grid(int x, int y, int z)
{
    layers[z].init_grid(x, y);
}

// Make n a grid note it is assumed that n is already node in
// the graph and that n->data is of type img_node_data.
// it returns a node which was previously at the position.
node* img_graph::make_grid_node(node *n)
{
    //if (typeid(n->data) != typeid(img_node_data)) throw;

    img_node_data* nd = (img_node_data*)n->data;
    node*& pn = r_node_at(nd->x, nd->y, nd->z);
    node* result = pn;
    node* x = pn;

    pn = n;
    n->set_attr(IMG_NODE_ATTR | GRID_NODE_ATTR);
    while (x != nullptr) {
        x->clear_attr(IMG_NODE_ATTR);
        x->clear_attr(GRID_NODE_ATTR);
        x = ((img_node_data*)x->data)->next;
    }
    return result;
}

node* img_graph::add_grid_node(img_node_data* nd, int i, int j, int k)
{
    nd->x = i; nd->y = j; nd->z = k;

    node** pn = &r_node_at(i, j, k);
    node* nx = *pn;

    if (nx == nullptr) {
        *pn = add_node(nd, IMG_NODE_ATTR | GRID_NODE_ATTR);
    } else {
        bool setattr = true;

        //while (nx != nullptr && data_greater((img_node_data*)nx->data, nd)) {
        while (nx != nullptr && !data_less((img_node_data*)nx->data, nd)) {
            pn = &((img_node_data*)nx->data)->next;
            nx = *pn;
            setattr = false;
        }
        if (!setattr) *pn = add_node(nd, GRID_NODE_ATTR);
        else {
            *pn = add_node(nd, IMG_NODE_ATTR | GRID_NODE_ATTR);
            nx->clear_attr(IMG_NODE_ATTR);
        }
        ((img_node_data*)(*pn)->data)->next = nx;
    }
    return *pn;
}

void img_graph::add_grid_node(node* n)
{
    if (n == nullptr) return;

    img_node_data* nd = (img_node_data*)n->data;
    node** pn = &r_node_at(nd->x, nd->y, nd->z);
    node* nx = *pn;

    if (nx == nullptr) {
        n->set_attr(IMG_NODE_ATTR | GRID_NODE_ATTR);
        *pn = n; 
    } else {
        bool setattr = true;

        while (nx != nullptr && data_greater((img_node_data*)nx->data, nd)) {
            pn = &((img_node_data*)nx->data)->next;
            nx = *pn;
            setattr = false;
        }
        if (setattr) {
            n->set_attr(IMG_NODE_ATTR | GRID_NODE_ATTR);
            nx->clear_attr(IMG_NODE_ATTR);
        } else {
            n->set_attr(GRID_NODE_ATTR);
            n->clear_attr(IMG_NODE_ATTR);
        }
        *pn = n;
        ((img_node_data*)n->data)->next = nx;
    }
}


void img_graph::copy_to(graph* dest, const map<node*, node*>& cmap)
{
    graph::copy_to(dest, cmap);

    img_graph* destd = (img_graph*)dest;

    for (unsigned k = 0; k < MAX_LAYER_NUMBER; ++k) {
        destd->layers[k].x_size = layers[k].x_size;
        destd->layers[k].y_size = layers[k].y_size;
    }

}

void img_graph::delete_nodes()
{
    for (list<node*>::iterator i = nodes.begin(); i != nodes.end(); ++i) {
        node* n = *i;
        img_node_data* d = (img_node_data*)n->data;

        if (n->is_attr_set(IMG_NODE_ATTR) && grid(d->z) != nullptr) {
            r_node_at(d->x, d->y, d->z) = nullptr;
        }
    }
    graph::delete_nodes();
}

void img_graph::delete_nodes(const set<node*>& ns) 
{
    for (set<node*>::const_iterator i = ns.begin(); i != ns.end(); ++i) {
        node* n = *i;
        img_node_data* d = (img_node_data*)n->data;

        if (grid(d->z) != nullptr) {
            if (n->is_attr_set(IMG_NODE_ATTR)) {
                r_node_at(d->x, d->y, d->z) = d->next;
                if (d->next != nullptr) d->next->set_attr(IMG_NODE_ATTR);
            } else if (n->is_attr_set(GRID_NODE_ATTR)) {
                node* m = r_node_at(d->x, d->y, d->z);
                img_node_data* md = (img_node_data*)m->data;

                while (m != nullptr && md->next != n) {
                    m = md->next;
                    md = (img_node_data*)m->data;
                }
                if (m != nullptr) {
                    md->next = d->next;
                }
            }
        }
    }
    graph::delete_nodes(ns);
}

void img_graph::delete_node(node* n)
{
    if (n->is_attr_set(IMG_NODE_ATTR)) {
        img_node_data* d = (img_node_data*)n->data;
        r_node_at(d->x, d->y, d->z) = nullptr;
    }
    graph::delete_node(n);
}

// Delete edges between nodes od grids <= z; edges, for example, between "layers"
// z + 1 and 0 are not removed
void img_graph::delete_edges_leq(int index, int z)
{
    for (auto iter = nodes.begin(); iter != nodes.end(); ++iter)
        if (z_coo(*iter) <= z)
            (*iter)->delete_edges([index, z](int name, node* n, edge_data*) { return name == index && z_coo(n) <= z; });

}

void img_graph::change_data(node* n, img_node_data* newd, bool dispose)
{
    img_node_data* oldnd = (img_node_data*)n->data;

    if (n->is_attr_set(GRID_NODE_ATTR)) 
        newd->next = ((img_node_data*)n->data)->next;
    newd->x = oldnd->x;
    newd->y = oldnd->y;
    newd->z = oldnd->z;
    graph::change_data(n, newd, dispose);
}

void img_graph::remove_grid_node(const pair<node*, node**>& dp)
{
    if (dp.first != nullptr) {
        node* newnx = ((img_node_data*)dp.first->data)->next;

        if (newnx != nullptr && dp.first->is_attr_set(IMG_NODE_ATTR)) newnx->set_attr(IMG_NODE_ATTR);
        *dp.second = newnx;
        dp.first->new_attr(0);
    }
}

pair<node*,node**> img_graph::find_at(img_node_data* d, int i, int j, int k /* = 0 */)
{
    node** pn = &r_node_at(i, j, k);
    node* n = *pn;

    while (n != nullptr) {
        if (data_equal((img_node_data*)n->data, d)) return pair<node*, node**>(n, pn);
        pn = &((img_node_data*)n->data)->next;
        n = *pn;
    }
    return pair<node*, node**>(n, pn);
}

bool img_graph::img_node_greater(img_graph* g, node* n1, node* n2)
{
    return g->data_greater((img_node_data*)n1->data, (img_node_data*)n2->data);
}

int img_graph::count_at(int i, int j, int k /* = 0 */)
{
    int result = 0;
    node* n = node_at(i, j, k);

    while (n != nullptr) {
        ++result;
        n = ((img_node_data*)n->data)->next;
    }
    return result;
}

void img_graph::sort_at(int i, int j, int k /* = 0 */) 
{
    vector<node*> nlist;
    img_node_comparer comp(this);
    node** np = &r_node_at(i, j, k);
    node* n = *np;

    while (n != nullptr) {
        nlist.push_back(n);
        n = ((img_node_data*)n->data)->next;
    }
    sort(nlist.begin(), nlist.end(), comp);

    for (vector<node*>::iterator iter = nlist.begin(); iter != nlist.end(); ++iter) {
        node* n = *iter;

        n->clear_attr(IMG_NODE_ATTR);
        *np = n;
        np = &((img_node_data*)n->data)->next;
    }
    *np = nullptr;
    if (!nlist.empty()) 
        nlist.front()->set_attr(IMG_NODE_ATTR);
}

void img_graph::sort_at(node* n)
{
    img_node_data* nd = (img_node_data*)n->data;
    sort_at(nd->x, nd->y, nd->z);
}

void img_graph::all_grid_nodes(vector<node*>& result, const vector<node*> src)
{
    result.clear();
    for (vector<node*>::const_iterator iter = src.begin(); iter != src.end(); ++iter) {
        node* n = *iter;

        while (n != nullptr) {
            result.push_back(n);
            n = ((img_node_data*)n->data)->next;
        }
    }
}

void img_graph::connect_region(node* n, matrix<bool>& r, int name)
{
    if (!n->is_attr_set(IMG_NODE_ATTR)) return;

    img_node_data* nd = (img_node_data*)n->data;
    int rwidth = (int)r.width, rheight = (int)r.height;
    int x_start, y_start, x_end, y_end;
    int ix_start, iy_start;
    int z = nd->z;
    node* m;

    x_start = max(nd->x - rwidth/2, 0); 
    x_end = min(x_start + rwidth, x_size(z));
    y_start = max(nd->y - rheight/2, 0); 
    y_end = min(y_start + rheight, y_size(z));
    ix_start = x_start - (nd->x - rwidth/2);
    iy_start = y_start - (nd->y - rheight/2);

    for (int i = ix_start, x = x_start; x < x_end; ++i, ++x) {
        for (int j = iy_start, y = y_start; y < y_end; ++j, ++y) {
            if (r(i, j)) {
                m = node_at(x, y, z);
                if (m != n) {
                    if (m == nullptr) m = add_grid_node(new img_node_data(), x, y, z);
                    add_edge(n, m, name);
                }
            }
        }
    }
}

void img_graph::get_region(node* n, matrix<bool>& r, set<node**>& result)
{
    if (!n->is_attr_set(IMG_NODE_ATTR)) return;

    img_node_data* nd = (img_node_data*)n->data;
    int rwidth = (int)r.width, rheight = (int)r.height;
    int x_start, y_start, x_end, y_end;
    int ix_start, iy_start;
    int z = nd->z;
    node** m;

    x_start = max(nd->x - rwidth/2, 0); 
    x_end = min(nd->x + rwidth/2 + 1, x_size(z));
    y_start = max(nd->y - rheight/2, 0); 
    y_end = min(nd->y + rheight/2 + 1, y_size(z));
    ix_start = x_start - (nd->x - rwidth/2);
    iy_start = y_start - (nd->y - rheight/2);

    result.clear();
    for (int i = ix_start, x = x_start; x < x_end; ++i, ++x) {
        for (int j = iy_start, y = y_start; y < y_end; ++j, ++y) {
            if (r(i, j)) {
                m = &r_node_at(x, y, z);
                if (*m != n) result.insert(m);
            }
        }
    }
}

int img_graph::get_neighbor_set_extents2(node* n, int name)
{
    rectangle2<int> box;
        
    foreach_neighbor(n, name, iter) {
        img_node_data* d = (img_node_data*)neighbor_node_data(iter);
        box.eat(d->x, d->y);
    }
    return box.size2();
}

int img_graph::get_nodes_in_rect(vector<node*>& result, int z, const irectangle2& rect)
{
    if (!grid(z)) init_grid(z);

    int minj = ::max<int>(0, rect.ll.y);
    int maxi = ::min<int>(rect.ur.x, x_size(z));
    int maxj = ::min<int>(rect.ur.y, y_size(z));
    int count = 0;
    
    for (int i = max(0, rect.ll.x); i < maxi; ++i) {
        for (int j = minj; j < maxj; ++j) {
            node* n = node_at(i, j, z);

            while (n != nullptr) {
                result.push_back(n);
                n = ((img_node_data*)n->data)->next;
                ++count;
            }
        }
    }
    return count;
}

vector<node*> img_graph::get_grid_nodes_circular(int x, int y, int z, int dx, int dy)
{
    if (!grid(z)) init_grid(z);

    int imin = max(x - dx, 0);
    int imax = min(x + dx + 1, x_size(z));
    int jmin = max(y - dy, 0);
    int jmax = max(y + dy + 1, y_size(z));
    vector<node*> result;

    for (int i = imin; i < imax; ++i)
        for (int j = jmin; j < jmax; ++j) {
            node* n = node_at(i, j, z);

            if (n != nullptr) {
                int d2 = (i - x)*(i - x) + (j - y)*(j - y);

                if (d2 < dx*dy) 
                    result.push_back(n);
            }
        }
    return result;
}

node* img_graph::get_closest_grid_node(int x, int y, int z, int dx, int dy)
{
    if (!grid(z)) init_grid(z);

    int imin = max(x - dx, 0);
    int imax = min(x + dx + 1, x_size(z));
    int jmin = max(y - dy, 0);
    int jmax = max(y + dy + 1, y_size(z));
    node* result = nullptr;
    int mindist2 = INT_MAX;


    for (int i = imin; i < imax; ++i)
        for (int j = jmin; j < jmax; ++j) {
            node* n = node_at(i, j, z);

            if (n != nullptr) {
                int d2 = (i - x)*(i - x) + (j - y)*(j - y);

                if (d2 < mindist2) {
                    result = n;
                    mindist2 = d2;
                }
            }
        }
    return result;
}


void img_graph::write_vgr_label(ostream& os, node* n, int count)
{
    os << '\"' << count << '\"';
}

// Writes graph in "Pajek" format
// os - output stream
// nodes - list of nodes to export
// edgename - only export edges with edgename; call edgename = -1 exports all edges
void img_graph::write_vgr(ostream& os, const vector<node*>& nodes, int edgename)
{
    vector<node*>::const_iterator iter;
    map<node*, int> enumeration;
    int count;
    int minx, maxx, miny, maxy;
    node* n;
    img_node_data* nd;

    if (nodes.empty()) return;

    // make enumeration & determine maxx, minx, miny, maxy

    nd = (img_node_data*)(nodes.front())->data;   
    maxx = minx = nd->x;
    maxy = miny = nd->y;
    count = 0;
    for (iter = nodes.begin(); iter != nodes.end(); ++iter) {
        n = *iter;
        nd = (img_node_data*)n->data;   

        enumeration.insert(pair<node*, int>(n, count++));
        if (nd->x > maxx) maxx = nd->x; else if (nd->x < minx) minx = nd->x;
        if (nd->y > maxy) maxy = nd->y; else if (nd->y < miny) miny = nd->y;
    }

    int xsize = maxx - minx;
    int ysize = maxy - miny;

    // write vertices
    os << "*Vertices " << nodes.size() << endl;
    count = 0;
    for (iter = nodes.begin(); iter != nodes.end(); ++iter) {
        n = *iter;
        nd = (img_node_data*)n->data;   

        ++count;
        os << count << ' ';
        write_vgr_label(os, n, count);
        os << ' ';
        os << (double)(nd->x - minx)/xsize << ' ' << (double)(nd->y - miny)/ysize << ' ' << 0.5 << endl;
    }

    // write arcs
    // os << "*Arcs" << endl;

    // write edges
    os << "*Edges" << endl;
    count = 0;
    for (iter = nodes.begin(); iter != nodes.end(); ++iter) {
        n = *iter;
        nd = (img_node_data*)n->data;   

        //for (multimap<int, edge_pair>::iterator _iter = (n)->neighbors.begin();
        //    _iter != (n)->neighbors.end(); ++_iter) {
        forall_neighbors(n, i) {
            if (edgename < 0 || i->first == edgename) {
                node* m = i->second.first;
                map<node*, int>::iterator pos = enumeration.find(m);

                if (pos != enumeration.end() && pos->second > count) 
                    os << count + 1 << ' ' << pos->second + 1 << ' ' << i->first << endl;
            }
        }
        ++count;
    }
}


void img_graph::copy_to(streamable* p, cloner& cl)
{
    graph::copy_to(p, cl);

    for (unsigned k = 0; k < MAX_LAYER_NUMBER; ++k) {
        ((img_graph*)p)->layers[k].x_size = layers[k].x_size;
        ((img_graph*)p)->layers[k].y_size = layers[k].y_size;
    }
}

void img_graph::read_from_stream(istreamer& is)
{
    graph::read_from_stream(is);

    unsigned maxk = (is.get_version() == 1.0) ? 8 : MAX_LAYER_NUMBER;

    for (unsigned k = 0; k < maxk; ++k) {
        is.read(layers[k].x_size);
        is.read(layers[k].y_size);
    }
}

void img_graph::write_to_stream(ostreamer& os)
{
    graph::write_to_stream(os);

    for (unsigned k = 0; k < MAX_LAYER_NUMBER; ++k) {
        os.write(x_size(k));
        os.write(y_size(k));
    }
}

// earth
///////////////////////////////////////////////////////////////////////////////

earth::earth(const matrix<double>& m, int cir /* = -1 */) : 
    img_graph((int)m.width, (int)m.height)
{
    for (int i = 0; i < (int)m.width; ++i) {
        for (int j = 0; j < (int)m.height; ++j) {
            add_grid_node(new img_node_data_t<double>(m(i, j)), i, j);
        }
    }

    if (cir < 0) return;

    int cir2 = cir*cir;
    int cx = x_size()/2, cy = y_size()/2;
    node* c = node_at(cx, cy);

    if (c == nullptr) return;

    c->set_attr(EARTH_CENTER_ATTR);

    int minx = max(0, cx - cir);
    int maxx = min(x_size() - 1, cx + cir);
    int miny = max(0, cy - cir);
    int maxy = min(y_size() - 1, cy + cir);
    int edgename = atom("island_center").get_index();

    for (int i = minx; i <= maxx; ++i) {
        for (int j = miny; j <= maxy; ++j) {
            node* n = node_at(i, j);

            if (n != nullptr && distance2(cx, cy, i, j) <= cir2) { 
                n->set_attr(ISLAND_ATTR);
                add_edge(n, c, edgename);
            }
        }
    }
}

bool earth::on_island(int sx, int sy)
{
    node* n = node_at_safe(sx, sy);
    if (n == nullptr || !n->is_attr_set(ISLAND_ATTR)) return false;
    else return true;
}

node* earth::make_island(int sx, int sy, double thresh, int min_radius /* = 0 */)
{
    node* start = node_at_safe(sx, sy);

    if (start == nullptr || t_img_data<double>(start) < thresh) return nullptr;

    int edgename = atom("island_center").get_index();

    if (start->is_attr_set(ISLAND_ATTR)) return start->get_neighbor(edgename);
    
    queue<node*> q;
    int width = x_size(), height = y_size();
    int x, y, newx, newy;
    int min_radius2 = min_radius*min_radius;
    
    q.push(start);
    while (!q.empty()) {
        node* n = q.front();
        
        x = x_coo(n); y = y_coo(n);
        if (n != nullptr && !n->is_attr_set(ISLAND_ATTR) && 
                (t_img_data<double>(n) >= thresh || distance2(sx, sy, x, y) <= min_radius2)) {
            n->set_attr(ISLAND_ATTR);
            add_edge(n, start, edgename);
            if ((newx = x - 1) >= 0) {
                q.push(node_at(newx, y));
                if ((newy = y - 1) >= 0) q.push(node_at(newx, newy));
                if ((newy = y + 1) < height) q.push(node_at(newx, newy));
            }
            if ((newx = x + 1) < width) {
                q.push(node_at(newx, y));
                if ((newy = y - 1) >= 0) q.push(node_at(newx, newy));
                if ((newy = y + 1) < height) q.push(node_at(newx, newy));
            }
            if ((newy = y - 1) >= 0) q.push(node_at(x, newy));
            if ((newy = y + 1) < height) q.push(node_at(x, newy));

        }
        q.pop();
    }
    return start;
}

node* earth::get_center(int sx, int sy)
{
    node* start = node_at_safe(sx, sy);

    if (start != nullptr) return start->get_neighbor(atom("island_center")); 
    else return nullptr;
}

void earth::print(ostream& os)
{
    int width = x_size();
    int height = y_size();

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            if (node_at(i, j)->is_attr_set(ISLAND_ATTR)) os << "1 "; else os << "0 ";
        }
        os << endl;
    }
}   


