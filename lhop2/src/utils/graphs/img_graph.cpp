// graph library includes definitions found in
//  - img_graph.h
////////////////////////////////////////////////

#include <algorithm>
#include <queue>
//#include "utils/atom.h"
#include "utils/serialization/streaming.h"
#include "utils/graphs/graph.h"
#include "img_graph.h"
#include "utils/structures.h"
#include "utils/utils.h"

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

void img_graph::nodes_at(vector<node*>& result, node* n)
{
    while (n != nullptr) {
        result.push_back(n);
        n = ((img_node_data*)n->data)->next;
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
// unused
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

// unused 
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

int img_graph::get_neighbor_set_extents2(node* n, int name)
{
    rectangle2<int> box;
        
    foreach_neighbor(n, name, iter) {
        img_node_data* d = (img_node_data*)neighbor_node_data(iter);
        box.eat(d->x, d->y);
    }
    return box.size2();
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
