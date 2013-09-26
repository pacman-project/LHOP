/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// layer_1 

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>

#if defined WIN32 | defined WIN64
#include <crtdbg.h>
#endif

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stack>
#include <math.h>

#include "../utils/utils.h"
#include "../graphs/graph_utils.h"
#include "part_lib.h"
#include "layer_1_result.h"
#include "layer_1_creators.h"
#include <highgui.h>

//#define DEBUG_OUTPUT

using namespace std;

const double response_map::empty_val = numeric_limits<double>::quiet_NaN();

bool debug_output = false;

void set_debug_output(bool val) { debug_output = val; }
bool get_debug_output() { return debug_output; }

int response_from_string(string str)
{
    typedef map<string, int> rmap_t;

    static rmap_t rmap;

    if (rmap.empty()) {
        rmap.insert(rmap_t::value_type("R_RESPONSE", 0));
        rmap.insert(rmap_t::value_type("G_RESPONSE", 2));
        rmap.insert(rmap_t::value_type("S_RESPONSE", 3));
        rmap.insert(rmap_t::value_type("RR_RESPONSE", 5));
        rmap.insert(rmap_t::value_type("X_RESPONSE", 6));
    }
    ::transform(str.begin(), str.end(), str.begin(), ::toupper);

    rmap_t::iterator miter = rmap.find(str);

    if (miter == rmap.end()) return -1;
    else return miter->second;
}

// -1 for "better-is-smaller" responses and +1 for "better-is-bigger" responses
int response_type(int r)
{
    if (r == S_RESPONSE) return -1; else return 1;
}

double reverse_s_response(double s, double maxs /* = 5.0 */)
{
    return 1 - min(s, maxs)/maxs;
}

void response_map::write_to_stream(ostreamer& os)
{ 
    os.write(size);

    int i = 0, j = 0;

    while (i < size) {
        if (!_isnan(v[j])) { 
            os.write(j); os.write(v[j]); 
            ++i;
        }
        ++j;
    }
}

void response_map::read_from_stream(istreamer& is)
{
    is.read(size);
    for (int i = 0; i < size; ++i) {
        int k;
        double val;

        is.read(k); is.read(val);
        v[k] = val;
    }
}



// layer1_result
///////////////////////////////////////////////////////////////////////////////

// static fields
//////////////////

r_response_spf layer1_result::rspf;
v_response_spf layer1_result::vspf;
simple_g_response_spf layer1_result::sgspf;
identity_spf layer1_result::idspf;


layer1_result::~layer1_result() {
#ifdef OPENCL
	for (int i = 0; i < ocl_shape_nodes.size(); i++) {
		if (ocl_shape_nodes[i].first != nullptr) delete[] ocl_shape_nodes[i].first;
	}
	for (int i = 0; i < ocl_edges.size(); i++)  {
		if (ocl_edges[i].first != nullptr) delete[] ocl_edges[i].first;
	}
	for (int i = 0; i < ocl_shape_nodes_coord.size(); i++)  {
		if (ocl_shape_nodes_coord[i].first != nullptr) delete[] ocl_shape_nodes_coord[i].first;
	}
	for (int i = 0; i < ocl_shape_nodes_inhib_coord.size(); i++)  {
		if (ocl_shape_nodes_inhib_coord[i].first != nullptr) delete[] ocl_shape_nodes_inhib_coord[i].first;
	}
#endif
}

void layer1_result::inhibit(int z)
{
    int w = x_size(z);
    int h = y_size(z);
    layer1_data* d;
    iimg forb(w, h, 0);
    int i, j;
    vector<node*>& s_nodes = shape_nodes[z];
    vector<node*>& s_nodes_inhib = shape_nodes_inhib[z];
    
    sort(s_nodes.begin(), s_nodes.end(), layer1_data::greater1n); 
    s_nodes_inhib.clear();

    vector<node*>::iterator iter, enditer;

    for (iter = s_nodes.begin(); iter != s_nodes.end(); ++iter)  {
        d = (layer1_data*)((*iter)->data);
        i = d->x; j = d->y;
        if (forb(i, j) == 0) {
            s_nodes_inhib.push_back(*iter);
            (*iter)->set_attr(NODE_REDUCED_ATTR);
            forb(i - 1, j - 1) = 1;
            forb(i - 1, j) = 1;
            forb(i - 1, j + 1) = 1;
            forb(i, j - 1) = 1;
            forb(i, j) = 1;
            forb(i, j + 1) = 1;
            forb(i + 1, j - 1) = 1;
            forb(i + 1, j) = 1;
            forb(i + 1, j + 1) = 1;
        }
    }

}

// "Inhibits" boxes using a greedy method
// ***Is there a better (fast, correct) way of doing this!?*** NP-hard?
// 'boxes' is a *sorted* list of boxes (pair.second) and their values (pair.first)
// The result is a subset of 'boxes' s.t. |b_i \cap b_j|/|b_i \cup b_j| < t for all i != j
void inhibit_boxes(list<layer1_result::box_data_t>& boxes, double t)
{
	typedef layer1_result::box_data_t list_item_t;
	typedef list<list_item_t> list_t;

	list_t result;

	for (list_t::const_iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
		double maxt = 0.0;
		const irectangle2& bbox = biter->box;

		for (list_t::iterator riter = result.begin(); riter != result.end(); ++riter) {
			irectangle2& rbox = riter->box;
            double r = (double)(rbox.intersection(bbox).area())/min<int>(bbox.area(), rbox.area());
            //double r = (double)(rbox.intersection(bbox).area())/rbox.union_area(bbox);

			if (r > maxt) maxt = r;
		}
		if (maxt <= t)
			result.push_back(*biter);
	}
	boxes = result;
}

// Checks for hits and misses of boxes, returns pointers to "hit" boxes
void check_hits(vector<const layer1_result::box_data_t*>& hitboxes, vector<const layer1_result::box_data_t*>& missboxes,
                vector<bool>& hits, int& misses, const list<pair<irectangle2, int> >& gtrs,
                const list<layer1_result::box_data_t>& boxes, double thresh)
{
    hitboxes.clear();
    missboxes.clear();

    misses = 0;
    hits.resize(gtrs.size(), false);
    for (list<layer1_result::box_data_t>::const_iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
        const irectangle2& r = biter->box;
        double bestf = 0.0;
        int rcount = 0;

        for (list<pair<irectangle2, int> >::const_iterator iter = gtrs.begin(); iter != gtrs.end(); ++iter) {
            const pair<irectangle2, int>& gtr = *iter;

            if (biter->m == gtr.second) {
                double f = (double)(r.intersection(gtr.first)).area()/r.union_area(gtr.first);

                if (f > bestf) bestf = f;
                if (f >= thresh) {
                    hits[rcount] = true;
                    hitboxes.push_back(&(*biter));
                }
                ++rcount;
            }
        }
        if (bestf < thresh) {
            ++misses;
            missboxes.push_back(&(*biter));
        }
    }
}

//
//void layer1_result::count_hits(vector<bool>& hits, int& misses, vector<pair<double, bool> >& confidence,
//     const list<pair<irectangle2, int> >& gtrs, const list<box_data_t>& boxes, double thresh)
//{
//    misses = 0;
//    hits.resize(gtrs.size(), false);
//    for (list<box_data_t>::const_iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
//        const irectangle2& r = biter->box;
//        double bestf = 0.0;
//        int count = 0;
//        
//        for (list<pair<irectangle2, int> >::const_iterator iter = gtrs.begin(); iter != gtrs.end(); ++iter) {
//            const pair<irectangle2, int>& gtr = *iter;
//
//			if (biter->m == gtr.second) {
//				double f = (double)(r.intersection(gtr.first)).area()/r.union_area(gtr.first);
//
//	            if (f > bestf) bestf = f;
//		        if (f >= thresh) hits[count] = true;
//			}
//            ++count;
//        }
//        // we did not hit anything?
//        if (bestf < thresh) ++misses;
//        confidence.push_back(pair<double, bool>(biter->get_value(), bestf >= thresh));
//    }
//}

                
                
                //void layer1_result::

void layer1_result::add_results(layer1_result* res1, layer1_result* res2, int z)
{
    vector<layer1_result*> results;

    results.push_back(res1);
    results.push_back(res2);
    add_results(results, z);
}


void layer1_result::add_results(const vector<layer1_result*>& results, int z)
{
    vector<node*>& s_nodes = shape_nodes[z];

    for (vector<layer1_result*>::const_iterator iter = results.begin(); iter != results.end(); ++iter) {
        vector<node*>& s_nodes_src = (*iter)->shape_nodes[z];
        vector<node*>::const_iterator riter;
        node* n;

        for (riter = s_nodes_src.begin(); riter != s_nodes_src.end(); ++riter) {
            layer1_data* srcd = (layer1_data*)((*riter)->data);

            if ((n = node_at(srcd->x, srcd->y, z)) == nullptr) {
                node* nn = add_grid_node(new layer1_data(srcd->r, srcd->m), srcd->x, srcd->y, z);
				s_nodes.push_back(nn);
            } else {
                layer1_data* d = (layer1_data*)n->data;

                if (d != nullptr && d->val() < srcd->val()) {
                    d->r = srcd->r;
                    d->m = srcd->m;
                }
            } 
        }
    }
    inhibit(z);
}

void layer1_result::add_reconstruction_edges_back(int z)
{
    if (z < 0 || z > max_layer_index()) return;

    vector<node*>& s_nodes = shape_nodes[z];
    int to_prev = atom("toPrevLayer").get_index();
    int to_next = atom("toNextLayer0").get_index() + z;

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        do {
            set<node*> nbs;
        
            recurse_from_node(n, to_prev, nbs);
            for (set<node*>::iterator siter = nbs.begin(); siter != nbs.end(); ++siter) {
                add_edge_unique(*siter, n, to_next);
            }
            n = ((layer1_data*)n->data)->next;
        } while (n != nullptr);
    }
}

void layer1_result::add_reconstruction_edges(int z)
{
    if (z < 0 || z > max_layer_index()) return;

    vector<node*>& s_nodes = shape_nodes[z];
    int to_prev = atom("toPrevLayer").get_index();
    int to_0 = atom("toLayer0").get_index();
    int to_next = atom("toNextLayer0").get_index() + z;

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        do {
            set<node*> nbs;
            
            recurse_from_node(n, to_prev, nbs);
            for (set<node*>::iterator siter = nbs.begin(); siter != nbs.end(); ++siter) {
                add_edge_unique(n, *siter, to_0, to_next);
            }
            n = ((layer1_data*)(n->data))->next;
        } while (n != nullptr);
    }
}

void layer1_result::add_reconstruction_edges_link(int z)
{
    if (z < 0 || z > max_layer_index()) return;

    vector<node*>& s_nodes = shape_nodes[z];
    int to_prev = atom("toPrevLayer").get_index();
    int to_0 = atom("toLayer0").get_index();
    int to_next = atom("toNextLayer0").get_index() + z;

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        do {
            set<node*> nbs;
            
            recurse_and_link(n, to_prev, to_0, nbs);
            for (set<node*>::iterator siter = nbs.begin(); siter != nbs.end(); ++siter) {
                add_edge_unique(*siter, n, to_next);
            }
            n = ((layer1_data*)(n->data))->next;
        } while (n != nullptr);
    }
}

void layer1_result::add_reconstruction_edges_fwd(int z)
{
    if (z < 0 || z > max_layer_index()) return;

    vector<node*>& s_nodes = shape_nodes[z];
    int to_prev = atom("toPrevLayer").get_index();
    int to_0 = atom("toLayer0").get_index();

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        do {
            if (!n->has_neighbor(to_0)) {
                set<node*> nbs;
                
                recurse_from_node(n, to_prev, nbs);
                for (set<node*>::iterator siter = nbs.begin(); siter != nbs.end(); ++siter) {
                    if (node_layer(*siter) == 0) 
                        add_edge_unique(n, *siter, to_0);
                }
            }
            n = ((layer1_data*)(n->data))->next;
        } while (n != nullptr);
    }
}

void layer1_result::add_reconstruction_edges_fwd_link(int z)
{
    if (z < 0 || z > max_layer_index()) return;

    vector<node*>& s_nodes = shape_nodes[z];
    int to_prev = atom("toPrevLayer").get_index();
    int to_0 = atom("toLayer0").get_index();

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        do {
            recurse_and_link(n, to_prev, to_0);
            n = ((layer1_data*)(n->data))->next;
        } while (n != nullptr);
    }
}

void layer1_result::add_reconstruction_edges_leq_fwd(int z)
{
    if (z < 0 || z > max_layer_index()) return;

    int to_prev = atom("toPrevLayer").get_index();
    int to_0 = atom("toLayer0").get_index();

    for (int ly = 1; ly <= z; ++ly) {
        vector<node*>& s_nodes = shape_nodes[ly];
        
        for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
            node* n = *iter;

            do {
                set<node*> nbs;
                
                if (!n->has_neighbor(to_0)) {
                    n->get_neighbor_set(to_prev, nbs);
                    if (ly == 1) {
                        for (set<node*>::iterator siter = nbs.begin(); siter != nbs.end(); ++siter) {
                            add_edge_unique(n, *siter, to_0);
                        }
                    } else {
                        set<node*> nnbs;

                        for (set<node*>::iterator siter = nbs.begin(); siter != nbs.end(); ++siter) 
                            (*siter)->add_to_neighbor_set(to_0, nnbs);
                        for (set<node*>::iterator nsiter = nnbs.begin(); nsiter != nnbs.end(); ++nsiter) 
                            add_edge_unique(n, *nsiter, to_0);
                    }
                }
                n = ((layer1_data*)(n->data))->next;
            } while (n != nullptr);
        }
    }
}

int layer1_result::max_layer_index()
{
    return (int)shape_nodes.size() - 1;
}

int layer1_result::get_max_part_index(int z)
{
    if (z < 0 || z > max_layer_index()) return -1;

    int result = -1;
    vector<node*>& s_nodes = shape_nodes[z];

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        int m = ((layer1_data*)(*iter)->data)->m;

        if (m > result) result = m;
    }
    return result;
}

void layer1_result::get_parts(set<int>& parts, int z)
{
    parts.clear();
    if (z < 0 || z > max_layer_index()) return;

    vector<node*>& s_nodes = shape_nodes[z];

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd;
        do { 
            nd = (layer1_data*)n->data;
            if (!n->is_attr_set(HIDDEN_NODE_ATTR))
                parts.insert(nd->m);
            n = nd->next;
        } while (n != nullptr);
    }
}

void layer1_result::shift_part_id(layer1_result* res, const int layer, const int offset)
{
	vector<node*>& s_nodes = res->shape_nodes[layer];

    for(vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
		node* n = *iter;

        while(n != nullptr) {
			layer1_data* nd = (layer1_data*)n->data;
			nd->m += offset;
			n = nd->next;
		}
	}
}

double layer1_result::get_q_quantile_val(int layer, double invq)
{
    if (layer < 0 || layer > max_layer_index()) return 0.0;

    vector<node*>& s_nodes = shape_nodes[layer];
    vector<double> values;

    values.reserve(3*s_nodes.size());
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        while (n != nullptr && !n->is_attr_set(HAS_NEXT_LAYER)) {
            layer1_data* nd = (layer1_data*)n->data;
            
            values.push_back(nd->vval());
            n = nd->next;
        } 
    }
    sort(values.begin(), values.end(), greater<double>());

    if (values.empty()) return 0.0;

    size_t pos = (size_t)(invq*values.size());

    if (pos < 0) pos = 0;
    else if (pos >= values.size()) pos = values.size() - 1;

    return values[pos];
}

// It only merges toPrevLayer edges!!!
void layer1_result::merge(layer1_result* res, int border /* = 100 */, double mfactor /* = 1.0 */) 
{
	double sborder = border; // scaled border for this
	double sborder1 = border;
	double rsborder = border; // scaled border for res
	double rsborder1 = border;
    int to_prev = atom("toPrevLayer");

    for (int layer = 0; layer <= res->max_layer_index(); ++layer) {
		sborder1 = sborder;
		rsborder1 = rsborder;

		if (layer != 0) {
			sborder = (sborder * x_size(layer)) / (double)x_size(layer - 1);	// border scales with layer
			rsborder = (rsborder * res->x_size(layer)) / (double)res->x_size(layer - 1);
		}

		double factor = (double)(x_size(layer) - 2*sborder) / (res->x_size(layer) - 2*rsborder);
        double factor1 = (layer == 0) ? 1.0 : (double)(x_size(layer - 1) - 2*sborder1)/(res->x_size(layer - 1) - 2*rsborder1);
        vector<node*>& s_nodes = res->shape_nodes[layer];

        while (max_layer_index() < layer) {
            new_grid(res->x_size(max_layer_index() + 1), res->y_size(max_layer_index() + 1), layer);
            shape_nodes.push_back(vector<node*>());
            shape_nodes_inhib.push_back(vector<node*>());
            info.push_back(layer_info());
        }

        if (!grid(layer)) init_grid(layer);
        if (layer > 0 && !grid(layer - 1)) init_grid(layer - 1);
        for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
            node* n = *iter;
            
            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;
                int new_x = int_round(mfactor*factor*(nd->x - rsborder) + sborder);
                int new_y = int_round(mfactor*factor*(nd->y - rsborder) + sborder);

				if (new_x < x_size(layer) && new_y < y_size(layer)) {
					pair<node*, node**> fpair = find_at(nd, new_x, new_y, layer);
	                
					if (fpair.first == nullptr) {
						node* newn = add_grid_node(new layer1_data(nd->r, nd->m), new_x, new_y, layer);
                        int count = 0;
	                    
						forall_neighbors(n, niter) {
							if (neighbor_index(niter) != to_prev)
								continue;
							layer1_data* nnd = (layer1_data*)neighbor_node_data(niter);
                            int nnlayer = nnd->z;

							if (nnlayer >= layer) continue;

                            double nnborder = border*x_size(nnlayer)/(double)x_size(0);
                            double rnnborder = border*res->x_size(nnlayer)/(double)res->x_size(0);
                            double nnfactor = (double)(x_size(nnlayer) - 2*nnborder)/(res->x_size(nnlayer) - 2*rnnborder);

							int x = int_round(nnfactor*(nnd->x - rnnborder) + nnborder);
							int y = int_round(nnfactor*(nnd->y - rnnborder) + nnborder);
							pair<node*, node**> fpair1 = find_at(nnd, x, y, nnlayer);

                            if (fpair1.first != nullptr) { // though it should never be nullptr
                                edge_data_name* ed = dynamic_cast<edge_data_name*>(neighbor_edge_data(niter));

                                if (ed == nullptr) add_edge(newn, fpair1.first, neighbor_index(niter));
                                else add_edge_2(newn, fpair1.first, new edge_data_name(*ed), neighbor_index(niter));
                            } else {
                                //throw new_libhop_exception("fpair1.first is nullptr.");
                            }
                            ++count;
						}
                        //if (count == 0  && layer > 0) throw new_libhop_exception("no edges added!");

					}					
				}
				n = nd->next;
            } 
        }
        update_and_inhibit(layer);
    }
/*    for (iter_t aiter = nodes.begin(); aiter != nodes.end(); ++aiter) {
        node* n = *aiter; 
        layer1_data* nd = (layer1_data*)n->data;

        foreach_neighbor(n, to_prev, naiter) {
            edge_data_name* ed = dynamic_cast<edge_data_name*>(neighbor_edge_data(naiter));
            layer1_data* nnd = (layer1_data*)neighbor_node_data(naiter);

            if (ed == nullptr) cout << "From " << nd->x << ", " << nd->y << " ---> " << nnd->x << ", " << nnd->y << endl;
        }
    }
*/
}

/**
 * Reduces number of parts in layer k by only to taking parts with max G_RESPONSE
 * in local neighbourhood 3x3. Effectively there should only be one part in 3x3
 * neighbourhood for each layer. Note: that part can stil point to next part on 
 * same location but with different type.
 */
void layer1_result::update_and_inhibit(int k)
{
    int w = x_size(k);
    int h = y_size(k);
    layer1_data* d;
    matrix<int> forb(w, h, 0);
    int i, j;
    int w1 = w - 1, h1 = h - 1;

    vector<node*>& s_nodes = shape_nodes[k];
    vector<node*>& s_nodes_inhib = shape_nodes_inhib[k];

    // fill s_nodes
    list<node*>::iterator niter;

    s_nodes.clear();
    for (niter = nodes.begin(); niter != nodes.end(); ++niter) {
        d = (layer1_data*)(*niter)->data;
        if (d->z == k && d->x > 0 && d->x < w1 && d->y > 0 && d->y < h1
                && (*niter)->is_attr_set(IMG_NODE_ATTR)) { 
            s_nodes.push_back(*niter);
        }
    }
    sort(s_nodes.begin(), s_nodes.end(), layer1_data::greater1n); 


    // inhibit s_nodes
    vector<node*>::iterator iter;

    s_nodes_inhib.clear();
    for (iter = s_nodes.begin(); iter != s_nodes.end(); ++iter)  {
        d = (layer1_data*)(*iter)->data;
        i = d->x; j = d->y;
        if (forb(i, j) == 0) {
            s_nodes_inhib.push_back(*iter);
            (*iter)->set_attr(NODE_REDUCED_ATTR);
            forb(i - 1, j - 1) = 1;
            forb(i - 1, j) = 1;
            forb(i - 1, j + 1) = 1;
            forb(i, j - 1) = 1;
            forb(i, j) = 1;
            forb(i, j + 1) = 1;
            forb(i + 1, j - 1) = 1;
            forb(i + 1, j) = 1;
            forb(i + 1, j + 1) = 1;
        }
    }
}

img* layer1_result::get_image(set<node*>& nodes, int x_size, int y_size,
    bool paintuncov, const color& uncovered, const color& defcol, 
    const vector<int>& parts, const vector<color>& colors, bool bigpoints /* = false*/)
{
    img* result = new img(x_size, y_size, 0.0, false);  // colored image
    int R, G, B, uR, uG, uB;

    if (paintuncov) uncovered.to_rgb(uR, uG, uB);
    for (set<node*>::iterator i = nodes.begin(); i != nodes.end(); ++i) {
        node* n = *i;
        layer1_data* d = (layer1_data*)n->data;

        if (paintuncov && !n->is_attr_set(HAS_NEXT_LAYER)) {
            result->set_color(d->x, d->y, d->val()*uR, d->val()*uG, d->val()*uB);
        } else {
            if (parts.empty()) defcol.to_rgb(R, G, B);
            else {
                vector<int>::const_iterator result;

                while (n != nullptr) {
                    d = (layer1_data*)n->data;
                    if ((result = find(parts.begin(), parts.end(), d->m)) != parts.end()) break;
                    n = d->next;
                }
                if (n != nullptr) {
                    size_t index = result - parts.begin();
                    if (index < colors.size()) colors[index].to_rgb(R, G, B); 
                    else defcol.to_rgb(R, G, B);
                } else {
                    n = *i;
                    d = (layer1_data*)n->data;
                    defcol.to_rgb(R, G, B);
                }
            } 
			HOP_REAL rcol(COL_TO_REAL(d->val()*R, d->val()*G, d->val()*B));
            if (bigpoints) result->draw_big_point(d->x, d->y, rcol);
			else result->set_color(d->x, d->y, d->val()*R, d->val()*G, d->val()*B);
		}
    }
    return result;
}

img* layer1_result::get_image_new(int z,
        const vector<int>& parts, const color& defcol, const vector<color>& colors, bool bigpoints /* = false*/)
{
	img* result = new img(x_size(0), y_size(0), 0.0, false);  // colored image
    int R, G, B, uR, uG, uB;
	
	if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);

	set<node*> end_nodes;
    layer_predicate pred(0);
	recurse(shape_nodes[z], atom("toPrevLayer").get_index(), pred, end_nodes);
	//color all with default color
	for (set<node*>::iterator i = end_nodes.begin(); i != end_nodes.end(); ++i) {
        node* n = *i;
        layer1_data* d = (layer1_data*)n->data;

        defcol.to_rgb(R, G, B);
		HOP_REAL rcol(COL_TO_REAL(d->val()*R, d->val()*G, d->val()*B));
		if (bigpoints) result->draw_big_point(d->x, d->y, rcol);
		else result->set_color(d->x, d->y, d->val()*R, d->val()*G, d->val()*B);
	}
	//color parts
	for (int i = 0; i < (int)parts.size(); ++i) {
		set<node*> nodes;
		set<node*> endnodes;

		get_layer_nodes(nodes, z, vector<int>(1, parts[i]), -1.0, -1.0, -1.0, 100.0);
		recurse(nodes, atom("toPrevLayer").get_index(), pred, endnodes);
		for(set<node*>::iterator niter = endnodes.begin(); niter != endnodes.end(); ++niter) {
			layer1_data* d = (layer1_data*)(*niter)->data;
			colors[i].to_rgb(R, G, B); 
			HOP_REAL rcol(COL_TO_REAL(d->val()*R, d->val()*G, d->val()*B));
			if (bigpoints) {
				result->draw_big_point(d->x, d->y, rcol);
			} else {
				result->set_color(d->x, d->y, d->val()*R, d->val()*G, d->val()*B);
			}
		}
	}
    return result;
}

img* layer1_result::get_image_inhib_simple(int z)
{
    img* result = new img(x_size(z), y_size(z), 0.0);
    vector<node*>& nodes = shape_nodes_inhib[z];

    for (unsigned i = 0; i < nodes.size(); ++i) {
        node* n = nodes[i];
        layer1_data* d = (layer1_data*)n->data;

        result->at(d->x, d->y) = d->val(); //set_color(d->x, d->y, d->val, d->val, d->val);
    }
    return result;
}

img* layer1_result::get_image_inhib_simple_01(int z)
{
    img* result = new img(x_size(z), y_size(z), 0.0);
    vector<node*>& nodes = shape_nodes_inhib[z];

    for (unsigned i = 0; i < nodes.size(); ++i) {
        node* n = nodes[i];
        layer1_data* d = (layer1_data*)n->data;

        if (d->val() > 0) 
            result->at(d->x, d->y) = 1.0; 
    }
    return result;
}

img* layer1_result::get_image(int z, int zt, bool paintuncov, const color& uncovered, const color& defcol, 
    const vector<int>& parts, const vector<color>& colors, bool bigpoints /* = false*/)
{
    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);
    if (zt < 0) zt = 0; else zt = min(zt, z);

    set<node*> end_nodes;
    layer_predicate pred(zt);

    recurse(shape_nodes[z], atom("toPrevLayer").get_index(), pred, end_nodes);
    return get_image(end_nodes, x_size(zt), y_size(zt), paintuncov, uncovered, defcol, parts, colors, bigpoints);
}

img* layer1_result::get_image_inhib(int z, int zt, bool paintuncov, const color& uncovered, const color& defcol, 
    const vector<int>& parts, const vector<color>& colors)
{
    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);
    if (zt < 0) zt = 0; else zt = min(zt, z);

    set<node*> end_nodes;
    layer_predicate pred(zt);

    recurse(shape_nodes_inhib[z], atom("toPrevLayer").get_index(), pred, end_nodes);
    return get_image(end_nodes, x_size(zt), y_size(zt), paintuncov, uncovered, defcol, parts, colors);
}        

img* layer1_result::get_part_image(int z, int zt, node* p)
{
    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);
    if (zt < 0) zt = 0; else zt = min(zt, z);

    set<node*> end_nodes;
    layer_predicate pred(zt);
    vector<node*> start_node;

    start_node.push_back(p);

    recurse(start_node, atom("toPrevLayer").get_index(), pred, end_nodes);
    return get_image(end_nodes, x_size(zt), y_size(zt), false, color(20, 20, 20), color(255, 255, 255), 
        vector<int>(), vector<color>());
}

// Get bounding rectangle of all "reconstruction" nodes of 'n'.
irectangle2 layer1_result::get_box_with_cached_link(node* n)
{
	set<node*> starting_nodes;
    set<node*> nodes;
    irectangle2 result;

    recurse_and_link(n, atom("toPrevLayer"), atom("toLayer0"), nodes);
    for (set<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;
        if (nd->z == 0) result.eat(nd->x, nd->y);
    }
    return result;
}



// Get bounding rectangle of all "reconstruction" nodes of 'n'.
irectangle2 layer1_result::get_box(node* n)
{
    set<node*> nodes;
    irectangle2 result;

    recurse_from_node(n, atom("toPrevLayer"), nodes);
    for (set<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;
        if (nd->z == 0) result.eat(nd->x, nd->y);
    }
    return result;
}

// Get bounding rectangle of all "reconstruction" nodes of 'n'.
// Here we assume that "toLayer0" edges are present.
irectangle2 layer1_result::get_box_0(node* n)
{
    set<node*> nodes;
    irectangle2 result;

	n->get_neighbor_set(atom("toLayer0"), nodes);
    for (set<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;
        if (nd->z == 0) result.eat(nd->x, nd->y);
    }
    return result;
}

// p: library part for type n (no checking)
// 
vector<ipoint2> robust_PCA(layer1_result* res, node* n, node* p)
{
    static map<int, vector<pair<int, ipoint2> > > libptsmap;

    vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);

    if (pd == nullptr) {
        cout << "Non-vs-part, skipping." << endl;
        return vector<ipoint2>();
    }

    vector<pair<int, ipoint2> >& ipts = libptsmap[pd->type];
    vector<pair<int, ipoint2> > ptsm;

    get_node_geo(ptsm, res, n);
    ptsm = inhibit_point_set(ptsm, 5);

    // DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG 
    //cv::Mat debugim = point_image(cast_vector<dpoint2, ipoint2>(extract_second<ipoint2>(ptsm.begin(), ptsm.end())));
    //cv::imwrite("c:\\work\\tmp\\orig.png", debugim);
    // DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG 

    pair<double, dpoint2> delta = translate_and_scale(ptsm);

    if (ipts.empty()) {
        ipts = get_library_geo_pieces(p, pd->layer);
        translate_and_scale(ipts);
    }

    vector<int> match = piecewise_point_matching_p(ptsm, ipts);

    set<int> matchedpcs, allpcs;
    vector<dpoint2> dptsm(ipts.size());
    vector<int> dptsmind;

    dptsmind.reserve(dptsm.size());
    for (int i = 0; i < (int)dptsm.size(); ++i) {
        int j = match[i];

        if (j >= 0 && j < (int)ptsm.size()) {
            dptsm[i] = (dpoint2)ptsm[j].second;
            matchedpcs.insert(ipts[i].first);
            dptsmind.push_back(i);
        } /* else {
            dptsm[i] = (dpoint2)ipts[i].second; // + dpoint2(random_real(-1, 1), random_real(-1, 1));
        }*/
        allpcs.insert(ipts[i].first);
    }
    dptsm = take(dptsm, dptsmind);

    translate_and_scale(dptsm);
    cv::Mat1d B = flatten(dptsm).t();
    cv::Mat1d A(B.rows, pd->pcad.eigenvalues.rows, 0.0);
    cv::Mat1d R;

    for (int j = 0; j < pd->pcad.eigenvalues.rows; ++j) {
        vector<dpoint2> v = partition(pd->pcad.eigenvectors.row(j));
        v = take(v, dptsmind);

        cv::Mat1d vf = flatten(v);

        for (int i = 0; i < vf.cols; ++i) 
            A.at<double>(i, j) = vf.at<double>(0, i);
    }
    cv::solve(A, B, R, cv::DECOMP_SVD);

    cv::Mat1d result = pd->pcad.mean.clone();

    for (int i = 0; i < R.rows; ++i) {
        result = result + pd->pcad.eigenvectors.row(i)*R.at<double>(i, 0);
    }

    vector<dpoint2> resultd = partition(result);

    for (auto riter = resultd.begin(); riter != resultd.end(); ++riter) {
        *riter *= delta.first;
        *riter += delta.second;
    }
    // DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG 
    //cv::Mat debugim2 = point_image(resultd);
    //cv::imwrite("c:\\work\\tmp\\trans.png", debugim2);
    // DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG 
    return cast_vector<ipoint2, dpoint2>(resultd);
}

// Returns all boxes on the z-th layer with their value set to certain response.
// Boxes are "resized" by 'factor', 'border0' is border0 of result resized to.
// Structure of type 'response_filter' is used to filter nodes.
// Note: new boxes are added to the list!
//	void get_boxes(list<box_data_t>& boxes, int z, int response, part_lib* library, const response_filter& rfilter, 
//        double factor, int border0);

//    void get_boxes(list<box_data_t>& boxes, part_lib* library, int z, int response, bool use_lib_responses, bool revPCA,
//        const response_filter& rfilter, double factor, int border0);
void layer1_result::get_boxes(list<box_data_t>& boxes, part_lib* library, int z, int response, bool use_lib_responses, 
    bool revPCA, const response_filter& rfilter, double factor, int border0)
{
	if (z < 0 || z > max_layer_index()) return;

	int edgename = atom("toPrevLayer");
    int to0 = atom("toLayer0");

	for (list<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
		node* n = *iter;
		layer1_data* nd = (layer1_data*)n->data;

		if (n->is_attr_set(NODE_DELETED_ATTR))
			continue;
        if (nd->z == z) {
            node* nn = n->get_neighbor(edgename);
            layer1_data* nnd = (layer1_data*)nn->data;

            if (rfilter.check(nnd, use_lib_responses ? library : nullptr)) {
            //if (rfilter.check(nd)) {
			    irectangle2 box;
			    set<node*> rec;

                recurse_and_link(n, edgename, to0, rec);
                if (revPCA && library != nullptr) {
                    node* p = library->parts[nnd->z][nnd->m];
                    vector<ipoint2> pts = robust_PCA(this, nn, p);

                    for (auto piter = pts.begin(); piter != pts.end(); ++piter)
                        box.eat(*piter);
                }

                if (box.invalid()) {
			        //recurse_from_node(n, edgename, rec);
			        for (set<node*>::iterator niter = rec.begin(); niter != rec.end(); ++niter) {
				        layer1_data* nnnd = (layer1_data*)(*niter)->data;

				        if (nnnd->z == 0) box.eat(nnnd->x, nnnd->y);
			        }
                }
                box -= ipoint2(border, border);
                box.resize(factor);
                box += ipoint2(border0, border0);
			    boxes.push_back(box_data_t(box, nd->r, nd->r(response), (int)rec.size(), nd->m, nnd->m));
            }
		}
	}		 
}

// Gets predicted bounding boxes of selected (filtered) nodes. See 'get_boxes' above. 
// 'pred_boxes' are predctions for each part. 
void layer1_result::get_predicted_boxes(list<box_data_t>& boxes, int z, int response, 
    const response_filter& rfilter, const vector<irectangle2>& pred_boxes, double factor, int border0)
{
	if (z < 0 || z > max_layer_index()) return;

	int edgename = atom("toPrevLayer");

    double layerfactor = (double)x_size(0)/x_size(z - 1);

	for (list<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
		node* n = *iter;
		layer1_data* nd = (layer1_data*)n->data;

		if (n->is_attr_set(NODE_DELETED_ATTR))
			continue;
        if (nd->z == z) {
            node* nn = n->get_neighbor(edgename);
            layer1_data* nnd = (layer1_data*)nn->data;

            if (rfilter.check(nd)) {
                if ((int)pred_boxes.size() <= nnd->m) {
                    cout << "There is no predicted bounding box for part #" << nnd->m << endl;
                    throw;
                }

			    irectangle2 box = pred_boxes[nnd->m];
                
                box += ipoint2(int_round(layerfactor*nnd->x), int_round(layerfactor*nnd->y));
                box -= ipoint2(border, border);
                box.resize(factor);
                box += ipoint2(border0, border0);
			    boxes.push_back(box_data_t(box, nd->r, nd->r(response), 0, nd->m, nnd->m));
            }
		}
	}		 
}

void layer1_result::get_boxes(list<box_data_t>& boxes, int z, const set<int>& types, double inhibit, int resp_type /*= G_RESPONSE */ )
{
    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);

    vector<node*>::iterator iter = shape_nodes[z].begin(), end = shape_nodes[z].end();

    boxes.clear();
    for (; iter != end; ++iter) {
        node* n = *iter;

        while (n != nullptr && n->data != nullptr) {
            layer1_data* nd = (layer1_data*)n->data;

            if (types.empty() || types.find(nd->m) != types.end()) {
                rectangle2<int> box;
                set<node*> end_nodes;
                node* nn = n->get_neighbor(atom("toPrevLayer"));
                int nnm = (nn == nullptr) ? -1 : node_type(nn);

                recurse_from_node(n, atom("toPrevLayer").get_index(), end_nodes);
                for (set<node*>::iterator siter = end_nodes.begin(); siter != end_nodes.end(); ++siter) {
                    layer1_data* d = (layer1_data*)(*siter)->data;
                    if (d->z == 0) box.eat(d->x, d->y);
                }
                if (!box.invalid()) boxes.push_back(box_data_t(box, nd->r, nd->r(resp_type), (int)end_nodes.size(), nd->m, nnm));
            }
            n = nd->next;
        }
    }
	
    if (inhibit > 0.0) 
        irectangle2::inhibit(boxes, inhibit, box_data_t::box_data_f());
}

void layer1_result::get_boxes_inhibited(list<box_data_t>& boxes, int z, 
    const set<int>& types, double inhibit, double thresh)
{
    typedef pair<double, node*> pair_t;

    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);

    vector<node*>::iterator iter, end = shape_nodes[z].end();
    vector<pair_t> order;
    set<node*> hit;

    for (iter = shape_nodes[z].begin(); iter != end; ++iter) {
        node* n = *iter;
        layer1_data* nd;

        while (n != nullptr && n->data != nullptr) {
            nd = (layer1_data*)n->data;
            order.push_back(pair_t(nd->val(), n));
            n = nd->next;
        }
    }
    sort(order.begin(), order.end(), greater<pair_t>());

    boxes.clear();
    for (vector<pair_t>::iterator iter = order.begin(); iter != order.end(); ++iter) {
        //    for (iter = shape_nodes[z].begin(); iter != end; ++iter) {
        node* n = iter->second;
        layer1_data* nd = (layer1_data*)n->data;
        set<node*> end_nodes;
        rectangle2<int> box;
        node* nn = n->get_neighbor(atom("toPrevLayer"));
        int nnm = (nn == nullptr) ? -1 : node_type(nn);

        recurse_from_node(n, atom("toPrevLayer").get_index(), end_nodes);
        if ((types.empty() || types.find(nd->m) != types.end()) &&
                (intersection_size(end_nodes, hit) <= (int)(thresh * end_nodes.size()))) {
            for (set<node*>::iterator siter = end_nodes.begin(); siter != end_nodes.end(); ++siter) {
                layer1_data* d = (layer1_data*)(*siter)->data;
                if (d->z == 0) box.eat(d->x, d->y);
            }
        }
        if (!box.invalid()) boxes.push_back(box_data_t(box, nd->r, nd->val(), (int)end_nodes.size(), nd->m, nnm));
        hit.insert(end_nodes.begin(), end_nodes.end());
    }
    if (inhibit > 0.0) irectangle2::inhibit(boxes, inhibit, box_data_t::box_data_f());
}

img* layer1_result::get_image_boxed(int z, const vector<int>& parts, const vector<color>& c, 
    const color& defcol, bool inhibit)
{
    typedef pair<irectangle2, color> box_data_t;

    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);

    img* res = get_image_reconstructed(0, 0, vector<int>(), true, false);
    set<node*> end_nodes;
    vector<int> sparts(parts.begin(), parts.end());

    sort(sparts.begin(), sparts.end());
    vector<node*>::iterator iter = shape_nodes[z].begin(), end = shape_nodes[z].end();
    list<box_data_t> boxes;

    for (; iter != end; ++iter) {
        node* n = *iter;

        while (n != nullptr && n->data != nullptr) {
            layer1_data* nd = (layer1_data*)n->data;

            if (sparts.empty() || binary_search(sparts.begin(), sparts.end(), nd->m)) {
                rectangle2<int> box;

                end_nodes.clear();
                recurse_from_node(*iter, atom("toPrevLayer").get_index(), end_nodes);
                for (set<node*>::iterator siter = end_nodes.begin(); siter != end_nodes.end(); ++siter) {
                    layer1_data* d = (layer1_data*)(*siter)->data;
                    if (d->z == 0) box.eat(d->x, d->y);
                }
                if (!box.invalid()) 
                    boxes.push_back(box_data_t(box, (nd->m < (int)c.size()) ? c[nd->m] : defcol)); 
            }
            n = nd->next;
        } 
    }
    if (inhibit) irectangle2::inhibit(boxes, 0.1);
    for (list<box_data_t>::iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
        res->draw_box(biter->first, biter->second);
    }
    return res;
}

img* layer1_result::get_part_reconstructed(int z, int zt, node* p)
{
    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);
    if (zt < 0) zt = 0; else zt = min(zt, z);

    img* res = new img(x_size(0), y_size(0), 0.0);

    vector<img*> masks;
    set<node*> end_nodes;
    layer_predicate pred(zt);
    vector<node*> start_nodes;

    start_nodes.push_back(p);
    get_real_masks(masks);

    //
    layer1_data* d = (layer1_data*)p->data;
    cout << '(' << d->x << ',' << d->y << ')' << '-' << d->val() << ' ';
    //

    recurse(start_nodes, atom("toPrevLayer").get_index(), pred, end_nodes);
    for (set<node*>::iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter) {
        layer1_data* d = (layer1_data*)(*iter)->data;

        if (d->z != 0) continue;
        img* m = masks[d->m];
        
        res->combine_max(*m, d->x - (int)(m->width)/2, d->y - (int)(m->height)/2, pow(d->val(), 1.0));
    }
    for (unsigned i = 0; i < masks.size(); i++) 
        if (masks[i]) delete masks[i];
    return res;
}

struct predicate_recurse3_best {
    int edge_name, lyr;

    predicate_recurse3_best(int name, int l = 0) : edge_name(name), lyr(l) { }

    int operator()(set<node*>& neighbors, node* n) const 
    { 
        typedef pair<node*, double> value_t;

        map<int, value_t> bestm;
        pair<map<int, value_t>::iterator, bool> ins_pair;
        double bestval = 0.0;
        

        foreach_neighbor(n, edge_name, i) {
            layer1_data* nd = (layer1_data*)neighbor_node_data(i);
            node* nn = neighbor_node(i);

            if (nd->z < lyr) continue;

            ins_pair = bestm.insert(pair<int, value_t>(nd->m, value_t(nn, nd->val())));
            if (!ins_pair.second && nd->val() > ins_pair.first->second.second) {
                ins_pair.first->second.first = nn;
                ins_pair.first->second.second = nd->val();
            }
        }
        for (map<int, value_t>::iterator iter = bestm.begin(); iter != bestm.end(); ++iter)
            neighbors.insert(iter->second.first);
        return (int)bestm.size();
    }
};


void layer1_result::part_reconstruction(map<int, img>& result, int z)
{
    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);

    vector<img*> masks;
    //plus<HOP_REAL> operation;
    max_function<HOP_REAL> operation;
    vector<node*> start_nodes(1, (node*)nullptr);
    vector<node*>& s_nodes = shape_nodes[z];
    set<node*> end_nodes;
    layer_predicate pred(0);
    int edge_name = atom("toPrevLayer").get_index();
    double factor = (double)x_size(0)/x_size(z);

    get_real_masks(masks);
    for (vector<node*>::iterator niter = s_nodes.begin(); niter != s_nodes.end(); ++niter) {
        node* p = *niter;
        start_nodes[0] = p;
        layer1_data* pd = (layer1_data*)p->data;
        map<int, img>::iterator miter = result.find(pd->m);
        int dx = (int)(factor * pd->x);
        int dy = (int)(factor * pd->y);

        if (miter == result.end()) 
            miter = result.insert(pair<int, img>(pd->m, img())).first;

        img& mat = miter->second;
        img temp;

        end_nodes.clear();
        recurse(start_nodes, edge_name, pred, end_nodes);
        for (set<node*>::iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter) {
            layer1_data* d = (layer1_data*)(*iter)->data;

            if (d->z != 0) continue;

            img* m = masks[d->m];

            temp.blt_central2_ar(*m, (int)m->width/2, (int)m->height/2, 
                d->x - dx, d->y - dy, operation);
        }
        if (temp.count_geq(0.1) > mat.count_geq(0.1)) 
            mat = temp;
        //mat.save("c:\\test.png");
        //int x;
        //cin >> x;
        //ipoint2 center = simple_shape_matching(mat, temp, ipoint2((int)mat.width/2, (int)mat.height/2), 
        //    ipoint2((int)temp.width/2, (int)temp.height/2), (int)(factor + 1.0));
        //cout << '.';
        //mat.blt_central2_ar(temp, center.x, center.y, (int)mat.width/2, (int)mat.height/2, operation);
        
    }
    for (unsigned i = 0; i < masks.size(); i++) 
        if (masks[i]) delete masks[i];
}

img* layer1_result::get_image_reconstructed(int z, int zt, const vector<int>& parts,
    bool colored /* = false */, bool drawbox /* = false */, double factorpow /* = 1.0 */, bool all /* = false */)
{
    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);
    if (zt < 0) zt = 0; else zt = min(zt, z);

    img* res = new img(x_size(0), y_size(0), 0.0, !colored);

    if (shape_nodes[z].empty()) return res;

    vector<img*> masks;
    set<node*> end_nodes;
    vector<int> sparts(parts.begin(), parts.end());
    layer_predicate pred(zt);

    get_real_masks(masks);
    sort(sparts.begin(), sparts.end());
    rectangle2<int> box;
    vector<node*> all_nodes;
    combine_t combinef = get_combine_function();

    if (all) all_grid_nodes(all_nodes, shape_nodes[z]); 
    else all_nodes.insert(all_nodes.begin(), shape_nodes[z].begin(), shape_nodes[z].end());
    if (sparts.empty()) {
        recurse(all_nodes, atom("toPrevLayer").get_index(), pred, end_nodes);
    } else {
        vector<node*> srcnodes;
        int type;

        for (vector<node*>::iterator iter = all_nodes.begin(); iter != all_nodes.end(); ++iter) {
            type = ((layer1_data*)(*iter)->data)->m;
            if (binary_search(sparts.begin(), sparts.end(), type)) srcnodes.push_back(*iter);
        }
        recurse(srcnodes, atom("toPrevLayer").get_index(), pred, end_nodes);
    }
    for (set<node*>::iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter) {
        layer1_data* d = (layer1_data*)(*iter)->data;

        if (d->z != 0) continue;
        img* m = masks[d->m];
        
        (res->*combinef)(*m, d->x - (int)(m->width)/2, d->y - (int)(m->height)/2, pow(d->vval(), factorpow));
        box.eat(d->x, d->y);
    }
    //if (drawbox) res->draw_box(box, 1.0);
    for (unsigned i = 0; i < masks.size(); i++) 
        if (masks[i]) delete masks[i];
    return res;
}   

void layer1_result::distance_clustering(vector<vector<ipoint2> >& clusters, int layer, int cluster_n)
{
    clusters.clear();

    if (layer < 0 || layer >= layer_count()) return;

    vector<node*>& s_nodes = shape_nodes[layer];

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;

        clusters.push_back(vector<ipoint2>());
        clusters.back().push_back(ipoint2(nd->x, nd->y));
    }
    while ((int)clusters.size() > cluster_n) {
        hierarchical_clustering_iter(clusters, ipoint2_set_distance);
    }
}

void layer1_result::get_layer_nodes(set<node*>& result, int z, const vector<int>& parts, 
    double thresh, double thresh2, double thresh3, double thresh4, dpoint2* within_bounds)
{
    //cout << "get_layer_nodes 1" << endl;

    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);

    result.clear();
    if (shape_nodes[z].empty()) return;

    vector<int> sparts(parts.begin(), parts.end());
    vector<node*>& s_nodes = shape_nodes[z];
    bool first_only = thresh < 0.0;

    sort(sparts.begin(), sparts.end());

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;
		double dist = within_bounds != nullptr ? dpoint2::distance2(*within_bounds, dpoint2(nd->x,nd->y)) : -1;

        do {
            layer1_data* nd = (layer1_data*)n->data;

            if (!n->is_attr_set(HIDDEN_NODE_ATTR)) {
                if (first_only) {            
                    if (binary_search(sparts.begin(), sparts.end(), nd->m) && (within_bounds == nullptr || dist < 10)) 
	                result.insert(n);
                    break;
                } 
                if (nd->r(R_RESPONSE) >= thresh && 
                        nd->r.get_response(G_RESPONSE, 1.0) >= thresh2 && 
                        nd->r.get_response(RR_RESPONSE, 1.0) >= thresh3 &&
                        nd->r.get_response(S_RESPONSE, 5.0) <= thresh4 &&
    				    binary_search(sparts.begin(), sparts.end(), nd->m) && 
                        (within_bounds == nullptr || dist < 10 ) )
                    result.insert(n);
            } //else cout << "hidden node not inserted" << endl;
            n = nd->next;
        } while (n != nullptr);
        //}
    }   
}   

void layer1_result::get_layer_nodes(map<int, set<node*> >& nodes, int z, const vector<int>& parts, 
    double thresh, double thresh2, double thresh3, double thresh4, dpoint2* within_bounds)
{
    //cout << "get_layer_nodes 2; z = " << z << endl;

	//get nodes in a part map
	typedef map<int,set<node*> > nmap;

	if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);
	nodes.clear();
	if (shape_nodes[z].empty()) return;

    vector<node*>& s_nodes = shape_nodes[z];
    bool first_only = thresh < 0.0;

	for(vector<int>::const_iterator iter = parts.begin(); iter != parts.end(); ++iter) {
		nodes.insert(pair<int, set<node*> >(*iter, set<node*>()));
	}
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;
		// calculate distance only for first one since all others should be on same location
		double dist = within_bounds != nullptr ? dpoint2::distance2(*within_bounds, dpoint2(nd->x,nd->y)) : -1;

        do {
			nd = (layer1_data*)n->data;
            if (!n->is_attr_set(HIDDEN_NODE_ATTR)) {
                if (first_only) {
			        if (within_bounds == nullptr || dist < 10) {
				        nmap::iterator el = nodes.find(nd->m);

				        if (el != nodes.end()) el->second.insert(n);
			        }
                    break;
                }
                if (nd->r(R_RESPONSE) >= thresh && 
                        nd->r.get_response(G_RESPONSE) >= thresh2 && 
                        nd->r.get_response(RR_RESPONSE) >= thresh3 && 
                        nd->r.get_response(S_RESPONSE) <= thresh4 &&
                        (within_bounds == nullptr || dist < 10 )) {
				    nmap::iterator el = nodes.find(nd->m);

				    if(el != nodes.end()) el->second.insert(n);
				}
			} // else cout << "hidden node not inserted" << endl;
            n = nd->next;				
        } while (n != nullptr);
    }
}

void layer1_result::get_reconstruction_nodes(set<node*>& result, int z, const vector<int>& parts)
{
    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);

    result.clear();
    if (shape_nodes[z].empty()) return;

    layer_predicate pred(0);
    vector<int> sparts(parts.begin(), parts.end());
    vector<node*> srcnodes;
    vector<node*>& s_nodes = shape_nodes[z];
    int type;

    sort(sparts.begin(), sparts.end());
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        type = ((layer1_data*)(*iter)->data)->m;
        if (binary_search(sparts.begin(), sparts.end(), type)) srcnodes.push_back(*iter);
    }
    recurse(srcnodes, atom("toPrevLayer").get_index(), pred, result);
}

img* layer1_result::get_image_reconstructed(part_lib* library, int z, const vector<int>& parts)
{
    z = min(z, (int)shape_nodes.size() - 1);

    img* res = new img(x_size(0), y_size(0), 0.0, false);

    vector<int> sparts(parts.begin(), parts.end());

    sort(sparts.begin(), sparts.end());

    vector<node*>& s_nodes = shape_nodes[z];

    if (s_nodes.empty()) return res;

    vector<node*>& lparts = library->parts[z];
    vector<node*>::iterator iter = s_nodes.begin();
    layer_predicate pred(0);
    node* n = *iter;
    HOP_REAL red = COL_TO_REAL(255, 0, 0);
    int to_prev = atom("toPrevLayer").get_index();

    while (n != nullptr) {
        layer1_data* nd = (layer1_data*)n->data;

        if (sparts.empty() || binary_search(sparts.begin(), sparts.end(), nd->m)) {
            part_data* pd = (part_data*)lparts[nd->m]->data;
            set<node*> end_nodes;
            ipoint2 p(0, 0);
            int count = 0;
            matrix<double> mask;
            ipoint2 maskc;

            maskc = pd->get_mask(mask, lparts[nd->m], library);
            end_nodes.clear();
            recurse_from_node(n, to_prev, pred, end_nodes);
            for (set<node*>::iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter) {
                layer1_data* nd = (layer1_data*)(*iter)->data;
                p.add(nd->x, nd->y);
                ++count;
            }
            int x = int_round(((double)p.x)/count);
            int y = int_round(((double)p.y)/count);

            res->combine_max(img(mask, false), x - maskc.x, y - maskc.y, 1.0);
            res->draw_big_point(x, y, red);
        }

        n = nd->next;
        if (n == nullptr && ++iter != s_nodes.end()) n = *iter;
    }
    return res;
}

class MinMaxNodeFilter {
private:
	node* min_x_node;
	node* min_y_node;
	node* max_x_node;
	node* max_y_node;

	set<node*> nset;
public:
	typedef set<node*>::iterator iterator;

	MinMaxNodeFilter() {
		this->min_x_node = nullptr; this->max_x_node = nullptr;
		this->min_y_node = nullptr; this->max_y_node = nullptr;
	}
	iterator begin() {
		if (nset.empty()) {
			nset.insert(min_x_node); nset.insert(min_y_node); nset.insert(max_x_node); nset.insert(max_y_node);
		}
		return nset.begin();
	}
	iterator end() {
		return nset.end();
	}
	void insert(iterator start, iterator end) {
		// insert only min or max value
		for (iterator it = start; it != end; it++) {
			node* n = *it;			
			layer1_data* nd = (layer1_data*)n->data;

			if (min_x_node == nullptr) min_x_node = n;
			if (min_y_node == nullptr) min_y_node = n;
			if (max_x_node == nullptr) max_x_node = n;
			if (max_y_node == nullptr) max_y_node = n;

			layer1_data* min_x_data = (layer1_data*)min_x_node->data;
			layer1_data* min_y_data = (layer1_data*)min_y_node->data;
			layer1_data* max_x_data = (layer1_data*)max_x_node->data;
			layer1_data* max_y_data = (layer1_data*)max_y_node->data;

			min_x_node = min_x_data->x > nd->x ? n : min_x_node;
			min_y_node = min_y_data->y > nd->y ? n : min_y_node;
			max_x_node = max_x_data->x < nd->x ? n : max_x_node;
			max_y_node = max_y_data->y < nd->y ? n : max_y_node;
		}
	}
	bool empty() {
		return min_x_node == nullptr || min_y_node == nullptr || max_x_node == nullptr || max_y_node == nullptr;
	}
};

void layer1_result::save_visview(const string& dir, const string& nameint, const string& fname, 
    part_lib* library, int z, int zmax, const map<node*,int>& prevlayermap, int save_mode)
{
    if (z > zmax || empty(z)) return;

    int edgename = atom("toPrevLayer").get_index();	
    //int edge0name = atom("toLayer0").get_index();
	int minmaxedge0name = atom("toMinMaxLayer0").get_index();
    vector<irectangle2> librectangles;

    if (library != nullptr) {
        vector<node*>& parts = library->parts[z];

        librectangles.resize(parts.size());
        for (int i = 0; i < (int)parts.size(); ++i) {
            part_data* pd = (part_data*)parts[i]->data;
            librectangles[i] = pd->get_mask_box(parts[i], library);
        }
    }

	string name, names;

	if ((save_mode & VVE_SIMPLE_FILENAMES) != 0) {
		stringstream scale_num(get_extension(fname, "_").substr(1));			
			
		if ((save_mode & VVE_MATLAB_SCALES) != 0) {
			int num;
			scale_num >> num;
			num++;
			
			stringstream matlab_scale_num;
			matlab_scale_num << num ;
			scale_num.str(matlab_scale_num.str());
		}
		
        name = dir + change_extension(fname, "", "_") + "_" + (z + 1) + "_" + scale_num.str() + ".txt";
        names = dir + change_extension(fname, "", "_") + "_" + (z + 1) + "_" + scale_num.str() + "_factor.txt";
    } else {
        name = dir + nameint + "_image_" + (z + 1) + "_1.txt";
        names = dir + nameint + "_image_" + (z + 1) + "_1_factor.txt";
    }

    ofstream oss(names.c_str());
    oss << (double)(x_size(0) - 2*border)/original_width << '\n';
    oss.close();

    ofstream os;

    os.open(name.c_str(), ios::trunc);
    if (os.fail()) return;

    vector<node*>& nodes = shape_nodes[z];
    double factor = 1.0;
    int count = (int)nodes.size();
    map<node*,int> layermap;
    int linei = 0;

    if ((save_mode & VVE_TRUE_COORDINATES) != 0) 
        factor = (double)x_size(0)/x_size(z);

    for (int i = 0; i < count; ++i) {
        node* n = nodes[i];
        layer1_data* d;

        do {
            d = (layer1_data*)(n->data);

            if (!n->is_attr_set(NODE_DELETED_ATTR)) {
                layermap.insert(pair<node*,int>(n, ++linei));

                int xpos = (int)(factor*d->x);
                int ypos = (int)(factor*d->y);

                if ((save_mode & VVE_BOUNDING_BOXES) == 0) {
                    os << xpos + 1 << ',' << ypos + 1 << ",1," << d->m + 1 << ',';
                } else {
                    xpos -= border;
                    ypos -= border;

			        set<node*> nset;								
                    irectangle2 box;

                    recurse_and_link<MinMaxNodeFilter>(n, edgename, minmaxedge0name, nset);
                    box = node_set_bounding_rectangle(nset.begin(), nset.end());
                    box.ll -= border;
                    box.ur -= border;

                    ipoint2 boxc = box.center();

                    os << boxc.x + 1 << ',' << boxc.y + 1 << ",1," << d->m + 1 << ',';                
                    os << (box.ll.x + 1) << ',' << (box.ll.y + 1) << ',' << (box.ur.x + 1) << ',' << (box.ur.y + 1) << ',';

                    if (librectangles.empty()) { // || d->p.is_zero()) {
                        os << "0,0,0,0,";
                    } else {
                        irectangle2 r = librectangles[d->m];

                        r += ipoint2(xpos, ypos);
                        os << (r.ll.x + 1) << ',' << (r.ll.y + 1) << ',' << (r.ur.x + 1) << ',' << (r.ur.y + 1) << ',';
                    }
                }
                os << d->r(R_RESPONSE) << ',' << d->r(G_RESPONSE) << ',' << d->r(RR_RESPONSE) << '\n';
            }
        } while ((n = d->next) != nullptr);
    }
    os.close();

    // save links
    if (z > 0) {
        string name2;

        if ((save_mode & VVE_SIMPLE_FILENAMES) != 0) {
			stringstream scale_num(get_extension(fname, "_").substr(1));			
			
			if ((save_mode & VVE_MATLAB_SCALES) != 0) {
				int num;
				scale_num >> num;
				num++;
				
				stringstream matlab_scale_num;
				matlab_scale_num << num ;
				scale_num.str(matlab_scale_num.str());
			}
			
			name = dir + change_extension(fname, "", "_") + "_" + (z + 1) + "_" + scale_num.str() + "_links.txt";
			name2 = dir + change_extension(fname, "", "_") + "_" + (z + 1) + "_" + scale_num.str() + "_subparts.txt";
           //name = dir + fname + "_" + (z + 1) + "_links.txt";
           //name2 = dir + fname + "_" + (z + 1) + "_subparts.txt";
        } else {
           name = dir + nameint + "_image_" + (z + 1) + "_1_links.txt";
           name2 = dir + nameint + "_image_" + (z + 1) + "_1_subparts.txt";
        }

        ofstream os2;

        os.open(name.c_str(), ios::trunc);
        if (os.fail()) return;
        os2.open(name2.c_str(), ios::trunc);
        if (os2.fail()) return;
    
        int edgename = atom("toPrevLayer").get_index();
        for (int i = 0; i < count; ++i) {
            node* n = nodes[i];
            
            do {
                if (!n->is_attr_set(NODE_DELETED_ATTR)) {
                    bool first = true;

                    foreach_neighbor(n, edgename, iter) {
                        node* m = neighbor_node(iter);
                        edge_data_name* ed = (edge_data_name*)neighbor_edge_data(iter);

                        map<node*, int>::const_iterator j = prevlayermap.find(m);
                        if (j != prevlayermap.end()) {
                            if (first) first = false; else { os << ','; os2 << ','; }
                            os << j->second;
                            if (ed == nullptr) os2 << 0; else os2 << ed->index;
                        }
                    }
                    os << '\n';
                    os2 << '\n';
                }
            } while ((n = ((layer1_data*)n->data)->next) != nullptr);
        }
        os.close();
        os2.close();
    }
    //layermap.clear();
    save_visview(dir, nameint, fname, library, z + 1, zmax, layermap, save_mode);
}

void layer1_result::save_visview(const string& dir, const string& nameint, const string& fname, 
    part_lib* library, int z, int save_mode)
{
    map<node*,int> layermap;

    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);

    if (library != nullptr) 
        part_data::set_contractions(library->contractions, MAX_LAYER_NUMBER);

    save_visview(dir, nameint, fname, library, 0, z, layermap, save_mode);
}

int layer1_result::get_closest_node_distance2(vector<node*>& nodes, const point2<int>& p)
{
    if (nodes.size() == 0) return 0;
    
    vector<node*>::iterator iter = nodes.begin();
    layer1_data* nd = (layer1_data*)(*iter)->data;
    int dist2, mindist2 = p.distance2(nd->x, nd->y);

    while (++iter != nodes.end()) {
        nd = (layer1_data*)(*iter)->data;
        dist2 = p.distance2(nd->x, nd->y);
        if (mindist2 > dist2) mindist2 = dist2;
    } 
    return mindist2;
}

node* layer1_result::find_node_in_region(const point2<int>& p, int z, int deltax, int deltay, 
    int m, unsigned attr, double valthresh /* = 0.0 */)
{
    if (!grid(z)) init_grid(z);

    int minx = std::max<int>(0, p.x - deltax);
    int maxx = std::min<int>(x_size(z), p.x + deltax + 1);
    int miny = std::max<int>(0, p.y - deltay);
    int maxy = std::min<int>(y_size(z), p.y + deltay + 1);

    for (int i = minx; i < maxx; ++i) {
        for (int j = miny; j < maxy; ++j) {
            node* n = node_at(i, j, z);

            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if (n->is_attr_set(attr) && nd->m == m && nd->val() >= valthresh) return n;
                n = nd->next;
            }
        }
    }
    return nullptr;           
}

node* layer1_result::find_node_in_region(node* n, int deltax, int deltay, 
    int m, unsigned attr, double valthresh /* = 0.0 */)
{
    layer1_data* nd = (layer1_data*)n->data;

    return find_node_in_region(ipoint2(nd->x, nd->y), nd->z, deltax, deltay, m, attr, valthresh);
}

node* layer1_result::m_node_at(int x, int y, int z, int m)
{
    if (x < 0 || x >= x_size(z) || y < 0 || y >= y_size(z)) return nullptr;
    node* n = node_at(x, y, z);

    while (n != nullptr) {
        layer1_data* nd = (layer1_data*)n->data;

        if (nd->m == m) return n;
        n = nd->next;
    }
    return nullptr;
}

void layer1_result::find_nodes_in_region(vector<node*>& result, int x, int y, int z, int dx, int dy,
    int m, unsigned attr, double valthresh /* = 0.0 */)
{
    result.clear();

    if (!grid(z)) init_grid(z);

    int minx = std::max<int>(0, x - dx);
    int maxx = std::min<int>(x_size(z), x + dx + 1);
    int miny = std::max<int>(0, y - dy);
    int maxy = std::min<int>(y_size(z), y + dy + 1);

    for (int i = minx; i < maxx; ++i) {
        for (int j = miny; j < maxy; ++j) {
            node* n = node_at(i, j, z);

            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if (n->is_attr_set(attr) && nd->m == m && nd->val() >= valthresh) 
                    result.push_back(n);
                n = nd->next;
            }
        }
    }
}

// Return matrix with types of projections of nodes along edges with edgename
// (in fact: region of radius r around these points).
void layer1_result::layer_projection_matrix(matrix<lpr_matrix_type>& result, const vector<node*>& nodes, int edgename, int r)
{
    bool init = true;
    int xsize, ysize;

    for (vector<node*>::const_iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        int m = ((layer1_data*)n->data)->m;

        foreach_neighbor(n, edgename, niter) {
            layer1_data* d = (layer1_data*)neighbor_node_data(niter);
    
            if (init) { 
                xsize = x_size(d->z); ysize = y_size(d->z);
                result.resize(xsize, ysize); 
                init = false; 
            }

            int minx = std::max<int>(0, d->x - r);
            int maxx = std::min<int>(xsize, d->x + r + 1);
            int miny = std::max<int>(0, d->y - r);
            int maxy = std::min<int>(ysize, d->y + r + 1);

            for (int i = minx; i < maxx; ++i) {
                for (int j = miny; j < maxy; ++j) 
                    result(i, j).add_m(m);
            }

        }
    }
}

double layer1_result::convolve(int x, int y, matrix<double>& m, int mt, int z)
{
    if (!grid(z)) init_grid(z);

    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    double result = 0.0;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return 0.0;
    for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {
        for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
            node* n = node_at(ii, jj, z);

            if (n != nullptr /* && n->is_attr_set(NODE_REDUCED_ATTR )*/) {
                layer1_data* nd = (layer1_data*)n->data;
                if (nd->m == mt) result += nd->val()*m(i, j); 
            }
        }
    }
    return result;
}

pair<node*, double> layer1_result::convolve_max(int x, int y, matrix<double>& m, int mt, double thresh, int z)
{
    if (!grid(z)) init_grid(z);

    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    double r, result = 0.0;
    node* nresult = nullptr;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return pair<node*, double>(nullptr, 0.0);
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) { 
            node* n = node_at(ii, jj, z);

            if (n != nullptr /* && n->is_attr_set(NODE_REDUCED_ATTR )*/) {
                layer1_data* nd = (layer1_data*)n->data;
                double rr = nd->r(R_RESPONSE);

                if (nd->m == mt && rr >= thresh) {
                    r = rr*m(i, j); 
                    if (r > result) { result = r; nresult = n; }
                }
            }
        }
    }
    return pair<node*, double>(nresult, result);
}

double layer1_result::convolve_all(int x, int y, matrix<double>& m, int mt, double thresh, int z)
{
    if (!grid(z)) init_grid(z);

    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    double result = 0.0;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return 0.0;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {        
            node* n = node_at(ii, jj, z);

            if (n != nullptr /*&& n->is_attr_set(NODE_REDUCED_ATTR)*/) {
                layer1_data* nd;

                do {
                    nd = (layer1_data*)n->data;
                    double rr = nd->r(R_RESPONSE);

                    if (nd->m == mt && rr >= thresh) { result += rr*m(i, j); break; }
                } while ((n = nd->next) != nullptr);
            }
        }
    }
    return result;
}

pair<node*, double> layer1_result::convolve_max_all(int x, int y, matrix<double>& m, int mt, double thresh, int z)
{
    if (!grid(z)) init_grid(z);

    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    double r, result = 0.0;
    node* nresult = nullptr;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return pair<node*, double>(nullptr, 0.0);
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {        
            node* n = node_at(ii, jj, z);

            if (n != nullptr /*&& n->is_attr_set(NODE_REDUCED_ATTR)*/) {
                layer1_data* nd;

                do {
                    nd = (layer1_data*)n->data;
                    double rr = nd->r(G_RESPONSE);
    
                    if (nd->m == mt && rr >= thresh) { 
                        r = rr*m(i, j); 
                        if (r > result) { result = r; nresult = n; }
                        break; 
                    }
                } while ((n = nd->next) != nullptr);
            }
        }
    }
    return pair<node*, double>(nresult, result);
}

double layer1_result::schur_product_max(int x, int y, const matrix<double>& m, int type, int z, 
    schur_product_function* f, double thresh, vector<node*>& res)
{
    if (!grid(z)) init_grid(z);

    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    double result = 0.0;
    double r;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return 0.0;
    for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {
        for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
            node* n = node_at(ii, jj, z);
            if (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if (nd->m == type) {
                    r = (*f)(nd)*m(i, j);
                    if (r > thresh) { 
                        res.push_back(n); 
                        if (r > result) result = r;
                    }
                }
            }
        }
    }
    return result;
}

double layer1_result::schur_product_max(int x, int y, const matrix<double>& m, int type, int z, 
        schur_product_function* f, double thresh, double typethresh, bin_function<int,int,double>& typecompare, 
        vector<node*>& res)
{
    if (!grid(z)) init_grid(z);

    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    double result = 0.0;
    double r;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return 0.0;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {        
            node* n = node_at(ii, jj, z);
            if (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;
                double typecmp = typecompare(nd->m, type);

                if (typecmp >= typethresh) {
                    r = (*f)(nd)*m(i, j);
                    if (r > thresh) { 
                        res.push_back(n); 
                        if (r > result) result = r;
                    }
                }
            }
        }
    }
    return result;
}

dnpair layer1_result::schur_product_max_all(int xo, int yo, int dx, int dy, int index, matrix<double>& m, const map<int, double>& appmap, 
        int z, double c, schur_product_function* f, double thresh, sp_result_t& res)
{
    if (!grid(z)) init_grid(z);

    int x = (int)((xo + dx)*c), y = (int)((yo + dy)*c);
    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    dnpair result(0.0, nullptr);
    double r;
    map<int, double>::const_iterator appiter;

	if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return result;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {
            node* n = node_at(ii, jj, z);
            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if ((appiter = appmap.find(nd->m)) != appmap.end()) {
                    r = (*f)(nd)*m(i, j)*appiter->second;

					if (r > thresh) { 
                        res.push_back(sp_result_data_t(n, index, r)); 
                        if (r > result.first) { result.first = r; result.second = n; }
                    }
                }
                n = nd->next;
            }
        }
    }
    return result;
}

dnpair layer1_result::schur_product_max_all_new(sp_result_t& res, int xo, int yo, int z, double c,
    schur_product_params& spd)
{
    if (!grid(z)) init_grid(z);

    matrix<double>& m = spd.pdata->distr;
    int x = (int)((xo + spd.pdata->x)*c), y = (int)((yo + spd.pdata->y)*c);
    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    dnpair result(0, nullptr);
    double r;
    map<int, pair<vector<int>, double> >::const_iterator appiter;

	if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return result;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {
            node* n = node_at(ii, jj, z);
            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if ((appiter = spd.pdata->app.find(nd->m)) != spd.pdata->app.end()) {
					r = (*spd.f)(nd)*m(i, j)*appiter->second.second;

					if (r > spd.thresh) {
                        //double benergy, scdistance;

                        //if (scmap.empty() || check_geometry(benergy, scdistance, scmap, n, 
                        //        spd.pdata->geo, spd.center, spd.dthresh, spd.scthresh)) {
                        res.push_back(sp_result_data_t(n, spd.pdata->index, r)); 
                        if (r > result.first) { result.first = r; result.second = n; }
                    }
                }
                n = nd->next;
            }
        }
    }
    return result;
}

dnpair layer1_result::schur_product_max_all_ol(int x, int y, int dx, int dy, int index, matrix<double>& m, 
		const map<int, pair<vector<int>, double> >& appmap, 
        int z, schur_product_function* f, double thresh, sp_result_t& res)
{
    if (!grid(z)) init_grid(z);

    //int x = (int)((xo + dx)*c), y = (int)((yo + dy)*c);
    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    dnpair result(0.0, nullptr);
    double r;
    map<int, pair<vector<int>, double> >::const_iterator appiter;

	if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return result;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {
            node* n = node_at(ii, jj, z);
            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if ((appiter = appmap.find(nd->m)) != appmap.end()) {
					r = (*f)(nd)*m(i, j)*appiter->second.second;

					if (r > thresh) { 
                        res.push_back(sp_result_data_t(n, index, r)); 
                        if (r > result.first) { result.first = r; result.second = n; }
                    }
                }
                n = nd->next;
            }
        }
    }
    return result;
}

dnpair layer1_result::schur_product_max_all_best(int xo, int yo, int dx, int dy, int index, matrix<double>& m, 
        const map<int, double>& appmap, int z, double c, schur_product_function* f, double thresh, sp_result_t& res)
{
    if (!grid(z)) init_grid(z);

    int x = (int)((xo + dx)*c), y = (int)((yo + dy)*c);
    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    dnpair result(0.0, nullptr);
    double r;
    map<int, double>::const_iterator appiter;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return result;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {        
            node* n = node_at(ii, jj, z);
            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if ((appiter = appmap.find(nd->m)) != appmap.end()) {
                    r = (*f)(nd)*m(i, j)*appiter->second;
                    if (r > thresh && r > result.first) { result.first = r; result.second = n; }
                }
                n = nd->next;
            }
        }
    }
    if (result.second != nullptr) res.push_back(sp_result_data_t(result.second, index, result.first));
    return result;
}

dnpair layer1_result::schur_product_max_all(int xo, int yo, int dx, int dy, int index, matrix<double>& m, int type, int z, double c,
        schur_product_function* f, double thresh, double typethresh, bin_function<int,int,double>& typecompare, sp_result_t& res)
{
    if (!grid(z)) init_grid(z);

    int x = (int)((xo + dx)*c), y = (int)((yo + dy)*c);
    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    dnpair result(0.0, nullptr);
    double r;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return result;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {        
            node* n = node_at(ii, jj, z);
            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;
                double typecmp = typecompare(nd->m, type);

                if (typecmp >= typethresh) {
                    r = (*f)(nd)*m(i, j);
                    if (r > thresh) { 
                        res.push_back(sp_result_data_t(n, index, r)); 
                        if (r > result.first) { result.first = r; result.second = n; }
                    }
                }
                n = nd->next;
            }
        }
    }
    return result;
}

dnpair layer1_result::schur_product_max_all_best(int xo, int yo, int dx, int dy, int index, matrix<double>& m, int type, int z, double c,
        schur_product_function* f, double thresh, double typethresh, bin_function<int,int,double>& typecompare, sp_result_t& res)
{
    if (!grid(z)) init_grid(z);

    int x = (int)((xo + dx)*c), y = (int)((yo + dy)*c);
    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    dnpair result(0.0, nullptr);
    double r;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return result;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {        
            node* n = node_at(ii, jj, z);
            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;
                double typecmp = typecompare(nd->m, type);

                if (typecmp >= typethresh) {
                    r = (*f)(nd)*m(i, j);
                    if (r > thresh && r > result.first) { result.first = r; result.second = n; }
                }
                n = nd->next;
            }
        }
    }
    if (result.second != nullptr) res.push_back(sp_result_data_t(result.second, index, result.first)); 
    return result;
}

dnpair layer1_result::schur_product_max_all(int xo, int yo, int dx, int dy, int index, matrix<double>& m, int type, int z, double c,
        schur_product_function* f, double thresh, sp_result_t& res)
{
    if (!grid(z)) init_grid(z);

    int x = (int)((xo + dx)*c), y = (int)((yo + dy)*c);
    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    dnpair result(0.0, nullptr);
    double r;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return result;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {        
            node* n = node_at(ii, jj, z);
            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if (nd->m == type) {
                    r = (*f)(nd)*m(i, j);
                    if (r > thresh) { 
                        res.push_back(sp_result_data_t(n, index, r)); 
                        if (r > result.first) { result.first = r; result.second = n; }
                    }
                }
                n = nd->next;
            }
        }
    }
    return result;
}

double layer1_result::schur_product_max_all_dir(int xo, int yo, int dx, int dy, int index, double minfac, double maxfac,
    matrix<double>& m, const set<int>& typeset, int z, schur_product_function* f, double thresh, sp_result_t& res, double& factor)
{
    if (!grid(z)) init_grid(z);

    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    double r;
    double newfac = 1.0;
    ipoint2 P(int_round(minfac*dx), int_round(minfac*dy));
    ipoint2 Q(int_round(maxfac*dx), int_round(maxfac*dy));
    int steps = (P - Q).abs().maximum();
    sp_result_t res1, res2;
    sp_result_t* tmpres = &res1, * bestres = &res2;
    double bestresult = 0.0;

    if (steps != 0) newfac = (maxfac - minfac)/steps;

    for (int step = 0; step <= steps; ++step) {
        double tmpfactor = minfac + step*newfac;
        int x = xo + int_round(tmpfactor*dx);
        int y = yo + int_round(tmpfactor*dy);
        double tmpresult = 0.0;

        tmpres->clear();
        if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) 
            break;
		for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
			for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {
                node* n = node_at(ii, jj, z);
                while (n != nullptr) {
                    layer1_data* nd = (layer1_data*)n->data;

                    if (typeset.find(nd->m) != typeset.end()) {
                        r = (*f)(nd)*m(i, j);
                        if (r > thresh) { 
                            tmpres->push_back(sp_result_data_t(n, index, r)); 
                            if (r > tmpresult) tmpresult = r;
                        }
                    }
                    n = nd->next;
                }
            }
        }
        if (tmpresult > bestresult) { 
            swap(tmpres, bestres);
            bestresult = tmpresult;
            factor = tmpfactor;
        }
    }
    res.insert(res.end(), bestres->begin(), bestres->end());
    return bestresult;
}

double layer1_result::schur_product_max_all_dir(int xo, int yo, int dx, int dy, int index, double minfac, double maxfac,
    matrix<double>& m, int type, int z, schur_product_function* f, double thresh, sp_result_t& res, double& factor)
{
    if (!grid(z)) init_grid(z);

    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    double r;
    double newfac = 1.0;
    ipoint2 P(int_round(minfac*dx), int_round(minfac*dy));
    ipoint2 Q(int_round(maxfac*dx), int_round(maxfac*dy));
    int steps = (P - Q).abs().maximum();
    sp_result_t res1, res2;
    sp_result_t* tmpres = &res1, * bestres = &res2;
    double bestresult = 0.0;

    if (steps != 0) newfac = (maxfac - minfac)/steps;
    
    for (int step = 0; step <= steps; ++step) {
        int x = xo + int_round(minfac + step*newfac*dx);
        int y = yo + int_round(minfac + step*newfac*dy);
        double tmpresult = 0.0;

        tmpres->clear();
        if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) 
            break;
		for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
			for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {            
                node* n = node_at(ii, jj, z);
                while (n != nullptr) {
                    layer1_data* nd = (layer1_data*)n->data;

                    if (nd->m == type) {
                        r = (*f)(nd)*m(i, j);
                        if (r > thresh) { 
                            tmpres->push_back(sp_result_data_t(n, index, r)); 
                            if (r > tmpresult) tmpresult = r;
                        }
                    }
                    n = nd->next;
                }
            }
        }
        if (tmpresult > bestresult) { 
            swap(tmpres, bestres);
            bestresult = tmpresult;
        }
    }
    res.insert(res.end(), bestres->begin(), bestres->end());
    return bestresult;
}

dnpair layer1_result::schur_product_max_all_best(int xo, int yo, int dx, int dy, int index, matrix<double>& m, int type, int z, double c,
        schur_product_function* f, double thresh, sp_result_t& res)
{
    if (!grid(z)) init_grid(z);

    int x = (int)((xo + dx)*c), y = (int)((yo + dy)*c);
    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    dnpair result(0.0, nullptr);
    double r;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return result;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {        
            node* n = node_at(ii, jj, z);
            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if (nd->m == type) {
                    r = (*f)(nd)*m(i, j);
                    if (r > thresh && r > result.first) { result.first = r; result.second = n; }
                }
                n = nd->next;
            }
        }
    }
    if (result.second != nullptr) res.push_back(sp_result_data_t(result.second, index, result.first)); 
    return result;
}

dnpair layer1_result::schur_product_max_all(int x, int y, matrix<double>& m, int type, int z, 
        schur_product_function* f, double thresh)
{
    if (!grid(z)) init_grid(z);

    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    dnpair result(0.0, nullptr);
    double r;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return result;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {        
            node* n = node_at(ii, jj, z);
            while (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;				
				
                if (nd->m == type) {
                    r = (*f)(nd)*m(i, j);
                    if (r > thresh && r > result.first) { result.first = r; result.second = n; }
                }
                n = nd->next;
            }
        }
    }
    return result;
}

dnpair layer1_result::schur_product_max(int x, int y, const matrix<double>& m, int type, int z, 
    schur_product_function* f, double thresh)
{
    if (!grid(z)) init_grid(z);

    int cx = (int)m.width/2, cy = (int)m.height/2;
    int i, j, ii, jj;
    dnpair result(0.0, nullptr);
    double r;

    if (x + cx >= x_size(z) || x - cx < 0 || y + cy >= y_size(z) || y - cy < 0) return result;
	for (jj = y - cy, j = 0; j < (int)m.height; ++j, ++jj) {
		for (ii = x - cx, i = 0; i < (int)m.width; ++i, ++ii) {        
            node* n = node_at(ii, jj, z);
            if (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if (nd->m == type) {
                    r = (*f)(nd)*m(i, j);
                    if (r > result.first) { result.first = r; result.second = n; }
                }
            }
        }
    }
    return result;
}

int layer1_result::find_nodes_in_rect(sp_result_t& result, int z, const irectangle2& rect, 
    const map<int, double>& appmap, int maxn /* = INT_MAX */)
{

    if (!grid(z)) init_grid(z);

    int minj = ::max<int>(0, rect.ll.y);
    int maxi = ::min<int>(rect.ur.x, x_size(z));
    int maxj = ::min<int>(rect.ur.y, y_size(z));
    int count = 0;
    
    for (int i = max(0, rect.ll.x); i < maxi; ++i) {
        for (int j = minj; j < maxj; ++j) {
            node* n = node_at(i, j, z);

            if (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if (appmap.find(nd->m) != appmap.end()) {
                    result.push_back(sp_result_data_t(n, INT_MAX, nd->val()));
                    if (++count >= maxn) return count;
                }
            }
        }
    }
    return count;
}

node* layer1_result::find_node_in_rect(int z, const irectangle2& rect, const map<int, double>& appmap)
{
    if (!grid(z)) init_grid(z);

    int minj = ::max<int>(0, rect.ll.y);
    int maxi = ::min<int>(rect.ur.x, x_size(z));
    int maxj = ::min<int>(rect.ur.y, y_size(z));
    
    for (int i = max(0, rect.ll.x); i < maxi; ++i) {
        for (int j = minj; j < maxj; ++j) {
            node* n = node_at(i, j, z);

            if (n != nullptr) {
                layer1_data* nd = (layer1_data*)n->data;

                if (appmap.find(nd->m) != appmap.end()) return n;
            }
        }
    }
    return nullptr;
}

// We assume that rect is "ur-open".
bool layer1_result::box_empty(int z, const irectangle2& rect)
{
    if (!grid(z)) init_grid(z);

    int minj = ::max<int>(0, rect.ll.y);
    int maxi = ::min<int>(rect.ur.x, x_size(z));
    int maxj = ::min<int>(rect.ur.y, y_size(z));
    
    for (int i = max(0, rect.ll.x); i < maxi; ++i) {
        for (int j = minj; j < maxj; ++j) {
            if (node_at(i, j, z) != nullptr) return false;
        }
    }
    return true;
}

void layer1_result::copy_to(streamable* p, cloner& cl)
{
    img_graph::copy_to(p, cl);

    layer1_result* dest = (layer1_result*)p;

    dest->to_neighbor = to_neighbor;
    dest->layer1_region_threshold = layer1_region_threshold;
    dest->layer1_threshold = layer1_threshold;
    dest->layer1_3x3bound = layer1_3x3bound;
    dest->layer1_neighb_radius = layer1_neighb_radius;
    dest->attr = attr;
    dest->shape_nodes.resize(shape_nodes.size());
    for (size_t i = 0; i < shape_nodes.size(); ++i) {
        vector<node*>& src_nodes = shape_nodes[i];
        vector<node*>& dest_nodes = dest->shape_nodes[i];

        for (vector<node*>::iterator iter = src_nodes.begin(); iter != src_nodes.end(); ++iter) {
            dest_nodes.push_back((node*)cl.get_copy(*iter));
        }
    }
    dest->shape_nodes_inhib.resize(shape_nodes_inhib.size());
    for (size_t i = 0; i < shape_nodes_inhib.size(); ++i) {
        vector<node*>& src_nodes = shape_nodes_inhib[i];
        vector<node*>& dest_nodes = dest->shape_nodes_inhib[i];

        for (vector<node*>::iterator iter = src_nodes.begin(); iter != src_nodes.end(); ++iter) {
            dest_nodes.push_back((node*)cl.get_copy(*iter));
        }
    }
    dest->info = info;

}

void layer1_result::read_from_stream(istreamer& is)
{
    img_graph::read_from_stream(is);

    unsigned size;

    is.read(size);
    shape_nodes.resize(size);
    for (unsigned i = 0; i < size; ++i) is.read_streamable_p(shape_nodes[i]);
    is.read(size);
    shape_nodes_inhib.resize(size);
    for (unsigned i = 0; i < size; ++i) is.read_streamable_p(shape_nodes_inhib[i]);
    is.read(size);
    info.resize(size);
    for (unsigned i = 0; i < size; ++i) info[i].read_from_stream(is);
    
    is.read(to_neighbor);
    is.read(layer1_region_threshold);
    is.read(layer1_threshold);
    is.read(layer1_3x3bound);
    is.read(layer1_neighb_radius);
    if (is.get_version() >= 2.8) {
        is.read(original_width);
        is.read(border);
    } else {
        original_width = 0;
        border = 0;
    }
}

void layer1_result::write_to_stream(ostreamer& os)
{
    img_graph::write_to_stream(os);
    
    size_t size;

    size = shape_nodes.size();
    os.write((unsigned)size);
    for (size_t i = 0; i < size; ++i) os.write_streamable_p(shape_nodes[i]);
    size = shape_nodes_inhib.size();
    os.write((unsigned)size);
    for (size_t i = 0; i < size; ++i) os.write_streamable_p(shape_nodes_inhib[i]);
    size = info.size();
    os.write((unsigned)size);
    for (size_t i = 0; i < size; ++i) info[i].write_to_stream(os);
    
    os.write(to_neighbor);
    os.write(layer1_region_threshold);
    os.write(layer1_threshold);
    os.write(layer1_3x3bound);
    os.write(layer1_neighb_radius);
    os.write(original_width);
    os.write(border);
}

void layer1_result::write_vgr_label(ostream& os, node* n, int count)
{
    int name = atom("toPrevLayer").get_index();

    os << '\"' << *((layer1_data*)n->data) << " - " << name << " - " <<
        ::sqrt((double)get_neighbor_set_extents2(n, name)) << " (" << count-1 << ")\"";
}

void layer1_result::write_clu(ostream& os, const vector<node*> nodes)
{
    vector<node*>::const_iterator iter;

    os << "*Vertices " << nodes.size() << endl;
    for (iter = nodes.begin(); iter != nodes.end(); ++iter) {
        layer1_data* d = (layer1_data*)(*iter)->data;
        os << d->m << endl;
    }
}

// Write suitable for importing into Mathematica:
// {list of types, list of edges, list of coordinates}
void layer1_result::write_mma2(ostream& os, const vector<node*>& nodes, int edgename, bool alltypes /* = false */ )
{
    map<node*, int> node_enum;
    int size = (int)nodes.size();

    os << "{{";
    for (int i = 0; i < size; ++i) {
        node* n = nodes[i];
        layer1_data* nd = (layer1_data*)n->data;

        if (i > 0) os << ',';
        if (!alltypes) os << nd->m;
        else {
            int ti = 0;

            os << '{';
            while (true) {
                if (ti > 0) os << ',';
                os << '{' << nd->m << ',' << nd->val() << '}';
                if ((n = nd->next) == nullptr) break;
                nd = (layer1_data*)n->data;
                ti++;
            }
            os << '}';
        }
        node_enum.insert(pair<node*, int>(n, i));
    }
    os << "}," << endl << '{';

    bool first = true;

    for (int i = 0; i < size; ++i) {
        node* n = nodes[i];

        foreach_neighbor(n, edgename, iter) {
            node* m = neighbor_node(iter);
            map<node*, int>::iterator f = node_enum.find(m);

            if (f == node_enum.end()) continue;

            if (!first) os << ','; else first = false;
            os << '{' << i + 1 << ',' << f->second + 1 << '}';
        }
    }
    os << "}," << endl << '{';
    for (int i = 0; i < size; ++i) {
        node* n = nodes[i];
        layer1_data* nd = (layer1_data*)n->data;

        if (i > 0) os << ',';
        os << '{' << nd->x << ',' << nd->y << '}';
    }
    os << "}}" << endl;
}

// Write the *whole* structure in Mathematica format.
// {{{{layer1_p1_x, layer1_p1_y, layer1_p1_type, layer1_p1_value}, {}}, ...},
//  {{{layer1_p2_x, layer1_p2_y, layer1_p2_type, layer1_p2_value}, {to_prevlayer, ...}}, ...}
//  ...
// }
//
void layer1_result::write_mma2(ostream& os)
{
    os << '{';
    for (int i = 0; i < (int)shape_nodes.size(); ++i) {
        vector<node*>& nodes = shape_nodes[i];
        map<node*, int> nodemap;

        if (i > 0) {
            os << ',' << endl;
            vector<node*>& nodesp = shape_nodes[i - 1];
            int count = 0;
            
            for (vector<node*>::iterator iter = nodesp.begin(); iter != nodesp.end(); ++iter) 
                nodemap.insert(pair<node*, int>(*iter, ++count));
        }

        os << '{';
        for (int j = 0; j < (int)nodes.size(); ++j) {
            node* n = nodes[j];
            layer1_data* nd = (layer1_data*)n->data;

            if (j > 0) os << ',' << endl;

            os << "{{"; 
            os << nd->x << ',' << nd->y << ',' << nd->m << ',' << nd->val() << "},{";
            if (!nodemap.empty()) {
                int edge = atom("toPrevLayer").get_index();
                int ncount = 0;
                foreach_neighbor(n, edge, niter) {
                    map<node*,int>::iterator f = nodemap.find(neighbor_node(niter));
                    
                    if (f != nodemap.end()) {
                        if (ncount > 0) os << ',';
                        os << f->second;
                        ++ncount;
                    }
                }
            }
            os << "}}";
        }
        os << '}';
    }
    os << '}' << endl;
}

// Write the layer
void layer1_result::write_mma3(ostream& os, int m)
{
    if (m < 0 || m >= layer_count()) {
        os << "{}" << endl;
        return;
    }

    vector<node*>& s_nodes = shape_nodes[m];
    bool first = true;

    os << '{';
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        
        if (first) first = false; else os << ',';
        os << '{';
        while (n != nullptr) {
            layer1_data* nd = (layer1_data*)n->data;

            os << '{' << nd->x << ',' << nd->y << ',' << nd->m << ',' << nd->val() << '}';
            n = nd->next;
            if (n != nullptr) os << ',';
        }
        os << '}' << endl;
    }
    os << '}' << endl;
}


// Write the reconstruction of each node on slayer to dlayer, dlayer <= slayer,
// in Mathematica format.
// Format of the output is:
//   {{x_coo_of_node1, y_coo_of_node1, type_of_node1, value_of_node1, attr_of_node1, 
//       {desc_type, {{rec_node1_type, rec_node1_val, rec_node1_attr}, 
//                    {rec_node2_type, rec_node2_val, rec_node2_attr}, ...}}},
//     ...
//   }
void layer1_result::write_mma4(ostream& os, int slayer, int dlayer) 
{
    if (slayer < 0 || slayer >= layer_count() || dlayer < 0 || dlayer > slayer) {
        os << "{}" << endl;
        return;
    }
    
    vector<node*>& s_nodes = shape_nodes[slayer];
    layer_predicate pred(dlayer);
    int name = atom("toPrevLayer");

    os << '{';
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;
        multimap<int, node*> cmap;
        multimap<int, node*>::iterator cmapiter;
        bool first1 = true;

        if (iter != s_nodes.begin()) os << ',';
        os << '{' << nd->x << ',' << nd->y << ',' << nd->m << ',' << nd->val() << ',' << n->attr << ',';
        os << '{';

        foreach_neighbor(n, name, niter) {
            edge_data_name* ed = (edge_data_name*)neighbor_edge_data(niter);
            if (ed != nullptr) cmap.insert(pair<int, node*>(ed->index, neighbor_node(niter)));
        }

        if (cmap.empty()) {
            foreach_neighbor(n, name, niter) {
                cmap.insert(pair<int, node*>(0, neighbor_node(niter)));
            }
        }

        for (cmapiter = cmap.begin(); cmapiter != cmap.end(); ) {
            vector<node*> srcnodes;
            set<node*> recnodes;
            bool first2 = true;

            if (first1) first1 = false; else os << ',';
            os << '{' << cmapiter->first << ',' << '{';

            for (int p = cmapiter->first; cmapiter != cmap.end() && cmapiter->first == p; ++cmapiter) {
                srcnodes.push_back(cmapiter->second);
            }
            recurse(srcnodes, name, pred, recnodes);

            for (set<node*>::iterator riter = recnodes.begin(); riter != recnodes.end(); ++riter) {
                node* rn = *riter;
                layer1_data* rnd = (layer1_data*)rn->data;

                if (rnd->z == dlayer) {
                    if (first2) first2 = false; else os << ',';
                    os << '{' << rnd->x << ',' << rnd->y << ',' << rnd->m << ',' << rnd->val() << ',' << rn->attr << '}';
                }
            }
            os << '}' << '}' << endl;
        }
        os << '}' << '}' << endl;
    }
    os << '}' << endl;
}

// Write the layer
void layer1_result::write_matlab(ostream& os, int m)
{
    if (m < 0 || m >= layer_count()) {
        os << "{}" << endl;
        return;
    }

    vector<node*>& s_nodes = shape_nodes[m];
    bool first = true;

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        
        while (n != nullptr) {
            layer1_data* nd = (layer1_data*)n->data;

            os << nd->x << ' ' << nd->y << ' ' << nd->m << ' ' << nd->val() << endl;
            n = nd->next;
        }
    }
}

void layer1_result::get_edge_info(set<itriple>& info)
{
    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;
        
        forall_neighbors(n, niter) {
            int name = neighbor_index(niter);
            layer1_data* nnd = (layer1_data*)neighbor_node_data(niter);

            info.insert(itriple(nd->z, nnd->z, name));
        }
    }
}
void layer1_result::write_node_info(ostream& os, const vector<node*>& nodes)
{
    bool firstn = true;

    os << '{';
    for (vector<node*>::const_iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;
        bool firstnb = true;

        if (firstn) firstn = false; else os << ',' << endl;
        os << '{' << nd->m << ',' << nd->val() << ',' << '{';
        forall_neighbors(n, niter) {
            if (firstnb) firstnb = false; else os << ',';
            os << '\"' << atom::get_name(neighbor_index(niter)) << '\"';
        }
        os << '}' << '}';
    }
    os << '}';
}

void layer1_result::write_edge_info(ostream& os, const set<itriple>& info)
{
    for (set<itriple>::const_iterator siter = info.begin(); siter != info.end(); ++siter) {
        const itriple& t = *siter;

        os << "layer " << t.first << " --> " << t.second << "; name \"" << 
            atom::get_name(t.third) << '\"' << endl;
    }
}

void layer1_result::write_edge_info(ostream& os)
{
    set<itriple> result;

    get_edge_info(result);
    write_edge_info(os, result);
}

void layer1_result::type_statistics(vector<int>& stat, const vector<node*>& nodes, bool next)
{
    for (vector<node*>::const_iterator i = nodes.begin(); i != nodes.end(); ++i) {
        node* n = *i;
        layer1_data* nd;
        
        do {
            nd = (layer1_data*)n->data;
            if (nd->m >= (int)stat.size()) stat.resize(nd->m + 1, 0);
            ++stat[nd->m];
        } while (next && (n = nd->next) != nullptr);
    }
}

void layer1_result::type_statistics(vector<int>& stat, int layer, bool next)
{
    if (layer >= 0 && layer < layer_count())
        type_statistics(stat, shape_nodes[layer], next);
}

void layer1_result::get_covering_statistics(const vector<node*>& nodes, const vector<int>& parts, vector<int>& stat)
{
    stat.clear();
    vector<int>::const_iterator piter = parts.begin();

    while (piter != parts.end()) {
        vector<int> sparts(parts.begin(), ++piter);
        vector<node*> srcnodes;
        set<node*> endnodes;

        sort(sparts.begin(), sparts.end());
        for (vector<node*>::const_iterator niter = nodes.begin(); niter != nodes.end(); ++niter) {
            node* n = *niter;
            layer1_data* nd;
            
            do {
                nd = (layer1_data*)n->data;
                if (binary_search(sparts.begin(), sparts.end(), nd->m)) srcnodes.push_back(n);
            } while ((n = nd->next) != nullptr);
        }
        recurse(srcnodes, atom("toPrevLayer").get_index(), endnodes);
        stat.push_back((int)endnodes.size());
    }
}

double layer1_result::size_factor(int srclayer, int destlayer)
{
    return (double)x_size(destlayer)/x_size(srclayer);
}

void layer1_result::adjust_rectangle(irectangle2& rect, int srclayer, int destlayer)
{
    rect.resize((double)x_size(destlayer)/x_size(srclayer));
}

void layer1_result::get_box(irectangle2& box, int z)
{
    box.invalidate();

    if (z < 0 || z > max_layer_index()) return;

    vector<node*>& nodes = shape_nodes[z];

    for (vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;
        box.eat(nd->x, nd->y);
    }
}

ipoint2 layer1_result::get_center_of_mass(int z)
{
    if (z < 0 || z > max_layer_index()) return ipoint2();

    vector<node*>& nodes = shape_nodes[z];
    int xsum = 0, ysum = 0;
    int count = 0;

    for (vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;
        xsum += nd->x; ysum += nd->y;
        ++count;
    }
    return ipoint2(int_round((double)xsum/count), int_round((double)ysum/count));
}

// result[i] is area(BB_i ^ gtr)/area(BB_i u gtr) for all nodes of type found in 'types'.
void layer1_result::check_with_groundtruth(vector<double>& result, const list<irectangle2>& gtrs, 
        int z, const set<int>& types, double inhibit)
{
    list<box_data_t> boxes;

    result.clear();
    get_boxes(boxes, z, types, inhibit);
    for (list<box_data_t>::iterator iter = boxes.begin(); iter != boxes.end(); ++iter) {
        irectangle2& r = iter->box;
        double bestf = 0.0;

        for (list<irectangle2>::const_iterator riter = gtrs.begin(); riter != gtrs.end(); ++riter) {
            const irectangle2& gtr = *riter;
            double f = gtr.invalid() ? 0.0 : (double)(r.intersection(gtr).area())/r.union_area(gtr);

            if (f > bestf) bestf = f;
        }
        result.push_back(bestf);
    }
}

// Keeps only nodes with area(projection_bb ^ gtr)/area(projection_bb u gtr) >= thresh for sone
// rectangle in gtr. Returns the nodes for which this is not true.
vector<node*> layer1_result::filter_with_groundtruth(vector<node*>& nodes, const list<irectangle2>& gtr, double thresh)
{
    vector<node*> result;

    if (gtr.empty()) return nodes;

    vector<node*>::iterator iter = nodes.begin();

    while (iter != nodes.end()) {
        irectangle2 box = get_box(*iter);

        if (box.invalid()) {
            result.push_back(*iter);
            iter = nodes.erase(iter);
        } else {
            bool inside = false;

            for (list<irectangle2>::const_iterator biter = gtr.begin(); !inside && biter != gtr.end(); ++biter) {
                double r = (double)box.intersection(*biter).area()/box.union_area(*biter);
                inside = r >= thresh;
            }
            if (inside) ++iter; 
            else {
                result.push_back(*iter);
                iter = nodes.erase(iter);
            }
        }
    }
    return result;
}

vector<node*> layer1_result::filter_with_groundtruth(vector<node*>& nodes, const list<irectangle2>& gtr, 
	double pthresh, double nthresh)
{
    vector<node*> result;

    if (gtr.empty()) return nodes;

    vector<node*>::iterator iter = nodes.begin();

    while (iter != nodes.end()) {
        irectangle2 box = get_box(*iter);

        if (box.invalid()) {
            result.push_back(*iter);
            iter = nodes.erase(iter);
        } else {
            bool inside = false;
			double minr = 1.0, maxr = 0.0;

            for (list<irectangle2>::const_iterator biter = gtr.begin(); !inside && biter != gtr.end(); ++biter) {
                double r = (double)box.intersection(*biter).area()/box.union_area(*biter);
                inside = r >= pthresh;
				if (r < minr) minr = r;
				if (r > maxr) maxr = r;
            }
            if (inside) ++iter; 
            else {
                if (maxr <= nthresh) result.push_back(*iter);
                iter = nodes.erase(iter);
            }
        }
    }
    return result;
}


void layer1_result::count_hits(map<int, pair<int, int> >& stat, const list<irectangle2>& gtrs, 
    const list<box_data_t>& boxes, double thresh)
{
    typedef map<int, pair<int, int> > stat_t;

    for (list<box_data_t>::const_iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
        const box_data_t& bd = *biter;

        for (list<irectangle2>::const_iterator iter = gtrs.begin(); iter != gtrs.end(); ++iter) {
            const irectangle2& gtr = *iter;
            double f = (double)(bd.box.intersection(gtr)).area()/bd.box.union_area(gtr);
            stat_t::iterator fiter = stat.find(bd.m);
    
            if (fiter == stat.end()) 
                stat.insert(stat_t::value_type(bd.m, (f >= thresh) ? pair<int, int>(1, 0) : pair<int, int>(0, 1)));
            else 
                if (f >= thresh) ++(fiter->second.first); else ++(fiter->second.second);
        }
    }
}

void layer1_result::count_hits(vector<bool>& hits, int& misses, vector<pair<double, bool> >& confidence,
     const list<irectangle2>& gtrs, const list<box_data_t>& boxes, double thresh)
{
    misses = 0;
    hits.resize(gtrs.size(), false);
    for (list<box_data_t>::const_iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
        const irectangle2& r = biter->box;
        double bestf = 0.0;
        int count = 0;
        
        for (list<irectangle2>::const_iterator iter = gtrs.begin(); iter != gtrs.end(); ++iter) {
            const irectangle2& gtr = *iter;
            double f = (double)(r.intersection(gtr)).area()/r.union_area(gtr);

            if (f > bestf) bestf = f;
            if (f >= thresh) hits[count] = true;
            ++count;
        }
        // we did not hit anything?
        if (bestf < thresh) ++misses;
        confidence.push_back(pair<double, bool>(biter->get_value(), bestf >= thresh));
    }
}

void layer1_result::count_hits(vector<bool>& hits, int& misses, vector<pair<double, bool> >& confidence,
     const list<pair<irectangle2, int> >& gtrs, const list<box_data_t>& boxes, double thresh)
{
    misses = 0;
    hits.resize(gtrs.size(), false);
    for (list<box_data_t>::const_iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
        const irectangle2& r = biter->box;
        double bestf = 0.0;
        int count = 0;
        
        for (list<pair<irectangle2, int> >::const_iterator iter = gtrs.begin(); iter != gtrs.end(); ++iter) {
            const pair<irectangle2, int>& gtr = *iter;

			if (biter->m == gtr.second) {
				double f = (double)(r.intersection(gtr.first)).area()/r.union_area(gtr.first);

	            if (f > bestf) bestf = f;
		        if (f >= thresh) hits[count] = true;
			}
            ++count;
        }
        // we did not hit anything?
        if (bestf < thresh) ++misses;
        confidence.push_back(pair<double, bool>(biter->get_value(), bestf >= thresh));
    }
}

void layer1_result::count_hits_inhibited(vector<bool>& hits, int& misses, vector<pair<double, bool> >& confidence,
    const list<irectangle2>& gtrs, int z, const set<int>& types, double thresh)
{
    list<box_data_t> boxes;
    get_boxes_inhibited(boxes, z, types, 0.8, 0.3);

    misses = 0;
    hits.resize(gtrs.size(), false);
    for (list<box_data_t>::iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
        irectangle2& r = biter->box;
        double bestf = 0.0;
        int count = 0;
        
        for (list<irectangle2>::const_iterator iter = gtrs.begin(); iter != gtrs.end(); ++iter) {
            const irectangle2& gtr = *iter;
            double f = (double)(r.intersection(gtr)).area()/r.union_area(gtr);

            if (f > bestf) bestf = f;
            if (f >= thresh) hits[count] = true;
            ++count;
        }
        // we did not hit anything?
        if (bestf < thresh) ++misses;
        confidence.push_back(pair<double, bool>(biter->get_value(), bestf >= thresh));
    }
}

// Get nodes on layer z with type (m) in mset.
// If mset is empty, the is no restriction on type.
// Note: it does not delete 'result'; new nodes are appended to the corresponding vector!
void layer1_result::get_nodes(vector<node*>& result, int z, const set<int>& mset)
{
    typedef map<int, vector<node*> > result_t;

    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;

        if (nd->z == z && (mset.empty() || mset.find(nd->m) != mset.end())) {
            result.push_back(*iter);
        }
    }
}

void layer1_result::get_contractions(vector<double>& contractions)
{
    contractions.clear();
    for (int i = 1; i <= max_layer_index(); ++i) 
        contractions.push_back((double)x_size(i - 1)/x_size(i));
}

// Get nodes on layers z0 <= z <= z1 with type (m) in mset.
// If mset is empty, the is no restriction on type.
// Note: it does not delete 'result'; new nodes are appended to the corresponding vector!
void layer1_result::get_nodes(vector<node*>& result, int z0, int z1, const set<int>& mset)
{
    typedef map<int, vector<node*> > result_t;

    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;

        if (z0 <= nd->z && nd->z <= z1 && (mset.empty() || mset.find(nd->m) != mset.end())) {
            result.push_back(*iter);
        }
    }
}

// Returns (mean, variance) of the responses of "top" (best R_RESPONSE) nodes.
ddpair layer1_result::top_response_distribution(int layer, int response)
{
    if (layer < 0 || layer > max_layer_index()) return ddpair(0.0, 0.0);

    online_distribution result;
    vector<node*>& s_nodes = shape_nodes[layer];

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;

        result.new_data(nd->r.get_response(response));
    }
    return ddpair(result.get_mean(), result.get_variance());
}

// Returns (mean, variance) of the responses of nodes on a specific layer.
ddpair layer1_result::response_distribution(int layer, int response)
{
    if (layer < 0 || layer > max_layer_index()) return ddpair(0.0, 0.0);

    online_distribution result;
    vector<node*>& s_nodes = shape_nodes[layer];

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;

        do {
            layer1_data* nd = (layer1_data*)n->data;
            result.new_data(nd->r.get_response(response));
            n = nd->next;
        } while (n != nullptr);
    }
    return ddpair(result.get_mean(), result.get_variance());
}

void layer1_result::copy_to(graph* dest, const map<node*, node*>& cmap)
{
    img_graph::copy_to(dest, cmap);

    layer1_result* ly1r = (layer1_result*)dest;

    ly1r->to_neighbor = to_neighbor;
    ly1r->layer1_region_threshold = layer1_region_threshold;
    ly1r->layer1_threshold = layer1_threshold;
    ly1r->layer1_3x3bound = layer1_3x3bound;      
    ly1r->layer1_neighb_radius = layer1_neighb_radius;
    ly1r->attr = attr; 

    ly1r->fill_shape_node_vectors();
}

struct response_ordering_predicate_gt {
    int response;

    response_ordering_predicate_gt(int r) : response(r) { }
    bool operator()(node* n, node* m) 
    {
        return ((layer1_data*)n->data)->r(response) > ((layer1_data*)m->data)->r(response);
    }
};

void layer1_result::sorted_nodes_at(vector<node*>& v, int x, int y, int z, int response)
{
    if (!grid(z)) init_grid(z);

    v.clear();
    node* n = node_at_safe(x, y, z);

    while (n != nullptr) {
        v.push_back(n);
        n = ((layer1_data*)n->data)->next;
    }
    
    response_ordering_predicate_gt sortp(response);

    sort(v.begin(), v.end(), sortp);
}

void layer1_result::fill_shape_node_vectors()
{
    shape_nodes.reserve(10);
    shape_nodes_inhib.reserve(10);
    info.reserve(10);
    for (list<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;

        if ((int)shape_nodes.size() <= nd->z) {
            shape_nodes.resize(nd->z + 1);
            shape_nodes_inhib.resize(nd->z + 1);
            info.resize(nd->z + 1);
        }
        
        if (n->is_attr_set(IMG_NODE_ATTR)) {
            shape_nodes[nd->z].push_back(n);
            if (n->is_attr_set(NODE_REDUCED_ATTR)) shape_nodes_inhib[nd->z].push_back(n);
        }
        
    }
    for (int i = 0; i < (int)shape_nodes.size(); ++i) {
        sort(shape_nodes[i].begin(), shape_nodes[i].end(), layer1_data::less1n);
        sort(shape_nodes_inhib[i].begin(), shape_nodes_inhib[i].begin(), layer1_data::less1n);
    }

}

// On layer z keep only the given parts.
//   z = layer index
//   parts = parts to keep
// Restrictions: layer z should be the last layer.
void layer1_result::dilute_layer(int z, const vector<int>& parts)
{
    if (z != max_layer_index()) return;

    if (!grid(z)) init_grid(z);

    vector<node*>& s_nodes = shape_nodes[z];
    vector<node*>& s_nodes_inhib = shape_nodes_inhib[z];
    vector<int> index_map(*max_element(parts.begin(), parts.end()) + 1, -1);
    int index_map_size = (int)index_map.size();
    int i = 0;

    for (vector<int>::const_iterator iter = parts.begin(); iter != parts.end(); ++iter) 
        index_map[*iter] = i++;

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;
        node** pn = &r_node_at(nd->x, nd->y, nd->z);
        node* nx = *pn;

        while (nx != nullptr) {
            layer1_data* nxd = (layer1_data*)nx->data;
            if (nxd->m < index_map_size && (i = index_map[nxd->m]) >= 0) {
                pn = &nxd->next;
                nxd->m = i;
            } else {
                *pn = nxd->next;
                nx->set_attr(NODE_DELETED_ATTR);
            }
            nx = *pn;
        }
    }
    remove_nodes(s_nodes, NODE_DELETED_ATTR);
    remove_nodes(s_nodes_inhib, NODE_DELETED_ATTR);
    graph::delete_nodes(NODE_DELETED_ATTR);
}

void layer1_result::delete_layers_geq(int z, bool set_has_next_layer_attr /* = true */)
{
    set<node*> to_delete;

    for (int l = z; l < (int)shape_nodes.size(); ++l) {
        to_delete.insert(shape_nodes[l].begin(), shape_nodes[l].end());
        shape_nodes[l].clear();
        shape_nodes_inhib[l].clear();
    }
    delete_nodes(to_delete);
    if (set_has_next_layer_attr && z > 0 && z - 1 < layer_count()) 
        clear_attr(shape_nodes[z - 1].begin(), shape_nodes[z - 1].end(), HAS_NEXT_LAYER);
}

void layer1_result::delete_nodes(const set<node*>& ns)
{
    for (int l = 0; l <= max_layer_index(); ++l) {
        vector<node*>& s_nodes = shape_nodes[l];
        vector<node*>& s_nodes_inhib = shape_nodes_inhib[l];
        vector<node*>::iterator iter;
        
        iter = s_nodes.begin();
        while (iter != s_nodes.end()) {
            if (ns.find(*iter) == ns.end()) ++iter;
            else iter = s_nodes.erase(iter);
        }

        iter = s_nodes_inhib.begin();
        while (iter != s_nodes_inhib.end()) {
            if (ns.find(*iter) == ns.end()) ++iter;
            else iter = s_nodes_inhib.erase(iter);
        }

    }
    img_graph::delete_nodes(ns);
}

void layer1_result::delete_nodes(unsigned attr)
{
    set<node*> to_delete;

    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;

        if (n->is_attr_set(HIDDEN_NODE_ATTR)) 
            to_delete.insert(n);
    }
    delete_nodes(to_delete);
}

void layer1_result::remove_grid_nodes(int x, int y, int z, unsigned attr, double valthresh)
{
    node** np = &r_node_at(x, y, z);
    node* n = *np;

    while (n != nullptr) {
        layer1_data* nd = (layer1_data*)n->data;
        if (n->is_attr_set(attr) && nd->val() < valthresh) {
            n->set_attr(NODE_DELETED_ATTR);
            *np = nd->next;
        } else {
            np = &nd->next;
        }
        n = *np;
    }
    if ((n = node_at(x, y, z)) != nullptr) n->set_attr(IMG_NODE_ATTR);
}

void layer1_result::filter_nodes(int z, const irectangle2& rect)
{
    for (list<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;

        if (nd->z == z && !rect.inside(ipoint2(nd->x, nd->y))) {
            n->set_attr(NODE_DELETED_ATTR);
            r_node_at(nd->x, nd->y, nd->z) = nullptr;  // f**** dangerous!
        }
    }
    remove_nodes(shape_nodes[z], NODE_DELETED_ATTR);
    remove_nodes(shape_nodes_inhib[z], NODE_DELETED_ATTR);
    graph::delete_nodes(NODE_DELETED_ATTR);
}


// layer1_result_struct
///////////////////////////////////////////////////////////////////////////////

void layer1_result_struct::copy_to(graph* dest, const map<node*, node*>& cmap)
{
    layer1_result::copy_to(dest, cmap);

    layer1_result_struct* destd = (layer1_result_struct*)dest;

    destd->gabor_lambda = gabor_lambda;
    destd->gabor_gamma = gabor_gamma;
    destd->gabor_bw = gabor_bw;
    destd->gabor_step = gabor_step;
}

void layer1_result_struct::get_masks(vector<img*>& masks) 
{ 
    layer1_creator_struct::get_filter_kernels(masks, gabor_lambda, gabor_gamma, gabor_bw, gabor_step);
}

void layer1_result_struct::get_real_masks(vector<img*>& masks)
{
    get_masks(masks);
    int size = (int)masks.size();
    for (int i = size/2; i < size; ++i) delete masks[i];
    masks.erase(masks.begin() + size/2 , masks.end());
}

void layer1_result_struct::get_regions(vector<matrix<bool>*>& regions, config_dictionary& cfg)
{
    vector<img*> masks;

    layer1_creator_struct::get_filter_kernels(masks,
        cfg.get_value_double("region_lambda", 7.0),
        cfg.get_value_double("region_gamma", 1.0),
        cfg.get_value_double("region_bw", 2.75),
        gabor_step); 

    size_t nmasks2 = masks.size()/2;

    regions.resize(nmasks2);
    for (size_t i = 0; i < nmasks2; ++i) 
        regions[i] = masks[i]->get_bool_matrix(layer1_region_threshold);

    for (size_t i = 0; i < masks.size(); ++i) 
        delete masks[i];
}

void layer1_result_struct::get_part_data(vector<part_data*>& pd, config_dictionary& cfg)
{
    vector<img*> masks;
    vector<img*> r_masks;
    matrix<bool>* region;

    get_masks(masks);
    layer1_creator_struct::get_filter_kernels(r_masks,
        cfg.get_value_double("region_lambda", 7.0),
        cfg.get_value_double("region_gamma", 1.0),
        cfg.get_value_double("region_bw", 2.75),
        gabor_step); 

    size_t nmasks2 = masks.size()/2;

    pd.clear();
    for (size_t i = 0; i < nmasks2; ++i) {
        region = r_masks[i]->get_bool_matrix(layer1_region_threshold);
        pd.push_back(new part_data(*masks[i], *region, 0, (int)i));
        delete region;
    }
    for (size_t i = 0; i < masks.size(); ++i) {
        delete masks[i];
        delete r_masks[i];
    }
}

void layer1_result_struct::copy_to(streamable* p, cloner& cl)
{
    layer1_result::copy_to(p, cl);

    ((layer1_result_struct*)p)->gabor_lambda = gabor_lambda;
    ((layer1_result_struct*)p)->gabor_gamma = gabor_gamma;
    ((layer1_result_struct*)p)->gabor_bw = gabor_bw;
    ((layer1_result_struct*)p)->gabor_step = gabor_step;
}

void layer1_result_struct::read_from_stream(istreamer& is)
{
    layer1_result::read_from_stream(is);

    is.read(gabor_lambda);
    is.read(gabor_gamma);
    is.read(gabor_bw);
    is.read(gabor_step);
}

void layer1_result_struct::write_to_stream(ostreamer& os)
{
    layer1_result::write_to_stream(os);
    
    os.write(gabor_lambda);
    os.write(gabor_gamma);
    os.write(gabor_bw);
    os.write(gabor_step);
}

// layer1_result_color
///////////////////////////////////////////////////////////////////////////////

void layer1_result_color::copy_to(graph* dest, const map<node*, node*>& cmap)
{
    layer1_result::copy_to(dest, cmap);

    layer1_result_color* destd = (layer1_result_color*)dest;

    destd->gabor_size = gabor_size;
    destd->n_rotations = n_rotations;
}

void layer1_result_color::get_masks(vector<img*>& masks) 
{ 
    //layer1_creator_color::get_masks(masks, 7, 6);
	layer1_creator_color::get_masks(masks, 0, 0);	// dummy args!
}

void layer1_result_color::get_real_masks(vector<img*>& masks)
{
    get_masks(masks);
    int size = (int)masks.size();
    for (int i = size/2; i < size; ++i) delete masks[i];
    masks.erase(masks.begin() + size/2 , masks.end());
}

void layer1_result_color::get_regions(vector<matrix<bool>*>& regions, config_dictionary& cfg)
{
    vector<img*> masks;

    layer1_creator_struct::get_filter_kernels(masks,
        cfg.get_value_double("struct_region_lambda", 7.0),
        cfg.get_value_double("struct_region_gamma", 1.0),
        cfg.get_value_double("struct_region_bw", 2.75),
        int_round(180.0/n_rotations)); 

    size_t nmasks2 = masks.size()/2;

    regions.resize(nmasks2);
    for (size_t i = 0; i < nmasks2; ++i) 
        regions[i] = masks[i]->get_bool_matrix(layer1_region_threshold);

    for (size_t i = 0; i < masks.size(); ++i) 
        delete masks[i];
}

void layer1_result_color::get_part_data(vector<part_data*>& pd, config_dictionary& cfg)
{
    vector<img*> masks;
    vector<img*> r_masks;
    matrix<bool>* region;

    get_masks(masks);
    layer1_creator_struct::get_filter_kernels(r_masks,
        cfg.get_value_double("region_lambda", 7.0),
        cfg.get_value_double("region_gamma", 1.0),
        cfg.get_value_double("region_bw", 2.75),
        int_round(180.0/n_rotations)); 

    size_t nmasks2 = masks.size()/2;

    pd.clear();
    for (size_t i = 0; i < nmasks2; ++i) {
        region = r_masks[i]->get_bool_matrix(layer1_region_threshold);
        pd.push_back(new part_data(*masks[i], *region, 0, (int)i));
        delete region;
    }
    for (size_t i = 0; i < masks.size(); ++i) {
        delete masks[i];
        delete r_masks[i];
    }
}

void layer1_result_color::copy_to(streamable* p, cloner& cl)
{
    layer1_result::copy_to(p, cl);

    ((layer1_result_color*)p)->gabor_size = gabor_size;
    ((layer1_result_color*)p)->n_rotations = n_rotations;
}

void layer1_result_color::read_from_stream(istreamer& is)
{
    layer1_result::read_from_stream(is);

    is.read(gabor_size);
    is.read(n_rotations);
}

void layer1_result_color::write_to_stream(ostreamer& os)
{
    layer1_result::write_to_stream(os);
    
    os.write(gabor_size);
    os.write(n_rotations);
}

// layer1_result_app
///////////////////////////////////////////////////////////////////////////////

void layer1_result_app::get_masks(vector<img*>& masks) 
{ 
    layer1_creator_app::get_filter_kernels(masks, gabor_lambda, gabor_gamma, gabor_bw, gabor_step); 
}

void layer1_result_app::copy_to(streamable* p, cloner& cl)
{
    layer1_result::copy_to(p, cl);

    ((layer1_result_app*)p)->gabor_lambda = gabor_lambda;
    ((layer1_result_app*)p)->gabor_gamma = gabor_gamma;
    ((layer1_result_app*)p)->gabor_bw = gabor_bw;
    ((layer1_result_app*)p)->gabor_step = gabor_step;
}

void layer1_result_app::read_from_stream(istreamer& is)
{
    layer1_result::read_from_stream(is);

    is.read(gabor_lambda);
    is.read(gabor_gamma);
    is.read(gabor_bw);
    is.read(gabor_step);
}

void layer1_result_app::write_to_stream(ostreamer& os)
{
    layer1_result::write_to_stream(os);
    
    os.write(gabor_lambda);
    os.write(gabor_gamma);
    os.write(gabor_bw);
    os.write(gabor_step);
}

// layer1_result_dog
///////////////////////////////////////////////////////////////////////////////

void layer1_result_dog::get_masks(vector<img*>& masks) 
{ 
    layer1_creator_dog::get_filter_kernels(masks, sigma_inner, sigma_outer, mask_size_factor); 
}

void layer1_result_dog::copy_to(streamable* p, cloner& cl)
{
    layer1_result::copy_to(p, cl);

    ((layer1_result_dog*)p)->sigma_inner = sigma_inner;
    ((layer1_result_dog*)p)->sigma_outer = sigma_outer;
    ((layer1_result_dog*)p)->mask_size_factor = mask_size_factor;
}

void layer1_result_dog::read_from_stream(istreamer& is)
{
    layer1_result::read_from_stream(is);

    is.read(sigma_inner);
    is.read(sigma_outer);
    is.read(mask_size_factor);
}


void layer1_result_dog::write_to_stream(ostreamer& os)
{
    layer1_result::write_to_stream(os);
    
    os.write(sigma_inner);
    os.write(sigma_outer);
    os.write(mask_size_factor);
}



// layer1_result_loggabor
///////////////////////////////////////////////////////////////////////////////

void layer1_result_loggabor::get_masks(vector<img*>& masks) 
{ 
    layer1_creator_loggabor::get_filter_kernels(masks); 
}


// global functions
///////////////////////////////////////////////////////////////////////////////


void read_layer1_result(layer1_result*& res, const string& fname)
{
    try {
        res = (layer1_result*)streamable::read(fname);
    } catch (...) {
        res = nullptr;
    }
}

layer1_result* from_visview(const string& fname)
{
    ifstream is(fname.c_str());


/*

    for (int i = 0; i < count; ++i) {
        node* n = nodes[i];
        layer1_data* d;

        do {
            layermap.insert(pair<node*,int>(n, ++linei));
            d = (layer1_data*)(n->data);
            os << d->x + 1 << ',' << d->y + 1 << ",1," << d->m + 1 << ',' 
                << d->r(R_RESPONSE) << ',' << d->r(G_RESPONSE) << ',' << d->r(RR_RESPONSE) << endl;
        } while ((n = d->next) != nullptr);
    }
    os.close();

    // save links
    if (z > 0) {
        string name2;

        if (save_mode == 1) {
           name = dir + fname + "_" + (z + 1) + "_links.txt";
           name2 = dir + fname + "_" + (z + 1) + "_subparts.txt";
        } else {
           name = dir + nameint + "_image_" + (z + 1) + "_1_links.txt";
           name2 = dir + nameint + "_image_" + (z + 1) + "_1_subparts.txt";
        }

        ofstream os2;

        os.open(name.c_str(), ios::trunc);
        if (os.fail()) return;
        os2.open(name2.c_str(), ios::trunc);
        if (os2.fail()) return;
    
        int edgename = atom("toPrevLayer").get_index();
        for (int i = 0; i < count; ++i) {
            node* n = nodes[i];
            
            do {
                bool first = true;

                foreach_neighbor(n, edgename, iter) {
                    node* m = neighbor_node(iter);
                    edge_data_name* ed = (edge_data_name*)neighbor_edge_data(iter);

                    map<node*, int>::const_iterator j = prevlayermap.find(m);
                    if (j != prevlayermap.end()) {
                        if (first) first = false; else { os << ','; os2 << ','; }
                        os << j->second;
                        if (ed == nullptr) os2 << -1; else os2 << 200*(ed->data.x + 100) + ed->data.y + 100;
                    }
                }
                os << endl;
                os2 << endl;
            } while ((n = ((layer1_data*)n->data)->next) != nullptr);
        }
*/
    
    is.close();
    return nullptr;
}

void read_layer1_result_visview(layer1_result*& res, const string& fname)
{
    try {
        
    } catch (...) {
        
    }
}

void save_layer1_result(layer1_result* res, const string& fname)
{
    res->save(fname, -1);
}


layer1_result* read_layer1_result(const string& fname)
{
    try {
        return (layer1_result*)streamable::read(fname);
    } catch (...) {
        return nullptr;
    }
}

void save_node_set_mma(const string& fname, const set<node*>& nodes)
{
    ofstream os(fname.c_str());

    os << '{';
    for (set<node*>::const_iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;

        if (iter != nodes.begin()) os << ',' << endl;
        os << '{' << nd->m << ',';
        os << '{' << nd->x << ',' << nd->y << ',' << nd->z << '}';
        os << '}';
    }
    os << '}' << endl;
    os.close();
}


bool compare_node_data(node* n1, node* n2)
{
	layer1_data* n1d = (layer1_data*) n1->data;
	layer1_data* n2d = (layer1_data*) n2->data;
	//location
	if(n1d->x < n2d->x) return true;
	else return false;
	if(n1d->y < n2d->y) return true;
	else return false;
	if(n1d->z < n1d->z) return true;
	else return false;
	//part
	if(n1d->m < n2d->m) return true;
	else return false;
	if(n1d->val() < n1d->val()) return true;
	return false;
}



bool compare_node_neighbours(node* n1, node* n2)
{
	//sort neighbours and compare them
	vector<node*> nn1;
	vector<node*> nn2;
	forall_neighbors(n1, ineigh)
	{
		nn1.push_back(neighbor_node(ineigh));
	}
	forall_neighbors(n2, ineigh)
	{
		nn2.push_back(neighbor_node(ineigh));
	}
	if(nn1.size() < nn2.size()) return true;
	sort(nn1.begin(), nn1.end(), compare_node_data);
	sort(nn2.begin(), nn2.end(), compare_node_data);
	vector<node*>::iterator i1 = nn1.begin();
	vector<node*>::iterator i2 = nn2.begin();
	for(; i1 != nn1.end();)
	{
		if(compare_node_data(*i1,*i2) == true) return true;
		else return false;
		++n1;
		++n2;
	}
	return false;
}

bool compare_node(node* n1, node* n2)
{
	if(compare_node_data(n1, n2) == true) return true;
	else return false;
	if(compare_node_neighbours(n1, n2) == true) return true;
	return false;
}

bool compare_node_issame(node* n1, node* n2)
{
	// if !n1<n2 && !n2<n1
	if(compare_node(n1,n2) == compare_node(n1,n2)) return true;
	return false;
}

void print_node_short(node* n)
{
	printf("\t");
	layer1_data* nd = (layer1_data*)n->data;
	printf("x:%d y:%d m:%d v:%f | ",nd->x, nd->y, nd->m, nd->val());
}

void print_node(node* n)
{
	print_node_short(n);
	layer1_data* nd = (layer1_data*)n->data;
	n = nd->next;
	while(n != nullptr)
	{
		layer1_data* nd = (layer1_data*)n->data;
		printf("m:%d v:%f | ",nd->m, nd->val());
		//also print next nodes
		n = nd->next;
	}
	printf("\n");
}

template <class InputIterator1, class InputIterator2, class OutputIterator, class Compare>
  void
    custom_set_symmetric_difference ( InputIterator1 first1, InputIterator1 last1,
                               InputIterator2 first2, InputIterator2 last2,
                               OutputIterator result1, OutputIterator result2, Compare comp )
{
  while (true)
  {
    if (comp(*first1,*first2)) *result1++ = *first1++;
    else if (comp(*first2,*first1)) *result2++ = *first2++;
    else { first1++; first2++; }
	if (first1==last1) 
	{
		copy(first2,last2,result2);
		return;
	}
	if (first2==last2)
	{
		copy(first1,last1,result1);
		return;
	}
  }
  
    
}


/* ************ REMOVED ********  
void layer1_result::printLyDiff(layer1_result* r2)
{
	bool diff = false;
	for(int i = 0; i <= max_layer_index() || i <= r2->max_layer_index(); ++i)
	{
		//each layer
		printf("layer %d ----",i);
		//check if layer exists
		if(i > max_layer_index())
		{
			printf("\nobj1: no layer");
			diff = true;
			continue;
		}
		if(i > r2->max_layer_index())
		{
			printf("\nobj2: no layer");
			diff = true;
			continue;
		}
		//compare dim.
		if(x_size(i) != r2->x_size(i) || y_size(i) != r2->y_size(i))
		{
			printf("Different grid dimensions\n");
			diff = true;
			continue;
		}
		else
		{
			printf("%d x %d\n",x_size(i),y_size(i));
		}

		//compare grid nodes and neighbours by sorting them by their attrib. and comparing them
		vector<node*> lynodes1;
		vector<node*> lynodes2;
		get_layer_nodes(lynodes1, i);
		r2->get_layer_nodes(lynodes2, i);
		if(lynodes1.size() != lynodes2.size())
		{
			printf("Different number of shape nodes. file1: %d  file2: %d.\n", lynodes1.size(), lynodes2.size());
			diff = true;
			//continue;
		}
		sort(lynodes1.begin(), lynodes1.end(), compare_node);
		sort(lynodes2.begin(), lynodes2.end(), compare_node);
		vector<node*> diffnodes1;
		vector<node*> diffnodes2;
		custom_set_symmetric_difference ( lynodes1.begin(), lynodes1.end(),
			lynodes2.begin(), lynodes2.end(), back_inserter(diffnodes1),  back_inserter(diffnodes2), compare_node);
		
		int diffsum = diffnodes1.size() + diffnodes2.size();
		if(diffsum) 
		{
			printf("%d nodes differ:\n",diffsum);
			printf("set1 %d nodes:\n",diffnodes1.size());
			for_each(diffnodes1.begin(), diffnodes1.end(), print_node);
			printf("set2 %d nodes:\n",diffnodes2.size());
			for_each(diffnodes2.begin(), diffnodes2.end(), print_node);
		}
	}
	if(!diff) printf("Files are indentical.");
}

*/

void layer1_result::get_layer_nodes(vector<node*>& result, int z)
{
	vector<node*>& s_nodes = shape_nodes[z];
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        do {
            layer1_data* nd = (layer1_data*)n->data;
			result.push_back(n);
            n = nd->next;
        } while (n != nullptr);
    }
}

void layer1_result::get_layer_nodes(vector<node*>& result, int z, double bestpercent)
{
	vector<node*>& s_nodes = shape_nodes[z];
    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        double thresh = 0.0;

        do {
            layer1_data* nd = (layer1_data*)n->data;
            double val = nd->vval();

            if (thresh == 0.0) thresh = bestpercent*val;
            else if (val < thresh) break;
			result.push_back(n);
            n = nd->next;
        } while (n != nullptr);
    }
}

void layer1_result::debug_print_node_value_less(const double val, int layer)
{
	vector<node*> nodes;
	get_layer_nodes(nodes, layer);
	for(vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter)
	{
		layer1_data* nd = (layer1_data*)(*iter);
		if(nd->val() < val)
		{
			print_node_short(*iter);
		}
	}
}

void check_link(node* n, int link_name)
{
    if (!n->has_neighbor(link_name)) {
        cout << "Neighbor '" << link_name << "' does not exist." << endl;
        throw new_libhop_exception("Neighbor does not exist. Probably toLayer0?");
    }
}

// Get path map 'pmap' from node 'n'.
// Return value: end-node for path (0, 0, ..., 0).
ipoint2 get_path_map(path_map_t& pmap, layer1_result* res, const scmap_t& scmap, node* n, bool link)
{
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");
    ipoint2 result;

    pmap.clear();
    if (node_layer(n) == 0) {
        layer1_data* nd = (layer1_data*)n->data;
        hpoint_t& hp = pmap[vector<int>()];
        ipoint2 p(nd->x, nd->y);

        hp.p = p;
        hp.h = scmap.find(p)->second;
        result = p;
    } else { 
        //clock_t start, end;
        //start = clock();
        if (link) link_path(res, n, toprev, to0);
        else check_link(n, to0);
        //
        //end = clock();
        //cout << (double)(end - start)/CLK_TCK << "s ";
        foreach_neighbor (n, to0, niter) {
            layer1_data* nd = (layer1_data*)neighbor_node_data(niter);
            edge_path_data_t* ed = (edge_path_data_t*)neighbor_edge_data(niter);
            vector<int> v = ed->data.edges();
            pair<path_map_t::iterator, bool> ibpair = pmap.insert(path_map_t::value_type(ed->data.edges(), hpoint_t()));
            ipoint2 p(nd->x, nd->y);

            if (v.empty() || *max_element(v.begin(), v.end()) == 0)
                result = p;
            if (ibpair.second) {
                ibpair.first->second.h = scmap.find(p)->second;
                ibpair.first->second.p = p;
            } // else average ???!!!
        }
    }
    return result;
}

node* follow_center_link(layer1_result* res, node* n)
{
    if (node_layer(n) == 0) return n;

    int name = atom("toPrevLayer");

    foreach_neighbor (n, name, niter) {
        node* n = neighbor_node(niter);
        edge_data_name* nned = (edge_data_name*)neighbor_edge_data(niter);

        if (nned == nullptr) 
            throw new_libhop_exception("Edge data not found; use add_edge_names = true");
        if (nned->index == 0)
            return follow_center_link(res, n);
    }
    return n;
}

// Returns geometry-, shape context matching values:
// 'benergy' is TPS bending energy.
// 'scdistance' is a sum (not average!) of sc histogram distances (chi^2)
// 'maxpts' is maximal number of points used when (if) calculating bending energy
// bending energy (TPS) is only calculated if 'calcbenergy' is true
// Returns the number of unmatched points.
int part_geometry_matching(double& benergy, vector<dpoint2>& dvector, double& scdistance, 
    const path_map_t& im, const path_map_t& jm, bool calcbenergy)
{
    int result = 0;

    benergy = 0.0; // bending energy
    scdistance = 0.0; // sum of shape context histogram distances
    vector<dpoint2> pts1, pts2;

    for (path_map_t::const_iterator iter = im.begin(); iter != im.end(); ++iter) {
        path_map_t::const_iterator jter = jm.find(iter->first);

        if (jter == jm.end()) {
            ++result;
            // throw new_libhop_exception("Unexpected error in ::part_geometry_matching.");
        } else { 
            pts1.push_back((dpoint2)iter->second.p);
            pts2.push_back((dpoint2)jter->second.p);
            
            scdistance += sc_histogram_distance(iter->second.h, jter->second.h);
            dvector.push_back((dpoint2)(iter->second.p - jter->second.p));
        } 
    }

    if (calcbenergy) {
        translate_and_scale(pts1);
        translate_and_scale(pts2);

        double a = optimal_rotation(pts1, pts2);

        rotate(pts1, a);
        benergy = TPS_energy(pts1, pts2);
    }
    return result;
}

// Returns geometry-, shape context matching values;
// 'p' and 'q' are library nodes with geometry to compare
// See overloaded function above for more info.
int part_geometry_matching(double& benergy, vector<dpoint2>& dvector, double& scdistance, 
    node* p, node* q, bool calcbenergy)
{
    path_map_t ppm, qpm;

    get_library_geo(ppm, p);
    get_library_geo(qpm, q);
    return part_geometry_matching(benergy, dvector, scdistance, ppm, qpm, calcbenergy);
}

// Matches 'im' to 'jm' and vica-verse; for each of the two distances returns 
// the largest of the two values
void part_geometry_matching_sym(double& benergy, vector<dpoint2>& dvector, double& scdistance, 
    const path_map_t& im, const path_map_t& jm, bool calcbenergy)
{
    double be1, be2, sc1, sc2;
    int unm;
    double f;

    unm = part_geometry_matching(be1, dvector, sc1, im, jm, calcbenergy);
    f = (double)im.size()/(im.size() - unm); // divide by zero?
    be1 *= f; sc1 *= f/im.size();
    //unm = part_geometry_matching(be2, sc2, jm, im);
    //f = (double)jm.size()/(jm.size() - unm);
    //be2 *= f; sc2 *= f/jm.size();
    //benergy = max(be1, be2);
    //scdistance = max(sc1, sc2);
    benergy = be1;
    scdistance = sc1;
}



// Checks geometry of node 'n' against learned geometry given by 'gmap'.
// 'center' is assumed center of matching; i.e. from the support of n, center should be subtracted before
// (distance) matching.
/*void layer1_result::check_geometry(double& benergy, vector<dpoint2>& dvector, double& scdistance, 
    const scmap_t& scmap, node* n, const path_map_t& gmap, const ipoint2& center) 
{
    path_map_t pm;
    int unmatched;
    double f;

    get_path_map(pm, this, scmap, n);
    for (path_map_t::iterator pmiter = pm.begin(); pmiter != pm.end(); ++pmiter)
        pmiter->second.p -= center;
    unmatched = part_geometry_matching(benergy, dvector, scdistance, gmap, pm);
    f = (double)gmap.size()/(gmap.size() - unmatched); // divide by zero?
    benergy *= f; 
    scdistance *= f/gmap.size();
}*/



// Public functions
///////////////////////////////////////////////////////////////////////////////

// "Keeps" (i.e. does not mark them with HIDDEN_NODE_ATTR) nodes on layer 'z';
// 'thresh' is the inhibition threshold which is performed greedy w/r to response 'response' 
// on best 'maxn' nodes.
void inhibit_layer(layer1_result* res, int z, int response, int maxn, double thresh)
{
    int to_prev = atom("toPrevLayer");
    int to_0 = atom("toLayer0");
    vector<node*> nvec;

    for (list<node*>::iterator iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
        if (node_layer(*iter) == z) {
            nvec.push_back(*iter);
            (*iter)->set_attr(HIDDEN_NODE_ATTR); // We hide all nodes by default and "unhide" them below
        }
    }
    sort(nvec.begin(), nvec.end(), response_sort_f(response));

    // 'Inhibit' nodes
    list<irectangle2> boxes;
    int maxi = min<int>(maxn, (int)nvec.size());

    for (int i = 0; i < maxi; ++i) {
        node* n = nvec[i];
        set<node*> nset;
        irectangle2 box;
		double maxt = 0.0;

        res->recurse_and_link(n, to_prev, to_0, nset);
        box = node_set_bounding_rectangle(nset.begin(), nset.end());

        for (list<irectangle2>::iterator riter = boxes.begin(); riter != boxes.end(); ++riter) {
			irectangle2& rbox = *riter;
            double r = (double)(rbox.intersection(box).area())/box.union_area(rbox);

			if (r > maxt) maxt = r;
		}
        if (maxt <= thresh) {
			boxes.push_back(box);
            n->clear_attr(HIDDEN_NODE_ATTR);
        }
    }

    int count = 0;
    for (vector<node*>::iterator viter = nvec.begin(); viter != nvec.end(); ++viter) {
        if ((*viter)->is_attr_set(HIDDEN_NODE_ATTR)) ++count;
    }
    cout << count << " hidden nodes at layer " << z << endl;
}




#ifdef OPENCL
// OpenCL functions
////////////////////////////////

static double cc_clock_timer_ocl_make_data = 0;
static int cc_clock_count_ocl_make_data = 0;

void layer1_result::ocl_make_data(int layer, bool overrride_existing_ocl_data) {

	CStopWatch clock_1;

	clock_1.startTimer();

	// if data for current layer already exists then do nothing unless overrride_existing_ocl_data is set to true (default is false)
	if (overrride_existing_ocl_data == false && ocl_shape_nodes.size() > layer && ocl_shape_nodes[layer].second > 0)  {
		return;
	}

	// verify that shape_nodes for this layer exist
	if (this->shape_nodes.size() < layer) {
		cout << "Unable to make OpenCL data from layer1_result if no shape nodes for this layer were found." << endl;	 
		throw std::exception();
	}
	// verify that we already have made opencl data for lower layers
	if (ocl_shape_nodes.size() < layer || ocl_edges.size() < layer ||
		ocl_shape_nodes_coord.size() < layer || ocl_shape_nodes_inhib_coord.size() < layer) {
		ocl_make_data(layer - 1, overrride_existing_ocl_data);
		//cout << "Unable to make OpenCL data from layer1_result if previous layer does not have OpenCL data." << endl;
		//throw std::exception();
	}

	// we need to copy all those connections too
	int edge_type_mapping[][2] = { 
		{ atom("lyrCenterBack").get_index(), OCL_LYR_CENTER_BACK_EDGE },
		{ atom("lyrSrc").get_index(), OCL_LYR_SRC_EDGE },
		{ atom("lyrForbidden").get_index(), OCL_LYR_FORBIDDEN_EDGE },
		{ atom("toPrevLayer").get_index(), OCL_TO_PREV_LAYER_EDGE },
		{ atom("toSimilar").get_index(), OCL_TO_SIMILAR_EDGE },
		{ atom("toLayer0").get_index(), OCL_TO_LAYER_0 }
	};

	// make room for new layer
	ocl_shape_nodes.resize(layer + 1);
	ocl_edges.resize(layer + 1);
	ocl_shape_nodes_coord.resize(layer + 1);
	ocl_shape_nodes_coord_non_zero_count.resize(layer + 1);
	ocl_shape_nodes_inhib_coord.resize(layer + 1);
	ocl_shape_nodes_inhib_coord_non_zero_count.resize(layer + 1);

	int width = x_size(layer);
	int height = y_size(layer);

	int connections_count_all = 0;
	int part_counter = 0;

	// get array of nodes saved based on their location (x,y) on specific layer
	img_graph_layer& pos_layer = this->layers[layer];

	if (grid(layer) == nullptr)
		init_grid(layer);

	// in first loop only count how many parts we have and how many edges we have
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {

			// get node at this position
			node* n = pos_layer.node_at(i,j);
			
			// iterate over all nodes in same position
			while (n != nullptr) {

				// get data in node
				layer1_data* n_data = dynamic_cast<layer1_data*>(n->data);

				// count how many types of edges (connections we have)
				for (int j = 0; j < OCL_EDGE_COUNT; j++) {
					// get index from atom
					int atom_index = edge_type_mapping[j][0];
					// count how many elements we have
					int count = n->count_neighbors(atom_index);
					// and save to variable
					connections_count_all += count;
				}

				// move to next part in same position (i,j)
				n = n_data->next;

				// increment part counter
				part_counter++;
			}
		}
	}

	// reserve enough space for all parts in this layer
	int parts_count = part_counter; 
	int size = width * height;

	//cout << "Making memory object for parts for layer " << layer << " with parts " << parts_count << " and size " << parts_count * sizeof(ocl_part_data) << endl;

	ocl_layer1_data* parts = new ocl_layer1_data[parts_count];
	memset(parts, 0, parts_count * sizeof(ocl_layer1_data));

	ocl_layer1_data_coordinates* parts_coord = new ocl_layer1_data_coordinates[size];
	ocl_layer1_data_coordinates* parts_inhib_coord = new ocl_layer1_data_coordinates[size];

	//cout << "	size of parts_coord and parts_inhib_coord is " << size << "or size " << size * sizeof(ocl_layer1_data_coordinates) << endl;
	memset(parts_coord, 0, size * sizeof(ocl_layer1_data_coordinates));
	memset(parts_inhib_coord, 0, size * sizeof(ocl_layer1_data_coordinates));

	//cout << "Making memory object for for connections layer " << layer << " with edges " << connections_count_all << " and size " << connections_count_all * sizeof(ocl_edge_data_ip2) << endl;
	// we need to create memory for all edges since we only now know correct size
	ocl_edge_data_ip2* edges = new ocl_edge_data_ip2[connections_count_all];
	memset(edges, 0, connections_count_all * sizeof(ocl_edge_data_ip2));

	int coord_part_counter = 0;
	int edge_counter = 0;
	part_counter = 0;
	
	// should copy all parts into this new memory - but make sure to copy them in correct order
	// and initialize all edges
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {

			// get node at this position
			node* n = pos_layer.node_at(i,j);
			
			// iterate over all nodes in same position
			while (n != nullptr) {
				layer1_data* n_data = dynamic_cast<layer1_data*>(n->data);
				
				/**
				 * Initializing parts (layer1_data).
				 */ 

				// set pointer from n_data to this opencl structure
				n_data->ocl_struct = &parts[part_counter];

				// now copy each individual value to new memory
				parts[part_counter].x = n_data->x;
				parts[part_counter].y = n_data->y;
				parts[part_counter].z = n_data->z;
				parts[part_counter].m = n_data->m;
				parts[part_counter].attr = n->attr;
				parts[part_counter].response.R = n_data->r(R_RESPONSE);
				parts[part_counter].response.G = n_data->r(G_RESPONSE);
				parts[part_counter].response.RR = n_data->r(RR_RESPONSE);

				// set coordinates in parts_coord
				int coord_index = width * j + i;

				// if parts_coord are not 0 then this one is not first one so do not set offset
				if (parts_coord[coord_index].size == 0) {
					parts_coord[coord_index].offset = part_counter;
					// also count how many parts non zero locations we have
					coord_part_counter++;
				}
				
				// increment number of parts for this position
				parts_coord[coord_index].size++;

				/**
				 * Initializing edges / connections (edge_data_name).
				 */ 

				// for each edge type go through all edges
				for (int j = 0; j < OCL_EDGE_COUNT; j++) {
					// get index from atom
					int atom_index = edge_type_mapping[j][0];
					int edge_index = edge_type_mapping[j][1];

					// then copy all neighbors to global memory (also count how many we will copy)
					int copied = 0;
					foreach_neighbor(n, atom_index, iter) {
						node* pc = neighbor_node(iter); 
						layer1_data* pc_data = dynamic_cast<layer1_data*>(pc->data); 
						edge_data_name* pc_edge = dynamic_cast<edge_data_name*>(neighbor_edge_data(iter)); 
		
						// pc_edge might be null for some atom_index
						if (pc_edge != nullptr) {

							// set pointer from pc_edge to this opencl structure
							pc_edge->ocl_struct = &edges[edge_counter + copied];

                            edges[edge_counter + copied].x = pc_edge->index;
							edges[edge_counter + copied].y = pc_edge->index;
							edges[edge_counter + copied].r = pc_edge->r;

							edges[edge_counter + copied].type = OCL_EDGE_DATA_IP2_TYPE;
						} else {
							edges[edge_counter + copied].type = OCL_EDGE_DATA_EMPTY_TYPE;
						}

						// we need to find position of part being pointed to by this edge
						ocl_layer1_data* layer_pointed_parts = this->ocl_shape_nodes[pc_data->z].first;
						
						// calculate offset using pointer arithmetic
						edges[edge_counter + copied].node.offset = (pc_data->ocl_struct - layer_pointed_parts);
						edges[edge_counter + copied].node.layer = pc_data->z;
						
						copied++;
					}
					
					// set offset i.e. position of first element in global memory (end of previous counter)
					// and its size				
					n_data->ocl_struct->edges_loc[edge_index].offset = copied > 0 ? edge_counter : 0;
					n_data->ocl_struct->edges_loc[edge_index].size = copied;

					// increment counter for number of copied edges
					edge_counter += copied;
				}				
				
				// move to next part in same position (i,j)
				n = n_data->next;
				
				// increment edge counter
				part_counter++;
			}			
		}
	}

	//cout << "counter = " << counter << endl;

	/**
	 * Initializing inhibited parts.
	 */ 

	// we must also set coordinates for shape_nodes_inhib
	vector<node*> parts_inhib = shape_nodes_inhib[layer];
	int parts_inhib_size = parts_inhib.size();

	for (int i = 0; i < parts_inhib_size; i++) {
		// get node
		node* n = parts_inhib[i];
		// and data 
		layer1_data* n_data = (layer1_data*)n->data;

		// set coordinates in parts_inhib_coord
		int coord_index = width * n_data->y + n_data->x;

		// set offset using pointer arithmetics
		parts_inhib_coord[coord_index].offset = (n_data->ocl_struct - parts);
		parts_inhib_coord[coord_index].size = 1;

	}


	// save opencl all pointers for OpenCL data
	ocl_shape_nodes[layer] = pair<ocl_layer1_data*, int>(parts, parts_count);
	ocl_edges[layer] = pair<ocl_edge_data_ip2*, int>(edges, connections_count_all);

	ocl_shape_nodes_coord[layer] = pair<ocl_layer1_data_coordinates*, int>(parts_coord, size);
	ocl_shape_nodes_coord_non_zero_count[layer] = coord_part_counter;
	ocl_shape_nodes_inhib_coord[layer] = pair<ocl_layer1_data_coordinates*, int>(parts_inhib_coord, size);
	ocl_shape_nodes_inhib_coord_non_zero_count[layer] = parts_inhib_size;

	clock_1.stopTimer();

	double time = clock_1.getElapsedTime();
	cc_clock_count_ocl_make_data++;
	cc_clock_timer_ocl_make_data = time;
	//cout << "time to copy one layer " << time << endl;
	//cout << "avgrage time per layer with " << cc_clock_count_ocl_make_data << " iterations is " << cc_clock_timer_ocl_make_data/cc_clock_count_ocl_make_data << endl;

}


void layer1_result::make_data_from_ocl(int layer, bool add_edge_names,  bool overrride_existing_data) {

	clock_make_data_from_ocl.startTimer();
	
	int edge_type_mapping[][2] = { 
		{ atom("lyrCenterBack").get_index(), OCL_LYR_CENTER_BACK_EDGE },
		{ atom("lyrSrc").get_index(), OCL_LYR_SRC_EDGE },
		{ atom("lyrForbidden").get_index(), OCL_LYR_FORBIDDEN_EDGE },
		{ atom("toPrevLayer").get_index(), OCL_TO_PREV_LAYER_EDGE },
		{ atom("toSimilar").get_index(), OCL_TO_SIMILAR_EDGE },
		{ atom("toLayer0").get_index(), OCL_TO_LAYER_0 }
	};

	// make enough room for this new layer
    while ((int)this->shape_nodes.size() <= layer) this->shape_nodes.push_back(vector<node*>());
    while ((int)this->shape_nodes_inhib.size() <= layer) this->shape_nodes_inhib.push_back(vector<node*>());
    while ((int)this->info.size() <= layer) this->info.push_back(layer_info());

	// find out size of this layer and set it (should be already set)

	// from each valid position create node and layer1_data (from ocl_layer1_data) 
		// also set neighboors at the same time
	ocl_layer1_data* s_nodes = ocl_shape_nodes[layer].first;
	ocl_edge_data_ip2* edges = ocl_edges[layer].first;

	ocl_layer1_data_coordinates* coords = ocl_shape_nodes_coord[layer].first;
	ocl_layer1_data_coordinates* coords_inhib = ocl_shape_nodes_inhib_coord[layer].first;

	int s_nodes_size = ocl_shape_nodes[layer].second;
	int edges_size = ocl_edges[layer].second;

	int coords_size = ocl_shape_nodes_coord[layer].second;

	vector<node*>& new_layer_shape_nodes = shape_nodes[layer];
	vector<node*>& new_layer_shape_nodes_inhib = shape_nodes_inhib[layer];

	int shape_nodes_save_offset = 0;
	new_layer_shape_nodes.resize(ocl_shape_nodes_coord_non_zero_count[layer]);

	int shape_nodes_inhib_save_offset = 0;
	new_layer_shape_nodes_inhib.resize(ocl_shape_nodes_inhib_coord_non_zero_count[layer]);

	//layer1_data* layer1_data_buffer = new layer1_data[s_nodes_size];

	vector<node_pair> edges_array;
	int j = 0;

	for (int data_offset = 0; data_offset < s_nodes_size; data_offset += j) {
		

		// now get position where ocl_layer1_data are saved and construct node and layer1_data for each one
		//int data_offset = coords[i].offset;

		int x = s_nodes[data_offset].x;
		int y = s_nodes[data_offset].y;
		int z = s_nodes[data_offset].z;

		// make next parts at same position and append them to linked list of first part
		j = 0;
		node* new_n = nullptr;
		while (data_offset + j < s_nodes_size && s_nodes[data_offset + j].x == x && s_nodes[data_offset + j].y == y) {
			ocl_layer1_data* ocl_part = &s_nodes[data_offset + j];

			layer1_data* org_part_data = new layer1_data(ocl_part);

			org_part_data->x = x;
			org_part_data->y = y;
			org_part_data->z = z;
			
			if (new_n == nullptr) {

				node** ptr_new_n = &r_node_at(x, y, z);

				//*ptr_new_n = add_node(org_part_data, ocl_part->attr | IMG_NODE_ATTR);

				*ptr_new_n = node::new_node(ocl_part->attr | IMG_NODE_ATTR);
				new_n = *ptr_new_n;				

				// make node with correct layer1_data and push it to the list
				new_layer_shape_nodes[shape_nodes_save_offset++] = new_n;

				// also push it to list of inhibited nodes if indicated so by coords_inhib[i] 
				if (coords_inhib[y * x_size(layer) + x].size > 0)  {
					new_layer_shape_nodes_inhib[shape_nodes_inhib_save_offset++] = new_n;
				}

			} else {
				//((layer1_data*)new_n->data)->next = add_node(org_part_data, ocl_part->attr);
				((layer1_data*)new_n->data)->next = node::new_node(ocl_part->attr);	
				new_n = ((layer1_data*)new_n->data)->next;
			}

			new_n->data = org_part_data;
			nodes.push_back(new_n);


			node::iterator last_inserted = new_n->neighbors.begin();

/*			int edges_count = 0;
			for (int k = 0; k < OCL_EDGE_COUNT; k++) {
				edges_count += ocl_part->edges_loc[k].size;
			}

			edges_array.resize(edges_count);
			
			edges_count = 0;
*/
			// connect neighbors				
			for (int k = 0; k < OCL_EDGE_COUNT; k++) {
				// skip if we do not have any connection of this type
				if (ocl_part->edges_loc[k].size > 0 && ocl_part->edges_loc[k].offset + ocl_part->edges_loc[k].size <= edges_size) {
					
					int atom_index = edge_type_mapping[k][0];

					for (int m = 0; m < ocl_part->edges_loc[k].size; m++) {
						if (ocl_part->edges_loc[k].offset + m >= edges_size || ocl_part->edges_loc[k].offset + m < 0)
							continue;

						// get edge data
						ocl_edge_data_ip2* edge_data = &edges[ocl_part->edges_loc[k].offset + m];

						if (edge_data == nullptr ||  edge_data->node.layer < 0 || edge_data->node.layer >= ocl_shape_nodes.size() || ocl_shape_nodes[edge_data->node.layer].first == nullptr)
							continue;
							
						ocl_layer1_data* ocl_target_par = &ocl_shape_nodes[edge_data->node.layer].first[edge_data->node.offset];

						if (ocl_target_par == nullptr ) 
							continue;
						
						// then find actual part
						node* target_node = node_at(ocl_target_par->x, ocl_target_par->y, ocl_target_par->z);

						// different parts can be at same location so find correct one within
						while (target_node != nullptr && ((layer1_data*)target_node->data)->ocl_struct != ocl_target_par)
							target_node = ((layer1_data*)target_node->data)->next;
						
						//add_edge_2(new_n, target_node, new edge_data_name(edge_data), atom_index);
						if (add_edge_names == false) {
							last_inserted = new_n->neighbors.insert(last_inserted, node_pair(atom_index, edge_pair(target_node, nullptr)));
						} else {
							last_inserted = new_n->neighbors.insert(last_inserted, node_pair(atom_index, edge_pair(target_node, new edge_data_name(edge_data))));
						}
//						edges_array[edges_count++] = node_pair(atom_index, edge_pair(target_node, new edge_data_name(edge_data)));
					}
				}
			}
			//new_n->neighbors.insert(edges_array.begin(),edges_array.begin() + edges_count);
			
			j++;
		}
	}	
	clock_make_data_from_ocl.stopTimer();
	//cout << " exported to org structure from ocl data in : " << clock_make_data_from_ocl.getElapsedTime() << endl;
}

#endif

void link_path(layer1_result* res, node* n, int edge_name, int link_name)
{
    if (n->has_neighbor(link_name)) 
        return;

    foreach_neighbor(n, edge_name, iter) {
        node* nn = neighbor_node(iter);
        edge_data_name* nned = (edge_data_name*)neighbor_edge_data(iter);
        layer1_data* nnd = (layer1_data*)nn->data;
        int ecount = 0;
        
        if (nned == nullptr)  {
            cout << "nn->z = " << node_layer(nn) << endl;
            throw new_libhop_exception("Edge data not found; use add_edge_names = true");
        }

        link_path(res, nn, edge_name, link_name);
        foreach_neighbor(nn, link_name, niter) {
            edge_path_data_t* ed = (edge_path_data_t*)neighbor_edge_data(niter);
            edge_path p = ed->data;
            
            p.push_back(nned->index, nnd->m);
            res->add_edge_2(n, neighbor_node(niter), new edge_path_data_t(p), link_name);
            ++ecount;
        }
        if (ecount == 0) { // nn is layer1 node !? :pray
            edge_path p;
            
            p.push_back(nned->index, nnd->m);
            res->add_edge_2(n, nn, new edge_path_data_t(p), link_name);
        }
    }
}

/*void link_path(layer1_result* res, node* n, int edge_name, int link_name)
{
    scope_lock(res->lock);

    if (!n->has_neighbor(link_name)) {

        link_path_rec(res, n, edge_name, link_name);
    }
}*/

void get_histogram(vector<int>& h, const K_bin& bin, img_graph* g, int layer, const ipoint2& center, int width, int height)
{
    if (!g->grid(layer)) 
        g->init_grid(layer);

    int mcx = width/2, mcy = height/2;
    int maxx = min<int>(g->x_size(layer), center.x + mcx + 1);
    int maxy = min<int>(g->y_size(layer), center.y + mcy + 1);
    int minx = max<int>(0, center.x - mcx);
    int miny = max<int>(0, center.y - mcy);
    ip2_vector pts;

    for (int x = minx; x < maxx; ++x) {
        for (int y = miny; y < maxy; ++y) {
            node* n = g->node_at(x, y, layer);

            if (n != nullptr) {
                img_node_data* nd = (img_node_data*)n->data;
                
                pts.push_back(ipoint2(nd->x, nd->y));
            }
        }
    }
    bin.get_histogram(h, pts, center);
}

void get_sc_map(scmap_t& scmap, layer1_result* res, const K_bin& bin, bool normalize)
{
    typedef map<ipoint2, vector<int> > map_t;

    vector<node*>& s_nodes = res->shape_nodes[0];

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;
        ipoint2 p(nd->x, nd->y);
        vector<int> h;

        get_histogram(h, bin, res, 0, p, bin.width(), bin.height());
        if (normalize) normalize_histogram(scmap[p], h);
        else cast_vector(scmap[p], h);
    }
}

// Keeps only 'k' best boxes of each type.
// Note: 'boxes' list is resurned sorted by box_data::val 
void keep_best_boxes(list<layer1_result::box_data_t>& boxes, int k)
{
    boxes.sort(greater<layer1_result::box_data_t>());

    map<int, int> mcount;
    list<layer1_result::box_data_t>::iterator iter = boxes.begin();

    while (iter != boxes.end()) {
        if (mcount[iter->m]++ < k) iter++; else iter = boxes.erase(iter);
    }
}

void get_node_geo_p(vector<pair<int, ipoint2> >& ptsm, layer1_result* res, node* n)
{
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");

    ptsm.clear();
    foreach_neighbor (n, toprev, iter) {
        node* m = neighbor_node(iter);
        edge_data_name* med = (edge_data_name*)neighbor_edge_data(iter);
        layer1_data* md = (layer1_data*)m->data;
        set<node*> nset;
        set<ipoint2> pset;

        res->recurse_and_link(m, toprev, to0, nset);
        node_set_to_point_set(pset, nset.begin(), nset.end());

        for (set<ipoint2>::iterator pseti = pset.begin(); pseti != pset.end(); ++pseti) {
            ptsm.push_back(pair<int, ipoint2>(med->index, *pseti));
        }
    }
}

//void get_boxes(list<layer1_result::box_data_t>& boxes, layer1_result* res, const map<node*, vector<node*> >& clustering, 
//    int cattype, int response, double factor, int border0)
//{
//	int edgename = atom("toPrevLayer");
//    int to0 = atom("toLayer0");
//
//    for (auto cliter = clustering.begin(); cliter != clustering.end(); ++cliter) {
//        layer1_data* nd = (layer1_data*)cliter->first->data;
//		irectangle2 box;
//		set<node*> rec;
//
//        res->recurse_and_link(cliter->second, edgename, to0, rec);
//			    
//		for (auto niter = rec.begin(); niter != rec.end(); ++niter) {
//			layer1_data* nnnd = (layer1_data*)(*niter)->data;
//
//			if (nnnd->z == 0) box.eat(nnnd->x, nnnd->y);
//		}
//        box -= ipoint2(res->border, res->border);
//        box.resize(factor);
//        box += ipoint2(border0, border0);
//		boxes.push_back(layer1_result::box_data_t(box, nd->r, nd->r(response), (int)rec.size(), nd->m, cattype));
//    }
//}

void filter_nodes(list<node*>& nodes, layer1_result* res, int layer, const response_filter& rsf)
{
    nodes.clear();
    for (auto niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;

        if (nd->z == layer && rsf.check(nd)) 
            nodes.push_back(n);
    }
}

vector<node*> filter_nodes(const vector<node*>& nv, const response_filter& rsf)
{
	vector<node*> result;

    for (auto niter = nv.begin(); niter != nv.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;

        if (rsf.check(nd)) 
            result.push_back(n);
    }
	return result;
}


void sample_tree(set<node*>& result, const set<node*>& ns)
{
    int name = atom("toPrevLayer");
    set<node*> neighbors;

    for (auto iter = ns.begin(); iter != ns.end(); ++iter) {
        node* n = *iter;
        int count = 0;
        map<int, pair<vector<node*>, vector<double> > > nnbs;

        foreach_neighbor (n, name, niter) {
            edge_data_name* ed = (edge_data_name*)neighbor_edge_data(niter);

            if (ed == nullptr) throw;
            nnbs[ed->index].first.push_back(neighbor_node(niter));
            nnbs[ed->index].second.push_back(ed->r);
        }
        if (nnbs.empty()) result.insert(n);
        else {
            for (auto nniter = nnbs.begin(); nniter != nnbs.end(); ++nniter) 
                neighbors.insert(nniter->second.first[random_discrete(nniter->second.second)]);
        }
    }

    if (!neighbors.empty()) 
        sample_tree(result, neighbors);
}
    
// Recursively samples tree from nodes 'ns' in 'res' following toPrevLayer edges. 
// Path is determined according to the values on edges (edge_data_name::r)
// If some edge does not have edge data, exception is thrown.
void sample_tree(set<node*>& result, node* n)
{
    set<node*> ns;

    ns.insert(n);
    sample_tree(result, ns);
}

void get_node_geo(vector<pair<int, ipoint2> >& ptsm, layer1_result* res, node* n)
{
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");

    res->recurse_and_link(n, toprev, to0);

    ptsm.clear();
    foreach_neighbor(n, toprev, niter) {
        node* n = neighbor_node(niter);
        edge_data_name* nned = (edge_data_name*)neighbor_edge_data(niter);
        set<node*> nset;
        set<ipoint2> pset;

        if (nned == nullptr) 
            throw new_libhop_exception("Edge data not found; use add_edge_names = true");
        
        res->recurse_and_link(n, toprev, to0, nset);
        node_set_to_point_set(pset, nset.begin(), nset.end());

        for (auto piter = pset.begin(); piter != pset.end(); ++piter) 
            ptsm.push_back(pair<int, ipoint2>(nned->index, *piter));
    }
}

// Fills feature vector in "hoc" manner.
void fill_hoc_feature_vector(svm2::vector_t& v, layer1_result* res, int layer, int layer_size, const K_bin& bin, const ipoint2& bcenter,
    const irectangle2& rect)
{
    int ename = atom("toPrevLayer");
    int linkname = atom("toLayer0");
    double factor = (double)res->x_size(0)/res->x_size(layer);

	if (!res->grid(layer))
		res->init_grid(layer);
    v.resize(bin.bin_count()*layer_size, 0.0f);

	int imin = std::max<int>(int_round(rect.ll.x/factor), 0);
	int imax = std::min<int>(int_round(rect.ur.x/factor), res->x_size(layer));
	int jmin = std::max<int>(int_round(rect.ll.y/factor), 0);
	int jmax = std::min<int>(int_round(rect.ur.y/factor), res->y_size(layer));

	for (int i = imin; i < imax; ++i) {
		for (int j = jmin; j < jmax; ++j) {
			node* n = res->node_at(i, j, layer);

			while (n != nullptr) {
				layer1_data* nd = (layer1_data*)n->data;
				ipoint2 p(int_round(nd->x*factor), int_round(nd->y*factor));
                int b = bin.get_bin(p, bcenter);

                if (b >= 0) 
                    v[layer_size*b + nd->m] += nd->vval();
				n = nd->next;
			}
		}
	}


    //for (auto iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
    //    node* n = *iter;
    //    layer1_data* nd = (layer1_data*)n->data;

    //    if (nd->z == layer) {
    //        /*set<node*> nset;
    //        vector<double> binh(bin.bin_count(), 0.0);

    //        res->recurse_and_link(n, ename, linkname, nset);
    //        for (auto niter = nset.begin(); niter != nset.end(); ++niter) {
    //            ipoint2 p = node_coordinates(*niter);

    //            if (rect.inside(p)) {
    //                int b = bin.get_bin(p, bcenter);   

    //                if (b >= 0) binh[b] += nd->vval();
    //            }
    //        }

    //        double binmax = *max_element(binh.begin(), binh.end());

    //        if (binmax > 0.0) {
    //            for (int i = 0; i < binh.size(); ++i)
    //                if (binh[i] != 0.0) 
    //                    v[layer_size*i + nd->m] += binh[i]/binmax;
    //        }
    //        */
    //        ipoint2 p(int_round(nd->x*factor), int_round(nd->y*factor));

    //        if (rect.inside(p)) {
    //            int b = bin.get_bin(p, bcenter);

    //            if (b >= 0) 
    //                v[layer_size*b + nd->m] += nd->vval();
    //        }

    //    }
    //}
}


// Return clusters of nodes  
void cluster_detections_ms(map<node*, node*>& clmap, layer1_result* res, 
    const list<node*>& nodes, double sigma, int steps)
{
    double eps = 1E-6;
    double c = -1.0/(sigma*sigma*2);
    double radius = 2.0*sigma;
    double radius2 = radius*radius;
    double eps2 = eps*eps;

	clmap.clear();
    res->set_attr(nodes.begin(), nodes.end(), ATTR_MARKED);

    // Make clustering
    for (auto niter = nodes.begin(); niter != nodes.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;
        //auto fiter = clmap.find(n);

        //if (fiter != clmap.end()) 
        //    continue;
		if (n->is_attr_set(ATTR_MARKED1))
			continue;

        vector<dpoint2> pv;

        for (auto nniter = nodes.begin(); nniter != nodes.end(); ++nniter) {
            layer1_data* nnd = (layer1_data*)(*nniter)->data;

            if (nnd->m == nd->m) 
                pv.push_back(dpoint2(nnd->x, nnd->y));
        }
        
        vector<dpoint2> path;

        mean_shift(path, pv, dpoint2(nd->x, nd->y), radius, c, eps, steps);

        node* cn = get_closest_node(res, nd->z, nd->m, path.back(), 0, radius, ATTR_MARKED);
        auto cfiter = clmap.find(cn);

        if (cfiter != clmap.end()) 
            cn = cfiter->second;

        if (cn != nullptr) {
            for (auto piter = path.begin(); piter != path.end(); ++piter) {
                node* m = get_closest_node(res, nd->z, nd->m, *piter, 0, radius, ATTR_MARKED);

				if (m != nullptr && !m->is_attr_set(ATTR_MARKED1)) {//clmap.find(m) == clmap.end()) 
                    clmap.insert(pair<node*, node*>(m, cn));
					m->set_attr(ATTR_MARKED1);
				}
            }
        }
    }

    res->clear_attr(nodes.begin(), nodes.end(), ATTR_MARKED);
	res->clear_attr(nodes.begin(), nodes.end(), ATTR_MARKED1);
}

void cluster_detections_ms(map<node*, vector<node*> >& clusters, layer1_result* res, 
    const list<node*>& nodes, double sigma, int steps)
{
	map<node*, node*> clmap;

	cluster_detections_ms(clmap, res, nodes, sigma, steps);
    clusters.clear();
    for (auto cliter = clmap.begin(); cliter != clmap.end(); ++cliter)
        clusters[cliter->second].push_back(cliter->first);
}


void cluster_detections_ms(map<node*, vector<node*> >& clusters, layer1_result* res, int layer, const response_filter& rsf,
    double sigma, int steps)
{
    list<node*> nodes;

    filter_nodes(nodes, res, layer, rsf);
    cluster_detections_ms(clusters, res, nodes, sigma, steps);
}


node* get_closest_node(layer1_result* res, int layer, int type, const dpoint2& p, int minr, int maxr,
    unsigned attr)
{
    if (layer < 0 || layer >= res->max_layer_index() || res->shape_nodes[layer].empty())
        return nullptr;

    if (!res->grid(layer))
        res->init_grid(layer);

    ipoint2 ip(int_round(p.x), int_round(p.y));
    int minr2 = minr*minr;
    int maxr2 = maxr*maxr;
    int minx = std::max<int>(0, ip.x - maxr);
    int maxx = std::min<int>(res->x_size(layer) + 1, ip.x + maxr + 1);
    int miny = std::max<int>(0, ip.y - maxr);
    int maxy = std::min<int>(res->y_size(layer) + 1, ip.y + maxr + 1);
    node* result = nullptr;
    int resultd2 = INT_MAX;

    for (int x = minx; x < maxx; ++x) {
        for (int y = miny; y < maxy; ++y) {
            ipoint2 q(x, y);
            int d2 = (ip - q).norm2();

            if (d2 >= minr2 && d2 <= maxr2 && d2 < resultd2) {
                node* n = res->node_at(x, y, layer);

                while (n != nullptr) {
                    layer1_data* nd = (layer1_data*)n->data;
                    
                    if (nd->m == type && n->is_attr_set(attr)) {
                        result = n;
                        resultd2 = d2;
                        break;
                    }
                    n = nd->next;
                }
            }
        }
    }
    //if (resultd2 < INT_MAX)
        return result;
    //else
    //    return get_closest_position(res, layer, p, maxr, 2*maxr);
}

irectangle2 bounding_rectangle_of_projection(layer1_result* res, const vector<node*>& nodes)
{
	int toprev = atom("toPrevLayer");
	int to0 = atom("toLayer0");
	irectangle2 result;

    for (auto niter = nodes.begin(); niter != nodes.end(); ++niter) {
        set<node*> nset;

        res->recurse_and_link(*niter, toprev, to0, nset);

        irectangle2 rect1 = bounding_rectangle_of_nodes(nset.begin(), nset.end());

		result.eat(rect1.ll);
		result.eat(rect1.ur);
	}
	return result;
}

double layer0_cover_ratio(layer1_result* res, int l)
{
	int toprev = atom("toPrevLayer");
	int to0 = atom("toLayer0");
    set<ipoint2> ly0cover;

    for (auto niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
        node* n = *niter;
        set<node*> rec;

        if (z_coo(n) == l) {
            res->recurse_and_link(n, toprev, to0, rec);
            node_set_to_point_set(ly0cover, rec.begin(), rec.end());
        }
    }
    
    return (double)ly0cover.size()/res->shape_nodes[0].size();
}

double check_svmt(const svm_data& svmd, int nchildren, node* n, bool dfvalue)
{
    int ename = atom("toPrevLayer");

    if (svmd.svm == nullptr) 
        return 1;

    layer1_data* nd = (layer1_data*)n->data;
    //int nfeatures = 3 + 3*nchildren;
    int nfeatures = 3 + nchildren;
    cv::Mat test(1, nfeatures, CV_32FC1, cv::Scalar_<float>(0.0));

    //test.at<float>(0, 0) = nd->r(R_RESPONSE);
    test.at<float>(0, 0) = nd->r(G_RESPONSE);
    test.at<float>(0, 1) = nd->r(RR_RESPONSE);
    test.at<float>(0, 2) = reverse_s_response(nd->r(S_RESPONSE));
    foreach_neighbor(n, ename, niter) {
        layer1_data* nnd = (layer1_data*)neighbor_node_data(niter);
        edge_data_name* nned = dynamic_cast<edge_data_name*>(neighbor_edge_data(niter));

        if (nned == nullptr) {
            cout << "Edge names are not present. Use 'add_edge_names = true' in inference configuration file." << endl;
            throw;
        }

        //if ((float)nnd->r(G_RESPONSE) > test.at<float>(0, 3 + 3*nned->index)) {
        //    test.at<float>(0, 3 + 3*nned->index) = (float)nnd->r(G_RESPONSE);
        //    test.at<float>(0, 3 + 3*nned->index + 1) = (float)(nnd->x - nd->x)/5.0;
        //    test.at<float>(0, 3 + 3*nned->index + 2) = (float)(nnd->y - nd->y)/5.0;
        //}
        if ((float)nnd->r(G_RESPONSE) > test.at<float>(0, 3 + nned->index)) {
            test.at<float>(0, 3 + nned->index) = (float)nnd->r(G_RESPONSE);
        }
    }

    return svmd.svm->predict(test, dfvalue);
}


