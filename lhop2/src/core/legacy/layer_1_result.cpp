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
#include <opencv2/highgui/highgui.hpp>

#include "utils/utils.h"
#include "utils/graphs/graph_utils.h"

#include "core/legacy/part_lib.h"
#include "core/legacy/inference/layer_1_creators.h"

#include "core/legacy/layers.h"

#include "layer_1_result.h"

using namespace std;

const double response_map::empty_val = numeric_limits<double>::quiet_NaN();

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

void layer1_result::add_reconstruction_edges(int z)
{
    if (z < 0 || z > max_layer_index()) return;

    vector<node*>& s_nodes = shape_nodes[z];
    int to_prev = EdgeConnection::TO_PREV_LAYER;
    int to_0 = EdgeConnection::TO_LAYER0;
    int to_next = EdgeConnection::TO_NEXT_LAYER0 + z;

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
    int to_prev = EdgeConnection::TO_PREV_LAYER;
    int to_0 = EdgeConnection::TO_LAYER0;
    int to_next = EdgeConnection::TO_NEXT_LAYER0 + z;

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
    int to_prev = EdgeConnection::TO_PREV_LAYER;
    int to_0 = EdgeConnection::TO_LAYER0;

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

int layer1_result::max_layer_index()
{
    return (int)shape_nodes.size() - 1;
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
    int to_prev = EdgeConnection::TO_PREV_LAYER;

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

					}					
				}
				n = nd->next;
            } 
        }
        update_and_inhibit(layer);
    }
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

img* layer1_result::get_image(int z, int zt, bool paintuncov, const color& uncovered, const color& defcol, 
    const vector<int>& parts, const vector<color>& colors, bool bigpoints /* = false*/)
{
    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);
    if (zt < 0) zt = 0; else zt = min(zt, z);

    set<node*> end_nodes;
    layer_predicate pred(zt);

    recurse(shape_nodes[z], EdgeConnection::TO_PREV_LAYER, pred, end_nodes);
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

    recurse(shape_nodes_inhib[z], EdgeConnection::TO_PREV_LAYER, pred, end_nodes);
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

    recurse(start_node, EdgeConnection::TO_PREV_LAYER, pred, end_nodes);
    return get_image(end_nodes, x_size(zt), y_size(zt), false, color(20, 20, 20), color(255, 255, 255), 
        vector<int>(), vector<color>());
}

// Get bounding rectangle of all "reconstruction" nodes of 'n'.
irectangle2 layer1_result::get_box_with_cached_link(node* n)
{
	set<node*> starting_nodes;
    set<node*> nodes;
    irectangle2 result;

    recurse_and_link(n, EdgeConnection::TO_PREV_LAYER, EdgeConnection::TO_LAYER0, nodes);
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

    recurse_from_node(n, EdgeConnection::TO_PREV_LAYER, nodes);
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
        } 
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

	int edgename = EdgeConnection::TO_PREV_LAYER;
    int to0 = EdgeConnection::TO_LAYER0;

	for (list<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
		node* n = *iter;
		layer1_data* nd = (layer1_data*)n->data;

		if (n->is_attr_set(NODE_DELETED_ATTR))
			continue;
        if (nd->z == z) {
            node* nn = n->get_neighbor(edgename);
            layer1_data* nnd = (layer1_data*)nn->data;

            if (rfilter.check(nnd, use_lib_responses ? library : nullptr)) {
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
                node* nn = n->get_neighbor(EdgeConnection::TO_PREV_LAYER);
                int nnm = (nn == nullptr) ? -1 : node_type(nn);

                recurse_from_node(n, EdgeConnection::TO_PREV_LAYER, end_nodes);
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

img* layer1_result::get_image_boxed(part_lib* library, int z, const vector<int>& parts, const vector<color>& c, 
    const color& defcol, bool inhibit)
{
    typedef pair<irectangle2, color> box_data_t;

    if (z < 0) z = (int)shape_nodes.size() - 1;
    else z = min(z, (int)shape_nodes.size() - 1);

    img* res = get_image_reconstructed(library, 0, 0, vector<int>(), true, false);
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
                recurse_from_node(*iter, EdgeConnection::TO_PREV_LAYER, end_nodes);
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

img* layer1_result::get_image_reconstructed(part_lib* library, int z, int zt, const vector<int>& parts,
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

    //get_real_masks(masks);
    sort(sparts.begin(), sparts.end());
    rectangle2<int> box;
    vector<node*> all_nodes;
    combine_t combinef = get_combine_function();

    if (all) all_grid_nodes(all_nodes, shape_nodes[z]); 
    else all_nodes.insert(all_nodes.begin(), shape_nodes[z].begin(), shape_nodes[z].end());
    if (sparts.empty()) {
        recurse(all_nodes, EdgeConnection::TO_PREV_LAYER, pred, end_nodes);
    } else {
        vector<node*> srcnodes;
        int type;

        for (vector<node*>::iterator iter = all_nodes.begin(); iter != all_nodes.end(); ++iter) {
            type = ((layer1_data*)(*iter)->data)->m;
            if (binary_search(sparts.begin(), sparts.end(), type)) srcnodes.push_back(*iter);
        }
        recurse(srcnodes, EdgeConnection::TO_PREV_LAYER, pred, end_nodes);
    }
    for (set<node*>::iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter) {
        layer1_data* d = (layer1_data*)(*iter)->data;

        if (d->z != 0) continue;
        
		matrix<double> mask;
        ipoint2 maskc;

		vector<node*>& lparts = library->parts[d->z];

		part_data* pd = (part_data*)lparts[d->m]->data;		
		maskc = pd->get_mask(mask, lparts[d->m], library);
        
		img m(mask);
        (res->*combinef)(m, d->x - (int)(m.width)/2, d->y - (int)(m.height)/2, pow(d->vval(), factorpow));
        box.eat(d->x, d->y);
    }

    for (unsigned i = 0; i < masks.size(); i++) 
        if (masks[i]) delete masks[i];
    return res;
}   
void layer1_result::get_layer_nodes(set<node*>& result, int z, const vector<int>& parts, 
    double thresh, double thresh2, double thresh3, double thresh4, dpoint2* within_bounds)
{
    
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
            } 
            n = nd->next;
        } while (n != nullptr);

    }   
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
    int to_prev = EdgeConnection::TO_PREV_LAYER;

    while (n != nullptr) {
        layer1_data* nd = (layer1_data*)n->data;

        if (sparts.empty() || binary_search(sparts.begin(), sparts.end(), nd->m)) {
            part_data* pd = (part_data*)lparts[nd->m]->data;
            set<node*> end_nodes;
            ipoint2 p(0, 0);
            int count = 0;
            
            end_nodes.clear();
            recurse_from_node(n, to_prev, pred, end_nodes);
            for (set<node*>::iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter) {
                layer1_data* nd = (layer1_data*)(*iter)->data;
                p.add(nd->x, nd->y);
                ++count;
            }
            int x = int_round(((double)p.x)/count);
            int y = int_round(((double)p.y)/count);

			matrix<double> mask;
            ipoint2 maskc;

            maskc = pd->get_mask(mask, lparts[nd->m], library);
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

    int edgename = EdgeConnection::TO_PREV_LAYER;	
    //int edge0name = EdgeConnection::TO_LAYER0;
	int minmaxedge0name = EdgeConnection::TO_MINMAX_LAYER0;
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
                    box = get_nodes_bounding_rectangle(nset.begin(), nset.end());
					//box = bounding_rectangle_of_nodes(nset.begin(), nset.end());
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

        } else {
           name = dir + nameint + "_image_" + (z + 1) + "_1_links.txt";
           name2 = dir + nameint + "_image_" + (z + 1) + "_1_subparts.txt";
        }

        ofstream os2;

        os.open(name.c_str(), ios::trunc);
        if (os.fail()) return;
        os2.open(name2.c_str(), ios::trunc);
        if (os2.fail()) return;
    
        int edgename = EdgeConnection::TO_PREV_LAYER;
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

void save_layer1_result(layer1_result* res, const string& fname)
{
    res->save(fname, -1);
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
    int toprev = EdgeConnection::TO_PREV_LAYER;
    int to0 = EdgeConnection::TO_LAYER0;
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
        if (link) link_path(res, n, toprev, to0);
        else check_link(n, to0);

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

    int name = EdgeConnection::TO_PREV_LAYER;

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
    
    benergy = be1;
    scdistance = sc1;
}



// Public functions
///////////////////////////////////////////////////////////////////////////////

// "Keeps" (i.e. does not mark them with HIDDEN_NODE_ATTR) nodes on layer 'z';
// 'thresh' is the inhibition threshold which is performed greedy w/r to response 'response' 
// on best 'maxn' nodes.
void inhibit_layer(layer1_result* res, int z, int response, int maxn, double thresh)
{
    int to_prev = EdgeConnection::TO_PREV_LAYER;
    int to_0 = EdgeConnection::TO_LAYER0;
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
        box = get_nodes_bounding_rectangle(nset.begin(), nset.end());
		//box = bounding_rectangle_of_nodes(nset.begin(), nset.end());
		

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
    int toprev = EdgeConnection::TO_PREV_LAYER;
    int to0 = EdgeConnection::TO_LAYER0;

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

void sample_tree(set<node*>& result, const set<node*>& ns)
{
    int name = EdgeConnection::TO_PREV_LAYER;
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
    int toprev = EdgeConnection::TO_PREV_LAYER;
    int to0 = EdgeConnection::TO_LAYER0;

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
	return result;
}

irectangle2 bounding_rectangle_of_projection(layer1_result* res, const vector<node*>& nodes)
{
	int toprev = EdgeConnection::TO_PREV_LAYER;
	int to0 = EdgeConnection::TO_LAYER0;
	irectangle2 result;

    for (auto niter = nodes.begin(); niter != nodes.end(); ++niter) {
        set<node*> nset;

        res->recurse_and_link(*niter, toprev, to0, nset);

        //irectangle2 rect1 = bounding_rectangle_of_nodes(nset.begin(), nset.end());
		irectangle2 rect1 = get_nodes_bounding_rectangle(nset.begin(), nset.end());
		

		result.eat(rect1.ll);
		result.eat(rect1.ur);
	}
	return result;
}

double layer0_cover_ratio(layer1_result* res, int l)
{
	int toprev = EdgeConnection::TO_PREV_LAYER;
	int to0 = EdgeConnection::TO_LAYER0;
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
    int ename = EdgeConnection::TO_PREV_LAYER;

    if (svmd.svm == nullptr) 
        return 1;

    layer1_data* nd = (layer1_data*)n->data;

	int nfeatures = 3 + nchildren;
    cv::Mat test(1, nfeatures, CV_32FC1, cv::Scalar_<float>(0.0));

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

        if ((float)nnd->r(G_RESPONSE) > test.at<float>(0, 3 + nned->index)) {
            test.at<float>(0, 3 + nned->index) = (float)nnd->r(G_RESPONSE);
        }
    }

    return svmd.svm->predict(test, dfvalue);
}


