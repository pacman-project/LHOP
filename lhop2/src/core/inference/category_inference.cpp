// category_layer_inference
///////////////////////////////////////////////////////////////////////////////

#include "category_inference.h"

void register_category_inference_classes() {
	/************************************************************************* 
	 *           !!!! DO NOT MODIFY OR REMOVE THIS FUNCTION   !!!!!
	 ************************************************************************/
	ClassRegister::registerFactory<AbstractLayerInference::IFactory, CategoryLayerInference::Factory>();
	ClassRegister::registerFactory<AbstractLayerInference::IFactory, ObjectLayerInference::Factory>();
}

LayerOutputObject* CategoryLayerInference::performInference(const AbstractLayerInputObject* input_object) {
	/*
	// Use this function to add category layer only!
	//void layern_creator::add_layer3(layer1_result* res, int layer)

    int k = layer - 1;
    int k1 = k - 1;

    if (k <= 0) return;

    int x_size = res->x_size(k1);
    int y_size = res->y_size(k1);
    
    res->new_grid(x_size, y_size, k);
    while ((int)res->shape_nodes.size() <= k) res->shape_nodes.push_back(vector<node*>());
    while ((int)res->shape_nodes_inhib.size() <= k) res->shape_nodes_inhib.push_back(vector<node*>());
    while ((int)res->info.size() <= k) res->info.push_back(layer_info());

    vector<node*>& s_nodes1 = res->shape_nodes[k1];
    vector<node*>& s_nodes = res->shape_nodes[k];
    vector<node*>& s_nodes_inhib = res->shape_nodes_inhib[k];
    vector<node*>& parts1 = library->parts[k1];

    int to_prev_layer = EdgeConnection::TO_PREV_LAYER;
    int to_center = EdgeConnection::TO_LYR_CENTER_BACK;
        
    for (vector<node*>::iterator niter = s_nodes1.begin(); niter != s_nodes1.end(); ++niter) {
        node* n = *niter;

        do {
            layer1_data* nd = (layer1_data*)n->data;
            node* cpart = parts1[nd->m];
            part_data* cpartd = (part_data*)cpart->data;

            foreach_neighbor(cpart, to_center, iter) {
			    node* pc = neighbor_node(iter);        // neighbor of cpart
                part_data* pcd = (part_data*)pc->data; // data of pc

                if (!pc->is_attr_set(CATEGORY_NODE_ATTR)) 
                    continue;

                layer1_data* ly1d = new layer1_data(nd->r, pcd->type);
                pair<node*, node**> fpair = res->find_at(ly1d, nd->x, nd->y, k);
                node* newn = fpair.first;

				newn = res->add_grid_node(ly1d, nd->x, nd->y, k);
				if (!n->is_attr_set(HAS_NEXT_LAYER)) {
                    res->inc_covered(k1);
                    n->set_attr(HAS_NEXT_LAYER);
                }
                
				if (newn != nullptr) {
                    res->add_edge(newn, n, to_prev_layer);
                }
            }
            n = nd->next;
        } while (n != nullptr);
    }
    inhibit_result(res, k);
	*/

	return nullptr;
}


LayerOutputObject* ObjectLayerInference::performInference(const AbstractLayerInputObject* input_object) {
/*
	// bool layern_creator::add_layer7(layer1_result* res, const scmap_t& scmap, int layer, const irectangle2& region)

	typedef vector<pair<int, proj_data_7*> > parallel_candidates_t;

    int to_center_back = EdgeConnection::TO_LYR_CENTER_BACK;
    int to_part = EdgeConnection::TO_LYR_SOURCE;
    int to_prev_layer = EdgeConnection::TO_PREV_LAYER;
    int to_layer_0 = EdgeConnection::TO_LAYER0;
    int to_prev_layerI = EdgeConnection::TO_PREV_LAYER_I;

    int k = layer - 1;
    int k1 = k - 1;

    if (k <= 0) return false;

    vector<node*>& parts1 = library->parts[k1];
    vector<node*>& parts = library->parts[k];

    if (res->max_layer_index() < k1 || res->shape_nodes[k1].empty() || 
        parts1.empty() || parts.empty()) return false;

    while ((int)res->shape_nodes.size() <= k) res->shape_nodes.push_back(vector<node*>());
    while ((int)res->shape_nodes_inhib.size() <= k) res->shape_nodes_inhib.push_back(vector<node*>());
    while ((int)res->info.size() <= k) res->info.push_back(layer_info());

    set_candidate_thresholds(res, k1);

    // Set the contraction; either use the value from cfg file or, in case of
    // "rlayer" sets the predetermined correct values.
    double contraction = layer_contraction;
    
    int x_size = int_round(res->x_size(k1)/contraction); // !floor!
    int y_size = int_round(res->y_size(k1)/contraction); // !floor!

    vector<ipoint2> directions;
    int dirnum = variation_dimension/2;

    for (int i = -dirnum; i <= dirnum; ++i)
        for (int j = -dirnum; j <= dirnum; ++j) 
            directions.push_back(ipoint2(variation_factor*i, variation_factor*j));


    if (!res->grid(k)) res->new_grid(x_size, y_size, k);
    if (!res->grid(k1)) res->init_grid(k1);

    vector<node*> s_nodes1;

    res->get_layer_nodes(s_nodes1, k1);

    if (s_nodes1.empty()) return false;

    // Add reconstruction edges in case of shape checking (since this is not thread-safe, it must be
    // done in advance)
    if (shape_check) {
        for (auto niter = s_nodes1.begin(); niter != s_nodes1.end(); ++niter) {
            //link_path(res, *niter, to_prev_layer, to_layer_0);
            res->recurse_and_link(*niter, to_prev_layer, to_layer_0);
        }
    }
    library->init_basic_part_numbers(k);

    bool result = false;
    vector<list<proj_data_7> > candidates(s_nodes1.size());

    cout << "|";

    #pragma omp parallel for
    for (int ni = 0; ni < (int)s_nodes1.size(); ++ni) {
        node* n = s_nodes1[ni];
        layer1_data* nd = (layer1_data*)n->data;
        int ntype = nd->m;

        if (nd->r(R_RESPONSE) < candidate_r_threshold || nd->r(G_RESPONSE) < candidate_g_threshold ||
                (!region.invalid() && !region.inside(nd->x, nd->y))) 
            continue;

        node* cpart = parts1[ntype];
        part_data* cpartd = (part_data*)cpart->data;
        node* nc = follow_center_link(res, n);
        layer1_data* ncd = (layer1_data*)nc->data;

        foreach_neighbor (cpart, to_center_back, cpn_iter) {
            node* pc = neighbor_node(cpn_iter);

            if (!allowed_parts.empty() && allowed_parts.find(((lib_data*)pc->data)->type) == allowed_parts.end())
                continue;

            part_data_2c* pced = (part_data_2c*)neighbor_edge_data(cpn_iter);
            part_data* pcd = (part_data*)pc->data;
            int bpn = library->get_basic_part_number(pc); // basic part number

            if (min_part_distance >= 0 && res->find_node_in_region(ipoint2(nd->x + pcd->cmx - pced->x, nd->y + pcd->cmy - pced->y), k, 
                min_part_distance, min_part_distance, pcd->type, 0) != nullptr) continue;

            int part_count = pc->count_neighbors(to_part);
            int part_size = part_count;
            int failnum = (int)((1.0 - get_thresh(pcd, RR_THRESH, realization_ratio_threshold))*part_count);
                    
            for (double factor = min_factor; factor < max_factor + 10E-6; factor += 0.2) { // make 0.2 parameter...
                for (int d = 0; d < (int)directions.size(); ++d) {
                    double bpr = 0; // count how many basic parts are realized
                    double r_sum = 0.0;
                    double g_val = (g_response_operation == 0) ? 0.0 : 1.0;
                    int part_count = 1.0;
                    double realization_ratio = 0.0;

                    layer1_result::sp_result_t* conn = new layer1_result::sp_result_t();
                    const ipoint2& dir = directions[d];
                    response_map r;
                    vector<inference_data> idata;
                    
                    foreach_neighbor(pc, to_part, iter2) {
                        node* p = neighbor_node(iter2);
                        lib_data* pd = (lib_data*)p->data;
                        part_data_2* ed = (part_data_2*)neighbor_edge_data(iter2);
                        layer1_result::sp_result_t newconn;
                        double c = library->contraction(pd->layer + 2, k1 + 1);
                        dnpair max;

                        max = res->schur_product_max_all_ol(
                            (int)(c*(nd->x + factor*(ed->x - pced->x + dir.x))), 
                            (int)(c*(nd->y + factor*(ed->y - pced->y + dir.y))),
                            ed->x, ed->y, ed->index, 
                            ed->distr, ed->app, pd->layer, &layer1_result::sgspf, convolution_threshold, newconn);
                        if (max.first >= convolution_threshold) {
                            layer1_data* nnd = (layer1_data*)max.second->data;
                            if (convolution_link_threshold >= 1.0) 
                                conn->push_back(layer1_result::sp_result_data_t(max.second, ed->index, max.first));
                            else
                                filter_links(*conn, newconn, convolution_link_threshold * max.first);
                            r_sum += nnd->r(R_RESPONSE);
                            bpr += library->get_basic_part_number(p)*nnd->r(RR_RESPONSE);
                            if (g_response_operation == 0) g_val += max.first; 
                            else g_val *= max.first;

                            if (shape_check && !scmap.empty()) {
                                idata.push_back(inference_data(ed->index, max.second));
                            }
                        }
                        ++part_count;
                    }
                    // Calculate RR_RESPONSE
					realization_ratio = bpr/bpn;

                    r_sum += nd->r(R_RESPONSE);
                    r_sum /= part_count;
                    r_sum = ::pow(r_sum, r_response_pow);
                    g_val = ::pow(g_val * realization_ratio, g_response_pow);

                    if (r_sum > get_thresh(pcd, R_THRESH, r_response_threshold) && 
                            realization_ratio >= get_thresh(pcd, RR_THRESH, realization_ratio_threshold) &&
                            g_val >= get_thresh(pcd, G_THRESH, g_response_threshold)) {
                        int type = pcd->type;
                        vs_part_data* vspcd = dynamic_cast<vs_part_data*>(pc->data);
                        bool ok = true;
                        double benergy = 0.0;
                        double scdistance = 0.0;
                        double maxdistance = 0.0;
                        vector<node*> missnodes;

                        if (shape_check && !scmap.empty()) {
                            if (vspcd != nullptr && !vspcd->pcad.mean.empty()) {
                                vector<pair<int, ipoint2> > ipts;
                                pair<double, dpoint2> tsdata = get_node_geo(ipts, res, idata);
                                vector<dpoint2> misspos;
                                        
                                ddpair result = check_pca_geometry_p(misspos, ipts, vspcd->pcad, pc);
                                scdistance = result.first;
                                ok = scdistance <= get_thresh(pcd, S_THRESH, s_response_threshold);
                                if (ok && link_missing_support && !misspos.empty()) {
                                    vector<ipoint2> missposi;

                                    for (auto mpiter = misspos.begin(); mpiter != misspos.end(); ++mpiter) {
                                        missposi.push_back(ipoint2(mpiter->x*tsdata.first/100.0 + tsdata.second.x,
                                            mpiter->y*tsdata.first/100.0 + tsdata.second.y));
                                    }
                                    get_support_nodes(missnodes, res, missposi, max(int_round(tsdata.first * 0.05), 3));
                                }
                            } else {
                                vector<pair<int, ipoint2> > ipts;
                                pair<double, dpoint2> tsdata = get_node_geo(ipts, res, idata);

                                scdistance = check_geometry2(ipts, pc);
                                ok = scdistance <= get_thresh(pcd, S_THRESH, s_response_threshold);
                            }
                        }

                        if (ok) {
                            r.set_response(R_RESPONSE, r_sum);
                            r.set_response(RR_RESPONSE, realization_ratio);
                            r.set_response(G_RESPONSE, g_val); 
                            r.set_response(S_RESPONSE, scdistance); 
                            r.set_response(X_RESPONSE, benergy);
                            candidates[ni].push_back(proj_data_7(r, type, pcd->cmx - pced->x, pcd->cmy - pced->y, conn, missnodes));
                        } else 
                            delete conn;
                    } else
                        delete conn;
                }
                    
            }
            
        } // end of "foreach 'to_center_back' neighbor"
    }

    cout << "|";

    // Prepare for "fast insertion"
    map<ipoint2, parallel_candidates_t> cooccurrent;   // new position -> (candidate index, proj_data*)
    vector<double> values;
    double value_threshold;

    if (new_positions_threshold < INT_MAX)
        values.reserve(candidates.size()*2);
    for (int ci = 0; ci < (int)candidates.size(); ++ci) {
        layer1_data* nd = (layer1_data*)(s_nodes1[ci]->data);

        for (auto citer = candidates[ci].begin(); citer != candidates[ci].end(); ++citer) {
            int new_x = int_round((nd->x + citer->x)/contraction); // !floor!
            int new_y = int_round((nd->y + citer->y)/contraction); // !floor!

            cooccurrent[ipoint2(new_x, new_y)].push_back(pair<int, proj_data_7*>(ci, &(*citer)));
            if (new_positions_threshold < INT_MAX) 
                values.push_back(citer->val());
        }
    }

    if (values.size() < new_positions_threshold) value_threshold = 0.0;
    else {
        sort(values.begin(), values.end(), greater<double>());
        value_threshold = values[new_positions_threshold - 1] - 1E-6;
    }

    vector<parallel_candidates_t> par_candidates;

    for (auto coociter = cooccurrent.begin(); coociter != cooccurrent.end(); ++coociter) 
        par_candidates.push_back(coociter->second);

    // Result of the following block (shape checking): 
    //  - vector of new 'layer1_result's and 
    //  - vector of new edges

    vector<list<new_edges_data_t> > new_edges(par_candidates.size());
    vector<map<int, pair<layer1_data*, new_edges_data_t*> > > new_data(par_candidates.size());
    
    for (int pci = 0; pci < (int)par_candidates.size(); ++pci) {

        for (auto citer = par_candidates[pci].begin(); citer != par_candidates[pci].end(); ++citer) {
            int ci = citer->first;
            node* n = s_nodes1[ci];
            layer1_data* nd = (layer1_data*)n->data;
            const proj_data_7& best = *citer->second;

            new_edges[pci].push_back(vector<new_edge_data_t>());

            itriple new_edges_info = get_new_edges(new_edges[pci].back(), k1, best.d->begin(), best.d->end());  

            get_new_support_edges(new_edges[pci].back(), best.msupp);
            if (citer->second->val() > value_threshold && (!add_reconstruction_edges || new_edges_info.second == 0 || 
                    (double)new_edges_info.third/new_edges_info.second <= hypo_ratio_threshold)) {
                int new_x = int_round((nd->x + best.x)/contraction); // !floor!
                int new_y = int_round((nd->y + best.y)/contraction); // !floor!

                if (new_x >= 0 && new_x < x_size && new_y >= 0 && new_y < y_size) {
                    layer1_data* ly1d = new layer1_data(best.r, best.type);

                    ly1d->x = new_x;
                    ly1d->y = new_y;
                    ly1d->z = k;

                    bool exists = new_data[pci].find(ly1d->m) != new_data[pci].end();
                    pair<layer1_data*, new_edges_data_t*>& depair = new_data[pci][ly1d->m];

                    if (!exists) {
                        depair.first = ly1d;
                        depair.second = &new_edges[pci].back();
                    } else if (depair.first->vval() >= ly1d->vval()) {
                        delete ly1d;
                    } else {
                        delete depair.first;
                        depair.first = ly1d;
                        depair.second = &new_edges[pci].back();
                    }

                    if (!n->is_attr_set(HAS_NEXT_LAYER)) {  // ?
                        res->inc_covered(k1);
                        n->set_attr(HAS_NEXT_LAYER);
                    }

                }
            }
            delete best.d;
        }
    }

    // We insert the candidates into the graph; at each position we insert nodes 
    // sorted by best value (G_RESPONSE)

    for (int ndi = 0; ndi < (int)new_data.size(); ++ndi) {
        vector<pair<layer1_data*, new_edges_data_t*> > v;

        v.reserve(new_data[ndi].size());
        for (auto nditer = new_data[ndi].begin(); nditer != new_data[ndi].end(); ++nditer) 
            v.push_back(nditer->second);
        sort(v.begin(), v.end(), [](const pair<layer1_data*, new_edges_data_t*>& a, const pair<layer1_data*, new_edges_data_t*>& b) 
            { return a.first->vval() < b.first->vval(); });
        for (auto viter = v.begin(); viter != v.end(); ++viter) {
            node* newn = res->add_grid_node(viter->first, viter->first->x, viter->first->y, k);

            add_edges(res, newn, *viter->second);
        }
    }
        

    double quot = res->cover_quotient(k1);

    count++;
    quot_sum += quot;
    if (quot < min_cover_quot) min_cover_quot = quot;

    if (copy_prev_layer) copy_layer(res, k);

    if (inhibit_layer_response >= 0) {
        inhibit_layer(res, k, inhibit_layer_response, inhibit_layer_max, inhibit_layer_thresh);
        if (inhibit_layer_delete) 
            res->delete_nodes(HIDDEN_NODE_ATTR);
    }

    inhibit_result(res, k);
    
	res->delete_edges(to_layer_0);
    return result;
	*/

	return nullptr;
}