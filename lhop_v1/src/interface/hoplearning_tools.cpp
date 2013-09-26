/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */

#ifdef WIN32
#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
#include <windows.h>
#endif

#include "hop.h"

#include <stdio.h>

#include <ctime>
#include <tuple>
#include "layers/layer_learning.h"
#include "layers/optimization.h"
#include "utils/convert.h"
#include "graphs/graph_utils.h"


void update_nbhoods(map_learning& slresult, const string& dir, const string& fpatt, int max)
{
    string fname = dir;
    layer1_result* res;
    list<string> files;
    int i = 0;

    list_directory(files, dir + fpatt);
    if (max < 0) max = 100000;
    for (list<string>::iterator iter = files.begin(); iter != files.end() && i < max; ++iter, ++i) {
        cout << "Processing: " << *iter << endl;
        read_layer1_result(res, fname + *iter);
        if (res != nullptr) { 
            slresult.prepare_for_update(res);
            slresult.update(res); 
            delete res; 
        }
    }
}

void display_EM_step_distributions(map<int, map<int, pair<distribution2, histogram> > >& result)
{
    typedef map<int, map<int, pair<distribution2, histogram> > > result_t; 

    for (result_t::iterator riter = result.begin(); riter != result.end(); ++riter) {
        string fname = string("c:\\temp\\part") + riter->first + "dist.png";
        result_t::mapped_type& m = riter->second;
        distribution2 d;
        histogram h;

        // Save distribution
        for (result_t::mapped_type::iterator miter = m.begin(); miter != m.end(); ++miter) {
            d.update(miter->second.first);
        }

        rmatrix mt = d.to_matrix();
        double r = mt.maximum();

        if (r > 0.0) mt /= r;

        img im(mt);
        im.save_jet_colormap(fname.c_str(), -1000);

        // Display histogram
        cout << "#" << riter->first << endl;
        for (result_t::mapped_type::iterator miter = m.begin(); miter != m.end(); ++miter) {
            cout << "  name: " << miter->first << " hist: " << miter->second.second << endl;
        }
        
    }

}

void EM_step(part_lib* library, int source_layer, const string& dir, const string& patt)
{
    typedef map<int, map<int, pair<distribution2, histogram> > > result_t;

    config_dictionary icfg; // inference configuration
    
    icfg.from_string(
        "realization_ratio_threshold = 0.9"
        ";add_edge_names = true"
        ";layer_contraction = 1.0"
        ";identity_g_response = true"
        ";candidate_g_threshold_percent = 0.9"
        ";candidate_r_threshold_percent = 0.7"
    );

    layern_creator creator(icfg);
    list<string> files;
    int layer = source_layer + 1;
    result_t result;

    // Gather statistics
    creator.set_library(library);
    list_directory(files, dir + patt);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;
    
        cout << "Processing: " << *file;
        read_layer1_result(res, dir + *file);
        if (res == nullptr) 
            cout << " error reading file" << endl;
        else {
            creator.add_layer(res, layer + 1, 0);

            ddpair d;
                
            d = res->top_response_distribution(layer, R_RESPONSE);
            double rthresh = std::max<double>(0.0, d.first/* - 2*d.second*/);
            d = res->response_distribution(layer, G_RESPONSE);
            double gthresh = std::max<double>(0.0, d.first/* - 2*d.second*/);

            update_mean_map(result, res, layer, library, rthresh, gthresh);

            delete res;
            cout << " done" << endl;
        }
    }
    // Display distribution (debug only)
    display_EM_step_distributions(result);

    // Set distributions
    update_means(result, layer, library);

    // Prevent library from being disposed by creator!
    creator.set_library(nullptr);
}

/***********************
       U N U S E D
void update_or_result(vector<int>& result, map<iipair, int>& result2, layer1_result* res, int layer,
    const nodes_at_sort_f& sortf)
{
    typedef vector<int> result_t;
    typedef map<iipair, int> result2_t;
    res->add_reconstruction_edges(layer);
    res->init_grid(layer);

    for (int i = 0; i < res->shape_nodes[layer].size(); ++i) {
        vector<node*> v;
        node* n = res->shape_nodes[layer][i];
        layer1_data* nd = (layer1_data*)n->data;

        //res->sorted_nodes_at(v, nd->x, nd->y, nd->z, R_RESPONSE);
        //cout << "Best r-response: " << ((layer1_data*)(v[0]->data))->r(R_RESPONSE) << endl;
        //res->sorted_nodes_at(v, nd->x, nd->y, nd->z, G_RESPONSE);
        //cout << "Best g-response: " << ((layer1_data*)(v[0]->data))->r(G_RESPONSE) << endl;

        vector<node*> nvec;

        res->nodes_at(nvec, nd->x, nd->y, nd->z);
        sort(nvec.begin(), nvec.end(), sortf);
        
        node* bestn = nvec.front();
        layer1_data* bestnd = (layer1_data*)bestn->data;
        set<node*> bestset;

        ++result[bestnd->m];

        bestn->get_neighbor_set(atom("toLayer0"), bestset);
        for (vector<node*>::iterator niter = ++nvec.begin(); niter != nvec.end(); ++niter) {
            node* n = *niter;
            layer1_data* nd = (layer1_data*)n->data;
            set<node*> nset;

            n->get_neighbor_set(atom("toLayer0"), nset);
            if (is_subset(nset.begin(), nset.end(), bestset.begin(), bestset.end())) {
                double rdif = bestnd->r(R_RESPONSE) - nd->r(R_RESPONSE);
                double gdif = bestnd->r(G_RESPONSE) - nd->r(G_RESPONSE);

                //cout << "m = " << nd->m << " issubset of bestm = " << bestnd->m << 
                //    "  r-dif = " << rdif << "  g-dif = " << gdif << endl;
                //if (rdif < 0.1 && gdif < 0.1) {
                //    result2_t::iterator r2iter = 
                //        result2.insert(result2_t::value_type(sort(iipair(nd->m, bestnd->m)), 0)).first;
                //    ++r2iter->second;
                //    ++result[nd->m];
                //}
            }
        }
                
        //for (vector<node*>::iterator iter = nvec.begin(); iter != nvec.end(); ++iter) {
        //    layer1_data* nnd = (layer1_data*)(*iter)->data;

        //    cout << "m = " << nnd->m << ", covering: " << (*iter)->count_neighbors(atom("toLayer0")) << 
        //        ", r-response: " << nnd->r(R_RESPONSE) << ", g-response: " <<
        //        nnd->r(G_RESPONSE) << endl;
        //}
        //cout << endl;
        //int dummy;
        //cin >> dummy;
    }
}
************************/

void update_or_library(vector<int>& result, map<iipair, int>& result2, part_lib* library, int layer)
{
    typedef vector<int> result_t;
    typedef map<iipair, int> result2_t;

    vector<int> final;
    result_t::iterator riter;

    //cout << "Result:" << "---------------" << endl;
    while (*(riter = max_element(result.begin(), result.end())) > 0) {
        int pos = riter - result.begin();

        //cout << pos << ": #" << *riter << endl;

        for (result2_t::iterator r2iter = result2.begin(); r2iter != result2.end(); ++r2iter) {
            if (r2iter->first.first == pos) { 
                result[r2iter->first.second] -= r2iter->second;
                r2iter->second = 0;
            } else if (r2iter->first.second == pos) {
                result[r2iter->first.first] -= r2iter->second;
                r2iter->second = 0;
            }
        }
        result[pos] = 0;
        final.push_back(pos);
    }

    //vector<iipair> result3;
    //for (result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
    //    result3.push_back(iipair(iter->second, iter->first));
    //}
    //sort(result3.begin(), result3.end(), greater<iipair>());

    //cout << endl;
    //cout << "Result:" << endl;
    //cout << "--------------------------------" << endl;
    //for (result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
    //    cout << iter->first << ": " << iter->second << endl;
    //}

    //vector<int> order;

    //for (vector<iipair>::iterator iter = result2.begin(); iter != result2.end(); ++iter) {
    //    cout << iter->second << ": " << iter->first << endl;
    //    order.push_back(iter->second);
    //}

    //library->save_all("outlib.png", layer + 1, final, false /*cfg.get_value_bool("show_labels", false)*/);
}

void apply_to_result(part_learning& slresult, const string& dir, const string& fpatt)
{
    string fname = dir;
    layer1_result* res;
    list<string> files;
    K_bin bin(12, 2, 4, 7);

    list_directory(files, fname + fpatt);
    for (list<string>::iterator iter = files.begin(); iter != files.end(); ++iter) {
        cout << "Processing: " << *iter << endl;

        read_layer1_result(res, fname + *iter);
        if (res != nullptr) {
            scmap_t scmap;

            get_sc_map(scmap, res, bin, true);
            slresult.update(res, scmap); 
            delete res; 
        }
    }
}



void optimize_layer(part_lib* library, int layer, int tokeep, const config_dictionary& srccfg, const string& nspace,
    const string& dir, const string& patt)
{
    config_dictionary icfg; // inference configuration + optimization-specific parameters
    
    icfg.from_namespace(srccfg, nspace); // fill from srccfg.namespace 
    icfg.from_string( // override some specific settings
        "realization_ratio_threshold = 0.9;"
        "layer_contraction = 1.0;" // !?
    );

    layern_creator creator(icfg);
    list<string> files;
    optimization_set_streamed workingset;

    // Load files and perform inference, store the results in workingset
    creator.set_library(library);
    list_directory(files, dir + patt);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;
    
        cout << "Processing: " << *file;
        read_layer1_result(res, dir + *file);
        if (res == nullptr) 
            cout << " error reading file" << endl;
        else {
            cout << " <I";

            creator.add_layer(res, layer + 1, 0);

            cout << "><R";

            res->add_reconstruction_edges(layer);

            cout << "><S";

            workingset.add(streamed_pointer(res));

            cout << '>';

            delete res;
            cout << " done" << endl;
        }
    }

    double cover_thresh = icfg.get_value_double("cover_threshold", 0.8);
    double int_thresh = icfg.get_value_double("intersection_threshold", 0.2);
    int bite_size = icfg.get_value_double("bite_size", 10);
    set<int> keepset = set_range(tokeep);

    // Perform optimization and keep selected parts.    
    set<int> sel_parts = optimize_layer(library, workingset, layer, keepset, 
        cover_thresh, int_thresh, bite_size /* angle_thresh, a_thresh, e_thresh */);
    vector<int> v(sel_parts.begin(), sel_parts.end());

    library->keep_parts(layer, v);

    // Prevent library from being disposed by creator!
    creator.set_library(nullptr);
}

void optimize_layer_locally(part_lib* library, int layer, int tokeep, const config_dictionary& srccfg, const string& nspace,
    const string& dir, const string& patt)
{
    config_dictionary icfg; // inference configuration

    icfg.from_namespace(srccfg, nspace); // fill from srccfg.namespace 
    icfg.from_string( // override some specific settings
        "realization_ratio_threshold = 0.9;"
        "layer_contraction = 1.0;" // !?
    );

    layern_creator creator(icfg);
    list<string> files;
    optimization_set_streamed workingset;

    // Load files and perform inference, store the results in workingset
    creator.set_library(library);
    list_directory(files, dir + patt);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;
    
        cout << "Processing: " << *file;
        read_layer1_result(res, dir + *file);
        if (res == nullptr) 
            cout << " error reading file" << endl;
        else {
            cout << " <I";

            creator.add_layer(res, layer + 1, 0);

            cout << "><R";

            //res->add_reconstruction_edges(layer);
            res->add_reconstruction_edges_leq_fwd(layer);

            cout << "><S";

            workingset.add(streamed_pointer(res));

            cout << '>';

            delete res;
            cout << " done" << endl;
        }
    }

    // Perform optimization and keep selected parts.    
    int nparts = icfg.get_value_int("sample_size", -1);
    double grthresh = icfg.get_value_double("update_threshold", 0.7);
    double sthresh = icfg.get_value_double("statistics_threshold", 0.01);
    int bitesize = icfg.get_value_int("bite_size", 10);
    set<int> keepset = set_range(tokeep);

    cout << "statistics threshold = " << sthresh << endl;

    set<int> sel_parts = optimize_layer_locally(library, workingset, layer, keepset, nparts, grthresh, 
        sthresh, bitesize);

    vector<int> v(sel_parts.begin(), sel_parts.end());

    library->keep_parts(layer, v);

    // Prevent library from being disposed by creator!
    creator.set_library(nullptr);
}

void do_find_parts(config_dictionary& cfg_learn, const string& patt)
{
    cout << "action = learn_parts" << endl;

	int source_layer_index = cfg_learn.get_value_int("source_layer_index", 0);

	// support legacy config name
	if (source_layer_index == 0) {
		source_layer_index = cfg_learn.get_value_int("layer", 0);
		if (source_layer_index > 0)
			cfg_learn.set_value("source_layer_index",source_layer_index);
	}


	string ly_namespace = string("ly") + (source_layer_index + 1);
	
	config_dictionary cfg;
	cfg.from_namespace_priority(cfg_learn, 1, ly_namespace.c_str());

    int part_max_number = cfg_learn.get_value_int("part_max_number", 20);
    string libname;

    part_learning slresult(cfg);
    part_lib* library;
    string fname;

    cfg.get_value(libname, "part_lib_name", true);
    read_library(libname, library);

    if (library == nullptr) 
        throw new_libhop_exception("Library can not be opened.");

    // Read neighborhood file (== maps == result of do_update)
    // and find maxima
    cout << "Reading maps and finding maxima";

    cfg.get_value(fname, "nbhoods_file", true);
    slresult.reset(fname);

    cout << endl;

    // Save maxima to files - optionally
    if (cfg.get_value_bool("save_maxima", false)) {
        cout << "Saving maxima images";

        slresult.display_maxima(cfg.get_value_string("image_name", "c:\\TEMP\\maxima_%d%d.bmp"));

        cout << endl;
    }

    // Find parts
    string src_dir = cfg.get_value_string("src_dir", "");

    end_dir(src_dir);

    cout << "Loading files" << endl;

    apply_to_result(slresult, src_dir, patt);

    cout << "Max # of updates: (unknown)" << endl;

    // Filter sequences (to satisfy thresholds on repeatability, ...?)

    // Printing sequences - optionally
    if (cfg.get_value_bool("print_sequences", false)) {
        cout << "Printing sequences (unsupported)" << endl;
    }

    // Add parts to library
    cout << "Adding parts to library";

    int tokeep = library->layer_size(slresult.get_source_layer() + 1);
    vector<int> indices;
     
    cfg.get_value(indices, "library_indices");
    slresult.add_to_library(library, part_max_number);
    
    cout << " (" << tokeep << " parts kept, ";
    cout << library->layer_size(slresult.get_source_layer() + 1) - tokeep << " parts added)" << endl;

    // Learn g distributions - optionally
    if (cfg.get_value_bool("learn_g_distributions", true)) {
        cout << endl << "Learning g-distribution" << endl;

        string gdnamespace = cfg.get_value_string("g_distribution_namespace", "g_distribution");

        learn_g_distributions(library, slresult.get_source_layer(), cfg, gdnamespace, src_dir, patt);

        cout << endl;
    }

    fname = cfg.get_value_string("library_image_file", "");
    library->drop_images();
    library->save_all(change_extension(fname, "-0.png").c_str(), slresult.get_source_layer() + 2, 0, -1, 
        cfg.get_value_bool("show_labels", false));

    // Perform local optimization
    if (cfg.get_value_bool("local_optimization_step", false)) {
        cout << endl << "Local optimization" << endl;

        string lonamespace = cfg.get_value_string("local_optimization_namespace", "l_optimization");

        optimize_layer_locally(library, slresult.get_source_layer() + 1, tokeep, cfg, lonamespace, src_dir, patt);

        cout << endl;
    }

    // Perform global optimization
    if (cfg.get_value_bool("global_optimization_step", false)) {
        cout << endl << "Global optimization" << endl;

        string gonamespace = cfg.get_value_string("global_optimization_namespace", "g_optimization");

        optimize_layer(library, slresult.get_source_layer() + 1, tokeep, cfg, gonamespace, src_dir, patt);

        cout << endl;
    }

    fname = cfg.get_value_string("library_image_file", "");
    library->drop_images();
    library->save_all(change_extension(fname, "-1.png").c_str(), slresult.get_source_layer() + 2, 0, -1, 
        cfg.get_value_bool("show_labels", false));

    // Perform one EM-step - optionally
    int emsteps = cfg.get_value_int("em_steps", 0);

    for (int ems = 0; ems < emsteps; ++ems) {
        cout << endl << "EM step #" << ems + 1 << endl;

        EM_step(library, slresult.get_source_layer(), src_dir, patt);
    
        library->drop_images();
        library->save_all(change_extension(fname, string("-em") + (ems + 1) + string(".png")).c_str(), 
            slresult.get_source_layer() + 2, 0, -1, cfg.get_value_bool("show_labels", false));
    }
    cout << endl;

    // Save part images - optionally
    if (fname != "") {
        library->drop_images();
        library->save_all(fname.c_str(), slresult.get_source_layer() + 2, 0, -1, 
            cfg.get_value_bool("show_labels", false));
    }

    fname = cfg.get_value_string("library_region_file", "");
    if (fname != "") 
        library->save_all_regions(fname.c_str(), slresult.get_source_layer() + 2);

    // Save library
    fname = cfg.get_value_string("save_library", "");
    if (fname != "") {
        cout << "Saving library";

        library->save(fname);

        cout << endl;
    }

    delete library;
}

void do_update(config_dictionary& cfg_learn, const string& patt)
{
	int source_layer_index = cfg_learn.get_value_int("source_layer_index", 0);
	// support legacy config name
	if (source_layer_index == 0) {
		source_layer_index = cfg_learn.get_value_int("layer", 0);
		if (source_layer_index > 0)
			cfg_learn.set_value("source_layer_index",source_layer_index);
	}
	string ly_namespace = string("ly") + (source_layer_index+1);

	
	config_dictionary cfg;
	cfg.from_namespace_priority(cfg_learn, 1, ly_namespace.c_str());

    map_learning slresult(cfg);
    part_lib* library;
    string libname;

    cfg.get_value(libname, "part_lib_name", true);
    read_library(libname, library);

    if (library == nullptr) {
        cout << "Error: Library can not be loaded." << endl;
        return;
    }

    vector<node*>& vv = library->parts[0]; 
    /*for (size_t i = 0; i < vv.size(); ++i){
        part_data* d = (part_data*)vv[i]->data;
        d->region.print();
    }*/
    
    string src_dir = cfg.get_value_string("src_dir", "");
    string src_patt = cfg.get_value_string("src_patt", "");
    string dest;

    end_dir(src_dir);
    cfg.get_value(dest, "nbhoods_file", true);

    cout << "Updating neighborhoods..." << endl;

    update_nbhoods(slresult, src_dir, (patt.size() == 0) ? src_patt : patt, 
        cfg.get_value_int("max_files", -1));

    cout << "Writing neighborhoods...";

    slresult.write_maxima_to_stream(dest);

    cout << endl << "Saving neighborhoods images...";

    if (cfg.get_value_bool("save_nbhoods", false))
        slresult.display_statistics(cfg.get_value_string("image_name", "c:\\TEMP\\nbhood_%d%d.bmp"));

    cout << endl;
}

void merge_library(config_dictionary& cfg)
{
    string src;
    string dest;
    int layer;
    double thresh;
    part_lib* srclib = nullptr;

    cfg.get_value(src, "part_lib_name", true);
    cfg.get_value(dest, "save_library", true);
    cfg.get_value(layer, "layer", true);
    cfg.get_value(thresh, "merge_threshold", true);

    int tol = cfg.get_value_int("merge_tolerance", 0);

    read_library(src, srclib);

    if (srclib == nullptr) return;

    int mergecount = srclib->merge_part_types(layer, thresh, tol);
    cout << "Number of merged parts: " << mergecount << endl;

    srclib->save(dest);

    delete srclib;
}

void keep_parts(config_dictionary& cfg)
{
    string src;
    string dest;
    int layer;
    part_lib* srclib = nullptr;
    vector<int> tokeep;

    cfg.get_value(src, "part_lib_name", true);
    cfg.get_value(dest, "save_library", true);
    cfg.get_value(layer, "layer", true);
    cfg.get_value(tokeep, "keep_parts", true);

    read_library(src, srclib);

    if (srclib == nullptr) return;

    srclib->keep_parts(layer, tokeep);
    cout << "Keeping parts #";
    for (int i = 0; i < (int)tokeep.size(); ++i) 
        cout << tokeep[i] << ' ';
    cout << endl;

    srclib->save(dest);

    delete srclib;
}

void delete_parts(config_dictionary& cfg)
{
    string src;
    string dest;
    string file;
    int layer;
    part_lib* srclib = nullptr;
    set<int> to_delete;
    vector<int> to_delete_v;
    int i;

    cfg.get_value(src, "library", true);
    cfg.get_value(dest, "out_library", true);
    cfg.get_value(layer, "layer", true);
    cfg.get_value(file, "file", true);

    read_library(src, srclib);

    if (srclib == nullptr) return;

    ifstream is(file.c_str());

    if (is.is_open()) {
        while (true) {
            is >> i;
            if (is.fail()) break;
            to_delete.insert(i);
        }
    }
    to_delete_v.assign(to_delete.begin(), to_delete.end());
    srclib->delete_parts(layer, to_delete_v);
    srclib->save(dest);

    delete srclib;
}

void optimization_initialization(optimization_data& odata, config_dictionary& cfg, const string& pattern)
{
    string dir = cfg.get_value_string("src_dir", "");
    list<string> files;
    layer1_result* res;
    int count = 0;
    int maxcount = cfg.get_value_int("file_limit", INT_MAX);
    int start_layer = cfg.get_value_int("start_layer", 1);
    int end_layer = cfg.get_value_int("end_layer", 2);

    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + pattern);
    cout << "from_file = " << (cfg.get_value_bool("from_file", false) ? "true" : "false") << endl;
    cout << "pattern = " << pattern << endl;
    cout << "src_dir = " << dir << endl;
    cout << "start_layer = " << start_layer << endl;
    cout << "end_layer = " << end_layer << endl;
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        if (++count > maxcount) break;

        cout << "Processing " << *file << endl;
        read_layer1_result(res, dir + *file);
        if (res) {
            odata.add_to_test_set(streamed_pointer(res));
            delete res;
        }
    }
}

// Optimizes layers from start_layer to min_layer
void default_optimization(config_dictionary& cfg, const string& srcfile)
{
    cfg.from_file(cfg.get_value_string("import_configuration", "").c_str());

    optimization_data odata(cfg);

    optimization_initialization(odata, cfg, srcfile);

    if (odata.test_set_empty()) {
        cout << "Can not perform optimization on an empty set" << endl;
        return;
    }

    int overall_steps = cfg.get_value_int("overall_steps", 2);
    vector<cr_optimization_base*> lopts;
    int final_layer;

	vector<int> optseq;
    if (cfg.is_defined("start_layer")) {
        int start_layer = cfg.get_value_int("start_layer", 1);
        int end_layer = cfg.get_value_int("end_layer", 2);

		for (int l = start_layer; l <= end_layer; ++l) {
			optseq.push_back(l);
		}
    } else {
        cfg.get_value(optseq, "optimization_sequence", true);		
    }

	for (int i = 0; i < (int)optseq.size(); ++i) {
		string ly_namespace = string("ly") + optseq[i];
		config_dictionary ly_cfg;
		ly_cfg.from_namespace_priority(cfg, 1, ly_namespace.c_str());
        lopts.push_back(new cr_layer_optimization(&odata, ly_cfg, optseq[i]));
	}
    final_layer = optseq.back();

    cr_set_optimization cropt(&odata, overall_steps, lopts);

    // Fill the optimizer 
    cropt.set_library(odata.get_library());
    cropt.add_to_init_test_set(&odata);

    // Execute 
    cropt.execute();   
    
    // Get the results
    part_lib* lresult = cropt.get_best_library();
    list<streamed_pointer> sresult;

    cropt.get_best_test_set(sresult);

    // Save results
    string lib_export_name = cfg.get_value_string("lib_export_name", "lib");
    string res_export_name = cfg.get_value_string("res_export_name", "res");
	string name = cfg.get_value_string("out_library", lib_export_name + (final_layer ) + ".plb");
    int i = 0;

    PRINT_INFO("");
    PRINT_INFO("SAVING");

	lresult->save(name);

    if (cfg.get_value_bool("save_test_set", false)) {
        for (list<streamed_pointer>::iterator iter = sresult.begin(); iter != sresult.end(); ++iter) {
            layer1_result* res = (layer1_result*)iter->get();

            name = res_export_name + (i++) + string(".ly") + (final_layer );
            PRINT_INFO(name);
            res->delete_edges_complement(atom("toPrevLayer").get_index());
            res->save(name);
            delete res;
        }
    }

}

void perform_optimizations(config_dictionary& cfg, const string& srcfile)
{
    string type = cfg.get_value_string("optimization_type", "default");
    clock_t start, end;

    start = clock();

    if (type == "default") default_optimization(cfg, srcfile);
    end = clock();

    cout << "Processing time: " << (double)(end - start)/CLOCKS_PER_SEC << endl;
}

void from_images(config_dictionary& cfg, const string& srcfile)
{
    cout << "From images" << endl;
}


void random_selection(list<string>& files, int n, const string& file)
{
    if (n < 0 || n >= (int)files.size()) return;

    vector<int> perm;

    random_permutation(perm, (int)files.size());
    
    vector<int> todrop(perm.begin() + n, perm.end());
    ofstream os(file.c_str());
    list<string>::iterator iter;
    int count = 0;

    sort(todrop.begin(), todrop.end());
    iter = files.begin();
    while (iter != files.end()) {
        if (!binary_search(todrop.begin(), todrop.end(), count++)) ++iter;
        else {
            if (!os.fail()) os << *iter << endl;
            iter = files.erase(iter);
        }
    }
}

void learn_g_distribution(const config_dictionary& cfg, const string& fname)
{
    typedef map<pair<int, int>, online_distribution> result_t;

    cout << "Doing learn_g_responses" << endl;

    int layer;
    string srcdir = cfg.get_value_string("src_dir", "");
    string libname, outlibname;
    list<string> files;

    cfg.get_value(layer, "layer", true);
    cfg.get_value(libname, "library", true);
    cfg.get_value(outlibname, "out_library", true);

    part_lib* library;

    read_library(libname, library);

    if (library == nullptr) {
        cout << "Library " << libname << " can not be opened." << endl;
        return;
    }

    // Update "result"
    result_t result;

    end_dir(srcdir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, fname, srcdir, cfg.get_value_string("pattern", "")))
        list_directory(files, srcdir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        cout << "Processing: " << *file;
        read_layer1_result(res, srcdir + *file);
        if (res == nullptr) 
            cout << " error reading file" << endl;
        else {
            update_g_distribution(result, res, layer, library);

            delete res;
            cout << " done" << endl;
        }
    }

    // Set distribution
    set_g_distribution(library, layer, result, cfg.get_value_double("min_variance", 0.1), false);

    // Save library
    library->save(outlibname);
    
    delete library;
}

void perform_EM(const config_dictionary& cfg, const string& fname)
{
    typedef map<int, map<int, pair<distribution2, histogram> > > result_t;

    cout << "Doing EM (to be merged into part learning)" << endl;
    cout << "==========================================" << endl;

    int layer;
    string srcdir = cfg.get_value_string("src_dir", "");
    //int steps = cfg.get_value_int("steps", 1);
    string libname, outlibname;
    list<string> files;
    double rthresh, gthresh;

    cfg.get_value(layer, "layer", true);
    cfg.get_value(libname, "library", true);
    cfg.get_value(outlibname, "out_library", true);
    cfg.get_value(rthresh, "r_response_threshold", true);
    cfg.get_value(gthresh, "g_response_threshold", true);
   

    part_lib* library;

    read_library(libname, library);

    if (library == nullptr) {
        cout << "Library " << libname << " can not be opened." << endl;
        return;
    }

    // Update "result"
    result_t result;

    end_dir(srcdir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, fname, srcdir, cfg.get_value_string("pattern", "")))
        list_directory(files, srcdir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        cout << "Processing: " << *file;
        read_layer1_result(res, srcdir + *file);
        if (res == nullptr) 
            cout << " error reading file" << endl;
        else {
            cout << "<U";
            update_mean_map(result, res, layer, library, rthresh, gthresh);
            cout << ">";

            delete res;
            cout << " done" << endl;
        }
    }

    
    // Print (DEBUG)
    //for (result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
    //    cout << iter->first << ":";
    //    for (result_t::mapped_type::iterator miter = iter->second.begin(); miter != iter->second.end(); ++miter) {
    //        cout << " (" << miter->first << ") (" << miter->second.get_mean() << ")  ";
    //    }
    //    cout << endl;
    //}

    // For debug purposes
    display_EM_step_distributions(result);

    // Update mean
    update_means(result, layer, library);

    // Save library
    library->save(outlibname);
    library->save_all(change_extension(outlibname, ".png").c_str(), layer + 1, 0, -1, false);
    
    delete library;
}

void perform_optimization_test(const config_dictionary& cfg, const string& fname)
{
    cout << "Doing optimization (test only)" << endl;
    cout << "==============================" << endl;

    int layer;
    string srcdir = cfg.get_value_string("src_dir", "");
    string libname, outlibname;
    list<string> files;
    double cthresh;

    cfg.get_value(layer, "layer", true);
    cfg.get_value(libname, "library", true);
    cfg.get_value(outlibname, "out_library", true);
    cfg.get_value(cthresh, "cover_thresh", true);
   
    part_lib* library;

    read_library(libname, library);

    if (library == nullptr) {
        cout << "Library " << libname << " can not be opened." << endl;
        return;
    }
    string gonamespace = cfg.get_value_string("global_optimization_namespace", "g_optimization");

    optimize_layer(library, layer, 0, cfg, gonamespace, srcdir, fname);

    // Save library
    library->save(outlibname);
    library->save_all(change_extension(outlibname, ".png").c_str(), layer + 1, 0, -1, false);
    
    delete library;
}

void perform_local_optimization_test(const config_dictionary& cfg, const string& fname)
{
    cout << "Doing LOCAL optimization (test only)" << endl;
    cout << "====================================" << endl;

    int layer;
    string srcdir = cfg.get_value_string("src_dir", "");
    string libname, outlibname;
    list<string> files;
    int nparts;
    double sthresh;

    cfg.get_value(layer, "layer", true);
    cfg.get_value(libname, "library", true);
    cfg.get_value(outlibname, "out_library", true);
  
    part_lib* library;

    read_library(libname, library);

    if (library == nullptr) {
        cout << "Library " << libname << " can not be opened." << endl;
        return;
    }
    optimize_layer_locally(library, layer, 0, cfg, "l_optimization", srcdir, fname);

    // Save library
    library->save(outlibname);
    library->save_all(change_extension(outlibname, ".png").c_str(), layer + 1, 0, -1, false);
    
    delete library;
}


void learn_object_parts(const config_dictionary& cfg, const string& pattern)
{
    obj_learning learner(cfg);

    string dir = cfg.get_value_string("src_dir", "");
    string vdir = cfg.get_value_string("validation_src_dir", "");
    string vfilepos = cfg.get_value_string("positive_validation_files", "");
    string vfileneg = cfg.get_value_string("negative_validation_files", "");
    string libname = cfg.get_value_string("out_library", "olib.plb");
    string libimgfile = cfg.get_value_string("library_image_file", "");
    string vcreatorcfg = cfg.get_value_string("validation_creator", "");
    int rndchoice = cfg.get_value_int("random_select", -1);
    string catname = cfg.get_value_string("category_name", "");
	int maxcount = cfg.get_value_int("file_limit", INT_MAX);

    list<string> files;
    list<string> vfilespos;
    list<string> vfilesneg;
    layer1_result* res = nullptr;
    K_bin bin(12, 2, 4, 7);

    end_dir(dir);
    end_dir(vdir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + pattern);
    if (rndchoice > 0) random_selection(files, rndchoice, cfg.get_value_string("dropped_files", ""));
    if (!vfilepos.empty() || !vfileneg.empty()) {
        file_list_from_file(vfilespos, vfilepos, vdir, cfg.get_value_string("positive_pattern", ""));
        file_list_from_file(vfilesneg, vfileneg, vdir, cfg.get_value_string("negative_pattern", ""));
        if ((!vfilespos.empty() || !vfilesneg.empty()) && !vcreatorcfg.empty()) {
            learner.set_creator(new layern_creator(vcreatorcfg.c_str()));
        }
    }
    cout << "from_file = " << (cfg.get_value_bool("from_file", false) ? "true" : "false") << endl;
    cout << "pattern = " << pattern << endl;
    cout << "src_dir = " << dir << endl;

    // read validation files
    for (list<string>::iterator vfile = vfilespos.begin(); vfile != vfilespos.end(); ++vfile) {
        cout << "Loading positive validation file " << *vfile << endl;
        list<irectangle2> gtrs;

        read_layer1_result(res, vdir + *vfile);
        read_groundtruth(gtrs, vdir + *vfile, catname);
        if (res != nullptr) {
            if (gtrs.empty()) delete res; 
            else learner.add_validation_data(res, gtrs.front()); 
        }
    }
    for (list<string>::iterator vfile = vfilesneg.begin(); vfile != vfilesneg.end(); ++vfile) {
        cout << "Loading negative validation file " << *vfile << endl;

        read_layer1_result(res, vdir + *vfile);
        if (res != nullptr) 
            learner.add_validation_data(res, irectangle2()); 
    }


    // learn objects
	int count = 0;

    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        if (++count > maxcount) break;

        cout << "Processing " << *file << endl;
        read_layer1_result(res, dir + *file);
        if (res) {
            list<irectangle2> gtrs;
            clock_t start, end;
            scmap_t scmap;

            get_sc_map(scmap, res, bin, true);
            read_groundtruth(gtrs, dir + *file, catname);
            start = clock();
            learner.object_from_result(res, scmap, gtrs.empty() ? irectangle2() : gtrs.front());
            end = clock();
            cout << "Processing time: " << (double)(end - start)/CLOCKS_PER_SEC << endl;

            //res_list.push_back(res);
            delete res;
        }
    }
    int howmany = learner.add_to_library(catname);

    cout << howmany << " object parts added." << endl;

    learner.library->save(libname);

	if(!libimgfile.empty()) {
		cout << "Saving lib image: layer " << learner.layer << endl;
		for(int il=1; il<=learner.layer+3; il++) {
			string str;
			str += (string)"-" + il + ".png";
		    cout << "Saving lib image (layer " << il << ") to " << change_extension(libimgfile, str).c_str() << endl;
			learner.library->save_all(change_extension(libimgfile, str).c_str(), il, 0, -1,
		        cfg.get_value_bool("show_labels", false));
		}
		//cout << "Saving lib image (layer " << learner.layer + 3 << ") to " << libimgfile << endl;
		//learner.library->save_all(libimgfile.c_str(), learner.layer + 3, 0, -1, cfg.get_value_bool("show_labels", false));
	}
}

void learn_object_parts2(const config_dictionary& cfg, const string& pattern)
{
    cout << "Doing object part learning (learn_objects2)" << endl << endl;
 //   string vdir = cfg.get_value_string("validation_src_dir", "");
 //   string vfiles = cfg.get_value_string("positive_validation_files", "");
 //   string vfileneg = cfg.get_value_string("negative_validation_files", "");
 //   
 //   string vcreatorcfg = cfg.get_value_string("validation_creator", "");
 //   int rndchoice = cfg.get_value_int("random_select", -1);
 //   string catname = cfg.get_value_string("category_name", "");
 //   int maxcount = cfg.get_value_int("file_limit", INT_MAX);


    o_learning olearner(cfg);
    string dir = cfg.get_value_string("src_dir", "");
    string vdir = cfg.get_value_string("validation_src_dir", "");
    string vfiles = cfg.get_value_string("validation_files", "");
    string libname = cfg.get_value_string("out_library", "olib.plb");
    string catname;
    list<string> files;

    cfg.get_value(catname, "category_name", true);

    // Load train files 1st time

    cout << "Duplet Statistics" << endl;

    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + pattern);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        read_layer1_result(res, dir + *file);
        if (res == nullptr) {
            cout << "File " << *file << " can not be opened!" << endl;
        } else {
            cout << "Processing file " << *file;

            olearner.update_duplet_statistics(res);

            delete res;
            cout << " done." << endl;
        }
    }

    // Validation files

    cout << "Validation sets" << endl;

    // "Self"-validation
    if (cfg.get_value_bool("self_validation", false)) {
        if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
            list_directory(files, dir + pattern);
        for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
            layer1_result* res;

            read_layer1_result(res, dir + *file);
            if (res == nullptr) {
                cout << "File " << *file << " can not be opened!" << endl;
            } else {
                cout << "Processing file " << *file;

                list<irectangle2> gtrs;

                read_groundtruth(gtrs, vdir + *file, catname, ".groundtruth");
                if (gtrs.empty()) {
                    irectangle2 r = node_set_bounding_rectangle(res->shape_nodes[0].begin(), res->shape_nodes[0].end());

                    gtrs.push_back(r);

                }
                olearner.add_to_validation_set(res, gtrs);

                delete res;
                cout << " done." << endl;
            }
        }

    }


    files.clear();
    file_list_from_file(files, vfiles, vdir, cfg.get_value_string("validation_pattern", ""));
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        read_layer1_result(res, vdir + *file);
        if (res == nullptr) {
            cout << "File " << *file << " can not be opened!" << endl;
        } else {
            cout << "Processing file " << *file;
            list<irectangle2> gtrs;

            read_groundtruth(gtrs, vdir + *file, catname, ".groundtruth");

            if (gtrs.empty())  // add a dummy rectangle which can never be hit...
                gtrs.push_back(irectangle2(2*res->x_size(0), 2*res->x_size(0), 3*res->x_size(0), 3*res->y_size(0)));

            olearner.add_to_validation_set(res, gtrs);

            delete res;
            cout << " done." << endl;
        }
    }
    
    // Make models

    cout << "Make models" << endl;

    K_bin bin(12, 2, 4, 7);
    
    files.clear();
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + pattern);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        read_layer1_result(res, dir + *file);
        if (res == nullptr) {
            cout << "File " << *file << " can not be opened!" << endl;
        } else {
            scmap_t scmap;

            cout << "Processing file " << *file;

            get_sc_map(scmap, res, bin, true);
            olearner.learn_models(res, scmap);

            delete res;
            cout << " done." << endl;
        }
        //break;
    }
    olearner.finalize();
    olearner.get_library()->save(libname);
}

void learn_ellipses(const config_dictionary& cfg, const string& pattern)
{
    typedef list<pair<irectangle2, string> > gtr_list_t;

    ell_learning learner(cfg);

    string dir = cfg.get_value_string("src_dir", "");
    string libname, svmdir;
    list<string> files;

    cfg.get_value(libname, "out_library", true);
    cfg.get_value(svmdir, "svm_dir", true);
    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + pattern);

    // read validation files
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        cout << "Processing " << *file;

        read_layer1_result(res, dir + *file);
        if (res == nullptr) 
            cout << " error reading file." << endl;
        else {
            gtr_list_t gtrs;

            read_groundtruth(gtrs, dir + *file, "", ".groundtruth", true);
            learner.update(res, gtrs);
            cout << " done." << endl;
            delete res;
        }
    }
	
	learner.learn_all(libname + ".xml");
    learner.get_library()->save(libname);

    cout << "End" << endl;
    
}

void learn_shapes(const config_dictionary& cfg, const string& pattern)
{
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");

    string dir = cfg.get_value_string("src_dir", "");
    string srclib, outlib;
    int category_layer;
    string category;
    double cover_threshold = cfg.get_value_double("cover_threshold", 0.9);
    list<string> files;

    cfg.get_value(srclib, "library", true);
    cfg.get_value(outlib, "out_library", true);
    cfg.get_value(category_layer, "category_layer", true);
    cfg.get_value(category, "category", true);
    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + pattern);

    map<int, shape_model_maker> makermap;
    part_lib* library;
    
    read_library(srclib, library);
    
    if (library == nullptr) {
        cout << "Can not find library '" << srclib << "'" << endl;
        return;
    }

    set<int> types;
    int layer = category_layer - 1;

    library->get_object_parts(types, category_layer, category);

    // read files
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        cout << "Processing " << *file;

        read_layer1_result(res, dir + *file);
        if (res == nullptr) 
            cout << " error reading file." << endl;
        else {
            if (layer < 0 || layer > res->max_layer_index()) 
                cout << " Layer " << layer << " does not exist!" << endl;
            else {
                set<ipoint2> p0;
            
                node_set_to_point_set(p0, res->shape_nodes[0].begin(), res->shape_nodes[0].end());
                for (vector<node*>::iterator niter = res->shape_nodes[layer].begin(); 
                        niter != res->shape_nodes[layer].end(); ++niter) {
                    node* n = *niter;
                    layer1_data* nd = (layer1_data*)n->data;
                    set<node*> supp;
                    set<ipoint2> psupp;

                    if (types.find(nd->m) == types.end())
                        continue;

                    //cout << "L<";
                    link_path(res, n, toprev, to0);
                    //cout << ">";
                    n->get_neighbor_set(to0, supp);
                    node_set_to_point_set(psupp, supp.begin(), supp.end());

                    double q = (double)psupp.size()/p0.size();

                    if (q >= cover_threshold) {
                        ipoint2 mean;
                        double scale;

                        ipoint2_set_normalization_parameters(mean, scale, psupp.begin(), psupp.end());
                        foreach_neighbor(n, to0, nniter) {
                            edge_path_data_t* ned = (edge_path_data_t*)neighbor_edge_data(nniter);
                            layer1_data* nnd = (layer1_data*)neighbor_node_data(nniter);
                            ipoint2 p((int)(100*(nnd->x - mean.x)/scale), (int)(100*(nnd->y - mean.y)/scale));

                            makermap[nd->m].update(ned->data.edges(), p);
                        }
                    }
                    
                }
            }
            cout << endl;
            delete res;
        }
    }

    vector<int> todelete;

    for (int pi = 0; pi < (int)library->parts[layer].size(); ++pi) {
        map<int, shape_model_maker>::iterator miter = makermap.find(pi);
        //cout << "***   ";

        if (miter == makermap.end()) {
            if (types.find(pi) != types.end()) 
                todelete.push_back(pi);
        } else {
            node* p = library->parts[layer][pi];
            //cout << "!!!!  ";
            part_data* pd = (part_data*)p->data;
            spart_data* spd = new spart_data(*pd, miter->second.get_model());
            string name = string("part") + pi + string(".m");

            cout << name << endl;
            spd->get_model().save_mathematica(name);
            p->data = spd;
            delete pd;
        }
    }
    library->delete_parts(layer, todelete);
    library->save(outlib);

    cout << todelete.size() << " models deleted." << endl;

    delete library;
    
}

void set_s_thresholds(const config_dictionary& cfg, const string& pattern)
{
    part_lib* library;
    int layer = cfg.get_value_int("layer", -1);

    read_library(pattern, library);

    if (library == nullptr) {
        cout << "Library " << pattern << " can not be opened." << endl;
        return;
    }

    if (layer < 0 || layer > library->max_layer_index()) {
        cout << "Layer " << layer << " does not exist in library" << endl;
        delete library;
        return;
    }

    for (int i = 0; i < library->layer_size(layer); ++i) {
        node* p = library->parts[layer][i];
        spart_data* pd = dynamic_cast<spart_data*>(p->data);
        
        if (pd == nullptr) {
            cout << "Part " << i << " is not spart_data, skipping" << endl;
            continue;
        }

        string catname = library->get_category(p);

        if (catname == "") {
            cout << "Category name empty, skipping" << endl;
            continue;
        }
                
        double val = cfg.get_value_double(catname, -1.0);

        if (val >= 0.0) {
            pd->td.set_thresh(S_RESPONSE, val);
        }
    }

    library->save(cfg.get_value_string("out_library", pattern));
    delete library;
}


void write_matching(const string& fname, vector<dpoint2> srcv, vector<dpoint2> destv)
{
    ofstream os(fname.c_str());
    
    translate_and_scale(srcv);
    translate_and_scale(destv);
    vector<int> perm = point_matching(srcv, destv);
    permute(srcv, perm);

    for (int i = 0; i < srcv.size(); ++i) {
        os << srcv[i].x << ',' << srcv[i].y << ',' << destv[i].x << ',' << destv[i].y << endl;
    }

    os.close();
}

cv::Mat back_projection(const cv::Mat& data, const clique_data& cd)
{
    cv::Mat coeffs;
    cv::Mat result;
    
    gemm(data - cd.mean, cd.eigenvectors, 1, cv::Mat(), 0, coeffs, cv::GEMM_2_T);
    gemm(coeffs, cd.eigenvectors, 1, cd.mean, 1, result, 0);
    return result;
}

// Greedy selects cliques: in each step select the biggest clique C with the smallest intersection
// percent with already selected parts P; if |P \cap C|/|C| < thresh then intersection is considered 0.
// When clique C is selected, parts which can be explained sufficiently well with PCA 
void reduce_clique_set(part_lib* library, int layer, vector<clique_data>& cliques, double thresh)
{
    typedef tuple<double, int, int> int_item_t;

    vector<vector<dpoint2> > geovec(library->layer_size(layer), vector<dpoint2>());

    // Get part geometry for each part in library 
    for (int pi = 0; pi < library->layer_size(layer); ++pi) { 
        node* p = library->parts[layer][pi];
        vs_part_data* vspd = dynamic_cast<vs_part_data*>(p->data);
        path_map_t pm;
        vector<ipoint2> pts;

        if (vspd != nullptr) 
            geovec[pi] = partition(vspd->pcad.mean);
        else {
            get_library_geo(pm, p);
            for (path_map_t::iterator iter = pm.begin(); iter != pm.end(); ++iter) {
                pts.push_back(iter->second.p);
            }
            geovec[pi] = cast_vector<dpoint2, ipoint2>(inhibit_point_set(pts, 2));
        }
    }
    // Calculate distances between parts and cliques; cliquedist(p, c) is a distance of part p to clique c.
    matrix<ddpair> cliquedist(library->layer_size(layer), (int)cliques.size());

    for (int pi = 0; pi < library->layer_size(layer); ++pi) {
        vector<dpoint2> gv = geovec[pi];

        for (int ci = 0; ci < (int)cliques.size(); ++ci) {
            //cout << pi << '-' << ci << ' ';
            vector<dpoint2> v = partition(cliques[ci].mean);
            vector<dpoint2> gvresized = get_resized_vector(gv, (int)v.size());

            translate_and_scale(gvresized);

            vector<int> perm = point_matching(gvresized, v);

            //write_matching(string("c:\\work\\tmp\\match") + pi + string("-") + ci + string(".m"), gvresized, v);
            
            permute(gvresized, perm);
            
            cv::Mat data = flatten(gvresized);
            cv::Mat bproj = back_projection(data, cliques[ci]);
            double erfv = erf_value(data, cliques[ci].mean, cliques[ci].eigenvectors, cliques[ci].eigenvalues);

            cliquedist(pi, ci).first = cv::norm(bproj, data, cv::NORM_L2);
            cliquedist(pi, ci).second = 1 - erfv;
        }
    }

    vector<ddpair> maxintrad(cliques.size(), ddpair(0.0, 0.0));

    for (int i = 0; i < (int)cliques.size(); ++i) {
        for (auto piter = cliques[i].parts.begin(); piter != cliques[i].parts.end(); ++piter) {
            if (cliquedist(*piter, i).first > maxintrad[i].first)
                maxintrad[i].first = cliquedist(*piter, i).first;
            if (cliquedist(*piter, i).second > maxintrad[i].second)
                maxintrad[i].second = cliquedist(*piter, i).second;
        }
    }

    set<int> tocover;   // remaining parts; initially a union of all cliques
    set<int> covered;   // already covered parts

    for (auto cliter = cliques.begin(); cliter != cliques.end(); ++cliter) 
        tocover.insert(cliter->parts.begin(), cliter->parts.end());

    set<int> tokeep;    // indices of cliques to keep

    while (!tocover.empty()) {
        vector<int_item_t> intersections;

        for (int i = 0; i < (int)cliques.size(); ++i) {
            if (tokeep.find(i) != tokeep.end()) 
                continue;

            int s = (int)cliques[i].parts.size();
            double f = (double)intersection_size(cliques[i].parts, covered)/s;

            if (f < thresh) f = 0.0;
            intersections.push_back(int_item_t(f, -s, i));
        }

        if (intersections.empty()) 
            break;

        sort(intersections.begin(), intersections.end());

        int selc = get<2>(intersections[0]);

        // Change 'selc' clique -- add parts to it.
        for (auto piter = tocover.begin(); piter != tocover.end(); ++piter) {
            if (cliquedist(*piter, selc).first < 1.1*maxintrad[selc].first && 
                    cliquedist(*piter, selc).second < 1.1*maxintrad[selc].second)
                cliques[selc].parts.insert(*piter);
        }

        tokeep.insert(selc);
        covered.insert(cliques[selc].parts.begin(), cliques[selc].parts.end());

        set<int> tmps;

        set_difference(tmps, tocover, cliques[selc].parts);
        tocover = tmps;
    }
 
    for (int i = (int)cliques.size() - 1; i >= 0; --i)
        if (tokeep.find(i) == tokeep.end()) {
            cliques.erase(cliques.begin() + i);
        }
}

// Greedy selection of cliques; cliques with all elements included in previous 
// cliques are removed. If 'sortvector' is true, cliques are first sorted by their
// size.
void reduce_clique_set_by_inclusion(vector<clique_data>& cliques, bool sortvector)
{
    if (sortvector) {
        sort(cliques.begin(), cliques.end(),
            [](const clique_data& c1, const clique_data& c2) { return c1.parts.size() > c2.parts.size(); }
        );
    }
    
    auto citer = cliques.begin();
    set<int> selected;

    while (citer != cliques.end()) {
        if (includes(selected.begin(), selected.end(), citer->parts.begin(), citer->parts.end())) 
            citer = cliques.erase(citer);
        else {
            selected.insert(citer->parts.begin(), citer->parts.end());
            ++citer;
        }
    }
}

void reduce_clique_set_pca(part_lib* library, int layer, vector<clique_data>& cliques, double norm)
{
    vector<vector<dpoint2> > geovec(library->layer_size(layer), vector<dpoint2>());

    // 1. Get part geometry for each part in library 
    for (int pi = 0; pi < library->layer_size(layer); ++pi) { 
        node* p = library->parts[layer][pi];
        vs_part_data* vspd = dynamic_cast<vs_part_data*>(p->data);
        path_map_t pm;
        vector<ipoint2> pts;

        if (vspd != nullptr) 
            geovec[pi] = partition(vspd->pcad.mean);
        else {
            get_library_geo(pm, p);
            for (path_map_t::iterator iter = pm.begin(); iter != pm.end(); ++iter) {
                pts.push_back(iter->second.p);
            }
            geovec[pi] = cast_vector<dpoint2, ipoint2>(inhibit_point_set(pts, 2));
        }
    }

    // 2. Calculate norm of part geometry to its "back"-projection from clique space -- for each part and each clique
    //    matrix(p, c)
    matrix<double> cliquedist(library->layer_size(layer), (int)cliques.size());

    for (int pi = 0; pi < library->layer_size(layer); ++pi) {
        vector<dpoint2> gv = geovec[pi];

        for (int ci = 0; ci < (int)cliques.size(); ++ci) {
            //cout << pi << '-' << ci << ' ';
            vector<dpoint2> v = partition(cliques[ci].mean);
            double sizef = max((double)v.size()/gv.size(), (double)gv.size()/v.size());
            vector<dpoint2> gvresized = get_resized_vector(gv, (int)v.size());

            sizef = sqr(sizef);
            translate_and_scale(gvresized);

            vector<int> perm = point_matching(gvresized, v);

            //write_matching(string("c:\\work\\tmp\\match") + pi + string("-") + ci + string(".m"), gvresized, v);
            
            permute(gvresized, perm);
            
            cv::Mat data = flatten(gvresized);

            cliquedist(pi, ci) = subspace_distance(data, cliques[ci].mean, cliques[ci].eigenvectors,
                cliques[ci].eigenvalues, 2.0)*sizef;
        }
    }

    // Remove redundant cliques based on difference on back-projection
    vector<int> toremove;

    while (true) {

        // 3. Clique D is absorbed in clique C if D can be expressed with basis of C with small error
        vector<vector<int> > nabsorbed(cliques.size(), vector<int>());   // nabsorbed[i] counts cliques absorbed in clique i

        for (int ci = 0; ci < (int)cliques.size(); ++ci) {
            if (cliques[ci].parts.empty())
                continue;

            for (int di = 0; di < (int)cliques.size(); ++di) {
                const set<int>& D = cliques[di].parts;
                double maxdiff = 0.0;

                if (D.empty())
                    continue;

                for (set<int>::const_iterator iter = D.begin(); iter != D.end(); ++iter) {
                    int e = *iter;
                    double diff = cliquedist(e, ci);

                    if (diff > maxdiff) maxdiff = diff;
                }
                //cout << di << " is absorbed in " << ci << " with " << maxdiff << endl;
                if (ci != di && maxdiff <= norm) 
                    nabsorbed[ci].push_back(di);
                
            }
        }

        // 4. Find clique which absorbes most cliques
        vector<vector<int> >::iterator miter = max_element(nabsorbed.begin(), nabsorbed.end(), size_less<int>);

        // 4a. If there is no clique which absorbes some other clique, quit
        if (miter == nabsorbed.end() || miter->size() == 0)
            break;

        // 5. Keep such clique and remove absorbed cliques
        int keepi = (int)(miter - nabsorbed.begin());
        vector<int>& torm = *miter;

        for (vector<int>::iterator riter = torm.begin(); riter != torm.end(); ++riter) {
            cliques[keepi].parts.insert(cliques[*riter].parts.begin(), cliques[*riter].parts.end());
            toremove.push_back(*riter);
            cliques[*riter].parts.clear();
        }
    }
    sort(toremove.begin(), toremove.end(), greater<int>());
    for (vector<int>::iterator riter = toremove.begin(); riter != toremove.end(); ++riter) 
        cliques.erase(cliques.begin() + *riter);
}

void add_vs_part(part_lib* library, int layer, const clique_data& clique)
{
    int lyrSimilar = atom("lyrSimilar");


    // Add new node representing clique (vs-node)
    vs_part_data* pd = new vs_part_data();

    pd->cmx = 0;
    pd->cmy = 0;
    pd->layer = layer;
    pd->type = (int)library->parts[layer].size();
    pd->pcad.mean = clique.mean;
    pd->pcad.eigenvalues = clique.eigenvalues;
    pd->pcad.eigenvectors = clique.eigenvectors;
    pd->pcad.sizefactor = clique.sizefactor;
    
    node* n = library->add_node(pd, VS_PART_ATTR);
    set<node*> members;

    library->parts[layer].push_back(n);

    // VS part  <--- lyrVSRoot ------ "ordinary part"
    //          ---- lyrVSMember ---> 

    for (set<int>::const_iterator iter = clique.parts.begin(); iter != clique.parts.end(); ++iter)  {
        node* nn = library->parts[layer][*iter];

        members.insert(nn);
        foreach_neighbor (nn, lyrSimilar, nniter) {
            members.insert(neighbor_node(nniter));
        }
    }
    for (set<node*>::iterator miter = members.begin(); miter != members.end(); ++miter) {
        library->add_edge_unique(n, *miter, atom("lyrVSMember"), atom("lyrVSRoot"));
    }

}

void add_vs_parts(part_lib* library, int layer, const vector<clique_data>& cliques)
{
    // Add cliques as new parts...
    set<int> cparts;
    int layersize = (int)library->parts[layer].size();

    cout << "Adding parts >= " << layersize << endl;

    for (int ci = 0; ci < (int)cliques.size(); ++ci) {
        add_vs_part(library, layer, cliques[ci]);
        cparts.insert(cliques[ci].parts.begin(), cliques[ci].parts.end());
    }

    // Add parts which are not members of any clique as "singleton" cliques...
    for (int i = 0; i < layersize; ++i) {
        node* p = library->parts[layer][i];
        vs_part_data* vspd = dynamic_cast<vs_part_data*>(p->data);

        if (cparts.find(i) != cparts.end() || !is_sim_root(p))
            continue;
        
        clique_data cd;

        cd.part0 = i;
        cd.parts.insert(i);

        if (vspd != nullptr) {
            cd.mean = vspd->pcad.mean;
            cd.eigenvectors = vspd->pcad.eigenvectors;
            cd.eigenvalues = vspd->pcad.eigenvalues;
            cd.sizefactor = vspd->pcad.sizefactor;
        } else {
            path_map_t pm;
            vector<ipoint2> pts;

            get_library_geo(pm, p);
            for (path_map_t::iterator iter = pm.begin(); iter != pm.end(); ++iter) {
                pts.push_back(iter->second.p);
            }
            vector<dpoint2> ptsd = cast_vector<dpoint2, ipoint2>(inhibit_point_set(pts, 2));
            double sizef = translate_and_scale(ptsd).first;

            cd.mean = flatten(ptsd);
            cd.eigenvectors = cv::Mat(1, cd.mean.cols, CV_64F, cv::Scalar(0.0));
            cd.eigenvalues = cv::Mat(1, 1, CV_64F, cv::Scalar(0.0));
            cd.sizefactor = sizef;
        }
        add_vs_part(library, layer, cd);
    }
}

template<class T, class I> void all_ordered_pairs(list<pair<T, T> >& ps, I begin, I end)
{
    ps.clear();

    for (I i1 = begin; i1 != end; ++i1) {
        for (I i2 = begin; i2 != end; ++i2) {
            ps.push_back(pair<T, T>(*i1, *i2));
        }
    }
}

void save_matching(const string& fname, const vector<int>& perm, const vector<dpoint2>& pts1, 
    const vector<dpoint2>& pts2)
{
    ofstream os(fname.c_str());

    for (int i = 0; i < perm.size(); ++i) {
        os << (perm[i] + 1) << ',' << pts1[i].x << ',' << pts1[i].y << ',' << pts2[i].x << ',' << pts2[i].y << endl;
    }
    os.close();
}


// Add "lyrSimVal" edges between parts of the same layer 
void do_add_similarity_edges(const config_dictionary& cfg, const char* fname)
{
    unique_ptr<part_lib> library(read_library(fname));

    if (library == nullptr) {
        cout << "Can not open library file '" << fname << "'" << endl;
        return;
    }

    vector<int> layers;
    int max_geo_size = 20;

    cfg.get_value(layers, "layers", true);
    for (int i = 0; i < layers.size(); ++i) {
        add_similarity_edges(library.get(), layers[i], max_geo_size, 0.0);
    }
    library->save(cfg.get_value_string("out_library", fname));
}

string format_string(const string& s, const string& a)
{
    char buf[1000];

    sprintf(buf, s.c_str(), a.c_str());
    return string(buf);
}

void test(const config_dictionary& cfg, const char* pattern)
{
    cout << "Testing library operations." << endl;
    part_lib* library;
    double benergy;
    vector<dpoint2> dvector;
    double scdistance;

    read_library(pattern, library);
    if (library == nullptr) {
        cout << "Library can not be opened!" << endl;
        return;
    }

    vector<node*>& parts = library->parts[3];

    for (int i = 0; i < 1 /*(int)parts.size()*/; ++i) {
        node* p = parts[i];
        vector<int> cluster;

        cluster.push_back(i);
        cout << "Part " << i << ": ";
        while ((int)cluster.size() < 5) { 
            double maxminscdistance = 0.0;
            int minmaxj = -1;

            for (int j = i + 1; j < (int)parts.size(); ++j) {
                node* q = parts[j];
                path_map_t qpm;
                double minscdistance = 1e5;

                get_library_geo(qpm, q);
                for (int c = 0; c < (int)cluster.size(); ++c) {
                    node* cp = parts[cluster[c]];
                    path_map_t cpm;

                    get_library_geo(cpm, cp);
                    part_geometry_matching(benergy, dvector, scdistance, qpm, cpm, false);
                    if (scdistance < minscdistance) minscdistance = scdistance;
                }
                if (minscdistance > maxminscdistance) {
                    maxminscdistance = minscdistance;
                    minmaxj = j;
                }
            }
            cout << "to cluster: " << minmaxj << endl;
            cluster.push_back(minmaxj);
        }
    }

    delete library;
}

void old_test(const config_dictionary& cfg, const char* pattern)
{
    map_learning ml(cfg);

    string dir = cfg.get_value_string("src_dir", "");
    int part_max_number = cfg.get_value_int("part_max_number", 20);
    string libname;

    list<string> files;
    part_lib* library;

    cfg.get_value(libname, "library", true);
    read_library(libname, library);

    if (library == nullptr) {
        cout << "Can not open library '" << libname << "'" << endl;
        return;
    }

    ml.reset(library);

    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + pattern);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        read_layer1_result(res, dir + *file);
        if (res != nullptr) {
            cout << "Processing " << *file << " ...";

            ml.prepare_for_update(res);
            ml.update(res);

            cout << " done." << endl;

            delete res;
        }
    }
    ml.display_statistics("c:\\temp\\map_%03d-%03d.png");

    // part learning
    part_learning pl(cfg);

    // find maxima
    pl.reset(ml);
    pl.display_maxima("c:\\temp\\maxima_%03d-%03d.png");

    // make "K-bin"
    K_bin bin(12, 2, 4, 7);

    // update part statistics
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        read_layer1_result(res, dir + *file);

        if (res != nullptr) {
            cout << "Processing " << *file << " ...";

            scmap_t scmap;
            clock_t start, end;

            start = clock();
            get_sc_map(scmap, res, bin, cfg.get_value_bool("normalize_histogram", false));
            end = clock();

            cout << " (sc " << (double)(end - start)/CLOCKS_PER_SEC << " sec)";

            pl.update(res, scmap);

            cout << " done." << endl;

            delete res;
        }
    }
    cout << "Adding to library..." << endl;
    pl.add_to_library(library, part_max_number, 0,
        cfg.get_value_int("cluster_size", 5),
        cfg.get_value_double("merge_distance_threshold", 2),
        cfg.get_value_double("merge_sc_threshold", 1));
    library->save_all_sc("c:\\temp\\lib.png", pl.get_source_layer() + 2, 0, -1, cfg.get_value_bool("show_labels", false));
    library->save(cfg.get_value_string("out_library", "lib2new.plb"));

    cout << "Performing matching..." << endl;
    matrix<double> m;

    //pl.similarity_matrix(m, part_max_number, 0.0 /* dummy */);
    //m.print();
    
    //cout << "Mathematica output..." << endl;
    //pl.export_stat("c:\\temp\\lib.m", part_max_number);

    
    cout << "Merge parts..." << endl;

    vector<vector<part_learning::similarity_item> > mresult;

    pl.merge_stat_sc(mresult, part_max_number, 
        cfg.get_value_int("cluster_size", 5),
        cfg.get_value_int("merge_space_size", INT_MAX), 0,
        cfg.get_value_double("merge_distance_threshold", 2),
        cfg.get_value_double("merge_sc_threshold", 1));

    for (int i = 0; i < (int)mresult.size(); ++i) {
        cout << "Cluster " << i << ":";
        for (vector<part_learning::similarity_item>::iterator siter = mresult[i].begin(); siter != mresult[i].end(); ++siter) {
            cout << " " << siter->index;
        }
        cout << endl;
    }


    delete library;
}

void do_finalize_library(const config_dictionary& cfg, const char* pattern)
{
    string outlibname;
    part_lib* library;

    cfg.get_value(outlibname, "out_library", true);

    read_library(pattern, library);

    if (library == nullptr) return;

    library->delete_unused_parts();
    library->save(outlibname);
    
    delete library;
}

void do_validate_parts(const config_dictionary& cfg, const char* pattern)
{

    string dir = cfg.get_value_string("src_dir", "");
    list<string> files;
    string libname, outlibname;
    layer1_result* res;
    part_lib* library;
    list<irectangle2> gtruths;
    int object_index, object_layer;
    bool keep_nonexistent;
    double thresh;
    double ratio;
    string extension;
    validation_result_t result;
    set<int> parts;

    keep_nonexistent = cfg.get_value_bool("keep_nonexistent", true);
    cfg.get_value(object_layer, "object_layer", true);
    cfg.get_value(object_index, "object_index", true);
    cfg.get_value(thresh, "threshold", true);
    cfg.get_value(ratio, "tf_ratio", true);
    cfg.get_value(libname, "library", true);
    outlibname = cfg.get_value_string("out_library", "");
    cfg.get_value(extension, "extension", true);
    read_library(libname, library);

    cout << "Validating object #" << object_index << " on layer " << object_layer << endl;
    library->get_parent_parts(parts, object_layer, object_index);
    end_dir(dir);  
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + pattern);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cout << "Processing " << *file << endl;
        read_layer1_result(res, dir + *file);
        read_groundtruth(gtruths, dir + *file);
        if (res != nullptr) {
            if (!gtruths.empty()) 
                validate_parts(result, res, parts, gtruths, object_layer - 1, thresh);
            delete res;
        }
    }

    vector<int> to_delete;

    cout << "Final results after the validation" << endl;
    for (validation_result_t::iterator riter = result.begin(); riter != result.end(); ++riter) {
        cout << "Part #" << riter->first << " T: " << riter->second.first << "  F: " 
            << riter->second.second;
    
        if (riter->second.first != 0 || riter->second.second != 0) {
            if (riter->second.first == 0 || riter->second.second != 0 &&
                (double)(riter->second.first)/riter->second.second < ratio) {
                    to_delete.push_back(riter->first);
                    cout << " (D)";
            }
        } else if (riter->second.first == 0 && riter->second.second == 0 && !keep_nonexistent) {
            to_delete.push_back(riter->first);
            cout << " (D)";
        }

        cout << endl;
    }

    if (outlibname != "") {
        library->delete_parts(object_layer - 1, to_delete);
        library->save(outlibname);
    } else {
        string out = cfg.get_value_string("out_file", "to_delete.txt");
        ofstream os(out.c_str(), ios::out | ios::app);

        for (vector<int>::iterator iter = to_delete.begin(); iter != to_delete.end(); ++iter) {
            os << *iter << ' ';
        }
        os << endl;
        os.close();
    }

    delete library;

}

void make_nbsets(map<set<int>, int>& result, const vector<int>& partsizes, layer1_result* res, int layer, int maxsize)
{
    typedef map<set<int>, int> result_t;

    vector<node*>& s_nodes = res->shape_nodes[layer];

    for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        node* n = *iter;
        int bestsize = -1;
        vector<int> bestm;

        while (n != nullptr) {
            layer1_data* nd = (layer1_data*)n->data;
           
            if (partsizes[nd->m] == bestsize)
                bestm.push_back(nd->m);
            else if (partsizes[nd->m] > bestsize) { 
                bestsize = partsizes[nd->m];
                bestm.clear();
                bestm.push_back(nd->m);
            }
            n = nd->next;
        }

        set<int> bestset(bestm.begin(), bestm.begin() + min<size_t>((size_t)maxsize, bestm.size()));
        pair<result_t::iterator, bool> ibpair = result.insert(result_t::value_type(bestset, 1)); 

        if (!ibpair.second) ++(ibpair.first->second);
    }
}

void greedy_min_intersection(vector<iipair>& result, map<set<int>, int>& nbsets)
{
    cout << "greedy_min_intersection does not work properly!" << endl;
    throw exception();

    typedef set<int> iset;
    typedef iset* iset_p;
    typedef triple<const iset*, int, bool> mapped_t;
    typedef map<iset, int> nbsets_t;
    typedef multimap<int, mapped_t*> tsmap_t;

    tsmap_t tsmap;
    set<int> active;
    list<mapped_t> mapped;

    // Initializations
    for (nbsets_t::iterator siter = nbsets.begin(); siter != nbsets.end(); ++siter) {
        mapped.push_back(mapped_t(&(siter->first), siter->second, true));
        for (iset::const_iterator iter = siter->first.begin(); iter != siter->first.end(); ++iter) {
            tsmap.insert(tsmap_t::value_type(*iter, &mapped.back()));
            active.insert(*iter);
        }
    }
    
    // Greedy algorithm
    while (!active.empty()) {
        iset::iterator aiter = active.begin(); 
        int bestcount = 0;
        int bestactive;
        
        while (aiter != active.end()) {
            int i = *aiter;
            tsmap_t::iterator tsiter = tsmap.find(i);
            int ucount = 0;

            while (tsiter != tsmap.end() && tsiter->first == i) {
                if (tsiter->second->third) ucount += tsiter->second->second;
                ++tsiter;
            } 
            if (ucount == 0) active.erase(aiter);
            else {
                if (ucount > bestcount) { bestcount = ucount; bestactive = i; }
                ++aiter;
            }
        }
        if (bestcount > 0) {
            result.push_back(iipair(bestactive, bestcount));
            for (tsmap_t::iterator tsiter = tsmap.find(bestactive); tsiter != tsmap.end() && tsiter->first == bestactive; ++tsiter) {
                tsiter->second->third = false;
            }
        }
    }
}

void do_optimize_library(config_dictionary& cfg, const char* pattern)
{
    typedef vector<iipair> result_t;
    typedef map<set<int>, int> nbsets_t;

    cout << "doing do_optimize_library (TESTING)" << endl;

    string libname;
    string dir = cfg.get_value_string("src_dir", "");
    int maxf = cfg.get_value_int("max_files", INT_MAX);
    int maxnbset = cfg.get_value_int("max_local_set", INT_MAX);
    int layer;
    list<string> files;
    int lyrsrc = atom("lyrSrc");

    cfg.get_value(libname, "library", true);
    cfg.get_value(layer, "layer", true);

    part_lib* library;

    read_library(libname, library);

    if (library == nullptr) {
        cout << "Library can not be loaded!" << endl;
        return;
    }
    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + pattern);

    // Count part sizes
    vector<int> partsizes(library->parts[layer].size(), 0);

    for (int i = 0; i < (int)partsizes.size(); ++i) {
        partsizes[i] = library->parts[layer][i]->count_neighbors(lyrsrc) + 1;
    }

    // Read layer1_result files and make "nbsets"
    nbsets_t nbsets;

    int fcount = 0;

    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        if (fcount >= maxf) break;

        layer1_result* res;

        cout << "Processing: " << *file;
        read_layer1_result(res, dir + *file);
        if (res != nullptr) {
            make_nbsets(nbsets, partsizes, res, layer, maxnbset);
            delete res;
            fcount++;
            cout << " done";
        }
        cout << endl;
    }
    
    // Greedy selection of parts 
    result_t result;

    greedy_min_intersection(result, nbsets);

    // Print the result
    cout << "Part selection result" << endl;
    cout << result << endl;

    // Keep only best parts & save
    if (!result.empty()) {
        double thresh = cfg.get_value_double("keep_parts_threshold", 0.1);
        vector<int> tokeep;
        int min = (int)(result[0].second * thresh);
        string fname;

        for (result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
            if (iter->second >= min) tokeep.push_back(iter->first);
        }
        library->keep_parts(layer, tokeep);
        
        fname = cfg.get_value_string("image_file", "");
        if (!fname.empty()) 
            library->save_all(fname.c_str(), layer + 1, 0, -1, cfg.get_value_bool("show_labels", true));
        fname = cfg.get_value_string("out_library", "");
        if (!fname.empty()) 
            library->save(fname);
    }
    
    // Inferrence
    //lc.allowed_parts

    delete library;
}


void do_learn_thresholds(const config_dictionary& cfg, const char* pattern)
{
    t_learning tl(cfg);

    string pdir = cfg.get_value_string("positive_src_dir", "");
    string ndir = cfg.get_value_string("negative_src_dir", "");
    string pfile, nfile, libfile;
    list<string> pfiles, nfiles;
    vector<string> categories;
	vector<string> threshold_types;

    cfg.get_value(libfile, "out_library", true);
    cfg.get_value(pfile, "positive_files", true);
    cfg.get_value(nfile, "negative_files", true);
    cfg.get_value(categories, "categories", true);
    cfg.get_value(threshold_types, "threshold_types", true);

    file_list_from_file(pfiles, pfile, pdir, cfg.get_value_string("positive_pattern", ""));
    file_list_from_file(nfiles, nfile, ndir, cfg.get_value_string("negative_pattern", ""));

	int thrtype;
    tl.begin();
	for (vector<string>::iterator iteT = threshold_types.begin(); iteT != threshold_types.end(); ++iteT) {
		if(!(*iteT).compare("r"))
			thrtype= R_THRESH;
		else if(!(*iteT).compare("g"))
			thrtype= G_THRESH;
		else if(!(*iteT).compare("rr"))
			thrtype= RR_THRESH;
        else if(!(*iteT).compare("s"))
            thrtype= S_THRESH;
		else {
			cout << "Invalid threshold type: " << (*iteT) << endl;
			cout << "Please use lowercase comma separated threshold type names (e.g. g,rr)." << endl;
			continue;
		}
		cout << "Learning threshold: " << (*iteT) << endl;
		for (vector<string>::iterator iter = categories.begin(); iter != categories.end(); ++iter) {
			// check if name is number (may be old version where numbers were used for categories)
			// so inform user to change it to string
			try {
				int i = Convert::ToInteger<int>((*iter).c_str(), (unsigned)(*iter).size());
				// display warning to user
				cout << "WARNING: Found integer for category name instead of string (integer value: " << i <<  ") in threshold learning !!" << endl;
				cout << "\t\tYou may be using old version of config where integers were used for category definition." << endl;
				cout << "\t\tYou should be using names instead of numbers !!" << endl;
			} catch (InvalidConversionException& ex) {
				// apparently it is not a number, continue
			}
			cout << "Learning category: " << (*iter) << endl;
            tl.learn(pdir, pfiles, ndir, nfiles, *iter, thrtype, thrtype == S_THRESH ? -1 : 1);
		}
    }
    tl.end();
    tl.get_library()->save(libfile);

	string libimgfile = cfg.get_value_string("library_image_file", "");
	bool showlabels = cfg.get_value_bool("show_labels", false);
	int nlayers = tl.get_library()->layer_count;
	if(!libimgfile.empty()) {
		cout << "Saving lib image to " << libimgfile << " layers " << nlayers << endl;
		for(int il=1; il<nlayers; il++) {
			string str;
			str += (string)"-" + il + ".png";
		    cout << "Saving lib image (layer " << il << ") to " << change_extension(libimgfile, str).c_str() << endl;
			tl.get_library()->save_all(change_extension(libimgfile, str).c_str(), il, 0, -1, showlabels);
		}
	}
}


HOPINTERFACE_API bool hop_learning_tools(const char* cfg_file, const char* patt, const char* params)
{
	bool ok = true;
    try {
		// all config values
        config_dictionary cfg_all(cfg_file);
		cfg_all.from_string(params);
		cfg_all.update_namespace_references();

		// get only 'learning' namespace (also with root values)
		config_dictionary cfg_learn;
		cfg_learn.from_namespace_priority(cfg_all, 1, "learning");

        //string action = cfg_learn.get_value_string("action", "update");
		string action = cfg_learn.get_value_string("action", "optimization");

		// then get only values needed for specific action
		config_dictionary cfg_action;
		cfg_action.from_namespace_priority(cfg_learn, 1, action.c_str());

        //if (action == "update") do_update(cfg_action, patt); 
        //else if (action == "learn_parts") do_find_parts(cfg_action, patt);
        //else if (action == "merge_parts") merge_library(cfg_action);
        //??else if (action == "keep_parts") keep_parts(cfg_action);
        //??else if (action == "delete_parts") delete_parts(cfg_action);
        if (action == "optimization") perform_optimizations(cfg_action, patt);
        //??else if (action == "optimization_from_images") from_images(cfg_action, patt);
        else if (action == "learn_objects") learn_object_parts(cfg_action, patt);
        else if (action == "learn_objects2") learn_object_parts2(cfg_action, patt);
        //else if (action == "validate_objects") do_validate_parts(cfg_action, patt);
        // ?? else if (action == "finalize_library") do_finalize_library(cfg_action, patt);
        // ?? else if (action == "optimize_library") do_optimize_library(cfg_action, patt);
        // ?? else if (action == "learn_thresholds") do_learn_thresholds(cfg_action, patt);
        //else if (action == "learn_g_distribution") learn_g_distribution(cfg_action, patt);
        //else if (action == "EM") perform_EM(cfg_action, patt);
        //else if (action == "optimization_test") perform_optimization_test(cfg_action, patt);
        //else if (action == "local_optimization_test") perform_local_optimization_test(cfg_action, patt);
        //else if (action == "learn_ellipses") learn_ellipses(cfg_action, patt);
        //else if (action == "learn_shapes") learn_shapes(cfg_action, patt);
        //else if (action == "set_s_thresholds") set_s_thresholds(cfg_action, patt);
        //??else if (action == "test") test(cfg_action, patt);
        else if (action == "add_similarity_edges") do_add_similarity_edges(cfg_action, patt);

    } catch (const libhop_exception& e) {
		cout << e.what() << endl;
		ok = false;
    } catch (const exception& e) {
		cout << "General exception '" << typeid(e).name() << "' with message: '" << e.what() << "'" << endl;
		ok = false;
	}
	return ok;
}
