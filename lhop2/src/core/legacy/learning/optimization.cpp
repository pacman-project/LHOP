// optimization classes - definitions

#include "optimization.h"


// global variables
///////////////////////////////////////////////////////////////////////////////



// optimization_data
///////////////////////////////////////////////////////////////////////////////

optimization_data::optimization_data(const ConfigDictionary& cfg, part_lib* lib) :
    optimizers(MAX_LAYER_NUMBER),
    creators(MAX_LAYER_NUMBER),
    mlearners(MAX_LAYER_NUMBER),
    plearners(MAX_LAYER_NUMBER),
	library(lib)
{
    initialize(cfg);
}

optimization_data::~optimization_data() 
{
	for (auto iter = optimizers.begin(); iter != optimizers.end(); ++iter)
		if (*iter != nullptr) delete *iter;
	for (auto iter = creators.begin(); iter != creators.end(); ++iter)
		if (*iter != nullptr) delete *iter;
	for (auto iter = mlearners.begin(); iter != mlearners.end(); ++iter)
		if (*iter != nullptr) delete *iter;
	for (auto iter = plearners.begin(); iter != plearners.end(); ++iter)
		if (*iter != nullptr) delete *iter;   
}
void optimization_data::add_to_test_set(streamed_pointer res)
{
    test_set.push_back(res);
}

void optimization_data::initialize(const ConfigDictionary& cfg)
{
    char key[200];
    int min_layer, max_layer;

    cfg.getValue(min_layer, "min_layer", true);
    cfg.getValue(max_layer, "max_layer", true);

    // initialize creators & learners
    for (int l = min_layer; l <= max_layer; ++l) {

		string ly_namespace = string("ly") + l;

		ConfigDictionary ly_cfg;
		ly_cfg.fromNamespacePriority(cfg, 1, ly_namespace.c_str());

        if (ly_cfg.doesNamespaceExists("infer")) {
            ConfigDictionary ccfg;
			ccfg.fromNamespacePriority(ly_cfg, 1, "infer");

			cout << "Reading creator configuration in namespace " << ly_namespace  << ".infer" << endl;
            creators[l] = new layern_creator(ccfg);            
		} else {
			throw new_libhop_exception("Warning: No '" + string(ly_namespace) + ".infer' namespace found. Please define namespace or adjust min_layer / max_layer.");
		}
        
        if (ly_cfg.doesNamespaceExists("learn")) {
            ConfigDictionary lcfg;
			lcfg.fromNamespacePriority(ly_cfg, 1, "learn");
			
			cout << "Reading learner configuration in namespace " << ly_namespace  << ".learn" << endl;
            mlearners[l] = new map_learning(lcfg);
            plearners[l] = new part_learning(lcfg);
		} else {
			throw new_libhop_exception("Warning: No '" + string(ly_namespace) + ".learn' namespace found. Please define namespace or adjust min_layer / max_layer.");
		}

        if (ly_cfg.doesNamespaceExists("optimize")) {
            ConfigDictionary ocfg;
			ocfg.fromNamespacePriority(ly_cfg, 1, "optimize");

			cout << "Reading optimizer configuration in namespace " << ly_namespace  << ".optimize" << endl;
            optimizers[l] = new layer1_optimization_new(ocfg, l-1);
            
		} else {
			throw new_libhop_exception("Warning: No '" + string(ly_namespace) + ".optimize' namespace found. Please define namespace (at least optimize.optimize = false should be set) or adjust min_layer / max_layer.");
		}

    }

	// read library of not yet provided by user
	if (library == nullptr) {
		string libname;
		cfg.getValue(libname, "library", true);

		library = part_lib::read(libname);

		// throw exception if we still do not have library
		if (library == nullptr)
			throw custom_libhop_exception(ConfigException, string("Library '" + libname + "' not found"));
	}
}



// cr_layer_optimization
///////////////////////////////////////////////////////////////////////////////

cr_layer_optimization::cr_layer_optimization(optimization_data* oopt, 
        const ConfigDictionary& cfg, int layer) :
    cr_optimization_base(),
    current_train_set()
{
    optimizer = oopt->get_optimizer(layer);
    creator = oopt->get_creator(layer);
    mlearner = oopt->get_map_learner(layer);
    plearner = oopt->get_part_learner(layer);
    current_library = nullptr;
    current_loop = 0;

    steps = cfg.getValueInt("steps", 100);
    init_types_threshold = cfg.getValueDouble("init_types_threshold", 0.0);
    loops = cfg.getValueInt("loops", 1);
    uncovered_radius = cfg.getValueInt("uncovered_radius", 4);
	em_steps = cfg.getValueInt("em_steps", 0);

    cfg.getValue(part_max_number, "part_max_number", 20);
    cfg.getValue(sorting_type, "sorting_type", 0);
    cluster_size = cfg.getValueInt("cluster_size", 5);
    merge_distance_threshold = cfg.getValueDouble("merge_distance_threshold", 2),
    merge_sc_threshold = cfg.getValueDouble("merge_sc_threshold", 1);

    maxima_images = cfg.getValueString("maxima_images", "");
    map_images = cfg.getValueString("map_images", "");
    library_image = cfg.getValueString("library_image", "");
    library_image_sc = cfg.getValueString("library_image_sc", "");
    lib_export_name = cfg.getValueString("lib_export_name", "");
    show_labels = cfg.getValueBool("show_labels", false);
    display_stat = cfg.getValueBool("display_stat", false);
    video = cfg.getValueBool("video", false);
}

void cr_layer_optimization::reset()
{
    current_train_set.clear();
    current_library = nullptr;
}

void cr_layer_optimization::set_library(part_lib* plibrary, bool keep_parts /* = true*/)
{
    current_library = plibrary;
    if (!keep_parts)
        current_library->delete_parts_geq(get_start_layer());
}

// Set HAS_NEXT_LAYER attributes to nodes on get_start_layer - 1 which are not
// covered with the highest layer (max_layer_index()).
// Node n is "covered" when its reconstruction is a subset of the reconstruction 
// all nodes on the highest layer.
void cr_layer_optimization::prepare_result(layer1_result* res, int start_layer)
{
    int topl = res->max_layer_index();
    int bottoml = start_layer - 1;
    set<node*> coverednt;
    vector<node*> tnodes;
    vector<node*> bnodes;
    set<ipoint2> covered;
    int covcount = 0, uncovcount = 0;

    res->get_nodes(tnodes, topl, set<int>());
    res->get_nodes(bnodes, bottoml, set<int>());
    res->add_reconstruction_edges_fwd(topl);
    res->add_reconstruction_edges_fwd(bottoml);

    res->get_neighbors(coverednt, tnodes.begin(), tnodes.end(), EdgeConnection::TO_LAYER0);
    node_set_to_point_set(covered, coverednt.begin(), coverednt.end());

    for (vector<node*>::iterator iter = bnodes.begin(); iter != bnodes.end(); ++iter) {
        node* n = *iter; 
        set<node*> nbset;
        set<ipoint2> nbpts;
        set<ipoint2> diff;

        n->get_neighbor_set(EdgeConnection::TO_LAYER0, nbset);
        node_set_to_point_set(nbpts, nbset.begin(), nbset.end());
        set_difference(diff, nbpts, covered);
        if (diff.empty()) {
            n->set_attr(HAS_NEXT_LAYER); // nbpts is subset of covered, do not update...
            ++covcount;
        } else {
            n->clear_attr(HAS_NEXT_LAYER);
            ++uncovcount;
        }
    }
    cout << '(' << covcount << '/' << uncovcount << ')';
    res->delete_layers_geq(bottoml + 1, false);
}

void cr_layer_optimization::add_to_test_set(streamed_pointer res) 
{
    if (res.is_null()) 
        return;
    
    if (current_library->max_layer_index() < get_start_layer()) {
        // Library was not learned up to the requested layer
        // and we assume that res does not contain get_start_layer
        current_train_set.push_back(res);
    } else {
        layer1_result* resp = (layer1_result*)res.get();

        if (resp->max_layer_index() >= get_start_layer()) {
            // We already have result which is inferred up to some higher layer 
            // Set HAS_NEXT_LAYER attributes to nodes on get_start_layer - 1 which are not
            // covered with the highest layer (max_layer_index()).
            // Node n is "covered" when its reconstruction is a subset of the reconstruction 
            // all nodes on the highest layer.
            prepare_result(resp, get_start_layer());
            res.set(resp);
            current_train_set.push_back(res);
        } else {
            current_train_set.push_back(res);
        }
        delete resp;

    }
}

double cr_layer_optimization::execute()
{
    double result = 1.0;

    PRINT_INFO("============ CREATION OF LAYER " << optimizer->get_layer() << " STARTED =============");
    for (current_loop = 0; current_loop < loops; ++current_loop, optimizer->inc_loop()) {
        int parts_to_keep = current_library->layer_size(get_start_layer());

        PRINT_INFO("__ Part learning __");
        part_learning();
        if (current_loop == 0)
            mlearner->dispose();
        PRINT_INFO("__ Optimization __");
        optimization(parts_to_keep);
        EM_step(parts_to_keep);
        if (!lib_export_name.empty()) {
            stringstream sstr;

            sstr << lib_export_name << (get_start_layer() + 1) << "-" << current_loop << ".plb";
            current_library->save(sstr.str());
        }
    }
    plearner->dispose();

    PRINT_INFO("============ OPTIMIZATION OF LAYER " << optimizer->get_layer() << " ENDED ==============");
    return result;
}

////////////////////////////////////////////////
// main functionality of part_learning mapreduce implementation
streamable* part_learning_mapreduce::map(streamable* item) {
	// we know item should be streamed_pointer*
	layer1_result* res = (layer1_result*)((streamed_pointer*)item)->get();
		
	part_learning* plearner_copy = (part_learning*)plearner->get_copy();

	scmap_t scmap;
	if (current_loop > 0) {
		res->delete_layers_geq(layer); // needed?
		creator->add_layer(res, layer + 1, 0);
		cr_layer_optimization::prepare_result(res, layer);
	}
    get_sc_map(scmap, res, *bin, true);
    plearner_copy->update(res, scmap);
        
	delete res;

	return plearner_copy;
}
streamable* part_learning_mapreduce::reduce(list<streamable*> &item_list){
	part_learning* plearner_copy = (part_learning*)plearner->get_copy();

	// we know each item in list should be part_learning*
	for (auto iter = item_list.begin(); iter != item_list.end(); ++iter) {
		plearner_copy->merge((part_learning*)*iter);
	}

	return plearner_copy;
}

////////////////////////////////////////////////
// main functionality of map_learning mapreduce implementation
streamable* map_learning_mapreduce::map(streamable* item) {
	// we know item should be streamed_pointer*
	layer1_result* res = (layer1_result*)((streamed_pointer*)item)->get();
		
	map_learning* mlearner_copy = (map_learning*)mlearner->get_copy();

	// prepare layer1_result for map_learning (Sets ALLOW_UPDATE_ATTR to all nodes.)
	mlearner_copy->prepare_for_update(res);
	// update statistics within map_learning with layer1_result
    mlearner_copy->update(res);

	delete res;

	return mlearner_copy;
}
streamable* map_learning_mapreduce::reduce(list<streamable*> &item_list) {
	map_learning* mlearner_copy = (map_learning*)mlearner->get_copy();

	// we know each item in list should be map_learning*
	for (auto iter = item_list.begin(); iter != item_list.end(); ++iter) {
		mlearner_copy->merge((map_learning*)*iter);
	}

	return mlearner_copy;
}

////////////////////////////////////////////////
// main functionality of layer_creator_mapreduce mapreduce implementation
streamable* layern_creator_mapreduce::map(streamable* item) {
	// we know item should be streamed_pointer*
	layer1_result* res = (layer1_result*)((streamed_pointer*)item)->get();
	
	// delete layers if requested
	if (delete_layers_geq)
		res->delete_layers_geq(layer-1);	

	// add layer
	creator->add_layer(res, layer, 0);

	// add links for reconstruction edges if requested 
	if (add_reconstruction_edges_link)
		res->add_reconstruction_edges_link(layer-1);	

	return res;
}

streamable* layern_creator_mapreduce::reduce(list<streamable*> &item_list) {
	return nullptr;
}

string layern_creator_mapreduce::map_get_key(streamable* item) { 
	// each item should be saved into its own group (must have distinct key)
	return ((streamed_pointer*)item)->get_name_only();
}

#include <time.h>
typedef part_learning part_learning_t;
void cr_layer_optimization::part_learning()
{
	clock_t start_t, end_t;
    int layer = get_start_layer();

	K_bin* bin = new K_bin(12, 2, 4, 7);
	start_t = clock();
    // Update neighborhoods only during the first loop
    if (current_loop == 0) {
        PRINT_INFO("Map update started");

        int count = 0;

        // Reset the learner
        mlearner->reset(current_library);
		
		if (mapreduce_deployer == nullptr) {
			// the same functionality as in true statements except it is implemneted using map/reduce system
			for (list<streamed_pointer>::iterator iter = current_train_set.begin(); iter != current_train_set.end(); ++iter) {
				layer1_result* res = (layer1_result*)iter->get();

				mlearner->prepare_for_update(res);
				mlearner->update(res);
				if (video) {
					stringstream sstr;

					end_dir(map_images);
					sstr << map_images << "map_ly-" << setfill('0') << setw(2) << layer;
					sstr << "_img-" << setfill('0') << setw(3) << count;
					sstr << "_map_%03d-%03d.png";
					mlearner->display_statistics(sstr.str());
				}
				++count;
				delete res;
			}			

		} else {
			/////////////////////////////////////////////////////////////////////////
			// CALL TO map_learning_mapreduce::map and map_learning_mapreduce::reduce
			// (check base_deployer_mapreudce for default calling implementation)

			// the same functionality as in true statements except it is implemneted using map/reduce system
			mapreduce_items* mr_process_list = new mapreduce_items();

			for (list<streamed_pointer>::iterator iter = current_train_set.begin(); iter != current_train_set.end(); ++iter) {
				mr_process_list->items.push_back(&*iter); // be careful with reference to item within current_train_set
			}

			map_learning_mapreduce* mlearning_mapreduce = new map_learning_mapreduce();

 			mlearning_mapreduce->set_map_learning(mlearner);
			
			mapreduce_result* result = mapreduce_deployer->submit(mr_process_list, mlearning_mapreduce);

			list<streamable*> stremable_results = result->get_result();

			if (stremable_results.size() == 1) {
				// results OK	
				streamable* mlearning_result = stremable_results.front();
				if (mlearning_result != mlearner) {
					// copy values to existing plearner
					*mlearner = *(map_learning*)mlearning_result;
					// delete mlearning_result
					delete mlearning_result;
				}
			} else {
				// results NOT OK				
				throw new_libhop_exception("cr_layer_optimization::part_learning() - mapreduce_deployer returned invalid results for map_learning_mapreduce");
			}

			delete result;
			delete mlearning_mapreduce;
		}
		
        if (!map_images.empty() && !video) {
            PRINT_INFO("Saving maps...");
            mlearner->display_statistics(map_images);
        }

        // Find maxima
        PRINT_INFO("Finding maxima...");
        plearner->reset(*mlearner);

        // Save maxima
        if (!maxima_images.empty()) {
            PRINT_INFO("Saving maxima...");
            plearner->display_maxima(maxima_images);
        }
    }	
    plearner->reset();

	end_t = clock();
	printf("TIME: map learning - %f\n",((float)(end_t-start_t))/CLOCKS_PER_SEC);

    // Learn parts
    PRINT_INFO("Learning parts...");
    if (current_loop > 0)
        creator->set_library(current_library);

	start_t = clock();
	if (mapreduce_deployer == nullptr) {
		for (list<streamed_pointer>::iterator iter = current_train_set.begin(); iter != current_train_set.end(); ++iter) {
			layer1_result* res = (layer1_result*)iter->get();
			scmap_t scmap;

			// when using multiple loops and have already done one pass (i.e. already have created parts for the same layer)
			// then we need to ignore statistics of parts/subparts coming from parts already covered by existing parts in library
			if (current_loop > 0) {
				// delete current layer just in case (since we might have created new library in previous loop the inferred parts could be different now)
				res->delete_layers_geq(layer); // needed?
				// infer parts with the current layer
				creator->add_layer(res, layer + 1, 0);
				// now flag all parts that will have to be ignored i.e. parts that are covered by the parts of already existing layer
				prepare_result(res, layer);
			}
			get_sc_map(scmap, res, *bin, true);
			plearner->update(res, scmap);
			delete res;
		}
		
	} else {
		/////////////////////////////////////////////////////////////////////////
		// CALL TO part_learning_mapreduce::map and part_learning_mapreduce::reduce 
		// (check base_deployer_mapreudce for default calling implementation)

		// the same functionality as in true statements except it is implemented using map/reduce system
		mapreduce_items* mr_process_list = new mapreduce_items();

		for (list<streamed_pointer>::iterator iter = current_train_set.begin(); iter != current_train_set.end(); ++iter) {
			mr_process_list->items.push_back(&*iter); // be careful with reference to item within current_train_set
		}

		part_learning_mapreduce* plearning_mapreduce = new part_learning_mapreduce();

		plearning_mapreduce->set_part_learning(plearner);
		plearning_mapreduce->set_creator(creator);
		plearning_mapreduce->set_bin(bin);
		plearning_mapreduce->set_current_loop(current_loop);
		plearning_mapreduce->set_layer(layer);

		mapreduce_result* result = mapreduce_deployer->submit(mr_process_list, plearning_mapreduce);

		list<streamable*> stremable_results = result->get_result();

		if (stremable_results.size() == 1) {
			// results OK	
			streamable* plearning_result = stremable_results.front();
			if (plearning_result != plearner) {
				// copy values to existing plearner
				*plearner = *(part_learning_t*)plearning_result;
				// delete plearning_result
				delete plearning_result;
			}
		} else {
			// results NOT OK				
			throw new_libhop_exception("cr_layer_optimization::part_learning() - mapreduce_deployer returned invalid results for part_learning_mapreduce");
		}		

		delete result;
		delete plearning_mapreduce;
	}
	end_t = clock();
	printf("TIME: part learning - %f\n",((float)(end_t-start_t))/CLOCKS_PER_SEC);
	
    PRINT_INFO("Total of " << plearner->get_stat_size() << " sequences found.");

    if (display_stat) plearner->save_stat();

	start_t = clock();
    plearner->add_to_library(current_library, part_max_number[current_loop], sorting_type[current_loop], cluster_size, 
        merge_distance_threshold, merge_sc_threshold);

	end_t = clock();
	printf("TIME: add to lib - %f\n",((float)(end_t-start_t))/CLOCKS_PER_SEC);
    PRINT_INFO("" << current_library->layer_size(layer) << " parts added to library (limit: " 
        << part_max_number[current_loop] << ")");

    if (!library_image.empty()) {
        current_library->save_all(library_image.c_str(), plearner->get_source_layer() + 2, 0, -1, show_labels);
    }
    if (!library_image_sc.empty()) {
        current_library->save_all_sc(library_image_sc.c_str(), plearner->get_source_layer() + 2, 0, -1, show_labels);
    }

    PRINT_INFO("Part learning ended...");

	delete bin;
}
void cr_layer_optimization::optimization(int parts_to_keep)
{
    if (!optimizer->on())
        return;

    // Change some settings of creator
    double mem_rrt = creator->realization_ratio_threshold;

    creator->realization_ratio_threshold = 0.9;

    int layer = get_start_layer();
    list<streamed_pointer> opt_train_set;

    // Apply new parts to images
    PRINT_INFO("Inference...");

	clock_t start_t = clock();

    creator->set_library(current_library);
	if (mapreduce_deployer == nullptr) {
		for (list<streamed_pointer>::iterator iter = current_train_set.begin(); iter != current_train_set.end(); ++iter) {
			layer1_result* res = (layer1_result*)iter->get();

			res->delete_layers_geq(layer);
			creator->add_layer(res, layer + 1, 0);
			res->add_reconstruction_edges_link(layer);
			opt_train_set.push_back(streamed_pointer(res));
			delete res;
		}
	} else {
		/////////////////////////////////////////////////////////////////////////
		// CALL TO layern_creator_mapreduce::map and layern_creator_mapreduce::reduce
		// (check base_deployer_mapreudce for default calling implementation)

		mapreduce_items* mr_process_list = new mapreduce_items();

		for (list<streamed_pointer>::iterator iter = current_train_set.begin(); iter != current_train_set.end(); ++iter) {
			mr_process_list->items.push_back(&*iter); // be careful with reference to item within current_train_set
		}

		layern_creator_mapreduce* creator_mapreduce = new layern_creator_mapreduce();

 		creator_mapreduce->set_creator(creator);
		creator_mapreduce->set_layer(layer + 1);
		
		creator_mapreduce->set_delete_layers_geq(true);
		creator_mapreduce->set_add_reconstruction_edges_link(true);

		mapreduce_result* result = mapreduce_deployer->submit(mr_process_list, creator_mapreduce);

		list<streamable*> stremable_results = result->get_result();

		if (stremable_results.size() == current_train_set.size()) {
			// results OK	
			for (auto iter = stremable_results.begin(); iter != stremable_results.end(); ++iter) {
				opt_train_set.push_back(*(static_cast<streamed_pointer*>(*iter)));
				delete *iter;
			}			
		} else {
			// results NOT OK				
			throw new_libhop_exception("cr_layer_optimization::get_best_test_set() - mapreduce_deployer returned invalid results for layern_creator_mapreduce (not all results from training set were returned)");
		}

		delete result;
		delete creator_mapreduce;
	}
	clock_t end_t = clock();
	printf("TIME: part optimization (inference) - %f\n",((float)(end_t-start_t))/CLOCKS_PER_SEC);

	start_t = clock();
    // Restore parameters of inference
    creator->realization_ratio_threshold = mem_rrt;

    // Update optimizer
    PRINT_INFO("Running optimization...");
    optimizer->set_steps(steps);
    optimizer->set_library(current_library, parts_to_keep);
    optimizer->add_to_test_set(opt_train_set.begin(), opt_train_set.end());
    optimizer->make_steps();
	// make sure to also keep similar parts (TO_LYR_SIMILAR) and to delete parts that
	// are root for similarity (i.e. TO_LYR_SIMROOT)
    optimizer->keep_final_parts(); 
    optimizer->print_final_parts();
    optimizer->display_final_parts();

	end_t = clock();
	printf("TIME: part optimization (finilizing) - %f\n",((float)(end_t-start_t))/CLOCKS_PER_SEC);
}

void cr_layer_optimization::get_best_test_set(list<streamed_pointer>& container)
{ 
	clock_t start_t = clock();
    container.clear();
    creator->set_library(current_library);

	if (mapreduce_deployer == nullptr) {
		for (list<streamed_pointer>::iterator iter = current_train_set.begin(); iter != current_train_set.end(); ++iter) {
			layer1_result* res = (layer1_result*)iter->get();
			creator->add_layer(res, get_start_layer() + 1, 0);
			container.push_back(streamed_pointer(res));
			delete res;
		}
	} else {
		/////////////////////////////////////////////////////////////////////////
		// CALL TO layern_creator_mapreduce::map and layern_creator_mapreduce::reduce
		// (check base_deployer_mapreudce for default calling implementation)

		mapreduce_items* mr_process_list = new mapreduce_items();

		for (list<streamed_pointer>::iterator iter = current_train_set.begin(); iter != current_train_set.end(); ++iter) {
			mr_process_list->items.push_back(&*iter); // be careful with reference to item within current_train_set
		}

		layern_creator_mapreduce* creator_mapreduce = new layern_creator_mapreduce();

 		creator_mapreduce->set_creator(creator);
		creator_mapreduce->set_layer(get_start_layer() + 1);

		mapreduce_result* result = mapreduce_deployer->submit(mr_process_list, creator_mapreduce);

		list<streamable*> stremable_results = result->get_result();

		if (stremable_results.size() == current_train_set.size()) {
			// results OK	
			for (auto iter = stremable_results.begin(); iter != stremable_results.end(); ++iter) {
				container.push_back(*(static_cast<streamed_pointer*>(*iter)));
				delete *iter;
			}			
		} else {
			// results NOT OK				
			throw new_libhop_exception("cr_layer_optimization::get_best_test_set() - mapreduce_deployer returned invalid results for layern_creator_mapreduce (not all results from training set were returned)");
		}

		delete result;
		delete creator_mapreduce;
	}
	clock_t end_t = clock();
	printf("TIME: get_best_test_set (inference) - %f\n",((float)(end_t-start_t))/CLOCKS_PER_SEC);
}

void cr_layer_optimization::EM_step(int parts_to_keep)
{
    if (em_steps <= 0) 
        return;

    typedef map<int, map<int, pair<distribution2, histogram> > > result_t;

    // Change some settings of creator
    double mem_rrt = creator->realization_ratio_threshold;
    bool mem_aen = creator->add_edge_names;
    double mem_contr = creator->layer_contraction;
    bool mem_igr = creator->identity_g_response;

    creator->realization_ratio_threshold = 0.9;
    creator->add_edge_names = true;
    creator->layer_contraction = 1.0;
    creator->identity_g_response = true;

    creator->set_library(current_library);

    for (int ems = 0; ems < em_steps; ++ems) {
        result_t result;
        int layer = get_start_layer();

        // Apply new parts to images
        PRINT_INFO("Running EM step " << ems + 1 << " of " << em_steps << "...");

        for (list<streamed_pointer>::iterator iter = current_train_set.begin(); iter != current_train_set.end(); ++iter) {
            layer1_result* res = (layer1_result*)iter->get();

            res->delete_layers_geq(layer);
            creator->add_layer(res, layer + 1, 0);
            res->add_reconstruction_edges(layer);

            update_mean_map(result, res, layer, current_library, 0.0, 0.0);
            delete res;
        }

        // Update library
        result_t fresult;

        for (result_t::iterator iter = result.begin(); iter != result.end(); ++iter) {
            if (iter->first >= parts_to_keep) 
                fresult.insert(*iter);
        }
        update_means(fresult, layer, current_library);
    }

    // Restore parameters of inference
    creator->realization_ratio_threshold = mem_rrt;
    creator->add_edge_names = mem_aen;
    creator->layer_contraction = mem_contr;
    creator->identity_g_response = mem_igr;

}

// cr_set_optimization
///////////////////////////////////////////////////////////////////////////////

cr_set_optimization::cr_set_optimization(optimization_data* oopt, int nsteps, const vector<cr_optimization_base*>& oset) :
    cr_optimization_base(),
    init_test_set(),
    current_train_set(),
    best_test_set(),
    optimizers(oset.begin(), oset.end())
{
    init_library = current_library = best_library = nullptr;
    steps = nsteps;
}

void cr_set_optimization::reset()
{
    init_library = current_library = best_library = nullptr;
    init_test_set.clear();
    current_train_set.clear();
    best_test_set.clear();
}

void cr_set_optimization::set_library(part_lib* plibrary, bool keep_parts /* = true */)
{
    init_library = plibrary;
    if (!keep_parts)
        init_library->delete_parts_geq(get_start_layer());
}

void cr_set_optimization::add_to_test_set(streamed_pointer resp) 
{
    if (!resp.is_null()) {
        layer1_result* res = (layer1_result*)resp.get();

        res->delete_layers_geq(get_start_layer());
        init_test_set.push_back(streamed_pointer(res));
        delete res;
    }
}

void cr_set_optimization::add_to_init_test_set(optimization_data* oopt)
{
    for (auto tsiter = oopt->get_test_set_begin(); tsiter != oopt->get_test_set_end(); ++tsiter) {
        if (tsiter->is_null())
            continue;

        layer1_result* resp = (layer1_result*)tsiter->get();

        if (resp->max_layer_index() >= get_start_layer())
            resp->delete_layers_geq(get_start_layer());
        else {
            cout << "Adding layer(s) to train example (from " << resp->max_layer_index() + 1 << 
                " to " << get_start_layer() - 1 << ").";
            for (int l = resp->max_layer_index() + 1; l < get_start_layer(); ++l) {
                layern_creator* creator = oopt->get_creator(l + 1);

                if (creator == nullptr) {
                    cout << "Creator for layer " << l << " doest not exist!. Set \"min_layer\" to correct value." << endl;
                    throw;
                }
                if (creator->get_library() == nullptr) 
                    creator->set_library(init_library);
                if (creator->get_library() == nullptr) {
                    cout << "Library is not set in cr_set_optimization." << endl;
                    throw;
                }

                oopt->get_creator(l + 1)->add_layer(resp, l + 1, 0);
            }
            cout << endl;
        }
        init_test_set.push_back(streamed_pointer(resp));
        delete resp;
    }
}


double cr_set_optimization::execute()
{
    double best = 0.0;
    double val = 0.0;
    list<streamed_pointer>::iterator iter;
    int scount = steps;

    best_test_set.clear();
    while (scount-- > 0) {
        PRINT_INFO("-----------< SET >-----------< " << scount << " >-----------< SET >-----------");
        current_train_set.clear();
        for (iter = init_test_set.begin(); iter != init_test_set.end(); ++iter) {
            layer1_result* res = (layer1_result*)iter->get();
            current_train_set.push_back(streamed_pointer(res));
        }

        current_library = (part_lib*)init_library->get_copy();

        // Loop through all optimizers
        for (list<cr_optimization_base*>::iterator oiter = optimizers.begin(); oiter != optimizers.end(); ++oiter) {
            cr_optimization_base* optimizer = *oiter;

            optimizer->reset();
            PRINT_INFO("__ SET LIBRARY");
            optimizer->set_library(current_library);
            PRINT_INFO("__ ADD RANGE");
            optimizer->add_range_to_test_set(current_train_set.begin(), current_train_set.end());
            PRINT_INFO("__ EXECUTE!");
            val = optimizer->execute();
            PRINT_INFO("__ GET BEST LIB");
            current_library = optimizer->get_best_library();
            optimizer->get_best_test_set(current_train_set);

			// reset in order to clear disk images of streamed_pointer saved in test set
			optimizer->reset();
        }

        if (val <= best) 
            delete current_library;
        else {
            PRINT_INFO("--- BEST FOUND");
            best_library = current_library;
            best_test_set.assign(current_train_set.begin(), current_train_set.end());
            best = val;
        }
    }
    return best;
}

// public functions definitions
///////////////////////////////////////////////////////////////////////////////
