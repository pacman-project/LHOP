
#include <opencv2/opencv.hpp>
#include <queue>

#include "optimization_learning.h"

#include "utils/graphs/graph_utils.h"


// layer1_optimization_new
///////////////////////////////////////////////////////////////////////////////

layer1_optimization_new::layer1_optimization_new(const ConfigDictionary& cfg, int player)
{
    layer = player;
    loop = 0;
    init(cfg);
}
 
void layer1_optimization_new::init(const ConfigDictionary& cfg)
{
    cfg.getValue(cover_thresh, "covered_threshold", 0.85);
    cfg.getValue(int_thresh, "intersection_threshold", 0.2);
    cfg.getValue(bite_size, "bite_size", 10);
    library_image_file = cfg.getValueString("library_image_file", "");
    show_labels = cfg.getValueBool("show_labels", true);
    if (cfg.isDefined("optimize")) optimization.set_val(cfg.getValueBool("optimize", true) ? 1 : 0);
    else cfg.getValue(optimization, "optimization", 1);
}

void layer1_optimization_new::set_library(part_lib* plibrary, int keep)
{
    library = plibrary;
    parts_to_keep = keep;
}

void layer1_optimization_new::inc_loop()
{
    loop++;
    cover_thresh.inc();
    int_thresh.inc();
    bite_size.inc();
    optimization.inc();
}

void layer1_optimization_new::add_to_test_set(streamed_pointer tres)
{
    workingset.push_back(tres);
}

void layer1_optimization_new::reset()
{
    workingset.clear();
}

double layer1_optimization_new::make_steps()
{
    if (optimization == 0) {
        for (int i = 0; i < library->layer_size(layer); ++i) 
            final_parts.insert(i);
        return 1.0;
    }

    set<int> init_parts;
    int maxi = min(library->layer_size(layer), parts_to_keep);

    for (int i = 0; i < maxi; ++i)
        init_parts.insert(i);

    cout << "  (Keeping " << init_parts.size() << " parts)" << endl;

    switch (optimization) {
        case 2: 
            final_parts = optimize_layer_2(library, workingset, layer, init_parts, cover_thresh, int_thresh, bite_size);
            break;
        default:
            final_parts = optimize_layer(library, workingset, layer, init_parts, cover_thresh, int_thresh, bite_size);
    }
    return 1.0;
}

void layer1_optimization_new::keep_final_parts()
{
    vector<int> fpv(final_parts.begin(), final_parts.end());

    library->keep_clusters(layer, final_parts);
}

void layer1_optimization_new::print_final_parts()
{
    cout << "Best parts:";
    for (set<int>::iterator iter = final_parts.begin(); iter != final_parts.end(); ++iter) {
        cout << ' ' << *iter;
    }
    cout << endl;
}

void layer1_optimization_new::display_final_parts(const string& file, bool show_labels)
{
    if (!file.empty()) {
        library->save_all(file.c_str(), layer + 1, 0, -1, show_labels);
        library->drop_images();    
    }
}

void layer1_optimization_new::display_final_parts()
{
    display_final_parts(library_image_file, show_labels);
}


struct nodes_at_sort_f {
    int nbname;

    nodes_at_sort_f() : nbname(EdgeConnection::TO_LAYER0) { }
    bool operator()(node* n, node* m) const
    {
        layer1_data* nd = (layer1_data*)n->data;
        layer1_data* md = (layer1_data*)m->data;

        return (0.9*nd->r(G_RESPONSE) + 0.1*nd->r(R_RESPONSE)) > (0.9*md->r(G_RESPONSE) + 0.1*md->r(R_RESPONSE));
    }
};



class optimize_layer_stat : public streamable {
public:
	vector<int> statistics;
	optimize_layer_stat(int num_parts) : statistics(vector<int>(num_parts)) {}

	// streamable implementations
	virtual streamable* make_instance() const { return new optimize_layer_stat(0); }
	virtual void read_from_stream(istreamer& is) { is.read(statistics); }
	virtual void write_to_stream(ostreamer& os) { os.write(statistics);	}
};

class optimize_layer_mapreduce : public base_mapreduce {
	set<int> current_set; vector<set<ipoint2>> ly0sets;
	int num_parts; int layer; double cover_thresh; double int_thresh;
public:
	optimize_layer_mapreduce() { }
	optimize_layer_mapreduce(const set<int>& current_set, const vector<set<ipoint2>>& ly0sets, int num_parts, int layer, double cover_thresh, double int_thresh) :
		current_set(current_set), ly0sets(ly0sets), num_parts(num_parts), layer(layer), cover_thresh(cover_thresh), int_thresh(int_thresh) {}
	virtual ~optimize_layer_mapreduce() { }

	////////////////////////////////////////////////
	// main functionality is written in map and reduce
	// this is part of optimize_layer() that gathres statistics for each layer1_result
	virtual streamable* map(streamable* item) {
		// we know item should be streamed_pointer*
		layer1_result* res = (layer1_result*)((streamed_pointer*)item)->get();
		
		optimize_layer_stat* result = new optimize_layer_stat(num_parts);		
		
		int tolayer0 = EdgeConnection::TO_LAYER0;

		vector<node*>& s_nodes = res->shape_nodes[layer];

		// Calc covering ratio
		vector<node*> visited;
		set<node*> coveredn;
		set<node*> maxcoveredn;
		set<ipoint2> covered;
		set<ipoint2> maxcovered;

		select_nodes_by_type(visited, res, s_nodes, current_set);
		res->recurse(visited, tolayer0, coveredn);
		node_set_to_region_set(covered, coveredn, ly0sets, 0);
		res->recurse(s_nodes, tolayer0, maxcoveredn);
		node_set_to_region_set(maxcovered, maxcoveredn, ly0sets, 0);

		double cover_ratio = (double)covered.size()/(double)maxcovered.size();

		if (cover_ratio >= cover_thresh) return result;

		// Update statistics
		for (vector<node*>::iterator sniter = s_nodes.begin(); sniter != s_nodes.end(); ++sniter) {
			node* n = *sniter;
			vector<node*> nvec;
			bool found = false;

			res->nodes_at(nvec, n);
			for (vector<node*>::iterator nviter = nvec.begin(); nviter != nvec.end(); ++nviter) {
				node* nn = *nviter;
				layer1_data* nnd = (layer1_data*)nn->data;
				if (current_set.find(nnd->m) == current_set.end()) {
					set<ipoint2> prset; 
					set<node*> nrset;

					res->recurse_from_node(nn, tolayer0, nrset);
					node_set_to_region_set(prset, nrset, ly0sets, 0);
					int is = intersection_size(prset, covered);

					if (is < prset.size()*int_thresh)
						++result->statistics[nnd->m];
				}
			}

		}
		
		delete res;

		return result;
	}
	virtual streamable* reduce(list<streamable*> &item_list) {

		optimize_layer_stat* result = new optimize_layer_stat(num_parts);		

		for (auto iter = item_list.begin(); iter != item_list.end(); ++iter) {
			optimize_layer_stat* partial_result = (optimize_layer_stat*)*iter;
			for (int i = 0; i < num_parts; i++)
				result->statistics[i] += partial_result->statistics[i];

		}

		return result;
	}

	virtual string map_get_key(streamable* item) { return "merged_result"; }

	virtual bool has_reduce() { return true; }
	virtual bool should_wrap_results() { return false; }

	// streamable implementations
	virtual streamable* make_instance() const { return new optimize_layer_mapreduce(); }
	virtual void read_from_stream(istreamer& is) {
		is.read(current_set);
		is.read(ly0sets);
		is.read(num_parts);
		is.read(layer);
		is.read(cover_thresh);
		is.read(int_thresh);
	}
	virtual void write_to_stream(ostreamer& os) {
		os.write(current_set);
		os.write(ly0sets);
		os.write(num_parts);
		os.write(layer);
		os.write(cover_thresh);
		os.write(int_thresh);
	}
};

/// Optimizes library on specific layer using workingset; algorithm parameters are 
/// * cover_thresh: threshold for covering ratio ( = |covered|/|maxcovered|)
///        if image is covered >= cover_threshold, skip it
/// * int_thresh: update part statistics if intersection of this part with already 
///        covered part of image >= int_thresh * size of part
/// * bite_size: how many parts are processed in one greedy step
/// * angle_thresh, a_thresh, e_thresh: thresholds for similarity of ellipses
///       (angle, major semi-axis, eccentricity)
/// Other
/// * init_set: initial parts (to be considered fixed)
/// Returns optimized set of parts (includes init_set).
set<int> optimize_layer(part_lib* library, list<streamed_pointer>& workingset, int layer, const set<int>& init_set,
    double cover_thresh, double int_thresh, int bite_size)
{
    typedef set<ipoint2> pset_t;

    if (library->layer_size(layer) == 0)
        return set<int>();

    nodes_at_sort_f sortf;
    set<int> current_set(init_set);
    vector<node*>& parts = library->parts[layer];
    vector<node*> bite;
    int tolayer0 = EdgeConnection::TO_LAYER0;
    vector<pset_t> ly0sets;

    library->get_regions(1, ly0sets);

	base_deployer_mapreudce* mapreduce_deployer = new openmp_deployer_mapreudce();

    do {        
		vector<int> statistics;
		int num_parts = library->layer_size(layer);		

		cout << '.';		
		
		//////////////////////////////////////////////////////////////////////////////
		// CALL TO optimize_layer_mapreduce::map and optimize_layer_mapreduce::reduce
		// (check base_deployer_mapreudce for default calling implementation)
		{	
			mapreduce_items* mr_process_list = new mapreduce_items();

			for (auto iter = workingset.begin(); iter != workingset.end(); ++iter) {
				mr_process_list->items.push_back(&*iter); 
			}

			optimize_layer_mapreduce* collect_stat_mapreduce = new optimize_layer_mapreduce(current_set, ly0sets, num_parts, layer, cover_thresh, int_thresh);

			mapreduce_result* result = mapreduce_deployer->submit(mr_process_list, collect_stat_mapreduce);

			list<streamable*> stremable_results = result->get_result();

			if (stremable_results.size() == 1) {
				// results OK	
				optimize_layer_stat* collect_stat_result = (optimize_layer_stat*)stremable_results.front();
			
				statistics = collect_stat_result->statistics;
			
				delete collect_stat_result;
			} else {
				// results NOT OK				
				throw new_libhop_exception("optimize_layer() - mapreduce_deployer returned invalid results for optimize_layer_mapreduce (not all results from training set were returned)");
			}

			delete mr_process_list;
			delete collect_stat_mapreduce;
		}

        // Selection process
        vector<int> statordering = ordering<int>(statistics.begin(), statistics.end(), greater<int>());
        vector<node*> candidateparts;

        cout << '|' << statistics[statordering[0]];
        for (int i = 0; i < (int)statordering.size() && (int)candidateparts.size() < bite_size*8 &&   // !!!!8
                statistics[statordering[i]] > 0; ++i) {
            candidateparts.push_back(parts[statordering[i]]);
        }
        get_dissimilar_cluster(bite, candidateparts, bite_size);
        for (vector<node*>::iterator biter = bite.begin(); biter != bite.end(); ++biter) {
            lib_data* pd = (lib_data*)(*biter)->data;

            current_set.insert(pd->type);

        }

    } while (!bite.empty());

	delete mapreduce_deployer;

    cout << endl;
    for (set<int>::iterator iter = current_set.begin(); iter != current_set.end(); ++iter) 
        cout << *iter << ' ';
    cout << endl;

    return current_set;
}

/// Optimizes library on specific layer using workingset; algorithm parameters are 
/// * cover_thresh: threshold for covering ratio ( = |covered|/|maxcovered|)
///        if image is covered >= cover_threshold, skip it
/// * int_thresh: update part statistics if maximal intersection of this part with some 
///        already chosen part <= int_thresh*(size of this part)
/// * bite_size: how many parts are processed in one greedy step
/// * angle_thresh, a_thresh, e_thresh: thresholds for similarity of ellipses
///       (angle, major semi-axis, eccentricity)
/// Other
/// * init_set: initial parts (to be considered fixed)
/// Returns optimized set of parts (includes init_set).
set<int> optimize_layer_2(part_lib* library, list<streamed_pointer>& workingset, int layer, const set<int>& init_set,
    double cover_thresh, double int_thresh, int bite_size)
{
    cout << "optimize_layer_2" << endl;
    //cout << "Int_thresh = " << int_thresh << endl;

    typedef set<ipoint2> pset_t;

    if (library->layer_size(layer) == 0)
        return set<int>();

    nodes_at_sort_f sortf;
    set<int> current_set(init_set);
    vector<node*>& parts = library->parts[layer];
    vector<node*> bite;
    int tolayer0 = EdgeConnection::TO_LAYER0;
    int tonext = EdgeConnection::TO_NEXT_LAYER;
    vector<pset_t> ly0sets;

    library->get_regions(1, ly0sets);
    do {
        vector<double> statistics(library->layer_size(layer), 0.0);

		cout << '.';
        for (auto iter = workingset.begin(); iter != workingset.end(); ++iter) {
			layer1_result* res = (layer1_result*)iter->get();
            vector<node*>& s_nodes = res->shape_nodes[layer];

            // Calc covering ratio
            vector<node*> visited;
           
            select_nodes_by_type(visited, res, s_nodes, current_set);
           
            // Link visited nodes "backwards"
            for (vector<node*>::iterator viter = visited.begin(); viter != visited.end(); ++viter) {
                node* vn = *viter;
                set<node*> vnsupp;
                
                res->recurse_from_node(vn, tolayer0, vnsupp);
                for (set<node*>::iterator siter = vnsupp.begin(); siter != vnsupp.end(); ++siter) {
                    res->add_edge(*siter, vn, tonext);  
                }
            }

            // Update statistics
            for (vector<node*>::iterator sniter = s_nodes.begin(); sniter != s_nodes.end(); ++sniter) {
                node* n = *sniter;
                vector<node*> nvec;

                res->nodes_at(nvec, n);
                for (vector<node*>::iterator nviter = nvec.begin(); nviter != nvec.end(); ++nviter) {
                    node* nn = *nviter;
                    layer1_data* nnd = (layer1_data*)nn->data;

                    if (current_set.find(nnd->m) == current_set.end()) {
                        set<node*> nnsupp;
                        set<ipoint2> nnsuppr; 
                        set<node*> bnodes;
                        bool update = true;

                        res->recurse_from_node(nn, tolayer0, nnsupp);
                        node_set_to_region_set(nnsuppr, nnsupp, ly0sets, 0);
                        res->get_neighbors(bnodes, nnsupp.begin(), nnsupp.end(), tonext);

                        for (set<node*>::iterator siter = bnodes.begin(); siter != bnodes.end(); ++siter) {
                            node* bn = *siter;

                            if (bn == nullptr || bn == nn) 
                                continue;

                            set<node*> bnsupp;
                            set<ipoint2> bnsuppr;

                            res->recurse_from_node(bn, tolayer0, bnsupp);
                            node_set_to_region_set(bnsuppr, bnsupp, ly0sets, 0);

                            if (intersection_size(bnsuppr, nnsuppr) >= nnsuppr.size()*int_thresh) {
                                update = false;
                                break;
                            }
                        }

                        if (update) {
                            statistics[nnd->m] += nnd->vval();
                        }

                    }
                }

            }

			delete res;
        }
        
        // Selection process
        vector<int> statordering = ordering<double>(statistics.begin(), statistics.end(), greater<double>());
        vector<node*> candidateparts;

        for (int i = 0; i < (int)statordering.size() && (int)candidateparts.size() < bite_size*8 &&   // !!!!8
                statistics[statordering[i]] > 0; ++i) {
            candidateparts.push_back(parts[statordering[i]]);
        }
        get_dissimilar_cluster(bite, candidateparts, bite_size);
        for (vector<node*>::iterator biter = bite.begin(); biter != bite.end(); ++biter) {
            lib_data* pd = (lib_data*)(*biter)->data;

            current_set.insert(pd->type);

        }

    } while (!bite.empty());

    cout << endl;
    for (set<int>::iterator iter = current_set.begin(); iter != current_set.end(); ++iter) 
        cout << *iter << ' ';
    cout << endl;

    return current_set;
}
