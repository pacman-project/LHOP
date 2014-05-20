// inference.cpp : Defines the entry point feature extraction functionality
//

#include "tools/main_toolset.h"

#include <stdio.h>
#include <string>
#include <cctype>

#include "core/legacy/layers.h"

#include "utils/utils.h"
#include "utils/configuration.h"

// Returns pair (#of nodes at layer, #of covered positions, i.e. shape_nodes.size())
iipair layer_node_count(layer1_result* res, int layer)
{
	if (layer < 0 || layer > res->max_layer_index()) 
		return iipair(0, 0);

	vector<node*>& snodes = res->shape_nodes[layer];
	int count = 0;

	for (vector<node*>::iterator iter = snodes.begin(); iter != snodes.end(); ++iter) {
		node* n = *iter;

		while (n != nullptr) {
			++count;
			n = ((layer1_data*)n->data)->next;
		}
	}
	return iipair(count, (int)snodes.size());
}


// Return the number of points (shape positions) where 'res1' is different
// from 'res' on layer 'layer'.
int layer_0_difference(layer1_result* res, layer1_result* res1, double tolerance)
{
	if (!res->grid(0)) res->init_grid(0);

	vector<node*>& s_nodes = res1->shape_nodes[0];
	int result = 0;

	for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
		node* n = *iter;
		layer1_data* nd = (layer1_data*)n->data;
		bool found = false;

		for (int dx = -2; dx < 3 && !found; ++dx) 
			for (int dy = -2; dy < 3 && !found; ++dy) {
				node* n0 = res->node_at_safe(nd->x + dx, nd->y + dy, 0);        

				while (n0 != nullptr) {
					layer1_data* n0d = (layer1_data*)n0->data;

					if (abs(n0d->r(R_RESPONSE) - nd->r(R_RESPONSE)) > tolerance)
						break;
					if (n0d->m == nd->m) {
						found = true;
						break;
					}
					n0 = n0d->next;
				}
			}
		if (!found) ++result;
	}
	return result;
}


void read_layer1_result_vector(vector<layer1_result*>& res, const list<string>& names, const string& srcdir)
{
	for (list<string>::const_iterator fiter = names.begin(); fiter != names.end(); ++fiter) {
		if (fiter == names.begin()) 
			cout << "Processing " << *fiter;
		else
			cout << '.';

		layer1_result* r;

		read_layer1_result(r, srcdir + *fiter); 
		if (r != nullptr) 
			res.push_back(r);
		else {
			cout << "error reading file." << endl;
			break;
		}
	}
	cout << ' ';
}

string LegacyInferenceToolset::getShortDescription() {
	return "Legacy version of inference - using only layern_creators code";
}
string LegacyInferenceToolset::getLongDescription() {
	return "";
}
string LegacyInferenceToolset::getUsageDescription() {
	return "";
}

bool LegacyInferenceToolset::areArgumentsValid(int argc, char* argv[]) {
	// TODO: do more inteligent verification
	return argc > 2  ? true : false;
}

/// Inference of *multiple layers* on files described by pattern -- using config cfg
void LegacyInferenceToolset::main(int argc, char* argv[]) {
	// parse input arguments
	const char* cfgfile = argv[2];						// path to configuration file
	const char* fpatt = argc > 3 ? argv[3] : "";		// optional pattern
	const char* params = argc > 4 ? argv[4] : "";		// optional string of additional configuration values

	ConfigDictionary cfg_all(cfgfile);
	cfg_all.fromString(params);
	cfg_all.updateNamespaceReferences();

	ConfigDictionary cfg;
	// copy values based on namespace hierarhy
	cfg.fromNamespacePriority(cfg_all,1,"inference");

	int layer_index = cfg.getValueInt("layer_index", 0);
	int start_layer = -1, end_layer = -1;

	cfg.getValue(start_layer, "start_layer", layer_index > 0 ? false : true);
	cfg.getValue(end_layer, "end_layer", layer_index > 0 ? false : true);	

	// either ('layer_index') or ('start_layer' and 'end_layer') can be defined but not both
	if (layer_index > 0) {
		// if both start_layer and end_layer are in config and at the same time layer_index is defined, then return error
		if (start_layer > 0 || end_layer > 0) {
			throw custom_libhop_exception(ConfigException, "Error: Ambiguous definition of layer. Found 'layer_index' and 'start_layer','end_layer' in config. Use either 'layer_index' or 'start_layer' and 'end_layer'.");
		} else {
			// if layer_index was define, then start_layer and end_layer should not be in config, so continue normaly
			start_layer = layer_index;
			end_layer = layer_index;	
		}
	}
		
	// validate start and end layer values
	if (start_layer <= 1) {			
		cout  << "Warning: Found start_layer that is <= 1 but create_layern cannot process 1. layer ... setting start_layer to 2" << endl;
		start_layer = 2;
	}
	if (end_layer < start_layer) {
		// DO NOT throw exception so no error is reporter (required by FileImageInferenceJob in rosette/leoparts)
		cout << "Found start_layer > end_layer. Unable to continue processing" << endl;
		return;
	}

	bool save_result = cfg.getValueBool("save_result", true);
	int time_count_split = cfg.getValueInt("time_count_split", end_layer);
	int maxcount = cfg.getValueInt("file_limit", INT_MAX);
	string globlibname = cfg.getValueString("part_lib_name", "");
	//int layer_merge = cfg.getValueInt("layer_merge", 0) - 1;
	bool video_mode = cfg.getValueBool("video_mode", false); // never used !!??
	double video_thresh = cfg.getValueDouble("video_mode_tolerance", 0.2); // never used !!??
	int adjustment_layer = cfg.getValueInt("adjustment_layer", -1) +1; // never used !!??
    bool shape_check = false;

	string srcdir = "";
	string outdir = "";
	string newext = ".lyx";

	vector<int> endlayers(MAX_LAYER_NUMBER, 0);
	vector<pair<bool,string> > save_intermediate_results(MAX_LAYER_NUMBER, make_pair<bool,string>(0,""));

	map<string, part_lib*> libmap;
	map<string, part_lib*>::iterator libmapiter;
	
    K_bin bin(12, 2, 4, 7); 
    bool merge_scales = false;

	part_lib* globlib = part_lib::read(globlibname);

	if (globlib != nullptr) 
		cout << "Warning: library \'" << globlibname << "\' will be used for all layers." << endl;

	vector<layern_creator*> layern_creator_list(end_layer+1, nullptr);

	for (int i = start_layer; i <= end_layer; ++i) {
		string ly_namespace = string("ly") + i;
		string key = string("end_layer") + i;

		// check if should save result for this intermediate layer (only for non-last layer)
		if (i != end_layer) {
			// use only values from lyN namespace e.g. ly3.save_result and ly3.result_extension for layer 3
			bool save_layer_i = cfg.getValueBool(ly_namespace + ".save_result",false);
			string save_layer_i_extention = cfg.getValueString(ly_namespace + ".result_extension", "", save_layer_i);
				
			save_intermediate_results[i].first = save_layer_i;
			save_intermediate_results[i].second = save_layer_i_extention;
		}
		endlayers[i] = cfg.getValueInt(key, i);

		ConfigDictionary cfgi;
		cfgi.fromNamespacePriority(cfg, 1, ly_namespace.c_str());
			
		part_lib* library;

		// use global library or find one we have already loaded and saved to libmap
		if (globlib != nullptr) {
			library = globlib;
		} else {
			string libname = cfgi.getValueString("part_lib_name", "");

			if (libname.empty()) {
				cout << "Error: Key \'part_lib_name\' in namespace \'" << ly_namespace << "\' expected but not found." << endl;
				return;
			}
			if ((libmapiter = libmap.find(libname)) != libmap.end()) 
				library = libmapiter->second;
			else {
				library = part_lib::read(libname);
				if (library == nullptr) {
					cout << "Error: Unable to find or load library \'" << libname  << "\'." << endl;
					return;
				}
				libmap.insert(pair<string, part_lib*>(libname, library));
			}
		}
            
		// construct creator and set proper library
		layern_creator_list[i] = new layern_creator(cfgi);
		layern_creator_list[i]->set_library(library, cfgi);

        if (!shape_check) 
            shape_check = layern_creator_list[i]->shape_check;

		// update src_dir, out_dir and result_extension if needed
		if (i == start_layer && !cfg.isDefined("src_dir")) 
			cfgi.getValue(srcdir, "src_dir");
		else if (i == end_layer && !cfg.isDefined("out_dir"))
			cfgi.getValue(outdir, "out_dir");
	        
		cfgi.getValue(newext, "result_extension");
        merge_scales = merge_scales || (cfgi.getValueInt("scale_merge", 1) > 1);
	}		
		
	if (cfg.isDefined("src_dir")) cfg.getValue(srcdir, "src_dir");
	if (cfg.isDefined("out_dir")) cfg.getValue(outdir, "out_dir");
	if (cfg.isDefined("result_extension")) cfg.getValue(newext, "result_extension");
		
	end_dir(srcdir);
	end_dir(outdir);

	list<list<string> > files;
	clock_t total_time = 0, total_time1 = 0;
	clock_t time, time1, start, end;
	int count = 0;

	if (merge_scales) {
        cout << "Warning: scale merging is enabled for one of the layers, please make sure that " 
                "from_file is set to true and pattern to an appropriate pattern." << endl;
		//file_lists_from_file(files, fpatt, srcdir, cfg.getValueString("pattern", "_*.ly1"));
		if (!cfg.getValueBool("from_file", false) || !file_lists_from_file(files, fpatt, srcdir, cfg.getValueString("pattern", "%s_*.ly1"))) {
			list<string> fl;
			set<string> fls;

			list_directory(fl, srcdir + fpatt);
			for (list<string>::iterator iter = fl.begin(); iter != fl.end(); ++iter) {
				fls.insert(change_extension(*iter, "", "_"));
			}
			file_lists_from_list(files, list<string>(fls.begin(), fls.end()), srcdir, cfg.getValueString("pattern", "%s_*.ly1"));
		}

	} else {
		list<string> fl;

		if (!cfg.getValueBool("from_file", false) || !file_list_from_file(fl, fpatt, srcdir, cfg.getValueString("pattern", "")))
			list_directory(fl, srcdir + fpatt);
		for (list<string>::iterator fliter = fl.begin(); fliter != fl.end(); ++fliter) {
			files.push_back(list<string>());
			files.back().push_back(*fliter);
		}
	}

	// Process files

	for (list<list<string> >::iterator fliter = files.begin(); fliter != files.end(); ++fliter) {
		if (++count > maxcount) break;

		vector<layer1_result*> res;

		read_layer1_result_vector(res, *fliter, srcdir);
		if (res.empty()) 
			continue;

		if (video_mode && video_thresh < 0) {
			int maxmerge = (int)(-video_thresh);
			int m = 1;

			while (m < maxmerge && ++fliter != files.end()) {
				vector<layer1_result*> newres;

				cout << endl;
				read_layer1_result_vector(newres, *fliter, srcdir);
				if (newres.empty()) 
					continue;

				for (int i = 0; i < (int)newres.size(); ++i) {
					if (i < (int)res.size()) 
						res[i]->merge(newres[i], 100);  // 100 !!!!
					delete newres[i];
				}
				++m;
			}
			while (--m > 0) --fliter;
		} else if (video_mode && video_thresh > 0) {
			layer1_result* refres = (layer1_result*)res[0]->get_copy_s();

			while (++fliter != files.end()) {
				vector<layer1_result*> newres;

				cout << endl;
				read_layer1_result_vector(newres, *fliter, srcdir);
				if (newres.empty()) 
					continue;

				int diff = layer_0_difference(refres, newres[0], 0.5); // 0.5 !!!!!

				if ((double)diff/refres->shape_nodes[0].size() > video_thresh) {
					--fliter;
					break;
				}
				cout << 'M';
				for (int i = 0; i < (int)newres.size(); ++i) {
					if (i < (int)res.size()) 
						res[i]->merge(newres[i], 100);  // 100 !!!!
					delete newres[i];
				}
			}
			delete refres;
		}
		if (fliter == files.end()) --fliter;

		// g-response adjustment section
		// parameters: 
		//    adjustment_layer: default 0, i.e. no adjustment
		//    adjustment_average: average number of states in the same position, default = 0.0
		//    adjustment_count: if adjustment_average is set then it has no effect, 
		//       otherwise it is a threshold for the number of hits.

		if (adjustment_layer >= start_layer && adjustment_layer <= end_layer) {
			double adjustment_average = cfg.getValueDouble("adjustment_average", 0.0);
			int adjustment_count = cfg.getValueInt("adjustment_count", 1000);

			start = clock();
			int i = start_layer; 
			while (i <= end_layer) {
				// bool mem = layern_creator_list[i].add_activation_edges;
				//
				// layern_creator_list[i].add_activation_edges = false;
				layern_creator_list[i]->add_layer(res, i, 0);
				//layern_creator_list[i].add_activation_edges = mem;
				cout << '[' << i << ']';

				if (i != adjustment_layer) {
					++i;
				} else {
					iipair stat = layer_node_count(res[0], i - 1);
	                    
					if (adjustment_average > 0 && (stat.second == 0 || (double)stat.first/stat.second < adjustment_average) ||
							stat.first < adjustment_count) {
						do {
							double& gthresh = layern_creator_list[i]->candidate_g_threshold_percent;

							gthresh *= 0.5;      // or something else
							if (gthresh > 0.01)  // ~0
								break;
							--i;
							gthresh /= 0.5;      // or something else ^^^
							//res->delete_layers_geq(i - 1);
						} while (i >= start_layer);
						if (i < start_layer)
							break;
					} else {
						++i;
					}
	                        
				}

				// save intermediate results
				if (save_intermediate_results[i].first == true) {
					list<string>::iterator fiter = fliter->begin();
					string ly_ext = save_intermediate_results[i].second;

					for (int i = 0; i < (int)res.size(); ++i) {
						save_layer1_result(res[i], outdir + change_extension(*fiter, ly_ext));
						++fiter;
					}
				}
			}
			end = clock();
			time = end - start;
			total_time += time;
			cout << " (time = " << (double)time/CLOCKS_PER_SEC << " sec)" << endl; 
		} else {
			time = 0; time1 = 0;
            vector<scmap_t> scmap(res.size());

            if (shape_check) {
                for (int i = 0; i < (int)res.size(); ++i) 
                    get_sc_map(scmap[i], res[i], bin, true);
            }

			for (int i = start_layer; i <= end_layer; ++i) {
				start = clock();
				layern_creator_list[i]->add_layer(res, i, endlayers[i]);
				end = clock();
				time += end - start;
                if (!res.empty())
					cout << '(' << res[0]->cover_quotient(i - 2) << ')' << '[' << (double)(end - start)/CLOCKS_PER_SEC << " sec]";
				if (i == time_count_split) time1 = time;  
				i = endlayers[i];
					
				// save intermediate results
				if (save_intermediate_results[i].first == true) {
					list<string>::iterator fiter = fliter->begin();
					string ly_ext = save_intermediate_results[i].second;

					for (int i = 0; i < (int)res.size(); ++i) {
						save_layer1_result(res[i], outdir + change_extension(*fiter, ly_ext));
						++fiter;
					}
				}
			}
			total_time += time;
			total_time1 += time1;
			cout << " (time = " << (double)time/CLOCKS_PER_SEC << " = " << 
				(double)time1/CLOCKS_PER_SEC << " + " << (double)(time - time1)/CLOCKS_PER_SEC << " sec)" << endl;
		}

		//cout << endl << "saving result" << endl;
		if (save_result) {
			list<string>::iterator fiter = fliter->begin();

			for (int i = 0; i < (int)res.size(); ++i) {
				save_layer1_result(res[i], outdir + change_extension(*fiter, newext));
				++fiter;
			}
		}
		//cout << "deleting result" << endl;
		for (int i = 0; i < (int)res.size(); ++i) 
			delete res[i];
		//cout << "done with this one" << endl;
	}
	    
	cout << "Total time = " << (double)total_time/CLOCKS_PER_SEC << " sec;";
	if (!files.empty()) {
		cout << " = " << (double)total_time/CLOCKS_PER_SEC/min<int>((int)files.size(), maxcount) << " sec/file";
		cout << " = " << (double)total_time1/CLOCKS_PER_SEC/min<int>((int)files.size(), maxcount) << " + ";
		cout << (double)(total_time - total_time1)/CLOCKS_PER_SEC/min<int>((int)files.size(), maxcount);
	}
	cout << endl;
	for (int i = start_layer; i <= end_layer; ++i) {
		layern_creator_list[i]->set_library(nullptr);
		delete layern_creator_list[i];
	}
	for (map<string, part_lib*>::iterator miter = libmap.begin(); miter != libmap.end(); ++miter) {
		delete miter->second;
	}
    if (globlib != nullptr) 
        delete globlib;
}