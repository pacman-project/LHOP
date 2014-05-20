// laydisplay.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cctype>
#include <tuple>

#include <opencv2/highgui/highgui.hpp>

#include "utils/hopmath.h"
#include "utils/graphs/img_graph.h"
#include "utils/graphs/graph_utils.h"

#include "utils/convert.h"
#include "utils/matrix.h"
#include "utils/utils.h"
#include "utils/ocv.h"

#include "core/legacy/layer_1_result.h"

#include "modules/histogram-of-compositions/hoc.h"

#include "tools/main_toolset.h"

#ifdef WIN32
#define TOLOWER std::tolower
#define TOUPPER std::toupper
#else
#define TOLOWER ::tolower
#define TOUPPER ::toupper
#endif


void use_library_thresholds(layer1_result* res, part_lib* library, int z, double svmdelta)
{	
    if (library == nullptr) 
        return;
	int edgename = EdgeConnection::TO_PREV_LAYER;

	for (list<node*>::iterator iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
		node* n = *iter;
		layer1_data* nd = (layer1_data*)n->data;

		if (n->is_attr_set(NODE_DELETED_ATTR))
			continue;
        if (nd->z == z) {
            node* nn = n->get_neighbor(edgename);
            layer1_data* nnd = (layer1_data*)nn->data;

            if (nnd->z > library->max_layer_index() || nnd->m >= library->layer_size(nnd->z)) {
                cout << "Library is not compatible with result!";
                throw;
            }

            bool ok = true;
            node* p = library->parts[nnd->z][nnd->m];
            vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);
            double svmdfval = 0.0;

            ok = p->is_attr_set(NODE_DELETED_ATTR) == false;

            if (ok && pd != nullptr && pd->svmt.svm != nullptr) {
                int nchildren = p->count_neighbors(EdgeConnection::TO_LYR_SOURCE);

                svmdfval = check_svmt(pd->svmt, nchildren, nn, true);

                ok = svmdfval > svmdelta;
                svmdfval = svmdfval - svmdelta;
                //cout << 'f' << f;
                //benergy = f + 0.51;
            }
            if (!ok) {
                n->set_attr(NODE_DELETED_ATTR);
                nn->set_attr(NODE_DELETED_ATTR);
            }

		}
	}		 
    
}

void save_boxes(layer1_result* res, int layer_index, double inhibitboxes, const string& out)
{
    list<layer1_result::box_data_t> boxes;
    ofstream os(out.c_str());

    res->get_boxes(boxes, layer_index, set<int>(), inhibitboxes);

    os << irectangle2(0, 0, res->x_size(0), res->y_size(0)) << endl;
    for (list<layer1_result::box_data_t>::iterator iter = boxes.begin(); iter != boxes.end(); ++iter) {
        os << iter->box << endl;
    }
                
}

void save_normal(ConfigDictionary& cfg, const string& fname, string out)
{
    layer1_result* res;
    img* im = nullptr;
    int layer_index = cfg.getValueInt("layer", -1);
    int end_layer_index = cfg.getValueInt("end_layer", 0);
    bool drawbox = cfg.getValueBool("draw_box", false);
    string mode = cfg.getValueString("mode", "r"); 
    double factorpow = cfg.getValueDouble("factor_power", 1.0);
    bool alltypes = cfg.getValueBool("all_types", false);
    bool inhibitboxes = cfg.getValueBool("inhibit_boxes", false);
    vector<int> parts, icolors;
    vector<color> colors;
    color defcolor, uncovered;
    bool paintuncov;

    cfg.getValue(parts, "parts");
    cfg.getValue(icolors, "part_colors");
    color::to_color_vector(colors, icolors, 0);

    icolors.clear();
    cfg.getValue(icolors, "default_color");
    color::to_color(defcolor, icolors, 255);

    icolors.clear();
    cfg.getValue(icolors, "uncovered_color");
    paintuncov = (icolors.size() != 0);
    color::to_color(uncovered, icolors, 20);

    read_layer1_result(res, fname);
    
    if (out == "") {
        string::size_type dot_pos;
    
        if ((dot_pos = fname.rfind(".")) == string::npos) out = fname + ".bmp";
        else out = fname.substr(0, dot_pos) + ".bmp";
    }
    
    switch (mode[0]) {
        case 'N' : 
            im = res->get_image(layer_index, end_layer_index, paintuncov, uncovered, defcolor, parts, colors); 
            break;
        case 'I' : 
            im = res->get_image_inhib(layer_index, end_layer_index, paintuncov, uncovered, defcolor, parts, colors); 
            break;
        case 'B' :
			{
				part_lib* library;
				string libname = cfg.getValueString("part_lib", "");
				library = part_lib::read(libname);

				im = res->get_image_boxed(library, layer_index, parts, colors, defcolor, inhibitboxes);
            
				if (cfg.getValueBool("save_boxes", false)) 
					save_boxes(res, layer_index, cfg.getValueDouble("box_inhibition_percent", 0.1),
						change_extension(out, ".bb"));

				if (library != nullptr) delete library;
			}
            break;
        case 'X' :
            if (layer_index < (int)res->shape_nodes_inhib.size()) {
                vector<node*>& nodes = res->shape_nodes_inhib[layer_index];
                char str[5];
                int i = 0;

                for (vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
                    im = res->get_part_image(layer_index, end_layer_index, *iter);
                    if (im) {
                        sprintf(str, "%d", i++);
                        string name = change_extension(out, string(str) + ".bmp");
                        im->save_normalized(name.c_str());
                        delete im;
                    }
                }
                cout << endl;
                im = nullptr;
            }
            break;
        case 'Y' :
            {
                part_lib* library;
                string libname = cfg.getValueString("part_lib", "");

                library = part_lib::read(libname);
                im = res->get_image_reconstructed(library, layer_index, parts);
                if (library) delete library;
            }
            break;
        default : 
			part_lib* library;
            string libname = cfg.getValueString("part_lib", "");
			library = part_lib::read(libname);

            if (!cfg.getValueBool("animation", false) || parts.empty()) {
                im = res->get_image_reconstructed(library, layer_index, end_layer_index, parts, false, drawbox, factorpow, alltypes);
                if (cfg.getValueBool("draw_groundtruth", false)) {
                    
					list<pair<irectangle2,string> > rect_list;
					read_groundtruth(rect_list, fname);

					for (list<pair<irectangle2,string> >::iterator iter = rect_list.begin();
							iter != rect_list.end(); iter++) {
						pair<irectangle2,string>& p = *iter;
						im->draw_box(p.first, 1.0);
					}
                }
            } else {
                vector<int> dispparts;
                char str[5];

                for (int i = 0; i < (int)parts.size(); ++i) {
                    dispparts.push_back(parts[i]);
                    im = res->get_image_reconstructed(library, layer_index, end_layer_index, dispparts, 
                        false, drawbox, factorpow);
                    sprintf(str, "%d", i);
                    if (im) {
                        string name = change_extension(out, string(str) + ".bmp");
                        cout << name << endl;
                        im->save_normalized(name.c_str());
                        delete im;
                    }
                }
                im = nullptr;
            }
			if (library != nullptr) delete library;
    }
    if (im) {
        im->save_normalized(out.c_str());
        delete im;
    }
}


void save_visview(ConfigDictionary& cfg, const string& fname, const string& out)
{
    layer1_result* res = nullptr;
    part_lib* library = nullptr;
    string outdir = cfg.getValueString("out_dir", "");
    string srcdir = cfg.getValueString("src_dir", "");
    string libname = cfg.getValueString("library", "");
    int save_mode = cfg.getValueInt("save_mode", 0);
    int layer_index = cfg.getValueInt("layer", -1);
    bool save_sc = cfg.getValueBool("save_sc", true);
    bool save_var = cfg.getValueBool("save_variations", false);
    double thresh = cfg.getValueDouble("thresh", 0.3);
    int resp = response_from_string(cfg.getValueString("response", "RR_RESPONSE"));
    double svmdelta = cfg.getValueDouble("svm_delta", 0.0);
    bool libthresh = cfg.getValueBool("use_library_thresholds", false);

    if (cfg.getValueBool("true_coordinates", false))
        save_mode |= VVE_TRUE_COORDINATES;
    if (cfg.getValueBool("bounding_boxes", false))
        save_mode |= VVE_BOUNDING_BOXES;
    if (cfg.getValueBool("simple_filenames", false))
        save_mode |= VVE_SIMPLE_FILENAMES;
	if (cfg.getValueBool("matlab_scale_indexes", false))
        save_mode |= VVE_MATLAB_SCALES;

    end_dir(outdir);
    end_dir(srcdir);

    if (libname.size() > 0) {
        library = part_lib::read(libname);

        if (library) {
            library->save_visview_filters(outdir);
            library->save_visview_centers(outdir);
            library->save_visview_parts(outdir); 
            library->save_visview_contractions(outdir);
            library->save_sc(outdir, false, save_sc, save_var);
        }
    }

    istringstream iss(out);
    int outindex;
    set<double> scales;
    set<string> bbgt;

    iss >> outindex;
    if (iss.fail()) outindex = 1;

    list<string> files;

    list_directory(files, srcdir + fname);

	string filename;
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        cout << "Processing: " << *fiter << endl;
        
        read_layer1_result(res, srcdir + *fiter);

        if (res) {
            ostringstream oss;
            string basefname = change_extension(*fiter, "", "_");
            double scale = round((double)(res->x_size() - 2*res->border)/res->original_width, 2);
            list<pair<irectangle2, string> > rectangles;

            read_groundtruth(rectangles, srcdir + *fiter);
            oss << outindex++;
            if (libthresh) 
                use_library_thresholds(res, library, library->max_layer_index(), svmdelta);
			filename = change_extension(*fiter, "");
            res->save_visview(outdir, oss.str(), filename, library, layer_index, save_mode);
            scales.insert(scale);
            if (bbgt.find(basefname) == bbgt.end()) {
                for (auto riter = rectangles.begin(); riter != rectangles.end(); ++riter) {
                    riter->first -= ipoint2(res->border, res->border);
                    riter->first.ll /= scale;
                    riter->first.ur /= scale;
                }
                save_groundtruth(rectangles, outdir, basefname + "_bbgt.txt", "", "", false);
                bbgt.insert(basefname);
            }
            delete res;
        }
    }

    // save_scales
    vector<double> sscales(scales.begin(), scales.end());
    ofstream sfs((outdir + "scales.txt").c_str());

    sort(sscales.begin(), sscales.end(), greater<double>());
    for (int i = 0; i < (int)sscales.size(); ++i) {
        if (i == 0 || sscales[i - 1]/sscales[i] > 1.1) 
            sfs << sscales[i] << endl;
    }
    sfs.close();

    // save classnames
    if (library) {
        set<string> names;

        for (auto piter = library->nodes.begin(); piter != library->nodes.end(); ++piter) {
            node* p = *piter;

            if (p->is_attr_set(CATEGORY_NODE_ATTR)) {
                cpart_data* pd = (cpart_data*)p->data;
                names.insert(pd->name);
            }
        }
        
        ofstream cfs((outdir + "classnames.txt").c_str());

        for (auto citer = names.begin(); citer != names.end(); ++citer) {
            cfs << *citer << endl;
        }
        cfs.close();
    }

    if (library) delete library;
}


void save_library(ConfigDictionary& cfg, const string& fname, string out)
{
    int layer = cfg.getValueInt("layer", 1);
	part_lib* library = part_lib::read(fname);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }
    if (layer >= (int)library->parts.size()) { 
        cout << "Layer " << layer << " not found in " << fname << "." << endl;
        return;
    }

    if (cfg.getValueBool("add_similarity_edges", false)) {
        part_lib::similarity_threshold = cfg.getValueDouble("threshold", 0.1);
        double contraction = cfg.getValueDouble("contraction", 0.5);

        cout << "Adding similarity edges with threshold " << part_lib::similarity_threshold << " and" << endl;
        cout << "contraction factor " << contraction << endl;
        //library->add_similarity_edges(layer - 1, contraction);
        cout << "Saving library..." << endl;
        library->save(out);
    } else if (cfg.isDefined("part")) {
        int part = cfg.getValueInt("part", -1);
        cout << "Saving part #" << part << " to " << out << "..." << endl;
        library->save_part(out.c_str(), layer, part, cfg.getValueDouble("threshold", 0.0));
    } else if (cfg.isDefined("indices")) {
        vector<int> l;

        cfg.getValue(l, "indices", true);
        cout << "Saving " << fname << " to " << out << endl;
        library->save_all(out.c_str(), layer, l, cfg.getValueBool("show_labels", true), cfg.getValueBool("mark_center", false));
    } else {
        cout << "Saving " << fname << " to " << out << "." << endl;
        vector<double> contr;

        cfg.getValue(contr, "layer_contractions");
        part_data::set_contractions(contr);
        library->save_all(out.c_str(), layer, 0, -1, cfg.getValueBool("show_labels", true), false, 
            cfg.getValueBool("mark_center", false), cfg.getValueBool("one_row", false));
    }
    delete library;
}

void save_library_sc(ConfigDictionary& cfg, const string& fname, string out)
{
    int layer = cfg.getValueInt("layer", 1);
    part_lib* library = part_lib::read(fname);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }
    if (layer < 1 || layer - 1 > library->max_layer_index()) { 
        cout << "Layer " << layer << " not found in " << fname << "." << endl;
    } else {
        cout << "Saving " << fname << " to " << out << "." << endl;
        library->save_all_sc(out.c_str(), layer, 0, -1, cfg.getValueBool("show_labels", true));
    }
    delete library;
}

void save_library_mean(ConfigDictionary& cfg, const string& fname, string out)
{
    int layer;
    part_lib* library = part_lib::read(fname);
    string pattern;
    bool savev = cfg.getValueBool("vectors", false);

    cfg.getValue(layer, "layer", true);
    pattern = cfg.getValueString("pattern", "part-%d-%03d.txt");
	
    end_dir(out);
    if (library == nullptr) {
        cout << "Can not open library." << endl;
        return;
    }
    if (layer < 0 || layer > library->max_layer_index()) {
        cout << "Layer " << layer << " does not exist." << endl;
    } else {
        for (int pi = 0; pi < library->parts[layer].size(); ++pi) {
            node* p = library->parts[layer][pi];
            part_data* pd = (part_data*)p->data;
            vs_part_data* vspd = dynamic_cast<vs_part_data*>(p->data);
            char s[255];

            sprintf(s, pattern.c_str(), layer, pi);

            ofstream os((out + s).c_str());
            vector<dpoint2> dpts;

            if (savev) {
                for (int i = 0; i < vspd->pcad.mean.cols; i += 2) {
                    os << vspd->pcad.mean.at<double>(0, i) << ',' << vspd->pcad.mean.at<double>(0, i + 1);
                    for (int j = 0; j < vspd->pcad.eigenvectors.rows; ++j) {
                        os << ',' << vspd->pcad.eigenvectors.at<double>(j, i) << ',' 
                            << vspd->pcad.eigenvectors.at<double>(j, i + 1);
                    }
                    os << '\n';
                }
            } else {
                if (vspd != nullptr) {
                    dpts = partition(vspd->pcad.mean);
                } else {
                    vector<ipoint2> ipts = get_library_geo(p);

                    dpts = cast_vector<dpoint2, ipoint2>(ipts);
                }
                translate_and_scale(dpts);
                for (auto piter = dpts.begin(); piter != dpts.end(); ++piter) {
                    dpoint2 p = (*piter)*50.0;

                    os << (int)(p.x) << ',' << (int)(p.y) << '\n';
                }
            }
        }
    }
    
    delete library;
}

void save_library_sc_mma(ConfigDictionary& cfg, const string& fname, string out)
{
    int layer = cfg.getValueInt("layer", 1);
    part_lib* library = part_lib::read(fname);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }
    if (layer < 1 || layer - 1 > library->max_layer_index()) { 
        cout << "Layer " << layer << " not found in " << fname << "." << endl;
    } else {
        cout << "Saving " << fname << " to " << out << "." << endl;
        library->save_all_sc_mma(out.c_str(), layer, 0, -1);
    }
    delete library;
}



void save_features(ConfigDictionary& cfg_root, const string& fpatt, string outname)
{
	ConfigDictionary cfg;
	cfg.fromNamespacePriority(cfg_root, 1, "hoc");
	
	string g_outdir = outname;
	string outdir = outname;
	string srcdir = cfg.getValueString("src_dir", "");

	bool read_windows_from_file = cfg.getValueBool("read_windows_from_file", false);	
	string gtext = cfg.getValueString("groundtruth_extension", ".groundtruth");

	string libname;
	cfg.getValue(libname, "library", true);

	part_lib* library = part_lib::read(libname);
    if (library == nullptr) {
        cout << "Can not open library (library is needed to extract the number of features)." << endl;
        return;
    }

	using namespace hoc;

	hoc_histogram_generator hoc_generator(library);

	hoc_generator.init_cfg(cfg);
	

	end_dir(g_outdir);
	end_dir(outdir);

    int count = 0;

	list<string> files; 

    end_dir(srcdir);
    if (!cfg.getValueBool("from_file", false) || !file_list_from_file(files, fpatt, srcdir, cfg.getValueString("pattern", "")))
        list_directory(files, srcdir + fpatt);

	// since we will be appending descriptors from different scale, we need to clear any existing values in all output files
	for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
		string::size_type scale_start = (*fiter).rfind("_");
		string name = (*fiter).substr(0, scale_start);
		
		stringstream ss;
		ss << outdir << name << "_det.txt";

		fclose(fopen(ss.str().c_str(),"w"));
	}

    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *fiter);
        if (res == nullptr) {
            cout << "Can not load '" << *fiter << "'" << endl;
        } else {
            cout << "Processing " << *fiter << "...";
			clock_t startt, endt;
			startt = clock();
			{

				string::size_type ext_start = (*fiter).rfind(".");
				string::size_type scale_start = (*fiter).rfind("_");

				string name = (*fiter).substr(0, scale_start);
				string scale_str = (*fiter).substr(scale_start+1);
				
				int scale;
				sscanf(scale_str.c_str(), "%d%*s", &scale);

				stringstream ss;
				ss << outdir << name << "_det.txt";

				ofstream os(ss.str().c_str(), ios_base::out | ios_base::app);

				list<pair<irectangle2, string> > rectangles;
				read_groundtruth(rectangles, srcdir + *fiter, "", gtext, false);

				int border = res->border;
				double scale_factor = res->original_width / ((double)(res->x_size(0) - 2*border));

				// Save (transformed) groundtruths
				if (!rectangles.empty()) {
					
					list<pair<irectangle2, string> > rectangles2(rectangles);
					ipoint2 bpt(border, border);

					stringstream gt_name;
					gt_name<< outdir << name << "_gt.txt";

					ofstream g_os(gt_name.str().c_str());

					for (list<std::pair<irectangle2,string> >::iterator iter = rectangles2.begin(); iter != rectangles2.end(); ++iter) {
						std::pair<irectangle2,string> p = *iter;
						irectangle2& r = p.first;
						string& cat_name = p.second;
				        
						r = (r - bpt) * scale_factor;					

						g_os << r.ll.x << ' ' << r.ll.y << ' ' << r.ur.x - r.ll.x << ' ' << r.ur.y - r.ll.y << ' ' << cat_name << endl;
					}
					g_os.close();

				}
				std::list<irectangle2> annotations_rects;

				if (read_windows_from_file) { 
						
					const string& annot_name =  srcdir + name + "_annot.txt";
					
					ifstream is(annot_name.c_str());

					if (is.is_open()) {
						while (true) {
							string line;
							// read and parse whole line as stream of string
							getline(is,line);
							stringstream sline(line);
							
							float x,y,w,h;
							sline >> x >> y >> w >> h;
							
							if (is.fail()) break;

							annotations_rects.push_back(irectangle2(x, y, x + w, y + h));
							
						}
					} else {
						cout << "ERROR: Unable to open annotation file '" << annot_name << "'. " << endl;
						return;
					}
				}
				// process image (res) for specific regions to produce list of descriptors (histograms of compositions)
				std::vector<histogram_descriptor>* result = hoc_generator.generate_descriptors(res, annotations_rects);

				for (std::vector<histogram_descriptor>::iterator desc = result->begin(); desc != result->end(); desc++) {
					
					int hist_size = (int)(*desc).hist.size();

					// skip any hitograms that have all zeros
					bool only_zeros = true;
					for (int f = 0; f < hist_size; ++f) {
						if ((*desc).hist[f] > 0) {
							only_zeros = false;
							break;
						}
					}

					if (only_zeros) 
						continue; // skiping since has all zeros (it is useless descriptor)

					os << scale << ' ' << (*desc).x << ' ' << (*desc).y << ' ' << (*desc).w << ' ' << (*desc).h << ' ';
					for (int f = 0; f < hist_size; ++f) {
						os << (*desc).hist[f] << ' ';
					}
					os << '\n';					
				}
				delete result;
				os.close();
            }
			endt = clock();
			cout << " (time = " << (double)(endt - startt)/CLOCKS_PER_SEC << ")" << endl;
            cout << endl;
            delete res;
            ++count;
			
        }
    }
	delete library;
	
}

void init_display(const char* file, const char* cfg_file, const char* outname, const char* params)
{

    cout << "laydisplay (" __DATE__ " / " __TIME__ ")" << endl;
    try {
        string cff = cfg_file;
        string mode;
    
        if (cff == "") cff = "c:\\programs\\laydisplay.cfg";
        ConfigDictionary cfg(cff);
        
        cfg.fromString(params);
        mode = cfg.getValueString("mode", "r");
        transform(mode.begin(), mode.end(), mode.begin(), TOUPPER);
        if (mode == "V") save_visview(cfg, file, outname);
        else if (mode == "L") save_library(cfg, file, outname);
        else if (mode == "LSC") save_library_sc(cfg, file, outname);
        else if (mode == "LSCM") save_library_sc_mma(cfg, file, outname);
        else if (mode == "LMEAN") save_library_mean(cfg, file, outname);
		else if (mode == "FEAT") save_features(cfg, file, outname);
        else save_normal(cfg, file, outname);

    } catch (const libhop_exception& e) {
		cout << e.what() << endl;
    } catch (const exception& e) {
		cout << "General exception '" << typeid(e).name() << "' with message: '" << e.what() << "'" << endl;
	} 	
}

void print_usage() 
{
    cout << "hopdisplay (" __DATE__ " / " __TIME__ ")" << endl;
    cout << "Usage: hopdisplay fname [config file] [out name]" << endl;
    cout << "config file options:" << endl;
    cout << "   mode            = n|i|r|v|s|l|m|p|c" << endl;
    cout << "       (normal, inhibited, reconstruction (default), visview, SVM, library, " << endl;
    cout << "           Mathematica, part statistics, part covering statistics)" << endl;
    cout << "   layer           = index of layer to display (-1 means the last one (default))" << endl;
    cout << "   parts           = list of parts to display (not specified means all)" << endl;
    cout << "   part_colors     = list of colors for parts" << endl;
    cout << "       (list of triples form [0,..,255]^3 specifying rgb color)" << endl;
    cout << "   default_color   = default color for parts" << endl;
    cout << "   uncovered_color = color for uncovered parts (if not specified, uses" << endl;
    cout << "       color for parts)" << endl;
    cout << "   out_dir         = output directory (visview)" << endl;
    cout << "   scales          = list of scales to include to SVM export" << endl;
    cout << " " << endl;
}



int main(int argc, char* argv[])
{
    ///_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
    ///_CrtSetBreakAlloc(326330);
    //init_streaming();   
    random_seed((int)time(nullptr));

    switch (argc) {
        case 2: init_display(argv[1], "", "", ""); break;
        case 3: init_display(argv[1], argv[2], "", ""); break;
        case 4: init_display(argv[1], argv[2], argv[3], ""); break;
        case 5: init_display(argv[1], argv[2], argv[3], argv[4]); break;
        default: print_usage();
    }

	return 0;
}

