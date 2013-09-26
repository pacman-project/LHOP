/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */

#ifdef WIN32
#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
#include <windows.h>
#endif

#include "hop.h"

#include <stdio.h>
#include <iostream>
#include <cctype>

#include "utils/utils.h"
#include "utils/img.h"
#include "layers/layer_1_creators.h"

namespace layer1 {

	#ifdef WIN32
	#define TOLOWER std::tolower
	#define TOUPPER std::toupper
	#else 
	#define TOLOWER ::tolower
	#define TOUPPER ::toupper
	#endif

	// Extract rectangular patches from the image (listed in groundtruth file).
	// Image is resized s.t. (for each gt. rectangle) the maximal dimension is 
	// equal to max_max_dimension, scales are created until the maximal dimension is
	// equal to min_max_dimension. By default min_max_dimension = max_max_dimension.
	void create_and_save_extraction(const config_dictionary& cfg, const char* pattern)
	{
		string srcdir = cfg.get_value_string("src_dir", "");
		string outdir = cfg.get_value_string("out_dir", "");
		string outprefix = cfg.get_value_string("out_prefix", "");
		string gtext = cfg.get_value_string("groundtruth_extension", ".groundtruth");
		string catname = cfg.get_value_string("category_name", "");
		int gtdelta = cfg.get_value_int("groundtruth_delta", 0);
		int border = cfg.get_value_int("border_size", 0);
		bool separate_colors = cfg.get_value_bool("separate_colors", false);
        string maskext = cfg.get_value_string("mask_extension", "");
        int maskthick = cfg.get_value_int("mask_thickening", 0);
		int maxdim;

		cfg.get_value(maxdim, "max_max_dimension", true); 
	    
		int mindim = cfg.get_value_int("min_max_dimension", maxdim);
	   
		list<string> files;

		end_dir(srcdir);
		end_dir(outdir);
		if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, pattern, srcdir))
			list_directory(files, srcdir + pattern);
		for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
			vector<img> srcimg;
            img maskimg;

			if (separate_colors) {
				img red, green, blue;
				img::read_colors(srcdir + *fiter, red, green, blue);
				
				srcimg.push_back(red);
				srcimg.push_back(green);
				srcimg.push_back(blue);
			} else {
				srcimg.push_back(img(srcdir + *fiter));
			}

			cout << "Processing " << *fiter;
			if (srcimg.size() < 0 || srcimg.front().width == 0 || srcimg.front().height == 0) {
				cout << " error reading file!" << endl;
				continue;
			}

            if (!maskext.empty()) {
                try {
                    img::read_grayscale(srcdir + change_extension(*fiter, maskext), maskimg);
                    if (maskthick > 0) maskimg.thicken_white(maskthick);
                } catch (...) {
                    cout << "Mask extennsion specified but mask file \"" << srcdir + change_extension(*fiter, maskext) << 
                        "\" can not be loaded." << endl;
                }
            }

			list<pair<irectangle2, string> > rectangles;
			int gtcount = 0;

			read_groundtruth(rectangles, srcdir + *fiter, catname, gtext);
			if (rectangles.empty()) {
				cout << " groundtruth file does not exist, skipping creation." << endl;
				continue;
			}
			for (list<pair<irectangle2, string> >::iterator riter = rectangles.begin(); riter != rectangles.end(); ++riter) {
				irectangle2 r = riter->first;

				r.grow(gtdelta);

				vector<img> srcimg2(srcimg.size());
                img maskimg2 = maskimg.cut(r);
				for (int i = 0; i < srcimg.size(); i++)
					srcimg2[i] = srcimg[i].cut(r);
                

				config_dictionary cfg2(cfg); // just change the dimensions...
				vector<layer1_result*> result;
				string maxdimstr = string("") + maxdim;
				string mindimstr = string("") + mindim;

				cfg2.set_value("init_size", maxdimstr);
				cfg2.set_value("scale_limit", mindimstr);

				lay1create::create_layer1(result, cfg2, srcimg2, maskimg2);

				// Set metadata and save images
				for (int i = 0; i < (int)result.size(); ++i) {
					string ext = string("_") + gtcount + string("_") + i + string(".ly1");

					result[i]->original_width = srcimg.front().width;
					result[i]->border = border;
					save_layer1_result(result[i], change_extension(outdir + outprefix + *fiter, ext));
					delete result[i];
				}
				gtcount++;
			}

			cout << endl;
		}
	}

	irectangle2 min_rectangle(const list<pair<irectangle2, string> >& rectangles)
	{
		irectangle2 result;

		for (list<pair<irectangle2, string> >::const_iterator riter = rectangles.begin(); riter != rectangles.end(); ++riter) {
			irectangle2 r = riter->first;

			if (result.invalid() || max(r.x_dim(), r.y_dim()) < max(result.x_dim(), result.y_dim()))
				result = r;
		}
		return result;
	}

	irectangle2 max_rectangle(const list<pair<irectangle2, string> >& rectangles)
	{
		irectangle2 result;

		for (list<pair<irectangle2, string> >::const_iterator riter = rectangles.begin(); riter != rectangles.end(); ++riter) {
			irectangle2 r = riter->first;

			if (result.invalid() || max(r.x_dim(), r.y_dim()) > max(result.x_dim(), result.y_dim()))
				result = r;
		}
		return result;
	}

	void split_image(list<vector<img> >& srcimglist, list<list<pair<irectangle2, string> > >& rectlist, 
		const vector<img>& srcimg, const list<pair<irectangle2, string> >& rectangles, int n, int overlap, double rthresh)
	{
		if (srcimg.size() <= 0) {
			cout << "Error in split_image(): Unable to split image when there is no valid srcimg." << endl;
			return;
		}

		if (n <= 1) {
			srcimglist.push_back(srcimg);
			rectlist.push_back(rectangles);
			return;
		}

		int img_width = srcimg.front().width;
		int img_height = srcimg.front().height;
		overlap /= 2;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				int xmin = max<int>(0, i*img_width/n - overlap);
				int xmax = min<int>(img_width, (i + 1)*img_width/n + overlap);
				int ymin = max<int>(0, j*img_height/n - overlap);
				int ymax = min<int>(img_height, (j + 1)*img_height/n + overlap);
				irectangle2 rect(xmin, ymin, xmax, ymax);

				vector<img> srcimg_cut(srcimg.size());
				for (int k = 0; k < srcimg.size(); k++)
					srcimg_cut[k] = srcimg[k].cut(rect);

				srcimglist.push_back(srcimg_cut);
				rectlist.push_back(list<pair<irectangle2, string> >());
				for (list<pair<irectangle2, string> >::const_iterator iter = rectangles.begin(); iter != rectangles.end(); ++iter) {
					irectangle2 isection = rect.intersection(iter->first);

					if (!isection.invalid() && (double)isection.area()/iter->first.area() >= rthresh) {
						rectlist.back().push_back(pair<irectangle2, string>(isection - ipoint2(xmin, ymin), iter->second));
					}
				}
			}
		}
	}

	string replaceChar(const string& str, const char *c, char r)
	{
		string ret=str;
		for(int i=ret.find_first_of(c); i!=string::npos; i=ret.find_first_of(c,i+1)) {
			ret = ret.replace(i,1,1,r);
		}
		return ret;
	}

	void save_results(vector<layer1_result*>& result, list<pair<irectangle2, string> >& rectangles,
		const string& name, const string& outdir, const string& outprefix, const string& suffix)
	{
		string fname = replaceChar(name,"\\/",'.');
		for (int i = 0; i < (int)result.size(); ++i) {
			string istr = string("") + i;

			if (i > 0) cout << ' ';
			cout << istr;

			// Update metadata save "result"
			save_layer1_result(result[i], change_extension(outdir + outprefix + fname, suffix + "_" + istr + ".ly1"));

			// Save (transformed) groundtruths
			if (!rectangles.empty()) {
				int border = result[i]->border;
				double factor = ((double)(result[i]->x_size(0) - 2*border))/result[i]->original_width;
				list<pair<irectangle2, string> > rectangles2(rectangles);
				ipoint2 bpt(border, border);

				for (list<pair<irectangle2, string> >::iterator riter = rectangles2.begin(); riter != rectangles2.end(); ++riter) {
					riter->first = (riter->first * factor) + bpt;
				}
				save_groundtruth(rectangles2, outdir, fname, outprefix, suffix + "_" + istr);
			}

			delete result[i];
		}

	}

	void flip_groundtruth(list<pair<irectangle2, string> >& rectangles, int width)
	{
		for (list<pair<irectangle2, string> >::iterator iter = rectangles.begin(); iter != rectangles.end(); ++iter) {
			irectangle2& r = iter->first;
			int tmp = r.ll.x;
	        
			r.ll.x = width - 1 - r.ur.x;
			r.ur.x = width - 1 - tmp;
		}
	}

	void create_and_save_basic(const config_dictionary& cfg, const char* pattern)
	{
		typedef list<pair<irectangle2, string> > rectangle_list_t;

		string srcdir = cfg.get_value_string("src_dir", "");
		string outdir = cfg.get_value_string("out_dir", "");
		string outprefix = cfg.get_value_string("out_prefix", "");
		string maskext = cfg.get_value_string("mask_extension", "");
		string gtext = cfg.get_value_string("groundtruth_extension", ".groundtruth");
		bool flip = cfg.get_value_bool("flip", false);
		bool separate_colors = cfg.get_value_bool("separate_colors", false);
		int merge_scales = cfg.get_value_int("merge_scales", 1);
		int border_size = cfg.get_value_int("border_size", 100);

		if (!maskext.empty()) {
			cout << "Warning: running creation using mask_extension = '" << maskext << "'" << endl;
		}

		list<string> files;

		end_dir(srcdir);
		end_dir(outdir);
		if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, pattern, srcdir))
			list_directory(files, srcdir + "/" + string(pattern));

        cout << "Creating layer 1 from '" << srcdir << "'" << endl;
        cout << "  to '" << outdir << "'" << endl;
        cout << files.size() << " files found (from file = " << (cfg.get_value_bool("from_file", false) ? "true" : "false") <<
            ", pattern = " << pattern << ")" << endl;

		for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {

			vector<img> srcimg;
			img im;

			if (separate_colors && false) {
				cout << "using separate colors " << endl;
				img red, green, blue;
				img::read_colors(srcdir + *fiter, red, green, blue);
				
				srcimg.push_back(red);			
				srcimg.push_back(green);
				srcimg.push_back(blue);
			} else {
				img::read_colors(srcdir + *fiter, im);
				srcimg.push_back(im);
			}

			cout << "Processing " << *fiter;		
			if (srcimg.size() < 0 || srcimg.front().width == 0 || srcimg.front().height == 0)  {
				cout << " error reading file!" << endl;
			} else {
				img maskimg;

				try {
                    if (!maskext.empty())
				        img::read_grayscale(srcdir + change_extension(*fiter, maskext), maskimg);
                    // Size of mask image
				} catch (...) {
                    cout << "Mask extension specified but mask file can not be loaded." << endl;
				}

				int img_width = srcimg.front().width;
				int img_height = srcimg.front().height;

				// Create
				list<vector<img> > srcimglist;
                list<vector<img> > srcmasklist;
				rectangle_list_t rectangles;
				list<rectangle_list_t> rectlist;
				int splitn = cfg.get_value_int("split", 1);
				int splito = cfg.get_value_int("split_overlap", 100);

				cout << " (";
				read_groundtruth(rectangles, srcdir + *fiter, "", gtext, false);
				if (flip) {
					for (vector<img>::iterator img_it = srcimg.begin(); img_it != srcimg.end(); img_it++) {
						img_it->flip_horizontal();
					}
					flip_groundtruth(rectangles, (int)img_width);
				}

				split_image(srcimglist, rectlist, srcimg, rectangles, splitn, splito, 0.75);
                if (!maskimg.empty() && splitn > 1)
                    cout << "Splitting is not supported when processing with mask image." << endl;
					
				int mings = cfg.get_value_int("min_groundtruth_size", -1);
				int mingl = cfg.get_value_int("min_groundtruth_limit", mings);
				list<vector<img> >::iterator imgiter;
				list<rectangle_list_t>::iterator rectiter;
				string suffix = "_a";

				for (imgiter = srcimglist.begin(), rectiter = rectlist.begin(); 
						imgiter != srcimglist.end() && rectiter != rectlist.end();
						++imgiter, ++rectiter) {

					vector<img>& srcimg = *imgiter;
					rectangle_list_t& rectangles = *rectiter;
					vector<layer1_result*> result;

					if (mings <= 0)  {

						rectangle_list_t rect_list;
						read_groundtruth(rect_list, srcdir + *fiter, "", ".regions", false);

						vector<irectangle2> mask_regions;
						for (rectangle_list_t::iterator iter = rect_list.begin(); iter != rect_list.end(); iter++) {
							irectangle2 rect(iter->first.ll.x,iter->first.ll.y,iter->first.ur.x,iter->first.ur.y);
							rect.grow(1.1);
							mask_regions.push_back(rect);
						}
						clock_t start = clock();
						if (mask_regions.size() > 0) 
							lay1create::create_layer1(result, cfg, srcimg, mask_regions);
						else
							lay1create::create_layer1(result, cfg, srcimg, maskimg);
						clock_t end = clock();
						//cout << "time it took: " << (double)(end-start)/CLOCKS_PER_SEC << " sec" << endl;
					} else {
						irectangle2 minr = min_rectangle(rectangles);
						int mins = max(minr.x_dim(), minr.y_dim());
						irectangle2 maxr = max_rectangle(rectangles);
						int maxs = max(maxr.x_dim(), maxr.y_dim());
						config_dictionary cfg2(cfg);

						if (!maxr.invalid() && !minr.invalid() && mins > 0 && maxs > 0) {
							cfg2.set_value("init_size", (int)(-100 * mings/(double)mins));

							double f = cfg.get_value_double("scale_factor", 1/::pow(2.0, 1.0/3.0)) * mingl/(double)maxs;

							cfg2.set_value("scale_limit", (int)(max(img_width, img_height)*f*0.9));
						}
						lay1create::create_layer1(result, cfg2, srcimg, maskimg);
					}

					for (vector<layer1_result*>::iterator res_it = result.begin(); res_it != result.end(); res_it++) {
						for (int i = 1; i < merge_scales && i < result.size(); i++) {
							(*res_it)->merge(result[i], border_size);
						}
					}
					for (int i = result.size() - 1; i > result.size() - merge_scales; i--) 
						delete result[i];
					result.resize(result.size() - merge_scales + 1 );

					save_results(result, rectangles, *fiter, outdir, outprefix, (srcimglist.size() > 1) ? suffix : ""); 
					++suffix[1];
				}
				//} else {
				//	img maskimg;

				//	try {
				//		img::read_grayscale(srcdir + change_extension(*fiter, maskext), maskimg);
				//	} catch (...) {
				//		maskimg = img();
				//	}
				//	if (maskimg.empty()) 
				//		cout << " Warning: mask file does not exist (skipping creation";
				//	else { 
				//		vector<layer1_result*> result;

				//		cout << " (";
				//		lay1create::create_layer1(result, cfg, srcimg, maskimg);

				//		for (vector<layer1_result*>::iterator res_it = result.begin(); res_it != result.end(); res_it++) {
				//			for (int i = 1; i < merge_scales && i < result.size(); i++) {
				//				(*res_it)->merge(result[i], border_size);
				//			}
				//		}
				//		for (int i = result.size() - 1; i > result.size() - merge_scales; i--) 
				//			delete result[i];

				//		result.resize(result.size() - merge_scales + 1 );
				//		
				//		rectangle_list_t tmp;
				//		save_results(result, tmp, *fiter, outdir, outprefix, ""); 
				//	}
				//}

				// Save

				//for (int i = 0; i < (int)result.size(); ++i) {
				//    string istr = string("") + i;

				//    if (i > 0) cout << ' ';
				//    cout << istr;

				//    // Update metadata save "result"
				//    result[i]->original_width = srcimg.width;
				//    result[i]->border = border;
				//    save_layer1_result(result[i], change_extension(outdir + outprefix + *fiter, "_" + istr + ".ly1"));

				//    // Save (transformed) groundtruths
				//    if (!rectangles.empty()) {
				//        double factor = ((double)(result[i]->x_size(0) - 2*border))/srcimg.width;
				//        list<pair<irectangle2, string> > rectangles2(rectangles);
				//        ipoint2 bpt(border, border);

				//        for (list<pair<irectangle2, string> >::iterator riter = rectangles2.begin(); riter != rectangles2.end(); ++riter) {
				//            riter->first = (riter->first * factor) + bpt;
				//        }
				//        save_groundtruth(rectangles2, outdir, *fiter, outprefix, "_" + istr);
				//    }

				//    delete result[i];
				//}
				cout << ')' << endl;
				
			}
		}
	}

}

HOPINTERFACE_API bool hop_1_inference(const char* cfgfile, const char* pattern, const char* params)
{
	bool ok = true;
	try {
		config_dictionary cfg_all(cfgfile);
		cfg_all.from_string(params);
		cfg_all.update_namespace_references();

		config_dictionary cfg; 

		// copy values based on namespace hierarhy
		cfg.from_namespace_priority(cfg_all,2,"inference","ly1");

		if (cfg.get_value_bool("opencl_enable",false)) {
			vector<string> kernel_paths;
			vector<string> used_devices;

			cfg.get_value(kernel_paths, "opencl_kernel_paths");
			cfg.get_value(used_devices, "opencl_devices");
			
			OpenCL::initializeContexts(kernel_paths, used_devices);
		}
		string mode = cfg.get_value_string("mode", "normal");
		transform(mode.begin(), mode.end(), mode.begin(), TOUPPER);
		if (mode == "NORMAL") layer1::create_and_save_basic(cfg, pattern);
		else if (mode == "EXTRACTION") layer1::create_and_save_extraction(cfg, pattern);
		else
			throw custom_libhop_exception(config_exception, "Invalid mode");

	} catch (const libhop_exception& e) {
		cout << e.what() << endl;
		ok = false;
    } catch (const exception& e) {
		cout << "General exception '" << typeid(e).name() << "' with message: '" << e.what() << "'" << endl;
		ok = false;
	}
	return ok;
}


namespace layern {

	/// Inference of *one layer* on files described by pattern -- using config cfg
	void create_and_save(const config_dictionary& cfg, const char* pattern)
	{
		int layer = cfg.get_value_int("layer_index", -1, true);

		string srcdir = cfg.get_value_string("src_dir", "");
		string outdir = cfg.get_value_string("out_dir", "");
		string newext = cfg.get_value_string("result_extension", ".lyx");
		string libname;

		part_lib* library;
		list<string> files;
		clock_t time = 0;

		end_dir(srcdir);
		end_dir(outdir);

		
		cfg.get_value(libname, "part_lib_name", true);
		read_library(libname, library);

		if (library == nullptr) {
			throw custom_libhop_exception(io_exception, "Can not open library file \'" + libname + "\'.");
		}

		if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, pattern, srcdir, cfg.get_value_string("pattern", "")))
			list_directory(files, srcdir + pattern);
		for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
			layer1_result* res;

			cout << "Processing " << *fiter;
			read_layer1_result(res, srcdir + *fiter);
			if (res == nullptr) cout << " error reading file!" << endl;       
			else {
				vector<layer1_result*> resv;
				clock_t t;
				double quot;

				resv.push_back(res);
				t = layncreate::create_layern(resv, cfg, library);
				quot = resv.front()->cover_quotient(layer - 2);
				save_layer1_result(resv.front(), outdir + change_extension(*fiter, newext));
				delete res;

				time += t;

				cout << " done with " << quot;
				cout << " (in " << (double)t/CLOCKS_PER_SEC << " sec)" << endl;
			}
		}
		cout << "Total time: " << (double)time/CLOCKS_PER_SEC << " sec" << endl;
		cout << "Average time: " << (double)time/CLOCKS_PER_SEC/files.size() << " sec/file" << endl;
		delete library;
	}


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


	/*
		adjustment_layer += 1;
		if (adjustment_layer >= start_layer && adjustment_layer <= end_layer) {
			int i = start_layer; 
			while (i <= end_layer) {
				// bool mem = layern_creator::creators[i].add_activation_edges;
				//
				// layern_creator::creators[i].add_activation_edges = false;
				layern_creator::creators[i]->add_layer(res, i);
				//layern_creator::creators[i].add_activation_edges = mem;

				if (i != adjustment_layer) {
					++i;
				} else {
					iipair stat = layer_node_count(res, i - 1);
	                
					if (adjustment_average > 0 && (stat.second == 0 || (double)stat.first/stat.last < adjustment_average) ||
							stat.first < adjustment_count) {
						do {
							double& gthresh = layer_creator::creators[i].g_response_threshold;

							gthresh *= 0.5;      // or something else
							if (gthresh > 0.001) // ~0
								break;
							--i;
							res->delete_layers_geq(i - 1);
						} while (i >= start_layer);
					} else {
						++i;
					}
	                    
				}
			}
		}







	*/

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

	/// Inference of *multiple layers* on files described by pattern -- using config cfg
	void create_and_save_multiple(const config_dictionary& cfg, const string& fpatt)
	{
		int layer_index = cfg.get_value_int("layer_index", 0);
		int start_layer = -1, end_layer = -1;

		cfg.get_value(start_layer, "start_layer", layer_index > 0 ? false : true);
		cfg.get_value(end_layer, "end_layer", layer_index > 0 ? false : true);	

		// either ('layer_index') or ('start_layer' and 'end_layer') can be defined but not both
		if (layer_index > 0) {
			// if both start_layer and end_layer are in config and at the same time layer_index is defined, then return error
			if (start_layer > 0 || end_layer > 0) {
				throw custom_libhop_exception(config_exception, "Error: Ambiguous definition of layer. Found 'layer_index' and 'start_layer','end_layer' in config. Use either 'layer_index' or 'start_layer' and 'end_layer'.");
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

		bool save_result = cfg.get_value_bool("save_result", true);
		int time_count_split = cfg.get_value_int("time_count_split", end_layer);
		int maxcount = cfg.get_value_int("file_limit", INT_MAX);
		string globlibname = cfg.get_value_string("part_lib_name", "");
		//int layer_merge = cfg.get_value_int("layer_merge", 0) - 1;
		bool video_mode = cfg.get_value_bool("video_mode", false);
		double video_thresh = cfg.get_value_double("video_mode_tolerance", 0.2);
		int adjustment_layer = cfg.get_value_int("adjustment_layer", -1) +1;
        bool shape_check = false;

		string srcdir = "";
		string outdir = "";
		string newext = ".lyx";

		vector<int> endlayers(MAX_LAYER_NUMBER, 0);
		vector<pair<bool,string> > save_intermediate_results(MAX_LAYER_NUMBER, make_pair<bool,string>(0,""));

		map<string, part_lib*> libmap;
		map<string, part_lib*>::iterator libmapiter;
		part_lib* globlib;
        K_bin bin(12, 2, 4, 7); 
        bool merge_scales = false;

		read_library(globlibname, globlib);

		if (globlib != nullptr) 
			cout << "Warning: library \'" << globlibname << "\' will be used for all layers." << endl;

		vector<layern_creator*> layern_creator_list(end_layer+1, nullptr);

		for (int i = start_layer; i <= end_layer; ++i) {
			string ly_namespace = string("ly") + i;
			string key = string("end_layer") + i;

			// check if should save result for this intermediate layer (only for non-last layer)
			if (i != end_layer) {
				// use only values from lyN namespace e.g. ly3.save_result and ly3.result_extension for layer 3
				bool save_layer_i = cfg.get_value_bool(ly_namespace + ".save_result",false);
				string save_layer_i_extention = cfg.get_value_string(ly_namespace + ".result_extension", "", save_layer_i);
				
				save_intermediate_results[i].first = save_layer_i;
				save_intermediate_results[i].second = save_layer_i_extention;
			}
			endlayers[i] = cfg.get_value_int(key, i);

			config_dictionary cfgi;
			cfgi.from_namespace_priority(cfg, 1, ly_namespace.c_str());
			
			part_lib* library;

			// use global library or find one we have already loaded and saved to libmap
			if (globlib != nullptr) {
				library = globlib;
			} else {
				string libname = cfgi.get_value_string("part_lib_name", "");

				if (libname.empty()) {
					cout << "Error: Key \'part_lib_name\' in namespace \'" << ly_namespace << "\' expected but not found." << endl;
					return;
				}
				if ((libmapiter = libmap.find(libname)) != libmap.end()) 
					library = libmapiter->second;
				else {
					read_library(libname, library);
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
			if (i == start_layer && !cfg.is_defined("src_dir")) 
				cfgi.get_value(srcdir, "src_dir");
			else if (i == end_layer && !cfg.is_defined("out_dir"))
				cfgi.get_value(outdir, "out_dir");
	        
			cfgi.get_value(newext, "result_extension");
            merge_scales = merge_scales || (cfgi.get_value_int("scale_merge", 1) > 1);
		}		
		
		if (cfg.is_defined("src_dir")) cfg.get_value(srcdir, "src_dir");
		if (cfg.is_defined("out_dir")) cfg.get_value(outdir, "out_dir");
		if (cfg.is_defined("result_extension")) cfg.get_value(newext, "result_extension");
		
		end_dir(srcdir);
		end_dir(outdir);

		list<list<string> > files;
		clock_t total_time = 0, total_time1 = 0;
		clock_t time, time1, start, end;
		int count = 0;

		if (merge_scales) {
            cout << "Warning: scale merging is enabled for one of the layers, please make sure that " 
                    "from_file is set to true and pattern to an appropriate pattern." << endl;
			//file_lists_from_file(files, fpatt, srcdir, cfg.get_value_string("pattern", "_*.ly1"));
			if (!cfg.get_value_bool("from_file", false) || !file_lists_from_file(files, fpatt, srcdir, cfg.get_value_string("pattern", "%s_*.ly1"))) {
				list<string> fl;
				set<string> fls;

				list_directory(fl, srcdir + fpatt);
				for (list<string>::iterator iter = fl.begin(); iter != fl.end(); ++iter) {
					fls.insert(change_extension(*iter, "", "_"));
				}
				file_lists_from_list(files, list<string>(fls.begin(), fls.end()), srcdir, cfg.get_value_string("pattern", "%s_*.ly1"));
			}

		} else {
			list<string> fl;

			if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(fl, fpatt, srcdir, cfg.get_value_string("pattern", "")))
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
				double adjustment_average = cfg.get_value_double("adjustment_average", 0.0);
				int adjustment_count = cfg.get_value_int("adjustment_count", 1000);

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

}

HOPINTERFACE_API bool hop_n_inference(const char* cfgfile, const char* pattern, const char* params)
{
	bool ok = true;
	try {
		config_dictionary cfg_all(cfgfile);
		cfg_all.from_string(params);
		cfg_all.update_namespace_references();

		config_dictionary cfg;
		// copy values based on namespace hierarhy
		cfg.from_namespace_priority(cfg_all,1,"inference");

		if (cfg.get_value_bool("opencl_enable",false)) {
			vector<string> kernel_paths;
			vector<string> used_devices;

			cfg.get_value(kernel_paths, "opencl_kernel_paths");
			cfg.get_value(used_devices, "opencl_devices");

			OpenCL::initializeContexts(kernel_paths, used_devices);
		}

		layern::create_and_save_multiple(cfg, pattern);
	} catch (const libhop_exception& e) {
		cout << e.what() << endl;
		ok = false;
    } catch (const exception& e) {
		cout << "General exception '" << typeid(e).name() << "' with message: '" << e.what() << "'" << endl;
		ok = false;
	}
	return ok;
}