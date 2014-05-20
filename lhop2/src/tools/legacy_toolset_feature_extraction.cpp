// feature_extraction.cpp : Defines the entry point feature extraction functionality
//

#include "tools/main_toolset.h"

#include <stdio.h>
#include <string>
#include <cctype>

#include "core/legacy/layers.h"

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
void create_and_save_extraction(const ConfigDictionary& cfg, const char* pattern)
{
	string srcdir = cfg.getValueString("src_dir", "");
	string outdir = cfg.getValueString("out_dir", "");
	string outprefix = cfg.getValueString("out_prefix", "");
	string gtext = cfg.getValueString("groundtruth_extension", ".groundtruth");
	string catname = cfg.getValueString("category_name", "");
	int gtdelta = cfg.getValueInt("groundtruth_delta", 0);
	int border = cfg.getValueInt("border_size", 0);
	bool separate_colors = cfg.getValueBool("separate_colors", false);
    string maskext = cfg.getValueString("mask_extension", "");
    int maskthick = cfg.getValueInt("mask_thickening", 0);
	int maxdim;

	cfg.getValue(maxdim, "max_max_dimension", true); 
	    
	int mindim = cfg.getValueInt("min_max_dimension", maxdim);
	   
	list<string> files;

	end_dir(srcdir);
	end_dir(outdir);
	if (!cfg.getValueBool("from_file", false) || !list_from_file(files, pattern, srcdir))
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
                

			ConfigDictionary cfg2(cfg); // just change the dimensions...
			vector<layer1_result*> result;
			string maxdimstr = string("") + maxdim;
			string mindimstr = string("") + mindim;

			cfg2.setValue("init_size", maxdimstr);
			cfg2.setValue("scale_limit", mindimstr);

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

void create_and_save_basic(const ConfigDictionary& cfg, const char* pattern)
{
	typedef list<pair<irectangle2, string> > rectangle_list_t;

	string srcdir = cfg.getValueString("src_dir", "");
	string outdir = cfg.getValueString("out_dir", "");
	string outprefix = cfg.getValueString("out_prefix", "");
	string maskext = cfg.getValueString("mask_extension", "");
	string gtext = cfg.getValueString("groundtruth_extension", ".groundtruth");
	bool flip = cfg.getValueBool("flip", false);
	bool separate_colors = cfg.getValueBool("separate_colors", false);
	int merge_scales = cfg.getValueInt("merge_scales", 1);
	int border_size = cfg.getValueInt("border_size", 100);

	if (!maskext.empty()) {
		cout << "Warning: running creation using mask_extension = '" << maskext << "'" << endl;
	}

	list<string> files;

	end_dir(srcdir);
	end_dir(outdir);
	if (!cfg.getValueBool("from_file", false) || !list_from_file(files, pattern, srcdir))
		list_directory(files, srcdir + "/" + string(pattern));

    cout << "Creating layer 1 from '" << srcdir << "'" << endl;
    cout << "  to '" << outdir << "'" << endl;
    cout << files.size() << " files found (from file = " << (cfg.getValueBool("from_file", false) ? "true" : "false") <<
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
			int splitn = cfg.getValueInt("split", 1);
			int splito = cfg.getValueInt("split_overlap", 100);

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
					
			int mings = cfg.getValueInt("min_groundtruth_size", -1);
			int mingl = cfg.getValueInt("min_groundtruth_limit", mings);
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
					ConfigDictionary cfg2(cfg);

					if (!maxr.invalid() && !minr.invalid() && mins > 0 && maxs > 0) {
						cfg2.setValue("init_size", (int)(-100 * mings/(double)mins));

						double f = cfg.getValueDouble("scale_factor", 1/::pow(2.0, 1.0/3.0)) * mingl/(double)maxs;

						cfg2.setValue("scale_limit", (int)(max(img_width, img_height)*f*0.9));
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

string LegacyFeatureExtractionToolset::getShortDescription() {
	return "Legacy version of feature extraction (first layer inference) - using only layer1_creators code";
}
string LegacyFeatureExtractionToolset::getLongDescription() {
	return "";
}
string LegacyFeatureExtractionToolset::getUsageDescription() {
	return "";
}

bool LegacyFeatureExtractionToolset::areArgumentsValid(int argc, char* argv[]) {
	// TODO: do more inteligent verification
	return argc > 2  ? true : false;
}

void LegacyFeatureExtractionToolset::main(int argc, char* argv[]) {

	// parse input arguments
	const char* cfgfile = argv[2];
	const char* patt = argc > 3 ? argv[3] : "";
	const char* params = argc > 4 ? argv[4] : "";	

	ConfigDictionary cfg_all(cfgfile);
	cfg_all.fromString(params);
	cfg_all.updateNamespaceReferences();

	ConfigDictionary cfg; 

	// copy values based on namespace hierarhy
	cfg.fromNamespacePriority(cfg_all,2,"inference","ly1");

	string mode = cfg.getValueString("mode", "normal");
	transform(mode.begin(), mode.end(), mode.begin(), TOUPPER);
	if (mode == "NORMAL") create_and_save_basic(cfg, patt);
	else if (mode == "EXTRACTION") create_and_save_extraction(cfg, patt);
	else
		throw custom_libhop_exception(ConfigException, "Invalid mode");
}