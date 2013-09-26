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
#include <cv.h>

#include "../utils/utils.h"
#include "part_lib.h"
#include "layer_1.h"


//#define DEBUG_OUTPUT

using namespace std;

// global functions
///////////////////////////////////////////////////////////////////////////////

// Fills edge_data_2::var field on 'lyrSrc' edges with the type of the dest 
// node and edge_data_2var::var field on 'lyrCenter' edges. If edge data is
// null on 'lyrCenter' edges it creates new object.
/*void update_var_fields(part_lib* library)
{
    int srcname = atom("lyrSrc");
    int centername = atom("lyrCenter");

    for (list<node*>::iterator niter = library->nodes.begin(); niter != library->nodes.end(); ++niter) {
        node* n = *niter;

        foreach_neighbor(n, srcname, iter) {
            part_data* pd = (part_data*)neighbor_node_data(iter);
            part_data_2* ed = dynamic_cast<part_data_2*>(neighbor_edge_data(iter));

            if (ed != nullptr) ed->app.insert(pair<int, double>(pd->type, 1.0));
        }

        // there is only one center, though, but it is convenient to use defines
        foreach_neighbor(n, centername, iter) {
            part_data* pd = (part_data*)neighbor_node_data(iter);
            part_data_2app* ed = dynamic_cast<part_data_2app*>(neighbor_edge_data(iter));

            if (ed != nullptr) ed->app.insert(pair<int, double>(pd->type, 1.0));
            else neighbor_edge_data(iter) = new part_data_2app(pd->type);
        }
    }
}*/

void read_library(const string& name, part_lib*& library)
{
    try {
		// !!!!! WARNING !!!!!!
		// anything added here MUST also be added to hop_library(hop_blob& blob) constructor in hop.cpp
		// !!!!! WARNING !!!!!!
        library = (part_lib*)streamable::read(name);
        library->update_var_fields();

        cv::FileStorage fs(name + ".yml", cv::FileStorage::READ);
        if (fs.isOpened()) {
            library->read_cv_storable_data(*fs);
            fs.release();
        }
    } catch (...) {
        library = nullptr;
    }
}

part_lib* read_library(const string& name)
{
    part_lib* result;

    read_library(name, result);
    return result;
}


// make set of integers 
//  (x, y) - center
//  width, height - dimensions of the picture
//  r - region
//  (cx, cy) - center of the region
//  result - note that it is NOT deleted, elements are added to this set!!
void get_region_set(int x, int y, int width, int height, const matrix<bool>& r, int cx, int cy, 
                    set<int>& result)
{
    int rwidth = (int)r.width, rheight = (int)r.height;
    int x_start, y_start, x_end, y_end;
    int ix_start, iy_start;

    x_start = max(x - cx, 0); 
    x_end = min(x + rwidth - cx, width);
    y_start = max(y - cy, 0); 
    y_end = min(y + rheight - cy, height);
    ix_start = x_start - (x - cx);
    iy_start = y_start - (y - cy);

    //result.clear();
    for (int i = ix_start, x = x_start; x < x_end; ++i, ++x) {
        for (int j = iy_start, y = y_start; y < y_end; ++j, ++y) {
            if (r(i, j)) result.insert(x + y*width);
        }
    }
}

img region_to_img(const matrix<bool>& m)
{
    img result((int)m.width, (int)m.height, 0.0);

    for_each_xy(m, i, j) {
        result(i, j) = m(i, j);
    }
    return result;
}

void read_groundtruth(list<irectangle2>& rectangles, const string& fname, 
					  const string& cat_name_only, const string& gt_extension)
{
	list<std::pair<irectangle2,string> > rectangles_with_cat;
	// read with category names
	read_groundtruth(rectangles_with_cat, fname, cat_name_only, gt_extension);
	// copy only rectangles and delete strings
	for (list<std::pair<irectangle2,string> >::iterator iter = rectangles_with_cat.begin(); 
		iter != rectangles_with_cat.end(); iter++) {
		std::pair<irectangle2,string> p = *iter;
		rectangles.push_back(p.first);
	}
	
}

void read_groundtruth(list<std::pair<irectangle2,string> >& rectangles, const string& fname, 
					  const string& cat_name_only, const string& gt_extension, bool display_warnings /* = true */)
//void read_groundtruth(list<irectangle2>& rectangles, const string& fname)
{
	//const string gt_extension = ".groundtruth";
	const string file = change_extension(fname, gt_extension);
    ifstream is(file.c_str());

	rectangles.clear();
	bool warning_displayed = !display_warnings;
    if (is.is_open()) {
        while (true) {
			string line;
			// read and parse whole line as stream of string
			getline(is,line);
			stringstream sline(line);

			rectangle2<double> rect;
			string cat_name;

			// first read rectangle values
			sline >> rect;
            if (is.fail()) break;
			
			// then read category name and trim it
			sline >> cat_name;
			trim(cat_name);
			
			if (cat_name.length() <= 0) {
				// category name is missing in groundtruth file
				// this may be old version so only ignore it and use all values
				if (warning_displayed == false) {
					// diplay warning only once for one file
					cout << "\nWARNING: category name is missing in groundtruth data (file: '" << file << "') - will use all bounding boxes in file" << endl;
					warning_displayed = true;
				}
			} else if (cat_name_only.length() > 0 && cat_name != cat_name_only) {
				// skip all ground truth values that do not belog to category cat_name_only
				continue;
			}
			rectangles.push_back(std::pair<irectangle2,string>(irectangle2((int)rect.ll.x, (int)rect.ll.y, (int)rect.ur.x, (int)rect.ur.y), cat_name));
        }
	} else {
		cout << "\nWARNING: Unable to open groundtruth file '" << file << "'" << endl;
	}
}

void save_groundtruth(const list<std::pair<irectangle2,string> >& rectangles, const string& outdir, const string& srcname, 
    const string& prefix, const string& suffix, bool changeext /* = true */)
{
    string fstr; 

    if (changeext) fstr = outdir + prefix + change_extension(srcname, suffix + ".groundtruth");
    else fstr = outdir + prefix + srcname + suffix;

    ofstream os(fstr.c_str());

    for (list<std::pair<irectangle2,string> >::const_iterator iter = rectangles.begin(); iter != rectangles.end(); ++iter) {
		const std::pair<irectangle2,string> p = *iter;
        const irectangle2& r = p.first;
		const string& cat_name = p.second;
        
        os << r.ll.x << ' ' << r.ll.y << ' ' << r.ur.x << ' ' << r.ur.y << ' ' << cat_name << endl;
    }
    os.close();
}

//  Get features from the nodes in set 'nodes'. Result maps: (angle sector, radial sector, type) |--> value
//  (Angle) sectors are defined as follows: if 'angles' = (a0, a1, a2, ..., an-1) then
//  sector si = {a : ai-1 <= a < ai}, i = 0, ..., n; and a-1 = 0, an = 360; angles are in degrees!
//  (Radial) sectors are defined as follows: if 'radii' = (r0, r1, r2, ..., rm-1) then
//  rsector ri = {r : ri-1 <= r < ri+1}, i = 0, ..., m-1.
void get_feature_map_radial(map<vector<int>, double>& result, vector<node*>& nodes, const ipoint2& center, 
    const vector<int>& angles, const vector<int>& radii)
{
    int maxr = 0;

    for (vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;
        int d = center.distance2(ipoint2(nd->x, nd->y));
        maxr = max<int>(maxr, (int)::sqrt((double)d));
    }

    vector<int> amap(360, (int)angles.size());
    vector<int> rmap(maxr + 1, (int)radii.size());
    int i, j;

    i = j = 0;
    while (i < (int)angles.size() && j < (int)amap.size()) {
        amap[j++] = i;
        if (j >= angles[i]) ++i;
    }
    i = j = 0;
    while (i < (int)radii.size() && j < (int)rmap.size()) {
        rmap[j++] = i;
        if (j >= radii[i]) ++i;
    }

    for (vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;

        while (n != nullptr) {
            layer1_data* nd = (layer1_data*)n->data;
            ipoint2 p(nd->x, nd->y);
            double angle = 180.0*atan2((double)(p.y - center.y), (double)(p.x - center.x))/M_PI;
            vector<int> tuple(3);

            if (angle < 0) angle += 360.0;
            tuple[0] = amap[(int)angle];
            int dbg = center.distance2(p);
            double dbg2 = sqrt((double)dbg);
            tuple[1] = rmap[(int)dbg2];
            tuple[2] = nd->m;

            map<vector<int>, double>::iterator miter = result.find(tuple);

            if (miter != result.end()) miter->second += nd->val();
            else result.insert(pair<vector<int>, double>(tuple, nd->val()));

            n = nd->next;
        }
    }

}

void get_feature_map_radial(vector<float>& result, vector<node*>& nodes, const ipoint2& center, 
    const vector<int>& amap, int acount, const vector<int>& rmap, int rcount, int mcount)
{
    result.assign(acount*rcount*mcount, 0);
    for (vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;

        while (n != nullptr) {
            layer1_data* nd = (layer1_data*)n->data;
            ipoint2 p(nd->x, nd->y);
            double angle = 180.0*atan2((double)(p.y - center.y), (double)(p.x - center.x))/M_PI;
            int a, r, m;

            if (angle < 0) angle += 360.0;
            a = amap[(int)angle];
            int dbg = center.distance2(p);
            double dbg2 = sqrt((double)dbg);
            r = rmap[(int)dbg2];
            m = nd->m;

            result[a*rcount*mcount + r*mcount + m] += (float)nd->val();

            n = nd->next;
        }
    }
}

void get_feature_map_radial(map<int, double>& result, vector<node*>& nodes, const ipoint2& center, 
    const vector<int>& angles, const vector<int>& radii, int mcount)
{
    map<vector<int>, double> m;

    get_feature_map_radial(m, nodes, center, angles, radii);
    for (map<vector<int>, double>::iterator iter = m.begin(); iter != m.end(); ++iter) {
        int index = iter->first[0]*((int)radii.size() + 1)*mcount + iter->first[1]*mcount + iter->first[2];
        result.insert(pair<int, double>(index, iter->second));
    }
}


void get_feature_maps(vector<int>& amap, int& acount, vector<int>& rmap, int& rcount, 
    const vector<int>& angles, const vector<int>& radii)
{
    const int maxr = 500;

    amap.assign(360, (int)angles.size());
    rmap.assign(maxr, (int)radii.size());
    acount = (int)angles.size() + 1;
    rcount = (int)radii.size() + 1;

    int i, j;

    i = j = 0;
    while (i < (int)angles.size() && j < (int)amap.size()) {
        amap[j++] = i;
        if (j >= angles[i]) ++i;
    }
    i = j = 0;
    while (i < (int)radii.size() && j < (int)rmap.size()) {
        rmap[j++] = i;
        if (j >= radii[i]) ++i;
    }
}

