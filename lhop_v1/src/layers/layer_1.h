/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// layer_1

#pragma once
#ifndef _LAYER_1_
#define _LAYER_1_

#include "../utils/img.h"
#include "../graphs/img_graph.h"
#include "layer_1_result.h"
#include "layer_1_creators.h"


// layer_1 utility functions
//////////////////////////////

void get_region_set(int x, int y, int width, int height,
                    const matrix<bool>& r, int cx, int cy, set<int>& result);

img region_to_img(const matrix<bool>& m);

template<class I> irectangle2 node_set_bounding_rectangle(I begin, I end)
{
    irectangle2 rect;

    for (; begin != end; ++begin) {
        layer1_data* nd = (layer1_data*)(*begin)->data;
        rect.eat(nd->x, nd->y);
    }
    return rect;
}

template<class I> node* increment_node_iter(node* n, I& iter, I end)
{
    if (iter == end || *iter == nullptr) return nullptr;

    node* result = ((layer1_data*)n->data)->next;

    if (result != nullptr) return result;
    ++iter;
    if (iter != end) return *iter;
    return nullptr;
}

template<class I> ipoint2 node_set_center(I begin, I end) 
{
    irectangle2 rect = node_set_bounding_rectangle(begin, end);
    if (rect.invalid()) return ipoint2::zero; else return rect.center();
}

/// fname must have an extension; it is automatically changed to .grundtruth
void read_groundtruth(list<irectangle2>& rectangles, const string& fname, 
					  const string& cat_name_only = "", const string& gt_extension = ".groundtruth");

void read_groundtruth(list<std::pair<irectangle2,string> >& rectangles, const string& fname,
					  const string& cat_name_only = "", const string& gt_extension = ".groundtruth",
                      bool display_warnings = true);

void save_groundtruth(const list<std::pair<irectangle2,string> >& rectangles, const string& outdir, const string& srcname, 
    const string& prefix, const string& suffix, bool changeext = true);

void get_feature_map_radial(map<int, double>& result, vector<node*>& nodes, const ipoint2& center,
    const vector<int>& angles, const vector<int>& radii, int mcount);

void get_feature_map_radial(vector<float>& result, vector<node*>& nodes, const ipoint2& center,
    const vector<int>& amap, int acount, const vector<int>& rmap, int rcount, int mcount);

void get_feature_maps(vector<int>& amap, int& acount, vector<int>& rmap, int& rcount,
    const vector<int>& angles, const vector<int>& radii);


#endif /* _LAYER_1_ */
