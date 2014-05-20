// Graph utilities
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _GRAPH_UTILS_H_
#define _GRAPH_UTILS_H_

#include "utils/graphs/graph.h"
#include "utils/structures.h"
#include "utils/img.h"


// public functions
///////////////////////////////////////////////////////////////////////////////

graph* point_graph(const vector<ipoint2>& pts);

vector<node*> dag_longest_path(graph* g, node* start);

vector<node*> tree_longest_path(graph* g);

vector<node*> dfs_closest_path(graph* g, const set<node*>& nset, node* start);

vector<int> point_matching_p(const vector<ipoint2>& pts1, const vector<ipoint2>& pts2);

vector<int> piecewise_point_matching_p(const vector<pair<int, ipoint2> >& pts1, const vector<pair<int, ipoint2> >& pts2);

vector<double> curve_gaps(graph* g, const vector<ipoint2>& pts);

bool continuous_point_set(const vector<ipoint2>& ipts, double factor = 3.0);

#endif