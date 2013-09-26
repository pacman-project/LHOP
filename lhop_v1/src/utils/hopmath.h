// Collection of "math" functions and algorithms used in LHOP

#pragma once

#include "structures.h"
#include <cv.h>

// misc math functions
///////////////////////////////////////////////////////////////////////////////

template<class T> inline T sqr(const T& v) 
{ 
    return v*v;
}


// Splines
///////////////////////////////////////////////////////////////////////////////

struct cubic_spline_data {
    double a0, a1, a2, a3;  // coefficients of cubic polynomial for x 
    double b0, b1, b2, b3;  // and y
    dpoint2 d0, d1; // tangents at t = 0 and t = 1
};

typedef vector<cubic_spline_data> cubic_spline;

void cubic_g1_spline(vector<cubic_spline_data>& spline, const vector<dpoint2>& T, double lambda = 0);

void Bspline(vector<cubic_spline_data>& spline, const vector<dpoint2>& T);

void Bspline_tangents(vector<dpoint2>& tangents, const vector<dpoint2>& t);

void spline_points(vector<dpoint2>& pts, const vector<cubic_spline_data>& spdv, int steps = 10);

void spline_points(vector<dpoint2>& pts, const vector<cubic_spline>& spdv, int steps = 10);

cv::Mat1d spline_image(const cubic_spline& splines, int thickness = 1, int imgborder = 50, double factor = 1.0);

cv::Mat1d point_image(const vector<dpoint2>& pts, int radius = 1, int imgborder = 50);

cv::Mat1d spline_image(const vector<cubic_spline>& splines, int thickness = 1, int imgborder = 50, double factor = 1.0);

double spline_length(cubic_spline& spline);

void keep_longest_splines(vector<cubic_spline>& splines, double thresh);

void save_spline(const string& s, const vector<cubic_spline_data>& spline, int thickness = 1);

void splines_from_points(vector<cubic_spline>& splines, vector<ipoint2> pts, int inhibition);

void tangents_from_points(vector<dpoint2>& tangents, const vector<ipoint2>& pts);

double subspace_distance(const cv::Mat& data, const cv::Mat& mean, const cv::Mat& evec, const cv::Mat& eval, 
    double sigmamul);
