// Collection of "math" functions and algorithms used in LHOP

#pragma once

#include "utils/structures.h"
#include <opencv2/opencv.hpp>

// misc math functions
///////////////////////////////////////////////////////////////////////////////

template<class T> inline T sqr(const T& v) 
{ 
    return v*v;
}

double subspace_distance(const cv::Mat& data, const cv::Mat& mean, const cv::Mat& evec, const cv::Mat& eval, 
    double sigmamul);
