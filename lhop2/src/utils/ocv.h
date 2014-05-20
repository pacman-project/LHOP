// OpenCV interface

#pragma once

#include <vector>
#include <list>
#include <map>
#include <string>
#include <algorithm>
#include "utils/structures.h"
#include "utils/serialization/streaming.h"

#if 0
#include <cv.h>
#include <ml.h>
#else
#include <opencv2/opencv.hpp>
#endif

// std::less specialization for cv::Point2 (required by the new code due for use in std::set!!)
namespace std
{
template <> 
struct less<cv::Point2i> {
	bool operator() (const cv::Point2i& a, const cv::Point2i& b) const {
		return a.y == b.y ? a.x < b.x : a.y < b.y ;
	}
};

// std::hash specialization for cv::Point2 (required by the new code due for use in std::unorderd_map!!)
template <>
struct hash<cv::Point2i> {
	std::hash<int> int_hash;
	size_t operator()(const cv::Point2i& obj) const {
		return int_hash(obj.x) ^ int_hash(obj.y);
	}
};
}

using namespace std;

// cv_serialization
/////////////////////

class cv_storable_data {
public:
    virtual void write_to_storage(CvFileStorage* fs) = 0;
    virtual void read_from_storage(CvFileStorage* fs) = 0;
};

void cv_write_to_stream(ostreamer& os, const cv::Mat& m);

void cv_read_from_stream(istreamer& is, cv::Mat& m);

void cv_write_to_stream_f(ostreamer& os, const cv::Mat& m);

void cv_read_from_stream_f(istreamer& is, cv::Mat& m);


// templates
///////////////////////////////////////////////////////////////////////////////

template<class T> void cv_write_csv(ostream& os, const cv::Mat& m)
{
    for (int r = 0; r < m.rows; ++r) {
        for (int c = 0; c < m.cols; ++c) {
            if (c > 0) os << ',';
            os << m.at<T>(r, c);
        }
        os << '\n';
    }
}

template<class T> cv::Mat flatten(const vector<point2<T> >& v)
{
    cv::Mat result(1, 2*(int)v.size(), CV_64F);

    for (int i = 0; i < (int)v.size(); ++i) {
        const point2<T>& p = v[i];

        result.at<double>(0, 2*i) = p.x;
        result.at<double>(0, 2*i + 1) = p.y;
    }
    return result;
}


///////////////////////////////////////////////////////////////////////////////

vector<dpoint2> partition(const cv::Mat& m, int i = 0, bool row = true);

double TPS_energy(const vector<dpoint2>& pts1, const vector<dpoint2>& pts2);
