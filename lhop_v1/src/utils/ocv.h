// OpenCV interface

#pragma once

#include <vector>
#include <list>
#include <map>
#include <string>
#include <algorithm>
#include "structures.h"
#include "streaming.h"

#include <cv.h>
#include <ml.h>

using namespace std;

// SVM (Support vector machine) interfaces
///////////////////////////////////////////////////////////////////////////////

// svm_learning class (used for learning)
/////////////////////////////////

class svm_learning {
protected:
    typedef float real_t;
    typedef vector<real_t> vector_t;
    typedef list<vector_t> list_t;

    int feature_count;

    list_t positive;
    list_t negative;

    vector_t mean;
    vector_t var;
public:
    svm_learning(int fcount);

    int add_positive_sample(const map<int, double>& fmap);
    int add_negative_sample(const map<int, double>& fmap);
    void train(const string& fname, const irectangle2& avgbox);

    static void predict(vector<int>& result, const string& fname, list_t& data);
protected:
    void calc_norm_data();
    void normalize(vector_t& v);
    void normalize();
};

// svm2
/////////

class svm2 {
public:
    typedef float real_t;
    typedef vector<real_t> vector_t;
    typedef list<vector_t> list_t;
protected:

public:
    svm2();

    void train(const string& fname, const list_t& positive, const list_t& negative, CvFileStorage* storage);

    //static void predict(vector<int>& result, const string& fname, list_t& data);
};

// svm_predictor
//////////////////

class svm_predictor {
protected:
    typedef float real_t;
    typedef vector<real_t> vector_t;
    typedef list<vector_t> list_t;

    int feature_count;
    irectangle2 avgbox;
    vector_t mean;
    vector_t var;

    void* predictor;
public:
	svm_predictor();
    svm_predictor(const string& fname);
    ~svm_predictor();

    irectangle2 get_average_box() { return avgbox; }
	bool valid() { return predictor != nullptr; }
	void init_predictor(const string& fname);
    real_t predict(const map<int, double>& fmap);
    real_t predict(vector<real_t>& fmap);
};

// svm2_predictor
///////////////////

class svm2_predictor {
    typedef float real_t;

    void* predictor;
public:
	svm2_predictor();
    svm2_predictor(CvFileStorage* storage, const string& fname);
    ~svm2_predictor();

	bool valid() { return predictor != nullptr; }
    void invalidate();
	void init_predictor(CvFileStorage* storage, const string& fname);
    real_t predict(vector<real_t>& fmap);
};

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
