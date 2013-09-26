// OpenCV interface

#include "ocv.h"
#include <cv.h>
#include <cv.hpp>
#include <cxmisc.h>
#include <fstream>
#include "structures.h"
#include "exceptions.h"

// local functions
///////////////////////////////////////////////////////////////////////////////

string change_extension(const string& fname, const string& ext);
//{
//    string result;
//    string::size_type dot_pos;
//
//    if ((dot_pos = fname.rfind(".")) == string::npos) result = fname + ext;
//    else result = fname.substr(0, dot_pos) + ext;
//    return result;
//}


// svm_learning
///////////////////////////////////////////////////////////////////////////////


svm_learning::svm_learning(int fcount) : 
    feature_count(fcount), 
    positive(), negative()
{
}

int svm_learning::add_positive_sample(const map<int, double>& fmap)
{
    positive.push_back(vector_t(feature_count, 0.0));
    for (map<int, double>::const_iterator iter = fmap.begin(); iter != fmap.end(); ++iter) {
        if (iter->first >= feature_count) throw new_libhop_exception("Feature out of range.");
        positive.back()[iter->first] = (real_t)iter->second;
    }
}

int svm_learning::add_negative_sample(const map<int, double>& fmap)
{
    negative.push_back(vector_t(feature_count, 0.0));
    for (map<int, double>::const_iterator iter = fmap.begin(); iter != fmap.end(); ++iter) {
        if (iter->first >= feature_count) throw new_libhop_exception("Feature out of range.");
        negative.back()[iter->first] = (real_t)iter->second;
    }
}

void svm_learning::train(const string& fname, const irectangle2& box)
{
    normalize();

    // Fill train data
    cv::Mat trainData((int)positive.size() + (int)negative.size(), feature_count, CV_32FC1);
    int sample = 0;

    for (list_t::iterator iter = positive.begin(); iter != positive.end(); ++iter) {
        memcpy(trainData.ptr(sample++), &(iter->at(0)), sizeof(real_t) * feature_count);
    }
    for (list_t::iterator iter = negative.begin(); iter != negative.end(); ++iter) {
        memcpy(trainData.ptr(sample++), &(iter->at(0)), sizeof(real_t) * feature_count);
    }

    // Fill other data needed to tell the svm
    cv::Mat responses(trainData.rows, 1, CV_32SC1);
    cv::Mat samples(trainData.rows, 1, CV_8UC1);
    cv::Mat selectedFeatures(feature_count, 1, CV_8UC1);

    samples.setTo(cv::Scalar::all(1));
    selectedFeatures.setTo(cv::Scalar::all(1));
    responses.rowRange(0, (int)positive.size()).setTo(cv::Scalar::all(0));
    responses.rowRange((int)positive.size(), trainData.rows).setTo(cv::Scalar::all(1));

    // Run SVM   :pray
    cv::SVMParams params;
    cv::SVM svm;

    params.svm_type = cv::SVM::C_SVC;
    params.kernel_type = cv::SVM::LINEAR;

    svm.train(trainData, responses, selectedFeatures, samples, params);
    svm.save(fname.c_str());
    
    // Save normalization data: mean and var
    ofstream os(change_extension(fname, ".norm").c_str());
    os << feature_count << endl;
    os << box << endl;
    os << mean << endl;
    os << var << endl;
    os.close();
}

// Calculate mean and variance.
void svm_learning::calc_norm_data()
{
    mean.assign(feature_count, 0.0);
    var.assign(feature_count, 1.0);

    if (positive.empty() && negative.empty()) return;

    for (list_t::iterator iter = positive.begin(); iter != positive.end(); ++iter) 
        mean += *iter;
    for (list_t::iterator iter = negative.begin(); iter != negative.end(); ++iter) 
        mean += *iter;
    mean /= (real_t)(positive.size() + negative.size());

    vector_t tmp(feature_count);

    for (list_t::iterator iter = positive.begin(); iter != positive.end(); ++iter) {
        tmp = *iter - mean;
        sqr(tmp);
        var += tmp;
    }
    for (list_t::iterator iter = negative.begin(); iter != negative.end(); ++iter) {
        tmp = *iter - mean;
        sqr(tmp);
        var += tmp;
    }
    var /= (real_t)(positive.size() + negative.size() - 1);
    sqrt(var);
}

void svm_learning::normalize(vector_t& v)
{
    v -= mean;
    v /= var;
}

void svm_learning::normalize()
{
    calc_norm_data();

    for (list_t::iterator iter = positive.begin(); iter != positive.end(); ++iter)
        normalize(*iter);
    for (list_t::iterator iter = negative.begin(); iter != negative.end(); ++iter) 
        normalize(*iter);
}

// svm_predictor
///////////////////////////////////////////////////////////////////////////////


svm_predictor::svm_predictor() : 
    mean(), var(), avgbox(), predictor(nullptr)
{
}

svm_predictor::svm_predictor(const string& fname) : 
    mean(), var(), avgbox(), predictor(nullptr)
{
	init_predictor(fname);
}

void svm_predictor::init_predictor(const string& fname)
{
	if (predictor != nullptr) {
		delete (cv::SVM*)predictor;
		predictor = nullptr;
	}

    cv::SVM* svm = new cv::SVM();

    // Load svm
    svm->load(fname.c_str());

    // Load normalization data
    ifstream is(change_extension(fname, ".norm").c_str());

	if (is.fail()) {
        cout << "svm::predict: Can not open normalization file." << endl;
		return;
	}

    is >> feature_count;

    mean.resize(feature_count, 0.0);
    var.resize(feature_count, 1.0);

    is >> avgbox;
    is >> mean;
    is >> var;
    is.close();

    predictor = svm;
}

svm_predictor::~svm_predictor()
{
	if (predictor != nullptr)
		delete (cv::SVM*)predictor;
}

svm_predictor::real_t svm_predictor::predict(const map<int, double>& fmap)
{
	if (predictor == nullptr) 
		return -1;

    cv::SVM* svm = (cv::SVM*)predictor;
    cv::Mat test(1, feature_count, CV_32FC1);
    vector<real_t> v(feature_count, 0.0);

    // Classify :pray
    for (map<int, double>::const_iterator iter = fmap.begin(); iter != fmap.end(); ++iter) {
        if (iter->first >= feature_count) return -1.0;
        v[iter->first] = (real_t)iter->second;
    }

    v -= mean;
    v /= var;
    memcpy(test.ptr(0), &v.at(0), feature_count * sizeof(real_t));
    return svm->predict(test);
}

svm_predictor::real_t svm_predictor::predict(vector<real_t>& v)
{
	if (predictor == nullptr) 
		return -1;

    cv::SVM* svm = (cv::SVM*)predictor;
    cv::Mat test(1, feature_count, CV_32FC1);

    v -= mean;
    v /= var;
    memcpy(test.ptr(0), &v.at(0), feature_count * sizeof(real_t));
    return svm->predict(test);
}

// svm2
///////////////////////////////////////////////////////////////////////////////

svm2::svm2()
{
}

void svm2::train(const string& fname, const list_t& positive, const list_t& negative, CvFileStorage* storage)
{
    if (positive.empty() && negative.empty()) return;

    int feature_count;
    
    if (!positive.empty()) feature_count = (int)positive.front().size();
    else feature_count = (int)negative.front().size();

    cv::Mat trainData((int)positive.size() + (int)negative.size(), feature_count, CV_32FC1);
    int sample = 0;

    for (list_t::const_iterator iter = positive.begin(); iter != positive.end(); ++iter) {
        memcpy(trainData.ptr(sample++), &(iter->at(0)), sizeof(real_t) * feature_count);
    }
    for (list_t::const_iterator iter = negative.begin(); iter != negative.end(); ++iter) {
        memcpy(trainData.ptr(sample++), &(iter->at(0)), sizeof(real_t) * feature_count);
    }

    // Fill other data needed for svm
    cv::Mat responses(trainData.rows, 1, CV_32SC1);
    cv::Mat samples(trainData.rows, 1, CV_8UC1);
    cv::Mat selectedFeatures(feature_count, 1, CV_8UC1);

    samples.setTo(cv::Scalar::all(1));
    selectedFeatures.setTo(cv::Scalar::all(1));
    responses.rowRange(0, (int)positive.size()).setTo(cv::Scalar::all(0));
    responses.rowRange((int)positive.size(), trainData.rows).setTo(cv::Scalar::all(1));

    // Run SVM   :pray
    cv::SVMParams params;
    cv::SVM svm;

    params.svm_type = cv::SVM::C_SVC;
    params.kernel_type = cv::SVM::LINEAR;

    svm.train(trainData, responses, selectedFeatures, samples, params);
    //svm.save(fname.c_str());
	svm.write(storage, fname.c_str());

	
    //cv::SVM svm2;

    //svm2.load(fname.c_str());
    //cout << "--Positives: ";
    //for (list_t::const_iterator iter = positive.begin(); iter != positive.end(); ++iter) {
    //    const vector_t& v = *iter;
    //    cv::Mat test(1, feature_count, CV_32FC1);

    //    memcpy(test.ptr(0), &v.at(0), feature_count * sizeof(real_t));
    //    real_t x = svm2.predict(test);
    //    cout << x << " ";
    //}
    //cout << endl << "--Negatives: ";
    //for (list_t::const_iterator iter = negative.begin(); iter != negative.end(); ++iter) {
    //    const vector_t& v = *iter;
    //    cv::Mat test(1, feature_count, CV_32FC1);

    //    memcpy(test.ptr(0), &v.at(0), feature_count * sizeof(real_t));
    //    real_t x = svm2.predict(test);
    //    cout << x << " ";
    //}
    //cout << endl;
    
}

// svm2_predictor
///////////////////////////////////////////////////////////////////////////////

svm2_predictor::svm2_predictor() 
    : predictor(nullptr)
{
}

svm2_predictor::svm2_predictor(CvFileStorage* storage, const string& fname)
    : predictor(nullptr)
{
	init_predictor(storage, fname);
}

void svm2_predictor::invalidate()
{
	if (predictor != nullptr) {
		delete (cv::SVM*)predictor;
		predictor = nullptr;
	}
}

void svm2_predictor::init_predictor(CvFileStorage* storage, const string& fname)
{
	if (predictor != nullptr) {
		delete (cv::SVM*)predictor;
		predictor = nullptr;
	}

    cv::SVM* svm = new cv::SVM();

    // Load svm
    //svm->load(fname.c_str());
	svm->read(storage, cvGetFileNodeByName(storage, nullptr, fname.c_str()));
    predictor = svm;
}

svm2_predictor::~svm2_predictor()
{
	if (predictor != nullptr)
		delete (cv::SVM*)predictor;
}

svm2_predictor::real_t svm2_predictor::predict(vector<real_t>& v)
{
	if (predictor == nullptr) 
		return -1;

    int feature_count = (int)v.size();
    cv::SVM* svm = (cv::SVM*)predictor;
    cv::Mat test(1, feature_count, CV_32FC1);

    memcpy(test.ptr(0), &v.at(0), feature_count * sizeof(real_t));
    return svm->predict(test);
}


// Global functions
///////////////////////////////////////////////////////////////////////////////

vector<dpoint2> partition(const cv::Mat& m, int i /* = 0 */, bool row /* = true */)
{
    vector<dpoint2> result;

    if (row) {
        result.resize(m.cols/2);
        for (int j = 0; j < m.cols/2; ++j) {
            result[j].x = m.at<double>(i, 2*j);
            result[j].y = m.at<double>(i, 2*j + 1);
        }
    } else {
        result.resize(m.rows/2);
        for (int j = 0; j < m.rows/2; ++j) {
            result[j].x = m.at<double>(i, 2*j);
            result[j].y = m.at<double>(i, 2*j + 1);
        }
    }
    return result;
}

double TPS_energy(const vector<dpoint2>& pts1, const vector<dpoint2>& pts2)
{
    if (pts1.size() != pts2.size()) 
        return -1.0;

    int pcount = (int)pts1.size();

    cv::Mat L(pcount + 3, pcount + 3, CV_64F, cv::Scalar(0.0));
    
    // Fill matrix L
    for (int i = 0; i < pcount; ++i) {
        double* irow = L.ptr<double>(i);

        for (int j = 0; j < i; ++j) {
            double d2 = pts1[i].distance2(pts1[j]);

            if (d2 < 1e-6) L.at<double>(j, i) = irow[j] = 0.0;
            else L.at<double>(j, i) = irow[j] = d2*log(d2);
        }
    }
    
    for (int i = 0; i < pcount; ++i) {
        double* irow = L.ptr<double>(i) + pcount;

        *irow = 1.0;
        *(++irow) = pts1[i].x;
        *(++irow) = pts1[i].y;
    }

    double* row = L.ptr<double>(pcount);

    for (int j = 0; j < pcount; ++j)
        row[j] = 1.0;
    row = L.ptr<double>(pcount + 1);
    for (int j = 0; j < pcount; ++j)
        row[j] = pts1[j].x;
    row = L.ptr<double>(pcount + 2);
    for (int j = 0; j < pcount; ++j)
        row[j] = pts1[j].y;

    // Make matrix V
    cv::Mat V(pcount + 3, 2, CV_64F, cv::Scalar(0.0));

    for (int i = 0; i < pcount; ++i) {
        V.at<double>(i, 0) = pts2[i].x;
        V.at<double>(i, 1) = pts2[i].y;
    }

    cv::Mat c = L.inv(cv::DECOMP_SVD)*V;
    cv::Mat K(L, cv::Rect(0, 0, pcount, pcount));

    // Compute and return bending energy
    cv::Mat cc(c, cv::Rect(0, 0, 2, pcount));
    cv::Mat Q = cc.t()*K*cc;
    cv::Scalar result = cv::mean(Q.diag());

    if (abs(result[0]) < 1E-10)
        return 0.0;

    if (result[0] < 0) { 
        //cout << "Negative energy!" << endl;
        /*
        cout << result[0] << endl;
        for (int i = 0; i < pcount; ++i) {
            cout << pts1[i] << ' ';
        }
        cout << endl;
        for (int i = 0; i < pcount; ++i) {
            cout << pts2[i] << ' ';
        }
        cout << endl;
        for (int i = 0; i < pcount; ++i) {
            for (int j = 0; j < pcount; ++j) {
                cout << K.at<double>(i, j) << ' ';
            }
            cout << endl;
        }
        
        throw;
        */
    }

    return result[0];
}

void cv_write_to_stream(ostreamer& os, const cv::Mat& m)
{
    os.write(m.rows);
    os.write(m.cols);
    for (int i = 0; i < m.rows; ++i)
        for (int j = 0; j < m.cols; ++j) 
            os.write(m.at<double>(i, j));
}

void cv_read_from_stream(istreamer& is, cv::Mat& m)
{
    int rows, cols; 

    is.read(rows);
    is.read(cols);
    m.create(rows, cols, CV_64F);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            is.read(m.at<double>(i, j));
}

void cv_write_to_stream_f(ostreamer& os, const cv::Mat& m)
{
    os.write(m.rows);
    os.write(m.cols);
    for (int i = 0; i < m.rows; ++i)
        for (int j = 0; j < m.cols; ++j) 
            os.write(m.at<float>(i, j));
}

void cv_read_from_stream_f(istreamer& is, cv::Mat& m)
{
    int rows, cols; 

    is.read(rows);
    is.read(cols);
    m.create(rows, cols, CV_32F);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            is.read(m.at<float>(i, j));
}


