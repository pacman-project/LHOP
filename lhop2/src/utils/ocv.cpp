// OpenCV interface

#include "ocv.h"
//#include <cv.h>
//#include <cv.hpp>
//#include <cxmisc.h>
#include <fstream>
#include "utils/structures.h"
#include "utils/exceptions.h"

// local functions
///////////////////////////////////////////////////////////////////////////////

string change_extension(const string& fname, const string& ext);


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


