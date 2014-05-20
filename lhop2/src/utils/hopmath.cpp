// Collection of "math" functions and algorithms used in LHOP

#include <algorithm>

#include "hopmath.h"
#include <opencv2/highgui/highgui.hpp>
#include "utils/graphs/graph_utils.h"
#include "utils/utils.h"

// Splines
///////////////////////////////////////////////////////////////////////////////
double subspace_distance(const cv::Mat& data, const cv::Mat& mean, const cv::Mat& evec, const cv::Mat& eval, 
    double sigmamul)
{
    cv::Mat coeffs;
    cv::Mat result;

    cv::gemm(data - mean, evec, 1, cv::Mat(), 0, coeffs, cv::GEMM_2_T);
    for (int i = 0; i < coeffs.rows; ++i) {
        double& d = coeffs.at<double>(i, 0);

        d = min(d, sigmamul*eval.at<double>(i, 0));
    }
    cv::gemm(coeffs, evec, 1, mean, 1, result, 0);
    return cv::norm(result, data, cv::NORM_L2);
}



