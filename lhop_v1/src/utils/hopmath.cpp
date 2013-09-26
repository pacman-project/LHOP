// Collection of "math" functions and algorithms used in LHOP

#include <algorithm>

#include "hopmath.h"
#include <highgui.h>
#include "graphs/graph_utils.h"
#include "utils/utils.h"

// Splines
///////////////////////////////////////////////////////////////////////////////

template<class T> inline T Delta(const vector<T>& v, int i)
{
    return v[i + 1] - v[i];
}

inline double angle(const dpoint2& a, const dpoint2& b)
{
    return acos(dot_product(a, b)/(norm(a)*norm(b)));
}

inline double sign(double x)
{
    return x < 0 ? -1 : (x > 0 ? 1 : 0);
}

void g1_tangents(vector<dpoint2>& d, vector<double>& a0, vector<double>& a1, const vector<dpoint2>& T, double lambda)
{
    int N = (int)T.size();
    int n = N - 1;

    if (N <= 1) return;

    vector<double> t(N);
    double x;
    dpoint2 X;

    // Init parametrization 't'
    t[0] = 0;
    for (int i = 1; i < N; ++i) {
        double dist = sqrt(norm(T[i] - T[i - 1]));
        t[i] = dist + t[i - 1];
    }
    x = t[n];
    for (int i = 0; i < N; ++i)
        t[i] /= x;

    // Calculate tangents
    d.resize(N);
    
    X = Delta(T, 0);
    d[0] = X/norm(X);
    for (int i = 1; i < n; ++i) {
        dpoint2 DTi = Delta(T, i);
        dpoint2 DTim1 = Delta(T, i - 1);
        double normDTi = norm(DTi);
        double normDTim1 = norm(DTim1);

        if (angle(DTim1, DTi) < 1E-6) 
            d[i] = DTi/normDTi;
        else {
            double z = sign(planar_cross_product(DTim1, DTi));
            dpoint2 u = DTim1.rotated90()*z;
            dpoint2 v = DTi.rotated90()*(-z);
            double dotuv = dot_product(u, v);
            double l;

            if (lambda > 0 && lambda < 1) 
                l = lambda;
            else if (dotuv == 0) 
                l = 0.5;
            else {
                double Dti = pow(Delta(t, i), 3);
                double Dtim1 = pow(Delta(t, i - 1), 3);
                double A = Dtim1*sqr(normDTi) + 2*Dti*dotuv - Dti*sqr(normDTim1);
                double B = sqr(Dti*sqr(normDTim1) - Dtim1*sqr(normDTi)) + 4*Dtim1*Dti*sqr(dotuv);
                
                l = (2*Dti/(A + sqrt(B)))*dotuv;
                if (l < 0 || l > 1) l = (2*Dti/(A - sqrt(B)))*dotuv;
            }
            
            dpoint2 w = u*l + v*(1 - l);
            d[i] = w/norm(w);
        }
    }
    X = Delta(T, n - 1);
    d[n] = X/norm(X);

    // Calculate "alphas"
    a0.resize(N);
    a1.resize(N);

    for (int i = 0; i < n; ++i) {
        a0[i] = dot_product(d[i], Delta(T, i));
        a1[i] = dot_product(d[i + 1], Delta(T, i));
    }
    a0[n] = a1[n] = 0.0;
}

void g1_spline_coefficients(cubic_spline_data& csd, const vector<dpoint2>& T, const vector<dpoint2>& d, 
    const vector<double>& a0, const vector<double>& a1, int i)
{
    csd.a0 = T[i].x; 
    csd.a1 = a0[i]*d[i].x;
    csd.a2 = -a1[i]*d[i + 1].x - 2*a0[i]*d[i].x + 3*T[i + 1].x - 3*T[i].x;
    csd.a3 = a1[i]*d[i + 1].x + a0[i]*d[i].x - 2*T[i + 1].x + 2*T[i].x;

    csd.b0 = T[i].y; 
    csd.b1 = a0[i]*d[i].y;
    csd.b2 = -a1[i]*d[i + 1].y - 2*a0[i]*d[i].y + 3*T[i + 1].y - 3*T[i].y;
    csd.b3 = a1[i]*d[i + 1].y + a0[i]*d[i].y - 2*T[i + 1].y + 2*T[i].y;

    csd.d0 = d[i];
    csd.d1 = d[i + 1];
}

void cubic_g1_spline(vector<cubic_spline_data>& spline, const vector<dpoint2>& T, double lambda /* = 0 */)
{
    vector<dpoint2> d;
    vector<double> a0; 
    vector<double> a1;
    int n = (int)T.size() - 1;

    g1_tangents(d, a0, a1, T, lambda);
    spline.resize(n);
    for (int i = 0; i < n; ++i) {
        g1_spline_coefficients(spline[i], T, d, a0, a1, i);
    }
}

void Bspline_coefficients(double& a, double& b, double& c, double&d, 
    double pim1, double pi, double pip1, double pip2, double f)
{
    a = (4*pi + pim1 + pip1)*f;
    b = 3*(pip1 - pim1)*f;
    c = 3*(pim1 + pip1 - 2*pi)*f;
    d = (3*pi - pim1 - 3*pip1 + pip2)*f;
}

void Bspline_coefficients(cubic_spline_data& spd, 
    const dpoint2& pim1, const dpoint2& pi, const dpoint2& pip1, const dpoint2& pip2, double f)
{
    Bspline_coefficients(spd.a0, spd.a1, spd.a2, spd.a3, pim1.x, pi.x, pip1.x, pip2.x, f);
    Bspline_coefficients(spd.b0, spd.b1, spd.b2, spd.b3, pim1.y, pi.y, pip1.y, pip2.y, f);
}

void Bspline_coefficients_d(double& a, double& b, double& c, 
    double pim1, double pi, double pip1, double pip2, double f)
{
    //{{-1, 3, -3, 1}, {3, -6, 3, 0}, {-3, 0, 3, 0}, {1, 4, 1, 0}};
    //{{-3, 9, -9, 3}, {6, -12, 6, 0}, {-3, 0, 3, 0}};
    //{3 pi - pim1 - 3 pip1 + pip2, -6 pi + 3 pim1 + 3 pip1, -3 pim1 + 3 pip1, 4 pi + pim1 + pip1}
    //{3 (3 pi - pim1 - 3 pip1 + pip2), 6 (-2 pi + pim1 + pip1), -3 pim1 + 3 pip1}
    a = 3*(pip1 - pim1)*f;
    b = 6*(pim1 + pip1 - 2*pi)*f;
    c = 3*(3*pi - pim1 - 3*pip1 + pip2)*f;
}

void Bspline_coefficients_d(cubic_spline_data& spd, 
    const dpoint2& pim1, const dpoint2& pi, const dpoint2& pip1, const dpoint2& pip2, double f)
{
    Bspline_coefficients_d(spd.a0, spd.a1, spd.a2, pim1.x, pi.x, pip1.x, pip2.x, f);
    Bspline_coefficients_d(spd.b0, spd.b1, spd.b2, pim1.y, pi.y, pip1.y, pip2.y, f);
    spd.a3 = 0;
    spd.b3 = 0;
}

// Returns point on line through 'A' and 'B'; if 't' = 0 it returns 
// 'A', if 't' = 1, it returns 'B'; If 'A' == 'B' it returns 'A' (= 'B')
// regardeless of 't'
dpoint2 line_point(const dpoint2& A, const dpoint2& B, double t)
{
    return B*t + A*(1 - t);
}

void Bspline(vector<cubic_spline_data>& spline, const vector<dpoint2>& T)
{
    if (T.empty())
        return;

    vector<dpoint2> newT(T.size() + 4);
    int n = T.size() - 1;

    newT[0] = newT[1] = T[0];
    newT[n + 4] = newT[n + 3] = T[n];
    for (int i = 0; i <= n; ++i) newT[i + 2] = T[i];
    
    spline.resize(n + 2);
    for (int i = 0; i < n + 2; ++i) {
        Bspline_coefficients(spline[i], newT[i], newT[i + 1], newT[i + 2], newT[i + 3], 1.0/6.0);
    }
}


void spline_points(vector<dpoint2>& pts, const cubic_spline_data& spd, int steps = 10)
{
    for (int ti = 0; ti <= steps; ++ti) {
        double t = (double)ti/steps;
        double t2 = t*t;
        double t3 = t2*t;

        pts.push_back(dpoint2(spd.a0 + t*spd.a1 + t2*spd.a2 + t3*spd.a3,
            spd.b0 + t*spd.b1 + t2*spd.b2 + t3*spd.b3));
    }
}

void spline_points(vector<dpoint2>& pts, const vector<cubic_spline_data>& spdv, int steps /* = 10 */)
{
    for (auto iter = spdv.begin(); iter != spdv.end(); ++iter) 
        spline_points(pts, *iter, steps);
}

void spline_points(vector<dpoint2>& pts, const vector<cubic_spline>& spdv, int steps /* = 10 */)
{
    for (auto iter = spdv.begin(); iter != spdv.end(); ++iter) 
        spline_points(pts, *iter, steps);
}

cv::Mat1d spline_image(const vector<cubic_spline>& splines, int thickness /* = 1 */, int imgborder /* = 50 */, double factor /* = 1.0 */)
{
    vector<dpoint2> allpts;

    for (auto iter = splines.begin(); iter != splines.end(); ++iter) 
        spline_points(allpts, *iter);
    for_each(allpts.begin(), allpts.end(), [factor](dpoint2& p) { p *= factor; });

    drectangle2 box = drectangle2::bounding_rectangle(allpts.begin(), allpts.end());

    cv::Mat1d img((int)box.y_dim() + 2*imgborder, (int)box.x_dim() + 2*imgborder, 0.0);

    for (auto iter = splines.begin(); iter != splines.end(); ++iter) {
        vector<dpoint2> pts;

        spline_points(pts, *iter);

        int imax = (int)pts.size() - 1;

        for (int i = 0; i < imax; ++i) {
            dpoint2 p = (pts[i]*factor - box.ll) + imgborder;
            dpoint2 q = (pts[i + 1]*factor - box.ll) + imgborder;

            cv::line(img, cv::Point(p.x, p.y), cv::Point(q.x, q.y), 255.0, thickness);
        }
    }
    return img;
}

double spline_length(cubic_spline& spline) 
{
    vector<dpoint2> pts;
    double length = 0.0;
    
    spline_points(pts, spline);
    for (int i = 0; i < (int)pts.size() - 1; ++i) {
        length += sqrt(pts[i].distance2(pts[i + 1]));
    }
    return length;
}

void keep_longest_splines(vector<cubic_spline>& splines, double thresh)
{
    if (splines.size() < 2) 
        return;

    vector<pair<double, int> > sizes(splines.size());
    vector<int> todelete;

    for (int i = 0; i < (int)splines.size(); ++i) {
        sizes[i].first = spline_length(splines[i]);
        sizes[i].second = i;
    }
    sort(sizes.begin(), sizes.end(), greater<pair<double, int> >());
    for (int i = 1; i < (int)sizes.size(); ++i) {
        if (sizes[i].first < thresh*sizes[0].first) todelete.push_back(sizes[i].second);
    }
    sort(todelete.begin(), todelete.end(), greater<int>());
    for (int i = 0; i < (int)todelete.size(); ++i) 
        splines.erase(splines.begin() + todelete[i]);

}

cv::Mat1d spline_image(const cubic_spline& splines, int thickness /* = 1 */, int imgborder /* = 50 */, double factor /* = 1.0 */)
{
    vector<cubic_spline> sspline;

    sspline.push_back(splines);
    return spline_image(sspline, thickness, imgborder, factor);
}

cv::Mat1d point_image(const vector<dpoint2>& pts, int radius /* = 1 */, int imgborder /* = 50 */)
{
    drectangle2 box = drectangle2::bounding_rectangle(pts.begin(), pts.end());
    cv::Mat1d img((int)box.y_dim() + 2*imgborder, (int)box.x_dim() + 2*imgborder, 0.0);

    for (auto piter = pts.begin(); piter != pts.end(); ++piter) {
        dpoint2 p = *piter - box.ll + dpoint2(imgborder, imgborder);

        cv::circle(img, cv::Point(p.x, p.y), radius, 255.0, -1);
    }
    return img;
}

void save_spline(const string& s, const vector<cubic_spline_data>& spline, int thickness /* = 1 */)
{
    const int imgborder = 50;

    cv::Mat img = spline_image(spline, thickness, imgborder, 100);

    cv::imwrite(s, img);
}

void splines_from_points(vector<cubic_spline>& splines, vector<ipoint2> pts, int inhibition)
{
    splines.clear();
    if (pts.empty()) 
        return;

    graph* pg = point_graph(pts);
    vector<node*> path = forest_longest_path(pg);

    while (path.size() > 2) { 
        vector<dpoint2> ppts;

        for (auto piter = path.begin(); piter != path.end(); ++piter) {
            node* n = *piter;
            int pi = ((node_data_t<int>*)n->data)->data;

            ppts.push_back((dpoint2)pts[pi]);
        }
        if (inhibition > 0 && pts.size() > 2) {
            ipoint2 p0 = ppts[0], p1 = ppts[ppts.size() - 1];
            
            ppts.erase(ppts.begin());
            ppts.erase(ppts.end() - 1);
            ppts = cast_vector<dpoint2, ipoint2>(inhibit_point_set(cast_vector<ipoint2, dpoint2>(ppts), inhibition));
            ppts.insert(ppts.begin(), p0);
            ppts.insert(ppts.end(), p1);
        }

        splines.push_back(cubic_spline());
        Bspline(splines.back(), ppts);

        set<node*> todelete; 

        for (auto piter = path.begin(); piter != path.end(); ++piter) {
            if ((*piter)->neighbors.size() < 3) todelete.insert(*piter);
        }

        pg->delete_nodes(todelete);
        path = forest_longest_path(pg);
    }

    delete pg;
}

void Bspline_tangents(vector<dpoint2>& tangents, const vector<dpoint2>& T)
{
    tangents.clear();

    if (T.empty())
        return;

    vector<dpoint2> newT(T.size() + 4);
    int n = T.size() - 1;

    newT[0] = line_point(T[0], T[1], -0.2);
    newT[1] = T[0];
    newT[2] = line_point(T[0], T[1], 0.2);
    newT[n + 2] = line_point(T[n - 1], T[n], 0.8);
    newT[n + 3] = T[n];
    newT[n + 4] = line_point(T[n - 1], T[n], 1.2);
    for (int i = 1; i < n; ++i) newT[i + 2] = T[i];

    tangents.resize(n + 1);
    for (int i = 1; i < n + 2; ++i) {
        cubic_spline_data spd;
        
        Bspline_coefficients_d(spd, newT[i], newT[i + 1], newT[i + 2], newT[i + 3], 1.0/6.0);
        tangents[i - 1].set(spd.a0, spd.b0);
        normalize(tangents[i - 1]);
    }

}

void tangents_from_points(vector<dpoint2>& tangents, const vector<ipoint2>& pts)
{
    tangents.resize(pts.size(), dpoint2::zero);
    if (pts.empty()) 
        return;

    graph* pg = point_graph(pts);
    vector<node*> path = tree_longest_path(pg);

    while (path.size() > 2) {
        vector<dpoint2> ppts;
        vector<int> indices;

        for (auto piter = path.begin(); piter != path.end(); ++piter) {
            node* n = *piter;
            int pi = ((node_data_t<int>*)n->data)->data;

            ppts.push_back((dpoint2)pts[pi]);
            indices.push_back(pi);
        }

        vector<dpoint2> bst;

        Bspline_tangents(bst, ppts);
        for (int i = 0; i < (int)ppts.size(); ++i) 
            tangents[indices[i]] = bst[i];

        set<node*> todelete; 

        for (auto piter = path.begin(); piter != path.end(); ++piter) {
            if ((*piter)->neighbors.size() < 3) todelete.insert(*piter);
        }

        pg->delete_nodes(todelete);
        path = tree_longest_path(pg);
    }

    delete pg;
}

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



