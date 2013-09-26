/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// miscellanous classes and functions

#include "platform.h"
#include "utils.h"
#include "misc.h"
#include "matching/matching.h"
#include "../graphs/graph_utils.h"

#include <string.h>
#include <cstdlib> 
#if defined WIN32 | defined WIN64
#include <windows.h>
#include <Psapi.h>
#else
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <fnmatch.h>
#endif


#include <iostream>
#include <cv.h>
#include <highgui.h>



//#pragma comment(linker, "/DEFAULTLIB:psapi.lib")


using namespace std;

// region
///////////////////////////////////////////////////////////////////////////////
/*
void region2::move(const ipoint2& p)
{
    for (set<ipoint2>::iterator iter = points.begin(); iter != points.end(); ++iter)
        *iter += p;
}

void region2::contract(double factor)
{
    if (factor >= 1.0) {
        for (iter_t iter = points.begin(); iter != points.end(); ++iter)
            *iter *= factor;
    } else {
        set<ipoint2> tmp;

        for (iter_t iter = points.begin(); iter != points.end(); ++iter)
            tmp.insert(*iter * factor);
        points = tmp;
    }

}
*/
void region2::add(const region2& r, const ipoint2& d, double factor)
{
    for (const_iter_t iter = r.points.begin(); iter != r.points.end(); ++iter) 
        points.insert((*iter + d) * factor);
}

void region2::add(const region2& r, const ipoint2& d)
{
    for (const_iter_t iter = r.points.begin(); iter != r.points.end(); ++iter) 
        points.insert(*iter + d);
}

void region2::assign(const matrix<double>& m, double thresh)
{
    points.clear();
    for_each_xy(m, i, j) {
        if (m(i, j) >= thresh) points.insert(ipoint2(i, j));
    }
}


// grid_points
////////////////

matrix<double> grid_point::default_filter(3, 3, 1);

grid_points::grid_points(int xs, int ys) : 
    xsize(max(0, xs)), ysize(max(0, ys)), points()
{
    if (xsize == 0 || ysize == 0) grid = 0;
    else {
        grid = new ppoint_t[xsize*ysize];
        memset(grid, 0, xsize*ysize*sizeof(ppoint_t));
    }
}

grid_points::grid_points(const grid_points& gp) : 
    xsize(gp.xsize), ysize(gp.ysize), points()
{
    if (!gp.grid) grid = 0;
    else {
        grid = new ppoint_t[xsize*ysize];
        memset(grid, 0, xsize*ysize*sizeof(ppoint_t));

        for (citer_t i = gp.begin(); i != gp.end(); ++i) {
            add(new point_t(**i));
        }
    }
}

grid_points::~grid_points()
{
    if (grid) {
        for (iter_t iter = points.begin(); iter != points.end(); ++iter) delete *iter;
        delete[] grid;
    }
}

int grid_points::intersection_size(const grid_points& gpts, double tol) const
{
    int result = 0;

    for (citer_t iter = points.begin(); iter != points.end(); ++iter) {
        if (gpts.convolve(**iter, point2<int>::zero) >= tol) ++result;
    }
    return result;
}

bool grid_points::match(const grid_points& gpts, double dtol, int ntol)
{
    int psize = gpts.size();
    int ssize = size();
    bool success = false;
    int nmatch;

    if (ntol < abs(psize - ssize)) return false;
    
    for (int i = min(ntol, psize - 1); i >= 0; --i) {
        ppoint_t mp = gpts.points[i];

        for (int j = ssize - 1; j >= 0; --j) {
            ppoint_t sp = points[j];
            double dmatch = 0.0;

            nmatch = 0;
            if (sp->type == mp->type) {
                point2<int> dp = *sp - *mp;

                for (citer_t iter = gpts.begin(); !success && iter != gpts.end(); ++iter) {
                    double conv = convolve(**iter, dp);

                    dmatch += conv;
                    if (conv >= dtol) ++nmatch;  
                    else if (psize - nmatch - ntol < j) break;
                }
            }
            if (psize - nmatch <= ntol) return true;
        }
    }
    return false;
}

int grid_points::best_match(const grid_points& gpts, double dtol, int ntol)
{
    int psize = gpts.size();
    int ssize = size();
    int result = 0;
    int nmatch;

    if (ntol < abs(psize - ssize)) return ntol + 1;
    
    for (int i = min(ntol, psize - 1); i >= 0; --i) {
        ppoint_t mp = gpts.points[i];

        for (int j = ssize - 1; j >= 0; --j) {
            ppoint_t sp = points[j];
            double dmatch = 0.0;

            nmatch = 0;
            if (sp->type == mp->type) {
                point2<int> dp = *sp - *mp;
                int k = psize;

                for (citer_t iter = gpts.begin(); iter != gpts.end(); ++iter) {
                    double conv = convolve(**iter, dp);

                    --k;
                    dmatch += conv;
                    if (conv >= dtol) ++nmatch;  
                    else if (k + nmatch + ntol < psize) break;
                }
            }
            if (nmatch > result) result = nmatch;
        }
    }
    return psize - result;
}

double grid_points::convolve(const grid_point& gp, const point2<int>& dp) const
{
    int gpx = gp.x + dp.x, gpy = gp.y + dp.y;
    int cx = (int)gp.filter->width/2;
    int cy = (int)gp.filter->height/2;
    int i0 = max(0, gpx - cx);
    int j0 = max(0, gpy - cy);
    int i1 = min(gpx + cx + 1, xsize);
    int j1 = min(gpy + cy + 1, ysize);
    int fi0 = i0 - (gpx - cx);
    int fj0 = j0 - (gpy - cy);
    matrix<double>& filter = *gp.filter;
    double result = 0.0;

    for (int i = i0, fi = fi0; i < i1; ++i, ++fi) {
        for (int j = j0, fj = fj0; j < j1; ++j, ++fj) {
            ppoint_t p = point_at(i, j);

            if (p != 0 && p->type == gp.type) { result += filter(fi, fj); }
        }
    }
    return result;
}

void grid_points::write_vgr(ostream& os) const
{
    if (points.empty()) return;

    int maxx, maxy;
    int minx, miny;
    ppoint_t p = points.front();

    maxx = minx = p->x;
    maxy = miny = p->y;
    for (citer_t iter = begin(); iter != end(); ++iter) {
        p = *iter;

        if (p->x > maxx) maxx = p->x; else if (p->x < minx) minx = p->x;
        if (p->y > maxy) maxy = p->y; else if (p->y < miny) miny = p->y;
    }

    int xsize = maxx - minx;
    int ysize = maxy - miny;

    // write vertices
    os << "*Vertices " << points.size() << endl;
    int count = 0;
    for (citer_t iter = begin(); iter != end(); ++iter) {
        p = *iter;

        ++count;
        os << count << ' '; // index
        os << "\"(" << p->x << ',' << p->y << ")-" << p->type << '\"'; // label
        os << ' ';
        os << (double)(p->x - minx)/xsize << ' ' << (double)(p->y - miny)/ysize << ' ' << 0.5 << endl;
    }

    // write edges
    os << "*Edges" << endl;
}

void grid_points::print(ostream& os) const
{
    for (citer_t iter = begin(); iter != end(); ++iter) {
        ppoint_t p = *iter;
        os << '(' << p->x << ',' << p->y << ")-" << p->type << ' ';
    }
}

// angular_list
///////////////////////////////////////////////////////////////////////////////

angle_list::angle_list() : angles(360, 0), minus(5), plus(5) { }

angle_list::angle_list(int m, int p) : angles(360, 0), minus(m), plus(p) { }

angle_list::angle_list(const angle_list& al) : angles(al.angles), minus(al.minus), plus(al.plus) { }

void angle_list::block_angle(int a, int m, int p, int value /* = 1 */)
{
    int end = a + p + 1;

    for (int i = a - m; i < end; ++i) {
        angles[mod(i, 360)] = value;
    }
}

// online_distribution
///////////////////////////////

void online_distribution::reset()
{
    n = 0; 
    mean = 0.0;
    M2 = 0.0;
}

void online_distribution::new_data(double x)
{
    double delta = x - mean;

    ++n;
    mean += delta/n;
    M2 += delta*(x - mean);
}

// normal_distribution2
/////////////////////////

void get_normal_distribution2(normal_distribution2& dist, const vector<ipoint2>& v)
{
    if (v.empty()) dist = normal_distribution2();

    dpoint2 mean(0.0, 0.0);
    matrix2x2 variance;

    for (vector<ipoint2>::const_iterator iter = v.begin(); iter != v.end(); ++iter) {
        const ipoint2& p = *iter;

        mean.x += p.x;    
        mean.y += p.y;
    }
    mean.x /= v.size();
    mean.y /= v.size();

    for (vector<ipoint2>::const_iterator iter = v.begin(); iter != v.end(); ++iter) {
        const ipoint2& p = *iter;
        dpoint2 v(p.x - mean.x, p.y - mean.y);

        variance.a += v.x*v.x;
        variance.b += v.x*v.y;
        variance.c += v.x*v.y;
        variance.d += v.y*v.y;
    }
    variance.a /= v.size();
    variance.b /= v.size();
    variance.c /= v.size();
    variance.d /= v.size();

    dist.reset(mean, variance);
}

// K_bin and shape context histogram
///////////////////////////////////////////////////////////////////////////////

// nangles: number of angles
// r: sorted vector of radii boundaries
K_bin::K_bin(int nangles, const vector<int>& rv)
{
    init_matrix(nangles, rv);
}

K_bin::K_bin()
{
    vector<int> rv;

    rv.push_back(20);
    init_matrix(6, rv);
}

K_bin::K_bin(int nangles, int r1)
{
    vector<int> rv;

    rv.push_back(r1);
    init_matrix(nangles, rv);
}

K_bin::K_bin(int nangles, int r1, int r2)
{
    vector<int> rv;

    rv.push_back(r1);
    rv.push_back(r2);
    init_matrix(nangles, rv);
}

K_bin::K_bin(int nangles, int r1, int r2, int r3)
{
    vector<int> rv;

    rv.push_back(r1);
    rv.push_back(r2);
    rv.push_back(r3);
    init_matrix(nangles, rv);
}

void K_bin::reset(int nangles, const vector<int>& rv)
{
    init_matrix(nangles, rv);
}

// nangles: number of angles
// r: sorted vector of radii boundaries
//
//   radii:  r[size - 1]  r[0]    0     r[0]   r[size - 1]
//              |.:........:......:......:..........:.|
//
//
void K_bin::init_matrix(int nangles, const vector<int>& rv)
{
    double angle = 2.0*M_PI/nangles;
    double adelta = angle/2.0;
    int msize = 2*rv.back() + 3;

    m.resize(msize, msize);
    nbins = (int)rv.size()*nangles;
    
    int mc = m.width/2;

    for (int i = 0; i < msize; ++i) {
        for (int j = 0; j < msize; ++j) {
            int x = i - mc, y = j - mc;

            if (x == 0 && y == 0) {
                m(i, j) = 0;
            } else {
                double a = atan2((double)y, (double)x) + adelta;

                if (a < 0) a += 2.0*M_PI;
                
                int aindex = (int)(a/angle);

                if (aindex >= nangles) aindex = nangles - 1; // just in case

                double r = sqrt((double)(x*x + y*y));
                int rindex = (int)rv.size();

                for (int k = 0; k < (int)rv.size(); ++k) {
                    if (r < rv[k]) { 
                        rindex = k; 
                        break; 
                    }
                }
                m(i, j) = rindex*nangles + aindex;
            }
        }
    }
}

// Return bin or -1
int K_bin::get_bin(const ipoint2& pt, const ipoint2& center /* = ipoint2::zero */) const
{
    ipoint2 p = pt + ipoint2((int)m.width/2 - center.x, (int)m.height/2 - center.y);

    if (p.x >= 0 && p.x < (int)m.width && p.y >= 0 && p.y < (int)m.height) {
        int b = m(p.x, p.y);

        if (b < nbins) return b;
    }
    return -1;
}

void K_bin::get_histogram(vector<int>& h, const ip2_vector& pts, const ipoint2& center /* = ipoint2::zero */) const
{
    ipoint2 pp((int)m.width/2 - center.x, (int)m.height/2 - center.y);

    h.resize(nbins, 0);
    for (ip2_vector::const_iterator iter = pts.begin(); iter != pts.end(); ++iter) {
        ipoint2 p = *iter + pp;

        if (p.x >= 0 && p.x < (int)m.width && p.y >= 0 && p.y < (int)m.height) {
            int bin = m(p.x, p.y);
            
            if (bin < nbins) ++h[bin];
        }
    }
}

//void K_bin::get_histogram(vector<int>& h, img_graph* g, int layer, const ipoint2& center) const
//{
//    if (!g->grid(layer)) 
//        g->init_grid(layer);
//
//    int mcx = (int)m.width/2, mcy = (int)m.height/2;
//    int dx = mcx - center.x, dy = mcy - center.y;
//    int maxx = min<int>(g->x_size(layer), center.x + mcx + 1);
//    int maxy = min<int>(g->y_size(layer), center.y + mcy + 1);
//    int minx = max<int>(0, center.x - mcx);
//    int miny = max<int>(0, center.y - mcy);
//
//    h.resize(nbins, 0);
//    for (int x = minx; x < maxx; ++x) {
//        for (int y = miny; y < maxy; ++y) {
//            node* n = g->node_at(x, y, layer);
//
//            if (n != nullptr) {
//                img_node_data* nd = (img_node_data*)n->data;
//                int bin = m(nd->x + dx, nd->y + dy);
//                
//                if (bin < nbins) ++h[bin];
//            }
//        }
//    }
//}

//vector<int> K_bin::get_histogram(const ip2_vector<int>& pts, const ipoint2& center /* = ipoint2::zero */)
//{
//    vector<int> result;
//
//    get_histogram(result, pts, center);
//    return result;
//}

void K_bin::print() const
{
    m.print_mathematica(cout);
}


void normalize_histogram(vector<double>& hn, const vector<int>& h)
{
    hn.resize(h.size(), 0.0);

    int count = 0;

    for (int i = 0; i < h.size(); ++i) 
        count += h[i];
    if (count == 0) 
        return;
    for (int i = 0; i < h.size(); ++i)
        hn[i] = (double)h[i]/count;

}

// Returns chi^2 distance between two histograms
double sc_histogram_distance(const vector<int>& h1, const vector<int>& h2)
{
    size_t maxk = min<int>(h1.size(), h2.size());
    double result = 0.0;

    for (size_t k = 0; k < maxk; ++k) {
        double den = h1[k] + h2[k];
        
        if (den != 0.0) {
            double x = h1[k] - h2[k];

            result += x*x/den;
        }
    }
    return result/2.0;
}

double sc_histogram_distance(const vector<double>& h1, const vector<double>& h2)
{
    size_t maxk = min<int>(h1.size(), h2.size());
    double result = 0.0;

    for (size_t k = 0; k < maxk; ++k) {
        double den = h1[k] + h2[k];
        
        if (den != 0.0) {
            double x = h1[k] - h2[k];

            result += x*x/den;
        }
    }
    return result/2.0;
}


// M_bin
///////////////////////////////////////////////////////////////////////////////

M_bin::M_bin(const string& fname) : m(), nbins(), centers(), areas(), max_distance2(), user_value(nullptr)
{
    matrix<int> fm = read_int_matrix(fname);

    if (!fm.empty()) init_matrix(fm);
}

// Get bin at (i, j), note that (0, 0) is center of "matrix".
int M_bin::get_bin(int i, int j) const
{
    int x = (int)m.width/2 + i, y = (int)m.height/2 + j;

    if (x < 0 || x >= (int)m.width || y < 0 || y >= (int)m.height)
        return -1;
    return m(x, y);
}

void M_bin::init_matrix(const matrix<int>& pm)
{
    m = pm;
    nbins = 0;
    for_each_element (m, i) {
        int b = m[i];

        if (b >= 0 && nbins <= b) nbins = b + 1;
    }
    centers.resize(nbins);
	areas.resize(nbins);
	max_distance2.resize(nbins);

    ipoint2 c((int)m.width/2, (int)m.height/2);

    for (int i = 0; i < nbins; ++i) {
        int count = 0;
        ipoint2& p = centers[i];
        
		int max_distance_px = 0;

        p = ipoint2::zero;
        for_each_xy_int (m, x, y) {
            int b = m(x, y);

            if (b == i) { p.add(x - c.x, y - c.y); ++count; }
        }

        if (count > 0) p /= count;

		float dx = 0, dy = 0;

		for_each_xy_int (m, x, y) {
            int b = m(x, y);
            if (b == i) {
				max_distance_px = max<int>(max_distance_px, p.distance2(ipoint2(x,y) - c));
			}
        }
		max_distance2[i] = max_distance_px;
		areas[i] = count;
    }
	
}

M_bin* M_bin::get_resized(int w, int h) const{
	printf("resizing bin from %d,%d to %d,%d\n",width(), height(), w,h);
	matrix<int> new_m(w,h);

	float width_scale = m.width / (float)new_m.width;
	float height_scale = m.height / (float)new_m.height;

	for (int y = 0; y < (int)(new_m).height; ++y) {
		for (int x = 0; x < (int)(new_m).width; ++x) {
			int old_x = max<int>((int)min<int>((int)ceil(x * width_scale - 0.5), (int)m.width-1), (int)0);
			int old_y = max<int>((int)min<int>((int)ceil(y * height_scale - 0.5), (int)m.height-1), (int)0);
			new_m(x,y) = m(old_x, old_y);
		}
	}
	M_bin* new_bin = new M_bin(new_m);
	new_bin->nbins = nbins;
	return new_bin;
}
// functions
///////////////////////////////////////////////////////////////////////////////

vector<ipoint2> inhibit_point_set(const vector<ipoint2>& pts, int radius)
{
    if (pts.empty()) return vector<ipoint2>();
    if (radius <= 0) return pts;

    irectangle2 box;
    vector<ipoint2> result;

    for (vector<ipoint2>::const_iterator piter = pts.begin(); piter != pts.end(); ++piter) 
        box.eat(*piter);

    matrix<bool> m(box.x_dim() + 1, box.y_dim() + 1);

    m.fill(true);
    for (vector<ipoint2>::const_iterator piter = pts.begin(); piter != pts.end(); ++piter) {
        ipoint2 p = *piter - box.ll;

        if (m(p.x, p.y)) {
            result.push_back(*piter);
            m.set_region_circ(p.x, p.y, radius, false);
        }
    }
    return result;
}

vector<pair<int, ipoint2> > inhibit_point_set(const vector<pair<int, ipoint2> >& pts, int radius)
{
    typedef pair<int, ipoint2> item_t;
    typedef vector<item_t> result_t;
    typedef map<int, vector<ipoint2> > map_t;

    if (pts.empty()) return result_t();
    if (radius <= 0) return pts;

    map_t ptsm;
    result_t result;

    result.reserve(pts.size());
    for (result_t::const_iterator piter = pts.begin(); piter != pts.end(); ++piter) {
        ptsm[piter->first].push_back(piter->second);
    }
    for (map_t::iterator miter = ptsm.begin(); miter != ptsm.end(); ++miter) {
        vector<ipoint2> ipts = inhibit_point_set(miter->second, radius);

        for (vector<ipoint2>::iterator viter = ipts.begin(); viter != ipts.end(); ++viter) 
            result.push_back(item_t(miter->first, *viter));
    }
    return result;
}

double ipoint2_set_distance(const vector<ipoint2>& v1, const vector<ipoint2>& v2)
{
    double result = numeric_limits<double>::max();

    for (vector<ipoint2>::const_iterator iter1 = v1.begin(); iter1 != v1.end(); ++iter1) {
        for (vector<ipoint2>::const_iterator iter2 = v2.begin(); iter2 != v2.end(); ++iter2) {
            int d = ipoint2::distance2(*iter1, *iter2);
            if (d < result) result = d;
        }
    }
    return result;
}

void gaussian_mask_n(int dimx, int dimy, double sigma, matrix<double>& mask)
{
    mask.resize(dimx, dimy);

    int centerx = dimx/2;
    int centery = dimy/2;
    double sum = 0.0;

    for_each_xy(mask, x, y) { 
        double d = (x - centerx)*(x - centerx) + (y - centery)*(y - centery);
        sum += mask(x, y) = ::exp(-d/2.0/sigma/sigma);
    }
    for_each_xy(mask, x, y) {
        mask(x, y) /= sum;
    }
}

void gaussian_mask(int dimx, int dimy, double sigma, matrix<double>& mask)
{
    mask.resize(dimx, dimy);

    int centerx = dimx/2;
    int centery = dimy/2;

    for_each_xy(mask, x, y) { 
        double d = (x - centerx)*(x - centerx) + (y - centery)*(y - centery);
        mask(x, y) = ::exp(-d/2.0/sigma/sigma);
    }
}

void expand_coordinates(iipair& p, double contraction)
{
    p.first = int_round(p.first*contraction); 
    p.second = int_round(p.second*contraction);
}

void expand_coordinates(const vector<iipair>& coorig, vector<iipair>& coo, double contraction)
{
    coo.assign(coorig.begin(), coorig.end());
    for (vector<iipair>::iterator iter = coo.begin(); iter != coo.end(); ++iter) {
        expand_coordinates(*iter, contraction);
    }
}

// Sum of dist(v1[i], v2[perm[i]]), i = 0, ..., v1.size()
// Returns -1 if perm[i] is not a valid index for v2 for some i
double distance_cost(const vector<int>& perm, const ip2_vector& v1, const ip2_vector& v2)
{
    double result = 0.0;

    if (perm.size() != v1.size()) return -1.0;
    for (int i = 0; i < (int)v1.size(); ++i) {
        int p = perm[i];

        if (p < 0 || p >= v2.size()) return -1.0;
        result += sqrt((double)v1[i].distance2(v2[p]));
    }
    return result;
}

// Calculates assignment (permutation) of points v1 to points v2 such that the sum of 
// distances is minimal
// Solves by examining all possibilities, so do not use it if v1 and v2 are large!
// Returns -1.0 if |v1| != |v2|
// Implement Hungarian method ;)
double min_distance_matching(vector<int>& perm, const ip2_vector& v1, const ip2_vector& v2)
{
    if (v1.size() != v2.size()) return -1.0;

    vector<int> tmpp = vector_range(0, (int)v1.size() - 1);
    double bestval = numeric_limits<double>::max();

    do {
        double val = distance_cost(tmpp, v1, v2);

        if (val < bestval) {
            bestval = val;
            perm = tmpp;
        }
    } while (next_permutation(tmpp.begin(), tmpp.end()));
    return bestval;
}

vector<int> inverse_permutation(vector<int>& perm)
{
    vector<int> result(perm.size(), 0.0);

    for (int i = 0; i < (int)perm.size(); ++i) {
        result[perm[i]] = i;
    }
    return result;
}

bool is_identity_permutation(const vector<int>& p)
{
    for (int i = 0; i < (int)p.size(); ++i) {
        if (p[i] != i) return false;
    }
    return true;
}

vector<int> identity_permutation(int n)
{
    vector<int> result(n, 0);

    for (int i = 0; i < n; ++i) 
        result[i] = i;
    return result;
}


// avector - result (vector of 360 values)
// angles - vector of sorted "boundary" angles; and angles[-1] > angles[0] must be true.
// avalues - values between angles; the first value is the value between angles[-1] and angles[0].
void make_angle_vector(vector<int>& avector, const vector<int>& angles, const vector<int>& values)
{
    avector.resize(360, -1);
    
    int size = (int)angles.size();
    int a = angles.back();
    int a1 = 360;
    int val = values.front();

    while (a < a1) avector[a++] = val;
    a = 0;
    for (int i = 0; i < size; ++i) {
        a1 = angles[i];
        val = values[i];
        while (a < a1) avector[a++] = val;
    }
 
}

void centralized_region_coordinates(set<iipair>& coo, const matrix<bool>& region)
{
    vector<iipair> vec;
    int avgi = 0, avgj = 0;
    int count = 0;

    for_each_xy(region, i, j) {
        if (region(i, j)) {
            avgi += (int)i; avgj += (int)j;
            ++count;
            vec.push_back(iipair((int)i, (int)j));
        }
    }
    iipair avg((int)((double)avgi/count), (int)((double)avgj/count));

    coo.clear();
    for (vector<iipair>::iterator iter = vec.begin(); iter != vec.end(); ++iter) {
        coo.insert(*iter - avg);
    }
}

void centralized_region_coordinates(set<wpoint>& coo, const matrix<double>& region)
{
    const double eps = 10E-6;
    vector<wpoint> vec;
    double avgi = 0.0, avgj = 0.0;
    double mass = 0;

    for_each_xy(region, i, j) {
        double r = region(i, j);
        if (r > eps || r < -eps) {
            avgi += i*r; avgj += j*r;
            mass += r;
            vec.push_back(wpoint((int)i, (int)j, r));
        }
    }
    iipair avg((int)(avgi/mass), (int)(avgj/mass));

    coo.clear();
    for (vector<wpoint>::iterator iter = vec.begin(); iter != vec.end(); ++iter) {
        coo.insert(*iter - avg);
    }
}
/*
void scale_region(set<wpoint>& coo, double factor)
{
    for (set<wpoint>::iterator siter = coo.begin(); siter != coo.end(); ++siter) {
        siter->i = (int)(siter->i*factor); 
        siter->j = (int)(siter->j*factor);
    }
}
*/
double intersection_weight(const set<wpoint>& s1, const set<wpoint>& s2)
{
    set<wpoint>::const_iterator i1 = s1.begin(), i2 = s2.begin();
    double result = 0.0;
	double sqsum = 0.0;

    while (i1 != s1.end() && i2 != s2.end()) {
        if (*i1 < *i2) ++i1; 
        else if (*i2 < *i1) ++i2;
        else {
			result += (i1->w)*(i2->w);
			sqsum += (i1->w)*(i1->w) + (i2->w)*(i2->w);
            ++i1; ++i2;
        }
    }
    return result/sqsum;
}


void random_seed(int seed)
{
	srand((unsigned)seed);
}

void random_permutation(vector<int>& perm, int n)
{
    int* arr = new int[n];
    int i, p, t;

    for (i = 0; i < n; ++i) arr[i] = i; 
    for (i = n - 1; i > 0; --i) {
        p = rand()%(i + 1);
        t = arr[p];
        arr[p] = arr[i];
        arr[i] = t;
    }
    array_to_vector(perm, arr, n);
    delete[] arr;
}

int random_discrete(const vector<double>& d)
{
    if (d.empty())
        return -1;

    double sum = 0;

    for (size_t i = 0; i < d.size(); ++i) sum += d[i];

    double rnd = random_real()*sum;
    int result = 0;

    sum = d[0];
    while (sum < rnd) sum += d[++result];
    return result;
}

// "Fits" an ellipse to a set of points.
// a is the larger of half-axes, b the smaller of half-axes and angle is from (-Pi/2, Pi/2]
void fit_ellipse(double& cx, double& cy, double& a, double& b, double& angle, const vector<dpoint2>& points)
{
    if (points.empty()) {
        cx = cy = 0;
        a = b = angle = 0.0;
        return;
    }

    double x, y;
    double mata = 0.0, matb = 0.0, matd = 0.0;
    dpoint2 center = dpoint2::center(points.begin(), points.end());
    double N = points.size();

    for (size_t i = 0; i < points.size(); ++i) {
        x = points[i].x - center.x;
        y = points[i].y - center.y;
        
        mata += x*x;
        matb += x*y;
        matd += y*y;
    }

    double l1, l2;
    dpoint2 v1, v2;

    eigensystem(l1, l2, v1, v2, mata, matb, matd);

    //cout << "l1 = " << l1 << endl;
    //cout << "l2 = " << l2 << endl;
    //cout << "v1 = " << v1 << endl;
    //cout << "v2 = " << v2 << endl;

    if (l1 >= l2) { 
        a = sqrt(l1/N); b = sqrt(l2/N); 
        angle = atan2(v1.y, v1.x);
    } else {
        a = sqrt(l2/N); b = sqrt(l1/N);
        angle = atan2(v2.y, v2.x);
    }
    if (angle < -M_PI_2) angle += M_PI;
    else if (angle > M_PI_2) angle -= M_PI;
    cx = center.x;
    cy = center.y;
}

void fit_ellipse(double& cx, double& cy, double& a, double& b, double& angle, const set<ipoint2>& points)
{
    vector<dpoint2> dpoints;

    dpoints.reserve(points.size());
    for (set<ipoint2>::const_iterator iter = points.begin(); iter != points.end(); ++iter) {
        dpoints.push_back(dpoint2(iter->x, iter->y));
    }
    fit_ellipse(cx, cy, a, b, angle, dpoints);
}

void fit_ellipse(ellipse& ell, const vector<dpoint2>& points)
{
    fit_ellipse(ell.cx, ell.cy, ell.a, ell.b, ell.angle, points);
}

void fit_ellipse(ellipse& ell, const set<ipoint2>& points)
{
    fit_ellipse(ell.cx, ell.cy, ell.a, ell.b, ell.angle, points);
}

// "Normalize" vector (translate & scale).
pair<double, dpoint2> translate_and_scale(vector<dpoint2>& v)
{
    if (v.empty()) return pair<double, dpoint2>();

    dpoint2 avg = dpoint2::center(v.begin(), v.end());
    double scale = 0.0;

    for (vector<dpoint2>::iterator iter = v.begin(); iter != v.end(); ++iter) {
        scale += iter->distance2(avg);
    }
    scale = ::sqrt(scale/(double)v.size());
    for (vector<dpoint2>::iterator iter = v.begin(); iter != v.end(); ++iter) {
        *iter -= avg;
        iter->div(scale);
    }
    return pair<double, dpoint2>(scale, avg);
}

pair<double, dpoint2> translate_and_scale(vector<pair<int, ipoint2> >& pts)
{
    vector<ipoint2> v = extract_second<ipoint2>(pts.begin(), pts.end());
    auto result = translate_and_scale(v, 100.0);

    replace_second(pts, v);
    return result;
}

pair<double, dpoint2> translate_and_scale(vector<ipoint2>& v, double factor /* = 100 */)
{
    vector<dpoint2> pts = cast_vector<dpoint2, ipoint2>(v);
    auto result = translate_and_scale(pts);

    v = cast_vector<ipoint2, dpoint2>(pts, factor);
    return result;
}

// Rotate vector 'v' for angle 'a' about origin.
void rotate(vector<dpoint2>& v, double a)
{
    double cosa = cos(a);
    double sina = sin(a);

    for (vector<dpoint2>::iterator iter = v.begin(); iter != v.end(); ++iter) {
        iter->set(cosa*iter->x - sina*iter->y, sina*iter->x + cosa*iter->y);
    }
}

// Optimal rotation of v -> ref.
double optimal_rotation(const vector<dpoint2>& v, const vector<dpoint2>& ref)
{
    if (v.empty() || v.size() != ref.size()) 
        return 0.0;

    double numer = 0.0, denom = 0.0;
    
    for (vector<dpoint2>::const_iterator viter = v.begin(), refiter = ref.begin(); viter != v.end(); ++viter, ++refiter) {
        numer += viter->x*refiter->y - viter->y*refiter->x;
        denom += viter->x*refiter->x + viter->y*refiter->y;
    }
    return atan2(numer, denom);
}

struct mp_plus {
    double d;
    mp_plus(double vd = 0.0) : d(vd) { }

    double operator()(const double& x, const double& y) const
    {
        return x + d*y;
    }
};

// place mask2 to each point of mask one
void mask_product(matrix<double>& result, matrix<double>& mask1, matrix<double>& mask2, int cx2, int cy2)
{
    const double eps = 10E-06;

    double d;
    mp_plus pf;

    result.resize(mask1.width + 2*mask2.width, mask1.height + 2*mask2.height);
    result.blt_central(mask1, (int)mask1.width/2, (int)mask1.height/2, (int)result.width/2, (int)result.height/2);
    for_each_xy(result, i, j) {
        d = mask1(i, j);
        if (d > eps || d < -eps) {
            pf.d = d;
            result.blt_central(mask2, cx2, cy2, (int)i, (int)j, pf);
        }
    }

    
}

bool delete_file(const string& path)
{
    return remove(path.c_str()) == 0;
}

unsigned long get_process_id()
{
#if defined WIN32 | defined WIN64
    return GetCurrentProcessId();
#else
	return getpid();
#endif
}

unsigned long get_working_set_size()
{
#if defined WIN32 | defined WIN64
    PROCESS_MEMORY_COUNTERS lMemInfo; 

    GetProcessMemoryInfo(GetCurrentProcess(), &lMemInfo, sizeof(lMemInfo));
    return (unsigned long)lMemInfo.WorkingSetSize;
#else
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    //http://www.gnu.org/s/libc/manual/html_node/Resource-Usage.html
    return (unsigned long) (usage.ru_idrss + usage.ru_ixrss + usage.ru_isrss) * 1024;
#endif
}

void to_unix_path(string& s)
{
    for_each(s.begin(), s.end(), [](char& c) { if (c == '\\') c = '/'; });
}

string wilcard_to_regex(const string& s)
{
    string result;

    for (string::const_iterator iter = s.begin(); iter != s.end(); ++iter) {
        if (*iter == '*') result += ".*";
#if defined WIN32 | defined WIN64
        else if (*iter == '\\') result += "/";
#endif
        else if (*iter == '?') result += '.';
        else if (*iter == '.') result += "\\.";
        else result += *iter;
    }
    return result;
}

#if defined WIN32 | defined WIN64

void list_directory_no_delete(list<string>& result, const string& pattern, const string& prefix)
{
    WIN32_FIND_DATA w32fd;
    HANDLE hFind = INVALID_HANDLE_VALUE;

    hFind = FindFirstFile(pattern.c_str(), &w32fd);
    if (hFind == INVALID_HANDLE_VALUE) 
        return;
    do {
        if (!(w32fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) 
            result.push_back(prefix + w32fd.cFileName);
    } while (FindNextFile(hFind, &w32fd) != 0);
    FindClose(hFind);
}

void list_directory(list<string>& names, const string& rpath)
{
    names.clear();
    list_directory_no_delete(names, rpath, "");
}
/*
void list_directory(list<string>& result, const string& srcdir, const string& rx, const string& prefix)
{
    WIN32_FIND_DATA w32fd;
    HANDLE hFind = INVALID_HANDLE_VALUE;

    hFind = FindFirstFile((srcdir + "*").c_str(), &w32fd);
    if (hFind == INVALID_HANDLE_VALUE) 
        return;
    do {
        if ((w32fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) && 
                strcmp(w32fd.cFileName, ".") != 0 && strcmp(w32fd.cFileName, "..") != 0) {
            list_directory(result, srcdir + w32fd.cFileName + "\\", rx, string(w32fd.cFileName) + "\\");
        }
    } while (FindNextFile(hFind, &w32fd) != 0);
    FindClose(hFind);
    list_directory_no_delete(result, srcdir + rx, prefix);
}

void list_directory(list<string>& names, const string& srcdir, const string& pattern)
{
    list_directory(names, srcdir, pattern, "");
}
*/

const char* get_current_path() {
	int buffer_length = GetCurrentDirectory(0,nullptr);

#ifdef _UNICODE
	LPTSTR buffer = (LPTSTR)malloc(sizeof(wchar_t) * buffer_length);
#else
	LPTSTR buffer = (LPTSTR)malloc(sizeof(char) * buffer_length);
#endif
	
	GetCurrentDirectory(buffer_length, buffer);
#ifdef _UNICODE
	// THIS MAY NOT HAVE BEEN TESTED !!
	int wBuffer_len = WideCharToMultiByte(CP_ACP, 0, buffer, -1, nullptr, 0, nullptr, nullptr);

	LPTSTR wBuffer = (LPTSTR)malloc(sizeof(char) * wBuffer_len);
	WideCharToMultiByte(CP_ACP, 0, buffer, -1, wBuffer, wBuffer_len, nullptr, nullptr);

	return (char*)wBuffer;
#else
	return (char*)buffer;
#endif
}

#else

void list_directory(list<string>& names, const string& rpath)
{
    DIR *dp;
    struct dirent *ep;

	size_t sep_pos = rpath.find_last_of("/\\");
	
	string folder, filename;
	if (sep_pos < 0) {
		folder = "./";
		filename = rpath;
	} else {
		folder = rpath.substr(0,sep_pos);
		filename = rpath.substr(sep_pos+1);
	}
    const char* pattern = filename.c_str();

    dp = opendir(folder.c_str());

    if (dp != nullptr)
	{
        while (ep = readdir (dp)) {
			if (fnmatch(pattern, ep->d_name, 0) == 0)
				names.push_back(string(ep->d_name));
		}
        (void) closedir (dp);
    } else
        perror ("Couldn't open the directory");
}

const char* get_current_path() {
	return (char*)get_current_dir_name();
}


#endif

bool create_directory(const string& folder) {	
	bool ok = false;
#if defined WIN32 | defined WIN64		
	ok = CreateDirectory(folder.c_str(), nullptr) > 0 ? true : false; // returns 0 on error
#else
	ok = mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0 ? true : false; // return 0 on success
#endif
	return ok;
}
bool create_directories(const string& folder) {	
	
	bool ok = false;
	// create dir
	if (folder.length() > 0 ) {
		// if parent dir does not exists then call create directory on parent
		ok = create_directory(folder);

		if (ok == false) {
			const string parent = get_parent_of_path(folder);		
			if (parent.length() > 0 && create_directories(parent))
				ok = create_directory(folder);
		}
	}
	return ok;
}

const string get_filename_of_path(const string& path) {
	size_t sep_pos = path.find_last_of("/\\");
	return (sep_pos == string::npos) ? path : path.substr(sep_pos+1);
}

const string get_parent_of_path(const string& path) {
	size_t sep_pos = path.find_last_of("/\\");
	return (sep_pos == string::npos) ? "" : path.substr(0,sep_pos);
}

bool file_list_from_file(list<string>& l, const string& fname, const string& dir, const string& ext /* = "" */)
{
    if (ext.empty())
        list_from_file(l, fname, "");
    else {
        list<string> pl;

        list_from_file(pl, fname, "");
        for (list<string>::iterator iter = pl.begin(); iter != pl.end(); ++iter) {
            char buf[1024];
            list<string> fnl;

            sprintf(buf, ext.c_str(), iter->c_str());
            list_directory(fnl, dir + buf);
            l.insert(l.end(), fnl.begin(), fnl.end());
        }
    }
    return !l.empty();
}

bool file_lists_from_list(list<list<string> >& l, const list<string>& pl, const string& dir, const string& ext)
{
    for (list<string>::const_iterator iter = pl.begin(); iter != pl.end(); ++iter) {
        char buf[1024];
        list<string> fnl;

        sprintf(buf, ext.c_str(), iter->c_str());
        list_directory(fnl, dir + buf);
        l.push_back(fnl);
    }
    return !l.empty();
}

bool file_lists_from_file(list<list<string> >& l, const string& fname, const string& dir, const string& ext)
{
    list<string> pl;

    list_from_file(pl, fname, "");
    return file_lists_from_list(l, pl, dir, ext);
}

int matrix_union(const rmatrix& m1, const rmatrix& m2, const ipoint2& sp, const ipoint2& dp, double tol)
{
    rmatrix tmp(m1);
    int result = 0;

    tmp.blt_central2_ar(m2, sp.x, sp.y, dp.x, dp.y, max_function<double>());
    for_each_iter(tmp, rmatrix::iterator, iter) {
        if (*iter >= tol) ++result;
    }
    return result;
}

matrix<int> read_int_matrix(const string& s)
{
    int mwidth = 0;
    vector<int> v;

    ifstream is(s.c_str());
    
    bool fail = is.fail();
    int wcount = 0;

    while (!fail && !is.eof()) {
        char c;
        bool newline = false;

        do {
            c = is.get();
            if (c == '\n') newline = true;
        } while (c == ' ' || c == '\t' || c == '\n');
        if (!is.eof()) {

            if (!newline) ++wcount;
            else {
                if (mwidth == 0) { 
                    if (wcount != 0) mwidth = wcount; 
                } else if (mwidth != wcount) fail = true; 
                wcount = 1;
            }

            if (!fail) {
                int elt;

                is.putback(c);
                is >> elt;

                fail = is.fail();
                if (!fail) v.push_back(elt);
            }
        }
    }
    is.close();

    if (fail || mwidth == 0 || ((int)v.size()) % mwidth != 0) 
        return matrix<int>();

    matrix<int> m(mwidth, (int)v.size()/mwidth);
    
    m.assign(v.begin(), v.end());
    return m;
}


matrix<float> read_float_matrix(const string& s)
{
    int mwidth = 0;
    vector<float> v;

    ifstream is(s.c_str());
    
    bool fail = is.fail();
    int wcount = 0;

    while (!fail && !is.eof()) {
        char c;
        bool newline = false;

        do {
            c = is.get();
            if (c == '\n') newline = true;
        } while (c == ' ' || c == '\t' || c == '\n');
        if (!is.eof()) {

            if (!newline) ++wcount;
            else {
                if (mwidth == 0) { 
                    if (wcount != 0) mwidth = wcount; 
                } else if (mwidth != wcount) fail = true; 
                wcount = 1;
            }

            if (!fail) {
                float elt;

                is.putback(c);
                is >> elt;

                fail = is.fail();
                if (!fail) v.push_back(elt);
            }
        }
    }
    is.close();

    if (fail || mwidth == 0 || ((int)v.size()) % mwidth != 0) 
        return matrix<float>();

    matrix<float> m(mwidth, (int)v.size()/mwidth);
    
    m.assign(v.begin(), v.end());
    return m;
}


ipoint2 simple_shape_matching(const rmatrix& m1, const rmatrix& m2, 
    const ipoint2& c1, const ipoint2& c2, int tol)
{
    ipoint2 result = c2;
    int val = matrix_union(m1, m2, c1, c2, 0.1);

    for (int i = 0; i < tol; ++i) {
        for (int j = 0; j < tol; ++j) {
            int v = matrix_union(m1, m2, c1, c2 + ipoint2(i, j), 0.1);
            if (v < val) {
                result = c2 + ipoint2(i, j);
                val = v;
            }
        }
    }
    return result;
}


void save_region_set(ostream& os, set<int>& sm, int x_size, int y_size)
{
    if (sm.empty()) return;
    for (int j = 0; j < y_size; ++j) {
        for (int i = 0; i < x_size; ++i) {
            if (sm.find(j*x_size + i) != sm.end()) os << '#'; else os << '.';
            os << ' ';
        }
        os << endl;
    }
}

void save_region_set(const string& fname, set<int>& sm, int x_size, int y_size)
{
    ofstream os(fname.c_str());
    save_region_set(os, sm, x_size, y_size);
    os.close();
}

void save_region_sets(ostream& os, matrix<set<int> >& sm, int x_size, int y_size)
{
    for (int j = 0; j < sm.height; ++j) 
        for (int i = 0; i < sm.width; ++i) {
            os << i << ' ' << j << endl;
            save_region_set(os, sm(i, j), x_size, y_size);
        }
}


// Returns a vector of integers (a 'permutation') such that 
// pts1[perm[i]] <-> pts2[i] is optimal assignment.
vector<int> point_matching(const vector<dpoint2>& pts1, const vector<dpoint2>& pts2)
{
	if (pts1.empty() || pts2.empty()) return vector<int>();

    matrix<double> dist(pts2.size(), pts1.size());
    matrix<int> match(pts2.size(), pts1.size());
    vector<int> result(pts2.size());

    for (int i = 0; i < (int)pts1.size(); ++i) 
        for (int j = 0; j < (int)pts2.size(); ++j) 
            dist(j, i) = /*sqrt(*/(double)pts1[i].distance2(pts2[j])/*)*/;
    best_matching(match.ptr(0, 0), dist.height, dist.width, dist.ptr(0, 0));
    for (int i = 0; i < (int)pts1.size(); ++i) {
        for (int j = 0; j < (int)pts2.size(); ++j) 
            if (match(j, i) == 1) {
                result[j] = i;
            }
    }
    return result;
}

vector<int> point_matching(const vector<ipoint2>& pts1, const vector<ipoint2>& pts2)
{
    return point_matching(cast_vector<dpoint2, ipoint2>(pts1), cast_vector<dpoint2, ipoint2>(pts2));
}


int filled_region_size(const vector<ipoint2>& ptsvec)
{
    if (ptsvec.empty()) 
        return 0;

    vector<ipoint2> pts = ptsvec;
    irectangle2 box = irectangle2::bounding_rectangle(pts.begin(), pts.end());

    for_each(pts.begin(), pts.end(), [box](ipoint2& p) { p -= box.ll; });

    graph* pg = point_graph(pts);
    vector<node*> path = tree_longest_path(pg);

    cv::Point* cpts = new cv::Point[path.size()];
    int npts = (int)path.size();

    for (int ni = 0; ni < (int)path.size(); ++ni) {
        node* pn = path[ni];
        int pi = ((node_data_t<int>*)pn->data)->data;

        cpts[ni].x = pts[pi].x;
        cpts[ni].y = pts[pi].y;
    }

    cv::Mat1b im(box.y_dim(), box.x_dim(), (uchar)0);
    int result = 0;

    cv::fillPoly(im, (const cv::Point**)(&cpts), &npts, 1, cv::Scalar(255));
    for (auto imiter = im.begin(); imiter != im.end(); ++imiter) {
        if (*imiter > 0) ++result;
    }
    //cv::imwrite("c:\\work\\tmp.png", im);

    delete cpts;
    delete pg;

    return result;
}

// Reads matching from comma separated file:
// Format of file:
// p[0], srcx[0], srcy[0], destx[0], desty[0], ...
// p[1], srcx[1], srcy[1], destx[1], desty[1], ...
// ...
// where p is a permutation of src that determines the matching between src and dest; 
// i.e. src[p[i]] matches with dest[i] for all i.
void read_matching(vector<int>& perm, vector<dpoint2>& src, vector<dpoint2>& dest, const string& fname)
{
    ifstream is;

    is.open(fname.c_str());

    perm.clear();
    src.clear();
    dest.clear();

    while (is.good()) {
        int p;
        double srcx, srcy, destx, desty;

        is >> p;
        is.ignore(INT_MAX, ',');
        is >> srcx;
        is.ignore(INT_MAX, ',');
        is >> srcy;
        is.ignore(INT_MAX, ',');
        is >> destx;
        is.ignore(INT_MAX, ',');
        is >> desty;
        if (is.good()) {
            perm.push_back(p - 1);
            src.push_back(dpoint2(srcx - 1, srcy - 1));
            dest.push_back(dpoint2(destx - 1, desty - 1));
        }
        is.ignore(INT_MAX, '\n');
    }
    is.close();
}

// Mean shift algorithm
// 'path': result; front() == init
// 'pv': vector of points
// 'init': initial value
// 'radius': neighborhood of current 'x'
// 'c': normal kernel parameter in e^{c t^2}
// 'eps', 'maxsteps': iteration limits
// returns last 'eps' i.e.: |x_n - x_{n-1}|
double mean_shift(vector<dpoint2>& path, const vector<dpoint2>& pv, const dpoint2& init, 
    double radius, double c, double eps, int maxsteps)
{
    double radius2 = radius*radius;
    double eps2 = eps*eps;
    double delta2 = numeric_limits<double>::infinity();
    dpoint2 x = init;

    path.clear();
    path.push_back(x);
    for (int s = 0; s < maxsteps; ++s) {
        dpoint2 num = dpoint2::zero;
        double denom = 0.0;

        for (vector<dpoint2>::const_iterator iter = pv.begin(); iter != pv.end(); ++iter) {
            const dpoint2& p = *iter;
            double dxn = (p - x).norm2();

            if (dxn <= radius2) {
                double k = ::exp(c*dxn);
                num += p*k;
                denom += k;
            }
        }
        num.x /= denom;
        num.y /= denom;
        delta2 = x.distance2(num);
        x = num;
        path.push_back(x);
        if (delta2 < eps2)
            break;
        //cout << x << endl;
    }
    return sqrt(delta2);
}

// streamed_pointer
///////////////////////////////////////////////////////////////////////////////

int streamed_pointer::pid = -1;
int streamed_pointer::nextid = 0;

string streamed_pointer::tmp_dir = TEMPORARY_DIRECTORY;

streamed_pointer::streamed_pointer(streamable* s)
{
    if (s == nullptr) namep = nullptr; 
    else create_name(s);
}

streamable* streamed_pointer::get()
{ 
    if (namep == nullptr) return nullptr; 
    else return streamable::read(namep->name); 
}

void streamed_pointer::set(streamable* p)
{
    if (namep != nullptr) p->save(namep->name, -1);
    else create_name(p);
}

streamed_pointer* streamed_pointer::from_file(const string& file, bool is_file_disposable /*= true*/)
{
	streamed_pointer* ptr = new streamed_pointer();
	ptr->create_from_name(file);
	if (ptr->is_null() == false)
		ptr->namep->disposable = is_file_disposable;
	return ptr;
}

void streamed_pointer::copy(const streamed_pointer& sp)
{
    if (sp.namep == namep) return;
    dispose(); 
    namep = sp.namep;
    if (namep != nullptr) ++namep->count;
}

void streamed_pointer::dispose()
{
    if (namep != nullptr && --namep->count <= 0) delete namep;
}

void streamed_pointer::create_name(streamable* p)
{
    this->create_from_name(generate_name(nextid));
    p->save(namep->name, -1);
}
void streamed_pointer::create_from_name(const string& name)
{
    namep = new counted_name(name);
    nextid++;
}
void streamed_pointer::read_from_stream(istreamer& is) {
	bool has_namep;
	is.read(has_namep);

	if (has_namep) {
		string name;
		is.read(name);

		namep = new counted_name(name);
	}
}
void streamed_pointer::write_to_stream(ostreamer& os) {
	os.write(namep != nullptr ? true : false);
	if (namep != nullptr)
		os.write(namep->name);

}