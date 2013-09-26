/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// layer_1 

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>

#if defined WIN32 | defined WIN64
#include <crtdbg.h>
#endif

#include <iostream>
#include <fstream>
#include <algorithm>
#include <stack>
#include <math.h>

#include "../utils/utils.h"
#include "../graphs/graph_utils.h"
#include "part_lib.h"
#include "layer_1.h"


//#define DEBUG_OUTPUT

using namespace std;

// local definitions
///////////////////////////////////////////////////////////////////////////////

/*inline bool logical_or(const bool& a, const bool& b) 
{
    return a || b;
}*/

struct part_str_comparer : public binary_function<part_str, part_str, bool> {
    int dist;
  
    part_str_comparer(int d = 0) : dist(d) { }

    bool operator()(const part_str& p, const part_str& q) const
        { return p.type == q.type && abs(p.x - q.x) <= dist && abs(p.y - q.y) <= dist; }
};

// thresh_data
///////////////////////////////////////////////////////////////////////////////

void thresh_data::save_mma(ostream& os)
{
    os << '{';
    for (container_t::iterator iter = tmap.begin(); iter != tmap.end(); ++iter) {
        if (iter != tmap.begin()) os << ',';
        os << *iter;
    }
    os << '}';
}

void thresh_data::write_to_stream(ostreamer& os)
{
    os.write((unsigned)tmap.size());
    for (container_t::iterator iter = tmap.begin(); iter != tmap.end(); ++iter) {
        os.write(iter->first); 
        os.write(iter->second);
    }
}

void thresh_data::read_from_stream(istreamer& is)
{
    unsigned size;

    is.read(size);
    for (unsigned i = 0; i < size; ++i) {
        container_t::key_type k;
        container_t::mapped_type v;

        is.read(k);
        is.read(v);
        tmap.insert(container_t::value_type(k, v));
    }
}

// lib_data
///////////////////////////////////////////////////////////////////////////////

int lib_data::get_basic_part_number(node* n)
{
    if (bpcount <= 0) {
        if (layer == 0) bpcount = 1;
        else if (n->is_attr_set(VS_PART_ATTR)) {
            int edge = atom("lyrVSMember");
            int count = 0;

            bpcount = 0;
            foreach_neighbor(n, edge, iter) {
                node* nn = neighbor_node(iter);
                lib_data* nnd = (lib_data*)neighbor_node_data(iter);

                bpcount += nnd->get_basic_part_number(nn);
                ++count;
            }
            bpcount /= count;
        } else {
            int cedge = atom("lyrCenter");
            int sedge = atom("lyrSrc");

            bpcount = 0;
            foreach_neighbor(n, sedge, iter) {
                lib_data* nnd = (lib_data*)neighbor_node_data(iter);
                bpcount += nnd->get_basic_part_number(neighbor_node(iter));
            }
            foreach_neighbor(n, cedge, iter) {
                lib_data* nnd = (lib_data*)neighbor_node_data(iter);
                bpcount += nnd->get_basic_part_number(neighbor_node(iter));
            }
        }
    }
    return bpcount;
}

// part_data
///////////////////////////////////////////////////////////////////////////////

vector<double> part_data::contractions(MAX_LAYER_NUMBER, 1.0);

part_data::part_data(const img& m, const matrix<bool>& r, int l, int t) 
    : cmx(0), cmy(0), mask(), region(), lib_data(l, t), ocl_struct(nullptr)
{
    m.to_point_map(mask, 0.0);
    r.to_point_set(region, false);
}

part_data::part_data(const part_data& pd) :
    lib_data(pd),
    cmx(pd.cmx), cmy(pd.cmy), mask(pd.mask), region(pd.region), ocl_struct(nullptr)
{ }


void part_data::set_contractions(const vector<double>& cv)
{
    size_t max = min(contractions.size(), cv.size());

    for (size_t i = 0; i < max; ++i) contractions[i] = cv[i];
}

void part_data::set_contractions(double carr[], int size)
{
    int max = min((int)contractions.size(), size);

    for (int i = 0; i < max; ++i) contractions[i] = carr[i];
}

double part_data::total_contraction(int layer) 
{
	double result = 1.0;

	for (int i = 1; i <= layer; ++i) result *= contractions[i];

	return result;
}

void part_data::read_from_stream(istreamer& is)
{
    lib_data::read_from_stream(is);
    is.read(cmx); is.read(cmy);
    if (layer == 0) {
        is.read(mask);
        is.read(region);
    } 
}

void part_data::write_to_stream(ostreamer& os)
{
    lib_data::write_to_stream(os);
    os.write(cmx); os.write(cmy);
    if (layer == 0) {
        os.write(mask);
        os.write(region);
    }
}

void part_data::merge_mask_maps(map<ipoint2, double>& result, const map<ipoint2, double>& src, 
    const ipoint2& delta, const double& factor, unsigned libtype)
{
    typedef map<ipoint2, double> mask_t;

    for (mask_t::const_iterator miter = src.begin(); miter != src.end(); ++miter) {
        ipoint2 p = (miter->first + delta) * factor;
        pair<mask_t::iterator, bool> iiter = result.insert(mask_t::value_type(p, miter->second));

        if (!iiter.second) 
            iiter.first->second = 
                (libtype == APP_LIB_ATTR) ? iiter.first->second + miter->second : std::max<double>(iiter.first->second, miter->second); 
    }
}

void part_data::make_mask(node* n, unsigned libtype, part_lib* library)
{
    if (!mask.empty()) return;

    if (n->is_attr_set(VS_PART_ATTR)) {
        int name = atom("lyrVSMember");

        foreach_neighbor (n, name, iter) {
            node* nn = neighbor_node(iter);
            part_data* nnd = (part_data*)neighbor_node_data(iter);

            nnd->make_mask(nn, libtype, library);
            merge_mask_maps(mask, nnd->mask, ipoint2::zero, 1.0, libtype);
        }
        return;
    }

    int srcname = atom("lyrSrc");
    double factor = total_contraction(layer);

    part_data* nd = (part_data*)n->data;
    node* c = n->get_neighbor(atom("lyrCenter"));
    ipoint2 cdelta(-nd->cmx, -nd->cmy);

    if (c != nullptr) {
        part_data* cd = (part_data*)c->data;

        // Add central part image
        cd->make_mask(c, libtype, library);
        merge_mask_maps(mask, cd->mask, cdelta * factor, 1.0, libtype);
    }

    // Add other parts
    foreach_neighbor(n, srcname, iter) {
        part_data* nd = (part_data*)neighbor_node_data(iter);
        part_data_2* ed = (part_data_2*)neighbor_edge_data(iter);

        for (part_data_2::app_map_t::iterator aiter = ed->app.begin(); aiter != ed->app.end(); ++aiter) {
            node* varn = library->parts[nd->layer][aiter->first];
            part_data* varnd = (part_data*)varn->data;

            varnd->make_mask(varn, libtype, library);
            merge_mask_maps(mask, varnd->mask, (cdelta + ipoint2(ed->x, ed->y)) * factor, 1.0, libtype);
        }
    }
}

void merge_region_sets(set<ipoint2>& result, const set<ipoint2>& src, 
    const ipoint2& delta, const double& factor)
{
    for (set<ipoint2>::const_iterator siter = src.begin(); siter != src.end(); ++siter) {
        result.insert((*siter + delta) * factor);
    }
}

void part_data::make_region(node* n)
{
    if (!region.empty()) return;

    int srcname = atom("lyrSrc");
	double factor = total_contraction(layer);

	part_data* nd = (part_data*)n->data;
    ipoint2 cdelta(-nd->cmx, -nd->cmy);
    node* c = n->get_neighbor(atom("lyrCenter"));

	if (c != nullptr) {
		part_data* cd = (part_data*)c->data;

	    // Add central part image
		cd->make_region(c);
		merge_region_sets(region, cd->region, cdelta * factor, 1.0);
	}

    // Add other parts
    foreach_neighbor(n, srcname, iter) {
        part_data* nd = (part_data*)neighbor_node_data(iter);
        part_data_2* ed = (part_data_2*)neighbor_edge_data(iter);

        nd->make_region(neighbor_node(iter));
        merge_region_sets(region, nd->region, (cdelta + ipoint2(ed->x, ed->y)) * factor, 1.0);
    }
}

img part_data::get_image(node* n, unsigned libtype, part_lib* library)
{ 
    make_mask(n, libtype, library);
    return img(mask, 0.0, true, libtype == STRUCT_LIB_ATTR);
}

matrix<double> part_data::get_mask(node* n, part_lib* library)
{
    make_mask(n, 0, library);
    return rmatrix(mask, 0.0);
}

irectangle2 part_data::get_mask_box(node* n, part_lib* library)
{
    make_mask(n, 0, library);
    
    irectangle2 box;

    box.eat(0, 0);
    for (mask_t::iterator iter = mask.begin(); iter != mask.end(); ++iter) 
        box.eat(iter->first);
    return box;
}

ipoint2 part_data::get_mask(matrix<double>& m, node* n, part_lib* library)
{
    make_mask(n, 0, library);
    
    irectangle2 box = get_mask_box(n, library);

    m.resize(box.x_dim() + 1, box.y_dim() + 1);
    m.fill(0.0);
    for (mask_t::iterator iter = mask.begin(); iter != mask.end(); ++iter) {
        m(iter->first.x - box.ll.x, iter->first.y - box.ll.y) = iter->second;
    }
    return -box.ll;
}

//matrix<bool> part_data::get_region_matrix(node* n)
//{
//    matrix<bool> result;
//    get_region_matrix(result, n);
//    return result;
//}

void part_data::reset_mask()
{
    if (layer != 0) mask.clear();
}

void part_data::reset_region()
{
    if (layer != 0) region.clear();
}

matrix<bool> part_data::get_region(node* n)
{
    matrix<bool> result;
    get_region_matrix(result, n);
    return result;
}

ipoint2 part_data::get_region_matrix(matrix<bool>& m, node* n)
{
    irectangle2 box;
    set<ipoint2>::iterator iter;

    make_region(n);
    box.eat(0, 0);
    for (iter = region.begin(); iter != region.end(); ++iter) {
        box.eat(*iter);
    }
    m.resize(box.x_dim() + 1, box.y_dim() + 1);
    m.fill(false);
    for (iter = region.begin(); iter != region.end(); ++iter) {
        m(iter->x - box.ll.x, iter->y - box.ll.y) = true;
    }
    return -box.ll;
}

void part_data::get_region_set(set<ipoint2>& s, node* n)
{
    ipoint2 center(0, 0);

    s.clear();
    make_region(n);

    if (region.empty()) return;
    for (region_t::iterator iter = region.begin(); iter != region.end(); ++iter) 
        center += *iter;
    center.set(int_round((double)center.x/region.size()), int_round((double)center.y/region.size()));
    for (region_t::iterator iter = region.begin(); iter != region.end(); ++iter) 
        s.insert(*iter - center);
}

// svm_data
///////////////////////////////////////////////////////////////////////////////

void svm_data::write_to_stream(ostreamer& os) const
{
    cv_write_to_stream_f(os, mean);
    cv_write_to_stream_f(os, sigma);
}

void svm_data::read_from_stream(istreamer& is)
{
    if (is.get_version() >= 3.3) {
        cv_read_from_stream_f(is, mean);
        cv_read_from_stream_f(is, sigma);
    }

}


// vs_part_data
///////////////////////////////////////////////////////////////////////////////

string vs_part_data::name() const
{
    stringstream s;

    s << "part_" << layer << "_" << type;
    return s.str();
}

void vs_part_data::read_from_stream(istreamer& is)
{
    part_data::read_from_stream(is);
    pcad.read_from_stream(is);
    svmt.read_from_stream(is);
    svmh.read_from_stream(is);
}

void vs_part_data::write_to_stream(ostreamer& os)
{
    part_data::write_to_stream(os);
    pcad.write_to_stream(os);
    svmt.write_to_stream(os);
    svmh.write_to_stream(os);
}

void vs_part_data::write_to_storage(CvFileStorage* fs)
{
    if (svmd.svm != nullptr)
        svmd.svm->write(fs, name().c_str());
    if (svmt.svm != nullptr)
        svmt.svm->write(fs, (name() + "t").c_str());
    if (svmh.svm != nullptr)
        svmh.svm->write(fs, (name() + "h").c_str());
}

void vs_part_data::read_from_storage(CvFileStorage* fs)
{
    CvFileNode* fn = cvGetFileNodeByName(fs, nullptr, name().c_str());

    if (fn != nullptr) {
        svmd.svm.reset(new cv::SVM());
        svmd.svm->read(fs, fn);
    }

    fn = cvGetFileNodeByName(fs, nullptr, (name() + "t").c_str());

    if (fn != nullptr) {
        svmt.svm.reset(new cv::SVM());
        svmt.svm->read(fs, fn);
    }

    fn = cvGetFileNodeByName(fs, nullptr, (name() + "h").c_str());

    if (fn != nullptr) {
        svmh.svm.reset(new cv::SVM());
        svmh.svm->read(fs, fn);
    }

}



// rpart_data
///////////////////////////////////////////////////////////////////////////////

img rpart_data::square[2] = { img(10, 10, COL_TO_REAL(255, 255, 255), false), 
    img(10, 10, COL_TO_REAL(255, 255, 0), false) };

img rpart_data::get_image(node* n)
{
    static int factor = 10; 

    if (layer == 0) return square[type % 2];

    if (!image.empty()) return image;

    int to_center = atom("lyrCenter");
    int to_src = atom("lyrSrc");
    irectangle2 box;
    node* c = n->get_neighbor(to_center);
    rpart_data* cd = (rpart_data*)c->data;
    img cimage = cd->get_image(c);

    box.eat(-(int)cimage.width/2 - factor*cx, -(int)cimage.height/2 - factor*cy);
    box.eat((int)cimage.width/2 - factor*cx, (int)cimage.height/2 - factor*cy);

    foreach_neighbor(n, to_src, niter) {
        node* p = neighbor_node(niter);
        rpart_data* pd = (rpart_data*)p->data;
        part_data_2c* ed = (part_data_2c*)neighbor_edge_data(niter);
        img im = pd->get_image(p);

        box.eat(factor*ed->x - (int)im.width/2 - factor*cx, factor*ed->y - (int)im.height/2 - factor*cy);
        box.eat(factor*ed->x + (int)im.width/2 - factor*cx, factor*ed->y + (int)im.height/2 - factor*cy);
    }

    image.grayscale = false;
    image.resize(box.x_dim() + 2*factor, box.y_dim() + 2*factor, 0.0);

    int x0 = (box.x_dim() + 2*factor)/2, y0 = (box.y_dim() + 2*factor)/2;

    image.blt_central(cimage, (int)cimage.width/2, (int)cimage.height/2, x0 - factor*cx, y0 - factor*cy);
    foreach_neighbor(n, to_src, niter) {
        node* p = neighbor_node(niter);
        rpart_data* pd = (rpart_data*)p->data;
        part_data_2c* ed = (part_data_2c*)neighbor_edge_data(niter);
        img im = pd->get_image(p);

        image.blt_central(im, (int)im.width/2, (int)im.height/2, x0 + factor*ed->x - factor*cx, y0 + factor*ed->y - factor*cy);
    }
    return image;
} 

// spart_data & support
///////////////////////////////////////////////////////////////////////////////

// shape_model and shape_model_maker
//////////////////////////////////////

void shape_model::write_to_stream(ostreamer& os)
{
    os.write((int)m.size());
    for (map_t::iterator iter = m.begin(); iter != m.end(); ++iter) {
        os.write(iter->first);
        iter->second.write_to_stream(os);
    }
    //os.write(center);
    //os.write(rcount);
}

void shape_model::read_from_stream(istreamer& is)
{
    int size;

    is.read(size);
    for (int i = 0; i < size; ++i) {
        map_t::key_type key;
        map_t::mapped_type val;
        
        is.read(key);
        val.read_from_stream(is);
        m.insert(map_t::value_type(key, val));
    }
}

void shape_model::save_mathematica(string name) const
{
    ofstream os(name.c_str());

    os << '{';
    for (map_t::const_iterator iter = m.begin(); iter != m.end(); ++iter) {
        if (iter != m.begin()) os << ',';
        os << '{';
        os << '{';
        for (vector<int>::const_iterator kiter = iter->first.begin(); kiter != iter->first.end(); ++kiter) {
            if (kiter != iter->first.begin()) os << ',';
            os << *kiter;
        }
        os << '}';

        dpoint2 mean = iter->second.get_mean();
        matrix2x2 var = iter->second.get_variance();

        os << ',';
        os << '{';
        os << '{' << mean.x << ',' << mean.y << '}';
        os << ',';
        os << "{{" << var.a << ',' << var.b << "},{" << var.c << ',' << var.d << "}}";
        os << '}';

        os << '}';
    }
    os << '}';
    os.close();
}

//void shape_model_maker::update(const set<ipoint2>& pts, const ipoint2 c)
//{
//    if (pts.empty()) 
//        return;
//
//    double scale = 0.0;
//    ipoint2 avg(0, 0);
//
//    for (set<ipoint2>::const_iterator iter = pts.begin(); iter != pts.end(); ++iter) {
//        ipoint2 p = *iter;
//        
//        avg += p;
//    }
//    avg /= (int)pts.size();
//    for (set<ipoint2>::const_iterator iter = pts.begin(); iter != pts.end(); ++iter) {
//        ipoint2 p = *iter - avg;
//        
//        scale += p.x*p.x + p.y*p.y;
//    }
//    scale = sqrt(scale/pts.size());
//    for (set<ipoint2>::const_iterator iter = pts.begin(); iter != pts.end(); ++iter) {
//        ipoint2 p = *iter - avg; 
//        
//        p.x = (int)(100.0*p.x/scale);
//        p.y = (int)(100.0*p.y/scale);
//        ptset.insert(p + center);
//    }
//
//}

void fill_circle(matrix<int>& m, const ipoint2& c, int r, int value)
{
    int mini = max(0, c.x - r), maxi = min(c.x + r + 1, (int)m.width);
    int minj = max(0, c.y - r), maxj = min(c.y + r + 1, (int)m.height);
    int r2 = r*r;

    for (int i = mini; i < maxi; ++i)
        for (int j = minj; j < maxj; ++j) {
            //if ((i - c.x)*(i - c.x) + (j - c.y)*(j - c.y) <= r2) 
            m(i, j) = value;
        }
}

int circle_intersection_size(const matrix<int>& m, const ipoint2& c, int r)
{
    int result = 0;
    int mini = max(0, c.x - r), maxi = min(c.x + r + 1, (int)m.width);
    int minj = max(0, c.y - r), maxj = min(c.y + r + 1, (int)m.height);
    int r2 = r*r;

    for (int i = mini; i < maxi; ++i)
        for (int j = minj; j < maxj; ++j) {
            //if ((i - c.x)*(i - c.x) + (j - c.y)*(j - c.y) <= r2) 
            if (m(i, j) != 0) ++result;
        }
    return result;
}

shape_model shape_model_maker::get_model() const
{
    shape_model result;

    for (stat_t::const_iterator iter = stat.begin(); iter != stat.end(); ++iter) {
        normal_distribution2 dist;

        get_normal_distribution2(dist, iter->second);
        result.m[iter->first] = dist;
    }
    return result;
}

double shape_model::convolve(const map<ip2_vector, ipoint2>& v) const
{
    return 1.0;
}

int shape_model::distance2(const vector<int>& v, const ipoint2& p) const
{
    map_t::const_iterator iter = m.find(v);

    if (iter == m.end()) return INT_MAX;
    else {
        const dpoint2& dp = iter->second.get_mean();   
        return p.distance2((int)dp.x, (int)dp.y);
    }
}

// part_data_2a
///////////////////////////////////////////////////////////////////////////////

void read_path_map(path_map_t& pm, istreamer& is)
{
    int size;

    is.read(size);
    for (int i = 0; i < size; ++i) {
        vector<int> key;
        
        is.read(key);

        hpoint_t& val = pm[key];

        is.read(val.p);
        is.read(val.h);
    }
}

void write_path_map(ostreamer& os, const path_map_t& pm)
{
    os.write((int)pm.size());
    for (path_map_t::const_iterator miter = pm.begin(); miter != pm.end(); ++miter) {
        os.write(miter->first); // vector<int>
        os.write(miter->second.p); // ipoint2
        os.write(miter->second.h); // vector<double>
    }
}

void part_data_2a::read_from_stream(istreamer& is)
{ 
    edge_data::read_from_stream(is); 

	if (is.get_version() >= 3.4) is.read(app); 
	else {
		map<int, double> app1;

		is.read(app1);
		app.clear();
		for (auto iter = app1.begin(); iter != app1.end(); ++iter) {
			app.insert(app_map_t::value_type(iter->first, app_map_t::mapped_type(vector<int>(), iter->second)));
		}
	}

    // Read 'geo' map
    if (is.get_version() <= 3.0) {
        path_map_t pm;

        read_path_map(pm, is);
        if (!pm.empty()) 
            geo[app.begin()->first] = pm;
    } else {
        int size;

        is.read(size);
        for (int i = 0; i < size; ++i) {
            int key;

            is.read(key);

            path_map_t& pm = geo[key];

            read_path_map(pm, is);
        }
    }
}

void part_data_2a::write_to_stream(ostreamer& os)
{
    edge_data::write_to_stream(os); 
    os.write(app); 
    os.write((int)geo.size());
    for (geo_map_t::const_iterator iter = geo.begin(); iter != geo.end(); ++iter) {
        os.write(iter->first);
        write_path_map(os, iter->second);
    }
}


// part_lib
///////////////////////////////////////////////////////////////////////////////


// static fields definition and initialization

double part_lib::similarity_threshold = 1E-6;

// constructors & destructors

part_lib::part_lib(unsigned libtype) : updated(false), attr(libtype), layer_count(0), parts(MAX_LAYER_NUMBER)
#ifdef OPENCL
	, ocl_made(false)
#endif
{  
    memset(layer_data, 0, MAX_LAYER_NUMBER * sizeof(streamable*));
	memset(contractions, 0, MAX_LAYER_NUMBER * sizeof(double));

}

part_lib::~part_lib()
{
    for (int i = 0; i < layer_count; ++i) {
        if (layer_data[i] != nullptr) delete layer_data[i];
    }
#ifdef OPENCL
	for (int i = 0; i < ocl_parts.size(); ++i) {
		if (ocl_parts[i].first != nullptr) delete[] ocl_parts[i].first;
	}
	for (int i = 0; i < ocl_edges.size(); ++i) {
		if (ocl_edges[i].first != nullptr) delete[] ocl_edges[i].first;
	}
	for (int i = 0; i < ocl_apps.size(); ++i) {
		if (ocl_apps[i].first != nullptr) delete[] ocl_apps[i].first;
	}	
#endif
}

// member functions definitions

void part_lib::get_structure_c(node* n, vector<part_str>& str)
{
    int edge_name = atom("lyrPrev").get_index();
    part_data* nd = (part_data*)n->data;

    str.clear();
    foreach_neighbor(n, edge_name, j) {
        part_data* pd = (part_data*)neighbor_node_data(j);
        part_data_2c* ed = (part_data_2c*)neighbor_edge_data(j);

        str.push_back(part_str(pd->layer - nd->layer, pd->type, ed->x, ed->y, nullptr));
    }
}

void part_lib::get_structure(node* n, vector<part_str>& str)
{
    int src_name = atom("lyrSrc").get_index();
    part_data* nd = (part_data*)n->data;

    str.clear();
    node* c = n->get_neighbor(atom("lyrCenter"));
    str.push_back(part_str(-1, ((part_data*)c->data)->type, 0, 0, nullptr));
    foreach_neighbor(n, src_name, j) {
        part_data* pd = (part_data*)neighbor_node_data(j);
        part_data_2* ed = (part_data_2*)neighbor_edge_data(j);

        str.push_back(part_str(pd->layer - nd->layer, pd->type, ed->x, ed->y, ed));
    }
}

void part_lib::get_structure(node* n, int name, vector<part_str>& str)
{
    part_data* nd = (part_data*)n->data;

    str.clear();
    foreach_neighbor(n, name, j) {
        part_data* pd = (part_data*)neighbor_node_data(j);
        part_data_2* ed = (part_data_2*)neighbor_edge_data(j);

        if (ed == nullptr) 
            str.push_back(part_str(pd->layer - nd->layer, pd->type, 0, 0, nullptr));
        else
            str.push_back(part_str(pd->layer - nd->layer, pd->type, ed->x, ed->y, ed));
        
    }
}


void part_lib::get_structure(const vector<iipair>& str, const vector<iipair>& coo,
    vector<part_str>& result)
{
    result.clear();
    result.push_back(part_str(str[0].first, str[0].second, 0, 0, nullptr));

    size_t max = min(str.size(), coo.size());
    for (size_t i = 1; i < max; ++i)
        result.push_back(part_str(str[i].first, str[i].second, coo[i].first, coo[i].second, nullptr));
}

// Returns true if left \subsetneq right (as parts)
// Center left[0] must be equal to right[0].
bool part_lib::is_subpart(const vector<part_str>& left, const vector<part_str>& right)
{
    if (left.size() == 0 || right.size() == 0) return false;
    if (left[0].type != right[0].type) return false;

    return is_subset(left.begin() + 1, left.end(), right.begin() + 1, right.end(), true);
}

// Returns true if left \subsetneq right (as parts)
// Center left[0] must be equal to right[0].
// difference is right \ left
bool part_lib::is_subpart(const vector<part_str>& left, const vector<part_str>& right,
    list<part_str>& diff)
{
    if (left.size() == 0 || right.size() == 0) return false;
    if (left[0].type != right[0].type) return false;

    bool subset = range_difference(left.begin() + 1, left.end(), right.begin() + 1, right.end(), diff);
    return subset && !diff.empty();
}

// Returns true if left \subsetneq right (as parts)
bool part_lib::is_subpart(node* leftn, const vector<part_str>& right)
{
    vector<part_str> left;

    get_structure(leftn, left);
    return is_subpart(left, right);
}

// Returns true if left \subsetneq right (as parts)
bool part_lib::is_subpart(node* leftn, node* rightn)
{
    vector<part_str> left, right;

    get_structure(leftn, left);
    get_structure(rightn, right);
    return is_subpart(left, right);
}

void part_lib::extend_subparts(int lyr, node* n, niterator begin, niterator end)
{
    vector<part_str> nstr, istr, forb;
    list<part_str> diff, diff2;
    list<part_str>::iterator j;
    int forbindex = atom("lyrForbidden").get_index();

    get_structure(n, nstr);
    for (niterator i = begin; i != end; ++i) {
        get_structure(*i, istr);
        if (is_subpart(nstr, istr, diff)) {
            // add new edge
            get_structure(n, forbindex, forb);
            range_difference(forb.begin(), forb.end(), diff.begin(), diff.end(), diff2);
            
            for (j = diff2.begin(); j != diff2.end(); ++j) {
                add_edge_2(n, 
                    parts[lyr - 1 + j->ldelta][j->type], 
                    new part_data_2(j->x, j->y, j->data->distr, j->type, -1), forbindex);
            }
        }
    }
}

// lyr = layer nodes we want to check
// str = structure of the part (indices to previous level parts); str[0] is center
// coo = coordinates of parts
// tol = tolerance for search
// returns the type of the part; < 0 if the part does not exist or data is invalid
int part_lib::part_exists_c(const vector<node*>& lyr, const vector<iipair>& str, const vector<iipair>& coo, 
    int tol /* = 0 */)
{
    if (str.size() != coo.size()) return -1;

    vector<node*>::const_iterator iter;
    part_str_comparer comp(tol);
    vector<part_str> srcstr;

    for (size_t i = 0; i < str.size(); ++i) 
        srcstr.push_back(part_str(str[i].first, str[i].second, coo[i].first, coo[i].second, nullptr));

    for (iter = lyr.begin(); iter != lyr.end(); ++iter) {
        vector<part_str> partstr;

        get_structure_c(*iter, partstr);
        if (predicate_equal(srcstr.begin(), srcstr.end(), partstr.begin(), partstr.end(), comp))
            return ((part_data*)(*iter)->data)->type;
    }
    return -1;
}

// lyr = layer nodes we want to check
// str = structure of the part (indices to previous level parts) -- (lyr. delta, type); str[0] is center
// coo = coordinates of parts
// returns the type of the part; < 0 if the part does not exist or data is invalid
int part_lib::part_exists(const vector<node*>& lyr, const vector<iipair>& str, const vector<iipair>& coo)
{
    if (str.size() != coo.size()) 
        throw new_libhop_exception("Invalid data in part_lib::part_exists.");

    vector<part_str> srcstr;

    for (size_t i = 0; i < str.size(); ++i) 
        srcstr.push_back(part_str(str[i].first, str[i].second, coo[i].first, coo[i].second, nullptr));

    for (vector<node*>::const_iterator iter = lyr.begin(); iter != lyr.end(); ++iter) {
        vector<part_str> partstr;

        get_structure(*iter, partstr);
        if (ranges_equal(srcstr.begin(), srcstr.end(), partstr.begin(), partstr.end()))
            return ((part_data*)(*iter)->data)->type;
    }
    return -1;
}

irectangle2 part_lib::subpart_box(node* p)
{
    irectangle2 box;
    int srcname = atom("lyrSrc");

    box.eat(0, 0);
    foreach_neighbor(p, srcname, siter) {
        part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
	    box.eat(ed->x, ed->y);
    }
    return box;
}

irectangle2 part_lib::subpart_positions(vector<ipoint2>& result, node* p)
{
    int srcname = atom("lyrSrc");
    irectangle2 box;

    result.clear();
    result.push_back(ipoint2::zero); // center
    box.eat(0, 0);
    foreach_neighbor(p, srcname, siter) {
        part_data_2* ed = (part_data_2*)neighbor_edge_data(siter);
        result.push_back(ipoint2(ed->x, ed->y));
        box.eat(ed->x, ed->y);
    }
    return box;
}


// lyr = layer we want to add part to
// pd = pointer to part_data structure (node data)
// str = structure of the part (indices to previous level parts) -- (lyr. diff, type); str[0] is center
//      Note: lyr. difference < 0
// coo = coordinates of parts
// distr = distrubution around the part
// tol = tolerance to check whether a part already exists
// the last three parametres are ignored when lay = 1
// returns type (index) of the part if it succedes -type-1 of the part if it already exists
int part_lib::add_part(int lyr, part_data* pd, const vector<iipair>& str, const vector<iipair>& coo, 
    const vector<matrix<double>*>& distr, double eqthresh /* = -1.0 */, double contraction /* = 1.0 */, 
    int tol /* = 0 */, int eqtol /* = 0 */)
{
    if (str.size() == 0 && lyr != 1) throw;
    if (contraction <= 0.0) throw new_libhop_exception("Contraction exception"); //Contraction not defined

    int part_index;

    if (lyr == 1) {
        node* n = add_node(pd);
    
        if (lyr > layer_count) ++layer_count;
        parts[lyr - 1].push_back(n);
        part_index = (int)parts[lyr - 1].size() - 1;
    } else {
        vector<node*>& lyrm1 = parts[lyr - 1 + str[0].first];

        int index = str[0].second;
        part_data* cdata = (part_data*)lyrm1[index]->data;
        vector<part_data*> pdata;
        iipair massc = center_of_mass(coo.begin(), coo.end());
        vector<iipair> cooc(coo.begin(), coo.end());
        vector<iipair> coon(coo.begin(), coo.end());
        iipair coo0(coon[0]);
        int extype;

        for (size_t i = 0; i < cooc.size(); ++i) cooc[i] -= massc;
        for (size_t i = 0; i < coon.size(); ++i) coon[i] -= coo0;

        // check whether a part already exists
        if ((extype = part_exists(parts[lyr - 1], str, coon)) >= 0)
            return -extype - 1;

        pd->cmx = massc.first - coo0.first;
        pd->cmy = massc.second - coo0.second;
        pd->layer = lyr - 1;
        pd->type = (int)parts[lyr - 1].size();

        // add node 
        node* n = add_node(pd);
    
        if (lyr > layer_count) ++layer_count;
        parts[lyr - 1].push_back(n);
        part_index = (int)parts[lyr - 1].size() - 1;
        
        // add an edge to the center
        if (index < (int)lyrm1.size()) {
            add_edge_2(n, lyrm1[index], new part_data_2a(index), nullptr, atom("lyrCenter"), atom("lyrCenterBack"));
            add_edge_2(n, lyrm1[index], new part_data_2c(cooc[0]), atom("lyrPrev"));
        }

        // add edges to other subparts
        pdata.push_back(cdata);
        for (size_t i = 1; i < str.size(); ++i) {
            vector<node*>& lyrm1 = parts[lyr - 1 + str[i].first];
            index = str[i].second;
            if (index < (int)lyrm1.size()) {
                add_edge_2(n, lyrm1[index], new part_data_2(coon[i], *distr[i], index, i), atom("lyrSrc"));
                add_edge_2(n, lyrm1[index], new part_data_2c(cooc[i]), atom("lyrPrev"));
                pdata.push_back((part_data*)lyrm1[index]->data);
            }
        }
        extend_subparts(lyr, n, parts[lyr - 1].begin(), parts[lyr - 1].end());
        //make_image(lyr, part_index, n, pdata, coo, massc, contraction, eqthresh, eqtol);
    }
    if (!updated || contractions[lyr - 1] <= 0.0) {
        contractions[lyr - 1] = contraction;
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);
        updated = true;
    } else if (contractions[lyr - 1] != contraction) 
		throw new_libhop_exception("Contraction exception"); // Contractions mismatch
    return part_index;
}

int part_lib::add_object_part(int lyr, part_data* pd, const vector<cluster_data_t>& objdata,
        const vector<matrix<double>*>& distr, double contraction /* = 1.0 */)
{
    typedef vector<pair<vector<ipoint2>, ipoint2> > objdata_t;

    objdata_t objdata2;

    for (vector<cluster_data_t>::const_iterator iter = objdata.begin(); iter != objdata.end(); ++iter) {
        objdata_t::value_type item(iter->types, iter->pos);

        objdata2.push_back(item);
    }
    
    int result = add_object_part(lyr, pd, objdata2, distr, contraction);

    if (result < 0) 
        return result;

    // Add geometry data to sub-parts
    node* p = parts[lyr - 1][result];
    int ename = atom("lyrSrc");

    foreach_neighbor(p, ename, niter) {
        part_data_2* ed = (part_data_2*)neighbor_edge_data(niter);

        ed->geo = objdata[ed->index].geo;
    }
    return result;
}

// Add part representing an object to the library.
// - 'lyr' layer (index of the layer + 1)
// - 'pd' part_data object; cmx and cmy are filled automatically to (0, 0)
// - 'objdata': vector of subparts; vector<ipoint2> is a set of equivalent subparts and
//   ipoint2 is a coordinate of the subpart *relative to the center*.
// - 'distr' a vector of distributions for each part.
// - 'contraction' is a contraction of the layer (default value is 1.0).
// Returns index of the part; -1 if an error occurred.
int part_lib::add_object_part(int lyr, part_data* pd, const vector<pair<vector<ipoint2>, ipoint2> >& objdata,
        const vector<matrix<double>*>& distr, double contraction /* = 1.0 */)
{
    typedef pair<vector<ipoint2>, ipoint2> cluster_data_t;
    typedef vector<cluster_data_t> obj_data_t;

    if (objdata.empty()) return -1;

    int layer = lyr - 1;

    if (!updated) {
        contractions[layer] = contraction;
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);
        updated = true;
    } else if (contractions[layer] != contraction) throw new_libhop_exception("Contraction exception");

    dpoint2 dmassc(0.0, 0.0);
 
    for (obj_data_t::const_iterator iter = objdata.begin(); iter != objdata.end(); ++iter) {
        double c = this->contraction(layer + iter->first.front().x + 2, layer);
        dmassc += dpoint2(iter->second.x, iter->second.y)/c;
    }
    dmassc = dmassc/(double)objdata.size();

    ipoint2 massc(int_round(dmassc.x), int_round(dmassc.y));

    pd->cmx = 0;
    pd->cmx = 0;
    pd->layer = layer;
    pd->type = (int)parts[layer].size();

    node* n = add_node(pd, OBJ_PART_ATTR);

    if (layer >= layer_count) ++layer_count;
    parts[layer].push_back(n);

    for (size_t i = 0; i < objdata.size(); ++i) {
        const cluster_data_t& cl = objdata[i];
        double c = this->contraction(layer + cl.first.front().x + 2, layer);
        ipoint2 coo = ipoint2((int)(cl.second.x/c), (int)(cl.second.y/c)) - massc;
        const ipoint2& cld = cl.first[0];
        node* nn = parts[layer + cld.x][cld.y];
        part_data_2* pd2 = new part_data_2(coo.x, coo.y, *distr[i], cld.y, i);
        
        for (size_t j = 1; j < cl.first.size(); ++j) {
            const ipoint2& cldj = cl.first[j];
            node* nnn = parts[layer + cldj.x][cldj.y];

			pd2->app.insert(pair<int, pair<int, double> >(cldj.y, pair<int, double>(i, 1.0)));
            if (cldj.x == -1) add_edge_2(nnn, n, new part_data_2c(coo.x, coo.y), atom("lyrCenterBack"));
        }
        add_edge_2(n, nn, pd2, atom("lyrSrc"));
        if (cld.x == -1) add_edge_2(nn, n, new part_data_2c(coo.x, coo.y), atom("lyrCenterBack"));
    }
    if (!updated || contractions[lyr - 1] <= 0.0) {
        contractions[lyr - 1] = contraction;
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);
        updated = true;
    } else if (contractions[lyr - 1] != contraction) 
		throw new_libhop_exception("Contraction exception"); // Contractions mismatch
    return pd->type;
}

string part_lib::get_category(node* p)
{
    node* cp = p->get_neighbor(atom("lyrCenterBack"));

    if (cp == nullptr || !cp->is_attr_set(CATEGORY_NODE_ATTR)) return string();
    return ((cpart_data*)cp->data)->name;
}


// See overloaded variant below.
int part_lib::find_equivalent_object_part(int lyr, const vector<cluster_data_t>& objdata, int masksize)
{
    typedef vector<pair<vector<ipoint2>, ipoint2> > objdata_t;

    objdata_t objdata2;

    for (vector<cluster_data_t>::const_iterator iter = objdata.begin(); iter != objdata.end(); ++iter) {
        objdata_t::value_type item(iter->types, iter->pos);

        objdata2.push_back(item);
    }
    return find_equivalent_object_part(lyr, objdata2, masksize);
}

// Returns the index of the equivalent object part or -1 if all parts are sufficiently 
// different from objdata
int part_lib::find_equivalent_object_part(int lyr, const vector<pair<vector<ipoint2>, ipoint2> >& objdata, 
    int masksize)
{
    typedef vector<pair<vector<ipoint2>, ipoint2> > objdata_t;
    typedef vector<pair<set<int>, ipoint2> > local_objdata_t;
    typedef set<int> iset;
    typedef matrix<iset*> objdata_matrix_t;

    const int msize = 200;

    int srcedge = atom("lyrSrc");
    vector<node*>& lyrparts = parts[lyr];

    // Calc mass center
    dpoint2 dmassc(0.0, 0.0);

    for (objdata_t::const_iterator iter = objdata.begin(); iter != objdata.end(); ++iter) {
        dmassc += dpoint2(iter->second.x, iter->second.y);
    }
    dmassc = dmassc/(double)objdata.size();

    // Make local structure
    local_objdata_t lobjdata;

    for (objdata_t::const_iterator iter = objdata.begin(); iter != objdata.end(); ++iter) {
        set<int> cluster;

        for (vector<ipoint2>::const_iterator viter = iter->first.begin(); viter != iter->first.end(); ++viter) {
            if (viter->x == -1) 
                cluster.insert(viter->y);
        }
        lobjdata.push_back(local_objdata_t::value_type(cluster, 
            ipoint2((int)(iter->second.x - dmassc.x), (int)(iter->second.y - dmassc.y))));
    }

    // Make matrix structure
    objdata_matrix_t obm(msize, msize, nullptr);

    for (local_objdata_t::iterator iter = lobjdata.begin(); iter != lobjdata.end(); ++iter) {
        obm.set_region_c(iter->second.x + msize/2, iter->second.y + msize/2, masksize/2, masksize/2, &iter->first);
    }

    // Check with parts
    for (vector<node*>::iterator piter = lyrparts.begin(); piter != lyrparts.end(); ++piter) {
        node* p = *piter;
        bool similar = true;
        
        if (!p->is_attr_set(OBJ_PART_ATTR)) continue; // or break since this is probably not the
                                                      // object layer
        foreach_neighbor (p, srcedge, spiter) {
            part_data_2* ed = (part_data_2*)neighbor_edge_data(spiter);
            if (obm(ed->x + msize/2, ed->y + msize/2) == nullptr) {
                similar = false;
                break;
            } else {
                set<int> appset;

                get_key_set(appset, ed->app);
                if (intersection_size(*obm(ed->x + msize/2, ed->y + msize/2), appset) < appset.size()/2.0) {
                    similar = false;
                    break;
                }
            }
        }
        if (similar) {
            return ((lib_data*)p->data)->type;
        }
    }
    return -1;
}

//void part_lib::make_image(int lyr, int index, node* n, vector<part_data*>& pdata, 
//    const vector<iipair>& coorig, iipair massc, double contraction, double eqthresh, int eqtol)
//{
//    if (lyr <= 1) return;
//
//    //contraction *= 1.0 + (lyr - 2) * 0.1;   // "magic" for nicer parts...
//
//    vector<iipair> coo;
//    expand_coordinates(coorig, coo, contraction);
//    expand_coordinates(massc, contraction);
//
//    part_data* nd = (part_data*)n->data;
//    int minx = massc.first, miny = massc.second;
//    int maxx = minx, maxy = miny;
//    size_t i;
//    
//    for (i = 0; i < pdata.size(); ++i) {
//        const iipair& p = coo[i];
//        if (p.first > maxx) maxx = p.first; else if (p.first < minx) minx = p.first;
//        if (p.second > maxy) maxy = p.second; else if (p.second < miny) miny = p.second;
//    }
//    int dimxl1 = sizes[lyr - 2].first;
//    int dimyl1 = sizes[lyr - 2].second;
//    int dimx = maxx - minx + 2*dimxl1;
//    int dimy = maxy - miny + 2*dimyl1;
//    int offx = minx - dimxl1, offy = miny - dimyl1;
//    std::logical_or<bool> logical_or;
//    
//    nd->cx = massc.first - offx; 
//    nd->cy = massc.second - offy;
//
//    //if (lyr < 5) {
//        //nd->mask.resize(dimx, dimy, 0);
//        nd->region.resize(dimx, dimy);
//    
//        for (i = 0; i < pdata.size(); ++i) {
//            const iipair& p = coo[i];
//            part_data* d = pdata[i];
//
//            //nd->mask.blt_central_max(d->mask, d->cx, d->cy, p.first - offx, p.second - offy);
//            nd->region.blt_central(d->region, d->cx, d->cy, p.first - offx, p.second - offy, logical_or);
//        }
//    //}
//    if (eqthresh < 0.0) {
//        nd->type = index;
//    } else {
//        int newtype = find_type_equal_by_region(parts[lyr - 1].begin(), parts[lyr - 1].end() - 1, nd, eqthresh, eqtol, eqtol);
//        nd->type = (newtype >= 0) ? newtype : index; 
//    }
//    //nd->mask.save_normalized("c:\\TEMP\\temp.bmp");
//}

int part_lib::add_category(int lyr, const vector<int>& members, const string& name)
{
    if (lyr < 2) return -1;

    if (lyr > layer_count) ++layer_count;

    int z = lyr - 1;
    cpart_data* pd = new cpart_data(name == "" ? (string("cat_") + parts[z].size()) : name);
    node* n = add_node(pd, CATEGORY_NODE_ATTR);
    int to_center = atom("lyrCenter").get_index();
    int to_center_back = atom("lyrCenterBack").get_index();
    int to_prev = atom("lyrPrev").get_index();

    pd->type = (int)parts[z].size();
    pd->layer = z;
    parts[z].push_back(n);

    for (vector<int>::const_iterator iter = members.begin(); iter != members.end(); ++iter) {
        int i = *iter;
        add_edge(n, parts[z - 1][i], to_center, to_center_back);
        add_edge_2(n, parts[z - 1][i], new part_data_2c(iipair(0, 0)), to_prev);
    }
    return pd->type;
}

int part_lib::merge_part_types(int lyr, double thresh, int tol)
{
    vector<node*> layer = parts[lyr - 1];
    int result = 0;
    
    for (vector<node*>::iterator iter = layer.begin(); iter != layer.end(); ++iter) {
        part_data* nd = (part_data*)(*iter)->data;
        int newtype = find_type_equal_by_region(layer.begin(), iter, *iter, thresh, tol, tol);

        if (newtype >= 0) { 
            cout << nd->type << " -> " << newtype << "  ";
            nd->type = newtype; 
            ++result; 
        }
    }
    return result;
}

void rotate_region(const set<iipair>& reg, set<iipair>& result, double deg)
{
	double c = ::cos(deg);
	double s = ::sin(deg);
	result.clear();
	for (set<iipair>::const_iterator iter = reg.begin(); iter != reg.end(); ++iter) {
		result.insert(iipair((int)(c*iter->first - s*iter->second), (int)(s*iter->first + c*iter->second)));
	}
}

void part_lib::rotational_equivalence(list<set<int> >& result, int lyr, int degnum, int thresh, int tolx, int toly)
{
	vector<node*>& layer = parts[lyr - 1];
	set<int> processed;
	int total = (int)layer.size();
	int count = 0;

    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);
	result.clear();
    for (vector<node*>::iterator iter = layer.begin(); iter != layer.end(); ++iter) {
        part_data* nd = (part_data*)(*iter)->data;
		int type = nd->type;

		cout << count++ << '/' << total << ' ';

		if (processed.find(type) != processed.end()) continue;
		
		result.push_back(set<int>());

		set<int>& orbit = result.back();

		for (vector<node*>::iterator iter2 = iter; iter2 != layer.end(); ++iter2) {
	        part_data* nd2 = (part_data*)(*iter2)->data;

			if (nd2->type != type) continue;

 			//cout << ';';

			set<iipair> cregionnd;
			set<iipair> cregionndr;
			double ddeg;
            matrix<bool> region;
		
			orbit.insert(type);
			processed.insert(type);
            nd2->get_region_matrix(region, *iter2);
			centralized_region_coordinates(cregionnd, region);
			for (int deg = 1; deg < degnum; ++deg) {
				int ftype;
				//cout << '.';
				ddeg = 2.0*M_PI*deg/(double)degnum;
				rotate_region(cregionnd, cregionndr, ddeg);
				ftype = find_type_equal_by_region(layer.begin(), layer.end(), cregionndr, thresh, tolx, toly);
				orbit.insert(ftype);
				processed.insert(ftype);
			}
		}
	}
}

struct ispair {
	double d;
	node* n1;
	node* n2;

	ispair(double vd, node* vn1, node* vn2) : d(vd), n1(vn1), n2(vn2) { }

	bool operator<(const ispair& q) const { return d < q.d; }
	bool operator>(const ispair& q) const { return d > q.d; }
	bool operator==(const ispair& q) const { return d == q.d; }
};

void centralized_region_coordinates(set<wpoint>& coo, const img& region)
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
int part_lib::make_similarity_clusters(const vector<node*> parts, vector<set<node*> > clusters, 
	double thresh)
{
	
	int size = (int)parts.size();
	matrix<double> similarities(size, size, 1000.0);
	priority_queue<ispair> queue;

	for (int i = 0; i < parts.size(); ++i) {
		node* n1 = parts[i];
		part_data* d1 = (part_data*)n1->data;
		set<wpoint> s1;

		centralized_region_coordinates(s1, d1->mask);
		for (int j = i + 1; j < parts.size(); ++j) {
			node* n2 = parts[j];
			part_data* d2 = (part_data*)n2->data;
			set<wpoint> s2;
			double wt;

			centralized_region_coordinates(s2, d2->mask);
			wt = intersection_weight(s1, s2);

			if (wt > thresh) {
				similarities(i, j) = similarities(j, i) = wd;
				queue.push(ispair(wt, n1, n2));
			}
		}
	}

}*/

void copy_region(const set<iipair>& src, vector<set<iipair> >& dest, int x0, int x1, int y0, int y1)
{
    dest.resize((x1 - x0 + 1)*(y1 - y0 + 1));
    int i = 0;
    set<iipair>::const_iterator iter;

    for (int x = x0; x <= x1; ++x) {
        for (int y = y0; y <= y1; ++y) {
            set<iipair>& s = dest[i++];
            iipair d(x, y);

            for (iter = src.begin(); iter != src.end(); ++iter) {
                s.insert(*iter + d);
            }
        }
    }
}


// Implementation issue: converting from set<ipoint> to matrix<bool> to set<iipair> !
int part_lib::find_type_equal_by_region(vector<node*>::iterator itera, vector<node*>::iterator iterb, 
                                        node* n, double thresh, int dx /* = 0 */, int dy /* = 0 */)
{
    set<iipair> cregionnd, cregion;
    vector<set<iipair> > cregionndv;
    vector<set<iipair> >::iterator riter;
    int isize, rsize, r2size;
    double quot, quot1, quot2, bestquot = -1.0;
    int result = -1;
    part_data* nd = (part_data*)n->data;
    matrix<bool> region;

    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);
    nd->get_region_matrix(region, n);
    centralized_region_coordinates(cregionnd, region);
    copy_region(cregionnd, cregionndv, -dx, dx, -dy, dy);
    rsize = (int)cregionnd.size();
    for (vector<node*>::iterator iter = itera; iter != iterb; ++iter) {
        part_data* ond = (part_data*)(*iter)->data;

        ond->get_region_matrix(region, *iter);
        centralized_region_coordinates(cregion, region);
        r2size = (int)cregion.size();
        for (riter = cregionndv.begin(); riter != cregionndv.end(); ++riter) {
            isize = intersection_size(cregion, *riter);
            if ((quot1 = (double)isize/rsize) >= thresh && (quot2 = (double)isize/r2size) >= thresh) { 
                quot = min(quot1, quot2);
                if (quot > bestquot) { bestquot = quot; result = ond->type; }
            }
        }
    } 
    return result;
}

// Implementation issue: converting from set<ipoint> to matrix<bool> to set<iipair> !
int part_lib::find_type_equal_by_region(vector<node*>::iterator itera, vector<node*>::iterator iterb, 
                                        set<iipair>& cregionnd, double thresh, int dx /* = 0 */, int dy /* = 0 */)
{
    set<iipair> cregion;
    vector<set<iipair> > cregionndv;
    vector<set<iipair> >::iterator riter;
    int isize, rsize, r2size;
    double quot, quot1, quot2, bestquot = -1.0;
    int result = -1;
    matrix<bool> region;

    copy_region(cregionnd, cregionndv, -dx, dx, -dy, dy);
    rsize = (int)cregionnd.size();

    for (vector<node*>::iterator iter = itera; iter != iterb; ++iter) {
        part_data* ond = (part_data*)(*iter)->data;

        ond->get_region_matrix(region, *iter);
        centralized_region_coordinates(cregion, region);
        r2size = (int)cregion.size();
        for (riter = cregionndv.begin(); riter != cregionndv.end(); ++riter) {
            isize = intersection_size(cregion, *riter);
            if ((quot1 = (double)isize/rsize) >= thresh && (quot2 = (double)isize/r2size) >= thresh) { 
                quot = min(quot1, quot2);
                if (quot > bestquot) { bestquot = quot; result = ond->type; }
            }
        }
    } 
    return result;
}

void part_lib::delete_parts(int lyr, vector<int> v)
{
    set<node*> todelete;
    vector<node*>& lyparts = parts[lyr];
    vector<node*>::iterator iter;

    sort(v.begin(), v.end());
    for (int i = (int)v.size() - 1; i >= 0; --i) {
        int index = v[i];

        if (index >= (int)lyparts.size()) continue;

        iter = lyparts.begin() + index;
        todelete.insert(*iter);
        lyparts.erase(iter);
    }
    for (int i = (int)lyparts.size() - 1; i >= 0; --i) {
        part_data* pd = (part_data*)lyparts[i]->data;

        pd->type = i;
    }
    delete_nodes(todelete);
}

void part_lib::delete_parts(int lyr)
{
    vector<node*>& lyparts = parts[lyr];

    delete_nodes(set<node*>(lyparts.begin(), lyparts.end()));
    lyparts.clear();
}

void part_lib::delete_parts_geq(int lyr)
{
    set<node*> to_delete;

    for (int i = lyr; i < layer_count; ++i) {
        to_delete.insert(parts[i].begin(), parts[i].end());
        parts[i].clear();
    }
    delete_nodes(to_delete);
}

void part_lib::delete_layers_geq(int lyr)
{
    delete_parts_geq(lyr);
    layer_count = lyr;
}


class recurse_collector {
public:
    set<node*> nodes;

    recurse_collector() : nodes() { }
    bool operator()(node* n, node* nn) { nodes.insert(n); nodes.insert(nn); return true; }
};

void part_lib::delete_unused_parts()
{
    if (layer_count < 2) return;

    set<node*> result;
    recurse_collector coll;
    
    recurse(parts[layer_count - 1], atom("lyrPrev").get_index(), coll, result);
    for (int l = 1; l < layer_count - 1; ++l) {
        vector<int> to_keep;

        for (set<node*>::iterator iter = coll.nodes.begin(); iter != coll.nodes.end(); ++iter) {
            lib_data* nd = (lib_data*)(*iter)->data;
            if (nd->layer == l) to_keep.push_back(nd->type);
        }
        keep_parts(l, to_keep);
    }
}

void part_lib::keep_parts(int lyr, const vector<int>& v)
{
    int size = (int)parts[lyr].size();
    vector<int> todelete;
    
    for (int i = 0; i < size; ++i) {
        if (find(v.begin(), v.end(), i) == v.end()) todelete.push_back(i);
    }
    delete_parts(lyr, todelete);
}

void part_lib::keep_clusters(int lyr, const set<int>& v)
{
    int edge = atom("lyrSimilar");
    int redge = atom("lyrSimRoot");
    set<int> tokeep;
    set<int> nonroots;

    for (set<int>::const_iterator viter = v.begin(); viter != v.end(); ++viter) {
        node* p = parts[lyr][*viter];
        lib_data* pd = (lib_data*)p->data;
                
        tokeep.insert(pd->type);
        foreach_neighbor(p, edge, niter) {
            lib_data* npd = (lib_data*)neighbor_node_data(niter);
            tokeep.insert(npd->type);
            if (neighbor_node(niter) != p) 
                nonroots.insert(npd->type);
        }
    }
    // Delete edges between non-root p and all similar to p
    for (set<int>::iterator nriter = nonroots.begin(); nriter != nonroots.end(); ++nriter) {
        node* p = parts[lyr][*nriter];
        set<node*> simp;
        
        p->get_neighbor_set(edge, simp);
        p->delete_edges(edge);
        for (set<node*>::iterator siter = simp.begin(); siter != simp.end(); ++siter) {
            (*siter)->delete_edges(redge, p);
        }
    }
    keep_parts(lyr, vector<int>(tokeep.begin(), tokeep.end()));
}

void part_lib::get_parent_parts(set<int>& pset, int layer_index, int object_index)
{
    if (object_index < 0 || object_index >= layer_size(layer_index)) return;
    
    node* src = parts[layer_index][object_index];
    int to_part = atom("lyrPrev").get_index();

    foreach_neighbor(src, to_part, iter) {
        lib_data* nd = (lib_data*)neighbor_node_data(iter);
        pset.insert(nd->type);
    }
}

void part_lib::get_part_types(set<int>& pset, int layer_index)
{
    if (layer_index < 0 || layer_index >= layer_count) return;
    
    vector<node*>& lparts = parts[layer_index];

    for (vector<node*>::iterator iter = lparts.begin(); iter != lparts.end(); ++iter) {
        lib_data* nd = (lib_data*)(*iter)->data;
        pset.insert(nd->type);
    }
}

double part_lib::similarity(node* part1, node* part2, double rfactor, int depth, double M(double[], int, int))
{
    typedef pair<part_data_2*, node*> perm_pair_t;
    typedef edge_data_t<double> sim_data_t;

    part_data* pd1 = (part_data*)part1->data;
    part_data* pd2 = (part_data*)part2->data;
    int to_similar = atom("toSimilar").get_index();
    
    if (pd1->layer == 0 && pd2->layer == 0) {
        matrix<double> mask1, mask2;
        ipoint2 c1, c2;

        c1 = pd1->get_mask(mask1, part1, this);
        c2 = pd2->get_mask(mask2, part2, this);
        return mask1.convolve_central(mask2, c2.x, c2.y, c1.x, c1.y);
    }

    int srcname = atom("lyrSrc").get_index();
    int centername = atom("lyrCenter").get_index();
    vector<perm_pair_t> perm1, perm2;

    forall_neighbors(part1, iter1) {
        int ename1 = neighbor_index(iter1);

        if (ename1 == srcname || ename1 == centername) 
            perm1.push_back(perm_pair_t((part_data_2*)neighbor_edge_data(iter1), neighbor_node(iter1)));
    }
    forall_neighbors(part2, iter2) {
        int ename2  = neighbor_index(iter2);

        if (ename2 == srcname || ename2 == centername) 
            perm2.push_back(perm_pair_t((part_data_2*)neighbor_edge_data(iter2), neighbor_node(iter2)));
    }

    ipoint2 c1(pd1->cmx, pd1->cmy);
    ipoint2 c2(pd2->cmx, pd2->cmy);

    // perm1 is fixed and always shorter (or equal) than perm2 !
    if (perm2.size() < perm1.size()) {
        perm2.swap(perm1);
        c2.swap(c1);
    }
    sort(perm2.begin(), perm2.end());

    int size1 = (int)perm1.size();
    int size2 = (int)perm2.size();
    ipoint2 vi1, vi2;
    double bestsim = 0.0;
    double Mparams[20];
    double Mresult;

    do {
        for (int i = 0; i < size1; ++i) {
            if (perm1[i].first == nullptr) vi1.set(0, 0); else vi1.set(perm1[i].first->x, perm1[i].first->y);
            if (perm2[i].first == nullptr) vi2.set(0, 0); else vi2.set(perm2[i].first->x, perm2[i].first->y);
          
            double sim;
            ipoint2 p1 = vi1 - c1, p2 = vi2 - c2;
            double n1 = sqrt((double)p1.norm2());
            double n2 = sqrt((double)p2.norm2());
            double corr;
            sim_data_t* sim_data = (sim_data_t*)get_edge_data(perm1[i].second, perm2[i].second, to_similar);

            if (sim_data != nullptr) sim = sim_data->data;
            else {
                sim = similarity(perm1[i].second, perm2[i].second, rfactor, depth + 1, M);
                add_edge_2(perm1[i].second, perm2[i].second, new sim_data_t(sim), to_similar);
            }
            sim = pow(sim, rfactor);
            
            if (n1 < 1.5 && n2 < 1.5) corr = 1.0;
            else {
                if (n1 < 1.5) corr = 1.0 - n2/5.0;
                else if (n2 < 1.5) corr = 1.0 - n1/5.0;
                else {
                    corr = (double)p1.inner_product(p2)/(n1*n2);
                    corr = max(0.0, corr);
                    corr *= corr;
                    corr *= (n1 > n2) ? sqrt(n2/n1) : sqrt(n1/n2);
                }
            }
            sim *= corr;
            Mparams[i] = sim;
        }
        for (int i = size1; i < size2; ++i) Mparams[i] = 0.0;
        Mresult = M(Mparams, size2, depth);
        if (Mresult > bestsim) bestsim = Mresult;
    } while (next_permutation(perm2));

    return bestsim;
}

// OLD CODE - DOES NOT CONSIDER ALL PERMUTATIONS
/*  forall_neighbors(part1, iter1) {
        int ename1 = neighbor_index(iter1);

        if (ename1 != srcname && ename1 != centername) continue;

        node* pi1 = neighbor_node(iter1);
        part_data_2* edi1 = (part_data_2*)neighbor_edge_data(iter1);
        ipoint2 vi1;

        if (edi1 == nullptr) vi1.set(0, 0); else vi1.set(edi1->x, edi1->y);

        ipoint2 p1 = vi1 - c1;
        double n1 = sqrt((double)p1.norm2());
        double bestsim = 0.0;
        node* best = nullptr;

        forall_neighbors(part2, iter2) {
            int ename2 = iter2->first;

            if (ename2 != srcname && ename2 != centername) continue;

            node* pi2 = neighbor_node(iter2);
            part_data_2* edi2 = (part_data_2*)neighbor_edge_data(iter2);
            ipoint2 vi2;

            if (edi2 == nullptr) vi2.set(0, 0); else vi2.set(edi2->x, edi2->y);

            double sim = similarity(pi1, pi2, rfactor*rfactor, depth + 1, M);
            ipoint2 p2 = vi2 - c2;
            double n2 = sqrt((double)p2.norm2());
            double corr;
            
            if (n1 < 1.5 && n2 < 1.5) corr = 1.0;
            else {
                if (n1 < 1.5) corr = 1.0 - n2/5.0;
                else if (n2 < 1.5) corr = 1.0 - n1/5.0;
                else corr = (double)p1.inner_product(p2)/(n1*n2);
            }
            // cout << corr << " ";
            //double corr = ed1->distr.convolve_central(ed2->distr, ed2->x, ed2->y, ed1->x, ed1->y);
            sim *= max(0.0, corr);
            if (sim > bestsim) bestsim = sim;
        }
        Mparams[Mparams_count++] = bestsim;
    }
    return M(Mparams, Mparams_count, depth);
}
*/

//////////////////////////////////////////////////
/*
        forall_neighbors(part1, iter1) {
            int ename1 = neighbor_index(iter1);
    
            if (ename1 == srcname || ename1 == centername) 
                perm1.push_back(pair<neighbor_edge_data(iter1), neighbor_node(iter1)));
        }
        forall_neighbors(part2, iter2) {
            int ename2  = neighbor_index(iter2);

            if (ename2 == srcname || ename2 == centername) 
                perm2.push_back(neighbor_edge_data(iter2));
        }
        //perm1 is fixed and always shorter (or equal)
        if (perm2.size() < perm1.size())
            perm2.swap(perm1);
        sort(perm2.begin(), perm2.end());

        int size1 = (int)perm1.size();
        ipoint2 c1(pd1->cmx, pd1->cmy);
        ipoint2 c2(pd2->cmx, pd2->cmy);
        ipoint2 vi1, vi2;
        double bestsim = 0.0;
        node* best = nullptr;

        do {
            for (int i = 0; i < size1; ++i) {
                if (perm1[i].first == nullptr) vi1.set(0, 0); else vi1.set(perm1[i]->x, perm1[i]->y);
                if (perm2[i].first == nullptr) vi2.set(0, 0); else vi2.set(perm2[i]->x, perm2[i]->y);
                
                
                double sim = similarity(perm1[i].second, perm2[i].second, rfactor*rfactor, depth + 1, M);
                ipoint2 p1 = vi1 - c1, p2 = vi2 - c2;
                double n1 = sqrt((double)p1.norm2());
                double n2 = sqrt((double)p2.norm2());
                double corr;
                
                if (n1 < 1.5 && n2 < 1.5) corr = 1.0;
                else {
                    if (n1 < 1.5) corr = 1.0 - n2/5.0;
                    else if (n2 < 1.5) corr = 1.0 - n1/5.0;
                    else corr = (double)p1.inner_product(p2)/(n1*n2);
                }
                sim *= max(0.0, corr);
            }
            if (sim > bestsim) { best = p2i; bestsim = sim; }
        } while (next_permutation(perm2));
        Mparams[Mparams_count++] = bestsim;


            






            node* p1i = (node*)iter1->second.first;
            part_data_2* ed1 = (part_data_2*)iter1->second.second;
            ipoint2 pivec;

            if (ed == nullptr) pivec.set(pd1->cmx, pd1->cmy); 
            else pivec.set(

            forall_neighbors(part2, iter2) {
                int ename2 = iter2->first;

                if (ename2 != srcname && ename2 != centername) continue;

                node* p2i = (node*)iter2->second.first;
                part_data_2* ed2 = (part_data_2*)iter2->second.second;
                double sim;

                if (ed1 == nullptr || ed2 == nullptr) {
                    if (ed1 == nullptr && ed2 == nullptr) sim = similarity(p1i, p2i, rfactor*rfactor, depth + 1, M);
                    else sim = 0.0;
                } else {
                    /* sim *=  
                    sim = similarity(p1i, p2i, rfactor*rfactor, depth + 1, M);
                    ipoint2 p1(ed1->x, ed1->y);
                    ipoint2 p2(ed2->x, ed2->y);
                    double n1 = sqrt((double)p1.norm2());
                    double n2 = sqrt((double)p2.norm2());
                    double corr;
                    
                    if (n1 < 1.5 && n2 < 1.5) corr = 1.0;
                    else {
                        if (n1 < 1.5) corr = 1.0 - n2/5.0;
                        else if (n2 < 1.5) corr = 1.0 - n1/5.0;
                        else corr = (double)p1.inner_product(p2)/(n1*n2);
                    }
                    // cout << corr << " ";
                    //double corr = ed1->distr.convolve_central(ed2->distr, ed2->x, ed2->y, ed1->x, ed1->y);
                    sim *= max(0.0, corr);
                }
                
                
            }
            Mparams[Mparams_count++] = bestsim;
        }

        return M(Mparams, Mparams_count, depth);
    }
}
*/

void part_lib::add_similarity_edge(node* p1, node* p2, const vector<int>& perm, double val)
{
    int name = atom("lyrSimRoot");

    if (p1 == nullptr || p2 == nullptr) return;
    if (is_sim_root(p2)) {
        //if (!p1->is_neighbor(p2, name)) add_edge_2(p1, p2, new part_data_sim(perm, val), name);
    } else {
    }
}

void part_lib::add_similarity_edge(int p1, int p2, int layer, const vector<int>& perm, double val)
{
    if (layer < 0 || layer > max_layer_index() || 
        p1 < 0 || p1 >= layer_size(layer) || 
        p2 < 0 || p2 >= layer_size(layer)) return;
    add_similarity_edge(parts[layer][p1], parts[layer][p2], perm, val);
}

double part_lib::defaultMfunction(double* arr, int n, int lyr)
{
    double result = 0.0;
    double thr = ::pow(similarity_threshold, lyr);

    for (int i = 0; i < n; ++i) {
        // cout << "M = " << arr[i] << " thr = " << thr << " lyr = " << lyr << endl;
        if (arr[i] < thr) return 0.0; 
        result += arr[i];
    }
    //cout << endl << endl;
    return result/n;
}

// Returns true if v1 is subset of v2.
bool subset(const vector<part_lib::rs_pair_t>& v1, const vector<part_lib::rs_pair_t>& v2)
{
    for (size_t i = 0; i < v1.size(); ++i) {
        bool ok = false;

        for (size_t j = 0; j < v2.size(); ++j) {
            if (v1[i].first.subset_of(v2[j].first) && intersection_size(v1[i].second, v2[j].second) == (int)v1[i].second.size()) {
                ok = true;
                break;
            }
        }
        if (!ok) return false;
    }
    return true;
}

// Merges coarse (this) library and fine library 'finelib'. Coarse layer is 'clayer' and fine 
// layer is 'flayer'. Fine layer must be a "subset" of coarse layer.
// 
part_lib* part_lib::c2f_merge(part_lib* finelib, int clayer, int flayer)
{
    typedef vector<rs_pair_t> vector_t;

    part_lib* newlib = (part_lib*)get_copy();
    int newlayer = clayer + 1;
    int cbackname = atom("lyrCenterBack");

    newlib->insert_empty_layer(newlayer);

    vector<node*>& oparts = newlib->parts[newlayer + 1];
    vector<node*>& newparts = newlib->parts[newlayer];
    vector<node*>& cparts = newlib->parts[clayer]; // clayer == newlayer - 1
    vector<node*>& fparts = finelib->parts[flayer];
    vector<vector_t> crvs(cparts.size());
    vector<vector_t> orvs(oparts.size());

    newlib->delete_edges(cparts.begin(), cparts.end(), cbackname);
    for (size_t i = 0; i < cparts.size(); ++i) {
        newlib->get_rsrc_neighbors(crvs[i], cparts[i]);
    }
    for (size_t i = 0; i < oparts.size(); ++i) {
        newlib->get_rsrc_neighbors(orvs[i], oparts[i]);
    }

    for (vector<node*>::iterator iter = fparts.begin(); iter != fparts.end(); ++iter) {
        node* fp = *iter;
        rpart_data* pd = new rpart_data(*((rpart_data*)fp->data));
        node* p = newlib->add_node(pd, fp->attr);

        pd->layer = newlayer;
        newparts.push_back(p);

        // Make copies of edges from 'p'. They (should;) go: 
        //   down to the last regular layer, up to the object layer
        forall_neighbors(fp, niter) {
            lib_data* nd = (lib_data*)neighbor_node_data(niter);
            edge_data* ed = neighbor_edge_data(niter);
            node* destn = (nd->layer <= flayer) ? newlib->parts[nd->layer][nd->type] : 
                newlib->parts[newlayer + nd->layer - flayer][nd->type];

            newlib->add_edge_2(p, destn, (edge_data*)ed->get_copy(), neighbor_index(niter));
        }

        vector_t prv;

        newlib->get_rsrc_neighbors(prv, p);

        // Add edges from coarse layer to 'p'; 'cp' (from 'newlayer' - 1) is 
        // connected to 'p' iff RN(p) is subset of RN(pc) with lyrCenterBack.
        // RN is a collection of "rectangles" given by the edge data on 
        // lyrSrc and lyrSrcM edges.
        for (size_t i = 0; i < crvs.size(); ++i) {
            if (subset(prv, crvs[i])) 
                newlib->add_edge_2(cparts[i], p, new part_data_2c(0, 0), cbackname);
        }

        // Add edges from 'p' to object layer 'newlayer' + 1. Connect 'p' to 
        // 'op' iff RN(op) is subset of RN(p) with edge lyrCenterBack
        //for (size_t i = 0; i < orvs.size(); ++i) {
        //    if (subset(orvs[i], prv))
        //        newlib->add_edge_2(p, oparts[i], new part_data_2c(0, 0), cbackname);
        //}

    }
    return newlib;
}

// Inserts new empty layer into the library. All layers >= 'layer' are moved "up".
// All lib_data.layer fields are updated accordingly.
void part_lib::insert_empty_layer(int layer)
{
    if (layer < 0 || layer > max_layer_index()) return;

    for (int l = max_layer_index(); l >= layer; --l) {
        parts[l + 1] = parts[l];

        vector<node*>& p = parts[l + 1];

        for (vector<node*>::iterator iter = p.begin(); iter != p.end(); ++iter) {
            lib_data* pd = (lib_data*)(*iter)->data;
            pd->layer += 1;
        }
    }
    parts[layer].clear();
    layer_count += 1;
}

void part_lib::get_src_neighbors(set<pair<ipoint2, int> >& nbset, node* n)
{
    int srcname = atom("lyrSrc");
    node* cn = n->get_neighbor(atom("lyrCenter"));
    lib_data* cnd = (lib_data*)cn->data;
    lib_data* nd = (lib_data*)n->data;

    foreach_neighbor(n, srcname, iter) {
        node* pn = neighbor_node(iter);
        lib_data* pnd = (lib_data*)neighbor_node_data(iter);
        part_data_2* ed = (part_data_2*)neighbor_edge_data(iter);
        set<int> typeset;

        get_similar_types(pnd->layer, pnd->type, nd->layer, nd->type, ed->x, ed->y, 0.0, typeset);
        typeset.insert(pnd->type);
        for (set<int>::iterator iter = typeset.begin(); iter != typeset.end(); ++iter) {
            nbset.insert(pair<ipoint2, int>(ipoint2(ed->x, ed->y), *iter));
        }
    }
}

void part_lib::get_src_neighbors(set<pair<ipoint2, int> >& nbset, int layer, int part)
{
    if (layer < 0 || layer > max_layer_index() || part < 0 || part > layer_size(layer) - 1)
        return;
    get_src_neighbors(nbset, parts[layer][part]);
}

// Get neighbors from node with attribute R_PART_ATTR.
// The result is a vector of pair<irectangle2, set<int> > where irectangle2 is a "map" range
// and set<int> is a set of types.
void part_lib::get_rsrc_neighbors(vector<rs_pair_t>& result, node* n)
{
    if (!n->is_attr_set(R_PART_ATTR)) return;

    int srcname = atom("lyrSrc");
    lib_data* nd = (lib_data*)n->data;
    set<int> typeset;

    // Note that there is no lyrCenter edge!
    foreach_neighbor(n, srcname, iter) {
        node* pn = neighbor_node(iter);
        lib_data* pnd = (lib_data*)neighbor_node_data(iter);
        part_data_2r* ed = (part_data_2r*)neighbor_edge_data(iter);
        
        typeset.clear();
        get_similar_types(pnd->layer, pnd->type, nd->layer, nd->type, ed->x, ed->y, 0.0, typeset);
        typeset.insert(pnd->type);
        result.push_back(rs_pair_t(ed->rect + ipoint2(ed->x, ed->y), typeset));
    }
}

// Return the type of the center of the part p.
// If the center part does not exist, it returns -1.
int part_lib::get_center_type(node* p)
{
    if (p == nullptr) return -1;
    
    node* pn = p->get_neighbor(atom("lyrCenter"));

    if (pn == nullptr) return -1;
    else return ((lib_data*)pn->data)->type;
}

pair<node*, part_data_2*> part_lib::get_neighbor_pair(node* p, int index)
{
    typedef pair<node*, part_data_2*> result_t;

    if (p == nullptr) return result_t(nullptr, nullptr);
    
    int srcname = atom("lyrSrc");

    foreach_neighbor(p, srcname, niter) {
        part_data_2* ped = (part_data_2*)neighbor_edge_data(niter);
        if (index == ped->index) return result_t(neighbor_node(niter), ped);
    }
    //if (pt.is_zero()) return p->get_neighbor(atom("lyrCenter"));
    return result_t(nullptr, nullptr);
}

bool part_lib::is_vs_layer(int layer)
{
    if (layer < 0 || layer > max_layer_index()) return false;

    vector<node*>& pvec = parts[layer];

    for (int i = 0; i < (int)pvec.size(); ++i) 
        if (pvec[i]->is_attr_set(VS_PART_ATTR) || dynamic_cast<vs_part_data*>(pvec[i]->data) != nullptr)
            return true;
    return false;
}

double part_lib::contraction(int l1, int l2)
{
    double result = 1.0;

    for (int l = l1; l <= l2; ++l) {
        result *= contractions[l];
    }
    return result;
}

void part_lib::add_similarity_edges(int layer, double contraction)
{
    if (layer >= layer_count) return;

    vector<node*>& lparts = parts[layer];
    vector<node*>::iterator iter1, iter2;
    int edgename = atom("toSimilar").get_index();

    for (iter1 = lparts.begin(); iter1 != lparts.end(); ++iter1) {
        for (iter2 = lparts.begin(); iter2 != lparts.end(); ++iter2) {
            double sim;

            if ((*iter1) == (*iter2)) sim = 1.0;
            else {
                sim = similarity(*iter1, *iter2, contraction, 1, defaultMfunction);
                // cout << sim << " ";
            }
            add_edge_2(*iter1, *iter2, new edge_data_t<double>(sim), edgename);
        }
    }
}

double part_lib::part_similarity(int layer, int t1, int t2)
{
    node* p1 = parts[layer][t1];
    node* p2 = parts[layer][t2];
    static int to_similar = atom("toSimilar").get_index();

    foreach_neighbor(p1, to_similar, iter) {
        if (neighbor_node(iter) == p2) 
            return ((edge_data_t<double>*)neighbor_edge_data(iter))->data;
    }
    return 0.0;
}

void part_lib::get_similar_types(int tlayer, int t, int clayer, int ct, int x, int y, double typethresh, map<int, double>& appmap)
{
    static int to_similar = atom("toSimilar").get_index();
    static int to_srcm = atom("lyrSrcM").get_index();
    node* p = parts[tlayer][t];

    foreach_neighbor(p, to_similar, iter) {
        if (((edge_data_t<double>*)neighbor_edge_data(iter))->data >= typethresh) 
            appmap.insert(pair<int, double>(((part_data*)neighbor_node_data(iter))->type, 1.0));
    }
    p = parts[clayer][ct];
    foreach_neighbor(p, to_srcm, iter) {
        part_data_2c* ed = (part_data_2c*)neighbor_edge_data(iter);
        if (ed->x == x && ed->y == y) 
            appmap.insert(pair<int, double>(((part_data*)neighbor_node_data(iter))->type, 1.0));
    }
}

void part_lib::get_similar_types(int tlayer, int t, int clayer, int ct, int x, int y, double typethresh, set<int>& typeset)
{
    static int to_similar = atom("toSimilar").get_index();
    static int to_srcm = atom("lyrSrcM").get_index();
    node* p = parts[tlayer][t];

    foreach_neighbor(p, to_similar, iter) {
        if (((edge_data_t<double>*)neighbor_edge_data(iter))->data >= typethresh) 
            typeset.insert(((part_data*)neighbor_node_data(iter))->type);
    }
    p = parts[clayer][ct];
    foreach_neighbor(p, to_srcm, iter) {
        part_data_2c* ed = (part_data_2c*)neighbor_edge_data(iter);
        if (ed->x == x && ed->y == y) 
            typeset.insert(((part_data*)neighbor_node_data(iter))->type);
    }
}

void part_lib::get_similar_types(int layer, int t, double thresh, vector<int>& result)
{
    static int to_similar = atom("toSimilar").get_index();
    node* p1 = parts[layer][t];

    result.clear();
    foreach_neighbor(p1, to_similar, iter) {
        if (((edge_data_t<double>*)neighbor_edge_data(iter))->data >= thresh) 
            result.push_back(((part_data*)neighbor_node_data(iter))->type);
    }
    if (result.empty()) result.push_back(t);
}

double part_lib::get_thresh(int name, int lyr, int type, double defval)
{
    if (lyr < 1 || lyr > max_layer_index() || type < 0 || type >= (int)parts[lyr].size())
        return defval;
    lib_data* pd = (lib_data*)(parts[lyr][type])->data;
    return pd->get_thresh(name, defval);
}

void part_lib::update_thresh(int name, int lyr, int type, double val)
{
    if (lyr < 1 || lyr > max_layer_index() || type < 0 || type >= (int)parts[lyr].size())
        return;

    node* p = parts[lyr][type];
    lib_data* pd = (lib_data*)p->data;

    pd->td.update_thresh(name, val);
}

void part_lib::set_thresh(int name, int lyr, int type, double val)
{
    if (lyr < 1 || lyr > max_layer_index() || type < 0 || type >= (int)parts[lyr].size())
        return;

    node* p = parts[lyr][type];
    lib_data* pd = (lib_data*)p->data;

    pd->td.set_thresh(name, val);
}

// Get path map from node p (i.e. merges path maps of subparts),
// 'tmap' is a map from index (subpart name) to type. (part_data_2a::app in not necessarily of size one!)
void part_lib::get_path_map(path_map_t& pm, node* p, const map<int, int>& tmap)
{
    int ename = atom("lyrSrc");

    pm.clear();
    foreach_neighbor(p, ename, iter) {
        part_data_2* ed = (part_data_2*)neighbor_edge_data(iter);

        if (ed != nullptr) {
            map<int, int>::const_iterator tmiter = tmap.find(ed->index);

            if (tmiter != tmap.end()) {
                path_map_t& pmi = ed->geo.find(tmiter->second)->second;
            
                for (path_map_t::iterator pmiter = pmi.begin(); pmiter != pmi.end(); ++pmiter) {
                    path_map_t::key_type key = pmiter->first;

                    key.insert(key.begin(), ed->index);
                    pm.insert(path_map_t::value_type(key, pmiter->second));
                }
            }
        }
    }

    part_data_2a* ced = (part_data_2a*)p->get_edge_data(atom("lyrCenter"));

    if (ced != nullptr) {
        map<int, int>::const_iterator tmiter = tmap.find(0);

        if (tmiter != tmap.end()) {
            path_map_t& pmi = ced->geo.find(tmiter->second)->second;

            for (path_map_t::iterator pmiter = pmi.begin(); pmiter != pmi.end(); ++pmiter) {
                path_map_t::key_type key = pmiter->first;

                key.insert(key.begin(), 0);
                pm.insert(path_map_t::value_type(key, pmiter->second));
            }
        }
    }
}

void part_lib::copy_to(streamable* p, cloner& cl)
{
    graph::copy_to(p, cl);

    ((part_lib*)p)->attr = attr;
    ((part_lib*)p)->layer_count = layer_count;
    for (int i = 0; i < MAX_LAYER_NUMBER; ++i) {
        ((part_lib*)p)->contractions[i] = contractions[i];
    }
    for (int i = 0; i < layer_count; ++i) {
        vector<node*>& src_parts = parts[i];
        vector<node*>& dest_parts = ((part_lib*)p)->parts[i];
        for (vector<node*>::iterator iter = src_parts.begin(); iter != src_parts.end(); ++iter) {
            dest_parts.push_back((node*)cl.get_copy(*iter));
        }
    }
}

void part_lib::read_from_stream(istreamer& is)
{
    graph::read_from_stream(is);

    if (is.get_version() < 3.0)
        throw new_libhop_exception("Invalid version of the library.");
    
    unsigned nlyr_parts;

    is.read(attr);
    is.read(layer_count);
    for (int p = 0; p < layer_count; ++p) {
        vector<node*>& lyr_parts = parts[p];
        is.read(nlyr_parts);
        lyr_parts.resize(nlyr_parts);

        for (unsigned lp = 0; lp < nlyr_parts; ++lp) {
            node*& n = lyr_parts[lp];
            is.read((streamable*&)n);
        }
    }
    is.read(contractions);
    is.read(layer_data);
}

void part_lib::write_to_stream(ostreamer& os)
{
    graph::write_to_stream(os);

    vector<node*>::iterator iiter;

    os.write(attr);
    os.write(layer_count);
    for (int i = 0; i < layer_count; ++i) {
        vector<node*>& lyr_parts = parts[i];

        os.write((unsigned)lyr_parts.size());
        for (iiter = lyr_parts.begin(); iiter != lyr_parts.end(); ++iiter) {
            os.write(*iiter);
        }

    }
    os.write(contractions, layer_count);
    os.write(layer_data, layer_count);
}

void part_lib::save(const string& file, int zlib_compression_level /* = Z_NO_COMPRESSION */)
{
    streamable::save(file, zlib_compression_level);
    
    cv::FileStorage fs(file + ".yml", cv::FileStorage::WRITE);

    if (fs.isOpened()) {
        write_cv_storable_data(*fs);
        fs.release();
    }
}

void part_lib::get_regions(int lyr, vector<matrix<bool> >& regions, 
                           vector<ipoint2>& region_centers)
{
    regions.clear();
    region_centers.clear();
    if (layer_count < lyr) return;

    vector<node*>& lyr_parts = parts[lyr - 1];
    vector<node*>::iterator iter;

    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);
    for (iter = lyr_parts.begin(); iter != lyr_parts.end(); ++iter) {
        part_data* d = (part_data*)(*iter)->data;
        regions.push_back(matrix<bool>());
        ipoint2 p = d->get_region_matrix(regions.back(), *iter);
        region_centers.push_back(p);
    }
}

void part_lib::get_regions(int lyr, vector<set<ipoint2> >& regions)
{
    regions.clear();
    if (layer_count < lyr) return;

    vector<node*>& lyr_parts = parts[lyr - 1];
    vector<node*>::iterator iter;

    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);
    for (iter = lyr_parts.begin(); iter != lyr_parts.end(); ++iter) {
        part_data* d = (part_data*)(*iter)->data;
        set<ipoint2> s;

        d->get_region_set(s, *iter);
        regions.push_back(s);
    }
}

void part_lib::drop_images()
{
    for (iter_t iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        
        if (typeid(*(n->data)) == typeid(part_data)) {
            part_data* pd = (part_data*)n->data;

            pd->reset_mask();
            pd->reset_region();
        }
    }
}

img part_lib::lib_image(int lyr, int min, int max, bool show_labels, bool make_border, 
        bool mark_center, bool one_row)
{
    if (layer_count < lyr) return img();
    if (max < 0) max = (int)parts[lyr - 1].size();
    max = std::min<int>(max, (int)parts[lyr - 1].size());
    min = std::max<int>(min, 0);

    vector<img> images;
    double minval = (std::numeric_limits<double>::max)();
    double maxval = (std::numeric_limits<double>::min)();

    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);

    for (int i = min; i < max; ++i) {
        node* n = parts[lyr - 1][i];
        lib_data* nd = (lib_data*)n->data;
        double immin, immax;

		img im = nd->get_image(n, attr, this);

		images.push_back(im);

		if (im.size() <= 0)
			continue;

        images.back().minmax(immin, immax);
        if (immin < minval) minval = immin;
        if (immax > maxval) maxval = immax;
		
    }

    for (int i = min; i < max; ++i) {
        node* n = parts[lyr - 1][i];
        img& im = images[i];

		if (im.size() <= 0)
			continue;

        // Normalize values & convert image to color image
        im -= minval;
        if (minval < maxval) im /= (maxval - minval);
        im = img(im, false, true);

        if (make_border) {
            im.draw_box(irectangle2(0, 0, im.width - 1, im.height - 1), COL_TO_REAL(192, 192, 192));
        }
        if (show_labels) {
            img number = img::number_to_img(i);
            int resizef = -100; //-100*lyr;

            number = number.get_resized(resizef, resizef);
            im.blt(number, 0, 0);
        }
        if (mark_center) {
            part_data* pd = dynamic_cast<part_data*>(n->data);

            if (pd != nullptr) {
                irectangle2 box = pd->get_mask_box(n, this);
                im.draw_big_point(-box.ll.x, -box.ll.y, COL_TO_REAL(255, 0, 0));
            }
        }
    }
    double zero = (minval < maxval) ? -minval/(maxval - minval) : 0.0;
    IMG_COLOR cl;

    cl.blue = cl.green = cl.red = COL_FROM_REAL(zero);
    return one_row ? img::concat_linear_fw(images) : img::concat(images, *((double*)(&cl)));
}

void part_lib::save_all(const char* fname, int lyr, int min /* = 0 */, int max /* = -1 */, 
    bool show_labels /* = true */, bool make_border /* = false */, bool mark_center /* = false */, bool one_row /* false */)
{
    lib_image(lyr, min, max, show_labels, make_border, mark_center, one_row).save_normalized(fname);
}

list<vector<dpoint2> > sample_from_vs_part(vs_part_data* pd, int nsamples)
{
    list<vector<dpoint2> > result;

    if (pd == nullptr || pd->pcad.eigenvalues.rows == 0) 
        return result;

    double sigma1 = sqrt(pd->pcad.eigenvalues.at<double>(0, 0));
    double sigma1d = 3*sigma1/(nsamples + 1);
    double s = -1.5*sigma1;
    cv::Mat mean = pd->pcad.mean;
    
    for (int n = 0; n < nsamples; ++n) {
        cv::Mat delta1 = pd->pcad.eigenvectors.row(0);
        

        if (pd->pcad.eigenvalues.rows == 1) {
            cv::Mat sample = mean + delta1*s;
            vector<dpoint2> pts = partition(sample);

            for (auto piter = pts.begin(); piter != pts.end(); ++piter) 
                (*piter) *= pd->pcad.sizefactor;
            result.push_back(pts);
        } else {
            double sigma2 = sqrt(pd->pcad.eigenvalues.at<double>(1, 0));
            double sigma2d = 3*sigma2/(nsamples + 1);
            double t = -1.5*sigma2;

            for (int m = 0; m < nsamples; ++m) {
                cv::Mat delta2 = pd->pcad.eigenvectors.row(1);
                cv::Mat sample = mean + delta1*s + delta2*t;
                vector<dpoint2> pts = partition(sample);

                for (auto piter = pts.begin(); piter != pts.end(); ++piter) 
                    (*piter) *= pd->pcad.sizefactor;
                result.push_back(pts);
                t += sigma2d;
            }
        }
        s += sigma1d;
    }
    return result; 
}

vector<dpoint2> random_sample_from_vs_part(vs_part_data* pd)
{
    if (pd == nullptr) return vector<dpoint2>();
    else return sample_from_pca(pd->pcad);
}

double part_scale_factor(node* p)
{
    if (p == nullptr) return 0.0;

    vector<dpoint2> pts = cast_vector<dpoint2, ipoint2>(get_library_geo(p));
    pair<double, dpoint2> sdpair = translate_and_scale(pts);

    return sdpair.first;
}

bool is_vs_part(node* p)
{
    int sname = atom("lyrSrc");
    int cname = atom("lyrCenter");

    if (p == nullptr) return false;
    if (dynamic_cast<vs_part_data*>(p->data) != nullptr) 
        return true;

    lib_data* pd = (lib_data*)p->data;

    if (pd->layer == 0) return false;
    forall_neighbors (p, piter) {
        if ((neighbor_index(piter) == sname || neighbor_index(piter) == cname) &&
            is_vs_part(neighbor_node(piter))) return true;
    }
    return false;
}

bool has_vs_children(node* p)
{
    int sname = atom("lyrSrc");
    int cname = atom("lyrCenter");

    if (p == nullptr) return false;

    lib_data* pd = (lib_data*)p->data;

    if (pd->layer == 0) return false;
    forall_neighbors (p, piter) {
        if ((neighbor_index(piter) == sname || neighbor_index(piter) == cname) &&
            is_vs_part(neighbor_node(piter))) return true;
    }
    return false;
}

ipoint2 sample_from_matrix_distribution(matrix<double> m)
{
    int x, y;
    double maxm = m.maximum();
    double val;
    double sum = 0.0;

    for_each_element(m, i) { sum += m[i]; }
    if (sum > 0) m /= sum;
    do {
        x = random_int(0, m.width);
        y = random_int(0, m.height);
        val = random_real()*maxm;
    } while (val > m(x, y));
    return ipoint2(x - m.width/2, y - m.height/2);
}

vector<pair<dpoint2, int> > random_sample_part(part_lib* library, node* n)
{
    typedef pair<dpoint2, int> result_item_t;

    int sname = atom("lyrSrc");
    int cname = atom("lyrCenter");
    int vsmname = atom("lyrVSMember");

    vector<result_item_t> result;
    set<node*> orset;

    if (n == nullptr) return result;

    n->get_neighbor_set(vsmname, orset);
    if (!orset.empty()) {
        vector<node*> v(orset.begin(), orset.end());
        return random_sample_part(library, v[random_int(0, (int)v.size())]);
    }

    part_data* nd = (part_data*)n->data;

    if (nd->layer == 0) {
        result.push_back(result_item_t(dpoint2::zero, nd->type));
        return result;
    }

    double factor = part_data::total_contraction(nd->layer);
    node* c = n->get_neighbor(cname);
    dpoint2 cdelta(-nd->cmx, -nd->cmy);

    if (c != nullptr) {
        part_data* cd = (part_data*)c->data;

        // Add central part pts
        vector<result_item_t> pts = random_sample_part(library, c);

        for (vector<result_item_t>::iterator ptiter = pts.begin(); ptiter != pts.end(); ++ptiter) {
            result_item_t p = *ptiter;

            p.first += (cdelta*factor);
            result.push_back(p);
        }
    }

    // Add other parts
    foreach_neighbor(n, sname, iter) {
        part_data* nd = (part_data*)neighbor_node_data(iter);
        part_data_2* ed = (part_data_2*)neighbor_edge_data(iter);
        vector<int> avec;
        //dpoint2 rnddelta = dpoint2::zero; 
        dpoint2 rnddelta = (dpoint2)sample_from_matrix_distribution(ed->distr);

        for (part_data_2::app_map_t::iterator aiter = ed->app.begin(); aiter != ed->app.end(); ++aiter) {
            avec.push_back(aiter->first);
        }

        node* varn = library->parts[nd->layer][avec[random_int(0, (int)avec.size())]];
        //part_data* varnd = (part_data*)varn->data;

        vector<result_item_t> pts = random_sample_part(library, varn);

        for (vector<result_item_t>::iterator ptiter = pts.begin(); ptiter != pts.end(); ++ptiter) {
            result_item_t p = *ptiter;
            dpoint2 delta = (cdelta + rnddelta + dpoint2(ed->x, ed->y))*factor;

            //merge_mask_maps(mask, varnd->mask, (cdelta + ipoint2(ed->x, ed->y)) * factor, 1.0, libtype);            
            p.first += delta;
            result.push_back(p);
        }
    }
    
    return result;

}

vector<dpoint2> random_sample_from_part(node* p)
{
    int sname = atom("lyrSrc");
    int cname = atom("lyrCenter");
    int vsmname = atom("lyrVSMember");

    vector<dpoint2> result;

    if (p == nullptr) return result;
    if (!is_vs_part(p)) return cast_vector<dpoint2, ipoint2>(get_library_geo(p));

    vector<pair<int, ipoint2> > gpieces = get_library_geo_pieces(p, 3);
    part_data* pd = (part_data*)p->data;

    if (!has_vs_children(p)) {
        vector<dpoint2> result = random_sample_from_vs_part(dynamic_cast<vs_part_data*>(p->data));
        
        //// normalize
        //if (gpieces.empty()) {
        /*double factor = 0.0;
        int count = 0;

        foreach_neighbor(p, vsmname, niter) {
            factor += part_scale_factor(neighbor_node(niter));
            ++count;
        }
        if (count > 0) factor /= count;
        if (factor == 0.0) factor = 4.5*part_data::total_contraction(pd->layer);
            
        mult_vector(result, factor);
        */
        //}
        //    vector<ipoint2> pts = extract_second<ipoint2>(gpieces.begin(), gpieces.end());
        //    irectangle2 ptsbox = irectangle2::bounding_rectangle(pts.begin(), pts.end());
        //    drectangle2 resultbox = drectangle2::bounding_rectangle(result.begin(), result.end());

        //    double f = sqrt((double)ptsbox.size2())/sqrt(resultbox.size2());
        //    mult_vector(result, f);
        //}
        return result;
    }

    forall_neighbors (p, piter) {
        node* pn = neighbor_node(piter);
        part_data_2* ed = (part_data_2*)neighbor_edge_data(piter);

        if (neighbor_index(piter) != sname && neighbor_index(piter) != cname) 
            continue;

        vector<dpoint2> partsample = random_sample_from_part(pn);

        // translate
        vector<ipoint2> pts;

        for (int i = 0; i < (int)gpieces.size(); ++i) 
            if (gpieces[i].first == ed->index) pts.push_back(gpieces[i].second);

        irectangle2 ptsbox = irectangle2::bounding_rectangle(pts.begin(), pts.end());
        drectangle2 samplebox = drectangle2::bounding_rectangle(partsample.begin(), partsample.end());
        dpoint2 cdelta = (dpoint2)ptsbox.center() - samplebox.center();
        double factor = sqrt((double)ptsbox.size2())/sqrt(samplebox.size2());

        for (int i = 0; i < (int)partsample.size(); ++i) {
            //partsample[i] = partsample[i]*factor;
            partsample[i] += cdelta;
        }

        vector<dpoint2> dpts = cast_vector<dpoint2, ipoint2>(pts);

        //result.insert(result.end(), dpts.begin(), dpts.end());
        result.insert(result.end(), partsample.begin(), partsample.end());
    }
    return result;
}

    //void save_all_sc(const char* fname, int lyr, int min = 0, int max = -1, bool show_labels = true, bool one_row = false);
    //void save_all_sc(const string& outdir, bool show_labels = true);

void part_lib::save_sc(const string& outdir, bool show_labels, bool sc, bool save_variations /* = true */)
{
    for (int l = 0; l < layer_count; ++l) {
        if (layer_size(l) > 0 && !parts[l][0]->is_attr_set(CATEGORY_NODE_ATTR)) {
            string fname = outdir + "L_" + (l + 1) + "_image.png";

            if (sc) save_all_sc(fname.c_str(), l + 1, 0, -1, show_labels, true, save_variations);
            else save_all(fname.c_str(), l + 1, 0, -1, false, false, false, true);
        }
    }
}

void part_lib::save_all_sc(const char* fname, int lyr, int min /* = 0 */, int max /* = -1 */, bool show_labels /* = true */,
    bool one_row /* = false */, bool save_variations /* = true */)
{
    const color colors[] = {color(255, 0, 0), color(0, 255, 0), color(0, 0, 255), 
        color(0, 255, 255), color(255, 255, 255)};

    int sname = atom("lyrSrc");
    int cname = atom("lyrCenter");

    if (layer_count < lyr) return;
    if (max < 0) max = (int)parts[lyr - 1].size();
    max = std::min<int>(max, (int)parts[lyr - 1].size());
    min = std::max<int>(min, 0);

    vector<img> images;
    double minval = (std::numeric_limits<double>::max)();
    double maxval = (std::numeric_limits<double>::min)();

    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);

    for (int i = min; i < max; ++i) {
        node* n = parts[lyr - 1][i];
        vs_part_data* vsnd = dynamic_cast<vs_part_data*>(n->data);
        double immin, immax;
        img im(0, 0, false);
        vector<pair<ipoint2, int> > pts;
        set<ipoint2> cloud;
        int maxx = INT_MIN;
        int minx = INT_MAX;
        int maxy = INT_MIN;
        int miny = INT_MAX;

        //if (!is_sim_root(n) && !save_variations) 
        //    continue;

        if (vsnd != nullptr && save_variations) {
            vector<dpoint2> dpts = partition(vsnd->pcad.mean);

            pts.resize(dpts.size());
            for (int i = 0; i < (int)dpts.size(); ++i) {
                pts[i].first = (ipoint2)(dpts[i]*vsnd->pcad.sizefactor);
                pts[i].second = 0;
                int a = pts[i].first.x, b = pts[i].first.y;

                if (a > maxx) maxx = a;
                if (a < minx) minx = a;
                if (b > maxy) maxy = b;
                if (b < miny) miny = b;
            }

            list<vector<dpoint2> > lpts = sample_from_vs_part(vsnd, 5);

            for (list<vector<dpoint2> >::iterator iter = lpts.begin(); iter != lpts.end(); ++iter) {
                for (vector<dpoint2>::iterator viter = iter->begin(); viter != iter->end(); ++viter) {
                    ipoint2 p = (ipoint2)(*viter);

                    cloud.insert(p);
                }
            }

        } else {
            forall_neighbors (n, niter) {
		        node* nn = neighbor_node(niter);
                part_data_2a* nned = dynamic_cast<part_data_2a*>(neighbor_edge_data(niter));

                if (nned != nullptr && !nned->geo.empty()) {
                    for (part_data_2a::geo_map_t::iterator giter = nned->geo.begin(); giter != nned->geo.end(); ++giter) {
                        for (path_map_t::iterator miter = giter->second.begin(); miter != giter->second.end(); ++miter) {
                            int a = miter->second.p.x, b = miter->second.p.y;
                            part_data_2* nned2 = dynamic_cast<part_data_2*>(neighbor_edge_data(niter));

                            if (a > maxx) maxx = a;
                            if (a < minx) minx = a;
                            if (b > maxy) maxy = b;
                            if (b < miny) miny = b;
                            pts.push_back(pair<ipoint2, int>(ipoint2(a, b), nned2 ? nned2->index : 0));
                        }

                    }
                }
            }
        }

        if (!pts.empty()) {
            vector<ipoint2> debugpts = extract_first<ipoint2>(pts.begin(), pts.end());
            double debugfac = translate_and_scale(debugpts).first;


            im.resize(maxx - minx + 20, maxy - miny + 20, 0.0);

            for (set<ipoint2>::iterator citer = cloud.begin(); citer != cloud.end(); ++citer) {
                im.set_at_safe(citer->x - minx + 10, citer->y - miny + 10, COL_TO_REAL(128, 128, 128));
            }
            for (vector<pair<ipoint2, int> >::iterator piter = pts.begin(); piter != pts.end(); ++piter) {
                color c = colors[piter->second % 5]; //color::from_hsv((piter->second % 6) * 60, 1.0, 1.0);

                im.draw_big_point(piter->first.x - minx + 10, piter->first.y - miny + 10, 
                    COL_TO_REAL(c.red(), c.green(), c.blue()));
            }

		    images.push_back(im);
        }

		if (im.size() <= 0)
			continue;

        images.back().minmax(immin, immax);
        if (immin < minval) minval = immin;
        if (immax > maxval) maxval = immax;
		
    }

    for (int i = min; i < max && i < images.size(); ++i) {
        node* n = parts[lyr - 1][i];
        img& im = images[i];

		if (im.size() <= 0)
			continue;

        // Normalize values & convert image to color image
        //im -= minval;
        //if (minval < maxval) im /= (maxval - minval);
        //im = img(im, false, true);

        if (show_labels) {
            img number = img::number_to_img(i);
            int resizef = -100; //-100*lyr;

            number = number.get_resized(resizef, resizef);

            if (is_sim_root(n)) {
                im.draw_rectangle(irectangle2(0, 0, number.width + 1, number.height + 1), COL_TO_REAL(196, 196, 196));
                number = number.to_color_img(COL_TO_REAL(0, 196, 196));
            } else {
                number = number.to_color_img(COL_TO_REAL(196, 196, 196));
            }
            im.blt(number, 1, 1);
            
        }
    }
    //double zero = (minval < maxval) ? -minval/(maxval - minval) : 0.0;
    //IMG_COLOR cl;

    //cl.blue = cl.green = cl.red = COL_FROM_REAL(zero);

    //img result = img::concat(images, *((double*)(&cl)));
    img result = one_row ? img::concat_linear_fw(images) : img::concat(images);
    result.save_normalized(fname);
}

void part_lib::save_all_sc_mma(const char* fname, int lyr, int min /* = 0 */, int max /* = -1 */)
{
    int sname = atom("lyrSrc");
    int cname = atom("lyrCenter");

    if (layer_count < lyr) return;
    if (max < 0) max = (int)parts[lyr - 1].size();
    max = std::min<int>(max, (int)parts[lyr - 1].size());
    min = std::max<int>(min, 0);

    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);

    ofstream os;

    os.open(fname);
    if (os.fail()) {
        cout << "Can not open file '" << fname << "." << endl;
        return;
    }

    os << '{';

    for (int i = min; i < max; ++i) {
        node* n = parts[lyr - 1][i];
        lib_data* nd = (lib_data*)n->data;
        bool first = true;

        if (i != min) os << ',';
        os << '{' << i << ',' << '{';

        forall_neighbors (n, niter) {
		    node* nn = neighbor_node(niter);
            part_data_2a* nned = dynamic_cast<part_data_2a*>(neighbor_edge_data(niter));


            if (nned != nullptr && !nned->geo.empty()) {
                for (part_data_2a::geo_map_t::iterator giter = nned->geo.begin(); giter != nned->geo.end(); ++giter) {
                    for (path_map_t::iterator miter = giter->second.begin(); miter != giter->second.end(); ++miter) {
                        vector<int> path = miter->first;
                        int a = miter->second.p.x, b = miter->second.p.y;
                        part_data_2* nned2 = dynamic_cast<part_data_2*>(neighbor_edge_data(niter));

                        path.insert(path.begin(), nned2 ? nned2->index : 0);
                        if (!first) os << ','; else first = false;
                        os << '{';
                        os << '{';
                        for (int i = 0; i < (int)path.size(); ++i) {
                            if (i != 0) os << ',';
                            os << path[i];
                        }
                        os << '}';
                        os << ',' << '{' << a << ',' << b << '}';
                        os << '}';
                    }
                }
            }
        }
        os << '}' << '}' << '\n';
    }
    os << '}' << endl;
    os.close();
}

void part_lib::save_all(const char* fname, int lyr, const vector<int>& l, bool show_labels, bool mark_center)
{
    if (layer_count < lyr) return;

    int max = (int)parts[lyr - 1].size();
    vector<img*> images;

    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);

    for (int i = 0; i < (int)l.size(); ++i) {
        int index = l[i];
        
        if (index >= 0 && index < max) {
            node* n = parts[lyr - 1][index];
            part_data* nd = (part_data*)n->data;

            images.push_back(new img(nd->get_image(n, attr, this)));
            if (show_labels) {
                img number = img::number_to_img(i);
                int resizef = -100; //-100*lyr;

                number = number.get_resized(resizef, resizef);
                images.back()->blt(number, 0, 0);
            }
            if (mark_center) {
                irectangle2 box = nd->get_mask_box(n, this);
                images.back()->draw_big_point(-box.ll.x, -box.ll.y, COL_TO_REAL(255, 0, 0));
            }

        }
    }
    img result = img::concat(images);
    result.save_normalized(fname);
    for (vector<img*>::iterator iter = images.begin(); iter != images.end(); ++iter)
        delete *iter;
}

void part_lib::save_all_regions(const char* fname, int lyr, int min, int max)
{
    if (layer_count < lyr) return;
    if (max < 0) max = (int)parts[lyr - 1].size();
    max = std::min<int>(max, (int)parts[lyr - 1].size());
    min = std::max<int>(min, 0);

    vector<img> images;
    for (int i = min; i < max; ++i) {
        node* p = parts[lyr - 1][i];
        images.push_back(region_to_img(((part_data*)p->data)->get_region(p)));
    }
    img result = img::concat(images);
    result.save_normalized(fname);
}

void part_lib::info(ostream& os)
{
    os << "Number of layers: " << layer_count << endl;
    for (int i = 0; i < layer_count; ++i) {
        vector<int> v = get_root_parts_map(i);
        int m1count = 0;
        bool thresholds = false;
        bool svmthresholds = false;

        for_each(v.begin(), v.end(), [&m1count](int i) { if (i == -1) ++m1count; });
        for (auto p = parts[i].begin(); p != parts[i].end(); ++p) {
            lib_data* pd = (lib_data*)(*p)->data;
            if (!thresholds) thresholds = !pd->td.empty();
            if (!svmthresholds) {
                vs_part_data* vspd = dynamic_cast<vs_part_data*>((*p)->data);
                if (vspd != nullptr) svmthresholds = vspd->svmt.svm != nullptr;
            }
        }

        os << "Layer: " << i << endl;
        os << "   - Number of parts: " << parts[i].size() << endl;
        os << "   - number of parts indexing to next layer: " << ((int)v.size() - m1count) << endl;
        os << "   - Thresholds are " << (thresholds ? "" : "not ") << "set." << endl;
        os << "   - SVM thresholds are " << (svmthresholds ? "" : "not ") << "set." << endl;
    }
}

void part_lib::save_part(const char* fname, int lyr, int part, double threshold /* = 0.0 */)
{
    if (lyr < 1 || lyr > layer_count || part < 0 || part >= layer_size(lyr - 1))
        return;

    vector<img> images;
    node* p = parts[lyr - 1][part];
    int edgename = atom("toSimilar").get_index();

    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);
    images.push_back(((part_data*)(p->data))->get_image(p, attr, this));
    foreach_neighbor(p, edgename, iter) {
        node* q = (node*)iter->second.first;
        edge_data_t<double>* qd = (edge_data_t<double>*)iter->second.second;

        //cout << " " << qd->data;
        if (qd->data > threshold) images.push_back(((part_data*)(q->data))->get_image(q, attr, this));
    }
    img result = img::concat(images);
    result.save_normalized(fname);
}

// save layer 1 masks
void part_lib::save_visview_filters(const string& dir)
{
    int count = (int)(parts[0].size());

	for (int i = 0; i < count; i++) {
        string name = dir + "filter_" + (i + 1) + ".txt";
        node* p = parts[0][i];
        ((part_data*)p->data)->get_image(p, attr, this).save_visview(name.c_str());
	}
}

// save centers
void part_lib::save_visview_centers(const string& dir, int lyr)
{
    if (lyr < 1 || lyr > layer_count) return;

    int lyrindex = lyr - 1;
    int count = (int)(parts[lyrindex].size());
    int edgename = atom("lyrSrc").get_index();
    int edgenameM = atom("lyrSrcM").get_index();

    string name = dir + "L_" + lyr + "_centers.txt";
    ofstream os;

    os.open(name.c_str(), ios::trunc);
    if (os.fail()) return;

    if (lyrindex == 0) {
        for (int i = 0; i < count; ++i)
            os << "0,0" << '\n';
    } else {
        for (int i = 0; i < count; i++) {
            node* n = parts[lyrindex][i];
            part_data* nd = dynamic_cast<part_data*>(n->data);
            int countM = 0;

            os << "0,0";  // center of central part = (0, 0)
            foreach_neighbor(n, edgename, iter) {
                part_data_2* ed = (part_data_2*)neighbor_edge_data(iter);
                os << ',' << ed->x << ',' << ed->y;
            }
            foreach_neighbor(n, edgenameM, iter) {
                part_data_2* ed = (part_data_2*)neighbor_edge_data(iter);
                os << ',' << ed->x << ',' << ed->y;
                if (++countM >= 1) break;
            }
            if (nd == nullptr) {
                os << ",0,0";
            } else {
                os << ',' << nd->cmx << ',' << nd->cmy;
            }
            os << '\n';
        }
	}
    os.close();
}

void part_lib::save_visview_contraction(const string& dir, int lyr)
{
    if (lyr < 1 || lyr > layer_count) return;

    string name = dir + "L_" + lyr + "_factor.txt";
    ofstream os(name.c_str(), ios::trunc);

    os << contractions[lyr - 1] << endl;
    os.close();
}

// save parts
void part_lib::save_visview_parts(const string& dir, int lyr)
{
    if (lyr < 1 || lyr > layer_count) return;

    int lyrindex = lyr - 1;
    int count = (int)(parts[lyrindex].size());
    int edgename = atom("lyrSrc");
    int cedgename = atom("lyrCenter");
    int c1edgename = atom("lyrCenter1");
    int vsmname = atom("lyrVSMember");

    string name;
    ofstream os, osm;

    name = dir + "L_" + lyr + "_parts.txt";
    os.open(name.c_str(), ios::trunc);
    name = dir + "L_" + lyr + "_maps.txt";
    osm.open(name.c_str(), ios::trunc);
    if (os.fail() || osm.fail()) return;

    if (lyrindex == 0) {
        for (int i = 1; i <= count; ++i)
            os << i << ',' << i << ",1," << i << ",1" << endl;
    } else {
        for (int i = 0; i < count; i++) {
            node* n = parts[lyrindex][i];
            vector<iipair> types;   // .first = type, .second = group
            vector<rmatrix*> maps;

            os << i + 1 << ','; // name
            if (n->is_attr_set(CATEGORY_NODE_ATTR)) {
                foreach_neighbor(n, cedgename, iter) {
                    part_data* pd = (part_data*)neighbor_node_data(iter);
                    types.push_back(iipair(pd->type + 1, 1));
                }
                os << (types.empty() ? 0 : types[0].first) << ',' << types.size();
            } else {
                if (n->is_attr_set(VS_PART_ATTR)) 
                    n = n->get_neighbor(vsmname);

                node* cn = n->get_neighbor(c1edgename);
                lib_data* pcd = nullptr;
                int centertype;
                int subpartcount = 0;

                if (cn == nullptr) cn = n->get_neighbor(cedgename);
                if (cn == nullptr) cn = n->get_neighbor(edgename);
                pcd = (part_data*)cn->data;
                centertype = pcd->type + 1;
                foreach_neighbor(n, edgename, iter) {
                    part_data* pd = (part_data*)neighbor_node_data(iter);
                    part_data_2* ed = (part_data_2*)neighbor_edge_data(iter);
                    set<int> typeset;

                    get_similar_types(pd->layer, pd->type, pcd->layer, pcd->type, ed->x, ed->y, 0.0, typeset);
                    subpartcount++;
                    types.push_back(iipair(pd->type + 1, subpartcount));
                    maps.push_back(&ed->distr);
                    for (set<int>::iterator siter = typeset.begin(); siter != typeset.end(); ++siter) {
                        if (*siter != pd->type) { 
                            types.push_back(iipair(*siter + 1, subpartcount));
                            maps.push_back(&ed->distr);
                        }
                    }
                    
                }
                os << centertype << ',' << types.size() + 1 << ',' << centertype << ",1";
            }
            for (size_t p = 0; p < types.size(); ++p) 
                os << ',' << types[p].first << ',' << types[p].second + 1;
            os << endl;
            for (size_t p = 0; p < maps.size(); ++p) {
                rmatrix* pmap = maps[p];

                if (p > 0) osm << ',';
                osm << pmap->width << ',' << pmap->height << ',';
                pmap->print_as_vector(osm);
            }
            osm << endl;
        }
    }
    os.close();
    osm.close();
}

void part_lib::save_visview_parts(const string& dir)
{
    for (int i = 0; i < layer_count; ++i)
        save_visview_parts(dir, i + 1);
}

void part_lib::save_visview_contractions(const string& dir)
{
    for (int i = 0; i < layer_count; ++i)
        save_visview_contraction(dir, i + 1);
}

void part_lib::save_visview_centers(const string& dir)
{
    for (int i = 0; i < layer_count; ++i)
        save_visview_centers(dir, i + 1);
}

void part_lib::write_vgr_label(ostream& os, node* n, int count)
{
    part_data* d = (part_data*)n->data;
    os << '\"' << d->type << '\"';
}

void part_lib::save_mma(ostream& os, int lyr)
{
    int z = lyr - 1;

    if (z < 0 || z >= max_layer()) {
        os << '{' << '}';
        return;
    }
    
    vector<node*>& zparts = parts[z];
    int to_prev = atom("lyrSrc");
    int to_center = atom("lyrCenter");
    bool first = true;
    int offsetz1 = 1, offsetz = 1;

    for (int i = 0; i < z - 1; ++i) offsetz1 += (int)parts[i].size();
    if (z > 0) offsetz = offsetz1 + (int)parts[z - 1].size();

    os << '{';
    for (vector<node*>::iterator iter = zparts.begin(); iter != zparts.end(); ++iter) {
        node* n = *iter;
        part_data* pd = (part_data*)n->data;
        node* c = n->get_neighbor(to_center);
        lib_data* cd = (lib_data*)n->data;

        if (first) first = false; else os << ',';
        os << '{' << offsetz + pd->type << ',' << offsetz1 + cd->type << '}';
        foreach_neighbor(n, to_prev, niter) {
            lib_data* npd = (lib_data*)neighbor_node_data(niter);
            
            if (first) first = false; else os << ',';
            os << '{' << offsetz + pd->type << ',' << offsetz1 + npd->type << '}';
        }
    }
    os << '}';
}

void part_lib::save_mma_back(ostream& os, int lyr)
{
    int z = lyr - 1;

    if (z < 0 || z >= max_layer()) {
        os << '{' << '}';
        return;
    }

    vector<node*>& zparts = parts[z];
    int backname = atom("lyrCenterBack");
    bool first = true;

    os << '{';
    for (vector<node*>::iterator iter = zparts.begin(); iter != zparts.end(); ++iter) {
        node* n = *iter;
        part_data* pd = (part_data*)n->data;

        foreach_neighbor(n, backname, niter) {
            lib_data* npd = (lib_data*)neighbor_node_data(niter);
            
            if (first) first = false; else os << ',';
            os << '{' << pd->type << ',' << npd->type << '}';
        }
    }
    os << '}';
}

void part_lib::save_mma2(ostream& os, int z)
{
    if (z < 0 || z >= max_layer()) {
        os << '{' << '}';
        return;
    }
    
    vector<node*>& zparts = parts[z];
    int tocenter = atom("lyrCenter");
    int tosrc = atom("lyrSrc");
    int tosrcm = atom("lyrSrcM");
    int toprev = atom("lyrPrev");
    os << '{';
    for (vector<node*>::iterator iter = zparts.begin(); iter != zparts.end(); ++iter) {
        node* n = *iter;
        part_data* pd = (part_data*)n->data;
        bool first = true;

        if (iter != zparts.begin()) os << ',';
        os << '{' << pd->type << ',' << '{' << pd->cmx << ',' << pd->cmy << '}' << ',';
        pd->td.save_mma(os);
        os << ',';
        os << '{';
        foreach_neighbor(n, tocenter, niter) {
            part_data* npd = (part_data*)neighbor_node_data(niter);
            part_data_2c* epd = (part_data_2c*)neighbor_edge_data(niter);
            
            if (first) first = false; else os << ',';
            os << '{' << npd->type << ",{0,0}}";
        }
        foreach_neighbor(n, tosrc, niter) {
            part_data* npd = (part_data*)neighbor_node_data(niter);
            part_data_2c* epd = (part_data_2c*)neighbor_edge_data(niter);
            
            if (first) first = false; else os << ',';
            os << '{' << npd->type << ',' << '{' << epd->x << ',' << epd->y << '}' << '}';
        }
        foreach_neighbor(n, tosrcm, niter) {
            part_data* npd = (part_data*)neighbor_node_data(niter);
            part_data_2c* epd = (part_data_2c*)neighbor_edge_data(niter);
            
            if (first) first = false; else os << ',';
            os << '{' << npd->type << ',' << '{' << epd->x << ',' << epd->y << '}' << '}';
        }
        os << '}' << '}';
    }
    os << '}';
}

void part_lib::save_mma(ostream& os)
{
    os << '{';
    for (int lyr = 1; lyr <= max_layer(); ++lyr) {
        if (lyr > 1) os << ',' << endl;
        save_mma(os, lyr);
    }
    os << '}' << endl;
}

void part_lib::save_mma_back(ostream& os)
{
    os << '{';
    for (int lyr = 1; lyr <= max_layer(); ++lyr) {
        if (lyr > 1) os << ',' << endl;
        save_mma_back(os, lyr);
    }
    os << '}' << endl;
}

// Save library as a graph in Mathematica format. It exports a list of two lists:
// - a list of edges: {Subscript[type,lyr]->Subscript[type,lyr],...}; 
//    if edgel == true: {{Subscript[type,lyr]->Subscript[type,lyr],"edge name"},...}
//    Only edges from the set 'edges' are exported. If 'edges' is empty all
//    edges are exported!
// - a list of vertex coordinates: {Subscript[type,lyr]->{x,y},...}
//   Where y == lyr and x is chosen s.t. that the vertices in the same layer are 
//   "centered" in 0.
// Only layers [0, maxl) are exported; in each layer types [0, maxt) are exported.
// Hint for the use: 
//   {edges, coordinates}=Get[file]; 
//   GraphPlot[edges,VertexCoordinateRules->coordinates];
//   
void part_lib::save_as_graph_mma(ostream& os, int maxl, int maxt, const set<int>& edges, bool edgel /* = true */)
{
    bool first = true;
    
    if (maxl < 0 || maxl > layer_count) 
        maxl = layer_count;

    os << '{';

    // Export the list of edges.
    os << '{';
    for (int l = 0; l < maxl; ++l) {
        vector<node*>& pv = parts[l];
        int maxi = std::min<int>((int)pv.size(), maxt);

        for (int i = 0; i < maxi; ++i) {
            node* p = pv[i];

            forall_neighbors(p, iter) {
                if (edges.empty() || edges.find(neighbor_index(iter)) != edges.end()) {
                    lib_data* pnd = (lib_data*)neighbor_node_data(iter);

                    if (pnd->layer < maxl && pnd->type < maxt) {
                        if (first) first = false; else os << ',';
                        if (edgel) os << '{';
                        os << "Subscript[" << i << ',' << l << "]->Subscript["
                            << pnd->type << ',' << pnd->layer << ']';
                        if (edgel) os << ',' << '\"' << atom::get_name(neighbor_index(iter)) << '\"' << '}';
                    }
                }
            }
        }
        os << endl;
    }

    // Export the coordinates.
    first = true;
    os << "},{";
    for (int l = 0; l < maxl; ++l) {
        vector<node*>& pv = parts[l];
        int maxi = std::min<int>((int)pv.size(), maxt);

        for (int i = 0; i < maxi; ++i) {
            if (first) first = false; else os << ',';
            os << "Subscript[" << i << ',' << l << "]->{" << i - (maxi - 1)/2.0 << ',' << 10*l << '}';
        }
        os << endl;
    }
    os << '}';
    os << '}' << endl;

}

void part_lib::display_thresholds(const vector<node*>& parts)
{
    cout << "Complete data of parts (for thresholds, -1 indicates that no threshold is set)" << endl;

    for (vector<node*>::const_iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        lib_data* pd = (lib_data*)(*iter)->data;

        cout << "  #" << pd->type;
        cout << " layer: " << pd->layer;
        cout << " R thresh: " << pd->get_thresh(R_RESPONSE, -1);
        cout << " G thresh: " << pd->get_thresh(G_RESPONSE, -1);
        cout << " RR thresh: " << pd->get_thresh(RR_RESPONSE, -1);
        cout << " S thresh: " << pd->get_thresh(S_RESPONSE, -1);
        cout << endl;
    }
}

void part_lib::display_layer_info(int layer) {

	if (layer >= (int)this->parts.size()) { 
		cout  << "Library does not have layer '" << layer << "'" << endl;
		return;
	}

	vector<node*>& parts = this->parts[layer];
    int to_part = atom("lyrSrc").get_index();
    int to_center = atom("lyrCenter").get_index();
    int to_forbidden = atom("lyrForbidden").get_index();
    int to_prev = atom("lyrPrev").get_index();

    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* n = *iter;
        node* m;
        part_data* pd = (part_data*)n->data;
        
        cout << pd->type << " [" << pd->cmx << ',' << pd->cmy << "] == ";
        if ((m = n->get_neighbor(to_center)) != nullptr) cout << ((part_data*)m->data)->type << "(C)";
        foreach_neighbor(n, to_part, niter) {
            part_data* pd = (part_data*)neighbor_node_data(niter);
            part_data_2* ed = (part_data_2*)neighbor_edge_data(niter);
            cout << ", " << pd->type << " [" << ed->x << ',' << ed->y << ']';
        }
        foreach_neighbor(n, to_forbidden, niter) {
            part_data* pd = (part_data*)neighbor_node_data(niter);
            part_data_2* ed = (part_data_2*)neighbor_edge_data(niter);
            cout << ", " << pd->type << "(F) [" << ed->x << ',' << ed->y << ']';
        }
        cout << " | ";
        foreach_neighbor(n, to_prev, niter) {
            part_data* pd = (part_data*)neighbor_node_data(niter);
            part_data_2c* ed = (part_data_2c*)neighbor_edge_data(niter);
            cout << ", " << pd->type << " [" << ed->x << ',' << ed->y << ']';
        }
        cout << endl;
    }
    cout << endl;
}

void part_lib::update_part_data_contractions()
{
    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);
}

vector<irectangle2> part_lib::get_predicted_boxes(int l)
{
	vector<irectangle2> result;

	if (l < 0 || l > max_layer_index()) return result;

	update_part_data_contractions();
	result.resize(layer_size(l), irectangle2());
	for (int i = 0; i < layer_size(l); ++i) {
        node* p = parts[l][i];
        vector<ipoint2> pts = get_library_geo(p);
        part_data* pd = (part_data*)p->data;
            
		//irectangle2 recbox = pd->get_mask_box(p, library);  // better calculation of center...
        irectangle2 geobox = irectangle2::bounding_rectangle(pts.begin(), pts.end());
            
        result[i] = geobox - geobox.center();// + recbox.center();
	}
	return result;
}


void part_lib::update_var_fields() {

    int srcname = atom("lyrSrc");
    int centername = atom("lyrCenter");

    for (list<node*>::iterator niter = this->nodes.begin(); niter != this->nodes.end(); ++niter) {
        node* n = *niter;

        foreach_neighbor(n, srcname, iter) {
            part_data* pd = (part_data*)neighbor_node_data(iter);
            part_data_2* ed = dynamic_cast<part_data_2*>(neighbor_edge_data(iter));

			if (ed != nullptr && ed->app.find(pd->type) == ed->app.end()) 
				ed->app.insert(pair<int, pair<int, double> >(pd->type, pair<int, double>(ed->index, 1.0)));
        }

        // there is only one center, though, but it is convenient to use defines
        foreach_neighbor(n, centername, iter) {
            part_data* pd = (part_data*)neighbor_node_data(iter);
            part_data_2a* ed = dynamic_cast<part_data_2a*>(neighbor_edge_data(iter));

            if (ed != nullptr) {
				if (ed->app.find(pd->type) == ed->app.end()) 
					ed->app.insert(pair<int, pair<int, double> >(pd->type, pair<int, double>(0, 1.0)));
			} else neighbor_edge_data(iter) = new part_data_2a(pd->type);
        }
    }
}

img part_lib::get_image(int layer, int part) 
{
	if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);

	img im;

	if (layer < 0 || layer >= (int)parts.size()) return im;
	if (part < 0 || part >= (int)parts[layer].size()) return im;

    node* n = parts[layer][part];
    im = ((part_data*)(n->data))->get_image(n, 0, this);

	return im;
}

int part_lib::get_basic_part_number(int layer, int type)
{
    return get_basic_part_number(parts[layer][type]);
}


int part_lib::get_basic_part_number(node* p)
{
    int result = ((part_data*)p->data)->bpcount;

    if (result > 0) return result;
    else {
        scope_lock slock(lock);

        return ((part_data*)p->data)->get_basic_part_number(p);
    }
}

void part_lib::init_basic_part_numbers(int layer)
{
    if (layer < 0 || layer > max_layer_index()) 
        return;
    for (auto iter = parts[layer].begin(); iter != parts[layer].end(); ++iter) {
        node* p = *iter;

        ((part_data*)p->data)->get_basic_part_number(p);
    }
}


void part_lib::get_object_parts(set<int>& result, int category_layer, const string& catname)
{
    int nbname = atom("lyrPrev");

    if (category_layer < 1 || category_layer > max_layer_index()) 
        return;

	// find category node with name sm
	node* n = nullptr;

	for (vector<node*>::iterator iter = parts[category_layer].begin(); 
		    iter != parts[category_layer].end(); ++iter) {
		n = *iter;
		cpart_data* nd = dynamic_cast<cpart_data*>(n->data);
		if (nd != nullptr && nd->name == catname) break;
	}
	if (n != nullptr) {
		foreach_neighbor(n, nbname, iter) {
			lib_data* pd = (lib_data*)neighbor_node_data(iter);
			result.insert(pd->type);
		}
	}
}

void part_lib::edge_paths(vector<edge_path>& paths, node* p)
{
    int tosrc = atom("lyrSrc");
    int tocenter = atom("lyrCenter");

    vector<edge_path> qpaths;
    node* q = p->get_neighbor(tocenter);

    paths.clear();
    if (q != nullptr) {
        lib_data* qd = (lib_data*)q->data;

        edge_paths(qpaths, q);
        for (vector<edge_path>::iterator iter = qpaths.begin(); iter != qpaths.end(); ++iter) {
            paths.push_back(*iter);
            paths.back().push_back(0, qd->type);
        }
    }
    foreach_neighbor (p, tosrc, iter) {
        node* q = neighbor_node(iter);
        lib_data* qd = (lib_data*)q->data;
        part_data_2* qed = (part_data_2*)neighbor_edge_data(iter);

        edge_paths(qpaths, q);
        for (vector<edge_path>::iterator iter = qpaths.begin(); iter != qpaths.end(); ++iter) {
            paths.push_back(*iter);
            paths.back().push_back(qed->index, qd->type);
        }
    }
}

node* get_neighbor_by_index(node* p, int index) 
{
    if (index == 0) return p->get_neighbor(atom("lyrCenter"));
    else {
        int name = atom("lyrSrc");

        foreach_neighbor (p, name, iter) {
            part_data_2* ed = (part_data_2*)neighbor_edge_data(iter);

            if (ed->index == index) return neighbor_node(iter);
        }
    }
    return nullptr;
}

void part_lib::read_cv_storable_data(CvFileStorage* fs)
{
    for (iter_t niter = nodes.begin(); niter != nodes.end(); ++niter) {
        node* n = *niter;
        cv_storable_data* cvsd = dynamic_cast<cv_storable_data*>(n->data);

        if (cvsd != nullptr) cvsd->read_from_storage(fs);
    }
}

void part_lib::write_cv_storable_data(CvFileStorage* fs)
{
    for (iter_t niter = nodes.begin(); niter != nodes.end(); ++niter) {
        node* n = *niter;
        cv_storable_data* cvsd = dynamic_cast<cv_storable_data*>(n->data);

        if (cvsd != nullptr) cvsd->write_to_storage(fs);
    }
}

template<class K, class V> void write_map_yml(cv::FileStorage& fs, const map<K, V>& m)
{
    fs << "[";
    for (auto iter = m.begin(); iter != m.end(); ++iter) {
        fs << "{";
        fs << "key" << iter->first;
        fs << "value" << iter->second;
        fs << "}";
    }
    fs << "]";
}

void write_app_map_yml(cv::FileStorage& fs, const map<int, pair<vector<int>, double> >& m)
{
    fs << "[";
    for (auto iter = m.begin(); iter != m.end(); ++iter) {
        fs << "{";
        fs << "key" << iter->first;
		fs << "value" << iter->second.second;
        fs << "}";
    }
    fs << "]";
}

void write_points_yml(cv::FileStorage& fs, const vector<ipoint2>& pts)
{
    fs << "[";
    for (auto iter = pts.begin(); iter != pts.end(); ++iter) {
        fs << "{";
        fs << "x" << iter->x;
        fs << "y" << iter->y;
        fs << "}";
    }
    fs << "]";
}

void part_lib::write_yml(cv::FileStorage& fs)
{
    // Contractions
    fs << "layers" << max_layer_index() + 1;

    fs << "contractions";
    fs << "[";
    for (int i = 1; i < MAX_LAYER_NUMBER; ++i) 
        fs << contractions[i];
    fs << "]";

    // Part nodes
    fs << "parts";
    fs << "[";
    for (int l = 0; l <= max_layer_index(); ++l) {
        for (int t = 0; t < layer_size(l); ++t) {
            node* p = parts[l][t];
            part_data* pd = (part_data*)p->data;

            // Write part data
            fs << "{";
            fs << "type" << pd->type;
            fs << "layer" << pd->layer;
            fs << "cmx" << pd->cmx; 
            fs << "cmy" << pd->cmy;
            if (p->is_attr_set(CATEGORY_NODE_ATTR)) 
                fs << "category" << ((cpart_data*)p->data)->name;
            fs << "geometry";
            write_points_yml(fs, get_library_geo(p));
            fs << "neighbors";
            fs << "[";

            // Write all edges
            forall_neighbors(p, niter) {
                int name = neighbor_index(niter);
                lib_data* pnd = (lib_data*)neighbor_node_data(niter);

                if (name == atom("lyrCenter")) {
                    fs << "{";
                    part_data_2a* ped = (part_data_2a*)neighbor_edge_data(niter);

                    fs << "name" << "lyrCenter";
                    fs << "target_type" << pnd->type;
                    fs << "target_layer" << pnd->layer;
                    fs << "app";
                    write_app_map_yml(fs, ped->app);
                    fs << "index" << 0;
                    fs << "}";
                } else if (name == atom("lyrCenterBack")) {
                    part_data_2c* pced = dynamic_cast<part_data_2c*>(neighbor_edge_data(niter));

                    fs << "{";
                    fs << "name" << "lyrCenterBack";
                    fs << "target_type" << pnd->type;
                    fs << "target_layer" << pnd->layer;
                    fs << "x" << ((pced == nullptr) ? 0 : pced->x);
                    fs << "y" << ((pced == nullptr) ? 0 : pced->y);
                    fs << "}";
                } else if (name == atom("lyrSrc")) {
                    part_data_2* ped = (part_data_2*)neighbor_edge_data(niter);

                    fs << "{";
                    fs << "name" << "lyrSrc";
                    fs << "target_type" << pnd->type;
                    fs << "target_layer" << pnd->layer;
                    fs << "app";
                    write_app_map_yml(fs, ped->app);
                    fs << "index" << ped->index;
                    fs << "x" << ped->x;
                    fs << "y" << ped->y;
                    fs << "distr" << cv::Mat1d(ped->distr.height, ped->distr.width, ped->distr.ptr(0, 0));
                    fs << "}";
                } else if (name == atom("lyrSimRoot")) {
                    fs << "{";
                    fs << "name" << "lyrSimRoot";
                    fs << "target_type" << pnd->type;
                    fs << "target_layer" << pnd->layer;
                    fs << "}";
                }
            }
            fs << "]";
            fs << "}";
        }
    }
    fs << "]";
}

vector<int> part_lib::get_root_parts_map(int layer)
{
    int rootedge = atom("lyrSimRoot");
    int simedge = atom("lyrSimilar");

    if (layer < 0 || layer > max_layer_index()) 
        return vector<int>();

    vector<node*>& parts = this->parts[layer];
    vector<int> result(parts.size(), -1);
    set<int> roots;

    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* p = *iter;
            
        foreach_neighbor(p, rootedge, niter) {
            lib_data* pnd = (lib_data*)neighbor_node_data(niter);
            roots.insert(pnd->type);
        }
    }

    int count = 0;

    for (auto riter = roots.begin(); riter != roots.end(); ++riter) {
        result[*riter] = count++;
    }
    return result;
}


// Maps path 'ppath' starting from 'p' to equivalent path starting from 'q'.
// Returns true if qpath is returned and false if:
// - p or q is null, - p !~ q, - ppath is too long
/*bool part_lib::path_matching(vector<int>& qpath, node* p, const vector<int>& ppath, node* q)
{
    int name = atom("lyrSimilar");

    qpath.clear();

    if (p == nullptr || q == nullptr) return false;

    lib_data* pd = (lib_data*)p->data;
    lib_data* qd = (lib_data*)q->data;
    part_data_sim* ed = (part_data_sim*)p->get_edge_data(name);

    if (ed == nullptr) return false;
    if (pd->layer == 0 && !ppath.empty()) return false;

    if (path_matching(qpath, get_neighbor_by_index(p)
    


}

// Maps paths 'ppaths' starting from 'p' to equivalent paths starting from 'q'.
bool part_lib::path_matching(vector<vector<int> >& qpaths, node* p, const vector<vector<int> >& ppaths, node* q)
{
    int name = atom("lyrSimilar");


    return true;
}
*/

#ifdef OPENCL
// function required by opencl implementation
////////////////////
void part_lib::ocl_make_data() {

	// we need to copy all those connections too
	int edge_type_mapping[][2] = { 
		{ atom("lyrCenterBack").get_index(), OCL_LYR_CENTER_BACK_EDGE },
		{ atom("lyrSrc").get_index(), OCL_LYR_SRC_EDGE },
		{ atom("lyrForbidden").get_index(), OCL_LYR_FORBIDDEN_EDGE },
		{ atom("toPrevLayer").get_index(), OCL_TO_PREV_LAYER_EDGE },
		{ atom("toSimilar").get_index(), OCL_TO_SIMILAR_EDGE },
		{ atom("toLayer0").get_index(), OCL_TO_LAYER_0 }
		
	};
	
	int num_layers = this->layer_count;

	// make room for all layers
	ocl_parts.resize(num_layers);
	ocl_edges.resize(num_layers);
	ocl_apps.resize(num_layers);

	/**
	 * Initializing parts (part_data).
	 */ 

	// now go for each layer
	for (int l = 0; l < num_layers; l++) {
		// reserve enough space for all parts in this layer
		vector<node*> layer_parts = parts[l];
		int parts_count = layer_parts.size();

		ocl_part_data* parts = new ocl_part_data[parts_count];
		memset(parts, 0, parts_count * sizeof(ocl_part_data));

		int connections_count_all = 0;
		int app_count_all = 0;

		// should copy all parts into this new memory - but make sure to copy them in correct order
		for (int i = 0; i < parts_count; i++) {
			// get node 
			node* n = layer_parts[i];
			// and data in node
			part_data* n_data = dynamic_cast<part_data*>(n->data);

			// set pointer from n_data to this opencl structure
			n_data->ocl_struct = &parts[i];

			// now copy values to new memory each at a time
			parts[i].cmx = n_data->cmx;
			parts[i].cmy = n_data->cmy;
			parts[i].layer = n_data->layer;
			parts[i].type = n_data->type;
			parts[i].attr = n->attr;
			parts[i].bpcount = n_data->get_basic_part_number(n);
			parts[i].thresh.R = n_data->td.get_thresh(R_RESPONSE, -1);
			parts[i].thresh.G = n_data->td.get_thresh(G_RESPONSE, -1);
			parts[i].thresh.RR = n_data->td.get_thresh(RR_RESPONSE, -1);			

			// and count how many types of edges (connections we have)
			for (int j = 0; j < OCL_EDGE_COUNT; j++) {
				// get index from atom
				int atom_index = edge_type_mapping[j][0];
				// count how many elements we have
				int count = n->count_neighbors(atom_index);
				// and save to variable
				connections_count_all += count;

				// count how many app values have each connection
				foreach_neighbor(n, atom_index, iter) {
					node* pc = neighbor_node(iter); 
                    part_data* pc_data = dynamic_cast<part_data*>(pc->data);
					edge_data* pc_edge_d = (edge_data*)neighbor_edge_data(iter); 

					if (pc_edge_d != nullptr && typeid(*pc_edge_d) == typeid(part_data_2)) {
						part_data_2* pc_edge = dynamic_cast<part_data_2*>(pc_edge_d);

						app_count_all += pc_edge->app.size();
					}					
				}
				
			}			
		}


		// we need to create memory for all edges since we only now know correct size
		ocl_part_data_2* edges = new ocl_part_data_2[connections_count_all];
		memset(edges, 0, connections_count_all * sizeof(ocl_part_data_2));


		ocl_app_data* apps = new ocl_app_data[app_count_all];
		memset(apps, 0, app_count_all * sizeof(ocl_app_data));

		// save opencl parts and edges pointer and continue with next layer (edges MUST be initialized in new loop)
		ocl_parts[l] = pair<ocl_part_data*,int>(parts, parts_count);
		ocl_edges[l] = pair<ocl_part_data_2*,int>(edges, connections_count_all);
		ocl_apps[l] = pair<ocl_app_data*,int>(apps, app_count_all);
	}

	/**
	 * Initializing edges / connections (part_data_2).
	 * We must first initialize parts for all layer, since edges may
	 * point to any layer and we will need to make opencl pointer (=offset)
	 * to data in global memory.
	 */ 
	for (int l = 0; l < num_layers; l++) {
		// reserve enough space for all parts in this layer
		vector<node*> layer_parts = parts[l];
		int parts_count = layer_parts.size();

		// get previusly created parts and edges for this layer
		ocl_part_data* parts = ocl_parts[l].first;
		ocl_part_data_2* edges = ocl_edges[l].first;
		ocl_app_data* apps = ocl_apps[l].first;

		int edge_counter = 0;
		int app_counter = 0;

		// and go again through all parts to initialize their edges
		for (int i = 0; i < parts_count; i++) {
			// get node
			node* n = layer_parts[i];
			// and data in node
			part_data* n_data = (part_data*)n->data;

			// for each edge type go through all edges
			for (int j = 0; j < OCL_EDGE_COUNT; j++) {
				// get index from atom
				int atom_index = edge_type_mapping[j][0];
				int edge_index = edge_type_mapping[j][1];

				// then copy all neighbors to global memory (also count how many we will copy)
				int copied = 0;
				foreach_neighbor(n, atom_index, iter) {
					node* pc = neighbor_node(iter); 
                    part_data* pc_data = dynamic_cast<part_data*>(pc->data);
					edge_data* pc_edge_d = (edge_data*)neighbor_edge_data(iter); 

					// pc_edge might be null for some atom_index 
					if (pc_edge_d != nullptr) {
						
						// set default app value (that is 0)
						edges[edge_counter + copied].app.offset = 0;
						edges[edge_counter + copied].app.size = 0;

						// encode edge_data based on type
						if (typeid(*pc_edge_d) == typeid(part_data_2)) {
							part_data_2* pc_edge = dynamic_cast<part_data_2*>(pc_edge_d);

							pc_edge->ocl_struct = &edges[edge_counter + copied];

							edges[edge_counter + copied].cast_type = OCL_EDGE_DATA_PART_DATA_2_TYPE;

							edges[edge_counter + copied].x = pc_edge->x;
							edges[edge_counter + copied].y = pc_edge->y;
							edges[edge_counter + copied].type = edge_index;
							edges[edge_counter + copied].gdistr.mean = pc_edge->gdistr.first;
							edges[edge_counter + copied].gdistr.variance = pc_edge->gdistr.second;
							
							// copy each individual values of distribution matrix since we need
							// also cast it from double to float

							if (pc_edge->distr.size() > OCL_EDGE_DATA_2_DISTRIBUTION_MAX_SIZE) {
								cout << "ERROR: Unable to convert library for OpenCL support. Using parts with edge distribution size '" << pc_edge->distr.size()  << "' while opencl supports only max distribution size '" << OCL_EDGE_DATA_2_DISTRIBUTION_MAX_SIZE << "' " << endl;
								cout << "\tAdjust OCL_EDGE_DATA_2_DISTRIBUTION_MAX_SIZE in utils.cl ADN cl_structures.h" << endl;
								throw new std::exception();
							}


							edges[edge_counter + copied].distr_size = pc_edge->distr.size();
							for (int k = 0; k < edges[edge_counter + copied].distr_size; k++) {
								edges[edge_counter + copied].distr[k] = (float)pc_edge->distr[k];
							}

							// copy each app value 
							int app_values_added = 0;
							for(part_data_2::app_map_t::iterator app_iter = pc_edge->app.begin(); 
									app_iter != pc_edge->app.end(); app_iter++) {
								apps[app_counter + app_values_added].type = app_iter->first;
                                apps[app_counter + app_values_added].value = app_iter->second.second;

								app_values_added++;
							}
							// set offset and size for app mapping
							edges[edge_counter + copied].app.offset = app_counter;
							edges[edge_counter + copied].app.size = app_values_added;

							// increment app_counter for next connection
							app_counter += app_values_added;

						} else if (typeid(*pc_edge_d) == typeid(part_data_2c)) {
							part_data_2c* pc_edge = dynamic_cast<part_data_2c*>(pc_edge_d);

							pc_edge->ocl_struct = &edges[edge_counter + copied];

							edges[edge_counter + copied].cast_type = OCL_EDGE_DATA_PART_DATA_2C_TYPE;

							edges[edge_counter + copied].x = pc_edge->x;
							edges[edge_counter + copied].y = pc_edge->y;						
						
						} else if (typeid(*pc_edge_d) == typeid(part_data_2r)) {
							part_data_2r* pc_edge = dynamic_cast<part_data_2r*>(pc_edge_d);

							pc_edge->ocl_struct = &edges[edge_counter + copied];

							edges[edge_counter + copied].cast_type = OCL_EDGE_DATA_PART_DATA_2R_TYPE;

							edges[edge_counter + copied].x = pc_edge->x;
							edges[edge_counter + copied].y = pc_edge->y;

							// rectange will be encoded into distr table
							edges[edge_counter + copied].distr[0] = pc_edge->rect.ll.x;
							edges[edge_counter + copied].distr[1] = pc_edge->rect.ll.y;
							edges[edge_counter + copied].distr[2] = pc_edge->rect.ur.x;
							edges[edge_counter + copied].distr[3] = pc_edge->rect.ur.y;

						} else {
							// do nothing
						}			
					}

					ocl_part_data* layer_pointed_parts = this->ocl_parts[pc_data->layer].first;

					// calculate offset using pointer arithmetic
					edges[edge_counter + copied].node.offset = (pc_data->ocl_struct - layer_pointed_parts);
					edges[edge_counter + copied].node.layer = pc_data->layer;
					
					copied++;
				}
				
				// set offset i.e. position of first element in global memory (end of previous edge_counter)
				// and its size				
				parts[i].edges_loc[edge_index].offset = copied > 0 ? edge_counter : 0;
				parts[i].edges_loc[edge_index].size = copied;

				// increment edge_counter for number of copied edges
				edge_counter += copied;
			}
		}

		// re-set sizes for ocl_edges and ocl_apps in case some memory slot were not used
		ocl_edges[l].second = edge_counter;
		ocl_apps[l].second = app_counter;
	}

	ocl_made = true;
}

#endif

// path_map_stat
///////////////////////////////////////////////////////////////////////////////

void path_map_stat::add(const path_map_t& m)
{
    for (path_map_t::const_iterator iter2 = m.begin(); iter2 != m.end(); ++iter2) {
        pair<path_map_t::iterator, bool> ibpair = pm.insert(path_map_t::value_type(iter2->first, iter2->second));

        if (!ibpair.second) {
            ibpair.first->second.h += iter2->second.h;
            ibpair.first->second.p += iter2->second.p;
        }
        ++cm[iter2->first];
    }
}

path_map_t path_map_stat::get() const
{
    path_map_t result;

    get(result);
    return result;
}

void path_map_stat::get(path_map_t& m) const
{
    m.clear();

    for (path_map_t::const_iterator iter = pm.begin(); iter != pm.end(); ++iter) {
        count_map_t::const_iterator cmiter = cm.find(iter->first);

        if (cmiter == cm.end()) {
            cout << "Error in path_map_stat::get(m)" << endl;
            throw;
        }

        path_map_t::mapped_type& val = m[iter->first];

        val = iter->second;
        divide_vector(val.h, cmiter->second);
        val.p /= cmiter->second;
    }
}

// pca_data
///////////////////////////////////////////////////////////////////////////////


void pca_data::write_to_stream(ostreamer& os) const
{
    cv_write_to_stream(os, mean);
    cv_write_to_stream(os, eigenvectors);
    cv_write_to_stream(os, eigenvalues);
    os.write(sizefactor);
}

void pca_data::read_from_stream(istreamer& is)
{
    cv_read_from_stream(is, mean);
    cv_read_from_stream(is, eigenvectors);
    cv_read_from_stream(is, eigenvalues);
    if (is.get_version() < 3.2) sizefactor = 1.0;
    else is.read(sizefactor);
}

// pca_learning
///////////////////////////////////////////////////////////////////////////////

int debugint = 0;

pca_data pca_learning::get_pca_data(int maxcomp /* = 5 */) const
{
    typedef list<vector<dpoint2> > dsamples_t;

    int minsample = INT_MAX;
    dsamples_t newsamples;

    for (samples_t::const_iterator siter = samples.begin(); siter != samples.end(); ++siter) {
        newsamples.push_back(cast_vector<dpoint2, ipoint2>(inhibit_point_set(*siter, 2)));
        
        if ((int)newsamples.back().size() < minsample) 
            minsample = (int)newsamples.back().size();
    }

    if (newsamples.size() == 1) {
        vector<dpoint2>& pts = newsamples.back();
        double factor = translate_and_scale(pts).first;
        cv::Mat mean = flatten(pts);
        cv::Mat eigenvectors(1, mean.cols, CV_64F, cv::Scalar(0.0));
        cv::Mat eigenvalues(1, 1, CV_64F, cv::Scalar(0.0));

        return pca_data(mean, eigenvectors, eigenvalues, factor);
    } else {
        cv::Mat pcadata((int)newsamples.size(), 2*minsample, CV_64F);
        int n = 0;
        double factorsum = 0.0;

        for (dsamples_t::iterator siter = newsamples.begin(); siter != newsamples.end(); ++siter) {
            vector<dpoint2>& pts = *siter;
            vector<int> perm;

            resize_vector(pts, minsample);
            factorsum += translate_and_scale(pts).first;
            if (siter == newsamples.begin()) perm = vector_range(0, (int)pts.size() - 1);
            else perm = point_matching(pts, newsamples.front());
            
            for (int i = 0; i < (int)perm.size(); ++i) {
                dpoint2 p = pts[perm[i]];

                pcadata.at<double>(n, 2*i) = p.x;
                pcadata.at<double>(n, 2*i + 1) = p.y;
            }
            ++n;
        }

        //string fname = string("c:\\temp\\file") + (debugint++) + ".csv";
        //ofstream os(fname.c_str());
        //cv_write_csv<double>(os, pcadata);
        //os.close();

        cv::PCA pca(pcadata, cv::Mat(), CV_PCA_DATA_AS_ROW, maxcomp);

        return pca_data(pca.mean, pca.eigenvectors, pca.eigenvalues, factorsum/n);
    }
}

// shape_learning_p
///////////////////////////////////////////////////////////////////////////////

/*void shape_learning_p::update(const vector<ipoint2>& pts)
{
    point_stat_t lpts;

    lpts.reserve(pts.size());
    for (vector<ipoint2>::const_iterator iter = pts.begin(); iter != pts.end(); ++iter) {
        lpts.push_back(labeled_point_t(0, *iter));
    }
    update(lpts);
}*/

// Note: 'pts' must (should) be inhibited and scaled!
void shape_learning_p::update(const point_stat_t& pts)
{
    samples.push_back(pts);
}

// Aligns samples to library model; returns premuted samples.
void shape_learning_p::align_pieces(list<vector<ipoint2> >& dsamples, const vector<pair<int, ipoint2> >& libmodel) const
{
    typedef map<int, vector<ipoint2> > map_t;

    map_t libmap;

    for (vector<pair<int, ipoint2> >::const_iterator piter = libmodel.begin(); piter != libmodel.end(); ++piter)
        libmap[piter->first].push_back(piter->second);

    for (samples_t::const_iterator siter = samples.begin(); siter != samples.end(); ++siter) {
        map_t smap;

        for (vector<pair<int, ipoint2> >::const_iterator piter = siter->begin(); piter != siter->end(); ++piter)
            smap[piter->first].push_back(piter->second);
        if (smap.size() != libmap.size())
            continue;
        
        vector<int> match = piecewise_point_matching_p(*siter, libmodel);
        vector<ipoint2> sample = extract_second<ipoint2>(siter->begin(), siter->end());

        permute_and_resize(sample, match);
        dsamples.push_back(sample);
    }
}

// Note: 'libmodel' must (should) be inhibited and scaled!
// 'maxcomp': max number of pca vectors (components)
// 'scalefactor' is factor of scaling...
pca_data shape_learning_p::get_pca_data(const vector<pair<int, ipoint2> >& libmodel, double scalefactor, int maxcomp /* = 5 */) const
{
    //vector<pair<int, ipoint2> > ipts = get_library_geo_pieces(p, 5);
    //translate_and_scale(ipts);
    typedef list<vector<dpoint2> > dsamples_t;

    dsamples_t dsamples;
    list<vector<ipoint2> > isamples;

    align_pieces(isamples, libmodel);
    for (auto isiter = isamples.begin(); isiter != isamples.end(); ++isiter) {
        dsamples.push_back(cast_vector<dpoint2, ipoint2>(*isiter));
        translate_and_scale(dsamples.back());
    }

    if (dsamples.empty()) {
        return pca_data();
    } else if (dsamples.size() == 1) {
        vector<dpoint2>& pts = dsamples.back();

        //translate_and_scale(pts);
        cv::Mat mean = flatten(pts);
        cv::Mat eigenvectors(1, mean.cols, CV_64F, cv::Scalar(0.0));
        cv::Mat eigenvalues(1, 1, CV_64F, cv::Scalar(0.0));

        return pca_data(mean, eigenvectors, eigenvalues, scalefactor);
    } else {
        cv::Mat pcadata((int)dsamples.size(), 2*(int)libmodel.size(), CV_64F);
        int n = 0;

        for (dsamples_t::iterator siter = dsamples.begin(); siter != dsamples.end(); ++siter) {
            vector<dpoint2>& pts = *siter;

            for (int i = 0; i < (int)pts.size(); ++i) {   // pts.size() == libmodel.size() !!!
                pcadata.at<double>(n, 2*i) = pts[i].x;
                pcadata.at<double>(n, 2*i + 1) = pts[i].y;
            }
            ++n;
        }

        //string fname = string("c:\\temp\\file") + (debugint++) + ".csv";
        //ofstream os(fname.c_str());
        //cv_write_csv<double>(os, pcadata);
        //os.close();

        cv::PCA pca(pcadata, cv::Mat(), CV_PCA_DATA_AS_ROW, maxcomp);

        return pca_data(pca.mean, pca.eigenvectors, pca.eigenvalues, scalefactor);
    }
}

svm_data shape_learning_p::get_svm_data(const vector<pair<int, ipoint2> >& libmodel) const
{
    typedef list<vector<ipoint2> > isamples_t;

    isamples_t isamples;
    vector<ipoint2> libpts = extract_second<ipoint2>(libmodel.begin(), libmodel.end());

    translate_and_scale(libpts);
    align_pieces(isamples, libmodel);
    for (auto isiter = isamples.begin(); isiter != isamples.end(); ++isiter) 
        translate_and_scale(*isiter);

    if (isamples.size() < 5) {  // Add artificial samples (resized libpts)
        isamples.clear();
        for (double factor = 0.8; factor < 1.21; factor += 0.05) {
            isamples.push_back(libpts);
            
            vector<ipoint2>& v = isamples.back();
            
            for_each(v.begin(), v.end(), [factor](ipoint2& p) { p.x = (int)(p.x*factor); p.y = (int)(p.y*factor); });
        }
    }

    // Train data vector: (x1, y1, squared size1, x2, y2, squared size2, ...)
    int nsamples = (int)libpts.size();
    int nfeatures = 3*nsamples;
    cv::Mat trainData(isamples.size(), nfeatures, CV_32FC1);
    int scount = 0;

    for (auto siter = isamples.begin(); siter != isamples.end(); ++siter) {
        vector<ipoint2>& samplev = *siter;
        float* ptr = trainData.ptr<float>(scount);

        for (int i = 0; i < (int)samplev.size(); ++i) {
            int x = samplev[i].x - libpts[i].x;
            int y = samplev[i].y - libpts[i].y;

            *(ptr++) = (float)(x);
            *(ptr++) = (float)(y);
            *(ptr++) = (float)(x*x + y*y);
        }
        ++scount;
    }

    cv::Mat responses(trainData.rows, 1, CV_32SC1, cv::Scalar_<int>(1));
    //cv::Mat samples(trainData.rows, 1, CV_8UC1, cv::Scalar_<uchar>(1));
    //cv::Mat selectedFeatures(nfeatures, 1, CV_8UC1, cv::Scalar_<uchar>(1));

    cv::SVMParams params;
    cv::SVM* svm = new cv::SVM();

    params.svm_type = cv::SVM::ONE_CLASS;
    params.kernel_type = cv::SVM::LINEAR;
    params.nu = 0.9;

    svm->train(trainData, responses, cv::Mat(), cv::Mat(), params);
   
    /* 
    DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG
    cout << "nu: " << params.nu << endl;
    for (int r = 0; r < trainData.rows; ++r) 
        cout << "r: " << r << ", " << svm->predict(trainData.row(r), true) << endl;
    cout << "z: " << svm->predict(cv::Mat(1, nfeatures, CV_32FC1, cv::Scalar_<float>(0.0)), true) << endl;

    int dummy;
    cin >> dummy;
    DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG DEBUG 
    */

    return svm_data(svm);
  
}

// Global functions
///////////////////////////////////////////////////////////////////////////////

path_map_t::const_iterator zero(const path_map_t& pm)
{
    if (pm.empty()) return pm.end();

    path_map_t::const_iterator iter = pm.begin();

    if (iter->first.empty()) return iter;
    else if (*max_element(iter->first.begin(), iter->first.end()) == 0) return iter;
    else return pm.end();
}

// Augments 'geo' maps with missing items in 'full'; returns the number of added paths.
// If 'normalize' is true, then added paths are normalized to path 0
int augment_path_map(path_map_t& geo, const path_map_t& full, bool normalize /* = false */)
{
    ipoint2 norm(0, 0);
    int result = 0;

    if (normalize) {
        path_map_t::const_iterator geozero = zero(geo);
        path_map_t::const_iterator fullzero = zero(full);

        if (geozero != geo.end() && fullzero != full.end())
            norm = geozero->second.p - fullzero->second.p;
    }
    for (path_map_t::const_iterator fiter = full.begin(); fiter != full.end(); ++fiter) {
        if (geo.find(fiter->first) == geo.end()) {
            path_map_t::iterator giter = geo.insert(path_map_t::value_type(fiter->first, fiter->second)).first;

            giter->second.p += norm;
            ++result;
        }
    }
    return result;
}

void add_to_path_map(path_map_t& pm, const path_map_t& pm2)
{
    for (path_map_t::const_iterator iter2 = pm2.begin(); iter2 != pm2.end(); ++iter2) {
        pair<path_map_t::iterator, bool> ibpair = pm.insert(path_map_t::value_type(iter2->first, iter2->second));

        if (!ibpair.second) {
            ibpair.first->second.h += iter2->second.h;
            ibpair.first->second.p += iter2->second.p;
        }
    }
}

void divide_path_map(path_map_t& pm, int n)
{
    for (path_map_t::iterator iter = pm.begin(); iter != pm.end(); ++iter) {
        divide_vector(iter->second.h, n);
        iter->second.p /= n;
    }
}

path_map_t reduce_path_map(const path_map_t& pm, int radius)
{
    path_map_t result;
    irectangle2 box;

    for (path_map_t::const_iterator iter = pm.begin(); iter != pm.end(); ++iter) {
        box.eat(iter->second.p);
    }

    matrix<int> m(box.x_dim() + 2*radius + 1, box.y_dim() + 2*radius + 1, 0);
    
    for (path_map_t::const_iterator iter = pm.begin(); iter != pm.end(); ++iter) {
        const ipoint2& p = iter->second.p;
        ipoint2 q(p.x - box.ll.x + radius, p.y - box.ll.y + radius);

        if (m(q.x, q.y) == 0) {
            m.set_region_c(q.x, q.y, radius, radius, 1);
            result.insert(path_map_t::value_type(iter->first, iter->second));
        } 
    }
    return result;
}

// Keeps only those keys in 'pm' which are also in 'srcpm'
path_map_t synchronize_path_map(const path_map_t& pm, const path_map_t& srcpm)
{
    path_map_t result;

    for (path_map_t::const_iterator iter = pm.begin(); iter != pm.end(); ++iter) {
        if (srcpm.find(iter->first) != srcpm.end()) {
            result.insert(path_map_t::value_type(iter->first, iter->second));
        }
    }
    return result;
}

// Extends path map 'pm' with 'pm1', with each path prepended 'ext'.
void path_map_union(path_map_t& pm, const path_map_t& pm1, int ext)
{
    for (path_map_t::const_iterator pm1iter = pm1.begin(); pm1iter != pm1.end(); ++pm1iter) {
        path_map_t::key_type key = pm1iter->first;

        key.insert(key.begin(), ext);
        pm.insert(path_map_t::value_type(key, pm1iter->second));
    }
}

// Prepend 'i' to path vectors in 'pm'.
path_map_t prepend(const path_map_t& pm, int i)
{
    path_map_t result;

    for (path_map_t::const_iterator iter = pm.begin(); iter != pm.end(); ++iter) {
        path_map_t::key_type key = iter->first;

        key.insert(key.begin(), i);
        result.insert(path_map_t::value_type(key, iter->second));
    }
    return result;
}

vector<ipoint2> get_path_map_points(const path_map_t& pm)
{
    vector<ipoint2> pts;
    
    for (path_map_t::const_iterator pmiter = pm.begin(); pmiter != pm.end(); ++pmiter) 
        pts.push_back(pmiter->second.p);
    return pts;
}

void inhibit_path_map(path_map_t& pm, int maxsize)
{
    /*
        typedef vector<double> histogram_t;

        struct hpoint_t {
        ipoint2 p;
        histogram_t h;
        };

        typedef map<vector<int>, hpoint_t> path_map_t;

    */
    typedef map<ipoint2, path_map_t::const_iterator> map_t;

    if (pm.size() <= maxsize) 
        return;

    path_map_t result;
    map_t keepmap;
    vector<ipoint2> pts;

    pts.reserve(pm.size());
    for (path_map_t::const_iterator pmiter = pm.cbegin(); pmiter != pm.cend(); ++pmiter) {
        if (keepmap.insert(map_t::value_type(pmiter->second.p, pmiter)).second)
            pts.push_back(pmiter->second.p);
    }

    int radius = 1;

    while (pts.size() > maxsize) {
        pts = inhibit_point_set(pts, radius);
        ++radius;
    }

    for (auto piter = pts.begin(); piter != pts.end(); ++piter) {
        auto fiter = keepmap.find(*piter);
        result.insert(path_map_t::value_type(fiter->second->first, fiter->second->second));
    }
    swap(pm, result);
}

// Returns geometry of the part node; of each subpart only geometry of one (main) "variant" is used
// Implement get_geo_avg where of each subpart average of "variants" will be calculated!
void get_library_geo(path_map_t& geo, node* pn)
{
    int ename = atom("lyrSrc");
    int cname = atom("lyrCenter");

    geo.clear();
    foreach_neighbor(pn, ename, iter) {
        lib_data* nnd = (lib_data*)neighbor_node_data(iter);
        part_data_2* ned = (part_data_2*)neighbor_edge_data(iter);
        part_data_2a::geo_map_t::iterator giter = ned->geo.find(nnd->type);

        if (giter != ned->geo.end()) 
            path_map_union(geo, giter->second, ned->index);
    }
    foreach_neighbor(pn, cname, iter) {
        lib_data* nnd = (lib_data*)neighbor_node_data(iter);
        part_data_2a* ned = (part_data_2*)neighbor_edge_data(iter);
        part_data_2a::geo_map_t::iterator giter = ned->geo.find(nnd->type);

        if (giter != ned->geo.end()) 
            path_map_union(geo, giter->second, 0);
    }
}

void get_library_geo(path_map_t& geo, node* pn, const map<int, int>& imap)
{
    int ename = atom("lyrSrc");
    int cname = atom("lyrCenter");

    geo.clear();
    foreach_neighbor(pn, ename, iter) {
        part_data_2* ned = (part_data_2*)neighbor_edge_data(iter);
        map<int, int>::const_iterator miter = imap.find(ned->index);

        if (miter != imap.end()) {
            part_data_2a::geo_map_t::iterator giter = ned->geo.find(miter->second);

            if (giter != ned->geo.end()) // should always be true!
                path_map_union(geo, giter->second, ned->index);
            else {
                cout << "get_library_geo error" << endl;
                throw;
            }
        }
    }
    foreach_neighbor(pn, cname, iter) {
        part_data_2a* ned = (part_data_2a*)neighbor_edge_data(iter);
        map<int, int>::const_iterator miter = imap.find(0);

        if (miter != imap.end()) {
            part_data_2a::geo_map_t::iterator giter = ned->geo.find(miter->second);

            if (giter != ned->geo.end()) // should always be true!
                path_map_union(geo, giter->second, 0);
            else {
                cout << "get_library_geo error (center)" << endl;
                throw;
            }
        }
    }

}

vector<ipoint2> get_library_geo(node* pn)
{
    path_map_t pm;

    get_library_geo(pm, pn);
    return get_path_map_points(pm);
}

vector<ipoint2> get_library_geo_all(node* pn)
{
    int ename = atom("lyrSrc");
    int cname = atom("lyrCenter");

    vector<ipoint2> result;

    foreach_neighbor(pn, ename, iter) {
        part_data_2* ned = (part_data_2*)neighbor_edge_data(iter);
        
        for (part_data_2a::geo_map_t::iterator giter = ned->geo.begin(); giter != ned->geo.end(); ++giter) {
            for (path_map_t::iterator piter = giter->second.begin(); piter != giter->second.end(); ++piter) {
                result.push_back(piter->second.p);
            }
        }
    }
    foreach_neighbor(pn, cname, iter) {
        part_data_2a* ned = (part_data_2*)neighbor_edge_data(iter);

        for (part_data_2a::geo_map_t::iterator giter = ned->geo.begin(); giter != ned->geo.end(); ++giter) {
            for (path_map_t::iterator piter = giter->second.begin(); piter != giter->second.end(); ++piter) {
                result.push_back(piter->second.p);
            }
        }
    }
    return result;
}

// Returns a labeled set of points of library part. Inhibition with 'radius' is performed
// if radius > 0.
vector<pair<int, ipoint2> > get_library_geo_pieces(node* pn, int radius)
{
    typedef pair<int, ipoint2> item_t;

    int ename = atom("lyrSrc");
    int cname = atom("lyrCenter");

    vector<item_t> result;

    foreach_neighbor(pn, ename, iter) {
        part_data_2* ned = (part_data_2*)neighbor_edge_data(iter);
        vector<item_t> lpts;
        
        for (part_data_2a::geo_map_t::iterator giter = ned->geo.begin(); giter != ned->geo.end(); ++giter) {
            for (path_map_t::iterator piter = giter->second.begin(); piter != giter->second.end(); ++piter) {
                lpts.push_back(item_t(ned->index, piter->second.p));
            }
        }
        lpts = inhibit_point_set(lpts, radius);
        result.insert(result.end(), lpts.begin(), lpts.end());
    }
    foreach_neighbor(pn, cname, iter) {
        part_data_2a* ned = (part_data_2*)neighbor_edge_data(iter);
        vector<item_t> lpts;

        for (part_data_2a::geo_map_t::iterator giter = ned->geo.begin(); giter != ned->geo.end(); ++giter) {
            for (path_map_t::iterator piter = giter->second.begin(); piter != giter->second.end(); ++piter) {
                lpts.push_back(item_t(0, piter->second.p));
            }
        }
        lpts = inhibit_point_set(lpts, radius);
        result.insert(result.end(), lpts.begin(), lpts.end());
    }
    return result;
}

void get_dissimilar_cluster(vector<node*>& cluster, node* initp, const vector<node*>& parts, int maxsize)
{
    double benergy;
    vector<dpoint2> dvector;
    double scdistance;

    cluster.clear();
    cluster.push_back(initp);
    while ((int)cluster.size() < maxsize) { 
        double maxmin_scdistance = 0.0;
        node* maxminp = nullptr;

        for (vector<node*>::const_iterator piter = parts.begin(); piter != parts.end(); ++piter) {
            node* q = *piter;
            path_map_t qpm;
            double min_scdistance = 1e5;

            get_library_geo(qpm, q);
            for (vector<node*>::iterator citer = cluster.begin(); citer != cluster.end(); ++citer) {
                node* cp = *citer;
                path_map_t cpm;

                get_library_geo(cpm, cp);
                part_geometry_matching(benergy, dvector, scdistance, qpm, cpm, false);
                if (scdistance < min_scdistance) min_scdistance = scdistance;
            }
            if (min_scdistance > maxmin_scdistance) {
                maxmin_scdistance = min_scdistance;
                maxminp = q;
            }
        }
        if (maxminp != nullptr) cluster.push_back(maxminp);
        else break;
    }
}

void get_dissimilar_cluster(vector<node*>& cluster, const vector<node*>& parts, int maxsize)
{
    cluster.clear();

    if (parts.empty()) return;
    get_dissimilar_cluster(cluster, parts.front(), parts, maxsize);
}

bool is_sim_root(node* p)
{
    if (p == nullptr) return false;

    int name = atom("lyrSimRoot");

    foreach_neighbor (p, name, niter) {
        node* pn = neighbor_node(niter);

        if (pn != p) return false;
    }
    return true;
}

cv::Mat back_projection(const cv::Mat& data, const pca_data& pcd)
{
    cv::Mat coeffs;
    cv::Mat result;
    
    gemm(data - pcd.mean, pcd.eigenvectors, 1, cv::Mat(), 0, coeffs, cv::GEMM_2_T);
    gemm(coeffs, pcd.eigenvectors, 1, pcd.mean, 1, result, 0);
    return result;
}

// When calculating coefficients; use only only positions where mask is 1; 'mask' should be 0-1 matrix (row) of
// the same dimension as 'data'.
cv::Mat back_projection(const cv::Mat& data, const cv::Mat& mask, pca_data& pcd)
{
    cv::Mat coeffs;
    cv::Mat result;

    gemm((data - pcd.mean).mul(mask), pcd.eigenvectors, 1, cv::Mat(), 0, coeffs, cv::GEMM_2_T);
    gemm(coeffs, pcd.eigenvectors, 1, pcd.mean, 1, result, 0);
    return result;
}


// Value of exp(-1/2 (x^T Sigma x)) for 
// 'x': coefficients of expansion of data - mean with eigenvectors
// 'Sigma': diag(eigenvalues)
// 'eps': skip dimension with sigma < eps
double erf_value(const cv::Mat& data, const cv::Mat& mean, const cv::Mat& eigenvectors, 
    const cv::Mat& eigenvalues, double eps /* = 0.01 */)
{
    cv::Mat coeffs;
    double result = 0.0;

    cv::gemm(data - mean, eigenvectors, 1, cv::Mat(), 0, coeffs, cv::GEMM_2_T);
    for (int i = 0; i < coeffs.cols; ++i) {
        double x = coeffs.at<double>(0, i);
        double s2 = eigenvalues.at<double>(i, 0);

        if (s2 > eps) result += x*x/s2;
    }
    return exp(-result/2.0);
}

double erf_value(const cv::Mat& data, const pca_data& pcd, double eps /* = 0.01 */)
{
    return erf_value(data, pcd.mean, pcd.eigenvectors, pcd.eigenvalues, eps);
}

vector<dpoint2> sample_from_pca(const pca_data& pcad)
{
    if (pcad.eigenvalues.rows == 0) 
        return vector<dpoint2>();

    cv::Mat mresult = pcad.mean.clone();

    for (int i = 0; i < pcad.eigenvalues.rows; ++i) {
        double sigma = sqrt(pcad.eigenvalues.at<double>(i, 0));
        double r = normal_random(0, sigma);
        cv::Mat delta = pcad.eigenvectors.row(i);

        mresult = mresult + delta*r;
    }

    vector<dpoint2> result = partition(mresult);

    for (auto piter = result.begin(); piter != result.end(); ++piter) 
        (*piter) *= pcad.sizefactor;
    return result;
}

vector<dpoint2> get_vs_part_geo(node* p)
{
    vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);

    if (pd == nullptr) return vector<dpoint2>();
    else return partition(pd->pcad.mean);
}

// Returns "covering" of 'pts' with splines

// Converts library part 'p' to a spline; if 'p' is a vs part then mean is used, otherwise, the usual
// reconstruction is used
void part_to_splines(vector<cubic_spline>& splines, node* p)
{
    vector<ipoint2> pts = cast_vector<ipoint2, dpoint2>(get_vs_part_geo(p));
    
    if (pts.empty())
        pts = get_library_geo(p);
    splines_from_points(splines, pts, 5);
}

cv::Mat part_to_image(node* p) 
{
    vector<cubic_spline> partsplines;

    part_to_splines(partsplines, p);
    if (partsplines.empty() || partsplines.front().empty()) return cv::Mat1d();
    else return spline_image(partsplines, 2);
}

void comb_gen_rec(list<vector<int> >& result, const vector<int>& nvec, int n)
{
    if (nvec.empty() || n >= (int)nvec.size())
        return;

    int nvecmax = nvec[n];
    auto riter = result.begin();

    while (riter != result.end()) {
        for (int i = 0; i < nvecmax; ++i) {
            vector<int> newitem = *riter;

            newitem.push_back(i);
            result.insert(riter, newitem);
        }
        riter = result.erase(riter);
    }
    comb_gen_rec(result, nvec, n + 1);
}

// Returns in 'result' all possible combinations of (a_0, a_1, ..., a_{n-1})
// with a_i = 0, ..., nvec_i - 1
void comb_gen(list<vector<int> >& result, const vector<int>& nvec)
{
    result.clear();
    if (nvec.empty())
        return;
    result.push_back(vector<int>());
    comb_gen_rec(result, nvec, 0);
}

cv::Mat pca_data_to_image(const pca_data& pcad, int eval, int evsamples, 
    double dsigma /* = 1.5 */, double ssizethresh /* = 0.5 */)
{
    const int imgborder = 50;

    double sizefactor = (pcad.sizefactor <= 1.0) ? 50 : pcad.sizefactor;

    if (pcad.eigenvalues.rows == 0) 
        return cv::Mat();

    double sigma0 = sqrt(pcad.eigenvalues.at<double>(0, 0));
    list<vector<int> > comblist;
    int maxe = min(eval, pcad.eigenvalues.rows);
    vector<int> nvec(maxe);
    vector<ipoint2> greypts;

    for (int i = 0; i < maxe; ++i) {
        double sigmai = sqrt(pcad.eigenvalues.at<double>(i, 0));
        int n = max(1, (int)(evsamples*sigmai/sigma0));

        nvec[i] = (n % 2 == 0) ? n + 1 : n;
    }

    comb_gen(comblist, nvec);

    for (auto citer = comblist.begin(); citer != comblist.end(); ++citer) {
        cv::Mat mresult = pcad.mean.clone();

        for (int i = 0; i < citer->size(); ++i) {
            int j = citer->at(i);
            int maxj = nvec[i];

            if (maxj == 1) continue;

            double sigmai = sqrt(pcad.eigenvalues.at<double>(i, 0));
            double s = dsigma*sigmai;
            cv::Mat delta = pcad.eigenvectors.row(i);

            mresult = mresult + delta*(s*(2*j - maxj + 1)/(maxj - 1));
        }
        vector<dpoint2> pts = partition(mresult);

        for (auto piter = pts.begin(); piter != pts.end(); ++piter) {
            dpoint2 dp = (*piter) * sizefactor;
            greypts.push_back(ipoint2(int_round(dp.x), int_round(dp.y)));
        }
    }

    irectangle2 box = irectangle2::bounding_rectangle(greypts.begin(), greypts.end());
    cv::Mat3b result((int)box.y_dim() + 2*imgborder, (int)box.x_dim() + 2*imgborder, cv::Vec3b(0, 0, 0));

    for (auto piter = greypts.begin(); piter != greypts.end(); ++piter) {
        ipoint2 p = *piter - box.ll + ipoint2(imgborder, imgborder);

        cv::circle(result, cv::Point(p.x, p.y), 1, cv::Scalar_<uchar>(128, 128, 128), -1);
    }

    vector<dpoint2> rpts = partition(pcad.mean);

    for (auto piter = rpts.begin(); piter != rpts.end(); ++piter) {
        ipoint2 p = ipoint2(int_round(sizefactor*piter->x), int_round(sizefactor*piter->y)) - box.ll + ipoint2(imgborder, imgborder);
        cv::circle(result, cv::Point(p.x, p.y), 2, cv::Scalar_<uchar>(0, 0, 255), -1);
    }
    return result;
}

// See "part_to_images" for the explanation.
vector<cv::Mat> pca_data_to_images(const pca_data& pcad, int eval, int evsamples, 
    double dsigma /* = 1.5 */, double ssizethresh /* = 0.5 */, bool usesplines /* = true */)
{
    vector<cv::Mat> result;
    vector<cubic_spline> splines;
    double sizefactor = (pcad.sizefactor <= 1.0) ? 50 : pcad.sizefactor;

    if (pcad.eigenvalues.rows == 0) 
        return result;

    double sigma0 = sqrt(pcad.eigenvalues.at<double>(0, 0));
    list<vector<int> > comblist;
    int maxe = min(eval, pcad.eigenvalues.rows);
    vector<int> nvec(maxe);

    for (int i = 0; i < maxe; ++i) {
        double sigmai = sqrt(pcad.eigenvalues.at<double>(i, 0));
        int n = max(1, (int)(evsamples*sigmai/sigma0));

        nvec[i] = (n % 2 == 0) ? n + 1 : n;
    }

    comb_gen(comblist, nvec);

    for (auto citer = comblist.begin(); citer != comblist.end(); ++citer) {
        cv::Mat mresult = pcad.mean.clone();

        for (int i = 0; i < citer->size(); ++i) {
            int j = citer->at(i);
            int maxj = nvec[i];

            if (maxj == 1) continue;

            double sigmai = sqrt(pcad.eigenvalues.at<double>(i, 0));
            double s = dsigma*sigmai;
            cv::Mat delta = pcad.eigenvectors.row(i);

            mresult = mresult + delta*(s*(2*j - maxj + 1)/(maxj - 1));
        }
        vector<dpoint2> pts = partition(mresult);

        for (auto piter = pts.begin(); piter != pts.end(); ++piter) 
            (*piter) *= sizefactor;

        drectangle2 box = drectangle2::bounding_rectangle(pts.begin(), pts.end());
        int inhibition = int_round(sqrt(box.size2())/4.0);

        if (usesplines) {
            splines_from_points(splines, cast_vector<ipoint2, dpoint2>(pts), inhibition);
            keep_longest_splines(splines, ssizethresh);
            result.push_back(spline_image(splines, 2));
        } else {
            vector<ipoint2> ipts = cast_vector<ipoint2, dpoint2>(pts);

            ipts = inhibit_point_set(ipts, 5);
            result.push_back(point_image(cast_vector<dpoint2, ipoint2>(ipts)));
        }
    }
    return result;

}

// Returns images of part 'p'. If 'p' is ordinary part, only one image is returned;
// If 'p' is vs part, for each of the first 'eval' pca eigenvectors, 'evsamples' "equidistant"
// samples are generated, from [-2 sigma, 2 sigma]
vector<cv::Mat> part_to_images(node* p, int eval, int evsamples, double dsigma /* = 1.5 */)
{
    vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);

    if (pd != nullptr) {
        return pca_data_to_images(pd->pcad, eval, evsamples, dsigma);
    } else {
        vector<cv::Mat> result;
        vector<cubic_spline> splines;

        splines_from_points(splines, get_library_geo(p), 5);
        result.push_back(spline_image(splines, 2));
        return result;
    } 

}


