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

#include "utils/utils.h"
#include "utils/graphs/graph_utils.h"

#include "core/legacy/constants.h"
#include "core/legacy/layer_1_result.h"

#include "part_lib.h"

using namespace std;

// local definitions
///////////////////////////////////////////////////////////////////////////////

struct part_str_comparer : public binary_function<part_str, part_str, bool> {
    int dist;
  
    part_str_comparer(int d = 0) : dist(d) { }

    bool operator()(const part_str& p, const part_str& q) const
        { return p.type == q.type && abs(p.x - q.x) <= dist && abs(p.y - q.y) <= dist; }
};

// thresh_data
///////////////////////////////////////////////////////////////////////////////

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
            int edge = EdgeConnection::TO_LYR_VS_MEMBER;
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
            int cedge = EdgeConnection::TO_LYR_CENTER;
            int sedge = EdgeConnection::TO_LYR_SOURCE;

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
    : cmx(0), cmy(0), mask(), region(), lib_data(l, t)
{
    m.to_point_map(mask, 0.0);
    r.to_point_set(region, false);
}

part_data::part_data(const part_data& pd) :
    lib_data(pd),
    cmx(pd.cmx), cmy(pd.cmy), mask(pd.mask), region(pd.region)
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
        int name = EdgeConnection::TO_LYR_VS_MEMBER;

        foreach_neighbor (n, name, iter) {
            node* nn = neighbor_node(iter);
            part_data* nnd = (part_data*)neighbor_node_data(iter);

            nnd->make_mask(nn, libtype, library);
            merge_mask_maps(mask, nnd->mask, ipoint2::zero, 1.0, libtype);
        }
        return;
    }

    int srcname = EdgeConnection::TO_LYR_SOURCE;
    double factor = total_contraction(layer);

    part_data* nd = (part_data*)n->data;
    node* c = n->get_neighbor(EdgeConnection::TO_LYR_CENTER);
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

    int srcname = EdgeConnection::TO_LYR_SOURCE;
	double factor = total_contraction(layer);

	part_data* nd = (part_data*)n->data;
    ipoint2 cdelta(-nd->cmx, -nd->cmy);
    node* c = n->get_neighbor(EdgeConnection::TO_LYR_CENTER);

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

void part_data::reset_mask()
{
    if (layer != 0) mask.clear();
}

void part_data::reset_region()
{
    if (layer != 0) region.clear();
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
{  
    memset(layer_data, 0, MAX_LAYER_NUMBER * sizeof(streamable*));
	memset(contractions, 0, MAX_LAYER_NUMBER * sizeof(double));

}

part_lib::~part_lib()
{
    for (int i = 0; i < layer_count; ++i) {
        if (layer_data[i] != nullptr) delete layer_data[i];
    }
}

// read functions

part_lib* part_lib::read(const string& name) {
	part_lib* library = nullptr;
	try {
		// !!!!! WARNING !!!!!!
		// anything added here MUST also be added to hop_library(hop_blob& blob) constructor in hop.cpp
		// !!!!! WARNING !!!!!!
        library = (part_lib*)streamable::read(name);
        library->update_var_fields();

        cv::FileStorage fs(name + ".yml", cv::FileStorage::READ);
        if (fs.isOpened()) {
            library->read_cv_storable_data(*fs);
            fs.release();
        }
    } catch (...) {
        library = nullptr;
    }
	return library;
}


// member functions definitions

// get list of all subparts (TO_LYR_CENTER + TO_LYR_SOURCE) as list of part_str classes
void part_lib::get_structure(node* n, vector<part_str>& str)
{
    int src_name = EdgeConnection::TO_LYR_SOURCE;
    part_data* nd = (part_data*)n->data;

    str.clear();
    node* c = n->get_neighbor(EdgeConnection::TO_LYR_CENTER);
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

// find out if any existing parts are exactly the same as new part n
// but with also added some other subparts; add difference of the parts
// as forbidden parts for this new n part
// (i.e when part n is subset of part m then if during inference part m 
// will fire then also part n will fire - but we should not allow for that
// and we therefore add forbidden parts to n to prevent such situations)
void part_lib::extend_subparts(int lyr, node* n, niterator begin, niterator end)
{
    vector<part_str> nstr, istr, forb;
    list<part_str> diff, diff2;
    list<part_str>::iterator j;
    int forbindex = EdgeConnection::TO_LYR_FORBIDDEN;

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

// lyr = layer we want to add part to
// pd = pointer to part_data structure (node data)
// str = structure of the part (indices to previous level parts) -- (lyr. diff, type); str[0] is center
//      Note: lyr. difference < 0
// coo = coordinates of parts
// distr = distribution around the part
// tol = tolerance to check whether a part already exists
// the last three parameters are ignored when lay = 1
// returns type (index) of the part if it succeeds -type-1 of the part if it already exists
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
            add_edge_2(n, lyrm1[index], new part_data_2a(index), nullptr, EdgeConnection::TO_LYR_CENTER, EdgeConnection::TO_LYR_CENTER_BACK);
            add_edge_2(n, lyrm1[index], new part_data_2c(cooc[0]), EdgeConnection::TO_LYR_PREV);
        }

        // add edges to other subparts
        pdata.push_back(cdata);
        for (size_t i = 1; i < str.size(); ++i) {
            vector<node*>& lyrm1 = parts[lyr - 1 + str[i].first];
            index = str[i].second;
            if (index < (int)lyrm1.size()) {
                add_edge_2(n, lyrm1[index], new part_data_2(coon[i], *distr[i], index, i), EdgeConnection::TO_LYR_SOURCE);
                add_edge_2(n, lyrm1[index], new part_data_2c(cooc[i]), EdgeConnection::TO_LYR_PREV);
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
    int ename = EdgeConnection::TO_LYR_SOURCE;

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
            if (cldj.x == -1) add_edge_2(nnn, n, new part_data_2c(coo.x, coo.y), EdgeConnection::TO_LYR_CENTER_BACK);
        }
        add_edge_2(n, nn, pd2, EdgeConnection::TO_LYR_SOURCE);
        if (cld.x == -1) add_edge_2(nn, n, new part_data_2c(coo.x, coo.y), EdgeConnection::TO_LYR_CENTER_BACK);
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
    node* cp = p->get_neighbor(EdgeConnection::TO_LYR_CENTER_BACK);

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

    int srcedge = EdgeConnection::TO_LYR_SOURCE;
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
int part_lib::add_category(int lyr, const vector<int>& members, const string& name)
{
    if (lyr < 2) return -1;

    if (lyr > layer_count) ++layer_count;

    int z = lyr - 1;
    cpart_data* pd = new cpart_data(name == "" ? (string("cat_") + parts[z].size()) : name);
    node* n = add_node(pd, CATEGORY_NODE_ATTR);
    int to_center = EdgeConnection::TO_LYR_CENTER;
    int to_center_back = EdgeConnection::TO_LYR_CENTER_BACK;
    int to_prev = EdgeConnection::TO_LYR_PREV;

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
struct ispair {
	double d;
	node* n1;
	node* n2;

	ispair(double vd, node* vn1, node* vn2) : d(vd), n1(vn1), n2(vn2) { }

	bool operator<(const ispair& q) const { return d < q.d; }
	bool operator>(const ispair& q) const { return d > q.d; }
	bool operator==(const ispair& q) const { return d == q.d; }
};
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
    
    recurse(parts[layer_count - 1], EdgeConnection::TO_LYR_PREV, coll, result);
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
    int edge = EdgeConnection::TO_LYR_SIMILAR;
    int redge = EdgeConnection::TO_LYR_SIMROOT;
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
double part_lib::similarity(node* part1, node* part2, double rfactor, int depth, double M(double[], int, int))
{
    typedef pair<part_data_2*, node*> perm_pair_t;
    typedef edge_data_t<double> sim_data_t;

    part_data* pd1 = (part_data*)part1->data;
    part_data* pd2 = (part_data*)part2->data;
    int to_similar = EdgeConnection::TO_SIMILAR;
    
    if (pd1->layer == 0 && pd2->layer == 0) {
        matrix<double> mask1, mask2;
        ipoint2 c1, c2;

        c1 = pd1->get_mask(mask1, part1, this);
        c2 = pd2->get_mask(mask2, part2, this);
        return mask1.convolve_central(mask2, c2.x, c2.y, c1.x, c1.y);
    }

    int srcname = EdgeConnection::TO_LYR_SOURCE;
    int centername = EdgeConnection::TO_LYR_CENTER;
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
pair<node*, part_data_2*> part_lib::get_neighbor_pair(node* p, int index)
{
    typedef pair<node*, part_data_2*> result_t;

    if (p == nullptr) return result_t(nullptr, nullptr);
    
    int srcname = EdgeConnection::TO_LYR_SOURCE;

    foreach_neighbor(p, srcname, niter) {
        part_data_2* ped = (part_data_2*)neighbor_edge_data(niter);
        if (index == ped->index) return result_t(neighbor_node(niter), ped);
    }
    //if (pt.is_zero()) return p->get_neighbor(EdgeConnection::TO_LYR_CENTER);
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
    int edgename = EdgeConnection::TO_SIMILAR;

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
void part_lib::get_similar_types(int tlayer, int t, int clayer, int ct, int x, int y, double typethresh, set<int>& typeset)
{
    static int to_similar = EdgeConnection::TO_SIMILAR;
    static int to_srcm = EdgeConnection::TO_LYR_SRC_M;
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

double part_lib::get_thresh(int name, int lyr, int type, double defval)
{
    if (lyr < 1 || lyr > max_layer_index() || type < 0 || type >= (int)parts[lyr].size())
        return defval;
    lib_data* pd = (lib_data*)(parts[lyr][type])->data;
    return pd->get_thresh(name, defval);
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

bool is_vs_part(node* p)
{
    int sname = EdgeConnection::TO_LYR_SOURCE;
    int cname = EdgeConnection::TO_LYR_CENTER;

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
    int sname = EdgeConnection::TO_LYR_SOURCE;
    int cname = EdgeConnection::TO_LYR_CENTER;

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

// NOT USED ANY MORE ??
vector<pair<dpoint2, int> > random_sample_part(part_lib* library, node* n)
{
    typedef pair<dpoint2, int> result_item_t;

    int sname = EdgeConnection::TO_LYR_SOURCE;
    int cname = EdgeConnection::TO_LYR_CENTER;
    int vsmname = EdgeConnection::TO_LYR_VS_MEMBER;

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

        dpoint2 rnddelta = (dpoint2)sample_from_matrix_distribution(ed->distr);

        for (part_data_2::app_map_t::iterator aiter = ed->app.begin(); aiter != ed->app.end(); ++aiter) {
            avec.push_back(aiter->first);
        }

        node* varn = library->parts[nd->layer][avec[random_int(0, (int)avec.size())]];

        vector<result_item_t> pts = random_sample_part(library, varn);

        for (vector<result_item_t>::iterator ptiter = pts.begin(); ptiter != pts.end(); ++ptiter) {
            result_item_t p = *ptiter;
            dpoint2 delta = (cdelta + rnddelta + dpoint2(ed->x, ed->y))*factor;

            p.first += delta;
            result.push_back(p);
        }
    }
    
    return result;

}

// NOT USED ANY MORE ??
vector<dpoint2> random_sample_from_part(node* p)
{
    int sname = EdgeConnection::TO_LYR_SOURCE;
    int cname = EdgeConnection::TO_LYR_CENTER;
    int vsmname = EdgeConnection::TO_LYR_VS_MEMBER;

    vector<dpoint2> result;

    if (p == nullptr) return result;
    if (!is_vs_part(p)) return cast_vector<dpoint2, ipoint2>(get_library_geo(p));

    vector<pair<int, ipoint2> > gpieces = get_library_geo_pieces(p, 3);
    part_data* pd = (part_data*)p->data;

    if (!has_vs_children(p)) {
        vector<dpoint2> result = random_sample_from_vs_part(dynamic_cast<vs_part_data*>(p->data));
        
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
            partsample[i] += cdelta;
        }

        vector<dpoint2> dpts = cast_vector<dpoint2, ipoint2>(pts);

        result.insert(result.end(), partsample.begin(), partsample.end());
    }
    return result;
}

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

    int sname = EdgeConnection::TO_LYR_SOURCE;
    int cname = EdgeConnection::TO_LYR_CENTER;

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
    img result = one_row ? img::concat_linear_fw(images) : img::concat(images);
    result.save_normalized(fname);
}

void part_lib::save_all_sc_mma(const char* fname, int lyr, int min /* = 0 */, int max /* = -1 */)
{
    int sname = EdgeConnection::TO_LYR_SOURCE;
    int cname = EdgeConnection::TO_LYR_CENTER;

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
void part_lib::save_part(const char* fname, int lyr, int part, double threshold /* = 0.0 */)
{
    if (lyr < 1 || lyr > layer_count || part < 0 || part >= layer_size(lyr - 1))
        return;

    vector<img> images;
    node* p = parts[lyr - 1][part];
    int edgename = EdgeConnection::TO_SIMILAR;

    if (!updated) 
        part_data::set_contractions(contractions, MAX_LAYER_NUMBER);
    images.push_back(((part_data*)(p->data))->get_image(p, attr, this));
    foreach_neighbor(p, edgename, iter) {
        node* q = (node*)iter->second.first;
        edge_data_t<double>* qd = (edge_data_t<double>*)iter->second.second;

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
    int edgename = EdgeConnection::TO_LYR_SOURCE;
    int edgenameM = EdgeConnection::TO_LYR_SRC_M;

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
    int edgename = EdgeConnection::TO_LYR_SOURCE;
    int cedgename = EdgeConnection::TO_LYR_CENTER;
    int c1edgename = EdgeConnection::TO_LYR_CENTER1;
    int vsmname = EdgeConnection::TO_LYR_VS_MEMBER;

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

void part_lib::display_layer_info(int layer) {

	if (layer >= (int)this->parts.size()) { 
		cout  << "Library does not have layer '" << layer << "'" << endl;
		return;
	}

	vector<node*>& parts = this->parts[layer];
    int to_part = EdgeConnection::TO_LYR_SOURCE;
    int to_center = EdgeConnection::TO_LYR_CENTER;
    int to_forbidden = EdgeConnection::TO_LYR_FORBIDDEN;
    int to_prev = EdgeConnection::TO_LYR_PREV;

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
void part_lib::update_var_fields() {

    int srcname = EdgeConnection::TO_LYR_SOURCE;
    int centername = EdgeConnection::TO_LYR_CENTER;

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

void part_lib::edge_paths(vector<edge_path>& paths, node* p)
{
    int tosrc = EdgeConnection::TO_LYR_SOURCE;
    int tocenter = EdgeConnection::TO_LYR_CENTER;

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

vector<int> part_lib::get_root_parts_map(int layer)
{
    int rootedge = EdgeConnection::TO_LYR_SIMROOT;
    int simedge = EdgeConnection::TO_LYR_SIMILAR;

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
	
		cv::PCA pca(pcadata, cv::Mat(), CV_PCA_DATA_AS_ROW, maxcomp);

        return pca_data(pca.mean, pca.eigenvectors, pca.eigenvalues, factorsum/n);
    }
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

vector<ipoint2> get_path_map_points(const path_map_t& pm)
{
    vector<ipoint2> pts;
    
    for (path_map_t::const_iterator pmiter = pm.begin(); pmiter != pm.end(); ++pmiter) 
        pts.push_back(pmiter->second.p);
    return pts;
}

void inhibit_path_map(path_map_t& pm, int maxsize)
{
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
    int ename = EdgeConnection::TO_LYR_SOURCE;
    int cname = EdgeConnection::TO_LYR_CENTER;

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
    int ename = EdgeConnection::TO_LYR_SOURCE;
    int cname = EdgeConnection::TO_LYR_CENTER;

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
// Returns a labeled set of points of library part. Inhibition with 'radius' is performed
// if radius > 0.
vector<pair<int, ipoint2> > get_library_geo_pieces(node* pn, int radius)
{
    typedef pair<int, ipoint2> item_t;

    int ename = EdgeConnection::TO_LYR_SOURCE;
    int cname = EdgeConnection::TO_LYR_CENTER;

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

    int name = EdgeConnection::TO_LYR_SIMROOT;

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