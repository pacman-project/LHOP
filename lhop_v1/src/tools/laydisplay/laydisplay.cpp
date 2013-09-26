// laydisplay.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cctype>
#include "../../utils/atom.h"
#include "../../utils/hopmath.h"
#include "../../graphs/img_graph.h"
#include "../../layers/layer_1.h"
#include "../../utils/convert.h"
#include "../../utils/matrix.h"
#include "../../utils/utils.h"
#include "../../layers/initialization.h"
#include "../../graphs/graph_utils.h"
#include "../../layers/layer_1_result.h"
#include <highgui.h>
#include "../../utils/ocv.h"
#include <tuple>

#ifdef WIN32
#define TOLOWER std::tolower
#define TOUPPER std::toupper
#else
#define TOLOWER ::tolower
#define TOUPPER ::toupper
#endif

// local structures
///////////////////////////////////////////////////////////////////////////////

struct annot_struct {
    int type;
    int x_size, y_size;
    int x0, y0;
    int x1, y1;
};

/* class svm_exception { 
public:
	string error;
	svm_exception(const string& e = "") : error(e) { }
	void print() { cout << e << endl; }
}; */

typedef pair<string, annot_struct> map_pair;

// local functions
////////////////////

void bin_write_double(ostream& os, double d)
{
    os.write((const char*)(&d), sizeof(d));
}


void bin_write_int(ostream& os, int i)
{
    os.write((const char*)(&i), sizeof(i));
}

void write_sparse_histogram(ostream& os, const vector<double>& h)
{
    typedef pair<int, double> item_t;
    //vector<item_t> sh;

    int zeros = count(h.begin(), h.end(), 0);
    int nonzeros = (int)h.size() - zeros;

    bin_write_int(os, (int)h.size());
    bin_write_int(os, nonzeros);    
    //sh.resize(nonzeros);
    for (int i = 0, c = 0; i < (int)h.size(); ++i) {
        if (h[i] != 0.0) {
            bin_write_int(os, i);
            bin_write_double(os, h[i]);
            //sh[c].first = i;
            //sh[c].second = h[i];
            //++c;
        }
    }
    //os.write((const char*)&(sh.at(0)), nonzeros*(sizeof(item_t)));
}

//void recurse_directory(set<string>& result, const string& pattern)
//{
//    WIN32_FIND_DATA FindFileData;
//    HANDLE hFind;
//
//    result.clear();
//    hFind = FindFirstFile(pattern.c_str(), &FindFileData);
//    if (hFind == INVALID_HANDLE_VALUE) return;
//    result.insert(FindFileData.cFileName);
//    while (FindNextFile(hFind, &FindFileData)) {
//        result.insert(FindFileData.cFileName);
//    }
//    FindClose(hFind);
//}


// local functions

void read_annotations(const string& annotname, map<string, annot_struct>& annotations)
{
    ifstream f;
    string s, s1, s2;
    annot_struct as;
    size_t i;
    
    f.open(annotname.c_str());

    annotations.clear();
    while (!f.fail()) {
        f >> s1 >> s2;
        f >> as.type;
        f >> as.y0 >> as.x0;
        f >> as.y1 >> as.x1;
        f >> as.y_size >> as.x_size;

		 --as.x0; --as.y0;
		 --as.x1; --as.y1;
        if ((i = s2.rfind(".")) != string::npos) s2 = s2.substr(0, i);

        transform(s.begin(), s.end(), s.begin(), TOLOWER);
        annotations.insert(map_pair(s1 + "_" + s2, as));
    }
    f.close();
}

void make_merge(vector<int>& typemap1, part_lib* library, int layer, double merget, int mergetol)
{
    vector<node*>& parts = library->parts[layer];
    int size = (int)parts.size();
    vector<int> origtypes(size, 0);
    vector<int> newtypes(size, 0);
    part_data* pd;
    int i;
    
    i = 0;
    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        pd = (part_data*)(*iter)->data;
        origtypes[i++] = pd->type;
    }

    cout << "Merging parts..." << endl;
    int mergecount = library->merge_part_types(layer + 1, merget, mergetol);
    cout << endl << "Number of merged parts: " << mergecount << endl;

    i = 0;
    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        pd = (part_data*)(*iter)->data;
        newtypes[i++] = pd->type;
    }
    typemap1.resize(parts.size(), 0);
    for (i = 0; i < size; ++i) {
        typemap1[origtypes[i]] = newtypes[i];
    }

}

int get_part_map(vector<int>& pmap, const vector<node*>& parts)
{
    set<int> types;

    for (vector<node*>::const_iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        part_data* pd = (part_data*)(*iter)->data;
        types.insert(pd->type);
    }

    int newi = 0;

    pmap.resize(parts.size(), 0);
    for (set<int>::iterator iter = types.begin(); iter != types.end(); ++iter) {
        pmap[*iter] = newi++;
    }
    return (int)types.size();
}
				
int append_prev_layer(vector<int>& typemap1, vector<int>& typemap, int layersize, part_lib* library, int layer)
{
    int prevlayersize = (int)library->parts[layer - 1].size();
    int origlayersize = (int)library->parts[layer].size();

    for (int i = 0; i < prevlayersize; ++i) {
        typemap1.push_back(origlayersize + i);
        typemap.push_back(layersize + i);
    }
    return layersize + prevlayersize;
}



int cell_index(const vector<int>& boundaries, int pos)
{
    int size = (int)boundaries.size() - 1;
    for (int i = 1; i < size; ++i) {
        if (pos <= boundaries[i]) return i - 1;
    }
    return size - 1;
}

void read_rectangles(vector<rectangle>& rectlist, const string& fname)
{
	if (fname.size() == 0) return;

    ifstream f;
    rectangle r;

    f.open(fname.c_str());
    while (!f.fail()) {
        f >> r;
        if (!f.fail()) rectlist.push_back(r);
    }
}

// method definitions



int index_from_name(const string& name)
{
    size_t pos = name.rfind("_");

    if (pos == string::npos || pos < 4) return 0;
    else {
        string s = name.substr(pos - 4, 4);
        return Convert::ToInteger<int>(s.c_str(), 4);
    }
}





struct dclassifier_data {
	vector<double> mu, vx;
	vector<double> w;
	vector<int> selection;
	vector<int> valid;
};


void part_paths(const list<pair<node*, vector<int> > >& plist, list<vector<int> >& result, int minlayer)
{
    int edge_name = atom("lyrSrc");
    list<pair<node*, vector<int> > > neighbors;
    
    for (auto iter = plist.begin(); iter != plist.end(); ++iter) {
		lib_data* pd = (lib_data*)iter->first;

		if (pd->layer < minlayer) 
			continue;

        node* cp = iter->first->get_neighbor(atom("lyrCenter"));

        if (cp != nullptr) {
            result.push_back(iter->second);
            result.back().push_back(0);
            neighbors.push_back(pair<node*, vector<int> >(cp, result.back()));
        }
                
        foreach_neighbor(iter->first, edge_name, niter) {
            part_data_2* pd = (part_data_2*)neighbor_edge_data(niter);

            result.push_back(iter->second);
            result.back().push_back(pd->index);
            neighbors.push_back(pair<node*, vector<int> >(neighbor_node(niter), result.back()));
        }
    }
    if (neighbors.size() > 0)
        part_paths(neighbors, result, minlayer);
}

void part_paths(node* p, list<vector<int> >& result, int minlayer)
{
    list<pair<node*, vector<int> > > plist;

    plist.push_back(pair<node*, vector<int> >(p, vector<int>()));
    part_paths(plist, result, minlayer);
}

void part_path_map(node* p, map<vector<int>, int>& result, int minlayer)
{
    typedef map<vector<int>, int> result_t;

    int i = 0;
    list<vector<int> > l;

    result.clear();
    part_paths(p, l, minlayer);
    for (auto iter = l.begin(); iter != l.end(); ++iter) {
        result.insert(result_t::value_type(*iter, i));
        ++i;
    }
}


vector<int> round_vector(const vector<double>& dv)
{
	vector<int> result(dv.size());

	for (size_t i = 0; i < dv.size(); ++i) {
		result[i] = int_round(dv[i]);
	}
	return result;
}

template<class T> vector<T> extract_from_binary(const vector<T>& v, const vector<int>& binv)
{
	vector<T> result;
	size_t maxsize = min(v.size(), binv.size());

	for (auto i = 0; i < maxsize; ++i) {
		if (binv[i] != 0) result.push_back(v[i]);
	}
	return result;
}

template<class T> vector<T> extract_from_indices(const vector<T>& v, const vector<int>& indv)
{
	vector<T> result;
	size_t maxv = v.size(), maxsize = indv.size();

	for (auto i = 0; i < maxsize; ++i) {
		const int& index = indv[i];

		if (index >= 0 && index < maxv) result.push_back(v[index]);
	}
	return result;
}

// Return neighbors of 'n' for edge 'name', 'n' is always first, if it is a nieighbor
vector<node*> neighbor_vector(node* n, int name)
{
	vector<node*> result;

	foreach_neighbor(n, name, niter) {
		if (neighbor_node(niter) == n) 
			result.insert(result.begin(), n);
		else
			result.push_back(neighbor_node(niter));
	}
	return result;
}

// Return neighbors of 'n' for edge 'name', 'n' is always added 
// and is added first
vector<node*> neighbor_vector_n(node* n, int name)
{
	vector<node*> result;

	result.push_back(n);
	foreach_neighbor(n, name, niter) {
		if (neighbor_node(niter) != n) 
			result.push_back(neighbor_node(niter));
	}
	return result;
}


vector<pair<node*, vector<int> > > app_vector(part_lib* library, node* p, int index)
{
	int esrc = atom("lyrSrc");
	vector<pair<node*, vector<int> > > result;

	foreach_neighbor(p, esrc, pniter) {
		node* pn = neighbor_node(pniter);
		lib_data* pnd = (lib_data*)pn->data;
		part_data_2* ed = (part_data_2*)neighbor_edge_data(pniter);
		
		if (ed->index == index) {
			for (auto aiter = ed->app.begin(); aiter != ed->app.end(); ++aiter) {
				if (aiter->first == pnd->type) {
					result.insert(result.begin(), pair<node*, vector<int> >(pn, aiter->second.first));
				} else {
					result.push_back(pair<node*, vector<int> >(library->parts[pnd->layer][aiter->first], aiter->second.first));
				}
			}
			return result;
		}
	}

	if (index == 0) {
		edge_pair pnp = p->get_neighbor_pair(atom("lyrCenter"));
		node* pn = pnp.first;
		lib_data* pnd = (lib_data*)pn->data;
		part_data_2* ed = (part_data_2*)pnp.second;
		
		for (auto aiter = ed->app.begin(); aiter != ed->app.end(); ++aiter) {
			if (aiter->first == pnd->type) {
				result.insert(result.begin(), pair<node*, vector<int> >(pn, aiter->second.first));
			} else {
				result.push_back(pair<node*, vector<int> >(library->parts[pnd->layer][aiter->first], aiter->second.first));
			}
		}
		return result;
	}
	return result;
}


// (p*q)[i] = p[q[i]]
vector<int> permutation_product(const vector<int>& p, const vector<int>& q)
{
	if (p.empty()) return q;
	if (q.empty()) return p;

	vector<int> result(q.size());

	for (int i = 0; i < q.size(); ++i) {
		if (q[i] < p.size()) 
			result[i] = p[q[i]];
		else {
			cout << "permutation_product: Permutation size mismatch." << endl;
			throw;
		}
	}
	return result;
}

typedef edge_data_t<vector<int> > perm_edge_data;


// Return value: if data is consistent with 
bool library_permutation(vector<int>& result, part_lib* library, int srclayer, int srcpart, int destpart, int index)
{
	result.clear();

	if (srclayer < 0 || srclayer > library->max_layer_index() || srcpart < 0 || srcpart >= library->layer_size(srclayer))
		return false;

    int src_name = atom("lyrSrc");
	int center_name = atom("lyrCenter");
	int sim_name = atom("lyrSimilar");
	node* p = library->parts[srclayer][srcpart];
	vector<node*> pnbhood = neighbor_vector_n(p, sim_name);

	/// DBG DBG DBG 
	//if (srclayer == 3 && srcpart == 0 && destpart == 20) {
	//	for (auto dbgiter = pnbhood.begin(); dbgiter != pnbhood.end(); ++dbgiter) 
	//		cout <<  ((lib_data*)(*dbgiter)->data)->type << ' ';
	//	cout << endl;
	//}
	/// DBG DBG DBG

	for (auto piter = pnbhood.begin(); piter != pnbhood.end(); ++piter) {
		node* sp = *piter;
		perm_edge_data* sped = (perm_edge_data*)p->get_edge_data(sp, sim_name);
		//vector<int> pperm = sped == nullptr ? vector<int>() : sped->data ;	// not needed
		vector<pair<node*, vector<int> > > apphood = app_vector(library, sp, index);

	/// DBG DBG DBG 
	//if (srclayer == 3 && srcpart == 0 && destpart == 20) {
	//	for (auto dbgiter = apphood.begin(); dbgiter != apphood.end(); ++dbgiter) 
	//		cout <<  ((lib_data*)(dbgiter->first)->data)->type << ' ';
	//	cout << endl;
	//}
	/// DBG DBG DBG

		for (auto aiter = apphood.begin(); aiter != apphood.end(); ++aiter) {
			node* spn = aiter->first;
			const vector<int>& aperm = aiter->second;
			vector<node*> spnnbhood = neighbor_vector_n(spn, sim_name);
			
	/// DBG DBG DBG 
	//if (srclayer == 3 && srcpart == 0 && destpart == 20) {
	//	for (auto dbgiter = spnnbhood.begin(); dbgiter != spnnbhood.end(); ++dbgiter) 
	//		cout <<  ((lib_data*)(*dbgiter)->data)->type << ' ';
	//	cout << endl;
	//}
	/// DBG DBG DBG

			for (auto spiter = spnnbhood.begin(); spiter != spnnbhood.end(); ++spiter) {
				node* sspn = *spiter;
				lib_data* sspnd = (lib_data*)sspn->data;
				perm_edge_data* spned = (perm_edge_data*)spn->get_edge_data(sspn, sim_name);
				vector<int> spperm = spned == nullptr ? vector<int>() : spned->data;

				if (sspnd->type == destpart) {
					result = permutation_product(aperm, spperm);
					return true;
				}
			}
		}
	}
	return false;
}

void subgraph_paths(part_lib* library, const list<tuple<node* , vector<int>, vector<int> > >& plist, map<node*, vector<int> >& result, int minlayer)
{
    typedef map<node*, vector<int> > result_t;
	typedef tuple<node* , vector<int>, vector<int> > list_item_t;

    int edge_name = atom("toPrevLayer");
    list<list_item_t> neighbors;
    
    for (auto iter = plist.begin(); iter != plist.end(); ++iter) {
		node* n = get<0>(*iter);
		layer1_data* nd = (layer1_data*)n->data;

        foreach_neighbor(n, edge_name, niter) {
			layer1_data* nnd = (layer1_data*)neighbor_node_data(niter);

			if (nnd->m < minlayer) continue;

            edge_data_name* ed = (edge_data_name*)neighbor_edge_data(niter);
            vector<int> r = get<2>(*iter);
			vector<int> permp = get<1>(*iter);

			if (permp.empty()) r.push_back(ed->index);
			else r.push_back(permp[ed->index]);
            result.insert(result_t::value_type(neighbor_node(niter), r));

			vector<int> newp;

			if (nd->z > 1) {	// newp is identity for layers < 2
				if (!library_permutation(newp, library, nd->z, nd->m, nnd->m, ed->index)) {
					cout << "library_permutation_ptr returned false." << endl;
					throw;
				}
			}

            neighbors.push_back(list_item_t(neighbor_node(niter), newp, r));
        }
    }
    if (neighbors.size() > 0) 
        subgraph_paths(library, neighbors, result, minlayer);
}

void subgraph_paths(part_lib* library, node* n, map<node*, vector<int> >& result, int minlayer)
{
	typedef tuple<node* , vector<int>, vector<int> > list_item_t;

    list<list_item_t> neighbors;

    neighbors.push_back(list_item_t(n, vector<int>(), vector<int>()));
    subgraph_paths(library, neighbors, result, minlayer);
}


void result_paths(const list<pair<node*, vector<int> > >& plist, list<vector<int> >& result)
{
    int edge_name = atom("lyrSrc");
    list<pair<node*, vector<int> > > neighbors;
    
    for (auto iter = plist.begin(); iter != plist.end(); ++iter) {
        node* cp = iter->first->get_neighbor(atom("lyrCenter"));

        if (cp != nullptr) {
            result.push_back(iter->second);
            result.back().push_back(0);
            neighbors.push_back(pair<node*, vector<int> >(cp, result.back()));
        }
                
        foreach_neighbor(iter->first, edge_name, niter) {
            part_data_2* pd = (part_data_2*)neighbor_edge_data(niter);

            result.push_back(iter->second);
            result.back().push_back(pd->index);
            neighbors.push_back(pair<node*, vector<int> >(neighbor_node(niter), result.back()));
        }
    }
    if (neighbors.size() > 0)
        result_paths(neighbors, result);
}

void make_path_dhistogram(vector<double>& hist, irectangle2& bbox, const map<vector<int>, int>& typemap, 
    part_lib* library, layer1_result* res, node* n, int depth)
{
    int edgename = atom("toPrevLayer");
    layer1_data* nd = (layer1_data*)n->data;
    map<node*, vector<int> > parsenodes;
	int histsize = 0;

	if (depth <= 0) 
		histsize = typemap.size();
	else {
		for (auto tmiter = typemap.begin(); tmiter != typemap.end(); ++tmiter) 
			if (tmiter->first.size() <= depth) histsize = max(histsize, tmiter->second + 1);
	}

    hist.resize(histsize);
    fill(hist.begin(), hist.end(), 0.0);
    bbox.invalidate();
    subgraph_paths(library, n, parsenodes, 0);
    for (auto piter = parsenodes.begin(); piter != parsenodes.end(); ++piter) {
        node* nn = piter->first;
        layer1_data* nnd = (layer1_data*)nn->data;
        auto fiter = typemap.find(piter->second);

        if (fiter != typemap.end()) {
			if (fiter->second < histsize)
				hist[fiter->second] += nnd->vval();
            if (nnd->z == 0)
                bbox.eat(nnd->x, nnd->y);
        }
    }
}


void make_dhistogram(vector<double>& hist, irectangle2& bbox, const vector<vector<int> >& typemap, 
    int hsize, layer1_result* res, node* n)
{
    int edgename = atom("toPrevLayer");
    layer1_data* nd = (layer1_data*)n->data;
    set<node*> parsenodes;

    hist.resize(hsize);
    fill(hist.begin(), hist.end(), 0.0);
    bbox.invalidate();
    res->subgraph_from_node(n, edgename, parsenodes);
    for (auto piter = parsenodes.begin(); piter != parsenodes.end(); ++piter) {
        node* nn = *piter;
        layer1_data* nnd = (layer1_data*)nn->data;

        hist[typemap[nnd->z][nnd->m]] += nnd->vval();
        if (nnd->z == 0)
            bbox.eat(nnd->x, nnd->y);
    }
}

pair<int, double> dhistogram_classification(irectangle2& box, layer1_result* res, node* n, const vector<vector<int> >& typemap, int hsize,
	const dclassifier_data& dcd)
{
	vector<double> h;

	make_dhistogram(h, box, typemap, hsize, res, n);
    //make_path_dhistogram(hd.h, box, typemap[nnd->m], res.get(), nn);

	// ---------------- DEBUG --------------
	//static int index = 0;
	//layer1_data* nd = (layer1_data*)n->data;
	//string s = "c:\\work\\tmp\\h-";
	//s = s + index;
	//s = s + ".hist";
	//ofstream os(s.c_str(), ios_base::binary | ios_base::out);
 //   
	//bin_write_int(os, 1);
	//bin_write_int(os, 0);
	//bin_write_double(os, nd->vval());
	//bin_write_double(os, nd->r(S_RESPONSE));
	//write_sparse_histogram(os, h);
	//
	//os.close();

 //   ofstream os1((s + ".model").c_str(), ios_base::binary | ios_base::out);
 //   bin_write_int(os1, 1);
	//bin_write_int(os1, nd->m);
	//bin_write_double(os1, 0.0);
	//bin_write_double(os1, 0.0);
 //   write_sparse_histogram(os1, dcd.w);

 //   os1.close();

	//---------------- DEBUG --------------
	//h = extract_from_indices(h, dcd.selection);
	h = extract_from_binary(h, dcd.valid);

	for (int i = 0; i < h.size(); ++i) {
		double& hi = h[i];

		hi = (hi - dcd.mu[i])/dcd.vx[i];
	}
	h.push_back(1.0);

	double P = 0.0;
	int label;

	for (int i = 0; i < h.size(); ++i) 
		P += h[i]*dcd.w[i];
	label = (P > 0) ? 1 : 0;
	P = 1.0/(1 + exp(-P));

	// ------- DEBUG --------
	//s = "c:\\work\\tmp\\res-";
	//s = s + index;
	//s = s + ".hist";
	//os.open(s.c_str(), ios_base::binary | ios_base::out);
	//if (label == 0) { bin_write_double(os, 0.0); bin_write_double(os, 1 - P); }
	//else { bin_write_double(os, 1.0); bin_write_double(os, P); }
	//os.close();
	//++index;
	// ------- DEBUG --------

    //if (index > 3) throw;

	if (label == 0) return pair<int, double>(0, 1.0 - P);
	else return pair<int, double>(1, P);
}

bool load_matlab(vector<double>& v, const string& file)
{
	ifstream is(file.c_str());
	double d;

	v.clear();
	if (!is.good()) return false;
	while (is.good()) {
		is >> d;
		if (is.good()) v.push_back(d);
	}
	return true;
}



bool load_dclassifier_data(dclassifier_data& cld, const string& directory, const string& category)
{
	vector<double> valid, selection;
	string dir = directory;

	end_dir(dir);
	load_matlab(cld.mu, dir + category + "-mu.txt");
	load_matlab(valid, dir + category + "-valid.txt");
	load_matlab(cld.vx, dir + category + "-vx.txt");
	load_matlab(cld.w, dir + category + "-w.txt");
	load_matlab(selection, dir + category + "-selection.txt");

	if (cld.mu.empty() || valid.empty() || cld.vx.empty() || cld.w.empty() || selection.empty()) 
		return false;

	cld.valid = round_vector(valid);
	cld.selection = round_vector(selection);
	cld.mu = extract_from_binary(cld.mu, cld.valid);
	cld.vx = extract_from_binary(cld.vx, cld.valid);

	for_each(cld.selection.begin(), cld.selection.end(), [](int& i) { --i; });	// base 1 -> base 0
    return true;

}

void load_dclassifier_data(vector<dclassifier_data>& cldv, const string& directory, const string& category)
{
    int mi = 1;

    while (true) {
        dclassifier_data cld;

        if (!load_dclassifier_data(cld, directory, category + mi))
            break;
        cldv.push_back(cld);
        ++mi;
    }
}



void write_hoc_vector(const svm2::vector_t& v, const string& fname)
{
	ofstream os(fname.c_str());
	
	for (int i = 0; i < v.size(); ++i) {
		if (i > 0) os << ',';
		os << v[i];
	}
	os.close();
}

void write_hoc_vector(const cv::Mat& m, const string& fname)
{
	ofstream os(fname.c_str());

	cv_write_csv<float>(os, m);
	os.close();
}

/*
double check_hoc(const svm_data& svmd, part_lib* library, layer1_result* res, node* n, const vector<int> layers,
    const K_bin& bin, int box_grow_factor)
{
    int ename = atom("toPrevLayer");
    int linkname = atom("toLayer0");
    layer1_data* nd = (layer1_data*)n->data;
    set<node*> nset;
    set<ipoint2> pset;
    svm2::vector_t v;
    
    res->recurse_and_link(n, ename, linkname, nset);
    irectangle2 bbox = bounding_rectangle_of_nodes(nset.begin(), nset.end());
    
    bbox.grow(box_grow_factor);
    
    v.clear();
    for (int i = 0; i < layers.size(); ++i) {
        int l = layers[i];
        svm2::vector_t v1;

        fill_hoc_feature_vector(v1, res, l, library->layer_size(l), bin, bbox.center(), bbox);
        //normalize_histogram(v1);
        v.insert(v.end(), v1.begin(), v1.end());
    }
#ifdef _DEBUG || DEBUG
	write_hoc_vector(v, "c:\\work\\tmp\\v0.txt");
#endif
    cv::Mat test(1, v.size(), CV_32FC1, cv::Scalar_<float>(0.0));

    for (int i = 0; i < v.size(); ++i) {
        float& f = test.at<float>(0, i) = (v[i] - svmd.mean.at<float>(0, i));
        
        if (svmd.sigma.at<float>(0, i) > 0) 
            f /= svmd.sigma.at<float>(0, i);
    }
#ifdef _DEBUG || DEBUG
	write_hoc_vector(test, "c:\\work\\tmp\\v.txt");
#endif

    return svmd.svm->predict(test, true);
}

*/

double check_hoc(const svm_data& svmd, part_lib* library, layer1_result* res, node* n, const vector<int> layers,
    const K_bin& bin, int box_grow_factor)
{
    int ename = atom("toPrevLayer");
    int linkname = atom("toLayer0");
    layer1_data* nd = (layer1_data*)n->data;
    set<node*> nset;
    set<ipoint2> pset;
    svm2::vector_t v;
    
    res->recurse_and_link(n, ename, linkname, nset);
    irectangle2 bbox = bounding_rectangle_of_nodes(nset.begin(), nset.end());
    
    bbox.grow(box_grow_factor);
    
    v.clear();
    for (int i = 0; i < layers.size(); ++i) {
        int l = layers[i];
        svm2::vector_t v1;

        fill_hoc_feature_vector(v1, res, l, library->layer_size(l), bin, bbox.center(), bbox);
        //normalize_histogram(v1);
        v.insert(v.end(), v1.begin(), v1.end());
    }
#ifdef _DEBUG || DEBUG
	write_hoc_vector(v, "c:\\work\\tmp\\v0.txt");
#endif
    cv::Mat test(1, v.size(), CV_32FC1, cv::Scalar_<float>(0.0));

    for (int i = 0; i < v.size(); ++i) {
        float& f = test.at<float>(0, i) = (v[i] - svmd.mean.at<float>(0, i));
        
        if (svmd.sigma.at<float>(0, i) > 0) 
            f /= svmd.sigma.at<float>(0, i);
    }
#ifdef _DEBUG || DEBUG
	write_hoc_vector(test, "c:\\work\\tmp\\v.txt");
#endif

    return svmd.svm->predict(test, true);
}

double check_hoc(const svm_data& svmd, part_lib* library, layer1_result* res, const vector<node*>& nv, const vector<int> layers,
    const K_bin& bin, int box_grow_factor)
{
    svm2::vector_t v;
	irectangle2 bbox = bounding_rectangle_of_projection(res, nv);

	bbox.grow(box_grow_factor);
    v.clear();
    for (int i = 0; i < layers.size(); ++i) {
        int l = layers[i];
        svm2::vector_t v1;

        fill_hoc_feature_vector(v1, res, l, library->layer_size(l), bin, bbox.center(), bbox);
        //normalize_histogram(v1);
        v.insert(v.end(), v1.begin(), v1.end());
    }
#ifdef _DEBUG || DEBUG
	write_hoc_vector(v, "c:\\work\\tmp\\v0.txt");
#endif
    cv::Mat test(1, v.size(), CV_32FC1, cv::Scalar_<float>(0.0));

    for (int i = 0; i < v.size(); ++i) {
        float& f = test.at<float>(0, i) = (v[i] - svmd.mean.at<float>(0, i));
        
        if (svmd.sigma.at<float>(0, i) > 0) 
            f /= svmd.sigma.at<float>(0, i);
    }
#ifdef _DEBUG || DEBUG
	write_hoc_vector(test, "c:\\work\\tmp\\v.txt");
#endif

    return svmd.svm->predict(test, true);
}


void use_library_thresholds(layer1_result* res, part_lib* library, int z, double svmdelta)
{
    if (library == nullptr) 
        return;
	int edgename = atom("toPrevLayer");

	for (list<node*>::iterator iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
		node* n = *iter;
		layer1_data* nd = (layer1_data*)n->data;

		if (n->is_attr_set(NODE_DELETED_ATTR))
			continue;
        if (nd->z == z) {
            node* nn = n->get_neighbor(edgename);
            layer1_data* nnd = (layer1_data*)nn->data;

            if (nnd->z > library->max_layer_index() || nnd->m >= library->layer_size(nnd->z)) {
                cout << "Library is not compatible with result!";
                throw;
            }

            bool ok = true;
            node* p = library->parts[nnd->z][nnd->m];
            vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);
            double svmdfval = 0.0;

            ok = p->is_attr_set(NODE_DELETED_ATTR) == false;

            if (ok && pd != nullptr && pd->svmt.svm != nullptr) {
                int nchildren = p->count_neighbors(atom("lyrSrc"));

                svmdfval = check_svmt(pd->svmt, nchildren, nn, true);

                ok = svmdfval > svmdelta;
                svmdfval = svmdfval - svmdelta;
                //cout << 'f' << f;
                //benergy = f + 0.51;
            }
            if (!ok) {
                n->set_attr(NODE_DELETED_ATTR);
                nn->set_attr(NODE_DELETED_ATTR);
            }

		}
	}		 
    
}


void update_hm_map(vector<int>& m, const vector<const layer1_result::box_data_t*>& pboxes, double prec)
{
    for (int i = 0; i < pboxes.size(); ++i) {
        int j = int_round(prec*pboxes[i]->val);

        if (j >= 0 && j < m.size()) ++m[j];
    }
}

void save_purge_info(const string& fname, layer1_result* res, const list<pair<irectangle2, string> >& gtruths, int layer, double ipercent, 
		   double resp)
{
	typedef list<pair<irectangle2, string> > gtr_t;

	set<node*> toremove;
	vector<pair<double, node*> > outside, inside;
    cout << '!';

	for (list<node*>::iterator niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
		node* n = *niter;
		layer1_data* nd = (layer1_data*)n->data;
		set<node*> rec;
		set<ipoint2> prec;

		if (nd->z != layer)
			continue;
        if (gtruths.empty()) 
            outside.push_back(pair<double, node*>(nd->r(resp), n));
        else {
		    res->recurse_from_node(n, atom("toPrevLayer"), rec);
		    node_set_to_point_set(prec, rec.begin(), rec.end());

		    irectangle2 box = irectangle2::bounding_rectangle(prec.begin(), prec.end());

		    for (gtr_t::const_iterator giter = gtruths.begin(); giter != gtruths.end(); ++giter) {
			    double x = ((double)(giter->first.intersection(box).area()))/giter->first.union_area(box);

			    if (x < ipercent) outside.push_back(pair<double, node*>(nd->r(resp), n));
                else inside.push_back(pair<double, node*>(nd->r(resp), n));
		    }
        }
	}

	sort(outside.begin(), outside.end());
    sort(inside.begin(), inside.end());

    ofstream os(fname.c_str());

    for (auto oiter = outside.begin(); oiter != outside.end(); ++oiter) {
        node* n = oiter->second;
        node* m = n->get_neighbor(atom("toPrevLayer"));

        if (m != nullptr) {
            layer1_data* nd = (layer1_data*)m->data;

            os << 0 << ',' << nd->m << ',' << nd->r(G_RESPONSE) << ',' << nd->r(RR_RESPONSE) << ',' << nd->r(S_RESPONSE) << endl;
        }
    }
    for (auto iiter = inside.begin(); iiter != inside.end(); ++iiter) {
        node* n = iiter->second;
        node* m = n->get_neighbor(atom("toPrevLayer"));

        if (m != nullptr) {
            layer1_data* nd = (layer1_data*)m->data;

            os << 1 << ',' << nd->m << ',' << nd->r(G_RESPONSE) << ',' << nd->r(RR_RESPONSE) << ',' << nd->r(S_RESPONSE) << endl;
        }
    }
}

void purge(layer1_result* res, const list<pair<irectangle2, string> >& gtruths, int layer, double ipercent, 
		   double resp, double percent)
{
	typedef list<pair<irectangle2, string> > gtr_t;

	if (percent <= 0.0) return;
	set<node*> toremove;
	vector<pair<double, node*> > outside;

	for (list<node*>::iterator niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
		node* n = *niter;
		layer1_data* nd = (layer1_data*)n->data;
		set<node*> rec;
		set<ipoint2> prec;

		if (nd->z != layer)
			continue;
        if (gtruths.empty()) 
            outside.push_back(pair<double, node*>(nd->r(resp), n));

        else {
		    res->recurse_from_node(n, atom("toPrevLayer"), rec);
		    node_set_to_point_set(prec, rec.begin(), rec.end());

		    irectangle2 box = irectangle2::bounding_rectangle(prec.begin(), prec.end());
			double maxx = 0.0;

		    for (gtr_t::const_iterator giter = gtruths.begin(); giter != gtruths.end(); ++giter) {
			    double x = ((double)(giter->first.intersection(box).area()))/giter->first.union_area(box);

				if (x > maxx) maxx = x;
		    }
			if (maxx < ipercent) outside.push_back(pair<double, node*>(nd->r(resp), n));
        }
	}

	sort(outside.begin(), outside.end());

    if (!outside.empty()) {
	
        int maxn = int_round(outside.size()*percent);

	    for (int i = 0; i < maxn; ++i) {
            int j = int_round((((double)outside.size())/maxn)*i);

            if (j < 0) j = 0; else if (j >= outside.size()) j = outside.size() - 1;

		    node* n = outside[j].second;

		    n->set_attr(NODE_DELETED_ATTR);
		    node* m = n->get_neighbor(atom("toPrevLayer"));
		    if (m != nullptr) m->set_attr(NODE_DELETED_ATTR);
	    }
    }
}

void save_boxes(layer1_result* res, int layer_index, double inhibitboxes, const string& out)
{
    list<layer1_result::box_data_t> boxes;
    ofstream os(out.c_str());

    res->get_boxes(boxes, layer_index, set<int>(), inhibitboxes);

    os << irectangle2(0, 0, res->x_size(0), res->y_size(0)) << endl;
    for (list<layer1_result::box_data_t>::iterator iter = boxes.begin(); iter != boxes.end(); ++iter) {
        os << iter->box << endl;
    }
                
}

//save_hm_table_parts(hmdir + "hm_table_parts.txt", hmmap, total_objects, total_images);
void save_hm_table_parts(const string& fname, map<int, pair<int, int> >&hmmap, int objnum, int imgnum)
{
    ofstream os(fname.c_str());

    os << objnum << ',' << imgnum << '\n';
    for (auto i = hmmap.begin(); i != hmmap.end(); ++i) {
		os << i->first << ',' << i->second.first << ',' << i->second.second << '\n';
    }
}


void save_hm_table(const string& fname, const vector<int>& hitmap, const vector<int>& missmap, double prec, int objnum, int imgnum)
{
    ofstream os(fname.c_str());

    os << objnum << ',' << imgnum << '\n';
    for (int i = 0; i < hitmap.size() && i < missmap.size(); ++i) {
        os << i/prec << ',';
        os << hitmap[i] << ',';
        os << missmap[i] << '\n';
    }
}

void remove_forbidden_part_boxes(list<layer1_result::box_data_t>& boxes, const set<int>& parts)
{
    typedef list<layer1_result::box_data_t> boxes_t;
    typedef list<pair<irectangle2, string> > gtruths_t;

	if (parts.empty())
		return;

    boxes_t::iterator biter = boxes.begin();

    while (biter != boxes.end()) {
		if (parts.find(biter->nm) != parts.end()) biter = boxes.erase(biter);
        else ++biter;
    }
}


void remove_ignore_region_boxes(list<layer1_result::box_data_t>& boxes, const list<pair<irectangle2, string> >& gtruths)
{
    typedef list<layer1_result::box_data_t> boxes_t;
    typedef list<pair<irectangle2, string> > gtruths_t;

    const string ignore_region_string = "Ignore-Region";

    boxes_t::iterator biter = boxes.begin();

    while (biter != boxes.end()) {
        bool remove = false;

        for (gtruths_t::const_iterator giter = gtruths.begin(); giter != gtruths.end(); ++giter) {
            if (giter->second == ignore_region_string && 
                    giter->first.intersection(biter->box).area() > 0) {
                remove = true;
                break;
            }
        }
        if (remove) biter = boxes.erase(biter);
        else ++biter;
    }
}


pair<int, double> path_dhistogram_classification(irectangle2& box, part_lib* library, layer1_result* res, node* n, const map<vector<int>, int>& typemap, 
	const dclassifier_data& dcd, int depth)
{
	vector<double> h;

    //void make_path_dhistogram(vector<double>& hist, irectangle2& bbox, const map<vector<int>, int>& typemap, 
    //  layer1_result* res, node* n)

    make_path_dhistogram(h, box, typemap, library, res, n, depth);
    //make_path_dhistogram(hd.h, box, typemap[nnd->m], res.get(), nn);

	// ---------------- DEBUG --------------
	//static int index = 0;
	//layer1_data* nd = (layer1_data*)n->data;
	//string s = "c:\\work\\tmp\\h-";
	//s = s + index;
	//s = s + ".hist";
	//ofstream os(s.c_str(), ios_base::binary | ios_base::out);
 //   
	//bin_write_int(os, 1);
	//bin_write_int(os, 0);
	//bin_write_double(os, nd->vval());
	//bin_write_double(os, nd->r(S_RESPONSE));
	//write_sparse_histogram(os, h);
	//
	//os.close();

 //   ofstream os1((s + ".model").c_str(), ios_base::binary | ios_base::out);
 //   bin_write_int(os1, 1);
	//bin_write_int(os1, nd->m);
	//bin_write_double(os1, 0.0);
	//bin_write_double(os1, 0.0);
 //   write_sparse_histogram(os1, dcd.w);

 //   os1.close();

	//---------------- DEBUG --------------
	//h = extract_from_indices(h, dcd.selection);
	h = extract_from_binary(h, dcd.valid);

	for (int i = 0; i < h.size(); ++i) {
		double& hi = h[i];

		hi = (hi - dcd.mu[i])/dcd.vx[i];
	}
	h.push_back(1.0);

	double P = 0.0;
	int label;

	for (int i = 0; i < h.size(); ++i) 
		P += h[i]*dcd.w[i];
	label = (P > 0) ? 1 : 0;
	P = 1.0/(1 + exp(-P));

	// ------- DEBUG --------
	//s = "c:\\work\\tmp\\res-";
	//s = s + index;
	//s = s + ".hist";
	//os.open(s.c_str(), ios_base::binary | ios_base::out);
	//if (label == 0) { bin_write_double(os, 0.0); bin_write_double(os, 1 - P); }
	//else { bin_write_double(os, 1.0); bin_write_double(os, P); }
	//os.close();
	//++index;
	// ------- DEBUG --------

    //if (index > 3) throw;

	if (label == 0) return pair<int, double>(0, 1.0 - P);
	else return pair<int, double>(1, P);
}



// local functions
///////////////////////////////////////////////////////////////////////////////




void test_export(config_dictionary& cfg, const string& fname, string out)
{
    layer1_result* res;
    img* im = nullptr;
    int layer_index = cfg.get_value_int("layer", -1);
    bool alltypes = cfg.get_value_bool("all_types", false);
    vector<int> parts;

    cfg.get_value(parts, "parts");
    read_layer1_result(res, fname);

    if (res != nullptr) {
        //layer1_result_export es;

        //res->export_structure(es, layer_index, parts, alltypes);
        //cout << "# of points: " << es.points.size() << endl;
        //cout << "# of lines: " << es.lines.size() << endl;
        //cout << "# of color_point: " << es.color_point.size() << endl;
        //cout << "# of point_color: " << es.color_point.size() << endl;
        //cout << "# of color_line: " << es.color_line.size() << endl;
    }
    delete res;
}

void save_video(config_dictionary& cfg, const string& fname, string out)
{
	layer1_result* res;
    img* im = nullptr;
	img* im_back = nullptr;
    int layer_index = cfg.get_value_int("layer", -1);
    int end_layer_index = cfg.get_value_int("end_layer", 0);
    bool drawbox = cfg.get_value_bool("draw_box", false);
    string mode = cfg.get_value_string("mode", "r"); 
    double factorpow = cfg.get_value_double("factor_power", 1.0);
    bool alltypes = cfg.get_value_bool("all_types", false);
    bool inhibitboxes = cfg.get_value_bool("inhibit_boxes", false);
	bool bigpoints = cfg.get_value_bool("draw_big_points", false);
    vector<int> parts, icolors;
    vector<color> colors;
    color defcolor, uncovered;
    bool paintuncov;

    cfg.get_value(parts, "parts");
    cfg.get_value(icolors, "part_colors");
    color::to_color_vector(colors, icolors, 0);

    icolors.clear();
    cfg.get_value(icolors, "default_color");
    color::to_color(defcolor, icolors, 255);

    icolors.clear();
    cfg.get_value(icolors, "uncovered_color");
    paintuncov = (icolors.size() != 0);
    color::to_color(uncovered, icolors, 20);

    read_layer1_result(res, fname);
	
	if (!cfg.get_value_bool("animation", false) || parts.empty()) {
		im_back = res->get_image_reconstructed(0, 0, vector<int>(), true, false, factorpow, true);
		//(int z, int zt, const vector<int>& parts,
  //  bool colored /* = false */, bool drawbox /* = false */, double factorpow /* = 1.0 */, bool all /* = false */)
		if (cfg.get_value_bool("draw_groundtruth", false)) {
			
			list<pair<irectangle2,string> > rect_list;
			read_groundtruth(rect_list, fname);

			for (list<pair<irectangle2,string> >::iterator iter = rect_list.begin();
					iter != rect_list.end(); iter++) {
				pair<irectangle2,string>& p = *iter;
				im_back->draw_box(p.first, 1.0);
			}

			/*ifstream is(change_extension(fname, ".groundtruth").c_str());
			rectangle2<double> rect;

			while (true) {
				is >> rect;
				if (is.fail()) break;
				irectangle2 irect((int)rect.ll.x, (int)rect.ll.y, (int)rect.ur.x, (int)rect.ur.y);
				im_back->draw_box(irect, 1.0);
			}*/
		}
	} else {
		vector<int> dispparts;
		char str[5];

		for (int i = 0; i < (int)parts.size(); ++i) {
			dispparts.push_back(parts[i]);
			im_back = res->get_image_reconstructed(0, 0, vector<int>(), true, false, factorpow, true);
			sprintf(str, "%d", i);
			if (im_back) {
				string name = change_extension(out, string(str) + ".bmp");
				cout << name << endl;
				//im_back->save_normalized(name.c_str());
				delete im_back;
			}
		}
		im_back = nullptr;
	}
	im = res->get_image_new(layer_index, parts, defcolor, colors, bigpoints);

	for(int x = 0; x < im->width; ++x)
	{
		for(int y = 0; y < im->height; ++y)
		{
			/*if(im->at(x,y) != 0)
			{
				im_back->at(x,y) = im->at(x,y);
			}*/
			if(im->at(x,y) == 0 && im_back->at(x,y) != 0)
			{
				im->at(x,y) = im_back->at(x,y);
			}
		}
	}
	if(im) im->save_normalized(out.c_str());
	//im_back->save_normalized(out.c_str());
	if (im_back) delete im_back;
	if(im) delete im;
}

void save_normal(config_dictionary& cfg, const string& fname, string out)
{
    layer1_result* res;
    img* im = nullptr;
    int layer_index = cfg.get_value_int("layer", -1);
    int end_layer_index = cfg.get_value_int("end_layer", 0);
    bool drawbox = cfg.get_value_bool("draw_box", false);
    string mode = cfg.get_value_string("mode", "r"); 
    double factorpow = cfg.get_value_double("factor_power", 1.0);
    bool alltypes = cfg.get_value_bool("all_types", false);
    bool inhibitboxes = cfg.get_value_bool("inhibit_boxes", false);
    vector<int> parts, icolors;
    vector<color> colors;
    color defcolor, uncovered;
    bool paintuncov;

    cfg.get_value(parts, "parts");
    cfg.get_value(icolors, "part_colors");
    color::to_color_vector(colors, icolors, 0);

    icolors.clear();
    cfg.get_value(icolors, "default_color");
    color::to_color(defcolor, icolors, 255);

    icolors.clear();
    cfg.get_value(icolors, "uncovered_color");
    paintuncov = (icolors.size() != 0);
    color::to_color(uncovered, icolors, 20);

    read_layer1_result(res, fname);
    
    if (out == "") {
        string::size_type dot_pos;
    
        if ((dot_pos = fname.rfind(".")) == string::npos) out = fname + ".bmp";
        else out = fname.substr(0, dot_pos) + ".bmp";
    }
    
    switch (mode[0]) {
        case 'N' : 
            im = res->get_image(layer_index, end_layer_index, paintuncov, uncovered, defcolor, parts, colors); 
            break;
        case 'I' : 
            im = res->get_image_inhib(layer_index, end_layer_index, paintuncov, uncovered, defcolor, parts, colors); 
            break;
        case 'B' :
            im = res->get_image_boxed(layer_index, parts, colors, defcolor, inhibitboxes);
            
            if (cfg.get_value_bool("save_boxes", false)) 
                save_boxes(res, layer_index, cfg.get_value_double("box_inhibition_percent", 0.1),
                    change_extension(out, ".bb"));
            break;
        case 'X' :
            if (layer_index < (int)res->shape_nodes_inhib.size()) {
                vector<node*>& nodes = res->shape_nodes_inhib[layer_index];
                char str[5];
                int i = 0;

                for (vector<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
                    //im = res->get_part_reconstructed(layer_index, end_layer_index, i);
                    im = res->get_part_image(layer_index, end_layer_index, *iter);
                    if (im) {
                        sprintf(str, "%d", i++);
                        string name = change_extension(out, string(str) + ".bmp");
                        im->save_normalized(name.c_str());
                        delete im;
                    }
                }
                cout << endl;
                im = nullptr;
            }
            break;
        case 'Y' :
            {
                part_lib* library;
                string libname = cfg.get_value_string("part_lib", "");

                read_library(libname, library);
                im = res->get_image_reconstructed(library, layer_index, parts);
                if (library) delete library;
            }
            break;
        default : 
            if (!cfg.get_value_bool("animation", false) || parts.empty()) {
                im = res->get_image_reconstructed(layer_index, end_layer_index, parts, false, drawbox, factorpow, alltypes);
                if (cfg.get_value_bool("draw_groundtruth", false)) {
                    
					list<pair<irectangle2,string> > rect_list;
					read_groundtruth(rect_list, fname);

					for (list<pair<irectangle2,string> >::iterator iter = rect_list.begin();
							iter != rect_list.end(); iter++) {
						pair<irectangle2,string>& p = *iter;
						im->draw_box(p.first, 1.0);
					}

					/*ifstream is(change_extension(fname, ".groundtruth").c_str());
                    rectangle2<double> rect;

                    while (true) {
                        is >> rect;
                        if (is.fail()) break;
                        irectangle2 irect((int)rect.ll.x, (int)rect.ll.y, (int)rect.ur.x, (int)rect.ur.y);
                        im->draw_box(irect, 1.0);
                    }*/

                }
            } else {
                vector<int> dispparts;
                char str[5];

                for (int i = 0; i < (int)parts.size(); ++i) {
                    dispparts.push_back(parts[i]);
                    im = res->get_image_reconstructed(layer_index, end_layer_index, dispparts, 
                        false, drawbox, factorpow);
                    sprintf(str, "%d", i);
                    if (im) {
                        string name = change_extension(out, string(str) + ".bmp");
                        cout << name << endl;
                        im->save_normalized(name.c_str());
                        delete im;
                    }
                }
                im = nullptr;
            }
    }
    if (im) {
        im->save_normalized(out.c_str());
        delete im;
    }
}


void save_visview(config_dictionary& cfg, const string& fname, const string& out)
{
    layer1_result* res = nullptr;
    part_lib* library = nullptr;
    string outdir = cfg.get_value_string("out_dir", "");
    string srcdir = cfg.get_value_string("src_dir", "");
    string libname = cfg.get_value_string("library", "");
    int save_mode = cfg.get_value_int("save_mode", 0);
    int layer_index = cfg.get_value_int("layer", -1);
    bool save_sc = cfg.get_value_bool("save_sc", true);
    bool save_var = cfg.get_value_bool("save_variations", false);
    double thresh = cfg.get_value_double("thresh", 0.3);
    int resp = response_from_string(cfg.get_value_string("response", "RR_RESPONSE"));
    double svmdelta = cfg.get_value_double("svm_delta", 0.0);
    bool libthresh = cfg.get_value_bool("use_library_thresholds", false);

    if (cfg.get_value_bool("true_coordinates", false))
        save_mode |= VVE_TRUE_COORDINATES;
    if (cfg.get_value_bool("bounding_boxes", false))
        save_mode |= VVE_BOUNDING_BOXES;
    if (cfg.get_value_bool("simple_filenames", false))
        save_mode |= VVE_SIMPLE_FILENAMES;
	if (cfg.get_value_bool("matlab_scale_indexes", false))
        save_mode |= VVE_MATLAB_SCALES;

    end_dir(outdir);
    end_dir(srcdir);

    if (libname.size() > 0) {
        read_library(libname, library);

        if (library) {
            library->save_visview_filters(outdir);
            library->save_visview_centers(outdir);
            library->save_visview_parts(outdir); 
            library->save_visview_contractions(outdir);
            library->save_sc(outdir, false, save_sc, save_var);
        }
    }

    istringstream iss(out);
    int outindex;
    set<double> scales;
    set<string> bbgt;

    iss >> outindex;
    if (iss.fail()) outindex = 1;

    list<string> files;

    list_directory(files, srcdir + fname);

	string filename;
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        cout << "Processing: " << *fiter << endl;
        
        read_layer1_result(res, srcdir + *fiter);

        if (res) {
            ostringstream oss;
            string basefname = change_extension(*fiter, "", "_");
            double scale = round((double)(res->x_size() - 2*res->border)/res->original_width, 2);
            list<pair<irectangle2, string> > rectangles;

            read_groundtruth(rectangles, srcdir + *fiter);
            oss << outindex++;
            if (libthresh) 
                use_library_thresholds(res, library, library->max_layer_index(), svmdelta);
			filename = change_extension(*fiter, "");
            res->save_visview(outdir, oss.str(), filename, library, layer_index, save_mode);
            scales.insert(scale);
            if (bbgt.find(basefname) == bbgt.end()) {
                for (auto riter = rectangles.begin(); riter != rectangles.end(); ++riter) {
                    riter->first -= ipoint2(res->border, res->border);
                    riter->first.ll /= scale;
                    riter->first.ur /= scale;
                }
                save_groundtruth(rectangles, outdir, basefname + "_bbgt.txt", "", "", false);
                bbgt.insert(basefname);
            }
            delete res;
        }
    }

    // save_scales
    vector<double> sscales(scales.begin(), scales.end());
    ofstream sfs((outdir + "scales.txt").c_str());

    sort(sscales.begin(), sscales.end(), greater<double>());
    for (int i = 0; i < (int)sscales.size(); ++i) {
        if (i == 0 || sscales[i - 1]/sscales[i] > 1.1) 
            sfs << sscales[i] << endl;
    }
    sfs.close();

    // save classnames
    if (library) {
        set<string> names;

        for (auto piter = library->nodes.begin(); piter != library->nodes.end(); ++piter) {
            node* p = *piter;

            if (p->is_attr_set(CATEGORY_NODE_ATTR)) {
                cpart_data* pd = (cpart_data*)p->data;
                names.insert(pd->name);
            }
        }
        
        ofstream cfs((outdir + "classnames.txt").c_str());

        for (auto citer = names.begin(); citer != names.end(); ++citer) {
            cfs << *citer << endl;
        }
        cfs.close();
    }

    if (library) delete library;
}

void save_graph(config_dictionary& cfg, const string& fname, string out)
{
    layer1_result* res;
    
    read_layer1_result(res, fname);

    if (res == nullptr) return;

    cout << "Processing " << fname << endl;

    int layer = cfg.get_value_int("layer", 0);
    int iradius = cfg.get_value_int("inner_radius", 2);
    int oradius = cfg.get_value_int("outer_radius", 5);
    vector<int> angles, namemap, values;
    string outcl;

    cfg.get_value(angles, "angles");
    if (angles.size() < 4 || angles.size() > 12) {
        int aarr[] = {15, 45, 75, 105, 135, 165, 195, 225, 255, 285, 315, 345};
        array_to_vector(angles, aarr, 12);
    }
    int anglescount = (int)angles.size();
    int* varr = new int[anglescount];
    char edgename[20];

    for (int i = 0; i < anglescount; ++i) {
        sprintf(edgename, "toDirection%d", i);
        varr[i] = atom(edgename).get_index();
    }
    array_to_vector(values, varr, anglescount);
    make_angle_vector(namemap, angles, values);
    delete[] varr;

    if (out == "") { 
        string::size_type dot_pos;
    
        if ((dot_pos = fname.rfind(".")) == string::npos) { out = fname + ".net"; outcl = fname + ".clu"; }
        else { out = fname.substr(0, dot_pos) + ".net"; outcl = fname.substr(0, dot_pos) + ".clu"; }
    }

    ofstream os(out.c_str());

    if (!os.fail()) {
        res->connect_neighbors_circular_1(res->shape_nodes_inhib[layer], iradius, oradius, NODE_REDUCED_ATTR, 
            namemap, true);
        if (cfg.get_value_string("node_type", "i") == "i") 
            res->write_vgr(os, res->shape_nodes_inhib[layer], -1);
        else
            res->write_vgr(os, res->shape_nodes[layer], -1); 
    }

    os.close();
    os.open(outcl.c_str());

    if (!os.fail()) {
        res->write_clu(os, res->shape_nodes_inhib[layer]);
    }

    delete res;
}

void save_graphs(config_dictionary& cfg, const string& fname, string out)
{
    string srcdir = cfg.get_value_string("src_dir", "");
    list<string> files;

    end_dir(srcdir);
    list_directory(files, srcdir + fname);
    for (list<string>::iterator iter = files.begin(); iter != files.end(); ++iter) {
        save_graph(cfg, srcdir + *iter, out);
    }
}

void drop_library_images(config_dictionary& cfg, const string& fname, string out)
{
    part_lib* library;

    read_library(fname, library);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }
    cout << "Dropping images." << endl;
    library->drop_images();
    cout << "Saving " << fname << " to " << out << "." << endl;
    library->save(out);
    delete library;
}

void save_library_mma(config_dictionary& cfg, const string& fname, string out)
{
    part_lib* library;

    cout << "Loading " << fname;
    read_library(fname, library);
    cout << " done." << endl;
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }

    ofstream os(out.c_str());

    if (os.is_open()) {
        cout << "Saving to: " << out;
        library->save_mma(os);
        cout << " done." << endl;
    }
    os.close();
    delete library;
}

void save_library_mma_back(config_dictionary& cfg, const string& fname, string out)
{
    part_lib* library;

    cout << "Loading " << fname;
    read_library(fname, library);
    cout << " done." << endl;
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }

    ofstream os(out.c_str());

    if (os.is_open()) {
        cout << "\"Back\" saving to: " << out;
        library->save_mma_back(os);
        cout << " done." << endl;
    }
    os.close();
    delete library;
}


void select_files(config_dictionary& cfg, const string& fname, string out)
{
    string dir = cfg.get_value_string("src_dir", "");
    int layer;
    vector<int> parts;
    list<string> files;
    layer1_result* res;
    ofstream os;

    cfg.get_value(layer, "layer", true);
    if (layer < 0) {
        cout << "parameter 'layer' must be >= 0.";
        return;
    }

    end_dir(dir);
    list_directory(files, dir + fname);
    sort(parts.begin(), parts.end());
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        //cout << "Processing " << *iter << endl;

        read_layer1_result(res, dir + *fiter);
        if (res == nullptr) continue;

        if (res->max_layer_index() >= layer) {
            vector<node*>& s_nodes = res->shape_nodes[layer];

            for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
                node* n = *iter;
                
                if (parts.empty()) {
                    cout << *fiter << endl;
                } else {
                    while (n != nullptr) {
                        layer1_data* nd = (layer1_data*)n->data;
                        if (binary_search(parts.begin(), parts.end(), nd->m)) {
                            cout << *iter;
                            break;
                        }
                        n = nd->next;
                    }
                }
            }
        }

        delete res;
    }
}

void display_c2f(config_dictionary& cfg, const string& file, string outname)
{
    part_lib* library;

    read_library(file, library);
    if (library == nullptr) {
        cout << "Can not open " << file << "." << endl;
        return;
    }
    int layer = cfg.get_value_int("layer", 5);
    ofstream os(outname.c_str());
    vector<node*>& parts = library->parts[layer];
    int lyrsrc = atom("lyrSrc");
    int lyrback = atom("lyrCenterBack");
    bool firstpart = true;

    os << '{';
    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* p = *iter;
        rpart_data* pd = (rpart_data*)p->data;
        bool firstnb = true;

        if (firstpart) firstpart = false; else os << ',';
        os << '{' << '{' << pd->cx << ',' << pd->cy << '}' << ',' << '{';
        foreach_neighbor (p, lyrsrc, niter) {
            part_data_2r* ed = (part_data_2r*)neighbor_edge_data(niter);
            irectangle2 rect = ed->rect + ipoint2(ed->x, ed->y);
            
            if (firstnb) firstnb = false; else os << ',';
            rect.mma_write(os);
        }
        os << '}';

        foreach_neighbor (p, lyrback, biter) {
            node* pp = neighbor_node(biter);
            part_data* ppd = (part_data*)pp->data;
            bool firstnb2 = true;

            os << ',' << '{' << ppd->cmx << ',' << ppd->cmy << '}' << ',' << '{';
            foreach_neighbor (pp, lyrsrc, niter) {
                part_data_2* ed = (part_data_2*)neighbor_edge_data(niter);
                irectangle2 rect(-(int)ed->distr.width/2, -(int)ed->distr.height/2, (int)ed->distr.width/2,
                    (int)ed->distr.height/2);
                
                rect += ipoint2(ed->x, ed->y);
                if (firstnb2) firstnb2 = false; else os << ',';
                rect.mma_write(os);
            }
            os << '}';
        }
        os << '}';
        os << endl;
    }
    os << '}' << endl;

    os.close();
    delete library;
}

void display_c2f_graph_layer(vector<node*>& parts, ostream& os)
{
    int lyrback = atom("lyrCenterBack");

    os << '{';
    for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
        node* p = *iter;
        lib_data* pd = (lib_data*)p->data;
        bool first = true;

        if (iter != parts.begin()) os << ',';
        os << '{' << pd->type + 1 << ',' << '{';
        foreach_neighbor(p, lyrback, niter) {
            lib_data* pnd = (lib_data*)neighbor_node_data(niter);

            if (first) first = false; else os << ',';
            os << pnd->type + 1;
        }
        os << '}' << '}';
    }
    os << '}';
}

void display_c2f_graph(config_dictionary& cfg, const string& file, string outname)
{
    part_lib* library;

    read_library(file, library);
    if (library == nullptr) {
        cout << "Can not open " << file << "." << endl;
        return;
    }
    int layer = cfg.get_value_int("layer", 5);
    ofstream os(outname.c_str());

    os << '{';
    for (int i = layer; i <= library->max_layer_index(); ++i) {
        if (i != layer) os << ',' << endl;
        display_c2f_graph_layer(library->parts[i], os);
    }
    os << '}';

    os.close();
    delete library;
}


void save_category_parts(config_dictionary& cfg, const string& fname, string out)
{
    part_lib* library;

    read_library(fname, library);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }

    ofstream os(out.c_str());
    int first = true;
    int toback = atom("lyrCenterBack");

    if (os.is_open()) {
        int layer = cfg.get_value_int("category_layer", 6);
        int part = cfg.get_value_int("part", 0);
        
        
        vector<node*>& parts = library->parts[layer - 1];
        for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
            node* p = *iter;
            lib_data* pd = (lib_data*)p->data;
            
            foreach_neighbor (p, toback, niter) {
                lib_data* pnd = (lib_data*)neighbor_node_data(niter);

                if (pnd->type == part) {
                    if (first) first = false; else os << ',';
                    os << pd->type;
                }
            }
        }
    }
    os.close();
    delete library;
}


void save_library_layer(config_dictionary& cfg, const string& fname, string out)
{
    part_lib* library;

    read_library(fname, library);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }

    ofstream os(out.c_str());

    if (os.is_open()) library->save_mma2(os, cfg.get_value_int("layer", 0));
    os.close();
    delete library;
}

void save_library_graph(config_dictionary& cfg, const string& fname, string out)
{
    part_lib* library;

    read_library(fname, library);
    if (library == nullptr) return;

    ofstream os(out.c_str());

    if (os.fail()) return;

    vector<string> parts;
    bool edgel = cfg.get_value_bool("edge_labels", true);
    int maxlayer = cfg.get_value_int("max_layer", -1);
    int maxtype = cfg.get_value_int("max_type", 100);

    cfg.get_value(parts, "edges");

    set<int> pset;

    for (size_t i = 0; i < parts.size(); ++i)
        pset.insert(atom(parts[i]));
    library->save_as_graph_mma(os, maxlayer, maxtype, pset, edgel);
    os.close();
    delete library;
}

void save_library_yml(config_dictionary& cfg, const string& fname, string out)
{
    unique_ptr<part_lib> library(read_library(fname));

    if (library == nullptr) {
        cout << "Can not open library '" << fname << "'" << endl;
    }
    cout << "Saving library '" << fname << "' to '" << out << "'" << endl;

    cv::FileStorage fs(out, cv::FileStorage::WRITE);

    library->write_yml(fs);
}


void save_library(config_dictionary& cfg, const string& fname, string out)
{
    int layer = cfg.get_value_int("layer", 1);
    part_lib* library;

    read_library(fname, library);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }
    if (layer >= (int)library->parts.size()) { 
        cout << "Layer " << layer << " not found in " << fname << "." << endl;
        return;
    }

    if (cfg.get_value_bool("add_similarity_edges", false)) {
        part_lib::similarity_threshold = cfg.get_value_double("threshold", 0.1);
        double contraction = cfg.get_value_double("contraction", 0.5);

        cout << "Adding similarity edges with threshold " << part_lib::similarity_threshold << " and" << endl;
        cout << "contraction factor " << contraction << endl;
        library->add_similarity_edges(layer - 1, contraction);
        cout << "Saving library..." << endl;
        library->save(out);
    } else if (cfg.is_defined("part")) {
        int part = cfg.get_value_int("part", -1);
        cout << "Saving part #" << part << " to " << out << "..." << endl;
        library->save_part(out.c_str(), layer, part, cfg.get_value_double("threshold", 0.0));
    } else if (cfg.is_defined("indices")) {
        vector<int> l;

        cfg.get_value(l, "indices", true);
        cout << "Saving " << fname << " to " << out << endl;
        library->save_all(out.c_str(), layer, l, cfg.get_value_bool("show_labels", true), cfg.get_value_bool("mark_center", false));
    } else {
        cout << "Saving " << fname << " to " << out << "." << endl;
        vector<double> contr;

        cfg.get_value(contr, "layer_contractions");
        part_data::set_contractions(contr);
        library->save_all(out.c_str(), layer, 0, -1, cfg.get_value_bool("show_labels", true), false, 
            cfg.get_value_bool("mark_center", false), cfg.get_value_bool("one_row", false));
    }
    delete library;
}

void save_library_sc(config_dictionary& cfg, const string& fname, string out)
{
    int layer = cfg.get_value_int("layer", 1);
    part_lib* library;

    read_library(fname, library);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }
    if (layer < 1 || layer - 1 > library->max_layer_index()) { 
        cout << "Layer " << layer << " not found in " << fname << "." << endl;
    } else {
        cout << "Saving " << fname << " to " << out << "." << endl;
        library->save_all_sc(out.c_str(), layer, 0, -1, cfg.get_value_bool("show_labels", true));
    }
    delete library;
}

void save_library_mean(config_dictionary& cfg, const string& fname, string out)
{
    int layer;
    part_lib* library;
    string pattern;
    bool savev = cfg.get_value_bool("vectors", false);

    cfg.get_value(layer, "layer", true);
    pattern = cfg.get_value_string("pattern", "part-%d-%03d.txt");

    read_library(fname, library);

    end_dir(out);
    if (library == nullptr) {
        cout << "Can not open library." << endl;
        return;
    }
    if (layer < 0 || layer > library->max_layer_index()) {
        cout << "Layer " << layer << " does not exist." << endl;
    } else {
        for (int pi = 0; pi < library->parts[layer].size(); ++pi) {
            node* p = library->parts[layer][pi];
            part_data* pd = (part_data*)p->data;
            vs_part_data* vspd = dynamic_cast<vs_part_data*>(p->data);
            char s[255];

            sprintf(s, pattern.c_str(), layer, pi);

            ofstream os((out + s).c_str());
            vector<dpoint2> dpts;

            if (savev) {
                for (int i = 0; i < vspd->pcad.mean.cols; i += 2) {
                    os << vspd->pcad.mean.at<double>(0, i) << ',' << vspd->pcad.mean.at<double>(0, i + 1);
                    for (int j = 0; j < vspd->pcad.eigenvectors.rows; ++j) {
                        os << ',' << vspd->pcad.eigenvectors.at<double>(j, i) << ',' 
                            << vspd->pcad.eigenvectors.at<double>(j, i + 1);
                    }
                    os << '\n';
                }
            } else {
                if (vspd != nullptr) {
                    dpts = partition(vspd->pcad.mean);
                } else {
                    vector<ipoint2> ipts = get_library_geo(p);

                    dpts = cast_vector<dpoint2, ipoint2>(ipts);
                }
                translate_and_scale(dpts);
                for (auto piter = dpts.begin(); piter != dpts.end(); ++piter) {
                    dpoint2 p = (*piter)*50.0;

                    os << (int)(p.x) << ',' << (int)(p.y) << '\n';
                }
            }
        }
    }
    
    delete library;
}

void save_library_sc_mma(config_dictionary& cfg, const string& fname, string out)
{
    int layer = cfg.get_value_int("layer", 1);
    part_lib* library;

    read_library(fname, library);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }
    if (layer < 1 || layer - 1 > library->max_layer_index()) { 
        cout << "Layer " << layer << " not found in " << fname << "." << endl;
    } else {
        cout << "Saving " << fname << " to " << out << "." << endl;
        library->save_all_sc_mma(out.c_str(), layer, 0, -1);
    }
    delete library;
}

void remove_library_layers(config_dictionary& cfg, const string& fname, string out)
{
    int layer;
    part_lib* library;

    cfg.get_value(layer, "layer", true);
    read_library(fname, library);
    if (library != nullptr) {
        library->delete_layers_geq(layer);
        library->save(out);
        delete library;
    }
}

void save_similarity_matrix(config_dictionary& cfg, const string& fname, string out)
{
    part_lib::similarity_threshold = cfg.get_value_double("threshold", 0.1);
    double contraction = cfg.get_value_double("contraction", 0.5);
    int layer = cfg.get_value_int("layer", 0);
    part_lib* library;

    read_library(fname, library);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }
    
    cout << "Creating similarity matrix with threshold " << part_lib::similarity_threshold << " and" << endl;
    cout << "contraction factor " << contraction << endl;

    vector<node*>& lparts = library->parts[layer];
    vector<node*>::iterator iter1, iter2;
    matrix<double> result(lparts.size(), lparts.size(), 0.0);
    int i, j;

    for (iter1 = lparts.begin(), i = 0; iter1 != lparts.end(); ++iter1, ++i) {
        for (iter2 = lparts.begin(), j = 0; iter2 != lparts.end(); ++iter2, ++j) {
            double& sim = result(i, j);

            if (i == j) sim = 1.0;
            else sim = library->similarity(*iter1, *iter2, contraction, 1, part_lib::defaultMfunction);
        }
    }

    cout << "Saving matrix ..." << endl;
    result.save_mathematica(out.c_str());

    delete library;
}

void print_parts(config_dictionary& cfg, const string& fname, string out)
{
    int layer = cfg.get_value_int("layer", 1);
    part_lib* library;

    read_library(fname, library);
    
    if (library == nullptr) { 
        cout << "Can not open " << fname << "." << endl; 
        return; 
    }

	library->display_layer_info(layer);

    delete library;
}

void save_mathematica2(config_dictionary& cfg, const string& fname, string out)
{
    string dir = cfg.get_value_string("src_dir", "");
    string suffix = cfg.get_value_string("suffix", "");
    string prefix = cfg.get_value_string("prefix", "");
    string outdir = cfg.get_value_string("out_dir", "");
    bool alltypes = cfg.get_value_bool("all_types", false);
    int layer = cfg.get_value_int("layer", -1);
    int edgename = atom(cfg.get_value_string("edge_name", "dummy")).get_index();
    list<string> files;
    layer1_result* res;
    ofstream os;

    end_dir(dir);
    end_dir(outdir);
    list_directory(files, dir + fname);
    for (list<string>::iterator iter = files.begin(); iter != files.end(); ++iter) {
        cout << "Processing " << *iter << endl;

        read_layer1_result(res, dir + *iter);
        if (res == nullptr) continue;

        os.open(change_extension(outdir + prefix + *iter, suffix + ".m").c_str());    
        if (layer < 0) res->write_mma2(os);
        else res->write_mma2(os, res->shape_nodes[layer], edgename, alltypes);
        os.close();
        delete res;
    }
}

void save_part_statistics(config_dictionary& cfg, const string& fname, string out)
{
    string dir = cfg.get_value_string("src_dir", "");
    string odir = cfg.get_value_string("out_dir", "");
    int layer = cfg.get_value_int("layer", 0);
    list<string> files;
    vector<int> ostat;

    end_dir(dir);
    end_dir(odir);
    list_directory(files, dir + fname);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        layer1_result* res;

        cout << "Processing " << *file; 
        read_layer1_result(res, dir + *file);

        if (res == nullptr) cout << " failed!" << endl;
        else {
            vector<int> fstat;

            res->type_statistics(fstat, layer, true);

            ofstream os((odir + change_extension(*file, ".txt")).c_str());

            if (!os.fail()) {
                int size = (int)fstat.size();

                if (ostat.size() < fstat.size()) ostat.resize(fstat.size(), 0);
                for (int i = 0; i < size; ++i) {
                    os << i << ',' << fstat[i] << endl;
                    ostat[i] += fstat[i];
                }
                os.close();
            }
        
            delete res;
            cout << " done." << endl;
        }
    }

    ofstream os(out.c_str());

    if (!os.fail()) {
        int size = (int)ostat.size();

        for (int i = 0; i < size; ++i)
            os << i << ',' << ostat[i] << endl;

        os.close();
    }
}

// Return value: the number of such boxes B1 in boxes1 s.t. there exists some box B2 in boxes2 
//   with the property that A = area(B1 ^ B2) >= area(B1)*percent and A >= area(B2)*percent.
int count_intersections(const list<irectangle2>& boxes1, const list<irectangle2>& boxes2, double percent)
{
    int result = 0;

    for (list<irectangle2>::const_iterator iter1 = boxes1.begin(); iter1 != boxes1.end(); ++iter1) {
        double area1 = iter1->area()*percent;

        for (list<irectangle2>::const_iterator iter2 = boxes2.begin(); iter2 != boxes2.end(); ++iter2) {
            double area = (iter2->intersection(*iter1)).area();

            if (area >= area1 && area >= iter2->area()*percent) ++result;
        }
    }
    return result;
}

/**
should be using read_groundtruth instead of this one
void read_boxes(list<irectangle2>& boxes, const string& file)
{
    ifstream is(file.c_str());
    rectangle2<double> rect;

    boxes.clear();
    while (true) {
        is >> rect;
        if (is.fail()) break;
        boxes.push_back(irectangle2((int)rect.ll.x, (int)rect.ll.y, (int)rect.ur.x, (int)rect.ur.y));
    }
}
*/
void compare_boxes(config_dictionary& cfg, const string& fname, string out)
{
    string dir = cfg.get_value_string("src_dir", "");
    string gtruthext = cfg.get_value_string("groundtruth_extension", ".groundtruth");
    double inhibpct = cfg.get_value_double("box_inhibition_percent", 0.1);
    double intersectionpct = cfg.get_value_double("box_intersection_percent", 0.3);
    int layer;
    list<string> files;
    layer1_result* res;
    ofstream os(out.c_str());

    cfg.get_value(layer, "layer", true);
    end_dir(dir);
    list_directory(files, dir + fname);
    for (list<string>::iterator iter = files.begin(); iter != files.end(); ++iter) {
        cout << "Processing " << *iter << endl;
        read_layer1_result(res, dir + *iter);
        if (res == nullptr || res->max_layer_index() < layer) continue;

        list<irectangle2> gtruths, pure_boxes;
        list<layer1_result::box_data_t> boxes;

		read_groundtruth(gtruths, dir + *iter, "", gtruthext);
        //read_boxes(gtruths, change_extension(dir + *iter, gtruthext));
        res->get_boxes(boxes, layer, set<int>(), inhibpct);

        for (list<layer1_result::box_data_t>::iterator biter = boxes.begin(); biter != boxes.end(); ++biter) 
            pure_boxes.push_back(biter->box);

        os << *iter << ' ' << count_intersections(gtruths, pure_boxes, intersectionpct) << 
            '/' << gtruths.size() << ' ' << pure_boxes.size() << endl;

        delete res;
    }
}

void save_part_covering_statistics(config_dictionary& cfg, const string& fname, string out)
{
    string dir = cfg.get_value_string("src_dir", "");
    string outdir = cfg.get_value_string("out_dir", "");
    list<string> files;
    layer1_result* res;
    vector<int> parts;
    ofstream os;
    int layer;
    
    cfg.get_value(layer, "layer", true);
    cfg.get_value(parts, "parts", true);

    if (layer < 0) return;

    end_dir(dir);
    end_dir(outdir);
    list_directory(files, dir + fname);
    os.open((outdir + out).c_str());    

    if (!os) {
        cout << "Can not open file for writing." << endl;
        return;
    }
    for (list<string>::iterator iter = files.begin(); iter != files.end(); ++iter) {
        vector<int> filestat;
    
        cout << "Processing " << *iter << endl;

        read_layer1_result(res, dir + *iter);
        if (res == nullptr || res->max_layer_index() < layer) continue;

        res->get_covering_statistics(res->shape_nodes[layer], parts, filestat);

        // write "filestat" to os
        os << *iter << ',';
        os << res->shape_nodes[0].size(); 
        for (int i = 0; i < (int)filestat.size(); ++i) {
            os << ',' << filestat[i];
        }
        os << endl;
        delete res;
    }
    os.close();
}

void average_part_normalization(double& d)
{
    d = log(d + 1.0);
}

void save_average_parts(config_dictionary& cfg, const string& fname, const string& out)
{
    layer1_result* res = nullptr;
    string outdir = cfg.get_value_string("out_dir", "");
    string srcdir = cfg.get_value_string("src_dir", "");
    int layer_index = cfg.get_value_int("layer", 0);
    map<int, img> result;
    list<string> files;

    end_dir(outdir);
    end_dir(srcdir);

    list_directory(files, srcdir + fname);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        cout << "Processing: " << *fiter << endl;
        read_layer1_result(res, srcdir + *fiter);
        if (res) {
            res->part_reconstruction(result, layer_index);
            delete res;
        }
    }
    vector<img*> images;

    for (map<int, img>::iterator iter = result.begin(); iter != result.end(); ++iter) {
        img& im = iter->second;
        img number = img::number_to_img(iter->first);

        if (layer_index > 3) 
            number = number.get_resized(-200, -200);
        
        //HOP_REAL im_min, im_max;

        //im.minmax(im_min, im_max);
        //im += im_min;
        //im /= (im_max - im_min);
        //map_through(im.begin(), im.end(), average_part_normalization);

        im.blt(number, 0, 0);
        images.push_back(&im);
    }

    img res_image = img::concat(images);

    res_image.save(out);
}

void save_edge_names(config_dictionary& cfg, const string& fname, const string& out)
{
    layer1_result* res = nullptr;
    string srcdir = cfg.get_value_string("src_dir", "");
    list<string> files;
    ofstream os(out.c_str());
    set<itriple> result;

    end_dir(srcdir);

    if (os.fail())
        return;

    list_directory(files, srcdir + fname);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        cout << "Processing: " << *fiter << endl;
        read_layer1_result(res, srcdir + *fiter);
        if (res) {
            res->get_edge_info(result);
            delete res;
        }
    }
    layer1_result::write_edge_info(os, result);
    os.close();
}

void save_mathematica(config_dictionary& cfg, const string& fname, string out)
{
    layer1_result* res = nullptr;
    string srcdir = cfg.get_value_string("src_dir", "");
    list<string> files;

    end_dir(srcdir);
    end_dir(out);
    list_directory(files, srcdir + fname);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        cout << "Processing: " << *fiter << endl;
        read_layer1_result(res, srcdir + *fiter);
        if (res) {
            ofstream os((out + change_extension(*fiter, ".m")).c_str());
            res->write_mma(os);
            os.close();
            delete res;
        }
    }
}

void sort_detections_by_size(vector<iipair>& result, layer1_result* res, int layer)
{
	result.clear();

	if (layer < 0 || layer >= res->layer_count()) return;

	vector<node*>& s_nodes = res->shape_nodes[layer];
	int toprev = atom("toPrevLayer");

	for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
		node* n = *iter;

		while (n != nullptr) {
			layer1_data* nd = (layer1_data*)n->data;
			set<node*> rec;
			irectangle2 r;

			res->recurse_from_node(n, toprev, rec);
			r = node_set_bounding_rectangle(rec.begin(), rec.end());
			result.push_back(iipair(r.area(), nd->m));
			n = nd->next;
		}
	}
	sort(result.begin(), result.end(), greater<iipair>());
}

void save_classification(config_dictionary& cfg, const string& fname, string out)
{
    layer1_result* res = nullptr;
    string srcdir = cfg.get_value_string("src_dir", "");
    int layer = cfg.get_value_int("layer", 0);
    list<string> files;

    end_dir(srcdir);
    end_dir(out);
    list_directory(files, srcdir + fname);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        cout << "Processing: " << *fiter << endl;
        read_layer1_result(res, srcdir + *fiter);
        if (res) {
			ofstream os((out + change_extension(*fiter, ".m")).c_str());
			vector<iipair> result;

			sort_detections_by_size(result, res, layer);
			os << '{';
			for (vector<iipair>::iterator viter = result.begin(); viter != result.end(); ++viter) {
				if (viter != result.begin()) os << ',';
				os << '{' << viter->first << ',' << viter->second << '}';
			}
			os << '}' << endl;
			os.close();
			delete res;
		}
	}
}

void save_mathematica3(config_dictionary& cfg, const string& fname, string out)
{
    layer1_result* res = nullptr;
    string srcdir = cfg.get_value_string("src_dir", "");
    int layer = cfg.get_value_int("layer", 0);
    list<string> files;

    end_dir(srcdir);
    end_dir(out);
    list_directory(files, srcdir + fname);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        cout << "Processing: " << *fiter << endl;
        read_layer1_result(res, srcdir + *fiter);
        if (res) {
            ofstream os((out + change_extension(*fiter, ".m")).c_str());
            res->write_mma3(os, layer);
            os.close();
            delete res;
        }
    }
}

void save_node_info(config_dictionary& cfg, const string& fname, string out)
{
    layer1_result* res = nullptr;
    string srcdir = cfg.get_value_string("src_dir", "");
    int layer = cfg.get_value_int("layer", 0);
    list<string> files;

    end_dir(srcdir);
    end_dir(out);
    list_directory(files, srcdir + fname);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        read_layer1_result(res, srcdir + *fiter);
        if (res) {
            cout << "Processing: " << *fiter << endl;
            ofstream os((out + change_extension(*fiter, ".m")).c_str());
            if (layer >= 0 && layer <= res->max_layer_index())
                res->write_node_info(os, res->shape_nodes[layer]);
            os.close();
            delete res;
        }
    }
}

void display_library_info(config_dictionary& cfg, const string& fname)
{
    part_lib* library;

    read_library(fname, library);
    if (library) {
        cout << "Library '" << fname << "' info" << endl << endl;
        library->info(cout);
        delete library;
    }
}

void save_objects(config_dictionary& cfg, const string& fname, string out)
{
    typedef pair<double, bool> confidence_t;

    string dir = cfg.get_value_string("src_dir", "");
    list<string> patterns;
    layer1_result* res;
    part_lib* library;
    list<irectangle2> gtruths;
    string extension;
    string libname;
    int layer;
    double thresh;
    int object_layer;
    set<int> parts;
    vector<confidence_t> confidence;
    ofstream os(out.c_str());

    if (!os.is_open())
        return;

    end_dir(dir);
    cfg.get_value(thresh, "threshold", true);
    cfg.get_value(object_layer, "object_layer", true);
    cfg.get_value(extension, "extension", true);
    cfg.get_value(libname, "library", true);
    layer = object_layer - 1;
    bool first = true;

    read_library(libname, library);

    // print object parts
    set<int> oparts;
    library->get_part_types(oparts, object_layer);
    for (set<int>::iterator oiter = oparts.begin(); oiter != oparts.end(); ++oiter) {
        set<int> nparts;
        int i = *oiter;

        cout << i << " --> ";
        library->get_parent_parts(nparts, object_layer, i);
        for (set<int>::iterator siter = nparts.begin(); siter != nparts.end(); ++siter) {
            cout << *siter << ' ';
        }
        cout << endl;
    }

    library->get_part_types(parts, layer);

    os << '{';

    if (!cfg.get_value_bool("from_file", false) || !list_from_file(patterns, fname, dir))
        patterns.push_back(fname);
    for (list<string>::iterator pattern = patterns.begin(); pattern != patterns.end(); ++pattern) {
        list<string> files;
        
        list_directory(files, dir + *pattern + string("_*") + extension);
        for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
            read_layer1_result(res, dir + *file);
            if (res != nullptr) {
                cout << '.';
            
                list<layer1_result::box_data_t> boxes;

                res->get_boxes_inhibited(boxes, layer, parts, 0.8, thresh);
                for (list<layer1_result::box_data_t>::iterator biter = boxes.begin(); biter != boxes.end(); ++biter) {
                    if (first) first = false; else os << ',';
                    os << '{' << '\"' << *pattern << '\"' << ',' << biter->val << ',' << biter->size << ',' << biter->m << '}';
                }

                delete res;
            }
        }
    }
    os << '}' << endl;
    cout << endl;

    os.close();
}

void dilute_layer(const config_dictionary& cfg, const char* pattern)
{

    string dir = cfg.get_value_string("src_dir", "");
    string out_dir = cfg.get_value_string("out_dir", "");
    list<string> files;
    layer1_result* res;
    int layer;
    vector<int> parts;

    cfg.get_value(layer, "layer", true);
    cfg.get_value(parts, "parts", true);

    cout << "Keeping parts " << parts << " on layer " << layer << endl;
    end_dir(dir);
    end_dir(out_dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, pattern, dir))
        list_directory(files, dir + pattern);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cout << "Processing " << *file;
        read_layer1_result(res, dir + *file);
        if (res != nullptr) {
            res->dilute_layer(layer, parts);
            save_layer1_result(res, out_dir + *file);
            delete res;
            cout << " done";
        }
        cout << endl;
    }
}

void save_bounding_box(const config_dictionary& cfg, const char* pattern)
{

    string dir = cfg.get_value_string("src_dir", "");
    list<string> files;
    layer1_result* res;

    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, pattern, dir))
        list_directory(files, dir + pattern);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cout << "Processing " << *file;
        read_layer1_result(res, dir + *file);
        if (res != nullptr) {
            set<node*> nodes0;
            irectangle2 result;
            
            res->recurse(res->shape_nodes[res->max_layer_index()], atom("toPrevLayer"), nodes0);
            for (set<node*>::iterator niter = nodes0.begin(); niter != nodes0.end(); ++niter) {
                layer1_data* nd = (layer1_data*)(*niter)->data;
                if (nd->z == 0) result.eat(nd->x, nd->y);
            }
            
            ofstream os(change_extension(dir + *file, ".groundtruth").c_str());
            os << result << endl;
            os.close();

            delete res;
            cout << " done";
        }
        cout << endl;
    }
}

void track_layer(const config_dictionary& cfg, const char* pattern, string outname)
{
    string dir = cfg.get_value_string("src_dir", "");
    list<string> files;
    layer1_result* res;
    int slayer = cfg.get_value_int("src_layer", -1);
    int dlayer = cfg.get_value_int("dest_layer", 0);

    end_dir(dir);
    if (outname.empty()) outname = "tx.m";
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, pattern, dir))
        list_directory(files, dir + pattern);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cout << "Processing " << *file;
        read_layer1_result(res, dir + *file);
        if (res != nullptr) {
            string fname = change_extension(dir + *file, outname);
            ofstream os(fname.c_str());

            res->write_mma4(os, slayer, dlayer);
            os.close();
            delete res;
            
        }
        cout << endl;
    }
    
}

void save_support_images(layer1_result* res, const string& basename, int layer, int part)
{

    if (layer < 0 || layer > res->max_layer_index())
        return;

    unique_ptr<img> im(res->get_image_reconstructed(0, 0, vector<int>(), false, false, 1.0, false));
    int count = 0;

    if (im == nullptr) 
        return;
    res->recurse_and_link(layer, atom("toPrevLayer"), atom("toLayer0"));
    res->rename_edges(atom("toMissing"), atom("toLayer0"));
    for (auto niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;

        if (nd->z == layer && nd->m == part) {
            stringstream sstr;
            set<node*> support;
            set<ipoint2> points;
            img im2 = im->to_color_img();
            
            res->recurse_and_link(n, atom("toPrevLayer"), atom("toLayer0"), support);
            node_set_to_point_set(points, support.begin(), support.end());
            
            for (auto piter = points.begin(); piter != points.end(); ++piter) {
                im2.draw_big_point(piter->x, piter->y, COL_TO_REAL(255, 0, 0));
            }

            sstr << basename << "-" << layer << "-" << setw(3) << setfill('0') << part << "-" << 
                setw(3) << setfill('0') << count++ << ".png";
            im2.save(sstr.str());
        }

    }
}

void domen_gt(const config_dictionary& cfg, const char* pattern, string outname)
{
    string dir = cfg.get_value_string("src_dir", "");
    string cat = cfg.get_value_string("category", "");
    list<string> files;

    end_dir(dir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, pattern, dir))
        list_directory(files, dir + pattern);
    for (auto fiter = files.begin(); fiter != files.end(); ++fiter) {
        list<string> files2;
        int x, y, w, h, spec;
        bool ok;

        cout << "Processing " << *fiter << " ";

        ifstream is((dir + *fiter).c_str());
        is >> x >> y >> w >> h;
        ok = is.good();
        is.close();

        if (ok) {
            list_directory(files2, dir + change_extension(*fiter, "_*.ly?"));
            for (auto fiter2 = files2.begin(); fiter2 != files2.end(); ++fiter2) {
                unique_ptr<layer1_result> res(read_layer1_result(dir + *fiter2));

                if (res == nullptr) 
                    continue;

                double factor = ((double)res->x_size(0) - 2*res->border)/res->original_width;

                ofstream os((dir + change_extension(*fiter2, ".groundtruth")).c_str());
                os << int_round(res->border + x*factor) << ' ';
                os << int_round(res->border + y*factor) << ' ';
                os << int_round(res->border + (w + x)*factor) << ' ';
                os << int_round(res->border + (h + y)*factor) << ' ';
                os << cat << endl;
                os.close();

                cout << '.';
            }
        }
        cout << endl;
    }

}


void save_support(const config_dictionary& cfg, const char* pattern, string outname)
{
    string dir = cfg.get_value_string("src_dir", "");
    string outdir;
    list<string> files;
    int layer;
    int part;

    outdir = cfg.get_value_string("out_dir", dir);
    cfg.get_value(layer, "layer", true);
    cfg.get_value(part, "part", true);
    end_dir(dir);
    end_dir(outdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, pattern, dir))
        list_directory(files, dir + pattern);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cout << "Processing " << *file;
        unique_ptr<layer1_result> res(read_layer1_result(dir + *file));

        if (res != nullptr) {
            string fname = change_extension(dir + *file, "");

            save_support_images(res.get(), fname, layer, part);
        }
        cout << endl;
    }
    
}

void save_part_reconstruction(config_dictionary& cfg, const string& fpatt, string outname)
{
    string srcdir = cfg.get_value_string("src_dir", "");
    int layer = cfg.get_value_int("layer", 0);
    string libname;

    part_lib* library;
    list<string> files; 

    cfg.get_value(libname, "library", true);
    read_library(libname, library);

    if (library == nullptr) {
        cout << "Can not open library file '" << libname << "'" << endl;
        return;
    }

    end_dir(srcdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fpatt, srcdir))
        list_directory(files, srcdir + fpatt);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *fiter);
        if (res == nullptr) {
            cout << "Can not load '" << *fiter << "'" << endl;
        } else {
            if (layer >= 0 && layer <= res->max_layer_index()) {
                vector<node*>& s_nodes = res->shape_nodes[layer];
                for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
                    node* n = *iter;
                    layer1_data* nd = (layer1_data*)n->data;
                    node* p = library->parts[layer][nd->m];
                    part_data* pd = (part_data*)p->data;
                    
                    
                }
            }
            delete res;
        }
    }


    delete library;
}

void remove_texture(config_dictionary& cfg, const string& fpatt, string outname)
{
    string srcdir = cfg.get_value_string("src_dir", "");
    string outdir = cfg.get_value_string("out_dir", "");
    int layer = cfg.get_value_int("layer", 0);
    int radius = cfg.get_value_int("radius", 3);
    int threshold = cfg.get_value_int("threshold", 2);
    string prefix = cfg.get_value_string("prefix", "");

    list<string> files; 

    end_dir(srcdir);
    end_dir(outdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fpatt, srcdir))
        list_directory(files, srcdir + fpatt);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *fiter);
        if (res == nullptr) {
            cout << "Can not load '" << *fiter << "'" << endl;
        } else {
            cout << "Processing " << *fiter << "...";

			int removed = laydisplay::remove_texture(res, layer, radius, threshold);
            save_layer1_result(res, outdir + prefix + *fiter);

            cout << " (removed: " << removed << ")";
            cout << " done" << endl;
            delete res;
        }
    }
}

#include "layers/hoc.h"

void save_features(config_dictionary& cfg_root, const string& fpatt, string outname)
{
	config_dictionary cfg;
	cfg.from_namespace_priority(cfg_root, 1, "hoc");
	
	string g_outdir = outname;
	string outdir = outname;
	string srcdir = cfg.get_value_string("src_dir", "");

	bool read_windows_from_file = cfg.get_value_bool("read_windows_from_file", false);	
	string gtext = cfg.get_value_string("groundtruth_extension", ".groundtruth");

	string libname;
	cfg.get_value(libname, "library", true);

	part_lib* library;

    read_library(libname, library);
    if (library == nullptr) {
        cout << "Can not open library (library is needed to extract the number of features)." << endl;
        return;
    }

	using namespace hoc;

	hoc_histogram_generator hoc_generator(library);

	hoc_generator.init_cfg(cfg);
	

	end_dir(g_outdir);
	end_dir(outdir);

    int count = 0;

	list<string> files; 

    end_dir(srcdir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, fpatt, srcdir, cfg.get_value_string("pattern", "")))
        list_directory(files, srcdir + fpatt);

	// since we will be appending descriptors from different scale, we need to clear any existing values in all output files
	for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
		string::size_type scale_start = (*fiter).rfind("_");
		string name = (*fiter).substr(0, scale_start);
		
		stringstream ss;
		ss << outdir << name << "_det.txt";

		fclose(fopen(ss.str().c_str(),"w"));
	}

    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *fiter);
        if (res == nullptr) {
            cout << "Can not load '" << *fiter << "'" << endl;
        } else {
            cout << "Processing " << *fiter << "...";
			clock_t startt, endt;
			startt = clock();
			{

				string::size_type ext_start = (*fiter).rfind(".");
				string::size_type scale_start = (*fiter).rfind("_");

				string name = (*fiter).substr(0, scale_start);
				string scale_str = (*fiter).substr(scale_start+1);
				
				int scale;
				sscanf(scale_str.c_str(), "%d%*s", &scale);

				stringstream ss;
				ss << outdir << name << "_det.txt";

				ofstream os(ss.str().c_str(), ios_base::out | ios_base::app);

				list<pair<irectangle2, string> > rectangles;
				read_groundtruth(rectangles, srcdir + *fiter, "", gtext, false);

				int border = res->border;
				double scale_factor = res->original_width / ((double)(res->x_size(0) - 2*border));

				// Save (transformed) groundtruths
				if (!rectangles.empty()) {
					
					list<pair<irectangle2, string> > rectangles2(rectangles);
					ipoint2 bpt(border, border);

					stringstream gt_name;
					gt_name<< outdir << name << "_gt.txt";

					ofstream g_os(gt_name.str().c_str());

					for (list<std::pair<irectangle2,string> >::iterator iter = rectangles2.begin(); iter != rectangles2.end(); ++iter) {
						std::pair<irectangle2,string> p = *iter;
						irectangle2& r = p.first;
						string& cat_name = p.second;
				        
						r = (r - bpt) * scale_factor;					

						g_os << r.ll.x << ' ' << r.ll.y << ' ' << r.ur.x - r.ll.x << ' ' << r.ur.y - r.ll.y << ' ' << cat_name << endl;
					}
					g_os.close();

				}
				std::list<irectangle2> annotations_rects;

				if (read_windows_from_file) { 
						
					const string& annot_name =  srcdir + name + "_annot.txt";
					
					ifstream is(annot_name.c_str());

					if (is.is_open()) {
						while (true) {
							string line;
							// read and parse whole line as stream of string
							getline(is,line);
							stringstream sline(line);
							
							float x,y,w,h;
							sline >> x >> y >> w >> h;
							
							if (is.fail()) break;

							annotations_rects.push_back(irectangle2(x, y, x + w, y + h));
							
						}
					} else {
						cout << "ERROR: Unable to open annotation file '" << annot_name << "'. " << endl;
						return;
					}
				}
				// process image (res) for specific regions to produce list of descriptors (histograms of compositions)
				std::vector<histogram_descriptor>* result = hoc_generator.generate_descriptors(res, annotations_rects);

				for (std::vector<histogram_descriptor>::iterator desc = result->begin(); desc != result->end(); desc++) {
					
					int hist_size = (int)(*desc).hist.size();

					// skip any hitograms that have all zeros
					bool only_zeros = true;
					for (int f = 0; f < hist_size; ++f) {
						if ((*desc).hist[f] > 0) {
							only_zeros = false;
							break;
						}
					}

					if (only_zeros) 
						continue; // skiping since has all zeros (it is useless descriptor)

					os << scale << ' ' << (*desc).x << ' ' << (*desc).y << ' ' << (*desc).w << ' ' << (*desc).h << ' ';
					for (int f = 0; f < hist_size; ++f) {
						os << (*desc).hist[f] << ' ';
					}
					os << '\n';					
				}
				delete result;
				os.close();
            }
			endt = clock();
			cout << " (time = " << (double)(endt - startt)/CLOCKS_PER_SEC << ")" << endl;
            cout << endl;
            delete res;
            ++count;
			
        }
    }
	delete library;
	/*
	config_dictionary cfg;
	cfg.from_namespace_priority(cfg_root, 1, "hoc");
	
    string srcdir = cfg.get_value_string("src_dir", "");
    vector<int> layer;
    string mname; // matrix name
	string cname; // dense histogram cluster filename
    string libname;
    	
	bool read_windows_from_file;
    bool use_sw;
	bool use_img_coordinates;
	bool interpolate_to_neighboors;
	int border;
	bool normalize;
	bool include_similar_parts;
	double similarity_sigma;
	bool output_headers_only;
	bool use_part_merging;

	bool create_spatial_pyramid_matching;

	string gtext;
	string g_outdir = outname;
	string outdir = outname;

    cfg.get_value(layer, "layer", true);
    cfg.get_value(mname, "matrix", false);
	cfg.get_value(cname, "cluster_centers", false);
    cfg.get_value(libname, "library", true);	

    int sw_xstep = cfg.get_value_int("sliding_window_xstep", 0);
    int sw_ystep = cfg.get_value_int("sliding_window_ystep", 0);
	int sw_width = cfg.get_value_int("sliding_window_width", 0);
    int sw_height = cfg.get_value_int("sliding_window_height", 0);

	include_similar_parts = cfg.get_value_bool("parts_similarity.enable", false);
	similarity_sigma = cfg.get_value_bool("parts_similarity.sigma", 0.25);

	use_img_coordinates = cfg.get_value_bool("use_img_coordinates", true);	
	interpolate_to_neighboors = cfg.get_value_bool("interpolate_to_neighboors", false);	
	normalize = cfg.get_value_bool("normalize", false);	

	output_headers_only = cfg.get_value_bool("output_headers_only", false);	
	use_part_merging = cfg.get_value_bool("use_part_merging", true);	

	create_spatial_pyramid_matching = cfg.get_value_bool("create_spatial_pyramid_matching", false);	

	int spm_pyramid_layers = cfg.get_value_int("spatial_pyramid_matching.pyramid_layers", 0);
	int spm_toppart_layer = cfg.get_value_int("spatial_pyramid_matching.top_part_layer", 0);

	read_windows_from_file = cfg.get_value_bool("read_windows_from_file", false);	
	
	gtext = cfg.get_value_string("groundtruth_extension", ".groundtruth");
	
    use_sw = sw_xstep > 0 && sw_ystep > 0;

    part_lib* library;

    read_library(libname, library);
    if (library == nullptr) {
        cout << "Can not open library (library is needed to extract the number of features)." << endl;
        return;
    }

	// find min and max layer
	int min_layer = library->max_layer_index();
	int max_layer = 0;
	for (int i = 0; i < layer.size(); i++) {
		min_layer = min(min_layer, layer[i]);
		max_layer = max(max_layer, layer[i]);
	}
	if (layer.size() <= 0 || min_layer < 0 || max_layer > library->max_layer_index()) {
        cout << "Invalid layer (does not exist in library)." << endl;
        delete library;
        return;
    }

	histogram_clusters* cluster_centers = nullptr;

	if (cname.empty() == false) {
		// read all values into array				
		matrix<float> values_mat = read_float_matrix(cname);
		
		cluster_centers = new histogram_clusters(&values_mat[0], values_mat.width, values_mat.height);
	}

    vector<M_bin*> bin;
	for (int i = 0; i < layer.size(); i++) {
		bin.push_back(mname.empty() ? nullptr : new M_bin(mname));
	}

    if (mname.empty() == false && bin[0]->bin_count() == 0) {
        cout << "Cannot read \"bin\" matrix or matrix is invalid (" << mname << ")." << endl;
        delete library;
		for (int i = 0; i < bin.size(); i++)
			delete bin[i];
		if (cluster_centers != nullptr)
			delete cluster_centers;
        return;
    }

	histogram_distance_function* dist_func = nullptr;
	string dist_func_str;
	if (cfg.get_value(dist_func_str, "bin_distance_function", false)) {
		config_dictionary cfg_hist;
		// copy values based on namespace hierarhy
		cfg_hist.from_namespace_priority(cfg,1,"bin_distance_function");

		if (dist_func_str.compare("gaussian") == 0 ) {
			dist_func = new gaussian_distance_f(cfg_hist.get_value_int("sigma", 1.0, true));
		} else if (dist_func_str.compare("weighted_centers") == 0 ) {
			dist_func = new all_cetners_distance_f();
		} else if (dist_func_str.compare("round_overlap") == 0 ) {
			dist_func = new rounded_overlap_distance_f(cfg_hist.get_value_double("sigma", 1, true), cfg_hist.get_value_double("region_extend", 1.2, true));
		} else if (dist_func_str.compare("multikernel_gaussian") == 0 ) {
			dist_func = new multikernel_gaussian_distance_f(cfg_hist.get_value_double("sigma", 0.7, true), cfg_hist.get_value_double("threshold", 0.1, true));		
		} else if (dist_func_str.compare("none") != 0) {
			cout << "Invalid bin_distance_function (supported only: 'none', 'gaussian', 'weighted_centers')." << endl;
			delete library;
			for (int i = 0; i < bin.size(); i++)
				if (bin[i] != nullptr)
					delete bin[i];
			if (cluster_centers != nullptr)
				delete cluster_centers;
			return;
		}
	}

	vector<float> sw_xsteps(layer.size(), sw_xstep); vector<float> sw_ysteps(layer.size(), sw_ystep);
	vector<float> sw_widths(layer.size(), sw_width); vector<float> sw_heights(layer.size(), sw_height);

	vector<int> type_counts(layer.size(), 0);
	vector<double> contractions_used(layer.size(), 0);

	histogram_descriptor_settings extraction_settings;
	extraction_settings.interpolate_to_neighboors = interpolate_to_neighboors;
	extraction_settings.normalize = normalize;
	extraction_settings.dist = dist_func;

	extraction_settings.include_similar_parts = include_similar_parts;
	extraction_settings.similarity_sigma = similarity_sigma;

	extraction_settings.output_headers_only = output_headers_only;
	extraction_settings.use_part_merging = use_part_merging;

	extraction_settings.use_spm = create_spatial_pyramid_matching;
	extraction_settings.spm_pyramid_layers = spm_pyramid_layers;
	extraction_settings.spm_toppart_layer = spm_toppart_layer;

	for (int i = 0; i < layer.size(); i++) {
		type_counts[i] = library->layer_size(layer[i]);

		contractions_used[i] = library->contraction(1, layer[i] + 1);

		// if using image coordinates then transform all input pixels to specific layer coordinates
		if (use_img_coordinates == true) {
			cout << "using image coordinatres - contracted to " << contractions_used[i] << endl;
			sw_xsteps[i] = sw_xstep / contractions_used[i] ;
			sw_ysteps[i] = sw_ystep / contractions_used[i] ;
			sw_widths[i] = sw_width / contractions_used[i] ;
			sw_heights[i] = sw_height / contractions_used[i] ;
		} else {
			cout << "using original layer coordinates" << endl;
		}

		// also adjust bin size to sliding window size or to 
		if (bin[i] != nullptr) {
			if (sw_widths[i] > 0 && sw_heights[i] > 0 && (sw_widths[i] != bin[i]->width() || sw_height != bin[i]->height())) {
				M_bin* new_bin = bin[i]->get_resized(sw_widths[i], sw_heights[i]);
				delete bin[i];
				bin[i] = new_bin;
			} else if (use_img_coordinates) {
				M_bin* new_bin = bin[i]->get_resized(bin[i]->width() / contractions_used[i], bin[i]->height() / contractions_used[i]);
				delete bin[i];
				bin[i] = new_bin;
			} else if (layer.size() > 1) {
				cout << "Cannot use multiply layers without using image coordinates or using sliding window" << endl;
				delete library;
				for (int j = 0; j < bin.size(); j++)
					delete bin[j];
				if (cluster_centers != nullptr)
					delete cluster_centers;
				return;
			}
		}
	}
	end_dir(g_outdir);
	end_dir(outdir);

    int count = 0;

	list<string> files; 

    end_dir(srcdir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, fpatt, srcdir, cfg.get_value_string("pattern", "")))
        list_directory(files, srcdir + fpatt);

	// since we will be appending descriptors from different scale, we need to clear any existing values in all output files
	for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
		string::size_type scale_start = (*fiter).rfind("_");
		string name = (*fiter).substr(0, scale_start);
		
		stringstream ss;
		ss << outdir << name << "_det.txt";

		fclose(fopen(ss.str().c_str(),"w"));
	}

    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *fiter);
        if (res == nullptr) {
            cout << "Can not load '" << *fiter << "'" << endl;
        } else {
            cout << "Processing " << *fiter << "...";
			clock_t startt, endt;
			startt = clock();
			bool has_layers = true;
			for (int i = 0; i < layer.size() && has_layers; i++)
				has_layers = layer[i] <= res->max_layer_index();

			if (has_layers == false) {
                cout << "Layer does not exist, skipping...";
			} else {
				

				string::size_type ext_start = (*fiter).rfind(".");
				string::size_type scale_start = (*fiter).rfind("_");

				string name = (*fiter).substr(0, scale_start);
				string scale_str = (*fiter).substr(scale_start+1);
				
				int scale;
				sscanf(scale_str.c_str(), "%d%*s", &scale);

				stringstream ss;
				ss << outdir << name << "_det.txt";

				ofstream os(ss.str().c_str(), ios_base::out | ios_base::app);

				list<pair<irectangle2, string> > rectangles;
				read_groundtruth(rectangles, srcdir + *fiter, "", gtext, false);

				int border = res->border;
				double scale_factor = res->original_width / ((double)(res->x_size(0) - 2*border));

				// Save (transformed) groundtruths
				if (!rectangles.empty()) {
					
					list<pair<irectangle2, string> > rectangles2(rectangles);
					ipoint2 bpt(border, border);

					stringstream gt_name;
					gt_name<< outdir << name << "_gt.txt";

					ofstream g_os(gt_name.str().c_str());

					for (list<std::pair<irectangle2,string> >::iterator iter = rectangles2.begin(); iter != rectangles2.end(); ++iter) {
						std::pair<irectangle2,string> p = *iter;
						irectangle2& r = p.first;
						string& cat_name = p.second;
				        
						r = (r - bpt) * scale_factor;					

						g_os << r.ll.x << ' ' << r.ll.y << ' ' << r.ur.x - r.ll.x << ' ' << r.ur.y - r.ll.y << ' ' << cat_name << endl;
					}
					g_os.close();

				}
				std::list<irectangle2> annotations_rects;

				if (read_windows_from_file) { 
						
					const string& annot_name =  srcdir + name + "_annot.txt";
					
					ifstream is(annot_name.c_str());

					if (is.is_open()) {
						while (true) {
							string line;
							// read and parse whole line as stream of string
							getline(is,line);
							stringstream sline(line);
							
							float x,y,w,h;
							sline >> x >> y >> w >> h;
							
							if (is.fail()) break;

							annotations_rects.push_back(irectangle2(x, y, x + w, y + h));
							
						}
					} else {
						cout << "ERROR: Unable to open annotation file '" << annot_name << "'. " << endl;
						return;
					}
				}

				std::vector<laydisplay::histogram_descriptor>* merged_descriptors = nullptr;

				for (int ly_index = 0; ly_index < layer.size(); ly_index++) {
					
					double contraction_factor = (double)res->x_size(0) / (double)res->x_size(layer[ly_index]);

					std::list<frectangle2> rects;
					if (read_windows_from_file) { 
						
						for (std::list<irectangle2>::iterator it = annotations_rects.begin(); it != annotations_rects.end(); it++) {
							// copy rectangle with adjusted coordinates (a bit of overhead when use_img_coordinates == false)
							double x = (*it).ll.x;
							double y = (*it).ll.y;
							double w = (*it).ur.x - x;
							double h = (*it).ur.y - y;

							if (use_img_coordinates) {
								x = (x / scale_factor + border) /contraction_factor;
								y = (y / scale_factor + border) /contraction_factor;
								w = (w / scale_factor) /contraction_factor;
								h = (h / scale_factor) /contraction_factor;
							}
							rects.push_back(frectangle2(x, y, x + w, y + h));
						}

					} else if (use_sw) {
						// SLIDING window; uses sw_xstep and sw_ystep
						// Row format (all items separated with spaces): 
						//   no. of image (count) 
						//   coordinates of center of histogram
						//   histogram values
						//fpoint2 bin_size(bin[ly_index]->width(), bin[ly_index]->height());
						fpoint2 bin_size(sw_widths[ly_index], sw_heights[ly_index]);
						
						for (float xpos = 0; xpos + bin_size.x < res->x_size(layer[ly_index]); xpos += sw_xsteps[ly_index]) {
							for (float ypos = 0; ypos + bin_size.y < res->y_size(layer[ly_index]); ypos += sw_ysteps[ly_index]) {
								fpoint2 ll(xpos, ypos);
								rects.push_back(frectangle2(ll, ll + bin_size));
							}
						}
					} else {
						// this is no longer supported as nobody used it 
						cout << "Cannot create HoC without supplying region windows or using sliding window" << endl;
						delete library;
						for (int j = 0; j < bin.size(); j++)
							if (bin[j] != nullptr)
								delete bin[j];
						return;
						
						//fpoint2 bin_size(bin[ly_index]->width(), bin[ly_index]->height());						
						//fpoint2 ll = fpoint2(res->x_size(layer[ly_index])/2.0, res->y_size(layer[ly_index])/2) - bin_size/2.0 ;
						//rects.push_back(frectangle2(ll, ll + bin_size));
						
					}
			
					// this function works with specific layer[ly_index] coordinates and also returns results in layer[ly_index] coordinates
					std::vector<laydisplay::histogram_descriptor>* descriptors = laydisplay::get_part_histograms(library, res, layer[ly_index], bin[ly_index], cluster_centers, rects, extraction_settings);

					// merge with descriptors on previous layer unless this layer was processed first, then just use existing descriptor list
					if (merged_descriptors == nullptr) {
						// just use this descriptor list
						merged_descriptors = descriptors;
						
						// update coordinates if needed
						if (use_img_coordinates == true)  {
							for (std::vector<laydisplay::histogram_descriptor>::iterator merged_desc = merged_descriptors->begin(); merged_desc != merged_descriptors->end(); merged_desc++) {
								(*merged_desc).x = ((*merged_desc).x * contraction_factor - border) * scale_factor;
								(*merged_desc).y = ((*merged_desc).y * contraction_factor - border) * scale_factor;
								(*merged_desc).w = ((*merged_desc).w * contraction_factor) * scale_factor;
								(*merged_desc).h = ((*merged_desc).h * contraction_factor) * scale_factor;
							}
						}
					} else {
						std::vector<laydisplay::histogram_descriptor>::iterator merged_desc = merged_descriptors->begin();
						std::vector<laydisplay::histogram_descriptor>::iterator desc = descriptors->begin();
						while (desc != descriptors->end() && merged_desc != merged_descriptors->end()) {

							// update coordinates if needed
							if (use_img_coordinates == true)  {
								(*desc).x = ((*desc).x * contraction_factor - border) * scale_factor;
								(*desc).y = ((*desc).y * contraction_factor - border) * scale_factor;
								(*desc).w = ((*desc).w * contraction_factor) * scale_factor;
								(*desc).h = ((*desc).h * contraction_factor) * scale_factor;
							}
							
							// verify that bbox coordinates do match (at least with a few px of diff)
							if (::abs((*desc).x - (*merged_desc).x) > 3 ||
								::abs((*desc).y - (*merged_desc).y) > 3 || 
								::abs((*desc).w - (*merged_desc).w) > 3 ||
								::abs((*desc).h - (*merged_desc).h) > 3) {
								// error occured if one on same location or size
								throw exception();
							}
							
							// merge both descriptors by inserting new ones at the end
							(*merged_desc).hist.insert((*merged_desc).hist.end(), (*desc).hist.begin(), (*desc).hist.end());

							// continue to next descriptor
							desc++; merged_desc++;
						}

						if (desc != descriptors->end() || merged_desc != merged_descriptors->end()) {
							// error occured if one did not come to end
							throw exception();
						}							 

						delete descriptors;
					}
				}
				for (std::vector<laydisplay::histogram_descriptor>::iterator desc = merged_descriptors->begin(); desc != merged_descriptors->end(); desc++) {
					
					int hist_size = (int)(*desc).hist.size();

					// skip any hitograms that have all zeros
					bool only_zeros = true;
					for (int f = 0; f < hist_size; ++f) {
						if ((*desc).hist[f] > 0) {
							only_zeros = false;
							break;
						}
					}

					if (only_zeros) 
						continue; // skiping since has all zeros (it is useless descriptor)

					os << scale << ' ' << (*desc).x << ' ' << (*desc).y << ' ' << (*desc).w << ' ' << (*desc).h << ' ';
					for (int f = 0; f < hist_size; ++f) {
						os << (*desc).hist[f] << ' ';
					}
					os << '\n';					
				}
				delete merged_descriptors;
				os.close();
            }
			endt = clock();
			cout << " (time = " << (double)(endt - startt)/CLOCKS_PER_SEC << ")" << endl;
            cout << endl;
            delete res;
            ++count;
			
        }
    }
    if (dist_func != nullptr)
		delete dist_func;
	for (int i = 0; i < bin.size(); i++)
		if (bin[i] != nullptr)
			delete bin[i];
    delete library;
	if (cluster_centers != nullptr)
		delete cluster_centers;
		*/
}

#include <functional>

struct max_function_class : public binary_function<double,double,double> {
  double operator() (double a, double b) const {return std::max(a,b);}
};

/*
void mean_shift(layer1_result* res, int layer, double sigma, double eps = 1.0E-6, int maxsteps = 100)
{
    if (layer < 0 || layer >= res->max_layer_index() || res->shape_nodes[layer].empty());
        return;

    if (!res->grid(layer)) 
        res->init_grid(layer);

    matrix<bool> m(res->x_size(layer), res->y_size(layer), false);
    vector<node*> nvec;

    for (list<node*>::iterator iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
        if (node_layer(*iter) == layer) {
            nvec.push_back(*iter);
            (*iter)->set_attr(HIDDEN_NODE_ATTR); // We hide all nodes by default and "unhide" them below
        }
    }
    sort(nvec.begin(), nvec.end(), response_sort_f(G_RESPONSE));

    double c = -1.0/(sigma*sigma*2);
    double radius = 2.0*sigma;
    double radius2 = radius*radius;
    double eps2 = eps*eps;

    for (vector<node*>::iterator nviter = nvec.begin(); nviter != nvec.end(); ++nviter) {
        node* n = *nviter;
        layer1_data* nd = (layer1_data*)n->data;

        if (m(nd->x, nd->y))
            continue;

        dpoint2 X(nd->x, nd->y);
        vector<dpoint2> track;
        double delta2 = numeric_limits<double>::infinity();

        track.push_back(X);
        for (int s = 0; s < maxsteps; ++s) {
            dpoint2 num = dpoint2::zero;
            double denom = 0.0;
            int minx = std::max<int>(0, nd->x - radius);
            int maxx = std::min<int>(nd->x + radius + 1, res->x_size(layer));
            int miny = std::max<int>(0, nd->y - radius);
            int maxy = std::min<int>(nd->y + radius + 1, res->y_size(layer));

            for (int x = minx; x < maxx; ++x) {
                for (int y = miny; y < maxy; ++y) {
                    dpoint2 P(x, y);
                    double dxn = (P - X).norm2();

                    if (dxn <= radius2) {
                        double k = ::exp(c*dxn);
                        num += P*k;
                        denom += k;
                    }
                }
            }
            num.x /= denom;
            num.y /= denom;
            delta2 = X.distance2(num);
            X = num;
            track.push_back(X);
            if (delta2 < eps2)
                break;
        }

        // Find closest point in res to each of the points in the "track".
        for (int t = 0; t < (int)track.size(); ++t) {
            ipoint2 p = get_closest_position(res, layer, track[t], 0, radius);

            m(p.x, p.y) = true;
        }

        vector<node*> cnodes;
        ipoint2 c = get_closest_position(res, layer, track.back(), 0, radius);

        res->nodes_at(cnodes, c.x, c.y, layer);
        res->clear_attr(cnodes.begin(), cnodes.end(), HIDDEN_NODE_ATTR);
    }
}

*/

vector<ipoint2> points_from_pca_data(const pca_data& pcad, double factor = 50)
{
    vector<dpoint2> pts = partition(pcad.mean);
    vector<ipoint2> result;
    double f = pcad.sizefactor == 1 ? factor : pcad.sizefactor;

    result.reserve(pts.size());
    for (auto piter = pts.begin(); piter != pts.end(); ++piter) {
        result.push_back(ipoint2((int)(f*piter->x), (int)(f*piter->y)));
    }
    return result;
}

dpoint2 transfer_vector(dpoint2 srct, dpoint2 srcv, dpoint2 destt)
{
    normalize(srct);
    normalize(destt);
    
    if (acos(dot_product(srct, destt)) > M_PI_2) 
        destt = -destt;

    double delta = asin(planar_cross_product(srct, destt));
    
    return rotate_point(srcv, delta);
}


map<int, ipoint2> geo_pieces_map(const vector<pair<int, ipoint2> >& pcs)
{
	map<int, vector<ipoint2> > ptsmap;

	for (auto giter = pcs.begin(); giter != pcs.end(); ++giter) 
		ptsmap[giter->first].push_back(giter->second);

	map<int, ipoint2> result;

	for (auto miter = ptsmap.begin(); miter != ptsmap.end(); ++miter) 
		result[miter->first] = ipoint2::center(miter->second.begin(), miter->second.end());
	return result;
}

bool identity_permutation(const vector<int>& perm)
{
	for (int i = 0; i < (int)perm.size(); ++i) {
		if (perm[i] != i) 
			return false;
	}
	return true;
}

vector<ipoint2> geo_pieces_vector(const vector<pair<int, ipoint2> >& pcs)
{
	if (pcs.empty()) 
		return vector<ipoint2>();

	map<int, vector<ipoint2> > ptsmap;
	int maxi = INT_MIN;

	for (auto giter = pcs.begin(); giter != pcs.end(); ++giter) {
		ptsmap[giter->first].push_back(giter->second);
		if (giter->first > maxi) 
			maxi = giter->first;
	}

	vector<ipoint2> result(maxi + 1);

	for (auto miter = ptsmap.begin(); miter != ptsmap.end(); ++miter) 
		result[miter->first] = ipoint2::center(miter->second.begin(), miter->second.end());
	return result;
}


// Return subpart coordinates
//
vector<ipoint2> subpart_coordinates(node* p)
{
	int ename = atom("lyrSrc");
	vector<ipoint2> result;

	foreach_neighbor(p, ename, niter) {
		part_data_2* pd = (part_data_2*)neighbor_edge_data(niter);
		if (pd->index >= result.size()) result.resize(pd->index + 1); 
		result[pd->index].set(pd->x, pd->y);
	}
	return result;
}

// Synchronize indices in app: int -> (*vector<int>*, double)
// i.e. indices in parse graph (synchonized due to or-ing)
void update_path_permutations(part_lib* library)
{
	if (library == nullptr) 
		return;

	int lysim = atom("lyrSimilar");

	for (int l = 1; l <= library->max_layer_index(); ++l) {
		for (int pi = 0; pi < library->layer_size(l); ++pi) {
			node* p = library->parts[l][pi];

			if (p->is_attr_set(CATEGORY_NODE_ATTR))
				continue;

			forall_neighbors(p, niter) {
				if (neighbor_index(niter) != atom("lyrCenter") && neighbor_index(niter) != atom("lyrSrc")) 
					continue;

				node* pn = neighbor_node(niter);
				lib_data* pnd = (lib_data*)neighbor_node_data(niter);
				part_data_2a* pned = (part_data_2a*)neighbor_edge_data(niter);
				vector<ipoint2> pngmap = geo_pieces_vector(get_library_geo_pieces(pn, pnd->layer));

				// Handle co-occurrence or-ing (done on object layer)
				for (auto giter = pned->geo.begin(); giter != pned->geo.end(); ++giter) {
					auto aiter = pned->app.find(giter->first);

					if (aiter == pned->app.end()) {
						cout << "Inconsistency between geo and app:" << endl;
						cout << "layer: " << l << ", part: " << pi << ", geo target: " << giter->first << endl;
						throw;
					}

					node* tp = library->parts[pnd->layer][giter->first];
					vector<ipoint2> tpgmap = geo_pieces_vector(get_library_geo_pieces(tp, pnd->layer));
					vector<int> perm = point_matching(pngmap, tpgmap);

					if (identity_permutation(perm)) aiter->second.first.clear();	// = vector<int>();
					else aiter->second.first = perm;
				}

				// handle similarity or-ing
				// it will be handled below
				/*foreach_neighbor(pn, tosim, sniter) {
					node* tp = neighbor_node(sniter);
					lib_data* tpd = (lib_data*)neighbor_node_data(sniter);

					if (pned->app.find(tpd->type) != pned->app.end()) 
						continue;

					vector<ipoint2> tpgmap = geo_pieces_vector(get_library_geo_pieces(tp, tpd->layer));
					vector<int> perm = point_matching(pngmap, tpgmap);
					if (identity_permutation(perm)) pned->app[tpd->type].first.clear();
					else pned->app[tpd->type].first  = perm;
					pned->app[tpd->type].second = 1.0;
				}*/

			}

			// handle similarity or-ing
			// add edge data with permutation on the lyrSimilar edge (which is supposed to be null)
			// hopefully all permutations will be identities...
			vector<ipoint2> pgmap = geo_pieces_vector(get_library_geo_pieces(p, l));

			foreach_neighbor(p, lysim, siter) {
				node* pn = neighbor_node(siter);
				lib_data* pnd = (lib_data*)neighbor_node_data(siter);
				vector<ipoint2> pngmap = geo_pieces_vector(get_library_geo_pieces(pn, pnd->layer));
				vector<int> perm = point_matching(pgmap, pngmap);

				if (identity_permutation(perm)) neighbor_edge_data(siter) = new perm_edge_data(vector<int>());
				else neighbor_edge_data(siter) = new perm_edge_data(perm);
			}
		}
	}
}

void test(config_dictionary& cfg, const string& fpatt, string out)
    //  parse model
{
	int layer = cfg.get_value_int("layer", 1);
	int maxpi = cfg.get_value_int("max_parts", 1);

	cout << "Exporting paths for layer: " << layer << endl << endl;

	unique_ptr<part_lib> library(read_library("C:\\work\\data\\multiclass\\Ferrari\\mug\\2\\olibm.plb"));
//=======    cout << "Files matching the pattern " << fpatt << ": " << endl;
//	find_files(fnames, "c:\\work\\data\\multiclass\\Pascal\\bike\\", fpatt);
//    for (auto fiter = fnames.begin(); fiter != fnames.end(); ++fiter) {
//        cout << *fiter << endl;
//    }
//>>>>>>> .theirs	*/
	for (int pi = 0; pi < library->layer_size(layer) && pi < maxpi; ++pi) {
		node* p = library->parts[layer][pi];

		cout << "Part #" << pi << endl;
		forall_neighbors(p, niter) {
			if (neighbor_index(niter) == atom("lyrCenter") || neighbor_index(niter) == atom("lyrSrc")) {
				part_data_2a* pd = (part_data_2a*)neighbor_edge_data(niter);
				int count = 0;
				for (auto giter = pd->geo.begin(); giter != pd->geo.end(); ++giter) {
					cout << "edge name:" << neighbor_index(niter);
					cout << "variant: " << count << " target:" << giter->first << " paths: ";
					for (auto piter = giter->second.begin(); piter != giter->second.end(); ++piter) {
						cout << "{{";
						for (auto viter = piter->first.begin(); viter != piter->first.end(); ++viter) {
							if (viter != piter->first.begin()) cout << ',';
							cout << *viter;
						}
						cout << "}, ";
						cout << "{" << piter->second.p.x << ',' << piter->second.p.y << "}},";
					} 
					cout << endl;
					++count;
				}
			}
		}
		
	}

	update_path_permutations(library.get());

	//unique_ptr<part_lib> library(read_library(fpatt));
 //   int part, layer;
 //   
 //   if (library == nullptr) {
 //       cout << "Can not read library " << fpatt << endl;
 //       return;
 //   }

 //   cfg.get_value(part, "part", true);
 //   cfg.get_value(layer, "layer", true);
 //   if (layer < 0 || layer > library->max_layer_index() || part < 0 || part >= library->layer_size(layer)) {
 //       cout << "Wrong layer/part." << endl;
 //       return;
 //   }
 //   
 //   node* p = library->parts[layer][part];
 //   map<int, vector<int> > pmap;
 //   

    //part_paths(p, result);
}

void save_gt_test(config_dictionary& cfg, const string& fpatt, string out)
{
    list<string> files;
    string srcdir = cfg.get_value_string("src_dir", "");
    string extension = cfg.get_value_string("extension", "_?.ly7");
    string catname = cfg.get_value_string("cat_name", "bike");

    end_dir(srcdir);
    list_directory(files, srcdir + fpatt);
    for (auto fiter = files.begin(); fiter != files.end(); ++fiter) {
        string fname = *fiter;
        string basename = change_extension(fname, "");
        list<string> pfiles;
        list<pair<irectangle2, string> > rectangles;

        list_directory(pfiles, srcdir + basename + extension);
        read_groundtruth(rectangles, srcdir + fname);
        for (auto pfiter = pfiles.begin(); pfiter != pfiles.end(); ++pfiter) {
            unique_ptr<layer1_result> res(read_layer1_result(srcdir + *pfiter));

            if (res != nullptr) {
                list<pair<irectangle2, string> > newrectangles;

                for (auto riter = rectangles.begin(); riter != rectangles.end(); ++riter) {
                    irectangle2 r = riter->first;
                    double f = ((double)(res->x_size(0) - 2*res->border))/res->original_width;
                    
                    r.ur.x = int_round(r.ur.x*f + res->border);
                    r.ur.y = int_round(r.ur.y*f + res->border);
                    r.ll.x = int_round(r.ll.x*f + res->border);
                    r.ll.y = int_round(r.ll.y*f + res->border);
                    if (riter->second == "") newrectangles.push_back(pair<irectangle2, string>(r, catname));
                    else newrectangles.push_back(pair<irectangle2, string>(r, riter->second));
                }
                save_groundtruth(newrectangles, srcdir, *pfiter, "", "", true);
            }
        }
    }
}

void testxxxx(config_dictionary& cfg, const string& fpatt, string out)
{
    list<int> l;

	for (int i = 0; i < 10; ++i) l.push_back(i);
	random_sublist(l, 4);
	for (auto liter = l.begin(); liter != l.end(); ++liter) {
		cout << *liter << ' ';
	}
}

void testx(config_dictionary& cfg, const string& fpatt, string out)
{
	string dir = fpatt;
	string category;
	
	cfg.get_value(category, "category", true);

	end_dir(dir);

	vector<double> mu, valid, vx, w, selection;

	load_matlab(mu, dir + category + "-mu.txt");
	load_matlab(valid, dir + category + "-valid.txt");
	load_matlab(vx, dir + category + "-vx.txt");
	load_matlab(w, dir + category + "-w.txt");
	load_matlab(selection, dir + category + "-selection.txt");

	if (mu.empty() || valid.empty() || vx.empty() || w.empty()) {
		cout << "One of vectors is empty or file(s) do(es) not exists." << endl;
		return;
	}
}


void test0(config_dictionary& cfg, const string& fpatt, string out)
{
    unique_ptr<layer1_result> res(read_layer1_result("D:\\work\\data\\multiclass\\Ferrari\\all\\layerx\\mug\\4\\test_m_muki_1.ly7"));

    if (res == nullptr)
        return;

    vector<node*> lnodes;

    res->get_layer_nodes(lnodes, 5);
    cout << "# of nodes: " << lnodes.size() << endl;

    map<node*, vector<node*> > clusters;
    int sigma = cfg.get_value_double("sigma", 2.0);

	clock_t start = clock();
    cluster_detections_ms(clusters, res.get(), list<node*>(lnodes.begin(), lnodes.end()), sigma, 100);
	clock_t end = clock();

    cout << "# of clusters: " << clusters.size() << endl;
	cout << " (time = " << (double)(end - start)/CLOCKS_PER_SEC << ")" << endl;

    ofstream os("test_clustering.m");
    int count = 0;

    for (auto cliter = clusters.begin(); cliter != clusters.end(); ++cliter) {
        set<node*> nset;

        res->recurse_and_link(cliter->first, atom("toPrevLayer"), atom("toLayer0"), nset);
        irectangle2 rect0 = bounding_rectangle_of_nodes(nset.begin(), nset.end());

        os << count << ' ' << rect0.ll.x << ' ' << rect0.ll.y << ' ' << rect0.ur.x << ' ' << rect0.ur.y << endl;

        for (auto ccliter = cliter->second.begin(); ccliter != cliter->second.end(); ++ccliter) {
            set<node*> nset1;

            res->recurse_and_link(*ccliter, atom("toPrevLayer"), atom("toLayer0"), nset1);
            irectangle2 rect1 = bounding_rectangle_of_nodes(nset1.begin(), nset1.end());

            os << count << ' ' << rect1.ll.x << ' ' << rect1.ll.y << ' ' << rect1.ur.x << ' ' << rect1.ur.y << endl;
        }

        ++count;
    }

    os.close();
    

}


// Matches points 'ptsm' to mean of pca in node 'p'. 
// If 'p' is not a "vs part" an empty vector is returned.
vector<ipoint2> match_to_model(const vector<ipoint2>& ptsm, node* p)
{
    vs_part_data* vspd = dynamic_cast<vs_part_data*>(p->data);

    if (vspd == nullptr || vspd->pcad.mean.empty()) 
        return vector<ipoint2>();

    vector<dpoint2> pts;
    vector<ipoint2> dptsm;
    vector<dpoint2> v = partition(vspd->pcad.mean);

    pts = get_resized_vector(cast_vector<dpoint2, ipoint2>(ptsm), (int)v.size());
    dptsm = cast_vector<ipoint2, dpoint2>(pts);
    translate_and_scale(pts);

    permute(dptsm, point_matching(pts, v));
    return dptsm;
}

vector<ipoint2> match_to_model(node* p, layer1_result* res, node* n)
{
    vs_part_data* vspd = dynamic_cast<vs_part_data*>(p->data);

    if (vspd == nullptr) 
        return vector<ipoint2>();

    vector<pair<int, ipoint2> > ipts = get_library_geo_pieces(p, 5);
    vector<pair<int, ipoint2> > ptsm;

    get_node_geo(ptsm, res, n);

    if (ipts.empty()) {
        vector<dpoint2> pts;
        vector<ipoint2> dptsm;
        vector<dpoint2> v = partition(vspd->pcad.mean);

        for (auto piter = ptsm.begin(); piter != ptsm.end(); ++piter) {
            pts.push_back((dpoint2)piter->second);
        }
        pts = get_resized_vector(pts, (int)v.size());
        dptsm = cast_vector<ipoint2, dpoint2>(pts);
        translate_and_scale(pts);

        permute(dptsm, point_matching(pts, v));
        return dptsm;
    } else {
        ptsm = inhibit_point_set(ptsm, 5);

        vector<int> match = piecewise_point_matching_p(ptsm, ipts);

        set<int> matchedpcs, allpcs;
        vector<ipoint2> dptsm(ipts.size());

        for (int i = 0; i < (int)dptsm.size(); ++i) {
            int j = match[i];

            if (j >= 0 && j < (int)ptsm.size()) {
                dptsm[i] = ptsm[j].second;
                matchedpcs.insert(ipts[i].first);
            } else {
                //dptsm[i] = ipts[i].second;    
                dptsm[i] = ipoint2(INT_MAX, INT_MAX);
            }
            allpcs.insert(ipts[i].first);
        }
        return dptsm;
    }
}

void merge_scales(config_dictionary& cfg, const string& fname, string outname)
{
    string dir = cfg.get_value_string("src_dir", "");
    string outdir = cfg.get_value_string("out_dir", "");
    vector<int> dscales;
    set<int> dropscales;

    cfg.get_value(dscales, "drop_scales");
    dropscales.insert(dropscales.begin(), dropscales.end());

    list<string> files;
    map<string, list<string> > groups;
    string ext;

    end_dir(dir);
    end_dir(outdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, dir))
        list_directory(files, dir + fname);
    for (auto fiter = files.begin(); fiter != files.end(); ++fiter) {
        string f = change_extension(*fiter, "");

        if (f.empty()) 
            continue;
            
        int scale = (int)(f.back() - '0');

        if (dropscales.find(scale) != dropscales.end())
            continue;

        string basename = change_extension(f, "", "_");
        groups[basename].push_back(*fiter);
        ext = get_extension(*fiter, ".");
    }

    for (auto giter = groups.begin(); giter != groups.end(); ++giter) {
        list<string>& group = giter->second;

        if (group.empty())
            continue;

        cout << "Merging: " << giter->first << "...";

        group.sort();
        
        layer1_result* res0 = nullptr;

        for (auto fiter = group.begin(); fiter != group.end(); ++fiter) {
            layer1_result* res = read_layer1_result(dir + *fiter);

            if (res0 == nullptr) res0 = res;
            else {
                res0->merge(res, res0->border);
                delete res;
            }
        }
        if (res0 != nullptr) {
            res0->save(outdir + giter->first + "_m" + ext);
            delete res0;
        }
        cout << endl;
    }

}

void eccv_export(config_dictionary& cfg, const string& fpatt, string outname)
{
    unique_ptr<part_lib> library(read_library(fpatt));
    int layer = cfg.get_value_int("object_layer", 5);
    int src = atom("lyrSrc");

    library->update_part_data_contractions();
	
    for (int pi = 0; pi < library->parts[layer].size(); ++pi) {
        node* p = library->parts[layer][pi];
        part_data* pd = (part_data*)p->data;
        double factor = pd->total_contraction(layer);
        ipoint2 cdelta(-pd->cmx, -pd->cmy);
        cv::Mat3b im(500, 500, cv::Vec3b(0, 0, 0));
        
        foreach_neighbor(p, src, iter) {
            part_data_2* ed = (part_data_2*)neighbor_edge_data(iter);
            vs_part_data* vspd = dynamic_cast<vs_part_data*>(neighbor_node_data(iter));

            if (vspd == nullptr) { 
                cout << "Not a vspart" << endl;
                throw;
            }
            cv::Mat3b imp = pca_data_to_image(vspd->pcad, 2, 5);
            ipoint2 delta = (cdelta + ipoint2(ed->x, ed->y)) * factor;

            cv::imwrite("A.png", imp);

            for (int r = 0; r < imp.rows; ++r) 
                for (int c = 0; c < imp.cols; ++c) {
                    cv::Vec3b v = imp.at<cv::Vec3b>(r, c);

                    if (v[0] != 0 || v[1] != 0 || v[2] != 0) {
                        cv::Vec3b& vv = im.at<cv::Vec3b>(im.rows/2 + delta.y - imp.rows/2 + r, im.cols/2 + delta.x - imp.cols/2 + c);
                        vv[0] = v[0]; vv[1] = v[1]; vv[2] = v[2];
                    }
                }
        }

        stringstream sstr;

        sstr << outname << "-" << setw(4) << setfill('0') << pi << ".png";
        cv::imwrite(sstr.str(), im);
    }

}

// Returns a vector of nodes on layer 'layer' of 'res' which covers 'mincover' fraction of
// the whole cover of layer 'layer' and satisfy the conditions of filter 'rf'.
// Intersections between parts must be < 'intthresh' of its size.
vector<node*> layer_cover(layer1_result* res, int layer, const response_filter& rf, double intthresh = 0.5,
    double mincover = 1.0)
{
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");
    vector<node*> nodes;
    set<ipoint2> ly0pts;
    map<node*, set<ipoint2> > recmap;
    vector<pair<double, node*> > valvec;
    vector<node*> result;

    res->get_layer_nodes(nodes, layer);
    for (auto niter = nodes.begin(); niter != nodes.end(); ++niter) {
        node* n = *niter;
        layer1_data* nd = (layer1_data*)n->data;

        if (!rf.check(nd))
            continue;

        set<node*> rec;
        set<ipoint2>& prec = recmap[n];

        res->recurse_and_link(n, toprev, to0, rec);
        node_set_to_point_set(ly0pts, rec.begin(), rec.end());
        node_set_to_point_set(prec, rec.begin(), rec.end());
        
        irectangle2 rect = irectangle2::bounding_rectangle(prec.begin(), prec.end());

        valvec.push_back(pair<double, node*>(rect.area(), n)); // rec.size() * nd->val() ????
    }
    sort(valvec.begin(), valvec.end(), greater<pair<double, node*> >());

    set<ipoint2> covpts;
    
    for (auto vviter = valvec.begin(); vviter != valvec.end(); ++vviter) {
        node* n = vviter->second;
        set<ipoint2>& rec = recmap[n];

        if (intersection_size(rec, covpts) < intthresh*rec.size()) {
            result.push_back(n);
            covpts.insert(rec.begin(), rec.end());
        }
        if (intersection_size(ly0pts, covpts) >= mincover*ly0pts.size()) 
            break;
    }
    return result;
}


void sample_part(config_dictionary& cfg, const string& fpatt, string outname)
{
    int layer, part;
    int samples = cfg.get_value_int("samples", 100);
    int picsize = cfg.get_value_int("image_size", 300);

    if (outname.empty() || outname == "dummy") {
        outname = "%s_part-%03d_%03d.png";
    }

    part_lib* library;

    cfg.get_value(part, "part");
    cfg.get_value(layer, "layer");
    read_library(fpatt, library);
    if (library == nullptr) {
        cout << "Can not read library." << endl;
        return;
    }
    library->update_part_data_contractions();
    if (layer < 0 || layer > library->max_layer_index() || part < 0 || part >= library->layer_size(layer))
        cout << "Invalid 'layer' or 'part' value." << endl;
    else {
        node* p = library->parts[layer][part];
        part_data* pd = (part_data*)p->data;
        vector<img> gabors(library->layer_size(0));

        for (int i = 0; i < library->layer_size(0); ++i)
            gabors[i] = library->get_image(0, i);

        max_function_class maxf;
	
        for (int n = 0; n < samples; ++n) {
            vector<pair<dpoint2, int> > lpts = random_sample_part(library, p);
            vector<dpoint2> pts = extract_first<dpoint2>(lpts.begin(), lpts.end());
            dpoint2 ptsc = dpoint2::center(pts.begin(), pts.end());
            img im(picsize, picsize, 0.0, true);

            for (int i = 0; i < (int)lpts.size(); ++i) {
                dpoint2 p = lpts[i].first - ptsc + dpoint2(picsize/2, picsize/2);
                img& gim = gabors[lpts[i].second];

                if (p.x >= 0 && p.y >= 0 && p.x < picsize && p.y < picsize) {
                    //im.draw_big_point(p.x, p.y, 1.0);
                    im.blt_central(gim, (int)gim.width/2, (int)gim.height/2, (int)p.x, (int)p.y, maxf);
                }
            }

            char fname[1000];

            sprintf(fname, outname.c_str(), change_extension(fpatt, "").c_str(), pd->type, n);

            im.save(fname, -150);
        }
    }
    delete library;
}

irectangle2 get_rectangle_at_layer(layer1_result* res, node* n, int layer)
{
    set<node*> nodes;
    irectangle2 r;

    res->subgraph_from_node(n, atom("toPrevLayer"), nodes);
    for (set<node*>::iterator iter = nodes.begin(); iter != nodes.end(); ++iter) {
        node* n = *iter;
        layer1_data* nd = (layer1_data*)n->data;

        if (nd->z == layer) 
            r.eat(nd->x, nd->y);
    }
    return r;
}


void write_ellipse(ostream& os, const ellipse& ell)
{
    os << '{';
    os << ell.cx << ',' << ell.cy;
    os << ',' << cos(2 * ell.angle);
    os << ',' << sin(2 * ell.angle);
    os << ',' << 3 * ell.a;
    os << ',' << ell.eccentricity();
    os << '}';
}

void save_ellipses(config_dictionary& cfg, const string& fpatt, string outname)
{
    typedef map<int, set<node*> > psmap_t;

    string srcdir = cfg.get_value_string("src_dir", "");
    string outdir = cfg.get_value_string("out_dir", "");
    int layer = cfg.get_value_int("layer", 5);
    int to_prev = atom("toPrevLayer");
    int to_0 = atom("toLayer0");

    list<string> files; 

    end_dir(srcdir);
    end_dir(outdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fpatt, srcdir))
        list_directory(files, srcdir + fpatt);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *fiter);
        if (res == nullptr) {
            cout << "Can not load '" << *fiter << "'" << endl;
        } else {
            cout << "Processing " << *fiter << "...";

            ofstream os(change_extension(srcdir + *fiter, "_ell.m").c_str());
            vector<node*>& s_nodes = res->shape_nodes[layer];

            cout << "<R";
            //res->add_reconstruction_edges_fwd(category_layer - 1);
            res->add_reconstruction_edges_leq_fwd(layer - 1);
            cout << ">";

            os << fixed;
            os << '{';
            for (vector<node*>::iterator niter = s_nodes.begin(); niter != s_nodes.end(); ++niter) {
                node* on = *niter;
                layer1_data* ond = (layer1_data*)on->data;
                set<node*> objnodes;
                ellipse ell;
                psmap_t psmap;
                
                // Update statistics for subparts
                foreach_neighbor(on, to_prev, iter) {
                    node* sn = neighbor_node(iter);
                    layer1_data* snd = (layer1_data*)sn->data;
                    edge_data_name* nned = (edge_data_name*)neighbor_edge_data(iter);
                    set<node*> nset;
                    psmap_t::iterator miter;
                    
                    if (nned == nullptr) {
                        cout << "No edge name, exiting.";
                        throw exception();
                    }
                    sn->get_neighbor_set(to_0, nset);
                    objnodes.insert(nset.begin(), nset.end());
                    miter = psmap.find(nned->index);
                    if (miter == psmap.end()) psmap.insert(psmap_t::value_type(nned->index, nset));
                    else miter->second.insert(nset.begin(), nset.end());
                }

                fit_ellipse_to_nodes(ell, objnodes.begin(), objnodes.end());
                //add_to_map(!isfalse, ond->m, index, ipoint2(INT_MIN, INT_MIN), ell);
                if (niter != s_nodes.begin()) os << ',';
                os << endl << '{';
                os << ond->m << ',';
                write_ellipse(os, ell);
                os << ',' << '{';
                
                for (psmap_t::iterator miter = psmap.begin(); miter != psmap.end(); ++miter) {
                    fit_ellipse_to_nodes(ell, miter->second.begin(), miter->second.end());
                    if (miter != psmap.begin()) os << ',';
                    write_ellipse(os, ell);
                    //add_to_map(!isfalse, ond->m, index, miter->first, ell);
                }
                os << '}' << '}';
            }
            os << '}' << endl;
            os.close();

            cout << " done" << endl;
            delete res;
        }
    }
}


void save_sanja(config_dictionary& cfg, const string& fpatt, string outname)
{
    string srcdir = cfg.get_value_string("src_dir", "");
    string outdir = cfg.get_value_string("out_dir", "");

    list<string> files; 

    end_dir(srcdir);
    end_dir(outdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fpatt, srcdir))
        list_directory(files, srcdir + fpatt);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *fiter);
        if (res == nullptr) {
            cout << "Can not load '" << *fiter << "'" << endl;
        } else {
            cout << "Processing " << *fiter << "...";

            ofstream os(change_extension(outdir + *fiter, ".txt").c_str());
            vector<node*>& s_nodes = res->shape_nodes[0];


            os << s_nodes.size() << endl; 
            os << 0 << endl; // no edges
            for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
                layer1_data* nd = (layer1_data*)(*iter)->data;

                os << nd->x << ' ' << nd->y << ' ' << nd->m + 1 << ' ' << nd->vval() << endl;
            }
            os.close();

            cout << " done" << endl;
            delete res;
        }
    }
}

void display_layer_info(const config_dictionary& cfg, const char* fname)
{
    int layer;
    layer1_result* res;

    cfg.get_value(layer, "layer", true);
    
    read_layer1_result(res, fname);
    if (res == nullptr) {
        cout << "Can not read '" << fname << "'." << endl;
        return;
    }
    if (layer < 0 || layer > res->max_layer_index()) {
        cout << "Layer '" << layer << "' does not exist." << endl;
        return;
    }

    int allnodes = 0;
    int topnodes = 0;
    int duplicates = 0;

    vector<node*>& snodes = res->shape_nodes[layer];

    for (vector<node*>::iterator iter = snodes.begin(); iter != snodes.end(); ++iter) {
        node* n = *iter;
        set<int> types;

        ++topnodes;

        do {
            layer1_data* nd = (layer1_data*)n->data;
            
            if (types.find(nd->m) != types.end()) ++duplicates;
            else types.insert(nd->m);
            n = nd->next;
            ++allnodes;
        } while (n != nullptr);
    }

    cout << "Statistics for '" << fname << "':" << endl << endl;
    cout << "All nodes: " << allnodes << endl;
    cout << "Top nodes: " << topnodes << endl;
    cout << "Duplicate types (should be 0!): " << duplicates << endl;

    cout << endl;
    delete res;
}

void save_pca(const config_dictionary& cfg, const char* fname)
{
    int layer;
    string outdir = cfg.get_value_string("out_dir", "");
    part_lib* library;

    end_dir(outdir);
    cfg.get_value(layer, "layer", true);
    
    read_library(fname, library);

    if (library == nullptr) {
        cout << "Can not open library." << endl;
        return;
    }
    if (layer < 0 || layer > library->max_layer_index()) {
        cout << "Layer " << layer << " does not exist in the library." << endl;
    } else {
        for (int i = 0; i < (int)library->parts[layer].size(); ++i) {
            node* p = library->parts[layer][i];
            vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);

            if (pd != nullptr) {
                ofstream os((outdir + "part_" + fill_left(pd->type, 3) + ".csv").c_str());
                
                cv_write_csv<double>(os, pd->pcad.mean);
                cv_write_csv<double>(os, pd->pcad.eigenvectors);
                cv_write_csv<double>(os, pd->pcad.eigenvalues);

                os.close();
            }
        }
    }

    delete library;
}

void save_csv(const config_dictionary& cfg, const char* fname)
{
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");

    string srcdir = cfg.get_value_string("src_dir", "");
    string outdir = cfg.get_value_string("out_dir", "");
    int layer;

    cfg.get_value(layer, "layer", true);

    list<string> files; 

    end_dir(srcdir);
    end_dir(outdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, srcdir))
        list_directory(files, srcdir + fname);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *fiter);
        if (res == nullptr) {
            cout << "Can not load '" << *fiter << "'" << endl;
        } else {
            cout << "Processing " << *fiter << "...";

            ofstream os(change_extension(outdir + *fiter, ".txt").c_str());

            if (!res->grid(0)) res->init_grid(0);
            if (layer >= 0 && layer <= res->max_layer_index()) {
                vector<node*>& s_nodes = res->shape_nodes[layer];
                double factor = 1.0; //(double)res->original_width/(res->x_size(0) - 2*res->border);

                for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
                    node* n = *iter;
                    layer1_data* nd = (layer1_data*)n->data;
                    set<node*> nset;
                    set<ipoint2> pset;
                    
                    os << nd->m << ',';
                    os << nd->r.get_response(R_RESPONSE, 0.0) << ',';
                    os << nd->r.get_response(G_RESPONSE, 0.0) << ',';
                    os << nd->r.get_response(RR_RESPONSE, 0.0) << ',';
                    os << nd->r.get_response(S_RESPONSE, 0.0) << ',';
                    os << nd->r.get_response(X_RESPONSE, 0.0) << ',';

                    res->recurse_and_link(n, toprev, to0, nset);
                    node_set_to_point_set(pset, nset.begin(), nset.end());
                    os << pset.size() << ',';
                    for (set<ipoint2>::iterator piter = pset.begin(); piter != pset.end(); ++piter) {
                        ipoint2 p = *piter;
                        node* nn = res->node_at(p.x, p.y, 0);
                        layer1_data* nnd = (layer1_data*)nn->data;

                        if (piter != pset.begin()) os << ',';
                        os << nnd->m << ',' << nnd->x << ',' << nnd->y;
                    }
                    os << '\n';
                }
            }

            os.close();
            delete res;
            cout << endl;
        }
    }
}

void display_thresholds(const config_dictionary& cfg, const char* fname)
{
    unique_ptr<part_lib> library(read_library(fname));

    if (library == nullptr) {
        cout << "Can not open library '" << fname << "'" << endl;
        return;
    }

    int layer = cfg.get_value_int("layer", 5);
    int thresh = response_from_string(cfg.get_value_string("response", "S_RESPONSE"));

    if (layer < 0 || layer > library->max_layer_index()) {
        cout << "Layer " << layer << " does not exist!" << endl;
        return;
    }
    
    for (int i = 0; i < library->parts[layer].size(); ++i) {
        cout << "Part #" << i << ", threshold: " << library->get_thresh(thresh, layer, i, -1)  << endl;
    }

}

void save_yml(const config_dictionary& cfg, const char* fname)
{
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");

    string srcdir = cfg.get_value_string("src_dir", "");
    string outdir = cfg.get_value_string("out_dir", "");
    int layer;

    cfg.get_value(layer, "layer", true);

    list<string> files; 

    end_dir(srcdir);
    end_dir(outdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, srcdir))
        list_directory(files, srcdir + fname);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *fiter);
        if (res == nullptr) {
            cout << "Can not load '" << *fiter << "'" << endl;
        } else {
            cout << "Processing " << *fiter << "...";

            cv::FileStorage os(change_extension(outdir + *fiter, ".yml"), cv::FileStorage::WRITE);

            os << "segments" << "[";
            if (!res->grid(0)) res->init_grid(0);
            if (layer >= 0 && layer <= res->max_layer_index()) {
                vector<node*>& s_nodes = res->shape_nodes[layer];
                double factor = (double)res->original_width/(res->x_size(0) - 2*res->border);

                for (vector<node*>::iterator iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
                    node* n = *iter;
                    layer1_data* nd = (layer1_data*)n->data;
                    set<node*> nset;
                    set<ipoint2> pset;
                    
                    os << "{";
                    os << "type" << nd->m;
                    os << "contrast" << nd->r.get_response(R_RESPONSE, 0.0);
                    os << "ratio" << nd->r.get_response(RR_RESPONSE, 0.0);
                    os << "precision" << nd->r.get_response(G_RESPONSE, 0.0);
                    //os << nd->r.get_response(S_RESPONSE, 0.0) << ' ';
                    //os << nd->r.get_response(X_RESPONSE, 0.0) << ' ';

                    res->recurse_and_link(n, toprev, to0, nset);
                    node_set_to_point_set(pset, nset.begin(), nset.end());
                    os << "border" << "[:";
                    for (set<ipoint2>::iterator piter = pset.begin(); piter != pset.end(); ++piter) {
                        os << cv::Point((int)(factor*(piter->x - res->border)), 
                            (int)(factor*(piter->y - res->border)));
                    }
                    os << "]";
                    os << "}";
                }
            }

            delete res;
            cout << endl;
        }
    }
}

void save_CSV(const config_dictionary& cfg, const char* fname)
{
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");

    string srcdir = cfg.get_value_string("src_dir", "");
    string outdir = cfg.get_value_string("out_dir", "");
    bool ts = cfg.get_value_bool("translate_and_scale", false);
    int inhibitr = cfg.get_value_int("inhibition_radius", 0);
    int layer;

    cfg.get_value(layer, "layer", true);

    list<string> files; 

    end_dir(srcdir);
    end_dir(outdir);
    if (!cfg.get_value_bool("from_file", false) || !list_from_file(files, fname, srcdir))
        list_directory(files, srcdir + fname);
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        layer1_result* res;

        read_layer1_result(res, srcdir + *fiter);
        if (res == nullptr) {
            cout << "Can not load '" << *fiter << "'" << endl;
        } else {
            cout << "Processing " << *fiter << "...";

            ofstream os(change_extension(outdir + *fiter, ".csv").c_str());

            if (!res->grid(0)) res->init_grid(0);
            if (layer >= 0 && layer <= res->max_layer_index()) {
                double factor = (double)res->original_width/(res->x_size(0) - 2*res->border);

                for (auto iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
                    node* n = *iter;
                    layer1_data* nd = (layer1_data*)n->data;

                    if (nd->z != layer) 
                        continue;

                    set<node*> nset;
                    set<ipoint2> pset;
                    
                    os << nd->m << ',';
                    os << nd->r.get_response(R_RESPONSE, 0.0) << ',';
                    os << nd->r.get_response(RR_RESPONSE, 0.0) << ',';
                    os << nd->r.get_response(G_RESPONSE, 0.0);
                    //os << nd->r.get_response(S_RESPONSE, 0.0) << ' ';
                    //os << nd->r.get_response(X_RESPONSE, 0.0) << ' ';

                    res->recurse_and_link(n, toprev, to0, nset);
                    node_set_to_point_set(pset, nset.begin(), nset.end());
                    if (!ts) {
                        for (set<ipoint2>::iterator piter = pset.begin(); piter != pset.end(); ++piter) {
                            os << ',' << (int)(factor*(piter->x - res->border)) << ',' << 
                                (int)(factor*(piter->y - res->border));
                        }
                    } else {
                        vector<ipoint2> ipv(pset.begin(), pset.end());

                        if (inhibitr > 1) ipv = inhibit_point_set(ipv, inhibitr);

                        vector<dpoint2> dpv = cast_vector<dpoint2, ipoint2>(ipv);
                        
                        translate_and_scale(dpv);
                        for (auto piter = dpv.begin(); piter != dpv.end(); ++piter) {
                            os << ',' << piter->x << ',' << piter->y;
                        }
                    }
                    os << '\n';
                }
            }

            delete res;
            cout << endl;
        }
    }
}

void back_statistics(const config_dictionary& cfg, const char* fname)
{
    cout << "Doing library_statistics" << endl;

    part_lib* library;
    int backname = atom("lyrCenterBack");
    int objectlayer = INT_MAX - 1;

    read_library(fname, library);
    if (library == nullptr) {
        cout << "Can not open library '" << fname << "'" << endl;
        return;
    }

	// Basic library statistics
	//
	cout << "Contractions vector: ";
	for (int i = 0; i < library->layer_count + 1; ++i) {
		if (i > 0) cout << ", ";
		cout << library->contractions[i];
	}
	cout << endl;

    // Back statistics 
    //
    for (int i = 0; i < library->layer_count; ++i) {
        vector<node*>& parts = library->parts[i];
        int max = 0, min = INT_MAX;
        double avg = 0.0;
        
        for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
            node* n = *iter;
            int count = 0;

            foreach_neighbor(n, backname, niter) {
                ++count;
            }
            if (count > max) max = count;
            if (count < min) min = count;
            avg += count;
        }
        avg /= parts.size();
        if (parts.size() > 0 && parts[0]->is_attr_set(OBJ_PART_ATTR)) 
            objectlayer = i;

        cout << "Layer " << i << ": size = " << parts.size() << "; min up = " << min;
        cout << "; max up = " << max << "; avg up = " << avg << endl;

    }

    // Category layer statistics
    //
    if (objectlayer + 1 <= library->max_layer_index()) {
        vector<node*>& parts = library->parts[objectlayer + 1];
        int prevname = atom("lyrPrev");
        int count = 0;

        cout << "Layer " << objectlayer + 1 << " is class layer: number of classes = " << parts.size() << endl;
        for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
            node* p = *iter;
            //set<vector<ipoint2>> posset;

            //foreach_neighbor(p, prevname, oiter) {
            //    vector<ipoint2> pos;
            //    irectangle2 box;

            //    box = library->subpart_positions(pos, neighbor_node(oiter));
            //    for (vector<ipoint2>::iterator i = pos.begin(); i != pos.end(); ++i) *i -= box.ll;
            //    sort(pos.begin(), pos.end());
            //    posset.insert(pos);
            //}
            cout << "  Class #" << ++count << " (" << ((cpart_data*)(p->data))->name << "): number of objects = " << p->count_neighbors(prevname) /*posset.size()*/ << endl;
        }
    }

    delete library;
}

void save_sim_matrix(const config_dictionary& cfg, const char* fname, const char* outname)
{
    int ename = atom("lyrSimVal");

    unique_ptr<part_lib> library(read_library(fname));
    int layer;
    double thresh;
    vector<int> subset;

    if (library == nullptr) 
        return;
    cfg.get_value(layer, "layer", true);
    cfg.get_value(subset, "subset");
    thresh = cfg.get_value_double("threshold", 0.0);

    if (layer < 0 || layer > library->max_layer_index())
        return;
    if (subset.empty()) 
        subset = vector_range(0, library->layer_size(layer) - 1);

    if (cfg.get_value_bool("save_image", false)) {
        img image = library->lib_image(layer + 1, 0, -1, false, false, false, true);
        int d1 = image.height;
        img result(d1*(subset.size() + 1), d1*(subset.size() + 1), 0.0, false);

        for (int si = 0; si < subset.size(); ++si) {
            int di = (si + 1)*d1;
            int dip1 = di + d1;
            int i = subset[si];

            result.blt(image.extract(i*d1, 0, (i + 1)*d1, d1), di, 0);
            result.blt(image.extract(i*d1, 0, (i + 1)*d1, d1), 0, di);
        }
        for (int si = 0; si < subset.size(); ++si) {
            int di = (si + 1)*d1;
            int dip1 = di + d1;
            int i = subset[si];
            node* p = library->parts[layer][i];

            for (int sj = 0; sj < subset.size(); ++sj) {
                int dj = (sj + 1)*d1;
                int djp1 = dj + d1;
                int j = subset[sj];
                node* q = library->parts[layer][j];
                edge_data_t<double>* ed = (edge_data_t<double>*)p->get_edge_data(q, ename);

                if (ed != nullptr) {
                    if (ed->data < thresh) {
                        double val = max(1.0 - ed->data, 0.0);
    
                        IMG_COLOR cl;

                        cl.blue = cl.green = cl.red = COL_FROM_REAL(val);
                        result.set_region(di, dip1, dj, djp1, *((double*)(&cl)));
                    }
                }
            }
        }

        result.save_normalized(change_extension(outname, ".png"));
    } 
    if (cfg.get_value_bool("save_matrix", true)) {
        ofstream os(change_extension(outname, ".m").c_str());

        for (int si = 0; si < subset.size(); ++si) {
            int i = subset[si];
            node* p = library->parts[layer][i];

            for (int sj = 0; sj < subset.size(); ++sj) {
                int j = subset[sj];
                node* q = library->parts[layer][j];
                edge_data_t<double>* ed = (edge_data_t<double>*)p->get_edge_data(q, ename);

                if (sj > 0) os << ',';
                os << ((ed == nullptr) ? 0 : ed->data);
            }
            os << endl;
        }
        os.close();
    }
}

void display_similarity(const config_dictionary& cfg, const char* fname)
{
    int layer;
    part_lib* library;
    int edge = atom("lyrSimRoot");

    cfg.get_value(layer, "layer", true);
    read_library(fname, library);
    if (library == nullptr) {
        cout << "Can not open library." << endl;
        return;
    }
    if (layer <= 0 || layer > library->max_layer_index()) {
        cout << "Layer does not exist." << endl;
    } else {
        vector<node*>& parts = library->parts[layer];

        for (vector<node*>::iterator iter = parts.begin(); iter != parts.end(); ++iter) {
            node* p = *iter;
            bool root = false;
            
            foreach_neighbor(p, edge, niter) {
                if (neighbor_node(niter) == p) {
                    root = true;
                    break;
                }
            }
            if (root) {
                lib_data* pd = (lib_data*)p->data;

                cout << "root: " << pd->type << " ~: ";
                for (vector<node*>::iterator jter = parts.begin(); jter != parts.end(); ++jter) {
                    node* pp = *jter;
                    lib_data* ppd = (lib_data*)pp->data;

                    foreach_neighbor(pp, edge, niter) {
                        if (neighbor_node(niter) == p) {
                            cout << ppd->type << " (";
                            cout << ((edge_data_t<double>*)neighbor_edge_data(niter))->data << ")  ";
                            break;
                        }
                    }
                }
                cout << endl;
            }
        }
    }

    delete library;

}

void read_bins(vector<irectangle2>& bins, const string& fname)
{
    ifstream is(fname.c_str());

    bins.clear();
    while (is.good()) {
        int a, b, c, d;

        is >> a >> b >> c >> d;

        bins.push_back(irectangle2(a, b, a + c, b + d));
    }
}

const int ODRAW_DEFAULT = 0;
const int ODRAW_CHULL = 1;


/*

1.  h = 0 
2.  for every point ai of A,
      2.1  shortest = Inf ;
      2.2  for every point bj of B
                    dij = d (ai , bj )
                    if dij < shortest then
                              shortest = dij
      2.3  if shortest > h then 
                    h = shortest 

*/

double node_set_distance_h(const vector<node*>& v1, const vector<node*>& v2)
{
    double result = 0.0;

    for (auto iter1 = v1.begin(); iter1 != v1.end(); ++iter1) {
        double mind = numeric_limits<double>::max();

        for (auto iter2 = v2.begin(); iter2 != v2.end(); ++iter2) {
            double dist = sqrt((double)(node_coordinates(*iter1).distance2(node_coordinates(*iter2))));
            if (dist < mind) mind = dist;
        }
        if (mind > result) result = mind;
    }
    return result;
}

double node_set_distance(const vector<node*>& v1, const vector<node*>& v2)
{
    return max(node_set_distance_h(v1, v2), node_set_distance_h(v2, v1));
}


/*
bool texture_box(layer1_result* res, int r, double emptythresh)
{
    int xdim = res->x_size(0);
    int ydim = res->y_size(0);
    int allbox = 0;
    int emptybox = 0;

    if (!res->grid(0)) 
        res->init_grid(0);
    for (int x = r; x + r < xdim; x += r) {
        for (int y = r; y + r < ydim; y += r) {
            int all = 0, count = 0;

            ++allbox;
            for (int i = x - r; i < x + r; ++i) {
                for (int j = y - r; j < y + r; ++j) {
                    if (i >= 0 && j >= 0 && i < xdim && j < ydim) {
                        ++all;
                        if (res->node_at(i, j, 0) != nullptr) 
                            ++count;
                    }
                }
            }
            if (count < emptythresh*all)
                ++emptybox;
        }
    }
    
    return emptybox < 0.1*allbox;
}
*/

bool is_texture(layer1_result* res, double contrast_thresh) 
{
    if (res->nodes.empty())
        return true;
    
    node* n = res->nodes.front();
    layer1_data* nd = (layer1_data*)n->data;

    return nd->vval() < contrast_thresh;
}

void write_yml_detections(const vector<vector<cv::Point> >& detpts, const string& fname)
{
    cv::FileStorage os(fname, cv::FileStorage::WRITE);

    os << "detections" << "[";
    for (int i = 0; i < (int)detpts.size(); ++i) {
        os << "{";
        os << "index" << i;
        os << "points" << "[:";
        for (auto piter = detpts[i].begin(); piter != detpts[i].end(); ++piter) 
            os << *piter;
        os << "]";
        os << "}";
    }
}



struct dhistogram_data {
    int p;              // +-1 true/false detection
    int model;          // type
    double r;           // entry 14 from visview export (r-response)
    double rr;          // additional (rr-response)
    double s;           // additional (s-response)
    vector<double> h;   // histogram
    irectangle2 bbox;   // bounding box
};

void save_dhistogram(const config_dictionary& cfg, const string& fpatt)
{
    int edgename = atom("toPrevLayer");
    int to0 = atom("toLayer0");


    string dir = cfg.get_value_string("src_dir", "");
    string outdir = cfg.get_value_string("out_dir", "");
    list<string> files;
    int layer;
    string catname;
    int catindex = -1;
    map<string, int> catmap;
    string libname;
    vector<vector<int> > typemap;
    int hsize = 0;
	string pattern;

    cfg.get_value(layer, "layer", true);
    cfg.get_value(catname, "category", true);
    cfg.get_value(libname, "library", true);

    unique_ptr<part_lib> library(read_library(libname));
    
	if (library == nullptr) {
		cout << "Can not open library (library is needed to extract category names)" << endl;
		return;
    }
	if (layer < 0 || layer > library->max_layer_index()) {
		cout << "Layer #" << layer << " does not exist." << endl;
		return;
	}
    for (vector<node*>::iterator iter = library->parts[layer].begin(); iter != library->parts[layer].end(); ++iter) {
		cpart_data* pd = dynamic_cast<cpart_data*>((*iter)->data);

		if (pd) {
			catmap.insert(pair<string, int>(pd->name, pd->type));
            if (pd->name == catname) 
                catindex = pd->type;
            
        }
	}
	if (catmap.empty()) {
		cout << "No category names can be extracted from layer " << layer << " in the library!" << endl;
		return;
	}
    if (catindex < 0) {
        cout << "Category " << catname << " is not present in the library!" << endl;
        return;
    }

    hsize = 0;
    typemap.resize(layer + 1);
    for (int l = 0; l <= layer; ++l) {
        typemap[l].resize(library->layer_size(l));
        for (int i = 0; i < typemap[l].size(); ++i) {
            typemap[l][i] = hsize++;
        }
    }

    cout << "Saving discrimination histograms from layer " << layer << endl;
    end_dir(dir);
    end_dir(outdir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, fpatt, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + fpatt);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cout << "Processing " << *file;

        unique_ptr<layer1_result> res(read_layer1_result(dir + *file));
    
        if (res != nullptr) {
            list<pair<irectangle2, string> > gtruths;
            list<dhistogram_data> histograms;

            read_groundtruth(gtruths, dir + *file, catname, ".groundtruth", false);
            for (auto niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
                node* n = *niter;
                layer1_data* nd = (layer1_data*)n->data;

                if (nd->z == layer) {
                    node* nn = n->get_neighbor(edgename);
                    layer1_data* nnd = (layer1_data*)nn->data;
                    node* p = library->parts[nnd->z][nnd->m];

                    if (!p->is_attr_set(NODE_DELETED_ATTR)) {
                        irectangle2 box;
                        bool positive = false;

                        histograms.push_back(dhistogram_data());

                        dhistogram_data& hd = histograms.back();

                        make_dhistogram(hd.h, box, typemap, hsize, res.get(), n);
                        box.grow(1.2);
                        for (auto gtiter = gtruths.begin(); gtiter != gtruths.end(); ++gtiter) {
                            if (gtiter->second == catname) {
                                double x = ((double)(gtiter->first.intersection(box).area()))/gtiter->first.union_area(box);

                                if (x >= 0.5) {
                                    positive = true;
                                }
                            }
                        }
                        hd.p = positive ? 1 : -1;
                        hd.r = nnd->r(G_RESPONSE);
                        hd.rr = nnd->r(RR_RESPONSE);
                        hd.s = nnd->r(S_RESPONSE);
                        hd.model = nnd->m;
                    }
                        
                }

            }
            
            string ofname = change_extension(outdir + *file, ".hist");
            ofstream os(ofname, ios_base::binary | ios_base::out);

            if (!histograms.empty()) {
                bin_write_int(os, (int)histograms.size());
                for (auto hliter = histograms.begin(); hliter != histograms.end(); ++hliter) {
                    bin_write_int(os, hliter->p*(hliter->model + 1));
                    bin_write_double(os, hliter->r);
                    bin_write_double(os, hliter->s);
                    //bin_write_double(os, hliter->h.size());
                    //os.write((const char*)(&hliter->h.at(0)), sizeof(double) * hliter->h.size());
                    write_sparse_histogram(os, hliter->h);
                }
            }
            os.close();
        }
        cout << endl;
    }


}

void save_path_dhistogram(const config_dictionary& cfg, const string& fpatt)
{
    typedef map<vector<int>, int> map_t;

    int edgename = atom("toPrevLayer");
    int to0 = atom("toLayer0");

    string dir = cfg.get_value_string("src_dir", "");
    string outdir = cfg.get_value_string("out_dir", "");
    list<string> files;
    int layer;
	int depth;
    string catname;
    int catindex = -1;
    map<string, int> catmap;
    string libname;
    vector<map_t> typemap;
    int hsize = 0;
	string pattern;

    cfg.get_value(layer, "layer", true);
    cfg.get_value(catname, "category", true);
    cfg.get_value(libname, "library", true);
	depth = cfg.get_value_int("histogram_depth", -1);	// -1 or big value means all layers

    unique_ptr<part_lib> library(read_library(libname));
    
	if (library == nullptr) {
		cout << "Can not open library." << endl;
		return;
    }
	if (layer < 0 || layer > library->max_layer_index()) {
		cout << "Layer #" << layer << " does not exist." << endl;
		return;
	}
    for (vector<node*>::iterator iter = library->parts[layer].begin(); iter != library->parts[layer].end(); ++iter) {
		cpart_data* pd = dynamic_cast<cpart_data*>((*iter)->data);

		if (pd) {
			catmap.insert(pair<string, int>(pd->name, pd->type));
            if (pd->name == catname) 
                catindex = pd->type;
            
        }
	}
	if (catmap.empty()) {
		cout << "No category names can be extracted from layer " << layer << " in the library!" << endl;
		return;
	}
    if (catindex < 0) {
        cout << "Category " << catname << " is not present in the library!" << endl;
        return;
    }

    typemap.resize(library->layer_size(layer - 1));
    for (int pi = 0; pi < library->layer_size(layer - 1); ++pi) {
        part_path_map(library->parts[layer - 1][pi], typemap[pi], 0);
    }

	update_path_permutations(library.get());

    cout << "Saving discrimination histograms from layer " << layer << endl;
    end_dir(dir);
    end_dir(outdir);
    if (!cfg.get_value_bool("from_file", false) || !file_list_from_file(files, fpatt, dir, cfg.get_value_string("pattern", "")))
        list_directory(files, dir + fpatt);
    for (list<string>::iterator file = files.begin(); file != files.end(); ++file) {
        cout << "Processing " << *file;

        unique_ptr<layer1_result> res(read_layer1_result(dir + *file));
    
        if (res != nullptr) {
            list<pair<irectangle2, string> > gtruths;
            list<dhistogram_data> histograms;

            read_groundtruth(gtruths, dir + *file, catname, ".groundtruth", false);
            for (auto niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
                node* n = *niter;
                layer1_data* nd = (layer1_data*)n->data;

                if (nd->z == layer) {
                    node* nn = n->get_neighbor(edgename);
                    layer1_data* nnd = (layer1_data*)nn->data;
                    node* p = library->parts[nnd->z][nnd->m];

                    if (!p->is_attr_set(NODE_DELETED_ATTR)) {
                        irectangle2 box;
                        bool positive = false;

                        histograms.push_back(dhistogram_data());

                        dhistogram_data& hd = histograms.back();

                        //make_path_dhistogram(vector<double>& hist, irectangle2& bbox, const map<vector<int>, int>& typemap, layer1_result* res, node* n)
						make_path_dhistogram(hd.h, box, typemap[nnd->m], library.get(), res.get(), nn, depth);
                        box.grow(1.2);
                        for (auto gtiter = gtruths.begin(); gtiter != gtruths.end(); ++gtiter) {
                            if (gtiter->second == catname) {
                                double x = ((double)(gtiter->first.intersection(box).area()))/gtiter->first.union_area(box);

                                if (x >= 0.5) {
                                    positive = true;
                                }
                            }
                        }
                        hd.p = positive ? 1 : -1;
                        hd.r = nnd->r(G_RESPONSE);
                        hd.rr = nnd->r(RR_RESPONSE);
                        hd.s = nnd->r(S_RESPONSE);
                        hd.model = nnd->m;
                    }
                        
                }

            }
            
            string ofname = change_extension(outdir + *file, ".hist");
            ofstream os(ofname, ios_base::binary | ios_base::out);

            if (!histograms.empty()) {
                bin_write_int(os, (int)histograms.size());
                for (auto hliter = histograms.begin(); hliter != histograms.end(); ++hliter) {
                    bin_write_int(os, hliter->p*(hliter->model + 1));
                    bin_write_double(os, hliter->r);
                    bin_write_double(os, hliter->s);
                    //bin_write_double(os, hliter->h.size());
                    //os.write((const char*)(&hliter->h.at(0)), sizeof(double) * hliter->h.size());
                    write_sparse_histogram(os, hliter->h);
                }
            }
            os.close();
        }
        cout << endl;
    }


}

void init_display(const char* file, const char* cfg_file, const char* outname, const char* params)
{

    cout << "laydisplay (" __DATE__ " / " __TIME__ ")" << endl;
    try {
        string cff = cfg_file;
        string mode;
    
        if (cff == "") cff = "c:\\programs\\laydisplay.cfg";
        config_dictionary cfg(cff);
        
        cfg.from_string(params);
        mode = cfg.get_value_string("mode", "r");
        transform(mode.begin(), mode.end(), mode.begin(), TOUPPER);
        if (mode == "M") save_mathematica2(cfg, file, outname);
        else if (mode == "V") save_visview(cfg, file, outname);
        else if (mode == "CSV") save_CSV(cfg, file);
        else if (mode == "G") save_graphs(cfg, file, outname);
        else if (mode == "L") save_library(cfg, file, outname);
        else if (mode == "LYML") save_library_yml(cfg, file, outname);
        else if (mode == "LSC") save_library_sc(cfg, file, outname);
        else if (mode == "LSCM") save_library_sc_mma(cfg, file, outname);
        else if (mode == "LMEAN") save_library_mean(cfg, file, outname);
        else if (mode == "GL") save_library_graph(cfg, file, outname);
        else if (mode == "DELLIB") remove_library_layers(cfg, file, outname);
        else if (mode == "LM") save_library_mma(cfg, file, outname);
        else if (mode == "LMB") save_library_mma_back(cfg, file, outname);
        else if (mode == "STR") print_parts(cfg, file, outname);
        else if (mode == "P") save_part_statistics(cfg, file, outname);
        else if (mode == "C") save_part_covering_statistics(cfg, file, outname);
        else if (mode == "BOX") compare_boxes(cfg, file, outname);
        else if (mode == "DROPIMG") drop_library_images(cfg, file, outname);
        else if (mode == "~") save_similarity_matrix(cfg, file, outname);
        else if (mode == "T") test_export(cfg, file, outname);
        else if (mode == "AP") save_average_parts(cfg, file, outname);
        else if (mode == "E") save_edge_names(cfg, file, outname);
        else if (mode == "M2") save_mathematica(cfg, file, outname);
        else if (mode == "M3") save_mathematica3(cfg, file, outname);
        else if (mode == "LI") display_library_info(cfg, file);
        else if (mode == "OBJECTS") save_objects(cfg, file, outname);
        else if (mode == "NI") save_node_info(cfg, file, outname);
        else if (mode == "D") dilute_layer(cfg, file);
        else if (mode == "TX") track_layer(cfg, file, outname);
        else if (mode == "SAVEBOX") save_bounding_box(cfg, file);
        else if (mode == "LL") save_library_layer(cfg, file, outname);
        else if (mode == "SELECT") select_files(cfg, file, outname);
        else if (mode == "C2F") display_c2f(cfg, file, outname);
        else if (mode == "C2F_GRAPH") display_c2f_graph(cfg, file, outname);
		else if (mode == "CLASS") save_classification(cfg, file, outname);
        else if (mode == "CATPARTS") save_category_parts(cfg, file, outname);
		else if (mode == "VIDEO") save_video(cfg, file, outname);
        else if (mode == "PR") save_part_reconstruction(cfg, file, outname);
        else if (mode == "TXT") remove_texture(cfg, file, outname);
        else if (mode == "ELLIPSES") save_ellipses(cfg, file, outname);
        else if (mode == "SANJA") save_sanja(cfg, file, outname);
        else if (mode == "SIM") display_similarity(cfg, file);
        else if (mode == "LAYERINFO") display_layer_info(cfg, file);
        else if (mode == "BACKSTAT") back_statistics(cfg, file);
        else if (mode == "PCA") save_pca(cfg, file);
        else if (mode == "CSV") save_csv(cfg, file);
        else if (mode == "YML") save_yml(cfg, file);
        else if (mode == "THRESHOLDS") display_thresholds(cfg, file);
		else if (mode == "FEAT") save_features(cfg, file, outname);
        else if (mode == "SAMPLE") sample_part(cfg, file, outname);
        else if (mode == "TEST") test(cfg, file, outname);
        else if (mode == "ECCV") eccv_export(cfg, file, outname);
        else if (mode == "MERGE") merge_scales(cfg, file, outname);
        else if (mode == "SIM_MATRIX") save_sim_matrix(cfg, file, outname);
        else if (mode == "SUPPORT") save_support(cfg, file, outname);
        else if (mode == "DHISTOGRAM") save_dhistogram(cfg, file); // Histogram for "Matej discrimination learning"
        else if (mode == "PATH_DHISTOGRAM") save_path_dhistogram(cfg, file); // Histogram for "Matej discrimination learning"
//        else if (mode == "FEAT2") save_features2(cfg, file, outname);
        else if (mode == "DOMENGT") domen_gt(cfg, file, outname);
        else save_normal(cfg, file, outname);

    } catch (const libhop_exception& e) {
		cout << e.what() << endl;
    } catch (const exception& e) {
		cout << "General exception '" << typeid(e).name() << "' with message: '" << e.what() << "'" << endl;
	} 

}

void print_usage() 
{
    cout << "hopdisplay (" __DATE__ " / " __TIME__ ")" << endl;
    cout << "Usage: hopdisplay fname [config file] [out name]" << endl;
    cout << "config file options:" << endl;
    cout << "   mode            = n|i|r|v|s|l|m|p|c" << endl;
    cout << "       (normal, inhibited, reconstruction (default), visview, SVM, library, " << endl;
    cout << "           Mathematica, part statistics, part covering statistics)" << endl;
    cout << "   layer           = index of layer to display (-1 means the last one (default))" << endl;
    cout << "   parts           = list of parts to display (not specified means all)" << endl;
    cout << "   part_colors     = list of colors for parts" << endl;
    cout << "       (list of triples form [0,..,255]^3 specifying rgb color)" << endl;
    cout << "   default_color   = default color for parts" << endl;
    cout << "   uncovered_color = color for uncovered parts (if not specified, uses" << endl;
    cout << "       color for parts)" << endl;
    cout << "   out_dir         = output directory (visview)" << endl;
    cout << "   scales          = list of scales to include to SVM export" << endl;
    cout << " " << endl;
}


int main(int argc, char* argv[])
{

    ///_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
    ///_CrtSetBreakAlloc(326330);
    init_atoms();
    init_streaming();   
    random_seed((int)time(nullptr));

    switch (argc) {
        case 2: init_display(argv[1], "", "", ""); break;
        case 3: init_display(argv[1], argv[2], "", ""); break;
        case 4: init_display(argv[1], argv[2], argv[3], ""); break;
        case 5: init_display(argv[1], argv[2], argv[3], argv[4]); break;
        default: print_usage();
    }

	return 0;
}

