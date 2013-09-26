// test.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <tchar.h>
#include <windows.h>
#include <psapi.h>

#include "layers/layer_1.h"


//#include "stdafx.h"
// #include <iostream>
// #include <ctime>
// #include "../img/img.h"
#include "layers/optimization.h"
#include "layers/initialization.h"


#include <ctime>

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <sys/stat.h>
#include <crtdbg.h>
#include <hash_map>
#include <queue>
#include "utils/ocv.h"
#include "utils/convert.h"
#include <cstdarg>
#include "utils/utils.h"
#include "boost/filesystem.hpp"
#include "matching.h"
#include "graphs/graph_utils.h"

//#define _SAVEMAPS
#define PRINT_INFO(_info) cerr << _info << endl;
#define PRINT_DOT() cerr << '.';
//#define PRINT_INFO(_info) 
//#define PRINT_DOT()
//#define LOG_SP(_info) log_stream << _info << endl; log_stream.flush();
#define LOG_SP(_info) cerr << _info << endl; 
//#define LOG_SP(_info)


using namespace std;

ostream& print_matrix(ostream& os, const cv::Mat m)
{
    for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.cols; ++j) {
            os << m.at<double>(i, j) << ' ';
        }
        cout << endl;
    }
    return os;
}


/*
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

    /////////
        //cout << endl;
        //for (int i = 0; i < pcount+3; ++i) {
        //    for (int j = 0; j < pcount+3; ++j) {
        //        cout << L.at<double>(i, j) << ' ';
        //    }
        //    cout << endl;
        //}
    /////////

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

    return result[0];
}

// "Normalize" vector (translate & scale).
void translate_and_scale(vector<dpoint2>& v)
{
    if (v.empty()) return;

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
}

// Rotate vector for angle a about origin.
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

*/



template<class T> void write_CSV(ostream& os, const cv::Mat& m)
{
    for (int r = 0; r < m.rows; ++r) {
        for (int c = 0; c < m.cols; ++c) {
            if (c > 0) os << ',';
            os << m.at<T>(r, c);
        }
        os << '\n';
    }
}





void test1()
{
    //vector<dpoint2> pts1 = make_vector<dpoint2>(4, dpoint2(1, 0), dpoint2(4, 0), dpoint2(5, 0), dpoint2(6, 0));
    //vector<dpoint2> pts2 = make_vector<dpoint2>(4, dpoint2(1, 0.5), dpoint2(4, 0.6), dpoint2(5, 0.4), dpoint2(6, 0.4));
    //vector<dpoint2> pts2 = make_vector<dpoint2>(4, dpoint2(1, 0.5), dpoint2(4, 0.6), dpoint2(5, -0.4), dpoint2(6, 0.4));
    //vector<dpoint2> pts1 = make_vector<dpoint2>(4, dpoint2(-2, 0), dpoint2(-1, 0), dpoint2(1, 0), dpoint2(2, 0));
    //vector<dpoint2> pts2 = make_vector<dpoint2>(4, dpoint2(-2, 1), dpoint2(-1, 0.5), dpoint2(1, -0.5), dpoint2(2, -1));

    //cout << "Bending energy is " << TPS_energy(pts1, pts2) << endl;
    //cout << "Optimal pts1 -> pts2 rotation is: " << optimal_rotation(pts1, pts2) << endl;
    //rotate(pts1, optimal_rotation(pts1, pts2));
    //cout << "Rotated points: " << pts2 << endl;

    /*
    string match_dir = "D:\\work\\vshapes\\data\\matchings\\match\\";
    string src_dir = "D:\\work\\vshapes\\data\\matchings\\layerx\\";
    list<string> files;

    list_directory(files, src_dir + "figA*_0.ly4");
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        boost::filesystem::path ipath(*fiter);
        layer1_result* res1;

        read_layer1_result(res1, *fiter);
        if (res1 == NULL)
            continue;

        for (list<string>::iterator fjter = files.begin(); fjter != files.end(); ++fjter) {
            if (*fiter == *fjter) 
                continue;

            boost::filesystem::path jpath(*fjter);
            layer1_result* res2;

            read_layer1_result(res2, *fjter);
            if (res2 == NULL) 
                continue;

            string match_name = match_dir + ipath.stem() + "-" + jpath.stem() + ".txt";
            


            
        }
            
    }*/

    config_dictionary cfg("c:\\work\\data\\multiclass\\vshapes\\test.cfg");
    sim_learning slearner(cfg);
    part_lib* library;
    vector<int> perm;
    vector<dpoint2> src, dest;
    layer1_result* res1;
    layer1_result* res2;

    read_library("c:\\work\\data\\multiclass\\vshapes\\lib4.plb", library);
    for (int i = 1; i <= 11/*11*/; ++i)
        for (int j = 1; j <= 11/*11*/; ++j) {
            string iname = string("figA") + i + string("_0");
            string jname = string("figA") + j + string("_0");

            //read_matching(perm, src, dest, string("D:\\work\\vshapes\\data\\matchings\\match\\") + 
//                iname + string("-") + jname + string(".txt"));
            read_layer1_result(res1, string("D:\\work\\vshapes\\data\\matchings\\layerx\\") + iname + string(".ly4"));
            read_layer1_result(res2, string("D:\\work\\vshapes\\data\\matchings\\layerx\\") + jname + string(".ly4"));

            if (res1 != NULL && res2 != NULL && !perm.empty()) {
                cout << i << " <-> " << j << "; length: " << perm.size() << endl;
                slearner.update_similarities(res1, res2, perm, src, dest);
            }

            delete res1;
            delete res2;
    }


    slearner.print();

    graph* g = slearner.similarity_graph(cfg.get_value_double("matching_threshold", 0.33));

 //   pca_test(library, cfg.get_value_int("layer", 3), g);
    library->save("c:\\work\\data\\multiclass\\vshapes\\lib4vs.plb");


    delete library;
    delete g;
    //res2->save_visview("c:\\work\\", "1", "figA2_0", 0, 1);
}

vector<ipoint2> sort_points(const vector<ipoint2>& pts)
{
    return vector<ipoint2>();
}






struct average_operator_int {
    void add_to(int& a, int b) const { a += b; }
    int div_int(int a, int n) const { return a / n; }
};

    
void test2() 
{
    graph* g = new graph();

    node* nodes[5];

    for (int i = 0; i < 5; ++i) 
        nodes[i] = g->add_node(new node_data_t<int>(i));

    g->add_edge(nodes[0], nodes[1], 0, 0);
    g->add_edge(nodes[0], nodes[2], 0, 0);
    g->add_edge(nodes[1], nodes[2], 0, 0);
    g->add_edge(nodes[2], nodes[3], 0, 0);
    //g->add_edge(nodes[2], nodes[4], 0, 0);

    list<set<node*> > cliques;

    get_maximal_cliques(cliques, g);

    for (list<set<node*> >::iterator iter = cliques.begin(); iter != cliques.end(); ++iter) {
        for (set<node*>::iterator niter = iter->begin(); niter != iter->end(); ++niter) {
            node* n = *niter;
            node_data_t<int>* nd = (node_data_t<int>*)n->data;

            cout << nd->data << ' ';
        }
        cout << endl;
    }
    delete g;

    part_lib* library;

    read_library("c:\\work\\data\\multiclass\\vshapes\\lib4.plb", library);

    node* p1 = library->parts[3][0];
    path_map_t pm1, pm2;
    vector<ipoint2> pts1, pts2;
    ofstream os;
    
    os.open("C:\\work\\testgraphpts1.m");
    get_library_geo(pm1, p1);
    for (path_map_t::iterator iter = pm1.begin(); iter != pm1.end(); ++iter) {
        pts1.push_back(iter->second.p);
        os << pts1.back().x << ',' << pts1.back().y << endl;
    }
    g = point_graph(pts1);
    os.close();

    os.open("c:\\work\\tmp\\tree.m");
    for (graph::iter_t iter = g->nodes.begin(); iter != g->nodes.end(); ++iter) {
        node* n = *iter;
        node_data_t<int>* nd = (node_data_t<int>*)n->data;

        forall_neighbors(n, niter) {
            node_data_t<int>* nnd = (node_data_t<int>*)neighbor_node_data(niter);

            os << pts1[nd->data].x << ',' << pts1[nd->data].y << ',' << pts1[nnd->data].x << ',' << pts1[nnd->data].y << '\n';
        }
    }
    os.close();
    point_graph_image(g, pts1).save("c:\\work\\tmp\\tree.png");

    os.open("C:\\work\\testgraphpts2.m");
    get_library_geo(pm2, library->parts[3][30]);
    for (path_map_t::iterator iter = pm2.begin(); iter != pm2.end(); ++iter) {
        pts2.push_back(iter->second.p);
        os << pts2.back().x << ',' << pts2.back().y << endl;
    }
    os.close();

    reverse(pts2.begin(), pts2.end());
    vector<int> perm = point_matching(pts1, pts2);
    os.open("c:\\work\\matching.m");
    for (int i = 0; i < perm.size(); ++i) {
        os << pts1[perm[i]].x << ',' << pts1[perm[i]].y << ',' << pts2[i].x << ',' << pts2[i].y << '\n';
    }
    os.close();

    os.open("C:\\work\\testgraph.m");
    g->write_mma(os);
    os.close();

    cout << pts2 << endl;
    resize_vector(pts2, 20);
    cout << pts2 << endl;

    average_learning<int, int, average_operator_int> avgl;

    avgl.update(1, 2);
    avgl.update(1, 4);
    avgl.update(2, 22);

    map<int, int> avgres;
    
    avgl.get_average_map(avgres);


    delete g;

    delete library;
    
}

void test3() 
{
    part_lib* library;

    read_library("c:\\work\\data\\multiclass\\vshapes\\lib5.plb", library);

    int layer = 4;
    vector<node*>& parts = library->parts[layer];

    for (int i = 0; i < (int)parts.size(); ++i) {
        node* p = parts[i];
        vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);

        if (pd == NULL) continue;

        stringstream s;

        s << "c:\\work\\tmp\\lib5part";
        s.width(3);
        s.fill('0');
        s << i << ".csv";

        ofstream os(s.str().c_str());

        cv_write_csv<double>(os, pd->pcad.mean);
        cv_write_csv<double>(os, pd->pcad.eigenvectors);
        cv_write_csv<double>(os, pd->pcad.eigenvalues);

        os.close();
    }
}

void get_matching_permutation(map<vector<int>, int>& result, const path_map_t& pm, const cv::Mat& target)
{
	typedef map<vector<int>, int> result_t;

	vector<dpoint2> pmpts;
	vector<vector<int> > ppmap;  // path-point index map (path map is map and has no indices!)
	vector<dpoint2> dpts = partition(target);
	
	// Fill point vector
	pmpts.reserve(pm.size());
	ppmap.reserve(pm.size());
    for (path_map_t::const_iterator pmiter = pm.begin(); pmiter != pm.end(); ++pmiter) {
		pmpts.push_back((dpoint2)pmiter->second.p);
		ppmap.push_back(pmiter->first);
	}
	translate_and_scale(pmpts);
	
	// Resize pmpts and ppmap (synchronously) to the size of 'target'
    resize_vector(pmpts, dpts.size());
    resize_vector(ppmap, dpts.size());
	
	// Permute
	vector<int> perm = point_matching(pmpts, dpts);
	
	permute(ppmap, perm);
	
	result.clear();
	for (int i = 0; i < (int)ppmap.size(); ++i) 
		result.insert(result_t::value_type(ppmap[i], i));
}


void update_permutations(part_lib* library, int layer)
{
    int name = atom("lyrSimRoot");
    vector<node*>& parts = library->parts[layer];

    for (vector<node*>::iterator piter = parts.begin(); piter != parts.end(); ++piter) {
        node* p = *piter;
        vs_part_data* pd = dynamic_cast<vs_part_data*>(p->data);
        path_map_t ppm;

        get_library_geo(ppm, p);
        foreach_neighbor(p, name, niter) {
            node* pn = neighbor_node(niter);
            part_data_sim* pned = dynamic_cast<part_data_sim*>(neighbor_edge_data(niter));

            if (pd != NULL && pned != NULL)
                get_matching_permutation(pned->perm, ppm, pd->pcad.mean);
        }
    }
}

void test4()
{
    string src_dir = "D:\\work\\data\\multiclass\\apple-ferrari-sc\\layerx\\";
    double mergenorm = 1.0;
    int layer = 5;

    list<string> files;
    part_lib* library;
    pca_merging pmg(layer, 5);

    read_library("C:\\work\\data\\multiclass\\apple-ferrari-sc\\olib.plb", library);
    
    list_directory(files, src_dir + "traincut_*_0.ly7");
    for (list<string>::iterator fiter = files.begin(); fiter != files.end(); ++fiter) {
        boost::filesystem::path ipath(*fiter);
        layer1_result* res;

        cout << "Processing: " << *fiter << "...";
        read_layer1_result(res, src_dir + *fiter);
        pmg.update_stat(res, response_filter());
        if (res != NULL) 
            delete res;
        cout << endl;
    }

    pmg.update_library(library, mergenorm);
    library->save("C:\\work\\data\\multiclass\\apple-ferrari-sc\\olibm.plb");
    //cerr << "read" << endl;
    //library->save_all_sc("c:\\work\\data\\multiclass\\vshapes\\lib5sc.png", 5); 

    delete library;


    //layer1_result* res; 

    //read_layer1_result(res, "D:\\work\\vshapes\\data\\matchings\\layerx\\figA1_0.ly4");
    //streamed_pointer sp(res);
    //vector<streamed_pointer> ptrs;

    //ptrs.push_back(sp);

    //streamed_pointer sq;

    //sq = ptrs.front();
    //sq = streamed_pointer(res);
    //sq = NULL;

    //
    //delete res;
}

// Mean shift algorithm
// 'x': result
// 'pv': vector of points
// 'init': initial value
// 'radius': neighborhood of current 'x'
// 'c': normal kernel parameter in e^{c t^2}
// 'eps', 'maxsteps': iteration limits
// returns last 'eps' i.e.: |x_n - x_{n-1}|
double mean_shift(dpoint2& x, const vector<dpoint2>& pv, const dpoint2& init, double radius, double c, double eps, int maxsteps)
{
    double radius2 = radius*radius;
    double eps2 = eps*eps;
    double delta2 = numeric_limits<double>::infinity();

    x = init;
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
        if (delta2 < eps2)
            break;
        //cerr << x << endl;
    }
    return sqrt(delta2);
}
/*

void inhibit_layer(layer1_result* res, int z, int response, int maxn, double thresh)
{
    int to_prev = atom("toPrevLayer");
    int to_0 = atom("toLayer0");
    vector<node*> nvec;

    for (list<node*>::iterator iter = res->nodes.begin(); iter != res->nodes.end(); ++iter) {
        if (node_layer(*iter) == z) {
            nvec.push_back(*iter);
            (*iter)->set_attr(HIDDEN_NODE_ATTR); // We hide all nodes by default and "unhide" them below
        }
    }
    sort(nvec.begin(), nvec.end(), response_sort_f(response));

    // 'Inhibit' nodes
    list<irectangle2> boxes;
    int maxi = min<int>(maxn, (int)nvec.size());

    for (int i = 0; i < maxi; ++i) {
        node* n = nvec[i];
        set<node*> nset;
        irectangle2 box;
		double maxt = 0.0;

        res->recurse_and_link(n, to_prev, to_0, nset);
        box = node_set_bounding_rectangle(nset.begin(), nset.end());

        for (list<irectangle2>::iterator riter = boxes.begin(); riter != boxes.end(); ++riter) {
			irectangle2& rbox = *riter;
            double r = (double)(rbox.intersection(box).area())/box.union_area(rbox);

			if (r > maxt) maxt = r;
		}
        if (maxt <= thresh) {
			boxes.push_back(box);
            n->clear_attr(HIDDEN_NODE_ATTR);
        }
    }

    int count = 0;
    for (vector<node*>::iterator viter = nvec.begin(); viter != nvec.end(); ++viter) {
        if ((*viter)->is_attr_set(HIDDEN_NODE_ATTR)) ++count;
    }
    cout << count << " hidden nodes at layer " << z << endl;
}

**/


ipoint2 get_closest_position(layer1_result* res, int layer, const dpoint2& p, int minr, int maxr)
{
    if (layer < 0 || layer >= res->max_layer_index() || res->shape_nodes[layer].empty())
        return ipoint2(-1, -1);

    if (!res->grid(layer))
        res->init_grid(layer);

    ipoint2 ip(int_round(p.x), int_round(p.y));
    int minr2 = minr*minr;
    int maxr2 = maxr*maxr;
    int minx = std::max<int>(0, ip.x - maxr);
    int maxx = std::min<int>(res->x_size(layer) + 1, ip.x + maxr + 1);
    int miny = std::max<int>(0, ip.y - maxr);
    int maxy = std::min<int>(res->y_size(layer) + 1, ip.y + maxr + 1);
    ipoint2 result(-1, -1);
    int resultd2 = INT_MAX;

    for (int x = minx; x < maxx; ++x) {
        for (int y = miny; y < maxy; ++y) {
            ipoint2 q(x, y);
            int d2 = (ip - q).norm2();

            if (d2 >= minr2 && d2 <= maxr2) 
                if (d2 < resultd2) { result = q; resultd2 = d2; }
        }
    }
    if (resultd2 < INT_MAX)
        return result;
    else
        return get_closest_position(res, layer, p, maxr, 2*maxr);
}


// mean_shift inhibition
// 
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

/*
vector<double> curve_gaps(graph* g, const vector<ipoint2>& pts)
{
    set<node*> visited;
    vector<double> result;

    for (graph::iter_t niter = g->nodes.begin(); niter != g->nodes.end(); ++niter) {
        node* n = *niter;
        node_data_t<int>* nd = (node_data_t<int>*)n->data;

        forall_neighbors (n, iter) {
            node* nn = neighbor_node(iter);
            node_data_t<int>* nnd = (node_data_t<int>*)nn->data;
            
            if (visited.find(nn) == visited.end()) 
                result.push_back(sqrt((double)(pts[nnd->data].distance2(pts[nd->data]))));
        }
        visited.insert(n);
    }
    return result;
}
*/

// Auxillary function for 'point_matching_p'
//void match_paths(vector<int>& result, double& val, 
//    graph* g1, graph* g2, const vector<node*> lpath1, const vector<node*> lpath2,
//    const vector<ipoint2>& pts1, const vector<ipoint2>& pts2)
//{
//    double f = (double)(lpath2.size() - 1)/(lpath1.size() - 1);
//    set<node*> matched; // Matched nodes of g2.
//
//    val = 0.0; // Sum of distances^2
//    result.resize(pts2.size());
//    
//    // Match lpath1 to lpath2
//    for (int i = 0; i < (int)lpath1.size(); ++i) {
//        int j = int_round(f*i);
//        node* n1 = lpath1[i];
//        node* n2 = lpath2[j];
//        int i1 = ((node_data_t<int>*)n1->data)->data;
//        int i2 = ((node_data_t<int>*)n2->data)->data;
//
//        val += pts1[i1].distance2(pts2[i2]);
//        result[i2] = i1;
//        matched.insert(n2);
//    }
//
//    // Match rest
//    set<node*> unmatched;
//
//    set_difference(unmatched, set<node*>(g2->nodes.begin(), g2->nodes.end()), matched);
//    while (!unmatched.empty()) {
//        node* n = *unmatched.begin();
//        vector<node*> p = dfs_closest_path(g2, matched, n);
//
//        for (int pi = 0; pi + 1 < (int)p.size(); ++pi) {
//            node* pn = p[pi];
//            int i2 = ((node_data_t<int>*)pn->data)->data;
//            int i1 = ((node_data_t<int>*)p.back()->data)->data;
//            
//            result[i2] = result[i1];
//            val += pts1[result[i1]].distance2(pts2[i2]);
//            unmatched.erase(pn);
//            matched.insert(pn);
//        }
//    }
//}
//
//void save_matching_dbg(const string& fname, const vector<int>& perm, const vector<dpoint2>& pts1, 
//    const vector<dpoint2>& pts2)
//{
//    ofstream os(fname.c_str());
//
//    for (int i = 0; i < perm.size(); ++i) {
//        os << (perm[i] + 1) << ',';
//        if (i < pts1.size()) os << pts1[i].x << ',' << pts1[i].y << ',';
//        else os << 0 << ',' << 0 << ',';
//        os << pts2[i].x << ',' << pts2[i].y << endl;
//    }
//    os.close();
//}
//
//
//// Match points 'pts1' to points 'pts2' using matching of longest paths in the 
//// minimum spanning trees; retruns mapping 'm' such that 
//// pts1[m[i]] <-> pts2[i]. Note that 'm' is not necessarily a permutation!
//vector<int> point_matching_p(const vector<ipoint2>& pts1, const vector<ipoint2>& pts2)
//{
//    graph* g1 = point_graph(pts1);
//    vector<node*> lpath1 = tree_longest_path(g1);
//    graph* g2 = point_graph(pts2);
//    vector<node*> lpath2 = tree_longest_path(g2);
//
//    vector<int> result;
//    double val;
//    vector<int> resultinv;
//    double valinv;
//
//    match_paths(result, val, g1, g2, lpath1, lpath2, pts1, pts2);
//    //save_matching_dbg("c:\\work\\match.csv", result, cast_vector<dpoint2, ipoint2>(pts1), 
//    //    cast_vector<dpoint2, ipoint2>(pts2));
//    reverse(lpath1.begin(), lpath1.end());
//    match_paths(resultinv, valinv, g1, g2, lpath1, lpath2, pts1, pts2);
//    //save_matching_dbg("c:\\work\\matchinv.csv", resultinv, cast_vector<dpoint2, ipoint2>(pts1), 
//    //    cast_vector<dpoint2, ipoint2>(pts2));
//
//    delete g2;
//    delete g1;
//
//    if (valinv < val) return resultinv;
//    else return result;
//}
//

template<class A, class B> void save(const string& fname, const vector<pair<A, point2<B> > >& v)
{
    ofstream os(fname.c_str());

    for (int i = 0; i < (int)v.size(); ++i) {
        os << v[i].first << ',' << v[i].second.x << ',' << v[i].second.y << '\n';
    }
    os.close();
}

void test5(const char* arg1, const char* arg2)
{
    //ifstream is("c:\\work\\pts.m");
    //vector<dpoint2> pts;

    //while (is.good()) {
    //    double x, y;
    //    is >> x; 
    //    is.ignore(INT_MAX, ',');
    //    is >> y; 
    //    if (is.good()) pts.push_back(dpoint2(x, y));
    //    is.ignore(INT_MAX, '\n');
    //}
    //is.close();
    //
    //dpoint2 x;
    //double diff = mean_shift(x, pts, pts[1], 2.5, -1/8.0, 1.0E-6, 100);

    //cout << "Diff: " << diff << endl;
    //cout << "x: " << x.x << ',' << x.y << endl;
    
    // _____________ LIBRARY ___________
    part_lib* library;
    int layer = 5;
    int type = 23;

    read_library(arg1, library);

    node* p = library->parts[layer][type];
    //vs_part_data* vspd = dynamic_cast<vs_part_data*>(p->data);

    //if (vspd == NULL) cerr << "0 does not exist!" << endl;

    //vector<dpoint2> pts = partition(vspd->pcad.mean);
    //vector<ipoint2> ipts = cast_vector<ipoint2, dpoint2, double>(pts, 100.0);
    vector<pair<int, ipoint2> > ipts = get_library_geo_pieces(p, 5);

    //ipts = inhibit_point_set(ipts, 5);
    translate_and_scale(ipts);
    save("c:\\work\\pts1.csv", ipts);

    vector<ipoint2> ppts = extract_second<ipoint2>(ipts.begin(), ipts.end());
    graph* g = point_graph(ppts);
    point_graph_image(g, ppts).save("c:\\work\\pts1.png");

    delete g;

    //vector<int> match = piecewise_point_matching_p(ipts, ipts2);

    //save_matching_dbg("c:\\work\\match.csv", match, ppts, ppts2);

    
    // ________________ IMAGE ________________
    int toprev = atom("toPrevLayer");
    int to0 = atom("toLayer0");

    typedef vector<pair<int, ipoint2> > point_vector_t;


    layer1_result* res;

    read_layer1_result(res, arg2);
    for (int i = 0; i < (int)res->shape_nodes[layer].size(); ++i) {
        node* n = res->shape_nodes[layer][i];

        while (n != NULL) {
            layer1_data* nd = (layer1_data*)n->data;

            if (nd->m == type) {
                point_vector_t ptsm;

                foreach_neighbor (n, toprev, iter) {
                    node* m = neighbor_node(iter);
                    edge_data_name* med = (edge_data_name*)neighbor_edge_data(iter);
                    layer1_data* md = (layer1_data*)m->data;
                    set<node*> nset;
                    set<ipoint2> pset;

                    res->recurse_and_link(m, toprev, to0, nset);
                    node_set_to_point_set(pset, nset.begin(), nset.end());

                    for (set<ipoint2>::iterator pseti = pset.begin(); pseti != pset.end(); ++pseti) {
                        ptsm.push_back(pair<int, ipoint2>(med->index, *pseti));
                    }
                }
                ptsm = inhibit_point_set(ptsm, 5);
                translate_and_scale(ptsm);

                save("c:\\work\\ptsm.csv", ptsm);
               
                vector<int> match = piecewise_point_matching_p(ipts, ptsm);
                vector<ipoint2> ptsm2 = extract_second<ipoint2>(ptsm.begin(), ptsm.end());

                save_matching_dbg("c:\\work\\match.csv", match, ppts, ptsm2);

                cout << "Check result and enter some int: ";
                int j;
                cin >> j;
            }
            n = nd->next;
        }
    }

    //node* p2 = library->parts[layer][4];
    //vs_part_data* vspd2 = dynamic_cast<vs_part_data*>(p2->data);

    //if (vspd2 == NULL) cerr << "4 does not exist!" << endl;

    //vector<dpoint2> pts2 = partition(vspd2->pcad.mean);
    //vector<ipoint2> ipts2 = cast_vector<ipoint2, dpoint2, double>(pts2, 100.0);

    //graph* g2 = point_graph(ipts2);
    //point_graph_image(g2, ipts2).save("c:\\work\\pts2.png");

    //point_matching_p(ipts, ipts2);


    
    //for (int i = library->parts[layer].size() - 1; i > 0; --i) {
    //    node* p = library->parts[layer][i];
    //    vs_part_data* vspd = dynamic_cast<vs_part_data*>(p->data);

    //    if (vspd == NULL) continue;

    //    vector<dpoint2> pts = partition(vspd->pcad.mean);
    //    vector<ipoint2> ipts = cast_vector<ipoint2, dpoint2, double>(pts, 100.0);

    //    graph* g = point_graph(ipts);
    //    point_graph_image



    //    /*vector<double> gaps = curve_gaps(g, ipts);
    //    double m = median_value(gaps, 0.8);

    //    cerr << i << "; Max: " << *max_element(gaps.begin(), gaps.end()) << "; med: " << m << endl;

    //    //if (!continuous_curve(g, ipts)) cerr << "Part " << i << " is not continuous." << endl;*/

    //    //point_graph_image(g, ipts).save("c:\\work\\tmp\\part" + fill_left(i, 3) + ".png");

    //    vector<node*> lpath = tree_longest_path(g);

    //    cout << "Longest path:";
    //    cout << " size: " << lpath.size() << "; total nodes: " << g->nodes.size() << endl;

    //    //for (int i = 0; i < lpath.size(); ++i) {
    //    //    node_data_t<int>* nd = (node_data_t<int>*)lpath[i]->data;
    //    //    cout << ',' << ipts[nd->data].x << ',' << ipts[nd->data].y;
    //    //}
    //    //cout << endl;
    //    //for (int i = 0; i < ipts.size(); ++i) {
    //    //    cout << ',' << ipts[i].x << ',' << ipts[i].y;
    //    //}

    //    set<node*> diff;

    //    set_difference(diff, set<node*>(g->nodes.begin(), g->nodes.end()),
    //        set<node*>(lpath.begin(), lpath.end()));

    //    cout << "Difference: " << diff.size() << endl;
    //    for (set<node*>::iterator diter = diff.begin(); diter != diff.end(); ++diter) {
    //        vector<node*> p = dfs_closest_path(g, set<node*>(lpath.begin(), lpath.end()), *diter);
    //        cout << p.size() << ' ';
    //    }
    //    
    //    cout << endl;

    //    int l;
    //    cin >> l;

    //    delete g;
    //}

    delete library;
}


int _tmain(int argc, _TCHAR* argv[])
{
    _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
    //_CrtSetBreakAlloc(338);
    init_atoms();
    init_streaming();
    
    //int* boo = new int;

    cerr << "test (" __DATE__ " / " __TIME__ ")" << endl;
    srand((unsigned)time(0)); 
    //srand(2009);

    cerr << argc << " arguments given" << endl;

    switch (argc) {
        case 1 : test1(); break;
        case 2 : test2(); break;
        case 3 : test3(); break;
        case 4 : test4(); break;
        case 5 : test5(argv[1], argv[2]); break;
        default: 
            cout << "test files ['m' | 'n']" << endl;
            cout << "  (m = Mathematica, n = Pajek)" << endl;
    }
    return 0;
}

