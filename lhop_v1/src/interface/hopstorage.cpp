/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */

#ifdef WIN32
#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN // Exclude rarely-used stuff from Windows headers
#include <windows.h>
#endif

#include <stdio.h>

#include "layers/layer_n_creators.h"
#include "layers/layer_1_creators.h"

#include "hopstorage.h"

using namespace std; 

int** _hop_child_indices(const hop_result* res, int layer)
{
    typedef map<node*, int> map_t;

    if (!res) 
        return 0;

    layer1_result* result = (layer1_result*)res->data();
    
    if (!result || layer < 1 || layer > result->max_layer_index())
        return 0;
    
    int edgename = atom("toPrevLayer").get_index();
    map_t indexmap;
    vector<node*>& snodes = result->shape_nodes[layer];
    vector<node*>& snodes1 = result->shape_nodes[layer - 1];
    int** carray = new int*[snodes.size()];
    vector<int> cvector;

    cvector.reserve(100);
    for (int i = 0; i < (int)snodes1.size(); ++i) {
        indexmap.insert(map_t::value_type(snodes1[i], i));
    }
    for (int i = 0; i < (int)snodes.size(); ++i) {
        carray[i] = nullptr;
        node* n = snodes[i];

        cvector.clear();
        foreach_neighbor(n, edgename, iter) {
            map_t::iterator miter = indexmap.find(neighbor_node(iter));

            if (miter != indexmap.end()) cvector.push_back(miter->second);
        }
        cvector.push_back(-1);
        carray[i] = new int[cvector.size()];
        memcpy(carray[i], &cvector.at(0), sizeof(int)*cvector.size());
    }
    return carray;
}

void _hop_store_result_blob(const hop_result* res, const char* fname, bool compress) {
    if (compress)
        ((streamable*)(res->data()))->save(fname, -1);
    else ((streamable*)(res->data()))->save(fname);
}

void _hop_store_result_protobuf(const hop_result* res, const char* fname, bool compress) {
#ifdef PROTOBUF_SUPPORT
	using namespace com::vicos::hop::proto;
	using namespace google::protobuf::io;

    using namespace com::vicos::hop::proto;
    using namespace google::protobuf::io;

    Pyramid p;

    for (int i = 0; i < res->get_layer_count(); i++) {

        Layer* l = p.add_layer();

        // inverted coordinates !!
        l->set_width(res->get_width(i));
        l->set_height(res->get_height(i));

        int** indices = _hop_child_indices(res, i);

        for (int j = 0; j < res->get_node_count(i); j++) {
            Layer_Node* n = l->add_node();
            hop_node hn = res->get_node(i, j);

            // inverted coordinates !!
            n->set_x(hn.x);
            n->set_y(hn.y);
            n->set_weight(hn.r_response);
            n->set_id(hn.type);

            if (indices && indices[j]) {
                for (int k = 0; indices[j][k] > -1; k++) {
                    n->add_link(indices[j][k]);
                }

                delete [] indices[j];
            }
        }
        if (indices)
            delete [] indices;
    }

    fstream out(fname, ios::out | ios::binary | ios::trunc); 
    ZeroCopyOutputStream *outs = nullptr;
    GzipOutputStream *zips = nullptr;

    OstreamOutputStream *osts = new OstreamOutputStream(&out);
    if (compress) {
        zips = new GzipOutputStream(osts);
        outs = zips;
    } else outs = osts;
    
    p.SerializeToZeroCopyStream(outs); 

    if (compress)
        delete zips;
    delete osts;

    out.close();

#else
    throw hop_exception("Protobuf support not enabled during compilation!");
#endif

}




