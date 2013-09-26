/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */

// hop-test.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <stdio.h>
#include "../../interface/hop.h"


const char* cfg1 = "type = struct; init_size = -150; scale_limit = 100;"
    "layer1_threshold = 0.1; power_correction = 0.6; scale_sigma = 1.0;"
    "scale_factor = 0.707107";

const char* cfg2 = "layer_index = 2; layer_contraction = 2.0;"
    "r_response_threshold = 0.3; g_response_threshold = 0.3; candidate_r_threshold_percent = 0.8;"
    "candidate_g_threshold_percent = 0.8; realization_ratio_threshold = 0.9";

const char* cfg3 = "layer_index = 3; layer_contraction = 2.0;"
    "r_response_threshold = 0.3; g_response_threshold = 0.3; candidate_r_threshold_percent = 0.8;"
    "candidate_g_threshold_percent = 0.8; realization_ratio_threshold = 0.9";

const char* dcfg = "mode = r; layer = 0; end_layer = 0";

int main(int argc, char* argv[])
{
    // Read image, read library
    hop_image im = hop_read_image(argv[1]);

    hop_library lib = hop_read_library(argv[2]);
    

    std::cout << "Image dimensions: " << im.get_width() << " x " << im.get_height() << std::endl;

    // Inference
    hop_result* result;
    int count;

    count = hop_inference(result, im, cfg1);     // layer1
    hop_inference(result[0], lib, cfg2);  // layer2
    hop_inference(result[0], lib, cfg3);  // layer3

    // Display rarray[0]
    hop_image resim = hop_display(result[0], dcfg);
    hop_save_image(resim, "result.png");


    // rarray[0] -> blob
    hop_blob blb = result[0].to_blob();

    // blob -> object
    hop_result res(blb);
    
    // Check if it is identical
    resim = hop_display(res, dcfg);
    hop_save_image(resim, "result1.png");

    // image -> blob
    hop_blob imblob = im.to_blob();
    hop_image im2(imblob);
    hop_save_image(im2, "sample1-1.png");

    // Saving to file
    result[0].to_file("result.ly3", HOP_FORMAT_BLOB | HOP_FORMAT_COMPRESS);

    // Get one node on 3rd layer and its nodes
    hop_node n = result[0].get_node(2, 0);  
    hop_nodes nc = result[0].get_child_nodes(n);

    std::cout << "node position: (" << n.x << ", " << n.y << ") node type: " << n.type << " response: " << n.r_response << std::endl;
    std::cout << "Children count: " << nc.size << std::endl;
    for (int i = 0; i < nc.size; ++i) {
        std::cout << "   #" << i + 1 << " position: (" << nc.items[i].x << ", " << nc.items[i].y 
            << ") node type: " << nc.items[i].type << " response: " << nc.items[i].r_response << std::endl;
    }

    // Get child node indices
    hop_indices indices = result[0].get_child_nodes(2, 0);

    std::cout << "Children count: " << indices.size << std::endl;
    for (int i = 0; i < indices.size; ++i) {
        std::cout << "   #" << i << ": " << indices.items[i] << std::endl;
    }
    delete[] indices.items;
    
    delete[] result;
    return 0;
}

