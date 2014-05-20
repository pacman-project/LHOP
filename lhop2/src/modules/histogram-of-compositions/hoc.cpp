/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// layer_1 

#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>

#if defined WIN32 | defined WIN64
#include <crtdbg.h>
#endif

#include "hoc.h"

#include "utils/utils.h"
#include "utils/graphs/graph_utils.h"


using namespace std;

namespace hoc {

histogram_distance_function* histogram_distance_function::new_instance(const string& type_str, const ConfigDictionary& cfg) {
	histogram_distance_function* distance_func = nullptr;
	
	if (type_str.compare("gaussian") == 0 ) {
		distance_func = new gaussian_distance_f(cfg.getValueInt("sigma", 1.0, true));
	} else if (type_str.compare("weighted_centers") == 0 ) {
		distance_func = new all_cetners_distance_f();
	} else if (type_str.compare("round_overlap") == 0 ) {
		distance_func = new rounded_overlap_distance_f(cfg.getValueDouble("sigma", 1, true), cfg.getValueDouble("region_extend", 1.2, true));
	} else if (type_str.compare("multikernel_gaussian") == 0 ) {
		distance_func = new multikernel_gaussian_distance_f(cfg.getValueDouble("sigma", 0.7, true), cfg.getValueDouble("threshold", 0.1, true));		
	} else if (type_str.compare("none") != 0) {
		throw new_libhop_exception("Invalid bin_distance_type (supported only: 0 (none), 1 (gaussian), 2 (rounded_overlap), 4 (multikernel gaussian)." );
	}	

	return distance_func;
}

// factory method for histograming implementations
histograming_impl* histograming_impl::new_instance(part_lib* library, const int current_layer, const string& type_str, const ConfigDictionary& cfg) {
	histograming_impl* hist_impl = nullptr;

	bool include_similar_parts = cfg.getValueBool("parts_similarity.enable", false);
	float similarity_sigma = cfg.getValueBool("parts_similarity.sigma", 0.25);

	bool use_part_merging = cfg.getValueBool("use_part_merging", true);

	int type_count = library->layer_size(current_layer);

	vector<int> merged_parts_indexes = library->get_root_parts_map(current_layer);		
			
	if (merged_parts_indexes.size() > 0 && use_part_merging) {
		// calculate how many different parts we will have by searching for highest index
		type_count = 0;
		for (int i = 0; i < merged_parts_indexes.size(); i++)
			type_count = max(type_count, merged_parts_indexes[i]);
		// increase by one as we only found max (zero based) index
		type_count++;
	} else {
		// fill mapping with identety (indexes maps to same value)
		merged_parts_indexes.resize(type_count);
		for (int i = 0; i < type_count; i++)
			merged_parts_indexes[i] = i;
	}

	if (type_str.compare("plain_hoc") == 0) {
		hist_impl = new plain_hoc_histograming(current_layer, type_count, merged_parts_indexes, include_similar_parts, similarity_sigma);
	} else if (type_str.compare("upper_layer_filtering") == 0) {
		hist_impl = new upper_layer_filter_histograming(cfg.getValueInt("top_part_layer", current_layer), current_layer, type_count, merged_parts_indexes);
	} else if (type_str.compare("dense_featureset") == 0) {

		// create plain histograming implementation for local features
		//histograming_impl* local_feature_hist_impl = new plain_hoc_histograming(current_layer, type_count, merged_parts_indexes, include_similar_parts, similarity_sigma);
		string feature_hist_str("plain_hoc");

		cfg.getValue(feature_hist_str, "feature_histograming", false);

		ConfigDictionary feature_hist_cfg;
		feature_hist_cfg.fromNamespacePriority(cfg, 1, "feature_histograming");

		histograming_impl* local_feature_hist_impl = histograming_impl::new_instance(library, current_layer, feature_hist_str, feature_hist_cfg);
		
		// create requested binning
		// use 'grid' as default binning for feature
		string feature_bin_str("grid_binning");

		cfg.getValue(feature_bin_str, "feature_binning", false);

		ConfigDictionary feature_bin_cfg;
		feature_bin_cfg.fromNamespacePriority(cfg, 1, "feature_binning");

		binning_impl* local_feature_bin_impl = binning_impl::new_instance(library, current_layer, feature_bin_str, feature_bin_cfg);
		
		int stride = cfg.getValueInt("stride", 2);
		int window_size_x = cfg.getValueInt("window_size_x", 16);
		int window_size_y = cfg.getValueInt("window_size_y", 16);

		string cluster_filename;
		cfg.getValue(cluster_filename, "cluster_centers", true);

		// read all values into array				
		matrix<float> values_mat = read_float_matrix(cluster_filename);
		
		histogram_clusters* cluster_centers = new histogram_clusters(&values_mat[0], values_mat.width, values_mat.height);

		// create dense featureset implementation of histograming with standard matrix binning
		hist_impl = new dense_featureset_histograming(current_layer, stride, window_size_x,window_size_y, cluster_centers, local_feature_hist_impl, local_feature_bin_impl);
	} else {
		throw new_libhop_exception("Invalid histogram_implementation (supported only: plain_hoc (default), upper_layer_filtering, dense_featureset." );
	}

	return hist_impl;
}

// factory method for binning implementations
binning_impl* binning_impl::new_instance(part_lib* library, const int current_layer, const string& type_str, const ConfigDictionary& cfg) {
	binning_impl* bin_impl = nullptr;

	if (type_str.compare("matrix_binning") == 0) {
		
		bool legacy_compatible = cfg.getValueBool("legacy_compatible", true);	
		bool interpolate_to_neighboors = cfg.getValueBool("interpolate_to_neighboors", false);	
		
		histogram_distance_function* distance_func = nullptr;
		string dist_func_str;

		if (cfg.getValue(dist_func_str, "bin_distance_function", false)) {
			ConfigDictionary cfg_dist_func;
			// copy values based on namespace hierarhy
			cfg_dist_func.fromNamespacePriority(cfg,1,"bin_distance_function");

			distance_func = histogram_distance_function::new_instance(dist_func_str, cfg_dist_func);
		}

		string matrix_filename;
		cfg.getValue(matrix_filename, "matrix", true);

		bin_impl = new matrix_binning(new M_bin(matrix_filename), distance_func, interpolate_to_neighboors, legacy_compatible);
	} else if (type_str.compare("spatial_pyramid_match") == 0) {
		bin_impl = new spm_binning(cfg.getValueInt("pyramid_layers", 3));
	} else if (type_str.compare("grid_binning") == 0) {
		bin_impl = new grid_binning(cfg.getValueInt("width", 0),cfg.getValueInt("height", 0));
	} else {
		throw new_libhop_exception("Invalid binning_implementation (supported only: matrix_binning (default), spatial_pyramid_match, grid_binning." );
	}
	bin_impl->set_current_layer(current_layer);

	return bin_impl;
}

void hoc_histogram_generator::init_cfg(const ConfigDictionary& cfg) {
	
	// get list of layers on which we will process 
    cfg.getValue(layers, "layer", true);    

	// sliding windows cfg
    sw_xstep = cfg.getValueInt("sliding_window_xstep", 0);
    sw_ystep = cfg.getValueInt("sliding_window_ystep", 0);
	sw_width = cfg.getValueInt("sliding_window_width", 0);
    sw_height = cfg.getValueInt("sliding_window_height", 0);

	padding_x = cfg.getValueInt("padding_x", 0);
	padding_y = cfg.getValueInt("padding_y", 0);

	// we should always be using img coordinates 
	bool use_img_coordinates = cfg.getValueBool("use_img_coordinates", true);	// make sure it is not FALSE as this is not supported any more !!!!
	
	normalize = cfg.getValueBool("normalize", false);	
	output_headers_only = cfg.getValueBool("output_headers_only", false);	
	
	// histograming and binning implementation (use "plain_hoc" and "matrix" by default as lagacy support)
	string histograming_str("plain_hoc");
	ConfigDictionary cfg_hist;

	if (cfg.getValue(histograming_str, "histogram_implementation", false))
		cfg_hist.fromNamespacePriority(cfg,1,(const char*)histograming_str.c_str());
	else 
		cfg_hist = cfg;

	string binning_str("matrix_binning");
	ConfigDictionary cfg_bin;
	if (cfg.getValue(binning_str, "binning_implementation", false))
		cfg_bin.fromNamespacePriority(cfg,1,(const char*)binning_str.c_str());
	else
		cfg_bin = cfg;

	// create seperate instance of histograming/binning for each layer due to better optimization (matrix_binning can use matrix with specific size for each layer - no need to resize it each time - also distance function caches its values on this matrix)
	for (int i = 0; i < layers.size(); i++) {
		hist_impl.push_back(histograming_impl::new_instance(library, layers[i], histograming_str, cfg_hist));
		bin_impl.push_back(binning_impl::new_instance(library, layers[i], binning_str, cfg_bin));
	}
}

// returns list of sliding windows generated for specific image size that has been rescaled using scale_factor
// sliding windows parameters are applied on original unscaled image coordinates except for (sw_xstep, sw_ystep) which is multipled by scale_factor
// regions that are returned are in original unscaled image coordinates
// 
std::list<frectangle2> hoc_histogram_generator::generate_sliding_windows(layer1_result* res, double scale_factor) {

	std::list<frectangle2> sliding_windows_regions;

	fpoint2 bin_size(sw_width* scale_factor, sw_height* scale_factor);
	float image_width = (res->x_size(0)- 2*res->border) * scale_factor;
	float image_height = (res->y_size(0)- 2*res->border) * scale_factor;

	//printf("xpos = %d; xpos + %f <= %f + %d; xpos += %d\n",-padding_x, bin_size.x, image_width, padding_x, sw_xstep);
	//printf("\typos = %d; ypos + %f <= %f + %d; ypos += %d\n",-padding_y, bin_size.y, image_height, padding_y, sw_ystep);
			
	for (float xpos = -padding_x; xpos + bin_size.x <= image_width + padding_x; xpos += sw_xstep*scale_factor) {
		for (float ypos = -padding_y; ypos + bin_size.y <= image_height + padding_y; ypos += sw_ystep*scale_factor) {
			fpoint2 ll(xpos, ypos);
			sliding_windows_regions.push_back(frectangle2(ll, ll + bin_size));
		}	  
	}

	return sliding_windows_regions;
}
vector<histogram_descriptor>* hoc_histogram_generator::generate_descriptors(layer1_result* res, const std::list<irectangle2> regions) {
	
	int border = res->border;
	double scale_factor = res->original_width / ((double)(res->x_size(0) - 2*border));

	std::vector<histogram_descriptor>* merged_descriptors = nullptr;
	
	int layers_count = layers.size();
	
	bool has_layers = true;
	for (int i = 0; i < layers.size() && has_layers; i++)
		has_layers = layers[i] <= res->max_layer_index();

	if (has_layers) {
		for (int i = 0; i < layers_count; i++) {
			double contraction_factor = (double)res->x_size(0) / (double)res->x_size(layers[i]);		
		
			std::list<frectangle2> layer_i_regions;

			// generate either sliding windows across whole image or use user supplied window regions
			if (regions.size() <= 0) {
				if (sw_xstep > 0 && sw_ystep > 0)
					layer_i_regions = generate_sliding_windows(res, scale_factor);
				else
					throw new_libhop_exception("Invalid arguments: user must specify either region windows or define sliding windows parameters");
			} else {
				// just copy all rectangles (they should be in image coordinate system)
				for (std::list<irectangle2>::const_iterator rect = regions.begin(); rect != regions.end(); rect++) {
					layer_i_regions.push_back(frectangle2(rect->ll.x, rect->ll.y, rect->ur.x, rect->ur.y));
				}
			}

			// adjust all rectangle coords to layer coordinates if they are in image coordinate system
			for (std::list<frectangle2>::iterator rect = layer_i_regions.begin(); rect != layer_i_regions.end(); rect++) {
				rect->ll.x = (rect->ll.x / scale_factor + border) /contraction_factor;
				rect->ll.y = (rect->ll.y / scale_factor + border) /contraction_factor;
				rect->ur.x = (rect->ur.x / scale_factor + border) /contraction_factor;
				rect->ur.y = (rect->ur.y / scale_factor + border) /contraction_factor;
			}
			std::vector<histogram_descriptor>* descriptors = new std::vector<histogram_descriptor>();
			descriptors->resize(layer_i_regions.size());
			int j = 0;
	
			for (std::list<frectangle2>::const_iterator rect = layer_i_regions.begin(); rect != layer_i_regions.end(); rect++) {

				(*descriptors)[j].x = rect->ll.x;
				(*descriptors)[j].y = rect->ll.y;
				(*descriptors)[j].w = rect->ur.x - rect->ll.x;
				(*descriptors)[j].h = rect->ur.y - rect->ll.y;
				if (output_headers_only == false)
					(*descriptors)[j].hist = get_histogram(hist_impl[i], bin_impl[i], res, i, irectangle2(rect->ll.x,rect->ll.y,rect->ur.x,rect->ur.y), normalize);

				double sum = 0;
				for (std::vector<float>::iterator it = (*descriptors)[j].hist.begin(); it != (*descriptors)[j].hist.end(); it++) {
					sum += *it;
				}
				j++;
			}

			// merge with descriptors on previous layer unless this layer was processed first, then just use existing descriptor list
			if (merged_descriptors == nullptr) {
				// just use this descriptor list
				merged_descriptors = descriptors;
			
				// update coordinates if needed
				//if (use_img_coordinates == true)  {
				for (std::vector<histogram_descriptor>::iterator merged_desc = merged_descriptors->begin(); merged_desc != merged_descriptors->end(); merged_desc++) {
					(*merged_desc).x = ((*merged_desc).x * contraction_factor - border) * scale_factor;
					(*merged_desc).y = ((*merged_desc).y * contraction_factor - border) * scale_factor;
					(*merged_desc).w = ((*merged_desc).w * contraction_factor) * scale_factor;
					(*merged_desc).h = ((*merged_desc).h * contraction_factor) * scale_factor;
				}
				//}
			} else {
				std::vector<histogram_descriptor>::iterator merged_desc = merged_descriptors->begin();
				std::vector<histogram_descriptor>::iterator desc = descriptors->begin();
				while (desc != descriptors->end() && merged_desc != merged_descriptors->end()) {

					// update coordinates if needed
					//if (use_img_coordinates == true)  {
					(*desc).x = ((*desc).x * contraction_factor - border) * scale_factor;
					(*desc).y = ((*desc).y * contraction_factor - border) * scale_factor;
					(*desc).w = ((*desc).w * contraction_factor) * scale_factor;
					(*desc).h = ((*desc).h * contraction_factor) * scale_factor;
					//}

					// verify that bbox coordinates do match (at least with a few px of diff)
					if (::abs((*desc).x - (*merged_desc).x) > 3 ||
						::abs((*desc).y - (*merged_desc).y) > 3 || 
						::abs((*desc).w - (*merged_desc).w) > 3 ||
						::abs((*desc).h - (*merged_desc).h) > 3) {
						// error occured if one on same location or size
							throw new_libhop_exception("Error while merging histogram of different layers: widnow sizes of histogram do not match");
					}
				
					// merge both descriptors by inserting new ones at the end
					(*merged_desc).hist.insert((*merged_desc).hist.end(), (*desc).hist.begin(), (*desc).hist.end());

					// continue to next descriptor
					desc++; merged_desc++;
				}

				if (desc != descriptors->end() || merged_desc != merged_descriptors->end()) {
					// error occured if one did not come to end
					throw new_libhop_exception("Error while merging histogram of different layers: count of histogram do not match");
				}							 

				delete descriptors;
			}
		}
	} else {
		// no layers found - just return empty list
		merged_descriptors = new std::vector<histogram_descriptor>();
	}
	return merged_descriptors;
}


// main method for generating a single histogram from layer1_result 'res' on specific 'layer' withing specific region 'r' using class specific implemention of histograming and binning
// this method works using abstract classes while specific feature extraction and region binning is implemented in other classes (see histograming_impl and binning_impl)
// in general this method iterates over all features retuned by histograming implemnetation and iterates over each bin returned by binning implementations to account for each feature type in final histogram
vector<float> hoc_histogram_generator::get_histogram(histograming_impl* hist_impl, binning_impl* bin_impl, layer1_result* res, int layer, const irectangle2& r, const bool normalize) {

	histograming_impl::iterator* features_iter = hist_impl->get_features_iterator(res, r);
	
	binning_impl::iterator* bin_iter = bin_impl->get_bin_iterator(res, r);

	int max_feature_types = hist_impl->get_max_features();
	int number_bins = bin_impl->get_num_bins();

	vector<float> result(number_bins*max_feature_types, 0.0);
	
	// do as initialization
	features_iter->reset();

	// ignore zero size region
	if (r.x_dim() > 0 && r.y_dim() > 0) {
		while (features_iter->move_next()) {

			ipoint2 feat_location = features_iter->location() - r.ll;
		
			int feat_type = features_iter->type();
			float feat_value = features_iter->value();

			bin_iter->reset(feat_location);

			//printf("(n)part loc: %d,%d, part type %d\n",feat_location.x,feat_location.y,feat_type);

			while (bin_iter->move_next()) {

				int hist_location = bin_iter->bin_index();
				float bin_weight = bin_iter->bin_weight();

				result[hist_location * max_feature_types + feat_type] += feat_value * bin_weight;
			}		
		}

		if (normalize) {
			// normalize to sum == 1
			float sum = 0;
			for (int i = 0; i < result.size(); i++) 
				sum += result[i];
			if (sum != 0) {
				for (int i = 0; i < result.size(); i++) 
					result[i] /= sum;
			}
		}
	}
	delete features_iter;
	delete bin_iter;

	return result;
}

}