/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// layer_1

#pragma once
#ifndef _HOC_1_
#define _HOC_1_

#include <algorithm>
#include <utility>
#include <bitset>
#include <vector>

#include "../utils/img.h"
#include "../graphs/img_graph.h"
#include "../utils/matrix.h"
#include "../utils/misc.h"
#include "../utils/utils.h"

#include "layer_1_result.h"

using namespace std;

namespace hoc {

struct histogram_distance_function {
	//virtual double get_normalization_value(M_bin& b, ipoint2 loc, int current_bin) = 0;
	virtual double operator()(M_bin& b, ipoint2 loc, int current_bin) = 0; 
	static histogram_distance_function* new_instance(const string& type_str, const config_dictionary& cfg);
};

struct no_distance_f : public histogram_distance_function {
	//virtual double get_normalization_value(M_bin& b, ipoint2 loc, int current_bin) { return 1;}
	virtual double operator()(M_bin& b, ipoint2 loc, int current_bin) { return 1.0; }
};

struct gaussian_distance_f : public histogram_distance_function {
	double sigma;
	gaussian_distance_f(double s) : sigma(s) {}
	/*
	virtual double get_normalization_value(M_bin& b, ipoint2 loc, int current_bin) {
		//int bin = b.get_bin(loc.x, loc.y);
		int bin = current_bin; // should use current bin as area normalization
		return 1 / ::sqrt((double)b.bin_area(bin));
	}*/
	virtual double operator()(M_bin& b, ipoint2 nd_loc, int current_bin) { 
		//int bin = b.get_bin(loc.x, loc.y);
		// use distance to current bin, but normalized to sqrt(area) of bin
		double distance = ::sqrt((double)b.bin_center(current_bin).distance2(nd_loc) / (double)b.bin_area(current_bin));
		return ::exp(- (distance * distance) / (2 * sigma * sigma)); // simple gaussian function
	}
};

class all_cetners_distance_user_object : public M_bin_user_object {
public:
    typedef pair<double, vector<int> > mapped_distance;
    mapped_distance* dist;
    
    virtual ~all_cetners_distance_user_object() {
      delete [] dist;
    }
    void allocate(size_t s) {
      dist = new mapped_distance[s];
      
      for (int i = 0; i <s; i++) {			    
	  dist[i] = mapped_distance();
	  dist[i].first = 0;
	  dist[i].second.resize(0);
      }
      
    }
    
    
};

struct histogram_sort_distances_method {
	bool operator() (pair<double,int> i,pair<double,int> j) { return (i.first<j.first);}
};

struct all_cetners_distance_f : public histogram_distance_function {
	
	virtual double operator()(M_bin& b, ipoint2 nd_loc, int current_bin) { 
		all_cetners_distance_user_object* saved_distances_obj = (all_cetners_distance_user_object*)b.get_user_value();
		
		if (saved_distances_obj == nullptr) {
			// create new map of distances if it does not exists
			saved_distances_obj = new all_cetners_distance_user_object();
			saved_distances_obj->allocate(b.width()*b.height());
			b.set_user_value(saved_distances_obj);

		}
		all_cetners_distance_user_object::mapped_distance* saved_distances = saved_distances_obj->dist;

		int saved_distance_index = (nd_loc.x + b.width()/2)  + (nd_loc.y + b.height()/2) * b.width();

		double sum = 0;

		all_cetners_distance_user_object::mapped_distance* loc_distance = saved_distances + saved_distance_index;

		if (loc_distance->second.size() <= 0)  {
			// calculate 4 nearest neighbors and save them to saved_distances

			// find 4 nearest neighbors and use only those
			vector<pair<double,int> > distances;
			
			// get all distances
			for (int i = 0; i < b.bin_count(); i++) {
				double d = nd_loc.distance2(b.bin_center(i));
				distances.push_back(make_pair(d,i));
			}

			// sort distances
			histogram_sort_distances_method sort_method;
			std::sort(distances.begin(),distances.end(), sort_method);
			
			loc_distance->second.insert(loc_distance->second.begin(),b.bin_count(), false);
			
			
			// calculate sum of first/closest 4 bins
			int closes_count = min<int>(4, distances.size());
			for (int i = 0; i < closes_count; i++) {
				pair<double,int>& d = distances[i];
				sum += 1/::sqrt(d.first);
				loc_distance->second[d.second] = true;
			}

			loc_distance->first = sum;
		}

		if (loc_distance->second[current_bin]) {
			double distance = ::sqrt((double)b.bin_center(current_bin).distance2(nd_loc)) *  loc_distance->first;
			if (distance != distance) // handle distance == inf case; (in this case this means very near i.e. at center and should return weight 1
				return 1; 
			else
				return 1/distance;
		} else {
			return 0; // handle zero normalized distance (in this case this means very far and should return weight 0
		}	
	}
};

struct rounded_overlap_distance_f : public histogram_distance_function {
	double sigma;
	double region_extention;
	rounded_overlap_distance_f(double s, double reg_ext) : sigma(s), region_extention(reg_ext) {}
	virtual double operator()(M_bin& b, ipoint2 nd_loc, int current_bin) { 
		
		// get center and max_distance for current bin
		float bin_max_dist2 = b.bin_max_distance2(current_bin) * region_extention;
		ipoint2 bin_center = b.bin_center(current_bin);

		// check if loc falls into this region 
		if (nd_loc.distance2(bin_center) > bin_max_dist2) {
			// not within region - return 0 as this location is not used
			return 0;
		} else {
			// assign region with proper normalization (i.e. gaussian weighting)
			double distance = ::sqrt((double)bin_center.distance2(nd_loc) / (double)b.bin_area(current_bin));

			return ::exp(- (distance * distance) / (2 * sigma * sigma)); // simple gaussian function
		}
	}
};

class multikernel_gaussian_distance_user_object : public M_bin_user_object {
public:
    vector<matrix<float> > gaussian_kernels;
    
	int width(int bin) { return gaussian_kernels[bin].width; }
	int height(int bin) { return gaussian_kernels[bin].height; }
    /*virtual ~multikernel_gaussian_distance_user_object() {
    }*/
    void calculate(M_bin& b, double sigma, double threshold) {      
		vector<float> kernel_dx,kernel_dy;
	
		kernel_dx.resize(b.bin_count());
		kernel_dy.resize(b.bin_count());
		
		float min_kernel_dx = INT_MAX;
		float min_kernel_dy = INT_MAX;

		ipoint2 c((int)b.width()/2, (int)b.height()/2);

		for (int i = 0; i < b.bin_count(); ++i) {
			int count = b.bin_area(i);
				
			ipoint2 p = b.bin_center(i);
			
			float dx = 0, dy = 0;

			for (int x = 0; x < b.width(); x++ ) {
				for ( int y = 0; y < b.height(); y++) {
					int bin = b.get_bin(x - c.x, y - c.y);
					if (bin == i) {

						dx += pow((float)(p.x - (x - c.x)),2);
						dy += pow((float)(p.y - (y - c.y)),2);
					}
				}
			}
			kernel_dx[i] = 1 / (float)(count - 1) * dx * sigma* sigma ;		
			kernel_dy[i] = 1 / (float)(count - 1) * dy * sigma* sigma;

			min_kernel_dx = min<float>(kernel_dx[i],min_kernel_dx);		
			min_kernel_dy = min<float>(kernel_dy[i],min_kernel_dy);

			// use same value for dy
			//kernel_dy[i] = kernel_dx[i];
			//min_kernel_dy = min_kernel_dx;

			//kernel_dx[i] = min<float>(kernel_dx[i],kernel_dy[i]);		
			//kernel_dy[i] = min<float>(kernel_dx[i],kernel_dy[i]);		
		}

		gaussian_kernels.resize(b.bin_count());

		// pre-calculate gaussian distribution for each bin
		for (int i = 0; i < b.bin_count(); ++i) {
			matrix<float> k(b.width(), b.height());
			for (int x = 0; x < b.width(); x++ ) {
				for ( int y = 0; y < b.height(); y++) {
					float xx = b.bin_center(i).x - (x - c.y);
					float yy = b.bin_center(i).y - (y - c.y);

					float e_val = (xx * xx) / (2* kernel_dx[i]*kernel_dx[i]) + (yy * yy) / (2* kernel_dy[i]*kernel_dy[i]);
					//float e_val = (xx * xx) / (2* min_kernel_dx*min_kernel_dx) + (yy * yy) / (2* min_kernel_dy*min_kernel_dy);
					k(x,y) = 1*exp(-e_val);
					if(k(x,y) != k(x,y)) {
						cout << "got nan in make gaussian_kernels: will set value to 1 if distance 0 otherwise to will set to 0"<< endl;
						k(x,y) = (xx == 0 && yy == 0) ? 1 : 0;
					}
				}
			}
			gaussian_kernels[i] = k;
		}
		
		// normalize gaussian distribution so on each location sum over all bins is 1
		for (int x = 0; x < b.width(); x++ ) {
			for ( int y = 0; y < b.height(); y++) {
				// get sum
				float sum_gauss = 0;
				for (int i = 0; i < b.bin_count(); ++i) {
					sum_gauss += gaussian_kernels[i](x,y);
					if(sum_gauss != sum_gauss)
						cout << "got nan in make sum of normalization of gaussian_kernels 1"<< endl;
				}
				// normalize and clip values lower then 0.1 (only if sum_gauss is valid number otherwise we will divide by zero (we get exception or NaN)
				if (sum_gauss <= 0) 
					continue;

				float norm_sum_gauss = 0;
				for (int i = 0; i < b.bin_count(); ++i) {
					float& v = gaussian_kernels[i](x,y);
					v  =  v/sum_gauss;

					if(v != v)
						cout << "got nan in make sum of normalization of gaussian_kernels 2"<< endl;

					if (v <= threshold)
						v = 0;
					norm_sum_gauss += v;
				}		
				// normalize again with clipped values (only if norm_sum_gauss is valid number otherwise we will divide by zero (we get exception or NaN)
				if (norm_sum_gauss <= 0)  {
					continue;
				}
				for (int i = 0; i < b.bin_count() ; ++i) {			
					gaussian_kernels[i](x,y) /= norm_sum_gauss;

					if(gaussian_kernels[i](x,y) != gaussian_kernels[i](x,y))
						cout << "got nan in make sum of normalization of gaussian_kernels 3"<< endl;

				}
			}
		}
    }

	float bin_gaussian_kernel(int b, int i, int j) {
		if  (b < 0 || b >= gaussian_kernels.size()) 
			return -1;
		
		// all input values to matrix have 0,0 coordinate for center (recalculate to matrix cooordinates)
		int x = (int)gaussian_kernels[0].width/2 + i, y = (int)gaussian_kernels[0].height/2 + j;

		if (x < 0 || x >= (int)gaussian_kernels[0].width || y < 0 || y >= (int)gaussian_kernels[0].height)
			return -1;

		return gaussian_kernels[b](x,y);
	}
    
};

struct multikernel_gaussian_distance_f : public histogram_distance_function {
	double sigma;
	double threshold;
	multikernel_gaussian_distance_f(double s, double t) : sigma(s), threshold(t) {}

	virtual double operator()(M_bin& b, ipoint2 nd_loc, int current_bin) { 

		multikernel_gaussian_distance_user_object* saved_mk_gaussian_obj = (multikernel_gaussian_distance_user_object*)b.get_user_value();
		
		if (saved_mk_gaussian_obj == nullptr) {
			// create new map of distances if it does not exists
			saved_mk_gaussian_obj = new multikernel_gaussian_distance_user_object();
			saved_mk_gaussian_obj->calculate(b, sigma, threshold);
			b.set_user_value(saved_mk_gaussian_obj);
		}
		if (saved_mk_gaussian_obj->width(current_bin) != b.width() || saved_mk_gaussian_obj->height(current_bin) != b.height())
			saved_mk_gaussian_obj->calculate(b, sigma, threshold);

		// use pre-computed gaussian kernel for this specific bin
		return saved_mk_gaussian_obj->bin_gaussian_kernel(current_bin, nd_loc.x, nd_loc.y);
	}
};



class histogram_descriptor {
public:
	float x, y, w, h;
	vector<float> hist;
};

class histograming_impl {
public:	
	class iterator {
	public:
		virtual ~iterator() {}
		virtual bool move_next() = 0;
  
		virtual void reset() = 0;
		virtual ipoint2 location() = 0; // must return location in image/layer coordinates (not relative to current rectangle region) but coordinates are in specific layer (they must be contracted etc.)
		virtual int type() = 0;
		virtual float value() = 0;
	};
	virtual int get_max_features() = 0;
	virtual iterator* get_features_iterator(layer1_result* res, const irectangle2& r) = 0; // rectangle region is in specific layer coordinates (they are already contracted) - implementation should make sure it knows which layer that is
	
	static histograming_impl* new_instance(part_lib* library, const int current_layer, const string& type_str, const config_dictionary& cfg);
};

class binning_impl {
public:	
	int current_layer;
	class iterator {
	public:
		virtual ~iterator() {}
		virtual bool move_next() = 0;
  
		virtual void reset(ipoint2 feat_location) = 0; // location returned by histograming_impl::iterator::location() but centered relative to current region rectangle location

		virtual int bin_index() = 0;
		virtual float bin_weight() = 0;
	};
	virtual int get_num_bins() = 0;
	virtual iterator* get_bin_iterator(layer1_result* res, const irectangle2& r) = 0; // rectangle region is in specific layer coordinates (they are already contracted) - implementation should make sure it knows which layer that is

	static binning_impl* new_instance(part_lib* library, const int current_layer, const string& type_str, const config_dictionary& cfg);

	void set_current_layer(int l) { current_layer = l; }
};

class hoc_histogram_generator {
	part_lib* library;

	vector<histograming_impl*> hist_impl;
	vector<binning_impl*> bin_impl;
	
	vector<int> layers;

	bool normalize;

	bool output_headers_only;

	int sw_xstep;
    int sw_ystep;
	int sw_width;
    int sw_height;

	int padding_x;
	int padding_y;

public:
	hoc_histogram_generator(part_lib* library) :  library(library) {}
	virtual ~hoc_histogram_generator() {
		// delete histograming and binning implementations as they were created by init_cfg
		for (int i = 0; i < hist_impl.size(); i++) {
			if (hist_impl[i] != nullptr)
				delete hist_impl[i];
		}
		for (int i = 0; i < bin_impl.size(); i++) {
			if (bin_impl[i] != nullptr)
				delete bin_impl[i];
		}
	}
	void init_cfg(const config_dictionary& cfg);

	// generates list of descriptors within specific regions of image 'res'
	vector<histogram_descriptor>* generate_descriptors(layer1_result* res, const std::list<irectangle2> regions);

	// static abstract implementation of main method that uses specific implementation of histograming and binning to produce histogram on 'res' on single 'layer' within region 'r'
	static vector<float> get_histogram(histograming_impl* hist_impl, binning_impl* bin_impl, layer1_result* res, int layer, const irectangle2& r, const bool normalize);

private:
	std::list<frectangle2> generate_sliding_windows(layer1_result* res, double scale_factor);
};



class plain_hoc_histograming : public histograming_impl {
	int layer;
	int maxfeat;
	vector<int> merged_parts_mapping;

	bool include_similar_parts;
	float similarity_sigma;
public:
	class iterator : public histograming_impl::iterator
	{
		layer1_result* res;

		plain_hoc_histograming* self;
		layer1_data* nd;
		
		// counters
		node* n;  
		int x,y;
		
		bool end_of_line;

		// start/stop conditions
		int minx,miny,maxx,maxy;
	public:
		iterator(plain_hoc_histograming* s, layer1_result* res, int minx, int miny, int maxx, int maxy) :
			self(s), res(res), minx(minx), miny(miny), maxx(maxx), maxy(maxy), n(nullptr), nd(nullptr), x(minx), y(miny-1), end_of_line(0) {}

		inline bool move_next() {
			// this should provide the same functionality as
			//for (int x = minx; x < maxx; ++x)
			//   for (int y = miny; y < maxy; ++y)
			//        node* n = node_at(x, y, layer);
			//        while (n != nullptr)
			//            n = nd->next;
	
			bool has_next = true;
			if (n != nullptr && nd->next != nullptr) {
				n = nd->next;
			} else {
				n = nullptr;
				while (n == nullptr) {
					if (y < maxy-1) {
						y++;
					} else {
						if (x < maxx-1) {
							x++;
						} else {
							// should stop iteration
							has_next = false;
							break;
						}
						// reset y
						y = miny;
					}
					// reset n
					n = res->node_at(x, y, self->layer);
				}
			}
			if (n != nullptr)
				nd = (layer1_data*)n->data;
			return has_next;
		}
		inline void reset() {
			if (!res->grid(self->layer)) res->init_grid(self->layer);
		}

		inline ipoint2 location() { return ipoint2(x,y); }
		inline int type() { return self->merged_parts_mapping[nd->m]; }
		inline float value() { return nd->vval(); }

	};

	plain_hoc_histograming(int layer, int maxfeat, vector<int> merged_parts_mapping, bool include_similar_parts, float similarity_sigma) : 
		layer(layer), maxfeat(maxfeat), merged_parts_mapping(merged_parts_mapping), include_similar_parts(include_similar_parts), similarity_sigma(similarity_sigma) {
		printf("using plain_hoc_histograming with settings:\n\tlayer: %d\n\tinclude_similar_parts: %d\n\tsimilarity_sigma: %f\n", layer, include_similar_parts,similarity_sigma); fflush(stdout);
	}
	
	inline virtual int get_max_features() { return maxfeat; }

	inline virtual histograming_impl::iterator* get_features_iterator(layer1_result* res, const irectangle2& r) {		
		ipoint2 c = r.center();
		int minx = max(0, c.x - r.x_dim()/2);
		int maxx = min(res->x_size(layer), c.x + r.x_dim()/2);
		int miny = max(0, c.y - r.y_dim()/2);
		int maxy = min(res->y_size(layer), c.y + r.y_dim()/2);

		return new iterator(this, res, minx, miny, maxx, maxy);
	}
};

class upper_layer_filter_histograming : public histograming_impl {
	int top_part_layer;
	int layer;
	int maxfeat;
	vector<int> merged_parts_mapping;

public:
	class layer_check
	{
	public:
		int layer;
		layer_check(int l) : layer(l) { }
		bool operator() (node* n, node* nn) const {
			return ((layer1_data*)nn->data)->z >= layer;
		}
		bool operator() (node* n) const {
			return ((layer1_data*)n->data)->z == layer;
		}
	};
	class iterator : public histograming_impl::iterator
	{
		layer1_result* res;
		upper_layer_filter_histograming* self;

		set<node*> to_layer_result;

		layer1_data* nd;
		layer1_data* top_part_nd;
		
		// counters
		node* top_part_n;  		
		int x,y;
		set<node*>::const_iterator	to_layer_iter;

		// start/stop conditions
		int minx,miny,maxx,maxy;
	public:
		iterator(upper_layer_filter_histograming* s, layer1_result* res, int minx, int miny, int maxx, int maxy) :
			self(s), res(res), minx(minx), miny(miny), maxx(maxx), maxy(maxy), top_part_n(nullptr), top_part_nd(nullptr), nd(nullptr), x(minx), y(miny-1) {}

		inline bool move_next() {
			// this should provide the same functionality as
			//for (int x = minx; x < maxx; ++x)
			//   for (int y = miny; y < maxy; ++y)
			//        node* n = node_at(x, y, layer);
			//        while (n != nullptr)
			//			  recurse05(start_node_list, atom("toPrevLayer"), layer_check_2(layer), to_layer_result);
			//			  for (set<node*>::const_iterator iter = to_layer_result.begin(); iter != to_layer_result.end(); iter++) {
			//				  count parts
			//            n = nd->next;
	
			bool has_next = true;
			if (to_layer_result.size() > 0 && ++to_layer_iter != to_layer_result.end()) {
				// move to next part in to_layer_result list
				//to_layer_iter++;
			} else {
				// move to next part at this location if any left				
				to_layer_result.clear();
				while (to_layer_result.size() <= 0 && has_next) {
					if (top_part_n != nullptr && top_part_nd->next != nullptr) {
						top_part_n = top_part_nd->next;
					} else {
						// otherwise continue moving to next location (x,y) until any part is found 
						top_part_n = nullptr;
						while (top_part_n == nullptr && has_next) {
							if (y < maxy-1) {
								y++;
							} else {
								if (x < maxx-1) {
									x++;
								} else {
									// should stop iteration
									has_next = false;
									break;
								}
								// reset y
								y = miny;
							}
							// reset n
							top_part_n = res->node_at(x, y, self->top_part_layer);
						}
					}
						// get layer1_data for next top layer part
					if (top_part_n != nullptr) {
						top_part_nd = (layer1_data*)top_part_n->data;

						// follow parts to lower layer and obtain list of parts to process next
						set<node*> start_node_list;
						start_node_list.insert(top_part_n);

						res->recurse05(start_node_list, atom("toPrevLayer"), upper_layer_filter_histograming::layer_check(self->layer), to_layer_result);

						// prepare iterator
						to_layer_iter = to_layer_result.begin();
					}
				}
			}

			// get layer1_data for next part
			if (to_layer_result.size() > 0 && *to_layer_iter != nullptr)
				nd = (layer1_data*)(*to_layer_iter)->data;			

			return has_next;
		}
		inline void reset() {
			if (!res->grid(self->top_part_layer)) res->init_grid(self->top_part_layer);
			//if (!res->grid(self->layer)) res->init_grid(self->layer);
		}
  
		inline ipoint2 location() { return ipoint2(nd->x,nd->y); }
		inline int type() { return self->merged_parts_mapping[nd->m]; }
		inline float value() { return nd->vval(); }

	};

	upper_layer_filter_histograming(int top_part_layer, int layer, int maxfeat, vector<int> merged_parts_mapping) : 
		top_part_layer(top_part_layer), layer(layer), maxfeat(maxfeat), merged_parts_mapping(merged_parts_mapping) {
		printf("using upper_layer_filter_histograming with settings:\n\ttop_part_layer: %d\n\tlayer: %d\n", top_part_layer, layer); fflush(stdout);
	}
	
	inline virtual int get_max_features() { return maxfeat; }

	inline virtual histograming_impl::iterator* get_features_iterator(layer1_result* res, const irectangle2& r) {		
		int minx = max(0, r.ll.x);
		int maxx = min(res->x_size(layer), r.ur.x);
		int miny = max(0, r.ll.y);
		int maxy = min(res->y_size(layer), r.ur.y);

		// adjust min/max x/y values to account for layer contraction
		double layer_diff_contraction = (double)res->x_size(layer) / (double)res->x_size(top_part_layer);

		minx /= layer_diff_contraction;
		maxx /= layer_diff_contraction;
		miny /= layer_diff_contraction;
		maxy /= layer_diff_contraction;

		return new iterator(this, res, minx, miny, maxx, maxy);
	}
};

class matrix_binning : public binning_impl {
	M_bin* original_bin;

	M_bin* bin;

	histogram_distance_function* dist_func;

	bool interpolate_to_neighboors;
	bool legacy_compatible;
	bool original_bin_already_resized;
public:

	class iterator : public binning_impl::iterator {
	public:
		matrix_binning* self;
		ipoint2 feat_location;

		int num_bins;
		int current_bin;
		bool processed;
		
		iterator(matrix_binning* s) : self(s), current_bin(0), processed(0), num_bins(s->bin->bin_count()) {
		}
		
		inline void reset(ipoint2 feat_location){
			// mirror feat_location over center for compatibility with previous implementation
			if (self->legacy_compatible)
				this->feat_location = -(feat_location - ipoint2(self->bin->width()/2,self->bin->height()/2)); // !!!!
			else
				this->feat_location = (feat_location - ipoint2(self->bin->width()/2,self->bin->height()/2));
			this->current_bin = -1;
			this->processed = 0;

		}
		inline bool move_next() {
			bool has_next = false;
			if (self->interpolate_to_neighboors) {
				// iterate through all bins
				if (current_bin++ < num_bins-1)
					has_next = true;
			} else {
				// use only one bin and finish
				current_bin = self->bin->get_bin(feat_location.x, feat_location.y);				
				has_next = true && processed == false;
				processed = true;
			}
			return has_next;
		}

		inline int bin_index() { return current_bin; }

		inline float bin_weight() {
			return self->dist_func != nullptr ? (*self->dist_func)(*self->bin, feat_location, current_bin) : 1;
		}
	};

	matrix_binning(M_bin* b, histogram_distance_function* d, bool i, bool l): bin(b), dist_func(d), interpolate_to_neighboors(i), legacy_compatible(l), original_bin_already_resized(false) {
		if (bin == nullptr)
			throw new_libhop_exception("Cannot create matrix_binning withut valid bin (M_bin is null)");
		// make a copy of bin that is used for resizing to proper rectangle size
		original_bin = new M_bin(*b);
		printf("using matrix_binning with settings:\n\tbin size x: %d\n\tbin size y: %d\n", original_bin->width(), original_bin->height()); fflush(stdout);
	}
	virtual ~matrix_binning() {
		if (bin != nullptr)
			delete bin;
		if (bin != nullptr)
			delete original_bin;
	}
	inline virtual int get_num_bins() {
		return bin->bin_count();
	}
	inline virtual binning_impl::iterator* get_bin_iterator(layer1_result* res, const irectangle2& r) {
		if (original_bin_already_resized == false) {
			double contraction_factor = (double)res->x_size(0) / (double)res->x_size(current_layer);
			original_bin = original_bin->get_resized(original_bin->width()/contraction_factor, original_bin->height()/contraction_factor);
			original_bin_already_resized  = true;
		}
		// resize bin to match size of region r
		int x_dim = floor(r.x_dim() + 0.5);
		int y_dim = floor(r.y_dim() + 0.5);
		if (x_dim != bin->width() || y_dim != bin->height()) {
			// delete current bin
			delete bin;
			// get resized bin
			bin = original_bin->get_resized(r.x_dim(), r.y_dim());
		} 
		return new iterator(this);
	}
};

class spm_binning : public binning_impl {

	int pyramid_layers;
	int number_bins;
public:

	class iterator : public binning_impl::iterator {
	public:
		spm_binning* self;
		ipoint2 feat_location;
	
		irectangle2 rect;

		int current_pyramid_layer;
		int current_pyramid_offset;
		int current_pyramid_bin;

		iterator(spm_binning* s, irectangle2 r) : self(s), current_pyramid_layer(-1), current_pyramid_offset(0), current_pyramid_bin(0), rect(r) {
			
		}
		
		inline void reset(ipoint2 feat_location){
			this->feat_location = feat_location;
			this->current_pyramid_layer = -1;
			this->current_pyramid_offset = 0;
			this->current_pyramid_bin = 0;
		}
		inline bool move_next() {
			bool has_next = false;
			if (current_pyramid_layer < self->pyramid_layers-1) {
				// add pyramid bin offset for previous layer
				current_pyramid_offset += pow(2.0,2*current_pyramid_layer);
				
				// move to next pyramid layer
				current_pyramid_layer++;
				
				// get number of splits for new pyramid layer
				int num_splits = pow(2.0,current_pyramid_layer);
				// get position of part within this pyramid layer
				int l_position_x = min(max((feat_location.x * num_splits) / rect.x_dim(), 0), num_splits-1);
				int l_position_y = min(max((feat_location.y * num_splits) / rect.y_dim(), 0), num_splits-1);
				
				// save index of bin
				current_pyramid_bin = l_position_x + num_splits * l_position_y;

				has_next = true;
			}
			return has_next;
		}

		inline int bin_index() { return current_pyramid_offset + current_pyramid_bin; }

		inline float bin_weight() {
			// use proper weights based on specific pyramid layer
			return current_pyramid_layer == 0 ? 1 / pow(2.0,self->pyramid_layers) : 1 / pow(2.0,self->pyramid_layers - current_pyramid_layer +1);
		}
	};

	spm_binning(int pyramid_layers): pyramid_layers(pyramid_layers)  {
		number_bins = 0;
		for (int l = 0; l < pyramid_layers; l++)
			number_bins += pow(2.0,l*2);
		printf("using spm_binning with settings:\n\tpyramid_layers: %d\n", pyramid_layers); fflush(stdout);
	}

	inline virtual int get_num_bins() {
		return number_bins;
	}
	inline virtual binning_impl::iterator* get_bin_iterator(layer1_result* res, const irectangle2& r) {
		return new iterator(this, r);
	}
};

class grid_binning : public binning_impl {

	int x_elements, y_elements;
public:

	class iterator : public binning_impl::iterator {
	public:
		grid_binning* self;
		ipoint2 feat_location;

		irectangle2 rect;

		bool processed;

		iterator(grid_binning* s, irectangle2 r) : self(s), rect(r), processed(false) {
		}
		
		inline void reset(ipoint2 feat_location){
			this->feat_location = feat_location;
			this->processed = false;
		}
		inline bool move_next() {			
			// this iterator should only iterate once (one element)
			
			// continue only if not yet processed
			bool has_next = processed == false;
			// set it to processed so next loop will stop
			processed = true;
			return has_next;
		}

		inline int bin_index() {
			int x_pos = min(max((feat_location.x * self->x_elements) / rect.x_dim(),0), self->x_elements-1);
			int y_pos = min(max((feat_location.y * self->y_elements)  / rect.y_dim(), 0), self->y_elements-1);
			return y_pos * self->x_elements + x_pos; 
		}

		inline float bin_weight() {
			return 1;
		}
	};

	grid_binning(int x_elements, int y_elements): x_elements(x_elements), y_elements(y_elements)  {
		printf("using grid_binning with settings:\n\tx_elements: %d\ty_elements: %d\n", x_elements, y_elements); fflush(stdout);
	}

	inline virtual int get_num_bins() {
		return x_elements*y_elements;
	}
	inline virtual binning_impl::iterator* get_bin_iterator(layer1_result* res, const irectangle2& r) {
		return new iterator(this, r);
	}
};

class histogram_clusters {
	vector<vector<float> > centers;
public:
	histogram_clusters(const float* values, const int width, const int height) {
		centers = vector<vector<float> >(height, vector<float>(width,0.0));
		for (int i = 0; i < height; i++) {
			vector<float>& c = centers[i];

			for (int j = 0; j < width; j++) {
				c[j] = values[i* width + j];
			}
		}
	}
	inline int size() const { return centers.size(); }

	inline int find_closest_center(vector<float> feature, float *return_best_distance = nullptr) const {
		float best_distance = FLT_MAX;
		int best_center_index = 0;
		for (int i = 0; i < centers.size(); i++) {
			const vector<float>& c = centers[i];

			int min_dimension = min(feature.size(), c.size());

			float current_distance = 0;			
			for (int j = 0; j < min_dimension; j++) {
				float vv = feature[j] - c[j];
				current_distance += vv*vv;
			}
			if (current_distance < best_distance){
				best_distance = current_distance;
				best_center_index = i;
			}
		}

		if (return_best_distance != nullptr)
			*return_best_distance = best_distance;
		return best_center_index;
	}
};

class dense_featureset_histograming : public histograming_impl {
	int layer;
	int number_featureset_centers;
	histogram_clusters* centers;;

	int stride;
	int window_size_x,window_size_y;

	histograming_impl* local_feature_hist_impl;
	binning_impl* local_feature_bin_impl;

public:	
	class iterator : public histograming_impl::iterator
	{
		layer1_result* res;
		dense_featureset_histograming* self;

		// counters
		int x,y;

		// start/stop conditions
		int minx,miny,maxx,maxy;

		float current_feature_center_distance;
		int current_feature_center;

		int adjusted_stride;
		int adjusted_window_size_x, adjusted_window_size_y;
	public:
		iterator(dense_featureset_histograming* s, layer1_result* res, int adjusted_stride, int adjusted_window_size_x, int adjusted_window_size_y, int minx, int miny, int maxx, int maxy) :
			self(s), res(res), adjusted_stride(adjusted_stride), adjusted_window_size_x(adjusted_window_size_x), adjusted_window_size_y(adjusted_window_size_y),
			minx(minx), miny(miny), maxx(maxx), maxy(maxy), x(minx), y(miny-adjusted_stride) {
		}

		inline bool move_next() {
			bool has_next = true;
			if (y < maxy-(adjusted_stride)) {
				y += adjusted_stride;
			} else {
				if (x < maxx-(adjusted_stride)) {
					x += adjusted_stride;
				} else {
					// should stop iteration
					has_next = false;
				}
				// reset y
				y = miny;
			}
			// calculate feature only if should continue
			if (has_next)
				calculate_feature();
			return has_next;
		}
		inline void reset() {
			x = minx;
			y = miny - adjusted_stride;
		}
		inline void calculate_feature() {
			// on location (x,y) with size (window_size_x,window_size_y) calculate small local HoC using plane histograming with simple grid binning
			vector<float> feature = hoc_histogram_generator::get_histogram(self->local_feature_hist_impl, self->local_feature_bin_impl, res, self->layer, irectangle2(x,y,x+adjusted_window_size_x,y+adjusted_window_size_y), false);
			
			// find closest center of this feature (assining it to closest known clusters)
			current_feature_center = self->centers->find_closest_center(feature, &current_feature_center_distance);
		}
	  
		inline ipoint2 location() { return ipoint2(x,y); }
		inline int type() { return current_feature_center; }
		inline float value() { return 1; } // maybe use current_feature_center_distance?

	};

	dense_featureset_histograming(int layer, int stride, int win_size_x, int win_size_y, histogram_clusters* centers, histograming_impl* local_feature_hist_impl, binning_impl* local_feature_bin_impl) : 
		layer(layer), stride(stride), window_size_x(win_size_x), window_size_y(win_size_y), centers(centers), number_featureset_centers(centers->size()), 
		local_feature_hist_impl(local_feature_hist_impl), local_feature_bin_impl(local_feature_bin_impl) {
		printf("using dense_featureset_histograming with settings:\n\tlayer: %d\n\tstride: %d\n\twin_size_x: %d\n\twin_size_x: %d\tnumber centers: %d\n", layer, stride, win_size_x, win_size_y, number_featureset_centers); fflush(stdout);
	}
	virtual ~dense_featureset_histograming() {
		if (local_feature_hist_impl != nullptr)
			delete local_feature_hist_impl;
		if (local_feature_bin_impl != nullptr)
			delete local_feature_bin_impl;
		if (centers != nullptr)
			delete centers;
	}
	inline virtual int get_max_features() { return number_featureset_centers; }

	inline virtual histograming_impl::iterator* get_features_iterator(layer1_result* res, const irectangle2& r) {
		int minx = max(0, r.ll.x);
		int maxx = min(res->x_size(layer), r.ur.x);
		int miny = max(0, r.ll.y);
		int maxy = min(res->y_size(layer), r.ur.y);

		// adjust windows size values to account for layer contraction
		double layer_diff_contraction = (double)res->x_size(0) / (double)res->x_size(layer);
		
		int adjusted_window_size_x = max<int>(1,this->window_size_x / layer_diff_contraction);
		int adjusted_window_size_y = max<int>(1,this->window_size_y / layer_diff_contraction);

		int adjusted_stride = max<int>(1,this->stride / layer_diff_contraction);

		return new iterator(this, res, adjusted_stride, adjusted_window_size_x, adjusted_window_size_y, minx, miny, maxx, maxy);
	}
};

}

#endif /* _HOC_1_ */
