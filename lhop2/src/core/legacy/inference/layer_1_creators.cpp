
// layer_1_creators
///////////////////////////////////////////////////////////////////////////////

//#include "utils/atom.h"
#include "utils/graphs/img_graph.h"
#include "utils/misc.h"
#include "utils/utils.h"
#include "layer_1_creators.h"

#define THREAD_TEST__
#ifdef THREAD_TEST
struct threadData
{
	char* arg1;
	char* arg2;
	char* arg3;
};
#endif

// layer1_creator
///////////////////////////////////////////////////////////////////////////////

layer1_creator::layer1_creator() : 
    filter_kernels(),
    layer1_region_threshold(0.1),
    layer1_threshold(0.1),
    layer1_3x3bound(3),
    layer1_neighb_radius(7),
    power_correction(1.0),
    normalization_percent(0.0),
    response_percent(0.8),
	use_opencl(false),
    border(0),	
	make_dummy_result(false)
{
    to_neighbor = EdgeConnection::TO_NEIGHBOOR;
}

layer1_creator::~layer1_creator()
{

    delete_filter_kernels();
}

void layer1_creator::delete_filter_kernels()
{
    for (unsigned i = 0; i < filter_kernels.size(); ++i)
        delete filter_kernels[i];
    for (unsigned i = 0; i < regions.size(); ++i)
        delete regions[i];
    filter_kernels.clear();
    regions.clear();
}

int layer1_creator::create_result_vector(vector<layer1_result*>& results, const img& srcim, const img& maskim) {
	
	int originalw = srcim.width;

	results.clear();

	if (srcim.empty()) return -1;

	img im;
    img mim;

	if (init_size != -100) srcim.get_resized(im, init_size);
	else im = srcim;

	if (!maskim.empty())
		if (init_size != -100) maskim.get_resized(mim, init_size);
		else mim = maskim;
	else
		mim = img(im.width, im.height, 1, true);

	if (scale_mask_size >= 3 && scale_sigma > 0.0) 
		im.blur(scale_mask_size, scale_sigma);

	int swidth = (int)im.width;
	int sheight = (int)im.height;

	int count = 0;
    bool process = true;

	while (process) {
		layer1_result* result = create_result(im, mim, originalw);

		results.push_back(result);
		result = nullptr;

		swidth = (int)(swidth*scale_factor);
		sheight = (int)(sheight*scale_factor);
        ++count;
			
		if (swidth < scale_limit && sheight < scale_limit || count >= max_scales) 
            process = false;
        else {
            im = im.get_resized(swidth, sheight);
            if (!mim.empty()) mim = mim.get_resized(swidth, sheight);
		}
	}
}
int layer1_creator::create_result_vector(vector<layer1_result*>& results, const img& srcim, const vector<irectangle2>& maskregions)
{
    int originalw = srcim.width;
	img maskim;
	
	// generate layer1 from list of bounding box regions where each region indicates possible object location at specific scale

	// make a histogram of all scales
	int num_scale_clusters = 40;
	int trained_model_size = 150;

	float max_scale = 0;
	float min_scale = FLT_MAX;
	vector<float> scales(maskregions.size(),0);

	for (int i = 0; i < maskregions.size(); i++) {
		float scale = trained_model_size/(float)min<int>(maskregions[i].x_dim(),maskregions[i].y_dim());
		max_scale = max<float>(max_scale, scale);
		min_scale = min<float>(min_scale, scale);
		scales[i] = scale;
	}

	// split scale range into num_scale_clusters of different clusters
	vector<float> scale_histogram_centers(num_scale_clusters,0);

	float scale_step = (max_scale - min_scale) / (float)num_scale_clusters;
	for (int i = 0; i < num_scale_clusters; i++) {
		scale_histogram_centers[i] = min_scale + i * scale_step;
	}

	// assign each region to closes histogram scale center
	vector<list<irectangle2> > scale_clusters(num_scale_clusters,list<irectangle2>());

	for (int i = 0; i < maskregions.size(); i++) {
		const irectangle2& reg = maskregions[i];
		float scale = scales[i];

		// find closes scale center
		float min_dist = FLT_MAX;
		int min_dist_index = -1;

		for (int j = 0; j < scale_histogram_centers.size(); j++) {
			float dist = scale - scale_histogram_centers[j];
			dist *= dist;
			if (dist < min_dist)  {
				min_dist = dist;
				min_dist_index = j;
			}
		}

		// add this region into cluster of closes scale
		scale_clusters[min_dist_index].push_back(maskregions[i]);
	}

	// for each scale find its assosiated regions
	for (int i = num_scale_clusters-1; i >= 0; i--) {
		float scale = scale_histogram_centers[i];
		list<irectangle2>& scale_regions = scale_clusters[i];

		// skip scales with no regions
		if (scale_regions.size() <= 0)
			continue;

		if (scale * srcim.width > 5000 || scale * srcim.height > 5000 ) {
			cout << "Image size bigger then 5000px: " << scale * srcim.width << " width and " << scale * srcim.height << " height " << endl;
			continue;
		}

		img im;
		img mim;

		img maskim(srcim.width, srcim.height, 0, true);
		// generate image mask based on all regions for this scale
		for (list<irectangle2>::const_iterator iter = scale_regions.begin(); iter != scale_regions.end(); iter++) {
			const irectangle2& reg = *iter;
			maskim.set_region(reg.ll.x, reg.ur.x, reg.ll.y, reg.ur.y , 1);
		}

		// prepare image for this scale
		srcim.get_resized(im, (int)(-scale*100));
		maskim.get_resized(mim, (int)(-scale*100));

		// now process image on this scale
		layer1_result* result = create_result(im, mim, originalw);
            

		results.push_back(result);
		result = nullptr;
	}
	return 0;	
}

layer1_result* layer1_creator::create_result(img& im, img& maskim, int originalw) {
	int swidth = (int)im.width;
	int sheight = (int)im.height;

	// convert to grayscale since all images will be in colors
	im.to_grayscale();

	layer1_result* result = init_result();
	
	result->new_grid(swidth + 2*border, sheight + 2*border, 0);

	if (make_dummy_result == true) {
		result->shape_nodes.resize(1);
		result->shape_nodes_inhib.resize(1);
	} else {

        vector<img*> images;
        img* max_image = nullptr;
        iimg* max_image_src = nullptr;

		make_images(images, im, maskim);

		make_max_image(max_image, max_image_src, images);

		make_result(result, max_image, max_image_src, maskim);

		finalize(result, images, max_image, max_image_src);

		inhibit_result(result);			

		for (int i = 0; i < (int)images.size(); ++i)
            if (images[i] != nullptr) delete images[i];
        if (max_image) delete max_image;
        if (max_image_src) delete max_image_src;
	}
	result->border = border;
	result->original_width = originalw;

	return result;
}

void layer1_creator::image_processing(img* im)
{
    
}

bool layer1_creator::init_default_from_result(layer1_result* res)
{
    bool changed = false;

    if (to_neighbor != res->to_neighbor) { changed = true; to_neighbor = res->to_neighbor; }
    if (layer1_region_threshold != res->layer1_region_threshold) 
        { changed = true; layer1_region_threshold = res->layer1_region_threshold; }
    if (layer1_threshold != res->layer1_threshold)
        { changed = true; layer1_threshold = res->layer1_threshold; }
    if (layer1_3x3bound != res->layer1_3x3bound) { changed = true; layer1_3x3bound = res->layer1_3x3bound; }
    if (layer1_neighb_radius != res->layer1_neighb_radius)
        { changed = true; layer1_neighb_radius = res->layer1_neighb_radius; }
    return changed;
}

void layer1_creator::init_result_with_default(layer1_result* res)
{
    res->to_neighbor = to_neighbor;
    res->layer1_region_threshold = layer1_region_threshold;
    res->layer1_threshold = layer1_threshold;
    res->layer1_3x3bound = layer1_3x3bound;
    res->layer1_neighb_radius = layer1_neighb_radius;
}


layer1_result* layer1_creator::init_result()
{
    return nullptr;
}

void layer1_creator::make_filter_kernels()
{

}

void layer1_creator::make_images(vector<img*>& images, img& im, img& immask)
{

}

int layer1_creator::make_images_get_count() 
{
	return 0;
}

void layer1_creator::make_max_image(img*& max_image, iimg*& max_image_src, const vector<img*>& images)
{
    int image_count = (int)images.size();

    if (image_count == 0) {
        max_image = nullptr;
        return;
    }

    HOP_REAL** ptrs = new HOP_REAL*[image_count];
    int w = (int)images[0]->width;
    int h = (int)images[0]->height;

    for (int i = 0; i < image_count; ++i) {
        ptrs[i] = images[i]->ptr(0, 0);
    }

    max_image = new img(w, h);
    max_image_src = new iimg(w, h);
    HOP_REAL* ptr = max_image->ptr(0, 0);
    int* iptr = max_image_src->ptr(0, 0);
    HOP_REAL tmpmax;
    int tmpi;

    for (int j = 0; j < h; ++j) {
        for (int i = 0; i < w; ++i) {
            tmpmax = 1.0e-10; 
            tmpi = 0;
            for (int k = 0; k < image_count; ++k) {
                if (*ptrs[k] > tmpmax) { tmpi = k; tmpmax = *ptrs[k]; }
                ++(ptrs[k]);
            } 
            *(ptr++) = tmpmax;
            *(iptr++) = tmpi;
        }
    }
    delete[] ptrs;
}

void layer1_creator::make_result(layer1_result* result, img* max_image, iimg* max_image_src, const img& maskimg)
{
    if (result->shape_nodes.size() == 0) {
        result->shape_nodes.resize(1);
        result->shape_nodes_inhib.resize(1);
    }

    int w = (int)max_image->width;
    int h = (int)max_image->height;
    double wfactor = maskimg.width == 0 ? 0.0 : (double)maskimg.width/w;
    double hfactor = maskimg.height == 0 ? 0.0 : (double)maskimg.height/h;
    bool mask = maskimg.width > 0 && maskimg.height > 0;

    if (w < 3 || h < 3) return;

    img* sum2x2;
    HOP_REAL sum2x2max = max_image->p_get_2x2sum(sum2x2);

    HOP_REAL* ptrs[9];
    HOP_REAL* ptrs2[9];
    HOP_REAL threshold = max_image->maximum()*layer1_threshold;
    HOP_REAL threshold2 = sum2x2max*layer1_threshold;
    
    //cout << endl << threshold << endl << layer1_threshold << endl;

    ptrs[0] = max_image->ptr(0, 0); ptrs2[0] = sum2x2->ptr(0, 0);
    ptrs[1] = max_image->ptr(1, 0); ptrs2[1] = sum2x2->ptr(1, 0);
    ptrs[2] = max_image->ptr(2, 0); ptrs2[2] = sum2x2->ptr(2, 0);
    ptrs[3] = max_image->ptr(0, 1); ptrs2[3] = sum2x2->ptr(0, 1);
    ptrs[4] = max_image->ptr(1, 1); ptrs2[4] = sum2x2->ptr(1, 1);
    ptrs[5] = max_image->ptr(2, 1); ptrs2[5] = sum2x2->ptr(2, 1);
    ptrs[6] = max_image->ptr(0, 2); ptrs2[6] = sum2x2->ptr(0, 2);
    ptrs[7] = max_image->ptr(1, 2); ptrs2[7] = sum2x2->ptr(1, 2);
    ptrs[8] = max_image->ptr(2, 2); ptrs2[8] = sum2x2->ptr(2, 2);

    HOP_REAL** mid = &ptrs[4];
    HOP_REAL** mid2 = &ptrs2[4];;
    int* iptr = max_image_src->ptr(1, 1);
    int gcount, gcount2;
    int k;
    node* n;
    vector<node*>& s_nodes = result->shape_nodes[0];

    for (int j = 1; j < h - 1; ++j) {
        for (int i = 1; i < w - 1; ++i) {
            if (!mask || maskimg.at(int_round(i*wfactor), int_round(j*hfactor), 0.0) > 1E-6) {
                gcount = gcount2 = 0;
                if (**mid >= threshold) {
                    for (k = 0; k < 9; ++k) { if (*ptrs[k] >= **mid) ++gcount; }
                    if (gcount <= layer1_3x3bound) {
                        n = result->add_grid_node(new layer1_data(**mid, *iptr), i + border, j + border);
                        s_nodes.push_back(n);

                    } else {
                        for (k = 0; k < 9; ++k) { if (*ptrs2[k] >= **mid2) ++gcount2; }
                        if (gcount2 <= layer1_3x3bound) {
                            n = result->add_grid_node(new layer1_data(**mid, *iptr), i + border, j + border);
                            s_nodes.push_back(n);
                        }
                    }
                }
            }
            for (k = 0; k < 9; ++k) { ++(ptrs[k]); ++(ptrs2[k]); }
            ++iptr;
        }
        for (int k = 0; k < 9; ++k) { ptrs[k] += 2; ptrs2[k] += 2; }
        iptr += 2;
    }
    if (sum2x2 != nullptr) delete sum2x2;            
}


void layer1_creator::inhibit_result(layer1_result* result)
{
    //! we could do (approximately) the same thing in the process of generating max_image?

    result->inhibit(0);
}

// Does the following:
// - normalizes shape nodes
// - adds all filter responses to nodes.
// - adds in/out flags if imb != nullptr
void layer1_creator::finalize(layer1_result* result, vector<img*>& images, img* max_image, iimg* max_image_src)
{

    vector<node*>& s_nodes = result->shape_nodes[0];
    vector<node*>::iterator iter;

    if (s_nodes.empty()) 
        return;

    sort(s_nodes.begin(), s_nodes.end(), response_sort_f(R_RESPONSE)); 

    double max;

    // normalize shape_nodes
    if (normalization_percent <= 0.0 || normalization_percent >= 1.0) {
        max = ((layer1_data*)(s_nodes.front()->data))->r(R_RESPONSE);
    } else {
        int i = (int)(s_nodes.size() * normalization_percent);

        if (i >= (int)s_nodes.size()) i = (int)s_nodes.size() - 1;
        max = ((layer1_data*)s_nodes[i]->data)->r(R_RESPONSE);
    }

    if (max <= 0.0) return; // !?

    // add all filter responses & normalize
    bool powercorr = (power_correction != 1.0);
    int image_count = (int)images.size();
	
    for (iter = s_nodes.begin(); iter != s_nodes.end(); ++iter) {
        layer1_data* nd = (layer1_data*)(*iter)->data;
        int x = nd->x - border, y = nd->y - border;
        int prohib = max_image_src->at(x, y);
        double origd = nd->r(R_RESPONSE);
        double d = origd;

        if (normalize) d /= max;
        if (d > 1.0) d = 1.0;
        if (powercorr) d = pow(d, power_correction);
        nd->r.set_response(R_RESPONSE, d);
		nd->r.set_response(G_RESPONSE, 1.0);
		nd->r.set_response(RR_RESPONSE, 1.0);
        nd->r.set_response(S_RESPONSE, 0.0);
		node* first_node = *iter;

        for (int i = 0; i < image_count; ++i) {
            if (i != prohib) {
                double val = images[i]->at(x, y);

                if (val >= origd*response_percent) {
                    if (normalize) val /= max;
                    if (val > 1.0) val = 1.0;
                    if (powercorr) val = pow(val, power_correction);

                    layer1_data* nnd = new layer1_data(val, i);

		            nnd->r.set_response(G_RESPONSE, 1.0);
		            nnd->r.set_response(RR_RESPONSE, 1.0);
                    nnd->r.set_response(S_RESPONSE, 0.0);
                    node* nn = result->add_grid_node(nnd, nd->x, nd->y, 0);					
					if (nn->is_attr_set(IMG_NODE_ATTR)) {
						first_node = nn;
					}
                }
            }
        }
		// if new node (nn) gets added as first in line then it must replace current node in s_nodes
		if (first_node != *iter)
			*iter = first_node;
    }

}

unsigned layer1_creator::lib_type() { return STRUCT_LIB_ATTR; }

part_lib* layer1_creator::get_part_lib(const ConfigDictionary& cfg)
{
    vector<part_data*> pd;
    
    get_part_data(pd, cfg);
    part_lib* result = new part_lib(lib_type());
    for (size_t i = 0; i < pd.size(); ++i) {
        result->add_part(1, pd[i], vector<iipair>(), vector<iipair>(), vector<matrix<double>*>());
    }
    return result;
}

void layer1_creator::cfg_init(const ConfigDictionary& cfg)
{
    cfg.getValue(layer1_region_threshold, "layer1_region_threshold");
    cfg.getValue(layer1_threshold, "layer1_threshold");
    cfg.getValue(power_correction, "power_correction");
    cfg.getValue(normalization_percent, "normalization_percent");
    cfg.getValue(response_percent, "response_percent");
    cfg.getValue(layer1_3x3bound, "layer1_3x3bound");
    cfg.getValue(layer1_neighb_radius, "layer1_neighb_radius");  
    cfg.getValue(border, "border_size");
    normalize = cfg.getValueBool("normalize", true);
	scale_limit = cfg.getValueInt("scale_limit", 200);
    max_scales = cfg.getValueInt("max_scales", 99);
    init_size = cfg.getValueInt("init_size", 0);
    scale_factor = cfg.getValueDouble("scale_factor", 1/::pow(2.0, 1.0/3.0));
    scale_sigma = cfg.getValueDouble("scale_sigma", 0.0);
    scale_mask_size = (int)(5.0 * scale_sigma);

	make_dummy_result = cfg.getValueBool("make_dummy_result", false);	

	use_opencl = cfg.getValueBool("use_opencl", false);	

	if (use_opencl) {
		cout << "Warning: Found 'use_opencl = true' in config but binary not build with OpenCL support. Will use original implementation !!" << endl;
		use_opencl = false;
	}
}


void layer1_creator::init_filters()
{
	make_filter_kernels();
}

void layer1_creator::clear_filters()
{

}


