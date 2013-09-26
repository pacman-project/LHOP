
// layer_1_creators
///////////////////////////////////////////////////////////////////////////////

#include "../utils/atom.h"
#include "../graphs/img_graph.h"
#include "../utils/misc.h"
#include "../utils/utils.h"
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
    //images(), 
    filter_kernels(),
    //result(0),
    layer1_region_threshold(0.1),
    layer1_threshold(0.1),
    layer1_3x3bound(3),
    layer1_neighb_radius(7),
    power_correction(1.0),
    normalization_percent(0.0),
    response_percent(0.8),
    //max_image(0),
    //max_image_src(0),
    border(0),
	use_opencl(false),
	opencl_verify_result(false),
	make_dummy_result(false)
{
    to_neighbor = atom("toNeighbor").get_index();
}

layer1_creator::~layer1_creator()
{

    delete_filter_kernels();
//    delete_images();
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

//void layer1_creator::delete_images()
//{
//    for (unsigned i = 0; i < images.size(); ++i)
//        delete images[i];
//    images.clear();
//    if (max_image) { delete max_image; max_image = nullptr; }
//    if (max_image_src) { delete max_image_src; max_image_src = nullptr; }
//}


//int layer1_creator::create_result(vector<layer1_result*>& results, bimg* src, img* maskimg)
//{
//	results.clear();
//	if(src == nullptr) return -1;
//	img sim(*src);
//	create_result(results, &sim, maskimg);
//}


//int layer1_creator::create_result(vector<layer1_result*>& results, img* im, const img* maskimg)
//{
//    img im(srcim);
//
//    return create_result(results, &im, &maskimg);
//}

int layer1_creator::create_result_vector(vector<layer1_result*>& results, const img& srcim, const img& maskim) {
	
	int originalw = srcim.width;

	results.clear();

	if (srcim.empty()) return -1;

	img im;
    img mim;

//	clock_t startt,endt;

//	startt = clock();
	if (init_size != -100) srcim.get_resized(im, init_size);
	else im = srcim;
//	endt = clock(); cout << "Resize img took " << ((double)endt-(double)startt)/CLOCKS_PER_SEC << " sec" << endl;
//	startt = clock();
    if (!maskim.empty())
		if (init_size != -100) maskim.get_resized(mim, init_size);
		else mim = maskim;
	else
		mim = img(im.width, im.height, 1, true);
//	endt = clock(); cout << "Resize mask took " << ((double)endt-(double)startt)/CLOCKS_PER_SEC << " sec" << endl;

    //cerr << "scale_mas_size = " << scale_mask_size << " ";
//	startt = clock();
	if (scale_mask_size >= 3 && scale_sigma > 0.0) 
		im.blur(scale_mask_size, scale_sigma);
//	endt = clock(); cout << "Blur took " << ((double)endt-(double)startt)/CLOCKS_PER_SEC << " sec" << endl;
	int swidth = (int)im.width;
	int sheight = (int)im.height;

    //cerr << "swidth = " << swidth << " sheight = " << sheight << " ";
		
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
    //cerr << srcim.width << 'x' << srcim.height << " border = " << border << " ";
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
		
		//maskim.save("c:\\work\\libhop_vs10\\res\\apple-test\\mask.png");
		//im.save("c:\\work\\libhop_vs10\\res\\apple-test\\im_resized.png");

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
	
//	clock_t startt,endt;
//	startt = clock();
	result->new_grid(swidth + 2*border, sheight + 2*border, 0);
//	endt = clock(); cout << "New grid took " << ((double)endt-(double)startt)/CLOCKS_PER_SEC << " sec" << endl;

	if (make_dummy_result == true) {
		result->shape_nodes.resize(1);
		result->shape_nodes_inhib.resize(1);
	} else {

		if (use_opencl) {
#ifdef OPENCL
			ocl_create_result(result, &im, &maskim);
#else
			cout << "Binary not build with OpenCL support !!" << endl;
			throw new std::exception();
#endif
		} else {
            vector<img*> images;
            img* max_image = nullptr;
            iimg* max_image_src = nullptr;

//			startt = clock();
			make_images(images, im, maskim);
//			endt = clock(); cout << "Make images took " << ((double)endt-(double)startt)/CLOCKS_PER_SEC << " sec" << endl;
//			startt = clock();
			make_max_image(max_image, max_image_src, images);
//			endt = clock(); cout << "Make max images took " << ((double)endt-(double)startt)/CLOCKS_PER_SEC << " sec" << endl;

//			startt = clock();
			make_result(result, max_image, max_image_src, maskim);
//			endt = clock(); cout << "Make result took " << ((double)endt-(double)startt)/CLOCKS_PER_SEC << " sec" << endl;

//			startt = clock();
			finalize(result, images, max_image, max_image_src);
//			endt = clock(); cout << "Finilaize took " << ((double)endt-(double)startt)/CLOCKS_PER_SEC << " sec" << endl;

//			startt = clock();
			inhibit_result(result);			
//			endt = clock(); cout << "Inhibition took " << ((double)endt-(double)startt)/CLOCKS_PER_SEC << " sec" << endl;

			for (int i = 0; i < (int)images.size(); ++i)
                if (images[i] != nullptr) delete images[i];
            if (max_image) delete max_image;
            if (max_image_src) delete max_image_src;
		}
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
    /*int max = (int)result->shape_nodes_inhib.size();

    for (int i = 0; i < max; ++i) {
        result->connect_neighbors_circular(result->shape_nodes_inhib[i], 
            layer1_neighb_radius, layer1_neighb_radius, NODE_REDUCED_ATTR, to_neighbor);
    }
    
    node_neighbor_iterator(iter);
    node* gn;
    layer1_data* d;
    for (int i = 0; i < max; ++i) {
        foreach_neighbor (result->shape_nodes_inhib[i], to_neighbor, iter) {
            d = (layer1_data*)(iter->second)->data;
            gn = (*result)(d->x, d->y);
            
        }
        
    }*/

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

    // D E B U G
    /*max_image->save_mathematica("c:\\temp\\maximage.m");
	max_image->save_matlab("c:\\work\\data\\test\\maximum.m");
    max_image_src->save_mathematica("c:\\temp\\maximagesrc.m");
	max_image_src->save_matlab("c:\\work\\data\\test\\origin.m");
    for (int i = 0; i < image_count; ++i) {
        string name;
        name = string("c:\\temp\\image") + i + string(".m");
        images[i]->save_mathematica(name.c_str());
    }
	*/
    // D E B U G

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
            //cout << i << ' ' << pow(images[i]->at(x, y)/max, power_correction) << ' ' << prohib << ' ' << x << ',' << y << ' ' << d << endl;
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

#ifdef DEBUG_OUTPUT
	std::ofstream osx("x.mat");
	std::ofstream osy("y.mat");
	std::ofstream origin("origin.mat");
	for (iter = s_nodes.begin(); iter != s_nodes.end(); ++iter)
	{
		layer1_data* nd = (layer1_data*)(*iter)->data;
		if(iter != s_nodes.begin())
		{
			osx << ",";
			osy << ",";
			origin << ",";
		}
		osx << nd->x;
		osy << nd->y;
		origin << nd->m;
	}
	osx << endl;
	osy << endl;
	origin << endl;
	osx.close();
	osy.close();
	origin.close();

	std::ofstream osvalue("value.mat");
	double* tmp = (double*)malloc(sizeof(double) * s_nodes.size() * image_count);
	memset(tmp, 0xBB, sizeof(double) * s_nodes.size() * image_count);
	double* p = tmp;
	for(iter = s_nodes.begin(); iter != s_nodes.end(); iter++)
	{
		layer1_data* nd = (layer1_data*)(*iter)->data;
		while(1)
		{
			*(p + s_nodes.size() * nd->m) = nd->val;
			if(nd->next == nullptr) break;
			else nd = (layer1_data*)nd->next->data;
		}
		p++;
	}
	p = tmp;
	for(int i=0; i<image_count; i++)
	{
		for(int j=0; j<(int)s_nodes.size(); j++)
		{
			if(j!=0) osvalue << ",";
			osvalue << *p;
			p++;
		}
		osvalue << endl;
	}
	osvalue.close();
	free(tmp);

	vector<node*>& s_nodes_inhib = result->shape_nodes_inhib[0];
	std::ofstream osinhib("isInhib.mat");
	for(iter = s_nodes_inhib.begin(); iter != s_nodes_inhib.end(); iter++)
	{
		layer1_data* nd = (layer1_data*)(*iter)->data;
		if(iter != s_nodes_inhib.begin()) osinhib << ",";
		osinhib << nd->x;
	}
	osinhib << endl;
	for(iter = s_nodes_inhib.begin(); iter != s_nodes_inhib.end(); iter++)
	{
		layer1_data* nd = (layer1_data*)(*iter)->data;
		if(iter != s_nodes_inhib.begin()) osinhib << ",";
		osinhib << nd->y;
	}
	osinhib << endl;
	osinhib.close();

#endif
}

unsigned layer1_creator::lib_type() { return STRUCT_LIB_ATTR; }

part_lib* layer1_creator::get_part_lib(const config_dictionary& cfg)
{
    vector<part_data*> pd;
    
    get_part_data(pd, cfg);
    part_lib* result = new part_lib(lib_type());
    for (size_t i = 0; i < pd.size(); ++i) {
        result->add_part(1, pd[i], vector<iipair>(), vector<iipair>(), vector<matrix<double>*>());
    }
    return result;
}

void layer1_creator::cfg_init(const config_dictionary& cfg)
{
    cfg.get_value(layer1_region_threshold, "layer1_region_threshold");
    cfg.get_value(layer1_threshold, "layer1_threshold");
    cfg.get_value(power_correction, "power_correction");
    cfg.get_value(normalization_percent, "normalization_percent");
    cfg.get_value(response_percent, "response_percent");
    cfg.get_value(layer1_3x3bound, "layer1_3x3bound");
    cfg.get_value(layer1_neighb_radius, "layer1_neighb_radius");  
    cfg.get_value(border, "border_size");
    normalize = cfg.get_value_bool("normalize", true);
	scale_limit = cfg.get_value_int("scale_limit", 200);
    max_scales = cfg.get_value_int("max_scales", 99);
    init_size = cfg.get_value_int("init_size", 0);
    scale_factor = cfg.get_value_double("scale_factor", 1/::pow(2.0, 1.0/3.0));
    //scale_mask_size = cfg.get_value_int("scale_mask_size", 1);
    scale_sigma = cfg.get_value_double("scale_sigma", 0.0);
    scale_mask_size = (int)(5.0 * scale_sigma);

	make_dummy_result = cfg.get_value_bool("make_dummy_result", false);	

	use_opencl = cfg.get_value_bool("use_opencl", false);	
	opencl_verify_result = cfg.get_value_bool("opencl_verify_result", false);
	
	cfg.get_value(use_opencl_devices, "use_opencl_devices");

#ifdef OPENCL

	if (use_opencl == true && OpenCL::context_manager->is_initialized == false) {
		cout << "Warning: Found 'use_opencl = true' in config but OpenCL Context IS NOT initialized (maybe missing 'opencl_enable = true' in config)." << endl << "Continuing without OpenCL acceleration." << endl;
		use_opencl = false;
	}

#else
	if (use_opencl) {
		cout << "Warning: Found 'use_opencl = true' in config but binary not build with OpenCL support. Will use original implementation !!" << endl;
		use_opencl = false;
	}
#endif
}


void layer1_creator::init_filters()
{
	make_filter_kernels();
}

void layer1_creator::clear_filters()
{

}



// layer1_creator_struct
///////////////////////////////////////////////////////////////////////////////

layer1_creator_struct::layer1_creator_struct() :
    layer1_creator()
{
    gabor_lambda = 6.0;
    gabor_gamma = 0.75;
    gabor_bw = 1.75;
    gabor_step = 30;
}

layer1_creator_struct::~layer1_creator_struct()
{ }

void layer1_creator_struct::init_from_result(layer1_result_struct* res)
{
    bool changed;

    changed = init_default_from_result(res);

    if (gabor_lambda != res->gabor_lambda) { changed = true; gabor_lambda = res->gabor_lambda; }
    if (gabor_gamma != res->gabor_gamma) { changed = true; gabor_gamma = res->gabor_gamma; }
    if (gabor_bw != res->gabor_bw) { changed = true; gabor_bw = res->gabor_bw; }
    if (gabor_step != res->gabor_step) { changed = true; gabor_step = res->gabor_step; }

    if (changed) delete_filter_kernels();
}

layer1_result* layer1_creator_struct::init_result()
{
    layer1_result_struct* result1 = new layer1_result_struct();
    
    init_result_with_default(result1);

    result1->gabor_lambda = gabor_lambda;
    result1->gabor_gamma = gabor_gamma;
    result1->gabor_bw = gabor_bw;
    result1->gabor_step = gabor_step;
    return result1;
}

/*void layer1_creator_struct::make_masks()
{
    img* mask;
    int nmasks2 = 180/gabor_step;

    masks.resize(2*nmasks2);
    regions.resize(nmasks2);
    for (int angle = 0, i = 0; angle < 180; angle += gabor_step, ++i) {
        masks[i] = mask = img::gabor_mask(gabor_lambda, angle, 0.0, gabor_gamma, gabor_bw);
        regions[i] = mask->get_01_img(layer1_region_threshold);
        masks[nmasks2 + i] = img::gabor_mask(gabor_lambda, angle, 90.0, gabor_gamma, gabor_bw);
    }

}*/

void layer1_creator_struct::make_filter_kernels()
{
    get_filter_kernels(filter_kernels, gabor_lambda, gabor_gamma, gabor_bw, gabor_step);
}

void layer1_creator_struct::get_filter_kernels(vector<img*>& filter_kernels, 
    double gabor_lambda, double gabor_gamma, double gabor_bw, int gabor_step)
{
    img* filter_kernel;
    int nfilter_kernels2 = 180/gabor_step;
	
    filter_kernels.resize(2*nfilter_kernels2);
    for (int angle = 0, i = 0; i < nfilter_kernels2; angle += gabor_step, ++i) {
        filter_kernels[i] = filter_kernel = img::gabor_mask(gabor_lambda, angle, 0.0, gabor_gamma, gabor_bw);
        filter_kernels[nfilter_kernels2 + i] = img::gabor_mask(gabor_lambda, angle, 90.0, gabor_gamma, gabor_bw);
    }

}

void layer1_creator_struct::get_filter_kernels(vector<img*>& filter_kernels)
{
    get_filter_kernels(filter_kernels, gabor_lambda, gabor_gamma, gabor_bw, gabor_step);
}

void layer1_creator_struct::get_part_data(vector<part_data*>& pd, const config_dictionary& cfg)
{
    vector<img*> filter_kernels;
    vector<img*> r_filter_kernels;
    matrix<bool>* region;

    get_filter_kernels(filter_kernels, gabor_lambda, gabor_gamma, gabor_bw, gabor_step);
    get_filter_kernels(r_filter_kernels, 
        cfg.get_value_double("region_lambda", 7.0),
        cfg.get_value_double("region_gamma", 1.0),
        cfg.get_value_double("region_bw", 2.75),
        gabor_step);

    size_t nfilter_kernels2 = filter_kernels.size()/2;

    pd.clear();
    for (size_t i = 0; i < nfilter_kernels2; ++i) {
        region = r_filter_kernels[i]->get_bool_matrix(layer1_region_threshold);
        pd.push_back(new part_data(*filter_kernels[i], *region, 0, (int)i));
        delete region;
    }
    for (size_t i = 0; i < filter_kernels.size(); ++i) {
        delete filter_kernels[i];
        delete r_filter_kernels[i];
    }
}

void layer1_creator_struct::make_images(vector<img*>& images, img& im, img& immask)
{
    int nmasks2 = (int)filter_kernels.size()/2;

    images.resize(nmasks2);
	
	#pragma omp parallel for
    for (int i = 0; i < nmasks2; ++i) {
        img* mask;
        img* cosim, * sinim;

		clock_t startt, endt;
        mask = filter_kernels[i];
        //cosim = im.get_convolve(*mask);
		cosim = im.get_convolve(*mask, immask);
        cosim->sqr();
        mask = filter_kernels[nmasks2 + i];
        //sinim = im.get_convolve(*mask);
		sinim = im.get_convolve(*mask, immask);
        sinim->sqr();
        *cosim += *sinim;
        cosim->sqrt(immask);
        images[i] = cosim;
        delete sinim;
    }
}

int layer1_creator_struct::make_images_get_count() {

	// layer1_creator_struct creates one image from two convolution filters
	return filter_kernels.size()/2;
}

void layer1_creator_struct::cfg_init(const config_dictionary& cfg)
{
    layer1_creator::cfg_init(cfg);
    cfg.get_value(gabor_lambda, "gabor_lambda");
    cfg.get_value(gabor_gamma, "gabor_gamma");
    cfg.get_value(gabor_bw, "gabor_bw");                 
    cfg.get_value(gabor_step, "gabor_step");
}


// layer1_creator_color
///////////////////////////////////////////////////////////////////////////////


// gsign: +1, -1, 0 (positive, negative, normal)
img gabor_mask_rot(int size, double lambda, double sigma, int gsign)
{
	double size2 = (double)size/2;
	img result(size, size);
	double mean = 0.0;

/*If[Sqrt[x*x + y*y] <= rfCount/2,
 e = Exp[-(x*x + y*y)/(2*sigma*sigma)];
 e = e*Cos[2*Pi/lambda Sqrt[x*x + y*y]],
 e = 0*/

	for (int j = 0; j < size; ++j) {
		for (int i = 0; i < size; ++i) {
			int x = j - size/2;
			int y = i - size/2;
			int xxyy = x*x + y*y;
			double e = 0.0;

			if (sqrt((double)xxyy) <= size2) {
				e = exp(-xxyy/(2*sigma*sigma));
				e = e*cos(2*sqrt((double)xxyy)*M_PI/lambda);
			}
			result(j, i) = e;
			mean += e;
		}
	}
	mean = mean/(size*size);

	double std = 0.0;

	for (auto iter = result.begin(); iter != result.end(); ++iter) {
		*iter -= mean;
		std += (*iter)*(*iter);
	}
	std = sqrt(std/(result.size() - 1));

	if (std > 0)
	for (auto iter = result.begin(); iter != result.end(); ++iter) {
		*iter /= std;
	}

	switch (gsign) {
		case 1:
			for (auto iter = result.begin(); iter != result.end(); ++iter) {
				if (*iter < 0) *iter = 0;
			}
			break;
		case -1:
			for (auto iter = result.begin(); iter != result.end(); ++iter) {
				if (*iter > 0) *iter = 0;
			}
			break;
	}

	return result;
}

img gabor_mask_rot(int sizef, int gsign)
{
  double lambda = 2*sizef/3.8;
  double sigma = lambda*0.8;

  return gabor_mask_rot(sizef, lambda, sigma, gsign);
}

img gauss_filter1(double sigma, bool column = false)
{
	int min = (int)(-3.0*sigma - 0.5);
	int max = (int)(3.0*sigma + 0.5);
	img result(column ? 1 : (max - min + 1), column ? (max - min + 1) : 1);
	double sum = 0.0;
	
	for (int x = min, i = 0; x <= max; ++x, ++i) {
		sum += result[i] = exp(-x*x/(2*sigma*sigma));
	}
	for (auto iter = result.begin(); iter != result.end(); ++iter) 
		*iter /= sum;
	return result;
}

img gaussdx_filter1(double sigma, bool column = false)
{
	int min = (int)(-3.0*sigma - 0.5);
	int max = (int)(3.0*sigma + 0.5);
	img result(column ? 1 : (max - min + 1), column ? (max - min + 1) : 1);
	double sum = 0.0;
	double sigma2 = sigma*sigma;
	double norm = sqrt(2*M_PI)*sigma2*sigma;
	
	for (int x = min, i = 0; x <= max; ++x, ++i) {
		double val = -2*x*exp(-x*x/(2*sigma2));

		val /= norm;
		result[i] = val;
		sum += fabs(val);
	}
	sum /= 2.0;
	for (auto iter = result.begin(); iter != result.end(); ++iter) 
		*iter /= sum;
	return result;
}

// Return pair (Dx, Dy)
pair<img*, img*> gauss_deriv(img& im, double sigma)
{
	img G = gauss_filter1(sigma);
	img Gt = gauss_filter1(sigma, true);
	img D = gaussdx_filter1(sigma);
	img Dt = gaussdx_filter1(sigma, true);

	img* img1 = im.get_convolve(D);
	img* imgDx = img1->get_convolve(Gt);
	img* img2 = im.get_convolve(G);
	img* imgDy = img2->get_convolve(Dt);
	
	delete img1;
	delete img2;
	return pair<img*, img*>(imgDx, imgDy);
}

img* nonmax_sup(img& imgMag, img& imgDir)
{
	static int offx[] = {-1, -1, 0, 1, 1, 1, 0, -1, -1};
	static int offy[] = {0, -1, -1, -1, 0, 1, 1, 1, 0};

	int h = imgMag.height;
	int w = imgMag.width;

	img* imgMax = new img(w, h);
	for (int y = 0; y < h; ++y) {
		for (int x = 0; x < w; ++x) {
			double dir = imgDir(x, y);
			int idx = int_round(((dir + M_PI)/M_PI)*4);
		    int y1 = y + offy[idx], x1 = x + offx[idx];
			int y2 = y - offy[idx], x2 = x - offx[idx];

			x1 = max(0, x1); x1 = min(w - 1, x1);
			y1 = max(0, y1); y1 = min(h - 1, y1);
		    x2 = max(0, x2); x2 = min(w - 1, x2);
			y2 = max(0, y2); y2 = min(h - 1, y2);
			if (imgMag(x, y) >= imgMag(x1, y1) && imgMag(x, y) >= imgMag(x2, y2)) 
				imgMax->at(x, y) = imgMag(x, y);
		}
	}
	return imgMax;
}

img* nonmax_sup2(img& imgMag, img& imgDir, double layer1_threshold, int layer1_3x3bound)
{
    int w = (int)imgMag.width;
    int h = (int)imgMag.height;

    if (w < 3 || h < 3) return nullptr;

    img* sum2x2;
    HOP_REAL sum2x2max = imgMag.p_get_2x2sum(sum2x2);

    HOP_REAL* ptrs[9];
    HOP_REAL* ptrs2[9];
    HOP_REAL threshold = imgMag.maximum()*layer1_threshold;
    HOP_REAL threshold2 = sum2x2max*layer1_threshold;
    
    //cout << endl << threshold << endl << layer1_threshold << endl;

    ptrs[0] = imgMag.ptr(0, 0); ptrs2[0] = sum2x2->ptr(0, 0);
    ptrs[1] = imgMag.ptr(1, 0); ptrs2[1] = sum2x2->ptr(1, 0);
    ptrs[2] = imgMag.ptr(2, 0); ptrs2[2] = sum2x2->ptr(2, 0);
    ptrs[3] = imgMag.ptr(0, 1); ptrs2[3] = sum2x2->ptr(0, 1);
    ptrs[4] = imgMag.ptr(1, 1); ptrs2[4] = sum2x2->ptr(1, 1);
    ptrs[5] = imgMag.ptr(2, 1); ptrs2[5] = sum2x2->ptr(2, 1);
    ptrs[6] = imgMag.ptr(0, 2); ptrs2[6] = sum2x2->ptr(0, 2);
    ptrs[7] = imgMag.ptr(1, 2); ptrs2[7] = sum2x2->ptr(1, 2);
    ptrs[8] = imgMag.ptr(2, 2); ptrs2[8] = sum2x2->ptr(2, 2);

    HOP_REAL** mid = &ptrs[4];
    HOP_REAL** mid2 = &ptrs2[4];;
    int gcount, gcount2;
    int k;
    node* n;
	img* imgMax = new img(w, h, 0.0);

    for (int j = 1; j < h - 1; ++j) {
        for (int i = 1; i < w - 1; ++i) {
            gcount = gcount2 = 0;
            if (**mid >= threshold) {
                for (k = 0; k < 9; ++k) { if (*ptrs[k] >= **mid) ++gcount; }
                if (gcount <= layer1_3x3bound) {
					imgMax->at(i, j) = **mid;

                } else {
                    for (k = 0; k < 9; ++k) { if (*ptrs2[k] >= **mid2) ++gcount2; }
                    if (gcount2 <= layer1_3x3bound) {
						imgMax->at(i, j) = **mid;
                    }
                }
            }
            
            for (k = 0; k < 9; ++k) { ++(ptrs[k]); ++(ptrs2[k]); }
        }
        for (int k = 0; k < 9; ++k) { ptrs[k] += 2; ptrs2[k] += 2; }
    }
    if (sum2x2 != nullptr) delete sum2x2;   
	return imgMax;
}

layer1_creator_color::layer1_creator_color() :
	layer1_creator()
{
	gabor_size = 7;
	n_rotations = 6;
	gauss_sigma = 5.0;
	vmd_kappa = 5.0;
	gray_weight = 3.0;
	hard_resize = false;
	resize_sigma = 1.0;
}

layer1_creator_color::~layer1_creator_color()
{
	delete gaussFilter;
	delete gaussFilterT;
	delete gaussDxFilter;
	delete gaussDxFilterT;
}


void layer1_creator_color::cfg_init(const config_dictionary& cfg)
{
    layer1_creator::cfg_init(cfg);
    cfg.get_value(gabor_size, "gabor_size");
    cfg.get_value(n_rotations, "n_rotations");
    cfg.get_value(gauss_sigma, "gauss_sigma");                 
	cfg.get_value(vmd_kappa, "vmd_kappa");
	cfg.get_value(gray_weight, "gray_weight");
	cfg.get_value(hard_resize, "hard_resize");
	cfg.get_value(resize_sigma, "resize_sigma");
}

img* resize_and_blur(img* im, int neww, int newh, double blursigma)
{
	img* result = new img(*im);
	if (blursigma > 0)
		result->blur(3*blursigma, blursigma, false);
	*result = result->get_resized(neww, newh, false);
	return result;
}

// layer1 only!
void get_resized_result(layer1_result* newres, layer1_result* res, int neww, int newh)
{
	newres->new_grid(neww + 2*res->border, newh + 2*res->border, 0);
	newres->shape_nodes.resize(1);
	newres->shape_nodes_inhib.resize(1);

	double f = (double)neww/(res->x_size(0) - 2*res->border);

	for (auto niter = res->nodes.begin(); niter != res->nodes.end(); ++niter) {
		node* n = *niter;
		layer1_data* nd = (layer1_data*)n->data;

		if (nd->z != 0) continue;
		
		int newx = int_round(f*(nd->x - res->border));
		int newy = int_round(f*(nd->y - res->border));

		layer1_data* nnd = new layer1_data(nd->val(), nd->m);

		nnd->r.set_response(G_RESPONSE, 1.0);
		nnd->r.set_response(RR_RESPONSE, 1.0);
        nnd->r.set_response(S_RESPONSE, 0.0);

		newres->add_grid_node(nnd, newx + res->border, newy + res->border, 0);
	}
	for (auto niter = newres->nodes.begin(); niter != newres->nodes.end(); ++niter) {
		node* nn = *niter;

		if (nn->is_attr_set(IMG_NODE_ATTR)) 
			newres->shape_nodes[0].push_back(nn);
	}
	newres->inhibit(0);
}

void resize(int& width, int& height, int newsize)
{
	int newh, neww;

	if (newsize < 0) {
		neww = (int)(-(double)newsize*width/100.0);
		newh = (int)(-(double)newsize*height/100.0);
	} else {
        if (width >= height) { // "landscape"
            neww = newsize;
            newh = (newsize*height)/width;
        } else {
            newh = newsize;
            neww = (newsize*width)/height;
        }
	}
	width = neww;
	height = newh;
}

int layer1_creator_color::create_result_vector(vector<layer1_result*>& results, const img&, const vector<irectangle2>&)
{
    cout << "Called unimplemented version of create_result_vector" << endl;
    cout << "layer1_creator_color::create_result_vector(vector<layer1_result*>& results, const img&, const vector<irectangle2>&)" << endl;
    throw new_libhop_exception("layer1_creator_color::create_result_vector not implemented.");
}


layer1_result* layer1_creator_color::create_result(img& im, img& maskim, int originalw)
{
    cout << "Calling layer1_creator_color::create_result." << endl;
    cout << "layer1_creator_color::create_result_vector should be called in \"color\" mode." << endl;
    throw new_libhop_exception("layer1_creator_color::create_result called");
}


int layer1_creator_color::create_result_vector(vector<layer1_result*>& results, const img& srcim, const img& mask) 
{
	results.clear();

	if (srcim.empty()) return -1;
	if (srcim.grayscale) {
		cout << "layer1_creator_color::create_result: color srcim expeced." << endl;
		return -1;
	}

	img im;
    img mim;
    int originalw = srcim.width;
	int swidth = (int)srcim.width;
	int sheight = (int)srcim.height;
	bool skipfirst;

	resize(swidth, sheight, init_size);

	if (hard_resize || swidth < (int)im.width || sheight < (int)im.height) {
		if (init_size != -100) srcim.get_resized(im, init_size); else im = srcim;
		if (!mask.empty() && init_size != -100) mask.get_resized(mim, init_size); else mim = mask;
		swidth = (int)im.width;
		sheight = (int)im.height;
		skipfirst = false;
	} else {
		im = srcim;
		mim = mask;
		skipfirst = true;
	}
	if (scale_mask_size >= 3 && scale_sigma > 0.0) 
		im.blur(scale_mask_size, scale_sigma);

	img R, G, B;
	img gray;
	vector<img*> M;
	
	im.get_colors(R, G, B);
	//R.save("d:\\work\\R.png");
	//G.save("d:\\work\\G.png");
	//B.save("d:\\work\\B.png");
	gray = im;
	gray.to_grayscale();

	SO(M, R, G, B);
	//for (int i = 0; i < M.size(); ++i) {
	//	M[i]->save(string("d:\\work\\M_") + i + string(".png"));
	//}

	img* imgMag, * imgDir, * imgMax;

	edges(imgMag, imgDir, imgMax, M, gray);
	//imgMax->save("d:\\work\\Max.png");
	//imgDir->save("d:\\work\\Dir.png");
	//imgMag->save("d:\\work\\Mag.png");

	img* imgMagCurr;
	img* imgDirCurr;
	img* imgMaxCurr;

	if (skipfirst) {
		imgMagCurr = resize_and_blur(imgMag, swidth, sheight, resize_sigma);
		imgDirCurr = resize_and_blur(imgDir, swidth, sheight, resize_sigma);
		imgMaxCurr = resize_and_blur(imgMax, swidth, sheight, resize_sigma);
	} else {
		imgMagCurr = new img(*imgMag);
		imgDirCurr = new img(*imgDir);
		imgMaxCurr = new img(*imgMax);
	}

	int count = 0;
	bool process = true;

	while (process) {
		layer1_result* res = get_result(imgMagCurr, imgDirCurr, imgMaxCurr);

		res->original_width = originalw;
		res->border = border;
		inhibit_result(res);			
		results.push_back(res);
		
		++count;

		swidth = (int)(swidth*scale_factor);
		sheight = (int)(sheight*scale_factor);

		if (imgMagCurr != nullptr) delete imgMagCurr;
		if (imgDirCurr != nullptr) delete imgDirCurr;
		if (imgMaxCurr != nullptr) delete imgMaxCurr;

		if (swidth < scale_limit && sheight < scale_limit || count >= max_scales) 
			process = false;
		else {
			imgMagCurr = resize_and_blur(imgMag, swidth, sheight, resize_sigma);
			imgDirCurr = resize_and_blur(imgDir, swidth, sheight, resize_sigma);
			imgMaxCurr = resize_and_blur(imgMax, swidth, sheight, resize_sigma);
            //im = im.get_resized(swidth, sheight);
            //if (!mim.empty()) mim = mim.get_resized(swidth, sheight);
        }

	}



	//layer1_result* res = get_result(imgMag, imgDir, imgMax);

 //   res->original_width = originalw;
 //   res->border = border;
	//inhibit_result(res);			
	//results.push_back(res);

	//int count = 1;
	//bool process = true;

	//while (process) {
 //       ++count;

	//	swidth = (int)(swidth*scale_factor);
	//	sheight = (int)(sheight*scale_factor);

	//	if (swidth < scale_limit && sheight < scale_limit || count >= max_scales) 
 //           process = false;
 //       else {
 //           //im = im.get_resized(swidth, sheight);
 //           //if (!mim.empty()) mim = mim.get_resized(swidth, sheight);
 //       }

	//}
	for (auto Miter = M.begin(); Miter != M.end(); ++Miter) {
		if (*Miter != nullptr) 
			delete *Miter;
	}
	if (imgMag != nullptr) delete imgMag;
	if (imgDir != nullptr) delete imgDir;
	if (imgMax != nullptr) delete imgMax;

	return 0;
}

void layer1_creator_color::init_from_result(layer1_result_color* res)
{
    bool changed;

    changed = init_default_from_result(res);

    if (gabor_size != res->gabor_size) { changed = true; gabor_size = res->gabor_size; }
    if (n_rotations != res->n_rotations) { changed = true; n_rotations = res->n_rotations; }

	if (changed) delete_filter_kernels();
}

layer1_result* layer1_creator_color::init_result()
{
    layer1_result_color* result1 = new layer1_result_color();
    
    init_result_with_default(result1);

    result1->gabor_size = gabor_size;
    result1->n_rotations = n_rotations;
    return result1;
}


double layer1_creator_color::Omega(int ch, int col)
{
	const double Omega_arr[][3] = {
		{1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0}, 
		{2.0/sqrt(6.0), -1.0/sqrt(6.0), -1.0/sqrt(6.0)}, 
		{1.0/sqrt(6.0), 1.0/sqrt(6.0), -2.0/sqrt(6.0)},
		{1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0)},
		{-1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0}, 
		{-2.0/sqrt(6.0), 1.0/sqrt(6.0), 1.0/sqrt(6.0)}, 
		{-1.0/sqrt(6.0), -1.0/sqrt(6.0), 2.0/sqrt(6.0)},
		{-1.0/sqrt(3.0), -1.0/sqrt(3.0), -1.0/sqrt(3.0)}
	};

	return Omega_arr[ch][col];
}

int layer1_creator_color::Omega_size() 
{
	return 8;
}

void layer1_creator_color::make_filter_kernels()
{
    //get_masks(masks, gabor_size, int_round(180/n_rotations));
	for (auto miter = filter_kernels.begin(); miter != filter_kernels.end(); ++miter)
		if (*miter != nullptr) delete *miter;
	filter_kernels.resize(3);
	masks_norm.resize(3);
	filter_kernels[0] = new img(gabor_mask_rot(gabor_size, -1));
	filter_kernels[0]->absolute_value();
	masks_norm[0] = filter_kernels[0]->norm();
	filter_kernels[1] = new img(gabor_mask_rot(gabor_size, 1));
	filter_kernels[1]->absolute_value();
	masks_norm[1] = filter_kernels[1]->norm();
	filter_kernels[2] = new img(gabor_mask_rot(gabor_size, 0));
	masks_norm[2] = filter_kernels[2]->norm();

	gaussFilter = new img(gauss_filter1(gauss_sigma));
	gaussFilterT = new img(gauss_filter1(gauss_sigma, true));
	gaussDxFilter = new img(gaussdx_filter1(gauss_sigma));
	gaussDxFilterT = new img(gaussdx_filter1(gauss_sigma, true));
}

void layer1_creator_color::get_part_data(vector<part_data*>& pd, const config_dictionary& cfg)
{
    vector<img*> masks;
    vector<img*> r_masks;
    matrix<bool>* region;

    //get_masks(masks, gabor_size, n_rotations);
    //get_masks(r_masks, gabor_size, n_rotations);
    get_masks(masks);
    get_masks(r_masks);

    size_t nmasks2 = masks.size()/2;

    pd.clear();
    for (size_t i = 0; i < nmasks2; ++i) {
        region = r_masks[i]->get_bool_matrix(layer1_region_threshold);
        pd.push_back(new part_data(*masks[i], *region, 0, (int)i));
        delete region;
    }
    for (size_t i = 0; i < masks.size(); ++i) {
        delete masks[i];
        delete r_masks[i];
    }
}

void layer1_creator_color::SO(vector<img*>& M, img& R, img& G, img& B)
{
	const double Efactor2 = 0.225*0.225;

	img* fp[3];
	img* fm[3];

	#pragma omp parallel sections
	{
		#pragma omp section
		{
			fp[0] = R.get_convolve(*filter_kernels[1]);
		}
		#pragma omp section
		{
			fp[1] = G.get_convolve(*filter_kernels[1]);
		}
		#pragma omp section
		{
			fp[2] = B.get_convolve(*filter_kernels[1]);
		}
		#pragma omp section
		{
			fm[0] = R.get_convolve(*filter_kernels[0]);
		} 
		#pragma omp section
		{
			fm[1] = G.get_convolve(*filter_kernels[0]);
		}
		#pragma omp section
		{
			fm[2] = B.get_convolve(*filter_kernels[0]);
		}
	}

	int w = fp[0]->width;
	int h = fp[0]->height;
	int channels = Omega_size();

	M.resize(channels);
	for (int channel = 0; channel < channels; ++channel) {
		img* Mc = M[channel] = new img(w, h, 0.0);
		
		for (int color = 0; color < 3; ++color) {
			double o = Omega(channel, color);

			if (o != 0) {
				img* convolution = (o > 0) ? fp[color] : fm[color];
				double norm = (o > 0) ? masks_norm[1] : masks_norm[0];
				//double f = (o > 0) ? 1 : -1;
				double f = o;

				if (norm > 0) f /= norm;
				for (int i = 0; i < convolution->size(); ++i) 
					(*Mc)[i] += f*(*convolution)[i];
			}
		}
		for (int i = 0; i < Mc->size(); ++i) 
			if ((*Mc)[i] < 0) (*Mc)[i] = 0.0;
	}

	//M[0]->save("d:\\work\\M0DBG.png");

	img E(w, h, 0.0);

	for (int channel = 0; channel < channels; ++channel) {
		img* Mc = M[channel];

		for (int i = 0; i < E.size(); ++i) 
			E[i] += sqr((*Mc)[i]);
	}
	for (int i = 0; i < E.size(); ++i)
		E[i] /= channels;
	for (int channel = 0; channel < channels; ++channel) {
		img* Mc = M[channel];

		for (int i = 0; i < E.size(); ++i) 
			(*Mc)[i] = ::sqrt(sqr((*Mc)[i])/(Efactor2 + E[i]));
	}

	for (int i = 0; i < 3; ++i) {
		delete fp[i];
		delete fm[i];
	}
}


void layer1_creator_color::edges(img*& imgMag, img*& imgDir, img*& imgMax, vector<img*>& M, img& gray)
{
	if (M.empty()) {
		imgMag = imgDir = imgMax = nullptr;
		return;
	}

	int w = M[0]->width;
	int h = M[0]->height;

	imgMag = new img(w, h, 0.0);
	imgDir = new img(w, h, 0.0);

	img* imgDx = new img(w, h, 0.0);
	img* imgDy = new img(w, h, 0.0);
	double cnt = 0.0;

	// color channels
	for (int channel = 0; channel < M.size(); ++channel) {
		if (M.size() == 8 && (channel == 3 || channel == 7))
			continue;

		double alpha = channel < M.size()/2 ? 1 : -1;
		pair<img*, img*> imgD = gauss_deriv(*M[channel], gauss_sigma);
		img* Mc = M[channel];
		
		for (int i = 0; i < Mc->size(); ++i) {
			(*imgMag)[i] += sqr((*imgD.first)[i]) + sqr((*imgD.second)[i]);
			(*imgDx)[i] += alpha * (*imgD.first)[i];
			(*imgDy)[i] += alpha * (*imgD.second)[i];
		}
		cnt += 1.0;
		delete imgD.first;
		delete imgD.second;
	}

	// gray channel
	pair<img*, img*> imgDg = gauss_deriv(gray, gauss_sigma);
	for (int i = 0; i < gray.size(); ++i) {
		(*imgMag)[i] += sqr((*imgDg.first)[i]) + sqr((*imgDg.second)[i])*gray_weight;
		(*imgDx)[i] += (*imgDg.first)[i];
		(*imgDy)[i] += (*imgDg.second)[i];
	}
	cnt += gray_weight;
	delete imgDg.first;
	delete imgDg.second;

	for (int i = 0; i < imgMag->size(); ++i) {
		(*imgMag)[i] /= cnt;
		(*imgMag)[i] = sqrt((*imgMag)[i]);
		(*imgDx)[i] /= cnt;
		(*imgDy)[i] /= cnt;
		(*imgDir)[i] = (fabs((*imgDy)[i]) < 1E-6 && fabs((*imgDx)[i]) < 1E-6 == 0.0) ? 0.0 : atan2((*imgDy)[i], (*imgDx)[i]);
	}
	imgMax = nonmax_sup(*imgMag, *imgDir); //, layer1_threshold, layer1_3x3bound);
	delete imgDx;
	delete imgDy;
}

layer1_result* layer1_creator_color::get_result(img* imgMag, img* imgDir, img* imgMax)
{
	const int operm[] = { 5, 4, 3, 2, 1, 0 };
	
	vector<von_mises_distribution> vmd(12);

	for (int i = 0; i < 12; ++i) {
		vmd[i].reset_kappa(vmd_kappa);
		vmd[i].reset_mu((-180.0 + (i + 1)*30.0)*M_PI/180.0);
	}

	int w = imgMax->width;
	int h = imgMax->height;
	layer1_result* result = init_result();
	double threshold = imgMax->maximum()*layer1_threshold;

	result->new_grid(w + 2*border, h + 2*border, 0);
	result->shape_nodes.resize(1);
	result->shape_nodes_inhib.resize(1);

	vector<node*>& s_nodes = result->shape_nodes[0];

	for (int i = 0; i < w; ++i) {
		for (int j = 0; j < h; ++j) {
			double val = imgMax->at(i, j);

			if (val < threshold) 
				continue;

			double angle = imgDir->at(i, j);
			vector<pair<double, int> > oval(6);

			for (int k = 0; k < 6; ++k) {
				oval[k].second = k;
				oval[k].first = max(vmd[k].pdf_val1(angle), vmd[k + 6].pdf_val1(angle));
			}
			sort(oval.begin(), oval.end(), greater<pair<double, int> >());

			double maxoval = oval[0].first;

			for (int k = 0; k < 6; ++k) {
				if (oval[k].first >= maxoval*response_percent) {
					double f = oval[k].first/maxoval;

					f = f*val;
					if (power_correction != 1) f = pow(f, power_correction);

					layer1_data* nnd = new layer1_data(f, operm[oval[k].second]);

		            nnd->r.set_response(G_RESPONSE, 1.0);
		            nnd->r.set_response(RR_RESPONSE, 1.0);
                    nnd->r.set_response(S_RESPONSE, 0.0);

					node* n = result->add_grid_node(nnd, i + border, j + border, 0);

					if (k == 0) s_nodes.push_back(n);
				}
			}

		}
	}
	sort(s_nodes.begin(), s_nodes.end(), response_sort_f(R_RESPONSE)); 
	return result;
}

void layer1_creator_color::get_masks(vector<img*>& masks, int, int)
{
    const double gabor_lambda = 6.0;
    const double gabor_gamma = 0.75;
    const double gabor_bw = 1.75;
    const double gabor_step = 30;

    int nmasks2 = 180/gabor_step;

    masks.resize(2*nmasks2);
    for (int angle = 0, i = 0; i < nmasks2; angle += gabor_step, ++i) {
        masks[i] = img::gabor_mask(gabor_lambda, angle, 0.0, gabor_gamma, gabor_bw);
        masks[nmasks2 + i] = img::gabor_mask(gabor_lambda, angle, 90.0, gabor_gamma, gabor_bw);
    }

}

void layer1_creator_color::get_masks(vector<img*>& masks)
{
	get_masks(masks, 0, 0);	// dummy args!
}


	/*

ArcTan2[x_, y_] := 
 If[Abs[x] < 10^-6 && Abs[y] < 10^-6, 0, ArcTan[x, y]]
SetAttributes[ArcTan2, Listable]

ColorEdges[M_] :=
 Module[{i, sigma, imgDx, imgDy, imgMag, cntr, alph, imgDxc, imgDyc, 
   imgDir, imgMax},
  sigma = 5;
  imgDx = imgDy = 0;
  imgMag = 0;
  cntr = 0;
  For[i = 1, i <= 8, i++,
   If[i != 4 && i != 88,
     alph = If[i > 4, -1, 1];
     {imgDxc, imgDyc} = GaussDeriv[M[[i]], sigma];
     imgMag = imgMag + imgDxc^2 + imgDyc^2;
     imgDx = imgDx + alph*imgDxc;
     imgDy = imgDy + alph*imgDyc;
     cntr++;
     ];
   ];
  imgMag = imgMag/cntr;
  imgDx = imgDx/cntr;
  imgDy = imgDy/cntr;
  imgMag = Sqrt[imgMag];
  imgDir = ArcTan2[imgDx, imgDy];
  imgMax = NonmaxSup[imgMag, imgDir];
  {imgMag, imgDir, imgMax}
  ]


  */


// layer1_creator_app
///////////////////////////////////////////////////////////////////////////////

layer1_creator_app::layer1_creator_app() :
    layer1_creator()
{
    gabor_lambda = 10.0;
    gabor_gamma = 0.8;
    gabor_bw = 8.0;
    gabor_step = 45;
}

layer1_creator_app::~layer1_creator_app()
{

}

void layer1_creator_app::init_from_result(layer1_result_app* res)
{
    bool changed;

    init_default_from_result(res);

    if (gabor_lambda != res->gabor_lambda) { changed = true; gabor_lambda = res->gabor_lambda; }
    if (gabor_gamma != res->gabor_gamma) { changed = true; gabor_gamma = res->gabor_gamma; }
    if (gabor_bw != res->gabor_bw) { changed = true; gabor_bw = res->gabor_bw; }
    if (gabor_step != res->gabor_step) { changed = true; gabor_step = res->gabor_step; }

    if (changed) delete_filter_kernels();
}

layer1_result* layer1_creator_app::init_result()
{
    layer1_result_app* result1 = new layer1_result_app();

    init_result_with_default(result1);

    result1->gabor_lambda = gabor_lambda;
    result1->gabor_gamma = gabor_gamma;
    result1->gabor_bw = gabor_bw;
    result1->gabor_step = gabor_step;
    return result1;
}

void layer1_creator_app::cfg_init(const config_dictionary& cfg)
{
    layer1_creator::cfg_init(cfg);
    cfg.get_value(gabor_lambda, "gabor_lambda");
    cfg.get_value(gabor_gamma, "gabor_gamma");
    cfg.get_value(gabor_bw, "gabor_bw");                 
    cfg.get_value(gabor_step, "gabor_step");
}


/*
void layer1_creator_app::make_masks()
{
    int nmasks = 360/gabor_step;
    int nmasks2 = nmasks/2;

    masks.resize(nmasks);
    regions.resize(nmasks);
    for (int angle = 0, i = 0; angle < 360; angle += gabor_step, ++i) {
        masks[i] = img::gabor_mask(gabor_lambda, angle, 90.0, gabor_gamma, gabor_bw);

        img* tmpmask = img::gabor_mask(gabor_lambda, angle, 0.0, gabor_gamma, gabor_bw);
        regions[i] = tmpmask->get_01_img(layer1_region_threshold);
        delete tmpmask;
    }
   
}
*/

void layer1_creator_app::make_filter_kernels()
{
    get_filter_kernels(filter_kernels, gabor_lambda, gabor_gamma, gabor_bw, gabor_step);
}

void layer1_creator_app::get_filter_kernels(vector<img*>& filter_kernels, 
    double gabor_lambda, double gabor_gamma, double gabor_bw, int gabor_step)
{
    int nfilter_kernels = 360/gabor_step;

	filter_kernels.resize(nfilter_kernels);
    for (int angle = 0, i = 0; i < nfilter_kernels; angle += gabor_step, ++i) {
        filter_kernels[i] = img::gabor_mask(gabor_lambda, angle, 90.0, gabor_gamma, gabor_bw);
    }
	//filter_kernels.resize(nfilter_kernels/2);
}

void layer1_creator_app::make_images(vector<img*>& images, img& im, img& immask)
{
    img* sinim, * sinim2;
    int nmasks = (int)filter_kernels.size();
	int nmasks2 = nmasks/2;
	
    images.resize(make_images_get_count());
    for (int i = 0; i < nmasks2; ++i) {
        sinim = im.get_convolve(*filter_kernels[i], immask);
        sinim2 = new img(*sinim);
        sinim->positive();
        sinim2->negative2();
        images[i] = sinim;
        images[nmasks2 + i] = sinim2;
    }   
}

int layer1_creator_app::make_images_get_count() {

	// layer1_creator_app creates two images from one convolution filter
	return filter_kernels.size();//*2;
}

void layer1_creator_app::get_part_data(vector<part_data*>& pd, const config_dictionary& cfg)
{
    vector<img*> filter_kernels;
    vector<img*> r_filter_kernels;
    matrix<bool>* region;

    get_filter_kernels(filter_kernels, gabor_lambda, gabor_gamma, gabor_bw, gabor_step);
    get_filter_kernels(r_filter_kernels, 
        cfg.get_value_double("region_lambda", 7.0),
        cfg.get_value_double("region_gamma", 1.0),
        cfg.get_value_double("region_bw", 2.75),
        gabor_step);

    pd.clear();
    for (size_t i = 0; i < filter_kernels.size(); ++i) {
        region = r_filter_kernels[i]->get_bool_matrix(layer1_region_threshold);
        pd.push_back(new part_data(*filter_kernels[i], *region, 0, (int)i));
        delete region;
    }
    for (size_t i = 0; i < filter_kernels.size(); ++i) {
        delete filter_kernels[i];
        delete r_filter_kernels[i];
    }
}



// layer1_creator_dog
///////////////////////////////////////////////////////////////////////////////

layer1_creator_dog::layer1_creator_dog() :
    layer1_creator()
{
    sigma_inner = sqrt(2.0);
    sigma_outer = 2;
    mask_size_factor = 2.8;
}

layer1_creator_dog::~layer1_creator_dog()
{

}

void layer1_creator_dog::init_from_result(layer1_result_dog* res)
{
    bool changed;

    changed = init_default_from_result(res);

    if (sigma_inner != res->sigma_inner) { changed = true; sigma_inner = res->sigma_inner; }
    if (sigma_outer != res->sigma_outer) { changed = true; sigma_outer = res->sigma_outer; }
    if (mask_size_factor != res->mask_size_factor) { changed = true; mask_size_factor = res->mask_size_factor; }

    if (changed) delete_filter_kernels();
}

layer1_result* layer1_creator_dog::init_result()
{
    layer1_result_dog* result1 = new layer1_result_dog();

    init_result_with_default(result1);

    result1->sigma_inner = sigma_inner;
    result1->sigma_outer = sigma_outer;
    result1->mask_size_factor = mask_size_factor;
    return result1;
}


/*void layer1_creator_dog::make_masks()
{
    int maskw = 2*(int)(mask_size_factor*sigma_outer + 0.5) + 1;
    img* maskin = img::gaussian_mask(maskw, maskw, sigma_inner);
    img* maskout = img::gaussian_mask(maskw, maskw, sigma_outer);

    masks.resize(2);
    regions.resize(2);

    regions[0] = maskin->get_01_img(layer1_region_threshold);
    regions[1] = new iimg(*(regions[0]));

    *maskin -= *maskout;
    masks[0] = maskin;
    masks[1] = new img(maskin);
    masks[1]->neg();

    delete maskout;
}*/

void layer1_creator_dog::make_filter_kernels()
{
    get_filter_kernels(filter_kernels, sigma_inner, sigma_outer, mask_size_factor);
}

void layer1_creator_dog::get_filter_kernels(vector<img*>& filter_kernels, 
    double sigma_inner, double sigma_outer, double filter_kernel_size_factor)
{
    int filter_kernel_w = 2*(int)(filter_kernel_size_factor*sigma_outer + 0.5) + 1;
    img* filter_kernel_in = img::gaussian_mask(filter_kernel_w, filter_kernel_w, sigma_inner);
    img* filter_kernel_out = img::gaussian_mask(filter_kernel_w, filter_kernel_w, sigma_outer);

    filter_kernels.resize(2);

    *filter_kernel_in -= *filter_kernel_out;
    filter_kernels[0] = filter_kernel_in;
    filter_kernels[1] = new img(*filter_kernel_in);
    filter_kernels[1]->neg();

    delete filter_kernel_out;
}

void layer1_creator_dog::make_images(vector<img*>& images, img& im, img& immask)
{
    images.resize(make_images_get_count());

    images[0] = im.get_convolve(*filter_kernels[0], immask);
    images[1] = new img(*images[0]);
    images[0]->positive();
    images[1]->negative2();
}


int layer1_creator_dog::make_images_get_count() {

	// layer1_creator_dog creates only two images
	return 2;
}



void layer1_creator_dog::cfg_init(const config_dictionary& cfg)
{
    layer1_creator::cfg_init(cfg);
    cfg.get_value(sigma_inner, "sigma_inner");
    cfg.get_value(sigma_outer, "sigma_outer"); 
    cfg.get_value(mask_size_factor, "mask_size_factor");                 
}

void layer1_creator_dog::get_part_data(vector<part_data*>& pd, const config_dictionary& cfg)
{
    vector<img*> filter_kernels;
    vector<img*> r_filter_kernels;
    matrix<bool>* region;

    get_filter_kernels(filter_kernels, sigma_inner, sigma_outer, mask_size_factor);
    //get_filter_kernels(filter_kernels, gabor_lambda, gabor_gamma, gabor_bw, gabor_step);
    get_filter_kernels(r_filter_kernels, 
        cfg.get_value_double("region_sigma_inner", 1.41),
        cfg.get_value_double("region_sigma_outer", 2.0),
        cfg.get_value_double("region_mask_size_factor", 2.8));

    pd.clear();
    for (size_t i = 0; i < filter_kernels.size(); ++i) {
        region = r_filter_kernels[i]->get_bool_matrix(layer1_region_threshold);
        pd.push_back(new part_data(*filter_kernels[i], *region, 0, (int)i));
        delete region;
    }
    for (size_t i = 0; i < filter_kernels.size(); ++i) {
        delete filter_kernels[i];
        delete r_filter_kernels[i];
    }
}



// layer1_creator_loggabor
///////////////////////////////////////////////////////////////////////////////

layer1_creator_loggabor::layer1_creator_loggabor() :
    layer1_creator()
{
    
}

layer1_creator_loggabor::~layer1_creator_loggabor()
{

}

void layer1_creator_loggabor::init_from_result(layer1_result_loggabor* res)
{
    bool changed;

    changed = init_default_from_result(res);

    if (changed) delete_filter_kernels();
}

layer1_result* layer1_creator_loggabor::init_result()
{
    layer1_result_loggabor* result1 = new layer1_result_loggabor();

    init_result_with_default(result1);

    return result1;
}

void layer1_creator_loggabor::make_filter_kernels()
{
    get_filter_kernels(filter_kernels);
}

void layer1_creator_loggabor::get_filter_kernels(vector<img*>& filter_kernels, double dummy)
{
    int nfilter_kernels = 6;

    filter_kernels.resize(nfilter_kernels);
    for (int i = 0; i < 6; ++i) {
        filter_kernels[i] = img::log_gabor_mask(5, i + 1);
    }
}

void layer1_creator_loggabor::make_images(vector<img*>& images, img& im, img& immask)
{
	int nmasks = (int)filter_kernels.size();
	
    images.resize(make_images_get_count());
    for (int i = 0; i < nmasks; ++i) {
        images[i] = im.get_convolve(*filter_kernels[i],immask);
    }   
}


int layer1_creator_loggabor::make_images_get_count() {

	return filter_kernels.size();
}


void layer1_creator_loggabor::get_part_data(vector<part_data*>& pd, const config_dictionary& cfg)
{
    vector<img*> masks;
    vector<img*> r_masks;
    matrix<bool>* region;

    get_masks(masks);
    get_masks(r_masks);

    size_t nmasks = masks.size();

    pd.clear();
    for (size_t i = 0; i < nmasks; ++i) {
        region = r_masks[i]->get_bool_matrix(layer1_region_threshold);
        pd.push_back(new part_data(*masks[i], *region, 0, (int)i));
        delete region;
    }
    for (size_t i = 0; i < masks.size(); ++i) {
        delete masks[i];
        delete r_masks[i];
    }
}


void layer1_creator_loggabor::cfg_init(const config_dictionary& cfg)
{
    layer1_creator::cfg_init(cfg);

}
