/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// feature_extraction

#include "color_feature_extraction.h"

#include "utils/class_register.h"

void register_shape_color_feature_extraction_classes() {
	/************************************************************************* 
	 *           !!!! DO NOT MODIFY OR REMOVE THIS FUNCTION   !!!!!
	 ************************************************************************/

	//ClassRegister::registerFactory<AbstractFeatureExtraction::IFactory, ColorEdgesFeatureExtraction::Factory>();	
}
/*
// gsign: +1, -1, 0 (positive, negative, normal)
img gabor_mask_rot(int size, double lambda, double sigma, int gsign)
{
	double size2 = (double)size/2;
	img result(size, size);
	double mean = 0.0;

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


void layer1_creator_color::cfg_init(const ConfigDictionary& cfg)
{
    layer1_creator::cfg_init(cfg);
    cfg.getValue(gabor_size, "gabor_size");
    cfg.getValue(n_rotations, "n_rotations");
    cfg.getValue(gauss_sigma, "gauss_sigma");                 
	cfg.getValue(vmd_kappa, "vmd_kappa");
	cfg.getValue(gray_weight, "gray_weight");
	cfg.getValue(hard_resize, "hard_resize");
	cfg.getValue(resize_sigma, "resize_sigma");
}

img* resize_and_blur(img* im, int neww, int newh, double blursigma)
{
	img* result = new img(*im);
	if (blursigma > 0)
		result->blur(3*blursigma, blursigma, false);
	*result = result->get_resized(neww, newh, false);
	return result;
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

	gray = im;
	gray.to_grayscale();

	SO(M, R, G, B);
	
	img* imgMag, * imgDir, * imgMax;

	edges(imgMag, imgDir, imgMax, M, gray);
	
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
        }

	}

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

void layer1_creator_color::get_part_data(vector<part_data*>& pd, const ConfigDictionary& cfg)
{
    vector<img*> masks;
    vector<img*> r_masks;
    matrix<bool>* region;

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
}*/