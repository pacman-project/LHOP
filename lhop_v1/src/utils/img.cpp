//
// Wrapper for image functions
///////////////////////////////////////////////////////////////////////////////

#include <cv.h>
#include <highgui.h>
#include <cxcore.h>
#include <cv.hpp>
#include <cxcore.hpp>

#include <iostream>
#include <string.h>
#include <fstream>
#include <limits>
#include "img.h"
#include "fonts.h"





#ifndef WIN32
typedef unsigned char UCHAR, *PUCHAR; 
#endif

typedef unsigned long DWORD;

// local variables
///////////////////////////////////////////////////////////////////////////////

vector<img> font7x11;

// functions
///////////////////////////////////////////////////////////////////////////////

img jet_map()
{
    img result(256, 3, 0.0);
    const int ulen = 191;
    img tmp(ulen, 1);
    int i, j;

    for (i = 0; i < 64; ++i) tmp[i] = (i+1)/64.0;
    for (i = 64; i < 127; ++i) tmp[i] = 1.0;
    for (i = 127; i < ulen; ++i) tmp[i] = (ulen - i)/64.0;
    
    for (i = 96, j = 0; i < 256; ++i, ++j) result(i, 0) = tmp[j];
    for (i = 32, j = 0; i < 223; ++i, ++j) result(i, 1) = tmp[j];
    for (i = 0, j = ulen - 159; i < 159; ++i, ++j) result(i, 2) = tmp[j];

    return result;
}

void make_font(vector<img>& font, const unsigned int* def, int w, int h, HOP_REAL oneval = 1.0)
{
    const unsigned int* ptr = def;
    unsigned int b = 0, val = 0;

    font.clear();
    font.resize(256, img(w, h, 0.0));
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < 256 * w; ++x) {
            b >>= 1; 
            if (!b) { 
                b = 0x80000000; 
                val = *(ptr++); 
            }
            font[x/w](x % w, y) = (val & b) ? oneval : 0;
        }
    }
}

img* get_char(char c)
{
    if (font7x11.empty()) make_font(font7x11, font7x11_definition, 7, 11, 1.0);
    return &font7x11[c];
}

cv::Mat to_cv_image(const img& im) 
{
	int rows = im.height; 
	int cols = im.width;

	if (im.grayscale) {
		int type = (sizeof(HOP_REAL) == sizeof(double)) ? CV_64F : CV_32F;
		void* data = (void*)im.ptr(0, 0);

		return cv::Mat(rows, cols, type, data);
	} else {
		cv::Mat3b result(rows, cols);

		cv::MatIterator_<cv::Vec3b> pixels = result.begin();
		img::const_iterator riter = im.cbegin();
    
	    for (int i = rows*cols; i > 0; --i) {
			IMG_COLOR& ic = (IMG_COLOR&)(*riter);
			(*pixels)[0] = ic.blue;
			(*pixels)[1] = ic.green;
			(*pixels)[2] = ic.red;
			++pixels;
		    ++riter; 
		}
		return result;
    }
}

img from_cv(cv::Mat& image)
{
	bool gscale = image.type() == CV_64F || image.type() == CV_32F;

	img result(image.cols, image.rows, gscale);

	if (gscale) {
		memcpy(result.ptr(0, 0), image.data, sizeof(HOP_REAL)*image.rows*image.cols);
	} else {
		cv::MatConstIterator_<cv::Vec3b> pixels = image.begin<cv::Vec3b>();
		img::iterator riter = result.begin();
    
		for (int i = image.rows*image.cols; i > 0; --i) {
			IMG_COLOR& ic = (IMG_COLOR&)(*riter);
	        ic.red = (*pixels)[2]; 
			ic.green = (*pixels)[1];
			ic.blue = (*pixels)[0];
			++pixels;
			++riter; 
		}
    }
	return result;
}


// color - definitions 
///////////////////////////////////////////////////////////////////////////////

void color::to_color_vector(vector<color>& result, const vector<int>& cvec, int def)
{
    color* colptr;
    int max = (int)cvec.size();
    
    result.clear();
    for (int i = 0; i < max; ++i) {
        if (i % 3 == 0) { 
            result.push_back(color(def, def, def));
            colptr = &result.back();
        }
        colptr->set(i % 3, cvec[i]);
    }

}

void color::to_color(color& result, const vector<int>& cvec, int def)
{
    int size = (int)cvec.size();

    result.red() = (size > 0) ? cvec[0] : def;
    result.green() = (size > 1) ? cvec[1] : def;
    result.blue() = (size > 2) ? cvec[2] : def;
}

// h in [0, 360)
// s in [0, 1]
// v in [0, 1]
// Formula from: http://en.wikipedia.org/wiki/HSL_and_HSV
color color::from_hsv(double h, double s, double v)
{
    double c = v*s;

    h /= 60.0;

    double x = c*(1 - fabs(fmod(h, 2.0) - 1.0));

    x *= 255;
    c *= 255;

    if (0.0 <= h && h < 1.0) return color(c, x, 0); 
    else if (1.0 <= h && h < 2.0) return color(x, c, 0);
    else if (2.0 <= h && h < 3.0) return color(0, c, x);
    else if (3.0 <= h && h < 4.0) return color(0, x, c);
    else if (4.0 <= h && h < 5.0) return color(x, 0, c);
    else if (5.0 <= h && h < 6.0) return color(c, 0, x);
    else return color(0, 0, 0);
}



// img - definitions
///////////////////////////////////////////////////////////////////////////////

img::img(const string& name)
{
	read_grayscale(name, *this);
}

//img::img(const bimg& bim)
//{
//    bim.toImg(*this);
//}

// CUDA code
img::img(matrix<float>* m)
{
	resize(m->width,m->height);
	if(sizeof(HOP_REAL)==sizeof(float)){
		memcpy(ptr(0,0),m->ptr(0,0),m->size()*sizeof(float));
	}
	else if(sizeof(HOP_REAL)==sizeof(double)){
		img::iterator i1=begin();
		matrix<float>::const_iterator i2=m->begin();
		for(;i1!=end();++i1,++i2){
			*i1=(HOP_REAL)*i2;
		}
	}
}
// CUDA code ends

img::img(const img& im, bool gscale, bool normalized)
{
    grayscale = gscale;
    if (im.grayscale == gscale) *this = im;
    else if (gscale) {
        resize(im.width, im.height);

        const HOP_REAL* src = im.ptr(0, 0);
        HOP_REAL* dest = ptr(0, 0);

        for (int i = (int)(width*height) - 1; i >= 0; --i) {
            IMG_COLOR& ic = (IMG_COLOR&)(*src);
            *dest = (0.299*ic.red + 0.587*ic.green + 0.114*ic.blue)/255.0;
            ++dest; ++src;
        }
    } else {
        rmatrix::operator =(im);
        if (!normalized) normalize(0.0, 1.0);

        HOP_REAL* p = ptr(0, 0);

        for (int i = (int)(width*height) - 1; i >= 0; --i) {
            IMG_COLOR& ic = (IMG_COLOR&)(*p);
            ic.red = ic.blue = ic.green = COL_FROM_REAL(*p);
            ++p;
        }
    }

}

void img::gaussian_mask(img& result, int dimx, int dimy, double sigma)
{
    result.resize(dimx, dimy);

    int centerx = dimx/2;
    int centery = dimy/2;
    double sum = 0.0;

    for_each_xy_int(result, x, y) { 
        sum += result(x, y) = ::exp(-((x - centerx)*(x - centerx) + (y - centery)*(y - centery))/2.0/sigma/sigma);
    }
    for_each_xy_int(result, x, y) {
        result(x, y) /= sum;
    }
}

img* img::gaussian_mask(int dimx, int dimy, double sigma)
{
    img* result = new img();

    gaussian_mask(*result, dimx, dimy, sigma);
    return result;
}

img img::gaussian_mask2(int dimx, int dimy, double sigma)
{
    img result;

    gaussian_mask(result, dimx, dimy, sigma);
    return result;
}


img* img::gabor_mask(double lambda, double theta, double phi, double gamma, double bw)
{
    double sigma = (lambda * ::sqrt(::log(2.0)) * (::pow(2.0, bw) + 1.0)) / 
        (M_PI * ::sqrt(2.0) * (::pow(2.0, bw) - 1.0));
    theta = M_PI * theta / 180.0;
    phi = M_PI * phi / 180.0;

    int sz = 2*(int)floor(0.5 + sigma/gamma*2.6) + 1;
    img* F = new img(sz, sz);

    double d7, d8, u, v;
    double d9 = -2.0 * M_PI*M_PI * sigma*sigma / sz / sz;
	double d10 = sz / lambda;
    double d11 = ::sin(-theta);
    double d12 = ::cos(-theta);

	int xymin = -(int)floor(sz / 2.0);
	int xymax = (int)floor(sz / 2.0);

    double sum = 0.0;
    for (int k = xymin; k <= xymax; ++k) {
        for (int l = xymin; l <= xymax; ++l) {
		    u = k * d12 + l * d11;
			v = -k * d11 + l * d12;

            d8 = ::exp(-(u*u + gamma*gamma * v*v) / (2.0 * sigma*sigma));
            d7 = d8 * ::cos((2*M_PI * u) / lambda + phi);

			(*F)(k + sz/2, l + sz/2) = d7;
            sum += d7;
        }
    }
    *F -= sum/(sz*sz);

    return F;
}

img* img::log_gabor_mask(int scale, int orientation) {
	cout << "ERROR: log_gabor_mask NOT IMPLEMENTED - use source code in img_log_gabor.txt" << endl;
	throw std::exception();
	return nullptr;
}

img* img::get_convolve_ocl(img& mask)
{
#ifdef OPENCL
	int ret = 0;
	
	const best_cmd_queue_info& cmd_queue_info = OpenCL::context_manager->getBestCommandQueue();

	cl_context context = OpenCL::context_manager->getContext(cmd_queue_info.context_number);

	cl_command_queue cmd_queue = cmd_queue_info.queue;
	cl_device_id device_id = cmd_queue_info.device_id;

	int context_number = cmd_queue_info.context_number;

	cl_int2 image_size, local_data_size;
	cl_int4 convol_mask_size;

	image_size.s[0] = this->width;
	image_size.s[1] = this->height;

	convol_mask_size.s[0] = mask.width;
	convol_mask_size.s[1] = mask.height;
	convol_mask_size.s[2] = 1; // opencl convolution supports computation of multiply filters on one image - but we will only use one filter

	int kernel_err_host = 0;
	
	cl_mem kernel_err = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, sizeof(cl_int), &kernel_err_host, &ret);
	ocl_check_error(ret, CL_SUCCESS);


	// use pinned buffer since it will make faster memory transfer when using GPU (and will not be slower when using CPU
	cl_mem input_output_pinned = clCreateBuffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_float) * this->size(), nullptr, &ret);
	ocl_check_error(ret, CL_SUCCESS);

	cl_mem convol_mask_pinned = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR, sizeof(cl_float) * mask.size(), nullptr, &ret);
	ocl_check_error(ret, CL_SUCCESS);

	cl_float* input_output_host = (cl_float*)clEnqueueMapBuffer(cmd_queue, input_output_pinned, CL_FALSE, CL_MAP_WRITE | CL_MAP_READ, 0, sizeof(cl_float) * this->size(), 0, nullptr, nullptr, &ret);
	ocl_check_error(ret, CL_SUCCESS);

	cl_float* convolution_mask_host = (cl_float*)clEnqueueMapBuffer(cmd_queue, convol_mask_pinned, CL_FALSE, CL_MAP_READ, 0, sizeof(cl_float) * mask.size(), 0, nullptr, nullptr, &ret);
	ocl_check_error(ret, CL_SUCCESS);

	clFinish(cmd_queue);
	ocl_check_error(ret, CL_SUCCESS);

//	cl_float* input_output_host = new cl_float[this->size()];

	int i_pos = 0;
	for (int j = 0; j < this->height; j++) {
		for (int i = 0; i < this->width; i++) {
			input_output_host[i_pos++] = this->at(i,j);
		}
	}
//	cl_float* convolution_mask_host = new cl_float[mask.size()];

	for (int i = 0; i < mask.size(); i++) {
		convolution_mask_host[i] = mask[i];
	}

	if (cmd_queue_info.has_img_support == true && false) {
		// image supported - use convolution with image

		cl_kernel convolution_img_kernel = OpenCL::context_manager->getCommonKernel("image_convolve_img")[context_number];

		size_t max_kernel_work_group_size;
		cl_ulong kernel_local_mem_use;
		ret = clGetKernelWorkGroupInfo(convolution_img_kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &max_kernel_work_group_size, nullptr);
		ret |= clGetKernelWorkGroupInfo(convolution_img_kernel, device_id, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(cl_ulong), &kernel_local_mem_use, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		// use precalculated values for best work-group size
		int best_local_work_size_x = cmd_queue_info.best_work_groups_dim2.x;
		int best_local_work_size_y = cmd_queue_info.best_work_groups_dim2.y;

		// select best work-group size by decresing size until device can satisfy work group size requirements for this kernel
		while (best_local_work_size_x * best_local_work_size_y > max_kernel_work_group_size) {			
			// decrese one that is bigger so that we evenly distribute reduction
			if (best_local_work_size_x >= best_local_work_size_y) {
				best_local_work_size_x /= 2;
			} else {
				best_local_work_size_y /= 2;
			}
		}

		int best_group_size_x = (this->width) / best_local_work_size_x + 1;
		int best_group_size_y = (this->height) / best_local_work_size_y + 1;

		int best_global_work_size_x = best_group_size_x * best_local_work_size_x;
		int best_global_work_size_y = best_group_size_y * best_local_work_size_y;

		size_t global_work_load[2] = {best_global_work_size_x, best_global_work_size_y};
		size_t local_work_load[2] = {best_local_work_size_x, best_local_work_size_y};
		
		cl_image_format format;
		format.image_channel_order = CL_R;
		format.image_channel_data_type = CL_FLOAT;

		//cl_mem input_cl_img = clCreateImage2D(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, &format, this->width, this->height, 0, input_output_host, &ret);
		cl_mem input_cl_img = clCreateImage2D(context, CL_MEM_READ_ONLY, &format, this->width, this->height, 0, nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem output_cl_img = clCreateImage2D(context, CL_MEM_WRITE_ONLY , &format, this->width, this->height, 0, nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		//cl_mem convolution_mask = clCreateImage2D(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, &format, mask.width, mask.height, 0, convolution_mask_host, &ret);
		cl_mem convolution_mask = clCreateImage2D(context, CL_MEM_READ_ONLY, &format, mask.width, mask.height, 0, nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);
	
		ret = clSetKernelArg(convolution_img_kernel, 0, sizeof(cl_mem), (void*)&input_cl_img);
		ret |= clSetKernelArg(convolution_img_kernel, 1, sizeof(cl_int2), (void*)&image_size);
		ret |= clSetKernelArg(convolution_img_kernel, 2, sizeof(cl_mem), (void*)&convolution_mask);
		ret |= clSetKernelArg(convolution_img_kernel, 3, sizeof(cl_int4), (void*)&convol_mask_size);
		ret |= clSetKernelArg(convolution_img_kernel, 4, sizeof(cl_mem), (void*)&output_cl_img);
		ret |= clSetKernelArg(convolution_img_kernel, 5, sizeof(cl_mem), (void*)&kernel_err);

		cl_event evnt_copy[2];
		size_t origin[3] = {0,0,0};
		size_t region[3] = {width, height, 1};

		ret = clEnqueueWriteImage(cmd_queue, input_cl_img, CL_FALSE, origin, region, 0, 0, input_output_host, 0, nullptr, evnt_copy + 0);
		ocl_check_error(ret, CL_SUCCESS);

		region[0] = mask.width;
		region[1] = mask.height;
		ret = clEnqueueWriteImage(cmd_queue, convolution_mask, CL_FALSE, origin, region, 0, 0, convolution_mask_host, 0, nullptr, evnt_copy + 1);
		ocl_check_error(ret, CL_SUCCESS);

		cl_event evnt;
		ret = clEnqueueNDRangeKernel(cmd_queue, convolution_img_kernel, 2, nullptr, global_work_load, local_work_load, 2, evnt_copy, &evnt);
		ocl_check_error(ret, CL_SUCCESS);

		//size_t origin[3] = {0,0,0};
		//size_t region[3] = {this->width, this->height, 1};
		region[0] = width;
		region[1] = height;
		ret = clEnqueueReadImage(cmd_queue, output_cl_img, CL_FALSE, origin, region, 0, 0, input_output_host, 1, &evnt, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clEnqueueReadBuffer(cmd_queue, kernel_err, CL_FALSE, 0, sizeof(cl_int), &kernel_err_host, 1, &evnt, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clFinish(cmd_queue);
		ocl_check_error(ret, CL_SUCCESS);

	} else {

		// no image support - use normal version of convolution (will use local mem if device has local memory)

		// find out how big is local memory and use convolution with appropriate vector and group size

		int vector_size = 4;

		stringstream convol_kernel_name;
		convol_kernel_name << "image_convolve_vector" << vector_size << "";

		cl_kernel convolution_kernel = OpenCL::context_manager->getCommonKernel(convol_kernel_name.str())[context_number];

		// get memory usage and best work group size from opencl driver
		size_t max_kernel_work_group_size;
		cl_ulong kernel_local_mem_use;
		
		ret = clGetKernelWorkGroupInfo(convolution_kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &max_kernel_work_group_size, nullptr);
		ret |= clGetKernelWorkGroupInfo(convolution_kernel, device_id, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(cl_ulong), &kernel_local_mem_use, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		// use integer division in case if even size where convol_mask_size.s[0]-1 would not be ok
		int added_border_size = 2*(convol_mask_size.s[0]/2);

		// use only pre-calculated values by default (this are max values that we can use)
		// take into account number of compute uints that will work in parallel		
		int max_local_work_size_x_from_compute_units = this->width / (cmd_queue_info.max_compute_units * vector_size/2);

		frexp((float)max_local_work_size_x_from_compute_units, &max_local_work_size_x_from_compute_units);
		max_local_work_size_x_from_compute_units--;

		// use either default value (pre-computed max value) or max value based on number of computing units
		//int pow_x = min(cmd_queue_info.best_work_groups_dim2.pow_x, (int)::log((float)max_local_work_size_x_from_compute_units)/::log(2.0f));
		int pow_x = min(cmd_queue_info.best_work_groups_dim2.pow_x, max_local_work_size_x_from_compute_units);
		int pow_y = cmd_queue_info.best_work_groups_dim2.pow_y;

		int best_local_work_size_x = ::pow(2.0, pow_x);
		int best_local_work_size_y = cmd_queue_info.best_work_groups_dim2.y;

		// if will be using local memory then calculate best work-group size
		if (cmd_queue_info.local_mem_type == CL_LOCAL) {

			// calculate how much local memory we will need
			int memory_required = 4 * (best_local_work_size_x * vector_size + added_border_size) * (best_local_work_size_y + added_border_size) + kernel_local_mem_use;
		
			// select best work-group size by decresing size until device can satisfy local memory and work group size requirements for this kernel
			while (memory_required > cmd_queue_info.local_mem_size || best_local_work_size_x * best_local_work_size_y > max_kernel_work_group_size) {
				
				if (best_local_work_size_x >= best_local_work_size_y) {
					pow_x--;
					best_local_work_size_x = ::pow(2.0, pow_x);
				} else {
					pow_y--;
					best_local_work_size_y = ::pow(2.0, pow_y);
				}

				memory_required = 4 * (best_local_work_size_x * vector_size + added_border_size) * (best_local_work_size_y + added_border_size) + kernel_local_mem_use;
			}

			local_data_size.s[0] = best_local_work_size_x * vector_size + added_border_size;
			local_data_size.s[1] = best_local_work_size_y + added_border_size;
		} else {
			// else just use default values and
			// use min local data size since local mem will not be used

			local_data_size.s[0] = 1;
			local_data_size.s[1] = 1;
		}
		int best_group_size_x = ((this->width - added_border_size) / vector_size) / best_local_work_size_x + 1;
		int best_group_size_y = (this->height - added_border_size) / best_local_work_size_y + 1;

		int best_global_work_size_x = best_group_size_x * best_local_work_size_x;
		int best_global_work_size_y = best_group_size_y * best_local_work_size_y;

		size_t global_work_load[2] = {best_global_work_size_x, best_global_work_size_y};
		size_t local_work_load[2] = {best_local_work_size_x, best_local_work_size_y};

		//cl_mem input_image = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_float) * this->size(), input_output_host, &ret);
		cl_mem input_image = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float) * this->size(), nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		//cl_mem convolution_mask = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(cl_float) * mask.size(), convolution_mask_host, &ret);
		cl_mem convolution_mask = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(cl_float) * mask.size(), nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_mem output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(cl_float) * this->size(), nullptr, &ret);
		ocl_check_error(ret, CL_SUCCESS);

		cl_uint save_max = 0;
		cl_uint layer1_creator_type = 0; // 0 = normal convolution

		ret = clSetKernelArg(convolution_kernel, 0, sizeof(cl_mem), (void*)&input_image);
		ret |= clSetKernelArg(convolution_kernel, 1, sizeof(cl_int2), (void*)&image_size);
		ret |= clSetKernelArg(convolution_kernel, 2, sizeof(cl_mem), (void*)&convolution_mask);
		ret |= clSetKernelArg(convolution_kernel, 3, sizeof(cl_int4), (void*)&convol_mask_size);
		ret |= clSetKernelArg(convolution_kernel, 4, sizeof(cl_float) * local_data_size.s[0] * local_data_size.s[1] , nullptr);
		ret |= clSetKernelArg(convolution_kernel, 5, sizeof(cl_int2), (void*)&local_data_size);
		ret |= clSetKernelArg(convolution_kernel, 6, sizeof(cl_uint), (void*)&save_max);
		// since we do not need max values we can skip argument 7 and 8
		//ret |= clSetKernelArg(convolution_kernel, 7, sizeof(cl_mem), nullptr);
		//ret |= clSetKernelArg(convolution_kernel, 8, sizeof(cl_mem), nullptr);
		ret |= clSetKernelArg(convolution_kernel, 9, sizeof(cl_mem), (void*)&output);
		ret |= clSetKernelArg(convolution_kernel, 10, sizeof(cl_uint), (void*)&layer1_creator_type);		
		ret |= clSetKernelArg(convolution_kernel, 11, sizeof(cl_mem), (void*)&kernel_err);
		ocl_check_error(ret, CL_SUCCESS);

		cl_event evnt_copy[2];
		ret = clEnqueueWriteBuffer(cmd_queue, input_image, CL_FALSE, 0, sizeof(cl_float) * this->size(), input_output_host, 0, nullptr, evnt_copy + 0);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clEnqueueWriteBuffer(cmd_queue, convolution_mask, CL_FALSE, 0, sizeof(cl_float) * mask.size(), convolution_mask_host, 0, nullptr, evnt_copy + 1);
		ocl_check_error(ret, CL_SUCCESS);


		cl_event evnt;
		ret = clEnqueueNDRangeKernel(cmd_queue, convolution_kernel, 2, nullptr, global_work_load, local_work_load, 2, evnt_copy, &evnt);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clEnqueueReadBuffer(cmd_queue, output, CL_FALSE, 0, sizeof(cl_float) * this->size(), input_output_host, 1, &evnt, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clEnqueueReadBuffer(cmd_queue, kernel_err, CL_FALSE, 0, sizeof(cl_int), &kernel_err_host, 1, &evnt, nullptr);
		ocl_check_error(ret, CL_SUCCESS);

		ret = clFinish(cmd_queue);
		ocl_check_error(ret, CL_SUCCESS);

	}

	img* result = new img((int)width, (int)height, true);

	for (int j = mask.height/2; j < result->height - mask.height/2; j++) {
		for (int i = mask.width/2; i < result->width - mask.width/2; i++) {
			result->at(i,j) = input_output_host[j * result->width + i];
		}
	}

	return result;
#else
	cout << "OpenCL convolution requested but binary not build with OpenCL support!!" << endl;
	throw new std::exception();
#endif
}


img* img::get_convolve(img& mask, bool use_opencl)
{
	if (use_opencl == true && width * height > 150*150 && false){
		// use opencl if specified to use it and if size of image > 150*150
		cout << "using opencl convolution" << endl;
		return get_convolve_ocl(mask);
	} else {
		img* result = new img((int)width, (int)height, 0.0);
		int maskw = (int)mask.width;
		int maskw1 = maskw - 1;
		int maskh = (int)mask.height;
		int maskh1 = maskh - 1;
		int masksize = maskw*maskh;
		int masksize1 = masksize - 1;
		int wstep = (int)width - maskw1;
		int wskip = (int)width - maskw;
		int hstep = (int)height - maskh1;

		HOP_REAL* maskorig = mask.ptr(maskw1, maskh1);
		HOP_REAL** maskrgn = new HOP_REAL*[masksize];
		HOP_REAL* destptr = result->ptr(maskw/2, maskh/2);
		HOP_REAL* p;
		HOP_REAL sum;
		int i, j, k;
		
		p = ptr(0, 0);
		k = 0;
		for (j = maskh; j > 0; --j) {
			for (i = maskw; i > 0; --i)  
				maskrgn[k++] = p++;
			p += wskip;
		}
		for (j = hstep; j > 0 ; --j) {
			for (i = wstep; i > 0; --i) {
				sum = 0.0;
				p = maskorig;
				for (k = masksize1; k >= 0; --k) { 
					sum += (*(p--))*(*(maskrgn[k]++));
				}
				*(destptr++) = sum;
			}
			for (k = masksize1; k >= 0; --k) maskrgn[k] += maskw1;
			destptr += maskw1;
		}
		delete[] maskrgn;
		return result;
	}    
}

img* img::get_convolve(img& kernel, img& mask)
{
	img* result = new img((int)width, (int)height, 0.0);
	int kernelw = (int)kernel.width;
	int kernelw1 = kernelw - 1;
	int kernelh = (int)kernel.height;
	int kernelh1 = kernelh - 1;
	int kernelsize = kernelw*kernelh;
	int kernelsize1 = kernelsize - 1;
	int wstep = (int)width - kernelw1;
	int wskip = (int)width - kernelw;
	int hstep = (int)height - kernelh1;

	HOP_REAL* kernelorig = kernel.ptr(kernelw1, kernelh1);
	HOP_REAL** kernelrgn = new HOP_REAL*[kernelsize];
	HOP_REAL* destptr = result->ptr(kernelw/2, kernelh/2);
	HOP_REAL* p;
	HOP_REAL* p_mask;
	HOP_REAL sum;
	int i, j, k;
		
	p = ptr(0, 0);
	p_mask = mask.ptr(kernelw/2, kernelh/2);
	k = 0;
	for (j = kernelh; j > 0; --j) {
		for (i = kernelw; i > 0; --i) {
			kernelrgn[k++] = p++;
		}
		p += wskip;
	}
	int count_pos = 0;
	int count_neg = 0;
	for (j = hstep; j > 0 ; --j) {
		for (i = wstep; i > 0; --i) {
			sum = 0.0;
			if (*(p_mask++) != 0) {
				p = kernelorig;
				for (k = kernelsize1; k >= 0; --k) { 
					kernelrgn[k] += count_neg;
					sum += (*(p--))*(*(kernelrgn[k]));					
				}
				count_neg = 1;
				//count_pos++;
			} else {
				//for (k = kernelsize1; k >= 0; --k) {
				//	kernelrgn[k]++;
				//}
				count_neg++;
			}
			
			*(destptr++) = sum;
		}
		for (k = kernelsize1; k >= 0; --k) kernelrgn[k] += kernelw1;
		destptr += kernelw1;
		p_mask += kernelw1;
	}
	delete[] kernelrgn;
	return result;
}

void img::blur(int mask_size, double sigma, bool normalization /* = true */)
{
	int rows = height; 
	int cols = width;

	if (mask_size % 2 == 0) mask_size += 1;

	if (!grayscale) {
		cv::Mat image = to_cv_image(*this);
		cv::Mat blured_image;

		cv::GaussianBlur(image, blured_image, cv::Size(mask_size,  mask_size), sigma, sigma);
		*this = from_cv(blured_image);
	} else {
		int type = (sizeof(HOP_REAL) == sizeof(double)) ? CV_64F : CV_32F;
		void* data = (void*)this->ptr(0, 0);
		cv::Mat image = cv::Mat(rows, cols, type, data);

		// blur image (will be copied to new cv::Mat)
		cv::Mat blured_image;
		cv::GaussianBlur(image, blured_image, cv::Size(mask_size,  mask_size), sigma, sigma);

		// copy result to this image
		img::iterator iter = this->begin();

		cv::MatConstIterator_<HOP_REAL> pixels = blured_image.begin<HOP_REAL>();
		for (unsigned i = rows*cols; i > 0 ; --i) {
			*iter = (double)*pixels;
			++pixels; ++iter;
		}
		if (normalization) normalize(0.0, 1.0);
	}
}

//
//void img::print_matrix()
//{
//    cimg_forX(*this, x) {
//        cimg_forY(*this, y) {
//            std::cout << (*this)(x, y) << " ";
//        }
//        std::cout << std::endl;
//    }
//}

HOP_REAL img::p_get_2x2sum(img*& result)
{
    result = new img(width, height);
    
    HOP_REAL* ptrs[4];
    HOP_REAL* resptr = result->ptr(0, 1);
#if defined WIN32 | defined WIN64
    memcpy_s(result->ptr(0, 0), sizeof(HOP_REAL)*width, ptr(0, 0), sizeof(HOP_REAL)*width);
#else
    memcpy(result->ptr(0, 0), ptr(0, 0), sizeof(HOP_REAL)*width);
#endif
    ptrs[0] = ptr(0, 0);
    ptrs[1] = ptr(1, 0);
    ptrs[2] = ptr(0, 1);
    ptrs[3] = ptr(1, 1);

    HOP_REAL sum;
    HOP_REAL max = 0.0;
    int k;
    for (unsigned j = 0; j < height - 1; ++j) {
        *(resptr++) = *ptrs[2];
        for (unsigned i = 0; i < width - 1; ++i) {
            for (k = 0, sum = 0.0; k < 4; ++k) sum += *ptrs[k];
            if (sum > max) max = sum;
            *resptr = sum;
            for (k = 0; k < 4; ++k) ++(ptrs[k]);
            ++resptr;
        }
        for (k = 0; k < 4; ++k) ++(ptrs[k]);
    }
    return max;       
}


void img::read_grayscale(const string& name, img& result)
{

	const cv::Mat image = cv::imread(name, CV_LOAD_IMAGE_GRAYSCALE);

	if (image.empty()) {
		result.resize(0, 0);
		return;
	}

	unsigned rows = image.rows, columns = image.cols;
	result.resize((size_t)columns, (size_t)rows);
	result.grayscale = true;

	img::iterator iter = result.begin();

	int type = image.type();

	
	cv::MatConstIterator_<unsigned char> pixels = image.begin<unsigned char>();
	for (unsigned i = rows*columns; i > 0 ; --i) {
		*iter = (double)*pixels;
		++pixels; ++iter;
	}
	result.normalize(0.0, 1.0);
}

void img::read_colors(const string& name, img& result, bool as_hsv /* = false */)
{
	cv::Mat image = cv::imread(name, CV_LOAD_IMAGE_COLOR);

	if (as_hsv == true) {
		cv::Mat new_image;
		cv::cvtColor(image, new_image, CV_BGR2HSV);
		image = new_image;
	}

	unsigned int columns = image.cols, rows = image.rows;

	result.resize((size_t)columns, (size_t)rows);
	result.grayscale = false;
	
	// data is encoded in 3 x 8 byte fields where each field is one channel
	cv::MatConstIterator_<cv::Vec3b> pixels = image.begin<cv::Vec3b>();
    img::iterator riter = result.begin();
    
    for (unsigned i = rows*columns; i > 0; --i) {
        *riter = COL_TO_REAL((int)(*pixels)[2], (int)(*pixels)[1], (int)(*pixels)[0]); 
		++pixels;
        ++riter;
    }
}

void img::read_colors(const string& name, img& red, img& green, img& blue, bool as_hsv)
{

	cv::Mat image = cv::imread(name, CV_LOAD_IMAGE_COLOR);

	if (as_hsv == true) {
		cv::Mat new_image;
		cv::cvtColor(image, new_image, CV_BGR2HSV);
		image = new_image;
	}

	unsigned int columns = image.cols, rows = image.rows;

    red.resize((size_t)columns, (size_t)rows);
    green.resize((size_t)columns, (size_t)rows);
    blue.resize((size_t)columns, (size_t)rows);
	
	// data is encoded in 3 x 8 byte fields where each field is one channel
	cv::MatConstIterator_<cv::Vec3b> pixels = image.begin<cv::Vec3b>();
    img::iterator riter = red.begin();
    img::iterator giter = green.begin();
    img::iterator biter = blue.begin();
    
    for (unsigned i = rows*columns; i > 0; --i) {
        *riter = (double)(*pixels)[2]; // blue and red channels are switched
        *giter = (double)(*pixels)[1];
        *biter = (double)(*pixels)[0];

		++pixels;
        ++riter; ++giter; ++biter;
    }
    red.normalize(0.0, 1.0);
    green.normalize(0.0, 1.0);
    blue.normalize(0.0, 1.0);
}


//img* img::get_grayscale()
//{
//    img* result = new img(width, height, depth, 1);
//    cimg_forXYZ(*this, x, y, z) {
//        const HOP_REAL 
//            R = (*this)(x, y, z, 0),
//            G = (*this)(x, y, z, 1),
//            B = (*this)(x, y, z, 2);
//        (*result)(x, y, z) = 0.299*R + 0.587*G + 0.114*B;
//    }
//    return result;
//}
//
//
//iimg* img::get_01_img(HOP_REAL thresh)
//{
//    iimg* result = new iimg(width, height);
//
//    int* rptr = result->data;
//    cimg_for(*this, ptr, HOP_REAL) {
//        *(rptr++) = *ptr < thresh ? 0 : 1;
//    }
//    return result;
//}

matrix<bool>* img::get_bool_matrix(HOP_REAL thresh)
{
    matrix<bool>* result = new matrix<bool>(width, height);

    std::vector<bool>::iterator iter = result->begin();
    iterator iter2 = begin();
    img_for(*this, ptr) {
        *(iter++) = *ptr >= thresh;
    }
    return result;
}

img& img::normalize(HOP_REAL min, HOP_REAL max)
{
    HOP_REAL minval, maxval;
    HOP_REAL k, n;
   
    minmax(minval, maxval);
    if (minval == maxval) return *this;
    k = (max - min)/(maxval - minval);
    n = min - minval*k;

    img_for(*this, ptr) {
        *ptr = (*ptr)*k + n;
    }
    return *this;
}

img img::get_resized(int neww, int newh, bool normalized /* = true */) const
{
    img result;

    get_resized(result, neww, newh, normalized);
    return result;
}

void img::get_resized(img& result, int newsize, bool normalized /* = true */) const
{
    if (newsize < 0) 
        get_resized(result, newsize, newsize, normalized);
    else {
        int neww, newh;

        if (width >= height) { // "landscape"
            neww = newsize;
            newh = (newsize*height)/width;
        } else {
            newh = newsize;
            neww = (newsize*width)/height;
        }
        get_resized(result, neww, newh, normalized);
    }
}

void img::get_resized(img& result, int neww, int newh, bool normalized /* = true */) const
{
    if (width == 0 || height == 0) return;

    if (neww < 0) neww = (int)(-(double)neww*width/100.0);
    if (newh < 0) newh = (int)(-(double)newh*height/100.0);

	// create space for new resized image

	if (grayscale) {
		result.resize((size_t)neww, (size_t)newh);
		result.grayscale = true;

		// make cv::Mat header for src
		int rows = height; 
		int cols = width;
		int type = (sizeof(HOP_REAL) == sizeof(double)) ? CV_64F : CV_32F;
		void* data = (void*)this->ptr(0, 0);
		cv::Mat image = cv::Mat(rows, cols, type, data);

		// make cv::Mat header for dst
		void* result_data = (void*)result.ptr(0, 0);
		cv::Mat resized_image(newh, neww, type, result_data);

		cv::Mat image_float;
		cv::Mat resized_image_float;

		// in case we have 64 byte we need to convert to 32
		if (type != CV_32F) {
			image.convertTo(image_float, CV_32F);		
		} else {
			resized_image_float = resized_image;
		}

		// resize image		
		cv::resize(image_float, resized_image_float, cv::Size(neww,newh), CV_INTER_NN);

		if (type != CV_32F) {
			// we neet to resize back to 64 bytes
			resized_image_float.convertTo(resized_image, type);
		}

		if (normalized) result.normalize(0.0, 1.0);
	} else {
		cv::Mat cvresult;
		cv::Mat m = to_cv_image(*this);

		//cv::imwrite("c:\\work\\dbg2.png", m);

		cv::resize(m, cvresult, cv::Size(neww, newh), CV_INTER_NN);
		
		result = from_cv(cvresult);
	}
}

HOP_REAL img::norm()
{
	int type = (sizeof(HOP_REAL) == sizeof(double)) ? CV_64F : CV_32F;
	void* data = (void*)this->ptr(0, 0);
	cv::Mat image = cv::Mat(height, width, type, data);
	cv::Mat eval;
	double minev, maxev;

	cv::SVD::compute(image*image.t(), eval);
	cv::minMaxIdx(eval, &minev, &maxev);

	return (HOP_REAL)::sqrt(maxev);
}

void img::crop(int w, int h)
{
    if (w == (int)width && h == (int)height) return;

    int neww = min((int)width, w);
    int newh = min((int)height, h);
    img temp(*this);

    resize(neww, newh);
    for_each_xy (*this, i, j) {
        at(i, j) = temp.at(i, j);
    }
}

img img::cut(const irectangle2& rect) const
{
    if (empty()) return img();

    irectangle2 rect2 = rect.intersection(irectangle2(0, 0, width, height));
    img result(rect.x_dim(), rect.y_dim(), 0.0, grayscale);
    
    for (int i = rect2.ll.x, ii = 0; i < rect2.ur.x; ++i, ++ii) {
        for (int j = rect2.ll.y, jj = 0; j < rect2.ur.y; ++j, ++jj) {
            result(ii, jj) = at(i, j);
        }
    }
    return result;
}

void img::flip_horizontal()
{
    size_t maxi = width/2;

    for (size_t j = 0; j < height; ++j) {
        double* left = &at((size_t)0, j);
        double* right = &at(width - 1, j);

        for (size_t i = 0; i < maxi; ++i) {
            double tmp = *left;

            *left = *right;
            *right = tmp;
            ++left; 
            --right;
        }
    }
}

void img::save(const string& name, int newxsize /* = -100 */) const
{
    if (grayscale) save_grayscale(name, newxsize);
    else save_colored(name, newxsize);
}

unsigned int img::save(void*& mem, int newxsize /*= -100*/) const
{
    if (grayscale) return save_grayscale(mem, newxsize);
    else return 0;
}

void img::save_grayscale(const string& name, int newxsize) const
{
	
    if (width == 0 || height == 0) 
        return;

	img tmp(*this);
	tmp.normalize(0.0, 1.0);

	int rows = height; 
	int cols = width;
	int type = (sizeof(HOP_REAL) == sizeof(double)) ? CV_64F : CV_32F;
	void* data = (void*)tmp.ptr(0, 0);		
	cv::Mat image = cv::Mat(rows, cols, type, data);

	cv::Mat gray_image;
	image.convertTo(gray_image, CV_8U, 255);
	image = gray_image;

	if (newxsize != -100) {
		unsigned newwidth, newheight;

		if (newxsize < 0) {
			newwidth = (unsigned)((-newxsize)*width/100.0);
			newheight = (unsigned)((-newxsize)*height/100.0);
		} else {
			 newwidth = (unsigned)newxsize;
			 newheight = (unsigned)(((double)newwidth/width)*height);
		}

		cv::Mat resized_img;
		cv::resize(image, resized_img, cv::Size(newwidth, newheight));		
		image = resized_img;
	}		
	cv::imwrite(name, image);	
}

unsigned int img::save_grayscale(void*& mem, int newxsize) const
{
    if (width == 0 || height == 0) 
        return 0;

    img tmp(*this);
	tmp.normalize(0.0, 1.0);

	int rows = height; 
	int cols = width;
	int type = (sizeof(HOP_REAL) == sizeof(double)) ? CV_64F : CV_32F;
	void* data = (void*)tmp.ptr(0, 0);		
	cv::Mat image = cv::Mat(rows, cols, type, data);

	// we need to convert it to grayscale first
	cv::Mat gray_image;
	image.convertTo(gray_image, CV_8U, 255);
	image = gray_image;

	if (newxsize != -100) {
		unsigned newwidth, newheight;

		if (newxsize < 0) {
			newwidth = (unsigned)((-newxsize)*width/100.0);
			newheight = (unsigned)((-newxsize)*height/100.0);
		} else {
			 newwidth = (unsigned)newxsize;
			 newheight = (unsigned)(((double)newwidth/width)*height);
		}

		cv::Mat resized_img;
		cv::resize(image, resized_img, cv::Size(newwidth, newheight));		
		image = resized_img;
	}		


	vector<uchar> buf;
	cv::imencode("*.bmp", image, buf);

	mem = new uchar[buf.size()];
	memcpy(mem, &buf[0], buf.size());

	return buf.size();	
}

void img::save_colored(const string& name, int newxsize) const
{
    if (width == 0 || height == 0) 
        return;

    int size = (int)(width*height);
    IMG_COLOR4* pixels = new IMG_COLOR4[size];
    const IMG_COLOR* data = (IMG_COLOR*)ptr(0, 0);

    for (int i = 0; i < size; ++i) {
        pixels[i].blue = ((IMG_COLOR3*)data)->color.red;  // red <--> blue!
		pixels[i].green = ((IMG_COLOR3*)data)->color.green;
		pixels[i].red = ((IMG_COLOR3*)data)->color.blue;  // red <--> blue!
		data++;
    }

	int rows = height; 
	int cols = width;
	cv::Mat image = cv::Mat(rows, cols, CV_8UC3, pixels);

	if (newxsize != -100) {
		unsigned newwidth, newheight;

		if (newxsize < 0) {
			newwidth = (unsigned)((-newxsize)*width/100.0);
			newheight = (unsigned)((-newxsize)*height/100.0);
		} else {
			 newwidth = (unsigned)newxsize;
			 newheight = (unsigned)(((double)newwidth/width)*height);
		}

		cv::Mat resized_img;
		cv::resize(image, resized_img, cv::Size(newwidth, newheight));		
		image = resized_img;
	}
	cv::imwrite(name, image);	
    delete pixels;
}

void img::save_visview(const char* name)
{
    img result(*this);
    ofstream os;

    os.open(name, ios::trunc);
    if (os.fail()) return;

    result.normalize(0.0, 255.0);
    int height1 = (int)height - 1;
    for (int x = 0; x < (int)width; ++x) {
        for (int y = 0; y < height1; ++y) {
            os << std::max<int>(1, (int)(result(x, y) + 0.5)) << ", ";
        }
        os << std::max<int>(1, (int)(result(x, height1) + 0.5)) << endl;
    }
    os.close();
}

void img::save_jet_colormap(const char* name, int sizex)
{
    img ni(*this);
    
    ni.normalize(0.0, 255.0);

    img result(width, height, false); // color image
    img jmap = jet_map();
    int i;

    for_each_xy(ni, x, y) {
        i = (int)(ni(x, y));
        result.set_color((int)x, (int)y, jmap(i, 0), jmap(i, 1), jmap(i, 2));
    }
    result.save(name, sizex);
}

img& img::combine_max(const img& src, int desti, int destj, HOP_REAL factor)
{
    int i0 = max(0, desti) - desti;
    int j0 = max(0, destj) - destj;
    int i1 = min((int)width, desti + (int)src.width) - desti;
    int j1 = min((int)height, destj + (int)src.height) - destj;
    
    if (grayscale) {
        HOP_REAL* p;
    
        for (int i = i0; i < i1; ++i) {
            for (int j = j0; j < j1; ++j) {
                HOP_REAL val = factor*src(i, j);

                p = ptr(i + desti, j + destj);
                *p = std::max<double>(*p, val);
            }
        }
    } else {
        for (int i = i0; i < i1; ++i) {
            for (int j = j0; j < j1; ++j) {
                //const IMG_COLOR& srccol = (IMG_COLOR&)src.color_at(i, j);
                HOP_REAL srccol = 255.0 * src.at(i, j);
                IMG_COLOR& destcol = (IMG_COLOR&)at(i + desti, j + destj);
                
                destcol.red = std::max<int>(std::min<int>(COL_TYPE_MAX, (int)(srccol*factor)), destcol.red);
                destcol.green = std::max<int>(std::min<int>(COL_TYPE_MAX, (int)(srccol*factor)), destcol.green);
                destcol.blue = std::max<int>(std::min<int>(COL_TYPE_MAX, (int)(srccol*factor)), destcol.blue);
            }
        }
                    
    }
    return *this;
}

img& img::combine_sum(const img& src, int desti, int destj, HOP_REAL factor)
{
    int i0 = max(0, desti) - desti;
    int j0 = max(0, destj) - destj;
    int i1 = min((int)width, desti + (int)src.width) - desti;
    int j1 = min((int)height, destj + (int)src.height) - destj;
    
    if (grayscale) {
        HOP_REAL* p;
    
        for (int i = i0; i < i1; ++i) {
            for (int j = j0; j < j1; ++j) {
                HOP_REAL val = factor*src(i, j);

                p = ptr(i + desti, j + destj);
                //*p = std::min<double>(1.0, *p + val);
                *p += val;
            }
        }
    } else {
        for (int i = i0; i < i1; ++i) {
            for (int j = j0; j < j1; ++j) {
                //const IMG_COLOR& srccol = (IMG_COLOR&)src.color_at(i, j);
                HOP_REAL srccol = 255.0 * src.at(i, j);
                IMG_COLOR& destcol = (IMG_COLOR&)at(i + desti, j + destj);
                
                destcol.red = std::max<int>(std::min<int>(COL_TYPE_MAX, (int)(srccol*factor)), destcol.red);
                destcol.green = std::max<int>(std::min<int>(COL_TYPE_MAX, (int)(srccol*factor)), destcol.green);
                destcol.blue = std::max<int>(std::min<int>(COL_TYPE_MAX, (int)(srccol*factor)), destcol.blue);
            }
        }
                    
    }
    return *this;
}

//img& img::combine_sum(const img& src, int desti, int destj, HOP_REAL factor)
//{
//    if (desti < 0 || destj < 0 || 
//        desti + src.width > width ||
//        destj + src.height > height) return *this;
//
//    HOP_REAL* p;
//
//    for (unsigned i = 0; i < src.width; ++i) {
//        for (unsigned j = 0; j < src.height; ++j) {
//            p = ptr(i + desti, j + destj);
//            *p += factor*src(i, j);
//        }
//    }
//    return *this;
//}
//
//img& img::draw_frame(int x0, int y0, int x1, int y1, HOP_REAL color)
//{
//    for (int i = x0; i < x1; ++i) { at(i, y0) = color; at(i, y1) = color; }
//    for (int j = y0; j < y1; ++j) { at(x0, j) = color; at(x1, j) = color; }
//    return *this;
//}

// Copies src image to the picture.
// Point (cx, cy) on the image is mapped to (cdx, cdy) on this image.
img& img::blt_central_max(const img& src, int csx, int csy, int cdx, int cdy)
{
    int minx = cdx - csx, miny = cdy - csy;
    int maxdx = minx + (int)src.width, maxdy = miny + (int)src.height;
    int mindx = max(minx, 0);
    int mindy = max(miny, 0);
    int minsx = mindx - minx, minsy = mindy - miny;
    maxdx = min(maxdx, (int)width);
    maxdy = min(maxdy, (int)height);

    for (int x = mindx, sx = minsx; x < maxdx; ++x, ++sx)
        for (int y = mindy, sy = minsy; y < maxdy; ++y, ++sy) {
            HOP_REAL& pt = at(x, y);
            pt = max(pt, src(sx, sy));
        }
    return *this;
}

// Copies src image to the picture.
// Point (cx, cy) on the image is mapped to (cdx, cdy) on this image.
// Returns the size of the x- and y-offsets of the "old" image in the "new" image
// Remark: this function automatically resizes the image.
ipoint2 img::blt_central_max_ar(const img& src, int csx, int csy, int cdx, int cdy)
{
    int mindx = cdx - csx, mindy = cdy - csy;
    int offx, offy;
    
    if (mindx >= 0) offx = 0; else { offx = -mindx; mindx = 0; } 
    if (mindy >= 0) offy = 0; else { offy = -mindy; mindy = 0; }
    
    int maxdx = mindx + (int)src.width, maxdy = mindy + (int)src.height;
    int plusx = max(0, maxdx - ((int)width + offx - 1));
    int plusy = max(0, maxdy - ((int)height + offy - 1));

    if (offx > 0 || offy > 0 || plusx > 0 || plusy > 0) {
        img tmp(*this);

        resize(width + offx + plusx, height + offy + plusy, 0);
        fill(0);
        blt(tmp, offx, offy);
        //save("c:\\tmpx.bmp");
    }

    for (int x = mindx, sx = 0; x < maxdx; ++x, ++sx)
        for (int y = mindy, sy = 0; y < maxdy; ++y, ++sy) {
            HOP_REAL& pt = at(x, y);
            pt = max(pt, src(sx, sy));
        }
    return ipoint2(offx, offy);
}

//img& img::blt_central_max(const img& srcn, int csx, int csy, double factor, int cdx, int cdy)
//{
//    CImg<HOP_REAL> src = srcn.get_resize(floor(factor*srcn.width), floor(factor*srcn.height));
//    csx = floor(factor*csx);
//    csy = floor(factor*csy);
//
//    int minx = cdx - csx, miny = cdy - csy;
//    int maxdx = minx + src.width, maxdy = miny + src.height;
//    int mindx = std::max(minx, 0);
//    int mindy = std::max(miny, 0);
//    int minsx = mindx - minx, minsy = mindy - miny;
//    maxdx = std::min(maxdx, (int)width);
//    maxdy = std::min(maxdy, (int)height);
//    HOP_REAL* pt;
//
//    for (int x = mindx, sx = minsx; x < maxdx; ++x, ++sx)
//        for (int y = mindy, sy = minsy; y < maxdy; ++y, ++sy) {
//            pt = ptr(x, y);
//            *pt = std::max(*pt, src(sx, sy));
//        }
//    return *this;
//}
//
//
//img& img::blt_max(const img& src, int dx, int dy)
//{
//    int mindx = std::max(dx, 0);
//    int mindy = std::max(dy, 0);
//    int maxdx = std::min((int)width, mindx + (int)src.width);
//    int maxdy = std::min((int)height, mindy + (int)src.height);
//    int minsx = mindx - dx, minsy = mindy - dy;
//
//    for (int x = mindx, sx = minsx; x < maxdx; ++x, ++sx)
//        for (int y = mindy, sy = minsy; y < maxdy; ++y, ++sy) {
//            (*this)(x, y) += src(sx, sy);
//        }
//    return *this;
//}


img img::concat_linear(const std::vector<img*>& images)
{
    unsigned w = 0, maxh = 0;
    std::vector<img*>::const_iterator iter;

    for (iter = images.begin(); iter != images.end(); ++iter) {
        img* im = *iter;

        w += im->width;
        if (maxh < im->height) maxh = im->height;
    }

    img result((size_t)w, (size_t)maxh, images[0]->grayscale);

    w = 0;
    for (iter = images.begin(); iter != images.end(); ++iter) {
        img* im = *iter;

        result.blt(*im, (int)w, (int)(maxh - im->height));
        w += im->width;
    }
    return result;
}

img img::concat_linear(std::vector<img>& images)
{
    std::vector<img*> tmp;

    for (size_t i = 0; i < images.size(); ++i) 
        tmp.push_back(&images[i]);
    return concat_linear(tmp);
}

img img::concat_linear_fw(const std::vector<img*>& images)
{
    if (images.empty())
        return img();

    unsigned maxw = 0, maxh = 0;
    std::vector<img*>::const_iterator iter;

    for (iter = images.begin(); iter != images.end(); ++iter) {
        img* im = *iter;

        if (maxw < im->width) maxw = im->width;
        if (maxh < im->height) maxh = im->height;
    }

    img result((size_t)(maxw*images.size()), (size_t)maxh, images[0]->grayscale);
    int w = 0;

    for (iter = images.begin(); iter != images.end(); ++iter) {
        img* im = *iter;

        result.blt_central(*im, (int)(im->width/2), (int)(im->height/2), w + maxw/2, maxh/2);
        w += maxw;
    }
    return result;
}

img img::concat_linear_fw(std::vector<img>& images)
{
    std::vector<img*> tmp;

    for (size_t i = 0; i < images.size(); ++i) 
        tmp.push_back(&images[i]);
    return concat_linear_fw(tmp);
}

img img::concat(const std::vector<img*>& images, double zero /* = 0.0 */)
{
    size_t size = images.size();

    if (size == 0) return img();
    
    int ynum = (int)::sqrt((double)size/1.62) + 1;
    int xnum = (int)(ynum*1.62);
    int maxw = 0, maxh = 0;

    for (size_t i = 0; i < size; ++i) {
        double min, max;
        if (maxw < (int)images[i]->width) maxw = (int)images[i]->width;
        if (maxh < (int)images[i]->height) maxh = (int)images[i]->height;
    }

    int maxw2 = maxw/2, maxh2 = maxh/2;

    img result(maxw*xnum, maxh*ynum, zero, images[0]->grayscale);
    if (result.width == 0 || result.height == 0) 
        return result;
    for (size_t i = 0; i < size; ++i) {
        img& im = *images[i];

        result.blt_central(im, (int)im.width/2, (int)im.height/2, 
            ((int)i % xnum)*maxw + maxw2, ((int)i/xnum)*maxh + maxh2);
    }
    return result;
}

img img::concat(std::vector<img>& images, double zero /* = 0.0 */)
{
    std::vector<img*> tmp;

    for (size_t i = 0; i < images.size(); ++i) 
        tmp.push_back(&images[i]);
    return concat(tmp, zero);
}

img img::string_to_img(const std::string& str)
{
    vector<img*> chars(str.length(), nullptr);

    for (unsigned i = 0; i < str.length(); ++i) {
        chars[i] = get_char(str[i]);
    }
    return concat_linear(chars);
}

img img::number_to_img(int i)
{
    char str[100];

    sprintf(str, "%d", i);
    return string_to_img(str);
}

img& img::positive()
{
    img_for(*this, ptr) { 
        if (*ptr < 0.0) *ptr = 0.0;
    }
    return *this;
}

img& img::negative()
{
    img_for(*this, ptr) { 
        if (*ptr > 0.0) *ptr = 0.0;
    }
    return *this;
}

img& img::negative2()
{
    img_for(*this, ptr) { 
        if (*ptr > 0.0) *ptr = 0.0; else *ptr = -(*ptr);
    }
    return *this;
}

img& img::neg()
{
    img_for(*this, ptr) { 
        *ptr = -(*ptr);
    }
    return *this;
}



void img::write_to_stream(ostreamer& os) const
{
    os.write((unsigned)width);
    os.write((unsigned)height);
    os.write((unsigned)1);  // for "historical" reasons - depth
    os.write((unsigned)1);  // for "historical" reasons - dim
    img_for_const(*this, ptr) {
        os.write(*ptr);
    }
}


void img::read_from_stream(istreamer& is)
{
    unsigned width, height, dummy;
    is.read(width);
    is.read(height);
    is.read(dummy); // depth
    is.read(dummy); // dim
    
    resize((size_t)width, (size_t)height);

    img_for(*this, ptr) {
        is.read(*ptr);
    }
}
//
//// interp (see cimg.h):
////  - 0 = no interpolation : additional space is filled with 0.
////  - 1 = bloc interpolation (nearest point).
////  - 2 = mosaic : image is repeated if necessary.
////  - 3 = linear interpolation.
////  - 4 = grid interpolation.
////  - 5 = bi-cubic interpolation.
//void img::blur_resize(HOP_REAL factor, int startsize, int slimito, int msize, double sigma, 
//                      int interp, std::vector<img*>& result)
//{
//    img* mask = img::gaussian_mask(msize, msize, sigma);
//    unsigned size = std::min(width, height);
//    //char fname[100];
//    int slimit;
//    bool and;
//
//    if (slimito < 0) { slimit = -slimito; and = false; } else { slimit = slimito; and = true; } 
//    result.clear();
//    if (startsize <= slimit) {
//        result.push_back(new img(*this));
//    } else {
//        double initfactor = (width >= height) ? (double)startsize/width : (double)startsize/height;
//        result.push_back(this->p_get_convolve(mask));
//        result.back()->resize((int)(width*initfactor), (int)(height*initfactor), -100, -100, interp);
//        //sprintf(fname, "%d.bmp", 0);
//        //result.back()->save_normalized(fname);
//    }
//    if (slimit < 0 || factor >= 1.0) return;
//
//    unsigned limit = (unsigned)(slimit/factor);
//    unsigned lastw, lasth;
//    if (and) {
//		lastw = result.back()->width;
//		lasth = result.back()->height;
//        int index = 1;
//        while ((int)(lastw*factor) > limit && (int)(lasth*factor) > limit) {
//            result.push_back(result.back()->p_get_convolve(mask));
//            result.back()->resize((int)(lastw*factor), (int)(lasth*factor), -100, -100, interp);
//            //sprintf(fname, "%d.bmp", index++);
//            //result.back()->save_normalized(fname);
//			lastw = result.back()->width;
//			lasth = result.back()->height;
//        }
//    } else {
//		lastw = result.back()->width;
//		lasth = result.back()->height;
//        int index = 1;
//        while ((int)(lastw*factor) > limit || (int)(lasth*factor) > limit) {
//            result.push_back(result.back()->p_get_convolve(mask));
//            result.back()->resize((int)(lastw*factor), (int)(lasth*factor), -100, -100, interp);
//            //sprintf(fname, "%d.bmp", index++);
//            //result.back()->save_normalized(fname);
//			lastw = result.back()->width;
//			lasth = result.back()->height;
//        }
//    }
//    delete mask;
//}
//
//double img::convolve_central(const img& src, int csx, int csy, int cdx, int cdy)
//{
//    int minx = cdx - csx, miny = cdy - csy;
//    int maxdx = minx + src.width, maxdy = miny + src.height;
//    int mindx = std::max(minx, 0);
//    int mindy = std::max(miny, 0);
//    int minsx = mindx - minx, minsy = mindy - miny;
//    maxdx = std::min(maxdx, (int)width);
//    maxdy = std::min(maxdy, (int)height);
//    double result = 0.0;
//    double d1, d2;
//    double n1 = 0.0, n2 = 0.0;
//
//    for (int x = mindx, sx = minsx; x < maxdx; ++x, ++sx)
//        for (int y = mindy, sy = minsy; y < maxdy; ++y, ++sy) {
//            d1 = (*this)(x, y); d2 = src(sx, sy);
//            n1 += d1*d1; n2 += d2*d2;
//            result += d1*d2;
//        }
//    return result/::sqrt(n1)/::sqrt(n2);
//}

void img::draw_line(int x1, int y1, int x2, int y2, HOP_REAL color)
{
    // Uses Bresenham algorithm:
    // http://www.gamedev.net/reference/articles/article1275.asp
    int deltax = abs(x2 - x1);        // The difference between the x's
    int deltay = abs(y2 - y1);        // The difference between the y's
    int x = x1;                       // Start x off at the first pixel
    int y = y1;                       // Start y off at the first pixel
    int xinc1, xinc2, yinc1, yinc2;
    int den, num;
    int numadd, numpixels;

    if (x2 >= x1) { xinc1 = 1; xinc2 = 1; }     // The x-values are increasing
    else { xinc1 = -1; xinc2 = -1; }            // The x-values are decreasing

    if (y2 >= y1) { yinc1 = 1; yinc2 = 1; }     // The y-values are increasing
    else { yinc1 = -1; yinc2 = -1; }            // The y-values are decreasing

    if (deltax >= deltay) {         // There is at least one x-value for every y-value
        xinc1 = 0;                  // Don't change the x when numerator >= denominator
        yinc2 = 0;                  // Don't change the y for every iteration
        den = deltax;
        num = deltax / 2;
        numadd = deltay;
        numpixels = deltax;         // There are more x-values than y-values
    } else {                        // There is at least one y-value for every x-value
        xinc2 = 0;                  // Don't change the x for every iteration
        yinc1 = 0;                  // Don't change the y when numerator >= denominator
        den = deltay;
        num = deltay / 2;
        numadd = deltax;
        numpixels = deltay;         // There are more y-values than x-values
    }

    for (int curpixel = 0; curpixel <= numpixels; curpixel++) {
        set_at_safe(x, y, color);   // Draw the current pixel
        num += numadd;              // Increase the numerator by the top of the fraction
        if (num >= den) {           // Check if numerator >= denominator
            num -= den;             // Calculate the new numerator value
            x += xinc1;             // Change the x as appropriate
            y += yinc1;             // Change the y as appropriate
        }
        x += xinc2;                 // Change the x as appropriate
        y += yinc2;                 // Change the y as appropriate
    }
}    

void img::draw_box(const rectangle2<int>& rect, const color& c)
{
    draw_box(rect, COL_TO_REAL(c.red(), c.green(), c.blue()));
}

void img::draw_box(const rectangle2<int>& rect, HOP_REAL col)
{
    draw_line(rect.ll.x, rect.ll.y, rect.ur.x, rect.ll.y, col);
    draw_line(rect.ur.x, rect.ll.y, rect.ur.x, rect.ur.y, col);
    draw_line(rect.ur.x, rect.ur.y, rect.ll.x, rect.ur.y, col);
    draw_line(rect.ll.x, rect.ur.y, rect.ll.x, rect.ll.y, col);
}

void img::draw_rectangle(const rectangle2<int>& rect, HOP_REAL c)
{
    if (rect.invalid()) 
        return;

    int i0 = max(rect.ll.x, 0), j0 = max(rect.ll.y, 0);
    int i1 = min<int>(rect.ur.x + 1, width), j1 = min<int>(rect.ur.y + 1, height);

    for (int i = i0; i < i1; ++i)
        for (int j = j0; j < j1; ++j) 
            at(i, j) = c;
}

void img::draw_box(const rectangle2<int>& rect, HOP_REAL col, int w)
{
    for (int i = -w/2; i < (w + 1)/2; ++i) {
        rectangle2<int> r = rect;

        r.grow(i);
        draw_line(r.ll.x, r.ll.y, r.ur.x, r.ll.y, col);
        draw_line(r.ur.x, r.ll.y, r.ur.x, r.ur.y, col);
        draw_line(r.ur.x, r.ur.y, r.ll.x, r.ur.y, col);
        draw_line(r.ll.x, r.ur.y, r.ll.x, r.ll.y, col);
    }
}

void img::draw_big_point(int x, int y, HOP_REAL* col)
{
    at(x, y) = *col;
}

void img::draw_big_point(int x, int y, HOP_REAL col)
{
    for (int i = -1; i < 2; ++i) 
        for (int j = -1; j < 2; ++j) at(x + i, y + j) = col;
}

void img::revert_white(double white /* = 1.0*/)
{
    img_for(*this, ptr) { 
        if (*ptr >= white) *ptr = 0.0; else *ptr = 1.0;
    }
}

void img::replace_color(HOP_REAL what, HOP_REAL with)
{
    img_for(*this, iter) {
        if (*iter == what) *iter = with;
    }
}

img img::to_color_img(HOP_REAL color) const
{
    if (!grayscale) return *this;

    img result(width, height, 0.0, false);
    
    for_each_xy(*this, i, j) {
        double d = at(i, j);

        if (d > 0.0) result(i, j) = color;
    }
    return result;
}

img img::to_color_img() const
{
    if (!grayscale) return *this;

    img result(width, height, 0.0, false);
    
    for_each_xy(*this, i, j) {
        double d = at(i, j);
        int col = max(0, min(255, (int)(d*255)));

        if (d > 0.0) result(i, j) = COL_TO_REAL(col, col, col);
    }
    return result;
}

void img::thicken_white(int r, double white /* = 1.0*/)
{
    img temp(*this);

    for (int i = 0; i < (int)width; ++i) {
        for (int j = 0; j < (int)height; ++j) {
            if (temp.at(i, j) < white) continue;

            int i0 = max(0, i - r), i1 = min((int)width, i + r + 1);
            int j0 = max(0, j - r), j1 = min((int)height, j + r + 1);

            for (int ii = i0; ii < i1; ++ii) {
                for (int jj = j0; jj < j1; ++jj) {
                    at(ii, jj) = 1.0;
                }
            }
        }
    }
}

//void img::draw_boxes(const list<rectangle2<int> >& rlist, const color& c)
//{
//    list<rectangle2<int> >::const_iterator iter;
//
//    for (iter = rlist.begin(); iter != rlist.end(); ++iter) {
//        draw_box(*iter, c);
//    }
//}

void img::from_bits_grayscale(img &result, void *data, int dimx, int dimy, 
    unsigned long *masks, img::bit_type type)
{
    result.resize(dimx, dimy);
    result.grayscale = true;

    switch (type) {
        case BIT_TYPE_32 : 
            {
                int redMask = *masks;
                int greenMask = *(masks + 1);
                int blueMask = *(masks + 2);
                int* pixels = (int*)data;

                HOP_REAL* p = result.ptr(0, 0);
                for (int j = dimy; j > 0; --j) {
                    for (int i = dimx; i > 0; --i) {
                        int color = *pixels;
                        int red = (color & redMask) >> 16;
                        int green = (color & greenMask) >> 8;
                        int blue = (color & blueMask);
                        //fprintf(file, "%u, %u, %u", red, green, blue);
                        *p = (0.299*red + 0.587*green + 0.114*blue)/255.0;
                        ++pixels; ++p;
                    }
                }
            }
            break;
        default : 
            throw new_libhop_exception("Unsupported bitmap type.");
            break;
    };

}

void img::to_grayscale()
{
	if (!grayscale) {
		img_for(*this, iter) { 
			IMG_COLOR& ic = (IMG_COLOR&)(*iter);
			double d = (0.299*ic.red + 0.587*ic.green + 0.114*ic.blue)/255.0;
			*iter = d;
		}
		grayscale = true;
	}
}

void img::get_colors(img& R, img& G, img& B) const
{
	R.resize(width, height);
	R.grayscale = true;
	G.resize(width, height);
	G.grayscale = true;
	B.resize(width, height);
	B.grayscale = true;

	for (int i = 0; i < size(); i++) {
		IMG_COLOR& ic = (IMG_COLOR&)((*this)[i]);
		R[i] = (double)ic.red/255;
		G[i] = (double)ic.green/255;
		B[i] = (double)ic.blue/255;
	}
}

void img::from_bits_color(img &result, void *data, int dimx, int dimy, 
    unsigned long *masks, img::bit_type type)
{
    result.resize(dimx, dimy);
    result.grayscale = true;

    switch (type) {
        case BIT_TYPE_32 : 
            {
                int redMask = *masks;
                int greenMask = *(masks + 1);
                int blueMask = *(masks + 2);
                int* pixels = (int*)data;

                HOP_REAL* p = result.ptr(0, 0);
                for (int j = dimy; j > 0; --j) {
                    for (int i = dimx; i > 0; --i) {
                        int color = *pixels;
                        int red = (color & redMask) >> 16;
                        int green = (color & greenMask) >> 8;
                        int blue = (color & blueMask);
                        //fprintf(file, "%u, %u, %u", red, green, blue);
                        *p = COL_TO_REAL(red, green, blue);
                        ++pixels; ++p;
                    }
                }
            }
            break;
        default : 
            throw new_libhop_exception("Unsupported bitmap type.");
            break;
    };

}

void img::get_bits(void*& data, img::bit_type type)
{
    if (!grayscale)
        return;

    switch (type) {
        case BIT_TYPE_32 : 
            {
                data = malloc(width*height*32/8);

                HOP_REAL* p = ptr(0, 0);
                DWORD* dp = (DWORD*)data;
                HOP_REAL maxval, minval;
                HOP_REAL k, n;
                    
                minmax(minval, maxval);
                k = 255.0/(maxval - minval);
                n = -minval*k;
                for (int i = (int)width*(int)height; i > 0; --i) {
                    DWORD g = (DWORD)((*p)*k + n);
                    //DWORD res = (g | (g << 8) | (g << 16));
                    *dp = (g | (g << 8) | (g << 16));
                    ++p; ++dp;
                }
            }
            break;
        default: 
            throw new_libhop_exception("Unsupported bitmap type.");
            break;
    }
}

void img::get_bits(void*& data, int w, int h, img::bit_type type)
{
    if (!grayscale)
        return;

    switch (type) {
        case BIT_TYPE_32 : 
            {
                int size = w*h*32/8;				

				if (size <= 0) {
					data = nullptr;
					break;
				}
                data = malloc(size);
                memset(data, 0, size);

                HOP_REAL* p = ptr(0, 0);
                DWORD* dp = (DWORD*)data;
                HOP_REAL maxval, minval;
                HOP_REAL k, n;
                //int maxw = std::min(w, (int)width);
				int maxw = min(w, (int)width);
                //int maxh = std::min(h, (int)height);
				int maxh = min(h, (int)height);
                int wstepsrc = (int)width - maxw;
                int wstepdest = w - maxw;
                    
                minmax(minval, maxval);
                k = 255.0/(maxval - minval);
                n = -minval*k;
                for (int j = maxh; j > 0; --j) {
                    for (int i = maxw; i > 0; --i) {
                        DWORD g = (DWORD)((*p)*k + n);
                        
                        *dp = (g | (g << 8) | (g << 16));
                        ++p; ++dp;
                    }
                    p += wstepsrc;
                    dp += wstepdest;
                }
            }
            break;
        default: 
            throw new_libhop_exception("Unsupported bitmap type.");
            break;
    }
}

void point_img(img& im, const vector<ipoint2>& pts, int border, bool bigpoint)
{
    irectangle2 box = irectangle2::bounding_rectangle(pts.begin(), pts.end());

    im.resize(box.x_dim() + 2*border + 1, box.y_dim() + 2*border + 1, 0.0);
    for (auto piter = pts.begin(); piter != pts.end(); ++piter) {
        if (bigpoint) im.draw_big_point(piter->x + border - box.ll.x, piter->y + border - box.ll.y, 1.0);
        else im(piter->x + border - box.ll.x, piter->y + border - box.ll.y) = 1.0;
    }
}
