//
// Wrapper for image functions
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _IMG_H_
#define _IMG_H_

//#ifdef IMG_EXPORTS
//#define IMG_API __declspec(dllexport)
//#else
//#define IMG_API __declspec(dllimport)
//#endif

#define _USE_MATH_DEFINES  // necessary for <math.h> to define M_PI,.... constants, etc. 

#include <math.h>
#pragma warning (disable: 4996)
#include "matrix.h"
#include "structures.h"
#include "streaming.h"

#include "../opencl/cl_utils.h"
#include "../opencl/clContextManager.h"

#define HOP_REAL double
#define COL_TYPE unsigned char
#define COL_TYPE_MAX 255

#define COL_FROM_REAL(_color) (COL_TYPE)(255*(_color))

using namespace std;

typedef matrix<int> iimg;

// 8 byte double -> color
struct IMG_COLOR {
    unsigned int dummy1;    // 4 bytes
    COL_TYPE red;           // 1 byte
    COL_TYPE green;         // 1 byte
    COL_TYPE blue;          // 1 byte
    unsigned char dummy2;   // 1 byte
};

struct IMG_COLOR2 { 
    COL_TYPE red;
    COL_TYPE green;
    COL_TYPE blue;
    unsigned char dummy2;
};

struct IMG_COLOR3 {
    unsigned int dummy1;
    IMG_COLOR2 color;
};

struct IMG_COLOR4 { 
    COL_TYPE red;
    COL_TYPE green;
    COL_TYPE blue;
};


// public functions
///////////////////////////////////////////////////////////////////////////////

inline HOP_REAL COL_TO_REAL(int red, int green, int blue)
{
    IMG_COLOR result;
    HOP_REAL* presult = (HOP_REAL*)(&result);
    
    result.red = red;
    result.green = green;
    result.blue = blue;
    return *presult;
}

// color
//////////

struct color {
protected:
    int r, g, b;
public:
    color(int vr = 0, int vg = 0, int vb = 0) : r(vr), g(vg), b(vb) { }
    color(const color& c) : r(c.r), g(c.g), b(c.b) { }

    int red() const { return r; }
    double red01() const { return r/255.0; }
    int green() const { return g; }
    double green01() const { return g/255.0; }
    int blue() const { return b; }
    double blue01() const { return b/255.0; }
    int& red() { return r; }
    int& green() { return g; }
    int& blue() { return b; }

    void set(int i, int val) { if (i == 0) r = val; else if (i == 1) g = val; else if (i == 2) b = val; }

    void to_rgb(int& R, int& G, int& B) const { R = r; G = g; B = b; }

    static void to_color_vector(vector<color>& result, const vector<int>& cvec, int def = 0);
    static void to_color(color& result, const vector<int>& cvec, int def = 0);
    static color from_hsv(double h, double s, double v);
};



// img
///////////////////////////////////////////////////////////////////////////////

#define img_for(_img, _iter) for (rmatrix::iterator _iter = (_img).begin(); _iter != (_img).end(); ++_iter)
#define img_for_const(_img, _iter) for (rmatrix::const_iterator _iter = (_img).begin(); _iter != (_img).end(); ++_iter)


class img : public rmatrix {
public:
    enum bit_type { BIT_TYPE_32 };

    bool grayscale;

    img(const string& name);  // same as read_grayscale
    img(bool gscale = true) : rmatrix(), grayscale(gscale) { }
    img(size_type x, size_type y, bool gscale = true) : rmatrix(x, y), grayscale(gscale) { }
    img(int x, int y, bool gscale = true) : rmatrix(x, y), grayscale(gscale) { }
    img(int x, int y, HOP_REAL f, bool gscale = true) : rmatrix(x, y, f), grayscale(gscale) {  }
    img(size_type x, size_type y, HOP_REAL f, bool gscale = true) : rmatrix(x, y, f), grayscale(gscale) { }
    img(const matrix<double>& m, bool gscale = true, bool nneg = true) 
        : rmatrix(m), grayscale(gscale) { if (nneg) set_max(0.0); }
    img(const img& im) : rmatrix(im), grayscale(im.grayscale) { }
    //img(const bimg& bim);
    img(const map<ipoint2, HOP_REAL>& m, const HOP_REAL& val = 0.0, bool gscale = true, bool nneg = true) 
        : rmatrix(m, val), grayscale(gscale) { if (nneg) set_max(0.0); }
    img(const set<ipoint2>& m, const HOP_REAL& val0 = 0.0, const HOP_REAL& val1 = 1.0, bool gscale = true) 
        : rmatrix(m, val0, val1), grayscale(gscale) { }
    
	//CUDA code
	img(matrix<float>*);
	//CUDE code ends
    img(const img& im, bool gscale, bool normalized);

	void to_grayscale();
	void get_colors(img& R, img& G, img& B) const;

    static void read_grayscale(const string& name, img& result);
	static void read_colors(const string& name, img& result, bool as_hsv = false);

	static void read_colors(const string& name, img& red, img& green, img& blue, bool as_hsv = false);

    static void from_bits_grayscale(img& result, void* data, int dimx, int dimy, 
        unsigned long* masks, bit_type type);
    static void from_bits_color(img& result, void* data, int dimx, int dimy, 
        unsigned long* masks, bit_type type);
    void get_bits(void*& data, bit_type type);
    void get_bits(void*& data, int w, int h, bit_type type);

    void set_color(int x, int y, int red, int green, int blue)
    { 
        IMG_COLOR& ic = (IMG_COLOR&)at(x, y);
        ic.red = (COL_TYPE)red;
        ic.green = (COL_TYPE)green;
        ic.blue = (COL_TYPE)blue;
    }

    void set_color(int x, int y, double red, double green, double blue)
    {
        set_color(x, y, COL_FROM_REAL(red), COL_FROM_REAL(green), COL_FROM_REAL(blue));
    }

    void replace_color(HOP_REAL what, HOP_REAL with);
    img to_color_img(HOP_REAL color) const;
    img to_color_img() const;

	// opencl version of convolution	
	img* get_convolve_ocl(img& mask);
#ifdef OPENCL
    img* get_convolve(img& mask, bool use_opencl = true); // if using opencl then by default set use_opencl to true
#else
	img* get_convolve(img& mask, bool use_opencl = false);
#endif
	img* get_convolve(img& kernel, img& mask);
    void blur(int mask_size, double sigma, bool normalize = true);
    HOP_REAL p_get_2x2sum(img*& result);
    //img* get_grayscale();
    //iimg* get_01_img(HOP_REAL thresh);
    matrix<bool>* get_bool_matrix(HOP_REAL thresh);
    //void blur_resize(HOP_REAL factor, int startsize, int slimito, int msize, double sigma, int interp,
    //    std::vector<img*>& result);

    img& sqr() { img_for(*this, ptr) { const HOP_REAL& val = *ptr; *ptr = val*val; } return *this; }
	img& sqr(img &mask) {
		rmatrix::iterator ptr = begin();
		rmatrix::iterator ptr_mask = mask.begin();
		for (; ptr != end(); ++ptr, ++ptr_mask) {
			if (*ptr_mask != 0) {
				const HOP_REAL& val = *ptr; *ptr = val*val; 
				*ptr = val*val;
			}
		} 
		return *this; 
	}
    img& sqrt() { img_for(*this, ptr) { *ptr = ::sqrt(*ptr); } return *this; }
	img& absolute_value() { img_for(*this, ptr) { *ptr = ::fabs(*ptr); } return *this; }
	img& sqrt(img &mask) { 
		rmatrix::iterator ptr = begin();
		rmatrix::iterator ptr_mask = mask.begin();
		for (; ptr != end(); ++ptr, ++ptr_mask) {
			if (*ptr_mask != 0)
				*ptr = ::sqrt(*ptr); 
		} 
		return *this; 
	}
    img& positive();
    img& negative();
    img& negative2();
    img& neg();
	HOP_REAL norm();
	img& normalize(HOP_REAL min, HOP_REAL max);
    img& combine_max(const img&, int, int, HOP_REAL);
    img& combine_sum(const img&, int, int, HOP_REAL);
    //img& draw_frame(int x0, int y0, int x1, int y1, HOP_REAL color);
    img get_resized(int neww, int newh, bool normalize = true) const;
    void get_resized(img& result, int neww, int newh, bool normalize = true) const;
    void get_resized(img& result, int newsize, bool normalize = true) const;
    void crop(int w, int h);
    img cut(const irectangle2& rect) const;
    void flip_horizontal();

    img& blt_central_max(const img&, int, int, int, int);
    ipoint2 blt_central_max_ar(const img&, int, int, int, int);
    //img& blt_central_max(const img&, int, int, double, int, int);
    //img& blt_max(const img& src, int dx, int dy);
    static img concat(const std::vector<img*>&, double zero = 0.0);
    static img concat(std::vector<img>&, double zero = 0.0);
    static img concat_linear(const std::vector<img*>&);
    static img concat_linear(std::vector<img>&);
    static img concat_linear_fw(const std::vector<img*>&);
    static img concat_linear_fw(std::vector<img>&);
    static img string_to_img(const std::string& str);
    static img number_to_img(int i);

    //void print_matrix();
    void save(const string& name, int newxsize = -100) const;
    unsigned int save(void*& mem, int newxsize = -100) const;
    void save_normalized(const string& name, int newxsize = -100) const { save(name, newxsize); }
    void save_jet_colormap(const char*, int = -100);
    void save_visview(const char* name);

    //img& operator+=(const img& im) 
    //{
    //    HOP_REAL* ptrend = data + std::min(size(), im.size());
    //    HOP_REAL* ptrs = im.data;
    //    for (HOP_REAL *ptrd = data; ptrd < ptrend; ++ptrd) *ptrd += *(ptrs++);
    //    return *this;
    //}
    //img& operator-=(const img& im) 
    //{
    //    HOP_REAL* ptrend = data + std::min(size(), im.size());
    //    HOP_REAL* ptrs = im.data;
    //    for (HOP_REAL *ptrd = data; ptrd < ptrend; ++ptrd) *ptrd -= *(ptrs++);
    //    return *this;
    //}
    //img& operator-=(HOP_REAL d) 
    //{
    //    cimg_for(*this, ptr, HOP_REAL) { 
    //        *ptr -= d;
    //    }
    //    return *this;
    //}

    void write_to_stream(ostreamer& os) const;
    //static img read_from_stream(istreamer& is);
    void read_from_stream(istreamer& is);
    

    static void gaussian_mask(img& result, int dimx, int dimy, double sigma);
    static img gaussian_mask2(int dimx, int dimy, double sigma);
    static img* gaussian_mask(int dimx, int dimy, double sigma);
    static img* gabor_mask(double lambda, double theta, double phi, double gamma, double bw);
	static img* log_gabor_mask(int scale, int orientation);
    //double convolve_central(const img&, int , int, int, int);

    void draw_line(int x1, int y1, int x2, int y2, HOP_REAL color);
    void draw_box(const rectangle2<int>& rect, HOP_REAL color);
    void draw_box(const rectangle2<int>& rect, HOP_REAL color, int w);
    void draw_box(const rectangle2<int>& rect, const color& c);
    void draw_rectangle(const rectangle2<int>& rect, HOP_REAL c);
    void draw_big_point(int x, int y, HOP_REAL* col);
    void draw_big_point(int x, int y, HOP_REAL col);
    void revert_white(double white = 1.0);
    void thicken_white(int r, double white = 1.0);
    //void draw_boxes(const list<rectangle2<int> >& rlist, const color& c);

protected:
    unsigned int save_grayscale(void*& mem, int newxsize) const;
    void save_grayscale(const string& name, int newxsize) const;
    void save_colored(const string& name, int newxsize) const;

};


void point_img(img& im, const vector<ipoint2>& pts, int border, bool bigpoint);

template <class I> void point_img(img& im, I begin, I end, int border, bool bigpoint)
{
    point_img(im, vector<ipoint2>(begin, end), border, bigpoint);
}

#endif /* _IMG_H_ */


//// This class is exported from the img.dll
//class IMG_API Cimg {
//public:
//	Cimg(void);
//	// TODO: add your methods here.
//};
//
//extern IMG_API int nimg;
//
//IMG_API int fnimg(void);
