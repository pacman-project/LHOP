/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// matrix.h -- simple matrix based on vector<T>

#pragma once
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <functional>
#include <stdexcept>
#if defined WIN32 | defined WIN64
#include <windows.h>
#endif
#include "utils/structures.h"

#include "opencv2/opencv.hpp"

using namespace std;

#define for_each_xy(_matrix, _i, _j) for (size_t _i = 0; _i < (_matrix).width; ++_i)\
    for (size_t _j = 0; _j < (_matrix).height; ++_j)

#define for_each_xy_int(_matrix, _i, _j) for (int _i = 0; _i < (int)(_matrix).width; ++_i)\
    for (int _j = 0; _j < (int)(_matrix).height; ++_j)

#define for_each_element(_matrix, _i) for (size_t _i = 0, _max = (_matrix).size(); _i < _max; ++i)
#define for_each_iter(_matrix, _itertype, _iter) for (_itertype _iter = (_matrix).begin(); _iter != (_matrix).end(); ++_iter)
#define for_each_iter2(_matrix, _iter) for (_matrix::iter_t _iter = (_matrix).begin(); _iter != (_matrix).end(); ++_iter)

template<class T> class matrix : public std::vector<T> {
public:
    typedef typename vector<T>::iterator iter_t;


    typename vector<T>::size_type width, height;

    matrix() : std::vector<T>(), width(0), height(0) { }
    matrix(typename vector<T>::size_type w, typename vector<T>::size_type h) 
        : std::vector<T>(w*h), width(w), height(h) { }
    matrix(int w, int h) 
        : std::vector<T>(w*h), width((typename vector<T>::size_type)w), height((typename vector<T>::size_type)h) { }
    matrix(typename vector<T>::size_type w, typename vector<T>::size_type h, const T& val) 
        : std::vector<T>(w*h, val), width(w), height(h) { }
    matrix(int w, int h, const T& val) 
        : std::vector<T>(w*h, val), width((typename vector<T>::size_type)w), height((typename vector<T>::size_type)h) { }
    matrix(const matrix<T>& m) 
        : std::vector<T>(m), width(m.width), height(m.height) { }
    matrix(const std::map<ipoint2, T>& m, const T& val)
        : std::vector<T>(), width(0), height(0)
    {
        if (!m.empty()) {
            irectangle2 box;
            for (typename std::map<ipoint2, T>::const_iterator iter = m.begin(); iter != m.end(); ++iter) 
                box.eat(iter->first);
            resize(box.x_dim() + 1, box.y_dim() + 1, val);
            for (typename std::map<ipoint2, T>::const_iterator iter = m.begin(); iter != m.end(); ++iter) {
                const ipoint2& p = iter->first;
                std::vector<T>::at((p.x - box.ll.x) + width*(p.y - box.ll.y)) = iter->second;
            }
        }
    }
    matrix(const std::set<ipoint2>& m, const T& val0, const T& val1)
        : std::vector<T>(), width(0), height(0)
    {
        if (!m.empty()) {
            irectangle2 box;
            for (std::set<ipoint2>::const_iterator iter = m.begin(); iter != m.end(); ++iter) 
                box.eat(*iter);
            resize(box.x_dim() + 1, box.y_dim() + 1, val0);
            for (std::set<ipoint2>::const_iterator iter = m.begin(); iter != m.end(); ++iter) {
                const ipoint2& p = *iter;
                std::vector<T>::at((p.x - box.ll.x) + width*(p.y - box.ll.y)) = val1;
            }
        }
    }
	
    void resize(int w, int h)
    {
        resize(
            (w < 0) ? (typename vector<T>::size_type)(-w*width/100) : (typename vector<T>::size_type)w,
            (h < 0) ? (typename vector<T>::size_type)(-h*height/100) : (typename vector<T>::size_type)h
        );
    }

    void resize(typename vector<T>::size_type w, typename vector<T>::size_type h) { width = w; height = h; std::vector<T>::resize(w*h); }
    void resize(typename vector<T>::size_type w, typename vector<T>::size_type h, const T& init) { width = w; height = h; std::vector<T>::resize(w*h, init); }
    void dispose() 
    { 
        vector<T> emptyvec;

        width = height = 0; 
        this->swap(emptyvec);
    }

    void fill(const T& init) { std::fill(this->begin(), this->end(), init); }

    typename vector<T>::reference at(int i, int j) { return std::vector<T>::at(i + width*j); }
    T at(int i, int j, const T& defval) const { return (i < 0 || j < 0 || i >= width || j >= height) ? defval : std::vector<T>::at(i + width*j); }
    typename vector<T>::reference at(typename vector<T>::size_type i, typename vector<T>::size_type j) { return std::vector<T>::at(i + width*j); }
    typename vector<T>::const_reference at(int i, int j) const { return std::vector<T>::at(i + width*j); }
    typename vector<T>::const_reference at(typename vector<T>::size_type i, typename vector<T>::size_type j) const { return std::vector<T>::at(i + this->width*j); }
    typename vector<T>::iterator iter_at(int i, int j) { return this->begin() + i + this->width*j; }
    T* ptr(int i, int j) { return &at(i, j); }
    const T* ptr(int i, int j) const { return &at(i, j); }
    void set_at_safe(int i, int j, const T& val) 
    { 
        if (i >= 0 && j >= 0 && i < width && j < height) std::vector<T>::at(i + width*j) = val;
    }

    typename vector<T>::reference operator()(int i, int j) { return at(i, j); }
    typename vector<T>::reference operator()(typename vector<T>::size_type i, typename vector<T>::size_type j) { return at(i, j); }
    typename vector<T>::const_reference operator()(typename vector<T>::size_type i, typename vector<T>::size_type j) const { return at(i, j); }
    typename vector<T>::reference operator[](int i) { return std::vector<T>::at(i); }
    typename vector<T>::const_reference operator[](int i) const { return std::vector<T>::at(i); }
    typename vector<T>::reference operator[](typename vector<T>::size_type i) { return std::vector<T>::at(i); }
    typename vector<T>::const_reference operator[](typename vector<T>::size_type i) const { return std::vector<T>::at(i); }

    matrix<T>& operator=(const matrix<T>& m)
    {
        vector<T>::operator=(m);
        width = m.width;
        height = m.height;
        return *this;
    }

    void set_max(const T& val) 
    {
        for_each_iter(*this, typename vector<T>::iterator, iter) 
            if (*iter < val) *iter = val;
    }

    void set_region(int i0, int i1, int j0, int j1, const T& val) 
    {
        int imin = max(i0, 0), imax = min(i1, (int)width);
        int jmin = max(j0, 0), jmax = min(j1, (int)height);

        for (int i = imin; i < imax; ++i) {
            for (int j = jmin; j < jmax; ++j) at(i, j) = val;
        }
    }

    void set_region_circ(int cx, int cy, int radius, const T& val)
    {
        int imin = max(cx - radius, 0), imax = min(cx + radius + 1, (int)width);
        int jmin = max(cy - radius, 0), jmax = min(cy + radius + 1, (int)height);
        int radius2 = radius*radius;

        for (int i = imin; i < imax; ++i) {
            int i2 = (i - cx)*(i - cx);

            for (int j = jmin; j < jmax; ++j) 
                if (i2 + (j - cy)*(j - cy) <= radius2) at(i, j) = val;
        }
    }

    void set_region_c(int x, int y, int dx, int dy, const T& val)
    {
        set_region(x - dx, x + dx + 1, y - dy, y + dy + 1, val);
    }

    int count_region(int i0, int i1, int j0, int j1, const T& val) 
    {
        int imin = max(i0, 0), imax = min(i1, (int)width);
        int jmin = max(j0, 0), jmax = min(j1, (int)height);
		int result = 0;

        for (int i = imin; i < imax; ++i) {
            for (int j = jmin; j < jmax; ++j) 
				if (at(i, j) == val) ++result;
        }
		return result;
    }

    T maximum() const
    {
        if (this->size() == 0) return T();
        typename vector<T>::const_iterator result = this->begin();
        for (typename vector<T>::const_iterator i = this->begin() + 1; i != this->end(); ++i) {
            if (*result < *i) result = i;
        }
        return *result;
    }

    void minmax(T& min, T& max)
    {
        if (this->empty()) {
            min = max = T();
            return;
        }

        typename vector<T>::iterator mini, maxi;
        
        mini = maxi = this->begin();
        for (typename vector<T>::iterator i = this->begin() + 1; i != this->end(); ++i) {
            if (*maxi < *i) maxi = i; 
            else if (*i < *mini) mini = i;
        }
        min = *mini; max = *maxi;
    }

    T maximum(int i0, int i1, int j0, int j1)
    {
        if (this->size() == 0) return T();

        int imin = max(i0, 0), imax = min(i1, (int)width);
        int jmin = max(j0, 0), jmax = min(j1, (int)height);

        typename vector<T>::iterator result = iter_at(imin, jmin);
        typename vector<T>::iterator iter;
        for (int j = jmin; j < jmax; ++j) {
            iter = iter_at(imin, j);
            for (int i = imin; i < imax; ++i) {
                if (*result < *iter) result = iter;
                ++iter;
            }
        }
        return *result;
    }

    void change_values_leq(const T& cond, const T& newval) 
    { 
        for_each_iter(*this, typename vector<T>::iterator, iter) {
            T& v = *iter;
            if (v <= cond) v = newval;
        }
    }

    T maximum_c(int x, int y, int dx, int dy)
    {
        return maximum(x - dx, x + dx + 1, y - dy, y + dy + 1);
    }

    void add_border(int bx, int by, const T& d)
    {
        if (bx < 0 || by < 0) return;

        int tmpw = (int)width, tmph = (int)height;
        vector<T> tmp = *this;
        typename vector<T>::iterator titer, iter;
        int i, j;

        resize(width + 2*bx, height + 2*by); 
        
        int hmax = (int)height - by, wmax = (int)width - bx;

        for (j = 0; j < by; ++j)
            for (iter = this->begin() + j*this->width, i = 0; i < (int)this->width; ++i, ++iter) *iter = d;
        for (j = 0; j < tmph; ++j) {
            iter = this->begin() + (j + by)*width;
            for (i = 0; i < bx; ++i) *(iter++) = d;
            titer = tmp.begin() + j*tmpw;
            for (i = 0; i < tmpw; ++i) *(iter++) = *(titer++);
            for (i = 0; i < bx; ++i) *(iter++) = d;
        }
        for (j = 0; j < by; ++j)
            for (iter = this->begin() + (this->height - by + j)*width, i = 0; i < (int)width; ++i, ++iter) *iter = d;
    }

    void remove_border(int bx, int by) 
    {
        int tmpw = (int)width;
        vector<T> tmp = *this;

        resize(width - 2*bx, height - 2*by);

        typename vector<T>::iterator titer, iter;
        int i, j;

        iter = this->begin();
        titer = tmp.begin() + by*tmpw + bx;
        for (j = 0; j < (int)height; ++j) {
            for (i = 0; i < (int)width; ++i) *(iter++) = *(titer++);
            titer += 2*bx;
        }
    }

    matrix extract(int minx, int miny, int maxx, int maxy)
    {
        minx = max(0, minx);
        miny = max(0, miny);
        maxx = min(maxx, (int)width);
        maxy = min(maxy, (int)height);

        matrix result(maxx - minx, maxy - miny);

        for (int i = minx; i < maxx; ++i)
            for (int j = miny; j < maxy; ++j) {
                result(i - minx, j - miny) = at(i, j);
            }
        return result;
    }

    void blt_central(const matrix& src, int csx, int csy, int cdx, int cdy)
    {
        int minx = cdx - csx, miny = cdy - csy;
        int maxdx = minx + (int)src.width, maxdy = miny + (int)src.height;
        int mindx = max(minx, 0);
        int mindy = max(miny, 0);
        int minsx = mindx - minx, minsy = mindy - miny;
        maxdx = min(maxdx, (int)width);
        maxdy = min(maxdy, (int)height);

        for (int x = mindx, sx = minsx; x < maxdx; ++x, ++sx)
            for (int y = mindy, sy = minsy; y < maxdy; ++y, ++sy) 
                (*this)(x, y) = src(sx, sy);
    }

    void blt(const matrix& src, int minx, int miny)
    {
        int maxdx = minx + (int)src.width, maxdy = miny + (int)src.height;
        int mindx = max(minx, 0);
        int mindy = max(miny, 0);
        int minsx = mindx - minx, minsy = mindy - miny;
        maxdx = min(maxdx, (int)width);
        maxdy = min(maxdy, (int)height);

        for (int x = mindx, sx = minsx; x < maxdx; ++x, ++sx)
            for (int y = mindy, sy = minsy; y < maxdy; ++y, ++sy) 
                (*this)(x, y) = src(sx, sy);
    }

    template<class B> void blt_central(const matrix& src, int csx, int csy, int cdx, int cdy, 
        const B& f)
    {
        int minx = cdx - csx, miny = cdy - csy;
        int maxdx = minx + (int)src.width, maxdy = miny + (int)src.height;
        int mindx = max(minx, 0);
        int mindy = max(miny, 0);
        int minsx = mindx - minx, minsy = mindy - miny;
        maxdx = min(maxdx, (int)width);
        maxdy = min(maxdy, (int)height);

        for (int x = mindx, sx = minsx; x < maxdx; ++x, ++sx)
            for (int y = mindy, sy = minsy; y < maxdy; ++y, ++sy)  {
                (*this)(x, y) = f((*this)(x, y),src(sx, sy));
	    }
    }

    template<class B> void blt_central_ar(const matrix& src, int csx, int csy, int cdx, int cdy, 
        const B& f)
    {
        int mindx = cdx - csx, mindy = cdy - csy;
        int offx, offy;
        
        if (mindx >= 0) offx = 0; else { offx = -mindx; mindx = 0; } 
        if (mindy >= 0) offy = 0; else { offy = -mindy; mindy = 0; }
        
        int maxdx = mindx + (int)src.width, maxdy = mindy + (int)src.height;
        int plusx = max(0, maxdx - ((int)width + offx - 1));
        int plusy = max(0, maxdy - ((int)height + offy - 1));

        if (offx > 0 || offy > 0 || plusx > 0 || plusy > 0) {
            matrix tmp(*this);

            resize(width + offx + plusx, height + offy + plusy, 0);
            fill(this->begin(), this->end(), 0);
            blt(tmp, offx, offy);
        }

        for (int x = mindx, sx = 0; x < maxdx; ++x, ++sx)
            for (int y = mindy, sy = 0; y < maxdy; ++y, ++sy) 
                (*this)(x, y) = f((*this)(x, y), src(sx, sy));
    }

    // blt src matrix such that (csx, csy) -> center_of_matrix + (cdx, cdy)
    // The value is determined by binary function f.
    template<class B> void blt_central2_ar(const matrix& src, int csx, int csy, 
        int cdx, int cdy, const B& f)
    {
        cdx += (int)width/2;
        cdy += (int)height/2;
        int mindx = cdx - csx, mindy = cdy - csy;
        int maxdx = mindx + (int)src.width, maxdy = mindy + (int)src.height;
        int offx = 0, offy = 0;
        
        if (mindx < 0) offx = -mindx;
        if (mindy < 0) offy = -mindy;
        if (maxdx > (int)width) offx = max(offx, maxdx - (int)width);
        if (maxdy > (int)height) offy = max(offy, maxdy - (int)height);

        if (offx > 0 || offy > 0) {
            matrix tmp(*this);

            resize(width + 2*offx, height + 2*offy);
            fill(0);
            blt(tmp, offx, offy);
            mindx += offx;
            maxdx += offx;
            mindy += offy;
            maxdy += offy;
        }

        for (int x = mindx, sx = 0; x < maxdx; ++x, ++sx)
            for (int y = mindy, sy = 0; y < maxdy; ++y, ++sy) 
                (*this)(x, y) = f((*this)(x, y), src(sx, sy));
    }

    int count(const T& t)
    {
        int result = 0;

        for_each_element(*this, i) { 
            if (std::vector<T>::at(i) == t) ++result; 
        }
        return result;
    }

    int count_geq(const T& t)
    {
        int result = 0;

        for_each_element(*this, i) { 
            if (std::vector<T>::at(i) >= t) ++result; 
        }
        return result;
    }

    void print() 
    {
        for (typename vector<T>::size_type i = 0; i < width; ++i) {
            for (typename vector<T>::size_type j = 0; j < height; ++j) cout << at(i, j) << " ";
            cout << endl;
        }
    }
    matrix& operator*=(const T& t) 
    {
        for_each_iter(*this, typename vector<T>::iterator, iter) { *iter *= t; }
        return *this;
    }

    matrix& operator/=(const T& t) 
    {
        for_each_iter(*this, typename vector<T>::iterator, iter) { *iter /= t; }
        return *this;
    }

    matrix& operator+=(const T& t)
    {
        for_each_iter(*this, typename vector<T>::iterator, iter) { *iter += t; }
        return *this;
    }

    matrix& operator-=(const T& t)
    {
        for_each_iter(*this, typename vector<T>::iterator, iter) { *iter -= t; }
        return *this;
    }

    void shape_circle(int r, const T& inside, const T& outside) 
    {
        int dim = 2*r + 1;
        int r2 = r*r;

        resize(dim, dim, outside);
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                if ((i - r)*(j - r) <= r2) at(i, j) = inside;
            }
        }
    }

    T convolve_central(const matrix& src, int csx, int csy, int cdx, int cdy) const
    {
        int minx = cdx - csx, miny = cdy - csy;
        int maxdx = minx + (int)src.width, maxdy = miny + (int)src.height;
        int mindx = max(minx, 0);
        int mindy = max(miny, 0);
        int minsx = mindx - minx, minsy = mindy - miny;
        maxdx = min(maxdx, (int)width);
        maxdy = min(maxdy, (int)height);
        T result = 0;
        T d1, d2;
        T n1 = 0, n2 = 0;

        for (int x = mindx, sx = minsx; x < maxdx; ++x, ++sx)
            for (int y = mindy, sy = minsy; y < maxdy; ++y, ++sy) {
                d1 = (*this)(x, y); d2 = src(sx, sy);
                n1 += d1*d1; n2 += d2*d2;
                result += d1*d2;
            }

        if (n1 < 1.0E-7 || n2 < 1.0E-7) return 0;
        else return result/::sqrt(n1)/::sqrt(n2);
    }

    void to_point_map(std::map<ipoint2, T>& result, const T& zero) const
    {
        result.clear();
        ipoint2 center((int)width/2, (int)height/2);

        for_each_xy_int (*this, i, j) {
            if (at(i, j) != zero) result.insert(pair<ipoint2, T>(ipoint2(i, j) - center, at(i, j)));
        }
    }

    void to_point_set(std::set<ipoint2>& result, const T& zero) const
    {
        result.clear();
        ipoint2 center((int)width/2, (int)height/2);

        for_each_xy_int (*this, i, j) {
            if (at(i, j) != zero) result.insert(ipoint2(i, j) - center);
        }
    }

    void save_mathematica(const char* fname)
    {
        std::ofstream os(fname);
        print_mathematica(os);
        os.close();
    }

	void save_matlab(const char* fname)
    {
        std::ofstream os(fname);
        print_matlab(os);
        os.close();
    }

    void print_mathematica(std::ostream& os) const
    {
        os << '{';
        for (typename vector<T>::size_type i = 0; i < width; ++i) {
            if (i > 0) os << ',' << endl;
            os << '{'; 
            for (typename vector<T>::size_type j = 0; j < height; ++j) {
                if (j > 0) os << ',';
                os << (*this)(i, j);
            }
            os << '}';
        }
        os << '}';
    }

	void print_matlab(std::ostream& os)
	{
        for (typename vector<T>::size_type i = 0; i < width; ++i) {
            for (typename vector<T>::size_type j = 0; j < height; ++j) {
                if (j > 0) os << ',';
                os << (*this)(i, j);
            }
            os << endl;
        }
		os << endl;
    }

    void print_as_vector(std::ostream& os)
    {
        for (typename vector<T>::size_type j = 0; j < height; ++j) {
            if (j > 0) os << ',';
            for (typename vector<T>::size_type i = 0; i < width; ++i) {
                if (i > 0) os << ',';
                os << (*this)(i, j);
            }
        }
    }

	void normalize(T minimum,T maximum)
	{
		T minval;
		T maxval;
		minmax(minval,maxval);
		if(minval==maxval)return;
		T k;
		T n;
		k=(maximum-minimum)/(maxval-minval);
		n=minimum-minval*k;
		for(typename matrix<T>::iterator it=this->begin();it!=this->end();++it){
			*it=(*it)*k+n;
		}
	}

	cv::Mat toCvMat() {
		return std::logic_error("Specialization of method required: must be implemented for each different template type");
	}

private:	
	void copyToCvMat(cv::Mat& mat) {
		for (int i = 0; i < height; ++i)
			for (int j = 0; j < width; ++j)
				mat.at<T>(i,j) = at(j,i);
	}
};

typedef matrix<double> rmatrix;
typedef matrix<float> fmatrix;
typedef matrix<int> imatrix;

template<> inline
cv::Mat matrix<double>::toCvMat() {
	cv::Mat result(height, width, CV_64F);
	copyToCvMat(result);
	return result;
}

template<> inline
cv::Mat matrix<float>::toCvMat() {
	cv::Mat result(height, width, CV_32F);
	copyToCvMat(result);
	return result;
}

template<> inline
cv::Mat matrix<int>::toCvMat() {
	cv::Mat result(height, width, CV_32S);
	copyToCvMat(result);
	return result;
}

#endif  /* _MATRIX_H_ */


