// structures.h

#pragma once
#ifndef _STRUCTURES_H_
#define _STRUCTURES_H_

#include <vector>
#include <list>
#include <iostream>
#include <string>
#include <list>
#include <cmath>

using namespace std;

// "function" classes
///////////////////////////////////////////////////////////////////////////////



template<class P, class R> struct unary_function_t {
    typedef P argument_t;
    typedef R result_t;

    virtual R operator()(const P& x) const = 0;
};

typedef unary_function_t<double, double> unary_double_function;

// point2<T>
///////////////////////////////////////////////////////////////////////////////

template<class T> struct point2 {
    T x, y;

    point2() : x(0), y(0) { }
    point2(const T& xy) : x(xy), y(xy) { }
    point2(const T& px, const T& py) : x(px), y(py) { }
    point2(const point2& p) : x(p.x), y(p.y) { }
    template<class U> point2(const point2<U>& p) : x((T)p.x), y((T)p.y) { }

    static point2 zero;

    pair<T,T> to_pair() const { return pair<T,T>(x, y); }
    void swap(point2& p2) 
    { 
        T tmp;
        tmp = x; x = p2.x; p2.x = tmp;
        tmp = y; y = p2.y; p2.y = tmp;
    }

    point2& operator=(const point2& p) { x = p.x; y = p.y; return *this; }
    bool operator==(const point2& p) const { return x == p.x && y == p.y; }
    bool operator!=(const point2& p) const { return x != p.x || y != p.y; }
    bool operator<(const point2& p) const { return (x == p.x) ? y < p.y : x < p.x; }
    bool operator>(const point2& p) const { return (x == p.x) ? y > p.y : y > p.y; }
    point2 operator+(const point2& p) const { return point2(x + p.x, y + p.y); }
    point2 operator-(const point2& p) const { return point2(x - p.x, y - p.y); }
    point2 operator-() const { return point2(-x, -y); }
    point2& operator+=(const point2& p) { x += p.x; y += p.y; return *this; }
    point2& operator-=(const point2& p) { x -= p.x; y -= p.y; return *this; }
    point2& operator*=(double f) { x = (T)(x*f); y = (T)(y*f); return *this; }
	point2& operator*=(const point2& p) { x *= p.x; y *= p.y; return *this; }
    point2 operator/(const T& d) { return point2(x/d, y/d); }
    point2& operator/=(const T& d) { x /= d; y /= d; return *this; }
    point2& operator/=(double f) { x = (T)(x/f); y = (T)(y/f); return *this; }
    point2 operator*(double f) const { return point2((T)(x*f), (T)(y*f)); }
    template<class S> operator point2<S>() const { return point2<S>((S)x, (S)y); } 
    pair<T, T> to_pair() { return pair<T, T>(x, y); }
    bool is_zero() const { return x == 0 && y == 0; }
    T maximum() const { return max<T>(x, y); }
    T minimum() const { return min<T>(x, y); }
    point2 abs() const { return point2( std::abs(x), std::abs(y)); }

    point2 get_rotated(const point2& center, double angle) const
    {
        point2 q = *this - center;

        return center + point2((T)(q.x * cos(angle) + q.y * sin(angle)), 
            (T)(q.y * cos(angle) - q.x * sin(angle)));
    }

    point2& add(const T& nx, const T& ny) { x += nx; y += ny; return *this; }
    point2& set(const T& nx, const T& ny) { x = nx; y = ny; return *this; }
    point2& div(const T& d) { x /= d; y /= d; return *this; }

    T inner_product(const point2& p) const { return x*p.x + y*p.y; }
    T norm2() const { return x*x + y*y; }
    T distance2(const point2& p) const { return (*this - p).norm2(); }
    T distance2(const T& px, const T& py) const { return norm2(x - px, y - py); }
    point2 rotated90() { return point2(-y, x); }

    static T norm2(const T& x, const T& y) { return x*x + y*y; }
    static T distance2(const point2& p, const point2& q) { return (p - q).norm2(); }
    static bool x_less(const point2& p, const point2& q) { return p.x < q.x; }
    static bool y_less(const point2& p, const point2& q) { return p.y < q.y; }

    template<class I> static point2 center(I begin, I end)
    {
        point2 result;
        int count = 0;

        while (begin != end) {
            result += *begin;
            ++count;
            ++begin;
        }
        return result/(T)count;
    }

    template<class I> static point2 center_p(I begin, I end)
    {
        point2 result;
        int count = 0;

        while (begin != end) {
            result += **begin;
            ++count;
            ++begin;
        }
        return result/(T)count;
    }

    template<class I> pair<T, I> closest_point(I begin, I end) const
    {
        if (begin == end) return pair<T,I>(0, end);

        I closest = begin;
        T mindist = distance2(*closest);
        T d;

        while (++begin != end) {
            if ((d = distance2(*begin)) < mindist) { 
                mindist = d; 
                closest = begin; 
            }
        }
        return pair<T, I>(mindist, closest);
    }

    template<class I> pair<T, I> closest_point_p(I begin, I end) const
    {
        if (begin == end) return pair<T,I>(0, end);

        I closest = begin;
        T mindist = distance2(**closest);
        T d;

        while (++begin != end) {
            if ((d = distance2(**begin)) < mindist) { 
                mindist = d; 
                closest = begin; 
            }
        }
        return pair<T, I>(mindist, closest);
    }

    template<class I> I closest_point2(I begin, I end) const
    {
        if (begin == end) return end;

        I closest = begin;
        T mindist = distance2(closest->second);
        T d;

        while (++begin != end) {
            if ((d = distance2(begin->second)) < mindist) { 
                mindist = d; 
                closest = begin; 
            }
        }
        return closest;
    }

    template<class I> static void centralize_list(const point2& p, I begin, I end)
    {
        while (begin != end) {
            *begin -= p;
            ++begin;
        }
    }

    template<class I> static void expand_list(I begin, I end, double factor)
    {
        while (begin != end) { *begin *= factor; ++begin; }
    }

    template<class I> static vector< point2<T> > split_list(I begin, I end) 
    {
        vector<point2<T> > result;

        while (begin != end) {
            T x = *begin;
            if (++begin != end) { result.push_back(point2<T>(x, *begin)); ++begin; }
        }
        return result;
    }

    //void read_from_stream(istreamer& is) 
    //{
    //    is.read(x);
    //    is.read(y);
    //}

    //void write_to_stream(ostreamer& os) const
    //{
    //    os.write(x);
    //    os.write(y);
    //}

    void mma_write(ostream& os)
    {
        os << '{' << x << ',' << y << '}';
    }

    friend ostream& operator<<(ostream& os, const point2& pt)
    {
        os << pt.x << ' ' << pt.y;
        return os;
    }

};

typedef point2<int> ipoint2;
typedef vector<ipoint2> ip2_vector;
typedef vector<double> dvector;
typedef point2<double> dpoint2;
typedef point2<float> fpoint2;

//template<> point2<int> point2<int>::zero;
//template<> point2<double> point2<double>::zero;

inline double norm(const dpoint2& p)
{
    return sqrt(p.norm2());
}

inline void normalize(dpoint2& p)
{
    double n = sqrt(p.norm2());
    if (n != 0) { p.x /= n; p.y /= n; }
}

inline dpoint2 rotate_point(dpoint2& p, double angle)
{
    double cangle = cos(angle), sangle = sin(angle);

    return dpoint2(p.x*cangle - p.y*sangle, p.y*cangle + p.x*sangle);
}

template<class T> inline T dot_product(const point2<T>& a, const point2<T>& b)
{
    return a.x*b.x + a.y*b.y;
}

template<class T> inline T planar_cross_product(const point2<T>& a, const point2<T>& b)
{
    return a.x*b.y - a.y*b.x;
}


// line2<T>
///////////////////////////////////////////////////////////////////////////////

template<class T> struct line2 {
    point2<T> a, b; 
};

typedef line2<int> iline2;


// rectangle2<T>
///////////////////////////////////////////////////////////////////////////////

template<class T> struct rectangle2 {
    point2<T> ll, ur;

    rectangle2(const T& lx = 1, const T& ly = 1, const T& ux = 0, const T& uy = 0) :
        ll(lx, ly), ur(ux, uy) { }
    rectangle2(const point2<T>& pll, const point2<T>& pur) : ll(pll), ur(pur) { }
    rectangle2(const rectangle2& r) : ll(r.ll), ur(r.ur) { }

    T x_dim() const { return ur.x - ll.x; }
    T y_dim() const { return ur.y - ll.y; }
    T size2() const { return ll.distance2(ur); }

    bool invalid() const { return ur.x < ll.x || ur.y < ll.y; }
    void invalidate() { ll.x = ll.y = 1; ur.x = ur.y = 0; }

    void eat(const T& x, const T& y)
    {
        if (invalid()) { ll.x = ur.x = x;  ll.y = ur.y = y; }
        else {
            if (x < ll.x) ll.x = x; else if (x > ur.x) ur.x = x;
            if (y < ll.y) ll.y = y; else if (y > ur.y) ur.y = y;
        }
    }

    void eat(const point2<T>& p) { eat(p.x, p.y); }

    rectangle2 intersection(const rectangle2& r) const
    {
        return rectangle2(max(ll.x, r.ll.x), max(ll.y, r.ll.y), min(ur.x, r.ur.x), min(ur.y, r.ur.y));
    }

    rectangle2 bounding_rectangle(const rectangle2& r) const
    {
        return rectangle2(min(ll.x, r.ll.x), min(ll.y, r.ll.y), max(ur.x, r.ur.x), max(ur.y, r.ur.y));
    }

    T union_area(const rectangle2& r) const
    {
        rectangle2 i = intersection(r);
        T result = area() + r.area();
        if (i.invalid()) return result; else return result - i.area();
    }

    rectangle2 get_rotated(const point2<T>& center, double angle) const
    {
        rectangle2 result;
        ipoint2 q;

        result.eat(ll.get_rotated(center, angle));
        result.eat(point2<T>(ll.x, ur.y).get_rotated(center, angle));
        result.eat(ur.get_rotated(center, angle));
        result.eat(point2<T>(ur.x, ll.y).get_rotated(center, angle));
        return result;
    }

    rectangle2 get_scaled(double factor) const 
    {
        return rectangle2((T)(ll.x*factor), (T)(ll.y*factor), (T)(ur.x*factor), (T)(ur.y*factor));
    }

    T area() const { return invalid() ? 0 : (ur.x - ll.x)*(ur.y - ll.y); }

    rectangle2& operator=(const rectangle2& r)
    {
        ll = r.ll; ur = r.ur; 
        return *this;
    }

    rectangle2 operator+(const point2<T>& p) 
    {
        return rectangle2(ll + p, ur + p);
    }

    rectangle2 operator-(const point2<T>& p) 
    {
        return rectangle2(ll - p, ur - p);
    }

    rectangle2 operator+(const T& t) 
    {
        return rectangle2(ll.x + t, ll.y + t, ur.x + t, ur.y + t);
    }

    template<class M> rectangle2 operator*(const M& t) 
    {
        return rectangle2((T)(ll.x * t), (T)(ll.y * t), (T)(ur.x * t), (T)(ur.y * t));
    }

    rectangle2& operator+=(const point2<T>& p)
    {
        ll += p; ur += p;
        return *this;
    }

    rectangle2& operator-=(const point2<T>& p)
    {
        ll -= p; ur -= p;
        return *this;
    }

    void resize(double factor) 
    {
        ll *= factor; ur *= factor;
    }

    rectangle2& grow(const T& g) 
    {
        ll.x -= g; ll.y -= g;
        ur.x += g; ur.y += g;
        return *this;
    }

    rectangle2& grow(const point2<T>& p) 
    {
        ll.x -= p.x; ll.y -= p.y;
        ur.x += p.x; ur.y += p.y;
        return *this;
    }

    rectangle2& grow(double f) 
    {
        point2<T> c = center();

        ll -= c; ur -= c;
        ll *= f; ur *= f;
        ll += c; ur += c;
        return *this;
    }


    bool inside(const point2<T>& p) const
    {
        return p.x >= ll.x && p.x <= ur.x && p.y >= ll.y && p.y <= ur.y;
    }

    bool inside(const T& x, const T& y) const
    {
        return x >= ll.x && x <= ur.x && y >= ll.y && y <= ur.y;
    }

    // Returns true if 'r' is inside 'this'.
    bool inside(const rectangle2& r) const
    {
        return inside(r.ll) && inside(r.ur);
    }

    bool subset_of(const rectangle2& r) const
    {
        return r.inside(*this);
    }

    bool inside_llopen(const T& x, const T& y) const
    {
        return x > ll.x && x <= ur.x && y > ll.y && y <= ur.y;
    }

    bool inside_uropen(const T& x, const T& y) const
    {
        return x >= ll.x && x < ur.x && y >= ll.y && y < ur.y;
    }

    point2<T> center() const { return (ur + ll)/2; }

    //void read_from_stream(istreamer& is) 
    //{
    //    ll.read_from_stream(is);
    //    ur.read_from_stream(is);
    //}

    //void write_to_stream(ostreamer& os) const
    //{
    //    ll.write_to_stream(os);
    //    ur.write_to_stream(os);
    //}

    template<class I> static rectangle2<T> bounding_rectangle(I begin, I end) 
    {
        if (begin == end) return rectangle2<T>();

        point2<T> ll = *begin;
        point2<T> ur = *begin;

        while (++begin != end) {
            const point2<T>& p = *begin;

            if (p.x < ll.x) ll.x = p.x; else if (p.x > ur.x) ur.x = p.x;
            if (p.y < ll.y) ll.y = p.y; else if (p.y > ur.y) ur.y = p.y;
        }
        return rectangle2(ll, ur);
    }

    static void inhibit(list<rectangle2<T> >& l, double thr) 
    {
        for (typename list<rectangle2<T> >::iterator i = l.begin(); i != l.end(); ++i) {
            
            typename list<rectangle2<T> >::iterator j = l.begin();

            while (j != l.end()) {
                if (j != i){
                    rectangle2<T> r = j->intersection(*i);
                    T jarea = j->area();

                    if (abs(r.area() - jarea) <= thr*jarea) {
                        j = l.erase(j);
                        continue;
                    }
                }
                ++j;
            }
        }
    }

    template<class S> static void inhibit(list<pair<rectangle2<T>, S> >& l, double thr) 
    {
        for (typename list<pair<rectangle2<T>, S> >::iterator i = l.begin(); i != l.end(); ++i) {
            
            typename list<pair<rectangle2<T>, S> >::iterator j = l.begin();

            while (j != l.end()) {
                if (j != i){
                    rectangle2<T> r = j->first.intersection(i->first);
                    T jarea = j->first.area();

                    if (abs(r.area() - jarea) <= thr*jarea) {
                        j = l.erase(j);
                        continue;
                    }
                }
                ++j;
            }
        }
    }

    //template<class S> static void inhibit(list<S>& l, double thr, unary_function<S, rectangle2<T> >& f_get_box) 
    template<class S, class F> static void inhibit(list<S> l, double thr, const F& f_get_box) 
    {
        for (typename list<S>::iterator i = l.begin(); i != l.end(); ++i) {
            
            typename list<S>::iterator j = l.begin();

            while (j != l.end()) {
                if (j != i){
                    rectangle2<T> r = f_get_box(*j).intersection(f_get_box(*i));
                    T jarea = f_get_box(*j).area();

                    if (abs(r.area() - jarea) <= thr*jarea) {
                        j = l.erase(j);
                        continue;
                    }
                }
                ++j;
            }
        }
    }

    void split(list<rectangle2>& result, int xn, int yn) const
    {
        result.clear();

        if (xn < 1 || yn < 1) return;

        T xstep = x_dim()/xn;
        T ystep = y_dim()/yn;
        
        if (xstep == 0) xstep = 1;
        if (ystep == 0) ystep = 1;
        for (int i = 0; i < xn - 1; ++i) {
            for (int j = 0; j < yn - 1; ++j) {
                result.push_back(rectangle2(ll.x + i*xstep, ll.y + j*ystep, 
                    ll.x + (i + 1)*xstep, ll.y + (j + 1)*ystep));
            }
            result.push_back(rectangle2(ll.x + i*xstep, ll.y + (yn - 1)*ystep, 
                ll.x + (i + 1)*xstep, ur.y));
        }
        for (int j = 0; j < yn - 1; ++j) {
            result.push_back(rectangle2(ll.x + (xn - 1)*xstep, ll.y + j*ystep, 
                ur.x, ll.y + (j + 1)*ystep));
        }
        result.push_back(rectangle2(ll.x + (xn - 1)*xstep, ll.y + (yn - 1)*ystep, 
            ur.x, ur.y));

    }

    void mma_write(ostream& os) 
    {
        os << '{';
        ll.mma_write(os);
        os << ',';
        ur.mma_write(os);
        os << '}';
    }

    friend ostream& operator<<(ostream& os, const rectangle2& rect)
    {
        os << rect.ll << ' ' << rect.ur;
        return os;
    }

    friend istream& operator>>(istream& is, rectangle2& rect)
    {
        is >> rect.ll.x >> rect.ll.y >> rect.ur.x >> rect.ur.y;
        return is;
    }

};

typedef rectangle2<int> irectangle2;
typedef rectangle2<double> drectangle2;
typedef rectangle2<float> frectangle2;

// wpoint2<T>
///////////////////////////////////////////////////////////////////////////////

template<class T, class W> struct wpoint2 : public point2<T> {
    W w;

    wpoint2() : point2<T>() { }
    wpoint2(const T& x, const T& y, const W& pw) : point2<T>(x, y), w(pw) { }
    wpoint2(const wpoint2& p) : point2<T>(p), w(p.w) { }

    wpoint2& operator=(const wpoint2& p) { this.x = p.x; this.y = p.y; w = p.w; return *this; }
    bool operator==(const wpoint2& p) const { return this.x == p.x && this.y == p.y && w == p.w; }

    //virtual void read_from_stream(istreamer& is) 
    //{
    //    point2<T>::read_from_stream(is);
    //    is.read(w);
    //}

    //virtual void write_to_stream(ostreamer& os) const
    //{
    //    point2<T>::write_to_stream(os);
    //    os.write(w);
    //}

};

typedef wpoint2<int, double> idpoint2;

// point3<T>
///////////////////////////////////////////////////////////////////////////////

template<class T> struct point3 {
    T x, y, z;

    point3() : x(0), y(0), z(0) { }
    point3(const T& px, const T& py, const T& pz) : x(px), y(py), z(pz) { }
    point3(const point3& p) : x(p.x), y(p.y), z(p.z) { }

    static point3 zero;

    void swap(point3& p2) 
    { 
        T tmp;
        tmp = x; x = p2.x; p2.x = tmp;
        tmp = y; y = p2.y; p2.y = tmp;
        tmp = z; z = p2.z; p2.z = tmp;
    }

    point3& operator=(const point3& p) { x = p.x; y = p.y; z = p.z; return *this; }
    bool operator==(const point3& p) const { return x == p.x && y == p.y && z == p.z; }
    point3 operator+(const point3& p) const { return point3(x + p.x, y + p.y, z + p.z); }
    point3 operator-(const point3& p) const { return point3(x - p.x, y - p.y, z - p.z); }
    point3 operator-() const { return point3(-x, -y, -z); }
    point3& operator+=(const point3& p) { x += p.x; y += p.y; z += p.z; return *this; }
    point3& operator-=(const point3& p) { x -= p.x; y -= p.y; z -= p.z; return *this; }
    point3& operator*=(double f) { x = (T)(x*f); y = (T)(y*f); z = (T)(z*f); return *this; }
    point3 operator/(const T& d) { return point3(x/d, y/d, z/d); }
    point3& operator/=(const T& d) { x/d; y/d; z/d; return *this; }

    T inner_product(const point3& p) const { return x*p.x + y*p.y + z*p.z; }
    T norm2() const { return x*x + y*y + z*z; }
    T distance2(const point3& p) const { return (*this - p).norm2(); }
    point3& set(const T& nx, const T& ny, const T& nz) { x = nx; y = ny; z = nz; return *this; }

    static T distance2(const point3& p, const point3& q) { return (p - q).norm2(); }

    template<class I> static point3 center(I begin, I end)
    {
        point3 result;
        int count = 0;

        while (begin != end) {
            result += *begin;
            ++count;
            ++begin;
        }
        return result/(T)count;
    }

    template<class I> static void centralize_list(const point3& p, I begin, I end)
    {
        while (begin != end) {
            *begin -= p;
            ++begin;
        }
    }

    //void read_from_stream(istreamer& is) 
    //{
    //    is.read(x);
    //    is.read(y);
    //    is.read(z);
    //}

    //void write_to_stream(ostreamer& os) const
    //{
    //    os.write(x);
    //    os.write(y);
    //    os.write(z);
    //}

    friend ostream& operator<<(ostream& os, const point3& pt)
    {
        os << pt.x << ' ' << pt.y << ' ' << pt.z;
        return os;
    }

};


typedef point3<int> ipoint3;
//template<> point3<int> point3<int>::zero;


// line3<T>
///////////////////////////////////////////////////////////////////////////////

template<class T> struct line3 {
    point3<T> a, b; 
};

typedef line3<int> iline3;

// int_set
///////////////////////////////////////////////////////////////////////////////

class int_set {
protected:
    vector<int> v;
    int s;
public:
    int_set() : v() { }
    int_set(int capacity) : v(capacity, 0) { }
    
    void set_capacity(int capacity) { v.resize(capacity, 0); }
    void clear() { fill(v.begin(), v.end(), 0); s = 0; }
    void insert(int i) { if (v[i] == 0) { v[i] = 1; ++s; } }
    void add_range(int min, int max) { while (min <= max) insert(min++); }
    void erase(int i) { if (v[i] > 0) { v[i] = 0; --s; } }
    bool find(int i) const { return v[i] > 0; }
    int size() const { return s; }
    int capacity() { return (int)v.size(); }
    template<class I> void assign(I begin, I end)
    {
        clear();
        while (begin != end) insert(*begin++);
    }
    template<class I> void insert(I begin, I end) 
    { 
        while (begin != end) insert(*begin++);
    }
    template<class C> void insert_to(C& container) const
    {
        int index;
        vector<int>::const_iterator iter;

        for (iter = v.begin(), index = 0; iter != v.end(); ++iter, ++index) {
            if (*iter > 0) container.push_back(index);
        }
    }
    template<class C> void insert_complement_to(C& container) const
    {
        int index;
        vector<int>::const_iterator iter;

        for (iter = v.begin(), index = 0; iter != v.end(); ++iter, ++index) {
            if (*iter == 0) container.push_back(index);
        }
    }


};


// subset_builder
///////////////////////////////////////////////////////////////////////////////

template<class T> class subset_creator {
public:
    typedef typename std::list<std::list<T> > container_t;
    typedef typename std::list<std::list<T> >::iterator iter_t;

    container_t subsets;

    subset_creator() : subsets() { subsets.push_back(list<T>()); }
    template<class Cond> void add(const T& e, const Cond& condition) 
    {
        iter_t iter = subsets.begin(), end = --subsets.end();

        do {
            list<T>& s = *iter;

            if (condition(s, e)) {
                subsets.push_back(s);
                subsets.back().push_back(e);
            }
        } while (iter++ != end);
    } 

};

template<class I, class S> void map_sum(I begin, I end, const S& s)
{
    while (begin != end) {
        *begin += s;
        ++begin;
    }
}

template<class I, class S> void map_diff(I begin, I end, const S& s)
{
    while (begin != end) {
        *begin -= s;
        ++begin;
    }
}


template<class T> struct matrix22
{
	T a;
	T b;
	T c;
	T d;

    matrix22() : a(0), b(0), c(0), d(0) { }
	matrix22(const T e11, const T e12, const T e21, const T e22) : a(e11), b(e12), c(e21), d(e22) { }
	matrix22(const matrix22& m) : a(m.a), b(m.b), c(m.c), d(m.d) { }
	
	T tr()
	{
		return a + d;
	}
	T det()
	{
		return a*d-b*c;
	}

	/*matrix22 inv()
	{
		return (matrix22<T>(d,-b,-c,a) / (det() + numeric_limits<T>::epsilon()));
	}*/

	matrix22 inv()
	{
		T f = det();
		return matrix22<T>(d/f,-b/f,-c/f,a/f);
	}
	matrix22 inv_safe(T add_to_det)
	{
		T f = det() + add_to_det;
		return matrix22<T>(d/f,-b/f,-c/f,a/f);
	}

	matrix22 t()
	{
		return matrix22<T>(a,c,b,d);
	}
	
	matrix22& operator*=(T f) { a*=f; b*=f; c*=f; d*=f; return *this; }

	void eig(point2<T>& v1, point2<T>& v2, T& l1, T& l2)
	{
		T t = tr();
		T tmp = sqrt(t*t/ 4 - det());
		//eigenvalues
		l1 = t/2 + tmp;
		l2 = t/2 - tmp;
		if(l1<l2)
		{
			T tmp = l1;
			l1 = l2;
			l2 = tmp;
		}
		//eigenvectors
		if(c>0)
		{
			v1.x = l1 - d;
			v1.y = c;
			v1 /= sqrt(v1.norm2());
			v2.x = l2 -d;
			v2.y = c;
			v2 /= sqrt(v2.norm2());
		}
		else if(b>0)
		{
			v1.x = b;
			v1.y = l1-a;
			v1 /= sqrt(v1.norm2());
			v2.x = b;
			v2.y = l2-a;
			v2 /= sqrt(v2.norm2());
		}
		else
		{
			v1.x = 1;
			v1.y = 0;
			v2.x = 0;
			v2.y = 1;
		}
	}
	matrix22 operator/(const T& f) const { return matrix22(a/f, b/f, c/f, d/f); }
	matrix22& operator/=(const T& f) { a/=f; b/=f; c/=f; d/=f; return *this; }
	matrix22 operator*(const matrix22& m) const { return matrix22(a*m.a+b*m.c, a*m.b+b*m.d, c*m.a+d*m.c, c*m.b+d*m.d); }
	
};

// ellipse<T>
///////////////////////////////////////////////////////////////////////////////
//enum axis_tip
//{
//	MAJOR_P,
//	MAJOR_N,
//	MINOR_P,
//	MINOR_N
//};
//
//template<class T> struct ellipse {
//	point2<T> c;
//	point2<T> eigv1;	//eigenvectors
//	point2<T> eigv2;
//	T l1;				//eigenvalues (l1 ~ eigv1) , l1 >= l2
//	T l2;
//	T aa;
//	T bb;
//	T cc;
//	int n;
//    /*ellipse(const point2<T>& pc = point2<T>(0,0), const point2<T>& peigv1 = point2<T>(0,0), const point2<T>& peigv2 = point2<T>(0,0),
//		const T& pl1 = 0, const T& pl2 = 0):
//        c(pc),  eigv1(peigv1), eigv2(peigv2), l1(pl1), l2(pl2){ }*/
//	ellipse(const point2<T>& pc = point2<T>(0,0), const point2<T>& peigv1 = point2<T>(0,0), const point2<T>& peigv2 = point2<T>(0,0),
//		const T& pl1 = 0, const T& pl2 = 0, const T& paa = 0, const T& pbb = 0, const T& pcc = 0, const int pn = 0):
//        c(pc),  eigv1(peigv1), eigv2(peigv2), l1(pl1), l2(pl2), aa(paa), bb(pbb), cc(pcc), n(pn) { }
//	template<class I> static ellipse<T> ellipse_from_pca(I begin, I end) 
//    {
//        if (begin == end) return ellipse<T>();
//		
//		point2<T> center = point2<T>::center(begin, end);
//        //covariance matrix
//		//T mat[2][2];
//		///*00 01
//		//10 11*/
//		//memset(mat, 0, sizeof(T) * 4);
//		//int n = 0;
//		//for(I iter = begin; iter != end; ++iter)
//		//{
//		//	point2<T> p = *iter - center;
//		//	mat[0][0] += p.x * p.x;
//		//	mat[0][1] += p.x * p.y;
//		//	mat[1][1] += p.y * p.y;
//		//	n++;
//		//}
//		//mat[0][0] /= n;
//		//mat[0][1] /= n;
//		//mat[1][1] /= n;
//		//// *** minimum = 1 pixel
//		//if(mat[0][0] < (T)1.0) mat[0][0] = (T)1.0;
//		//if(mat[0][1] < (T)1.0) mat[0][1] = (T)1.0;
//		//if(mat[1][1] < (T)1.0) mat[1][1] = (T)1.0;
//		//T tr = mat[0][0] + mat[1][1];
//		//T det = mat[0][0]*mat[1][1] - mat[0][1]*mat[0][1];
//		//T tmp = sqrt(tr*tr / 4 - det);
//		////eigenvalues
//		//T l1 = tr/2 + tmp;
//		//T l2 = tr/2 - tmp;	
//		////eigenvectors
//		//point2<T> eigv1;
//		//point2<T> eigv2;
//		//if(mat[0][1])
//		//{
//		//	eigv1.x = l1 - mat[1][1];
//		//	eigv1.y = mat[0][1];
//		//	eigv1 /= sqrt(eigv1.norm2());
//		//	eigv2.x = l2 - mat[1][1];
//		//	eigv2.y = mat[0][1];
//		//	eigv2 /= sqrt(eigv2.norm2());
//		//}
//		//else
//		//{
//		//	eigv1.x = 1;
//		//	eigv1.y = 0;
//		//	eigv2.x = 0;
//		//	eigv2.y = 1;
//		//}
//		//l1 = sqrt(l1);
//		//l2 = sqrt(l2);
//		//if(l2 > l1)
//		//{
//		//	T tmpt = l1;
//		//	l1 = l2;
//		//	l2 = tmpt;
//		//	point2<T> tmpv = eigv1;
//		//	eigv1 = eigv2;
//		//	eigv2 = tmpv;
//		//}
//		////return ellipse(center, eigv1, eigv2, l1, l2, mat[0][0], mat[0][1], mat[1][1], n);
//		matrix22<float> m;
//		int n = 0;
//		for(I iter = begin; iter != end; ++iter)
//		{
//			point2<T> p = *iter - center;
//			m.a += p.x * p.x;
//			m.b += p.x * p.y;
//			m.d += p.y * p.y;
//			n++;
//		}
//		m /= (T)n;
//		m.c = m.b;
//		// *** minimum = 1 pixel
//		point2<float> v1;
//		point2<float> v2;
//		T l1;
//		T l2;
//		m.eig(v1, v2, l1, l2);
//		return ellipse(center, v1, v2, sqrt(l1), sqrt(l2), m.a ,m.b, m.d, n);
//    }
//	
//	void scale_transf(point2<float> scale)
//	{
//		c *= scale;
//		eigv1 *= scale;
//		eigv2 *= scale;
//		float se1 = sqrt(eigv1.norm2());
//		float se2 = sqrt(eigv2.norm2());
//		l1 *= se1;
//		l2 *= se2;
//		eigv1 *= se1;
//		eigv2 *= se2;
//		aa *= scale.x*scale.x;
//		bb *= scale.x*scale.y;
//		cc *= scale.y*scale.y;
//	}
//
//	point2<T> axis_tip(axis_tip type) const
//	{
//		switch(type)
//		{
//		case MAJOR_P:
//			return c + eigv1 * l1;
//		case MAJOR_N:
//			return c - eigv1 * l1;
//		case MINOR_P:
//			return c + eigv2 * l2;
//		case MINOR_N:
//			return c - eigv2 * l2;
//		}
//		return point2<T>();
//	}
//	T a() const
//	{
//		return l1;
//	}
//	T b() const
//	{
//		return l2;
//	}
//
//	//bool in_bounding_box(rectangle2<T> r)
//	//{
//	//	T angx = atan(bb/aa);
//	//	T angy = atan(cc/bb);
//	//	point2<T> m(cos(angx) * aa + sin(angx) * bb,
//	//		cos(angy) * bb + sin(angy) * cc);
//	//	return c.x - m.x > r.ll.x && c.x + m.x < r.ur.x && c.y - m.y > r.ll.y && c.y + m.y < r.ur.y;
//	//}
//	rectangle2<T> bounding_box()
//	{
//		//T angx = atan(bb/aa);
//		//T angy = atan(cc/bb);
//		//point2<T> m(sqrt(abs(cos(angx) * aa + sin(angx) *bb)),
//		//	sqrt(abs(cos(angy) * bb + sin(angy) * cc)));
//		//return rectangle2<T>(c - m, c + m);
//		matrix22<T>m(aa,bb,bb,cc);
//		m = m.inv();
//		point2<T> di(sqrt(m.d / (m.a * m.d - m.b * m.b)), sqrt(m.a / (m.a * m.d - m.b * m.b)));
//		return rectangle2<T>(c - di, c + di);
//	}
//
//	bool in_bounding_box(rectangle2<T>& r)
//	{
//		matrix22<T>m(aa,bb,bb,cc);
//		m = m.inv();
//		point2<T> di(sqrt(m.d / (m.a * m.d - m.b * m.b)), sqrt(m.a / (m.a * m.d - m.b * m.b)));
//		//find((feat1(:,1)+feat1(:,8))<im1x & (feat1(:,1)-feat1(:,8))>0 & (feat1(:,2)+feat1(:,9))<im1y & (feat1(:,2)-feat1(:,9))>0);
//		return c.x - di.x > r.ll.x && c.x + di.x < r.ur.x && c.y - di.y > r.ll.y && c.y + di.y < r.ur.y;
//	}
//
//	void merge(ellipse ell, T weight)
//	{
//		c = c + (ell.c - c) * weight;
//		//weight for eigenvectors is approx.
//		eigv1 = eigv1 + (ell.eigv1 - eigv1) * weight;
//		eigv1 /= sqrt(eigv1.norm2());
//		eigv2 = eigv2 + (ell.eigv2 - eigv2) * weight;
//		eigv2 /= sqrt(eigv2.norm2());
//		l1 = l1 + (ell.l1 - l1) * weight;
//		l2 = l2 + (ell.l2 - l2) * weight;
//	}
//
//};

// operations on vectors
//////////////////////////

template<class T> vector<T>& operator+=(vector<T>& v1, const vector<T>& v2)
{
	typename vector<T>::iterator i1 = v1.begin();
	typename vector<T>::const_iterator i2 = v2.begin();
	size_t max = min(v1.size(), v2.size());

	for (size_t i = 0; i < max; ++i) {
		*i1 += *i2;
		++i1; ++i2;
	}
	return v1;
}

template<class T> vector<T>& operator+=(vector<T>& v, const T& t)
{
    for (typename vector<T>::iterator i = v.begin(); i != v.end(); ++i) {
		*i += t;
	}
	return v;
}

template<class T> vector<T>& operator-=(vector<T>& v1, const vector<T>& v2)
{
	typename vector<T>::iterator i1 = v1.begin();
	typename vector<T>::const_iterator i2 = v2.begin();
	size_t max = min(v1.size(), v2.size());

	for (size_t i = 0; i < max; ++i) {
		*i1 -= *i2;
		++i1; ++i2;
	}
	return v1;
}

template<class T> vector<T>& operator-=(vector<T>& v, const T& t)
{
    for (typename vector<T>::iterator i = v.begin(); i != v.end(); ++i) {
		*i -= t;
	}
	return v;
}

template<class T> vector<T>& operator/=(vector<T>& v1, const vector<T>& v2)
{
	typename vector<T>::iterator i1 = v1.begin();
	typename vector<T>::const_iterator i2 = v2.begin();
	size_t max = min(v1.size(), v2.size());

	for (size_t i = 0; i < max; ++i) {
		*i1 /= *i2;
		++i1; ++i2;
	}
	return v1;
}

template<class T> vector<T>& operator/=(vector<T>& v, const T& t)
{
	for (typename vector<T>::iterator i = v.begin(); i != v.end(); ++i) {
		*i /= t;
	}
	return v;
}

template<class T> vector<T> operator+(const vector<T>& v1, const vector<T>& v2)
{
	typename vector<T>::const_iterator i1 = v1.begin();
	typename vector<T>::const_iterator i2 = v2.begin();
	size_t max = min(v1.size(), v2.size());
	vector<T> result(max);
	typename vector<T>::iterator i3 = result.begin();

	for (size_t i = 0; i < max; ++i) {
		*i3 = *i1 + *i2;
		++i1; ++i2; ++i3;
	}
	return result;
}

template<class T> vector<T> operator-(const vector<T>& v1, const vector<T>& v2)
{
	typename vector<T>::const_iterator i1 = v1.begin();
	typename vector<T>::const_iterator i2 = v2.begin();
	size_t max = min(v1.size(), v2.size());
	vector<T> result(max);
	typename vector<T>::iterator i3 = result.begin();

	for (size_t i = 0; i < max; ++i) {
		*i3 = *i1 - *i2;
		++i1; ++i2; ++i3;
	}
	return result;
}

template<class T> void sqr(vector<T>& v)
{
	for (typename vector<T>::iterator i = v.begin(); i != v.end(); ++i) {
		*i *= *i;
	}
}

template<class T> void sqrt(vector<T>& v)
{
	for (typename vector<T>::iterator i = v.begin(); i != v.end(); ++i) {
		*i = sqrt(*i);
	}
}

template<class First, class Second> ostream& operator<<(ostream& os, const pair<First,Second>& p)
{
    os << '(' << p.first << ',' << p.second << ')';
    return os;
}

template<class T> ostream& operator<<(ostream& os, const vector<T>& v)
{
    os << (int)v.size();
    for (typename vector<T>::const_iterator i = v.begin(); i != v.end(); ++i) {
        os << ' ' << *i;
    }
    return os;
}

template<class T> istream& operator>>(istream& is, vector<T>& v)
{
    int size;

    is >> size;
    v.resize(size);
    for (typename vector<T>::iterator i = v.begin(); i != v.end(); ++i) {
        is >> *i;
    }
    return is;
}


// various operations on containers
/////////////////////////////////////

template<class T> T& at(list<T>& l, int i) 
{
    for (typename list<T>::iterator iter = l.begin(); iter != l.end() && i >= 0; ++iter, --i) 
        if (i == 0) return *iter;
    throw exception();
}

template <typename I, class F> void map_through(I begin, I end, const F& f)
{
    while (begin != end) {
        f(*begin);
        ++begin;
    }
}


#endif /* _STRUCTURES_H_ */




