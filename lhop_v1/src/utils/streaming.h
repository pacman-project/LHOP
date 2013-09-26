// support for storing/restoring classes to/from stream

#pragma once
#ifndef _STREAMING_H_
#define _STREAMING_H_


#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <stack>
#if defined WIN32 | defined WIN64
#include <typeinfo.h>
#else
#include <typeinfo>
#endif
#include "structures.h"
#include "matrix.h"
#include "zlib_compress.h"
#include "exceptions.h"
#include <iomanip>

using namespace std;

// classes
///////////////////////////////////////////////////////////////////////////////

class ostreamer;
class istreamer;
class streamable;
class cloner;

class streaming_exception : public libhop_exception {
public:
    enum TYPE { UNKNOWN_ERROR = 0, CANNOT_OPEN_FILE, CLASS_NOT_REGISTERED, POINTER_ID_NOT_FOUND, NONSTREAMABLE_CLASS };

    TYPE type;

	streaming_exception(const string& file, const size_t line, TYPE t, const string s);

	streaming_exception(const string& file, const size_t line, streamable* c);
	streaming_exception(const string& file, const size_t line, TYPE t);
	streaming_exception(const string& file, const size_t line);

    void write();
private:
	string parse_error(TYPE t, const string s);
};

class class_id_store {
private:
    vector<streamable*> id_map;
    map<string, int> name_map;

public:    
    class_id_store() : id_map(), name_map() { }
    ~class_id_store();

    int get_id(streamable* c);
    streamable* get_instance(int id) { return (id >= 0 && id < (int)id_map.size()) ? id_map[id] : nullptr; }
    void register_class(streamable* c);
};

class streamable {
	
protected:
    static class_id_store _class_id_storage;

    streamable() { }
public:
    virtual ~streamable() { }

	// desended classes MUST implement following methods
    virtual streamable* make_instance() const = 0;
	virtual void read_from_stream(istreamer& is) {}
	virtual void write_to_stream(ostreamer& os) {}
    virtual void copy_to(streamable* p, cloner& cl) 
    { 
        throw streaming_exception(__FILE__, __LINE__, streaming_exception::UNKNOWN_ERROR, "Method 'copy_to' not implemented"); 
    }

	virtual void save(const string& file, int zlib_compression_level = Z_NO_COMPRESSION);

    static streamable* read(const string& file);

    static streamable* read_instance(istreamer& is);

    // Makes a copy by writing to - and reading from a stream
    // Should not be used if performance is critical.
    streamable* get_copy_s();

    streamable* get_copy();

    static class_id_store& get_class_id_store() { static class_id_store _class_id_storage; return _class_id_storage; }

    int get_type() { return get_class_id_store().get_id(this); }
};

// cloner
///////////////////////////////////////////////////////////////////////////////

class cloner {
protected:
    map<streamable*, streamable*> pointer_dict;

public:
    cloner();
    streamable* get_copy(streamable* p);
    static streamable* clone_streamable(streamable* p);
};

// ostreamers
///////////////////////////////////////////////////////////////////////////////

// ostreamer
//////////////

class ostreamer : public ostream {
private:
	bool do_compression;
	zlib_compress* z_compress_stream;
	

protected:
    map<streamable*, int> pointer_dict;
    stack<streamable*> obj_stack;

public:
    ostreamer(basic_streambuf<char>* buf, int zlib_compress_level = Z_NO_COMPRESSION);
    virtual ~ostreamer() {
		if (z_compress_stream != nullptr)
			delete z_compress_stream;
	}
    
    virtual void write_header();

    // should be used to save the whole "streamable" structure
    void write_structure(streamable*);

    // write_* methods should be called within "write_to_stream" methods only
    // to save the whole structure use the "write_structure" method.
    void write(streamable*);
    void write(const char* ptr, streamsize size) { 
		// should check if compression is needed and modify data acordningly 
		if (do_compression == false) {
			// write directly to stream
			ostream::write(ptr, size); 
			//ostream::flush();
		} else {
			// write to buffer and compress 
			z_compress_stream->write_to_buf_and_compress(ptr, size);
		}
		
	}
	/*void write(const bool& b) { 
		stringstream ss;
		ss << b << endl;
		write(ss.str().c_str(), sizeof(char) * ss.str().size()); 
	}
    void write(const int& i)  { 
		stringstream ss;
		ss << i << endl;
		write(ss.str().c_str(), sizeof(char) * ss.str().size()); 
	}
    void write(const unsigned& u)  { 
		stringstream ss;
		ss << u << endl;
		write(ss.str().c_str(), sizeof(char) * ss.str().size()); 
	}
    void write(const char& c)  { 
		stringstream ss;
		ss << c << endl;
		write(ss.str().c_str(), sizeof(char) * ss.str().size()); 
	}
    void write(const double& d)  { 
		stringstream ss;
		ss << setprecision(3) << d << endl;
		write(ss.str().c_str(), sizeof(char) * ss.str().size()); 
	}
    void write(const string& s) 
    { 
		write(s.data(), s.size()*sizeof(char));
    }*/
    void write(const bool& b) { write((char*)&b, sizeof(bool)); }
    void write(const int& i) { write((char*)&i, sizeof(int)); }
    void write(const unsigned& u) { write((char*)&u, sizeof(unsigned)); }
    void write(const char& c) { write((char*)&c, sizeof(char)); }
    void write(const double& d) { write((char*)&d, sizeof(double)); }
    void write(const float& d) { write((char*)&d, sizeof(float)); }
    void write(const string& s) 
    { 
        unsigned len = (unsigned)s.length();
        write(len);
        write(s.data(), len*sizeof(char));
    }
    template<class T> void write(const point2<T>& p) { write(p.x); write(p.y); }
    template<class T, class W> void write(const wpoint2<T, W>& p) { write(p.x); write(p.y); write(p.w); }
    template<class T> void write(const rectangle2<T>& r) { write(r.ll); write(r.ur); }
    //void write(const vector<streamable*>& v);
    template<class T> void write(const vector<T>& v);
    template<class T> void write(const vector<T*>& v);
    template<class T> void write(const list<T>& l);
    template<class T> void write(const matrix<T>& m);
    template<class T> void write(const matrix<T*>& m);
    template<class T> void write(const set<T>& s);
    template<class T, class S> void write(const map<T, S>& m);
    template<class T> void write(const T[], int);
    template<class T, class S> void write(const pair<T, S>& p) 
        { write(p.first); write(p.second); }
    template<class T> void write_streamable_p(const T[], int);
    template<class T> void write_streamable_p(const vector<T>& v);
    template<class T, class S> void write_streamable_p(const set<T, S>& s);
    template<class T> void write_streamable(const vector<T>& v);
    template<class T> void write_streamable(const list<T>& l);

};

// ofstreamer
///////////////

class ofstreamer : public ostreamer {
public:
    ofstreamer(int zlib_compression_level = Z_NO_COMPRESSION);
    virtual ~ofstreamer() { }
 
    void open(const char*);
    void open(const string& s) { open(s.c_str()); }
    void close();

protected:
    basic_filebuf<char> fbuf;
};

// ofzstreamer
////////////////
//
//class ofzstreamer : public ostreamer {
//public:
//    ofzstreamer();
// 
//    void open(const char*);
//    void open(const string& s) { open(s.c_str()); }
//    void close();
//
//protected:
//    gzfilebuf fbuf;
//};

// ostrstreamer
/////////////////

class ostrstreamer : public ostreamer {
public:
    ostrstreamer(int zlib_compress_level = Z_NO_COMPRESSION);
    virtual ~ostrstreamer() { }
    
    //streamsize pcount() const {	return sbuf.pcount(); }
    stringbuf* rdbuf() const { return (stringbuf*)&sbuf; }
protected:
    stringbuf sbuf;
};

// istreamers
///////////////////////////////////////////////////////////////////////////////

// istreamer
//////////////

class istreamer : public istream {
private:
	bool do_decompression;
	zlib_compress* z_decompress_stream;

protected:
    double version;

    map<int, streamable*> ptr_store;
    stack<streamable*> obj_stack;
public:
    istreamer(basic_streambuf<char>* buf);
	virtual ~istreamer();

    double get_version() const { return version; }
    virtual void read_header();

    // should be used to read the whole "streamable" structure
    streamable* read_structure();

    // read_* methods should be called within "read_from_stream" methods only
    // to read the whole structure use the "read_structure" method.
    void read(streamable*&);
    void read(char* ptr, streamsize size) { 
		// should compress if needed 
		if (do_decompression == false) { 
			// use normal compress
			istream::read(ptr, size); 
		} else {
			// decompress data, output data needed is "size"
			z_decompress_stream->read_to_buf_and_decompress(ptr, size);
		}
	}
    void read(bool& b) { read((char*)&b, sizeof(bool)); }
    void read(int& i) { read((char*)&i, sizeof(int)); }
    void read(unsigned& u) { read((char*)&u, sizeof(unsigned)); }
    void read(char& c) { read((char*)&c, sizeof(char)); }
    void read(double& d) { read((char*)&d, sizeof(double)); }
    void read(float& d) { read((char*)&d, sizeof(float)); }
    void read(string& s);
    template<class T> void read(point2<T>& p) { read(p.x); read(p.y); }
    template<class T, class W> void read(wpoint2<T, W>& p) { read(p.x); read(p.y); read(p.w); }
    template<class T> void read(rectangle2<T>& r) { read(r.ll); read(r.ur); }
    //void read(vector<streamable*>& v);
    template<class T> void read(vector<T>& v);
    template<class T> void read(list<T>& l);
    template<class T> void read(T v[]);
    template<class T> void read(vector<T*>& v);
    template<class T> void read(matrix<T>& m);
    template<class T> void read(matrix<T*>& m);
    template<class T> void read(set<T>& s);
    template<class T, class S> void read(map<T, S>& m);
    template<class T, class S> void read(pair<T,S>& p)
        { read(p.first); read(p.second); }
    template<class T> void read_streamable_p(T v[]);
    template<class T> void read_streamable_p(vector<T>& v);
    template<class T, class S> void read_streamable_p(set<T, S>& s);
    template<class T> void read_streamable(vector<T>& v);
    template<class T> void read_streamable(list<T>& l);

};

// ifstreamer
///////////////

class ifstreamer : public istreamer {
public:
    ifstreamer();

    void open(const char*);
    void open(const string& s) { open(s.c_str()); }

    void close();

protected:
    basic_filebuf<char> fbuf;
};

// ifzstreamer
////////////////
//
//class ifzstreamer : public istreamer {
//public:
//    ifzstreamer();
//
//    void open(const char*);
//    void open(const string& s) { open(s.c_str()); }
//
//    void close();
//
//protected:
//    gzfilebuf fbuf;
//};

class istrstreamer : public istreamer {
public: 
    istrstreamer(stringbuf* buf);

    //streamsize pcount() const {	return sbuf->pcount(); }
    stringbuf* rdbuf() const { return sbuf; }
protected:
    stringbuf* sbuf;
};


// ostreamer templates
///////////////////////////////////////////////////////////////////////////////

template<class T> void ostreamer::write(const std::vector<T>& v)
{
    typename std::vector<T>::const_iterator iter;

    write((unsigned)v.size());
    for (iter = v.begin(); iter != v.end(); ++iter) {
        write(*iter);
    }
}

template<class T> void ostreamer::write(const std::list<T>& l)
{
    typename std::list<T>::const_iterator iter;

    write((unsigned)l.size());
    for (iter = l.begin(); iter != l.end(); ++iter) {
        write(*iter);
    }
}

template<class T> void ostreamer::write(const std::set<T>& s)
{
    write((unsigned)s.size());
    for (typename set<T>::const_iterator iter = s.begin(); iter != s.end(); ++iter) {
        write(*iter);
    }
}

template<class T, class S> void ostreamer::write(const map<T, S>& m)
{
    write((unsigned)m.size());
    for (typename map<T, S>::const_iterator iter = m.begin(); iter != m.end(); ++iter) {
        write(iter->first);
        write(iter->second);
    }
}

template<class T> void ostreamer::write(const T v[], int count)
{
    write(count);
    for (int i = 0; i < count; ++i) {
        write(v[i]);
    }
}

template<class T> void ostreamer::write(const std::vector<T*>& v)
{
    typename std::vector<T*>::const_iterator iter;

    write((unsigned)v.size());
    for (iter = v.begin(); iter != v.end(); ++iter) {
        write(**iter);
    }
}

template<class T> void ostreamer::write(const matrix<T>& m)
{
    write((unsigned)m.width);
    write((unsigned)m.height);
    for_each_element(m, i) {
        write(m[i]);
    }
}

template<class T> void ostreamer::write(const matrix<T*>& m)
{
    write((unsigned)m.width);
    write((unsigned)m.height);
    for_each_element(m, i) {
        write(*m[i]);
    }
}

template<class T> void ostreamer::write_streamable_p(const T v[], int count)
{
    write(count);
    for (int i = 0; i < count; ++i) {
        write((streamable*)(v[i]));
    }
}

template<class T> void ostreamer::write_streamable_p(const vector<T>& v)
{
    write((unsigned)(v.size()));
    for (typename vector<T>::const_iterator iter = v.begin(); iter != v.end(); ++iter)
        write((streamable*)(*iter));
}

template<class T, class S> void ostreamer::write_streamable_p(const set<T, S>& s)
{
    write((unsigned)(s.size()));
    for (typename set<T, S>::const_iterator iter = s.begin(); iter != s.end(); ++iter)
        write((streamable*)(*iter));
}

template<class T> void ostreamer::write_streamable(const vector<T>& v)
{
    write((unsigned)(v.size()));
    for (typename vector<T>::const_iterator iter = v.begin(); iter != v.end(); ++iter)
        iter->write_to_stream(*this);
}

template<class T> void ostreamer::write_streamable(const list<T>& l)
{
    write((unsigned)(l.size()));
    for (typename list<T>::const_iterator iter = l.begin(); iter != l.end(); ++iter)
        iter->write_to_stream(*this);
}


// istreamer templates
///////////////////////////////////////////////////////////////////////////////

template<class T> void istreamer::read(T v[])
{
    int size;
    
    read(size);
    for (int i = 0; i < size; ++i) {
        read(v[i]);
    }
}

template<class T> void istreamer::read(std::vector<T>& v)
{
    unsigned size;
    
    read(size);
    v.resize(size);
    for (unsigned i = 0; i < size; ++i) {
        read(v[i]);
    }
}

template<class T> void istreamer::read(std::list<T>& l)
{
    unsigned size;
    
    read(size);
    for (unsigned i = 0; i < size; ++i) {
        T item;

        read(item);
        l.push_back(item);
    }
}

template<class T> void istreamer::read(set<T>& s)
{
    unsigned size;
    T item;

    read(size);
    for (unsigned i = 0; i < size; ++i) {
        read(item);
        s.insert(item);
    }
}
 
template<class T, class S> void istreamer::read(map<T, S>& m)
{
    unsigned size;
    T key;
    S value;

    read(size);
    for (unsigned i = 0; i < size; ++i) {
        read(key); read(value);
        m.insert(pair<T, S>(key, value));
    }
}


template<class T> void istreamer::read(std::vector<T*>& v)
{
    unsigned size;
    T* item;
    
    read(size);
    v.resize(size);
    for (unsigned i = 0; i < v.size(); ++i) {
        item = new T;
        read(*item);
        v[i] = item;
    }
}

template<class T> void istreamer::read(matrix<T>& m)
{
    unsigned width, height;
    T item;
    
    read(width);
    read(height);
    m.resize((size_t)width, (size_t)height);
    for_each_element(m, i) { 
        read(item);
        m[i] = item;
        //read((T&)(m[i]));
    }
}

template<class T> void istreamer::read(matrix<T*>& m)
{
    unsigned width, height;
    T* item;
    
    read(width);
    read(height);
    m.resize(width, height);
    for_each_element(m, i) { 
        item = new T;
        read(*item);
        m[i] = item;
    }
}

template<class T> void istreamer::read_streamable_p(T v[])
{
    int size;
    
    read(size);
    for (int i = 0; i < size; ++i) {
        read(v[i]);
    }
}


template<class T> void istreamer::read_streamable_p(vector<T>& v)
{
    unsigned size;
    streamable* s;

    read(size);
    for (; size > 0; --size) {
        read(s);
        v.push_back((T)s);
    }
}

template<class T, class S> void istreamer::read_streamable_p(set<T, S>& s)
{
    unsigned size;
    streamable* p;

    read(size);
    for (; size > 0; --size) {
        read(p);
        s.insert((T)p);
    }
}

template<class T> void istreamer::read_streamable(vector<T>& v)
{
    unsigned size;

    read(size);
    v.resize(size);
    for (unsigned i = 0; i < size; ++i) {
        v[i].read_from_stream(*this);
    }
}

template<class T> void istreamer::read_streamable(list<T>& l)
{
    unsigned size;

    read(size);
    for (; size > 0; --size) {
        l.push_back(T());
        l.back().read_from_stream(*this);
    }
}


#endif /* _STREAMING_H_ */
