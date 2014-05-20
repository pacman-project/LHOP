// graph edge classes
///////////////////////////////////////////////////////////////////////////////

#pragma once
#ifndef _GRAPH_EDGE_H_
#define _GRAPH_EDGE_H_

#include "utils/serialization/streaming.h"
#include "utils/smemory.h"

using namespace std;

class node;

class edge_data : public streamable, public custom_alloc {
public:
    edge_data() { }

    virtual void clone(const edge_data* d) { }
    edge_data* get_clone() const
    { 
        edge_data* res = new edge_data(); 
        res->clone(this); 
        return res;
    }
    virtual void gc_clone(const edge_data* d, const map<node*,node*>& cmap) { }
    edge_data* get_gc_clone(const map<node*, node*>& cmap) const
    {
        edge_data* res = (edge_data*)make_instance();
        res->gc_clone(this, cmap);
        return res;
    }

    virtual void copy_to(streamable* p, cloner& cl) { /* streamable::copy_to(p, cl); */ }
    virtual streamable* make_instance() const { return new edge_data(); }
    virtual void read_from_stream(istreamer& is) { streamable::read_from_stream(is); }
    virtual void write_to_stream(ostreamer& os) { streamable::write_to_stream(os); }

};

template<class T> class edge_data_t : public edge_data {
public:
    T data;

    edge_data_t() { }
    edge_data_t(const T& d) : data(d) { }

    virtual void clone(const edge_data_t* d) { data = ((edge_data_t*)d)->data; }

    virtual void copy_to(streamable* p, cloner& cl) 
    {
        edge_data::copy_to(p, cl);
        ((edge_data_t*)p)->data = data;
    }

    virtual streamable* make_instance() const { return new edge_data_t<T>(); }
    virtual void read_from_stream(istreamer& is) 
    { 
        edge_data::read_from_stream(is);
        is.read(data);
    }

    virtual void write_to_stream(ostreamer& os) 
    { 
        edge_data::write_to_stream(os); 
        os.write(data);
    }

    virtual ostream& print(ostream& os) const
    {
        os << data << endl;
        return os;
    }

};

template<class T> class edge_data_ts : public edge_data { 
public:
    T data;

    edge_data_ts() { }
    edge_data_ts(const T& d) : data(d) { }

    virtual void clone(const edge_data_ts* d) { data = ((edge_data_ts*)d)->data; }

    virtual void copy_to(streamable* p, cloner& cl) 
    {
        edge_data::copy_to(p, cl);
        ((edge_data_ts*)p)->data = (T)cl.get_copy((streamable*)data);
    }

    virtual streamable* make_instance() const { return new edge_data_ts<T>(); }
    virtual void read_from_stream(istreamer& is) 
    { 
        edge_data::read_from_stream(is);
        is.read((streamable*&)data);
    }

    virtual void write_to_stream(ostreamer& os) 
    { 
        edge_data::write_to_stream(os); 
        os.write((streamable*)data);
    }
};

template<class T> class edge_data_tw : public edge_data { 
public:
    T data;

    edge_data_tw() { }
    edge_data_tw(const T& d) : data(d) { }

    virtual void clone(const edge_data_tw* d) { data = ((edge_data_tw*)d)->data; }

    virtual void copy_to(streamable* p, cloner& cl) 
    {
        edge_data::copy_to(p, cl);
        ((edge_data_tw*)p)->data = data;
    }

    virtual streamable* make_instance() const { return new edge_data_tw<T>(); }
    virtual void read_from_stream(istreamer& is) 
    { 
        edge_data::read_from_stream(is);
        data.read_from_stream(is);
    }

    virtual void write_to_stream(ostreamer& os) 
    { 
        edge_data::write_to_stream(os); 
        data.write_to_stream(os);
    }
};


#endif /* _GRAPH_EDGE_H_ */

