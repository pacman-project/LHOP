// class streamable definitions

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <typeinfo>
#include "streaming.h"
#include "utils/utils.h"
/*
set<string> & types() {
        static set<string> types; 
        return types;
} 
*/

// streaming_exception
///////////////////////////////////////////////////////////////////////////////

/// constructors should be here due to forward declaration of streamable
streaming_exception::streaming_exception(const string& file, const size_t line, TYPE t, const string s): libhop_exception(file, line, parse_error(t, s),"streaming_exception"), type(t) {}
streaming_exception::streaming_exception(const string& file, const size_t line, streamable* c):  libhop_exception(file, line, parse_error(CLASS_NOT_REGISTERED, typeid(*c).name()),"streaming_exception"), type(CLASS_NOT_REGISTERED) {}
streaming_exception::streaming_exception(const string& file, const size_t line, TYPE t): libhop_exception(file, line, parse_error(t, ""),"streaming_exception"), type(t) {}
streaming_exception::streaming_exception(const string& file, const size_t line): libhop_exception(file, line, parse_error(UNKNOWN_ERROR, ""),"streaming_exception"), type(UNKNOWN_ERROR) {}


string streaming_exception::parse_error(TYPE t, const string s) {
    switch (t) {
        case POINTER_ID_NOT_FOUND :
            return "Pointer id too big while reading from stream."; 
        case CANNOT_OPEN_FILE :
            return "Can not open file " + s;
        case CLASS_NOT_REGISTERED :
            return "Class type not registered for streaming: " + s;
        case NONSTREAMABLE_CLASS : 
            return "Class '" + s + "' can not be streamed.";
        default:
            return "Streaming exception." + s;
    }
}
void streaming_exception::write()
{
    cout << what() << endl;
}

// cloner
///////////////////////////////////////////////////////////////////////////////

cloner::cloner() : pointer_dict()
{
}

streamable* cloner::get_copy(streamable* p)
{
    if (p == nullptr) return nullptr;

    map<streamable*, streamable*>::iterator fiter = pointer_dict.find(p);

    if (fiter != pointer_dict.end()) return fiter->second;
    
    streamable* result = p->make_instance();

    pointer_dict.insert(pair<streamable*, streamable*>(p, result));
    p->copy_to(result, *this);
    return result;
}

streamable* cloner::clone_streamable(streamable *p)
{
    cloner cl;
    
    return cl.get_copy(p);
}



// streamable
///////////////////////////////////////////////////////////////////////////////

streamable* streamable::read(const string& file)
{
    streamable* result = nullptr;
    ifstreamer is;

    is.open(file);
    result = is.read_structure();
    is.close();
    return result;
}

streamable* streamable::read_instance(istreamer& is)
{
    int id;

    is.read(id);
    streamable* s = get_class_id_store().get_instance(id);

    return (s == nullptr) ? nullptr : s->make_instance();
}

streamable* streamable::get_copy_s()
{
    ostrstreamer os;
    os.write_structure(this);
    istrstreamer is(os.rdbuf());
    streamable* result = is.read_structure();
    return result;
}

streamable* streamable::get_copy()
{
    cloner cl;
    return cl.get_copy(this);
}

void streamable::save(const string& file, int zlib_compression_level)
{
    ofstreamer os(zlib_compression_level);
		
    os.open(file);
    os.write_structure(this);
    os.close();
}


// class_id_store
///////////////////////////////////////////////////////////////////////////////

class_id_store::~class_id_store()
{
    vector<streamable*>::iterator iter;
    for (iter = id_map.begin(); iter != id_map.end(); ++iter) 
        delete *iter;
}


int class_id_store::get_id(streamable* c)
{
    map<string, int>::iterator iter = name_map.find(typeid(*c).name());
    if (iter == name_map.end()) throw streaming_exception(__FILE__,  __LINE__, c);
    else return iter->second;
}

// o(f,str,fz)streamer
///////////////////////////////////////////////////////////////////////////////

// ofstreamer
///////////////

ofstreamer::ofstreamer(int zlib_compression_level) :
    ostreamer(&fbuf, zlib_compression_level)
{
}

void ofstreamer::open(const char* fname)
{
	// get folder name and create it
	create_directories(get_parent_of_path(fname));

    fbuf.open(fname, ios::out | ios::binary | ios::trunc);
    if (fail()) throw streaming_exception(__FILE__, __LINE__, streaming_exception::CANNOT_OPEN_FILE, fname);
}

void ofstreamer::close()
{
    fbuf.close();
}

// ostrstreamer
/////////////////

ostrstreamer::ostrstreamer(int zlib_compress_level /* = Z_NO_COMPRESSION */) :
    ostreamer(&sbuf, zlib_compress_level), sbuf(ios_base::in | ios_base::out)
{
}


// ostreamer
//////////////

ostreamer::ostreamer(basic_streambuf<char>* buf, int zlib_compress_level) : 
    ostream(buf),
    pointer_dict(),
    obj_stack(),
	do_compression(false)
{
	// construct compress stream even if it is not needed
	z_compress_stream = new zlib_compress(zlib_compress_level);
	z_compress_stream->set_ostream(this);
}

void ostreamer::write_header()
{
    // 1.0 -- MAX_LAYER_NUMBER = 8
    // 1.1 -- MAX_LAYER_NUMBER = 16
    // 1.2 -- field part_lib::layer_info added
    // 1.3 -- part_data serialization changed
    // 1.4 -- field part_lib::attr added
    // 2.0 -- field lib_data::td added (thresholds)
    // 2.1 -- field layer1_data::r (response_map)
    // 2.2 -- field part_data_2::gdistr (distribution for G_RESPONSE) added
	// 2.3 -- added compression of all versions 
	// header must be uncompressed so that it can be read 
    // 2.4 -- val removed from layer1_data
    // 2.5 -- structure edge_data_ip2 changed
    // 2.6 -- field var added to part_data
    // 2.7 -- field part_data_2::td added
    // 2.8 -- layer1_data::border and layer1_data::original size added
    // 3.0 -- "new" learning with integrated shape learning
    // 3.1 -- part_data_2a::geo changed
    // 3.2 -- pca_data changed (added sizefactor and normal)
    // 3.3 -- svm_data changed (added mean and sigma)
	// 3.4 -- app: int -> double ==> int -> (int, double)
    write((double)3.4); 
	if (z_compress_stream != nullptr) {
		write(z_compress_stream->get_compression_level());
	} else {		
		// just write zero if z_compress_stream is not define (probably something wrong so send info to err)
		write((int)0); 
		std::cout << "Warning durring streamable header writing. Trying to write compression level to header but z_compress_stream == nullptr. \n\tWill not use compression !!" << endl;
	}
}

void ostreamer::write_structure(streamable* c)
{
    streamable* o;

    pointer_dict.clear();
	
	// header MUST be write without compression
	do_compression = false; // just in case
    
	write_header();

	// if should use compression then set do_compression flag
	// so it can be used by basic write method
	if (z_compress_stream != nullptr && z_compress_stream->get_compression_level() != 0)
		do_compression = true;
	else
		do_compression = false;

	// initialize compression
	if (do_compression == true)
		z_compress_stream->begin_compression();

    write(c);

    while (!obj_stack.empty()) {
        o = obj_stack.top();
        obj_stack.pop();
        o->write_to_stream(*this);
    }
	// finalize compression
	if (do_compression == true)
		z_compress_stream->end_compression();
}

void ostreamer::write(streamable* c)
{
    if (c == nullptr) {
        write((int)-1);
        return;
    }

    pair<map<streamable*, int>::iterator, bool> pr;
    int id = (int)pointer_dict.size();

    pr = pointer_dict.insert(pair<streamable*, int>(c, id));
    if (pr.second) {
        write(id);
        write(c->get_type());
        obj_stack.push(c);        
    } else {
        write((pr.first)->second);
    }
}

// i(f,fz,str)streamer
///////////////////////////////////////////////////////////////////////////////

// ifstreamer
//////////////

ifstreamer::ifstreamer() :
    istreamer(&fbuf)
{ }    

void ifstreamer::open(const char* file)
{
    if (fbuf.open(file, ios::in | ios::binary) == 0)
        throw streaming_exception(__FILE__, __LINE__, streaming_exception::CANNOT_OPEN_FILE, file);
}

void ifstreamer::close()
{
    fbuf.close();
}

// istrstreamer
/////////////////

istrstreamer::istrstreamer(stringbuf* buf) :
    istreamer(buf), 
    sbuf(buf)
{
}

// istreamer
//////////////

istreamer::istreamer(basic_streambuf<char>* buf) :
    istream(buf),
    ptr_store(),
    obj_stack(),	
	do_decompression(false)
{ 
	// create z_decompress_stream since we will not 
	// know if we might need it until we read header (where compression info is stored)
	z_decompress_stream = new zlib_compress();
	z_decompress_stream->set_istream(this);
}

istreamer::~istreamer()
{
	if (z_decompress_stream != nullptr)
		delete z_decompress_stream;
}
void istreamer::read_header()
{
	int compression = 0;
	// reading header should always be uncompressed
    read(version);
	// read compression only if version 2.3 or newer
	if (version >= 2.3) {
		read(compression);
	}


	// after all headers have been read, set do_compression that is used in basic write method
	if (compression != 0) 
		do_decompression = true;
}

streamable* istreamer::read_structure()
{
    streamable* result = nullptr;
    streamable* o;

    ptr_store.clear();

	// header must be read without compression
	do_decompression = false; // just in case
    
	read_header();

	// initialize decompression
	if (do_decompression == true)
		z_decompress_stream->begin_compression(true);

    read((streamable*&)result);

    while (!obj_stack.empty()) {
        o = obj_stack.top();
        obj_stack.pop();
        if (o == nullptr) throw streaming_exception(__FILE__,  __LINE__, streaming_exception::POINTER_ID_NOT_FOUND);
        o->read_from_stream(*this);
    }

	// finalize decompression
	if (do_decompression == true)
		z_decompress_stream->end_compression(true);

    return result;
}

void istreamer::read(streamable*& p)
{
    int id;

    read(id);
    if (id < 0) p = nullptr;
    else {
        map<int, streamable*>::iterator iter = ptr_store.find(id);

        if (iter != ptr_store.end()) {
            p = iter->second;
        } else {
            p = streamable::read_instance(*this);
			ptr_store.insert(pair<int, streamable*>(id, p));
            obj_stack.push(p);
        }
    } 
}


void istreamer::read(string& s) 
{ 
    unsigned len;
    
    read(len);
    
    char* buf = new char[len + 1];
    
    read(buf, sizeof(char)*len);
    buf[len] = '\0';
    s.assign(buf);
    delete[] buf;
}

