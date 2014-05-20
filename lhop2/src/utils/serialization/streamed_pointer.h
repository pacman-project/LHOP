// Pointer to streamable object which is kept on disk

#pragma once

#ifndef __STREAMED_POINTER_H__
#define __STREAMED_POINTER_H__

#include <string>
#include "utils/utils.h"
#include "utils/serialization/streaming.h"

using namespace std;


// streamed_pointer
/////////////////////

/// Pointer to streamable object which is kept on disk!
/// Used by optimization process.
/// Pointer is also streamable so it can be send/recived to mapreduce implementation
struct streamed_pointer : public streamable {
	friend class hop_streamed_pointer;
protected:
    struct counted_name {
        int count;
        string name;
		bool disposable;
        counted_name(const string& n) : name(n), count(1), disposable(1) { }
        ~counted_name() { if (disposable) delete_file(name); }
    };

    counted_name* namep;
    static int nextid;
    static int pid; // process id.
    static string tmp_dir;	
public:
    streamed_pointer() : namep(nullptr){ }
    streamed_pointer(streamable* s);
    streamed_pointer(const streamed_pointer& sp) : namep(nullptr) { copy(sp); }
    ~streamed_pointer() { dispose(); }
    
	string get_name_only() { 
		string name = (namep != nullptr) ? namep->name : ""; 
		size_t pos = name.find_last_of("/\\");
		return pos != string::npos ? name.substr(pos+1) : name;
	}
    bool is_null() { return namep == nullptr; }
    virtual streamable* get();
    virtual void set(streamable* p);
    virtual streamed_pointer& operator=(const streamed_pointer& sp) { copy(sp); return *this; }

	// streamable implementations
	virtual streamable* make_instance() const { return new streamed_pointer(); };
	virtual void read_from_stream(istreamer& is);
	virtual void write_to_stream(ostreamer& os);


	static streamed_pointer* from_file(const string& file, bool is_file_disposable = true) ;
protected:
    static int get_pid() { if (pid == -1) { pid = (int)get_process_id(); } return pid; }
    static string generate_name(int id) { return tmp_dir + "sobject_" + get_pid() + "_" + id + ".dat"; }
    void copy(const streamed_pointer& sp);
    void dispose();
    void create_name(streamable* s);
	void create_from_name(const string& name);
};

#endif /* __STREAMED_POINTER_H__*/