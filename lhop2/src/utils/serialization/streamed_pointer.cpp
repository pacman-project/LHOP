/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
// streamed_pointer
///////////////////////////////////////////////////////////////////////////////

#include "streamed_pointer.h"
#include "utils/platform.h"

int streamed_pointer::pid = -1;
int streamed_pointer::nextid = 0;

string streamed_pointer::tmp_dir = TEMPORARY_DIRECTORY;

streamed_pointer::streamed_pointer(streamable* s)
{
    if (s == nullptr) namep = nullptr; 
    else create_name(s);
}

streamable* streamed_pointer::get()
{ 
    if (namep == nullptr) return nullptr; 
    else return streamable::read(namep->name); 
}

void streamed_pointer::set(streamable* p)
{
    if (namep != nullptr) p->save(namep->name, -1);
    else create_name(p);
}

streamed_pointer* streamed_pointer::from_file(const string& file, bool is_file_disposable /*= true*/)
{
	streamed_pointer* ptr = new streamed_pointer();
	ptr->create_from_name(file);
	if (ptr->is_null() == false)
		ptr->namep->disposable = is_file_disposable;
	return ptr;
}

void streamed_pointer::copy(const streamed_pointer& sp)
{
    if (sp.namep == namep) return;
    dispose(); 
    namep = sp.namep;
    if (namep != nullptr) ++namep->count;
}

void streamed_pointer::dispose()
{
    if (namep != nullptr && --namep->count <= 0) delete namep;
}

void streamed_pointer::create_name(streamable* p)
{
    this->create_from_name(generate_name(nextid));
    p->save(namep->name, -1);
}
void streamed_pointer::create_from_name(const string& name)
{
    namep = new counted_name(name);
    nextid++;
}
void streamed_pointer::read_from_stream(istreamer& is) {
	bool has_namep;
	is.read(has_namep);

	if (has_namep) {
		string name;
		is.read(name);

		namep = new counted_name(name);
	}
}
void streamed_pointer::write_to_stream(ostreamer& os) {
	os.write(namep != nullptr ? true : false);
	if (namep != nullptr)
		os.write(namep->name);

}