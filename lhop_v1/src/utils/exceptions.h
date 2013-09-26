
#pragma once
#ifndef _EXCEPTIONS_
#define _EXCEPTIONS_

#include <exception>
#include <string>
#include <sstream>
#include <typeinfo>
#ifdef __linux__
#include <execinfo.h>
#endif

#define new_libhop_exception(msg) libhop_exception(__FILE__, __LINE__, string(msg))
#define custom_libhop_exception(exception_name, msg) exception_name(__FILE__, __LINE__, string(msg))

class libhop_exception : public std::exception
{
private:
	string msg;
	string file;
	size_t line;

	string what_msg;
public:

	libhop_exception(const string& m_file, const size_t m_line, const string& m_msg = string(), const string& ex_name = "libhop_exception"): msg(m_msg), file(m_file), line(m_line) { create_what_msg(ex_name); }
	libhop_exception(const string& m_file, const size_t m_line, const std::exception& from_ex, const string& ex_name = "libhop_exception"): file(m_file), line(m_line) {
		msg = "wrapped exception class '" + string(typeid(from_ex).name()) + "': " + string(from_ex.what());
		create_what_msg(ex_name);
	}
	virtual ~libhop_exception() throw() {};

	const char* what() const throw() { return what_msg.c_str(); }
	const char* get_file() const throw() { return file.c_str(); }
	size_t get_line() const throw() { return line; }

	virtual const char* get_message() {
		return msg.c_str();
	}

private:
	void create_what_msg(const string& ex_name){
		stringstream ss;
		ss << "Exception '" << ex_name << "' thrown in file '" << file << "' line '" << line << "' with message:\n" << get_message();

#ifdef __linux__
		void *buffer[100];
		char **strings;

		int nptrs = backtrace(buffer, 100);		

		strings = backtrace_symbols(buffer, nptrs);
		if (strings != nullptr) {
			ss << endl << endl << "Stack trace:" << endl;
			for (int j = 0; j < nptrs; j++)
				ss << strings[j] << endl;
			free(strings);
		}
#endif
		what_msg = ss.str();
	}
};

// DO NOT call explicitly use custom_libhop_exception macro (e.g. throw custom_libhop_exception(opencl_exception,"msg"))
class opencl_exception : public libhop_exception
{
public:
	opencl_exception(const string& file, const size_t line, const string& msg = string()): libhop_exception(file, line, msg, "opencl_exception") {}
	opencl_exception(const string& file, const size_t line, const std::exception& from_ex): libhop_exception(file, line, from_ex, "opencl_exception") {}
};

// DO NOT call explicitly use custom_libhop_exception macro (e.g. throw custom_libhop_exception(io_exception,"msg"))
class io_exception : public libhop_exception
{
public:
	io_exception(const string& file, const size_t line, const string& msg = string()): libhop_exception(file, line, msg, "io_exception") {}
	io_exception(const string& file, const size_t line, const std::exception& from_ex): libhop_exception(file, line, from_ex, "io_exception") {}
};

#endif


