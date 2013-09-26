
#ifndef _ZLIB_COMPRESS_H_
#define _ZLIB_COMPRESS_H_


#include <string>
#include <iostream>
#if defined WIN32 | defined WIN64
#include <typeinfo.h>
#else
#include <typeinfo>
#endif

#include <zlib.h>

using namespace std;

// zlib (de)compression class
///////////////////////////////////////////////////////////////////////////////

class zlib_compress {
private:
	static const int CHUNK = 512*1024; // 512 KB (bigger chunk better compression)

	// buffer for compression
	char* buffer_in;
	char* buffer_out;
	// current position of empty data on buffer_in (starts with 0) if using for write_to_buf_and_compress
	// else current possition of data to be read in buffer_out if using read_to_buf_and_decompress
	int buffer_pos;
	// zlib compression level 
	int compression_level;
	z_stream zlib_stream;

	istream* in;
	ostream* out;

	// compresses chunk from buffer_in and writes it to this->out stream
	void compress_chunk(int flush_flag);
	// reades chunk data from this->in stream to buffer_in and 
	// returns decompressed in buffer_out
	void decompress_chunk(int flush_flag);

public:

	/**
	 * This method must be called before any (de)compression begins
	 * It will call inflateInit or deflateInit and prepare all necessary buffers.
	 * Same function is used for compression and decompression.
	 * Set argument decompress = true for decompression otherwise compression is used (default).
	 */
	void begin_compression(bool decompress = false);

	/**
	 * Writes data of 'size' to intermediate buffer (buffer_in) and
	 * when buffer is full it compresses it to buffer_out and writes whole buffer to underlying stream.
	 */
	void write_to_buf_and_compress(const char* ptr, streamsize size);
	/**
	 * Reads data from underlying stream to intermedeiate buffer (buffer_out) and
	 * decompresses its to buffer_in. Data of size 'size' is then copied from buffer_in to ptr.
	 */
	void read_to_buf_and_decompress(const char* ptr, streamsize size);
	/**
	 * Method must be called on end of (de)compression so it will deallocate memory 
	 * and make call to inflateEnd or deflateEnd. 
	 * Same function is used for compression and decompression.
	 * Set argument decompress = true for decompression otherwise compression is used (default).
	 */
	void end_compression(bool decompress = false);
public:
	zlib_compress(int zlib_compression_level = Z_NO_COMPRESSION);
	virtual ~zlib_compress();

	// MUST set input stream from where data will be read to intermediate buffer
	void set_istream(istream* pin) { this->in = pin; }
	// MUST set output stream from where data will be writen from compressed data buffer
	void set_ostream(ostream* pout) { this->out = pout; }

	// gets level of compression (uses default z_lib constants for compression level)
	int get_compression_level() { return this->compression_level; }
};

#endif
