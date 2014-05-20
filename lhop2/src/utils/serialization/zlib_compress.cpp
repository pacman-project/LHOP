
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "zlib_compress.h"

zlib_compress::zlib_compress(int zlib_compression_level) :
buffer_in(nullptr), buffer_out(nullptr), buffer_pos(0), 
compression_level(zlib_compression_level), in(nullptr), out(nullptr)
{


}

zlib_compress::~zlib_compress() {
	// detele buffer is not empty
	if (buffer_in != nullptr)
		delete[] buffer_in;

	if (buffer_out != nullptr)
		delete[] buffer_out;
}

void zlib_compress::begin_compression(bool decompress) {
	// init structure (allocate deflate state)
	zlib_stream.zalloc = Z_NULL;
    zlib_stream.zfree = Z_NULL;
    zlib_stream.opaque = Z_NULL;

	// must start with avail_out = 0 since read_to_buf_and_decompress will use it to read first chunk from stream
	zlib_stream.avail_out = 0; 
	if (decompress == false) {
		deflateInit(&zlib_stream, compression_level);
	} else {
		inflateInit(&zlib_stream);

		// in case of inlation buffer_pos must be set to CHUNK, so new chunk will be read from stream
		buffer_pos = CHUNK;
		// also set zlib_stream.avil_in/out to 0, so that it will read from stream on first iteration
		zlib_stream.avail_in = 0;
		zlib_stream.avail_out = 0;
	}

	// also create buffer (and free it if already exists)
	if (buffer_in != nullptr)
		delete[] buffer_in;

	if (buffer_out != nullptr)
		delete[] buffer_out;

	buffer_in = new char[CHUNK];
	buffer_out = new char[CHUNK];
	
}

void zlib_compress::write_to_buf_and_compress(const char* ptr, streamsize size) {
	if (buffer_out == nullptr)
		// error
		return;
	// save data to buffer and write it to stream if already full

	// set start position and size of src data 
	int src_pos = 0;
	int src_size_left = size;

	do {
		// write as many chars as possible to buffer
		int n_buffer_avail_space = CHUNK - buffer_pos;

		// if no new space then compress this buffer and write it to underlying stream
		if (n_buffer_avail_space <= 0) {
			compress_chunk(Z_NO_FLUSH);

			// MUST recalculate available space since buffer_pos has changed
			n_buffer_avail_space = CHUNK - buffer_pos;
		}

		// now copy min(n_buffer_avail_space, src_size_left) data to buffer
		int data_copied = min(n_buffer_avail_space, src_size_left);
		memcpy(buffer_in + buffer_pos, ptr + src_pos, data_copied);

		// update position of buffer
		buffer_pos += data_copied;

		// update size and pos of src
		src_size_left -= data_copied;
		src_pos += data_copied;

		// if all data copied the finish else do another loop
	} while (src_size_left > 0 && src_pos > 0 && src_pos < size);
}

void zlib_compress::compress_chunk(int flush_flag) {
	// prepare input data
	zlib_stream.next_in = (Bytef*)buffer_in;
	zlib_stream.avail_in = buffer_pos;
	
	do {
		// prepare output buffer for compression
		zlib_stream.next_out = (Bytef*)buffer_out;
		zlib_stream.avail_out = CHUNK;

		// compress using DEFLATE
		int err = deflate(&zlib_stream, flush_flag);

		// check how much data has been returned
		int size_compressed_out = CHUNK - zlib_stream.avail_out;
		// write returned data to underlying stream
		out->write(buffer_out, size_compressed_out);
		
		// if still some data left in then do another loop
	} while (zlib_stream.avail_out == 0);

	// also prepare new buffer (just reset buffer_pos to 0)
	buffer_pos = 0;
}

void zlib_compress::read_to_buf_and_decompress(const char* ptr, streamsize size) {

	// read from stream and decompress until size of decompressed data is "size"

	int dst_pos = 0;
	int dst_size_left = size;

	do {
		int n_buffer_avail_data = CHUNK - buffer_pos;

		// check if buffer_out is full
		if (n_buffer_avail_data <= 0) { 
			// get new chunk from stream and compress it
			decompress_chunk(Z_NO_FLUSH);

			// MUST recalculate available space since buffer_pos has changed
			n_buffer_avail_data = CHUNK - buffer_pos;
		}

		// read max(size or CHUNK) of data and copy it to ptr
		int data_copied = min(n_buffer_avail_data, dst_size_left);
		memcpy((void*)(ptr + dst_pos), buffer_out + buffer_pos, data_copied);

		// update position of buffer
		buffer_pos += data_copied;

		// update size and pos of dst
		dst_size_left -= data_copied;
		dst_pos += data_copied;

		// if not all data read, then do another loop
	} while (dst_size_left > 0);
	
}

void zlib_compress::decompress_chunk(int flush_flag) {
	int ret = 0;

	do {
		// if buffer_in is already fully read then read new chunk from stream
		if (zlib_stream.avail_in == 0) {

			// prepare new input data
			zlib_stream.next_in = (Bytef*)buffer_in;
			zlib_stream.avail_in = CHUNK;

			// then read new chunk from underlying stream
			in->read(buffer_in, CHUNK);
		}

		// update next_out to buffer_out only if all buffer already used
		// (i.e. zlib_stream.next_out is 0 when buffer full)
		if (zlib_stream.avail_out == 0) {
			zlib_stream.next_out = (Bytef*)buffer_out;
			zlib_stream.avail_out = CHUNK;
		}

		// inflate
		ret = inflate(&zlib_stream, flush_flag);

		// repeat loop if buffer_out is not yet full or is end of stream
	} while (zlib_stream.avail_out != 0 && ret != Z_STREAM_END);


	buffer_pos = 0;
}
void zlib_compress::end_compression(bool decompress) {

	if (decompress == false) {
		// MUST deflate buffer if any data left

		// if buffer is empty then should still call deflate with Z_FINISH
		// in that case just write one byte of dummy data (maybe that is not neccessary)
		if (buffer_pos <= 0 || buffer_pos > CHUNK) {
			// just write one byte (actualy just set buffer_pos to 1, so only one byte will be sent)
			// TEST IT IF IT IS NEEDED !!!
			buffer_pos = 1;
		}

		// compress any chunk left and tell z_lib to flush it
		compress_chunk(Z_FINISH);

		deflateEnd(&zlib_stream);
	} else {
		// no need to inflate any more input since if it was not requested
		// by call to read_to_buf_and_decompress then it is not needed
		inflateEnd(&zlib_stream);
	}


	// finally free buffer if present (try to delete it even if no compression just in case)
	if (buffer_in != nullptr) {
		delete[] buffer_in;
		buffer_in = nullptr;
	}

	if (buffer_out != nullptr) {
		delete[] buffer_out;
		buffer_out = nullptr;
	}

}
