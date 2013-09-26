/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 4; tab-width: 4 -*- */
#ifndef __HOPSTORAGE
#define __HOPSTORAGE

#include "hop.h"

#ifdef PROTOBUF_SUPPORT
#include "hop.pb.h"
#include <iostream>
#include <fstream>
#include <google/protobuf/stubs/common.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/gzip_stream.h>

#endif

void _hop_store_result_blob(const hop_result* res, const char* fname, bool compress);
void _hop_store_result_protobuf(const hop_result* res, const char* fname, bool compress);

#endif

