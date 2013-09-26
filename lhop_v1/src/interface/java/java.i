%pragma(java) jniclassclassmodifiers="class"
%pragma(java) moduleclassmodifiers="public class"
%typemap("javapackage") SWIGTYPE & "com.vicos.hop"
%feature("autodoc", "1");
%include "javahead.swg"


%module(docstring="Java bindings for the libhop library",directors="1") JHOP

%feature("director") hop_deployer_mapreudce;
%feature("director") hop_learning_controler; 

%include "carrays.i"
%array_class(unsigned char, HopByteArray);

%{

#include <stdio.h>
#include <vector>
#include "hop.h"

using namespace std;

%}

%typemap(out) const char* {
	if ($1) { $result = jenv->NewStringUTF((const char *)$1); delete $1; }
}

/* void_data_bytes is byte buffer that should get mapped to java byte[] */ 
struct void_data_bytes;

%{
struct void_data_bytes {
	void* bytes;
	int size;
};

/* special type that is used only for DataBufferByte */ 
typedef void_data_bytes java_data_buffer_byte;
typedef void_data_bytes java_data_buffer_int;
%}

/* add methods for obtaining/creating java byte[] arrays from/to hop_blob */
%extend hop_blob {
	void_data_bytes getBytes(int offset = 0, int length = 0) {
		void_data_bytes result;		
		if (offset >= 0 && length > 0) {
			if (offset < $self->get_size()) {
				result.bytes = $self->data() + offset;
				result.size = min($self->get_size() - offset, length);
			}
		} else {
			result.bytes = $self->data();
			result.size = $self->get_size();
		}
		return result;
	}
	static hop_blob fromBytes(void_data_bytes data) {
		return hop_blob(data.bytes, data.size);
	}
};


%typemap(in) void_data_bytes {
	// $input is jByteArray
	// $1.bytes is void* (i.e. pointer to data)
	// $1.size is int (i.e. size of data)	
	
	$1.size = JCALL1(GetArrayLength, jenv, $input);
	
	if ($1.size > 0) {
		$1.bytes = malloc($1.size);
		
		void* arr = JCALL2(GetPrimitiveArrayCritical, jenv, $input, 0);
		memcpy($1.bytes, arr, $1.size);
		JCALL3(ReleasePrimitiveArrayCritical, jenv, $input, arr, 0);
	} else {
		$1.bytes = NULL;
	}
}

%typemap(out) void_data_bytes {
	// $result is jByteArray
	// $1.bytes is void* (i.e. pointer to data)
	// $1.size is int (i.e. size of data)
	//printf("creating new array\n");
	$result = JCALL1(NewByteArray, jenv, $1.size);
	
	//printf("size from byteArray %d\n", JCALL1(GetArrayLength, jenv, $result));
	
	//printf("if size > 0\n");
	if ($1.size > 0) {
		//printf("size is > 0, gettng primitive array\n");
		void* arr = JCALL2(GetPrimitiveArrayCritical, jenv, $result, 0);
		//printf("coping from %p to %p of size %d\n", $1.bytes, arr, $1.size);
		memcpy(arr, $1.bytes, $1.size);
		//printf("releasing\n");
		JCALL3(ReleasePrimitiveArrayCritical, jenv, $result, arr, 0);
	}
	//printf("end\n");
}

%typemap(jni) void_data_bytes "jbyteArray"
%typemap(jtype) void_data_bytes "byte[]"
%typemap(jstype) void_data_bytes "byte[]" 

%typemap(javain) void_data_bytes "$javainput"
%typemap(javaout) void_data_bytes {
	return $jnicall;
}

/* java_data_buffer_byte represents void* constructed from DataBufferByte  */
%typemap(in) java_data_buffer_byte {
	// $input is jObject of type DataBufferByte
	// $1 is void_data_bytes
		// $1.size is size of buffer
		// $1.bytes is buffer
		
	// verify input class instance type (should be DataBufferByte)
	jclass correct_class = JCALL1(FindClass, jenv, "java/awt/image/DataBufferByte");
	jboolean is_correct_instance = JCALL2(IsInstanceOf, jenv, $input, correct_class);
	
	if (is_correct_instance == true) {
		// prepare methods
		jmethodID method_getSize = JCALL3(GetMethodID, jenv, correct_class, "getSize","()I");
		jmethodID method_getNumBanks = JCALL3(GetMethodID, jenv, correct_class, "getNumBanks","()I");
		jmethodID method_getData = JCALL3(GetMethodID, jenv, correct_class, "getData","(I)[B");
		
		if (method_getSize == 0 || method_getNumBanks == 0 || method_getData == 0) {
			jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
		    JCALL2(ThrowNew, jenv, clazz, "No such method for class DataBufferByte: getSize, getNumBanks, getData");
		    return $null;
		}
		
		jint num_banks = JCALL2(CallIntMethod, jenv, $input, method_getNumBanks);
		$1.size = JCALL2(CallIntMethod, jenv, $input, method_getSize);
		
		$1.bytes = malloc($1.size); 
		
		int offset = 0;
		for (jint i = 0; i < num_banks; i++) {			
			jbyteArray byte_array = (jbyteArray)JCALL3(CallObjectMethod, jenv, $input, method_getData, i);
			int byte_array_size = JCALL1(GetArrayLength, jenv, byte_array);
			
			void* arr = JCALL2(GetPrimitiveArrayCritical, jenv, byte_array, 0);
			memcpy($1.bytes + offset, arr, byte_array_size);
			JCALL3(ReleasePrimitiveArrayCritical, jenv, byte_array, arr, 0);	
			
			offset += byte_array_size;
		}
	}
}

%typemap(jni) java_data_buffer_byte "jobject"
%typemap(jtype) java_data_buffer_byte "java.awt.image.DataBufferByte"
%typemap(jstype) java_data_buffer_byte "java.awt.image.DataBufferByte" 

%typemap(javain) java_data_buffer_byte "$javainput"
%typemap(javaout) java_data_buffer_byte {
	return $jnicall;
}

/* java_data_buffer_byte represents void* constructed from DataBufferInt  */
%typemap(in) java_data_buffer_int {
	// $input is jObject of type DataBufferInt
	// $1 is void_data_bytes
		// $1.size is size of buffer
		// $1.bytes is buffer
	
	// verify input class instance type (should be DataBufferInt)
	jclass correct_class = JCALL1(FindClass, jenv, "java/awt/image/DataBufferInt");
	jboolean is_correct_instance = JCALL2(IsInstanceOf, jenv, $input, correct_class);
	//printf("checking instance\n");
	if (is_correct_instance == true) {
		// prepare methods
		jmethodID method_getSize = JCALL3(GetMethodID, jenv, correct_class, "getSize","()I");
		jmethodID method_getNumBanks = JCALL3(GetMethodID, jenv, correct_class, "getNumBanks","()I");
		jmethodID method_getData = JCALL3(GetMethodID, jenv, correct_class, "getData","(I)[I");
		
		if (method_getSize == 0 || method_getNumBanks == 0 || method_getData == 0) {
			jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
		    JCALL2(ThrowNew, jenv, clazz, "No such method for class DataBufferInt: getSize, getNumBanks, getData");
		    return $null;
		}
		//printf("methods checked\n");
		jint num_banks = JCALL2(CallIntMethod, jenv, $input, method_getNumBanks);
		$1.size = JCALL2(CallIntMethod, jenv, $input, method_getSize) * 4;
		
		$1.bytes = malloc($1.size); 
		//printf("allocated bytes %d\n",$1.size);
		int offset = 0;
		for (jint i = 0; i < num_banks; i++) {
			//printf("allocating int array ... \n");
			jintArray int_array = (jintArray)JCALL3(CallObjectMethod, jenv, $input, method_getData, i);
			//printf("getting int array size \n");
			int byte_array_size = JCALL1(GetArrayLength, jenv, int_array) * 4;
			//printf("allocated int array of size %d \n", byte_array_size);
			
			void* arr = JCALL2(GetPrimitiveArrayCritical, jenv, int_array, 0);
			memcpy($1.bytes + offset, arr, byte_array_size);
			JCALL3(ReleasePrimitiveArrayCritical, jenv, int_array, arr, 0);	
			
			offset += byte_array_size;
		}
	}
	
}

%typemap(jni) java_data_buffer_int "jobject"
%typemap(jtype) java_data_buffer_int "java.awt.image.DataBufferInt"
%typemap(jstype) java_data_buffer_int "java.awt.image.DataBufferInt" 

%typemap(javain) java_data_buffer_int "$javainput"
%typemap(javaout) java_data_buffer_int {
	return $jnicall;
}

/* Expose method only made from bytes[] or DataBufferByte */ 

%ignore from_bytes_rgb32(int width, int height, void* bytes);

%extend hop_image {
	static hop_image from_bytes_rgb32(int width, int height, void_data_bytes data) {
		hop_image im = hop_image::from_bytes_rgb32(width, height, data.bytes);
		free(data.bytes);
		return im;
	}
	static hop_image from_bytes_rgb32(int width, int height, java_data_buffer_byte data) {
		hop_image im = hop_image::from_bytes_rgb32(width, height, data.bytes);
		free(data.bytes);
		return im;
	}
	static hop_image from_bytes_rgb32(int width, int height, java_data_buffer_int data) {
		hop_image im = hop_image::from_bytes_rgb32(width, height, data.bytes);
		free(data.bytes);
		return im;
	}
}

/* Wrapper for array of hop_results */

%typemap(out) hop_result_wrapper {
  jlong *arr;
  int i = 0;
  if (!$1.size) {
    return $null;
  }
  $result = JCALL1(NewLongArray, jenv, $1.size);
  if (!$result) {
    return $null;
  }
  arr = JCALL2(GetLongArrayElements, jenv, $result, 0);
  if (!arr) {
    return $null;
  }
  for (i=0; i<$1.size; i++) {
    arr[i] = 0; 
    hop_result t1 = $1.results[i];
    arr[i] = (jlong) new hop_result(t1);
  }
  JCALL3(ReleaseLongArrayElements, jenv, $result, arr, 0);
  delete [] $1.results;
}

/* These 3 typemaps tell SWIG what JNI and Java types to use */
%typemap(jni) hop_result_wrapper "jlongArray"
%typemap(jtype) hop_result_wrapper "long[]"
%typemap(jstype) hop_result_wrapper "hop_result[]"

%typemap(javaout) hop_result_wrapper {
    return $javaclassname.cArrayWrap($jnicall, true);
  }

%typemap(javacode) hop_result_wrapper %{
  protected static hop_result[] cArrayWrap(long[] cArray, boolean cMemoryOwn) {
    if (cArray == null) return null;
    hop_result[] arrayWrapper = new hop_result[cArray.length];
    for (int i=0; i<cArray.length; i++)
      arrayWrapper[i] = cArray[i] == 0 ? null : new hop_result(cArray[i], cMemoryOwn);
    return arrayWrapper;
  }
%}

%ignore hop_inference(hop_result*&, hop_image& , const char* );
%ignore hop_inference(hop_result*&, hop_image& , const std::list<irectangle2>, const char* );

%typemap(in) std::list<irectangle2> {
	// $input is jobject of type java.util.LinkedList<java.awt.Rectangle>
	// $1 is std::list<irectangle2>
	
	jclass linkedlist_class = JCALL1(FindClass, jenv, "java/util/LinkedList");
	
	jboolean is_correct_instance = JCALL2(IsInstanceOf, jenv, $input, linkedlist_class);
	//printf("is correct instance\n");
	if (is_correct_instance == true) {
		jclass iterator_class = JCALL1(FindClass, jenv, "java/util/Iterator");
		jclass rectangle_class = JCALL1(FindClass, jenv, "java/awt/Rectangle");
		if (linkedlist_class  == 0) printf("linkedlist_class is null\n");
		if (iterator_class  == 0) printf("iterator_class is null\n");
		if (rectangle_class  == 0) printf("rectangle_class is null\n");
		//printf("checking methods now\n");
		// prepare methods
		jmethodID method_iterator = JCALL3(GetMethodID, jenv, linkedlist_class, "iterator","()Ljava/util/Iterator;");

		jmethodID method_hasnext = JCALL3(GetMethodID, jenv, iterator_class, "hasNext","()Z");
		jmethodID method_next = JCALL3(GetMethodID, jenv, iterator_class, "next","()Ljava/lang/Object;");

		jmethodID method_getX = JCALL3(GetMethodID, jenv, rectangle_class, "getX","()D");
		jmethodID method_getY = JCALL3(GetMethodID, jenv, rectangle_class, "getY","()D");
		jmethodID method_getW = JCALL3(GetMethodID, jenv, rectangle_class, "getWidth","()D");
		jmethodID method_getH = JCALL3(GetMethodID, jenv, rectangle_class, "getHeight","()D");

		if (method_getX == 0 || method_getY == 0 || method_getW == 0 || method_getH == 0) {
			jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
			JCALL2(ThrowNew, jenv, clazz, "No such method for class java.awt.Rectangle: getX, getY, getWidth, getHeight");
			return $null;
		}
		if (method_hasnext == 0 || method_next == 0) {
			jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
			JCALL2(ThrowNew, jenv, clazz, "No such method for class java.util.Iterator: hasNext, next");
			return $null;
		}
		if (method_iterator == 0) {
			jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
			JCALL2(ThrowNew, jenv, clazz, "No such method for class java.util.LinkedList: iterator");
			return $null;
		}
		
		//printf("methods checked\n");
		jobject iterator = JCALL2(CallObjectMethod, jenv, $input, method_iterator);
		if (iterator == NULL) {
			jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
			JCALL2(ThrowNew, jenv, clazz, "Got null iterator from java/util/LinkedList object");
			return $null;
		}

		while (JCALL2(CallBooleanMethod, jenv, iterator, method_hasnext) == true) {
			jobject rect = JCALL2(CallObjectMethod, jenv, iterator, method_next);
			if (JCALL2(IsInstanceOf, jenv, rect, rectangle_class) == false) {
				jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
				JCALL2(ThrowNew, jenv, clazz, "Found invalid object type in LinkedList. Should be java.awt.Rectangle, but is not !!");
				return $null;
			}
			int x = (int)JCALL2(CallDoubleMethod, jenv, rect, method_getX);
			int y = (int)JCALL2(CallDoubleMethod, jenv, rect, method_getY);
			int w = (int)JCALL2(CallDoubleMethod, jenv, rect, method_getW);
			int h = (int)JCALL2(CallDoubleMethod, jenv, rect, method_getH);
			$1.push_back(irectangle2(x,y,x+w,y+h));
		}
	}
}

%typemap(out) std::list<irectangle2> {
	// $result is jobject of type java.util.LinkedList<java.awt.Rectangle>
	// $1 is std::list<irectangle2>
	jclass linkedlist_class = JCALL1(FindClass, jenv, "java/util/LinkedList");
	jclass rectangle_class = JCALL1(FindClass, jenv, "java/awt/Rectangle");

	// prepare methods
	jmethodID method_linkedlist_construct = JCALL3(GetMethodID, jenv, linkedlist_class, "<init>","()V");
	jmethodID method_add = JCALL3(GetMethodID, jenv, linkedlist_class, "add","(Ljava/lang/Object;)Z");
	jmethodID method_rectangle_construct = JCALL3(GetMethodID, jenv, rectangle_class, "<init>","(IIII)V");

	if (method_linkedlist_construct == 0) {
		jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
		JCALL2(ThrowNew, jenv, clazz, "No such method for class java.util.LinkedList: constructor");
		return $null;
	}
	if (method_rectangle_construct == 0) {
		jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
		JCALL2(ThrowNew, jenv, clazz, "No such method for class java.awt.Rectangle: constructor");
		return $null;
	}	

	$result = JCALL2(NewObject, jenv, linkedlist_class, method_linkedlist_construct);
	//printf("methods checked\n");
	for (std::list<irectangle2>::const_iterator it = $1.begin(); it != $1.end(); it++) {
		jobject rect = JCALL6(NewObject, jenv, rectangle_class, method_rectangle_construct, it->ll.x, it->ll.y, it->ur.x - it->ll.x, it->ur.y - it->ll.y);
		//printf("rect constructed\n");
		JCALL3(CallBooleanMethod, jenv, $result, method_add, rect);
		//printf("rect added\n");
	}
}

%typemap(jni) std::list<irectangle2> "jobject"
%typemap(jtype) std::list<irectangle2> "java.util.LinkedList<java.awt.Rectangle>"
%typemap(jstype) std::list<irectangle2> "java.util.LinkedList<java.awt.Rectangle>"
%typemap(javain) std::list<irectangle2> "$javainput"
%typemap(javaout) std::list<irectangle2> {
	return $jnicall;
}


%typemap(javacode) hop_histogram_descriptor %{
  // x,y,w,h,hist added manualy from java.i
  public float x,y,w,h;
  public float[] hist;
  public hop_histogram_descriptor(float x, float y, float w, float h, float[] hist) {
    this.x = x; this.y = y;
    this.w = w; this.h = h;
    this.hist = hist;
  }
%}

%typemap(out) hop_histogram_descriptor_wraper {
	// $result is jobjectArray (hop_histogram_descriptor[])
	// $1 is hop_histogram_descriptor_wraper

	jclass descriptor_class = JCALL1(FindClass, jenv, "com/vicos/hop/hop_histogram_descriptor");
	 
	if (descriptor_class == 0) {
		jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
		JCALL2(ThrowNew, jenv, clazz, "No such class found: missing class com/vicos/hop/hop_histogram_descriptor");
		return $null;
	}

	jmethodID descriptor_construct = JCALL3(GetMethodID, jenv, descriptor_class, "<init>","(FFFF[F)V");

	if (descriptor_construct == 0) {
		jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
		JCALL2(ThrowNew, jenv, clazz, "No such method for class com.vicos.hop.hop_histogram_descriptor: constructor");
		return $null;
	}

	//printf("got class and constructor\n");
	//printf("creting array of size %d\n",$1.size);
	
	$result = JCALL3(NewObjectArray,jenv, $1.size, descriptor_class, 0);
	//printf("created new array of size %d\n", $1.size);
	for (int i = 0; i < $1.size; i++) {
		//printf("converting index %d\n", i);
		int histSize = (*($1.results))[i].hist.size() ;
		float* hist = &((*($1.results))[i].hist)[0];
		jfloatArray jhist = JCALL1(NewFloatArray, jenv, histSize);
		
		//printf("float array has size %d\n", histSize);
		void* arr = JCALL2(GetPrimitiveArrayCritical, jenv, jhist, 0);
		//printf("coping from %p to %p of size %d\n", hist, arr, histSize * sizeof(float));
		memcpy(arr, hist, histSize * sizeof(float));
		//printf("releasing\n");
		JCALL3(ReleasePrimitiveArrayCritical, jenv, jhist, arr, 0);
		
		jobject element = JCALL7(NewObject, jenv, descriptor_class, descriptor_construct, (*($1.results))[i].x,(*($1.results))[i].y,(*($1.results))[i].w,(*($1.results))[i].h, jhist);
		JCALL3(SetObjectArrayElement, jenv, $result, i, element);

		JCALL1(DeleteLocalRef, jenv, element);
		JCALL1(DeleteLocalRef, jenv, jhist);
		
	}

	delete $1.results;
	//printf("done\n");
}

%typemap(jni) hop_histogram_descriptor_wraper "jobjectArray"
%typemap(jtype) hop_histogram_descriptor_wraper "hop_histogram_descriptor[]"
%typemap(jstype) hop_histogram_descriptor_wraper "hop_histogram_descriptor[]"
%typemap(javain) hop_histogram_descriptor_wraper "$javainput"
%typemap(javaout) hop_histogram_descriptor_wraper {
    return $jnicall;
}

%typemap(in) const int* bin {
    	// $input is jintArray
	// $1 is int*
	if ($input != 0) {
		int size = JCALL1(GetArrayLength, jenv, $input);
		$1 = (int*)malloc(size * sizeof(int));
		
		void* arr = JCALL2(GetPrimitiveArrayCritical, jenv, $input, 0);
		//printf("coping from %p to %p of size %d\n", arr, $1, size * sizeof(int));
		memcpy($1, arr, size * sizeof(int));
		//printf("releasing\n");
		JCALL3(ReleasePrimitiveArrayCritical, jenv, $input, arr, 0);
	}
}

%typemap(freearg) const int* bin {
	if ($1 != NULL)
		free($1);
}

%typemap(jni) const int* bin "jintArray"
%typemap(jtype) const int* bin "int[]"
%typemap(jstype) const int* bin "int[]"
%typemap(javain) const int* bin "$javainput"
%typemap(javaout) const int* bin {
    return $jnicall;
}

%typemap(in) const float* centers {
    	// $input is jfloatArray
	// $1 is float*
	if ($input != NULL) {
		int size = JCALL1(GetArrayLength, jenv, $input);
		$1 = (float*)malloc(size * sizeof(float));
		
		void* arr = JCALL2(GetPrimitiveArrayCritical, jenv, $input, 0);
		//printf("coping from %p to %p of size %d\n", arr, $1, size * sizeof(int));
		memcpy($1, arr, size * sizeof(float));
		//printf("releasing\n");
		JCALL3(ReleasePrimitiveArrayCritical, jenv, $input, arr, 0);
	}
}


%typemap(freearg) const float* centers {
	if ($1 != NULL)
		free($1);
}

%typemap(jni) const float* centers "jfloatArray"
%typemap(jtype) const float* centers "float[]"
%typemap(jstype) const float* centers "float[]"
%typemap(javain) const float* centers "$javainput"
%typemap(javaout) const float* centers {
    return $jnicall;
}

%typemap(in) (const int* layers, const int layers_count) {
    	// $input is jintArray
	// $1 is int*
	// $2 is int
	$2 = JCALL1(GetArrayLength, jenv, $input);
	$1 = (int*)malloc($2 * sizeof(int));
	
	void* arr = JCALL2(GetPrimitiveArrayCritical, jenv, $input, 0);
	//printf("coping from %p to %p of size %d\n", arr, $1, $2 * sizeof(int));
	memcpy($1, arr, $2 * sizeof(int));
	//printf("releasing\n");
	JCALL3(ReleasePrimitiveArrayCritical, jenv, $input, arr, 0);
}

%typemap(freearg) (const int* layers, const int layers_count) {
      free($1);
}

%typemap(jni) (const int* layers, const int layers_count) "jintArray"
%typemap(jtype) (const int* layers, const int layers_count) "int[]"
%typemap(jstype) (const int* layers, const int layers_count) "int[]"
%typemap(javain) (const int* layers, const int layers_count) "$javainput"
%typemap(javaout) (const int* layers, const int layers_count) {
    return $jnicall;
}

################################################


%typemap(javacode) hop_detection %{
  // x,y,w,h,responses added manualy from java.i
  public float x,y,w,h;
  public int category_id;
  public float[] responses;
  public hop_detection(float x, float y, float w, float h, int category_id, float[] responses) {
    this.x = x; this.y = y;
    this.w = w; this.h = h;
    this.category_id = category_id;
    this.responses = responses;
  }
%}

%typemap(out) hop_detections_wraper {
	// $result is jobject of type java.util.LinkedList<com.vicos.hop.hop_detection>
	// $1 is hop_detections_wraper	
	jclass linkedlist_class = JCALL1(FindClass, jenv, "java/util/LinkedList");
	jclass rectangle_class = JCALL1(FindClass, jenv, "java/awt/Rectangle");

	jclass detection_class = JCALL1(FindClass, jenv, "com/vicos/hop/hop_detection");
	 
	if (detection_class == 0) {
		jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
		JCALL2(ThrowNew, jenv, clazz, "No such class found: missing class com/vicos/hop/hop_detection");
		return $null;
	}
	
	// prepare methods
	jmethodID method_linkedlist_construct = JCALL3(GetMethodID, jenv, linkedlist_class, "<init>","()V");
	jmethodID method_add = JCALL3(GetMethodID, jenv, linkedlist_class, "add","(Ljava/lang/Object;)Z");
	jmethodID method_rectangle_construct = JCALL3(GetMethodID, jenv, rectangle_class, "<init>","(IIII)V");
	jmethodID method_detection_construct = JCALL3(GetMethodID, jenv, detection_class, "<init>","(FFFFI[F)V");
	
	if (method_linkedlist_construct == 0) {
		jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
		JCALL2(ThrowNew, jenv, clazz, "No such method for class java.util.LinkedList: constructor");
		return $null;
	}

	if (method_detection_construct == 0) {
		jclass clazz = JCALL1(FindClass, jenv, "java/lang/Exception");
		JCALL2(ThrowNew, jenv, clazz, "No such method for class com.vicos.hop.hop_detection: constructor");
		return $null;
	}

	$result = JCALL2(NewObject, jenv, linkedlist_class, method_linkedlist_construct);
	//printf("methods checked\n"); fflush(stdout);
	if ($1.bbox_list != NULL) {
		for (std::list<hop_detection>::const_iterator it = $1.bbox_list->begin(); it != $1.bbox_list->end(); it++) {
			int responsesSize = (*it).responses.size() ;
			const float* responses = &((*it).responses)[0];
			jfloatArray jresponses = JCALL1(NewFloatArray, jenv, responsesSize);
			
			//printf("float array has size %d\n", responsesSize); fflush(stdout);
			void* arr = JCALL2(GetPrimitiveArrayCritical, jenv, jresponses, 0);
			//printf("coping from %p to %p of size %d\n", responses, arr, responsesSize * sizeof(float));
			memcpy(arr, responses, responsesSize * sizeof(float)); fflush(stdout);
			//printf("releasing\n"); fflush(stdout);
			JCALL3(ReleasePrimitiveArrayCritical, jenv, jresponses, arr, 0);
			
			jobject detection = jenv->NewObject(detection_class, method_detection_construct, (float)(*it).box.ll.x, (float)(*it).box.ll.y, (float)(*it).box.ur.x - (*it).box.ll.x, (float)(*it).box.ur.y - (*it).box.ll.y, (*it).category_id, jresponses);
			//printf("detection constructed\n"); fflush(stdout);
			JCALL3(CallBooleanMethod, jenv, $result, method_add, detection);
			//printf("detection added\n"); fflush(stdout);
			
			JCALL1(DeleteLocalRef, jenv, jresponses);
		}
		
		delete $1.bbox_list;
	}
	
}

%typemap(jni) hop_detections_wraper "jobject"
%typemap(jtype) hop_detections_wraper "java.util.LinkedList<hop_detection>"
%typemap(jstype) hop_detections_wraper "java.util.LinkedList<hop_detection>"
%typemap(javain) hop_detections_wraper "$javainput"
%typemap(javaout) hop_detections_wraper {
    return $jnicall;
}

################################################

%typemap(in) std::vector<hop_result>& {
	// $input is jlongArray
	// $1 is std::vector<hop_result>
	int size = JCALL1(GetArrayLength, jenv, $input);
	long* arr = (long*)JCALL2(GetPrimitiveArrayCritical, jenv, $input, 0);
	$1 = new std::vector<hop_result>(size);
	for (int i = 0; i < size; i++)
	  (*$1)[i] = hop_result(*(hop_result*)arr[i]);
	JCALL3(ReleasePrimitiveArrayCritical, jenv, $input, arr, 0);
	
}
%typemap(freearg) std::vector<hop_result>& {
      delete $1;
}

%typemap(jni) std::vector<hop_result>& "jlongArray"
%typemap(jtype) std::vector<hop_result>& "long[]"
%typemap(jstype) std::vector<hop_result>& "hop_result[]"
%typemap(javain) std::vector<hop_result>& "$javaclassname.cArrayWrap($javainput)"

%typemap(javacode) std::vector<hop_result>& %{
protected static long[] cArrayWrap(hop_result[] list) {
  long[] array = new long[list.length];
  for (int i = 0; i < list.length; i++) {
    array[i] = hop_result.getCPtr(list[i]);
  }
  return array;
}

%}

##########################################################

%typemap(in) std::vector<hop_streamable>& {
	// $input is jlongArray
	// $1 is std::vector<hop_streamable>
	int size = JCALL1(GetArrayLength, jenv, $input);
	long* arr = (long*)JCALL2(GetPrimitiveArrayCritical, jenv, $input, 0);
	$1 = new std::vector<hop_streamable>(size);
	for (int i = 0; i < size; i++)
	  (*$1)[i] = hop_streamable(*(hop_streamable*)arr[i]);
	JCALL3(ReleasePrimitiveArrayCritical, jenv, $input, arr, 0);
	
}

%typemap(freearg) std::vector<hop_streamable> {
      delete $1;
}

%typemap(jni) std::vector<hop_streamable>& "jlongArray"
%typemap(jtype) std::vector<hop_streamable>& "long[]"
%typemap(jstype) std::vector<hop_streamable>& "hop_streamable[]"
%typemap(javain) std::vector<hop_streamable>& "$javaclassname.cArrayWrap($javainput)"

%typemap(javacode) std::vector<hop_streamable>& %{
protected static long[] cArrayWrap(hop_streamable[] list) {
  long[] array = new long[list.length];
  for (int i = 0; i < list.length; i++) {
    array[i] = hop_streamable.getCPtr(list[i]);
  }
  return array;
}

%}
##########################################################

%typemap(out) std::vector<hop_streamable> {
	// $1 is std::vector<hop_streamable>
	// $result is jlongArray
	jlong *arr;
	int i = 0;
	if (!$1.size()) {
	  return $null;
	}
	$result = JCALL1(NewLongArray, jenv, $1.size());
	if (!$result) {
	  return $null;
	}
	arr = JCALL2(GetLongArrayElements, jenv, $result, 0);
	if (!arr) {
	  return $null;
	}
	for (i=0; i<$1.size(); i++) {
	  arr[i] = 0; 
	  hop_streamable t1 = $1.at(i);
	  arr[i] = (jlong) new hop_streamable(t1);
	}
	JCALL3(ReleaseLongArrayElements, jenv, $result, arr, 0);
}

%typemap(jni) std::vector<hop_streamable> "jlongArray"
%typemap(jtype) std::vector<hop_streamable> "long[]"
%typemap(jstype) std::vector<hop_streamable> "hop_streamable[]"
%typemap(javaout) std::vector<hop_streamable> {
	return $javaclassname.cArrayWrap($jnicall, true);
}

%typemap(javacode) std::vector<hop_streamable> %{
protected static hop_streamable[] cArrayWrap(long[] cArray, boolean cMemoryOwn) {
    if (cArray == null) return null;
    hop_streamable[] arrayWrapper = new hop_streamable[cArray.length];
    for (int i=0; i<cArray.length; i++)
      arrayWrapper[i] = cArray[i] == 0 ? null : new hop_streamable(cArray[i], cMemoryOwn);
    return arrayWrapper;
  }
%}

%typemap(javacode) hop_streamable %{  
  public long getptr() {
		return swigCPtr;
  }
%}


%ignore hop_create_histograms(std::vector<laydisplay::histogram_descriptor>*&, hop_result&, std::list<irectangle2>, hop_histogram_generator&);
%ignore hop_get_detections(std::list<hop_detection>*&, hop_result&, hop_library&, const int);

%include "hop.h"


hop_result_wrapper hop_inference(hop_image& image, const char* params);
hop_result_wrapper hop_inference(hop_image& image, const std::list<irectangle2> mask_regions, const char* params);

hop_histogram_descriptor_wraper hop_create_histograms(hop_result& image, std::list<irectangle2> export_windows, hop_histogram_generator& hoc_generator);
hop_detections_wraper hop_get_detections(hop_result& res_wraper, hop_library& lib_wraper);
hop_detections_wraper hop_get_detections(hop_result& res_wraper, hop_library& lib_wraper, const int only_category);

struct hop_result_wrapper;
struct hop_histogram_descriptor_wraper;
struct hop_detections_wraper;

%{

struct hop_result_wrapper {
    hop_result* results;
    int size;
};

hop_result_wrapper hop_inference(hop_image& image, const char* params)
{
    hop_result_wrapper result;
    result.size = hop_inference(result.results, image, params);   
    return result;
}

hop_result_wrapper hop_inference(hop_image& image, const std::list<irectangle2> mask_regions, const char* params)
{
    hop_result_wrapper result;
    result.size = hop_inference(result.results, image, mask_regions, params);   
    return result;
}

struct hop_histogram_descriptor_wraper {
    std::vector<hop_histogram_descriptor>* results;
    int size;
};
hop_histogram_descriptor_wraper hop_create_histograms(hop_result& image, std::list<irectangle2> export_windows, hop_histogram_generator& hoc_generator) {
  hop_histogram_descriptor_wraper result;
  result.size = hop_create_histograms(result.results, image, export_windows, hoc_generator);
  return result;
}


struct hop_detections_wraper {
   std::list<hop_detection>* bbox_list;
   int size;
};

hop_detections_wraper hop_get_detections(hop_result& res_wraper, hop_library& lib_wraper, const int only_category) {
  hop_detections_wraper result;
  result.size = hop_get_detections(result.bbox_list, res_wraper, lib_wraper, only_category);
  return result;
}

hop_detections_wraper hop_get_detections(hop_result& res_wraper, hop_library& lib_wraper) {
  hop_detections_wraper result;
  result.size = hop_get_detections(result.bbox_list, res_wraper, lib_wraper);
  return result;
}

%}


