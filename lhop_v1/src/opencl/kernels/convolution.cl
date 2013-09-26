// include files
#include <utils.cl>

#if __DEVICE_LOCAL_MEM_TYPE__ > 0
	#define HAS_FAST_LOCAL_MEM 1
#endif



// functions image_convolve_vector16, image_convolve_vector8 and void image_convolve_vector4 MUST be defined in seperate file due to bug in NVIDIA compiler (__constant buffer gets corrupted)
__attribute__((vec_type_hint(float2)))
__kernel
void image_convolve_vector2(__global float* image, int2 image_size,
							__constant float* convol_mask, int4 convol_mask_size,
							__local float* local_data, int2 local_data_size,
							uint save_max, __global float* max_image, __global float* max_image_index,
							__global float* output, uint layer1_creator_type, __global int* error) {

	// define vector size 2 and include body of actualy function for convolution
	
	#define VECTOR_SIZE 2
	
	#include <core_convolution.cl>
	
	#undef 	VECTOR_SIZE
	
}


#ifdef __DEVICE_IMAGE_SUPPORTED__

const sampler_t sampler_my = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE;

__kernel
void image_convolve_img(__read_only image2d_t image, int2 image_size,
						__read_only image2d_t convol_mask, int4 convol_mask_size,
						__write_only image2d_t output, __global int* error) {
	
	// only two dimensions are allowed
	if (get_work_dim() != 2) {
		*error = 1;
		return;
	}

	int2 convol_mask_size_half = convol_mask_size.xy/(int2)2;

	if ((int)get_global_id(1) >= convol_mask_size_half.y && (int)get_global_id(1) < image_size.y - convol_mask_size_half.y &&
		(int)get_global_id(0) >= convol_mask_size_half.x && (int)get_global_id(0) < image_size.x - convol_mask_size_half.x) {
	
		int2 img_pos_start = (int2)( (int)get_global_id(0) - convol_mask_size_half.x, (int)get_global_id(1) - convol_mask_size_half.y);
		int2 mask_pos_start = convol_mask_size.xy - (int2)(1, 1);
		
		for (int convol_i = 0; convol_i < convol_mask_size.z; convol_i++) { 
		
			int2 img_pos = img_pos_start;
			int2 mask_pos = mask_pos_start;

			float sum = 0;
						
			for (int i = 0; i < convol_mask_size.y; i++) {
				
				img_pos.x = img_pos_start.x;
				mask_pos.x = mask_pos_start.x;
				
				for (int j = 0; j < convol_mask_size.x; j++) {
					sum += read_imagef(convol_mask, sampler_my, mask_pos).x * read_imagef(image, sampler_my, img_pos).x;
					img_pos.x++;
					mask_pos.x--;
				}
				
				img_pos.y++;
				mask_pos.y--;
			}			
			
			// save value back to global memory
			write_imagef(output, (int2)(get_global_id(0), get_global_id(1) + convol_i * image_size.y), (float4)(sum,0,0,1) );
			
			// move to next convolution kernel
			mask_pos_start.y += convol_mask_size.y;
		}

	}	
}


#else

__kernel
void image_convolve_img() {

}

#endif

