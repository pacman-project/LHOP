

#ifdef HAS_FAST_LOCAL_MEM
	int img_local_pos_offset = max(0, img_local_pos_offset_start);
#else
	int img_global_pos_offset = img_global_pos_offset_start;
#endif
	
	int mask_offset = convol_mask_size.x * convol_mask_size.y * (convol_i+1) - 1 ;
	
	vector_float sum = (vector_float)0;
	
	for (int i = 0; i < convol_mask_size.y; i++) {
		for (int j = 0; j < convol_mask_size.x; j++) {

#ifdef HAS_FAST_LOCAL_MEM
			sum = mad((vector_float)(convol_mask[mask_offset]), vector_load(0, local_data + img_local_pos_offset), sum);
			img_local_pos_offset++;
#else
			sum = mad((vector_float)(convol_mask[mask_offset]), vector_load(0, image + img_global_pos_offset), sum);			
			img_global_pos_offset++;
#endif
			mask_offset--;
		}
#ifdef HAS_FAST_LOCAL_MEM					
		img_local_pos_offset += skipped_local_pos_x;
#else
		img_global_pos_offset += skipped_global_pos_x;
#endif
	}
