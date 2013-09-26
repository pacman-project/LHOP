//
// ImageMagick <-> img, matrix objects
///////////////////////////////////////////////////////////////////////////////

#include <Magick++.h>
#include "img.h"

#if defined WIN32 | defined WIN64
void dib32_to_Magick(Magick::Image& result, void* data, int dimx, int dimy, unsigned long* masks);
#endif

void Magick_to_img(Magick::Image& image, img& dest);

void Magick_to_matrix(Magick::Image& image, matrix<float>& dest);

void Magick_to_img(Magick::Image& image, img& red, img& green, img& blue);


