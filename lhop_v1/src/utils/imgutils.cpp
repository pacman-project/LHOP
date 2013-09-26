//
// ImageMagick <-> img, matrix objects
///////////////////////////////////////////////////////////////////////////////

#include "imgutils.h"

#if defined WIN32 | defined WIN64
void dib32_to_Magick(Magick::Image& result, void* data, int dimx, int dimy, unsigned long* masks)
{
    UCHAR* buffer = new UCHAR[dimx*dimy*3];
    DWORD redMask = *masks;
    DWORD greenMask = *(masks + 1);
    DWORD blueMask = *(masks + 2);
    DWORD* pixels = (DWORD*)data;
    UCHAR* p = buffer + 3 * dimx * (dimy - 1);

    for (int j = dimy; j > 0; --j) {
        for (int i = dimx; i > 0; --i) {
            DWORD color = *pixels;

            *(p++) = (UCHAR)((color & redMask) >> 16);
            *(p++) = (UCHAR)((color & greenMask) >> 8);
            *(p++) = (UCHAR)(color & blueMask);
            ++pixels;
        }
        p -= 6*dimx;
    }
    // XXX: before "Magick::StorageType::CharPixel" was used here but that did not work using GCC
    result.read(dimx, dimy, "RGB", Magick::CharPixel, buffer);
    delete buffer;
}
#endif

void Magick_to_img(Magick::Image& image, img& dest)
{
    unsigned columns = image.columns(), rows = image.rows();

    dest.resize((size_t)columns, (size_t)rows);

    const Magick::PixelPacket* pixels = image.getConstPixels(0, 0, columns, rows);
    img::iterator iter = dest.begin();
    
    for (unsigned i = rows*columns; i > 0; --i) {
        *iter = (double)pixels->red;
        ++pixels; ++iter;
    }
    dest.normalize(0.0, 1.0);
}

void Magick_to_matrix(Magick::Image& image, matrix<float>& dest)
{
     unsigned columns = image.columns(), rows = image.rows();

    dest.resize((size_t)columns, (size_t)rows);

    const Magick::PixelPacket* pixels = image.getConstPixels(0, 0, columns, rows);
	matrix<float>::iterator iter = dest.begin();
    
    for (unsigned i = rows*columns; i > 0; --i) {
        *iter = (double)pixels->red;
        ++pixels; ++iter;
    }
}

void Magick_to_img(Magick::Image& image, img& red, img& green, img& blue)
{
    unsigned int columns = image.columns(), rows = image.rows();

    red.resize((size_t)columns, (size_t)rows);
    green.resize((size_t)columns, (size_t)rows);
    blue.resize((size_t)columns, (size_t)rows);

    const Magick::PixelPacket* pixels = image.getConstPixels(0, 0, columns, rows);
    img::iterator riter = red.begin();
    img::iterator giter = green.begin();
    img::iterator biter = blue.begin();
    
    for (unsigned i = rows*columns; i > 0; --i) {
        *riter = (double)pixels->red;
        *giter = (double)pixels->green;
        *biter = (double)pixels->blue;

        ++pixels; 
        ++riter; ++giter; ++biter;
    }
    red.normalize(0.0, 1.0);
    green.normalize(0.0, 1.0);
    blue.normalize(0.0, 1.0);
}
