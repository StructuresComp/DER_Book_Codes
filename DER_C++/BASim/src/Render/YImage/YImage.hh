/**
 * \file YImage.hh
 */

#ifndef YIMAGE_HH
#define YIMAGE_HH

namespace BASim {

/// file loading/saving automatically picks up changes to this
/// struct. The possibilities are: ARGB, ABGR, RGBA, BGRA.
struct YPixel {
  unsigned char r ;
  unsigned char g ;
  unsigned char b ;
  unsigned char a ;
};

/** Class for saving PNG images written by Yotam Gingold. */
class YImage
{
public:
  YImage() ;
  YImage( const YImage& ) ;
  virtual ~YImage() ;

  YImage& operator=( const YImage& ) ;

  bool save( const char* fname ) const ;
  bool load( const char* fname ) ;

  YPixel* data() ;
  const YPixel* data() const ;

  YPixel& at( int i, int j ) ;
  const YPixel& at( int i, int j ) const ;

  int width() const ;
  int height() const ;
  void resize( int width, int height ) ;

  /// flip vertically
  void flip() ;
  /// flip horizontally
  void mirror() ;
  /// average rgb
  void greyscale() ;
  /// sets the alpha channel of each pixel to alpha
  void setAllAlpha( unsigned char alpha );

protected:
  int m_width ;
  int m_height ;
  YPixel* m_data ; ///< raw image data
};

} // namespace BASim

#endif // YIMAGE_HH
