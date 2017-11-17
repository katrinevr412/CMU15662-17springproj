#include "texture.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE(sky): 
  // The starter code allocates the mip levels and generates a level 
  // map simply fills each level with a color that differs from its
  // neighbours'. The reference solution uses trilinear filtering
  // and it will only work when you have mipmaps.

  // Task 7: Implement this

  // check start level
	
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }
	//std::cout << "mipmap size: " << tex.mipmap[startLevel].texels.size() << std::endl;
  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;

  	while(1){
		baseWidth >>= 1;
		baseHeight >>= 1;
		if(!baseWidth || !baseHeight) break;
		MipLevel newML;
		newML.width = baseWidth;
		newML.height = baseHeight;
		newML.texels.clear();
		for(int i = 0; i < baseWidth * baseHeight; i++){
			vector<unsigned char> &tmp = tex.mipmap[startLevel].texels;
			int thisW = i % baseWidth, thisH = i / baseWidth;
			int buf[4] = {0};
			for(int q = 0; q < 4; q++){
				for(int r = 0; r < 2; r++){
					for(int s = 0; s < 2; s++){
						buf[q] += (int)tmp[4 * (2 * thisW + r + (2 * thisH + s) * 2 * baseWidth) + q];
					}
				}
			}
			for(int q = 0; q < 4; q++){
				newML.texels.push_back((unsigned char)(buf[q]/4.0));
			}
		}
		startLevel++;
		tex.mipmap.push_back(newML);
	}
	/*
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }
   //std::cout << "level size: " << tex.mipmap.size() << std::endl;

  */
}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {
  //cout << "level=" << level << endl;
  // Task ?: Implement nearest neighbour interpolation
  
  // return magenta for invalid level
  //cout << "u=" << u << ",v=" << v << endl;
  //return Color(150/255.0, 150/255.0, 150/255.0);
  if(level < 0 || level >= tex.mipmap.size()) return nearestAlertColor;
  //int tw = tex.width, th = tex.height;
  int mw = tex.mipmap[level].width, mh = tex.mipmap[level].height;
  int posW = (int) floor(u * mw), posH = (int) floor(v * mh);
  //cout << "posW=" << posW << ",posH=" << posH << endl;
  int offset = posW + posH * mw;
	//cout << "offset=" << offset << endl;
  std::vector<unsigned char> &tmp = tex.mipmap[level].texels;
  Color eh = Color(tmp[4*offset]/255.0, tmp[4*offset+1]/255.0, tmp[4*offset+2]/255.0, tmp[4*offset+3]/255.0);
	//cout << eh << endl;
	return eh;
}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
  // Task ?: Implement bilinear filtering

  // return magenta for invalid level
  if(level < 0 || level >= tex.mipmap.size()) return bilinearAlertColor;

	int mw = tex.mipmap[level].width, mh = tex.mipmap[level].height;
	u *= mw, v *= mh;
	int posW1 = (int) floor(u), posH1 = (int) floor(v);
	int posW2 = posW1 + 1, posH2 = posH1 + 1;
	int offset0 = posW1 + posH1 * mw, offset1 = posW2 + posH1 * mw, offset2 = posW1 + posH2 * mw, offset3 = posW2 + posH2 * mw;
	double wRatio = u - posW1, hRatio = v - posH1;
	std::vector<unsigned char> &tmp = tex.mipmap[level].texels;
	Color c1 = Color((tmp[4*offset0+0] * (1-wRatio) + tmp[4*offset1+0] * wRatio)/255.0,
			 (tmp[4*offset0+1] * (1-wRatio) + tmp[4*offset1+1] * wRatio)/255.0,
			 (tmp[4*offset0+2] * (1-wRatio) + tmp[4*offset1+2] * wRatio)/255.0,
			 (tmp[4*offset0+3] * (1-wRatio) + tmp[4*offset1+3] * wRatio)/255.0);
	Color c2 = Color((tmp[4*offset2+0] * (1-wRatio) + tmp[4*offset3+0] * wRatio)/255.0,
			 (tmp[4*offset2+1] * (1-wRatio) + tmp[4*offset3+1] * wRatio)/255.0,
			 (tmp[4*offset2+2] * (1-wRatio) + tmp[4*offset3+2] * wRatio)/255.0,
			 (tmp[4*offset2+3] * (1-wRatio) + tmp[4*offset3+3] * wRatio)/255.0);
	return Color(c1.r * (1-hRatio) + c2.r * hRatio,
		     c1.g * (1-hRatio) + c2.g * hRatio,
		     c1.b * (1-hRatio) + c2.b * hRatio,
		     c1.a * (1-hRatio) + c2.a * hRatio);
}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

  // Task 8: Implement trilinear filtering
	if(u_scale <= 1.0/tex.width || v_scale <= 1.0/tex.height){
		return sample_bilinear(tex, u, v);
	}
  // return magenta for invalid level
        if(u_scale > 1 || v_scale > 1) return trilinearAlertColor;
  	int level = 0;
	for(; level < tex.mipmap.size(); level++){
		if(u_scale <= (1 << level) * 1.0 / tex.width || v_scale <= (1 << level) * 1.0 / tex.height){
			break;
		}
	}
	//cout << "mipmap size: " << tex.mipmap.size() << "level: " << level << endl;
	double ratio1 = (tex.width * 1.0 / (1 << (level-1)) - 1.0 / u_scale) / (tex.width * 1.0 / (1 << level));
	double ratio2 = (tex.height * 1.0 / (1 << (level-1)) - 1.0 / v_scale) / (tex.height * 1.0 / (1 << level));
	double ratio = (ratio1 + ratio2) / 2;
	Color c1 = sample_bilinear(tex, u, v, level-1), c2 = sample_bilinear(tex, u, v, level);
	Color res(c1.r * (1-ratio) + c2.r * ratio,
		     c1.g * (1-ratio) + c2.g * ratio,
		     c1.b * (1-ratio) + c2.b * ratio,
		     c1.a * (1-ratio) + c2.a * ratio);
	//cout << res << endl;
	return res;
}

} // namespace CMU462
