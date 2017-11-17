#ifndef CMU462_TEXTURE_H
#define CMU462_TEXTURE_H

#include "CMU462.h"
#include <iostream>

namespace CMU462 {

static const int kMaxMipLevels = 14;

typedef enum SampleMethod{
  NEAREST,
  BILINEAR,
  TRILINEAR
} SampleMethod;

struct MipLevel {
  size_t width; 
  size_t height;
  std::vector<unsigned char> texels;
};

struct Texture {
  size_t width;
  size_t height;
  std::vector<MipLevel> mipmap;
	void print(){
		std::cout << "texture width: " << width << ", texture height: " << height << std::endl;
		int sz = mipmap.size();
		std::cout << "contains " << sz << " mips" << std::endl;
		for(int i = 0; i < sz; i++){
			std::cout << "	mip No. " << i+1 << ":" << std::endl;
			std::cout << "		width: " << mipmap[i].width << ", height:" << mipmap[i].height << std::endl;
			
			int szz = mipmap[i].texels.size();
			std::cout << "		there are " << szz << " texels." << std::endl;
			//if(szz > 1024) continue;
			//std::cout << "			";
			int buf[4] = {-1,-1,-1,-1};
			int cnt = 0;
			for(int j = 0; j < szz/4; j++){
				bool diff = false;
				
				for(int q = 0; q < 4; q++){
					if((int)mipmap[i].texels[4*j+q] != buf[q]){
						diff = true;
						//buf[q] = (int)mipmap[i].texels[4*j+q];
					}
				}
				if(diff){
					if(buf[0] != -1) std::cout << "			" << j-cnt << "-th: (" << buf[0] << "," << buf[1]
					<< "," << buf[2] << "," << buf[3] << "), "  << cnt << "times" << std::endl;
					cnt = 1;
					for(int q = 0; q < 4; q++){
						buf[q] = (int)mipmap[i].texels[4*j+q];
					}
				}
				else cnt++;
				
			}
			std::cout << "			" << szz/4-cnt << "-th: (" << buf[0] << "," << buf[1]
					<< "," << buf[2] << "," << buf[3] << "), "  << cnt << "times" << std::endl;
			std::cout << std::endl << std::endl;
		}
	}
};

class Sampler2D {
 public:

  Sampler2D( SampleMethod method ) : method ( method ) { }

  ~Sampler2D(){}

  virtual void generate_mips( Texture& tex, int startLevel ) = 0;

  virtual Color sample_nearest(Texture& tex, 
                               float u, float v, 
                               int level = 0) = 0;

  virtual Color sample_bilinear(Texture& tex, 
                                float u, float v, 
                                int level = 0) = 0;

  virtual Color sample_trilinear(Texture& tex, 
                                 float u, float v, 
                                 float u_scale, float v_scale) = 0;
  

  inline SampleMethod get_sample_method() const {
    return method;
  }
 
 protected:

  SampleMethod method;

}; // class Sampler2D

class Sampler2DImp : public Sampler2D {
 public:

  Sampler2DImp( SampleMethod method = TRILINEAR ) : Sampler2D ( method ) { }
  
  void generate_mips( Texture& tex, int startLevel );
  const Color nearestAlertColor = Color(0,0,1,1);
  const Color bilinearAlertColor = Color(0,1,0,1);
  const Color trilinearAlertColor = Color(1,0,0,1);

  Color sample_nearest(Texture& tex, 
                       float u, float v, 
                       int level = 0);

  Color sample_bilinear(Texture& tex, 
                        float u, float v, 
                        int level = 0);

  Color sample_trilinear(Texture& tex, 
                         float u, float v, 
                         float u_scale, float v_scale);
  
  Color sample(Texture& tex, float u, float v, int level = 0, float u_scale = 0, float v_scale = 0){
	//std::cout << "u-scale: " << u_scale << ", v-scale: " << v_scale << std::endl;
	switch(method){
		case NEAREST: return sample_nearest(tex, u, v, level);
		case BILINEAR: return sample_bilinear(tex, u, v, level);
		case TRILINEAR: return sample_trilinear(tex, u, v, u_scale, v_scale);
	}
  }
    

}; // class sampler2DImp

class Sampler2DRef : public Sampler2D {
 public:

  Sampler2DRef( SampleMethod method = TRILINEAR ) : Sampler2D ( method ) { }
  
  void generate_mips( Texture& tex, int startLevel );

  Color sample_nearest(Texture& tex, 
                       float u, float v, 
                       int level = 0);

  Color sample_bilinear(Texture& tex, 
                        float u, float v, 
                        int level = 0);

  Color sample_trilinear(Texture& tex, 
                         float u, float v, 
                         float u_scale, float v_scale);
  
}; // class sampler2DRef

} // namespace CMU462

#endif // CMU462_TEXTURE_H
