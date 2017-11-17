#ifndef CMU462_MATRIX3X3_H
#define CMU462_MATRIX3X3_H

#include <iosfwd>

#include "vector3D.h"

namespace CMU462 {  

/**
 * Defines a 3x3 matrix.
 * 3x3 matrices are extremely useful in computer graphics. 
 */
class Matrix3x3 {
  
  public:

  /** 
   * Sets all elements to val.
   */
  void zero(double val = 0.0 );

  /** 
   * Returns the determinant of A.
   */
  double det( void ) const;

  /** 
   * Returns the Frobenius norm of A.
   */
  double norm( void ) const;

  /** 
   * Returns the 3x3 identity matrix.
   */
  static Matrix3x3 identity( void );

  /** 
   * Returns a matrix representing the (left) cross product with u.
   */
  static Matrix3x3 crossProduct( const Vector3D& u );

  /**
   * Returns the ith column.
   */
        Vector3D& column( int i );
  const Vector3D& column( int i ) const;

  /** 
   * Returns the transpose of A.
   */
  Matrix3x3 T( void ) const;

  /** 
   * Returns the inverse of A.
   */
  Matrix3x3 inv( void ) const;

  // accesses element (i,j) of A using 0-based indexing
        double& operator()( int i, int j );
  const double& operator()( int i, int j ) const;

  // accesses the ith column of A   
        Vector3D& operator[]( int i );
  const Vector3D& operator[]( int i ) const;

  // increments by B
  void operator+=( const Matrix3x3& B );

  // returns -A
  Matrix3x3 operator-( void ) const;

  // returns A-B
  Matrix3x3 operator-( const Matrix3x3& B ) const;

  // returns c*A
  Matrix3x3 operator*( double c ) const;

  // returns A*B
  Matrix3x3 operator*( const Matrix3x3& B ) const;

  // returns A*x
  Vector3D operator*( const Vector3D& x ) const;

  // divides each element by x
  void operator/=( double x );   

  void setEntry(int i, int j, double val){
	if(j < 0 || j >= 3) return;
	switch(i){
		case 0:
			entries[j].x = val;
			break;
		case 1:
			entries[j].y = val;
			break;
		case 2:
			entries[j].z = val;
			break;
		default:
			return;	
	}
  }

  protected:

  // column vectors
  Vector3D entries[3];

}; // class Matrix3x3

// returns the outer product of u and v
Matrix3x3 outer( const Vector3D& u, const Vector3D& v );

// returns c*A
Matrix3x3 operator*( double c, const Matrix3x3& A );

// prints entries
std::ostream& operator<<( std::ostream& os, const Matrix3x3& A );
    
} // namespace CMU462

#endif // CMU462_MATRIX3X3_H