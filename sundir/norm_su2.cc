#ifndef __norm_su2_h__
#define __norm_su2_h__

inline void norm_su2(dc *vsu2 ) {

// Takes the first 1 rows of a 2x2 complex matrix,
// and makes an SU(2) matrix out of them.

  double xn;

// Normalize the row `0' of vsu2

  xn = 1./sqrt( norm(vsu2[0])
    + norm(vsu2[1]) );

  vsu2[0]*=xn;
  vsu2[1]*=xn;

// Build row `1' by taking +/- the complex conjugate of the elements of row `0'

  vsu2[2] = -conj( vsu2[1]);
  vsu2[3] =  conj( vsu2[0]);

}

#endif

