#ifndef __norm_su3_h__
#define __norm_su3_h__

void norm_su3(dc *vsu3 ) {

// Takes the first 2 rows of a 3x3 complex matrix,
// and makes an SU(3) matrix out of them.

  dc scalprod=dc(0.,0.);

  double xn;

// Normalize the row `0' of vsu3

  xn = 1./sqrt( norm(vsu3[0])
    + norm(vsu3[1]) 
    + norm(vsu3[2]) );

  vsu3[0]*=xn;
  vsu3[1]*=xn;
  vsu3[2]*=xn;

// Make the row `1' orthogonal to the row `0':

  scalprod = conj( vsu3[0] )*vsu3[3]
    + conj( vsu3[1] )*vsu3[4]
    + conj( vsu3[2] )*vsu3[5];

  vsu3[3] -= scalprod*vsu3[0];
  vsu3[4] -= scalprod*vsu3[1];
  vsu3[5] -= scalprod*vsu3[2];

// Normalize row `1':

  xn = 1./sqrt(  norm( vsu3[3] )
    + norm( vsu3[4] )
    + norm( vsu3[5] ) );

  vsu3[3] *= xn;
  vsu3[4] *= xn;
  vsu3[5] *= xn;

// Construct row `2' as the complex conjugate of the cross-product of the first two.
// (Automatically normalised, and orthogonal to first two rows.

  vsu3[6] = conj( vsu3[1]*vsu3[5] - vsu3[2]*vsu3[4] );
  vsu3[7] = conj( vsu3[2]*vsu3[3] - vsu3[0]*vsu3[5] );
  vsu3[8] = conj( vsu3[0]*vsu3[4] - vsu3[1]*vsu3[3] );



}

#endif

