#ifndef __norm_su4_h__
#define __norm_su4_h__

void norm_su4(dc *vsu4 ) {

// Takes the first 3 rows of a 4x4 complex matrix,
// and makes an SU(4) matrix out of them.

  dc determinant_conjugate_phase=dc(1.,0.);

  struct complex_double b[Ncol], DUMMY[1][1], WORK[2*Ncol];
  double AT[2*Ncolsquare];
  int i, j, ok, c1, c2, c3;
  char c4;

  dc scalprod=dc(0.,0.);

  double xn;

// Normalize the row `0' of vsu4

  xn = 1./sqrt( norm(vsu4[0])
    + norm(vsu4[1]) 
    + norm(vsu4[2]) 
    + norm(vsu4[3]) );

  vsu4[0]*=xn;
  vsu4[1]*=xn;
  vsu4[2]*=xn;
  vsu4[3]*=xn;

// Make the row `1' orthogonal to the row `0':

  scalprod = conj( vsu4[0] )*vsu4[4]
    + conj( vsu4[1] )*vsu4[5]
    + conj( vsu4[2] )*vsu4[6]
    + conj( vsu4[3] )*vsu4[7];

  vsu4[4] -= scalprod*vsu4[0];
  vsu4[5] -= scalprod*vsu4[1];
  vsu4[6] -= scalprod*vsu4[2];
  vsu4[7] -= scalprod*vsu4[3];

// Normalize row `1':

  xn = 1./sqrt(  norm( vsu4[4] )
    + norm( vsu4[5] )
    + norm( vsu4[6] )
    + norm( vsu4[7] ) );

  vsu4[4] *= xn;
  vsu4[5] *= xn;
  vsu4[6] *= xn;
  vsu4[7] *= xn;



// Make the row `2' orthogonal to the row `0':

  scalprod = conj( vsu4[0] )*vsu4[8]
    + conj( vsu4[1] )*vsu4[9]
    + conj( vsu4[2] )*vsu4[10]
    + conj( vsu4[3] )*vsu4[11];

  vsu4[8] -= scalprod*vsu4[0];
  vsu4[9] -= scalprod*vsu4[1];
  vsu4[10] -= scalprod*vsu4[2];
  vsu4[11] -= scalprod*vsu4[3];


// Make the row `2' orthogonal to the row `1':

  scalprod = conj( vsu4[4] )*vsu4[8]
    + conj( vsu4[5] )*vsu4[9]
    + conj( vsu4[6] )*vsu4[10]
    + conj( vsu4[7] )*vsu4[11];

  vsu4[8] -= scalprod*vsu4[4];
  vsu4[9] -= scalprod*vsu4[5];
  vsu4[10] -= scalprod*vsu4[6];
  vsu4[11] -= scalprod*vsu4[7];

// Normalize row `2':

  xn = 1./sqrt(  norm( vsu4[8] )
    + norm( vsu4[9] )
    + norm( vsu4[10] )
    + norm( vsu4[11] ) );

  vsu4[8] *= xn;
  vsu4[9] *= xn;
  vsu4[10] *= xn;
  vsu4[11] *= xn;



// Make the row `3' orthogonal to the row `0':

  scalprod = conj( vsu4[0] )*vsu4[12]
    + conj( vsu4[1] )*vsu4[13]
    + conj( vsu4[2] )*vsu4[14]
    + conj( vsu4[3] )*vsu4[15];

  vsu4[12] -= scalprod*vsu4[0];
  vsu4[13] -= scalprod*vsu4[1];
  vsu4[14] -= scalprod*vsu4[2];
  vsu4[15] -= scalprod*vsu4[3];


// Make the row `3' orthogonal to the row `1':

  scalprod = conj( vsu4[4] )*vsu4[12]
    + conj( vsu4[5] )*vsu4[13]
    + conj( vsu4[6] )*vsu4[14]
    + conj( vsu4[7] )*vsu4[15];

  vsu4[12] -= scalprod*vsu4[4];
  vsu4[13] -= scalprod*vsu4[5];
  vsu4[14] -= scalprod*vsu4[6];
  vsu4[15] -= scalprod*vsu4[7];


// Make the row `3' orthogonal to the row `2':

  scalprod = conj( vsu4[8] )*vsu4[12]
    + conj( vsu4[9] )*vsu4[13]
    + conj( vsu4[10] )*vsu4[14]
    + conj( vsu4[11] )*vsu4[15];

  vsu4[12] -= scalprod*vsu4[8];
  vsu4[13] -= scalprod*vsu4[9];
  vsu4[14] -= scalprod*vsu4[10];
  vsu4[15] -= scalprod*vsu4[11];

// Normalize row `3':

  xn = 1./sqrt(  norm( vsu4[12] )
    + norm( vsu4[13] )
    + norm( vsu4[14] )
    + norm( vsu4[15] ) );

  vsu4[12] *= xn;
  vsu4[13] *= xn;
  vsu4[14] *= xn;
  vsu4[15] *= xn;

// Impose unimodularity by arranging the phase of the last row
// to compensate for the phase of the determinant



  for (i=0; i<Ncol; i++)
  for(j=0; j<Ncol; j++) {
    AT[2*(j+Ncol*i)]=real(vsu4[j*Ncol+i]);
    AT[2*(j+Ncol*i)+1]=imag(vsu4[j*Ncol+i]);
  }
  c1=Ncol;
  c2=2*Ncol;
  c3=1;
  c4='N';

  zgeev_(&c4, &c4,&c1, AT, &c1, b, DUMMY, &c3, DUMMY, &c3, WORK, &c2, WORK, &ok);

  dc determinante=dc(1.,0.);
  if (ok==0) {
    for (i=0; i<Ncol; i++) {
      determinante*=dc(b[i].re, b[i].im);
    }
  }
  else { printf("An error occured"); exit(0); }

  determinant_conjugate_phase=conj(determinante);


  vsu4[12] *= determinant_conjugate_phase;
  vsu4[13] *= determinant_conjugate_phase;
  vsu4[14] *= determinant_conjugate_phase;
  vsu4[15] *= determinant_conjugate_phase;

}

#endif

