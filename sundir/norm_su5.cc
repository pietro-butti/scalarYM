#ifndef __norm_su5_h__
#define __norm_su5_h__

void norm_su5(dc *vsu5 ) {

// Takes the first 4 rows of a 5x5 complex matrix,
// and makes an SU(5) matrix out of them.

  dc determinant_conjugate_phase=dc(1.,0.);
  struct complex_double b[Ncol], DUMMY[1][1], WORK[2*Ncol];
  double AT[2*Ncolsquare];
  int i, j, ok, c1, c2, c3;
  char c4;

  dc scalprod=dc(0.,0.);
  double xn;

// Normalize the row `0' of vsu5

  xn = 1./sqrt( norm(vsu5[0])
    + norm(vsu5[1])
    + norm(vsu5[2])
    + norm(vsu5[3])
    + norm(vsu5[4]) );

  vsu5[0]*=xn;
  vsu5[1]*=xn;
  vsu5[2]*=xn;
  vsu5[3]*=xn;
  vsu5[4]*=xn;

// Make the row `1' orthogonal to the row `0':

  scalprod = conj( vsu5[0] )*vsu5[5]
    + conj( vsu5[1] )*vsu5[6]
    + conj( vsu5[2] )*vsu5[7]
    + conj( vsu5[3] )*vsu5[8]
    + conj( vsu5[4] )*vsu5[9];

  vsu5[5] -= scalprod*vsu5[0];
  vsu5[6] -= scalprod*vsu5[1];
  vsu5[7] -= scalprod*vsu5[2];
  vsu5[8] -= scalprod*vsu5[3];
  vsu5[9] -= scalprod*vsu5[4];


// Normalize the row `1'

  xn = 1./sqrt( norm(vsu5[5])
    + norm(vsu5[6])
    + norm(vsu5[7])
    + norm(vsu5[8])
    + norm(vsu5[9]) );

    vsu5[5]*=xn;
    vsu5[6]*=xn;
    vsu5[7]*=xn;
    vsu5[8]*=xn;
    vsu5[9]*=xn;


// Make the row `2' orthogonal to the row `0':

  scalprod = conj( vsu5[0] )*vsu5[10]
    + conj( vsu5[1] )*vsu5[11]
    + conj( vsu5[2] )*vsu5[12]
    + conj( vsu5[3] )*vsu5[13]
    + conj( vsu5[4] )*vsu5[14];

  vsu5[10] -= scalprod*vsu5[0];
  vsu5[11] -= scalprod*vsu5[1];
  vsu5[12] -= scalprod*vsu5[2];
  vsu5[13] -= scalprod*vsu5[3];
  vsu5[14] -= scalprod*vsu5[4];


// Make the row `2' orthogonal to the row `1':

  scalprod = conj( vsu5[5] )*vsu5[10]
    + conj( vsu5[6] )*vsu5[11]
    + conj( vsu5[7] )*vsu5[12]
    + conj( vsu5[8] )*vsu5[13]
    + conj( vsu5[9] )*vsu5[14];

  vsu5[10] -= scalprod*vsu5[5];
  vsu5[11] -= scalprod*vsu5[6];
  vsu5[12] -= scalprod*vsu5[7];
  vsu5[13] -= scalprod*vsu5[8];
  vsu5[14] -= scalprod*vsu5[9];


// Normalize the row `2'

  xn = 1./sqrt( norm(vsu5[10])
    + norm(vsu5[11])
    + norm(vsu5[12])
    + norm(vsu5[13])
    + norm(vsu5[14]) );

    vsu5[10]*=xn;
    vsu5[11]*=xn;
    vsu5[12]*=xn;
    vsu5[13]*=xn;
    vsu5[14]*=xn;


// Make the row `3' orthogonal to the row `0':

  scalprod = conj( vsu5[0] )*vsu5[15]
    + conj( vsu5[1] )*vsu5[16]
    + conj( vsu5[2] )*vsu5[17]
    + conj( vsu5[3] )*vsu5[18]
    + conj( vsu5[4] )*vsu5[19];

  vsu5[15] -= scalprod*vsu5[0];
  vsu5[16] -= scalprod*vsu5[1];
  vsu5[17] -= scalprod*vsu5[2];
  vsu5[18] -= scalprod*vsu5[3];
  vsu5[19] -= scalprod*vsu5[4];


// Make the row `3' orthogonal to the row `1':

  scalprod = conj( vsu5[5] )*vsu5[15]
    + conj( vsu5[6] )*vsu5[16]
    + conj( vsu5[7] )*vsu5[17]
    + conj( vsu5[8] )*vsu5[18]
    + conj( vsu5[9] )*vsu5[19];

  vsu5[15] -= scalprod*vsu5[5];
  vsu5[16] -= scalprod*vsu5[6];
  vsu5[17] -= scalprod*vsu5[7];
  vsu5[18] -= scalprod*vsu5[8];
  vsu5[19] -= scalprod*vsu5[9];


// Make the row `3' orthogonal to the row `2':

  scalprod = conj( vsu5[10] )*vsu5[15]
    + conj( vsu5[11] )*vsu5[16]
    + conj( vsu5[12] )*vsu5[17]
    + conj( vsu5[13] )*vsu5[18]
    + conj( vsu5[14] )*vsu5[19];

  vsu5[15] -= scalprod*vsu5[10];
  vsu5[16] -= scalprod*vsu5[11];
  vsu5[17] -= scalprod*vsu5[12];
  vsu5[18] -= scalprod*vsu5[13];
  vsu5[19] -= scalprod*vsu5[14];


// Normalize the row `3'

  xn = 1./sqrt( norm(vsu5[15])
    + norm(vsu5[16])
    + norm(vsu5[17])
    + norm(vsu5[18])
    + norm(vsu5[19]) );

    vsu5[15]*=xn;
    vsu5[16]*=xn;
    vsu5[17]*=xn;
    vsu5[18]*=xn;
    vsu5[19]*=xn;


// Make the row `4' orthogonal to the row `0':

  scalprod = conj( vsu5[0] )*vsu5[20]
    + conj( vsu5[1] )*vsu5[21]
    + conj( vsu5[2] )*vsu5[22]
    + conj( vsu5[3] )*vsu5[23]
    + conj( vsu5[4] )*vsu5[24];

  vsu5[20] -= scalprod*vsu5[0];
  vsu5[21] -= scalprod*vsu5[1];
  vsu5[22] -= scalprod*vsu5[2];
  vsu5[23] -= scalprod*vsu5[3];
  vsu5[24] -= scalprod*vsu5[4];


// Make the row `4' orthogonal to the row `1':

  scalprod = conj( vsu5[5] )*vsu5[20]
    + conj( vsu5[6] )*vsu5[21]
    + conj( vsu5[7] )*vsu5[22]
    + conj( vsu5[8] )*vsu5[23]
    + conj( vsu5[9] )*vsu5[24];

  vsu5[20] -= scalprod*vsu5[5];
  vsu5[21] -= scalprod*vsu5[6];
  vsu5[22] -= scalprod*vsu5[7];
  vsu5[23] -= scalprod*vsu5[8];
  vsu5[24] -= scalprod*vsu5[9];


// Make the row `4' orthogonal to the row `2':

  scalprod = conj( vsu5[10] )*vsu5[20]
    + conj( vsu5[11] )*vsu5[21]
    + conj( vsu5[12] )*vsu5[22]
    + conj( vsu5[13] )*vsu5[23]
    + conj( vsu5[14] )*vsu5[24];

  vsu5[20] -= scalprod*vsu5[10];
  vsu5[21] -= scalprod*vsu5[11];
  vsu5[22] -= scalprod*vsu5[12];
  vsu5[23] -= scalprod*vsu5[13];
  vsu5[24] -= scalprod*vsu5[14];


// Make the row `4' orthogonal to the row `3':

  scalprod = conj( vsu5[15] )*vsu5[20]
    + conj( vsu5[16] )*vsu5[21]
    + conj( vsu5[17] )*vsu5[22]
    + conj( vsu5[18] )*vsu5[23]
    + conj( vsu5[19] )*vsu5[24];

  vsu5[20] -= scalprod*vsu5[15];
  vsu5[21] -= scalprod*vsu5[16];
  vsu5[22] -= scalprod*vsu5[17];
  vsu5[23] -= scalprod*vsu5[18];
  vsu5[24] -= scalprod*vsu5[19];


// Normalize the row `4'

  xn = 1./sqrt( norm(vsu5[20])
    + norm(vsu5[21])
    + norm(vsu5[22])
    + norm(vsu5[23])
    + norm(vsu5[24]) );

    vsu5[20]*=xn;
    vsu5[21]*=xn;
    vsu5[22]*=xn;
    vsu5[23]*=xn;
    vsu5[24]*=xn;


// Impose unimodularity by arranging the phase of the last row
// to compensate for the phase of the determinant

  for (i=0; i<Ncol; i++)
  for (j=0; j<Ncol; j++) {
    AT[2*(j+Ncol*i)]=real(vsu5[j*Ncol+i]);
    AT[2*(j+Ncol*i)+1]=imag(vsu5[j*Ncol+i]);
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

  vsu5[20] *= determinant_conjugate_phase;
  vsu5[21] *= determinant_conjugate_phase;
  vsu5[22] *= determinant_conjugate_phase;
  vsu5[23] *= determinant_conjugate_phase;
  vsu5[24] *= determinant_conjugate_phase;


}

#endif
