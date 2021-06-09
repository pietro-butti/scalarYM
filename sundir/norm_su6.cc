#ifndef __norm_su5_h__
#define __norm_su5_h__

void norm_su6(dc *vsu6 ) {

// Takes the first 5 rows of a 6x6 complex matrix,
// and makes an SU(6) matrix out of them.

  dc determinant_conjugate_phase=dc(1.,0.);
  struct complex_double b[Ncol], DUMMY[1][1], WORK[2*Ncol];
  double AT[2*Ncolsquare];
  int i, j, ok, c1, c2, c3;
  char c4;

  dc scalprod=dc(0.,0.);
  double xn;

// Normalize the row `0' of vsu6

  xn = 1./sqrt( norm(vsu6[0])
    + norm(vsu6[1])
    + norm(vsu6[2])
    + norm(vsu6[3])
    + norm(vsu6[4])
    + norm(vsu6[5]) );

  vsu6[0]*=xn;
  vsu6[1]*=xn;
  vsu6[2]*=xn;
  vsu6[3]*=xn;
  vsu6[4]*=xn;
  vsu6[5]*=xn;

// Make the row `1' orthogonal to the row `0':

  scalprod = conj( vsu6[0] )*vsu6[6]
    + conj( vsu6[1] )*vsu6[7]
    + conj( vsu6[2] )*vsu6[8]
    + conj( vsu6[3] )*vsu6[9]
    + conj( vsu6[4] )*vsu6[10]
    + conj( vsu6[5] )*vsu6[11];

  vsu6[6] -= scalprod*vsu6[0];
  vsu6[7] -= scalprod*vsu6[1];
  vsu6[8] -= scalprod*vsu6[2];
  vsu6[9] -= scalprod*vsu6[3];
  vsu6[10] -= scalprod*vsu6[4];
  vsu6[11] -= scalprod*vsu6[5];


// Normalize the row `1'

  xn = 1./sqrt( norm(vsu6[6])
    + norm(vsu6[7])
    + norm(vsu6[8])
    + norm(vsu6[9])
    + norm(vsu6[10])
    + norm(vsu6[11]) );

    vsu6[6]*=xn;
    vsu6[7]*=xn;
    vsu6[8]*=xn;
    vsu6[9]*=xn;
    vsu6[10]*=xn;
    vsu6[11]*=xn;


// Make the row `2' orthogonal to the row `0':

  scalprod = conj( vsu6[0] )*vsu6[12]
    + conj( vsu6[1] )*vsu6[13]
    + conj( vsu6[2] )*vsu6[14]
    + conj( vsu6[3] )*vsu6[15]
    + conj( vsu6[4] )*vsu6[16]
    + conj( vsu6[5] )*vsu6[17];

  vsu6[12] -= scalprod*vsu6[0];
  vsu6[13] -= scalprod*vsu6[1];
  vsu6[14] -= scalprod*vsu6[2];
  vsu6[15] -= scalprod*vsu6[3];
  vsu6[16] -= scalprod*vsu6[4];
  vsu6[17] -= scalprod*vsu6[5];


// Make the row `2' orthogonal to the row `1':

  scalprod = conj( vsu6[6] )*vsu6[12]
    + conj( vsu6[7] )*vsu6[13]
    + conj( vsu6[8] )*vsu6[14]
    + conj( vsu6[9] )*vsu6[15]
    + conj( vsu6[10] )*vsu6[16]
    + conj( vsu6[11] )*vsu6[17];

  vsu6[12] -= scalprod*vsu6[6];
  vsu6[13] -= scalprod*vsu6[7];
  vsu6[14] -= scalprod*vsu6[8];
  vsu6[15] -= scalprod*vsu6[9];
  vsu6[16] -= scalprod*vsu6[10];
  vsu6[17] -= scalprod*vsu6[11];


// Normalize the row `2'

  xn = 1./sqrt( norm(vsu6[12])
    + norm(vsu6[13])
    + norm(vsu6[14])
    + norm(vsu6[15])
    + norm(vsu6[16])
    + norm(vsu6[17]) );

    vsu6[12]*=xn;
    vsu6[13]*=xn;
    vsu6[14]*=xn;
    vsu6[15]*=xn;
    vsu6[16]*=xn;
    vsu6[17]*=xn;


// Make the row `3' orthogonal to the row `0':

  scalprod = conj( vsu6[0] )*vsu6[18]
    + conj( vsu6[1] )*vsu6[19]
    + conj( vsu6[2] )*vsu6[20]
    + conj( vsu6[3] )*vsu6[21]
    + conj( vsu6[4] )*vsu6[22]
    + conj( vsu6[5] )*vsu6[23];

  vsu6[18] -= scalprod*vsu6[0];
  vsu6[19] -= scalprod*vsu6[1];
  vsu6[20] -= scalprod*vsu6[2];
  vsu6[21] -= scalprod*vsu6[3];
  vsu6[22] -= scalprod*vsu6[4];
  vsu6[23] -= scalprod*vsu6[5];


// Make the row `3' orthogonal to the row `1':

  scalprod = conj( vsu6[6] )*vsu6[18]
    + conj( vsu6[7] )*vsu6[19]
    + conj( vsu6[8] )*vsu6[20]
    + conj( vsu6[9] )*vsu6[21]
    + conj( vsu6[10] )*vsu6[22]
    + conj( vsu6[11] )*vsu6[23];

  vsu6[18] -= scalprod*vsu6[6];
  vsu6[19] -= scalprod*vsu6[7];
  vsu6[20] -= scalprod*vsu6[8];
  vsu6[21] -= scalprod*vsu6[9];
  vsu6[22] -= scalprod*vsu6[10];
  vsu6[23] -= scalprod*vsu6[11];


// Make the row `3' orthogonal to the row `2':

  scalprod = conj( vsu6[12] )*vsu6[18]
    + conj( vsu6[13] )*vsu6[19]
    + conj( vsu6[14] )*vsu6[20]
    + conj( vsu6[15] )*vsu6[21]
    + conj( vsu6[16] )*vsu6[22]
    + conj( vsu6[17] )*vsu6[23];

  vsu6[18] -= scalprod*vsu6[12];
  vsu6[19] -= scalprod*vsu6[13];
  vsu6[20] -= scalprod*vsu6[14];
  vsu6[21] -= scalprod*vsu6[15];
  vsu6[22] -= scalprod*vsu6[16];
  vsu6[23] -= scalprod*vsu6[17];


// Normalize the row `3'

  xn = 1./sqrt( norm(vsu6[18])
    + norm(vsu6[19])
    + norm(vsu6[20])
    + norm(vsu6[21])
    + norm(vsu6[22])
    + norm(vsu6[23]) );

    vsu6[18]*=xn;
    vsu6[19]*=xn;
    vsu6[20]*=xn;
    vsu6[21]*=xn;
    vsu6[22]*=xn;
    vsu6[23]*=xn;


// Make the row `4' orthogonal to the row `0':

  scalprod = conj( vsu6[0] )*vsu6[24]
    + conj( vsu6[1] )*vsu6[25]
    + conj( vsu6[2] )*vsu6[26]
    + conj( vsu6[3] )*vsu6[27]
    + conj( vsu6[4] )*vsu6[28]
    + conj( vsu6[5] )*vsu6[29];

  vsu6[24] -= scalprod*vsu6[0];
  vsu6[25] -= scalprod*vsu6[1];
  vsu6[26] -= scalprod*vsu6[2];
  vsu6[27] -= scalprod*vsu6[3];
  vsu6[28] -= scalprod*vsu6[4];
  vsu6[29] -= scalprod*vsu6[5];


// Make the row `4' orthogonal to the row `1':

  scalprod = conj( vsu6[6] )*vsu6[24]
    + conj( vsu6[7] )*vsu6[25]
    + conj( vsu6[8] )*vsu6[26]
    + conj( vsu6[9] )*vsu6[27]
    + conj( vsu6[10] )*vsu6[28]
    + conj( vsu6[11] )*vsu6[29];

  vsu6[24] -= scalprod*vsu6[6];
  vsu6[25] -= scalprod*vsu6[7];
  vsu6[26] -= scalprod*vsu6[8];
  vsu6[27] -= scalprod*vsu6[9];
  vsu6[28] -= scalprod*vsu6[10];
  vsu6[29] -= scalprod*vsu6[11];


// Make the row `4' orthogonal to the row `2':

  scalprod = conj( vsu6[12] )*vsu6[24]
    + conj( vsu6[13] )*vsu6[25]
    + conj( vsu6[14] )*vsu6[26]
    + conj( vsu6[15] )*vsu6[27]
    + conj( vsu6[16] )*vsu6[28]
    + conj( vsu6[17] )*vsu6[29];

  vsu6[24] -= scalprod*vsu6[12];
  vsu6[25] -= scalprod*vsu6[13];
  vsu6[26] -= scalprod*vsu6[14];
  vsu6[27] -= scalprod*vsu6[15];
  vsu6[28] -= scalprod*vsu6[16];
  vsu6[29] -= scalprod*vsu6[17];


// Make the row `4' orthogonal to the row `3':

  scalprod = conj( vsu6[18] )*vsu6[24]
    + conj( vsu6[19] )*vsu6[25]
    + conj( vsu6[20] )*vsu6[26]
    + conj( vsu6[21] )*vsu6[27]
    + conj( vsu6[22] )*vsu6[28]
    + conj( vsu6[23] )*vsu6[29];

  vsu6[24] -= scalprod*vsu6[18];
  vsu6[25] -= scalprod*vsu6[19];
  vsu6[26] -= scalprod*vsu6[20];
  vsu6[27] -= scalprod*vsu6[21];
  vsu6[28] -= scalprod*vsu6[22];
  vsu6[29] -= scalprod*vsu6[23];


// Normalize the row `4'

  xn = 1./sqrt( norm(vsu6[24])
    + norm(vsu6[25])
    + norm(vsu6[26])
    + norm(vsu6[27])
    + norm(vsu6[28])
    + norm(vsu6[29]) );

    vsu6[24]*=xn;
    vsu6[25]*=xn;
    vsu6[26]*=xn;
    vsu6[27]*=xn;
    vsu6[28]*=xn;
    vsu6[29]*=xn;


// Make the row `5' orthogonal to the row `0':

  scalprod = conj( vsu6[0] )*vsu6[30]
    + conj( vsu6[1] )*vsu6[31]
    + conj( vsu6[2] )*vsu6[32]
    + conj( vsu6[3] )*vsu6[33]
    + conj( vsu6[4] )*vsu6[34]
    + conj( vsu6[5] )*vsu6[35];

  vsu6[30] -= scalprod*vsu6[0];
  vsu6[31] -= scalprod*vsu6[1];
  vsu6[32] -= scalprod*vsu6[2];
  vsu6[33] -= scalprod*vsu6[3];
  vsu6[34] -= scalprod*vsu6[4];
  vsu6[35] -= scalprod*vsu6[5];


// Make the row `5' orthogonal to the row `1':

  scalprod = conj( vsu6[6] )*vsu6[30]
    + conj( vsu6[7] )*vsu6[31]
    + conj( vsu6[8] )*vsu6[32]
    + conj( vsu6[9] )*vsu6[33]
    + conj( vsu6[10] )*vsu6[34]
    + conj( vsu6[11] )*vsu6[35];

  vsu6[30] -= scalprod*vsu6[6];
  vsu6[31] -= scalprod*vsu6[7];
  vsu6[32] -= scalprod*vsu6[8];
  vsu6[33] -= scalprod*vsu6[9];
  vsu6[34] -= scalprod*vsu6[10];
  vsu6[35] -= scalprod*vsu6[11];


// Make the row `5' orthogonal to the row `2':

  scalprod = conj( vsu6[12] )*vsu6[30]
    + conj( vsu6[13] )*vsu6[31]
    + conj( vsu6[14] )*vsu6[32]
    + conj( vsu6[15] )*vsu6[33]
    + conj( vsu6[16] )*vsu6[34]
    + conj( vsu6[17] )*vsu6[35];

  vsu6[30] -= scalprod*vsu6[12];
  vsu6[31] -= scalprod*vsu6[13];
  vsu6[32] -= scalprod*vsu6[14];
  vsu6[33] -= scalprod*vsu6[15];
  vsu6[34] -= scalprod*vsu6[16];
  vsu6[35] -= scalprod*vsu6[17];


// Make the row `5' orthogonal to the row `3':

  scalprod = conj( vsu6[18] )*vsu6[30]
    + conj( vsu6[19] )*vsu6[31]
    + conj( vsu6[20] )*vsu6[32]
    + conj( vsu6[21] )*vsu6[33]
    + conj( vsu6[22] )*vsu6[34]
    + conj( vsu6[23] )*vsu6[35];

  vsu6[30] -= scalprod*vsu6[18];
  vsu6[31] -= scalprod*vsu6[19];
  vsu6[32] -= scalprod*vsu6[20];
  vsu6[33] -= scalprod*vsu6[21];
  vsu6[34] -= scalprod*vsu6[22];
  vsu6[35] -= scalprod*vsu6[23];


// Make the row `5' orthogonal to the row `4':

  scalprod = conj( vsu6[24] )*vsu6[30]
    + conj( vsu6[25] )*vsu6[31]
    + conj( vsu6[26] )*vsu6[32]
    + conj( vsu6[27] )*vsu6[33]
    + conj( vsu6[28] )*vsu6[34]
    + conj( vsu6[29] )*vsu6[35];

  vsu6[30] -= scalprod*vsu6[24];
  vsu6[31] -= scalprod*vsu6[25];
  vsu6[32] -= scalprod*vsu6[26];
  vsu6[33] -= scalprod*vsu6[27];
  vsu6[34] -= scalprod*vsu6[28];
  vsu6[35] -= scalprod*vsu6[29];


// Normalize the row `5'

  xn = 1./sqrt( norm(vsu6[30])
    + norm(vsu6[31])
    + norm(vsu6[32])
    + norm(vsu6[33])
    + norm(vsu6[34])
    + norm(vsu6[35]) );

    vsu6[30]*=xn;
    vsu6[31]*=xn;
    vsu6[32]*=xn;
    vsu6[33]*=xn;
    vsu6[34]*=xn;
    vsu6[35]*=xn;


// Impose unimodularity by arranging the phase of the last row
// to compensate for the phase of the determinant

  for (i=0; i<Ncol; i++)
  for (j=0; j<Ncol; j++) {
    AT[2*(j+Ncol*i)]=real(vsu6[j*Ncol+i]);
    AT[2*(j+Ncol*i)+1]=imag(vsu6[j*Ncol+i]);
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

  vsu6[30] *= determinant_conjugate_phase;
  vsu6[31] *= determinant_conjugate_phase;
  vsu6[32] *= determinant_conjugate_phase;
  vsu6[33] *= determinant_conjugate_phase;
  vsu6[34] *= determinant_conjugate_phase;
  vsu6[35] *= determinant_conjugate_phase;


}

#endif
