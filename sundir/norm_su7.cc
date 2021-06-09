#ifndef __norm_su5_h__
#define __norm_su5_h__

void norm_su7(dc *vsu7 ) {

// Takes the first 6 rows of a 7x7 complex matrix,
// and makes an SU(7) matrix out of them.

  dc determinant_conjugate_phase=dc(1.,0.);
  struct complex_double b[Ncol], DUMMY[1][1], WORK[2*Ncol];
  double AT[2*Ncolsquare];
  int i, j, ok, c1, c2, c3;
  char c4;

  dc scalprod=dc(0.,0.);
  double xn;

// Normalize the row `0' of vsu7

  xn = 1./sqrt( norm(vsu7[0])
    + norm(vsu7[1])
    + norm(vsu7[2])
    + norm(vsu7[3])
    + norm(vsu7[4])
    + norm(vsu7[5])
    + norm(vsu7[6]) );

  vsu7[0]*=xn;
  vsu7[1]*=xn;
  vsu7[2]*=xn;
  vsu7[3]*=xn;
  vsu7[4]*=xn;
  vsu7[5]*=xn;
  vsu7[6]*=xn;

// Make the row `1' orthogonal to the row `0':

  scalprod = conj( vsu7[0] )*vsu7[7]
    + conj( vsu7[1] )*vsu7[8]
    + conj( vsu7[2] )*vsu7[9]
    + conj( vsu7[3] )*vsu7[10]
    + conj( vsu7[4] )*vsu7[11]
    + conj( vsu7[5] )*vsu7[12]
    + conj( vsu7[6] )*vsu7[13];

  vsu7[7] -= scalprod*vsu7[0];
  vsu7[8] -= scalprod*vsu7[1];
  vsu7[9] -= scalprod*vsu7[2];
  vsu7[10] -= scalprod*vsu7[3];
  vsu7[11] -= scalprod*vsu7[4];
  vsu7[12] -= scalprod*vsu7[5];
  vsu7[13] -= scalprod*vsu7[6];


// Normalize the row `1'

  xn = 1./sqrt( norm(vsu7[7])
    + norm(vsu7[8])
    + norm(vsu7[9])
    + norm(vsu7[10])
    + norm(vsu7[11])
    + norm(vsu7[12])
    + norm(vsu7[13]) );

    vsu7[7]*=xn;
    vsu7[8]*=xn;
    vsu7[9]*=xn;
    vsu7[10]*=xn;
    vsu7[11]*=xn;
    vsu7[12]*=xn;
    vsu7[13]*=xn;


// Make the row `2' orthogonal to the row `0':

  scalprod = conj( vsu7[0] )*vsu7[14]
    + conj( vsu7[1] )*vsu7[15]
    + conj( vsu7[2] )*vsu7[16]
    + conj( vsu7[3] )*vsu7[17]
    + conj( vsu7[4] )*vsu7[18]
    + conj( vsu7[5] )*vsu7[19]
    + conj( vsu7[6] )*vsu7[20];

  vsu7[14] -= scalprod*vsu7[0];
  vsu7[15] -= scalprod*vsu7[1];
  vsu7[16] -= scalprod*vsu7[2];
  vsu7[17] -= scalprod*vsu7[3];
  vsu7[18] -= scalprod*vsu7[4];
  vsu7[19] -= scalprod*vsu7[5];
  vsu7[20] -= scalprod*vsu7[6];


// Make the row `2' orthogonal to the row `1':

  scalprod = conj( vsu7[7] )*vsu7[14]
    + conj( vsu7[8] )*vsu7[15]
    + conj( vsu7[9] )*vsu7[16]
    + conj( vsu7[10] )*vsu7[17]
    + conj( vsu7[11] )*vsu7[18]
    + conj( vsu7[12] )*vsu7[19]
    + conj( vsu7[13] )*vsu7[20];

  vsu7[14] -= scalprod*vsu7[7];
  vsu7[15] -= scalprod*vsu7[8];
  vsu7[16] -= scalprod*vsu7[9];
  vsu7[17] -= scalprod*vsu7[10];
  vsu7[18] -= scalprod*vsu7[11];
  vsu7[19] -= scalprod*vsu7[12];
  vsu7[20] -= scalprod*vsu7[13];


// Normalize the row `2'

  xn = 1./sqrt( norm(vsu7[14])
    + norm(vsu7[15])
    + norm(vsu7[16])
    + norm(vsu7[17])
    + norm(vsu7[18])
    + norm(vsu7[19])
    + norm(vsu7[20]) );

    vsu7[14]*=xn;
    vsu7[15]*=xn;
    vsu7[16]*=xn;
    vsu7[17]*=xn;
    vsu7[18]*=xn;
    vsu7[19]*=xn;
    vsu7[20]*=xn;


// Make the row `3' orthogonal to the row `0':

  scalprod = conj( vsu7[0] )*vsu7[21]
    + conj( vsu7[1] )*vsu7[22]
    + conj( vsu7[2] )*vsu7[23]
    + conj( vsu7[3] )*vsu7[24]
    + conj( vsu7[4] )*vsu7[25]
    + conj( vsu7[5] )*vsu7[26]
    + conj( vsu7[6] )*vsu7[27];

  vsu7[21] -= scalprod*vsu7[0];
  vsu7[22] -= scalprod*vsu7[1];
  vsu7[23] -= scalprod*vsu7[2];
  vsu7[24] -= scalprod*vsu7[3];
  vsu7[25] -= scalprod*vsu7[4];
  vsu7[26] -= scalprod*vsu7[5];
  vsu7[27] -= scalprod*vsu7[6];


// Make the row `3' orthogonal to the row `1':

  scalprod = conj( vsu7[7] )*vsu7[21]
    + conj( vsu7[8] )*vsu7[22]
    + conj( vsu7[9] )*vsu7[23]
    + conj( vsu7[10] )*vsu7[24]
    + conj( vsu7[11] )*vsu7[25]
    + conj( vsu7[12] )*vsu7[26]
    + conj( vsu7[13] )*vsu7[27];

  vsu7[21] -= scalprod*vsu7[7];
  vsu7[22] -= scalprod*vsu7[8];
  vsu7[23] -= scalprod*vsu7[9];
  vsu7[24] -= scalprod*vsu7[10];
  vsu7[25] -= scalprod*vsu7[11];
  vsu7[26] -= scalprod*vsu7[12];
  vsu7[27] -= scalprod*vsu7[13];


// Make the row `3' orthogonal to the row `2':

  scalprod = conj( vsu7[14] )*vsu7[21]
    + conj( vsu7[15] )*vsu7[22]
    + conj( vsu7[16] )*vsu7[23]
    + conj( vsu7[17] )*vsu7[24]
    + conj( vsu7[18] )*vsu7[25]
    + conj( vsu7[19] )*vsu7[26]
    + conj( vsu7[20] )*vsu7[27];

  vsu7[21] -= scalprod*vsu7[14];
  vsu7[22] -= scalprod*vsu7[15];
  vsu7[23] -= scalprod*vsu7[16];
  vsu7[24] -= scalprod*vsu7[17];
  vsu7[25] -= scalprod*vsu7[18];
  vsu7[26] -= scalprod*vsu7[19];
  vsu7[27] -= scalprod*vsu7[20];


// Normalize the row `3'

  xn = 1./sqrt( norm(vsu7[21])
    + norm(vsu7[22])
    + norm(vsu7[23])
    + norm(vsu7[24])
    + norm(vsu7[25])
    + norm(vsu7[26])
    + norm(vsu7[27]) );

    vsu7[21]*=xn;
    vsu7[22]*=xn;
    vsu7[23]*=xn;
    vsu7[24]*=xn;
    vsu7[25]*=xn;
    vsu7[26]*=xn;
    vsu7[27]*=xn;


// Make the row `4' orthogonal to the row `0':

  scalprod = conj( vsu7[0] )*vsu7[28]
    + conj( vsu7[1] )*vsu7[29]
    + conj( vsu7[2] )*vsu7[30]
    + conj( vsu7[3] )*vsu7[31]
    + conj( vsu7[4] )*vsu7[32]
    + conj( vsu7[5] )*vsu7[33]
    + conj( vsu7[6] )*vsu7[34];

  vsu7[28] -= scalprod*vsu7[0];
  vsu7[29] -= scalprod*vsu7[1];
  vsu7[30] -= scalprod*vsu7[2];
  vsu7[31] -= scalprod*vsu7[3];
  vsu7[32] -= scalprod*vsu7[4];
  vsu7[33] -= scalprod*vsu7[5];
  vsu7[34] -= scalprod*vsu7[6];


// Make the row `4' orthogonal to the row `1':

  scalprod = conj( vsu7[7] )*vsu7[28]
    + conj( vsu7[8] )*vsu7[29]
    + conj( vsu7[9] )*vsu7[30]
    + conj( vsu7[10] )*vsu7[31]
    + conj( vsu7[11] )*vsu7[32]
    + conj( vsu7[12] )*vsu7[33]
    + conj( vsu7[13] )*vsu7[34];

  vsu7[28] -= scalprod*vsu7[7];
  vsu7[29] -= scalprod*vsu7[8];
  vsu7[30] -= scalprod*vsu7[9];
  vsu7[31] -= scalprod*vsu7[10];
  vsu7[32] -= scalprod*vsu7[11];
  vsu7[33] -= scalprod*vsu7[12];
  vsu7[34] -= scalprod*vsu7[13];


// Make the row `4' orthogonal to the row `2':

  scalprod = conj( vsu7[14] )*vsu7[28]
    + conj( vsu7[15] )*vsu7[29]
    + conj( vsu7[16] )*vsu7[30]
    + conj( vsu7[17] )*vsu7[31]
    + conj( vsu7[18] )*vsu7[32]
    + conj( vsu7[19] )*vsu7[33]
    + conj( vsu7[20] )*vsu7[34];

  vsu7[28] -= scalprod*vsu7[14];
  vsu7[29] -= scalprod*vsu7[15];
  vsu7[30] -= scalprod*vsu7[16];
  vsu7[31] -= scalprod*vsu7[17];
  vsu7[32] -= scalprod*vsu7[18];
  vsu7[33] -= scalprod*vsu7[19];
  vsu7[34] -= scalprod*vsu7[20];


// Make the row `4' orthogonal to the row `3':

  scalprod = conj( vsu7[21] )*vsu7[28]
    + conj( vsu7[22] )*vsu7[29]
    + conj( vsu7[23] )*vsu7[30]
    + conj( vsu7[24] )*vsu7[31]
    + conj( vsu7[25] )*vsu7[32]
    + conj( vsu7[26] )*vsu7[33]
    + conj( vsu7[27] )*vsu7[34];

  vsu7[28] -= scalprod*vsu7[21];
  vsu7[29] -= scalprod*vsu7[22];
  vsu7[30] -= scalprod*vsu7[23];
  vsu7[31] -= scalprod*vsu7[24];
  vsu7[32] -= scalprod*vsu7[25];
  vsu7[33] -= scalprod*vsu7[26];
  vsu7[34] -= scalprod*vsu7[27];


// Normalize the row `4'

  xn = 1./sqrt( norm(vsu7[28])
    + norm(vsu7[29])
    + norm(vsu7[30])
    + norm(vsu7[31])
    + norm(vsu7[32])
    + norm(vsu7[33])
    + norm(vsu7[34]) );

    vsu7[28]*=xn;
    vsu7[29]*=xn;
    vsu7[30]*=xn;
    vsu7[31]*=xn;
    vsu7[32]*=xn;
    vsu7[33]*=xn;
    vsu7[34]*=xn;


// Make the row `5' orthogonal to the row `0':

  scalprod = conj( vsu7[0] )*vsu7[35]
    + conj( vsu7[1] )*vsu7[36]
    + conj( vsu7[2] )*vsu7[37]
    + conj( vsu7[3] )*vsu7[38]
    + conj( vsu7[4] )*vsu7[39]
    + conj( vsu7[5] )*vsu7[40]
    + conj( vsu7[6] )*vsu7[41];

  vsu7[35] -= scalprod*vsu7[0];
  vsu7[36] -= scalprod*vsu7[1];
  vsu7[37] -= scalprod*vsu7[2];
  vsu7[38] -= scalprod*vsu7[3];
  vsu7[39] -= scalprod*vsu7[4];
  vsu7[40] -= scalprod*vsu7[5];
  vsu7[41] -= scalprod*vsu7[6];


// Make the row `5' orthogonal to the row `1':

  scalprod = conj( vsu7[7] )*vsu7[35]
    + conj( vsu7[8] )*vsu7[36]
    + conj( vsu7[9] )*vsu7[37]
    + conj( vsu7[10] )*vsu7[38]
    + conj( vsu7[11] )*vsu7[39]
    + conj( vsu7[12] )*vsu7[40]
    + conj( vsu7[13] )*vsu7[41];

  vsu7[35] -= scalprod*vsu7[7];
  vsu7[36] -= scalprod*vsu7[8];
  vsu7[37] -= scalprod*vsu7[9];
  vsu7[38] -= scalprod*vsu7[10];
  vsu7[39] -= scalprod*vsu7[11];
  vsu7[40] -= scalprod*vsu7[12];
  vsu7[41] -= scalprod*vsu7[13];


// Make the row `5' orthogonal to the row `2':

  scalprod = conj( vsu7[14] )*vsu7[35]
    + conj( vsu7[15] )*vsu7[36]
    + conj( vsu7[16] )*vsu7[37]
    + conj( vsu7[17] )*vsu7[38]
    + conj( vsu7[18] )*vsu7[39]
    + conj( vsu7[19] )*vsu7[40]
    + conj( vsu7[20] )*vsu7[41];

  vsu7[35] -= scalprod*vsu7[14];
  vsu7[36] -= scalprod*vsu7[15];
  vsu7[37] -= scalprod*vsu7[16];
  vsu7[38] -= scalprod*vsu7[17];
  vsu7[39] -= scalprod*vsu7[18];
  vsu7[40] -= scalprod*vsu7[19];
  vsu7[41] -= scalprod*vsu7[20];


// Make the row `5' orthogonal to the row `3':

  scalprod = conj( vsu7[21] )*vsu7[35]
    + conj( vsu7[22] )*vsu7[36]
    + conj( vsu7[23] )*vsu7[37]
    + conj( vsu7[24] )*vsu7[38]
    + conj( vsu7[25] )*vsu7[39]
    + conj( vsu7[26] )*vsu7[40]
    + conj( vsu7[27] )*vsu7[41];

  vsu7[35] -= scalprod*vsu7[21];
  vsu7[36] -= scalprod*vsu7[22];
  vsu7[37] -= scalprod*vsu7[23];
  vsu7[38] -= scalprod*vsu7[24];
  vsu7[39] -= scalprod*vsu7[25];
  vsu7[40] -= scalprod*vsu7[26];
  vsu7[41] -= scalprod*vsu7[27];


// Make the row `5' orthogonal to the row `4':

  scalprod = conj( vsu7[28] )*vsu7[35]
    + conj( vsu7[29] )*vsu7[36]
    + conj( vsu7[30] )*vsu7[37]
    + conj( vsu7[31] )*vsu7[38]
    + conj( vsu7[32] )*vsu7[39]
    + conj( vsu7[33] )*vsu7[40]
    + conj( vsu7[34] )*vsu7[41];

  vsu7[35] -= scalprod*vsu7[28];
  vsu7[36] -= scalprod*vsu7[29];
  vsu7[37] -= scalprod*vsu7[30];
  vsu7[38] -= scalprod*vsu7[31];
  vsu7[39] -= scalprod*vsu7[32];
  vsu7[40] -= scalprod*vsu7[33];
  vsu7[41] -= scalprod*vsu7[34];


// Normalize the row `5'

  xn = 1./sqrt( norm(vsu7[35])
    + norm(vsu7[36])
    + norm(vsu7[37])
    + norm(vsu7[38])
    + norm(vsu7[39])
    + norm(vsu7[40])
    + norm(vsu7[41]) );

    vsu7[35]*=xn;
    vsu7[36]*=xn;
    vsu7[37]*=xn;
    vsu7[38]*=xn;
    vsu7[39]*=xn;
    vsu7[40]*=xn;
    vsu7[41]*=xn;


// Make the row `6' orthogonal to the row `0':

  scalprod = conj( vsu7[0] )*vsu7[42]
    + conj( vsu7[1] )*vsu7[43]
    + conj( vsu7[2] )*vsu7[44]
    + conj( vsu7[3] )*vsu7[45]
    + conj( vsu7[4] )*vsu7[46]
    + conj( vsu7[5] )*vsu7[47]
    + conj( vsu7[6] )*vsu7[48];

  vsu7[42] -= scalprod*vsu7[0];
  vsu7[43] -= scalprod*vsu7[1];
  vsu7[44] -= scalprod*vsu7[2];
  vsu7[45] -= scalprod*vsu7[3];
  vsu7[46] -= scalprod*vsu7[4];
  vsu7[47] -= scalprod*vsu7[5];
  vsu7[48] -= scalprod*vsu7[6];


// Make the row `6' orthogonal to the row `1':

  scalprod = conj( vsu7[7] )*vsu7[42]
    + conj( vsu7[8] )*vsu7[43]
    + conj( vsu7[9] )*vsu7[44]
    + conj( vsu7[10] )*vsu7[45]
    + conj( vsu7[11] )*vsu7[46]
    + conj( vsu7[12] )*vsu7[47]
    + conj( vsu7[13] )*vsu7[48];

  vsu7[42] -= scalprod*vsu7[7];
  vsu7[43] -= scalprod*vsu7[8];
  vsu7[44] -= scalprod*vsu7[9];
  vsu7[45] -= scalprod*vsu7[10];
  vsu7[46] -= scalprod*vsu7[11];
  vsu7[47] -= scalprod*vsu7[12];
  vsu7[48] -= scalprod*vsu7[13];


// Make the row `6' orthogonal to the row `2':

  scalprod = conj( vsu7[14] )*vsu7[42]
    + conj( vsu7[15] )*vsu7[43]
    + conj( vsu7[16] )*vsu7[44]
    + conj( vsu7[17] )*vsu7[45]
    + conj( vsu7[18] )*vsu7[46]
    + conj( vsu7[19] )*vsu7[47]
    + conj( vsu7[20] )*vsu7[48];

  vsu7[42] -= scalprod*vsu7[14];
  vsu7[43] -= scalprod*vsu7[15];
  vsu7[44] -= scalprod*vsu7[16];
  vsu7[45] -= scalprod*vsu7[17];
  vsu7[46] -= scalprod*vsu7[18];
  vsu7[47] -= scalprod*vsu7[19];
  vsu7[48] -= scalprod*vsu7[20];


// Make the row `6' orthogonal to the row `3':

  scalprod = conj( vsu7[21] )*vsu7[42]
    + conj( vsu7[22] )*vsu7[43]
    + conj( vsu7[23] )*vsu7[44]
    + conj( vsu7[24] )*vsu7[45]
    + conj( vsu7[25] )*vsu7[46]
    + conj( vsu7[26] )*vsu7[47]
    + conj( vsu7[27] )*vsu7[48];

  vsu7[42] -= scalprod*vsu7[21];
  vsu7[43] -= scalprod*vsu7[22];
  vsu7[44] -= scalprod*vsu7[23];
  vsu7[45] -= scalprod*vsu7[24];
  vsu7[46] -= scalprod*vsu7[25];
  vsu7[47] -= scalprod*vsu7[26];
  vsu7[48] -= scalprod*vsu7[27];


// Make the row `6' orthogonal to the row `4':

  scalprod = conj( vsu7[28] )*vsu7[42]
    + conj( vsu7[29] )*vsu7[43]
    + conj( vsu7[30] )*vsu7[44]
    + conj( vsu7[31] )*vsu7[45]
    + conj( vsu7[32] )*vsu7[46]
    + conj( vsu7[33] )*vsu7[47]
    + conj( vsu7[34] )*vsu7[48];

  vsu7[42] -= scalprod*vsu7[28];
  vsu7[43] -= scalprod*vsu7[29];
  vsu7[44] -= scalprod*vsu7[30];
  vsu7[45] -= scalprod*vsu7[31];
  vsu7[46] -= scalprod*vsu7[32];
  vsu7[47] -= scalprod*vsu7[33];
  vsu7[48] -= scalprod*vsu7[34];


// Make the row `6' orthogonal to the row `5':

  scalprod = conj( vsu7[35] )*vsu7[42]
    + conj( vsu7[36] )*vsu7[43]
    + conj( vsu7[37] )*vsu7[44]
    + conj( vsu7[38] )*vsu7[45]
    + conj( vsu7[39] )*vsu7[46]
    + conj( vsu7[40] )*vsu7[47]
    + conj( vsu7[41] )*vsu7[48];

  vsu7[42] -= scalprod*vsu7[35];
  vsu7[43] -= scalprod*vsu7[36];
  vsu7[44] -= scalprod*vsu7[37];
  vsu7[45] -= scalprod*vsu7[38];
  vsu7[46] -= scalprod*vsu7[39];
  vsu7[47] -= scalprod*vsu7[40];
  vsu7[48] -= scalprod*vsu7[41];


// Normalize the row `6'

  xn = 1./sqrt( norm(vsu7[42])
    + norm(vsu7[43])
    + norm(vsu7[44])
    + norm(vsu7[45])
    + norm(vsu7[46])
    + norm(vsu7[47])
    + norm(vsu7[48]) );

    vsu7[42]*=xn;
    vsu7[43]*=xn;
    vsu7[44]*=xn;
    vsu7[45]*=xn;
    vsu7[46]*=xn;
    vsu7[47]*=xn;
    vsu7[48]*=xn;


// Impose unimodularity by arranging the phase of the last row
// to compensate for the phase of the determinant

  for (i=0; i<Ncol; i++)
  for (j=0; j<Ncol; j++) {
    AT[2*(j+Ncol*i)]=real(vsu7[j*Ncol+i]);
    AT[2*(j+Ncol*i)+1]=imag(vsu7[j*Ncol+i]);
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

  vsu7[42] *= determinant_conjugate_phase;
  vsu7[43] *= determinant_conjugate_phase;
  vsu7[44] *= determinant_conjugate_phase;
  vsu7[45] *= determinant_conjugate_phase;
  vsu7[46] *= determinant_conjugate_phase;
  vsu7[47] *= determinant_conjugate_phase;
  vsu7[48] *= determinant_conjugate_phase;


}

#endif
