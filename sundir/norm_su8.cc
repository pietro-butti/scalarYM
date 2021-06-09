#ifndef __norm_su5_h__
#define __norm_su5_h__

void norm_su8(dc *vsu8 ) {

// Takes the first 7 rows of a 8x8 complex matrix,
// and makes an SU(8) matrix out of them.

  dc determinant_conjugate_phase=dc(1.,0.);
  struct complex_double b[Ncol], DUMMY[1][1], WORK[2*Ncol];
  double AT[2*Ncolsquare];
  int i, j, ok, c1, c2, c3;
  char c4;

  dc scalprod=dc(0.,0.);
  double xn;

// Normalize the row `0' of vsu8

  xn = 1./sqrt( norm(vsu8[0])
    + norm(vsu8[1])
    + norm(vsu8[2])
    + norm(vsu8[3])
    + norm(vsu8[4])
    + norm(vsu8[5])
    + norm(vsu8[6])
    + norm(vsu8[7]) );

  vsu8[0]*=xn;
  vsu8[1]*=xn;
  vsu8[2]*=xn;
  vsu8[3]*=xn;
  vsu8[4]*=xn;
  vsu8[5]*=xn;
  vsu8[6]*=xn;
  vsu8[7]*=xn;

// Make the row `1' orthogonal to the row `0':

  scalprod = conj( vsu8[0] )*vsu8[8]
    + conj( vsu8[1] )*vsu8[9]
    + conj( vsu8[2] )*vsu8[10]
    + conj( vsu8[3] )*vsu8[11]
    + conj( vsu8[4] )*vsu8[12]
    + conj( vsu8[5] )*vsu8[13]
    + conj( vsu8[6] )*vsu8[14]
    + conj( vsu8[7] )*vsu8[15];

  vsu8[8] -= scalprod*vsu8[0];
  vsu8[9] -= scalprod*vsu8[1];
  vsu8[10] -= scalprod*vsu8[2];
  vsu8[11] -= scalprod*vsu8[3];
  vsu8[12] -= scalprod*vsu8[4];
  vsu8[13] -= scalprod*vsu8[5];
  vsu8[14] -= scalprod*vsu8[6];
  vsu8[15] -= scalprod*vsu8[7];


// Normalize the row `1'

  xn = 1./sqrt( norm(vsu8[8])
    + norm(vsu8[9])
    + norm(vsu8[10])
    + norm(vsu8[11])
    + norm(vsu8[12])
    + norm(vsu8[13])
    + norm(vsu8[14])
    + norm(vsu8[15]) );

    vsu8[8]*=xn;
    vsu8[9]*=xn;
    vsu8[10]*=xn;
    vsu8[11]*=xn;
    vsu8[12]*=xn;
    vsu8[13]*=xn;
    vsu8[14]*=xn;
    vsu8[15]*=xn;


// Make the row `2' orthogonal to the row `0':

  scalprod = conj( vsu8[0] )*vsu8[16]
    + conj( vsu8[1] )*vsu8[17]
    + conj( vsu8[2] )*vsu8[18]
    + conj( vsu8[3] )*vsu8[19]
    + conj( vsu8[4] )*vsu8[20]
    + conj( vsu8[5] )*vsu8[21]
    + conj( vsu8[6] )*vsu8[22]
    + conj( vsu8[7] )*vsu8[23];

  vsu8[16] -= scalprod*vsu8[0];
  vsu8[17] -= scalprod*vsu8[1];
  vsu8[18] -= scalprod*vsu8[2];
  vsu8[19] -= scalprod*vsu8[3];
  vsu8[20] -= scalprod*vsu8[4];
  vsu8[21] -= scalprod*vsu8[5];
  vsu8[22] -= scalprod*vsu8[6];
  vsu8[23] -= scalprod*vsu8[7];


// Make the row `2' orthogonal to the row `1':

  scalprod = conj( vsu8[8] )*vsu8[16]
    + conj( vsu8[9] )*vsu8[17]
    + conj( vsu8[10] )*vsu8[18]
    + conj( vsu8[11] )*vsu8[19]
    + conj( vsu8[12] )*vsu8[20]
    + conj( vsu8[13] )*vsu8[21]
    + conj( vsu8[14] )*vsu8[22]
    + conj( vsu8[15] )*vsu8[23];

  vsu8[16] -= scalprod*vsu8[8];
  vsu8[17] -= scalprod*vsu8[9];
  vsu8[18] -= scalprod*vsu8[10];
  vsu8[19] -= scalprod*vsu8[11];
  vsu8[20] -= scalprod*vsu8[12];
  vsu8[21] -= scalprod*vsu8[13];
  vsu8[22] -= scalprod*vsu8[14];
  vsu8[23] -= scalprod*vsu8[15];


// Normalize the row `2'

  xn = 1./sqrt( norm(vsu8[16])
    + norm(vsu8[17])
    + norm(vsu8[18])
    + norm(vsu8[19])
    + norm(vsu8[20])
    + norm(vsu8[21])
    + norm(vsu8[22])
    + norm(vsu8[23]) );

    vsu8[16]*=xn;
    vsu8[17]*=xn;
    vsu8[18]*=xn;
    vsu8[19]*=xn;
    vsu8[20]*=xn;
    vsu8[21]*=xn;
    vsu8[22]*=xn;
    vsu8[23]*=xn;


// Make the row `3' orthogonal to the row `0':

  scalprod = conj( vsu8[0] )*vsu8[24]
    + conj( vsu8[1] )*vsu8[25]
    + conj( vsu8[2] )*vsu8[26]
    + conj( vsu8[3] )*vsu8[27]
    + conj( vsu8[4] )*vsu8[28]
    + conj( vsu8[5] )*vsu8[29]
    + conj( vsu8[6] )*vsu8[30]
    + conj( vsu8[7] )*vsu8[31];

  vsu8[24] -= scalprod*vsu8[0];
  vsu8[25] -= scalprod*vsu8[1];
  vsu8[26] -= scalprod*vsu8[2];
  vsu8[27] -= scalprod*vsu8[3];
  vsu8[28] -= scalprod*vsu8[4];
  vsu8[29] -= scalprod*vsu8[5];
  vsu8[30] -= scalprod*vsu8[6];
  vsu8[31] -= scalprod*vsu8[7];


// Make the row `3' orthogonal to the row `1':

  scalprod = conj( vsu8[8] )*vsu8[24]
    + conj( vsu8[9] )*vsu8[25]
    + conj( vsu8[10] )*vsu8[26]
    + conj( vsu8[11] )*vsu8[27]
    + conj( vsu8[12] )*vsu8[28]
    + conj( vsu8[13] )*vsu8[29]
    + conj( vsu8[14] )*vsu8[30]
    + conj( vsu8[15] )*vsu8[31];

  vsu8[24] -= scalprod*vsu8[8];
  vsu8[25] -= scalprod*vsu8[9];
  vsu8[26] -= scalprod*vsu8[10];
  vsu8[27] -= scalprod*vsu8[11];
  vsu8[28] -= scalprod*vsu8[12];
  vsu8[29] -= scalprod*vsu8[13];
  vsu8[30] -= scalprod*vsu8[14];
  vsu8[31] -= scalprod*vsu8[15];


// Make the row `3' orthogonal to the row `2':

  scalprod = conj( vsu8[16] )*vsu8[24]
    + conj( vsu8[17] )*vsu8[25]
    + conj( vsu8[18] )*vsu8[26]
    + conj( vsu8[19] )*vsu8[27]
    + conj( vsu8[20] )*vsu8[28]
    + conj( vsu8[21] )*vsu8[29]
    + conj( vsu8[22] )*vsu8[30]
    + conj( vsu8[23] )*vsu8[31];

  vsu8[24] -= scalprod*vsu8[16];
  vsu8[25] -= scalprod*vsu8[17];
  vsu8[26] -= scalprod*vsu8[18];
  vsu8[27] -= scalprod*vsu8[19];
  vsu8[28] -= scalprod*vsu8[20];
  vsu8[29] -= scalprod*vsu8[21];
  vsu8[30] -= scalprod*vsu8[22];
  vsu8[31] -= scalprod*vsu8[23];


// Normalize the row `3'

  xn = 1./sqrt( norm(vsu8[24])
    + norm(vsu8[25])
    + norm(vsu8[26])
    + norm(vsu8[27])
    + norm(vsu8[28])
    + norm(vsu8[29])
    + norm(vsu8[30])
    + norm(vsu8[31]) );

    vsu8[24]*=xn;
    vsu8[25]*=xn;
    vsu8[26]*=xn;
    vsu8[27]*=xn;
    vsu8[28]*=xn;
    vsu8[29]*=xn;
    vsu8[30]*=xn;
    vsu8[31]*=xn;


// Make the row `4' orthogonal to the row `0':

  scalprod = conj( vsu8[0] )*vsu8[32]
    + conj( vsu8[1] )*vsu8[33]
    + conj( vsu8[2] )*vsu8[34]
    + conj( vsu8[3] )*vsu8[35]
    + conj( vsu8[4] )*vsu8[36]
    + conj( vsu8[5] )*vsu8[37]
    + conj( vsu8[6] )*vsu8[38]
    + conj( vsu8[7] )*vsu8[39];

  vsu8[32] -= scalprod*vsu8[0];
  vsu8[33] -= scalprod*vsu8[1];
  vsu8[34] -= scalprod*vsu8[2];
  vsu8[35] -= scalprod*vsu8[3];
  vsu8[36] -= scalprod*vsu8[4];
  vsu8[37] -= scalprod*vsu8[5];
  vsu8[38] -= scalprod*vsu8[6];
  vsu8[39] -= scalprod*vsu8[7];


// Make the row `4' orthogonal to the row `1':

  scalprod = conj( vsu8[8] )*vsu8[32]
    + conj( vsu8[9] )*vsu8[33]
    + conj( vsu8[10] )*vsu8[34]
    + conj( vsu8[11] )*vsu8[35]
    + conj( vsu8[12] )*vsu8[36]
    + conj( vsu8[13] )*vsu8[37]
    + conj( vsu8[14] )*vsu8[38]
    + conj( vsu8[15] )*vsu8[39];

  vsu8[32] -= scalprod*vsu8[8];
  vsu8[33] -= scalprod*vsu8[9];
  vsu8[34] -= scalprod*vsu8[10];
  vsu8[35] -= scalprod*vsu8[11];
  vsu8[36] -= scalprod*vsu8[12];
  vsu8[37] -= scalprod*vsu8[13];
  vsu8[38] -= scalprod*vsu8[14];
  vsu8[39] -= scalprod*vsu8[15];


// Make the row `4' orthogonal to the row `2':

  scalprod = conj( vsu8[16] )*vsu8[32]
    + conj( vsu8[17] )*vsu8[33]
    + conj( vsu8[18] )*vsu8[34]
    + conj( vsu8[19] )*vsu8[35]
    + conj( vsu8[20] )*vsu8[36]
    + conj( vsu8[21] )*vsu8[37]
    + conj( vsu8[22] )*vsu8[38]
    + conj( vsu8[23] )*vsu8[39];

  vsu8[32] -= scalprod*vsu8[16];
  vsu8[33] -= scalprod*vsu8[17];
  vsu8[34] -= scalprod*vsu8[18];
  vsu8[35] -= scalprod*vsu8[19];
  vsu8[36] -= scalprod*vsu8[20];
  vsu8[37] -= scalprod*vsu8[21];
  vsu8[38] -= scalprod*vsu8[22];
  vsu8[39] -= scalprod*vsu8[23];


// Make the row `4' orthogonal to the row `3':

  scalprod = conj( vsu8[24] )*vsu8[32]
    + conj( vsu8[25] )*vsu8[33]
    + conj( vsu8[26] )*vsu8[34]
    + conj( vsu8[27] )*vsu8[35]
    + conj( vsu8[28] )*vsu8[36]
    + conj( vsu8[29] )*vsu8[37]
    + conj( vsu8[30] )*vsu8[38]
    + conj( vsu8[31] )*vsu8[39];

  vsu8[32] -= scalprod*vsu8[24];
  vsu8[33] -= scalprod*vsu8[25];
  vsu8[34] -= scalprod*vsu8[26];
  vsu8[35] -= scalprod*vsu8[27];
  vsu8[36] -= scalprod*vsu8[28];
  vsu8[37] -= scalprod*vsu8[29];
  vsu8[38] -= scalprod*vsu8[30];
  vsu8[39] -= scalprod*vsu8[31];


// Normalize the row `4'

  xn = 1./sqrt( norm(vsu8[32])
    + norm(vsu8[33])
    + norm(vsu8[34])
    + norm(vsu8[35])
    + norm(vsu8[36])
    + norm(vsu8[37])
    + norm(vsu8[38])
    + norm(vsu8[39]) );

    vsu8[32]*=xn;
    vsu8[33]*=xn;
    vsu8[34]*=xn;
    vsu8[35]*=xn;
    vsu8[36]*=xn;
    vsu8[37]*=xn;
    vsu8[38]*=xn;
    vsu8[39]*=xn;


// Make the row `5' orthogonal to the row `0':

  scalprod = conj( vsu8[0] )*vsu8[40]
    + conj( vsu8[1] )*vsu8[41]
    + conj( vsu8[2] )*vsu8[42]
    + conj( vsu8[3] )*vsu8[43]
    + conj( vsu8[4] )*vsu8[44]
    + conj( vsu8[5] )*vsu8[45]
    + conj( vsu8[6] )*vsu8[46]
    + conj( vsu8[7] )*vsu8[47];

  vsu8[40] -= scalprod*vsu8[0];
  vsu8[41] -= scalprod*vsu8[1];
  vsu8[42] -= scalprod*vsu8[2];
  vsu8[43] -= scalprod*vsu8[3];
  vsu8[44] -= scalprod*vsu8[4];
  vsu8[45] -= scalprod*vsu8[5];
  vsu8[46] -= scalprod*vsu8[6];
  vsu8[47] -= scalprod*vsu8[7];


// Make the row `5' orthogonal to the row `1':

  scalprod = conj( vsu8[8] )*vsu8[40]
    + conj( vsu8[9] )*vsu8[41]
    + conj( vsu8[10] )*vsu8[42]
    + conj( vsu8[11] )*vsu8[43]
    + conj( vsu8[12] )*vsu8[44]
    + conj( vsu8[13] )*vsu8[45]
    + conj( vsu8[14] )*vsu8[46]
    + conj( vsu8[15] )*vsu8[47];

  vsu8[40] -= scalprod*vsu8[8];
  vsu8[41] -= scalprod*vsu8[9];
  vsu8[42] -= scalprod*vsu8[10];
  vsu8[43] -= scalprod*vsu8[11];
  vsu8[44] -= scalprod*vsu8[12];
  vsu8[45] -= scalprod*vsu8[13];
  vsu8[46] -= scalprod*vsu8[14];
  vsu8[47] -= scalprod*vsu8[15];


// Make the row `5' orthogonal to the row `2':

  scalprod = conj( vsu8[16] )*vsu8[40]
    + conj( vsu8[17] )*vsu8[41]
    + conj( vsu8[18] )*vsu8[42]
    + conj( vsu8[19] )*vsu8[43]
    + conj( vsu8[20] )*vsu8[44]
    + conj( vsu8[21] )*vsu8[45]
    + conj( vsu8[22] )*vsu8[46]
    + conj( vsu8[23] )*vsu8[47];

  vsu8[40] -= scalprod*vsu8[16];
  vsu8[41] -= scalprod*vsu8[17];
  vsu8[42] -= scalprod*vsu8[18];
  vsu8[43] -= scalprod*vsu8[19];
  vsu8[44] -= scalprod*vsu8[20];
  vsu8[45] -= scalprod*vsu8[21];
  vsu8[46] -= scalprod*vsu8[22];
  vsu8[47] -= scalprod*vsu8[23];


// Make the row `5' orthogonal to the row `3':

  scalprod = conj( vsu8[24] )*vsu8[40]
    + conj( vsu8[25] )*vsu8[41]
    + conj( vsu8[26] )*vsu8[42]
    + conj( vsu8[27] )*vsu8[43]
    + conj( vsu8[28] )*vsu8[44]
    + conj( vsu8[29] )*vsu8[45]
    + conj( vsu8[30] )*vsu8[46]
    + conj( vsu8[31] )*vsu8[47];

  vsu8[40] -= scalprod*vsu8[24];
  vsu8[41] -= scalprod*vsu8[25];
  vsu8[42] -= scalprod*vsu8[26];
  vsu8[43] -= scalprod*vsu8[27];
  vsu8[44] -= scalprod*vsu8[28];
  vsu8[45] -= scalprod*vsu8[29];
  vsu8[46] -= scalprod*vsu8[30];
  vsu8[47] -= scalprod*vsu8[31];


// Make the row `5' orthogonal to the row `4':

  scalprod = conj( vsu8[32] )*vsu8[40]
    + conj( vsu8[33] )*vsu8[41]
    + conj( vsu8[34] )*vsu8[42]
    + conj( vsu8[35] )*vsu8[43]
    + conj( vsu8[36] )*vsu8[44]
    + conj( vsu8[37] )*vsu8[45]
    + conj( vsu8[38] )*vsu8[46]
    + conj( vsu8[39] )*vsu8[47];

  vsu8[40] -= scalprod*vsu8[32];
  vsu8[41] -= scalprod*vsu8[33];
  vsu8[42] -= scalprod*vsu8[34];
  vsu8[43] -= scalprod*vsu8[35];
  vsu8[44] -= scalprod*vsu8[36];
  vsu8[45] -= scalprod*vsu8[37];
  vsu8[46] -= scalprod*vsu8[38];
  vsu8[47] -= scalprod*vsu8[39];


// Normalize the row `5'

  xn = 1./sqrt( norm(vsu8[40])
    + norm(vsu8[41])
    + norm(vsu8[42])
    + norm(vsu8[43])
    + norm(vsu8[44])
    + norm(vsu8[45])
    + norm(vsu8[46])
    + norm(vsu8[47]) );

    vsu8[40]*=xn;
    vsu8[41]*=xn;
    vsu8[42]*=xn;
    vsu8[43]*=xn;
    vsu8[44]*=xn;
    vsu8[45]*=xn;
    vsu8[46]*=xn;
    vsu8[47]*=xn;


// Make the row `6' orthogonal to the row `0':

  scalprod = conj( vsu8[0] )*vsu8[48]
    + conj( vsu8[1] )*vsu8[49]
    + conj( vsu8[2] )*vsu8[50]
    + conj( vsu8[3] )*vsu8[51]
    + conj( vsu8[4] )*vsu8[52]
    + conj( vsu8[5] )*vsu8[53]
    + conj( vsu8[6] )*vsu8[54]
    + conj( vsu8[7] )*vsu8[55];

  vsu8[48] -= scalprod*vsu8[0];
  vsu8[49] -= scalprod*vsu8[1];
  vsu8[50] -= scalprod*vsu8[2];
  vsu8[51] -= scalprod*vsu8[3];
  vsu8[52] -= scalprod*vsu8[4];
  vsu8[53] -= scalprod*vsu8[5];
  vsu8[54] -= scalprod*vsu8[6];
  vsu8[55] -= scalprod*vsu8[7];


// Make the row `6' orthogonal to the row `1':

  scalprod = conj( vsu8[8] )*vsu8[48]
    + conj( vsu8[9] )*vsu8[49]
    + conj( vsu8[10] )*vsu8[50]
    + conj( vsu8[11] )*vsu8[51]
    + conj( vsu8[12] )*vsu8[52]
    + conj( vsu8[13] )*vsu8[53]
    + conj( vsu8[14] )*vsu8[54]
    + conj( vsu8[15] )*vsu8[55];

  vsu8[48] -= scalprod*vsu8[8];
  vsu8[49] -= scalprod*vsu8[9];
  vsu8[50] -= scalprod*vsu8[10];
  vsu8[51] -= scalprod*vsu8[11];
  vsu8[52] -= scalprod*vsu8[12];
  vsu8[53] -= scalprod*vsu8[13];
  vsu8[54] -= scalprod*vsu8[14];
  vsu8[55] -= scalprod*vsu8[15];


// Make the row `6' orthogonal to the row `2':

  scalprod = conj( vsu8[16] )*vsu8[48]
    + conj( vsu8[17] )*vsu8[49]
    + conj( vsu8[18] )*vsu8[50]
    + conj( vsu8[19] )*vsu8[51]
    + conj( vsu8[20] )*vsu8[52]
    + conj( vsu8[21] )*vsu8[53]
    + conj( vsu8[22] )*vsu8[54]
    + conj( vsu8[23] )*vsu8[55];

  vsu8[48] -= scalprod*vsu8[16];
  vsu8[49] -= scalprod*vsu8[17];
  vsu8[50] -= scalprod*vsu8[18];
  vsu8[51] -= scalprod*vsu8[19];
  vsu8[52] -= scalprod*vsu8[20];
  vsu8[53] -= scalprod*vsu8[21];
  vsu8[54] -= scalprod*vsu8[22];
  vsu8[55] -= scalprod*vsu8[23];


// Make the row `6' orthogonal to the row `3':

  scalprod = conj( vsu8[24] )*vsu8[48]
    + conj( vsu8[25] )*vsu8[49]
    + conj( vsu8[26] )*vsu8[50]
    + conj( vsu8[27] )*vsu8[51]
    + conj( vsu8[28] )*vsu8[52]
    + conj( vsu8[29] )*vsu8[53]
    + conj( vsu8[30] )*vsu8[54]
    + conj( vsu8[31] )*vsu8[55];

  vsu8[48] -= scalprod*vsu8[24];
  vsu8[49] -= scalprod*vsu8[25];
  vsu8[50] -= scalprod*vsu8[26];
  vsu8[51] -= scalprod*vsu8[27];
  vsu8[52] -= scalprod*vsu8[28];
  vsu8[53] -= scalprod*vsu8[29];
  vsu8[54] -= scalprod*vsu8[30];
  vsu8[55] -= scalprod*vsu8[31];


// Make the row `6' orthogonal to the row `4':

  scalprod = conj( vsu8[32] )*vsu8[48]
    + conj( vsu8[33] )*vsu8[49]
    + conj( vsu8[34] )*vsu8[50]
    + conj( vsu8[35] )*vsu8[51]
    + conj( vsu8[36] )*vsu8[52]
    + conj( vsu8[37] )*vsu8[53]
    + conj( vsu8[38] )*vsu8[54]
    + conj( vsu8[39] )*vsu8[55];

  vsu8[48] -= scalprod*vsu8[32];
  vsu8[49] -= scalprod*vsu8[33];
  vsu8[50] -= scalprod*vsu8[34];
  vsu8[51] -= scalprod*vsu8[35];
  vsu8[52] -= scalprod*vsu8[36];
  vsu8[53] -= scalprod*vsu8[37];
  vsu8[54] -= scalprod*vsu8[38];
  vsu8[55] -= scalprod*vsu8[39];


// Make the row `6' orthogonal to the row `5':

  scalprod = conj( vsu8[40] )*vsu8[48]
    + conj( vsu8[41] )*vsu8[49]
    + conj( vsu8[42] )*vsu8[50]
    + conj( vsu8[43] )*vsu8[51]
    + conj( vsu8[44] )*vsu8[52]
    + conj( vsu8[45] )*vsu8[53]
    + conj( vsu8[46] )*vsu8[54]
    + conj( vsu8[47] )*vsu8[55];

  vsu8[48] -= scalprod*vsu8[40];
  vsu8[49] -= scalprod*vsu8[41];
  vsu8[50] -= scalprod*vsu8[42];
  vsu8[51] -= scalprod*vsu8[43];
  vsu8[52] -= scalprod*vsu8[44];
  vsu8[53] -= scalprod*vsu8[45];
  vsu8[54] -= scalprod*vsu8[46];
  vsu8[55] -= scalprod*vsu8[47];


// Normalize the row `6'

  xn = 1./sqrt( norm(vsu8[48])
    + norm(vsu8[49])
    + norm(vsu8[50])
    + norm(vsu8[51])
    + norm(vsu8[52])
    + norm(vsu8[53])
    + norm(vsu8[54])
    + norm(vsu8[55]) );

    vsu8[48]*=xn;
    vsu8[49]*=xn;
    vsu8[50]*=xn;
    vsu8[51]*=xn;
    vsu8[52]*=xn;
    vsu8[53]*=xn;
    vsu8[54]*=xn;
    vsu8[55]*=xn;


// Make the row `7' orthogonal to the row `0':

  scalprod = conj( vsu8[0] )*vsu8[56]
    + conj( vsu8[1] )*vsu8[57]
    + conj( vsu8[2] )*vsu8[58]
    + conj( vsu8[3] )*vsu8[59]
    + conj( vsu8[4] )*vsu8[60]
    + conj( vsu8[5] )*vsu8[61]
    + conj( vsu8[6] )*vsu8[62]
    + conj( vsu8[7] )*vsu8[63];

  vsu8[56] -= scalprod*vsu8[0];
  vsu8[57] -= scalprod*vsu8[1];
  vsu8[58] -= scalprod*vsu8[2];
  vsu8[59] -= scalprod*vsu8[3];
  vsu8[60] -= scalprod*vsu8[4];
  vsu8[61] -= scalprod*vsu8[5];
  vsu8[62] -= scalprod*vsu8[6];
  vsu8[63] -= scalprod*vsu8[7];


// Make the row `7' orthogonal to the row `1':

  scalprod = conj( vsu8[8] )*vsu8[56]
    + conj( vsu8[9] )*vsu8[57]
    + conj( vsu8[10] )*vsu8[58]
    + conj( vsu8[11] )*vsu8[59]
    + conj( vsu8[12] )*vsu8[60]
    + conj( vsu8[13] )*vsu8[61]
    + conj( vsu8[14] )*vsu8[62]
    + conj( vsu8[15] )*vsu8[63];

  vsu8[56] -= scalprod*vsu8[8];
  vsu8[57] -= scalprod*vsu8[9];
  vsu8[58] -= scalprod*vsu8[10];
  vsu8[59] -= scalprod*vsu8[11];
  vsu8[60] -= scalprod*vsu8[12];
  vsu8[61] -= scalprod*vsu8[13];
  vsu8[62] -= scalprod*vsu8[14];
  vsu8[63] -= scalprod*vsu8[15];


// Make the row `7' orthogonal to the row `2':

  scalprod = conj( vsu8[16] )*vsu8[56]
    + conj( vsu8[17] )*vsu8[57]
    + conj( vsu8[18] )*vsu8[58]
    + conj( vsu8[19] )*vsu8[59]
    + conj( vsu8[20] )*vsu8[60]
    + conj( vsu8[21] )*vsu8[61]
    + conj( vsu8[22] )*vsu8[62]
    + conj( vsu8[23] )*vsu8[63];

  vsu8[56] -= scalprod*vsu8[16];
  vsu8[57] -= scalprod*vsu8[17];
  vsu8[58] -= scalprod*vsu8[18];
  vsu8[59] -= scalprod*vsu8[19];
  vsu8[60] -= scalprod*vsu8[20];
  vsu8[61] -= scalprod*vsu8[21];
  vsu8[62] -= scalprod*vsu8[22];
  vsu8[63] -= scalprod*vsu8[23];


// Make the row `7' orthogonal to the row `3':

  scalprod = conj( vsu8[24] )*vsu8[56]
    + conj( vsu8[25] )*vsu8[57]
    + conj( vsu8[26] )*vsu8[58]
    + conj( vsu8[27] )*vsu8[59]
    + conj( vsu8[28] )*vsu8[60]
    + conj( vsu8[29] )*vsu8[61]
    + conj( vsu8[30] )*vsu8[62]
    + conj( vsu8[31] )*vsu8[63];

  vsu8[56] -= scalprod*vsu8[24];
  vsu8[57] -= scalprod*vsu8[25];
  vsu8[58] -= scalprod*vsu8[26];
  vsu8[59] -= scalprod*vsu8[27];
  vsu8[60] -= scalprod*vsu8[28];
  vsu8[61] -= scalprod*vsu8[29];
  vsu8[62] -= scalprod*vsu8[30];
  vsu8[63] -= scalprod*vsu8[31];


// Make the row `7' orthogonal to the row `4':

  scalprod = conj( vsu8[32] )*vsu8[56]
    + conj( vsu8[33] )*vsu8[57]
    + conj( vsu8[34] )*vsu8[58]
    + conj( vsu8[35] )*vsu8[59]
    + conj( vsu8[36] )*vsu8[60]
    + conj( vsu8[37] )*vsu8[61]
    + conj( vsu8[38] )*vsu8[62]
    + conj( vsu8[39] )*vsu8[63];

  vsu8[56] -= scalprod*vsu8[32];
  vsu8[57] -= scalprod*vsu8[33];
  vsu8[58] -= scalprod*vsu8[34];
  vsu8[59] -= scalprod*vsu8[35];
  vsu8[60] -= scalprod*vsu8[36];
  vsu8[61] -= scalprod*vsu8[37];
  vsu8[62] -= scalprod*vsu8[38];
  vsu8[63] -= scalprod*vsu8[39];


// Make the row `7' orthogonal to the row `5':

  scalprod = conj( vsu8[40] )*vsu8[56]
    + conj( vsu8[41] )*vsu8[57]
    + conj( vsu8[42] )*vsu8[58]
    + conj( vsu8[43] )*vsu8[59]
    + conj( vsu8[44] )*vsu8[60]
    + conj( vsu8[45] )*vsu8[61]
    + conj( vsu8[46] )*vsu8[62]
    + conj( vsu8[47] )*vsu8[63];

  vsu8[56] -= scalprod*vsu8[40];
  vsu8[57] -= scalprod*vsu8[41];
  vsu8[58] -= scalprod*vsu8[42];
  vsu8[59] -= scalprod*vsu8[43];
  vsu8[60] -= scalprod*vsu8[44];
  vsu8[61] -= scalprod*vsu8[45];
  vsu8[62] -= scalprod*vsu8[46];
  vsu8[63] -= scalprod*vsu8[47];


// Make the row `7' orthogonal to the row `6':

  scalprod = conj( vsu8[48] )*vsu8[56]
    + conj( vsu8[49] )*vsu8[57]
    + conj( vsu8[50] )*vsu8[58]
    + conj( vsu8[51] )*vsu8[59]
    + conj( vsu8[52] )*vsu8[60]
    + conj( vsu8[53] )*vsu8[61]
    + conj( vsu8[54] )*vsu8[62]
    + conj( vsu8[55] )*vsu8[63];

  vsu8[56] -= scalprod*vsu8[48];
  vsu8[57] -= scalprod*vsu8[49];
  vsu8[58] -= scalprod*vsu8[50];
  vsu8[59] -= scalprod*vsu8[51];
  vsu8[60] -= scalprod*vsu8[52];
  vsu8[61] -= scalprod*vsu8[53];
  vsu8[62] -= scalprod*vsu8[54];
  vsu8[63] -= scalprod*vsu8[55];


// Normalize the row `7'

  xn = 1./sqrt( norm(vsu8[56])
    + norm(vsu8[57])
    + norm(vsu8[58])
    + norm(vsu8[59])
    + norm(vsu8[60])
    + norm(vsu8[61])
    + norm(vsu8[62])
    + norm(vsu8[63]) );

    vsu8[56]*=xn;
    vsu8[57]*=xn;
    vsu8[58]*=xn;
    vsu8[59]*=xn;
    vsu8[60]*=xn;
    vsu8[61]*=xn;
    vsu8[62]*=xn;
    vsu8[63]*=xn;


// Impose unimodularity by arranging the phase of the last row
// to compensate for the phase of the determinant

  for (i=0; i<Ncol; i++)
  for (j=0; j<Ncol; j++) {
    AT[2*(j+Ncol*i)]=real(vsu8[j*Ncol+i]);
    AT[2*(j+Ncol*i)+1]=imag(vsu8[j*Ncol+i]);
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

  vsu8[56] *= determinant_conjugate_phase;
  vsu8[57] *= determinant_conjugate_phase;
  vsu8[58] *= determinant_conjugate_phase;
  vsu8[59] *= determinant_conjugate_phase;
  vsu8[60] *= determinant_conjugate_phase;
  vsu8[61] *= determinant_conjugate_phase;
  vsu8[62] *= determinant_conjugate_phase;
  vsu8[63] *= determinant_conjugate_phase;


}

#endif
