#ifndef __norm_su5_h__
#define __norm_su5_h__

void norm_su9(dc *vsu9 ) {

// Takes the first 8 rows of a 9x9 complex matrix,
// and makes an SU(9) matrix out of them.

  dc determinant_conjugate_phase=dc(1.,0.);
  struct complex_double b[Ncol], DUMMY[1][1], WORK[2*Ncol];
  double AT[2*Ncolsquare];
  int i, j, ok, c1, c2, c3;
  char c4;

  dc scalprod=dc(0.,0.);
  double xn;

// Normalize the row `0' of vsu9

  xn = 1./sqrt( norm(vsu9[0])
    + norm(vsu9[1])
    + norm(vsu9[2])
    + norm(vsu9[3])
    + norm(vsu9[4])
    + norm(vsu9[5])
    + norm(vsu9[6])
    + norm(vsu9[7])
    + norm(vsu9[8]) );

  vsu9[0]*=xn;
  vsu9[1]*=xn;
  vsu9[2]*=xn;
  vsu9[3]*=xn;
  vsu9[4]*=xn;
  vsu9[5]*=xn;
  vsu9[6]*=xn;
  vsu9[7]*=xn;
  vsu9[8]*=xn;

// Make the row `1' orthogonal to the row `0':

  scalprod = conj( vsu9[0] )*vsu9[9]
    + conj( vsu9[1] )*vsu9[10]
    + conj( vsu9[2] )*vsu9[11]
    + conj( vsu9[3] )*vsu9[12]
    + conj( vsu9[4] )*vsu9[13]
    + conj( vsu9[5] )*vsu9[14]
    + conj( vsu9[6] )*vsu9[15]
    + conj( vsu9[7] )*vsu9[16]
    + conj( vsu9[8] )*vsu9[17];

  vsu9[9] -= scalprod*vsu9[0];
  vsu9[10] -= scalprod*vsu9[1];
  vsu9[11] -= scalprod*vsu9[2];
  vsu9[12] -= scalprod*vsu9[3];
  vsu9[13] -= scalprod*vsu9[4];
  vsu9[14] -= scalprod*vsu9[5];
  vsu9[15] -= scalprod*vsu9[6];
  vsu9[16] -= scalprod*vsu9[7];
  vsu9[17] -= scalprod*vsu9[8];


// Normalize the row `1'

  xn = 1./sqrt( norm(vsu9[9])
    + norm(vsu9[10])
    + norm(vsu9[11])
    + norm(vsu9[12])
    + norm(vsu9[13])
    + norm(vsu9[14])
    + norm(vsu9[15])
    + norm(vsu9[16])
    + norm(vsu9[17]) );

    vsu9[9]*=xn;
    vsu9[10]*=xn;
    vsu9[11]*=xn;
    vsu9[12]*=xn;
    vsu9[13]*=xn;
    vsu9[14]*=xn;
    vsu9[15]*=xn;
    vsu9[16]*=xn;
    vsu9[17]*=xn;


// Make the row `2' orthogonal to the row `0':

  scalprod = conj( vsu9[0] )*vsu9[18]
    + conj( vsu9[1] )*vsu9[19]
    + conj( vsu9[2] )*vsu9[20]
    + conj( vsu9[3] )*vsu9[21]
    + conj( vsu9[4] )*vsu9[22]
    + conj( vsu9[5] )*vsu9[23]
    + conj( vsu9[6] )*vsu9[24]
    + conj( vsu9[7] )*vsu9[25]
    + conj( vsu9[8] )*vsu9[26];

  vsu9[18] -= scalprod*vsu9[0];
  vsu9[19] -= scalprod*vsu9[1];
  vsu9[20] -= scalprod*vsu9[2];
  vsu9[21] -= scalprod*vsu9[3];
  vsu9[22] -= scalprod*vsu9[4];
  vsu9[23] -= scalprod*vsu9[5];
  vsu9[24] -= scalprod*vsu9[6];
  vsu9[25] -= scalprod*vsu9[7];
  vsu9[26] -= scalprod*vsu9[8];


// Make the row `2' orthogonal to the row `1':

  scalprod = conj( vsu9[9] )*vsu9[18]
    + conj( vsu9[10] )*vsu9[19]
    + conj( vsu9[11] )*vsu9[20]
    + conj( vsu9[12] )*vsu9[21]
    + conj( vsu9[13] )*vsu9[22]
    + conj( vsu9[14] )*vsu9[23]
    + conj( vsu9[15] )*vsu9[24]
    + conj( vsu9[16] )*vsu9[25]
    + conj( vsu9[17] )*vsu9[26];

  vsu9[18] -= scalprod*vsu9[9];
  vsu9[19] -= scalprod*vsu9[10];
  vsu9[20] -= scalprod*vsu9[11];
  vsu9[21] -= scalprod*vsu9[12];
  vsu9[22] -= scalprod*vsu9[13];
  vsu9[23] -= scalprod*vsu9[14];
  vsu9[24] -= scalprod*vsu9[15];
  vsu9[25] -= scalprod*vsu9[16];
  vsu9[26] -= scalprod*vsu9[17];


// Normalize the row `2'

  xn = 1./sqrt( norm(vsu9[18])
    + norm(vsu9[19])
    + norm(vsu9[20])
    + norm(vsu9[21])
    + norm(vsu9[22])
    + norm(vsu9[23])
    + norm(vsu9[24])
    + norm(vsu9[25])
    + norm(vsu9[26]) );

    vsu9[18]*=xn;
    vsu9[19]*=xn;
    vsu9[20]*=xn;
    vsu9[21]*=xn;
    vsu9[22]*=xn;
    vsu9[23]*=xn;
    vsu9[24]*=xn;
    vsu9[25]*=xn;
    vsu9[26]*=xn;


// Make the row `3' orthogonal to the row `0':

  scalprod = conj( vsu9[0] )*vsu9[27]
    + conj( vsu9[1] )*vsu9[28]
    + conj( vsu9[2] )*vsu9[29]
    + conj( vsu9[3] )*vsu9[30]
    + conj( vsu9[4] )*vsu9[31]
    + conj( vsu9[5] )*vsu9[32]
    + conj( vsu9[6] )*vsu9[33]
    + conj( vsu9[7] )*vsu9[34]
    + conj( vsu9[8] )*vsu9[35];

  vsu9[27] -= scalprod*vsu9[0];
  vsu9[28] -= scalprod*vsu9[1];
  vsu9[29] -= scalprod*vsu9[2];
  vsu9[30] -= scalprod*vsu9[3];
  vsu9[31] -= scalprod*vsu9[4];
  vsu9[32] -= scalprod*vsu9[5];
  vsu9[33] -= scalprod*vsu9[6];
  vsu9[34] -= scalprod*vsu9[7];
  vsu9[35] -= scalprod*vsu9[8];


// Make the row `3' orthogonal to the row `1':

  scalprod = conj( vsu9[9] )*vsu9[27]
    + conj( vsu9[10] )*vsu9[28]
    + conj( vsu9[11] )*vsu9[29]
    + conj( vsu9[12] )*vsu9[30]
    + conj( vsu9[13] )*vsu9[31]
    + conj( vsu9[14] )*vsu9[32]
    + conj( vsu9[15] )*vsu9[33]
    + conj( vsu9[16] )*vsu9[34]
    + conj( vsu9[17] )*vsu9[35];

  vsu9[27] -= scalprod*vsu9[9];
  vsu9[28] -= scalprod*vsu9[10];
  vsu9[29] -= scalprod*vsu9[11];
  vsu9[30] -= scalprod*vsu9[12];
  vsu9[31] -= scalprod*vsu9[13];
  vsu9[32] -= scalprod*vsu9[14];
  vsu9[33] -= scalprod*vsu9[15];
  vsu9[34] -= scalprod*vsu9[16];
  vsu9[35] -= scalprod*vsu9[17];


// Make the row `3' orthogonal to the row `2':

  scalprod = conj( vsu9[18] )*vsu9[27]
    + conj( vsu9[19] )*vsu9[28]
    + conj( vsu9[20] )*vsu9[29]
    + conj( vsu9[21] )*vsu9[30]
    + conj( vsu9[22] )*vsu9[31]
    + conj( vsu9[23] )*vsu9[32]
    + conj( vsu9[24] )*vsu9[33]
    + conj( vsu9[25] )*vsu9[34]
    + conj( vsu9[26] )*vsu9[35];

  vsu9[27] -= scalprod*vsu9[18];
  vsu9[28] -= scalprod*vsu9[19];
  vsu9[29] -= scalprod*vsu9[20];
  vsu9[30] -= scalprod*vsu9[21];
  vsu9[31] -= scalprod*vsu9[22];
  vsu9[32] -= scalprod*vsu9[23];
  vsu9[33] -= scalprod*vsu9[24];
  vsu9[34] -= scalprod*vsu9[25];
  vsu9[35] -= scalprod*vsu9[26];


// Normalize the row `3'

  xn = 1./sqrt( norm(vsu9[27])
    + norm(vsu9[28])
    + norm(vsu9[29])
    + norm(vsu9[30])
    + norm(vsu9[31])
    + norm(vsu9[32])
    + norm(vsu9[33])
    + norm(vsu9[34])
    + norm(vsu9[35]) );

    vsu9[27]*=xn;
    vsu9[28]*=xn;
    vsu9[29]*=xn;
    vsu9[30]*=xn;
    vsu9[31]*=xn;
    vsu9[32]*=xn;
    vsu9[33]*=xn;
    vsu9[34]*=xn;
    vsu9[35]*=xn;


// Make the row `4' orthogonal to the row `0':

  scalprod = conj( vsu9[0] )*vsu9[36]
    + conj( vsu9[1] )*vsu9[37]
    + conj( vsu9[2] )*vsu9[38]
    + conj( vsu9[3] )*vsu9[39]
    + conj( vsu9[4] )*vsu9[40]
    + conj( vsu9[5] )*vsu9[41]
    + conj( vsu9[6] )*vsu9[42]
    + conj( vsu9[7] )*vsu9[43]
    + conj( vsu9[8] )*vsu9[44];

  vsu9[36] -= scalprod*vsu9[0];
  vsu9[37] -= scalprod*vsu9[1];
  vsu9[38] -= scalprod*vsu9[2];
  vsu9[39] -= scalprod*vsu9[3];
  vsu9[40] -= scalprod*vsu9[4];
  vsu9[41] -= scalprod*vsu9[5];
  vsu9[42] -= scalprod*vsu9[6];
  vsu9[43] -= scalprod*vsu9[7];
  vsu9[44] -= scalprod*vsu9[8];


// Make the row `4' orthogonal to the row `1':

  scalprod = conj( vsu9[9] )*vsu9[36]
    + conj( vsu9[10] )*vsu9[37]
    + conj( vsu9[11] )*vsu9[38]
    + conj( vsu9[12] )*vsu9[39]
    + conj( vsu9[13] )*vsu9[40]
    + conj( vsu9[14] )*vsu9[41]
    + conj( vsu9[15] )*vsu9[42]
    + conj( vsu9[16] )*vsu9[43]
    + conj( vsu9[17] )*vsu9[44];

  vsu9[36] -= scalprod*vsu9[9];
  vsu9[37] -= scalprod*vsu9[10];
  vsu9[38] -= scalprod*vsu9[11];
  vsu9[39] -= scalprod*vsu9[12];
  vsu9[40] -= scalprod*vsu9[13];
  vsu9[41] -= scalprod*vsu9[14];
  vsu9[42] -= scalprod*vsu9[15];
  vsu9[43] -= scalprod*vsu9[16];
  vsu9[44] -= scalprod*vsu9[17];


// Make the row `4' orthogonal to the row `2':

  scalprod = conj( vsu9[18] )*vsu9[36]
    + conj( vsu9[19] )*vsu9[37]
    + conj( vsu9[20] )*vsu9[38]
    + conj( vsu9[21] )*vsu9[39]
    + conj( vsu9[22] )*vsu9[40]
    + conj( vsu9[23] )*vsu9[41]
    + conj( vsu9[24] )*vsu9[42]
    + conj( vsu9[25] )*vsu9[43]
    + conj( vsu9[26] )*vsu9[44];

  vsu9[36] -= scalprod*vsu9[18];
  vsu9[37] -= scalprod*vsu9[19];
  vsu9[38] -= scalprod*vsu9[20];
  vsu9[39] -= scalprod*vsu9[21];
  vsu9[40] -= scalprod*vsu9[22];
  vsu9[41] -= scalprod*vsu9[23];
  vsu9[42] -= scalprod*vsu9[24];
  vsu9[43] -= scalprod*vsu9[25];
  vsu9[44] -= scalprod*vsu9[26];


// Make the row `4' orthogonal to the row `3':

  scalprod = conj( vsu9[27] )*vsu9[36]
    + conj( vsu9[28] )*vsu9[37]
    + conj( vsu9[29] )*vsu9[38]
    + conj( vsu9[30] )*vsu9[39]
    + conj( vsu9[31] )*vsu9[40]
    + conj( vsu9[32] )*vsu9[41]
    + conj( vsu9[33] )*vsu9[42]
    + conj( vsu9[34] )*vsu9[43]
    + conj( vsu9[35] )*vsu9[44];

  vsu9[36] -= scalprod*vsu9[27];
  vsu9[37] -= scalprod*vsu9[28];
  vsu9[38] -= scalprod*vsu9[29];
  vsu9[39] -= scalprod*vsu9[30];
  vsu9[40] -= scalprod*vsu9[31];
  vsu9[41] -= scalprod*vsu9[32];
  vsu9[42] -= scalprod*vsu9[33];
  vsu9[43] -= scalprod*vsu9[34];
  vsu9[44] -= scalprod*vsu9[35];


// Normalize the row `4'

  xn = 1./sqrt( norm(vsu9[36])
    + norm(vsu9[37])
    + norm(vsu9[38])
    + norm(vsu9[39])
    + norm(vsu9[40])
    + norm(vsu9[41])
    + norm(vsu9[42])
    + norm(vsu9[43])
    + norm(vsu9[44]) );

    vsu9[36]*=xn;
    vsu9[37]*=xn;
    vsu9[38]*=xn;
    vsu9[39]*=xn;
    vsu9[40]*=xn;
    vsu9[41]*=xn;
    vsu9[42]*=xn;
    vsu9[43]*=xn;
    vsu9[44]*=xn;


// Make the row `5' orthogonal to the row `0':

  scalprod = conj( vsu9[0] )*vsu9[45]
    + conj( vsu9[1] )*vsu9[46]
    + conj( vsu9[2] )*vsu9[47]
    + conj( vsu9[3] )*vsu9[48]
    + conj( vsu9[4] )*vsu9[49]
    + conj( vsu9[5] )*vsu9[50]
    + conj( vsu9[6] )*vsu9[51]
    + conj( vsu9[7] )*vsu9[52]
    + conj( vsu9[8] )*vsu9[53];

  vsu9[45] -= scalprod*vsu9[0];
  vsu9[46] -= scalprod*vsu9[1];
  vsu9[47] -= scalprod*vsu9[2];
  vsu9[48] -= scalprod*vsu9[3];
  vsu9[49] -= scalprod*vsu9[4];
  vsu9[50] -= scalprod*vsu9[5];
  vsu9[51] -= scalprod*vsu9[6];
  vsu9[52] -= scalprod*vsu9[7];
  vsu9[53] -= scalprod*vsu9[8];


// Make the row `5' orthogonal to the row `1':

  scalprod = conj( vsu9[9] )*vsu9[45]
    + conj( vsu9[10] )*vsu9[46]
    + conj( vsu9[11] )*vsu9[47]
    + conj( vsu9[12] )*vsu9[48]
    + conj( vsu9[13] )*vsu9[49]
    + conj( vsu9[14] )*vsu9[50]
    + conj( vsu9[15] )*vsu9[51]
    + conj( vsu9[16] )*vsu9[52]
    + conj( vsu9[17] )*vsu9[53];

  vsu9[45] -= scalprod*vsu9[9];
  vsu9[46] -= scalprod*vsu9[10];
  vsu9[47] -= scalprod*vsu9[11];
  vsu9[48] -= scalprod*vsu9[12];
  vsu9[49] -= scalprod*vsu9[13];
  vsu9[50] -= scalprod*vsu9[14];
  vsu9[51] -= scalprod*vsu9[15];
  vsu9[52] -= scalprod*vsu9[16];
  vsu9[53] -= scalprod*vsu9[17];


// Make the row `5' orthogonal to the row `2':

  scalprod = conj( vsu9[18] )*vsu9[45]
    + conj( vsu9[19] )*vsu9[46]
    + conj( vsu9[20] )*vsu9[47]
    + conj( vsu9[21] )*vsu9[48]
    + conj( vsu9[22] )*vsu9[49]
    + conj( vsu9[23] )*vsu9[50]
    + conj( vsu9[24] )*vsu9[51]
    + conj( vsu9[25] )*vsu9[52]
    + conj( vsu9[26] )*vsu9[53];

  vsu9[45] -= scalprod*vsu9[18];
  vsu9[46] -= scalprod*vsu9[19];
  vsu9[47] -= scalprod*vsu9[20];
  vsu9[48] -= scalprod*vsu9[21];
  vsu9[49] -= scalprod*vsu9[22];
  vsu9[50] -= scalprod*vsu9[23];
  vsu9[51] -= scalprod*vsu9[24];
  vsu9[52] -= scalprod*vsu9[25];
  vsu9[53] -= scalprod*vsu9[26];


// Make the row `5' orthogonal to the row `3':

  scalprod = conj( vsu9[27] )*vsu9[45]
    + conj( vsu9[28] )*vsu9[46]
    + conj( vsu9[29] )*vsu9[47]
    + conj( vsu9[30] )*vsu9[48]
    + conj( vsu9[31] )*vsu9[49]
    + conj( vsu9[32] )*vsu9[50]
    + conj( vsu9[33] )*vsu9[51]
    + conj( vsu9[34] )*vsu9[52]
    + conj( vsu9[35] )*vsu9[53];

  vsu9[45] -= scalprod*vsu9[27];
  vsu9[46] -= scalprod*vsu9[28];
  vsu9[47] -= scalprod*vsu9[29];
  vsu9[48] -= scalprod*vsu9[30];
  vsu9[49] -= scalprod*vsu9[31];
  vsu9[50] -= scalprod*vsu9[32];
  vsu9[51] -= scalprod*vsu9[33];
  vsu9[52] -= scalprod*vsu9[34];
  vsu9[53] -= scalprod*vsu9[35];


// Make the row `5' orthogonal to the row `4':

  scalprod = conj( vsu9[36] )*vsu9[45]
    + conj( vsu9[37] )*vsu9[46]
    + conj( vsu9[38] )*vsu9[47]
    + conj( vsu9[39] )*vsu9[48]
    + conj( vsu9[40] )*vsu9[49]
    + conj( vsu9[41] )*vsu9[50]
    + conj( vsu9[42] )*vsu9[51]
    + conj( vsu9[43] )*vsu9[52]
    + conj( vsu9[44] )*vsu9[53];

  vsu9[45] -= scalprod*vsu9[36];
  vsu9[46] -= scalprod*vsu9[37];
  vsu9[47] -= scalprod*vsu9[38];
  vsu9[48] -= scalprod*vsu9[39];
  vsu9[49] -= scalprod*vsu9[40];
  vsu9[50] -= scalprod*vsu9[41];
  vsu9[51] -= scalprod*vsu9[42];
  vsu9[52] -= scalprod*vsu9[43];
  vsu9[53] -= scalprod*vsu9[44];


// Normalize the row `5'

  xn = 1./sqrt( norm(vsu9[45])
    + norm(vsu9[46])
    + norm(vsu9[47])
    + norm(vsu9[48])
    + norm(vsu9[49])
    + norm(vsu9[50])
    + norm(vsu9[51])
    + norm(vsu9[52])
    + norm(vsu9[53]) );

    vsu9[45]*=xn;
    vsu9[46]*=xn;
    vsu9[47]*=xn;
    vsu9[48]*=xn;
    vsu9[49]*=xn;
    vsu9[50]*=xn;
    vsu9[51]*=xn;
    vsu9[52]*=xn;
    vsu9[53]*=xn;


// Make the row `6' orthogonal to the row `0':

  scalprod = conj( vsu9[0] )*vsu9[54]
    + conj( vsu9[1] )*vsu9[55]
    + conj( vsu9[2] )*vsu9[56]
    + conj( vsu9[3] )*vsu9[57]
    + conj( vsu9[4] )*vsu9[58]
    + conj( vsu9[5] )*vsu9[59]
    + conj( vsu9[6] )*vsu9[60]
    + conj( vsu9[7] )*vsu9[61]
    + conj( vsu9[8] )*vsu9[62];

  vsu9[54] -= scalprod*vsu9[0];
  vsu9[55] -= scalprod*vsu9[1];
  vsu9[56] -= scalprod*vsu9[2];
  vsu9[57] -= scalprod*vsu9[3];
  vsu9[58] -= scalprod*vsu9[4];
  vsu9[59] -= scalprod*vsu9[5];
  vsu9[60] -= scalprod*vsu9[6];
  vsu9[61] -= scalprod*vsu9[7];
  vsu9[62] -= scalprod*vsu9[8];


// Make the row `6' orthogonal to the row `1':

  scalprod = conj( vsu9[9] )*vsu9[54]
    + conj( vsu9[10] )*vsu9[55]
    + conj( vsu9[11] )*vsu9[56]
    + conj( vsu9[12] )*vsu9[57]
    + conj( vsu9[13] )*vsu9[58]
    + conj( vsu9[14] )*vsu9[59]
    + conj( vsu9[15] )*vsu9[60]
    + conj( vsu9[16] )*vsu9[61]
    + conj( vsu9[17] )*vsu9[62];

  vsu9[54] -= scalprod*vsu9[9];
  vsu9[55] -= scalprod*vsu9[10];
  vsu9[56] -= scalprod*vsu9[11];
  vsu9[57] -= scalprod*vsu9[12];
  vsu9[58] -= scalprod*vsu9[13];
  vsu9[59] -= scalprod*vsu9[14];
  vsu9[60] -= scalprod*vsu9[15];
  vsu9[61] -= scalprod*vsu9[16];
  vsu9[62] -= scalprod*vsu9[17];


// Make the row `6' orthogonal to the row `2':

  scalprod = conj( vsu9[18] )*vsu9[54]
    + conj( vsu9[19] )*vsu9[55]
    + conj( vsu9[20] )*vsu9[56]
    + conj( vsu9[21] )*vsu9[57]
    + conj( vsu9[22] )*vsu9[58]
    + conj( vsu9[23] )*vsu9[59]
    + conj( vsu9[24] )*vsu9[60]
    + conj( vsu9[25] )*vsu9[61]
    + conj( vsu9[26] )*vsu9[62];

  vsu9[54] -= scalprod*vsu9[18];
  vsu9[55] -= scalprod*vsu9[19];
  vsu9[56] -= scalprod*vsu9[20];
  vsu9[57] -= scalprod*vsu9[21];
  vsu9[58] -= scalprod*vsu9[22];
  vsu9[59] -= scalprod*vsu9[23];
  vsu9[60] -= scalprod*vsu9[24];
  vsu9[61] -= scalprod*vsu9[25];
  vsu9[62] -= scalprod*vsu9[26];


// Make the row `6' orthogonal to the row `3':

  scalprod = conj( vsu9[27] )*vsu9[54]
    + conj( vsu9[28] )*vsu9[55]
    + conj( vsu9[29] )*vsu9[56]
    + conj( vsu9[30] )*vsu9[57]
    + conj( vsu9[31] )*vsu9[58]
    + conj( vsu9[32] )*vsu9[59]
    + conj( vsu9[33] )*vsu9[60]
    + conj( vsu9[34] )*vsu9[61]
    + conj( vsu9[35] )*vsu9[62];

  vsu9[54] -= scalprod*vsu9[27];
  vsu9[55] -= scalprod*vsu9[28];
  vsu9[56] -= scalprod*vsu9[29];
  vsu9[57] -= scalprod*vsu9[30];
  vsu9[58] -= scalprod*vsu9[31];
  vsu9[59] -= scalprod*vsu9[32];
  vsu9[60] -= scalprod*vsu9[33];
  vsu9[61] -= scalprod*vsu9[34];
  vsu9[62] -= scalprod*vsu9[35];


// Make the row `6' orthogonal to the row `4':

  scalprod = conj( vsu9[36] )*vsu9[54]
    + conj( vsu9[37] )*vsu9[55]
    + conj( vsu9[38] )*vsu9[56]
    + conj( vsu9[39] )*vsu9[57]
    + conj( vsu9[40] )*vsu9[58]
    + conj( vsu9[41] )*vsu9[59]
    + conj( vsu9[42] )*vsu9[60]
    + conj( vsu9[43] )*vsu9[61]
    + conj( vsu9[44] )*vsu9[62];

  vsu9[54] -= scalprod*vsu9[36];
  vsu9[55] -= scalprod*vsu9[37];
  vsu9[56] -= scalprod*vsu9[38];
  vsu9[57] -= scalprod*vsu9[39];
  vsu9[58] -= scalprod*vsu9[40];
  vsu9[59] -= scalprod*vsu9[41];
  vsu9[60] -= scalprod*vsu9[42];
  vsu9[61] -= scalprod*vsu9[43];
  vsu9[62] -= scalprod*vsu9[44];


// Make the row `6' orthogonal to the row `5':

  scalprod = conj( vsu9[45] )*vsu9[54]
    + conj( vsu9[46] )*vsu9[55]
    + conj( vsu9[47] )*vsu9[56]
    + conj( vsu9[48] )*vsu9[57]
    + conj( vsu9[49] )*vsu9[58]
    + conj( vsu9[50] )*vsu9[59]
    + conj( vsu9[51] )*vsu9[60]
    + conj( vsu9[52] )*vsu9[61]
    + conj( vsu9[53] )*vsu9[62];

  vsu9[54] -= scalprod*vsu9[45];
  vsu9[55] -= scalprod*vsu9[46];
  vsu9[56] -= scalprod*vsu9[47];
  vsu9[57] -= scalprod*vsu9[48];
  vsu9[58] -= scalprod*vsu9[49];
  vsu9[59] -= scalprod*vsu9[50];
  vsu9[60] -= scalprod*vsu9[51];
  vsu9[61] -= scalprod*vsu9[52];
  vsu9[62] -= scalprod*vsu9[53];


// Normalize the row `6'

  xn = 1./sqrt( norm(vsu9[54])
    + norm(vsu9[55])
    + norm(vsu9[56])
    + norm(vsu9[57])
    + norm(vsu9[58])
    + norm(vsu9[59])
    + norm(vsu9[60])
    + norm(vsu9[61])
    + norm(vsu9[62]) );

    vsu9[54]*=xn;
    vsu9[55]*=xn;
    vsu9[56]*=xn;
    vsu9[57]*=xn;
    vsu9[58]*=xn;
    vsu9[59]*=xn;
    vsu9[60]*=xn;
    vsu9[61]*=xn;
    vsu9[62]*=xn;


// Make the row `7' orthogonal to the row `0':

  scalprod = conj( vsu9[0] )*vsu9[63]
    + conj( vsu9[1] )*vsu9[64]
    + conj( vsu9[2] )*vsu9[65]
    + conj( vsu9[3] )*vsu9[66]
    + conj( vsu9[4] )*vsu9[67]
    + conj( vsu9[5] )*vsu9[68]
    + conj( vsu9[6] )*vsu9[69]
    + conj( vsu9[7] )*vsu9[70]
    + conj( vsu9[8] )*vsu9[71];

  vsu9[63] -= scalprod*vsu9[0];
  vsu9[64] -= scalprod*vsu9[1];
  vsu9[65] -= scalprod*vsu9[2];
  vsu9[66] -= scalprod*vsu9[3];
  vsu9[67] -= scalprod*vsu9[4];
  vsu9[68] -= scalprod*vsu9[5];
  vsu9[69] -= scalprod*vsu9[6];
  vsu9[70] -= scalprod*vsu9[7];
  vsu9[71] -= scalprod*vsu9[8];


// Make the row `7' orthogonal to the row `1':

  scalprod = conj( vsu9[9] )*vsu9[63]
    + conj( vsu9[10] )*vsu9[64]
    + conj( vsu9[11] )*vsu9[65]
    + conj( vsu9[12] )*vsu9[66]
    + conj( vsu9[13] )*vsu9[67]
    + conj( vsu9[14] )*vsu9[68]
    + conj( vsu9[15] )*vsu9[69]
    + conj( vsu9[16] )*vsu9[70]
    + conj( vsu9[17] )*vsu9[71];

  vsu9[63] -= scalprod*vsu9[9];
  vsu9[64] -= scalprod*vsu9[10];
  vsu9[65] -= scalprod*vsu9[11];
  vsu9[66] -= scalprod*vsu9[12];
  vsu9[67] -= scalprod*vsu9[13];
  vsu9[68] -= scalprod*vsu9[14];
  vsu9[69] -= scalprod*vsu9[15];
  vsu9[70] -= scalprod*vsu9[16];
  vsu9[71] -= scalprod*vsu9[17];


// Make the row `7' orthogonal to the row `2':

  scalprod = conj( vsu9[18] )*vsu9[63]
    + conj( vsu9[19] )*vsu9[64]
    + conj( vsu9[20] )*vsu9[65]
    + conj( vsu9[21] )*vsu9[66]
    + conj( vsu9[22] )*vsu9[67]
    + conj( vsu9[23] )*vsu9[68]
    + conj( vsu9[24] )*vsu9[69]
    + conj( vsu9[25] )*vsu9[70]
    + conj( vsu9[26] )*vsu9[71];

  vsu9[63] -= scalprod*vsu9[18];
  vsu9[64] -= scalprod*vsu9[19];
  vsu9[65] -= scalprod*vsu9[20];
  vsu9[66] -= scalprod*vsu9[21];
  vsu9[67] -= scalprod*vsu9[22];
  vsu9[68] -= scalprod*vsu9[23];
  vsu9[69] -= scalprod*vsu9[24];
  vsu9[70] -= scalprod*vsu9[25];
  vsu9[71] -= scalprod*vsu9[26];


// Make the row `7' orthogonal to the row `3':

  scalprod = conj( vsu9[27] )*vsu9[63]
    + conj( vsu9[28] )*vsu9[64]
    + conj( vsu9[29] )*vsu9[65]
    + conj( vsu9[30] )*vsu9[66]
    + conj( vsu9[31] )*vsu9[67]
    + conj( vsu9[32] )*vsu9[68]
    + conj( vsu9[33] )*vsu9[69]
    + conj( vsu9[34] )*vsu9[70]
    + conj( vsu9[35] )*vsu9[71];

  vsu9[63] -= scalprod*vsu9[27];
  vsu9[64] -= scalprod*vsu9[28];
  vsu9[65] -= scalprod*vsu9[29];
  vsu9[66] -= scalprod*vsu9[30];
  vsu9[67] -= scalprod*vsu9[31];
  vsu9[68] -= scalprod*vsu9[32];
  vsu9[69] -= scalprod*vsu9[33];
  vsu9[70] -= scalprod*vsu9[34];
  vsu9[71] -= scalprod*vsu9[35];


// Make the row `7' orthogonal to the row `4':

  scalprod = conj( vsu9[36] )*vsu9[63]
    + conj( vsu9[37] )*vsu9[64]
    + conj( vsu9[38] )*vsu9[65]
    + conj( vsu9[39] )*vsu9[66]
    + conj( vsu9[40] )*vsu9[67]
    + conj( vsu9[41] )*vsu9[68]
    + conj( vsu9[42] )*vsu9[69]
    + conj( vsu9[43] )*vsu9[70]
    + conj( vsu9[44] )*vsu9[71];

  vsu9[63] -= scalprod*vsu9[36];
  vsu9[64] -= scalprod*vsu9[37];
  vsu9[65] -= scalprod*vsu9[38];
  vsu9[66] -= scalprod*vsu9[39];
  vsu9[67] -= scalprod*vsu9[40];
  vsu9[68] -= scalprod*vsu9[41];
  vsu9[69] -= scalprod*vsu9[42];
  vsu9[70] -= scalprod*vsu9[43];
  vsu9[71] -= scalprod*vsu9[44];


// Make the row `7' orthogonal to the row `5':

  scalprod = conj( vsu9[45] )*vsu9[63]
    + conj( vsu9[46] )*vsu9[64]
    + conj( vsu9[47] )*vsu9[65]
    + conj( vsu9[48] )*vsu9[66]
    + conj( vsu9[49] )*vsu9[67]
    + conj( vsu9[50] )*vsu9[68]
    + conj( vsu9[51] )*vsu9[69]
    + conj( vsu9[52] )*vsu9[70]
    + conj( vsu9[53] )*vsu9[71];

  vsu9[63] -= scalprod*vsu9[45];
  vsu9[64] -= scalprod*vsu9[46];
  vsu9[65] -= scalprod*vsu9[47];
  vsu9[66] -= scalprod*vsu9[48];
  vsu9[67] -= scalprod*vsu9[49];
  vsu9[68] -= scalprod*vsu9[50];
  vsu9[69] -= scalprod*vsu9[51];
  vsu9[70] -= scalprod*vsu9[52];
  vsu9[71] -= scalprod*vsu9[53];


// Make the row `7' orthogonal to the row `6':

  scalprod = conj( vsu9[54] )*vsu9[63]
    + conj( vsu9[55] )*vsu9[64]
    + conj( vsu9[56] )*vsu9[65]
    + conj( vsu9[57] )*vsu9[66]
    + conj( vsu9[58] )*vsu9[67]
    + conj( vsu9[59] )*vsu9[68]
    + conj( vsu9[60] )*vsu9[69]
    + conj( vsu9[61] )*vsu9[70]
    + conj( vsu9[62] )*vsu9[71];

  vsu9[63] -= scalprod*vsu9[54];
  vsu9[64] -= scalprod*vsu9[55];
  vsu9[65] -= scalprod*vsu9[56];
  vsu9[66] -= scalprod*vsu9[57];
  vsu9[67] -= scalprod*vsu9[58];
  vsu9[68] -= scalprod*vsu9[59];
  vsu9[69] -= scalprod*vsu9[60];
  vsu9[70] -= scalprod*vsu9[61];
  vsu9[71] -= scalprod*vsu9[62];


// Normalize the row `7'

  xn = 1./sqrt( norm(vsu9[63])
    + norm(vsu9[64])
    + norm(vsu9[65])
    + norm(vsu9[66])
    + norm(vsu9[67])
    + norm(vsu9[68])
    + norm(vsu9[69])
    + norm(vsu9[70])
    + norm(vsu9[71]) );

    vsu9[63]*=xn;
    vsu9[64]*=xn;
    vsu9[65]*=xn;
    vsu9[66]*=xn;
    vsu9[67]*=xn;
    vsu9[68]*=xn;
    vsu9[69]*=xn;
    vsu9[70]*=xn;
    vsu9[71]*=xn;


// Make the row `8' orthogonal to the row `0':

  scalprod = conj( vsu9[0] )*vsu9[72]
    + conj( vsu9[1] )*vsu9[73]
    + conj( vsu9[2] )*vsu9[74]
    + conj( vsu9[3] )*vsu9[75]
    + conj( vsu9[4] )*vsu9[76]
    + conj( vsu9[5] )*vsu9[77]
    + conj( vsu9[6] )*vsu9[78]
    + conj( vsu9[7] )*vsu9[79]
    + conj( vsu9[8] )*vsu9[80];

  vsu9[72] -= scalprod*vsu9[0];
  vsu9[73] -= scalprod*vsu9[1];
  vsu9[74] -= scalprod*vsu9[2];
  vsu9[75] -= scalprod*vsu9[3];
  vsu9[76] -= scalprod*vsu9[4];
  vsu9[77] -= scalprod*vsu9[5];
  vsu9[78] -= scalprod*vsu9[6];
  vsu9[79] -= scalprod*vsu9[7];
  vsu9[80] -= scalprod*vsu9[8];


// Make the row `8' orthogonal to the row `1':

  scalprod = conj( vsu9[9] )*vsu9[72]
    + conj( vsu9[10] )*vsu9[73]
    + conj( vsu9[11] )*vsu9[74]
    + conj( vsu9[12] )*vsu9[75]
    + conj( vsu9[13] )*vsu9[76]
    + conj( vsu9[14] )*vsu9[77]
    + conj( vsu9[15] )*vsu9[78]
    + conj( vsu9[16] )*vsu9[79]
    + conj( vsu9[17] )*vsu9[80];

  vsu9[72] -= scalprod*vsu9[9];
  vsu9[73] -= scalprod*vsu9[10];
  vsu9[74] -= scalprod*vsu9[11];
  vsu9[75] -= scalprod*vsu9[12];
  vsu9[76] -= scalprod*vsu9[13];
  vsu9[77] -= scalprod*vsu9[14];
  vsu9[78] -= scalprod*vsu9[15];
  vsu9[79] -= scalprod*vsu9[16];
  vsu9[80] -= scalprod*vsu9[17];


// Make the row `8' orthogonal to the row `2':

  scalprod = conj( vsu9[18] )*vsu9[72]
    + conj( vsu9[19] )*vsu9[73]
    + conj( vsu9[20] )*vsu9[74]
    + conj( vsu9[21] )*vsu9[75]
    + conj( vsu9[22] )*vsu9[76]
    + conj( vsu9[23] )*vsu9[77]
    + conj( vsu9[24] )*vsu9[78]
    + conj( vsu9[25] )*vsu9[79]
    + conj( vsu9[26] )*vsu9[80];

  vsu9[72] -= scalprod*vsu9[18];
  vsu9[73] -= scalprod*vsu9[19];
  vsu9[74] -= scalprod*vsu9[20];
  vsu9[75] -= scalprod*vsu9[21];
  vsu9[76] -= scalprod*vsu9[22];
  vsu9[77] -= scalprod*vsu9[23];
  vsu9[78] -= scalprod*vsu9[24];
  vsu9[79] -= scalprod*vsu9[25];
  vsu9[80] -= scalprod*vsu9[26];


// Make the row `8' orthogonal to the row `3':

  scalprod = conj( vsu9[27] )*vsu9[72]
    + conj( vsu9[28] )*vsu9[73]
    + conj( vsu9[29] )*vsu9[74]
    + conj( vsu9[30] )*vsu9[75]
    + conj( vsu9[31] )*vsu9[76]
    + conj( vsu9[32] )*vsu9[77]
    + conj( vsu9[33] )*vsu9[78]
    + conj( vsu9[34] )*vsu9[79]
    + conj( vsu9[35] )*vsu9[80];

  vsu9[72] -= scalprod*vsu9[27];
  vsu9[73] -= scalprod*vsu9[28];
  vsu9[74] -= scalprod*vsu9[29];
  vsu9[75] -= scalprod*vsu9[30];
  vsu9[76] -= scalprod*vsu9[31];
  vsu9[77] -= scalprod*vsu9[32];
  vsu9[78] -= scalprod*vsu9[33];
  vsu9[79] -= scalprod*vsu9[34];
  vsu9[80] -= scalprod*vsu9[35];


// Make the row `8' orthogonal to the row `4':

  scalprod = conj( vsu9[36] )*vsu9[72]
    + conj( vsu9[37] )*vsu9[73]
    + conj( vsu9[38] )*vsu9[74]
    + conj( vsu9[39] )*vsu9[75]
    + conj( vsu9[40] )*vsu9[76]
    + conj( vsu9[41] )*vsu9[77]
    + conj( vsu9[42] )*vsu9[78]
    + conj( vsu9[43] )*vsu9[79]
    + conj( vsu9[44] )*vsu9[80];

  vsu9[72] -= scalprod*vsu9[36];
  vsu9[73] -= scalprod*vsu9[37];
  vsu9[74] -= scalprod*vsu9[38];
  vsu9[75] -= scalprod*vsu9[39];
  vsu9[76] -= scalprod*vsu9[40];
  vsu9[77] -= scalprod*vsu9[41];
  vsu9[78] -= scalprod*vsu9[42];
  vsu9[79] -= scalprod*vsu9[43];
  vsu9[80] -= scalprod*vsu9[44];


// Make the row `8' orthogonal to the row `5':

  scalprod = conj( vsu9[45] )*vsu9[72]
    + conj( vsu9[46] )*vsu9[73]
    + conj( vsu9[47] )*vsu9[74]
    + conj( vsu9[48] )*vsu9[75]
    + conj( vsu9[49] )*vsu9[76]
    + conj( vsu9[50] )*vsu9[77]
    + conj( vsu9[51] )*vsu9[78]
    + conj( vsu9[52] )*vsu9[79]
    + conj( vsu9[53] )*vsu9[80];

  vsu9[72] -= scalprod*vsu9[45];
  vsu9[73] -= scalprod*vsu9[46];
  vsu9[74] -= scalprod*vsu9[47];
  vsu9[75] -= scalprod*vsu9[48];
  vsu9[76] -= scalprod*vsu9[49];
  vsu9[77] -= scalprod*vsu9[50];
  vsu9[78] -= scalprod*vsu9[51];
  vsu9[79] -= scalprod*vsu9[52];
  vsu9[80] -= scalprod*vsu9[53];


// Make the row `8' orthogonal to the row `6':

  scalprod = conj( vsu9[54] )*vsu9[72]
    + conj( vsu9[55] )*vsu9[73]
    + conj( vsu9[56] )*vsu9[74]
    + conj( vsu9[57] )*vsu9[75]
    + conj( vsu9[58] )*vsu9[76]
    + conj( vsu9[59] )*vsu9[77]
    + conj( vsu9[60] )*vsu9[78]
    + conj( vsu9[61] )*vsu9[79]
    + conj( vsu9[62] )*vsu9[80];

  vsu9[72] -= scalprod*vsu9[54];
  vsu9[73] -= scalprod*vsu9[55];
  vsu9[74] -= scalprod*vsu9[56];
  vsu9[75] -= scalprod*vsu9[57];
  vsu9[76] -= scalprod*vsu9[58];
  vsu9[77] -= scalprod*vsu9[59];
  vsu9[78] -= scalprod*vsu9[60];
  vsu9[79] -= scalprod*vsu9[61];
  vsu9[80] -= scalprod*vsu9[62];


// Make the row `8' orthogonal to the row `7':

  scalprod = conj( vsu9[63] )*vsu9[72]
    + conj( vsu9[64] )*vsu9[73]
    + conj( vsu9[65] )*vsu9[74]
    + conj( vsu9[66] )*vsu9[75]
    + conj( vsu9[67] )*vsu9[76]
    + conj( vsu9[68] )*vsu9[77]
    + conj( vsu9[69] )*vsu9[78]
    + conj( vsu9[70] )*vsu9[79]
    + conj( vsu9[71] )*vsu9[80];

  vsu9[72] -= scalprod*vsu9[63];
  vsu9[73] -= scalprod*vsu9[64];
  vsu9[74] -= scalprod*vsu9[65];
  vsu9[75] -= scalprod*vsu9[66];
  vsu9[76] -= scalprod*vsu9[67];
  vsu9[77] -= scalprod*vsu9[68];
  vsu9[78] -= scalprod*vsu9[69];
  vsu9[79] -= scalprod*vsu9[70];
  vsu9[80] -= scalprod*vsu9[71];


// Normalize the row `8'

  xn = 1./sqrt( norm(vsu9[72])
    + norm(vsu9[73])
    + norm(vsu9[74])
    + norm(vsu9[75])
    + norm(vsu9[76])
    + norm(vsu9[77])
    + norm(vsu9[78])
    + norm(vsu9[79])
    + norm(vsu9[80]) );

    vsu9[72]*=xn;
    vsu9[73]*=xn;
    vsu9[74]*=xn;
    vsu9[75]*=xn;
    vsu9[76]*=xn;
    vsu9[77]*=xn;
    vsu9[78]*=xn;
    vsu9[79]*=xn;
    vsu9[80]*=xn;


// Impose unimodularity by arranging the phase of the last row
// to compensate for the phase of the determinant

  for (i=0; i<Ncol; i++)
  for (j=0; j<Ncol; j++) {
    AT[2*(j+Ncol*i)]=real(vsu9[j*Ncol+i]);
    AT[2*(j+Ncol*i)+1]=imag(vsu9[j*Ncol+i]);
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

  vsu9[72] *= determinant_conjugate_phase;
  vsu9[73] *= determinant_conjugate_phase;
  vsu9[74] *= determinant_conjugate_phase;
  vsu9[75] *= determinant_conjugate_phase;
  vsu9[76] *= determinant_conjugate_phase;
  vsu9[77] *= determinant_conjugate_phase;
  vsu9[78] *= determinant_conjugate_phase;
  vsu9[79] *= determinant_conjugate_phase;
  vsu9[80] *= determinant_conjugate_phase;


}

#endif
