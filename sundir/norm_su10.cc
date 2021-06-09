#ifndef __norm_su5_h__
#define __norm_su5_h__

void norm_su10(dc *vsu10 ) {

// Takes the first 9 rows of a 10x10 complex matrix,
// and makes an SU(10) matrix out of them.

  dc determinant_conjugate_phase=dc(1.,0.);
  struct complex_double b[Ncol], DUMMY[1][1], WORK[2*Ncol];
  double AT[2*Ncolsquare];
  int i, j, ok, c1, c2, c3;
  char c4;

  dc scalprod=dc(0.,0.);
  double xn;

// Normalize the row `0' of vsu10

  xn = 1./sqrt( norm(vsu10[0])
    + norm(vsu10[1])
    + norm(vsu10[2])
    + norm(vsu10[3])
    + norm(vsu10[4])
    + norm(vsu10[5])
    + norm(vsu10[6])
    + norm(vsu10[7])
    + norm(vsu10[8])
    + norm(vsu10[9]) );

  vsu10[0]*=xn;
  vsu10[1]*=xn;
  vsu10[2]*=xn;
  vsu10[3]*=xn;
  vsu10[4]*=xn;
  vsu10[5]*=xn;
  vsu10[6]*=xn;
  vsu10[7]*=xn;
  vsu10[8]*=xn;
  vsu10[9]*=xn;

// Make the row `1' orthogonal to the row `0':

  scalprod = conj( vsu10[0] )*vsu10[10]
    + conj( vsu10[1] )*vsu10[11]
    + conj( vsu10[2] )*vsu10[12]
    + conj( vsu10[3] )*vsu10[13]
    + conj( vsu10[4] )*vsu10[14]
    + conj( vsu10[5] )*vsu10[15]
    + conj( vsu10[6] )*vsu10[16]
    + conj( vsu10[7] )*vsu10[17]
    + conj( vsu10[8] )*vsu10[18]
    + conj( vsu10[9] )*vsu10[19];

  vsu10[10] -= scalprod*vsu10[0];
  vsu10[11] -= scalprod*vsu10[1];
  vsu10[12] -= scalprod*vsu10[2];
  vsu10[13] -= scalprod*vsu10[3];
  vsu10[14] -= scalprod*vsu10[4];
  vsu10[15] -= scalprod*vsu10[5];
  vsu10[16] -= scalprod*vsu10[6];
  vsu10[17] -= scalprod*vsu10[7];
  vsu10[18] -= scalprod*vsu10[8];
  vsu10[19] -= scalprod*vsu10[9];


// Normalize the row `1'

  xn = 1./sqrt( norm(vsu10[10])
    + norm(vsu10[11])
    + norm(vsu10[12])
    + norm(vsu10[13])
    + norm(vsu10[14])
    + norm(vsu10[15])
    + norm(vsu10[16])
    + norm(vsu10[17])
    + norm(vsu10[18])
    + norm(vsu10[19]) );

    vsu10[10]*=xn;
    vsu10[11]*=xn;
    vsu10[12]*=xn;
    vsu10[13]*=xn;
    vsu10[14]*=xn;
    vsu10[15]*=xn;
    vsu10[16]*=xn;
    vsu10[17]*=xn;
    vsu10[18]*=xn;
    vsu10[19]*=xn;


// Make the row `2' orthogonal to the row `0':

  scalprod = conj( vsu10[0] )*vsu10[20]
    + conj( vsu10[1] )*vsu10[21]
    + conj( vsu10[2] )*vsu10[22]
    + conj( vsu10[3] )*vsu10[23]
    + conj( vsu10[4] )*vsu10[24]
    + conj( vsu10[5] )*vsu10[25]
    + conj( vsu10[6] )*vsu10[26]
    + conj( vsu10[7] )*vsu10[27]
    + conj( vsu10[8] )*vsu10[28]
    + conj( vsu10[9] )*vsu10[29];

  vsu10[20] -= scalprod*vsu10[0];
  vsu10[21] -= scalprod*vsu10[1];
  vsu10[22] -= scalprod*vsu10[2];
  vsu10[23] -= scalprod*vsu10[3];
  vsu10[24] -= scalprod*vsu10[4];
  vsu10[25] -= scalprod*vsu10[5];
  vsu10[26] -= scalprod*vsu10[6];
  vsu10[27] -= scalprod*vsu10[7];
  vsu10[28] -= scalprod*vsu10[8];
  vsu10[29] -= scalprod*vsu10[9];


// Make the row `2' orthogonal to the row `1':

  scalprod = conj( vsu10[10] )*vsu10[20]
    + conj( vsu10[11] )*vsu10[21]
    + conj( vsu10[12] )*vsu10[22]
    + conj( vsu10[13] )*vsu10[23]
    + conj( vsu10[14] )*vsu10[24]
    + conj( vsu10[15] )*vsu10[25]
    + conj( vsu10[16] )*vsu10[26]
    + conj( vsu10[17] )*vsu10[27]
    + conj( vsu10[18] )*vsu10[28]
    + conj( vsu10[19] )*vsu10[29];

  vsu10[20] -= scalprod*vsu10[10];
  vsu10[21] -= scalprod*vsu10[11];
  vsu10[22] -= scalprod*vsu10[12];
  vsu10[23] -= scalprod*vsu10[13];
  vsu10[24] -= scalprod*vsu10[14];
  vsu10[25] -= scalprod*vsu10[15];
  vsu10[26] -= scalprod*vsu10[16];
  vsu10[27] -= scalprod*vsu10[17];
  vsu10[28] -= scalprod*vsu10[18];
  vsu10[29] -= scalprod*vsu10[19];


// Normalize the row `2'

  xn = 1./sqrt( norm(vsu10[20])
    + norm(vsu10[21])
    + norm(vsu10[22])
    + norm(vsu10[23])
    + norm(vsu10[24])
    + norm(vsu10[25])
    + norm(vsu10[26])
    + norm(vsu10[27])
    + norm(vsu10[28])
    + norm(vsu10[29]) );

    vsu10[20]*=xn;
    vsu10[21]*=xn;
    vsu10[22]*=xn;
    vsu10[23]*=xn;
    vsu10[24]*=xn;
    vsu10[25]*=xn;
    vsu10[26]*=xn;
    vsu10[27]*=xn;
    vsu10[28]*=xn;
    vsu10[29]*=xn;


// Make the row `3' orthogonal to the row `0':

  scalprod = conj( vsu10[0] )*vsu10[30]
    + conj( vsu10[1] )*vsu10[31]
    + conj( vsu10[2] )*vsu10[32]
    + conj( vsu10[3] )*vsu10[33]
    + conj( vsu10[4] )*vsu10[34]
    + conj( vsu10[5] )*vsu10[35]
    + conj( vsu10[6] )*vsu10[36]
    + conj( vsu10[7] )*vsu10[37]
    + conj( vsu10[8] )*vsu10[38]
    + conj( vsu10[9] )*vsu10[39];

  vsu10[30] -= scalprod*vsu10[0];
  vsu10[31] -= scalprod*vsu10[1];
  vsu10[32] -= scalprod*vsu10[2];
  vsu10[33] -= scalprod*vsu10[3];
  vsu10[34] -= scalprod*vsu10[4];
  vsu10[35] -= scalprod*vsu10[5];
  vsu10[36] -= scalprod*vsu10[6];
  vsu10[37] -= scalprod*vsu10[7];
  vsu10[38] -= scalprod*vsu10[8];
  vsu10[39] -= scalprod*vsu10[9];


// Make the row `3' orthogonal to the row `1':

  scalprod = conj( vsu10[10] )*vsu10[30]
    + conj( vsu10[11] )*vsu10[31]
    + conj( vsu10[12] )*vsu10[32]
    + conj( vsu10[13] )*vsu10[33]
    + conj( vsu10[14] )*vsu10[34]
    + conj( vsu10[15] )*vsu10[35]
    + conj( vsu10[16] )*vsu10[36]
    + conj( vsu10[17] )*vsu10[37]
    + conj( vsu10[18] )*vsu10[38]
    + conj( vsu10[19] )*vsu10[39];

  vsu10[30] -= scalprod*vsu10[10];
  vsu10[31] -= scalprod*vsu10[11];
  vsu10[32] -= scalprod*vsu10[12];
  vsu10[33] -= scalprod*vsu10[13];
  vsu10[34] -= scalprod*vsu10[14];
  vsu10[35] -= scalprod*vsu10[15];
  vsu10[36] -= scalprod*vsu10[16];
  vsu10[37] -= scalprod*vsu10[17];
  vsu10[38] -= scalprod*vsu10[18];
  vsu10[39] -= scalprod*vsu10[19];


// Make the row `3' orthogonal to the row `2':

  scalprod = conj( vsu10[20] )*vsu10[30]
    + conj( vsu10[21] )*vsu10[31]
    + conj( vsu10[22] )*vsu10[32]
    + conj( vsu10[23] )*vsu10[33]
    + conj( vsu10[24] )*vsu10[34]
    + conj( vsu10[25] )*vsu10[35]
    + conj( vsu10[26] )*vsu10[36]
    + conj( vsu10[27] )*vsu10[37]
    + conj( vsu10[28] )*vsu10[38]
    + conj( vsu10[29] )*vsu10[39];

  vsu10[30] -= scalprod*vsu10[20];
  vsu10[31] -= scalprod*vsu10[21];
  vsu10[32] -= scalprod*vsu10[22];
  vsu10[33] -= scalprod*vsu10[23];
  vsu10[34] -= scalprod*vsu10[24];
  vsu10[35] -= scalprod*vsu10[25];
  vsu10[36] -= scalprod*vsu10[26];
  vsu10[37] -= scalprod*vsu10[27];
  vsu10[38] -= scalprod*vsu10[28];
  vsu10[39] -= scalprod*vsu10[29];


// Normalize the row `3'

  xn = 1./sqrt( norm(vsu10[30])
    + norm(vsu10[31])
    + norm(vsu10[32])
    + norm(vsu10[33])
    + norm(vsu10[34])
    + norm(vsu10[35])
    + norm(vsu10[36])
    + norm(vsu10[37])
    + norm(vsu10[38])
    + norm(vsu10[39]) );

    vsu10[30]*=xn;
    vsu10[31]*=xn;
    vsu10[32]*=xn;
    vsu10[33]*=xn;
    vsu10[34]*=xn;
    vsu10[35]*=xn;
    vsu10[36]*=xn;
    vsu10[37]*=xn;
    vsu10[38]*=xn;
    vsu10[39]*=xn;


// Make the row `4' orthogonal to the row `0':

  scalprod = conj( vsu10[0] )*vsu10[40]
    + conj( vsu10[1] )*vsu10[41]
    + conj( vsu10[2] )*vsu10[42]
    + conj( vsu10[3] )*vsu10[43]
    + conj( vsu10[4] )*vsu10[44]
    + conj( vsu10[5] )*vsu10[45]
    + conj( vsu10[6] )*vsu10[46]
    + conj( vsu10[7] )*vsu10[47]
    + conj( vsu10[8] )*vsu10[48]
    + conj( vsu10[9] )*vsu10[49];

  vsu10[40] -= scalprod*vsu10[0];
  vsu10[41] -= scalprod*vsu10[1];
  vsu10[42] -= scalprod*vsu10[2];
  vsu10[43] -= scalprod*vsu10[3];
  vsu10[44] -= scalprod*vsu10[4];
  vsu10[45] -= scalprod*vsu10[5];
  vsu10[46] -= scalprod*vsu10[6];
  vsu10[47] -= scalprod*vsu10[7];
  vsu10[48] -= scalprod*vsu10[8];
  vsu10[49] -= scalprod*vsu10[9];


// Make the row `4' orthogonal to the row `1':

  scalprod = conj( vsu10[10] )*vsu10[40]
    + conj( vsu10[11] )*vsu10[41]
    + conj( vsu10[12] )*vsu10[42]
    + conj( vsu10[13] )*vsu10[43]
    + conj( vsu10[14] )*vsu10[44]
    + conj( vsu10[15] )*vsu10[45]
    + conj( vsu10[16] )*vsu10[46]
    + conj( vsu10[17] )*vsu10[47]
    + conj( vsu10[18] )*vsu10[48]
    + conj( vsu10[19] )*vsu10[49];

  vsu10[40] -= scalprod*vsu10[10];
  vsu10[41] -= scalprod*vsu10[11];
  vsu10[42] -= scalprod*vsu10[12];
  vsu10[43] -= scalprod*vsu10[13];
  vsu10[44] -= scalprod*vsu10[14];
  vsu10[45] -= scalprod*vsu10[15];
  vsu10[46] -= scalprod*vsu10[16];
  vsu10[47] -= scalprod*vsu10[17];
  vsu10[48] -= scalprod*vsu10[18];
  vsu10[49] -= scalprod*vsu10[19];


// Make the row `4' orthogonal to the row `2':

  scalprod = conj( vsu10[20] )*vsu10[40]
    + conj( vsu10[21] )*vsu10[41]
    + conj( vsu10[22] )*vsu10[42]
    + conj( vsu10[23] )*vsu10[43]
    + conj( vsu10[24] )*vsu10[44]
    + conj( vsu10[25] )*vsu10[45]
    + conj( vsu10[26] )*vsu10[46]
    + conj( vsu10[27] )*vsu10[47]
    + conj( vsu10[28] )*vsu10[48]
    + conj( vsu10[29] )*vsu10[49];

  vsu10[40] -= scalprod*vsu10[20];
  vsu10[41] -= scalprod*vsu10[21];
  vsu10[42] -= scalprod*vsu10[22];
  vsu10[43] -= scalprod*vsu10[23];
  vsu10[44] -= scalprod*vsu10[24];
  vsu10[45] -= scalprod*vsu10[25];
  vsu10[46] -= scalprod*vsu10[26];
  vsu10[47] -= scalprod*vsu10[27];
  vsu10[48] -= scalprod*vsu10[28];
  vsu10[49] -= scalprod*vsu10[29];


// Make the row `4' orthogonal to the row `3':

  scalprod = conj( vsu10[30] )*vsu10[40]
    + conj( vsu10[31] )*vsu10[41]
    + conj( vsu10[32] )*vsu10[42]
    + conj( vsu10[33] )*vsu10[43]
    + conj( vsu10[34] )*vsu10[44]
    + conj( vsu10[35] )*vsu10[45]
    + conj( vsu10[36] )*vsu10[46]
    + conj( vsu10[37] )*vsu10[47]
    + conj( vsu10[38] )*vsu10[48]
    + conj( vsu10[39] )*vsu10[49];

  vsu10[40] -= scalprod*vsu10[30];
  vsu10[41] -= scalprod*vsu10[31];
  vsu10[42] -= scalprod*vsu10[32];
  vsu10[43] -= scalprod*vsu10[33];
  vsu10[44] -= scalprod*vsu10[34];
  vsu10[45] -= scalprod*vsu10[35];
  vsu10[46] -= scalprod*vsu10[36];
  vsu10[47] -= scalprod*vsu10[37];
  vsu10[48] -= scalprod*vsu10[38];
  vsu10[49] -= scalprod*vsu10[39];


// Normalize the row `4'

  xn = 1./sqrt( norm(vsu10[40])
    + norm(vsu10[41])
    + norm(vsu10[42])
    + norm(vsu10[43])
    + norm(vsu10[44])
    + norm(vsu10[45])
    + norm(vsu10[46])
    + norm(vsu10[47])
    + norm(vsu10[48])
    + norm(vsu10[49]) );

    vsu10[40]*=xn;
    vsu10[41]*=xn;
    vsu10[42]*=xn;
    vsu10[43]*=xn;
    vsu10[44]*=xn;
    vsu10[45]*=xn;
    vsu10[46]*=xn;
    vsu10[47]*=xn;
    vsu10[48]*=xn;
    vsu10[49]*=xn;


// Make the row `5' orthogonal to the row `0':

  scalprod = conj( vsu10[0] )*vsu10[50]
    + conj( vsu10[1] )*vsu10[51]
    + conj( vsu10[2] )*vsu10[52]
    + conj( vsu10[3] )*vsu10[53]
    + conj( vsu10[4] )*vsu10[54]
    + conj( vsu10[5] )*vsu10[55]
    + conj( vsu10[6] )*vsu10[56]
    + conj( vsu10[7] )*vsu10[57]
    + conj( vsu10[8] )*vsu10[58]
    + conj( vsu10[9] )*vsu10[59];

  vsu10[50] -= scalprod*vsu10[0];
  vsu10[51] -= scalprod*vsu10[1];
  vsu10[52] -= scalprod*vsu10[2];
  vsu10[53] -= scalprod*vsu10[3];
  vsu10[54] -= scalprod*vsu10[4];
  vsu10[55] -= scalprod*vsu10[5];
  vsu10[56] -= scalprod*vsu10[6];
  vsu10[57] -= scalprod*vsu10[7];
  vsu10[58] -= scalprod*vsu10[8];
  vsu10[59] -= scalprod*vsu10[9];


// Make the row `5' orthogonal to the row `1':

  scalprod = conj( vsu10[10] )*vsu10[50]
    + conj( vsu10[11] )*vsu10[51]
    + conj( vsu10[12] )*vsu10[52]
    + conj( vsu10[13] )*vsu10[53]
    + conj( vsu10[14] )*vsu10[54]
    + conj( vsu10[15] )*vsu10[55]
    + conj( vsu10[16] )*vsu10[56]
    + conj( vsu10[17] )*vsu10[57]
    + conj( vsu10[18] )*vsu10[58]
    + conj( vsu10[19] )*vsu10[59];

  vsu10[50] -= scalprod*vsu10[10];
  vsu10[51] -= scalprod*vsu10[11];
  vsu10[52] -= scalprod*vsu10[12];
  vsu10[53] -= scalprod*vsu10[13];
  vsu10[54] -= scalprod*vsu10[14];
  vsu10[55] -= scalprod*vsu10[15];
  vsu10[56] -= scalprod*vsu10[16];
  vsu10[57] -= scalprod*vsu10[17];
  vsu10[58] -= scalprod*vsu10[18];
  vsu10[59] -= scalprod*vsu10[19];


// Make the row `5' orthogonal to the row `2':

  scalprod = conj( vsu10[20] )*vsu10[50]
    + conj( vsu10[21] )*vsu10[51]
    + conj( vsu10[22] )*vsu10[52]
    + conj( vsu10[23] )*vsu10[53]
    + conj( vsu10[24] )*vsu10[54]
    + conj( vsu10[25] )*vsu10[55]
    + conj( vsu10[26] )*vsu10[56]
    + conj( vsu10[27] )*vsu10[57]
    + conj( vsu10[28] )*vsu10[58]
    + conj( vsu10[29] )*vsu10[59];

  vsu10[50] -= scalprod*vsu10[20];
  vsu10[51] -= scalprod*vsu10[21];
  vsu10[52] -= scalprod*vsu10[22];
  vsu10[53] -= scalprod*vsu10[23];
  vsu10[54] -= scalprod*vsu10[24];
  vsu10[55] -= scalprod*vsu10[25];
  vsu10[56] -= scalprod*vsu10[26];
  vsu10[57] -= scalprod*vsu10[27];
  vsu10[58] -= scalprod*vsu10[28];
  vsu10[59] -= scalprod*vsu10[29];


// Make the row `5' orthogonal to the row `3':

  scalprod = conj( vsu10[30] )*vsu10[50]
    + conj( vsu10[31] )*vsu10[51]
    + conj( vsu10[32] )*vsu10[52]
    + conj( vsu10[33] )*vsu10[53]
    + conj( vsu10[34] )*vsu10[54]
    + conj( vsu10[35] )*vsu10[55]
    + conj( vsu10[36] )*vsu10[56]
    + conj( vsu10[37] )*vsu10[57]
    + conj( vsu10[38] )*vsu10[58]
    + conj( vsu10[39] )*vsu10[59];

  vsu10[50] -= scalprod*vsu10[30];
  vsu10[51] -= scalprod*vsu10[31];
  vsu10[52] -= scalprod*vsu10[32];
  vsu10[53] -= scalprod*vsu10[33];
  vsu10[54] -= scalprod*vsu10[34];
  vsu10[55] -= scalprod*vsu10[35];
  vsu10[56] -= scalprod*vsu10[36];
  vsu10[57] -= scalprod*vsu10[37];
  vsu10[58] -= scalprod*vsu10[38];
  vsu10[59] -= scalprod*vsu10[39];


// Make the row `5' orthogonal to the row `4':

  scalprod = conj( vsu10[40] )*vsu10[50]
    + conj( vsu10[41] )*vsu10[51]
    + conj( vsu10[42] )*vsu10[52]
    + conj( vsu10[43] )*vsu10[53]
    + conj( vsu10[44] )*vsu10[54]
    + conj( vsu10[45] )*vsu10[55]
    + conj( vsu10[46] )*vsu10[56]
    + conj( vsu10[47] )*vsu10[57]
    + conj( vsu10[48] )*vsu10[58]
    + conj( vsu10[49] )*vsu10[59];

  vsu10[50] -= scalprod*vsu10[40];
  vsu10[51] -= scalprod*vsu10[41];
  vsu10[52] -= scalprod*vsu10[42];
  vsu10[53] -= scalprod*vsu10[43];
  vsu10[54] -= scalprod*vsu10[44];
  vsu10[55] -= scalprod*vsu10[45];
  vsu10[56] -= scalprod*vsu10[46];
  vsu10[57] -= scalprod*vsu10[47];
  vsu10[58] -= scalprod*vsu10[48];
  vsu10[59] -= scalprod*vsu10[49];


// Normalize the row `5'

  xn = 1./sqrt( norm(vsu10[50])
    + norm(vsu10[51])
    + norm(vsu10[52])
    + norm(vsu10[53])
    + norm(vsu10[54])
    + norm(vsu10[55])
    + norm(vsu10[56])
    + norm(vsu10[57])
    + norm(vsu10[58])
    + norm(vsu10[59]) );

    vsu10[50]*=xn;
    vsu10[51]*=xn;
    vsu10[52]*=xn;
    vsu10[53]*=xn;
    vsu10[54]*=xn;
    vsu10[55]*=xn;
    vsu10[56]*=xn;
    vsu10[57]*=xn;
    vsu10[58]*=xn;
    vsu10[59]*=xn;


// Make the row `6' orthogonal to the row `0':

  scalprod = conj( vsu10[0] )*vsu10[60]
    + conj( vsu10[1] )*vsu10[61]
    + conj( vsu10[2] )*vsu10[62]
    + conj( vsu10[3] )*vsu10[63]
    + conj( vsu10[4] )*vsu10[64]
    + conj( vsu10[5] )*vsu10[65]
    + conj( vsu10[6] )*vsu10[66]
    + conj( vsu10[7] )*vsu10[67]
    + conj( vsu10[8] )*vsu10[68]
    + conj( vsu10[9] )*vsu10[69];

  vsu10[60] -= scalprod*vsu10[0];
  vsu10[61] -= scalprod*vsu10[1];
  vsu10[62] -= scalprod*vsu10[2];
  vsu10[63] -= scalprod*vsu10[3];
  vsu10[64] -= scalprod*vsu10[4];
  vsu10[65] -= scalprod*vsu10[5];
  vsu10[66] -= scalprod*vsu10[6];
  vsu10[67] -= scalprod*vsu10[7];
  vsu10[68] -= scalprod*vsu10[8];
  vsu10[69] -= scalprod*vsu10[9];


// Make the row `6' orthogonal to the row `1':

  scalprod = conj( vsu10[10] )*vsu10[60]
    + conj( vsu10[11] )*vsu10[61]
    + conj( vsu10[12] )*vsu10[62]
    + conj( vsu10[13] )*vsu10[63]
    + conj( vsu10[14] )*vsu10[64]
    + conj( vsu10[15] )*vsu10[65]
    + conj( vsu10[16] )*vsu10[66]
    + conj( vsu10[17] )*vsu10[67]
    + conj( vsu10[18] )*vsu10[68]
    + conj( vsu10[19] )*vsu10[69];

  vsu10[60] -= scalprod*vsu10[10];
  vsu10[61] -= scalprod*vsu10[11];
  vsu10[62] -= scalprod*vsu10[12];
  vsu10[63] -= scalprod*vsu10[13];
  vsu10[64] -= scalprod*vsu10[14];
  vsu10[65] -= scalprod*vsu10[15];
  vsu10[66] -= scalprod*vsu10[16];
  vsu10[67] -= scalprod*vsu10[17];
  vsu10[68] -= scalprod*vsu10[18];
  vsu10[69] -= scalprod*vsu10[19];


// Make the row `6' orthogonal to the row `2':

  scalprod = conj( vsu10[20] )*vsu10[60]
    + conj( vsu10[21] )*vsu10[61]
    + conj( vsu10[22] )*vsu10[62]
    + conj( vsu10[23] )*vsu10[63]
    + conj( vsu10[24] )*vsu10[64]
    + conj( vsu10[25] )*vsu10[65]
    + conj( vsu10[26] )*vsu10[66]
    + conj( vsu10[27] )*vsu10[67]
    + conj( vsu10[28] )*vsu10[68]
    + conj( vsu10[29] )*vsu10[69];

  vsu10[60] -= scalprod*vsu10[20];
  vsu10[61] -= scalprod*vsu10[21];
  vsu10[62] -= scalprod*vsu10[22];
  vsu10[63] -= scalprod*vsu10[23];
  vsu10[64] -= scalprod*vsu10[24];
  vsu10[65] -= scalprod*vsu10[25];
  vsu10[66] -= scalprod*vsu10[26];
  vsu10[67] -= scalprod*vsu10[27];
  vsu10[68] -= scalprod*vsu10[28];
  vsu10[69] -= scalprod*vsu10[29];


// Make the row `6' orthogonal to the row `3':

  scalprod = conj( vsu10[30] )*vsu10[60]
    + conj( vsu10[31] )*vsu10[61]
    + conj( vsu10[32] )*vsu10[62]
    + conj( vsu10[33] )*vsu10[63]
    + conj( vsu10[34] )*vsu10[64]
    + conj( vsu10[35] )*vsu10[65]
    + conj( vsu10[36] )*vsu10[66]
    + conj( vsu10[37] )*vsu10[67]
    + conj( vsu10[38] )*vsu10[68]
    + conj( vsu10[39] )*vsu10[69];

  vsu10[60] -= scalprod*vsu10[30];
  vsu10[61] -= scalprod*vsu10[31];
  vsu10[62] -= scalprod*vsu10[32];
  vsu10[63] -= scalprod*vsu10[33];
  vsu10[64] -= scalprod*vsu10[34];
  vsu10[65] -= scalprod*vsu10[35];
  vsu10[66] -= scalprod*vsu10[36];
  vsu10[67] -= scalprod*vsu10[37];
  vsu10[68] -= scalprod*vsu10[38];
  vsu10[69] -= scalprod*vsu10[39];


// Make the row `6' orthogonal to the row `4':

  scalprod = conj( vsu10[40] )*vsu10[60]
    + conj( vsu10[41] )*vsu10[61]
    + conj( vsu10[42] )*vsu10[62]
    + conj( vsu10[43] )*vsu10[63]
    + conj( vsu10[44] )*vsu10[64]
    + conj( vsu10[45] )*vsu10[65]
    + conj( vsu10[46] )*vsu10[66]
    + conj( vsu10[47] )*vsu10[67]
    + conj( vsu10[48] )*vsu10[68]
    + conj( vsu10[49] )*vsu10[69];

  vsu10[60] -= scalprod*vsu10[40];
  vsu10[61] -= scalprod*vsu10[41];
  vsu10[62] -= scalprod*vsu10[42];
  vsu10[63] -= scalprod*vsu10[43];
  vsu10[64] -= scalprod*vsu10[44];
  vsu10[65] -= scalprod*vsu10[45];
  vsu10[66] -= scalprod*vsu10[46];
  vsu10[67] -= scalprod*vsu10[47];
  vsu10[68] -= scalprod*vsu10[48];
  vsu10[69] -= scalprod*vsu10[49];


// Make the row `6' orthogonal to the row `5':

  scalprod = conj( vsu10[50] )*vsu10[60]
    + conj( vsu10[51] )*vsu10[61]
    + conj( vsu10[52] )*vsu10[62]
    + conj( vsu10[53] )*vsu10[63]
    + conj( vsu10[54] )*vsu10[64]
    + conj( vsu10[55] )*vsu10[65]
    + conj( vsu10[56] )*vsu10[66]
    + conj( vsu10[57] )*vsu10[67]
    + conj( vsu10[58] )*vsu10[68]
    + conj( vsu10[59] )*vsu10[69];

  vsu10[60] -= scalprod*vsu10[50];
  vsu10[61] -= scalprod*vsu10[51];
  vsu10[62] -= scalprod*vsu10[52];
  vsu10[63] -= scalprod*vsu10[53];
  vsu10[64] -= scalprod*vsu10[54];
  vsu10[65] -= scalprod*vsu10[55];
  vsu10[66] -= scalprod*vsu10[56];
  vsu10[67] -= scalprod*vsu10[57];
  vsu10[68] -= scalprod*vsu10[58];
  vsu10[69] -= scalprod*vsu10[59];


// Normalize the row `6'

  xn = 1./sqrt( norm(vsu10[60])
    + norm(vsu10[61])
    + norm(vsu10[62])
    + norm(vsu10[63])
    + norm(vsu10[64])
    + norm(vsu10[65])
    + norm(vsu10[66])
    + norm(vsu10[67])
    + norm(vsu10[68])
    + norm(vsu10[69]) );

    vsu10[60]*=xn;
    vsu10[61]*=xn;
    vsu10[62]*=xn;
    vsu10[63]*=xn;
    vsu10[64]*=xn;
    vsu10[65]*=xn;
    vsu10[66]*=xn;
    vsu10[67]*=xn;
    vsu10[68]*=xn;
    vsu10[69]*=xn;


// Make the row `7' orthogonal to the row `0':

  scalprod = conj( vsu10[0] )*vsu10[70]
    + conj( vsu10[1] )*vsu10[71]
    + conj( vsu10[2] )*vsu10[72]
    + conj( vsu10[3] )*vsu10[73]
    + conj( vsu10[4] )*vsu10[74]
    + conj( vsu10[5] )*vsu10[75]
    + conj( vsu10[6] )*vsu10[76]
    + conj( vsu10[7] )*vsu10[77]
    + conj( vsu10[8] )*vsu10[78]
    + conj( vsu10[9] )*vsu10[79];

  vsu10[70] -= scalprod*vsu10[0];
  vsu10[71] -= scalprod*vsu10[1];
  vsu10[72] -= scalprod*vsu10[2];
  vsu10[73] -= scalprod*vsu10[3];
  vsu10[74] -= scalprod*vsu10[4];
  vsu10[75] -= scalprod*vsu10[5];
  vsu10[76] -= scalprod*vsu10[6];
  vsu10[77] -= scalprod*vsu10[7];
  vsu10[78] -= scalprod*vsu10[8];
  vsu10[79] -= scalprod*vsu10[9];


// Make the row `7' orthogonal to the row `1':

  scalprod = conj( vsu10[10] )*vsu10[70]
    + conj( vsu10[11] )*vsu10[71]
    + conj( vsu10[12] )*vsu10[72]
    + conj( vsu10[13] )*vsu10[73]
    + conj( vsu10[14] )*vsu10[74]
    + conj( vsu10[15] )*vsu10[75]
    + conj( vsu10[16] )*vsu10[76]
    + conj( vsu10[17] )*vsu10[77]
    + conj( vsu10[18] )*vsu10[78]
    + conj( vsu10[19] )*vsu10[79];

  vsu10[70] -= scalprod*vsu10[10];
  vsu10[71] -= scalprod*vsu10[11];
  vsu10[72] -= scalprod*vsu10[12];
  vsu10[73] -= scalprod*vsu10[13];
  vsu10[74] -= scalprod*vsu10[14];
  vsu10[75] -= scalprod*vsu10[15];
  vsu10[76] -= scalprod*vsu10[16];
  vsu10[77] -= scalprod*vsu10[17];
  vsu10[78] -= scalprod*vsu10[18];
  vsu10[79] -= scalprod*vsu10[19];


// Make the row `7' orthogonal to the row `2':

  scalprod = conj( vsu10[20] )*vsu10[70]
    + conj( vsu10[21] )*vsu10[71]
    + conj( vsu10[22] )*vsu10[72]
    + conj( vsu10[23] )*vsu10[73]
    + conj( vsu10[24] )*vsu10[74]
    + conj( vsu10[25] )*vsu10[75]
    + conj( vsu10[26] )*vsu10[76]
    + conj( vsu10[27] )*vsu10[77]
    + conj( vsu10[28] )*vsu10[78]
    + conj( vsu10[29] )*vsu10[79];

  vsu10[70] -= scalprod*vsu10[20];
  vsu10[71] -= scalprod*vsu10[21];
  vsu10[72] -= scalprod*vsu10[22];
  vsu10[73] -= scalprod*vsu10[23];
  vsu10[74] -= scalprod*vsu10[24];
  vsu10[75] -= scalprod*vsu10[25];
  vsu10[76] -= scalprod*vsu10[26];
  vsu10[77] -= scalprod*vsu10[27];
  vsu10[78] -= scalprod*vsu10[28];
  vsu10[79] -= scalprod*vsu10[29];


// Make the row `7' orthogonal to the row `3':

  scalprod = conj( vsu10[30] )*vsu10[70]
    + conj( vsu10[31] )*vsu10[71]
    + conj( vsu10[32] )*vsu10[72]
    + conj( vsu10[33] )*vsu10[73]
    + conj( vsu10[34] )*vsu10[74]
    + conj( vsu10[35] )*vsu10[75]
    + conj( vsu10[36] )*vsu10[76]
    + conj( vsu10[37] )*vsu10[77]
    + conj( vsu10[38] )*vsu10[78]
    + conj( vsu10[39] )*vsu10[79];

  vsu10[70] -= scalprod*vsu10[30];
  vsu10[71] -= scalprod*vsu10[31];
  vsu10[72] -= scalprod*vsu10[32];
  vsu10[73] -= scalprod*vsu10[33];
  vsu10[74] -= scalprod*vsu10[34];
  vsu10[75] -= scalprod*vsu10[35];
  vsu10[76] -= scalprod*vsu10[36];
  vsu10[77] -= scalprod*vsu10[37];
  vsu10[78] -= scalprod*vsu10[38];
  vsu10[79] -= scalprod*vsu10[39];


// Make the row `7' orthogonal to the row `4':

  scalprod = conj( vsu10[40] )*vsu10[70]
    + conj( vsu10[41] )*vsu10[71]
    + conj( vsu10[42] )*vsu10[72]
    + conj( vsu10[43] )*vsu10[73]
    + conj( vsu10[44] )*vsu10[74]
    + conj( vsu10[45] )*vsu10[75]
    + conj( vsu10[46] )*vsu10[76]
    + conj( vsu10[47] )*vsu10[77]
    + conj( vsu10[48] )*vsu10[78]
    + conj( vsu10[49] )*vsu10[79];

  vsu10[70] -= scalprod*vsu10[40];
  vsu10[71] -= scalprod*vsu10[41];
  vsu10[72] -= scalprod*vsu10[42];
  vsu10[73] -= scalprod*vsu10[43];
  vsu10[74] -= scalprod*vsu10[44];
  vsu10[75] -= scalprod*vsu10[45];
  vsu10[76] -= scalprod*vsu10[46];
  vsu10[77] -= scalprod*vsu10[47];
  vsu10[78] -= scalprod*vsu10[48];
  vsu10[79] -= scalprod*vsu10[49];


// Make the row `7' orthogonal to the row `5':

  scalprod = conj( vsu10[50] )*vsu10[70]
    + conj( vsu10[51] )*vsu10[71]
    + conj( vsu10[52] )*vsu10[72]
    + conj( vsu10[53] )*vsu10[73]
    + conj( vsu10[54] )*vsu10[74]
    + conj( vsu10[55] )*vsu10[75]
    + conj( vsu10[56] )*vsu10[76]
    + conj( vsu10[57] )*vsu10[77]
    + conj( vsu10[58] )*vsu10[78]
    + conj( vsu10[59] )*vsu10[79];

  vsu10[70] -= scalprod*vsu10[50];
  vsu10[71] -= scalprod*vsu10[51];
  vsu10[72] -= scalprod*vsu10[52];
  vsu10[73] -= scalprod*vsu10[53];
  vsu10[74] -= scalprod*vsu10[54];
  vsu10[75] -= scalprod*vsu10[55];
  vsu10[76] -= scalprod*vsu10[56];
  vsu10[77] -= scalprod*vsu10[57];
  vsu10[78] -= scalprod*vsu10[58];
  vsu10[79] -= scalprod*vsu10[59];


// Make the row `7' orthogonal to the row `6':

  scalprod = conj( vsu10[60] )*vsu10[70]
    + conj( vsu10[61] )*vsu10[71]
    + conj( vsu10[62] )*vsu10[72]
    + conj( vsu10[63] )*vsu10[73]
    + conj( vsu10[64] )*vsu10[74]
    + conj( vsu10[65] )*vsu10[75]
    + conj( vsu10[66] )*vsu10[76]
    + conj( vsu10[67] )*vsu10[77]
    + conj( vsu10[68] )*vsu10[78]
    + conj( vsu10[69] )*vsu10[79];

  vsu10[70] -= scalprod*vsu10[60];
  vsu10[71] -= scalprod*vsu10[61];
  vsu10[72] -= scalprod*vsu10[62];
  vsu10[73] -= scalprod*vsu10[63];
  vsu10[74] -= scalprod*vsu10[64];
  vsu10[75] -= scalprod*vsu10[65];
  vsu10[76] -= scalprod*vsu10[66];
  vsu10[77] -= scalprod*vsu10[67];
  vsu10[78] -= scalprod*vsu10[68];
  vsu10[79] -= scalprod*vsu10[69];


// Normalize the row `7'

  xn = 1./sqrt( norm(vsu10[70])
    + norm(vsu10[71])
    + norm(vsu10[72])
    + norm(vsu10[73])
    + norm(vsu10[74])
    + norm(vsu10[75])
    + norm(vsu10[76])
    + norm(vsu10[77])
    + norm(vsu10[78])
    + norm(vsu10[79]) );

    vsu10[70]*=xn;
    vsu10[71]*=xn;
    vsu10[72]*=xn;
    vsu10[73]*=xn;
    vsu10[74]*=xn;
    vsu10[75]*=xn;
    vsu10[76]*=xn;
    vsu10[77]*=xn;
    vsu10[78]*=xn;
    vsu10[79]*=xn;


// Make the row `8' orthogonal to the row `0':

  scalprod = conj( vsu10[0] )*vsu10[80]
    + conj( vsu10[1] )*vsu10[81]
    + conj( vsu10[2] )*vsu10[82]
    + conj( vsu10[3] )*vsu10[83]
    + conj( vsu10[4] )*vsu10[84]
    + conj( vsu10[5] )*vsu10[85]
    + conj( vsu10[6] )*vsu10[86]
    + conj( vsu10[7] )*vsu10[87]
    + conj( vsu10[8] )*vsu10[88]
    + conj( vsu10[9] )*vsu10[89];

  vsu10[80] -= scalprod*vsu10[0];
  vsu10[81] -= scalprod*vsu10[1];
  vsu10[82] -= scalprod*vsu10[2];
  vsu10[83] -= scalprod*vsu10[3];
  vsu10[84] -= scalprod*vsu10[4];
  vsu10[85] -= scalprod*vsu10[5];
  vsu10[86] -= scalprod*vsu10[6];
  vsu10[87] -= scalprod*vsu10[7];
  vsu10[88] -= scalprod*vsu10[8];
  vsu10[89] -= scalprod*vsu10[9];


// Make the row `8' orthogonal to the row `1':

  scalprod = conj( vsu10[10] )*vsu10[80]
    + conj( vsu10[11] )*vsu10[81]
    + conj( vsu10[12] )*vsu10[82]
    + conj( vsu10[13] )*vsu10[83]
    + conj( vsu10[14] )*vsu10[84]
    + conj( vsu10[15] )*vsu10[85]
    + conj( vsu10[16] )*vsu10[86]
    + conj( vsu10[17] )*vsu10[87]
    + conj( vsu10[18] )*vsu10[88]
    + conj( vsu10[19] )*vsu10[89];

  vsu10[80] -= scalprod*vsu10[10];
  vsu10[81] -= scalprod*vsu10[11];
  vsu10[82] -= scalprod*vsu10[12];
  vsu10[83] -= scalprod*vsu10[13];
  vsu10[84] -= scalprod*vsu10[14];
  vsu10[85] -= scalprod*vsu10[15];
  vsu10[86] -= scalprod*vsu10[16];
  vsu10[87] -= scalprod*vsu10[17];
  vsu10[88] -= scalprod*vsu10[18];
  vsu10[89] -= scalprod*vsu10[19];


// Make the row `8' orthogonal to the row `2':

  scalprod = conj( vsu10[20] )*vsu10[80]
    + conj( vsu10[21] )*vsu10[81]
    + conj( vsu10[22] )*vsu10[82]
    + conj( vsu10[23] )*vsu10[83]
    + conj( vsu10[24] )*vsu10[84]
    + conj( vsu10[25] )*vsu10[85]
    + conj( vsu10[26] )*vsu10[86]
    + conj( vsu10[27] )*vsu10[87]
    + conj( vsu10[28] )*vsu10[88]
    + conj( vsu10[29] )*vsu10[89];

  vsu10[80] -= scalprod*vsu10[20];
  vsu10[81] -= scalprod*vsu10[21];
  vsu10[82] -= scalprod*vsu10[22];
  vsu10[83] -= scalprod*vsu10[23];
  vsu10[84] -= scalprod*vsu10[24];
  vsu10[85] -= scalprod*vsu10[25];
  vsu10[86] -= scalprod*vsu10[26];
  vsu10[87] -= scalprod*vsu10[27];
  vsu10[88] -= scalprod*vsu10[28];
  vsu10[89] -= scalprod*vsu10[29];


// Make the row `8' orthogonal to the row `3':

  scalprod = conj( vsu10[30] )*vsu10[80]
    + conj( vsu10[31] )*vsu10[81]
    + conj( vsu10[32] )*vsu10[82]
    + conj( vsu10[33] )*vsu10[83]
    + conj( vsu10[34] )*vsu10[84]
    + conj( vsu10[35] )*vsu10[85]
    + conj( vsu10[36] )*vsu10[86]
    + conj( vsu10[37] )*vsu10[87]
    + conj( vsu10[38] )*vsu10[88]
    + conj( vsu10[39] )*vsu10[89];

  vsu10[80] -= scalprod*vsu10[30];
  vsu10[81] -= scalprod*vsu10[31];
  vsu10[82] -= scalprod*vsu10[32];
  vsu10[83] -= scalprod*vsu10[33];
  vsu10[84] -= scalprod*vsu10[34];
  vsu10[85] -= scalprod*vsu10[35];
  vsu10[86] -= scalprod*vsu10[36];
  vsu10[87] -= scalprod*vsu10[37];
  vsu10[88] -= scalprod*vsu10[38];
  vsu10[89] -= scalprod*vsu10[39];


// Make the row `8' orthogonal to the row `4':

  scalprod = conj( vsu10[40] )*vsu10[80]
    + conj( vsu10[41] )*vsu10[81]
    + conj( vsu10[42] )*vsu10[82]
    + conj( vsu10[43] )*vsu10[83]
    + conj( vsu10[44] )*vsu10[84]
    + conj( vsu10[45] )*vsu10[85]
    + conj( vsu10[46] )*vsu10[86]
    + conj( vsu10[47] )*vsu10[87]
    + conj( vsu10[48] )*vsu10[88]
    + conj( vsu10[49] )*vsu10[89];

  vsu10[80] -= scalprod*vsu10[40];
  vsu10[81] -= scalprod*vsu10[41];
  vsu10[82] -= scalprod*vsu10[42];
  vsu10[83] -= scalprod*vsu10[43];
  vsu10[84] -= scalprod*vsu10[44];
  vsu10[85] -= scalprod*vsu10[45];
  vsu10[86] -= scalprod*vsu10[46];
  vsu10[87] -= scalprod*vsu10[47];
  vsu10[88] -= scalprod*vsu10[48];
  vsu10[89] -= scalprod*vsu10[49];


// Make the row `8' orthogonal to the row `5':

  scalprod = conj( vsu10[50] )*vsu10[80]
    + conj( vsu10[51] )*vsu10[81]
    + conj( vsu10[52] )*vsu10[82]
    + conj( vsu10[53] )*vsu10[83]
    + conj( vsu10[54] )*vsu10[84]
    + conj( vsu10[55] )*vsu10[85]
    + conj( vsu10[56] )*vsu10[86]
    + conj( vsu10[57] )*vsu10[87]
    + conj( vsu10[58] )*vsu10[88]
    + conj( vsu10[59] )*vsu10[89];

  vsu10[80] -= scalprod*vsu10[50];
  vsu10[81] -= scalprod*vsu10[51];
  vsu10[82] -= scalprod*vsu10[52];
  vsu10[83] -= scalprod*vsu10[53];
  vsu10[84] -= scalprod*vsu10[54];
  vsu10[85] -= scalprod*vsu10[55];
  vsu10[86] -= scalprod*vsu10[56];
  vsu10[87] -= scalprod*vsu10[57];
  vsu10[88] -= scalprod*vsu10[58];
  vsu10[89] -= scalprod*vsu10[59];


// Make the row `8' orthogonal to the row `6':

  scalprod = conj( vsu10[60] )*vsu10[80]
    + conj( vsu10[61] )*vsu10[81]
    + conj( vsu10[62] )*vsu10[82]
    + conj( vsu10[63] )*vsu10[83]
    + conj( vsu10[64] )*vsu10[84]
    + conj( vsu10[65] )*vsu10[85]
    + conj( vsu10[66] )*vsu10[86]
    + conj( vsu10[67] )*vsu10[87]
    + conj( vsu10[68] )*vsu10[88]
    + conj( vsu10[69] )*vsu10[89];

  vsu10[80] -= scalprod*vsu10[60];
  vsu10[81] -= scalprod*vsu10[61];
  vsu10[82] -= scalprod*vsu10[62];
  vsu10[83] -= scalprod*vsu10[63];
  vsu10[84] -= scalprod*vsu10[64];
  vsu10[85] -= scalprod*vsu10[65];
  vsu10[86] -= scalprod*vsu10[66];
  vsu10[87] -= scalprod*vsu10[67];
  vsu10[88] -= scalprod*vsu10[68];
  vsu10[89] -= scalprod*vsu10[69];


// Make the row `8' orthogonal to the row `7':

  scalprod = conj( vsu10[70] )*vsu10[80]
    + conj( vsu10[71] )*vsu10[81]
    + conj( vsu10[72] )*vsu10[82]
    + conj( vsu10[73] )*vsu10[83]
    + conj( vsu10[74] )*vsu10[84]
    + conj( vsu10[75] )*vsu10[85]
    + conj( vsu10[76] )*vsu10[86]
    + conj( vsu10[77] )*vsu10[87]
    + conj( vsu10[78] )*vsu10[88]
    + conj( vsu10[79] )*vsu10[89];

  vsu10[80] -= scalprod*vsu10[70];
  vsu10[81] -= scalprod*vsu10[71];
  vsu10[82] -= scalprod*vsu10[72];
  vsu10[83] -= scalprod*vsu10[73];
  vsu10[84] -= scalprod*vsu10[74];
  vsu10[85] -= scalprod*vsu10[75];
  vsu10[86] -= scalprod*vsu10[76];
  vsu10[87] -= scalprod*vsu10[77];
  vsu10[88] -= scalprod*vsu10[78];
  vsu10[89] -= scalprod*vsu10[79];


// Normalize the row `8'

  xn = 1./sqrt( norm(vsu10[80])
    + norm(vsu10[81])
    + norm(vsu10[82])
    + norm(vsu10[83])
    + norm(vsu10[84])
    + norm(vsu10[85])
    + norm(vsu10[86])
    + norm(vsu10[87])
    + norm(vsu10[88])
    + norm(vsu10[89]) );

    vsu10[80]*=xn;
    vsu10[81]*=xn;
    vsu10[82]*=xn;
    vsu10[83]*=xn;
    vsu10[84]*=xn;
    vsu10[85]*=xn;
    vsu10[86]*=xn;
    vsu10[87]*=xn;
    vsu10[88]*=xn;
    vsu10[89]*=xn;


// Make the row `9' orthogonal to the row `0':

  scalprod = conj( vsu10[0] )*vsu10[90]
    + conj( vsu10[1] )*vsu10[91]
    + conj( vsu10[2] )*vsu10[92]
    + conj( vsu10[3] )*vsu10[93]
    + conj( vsu10[4] )*vsu10[94]
    + conj( vsu10[5] )*vsu10[95]
    + conj( vsu10[6] )*vsu10[96]
    + conj( vsu10[7] )*vsu10[97]
    + conj( vsu10[8] )*vsu10[98]
    + conj( vsu10[9] )*vsu10[99];

  vsu10[90] -= scalprod*vsu10[0];
  vsu10[91] -= scalprod*vsu10[1];
  vsu10[92] -= scalprod*vsu10[2];
  vsu10[93] -= scalprod*vsu10[3];
  vsu10[94] -= scalprod*vsu10[4];
  vsu10[95] -= scalprod*vsu10[5];
  vsu10[96] -= scalprod*vsu10[6];
  vsu10[97] -= scalprod*vsu10[7];
  vsu10[98] -= scalprod*vsu10[8];
  vsu10[99] -= scalprod*vsu10[9];


// Make the row `9' orthogonal to the row `1':

  scalprod = conj( vsu10[10] )*vsu10[90]
    + conj( vsu10[11] )*vsu10[91]
    + conj( vsu10[12] )*vsu10[92]
    + conj( vsu10[13] )*vsu10[93]
    + conj( vsu10[14] )*vsu10[94]
    + conj( vsu10[15] )*vsu10[95]
    + conj( vsu10[16] )*vsu10[96]
    + conj( vsu10[17] )*vsu10[97]
    + conj( vsu10[18] )*vsu10[98]
    + conj( vsu10[19] )*vsu10[99];

  vsu10[90] -= scalprod*vsu10[10];
  vsu10[91] -= scalprod*vsu10[11];
  vsu10[92] -= scalprod*vsu10[12];
  vsu10[93] -= scalprod*vsu10[13];
  vsu10[94] -= scalprod*vsu10[14];
  vsu10[95] -= scalprod*vsu10[15];
  vsu10[96] -= scalprod*vsu10[16];
  vsu10[97] -= scalprod*vsu10[17];
  vsu10[98] -= scalprod*vsu10[18];
  vsu10[99] -= scalprod*vsu10[19];


// Make the row `9' orthogonal to the row `2':

  scalprod = conj( vsu10[20] )*vsu10[90]
    + conj( vsu10[21] )*vsu10[91]
    + conj( vsu10[22] )*vsu10[92]
    + conj( vsu10[23] )*vsu10[93]
    + conj( vsu10[24] )*vsu10[94]
    + conj( vsu10[25] )*vsu10[95]
    + conj( vsu10[26] )*vsu10[96]
    + conj( vsu10[27] )*vsu10[97]
    + conj( vsu10[28] )*vsu10[98]
    + conj( vsu10[29] )*vsu10[99];

  vsu10[90] -= scalprod*vsu10[20];
  vsu10[91] -= scalprod*vsu10[21];
  vsu10[92] -= scalprod*vsu10[22];
  vsu10[93] -= scalprod*vsu10[23];
  vsu10[94] -= scalprod*vsu10[24];
  vsu10[95] -= scalprod*vsu10[25];
  vsu10[96] -= scalprod*vsu10[26];
  vsu10[97] -= scalprod*vsu10[27];
  vsu10[98] -= scalprod*vsu10[28];
  vsu10[99] -= scalprod*vsu10[29];


// Make the row `9' orthogonal to the row `3':

  scalprod = conj( vsu10[30] )*vsu10[90]
    + conj( vsu10[31] )*vsu10[91]
    + conj( vsu10[32] )*vsu10[92]
    + conj( vsu10[33] )*vsu10[93]
    + conj( vsu10[34] )*vsu10[94]
    + conj( vsu10[35] )*vsu10[95]
    + conj( vsu10[36] )*vsu10[96]
    + conj( vsu10[37] )*vsu10[97]
    + conj( vsu10[38] )*vsu10[98]
    + conj( vsu10[39] )*vsu10[99];

  vsu10[90] -= scalprod*vsu10[30];
  vsu10[91] -= scalprod*vsu10[31];
  vsu10[92] -= scalprod*vsu10[32];
  vsu10[93] -= scalprod*vsu10[33];
  vsu10[94] -= scalprod*vsu10[34];
  vsu10[95] -= scalprod*vsu10[35];
  vsu10[96] -= scalprod*vsu10[36];
  vsu10[97] -= scalprod*vsu10[37];
  vsu10[98] -= scalprod*vsu10[38];
  vsu10[99] -= scalprod*vsu10[39];


// Make the row `9' orthogonal to the row `4':

  scalprod = conj( vsu10[40] )*vsu10[90]
    + conj( vsu10[41] )*vsu10[91]
    + conj( vsu10[42] )*vsu10[92]
    + conj( vsu10[43] )*vsu10[93]
    + conj( vsu10[44] )*vsu10[94]
    + conj( vsu10[45] )*vsu10[95]
    + conj( vsu10[46] )*vsu10[96]
    + conj( vsu10[47] )*vsu10[97]
    + conj( vsu10[48] )*vsu10[98]
    + conj( vsu10[49] )*vsu10[99];

  vsu10[90] -= scalprod*vsu10[40];
  vsu10[91] -= scalprod*vsu10[41];
  vsu10[92] -= scalprod*vsu10[42];
  vsu10[93] -= scalprod*vsu10[43];
  vsu10[94] -= scalprod*vsu10[44];
  vsu10[95] -= scalprod*vsu10[45];
  vsu10[96] -= scalprod*vsu10[46];
  vsu10[97] -= scalprod*vsu10[47];
  vsu10[98] -= scalprod*vsu10[48];
  vsu10[99] -= scalprod*vsu10[49];


// Make the row `9' orthogonal to the row `5':

  scalprod = conj( vsu10[50] )*vsu10[90]
    + conj( vsu10[51] )*vsu10[91]
    + conj( vsu10[52] )*vsu10[92]
    + conj( vsu10[53] )*vsu10[93]
    + conj( vsu10[54] )*vsu10[94]
    + conj( vsu10[55] )*vsu10[95]
    + conj( vsu10[56] )*vsu10[96]
    + conj( vsu10[57] )*vsu10[97]
    + conj( vsu10[58] )*vsu10[98]
    + conj( vsu10[59] )*vsu10[99];

  vsu10[90] -= scalprod*vsu10[50];
  vsu10[91] -= scalprod*vsu10[51];
  vsu10[92] -= scalprod*vsu10[52];
  vsu10[93] -= scalprod*vsu10[53];
  vsu10[94] -= scalprod*vsu10[54];
  vsu10[95] -= scalprod*vsu10[55];
  vsu10[96] -= scalprod*vsu10[56];
  vsu10[97] -= scalprod*vsu10[57];
  vsu10[98] -= scalprod*vsu10[58];
  vsu10[99] -= scalprod*vsu10[59];


// Make the row `9' orthogonal to the row `6':

  scalprod = conj( vsu10[60] )*vsu10[90]
    + conj( vsu10[61] )*vsu10[91]
    + conj( vsu10[62] )*vsu10[92]
    + conj( vsu10[63] )*vsu10[93]
    + conj( vsu10[64] )*vsu10[94]
    + conj( vsu10[65] )*vsu10[95]
    + conj( vsu10[66] )*vsu10[96]
    + conj( vsu10[67] )*vsu10[97]
    + conj( vsu10[68] )*vsu10[98]
    + conj( vsu10[69] )*vsu10[99];

  vsu10[90] -= scalprod*vsu10[60];
  vsu10[91] -= scalprod*vsu10[61];
  vsu10[92] -= scalprod*vsu10[62];
  vsu10[93] -= scalprod*vsu10[63];
  vsu10[94] -= scalprod*vsu10[64];
  vsu10[95] -= scalprod*vsu10[65];
  vsu10[96] -= scalprod*vsu10[66];
  vsu10[97] -= scalprod*vsu10[67];
  vsu10[98] -= scalprod*vsu10[68];
  vsu10[99] -= scalprod*vsu10[69];


// Make the row `9' orthogonal to the row `7':

  scalprod = conj( vsu10[70] )*vsu10[90]
    + conj( vsu10[71] )*vsu10[91]
    + conj( vsu10[72] )*vsu10[92]
    + conj( vsu10[73] )*vsu10[93]
    + conj( vsu10[74] )*vsu10[94]
    + conj( vsu10[75] )*vsu10[95]
    + conj( vsu10[76] )*vsu10[96]
    + conj( vsu10[77] )*vsu10[97]
    + conj( vsu10[78] )*vsu10[98]
    + conj( vsu10[79] )*vsu10[99];

  vsu10[90] -= scalprod*vsu10[70];
  vsu10[91] -= scalprod*vsu10[71];
  vsu10[92] -= scalprod*vsu10[72];
  vsu10[93] -= scalprod*vsu10[73];
  vsu10[94] -= scalprod*vsu10[74];
  vsu10[95] -= scalprod*vsu10[75];
  vsu10[96] -= scalprod*vsu10[76];
  vsu10[97] -= scalprod*vsu10[77];
  vsu10[98] -= scalprod*vsu10[78];
  vsu10[99] -= scalprod*vsu10[79];


// Make the row `9' orthogonal to the row `8':

  scalprod = conj( vsu10[80] )*vsu10[90]
    + conj( vsu10[81] )*vsu10[91]
    + conj( vsu10[82] )*vsu10[92]
    + conj( vsu10[83] )*vsu10[93]
    + conj( vsu10[84] )*vsu10[94]
    + conj( vsu10[85] )*vsu10[95]
    + conj( vsu10[86] )*vsu10[96]
    + conj( vsu10[87] )*vsu10[97]
    + conj( vsu10[88] )*vsu10[98]
    + conj( vsu10[89] )*vsu10[99];

  vsu10[90] -= scalprod*vsu10[80];
  vsu10[91] -= scalprod*vsu10[81];
  vsu10[92] -= scalprod*vsu10[82];
  vsu10[93] -= scalprod*vsu10[83];
  vsu10[94] -= scalprod*vsu10[84];
  vsu10[95] -= scalprod*vsu10[85];
  vsu10[96] -= scalprod*vsu10[86];
  vsu10[97] -= scalprod*vsu10[87];
  vsu10[98] -= scalprod*vsu10[88];
  vsu10[99] -= scalprod*vsu10[89];


// Normalize the row `9'

  xn = 1./sqrt( norm(vsu10[90])
    + norm(vsu10[91])
    + norm(vsu10[92])
    + norm(vsu10[93])
    + norm(vsu10[94])
    + norm(vsu10[95])
    + norm(vsu10[96])
    + norm(vsu10[97])
    + norm(vsu10[98])
    + norm(vsu10[99]) );

    vsu10[90]*=xn;
    vsu10[91]*=xn;
    vsu10[92]*=xn;
    vsu10[93]*=xn;
    vsu10[94]*=xn;
    vsu10[95]*=xn;
    vsu10[96]*=xn;
    vsu10[97]*=xn;
    vsu10[98]*=xn;
    vsu10[99]*=xn;


// Impose unimodularity by arranging the phase of the last row
// to compensate for the phase of the determinant

  for (i=0; i<Ncol; i++)
  for (j=0; j<Ncol; j++) {
    AT[2*(j+Ncol*i)]=real(vsu10[j*Ncol+i]);
    AT[2*(j+Ncol*i)+1]=imag(vsu10[j*Ncol+i]);
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

  vsu10[90] *= determinant_conjugate_phase;
  vsu10[91] *= determinant_conjugate_phase;
  vsu10[92] *= determinant_conjugate_phase;
  vsu10[93] *= determinant_conjugate_phase;
  vsu10[94] *= determinant_conjugate_phase;
  vsu10[95] *= determinant_conjugate_phase;
  vsu10[96] *= determinant_conjugate_phase;
  vsu10[97] *= determinant_conjugate_phase;
  vsu10[98] *= determinant_conjugate_phase;
  vsu10[99] *= determinant_conjugate_phase;


}

#endif
