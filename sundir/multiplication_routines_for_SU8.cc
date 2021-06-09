inline void mult_C_equals_AB_for_SU8(dc *C, dc *A, dc *B) {

  C[0] = A[0]*B[0]
        +A[1]*B[8]
        +A[2]*B[16]
        +A[3]*B[24]
        +A[4]*B[32]
        +A[5]*B[40]
        +A[6]*B[48]
        +A[7]*B[56];

  C[1] = A[0]*B[1]
        +A[1]*B[9]
        +A[2]*B[17]
        +A[3]*B[25]
        +A[4]*B[33]
        +A[5]*B[41]
        +A[6]*B[49]
        +A[7]*B[57];

  C[2] = A[0]*B[2]
        +A[1]*B[10]
        +A[2]*B[18]
        +A[3]*B[26]
        +A[4]*B[34]
        +A[5]*B[42]
        +A[6]*B[50]
        +A[7]*B[58];

  C[3] = A[0]*B[3]
        +A[1]*B[11]
        +A[2]*B[19]
        +A[3]*B[27]
        +A[4]*B[35]
        +A[5]*B[43]
        +A[6]*B[51]
        +A[7]*B[59];

  C[4] = A[0]*B[4]
        +A[1]*B[12]
        +A[2]*B[20]
        +A[3]*B[28]
        +A[4]*B[36]
        +A[5]*B[44]
        +A[6]*B[52]
        +A[7]*B[60];

  C[5] = A[0]*B[5]
        +A[1]*B[13]
        +A[2]*B[21]
        +A[3]*B[29]
        +A[4]*B[37]
        +A[5]*B[45]
        +A[6]*B[53]
        +A[7]*B[61];

  C[6] = A[0]*B[6]
        +A[1]*B[14]
        +A[2]*B[22]
        +A[3]*B[30]
        +A[4]*B[38]
        +A[5]*B[46]
        +A[6]*B[54]
        +A[7]*B[62];

  C[7] = A[0]*B[7]
        +A[1]*B[15]
        +A[2]*B[23]
        +A[3]*B[31]
        +A[4]*B[39]
        +A[5]*B[47]
        +A[6]*B[55]
        +A[7]*B[63];

  C[8] = A[8]*B[0]
        +A[9]*B[8]
        +A[10]*B[16]
        +A[11]*B[24]
        +A[12]*B[32]
        +A[13]*B[40]
        +A[14]*B[48]
        +A[15]*B[56];

  C[9] = A[8]*B[1]
        +A[9]*B[9]
        +A[10]*B[17]
        +A[11]*B[25]
        +A[12]*B[33]
        +A[13]*B[41]
        +A[14]*B[49]
        +A[15]*B[57];

  C[10] = A[8]*B[2]
        +A[9]*B[10]
        +A[10]*B[18]
        +A[11]*B[26]
        +A[12]*B[34]
        +A[13]*B[42]
        +A[14]*B[50]
        +A[15]*B[58];

  C[11] = A[8]*B[3]
        +A[9]*B[11]
        +A[10]*B[19]
        +A[11]*B[27]
        +A[12]*B[35]
        +A[13]*B[43]
        +A[14]*B[51]
        +A[15]*B[59];

  C[12] = A[8]*B[4]
        +A[9]*B[12]
        +A[10]*B[20]
        +A[11]*B[28]
        +A[12]*B[36]
        +A[13]*B[44]
        +A[14]*B[52]
        +A[15]*B[60];

  C[13] = A[8]*B[5]
        +A[9]*B[13]
        +A[10]*B[21]
        +A[11]*B[29]
        +A[12]*B[37]
        +A[13]*B[45]
        +A[14]*B[53]
        +A[15]*B[61];

  C[14] = A[8]*B[6]
        +A[9]*B[14]
        +A[10]*B[22]
        +A[11]*B[30]
        +A[12]*B[38]
        +A[13]*B[46]
        +A[14]*B[54]
        +A[15]*B[62];

  C[15] = A[8]*B[7]
        +A[9]*B[15]
        +A[10]*B[23]
        +A[11]*B[31]
        +A[12]*B[39]
        +A[13]*B[47]
        +A[14]*B[55]
        +A[15]*B[63];

  C[16] = A[16]*B[0]
        +A[17]*B[8]
        +A[18]*B[16]
        +A[19]*B[24]
        +A[20]*B[32]
        +A[21]*B[40]
        +A[22]*B[48]
        +A[23]*B[56];

  C[17] = A[16]*B[1]
        +A[17]*B[9]
        +A[18]*B[17]
        +A[19]*B[25]
        +A[20]*B[33]
        +A[21]*B[41]
        +A[22]*B[49]
        +A[23]*B[57];

  C[18] = A[16]*B[2]
        +A[17]*B[10]
        +A[18]*B[18]
        +A[19]*B[26]
        +A[20]*B[34]
        +A[21]*B[42]
        +A[22]*B[50]
        +A[23]*B[58];

  C[19] = A[16]*B[3]
        +A[17]*B[11]
        +A[18]*B[19]
        +A[19]*B[27]
        +A[20]*B[35]
        +A[21]*B[43]
        +A[22]*B[51]
        +A[23]*B[59];

  C[20] = A[16]*B[4]
        +A[17]*B[12]
        +A[18]*B[20]
        +A[19]*B[28]
        +A[20]*B[36]
        +A[21]*B[44]
        +A[22]*B[52]
        +A[23]*B[60];

  C[21] = A[16]*B[5]
        +A[17]*B[13]
        +A[18]*B[21]
        +A[19]*B[29]
        +A[20]*B[37]
        +A[21]*B[45]
        +A[22]*B[53]
        +A[23]*B[61];

  C[22] = A[16]*B[6]
        +A[17]*B[14]
        +A[18]*B[22]
        +A[19]*B[30]
        +A[20]*B[38]
        +A[21]*B[46]
        +A[22]*B[54]
        +A[23]*B[62];

  C[23] = A[16]*B[7]
        +A[17]*B[15]
        +A[18]*B[23]
        +A[19]*B[31]
        +A[20]*B[39]
        +A[21]*B[47]
        +A[22]*B[55]
        +A[23]*B[63];

  C[24] = A[24]*B[0]
        +A[25]*B[8]
        +A[26]*B[16]
        +A[27]*B[24]
        +A[28]*B[32]
        +A[29]*B[40]
        +A[30]*B[48]
        +A[31]*B[56];

  C[25] = A[24]*B[1]
        +A[25]*B[9]
        +A[26]*B[17]
        +A[27]*B[25]
        +A[28]*B[33]
        +A[29]*B[41]
        +A[30]*B[49]
        +A[31]*B[57];

  C[26] = A[24]*B[2]
        +A[25]*B[10]
        +A[26]*B[18]
        +A[27]*B[26]
        +A[28]*B[34]
        +A[29]*B[42]
        +A[30]*B[50]
        +A[31]*B[58];

  C[27] = A[24]*B[3]
        +A[25]*B[11]
        +A[26]*B[19]
        +A[27]*B[27]
        +A[28]*B[35]
        +A[29]*B[43]
        +A[30]*B[51]
        +A[31]*B[59];

  C[28] = A[24]*B[4]
        +A[25]*B[12]
        +A[26]*B[20]
        +A[27]*B[28]
        +A[28]*B[36]
        +A[29]*B[44]
        +A[30]*B[52]
        +A[31]*B[60];

  C[29] = A[24]*B[5]
        +A[25]*B[13]
        +A[26]*B[21]
        +A[27]*B[29]
        +A[28]*B[37]
        +A[29]*B[45]
        +A[30]*B[53]
        +A[31]*B[61];

  C[30] = A[24]*B[6]
        +A[25]*B[14]
        +A[26]*B[22]
        +A[27]*B[30]
        +A[28]*B[38]
        +A[29]*B[46]
        +A[30]*B[54]
        +A[31]*B[62];

  C[31] = A[24]*B[7]
        +A[25]*B[15]
        +A[26]*B[23]
        +A[27]*B[31]
        +A[28]*B[39]
        +A[29]*B[47]
        +A[30]*B[55]
        +A[31]*B[63];

  C[32] = A[32]*B[0]
        +A[33]*B[8]
        +A[34]*B[16]
        +A[35]*B[24]
        +A[36]*B[32]
        +A[37]*B[40]
        +A[38]*B[48]
        +A[39]*B[56];

  C[33] = A[32]*B[1]
        +A[33]*B[9]
        +A[34]*B[17]
        +A[35]*B[25]
        +A[36]*B[33]
        +A[37]*B[41]
        +A[38]*B[49]
        +A[39]*B[57];

  C[34] = A[32]*B[2]
        +A[33]*B[10]
        +A[34]*B[18]
        +A[35]*B[26]
        +A[36]*B[34]
        +A[37]*B[42]
        +A[38]*B[50]
        +A[39]*B[58];

  C[35] = A[32]*B[3]
        +A[33]*B[11]
        +A[34]*B[19]
        +A[35]*B[27]
        +A[36]*B[35]
        +A[37]*B[43]
        +A[38]*B[51]
        +A[39]*B[59];

  C[36] = A[32]*B[4]
        +A[33]*B[12]
        +A[34]*B[20]
        +A[35]*B[28]
        +A[36]*B[36]
        +A[37]*B[44]
        +A[38]*B[52]
        +A[39]*B[60];

  C[37] = A[32]*B[5]
        +A[33]*B[13]
        +A[34]*B[21]
        +A[35]*B[29]
        +A[36]*B[37]
        +A[37]*B[45]
        +A[38]*B[53]
        +A[39]*B[61];

  C[38] = A[32]*B[6]
        +A[33]*B[14]
        +A[34]*B[22]
        +A[35]*B[30]
        +A[36]*B[38]
        +A[37]*B[46]
        +A[38]*B[54]
        +A[39]*B[62];

  C[39] = A[32]*B[7]
        +A[33]*B[15]
        +A[34]*B[23]
        +A[35]*B[31]
        +A[36]*B[39]
        +A[37]*B[47]
        +A[38]*B[55]
        +A[39]*B[63];

  C[40] = A[40]*B[0]
        +A[41]*B[8]
        +A[42]*B[16]
        +A[43]*B[24]
        +A[44]*B[32]
        +A[45]*B[40]
        +A[46]*B[48]
        +A[47]*B[56];

  C[41] = A[40]*B[1]
        +A[41]*B[9]
        +A[42]*B[17]
        +A[43]*B[25]
        +A[44]*B[33]
        +A[45]*B[41]
        +A[46]*B[49]
        +A[47]*B[57];

  C[42] = A[40]*B[2]
        +A[41]*B[10]
        +A[42]*B[18]
        +A[43]*B[26]
        +A[44]*B[34]
        +A[45]*B[42]
        +A[46]*B[50]
        +A[47]*B[58];

  C[43] = A[40]*B[3]
        +A[41]*B[11]
        +A[42]*B[19]
        +A[43]*B[27]
        +A[44]*B[35]
        +A[45]*B[43]
        +A[46]*B[51]
        +A[47]*B[59];

  C[44] = A[40]*B[4]
        +A[41]*B[12]
        +A[42]*B[20]
        +A[43]*B[28]
        +A[44]*B[36]
        +A[45]*B[44]
        +A[46]*B[52]
        +A[47]*B[60];

  C[45] = A[40]*B[5]
        +A[41]*B[13]
        +A[42]*B[21]
        +A[43]*B[29]
        +A[44]*B[37]
        +A[45]*B[45]
        +A[46]*B[53]
        +A[47]*B[61];

  C[46] = A[40]*B[6]
        +A[41]*B[14]
        +A[42]*B[22]
        +A[43]*B[30]
        +A[44]*B[38]
        +A[45]*B[46]
        +A[46]*B[54]
        +A[47]*B[62];

  C[47] = A[40]*B[7]
        +A[41]*B[15]
        +A[42]*B[23]
        +A[43]*B[31]
        +A[44]*B[39]
        +A[45]*B[47]
        +A[46]*B[55]
        +A[47]*B[63];

  C[48] = A[48]*B[0]
        +A[49]*B[8]
        +A[50]*B[16]
        +A[51]*B[24]
        +A[52]*B[32]
        +A[53]*B[40]
        +A[54]*B[48]
        +A[55]*B[56];

  C[49] = A[48]*B[1]
        +A[49]*B[9]
        +A[50]*B[17]
        +A[51]*B[25]
        +A[52]*B[33]
        +A[53]*B[41]
        +A[54]*B[49]
        +A[55]*B[57];

  C[50] = A[48]*B[2]
        +A[49]*B[10]
        +A[50]*B[18]
        +A[51]*B[26]
        +A[52]*B[34]
        +A[53]*B[42]
        +A[54]*B[50]
        +A[55]*B[58];

  C[51] = A[48]*B[3]
        +A[49]*B[11]
        +A[50]*B[19]
        +A[51]*B[27]
        +A[52]*B[35]
        +A[53]*B[43]
        +A[54]*B[51]
        +A[55]*B[59];

  C[52] = A[48]*B[4]
        +A[49]*B[12]
        +A[50]*B[20]
        +A[51]*B[28]
        +A[52]*B[36]
        +A[53]*B[44]
        +A[54]*B[52]
        +A[55]*B[60];

  C[53] = A[48]*B[5]
        +A[49]*B[13]
        +A[50]*B[21]
        +A[51]*B[29]
        +A[52]*B[37]
        +A[53]*B[45]
        +A[54]*B[53]
        +A[55]*B[61];

  C[54] = A[48]*B[6]
        +A[49]*B[14]
        +A[50]*B[22]
        +A[51]*B[30]
        +A[52]*B[38]
        +A[53]*B[46]
        +A[54]*B[54]
        +A[55]*B[62];

  C[55] = A[48]*B[7]
        +A[49]*B[15]
        +A[50]*B[23]
        +A[51]*B[31]
        +A[52]*B[39]
        +A[53]*B[47]
        +A[54]*B[55]
        +A[55]*B[63];

  C[56] = A[56]*B[0]
        +A[57]*B[8]
        +A[58]*B[16]
        +A[59]*B[24]
        +A[60]*B[32]
        +A[61]*B[40]
        +A[62]*B[48]
        +A[63]*B[56];

  C[57] = A[56]*B[1]
        +A[57]*B[9]
        +A[58]*B[17]
        +A[59]*B[25]
        +A[60]*B[33]
        +A[61]*B[41]
        +A[62]*B[49]
        +A[63]*B[57];

  C[58] = A[56]*B[2]
        +A[57]*B[10]
        +A[58]*B[18]
        +A[59]*B[26]
        +A[60]*B[34]
        +A[61]*B[42]
        +A[62]*B[50]
        +A[63]*B[58];

  C[59] = A[56]*B[3]
        +A[57]*B[11]
        +A[58]*B[19]
        +A[59]*B[27]
        +A[60]*B[35]
        +A[61]*B[43]
        +A[62]*B[51]
        +A[63]*B[59];

  C[60] = A[56]*B[4]
        +A[57]*B[12]
        +A[58]*B[20]
        +A[59]*B[28]
        +A[60]*B[36]
        +A[61]*B[44]
        +A[62]*B[52]
        +A[63]*B[60];

  C[61] = A[56]*B[5]
        +A[57]*B[13]
        +A[58]*B[21]
        +A[59]*B[29]
        +A[60]*B[37]
        +A[61]*B[45]
        +A[62]*B[53]
        +A[63]*B[61];

  C[62] = A[56]*B[6]
        +A[57]*B[14]
        +A[58]*B[22]
        +A[59]*B[30]
        +A[60]*B[38]
        +A[61]*B[46]
        +A[62]*B[54]
        +A[63]*B[62];

  C[63] = A[56]*B[7]
        +A[57]*B[15]
        +A[58]*B[23]
        +A[59]*B[31]
        +A[60]*B[39]
        +A[61]*B[47]
        +A[62]*B[55]
        +A[63]*B[63];

};


inline void mult_C_equals_ABdagger_for_SU8(dc *C, dc *A, dc *B) {

  C[0] = A[0]*conj(B[0])
        +A[1]*conj(B[1])
        +A[2]*conj(B[2])
        +A[3]*conj(B[3])
        +A[4]*conj(B[4])
        +A[5]*conj(B[5])
        +A[6]*conj(B[6])
        +A[7]*conj(B[7]);

  C[1] = A[0]*conj(B[8])
        +A[1]*conj(B[9])
        +A[2]*conj(B[10])
        +A[3]*conj(B[11])
        +A[4]*conj(B[12])
        +A[5]*conj(B[13])
        +A[6]*conj(B[14])
        +A[7]*conj(B[15]);

  C[2] = A[0]*conj(B[16])
        +A[1]*conj(B[17])
        +A[2]*conj(B[18])
        +A[3]*conj(B[19])
        +A[4]*conj(B[20])
        +A[5]*conj(B[21])
        +A[6]*conj(B[22])
        +A[7]*conj(B[23]);

  C[3] = A[0]*conj(B[24])
        +A[1]*conj(B[25])
        +A[2]*conj(B[26])
        +A[3]*conj(B[27])
        +A[4]*conj(B[28])
        +A[5]*conj(B[29])
        +A[6]*conj(B[30])
        +A[7]*conj(B[31]);

  C[4] = A[0]*conj(B[32])
        +A[1]*conj(B[33])
        +A[2]*conj(B[34])
        +A[3]*conj(B[35])
        +A[4]*conj(B[36])
        +A[5]*conj(B[37])
        +A[6]*conj(B[38])
        +A[7]*conj(B[39]);

  C[5] = A[0]*conj(B[40])
        +A[1]*conj(B[41])
        +A[2]*conj(B[42])
        +A[3]*conj(B[43])
        +A[4]*conj(B[44])
        +A[5]*conj(B[45])
        +A[6]*conj(B[46])
        +A[7]*conj(B[47]);

  C[6] = A[0]*conj(B[48])
        +A[1]*conj(B[49])
        +A[2]*conj(B[50])
        +A[3]*conj(B[51])
        +A[4]*conj(B[52])
        +A[5]*conj(B[53])
        +A[6]*conj(B[54])
        +A[7]*conj(B[55]);

  C[7] = A[0]*conj(B[56])
        +A[1]*conj(B[57])
        +A[2]*conj(B[58])
        +A[3]*conj(B[59])
        +A[4]*conj(B[60])
        +A[5]*conj(B[61])
        +A[6]*conj(B[62])
        +A[7]*conj(B[63]);

  C[8] = A[8]*conj(B[0])
        +A[9]*conj(B[1])
        +A[10]*conj(B[2])
        +A[11]*conj(B[3])
        +A[12]*conj(B[4])
        +A[13]*conj(B[5])
        +A[14]*conj(B[6])
        +A[15]*conj(B[7]);

  C[9] = A[8]*conj(B[8])
        +A[9]*conj(B[9])
        +A[10]*conj(B[10])
        +A[11]*conj(B[11])
        +A[12]*conj(B[12])
        +A[13]*conj(B[13])
        +A[14]*conj(B[14])
        +A[15]*conj(B[15]);

  C[10] = A[8]*conj(B[16])
        +A[9]*conj(B[17])
        +A[10]*conj(B[18])
        +A[11]*conj(B[19])
        +A[12]*conj(B[20])
        +A[13]*conj(B[21])
        +A[14]*conj(B[22])
        +A[15]*conj(B[23]);

  C[11] = A[8]*conj(B[24])
        +A[9]*conj(B[25])
        +A[10]*conj(B[26])
        +A[11]*conj(B[27])
        +A[12]*conj(B[28])
        +A[13]*conj(B[29])
        +A[14]*conj(B[30])
        +A[15]*conj(B[31]);

  C[12] = A[8]*conj(B[32])
        +A[9]*conj(B[33])
        +A[10]*conj(B[34])
        +A[11]*conj(B[35])
        +A[12]*conj(B[36])
        +A[13]*conj(B[37])
        +A[14]*conj(B[38])
        +A[15]*conj(B[39]);

  C[13] = A[8]*conj(B[40])
        +A[9]*conj(B[41])
        +A[10]*conj(B[42])
        +A[11]*conj(B[43])
        +A[12]*conj(B[44])
        +A[13]*conj(B[45])
        +A[14]*conj(B[46])
        +A[15]*conj(B[47]);

  C[14] = A[8]*conj(B[48])
        +A[9]*conj(B[49])
        +A[10]*conj(B[50])
        +A[11]*conj(B[51])
        +A[12]*conj(B[52])
        +A[13]*conj(B[53])
        +A[14]*conj(B[54])
        +A[15]*conj(B[55]);

  C[15] = A[8]*conj(B[56])
        +A[9]*conj(B[57])
        +A[10]*conj(B[58])
        +A[11]*conj(B[59])
        +A[12]*conj(B[60])
        +A[13]*conj(B[61])
        +A[14]*conj(B[62])
        +A[15]*conj(B[63]);

  C[16] = A[16]*conj(B[0])
        +A[17]*conj(B[1])
        +A[18]*conj(B[2])
        +A[19]*conj(B[3])
        +A[20]*conj(B[4])
        +A[21]*conj(B[5])
        +A[22]*conj(B[6])
        +A[23]*conj(B[7]);

  C[17] = A[16]*conj(B[8])
        +A[17]*conj(B[9])
        +A[18]*conj(B[10])
        +A[19]*conj(B[11])
        +A[20]*conj(B[12])
        +A[21]*conj(B[13])
        +A[22]*conj(B[14])
        +A[23]*conj(B[15]);

  C[18] = A[16]*conj(B[16])
        +A[17]*conj(B[17])
        +A[18]*conj(B[18])
        +A[19]*conj(B[19])
        +A[20]*conj(B[20])
        +A[21]*conj(B[21])
        +A[22]*conj(B[22])
        +A[23]*conj(B[23]);

  C[19] = A[16]*conj(B[24])
        +A[17]*conj(B[25])
        +A[18]*conj(B[26])
        +A[19]*conj(B[27])
        +A[20]*conj(B[28])
        +A[21]*conj(B[29])
        +A[22]*conj(B[30])
        +A[23]*conj(B[31]);

  C[20] = A[16]*conj(B[32])
        +A[17]*conj(B[33])
        +A[18]*conj(B[34])
        +A[19]*conj(B[35])
        +A[20]*conj(B[36])
        +A[21]*conj(B[37])
        +A[22]*conj(B[38])
        +A[23]*conj(B[39]);

  C[21] = A[16]*conj(B[40])
        +A[17]*conj(B[41])
        +A[18]*conj(B[42])
        +A[19]*conj(B[43])
        +A[20]*conj(B[44])
        +A[21]*conj(B[45])
        +A[22]*conj(B[46])
        +A[23]*conj(B[47]);

  C[22] = A[16]*conj(B[48])
        +A[17]*conj(B[49])
        +A[18]*conj(B[50])
        +A[19]*conj(B[51])
        +A[20]*conj(B[52])
        +A[21]*conj(B[53])
        +A[22]*conj(B[54])
        +A[23]*conj(B[55]);

  C[23] = A[16]*conj(B[56])
        +A[17]*conj(B[57])
        +A[18]*conj(B[58])
        +A[19]*conj(B[59])
        +A[20]*conj(B[60])
        +A[21]*conj(B[61])
        +A[22]*conj(B[62])
        +A[23]*conj(B[63]);

  C[24] = A[24]*conj(B[0])
        +A[25]*conj(B[1])
        +A[26]*conj(B[2])
        +A[27]*conj(B[3])
        +A[28]*conj(B[4])
        +A[29]*conj(B[5])
        +A[30]*conj(B[6])
        +A[31]*conj(B[7]);

  C[25] = A[24]*conj(B[8])
        +A[25]*conj(B[9])
        +A[26]*conj(B[10])
        +A[27]*conj(B[11])
        +A[28]*conj(B[12])
        +A[29]*conj(B[13])
        +A[30]*conj(B[14])
        +A[31]*conj(B[15]);

  C[26] = A[24]*conj(B[16])
        +A[25]*conj(B[17])
        +A[26]*conj(B[18])
        +A[27]*conj(B[19])
        +A[28]*conj(B[20])
        +A[29]*conj(B[21])
        +A[30]*conj(B[22])
        +A[31]*conj(B[23]);

  C[27] = A[24]*conj(B[24])
        +A[25]*conj(B[25])
        +A[26]*conj(B[26])
        +A[27]*conj(B[27])
        +A[28]*conj(B[28])
        +A[29]*conj(B[29])
        +A[30]*conj(B[30])
        +A[31]*conj(B[31]);

  C[28] = A[24]*conj(B[32])
        +A[25]*conj(B[33])
        +A[26]*conj(B[34])
        +A[27]*conj(B[35])
        +A[28]*conj(B[36])
        +A[29]*conj(B[37])
        +A[30]*conj(B[38])
        +A[31]*conj(B[39]);

  C[29] = A[24]*conj(B[40])
        +A[25]*conj(B[41])
        +A[26]*conj(B[42])
        +A[27]*conj(B[43])
        +A[28]*conj(B[44])
        +A[29]*conj(B[45])
        +A[30]*conj(B[46])
        +A[31]*conj(B[47]);

  C[30] = A[24]*conj(B[48])
        +A[25]*conj(B[49])
        +A[26]*conj(B[50])
        +A[27]*conj(B[51])
        +A[28]*conj(B[52])
        +A[29]*conj(B[53])
        +A[30]*conj(B[54])
        +A[31]*conj(B[55]);

  C[31] = A[24]*conj(B[56])
        +A[25]*conj(B[57])
        +A[26]*conj(B[58])
        +A[27]*conj(B[59])
        +A[28]*conj(B[60])
        +A[29]*conj(B[61])
        +A[30]*conj(B[62])
        +A[31]*conj(B[63]);

  C[32] = A[32]*conj(B[0])
        +A[33]*conj(B[1])
        +A[34]*conj(B[2])
        +A[35]*conj(B[3])
        +A[36]*conj(B[4])
        +A[37]*conj(B[5])
        +A[38]*conj(B[6])
        +A[39]*conj(B[7]);

  C[33] = A[32]*conj(B[8])
        +A[33]*conj(B[9])
        +A[34]*conj(B[10])
        +A[35]*conj(B[11])
        +A[36]*conj(B[12])
        +A[37]*conj(B[13])
        +A[38]*conj(B[14])
        +A[39]*conj(B[15]);

  C[34] = A[32]*conj(B[16])
        +A[33]*conj(B[17])
        +A[34]*conj(B[18])
        +A[35]*conj(B[19])
        +A[36]*conj(B[20])
        +A[37]*conj(B[21])
        +A[38]*conj(B[22])
        +A[39]*conj(B[23]);

  C[35] = A[32]*conj(B[24])
        +A[33]*conj(B[25])
        +A[34]*conj(B[26])
        +A[35]*conj(B[27])
        +A[36]*conj(B[28])
        +A[37]*conj(B[29])
        +A[38]*conj(B[30])
        +A[39]*conj(B[31]);

  C[36] = A[32]*conj(B[32])
        +A[33]*conj(B[33])
        +A[34]*conj(B[34])
        +A[35]*conj(B[35])
        +A[36]*conj(B[36])
        +A[37]*conj(B[37])
        +A[38]*conj(B[38])
        +A[39]*conj(B[39]);

  C[37] = A[32]*conj(B[40])
        +A[33]*conj(B[41])
        +A[34]*conj(B[42])
        +A[35]*conj(B[43])
        +A[36]*conj(B[44])
        +A[37]*conj(B[45])
        +A[38]*conj(B[46])
        +A[39]*conj(B[47]);

  C[38] = A[32]*conj(B[48])
        +A[33]*conj(B[49])
        +A[34]*conj(B[50])
        +A[35]*conj(B[51])
        +A[36]*conj(B[52])
        +A[37]*conj(B[53])
        +A[38]*conj(B[54])
        +A[39]*conj(B[55]);

  C[39] = A[32]*conj(B[56])
        +A[33]*conj(B[57])
        +A[34]*conj(B[58])
        +A[35]*conj(B[59])
        +A[36]*conj(B[60])
        +A[37]*conj(B[61])
        +A[38]*conj(B[62])
        +A[39]*conj(B[63]);

  C[40] = A[40]*conj(B[0])
        +A[41]*conj(B[1])
        +A[42]*conj(B[2])
        +A[43]*conj(B[3])
        +A[44]*conj(B[4])
        +A[45]*conj(B[5])
        +A[46]*conj(B[6])
        +A[47]*conj(B[7]);

  C[41] = A[40]*conj(B[8])
        +A[41]*conj(B[9])
        +A[42]*conj(B[10])
        +A[43]*conj(B[11])
        +A[44]*conj(B[12])
        +A[45]*conj(B[13])
        +A[46]*conj(B[14])
        +A[47]*conj(B[15]);

  C[42] = A[40]*conj(B[16])
        +A[41]*conj(B[17])
        +A[42]*conj(B[18])
        +A[43]*conj(B[19])
        +A[44]*conj(B[20])
        +A[45]*conj(B[21])
        +A[46]*conj(B[22])
        +A[47]*conj(B[23]);

  C[43] = A[40]*conj(B[24])
        +A[41]*conj(B[25])
        +A[42]*conj(B[26])
        +A[43]*conj(B[27])
        +A[44]*conj(B[28])
        +A[45]*conj(B[29])
        +A[46]*conj(B[30])
        +A[47]*conj(B[31]);

  C[44] = A[40]*conj(B[32])
        +A[41]*conj(B[33])
        +A[42]*conj(B[34])
        +A[43]*conj(B[35])
        +A[44]*conj(B[36])
        +A[45]*conj(B[37])
        +A[46]*conj(B[38])
        +A[47]*conj(B[39]);

  C[45] = A[40]*conj(B[40])
        +A[41]*conj(B[41])
        +A[42]*conj(B[42])
        +A[43]*conj(B[43])
        +A[44]*conj(B[44])
        +A[45]*conj(B[45])
        +A[46]*conj(B[46])
        +A[47]*conj(B[47]);

  C[46] = A[40]*conj(B[48])
        +A[41]*conj(B[49])
        +A[42]*conj(B[50])
        +A[43]*conj(B[51])
        +A[44]*conj(B[52])
        +A[45]*conj(B[53])
        +A[46]*conj(B[54])
        +A[47]*conj(B[55]);

  C[47] = A[40]*conj(B[56])
        +A[41]*conj(B[57])
        +A[42]*conj(B[58])
        +A[43]*conj(B[59])
        +A[44]*conj(B[60])
        +A[45]*conj(B[61])
        +A[46]*conj(B[62])
        +A[47]*conj(B[63]);

  C[48] = A[48]*conj(B[0])
        +A[49]*conj(B[1])
        +A[50]*conj(B[2])
        +A[51]*conj(B[3])
        +A[52]*conj(B[4])
        +A[53]*conj(B[5])
        +A[54]*conj(B[6])
        +A[55]*conj(B[7]);

  C[49] = A[48]*conj(B[8])
        +A[49]*conj(B[9])
        +A[50]*conj(B[10])
        +A[51]*conj(B[11])
        +A[52]*conj(B[12])
        +A[53]*conj(B[13])
        +A[54]*conj(B[14])
        +A[55]*conj(B[15]);

  C[50] = A[48]*conj(B[16])
        +A[49]*conj(B[17])
        +A[50]*conj(B[18])
        +A[51]*conj(B[19])
        +A[52]*conj(B[20])
        +A[53]*conj(B[21])
        +A[54]*conj(B[22])
        +A[55]*conj(B[23]);

  C[51] = A[48]*conj(B[24])
        +A[49]*conj(B[25])
        +A[50]*conj(B[26])
        +A[51]*conj(B[27])
        +A[52]*conj(B[28])
        +A[53]*conj(B[29])
        +A[54]*conj(B[30])
        +A[55]*conj(B[31]);

  C[52] = A[48]*conj(B[32])
        +A[49]*conj(B[33])
        +A[50]*conj(B[34])
        +A[51]*conj(B[35])
        +A[52]*conj(B[36])
        +A[53]*conj(B[37])
        +A[54]*conj(B[38])
        +A[55]*conj(B[39]);

  C[53] = A[48]*conj(B[40])
        +A[49]*conj(B[41])
        +A[50]*conj(B[42])
        +A[51]*conj(B[43])
        +A[52]*conj(B[44])
        +A[53]*conj(B[45])
        +A[54]*conj(B[46])
        +A[55]*conj(B[47]);

  C[54] = A[48]*conj(B[48])
        +A[49]*conj(B[49])
        +A[50]*conj(B[50])
        +A[51]*conj(B[51])
        +A[52]*conj(B[52])
        +A[53]*conj(B[53])
        +A[54]*conj(B[54])
        +A[55]*conj(B[55]);

  C[55] = A[48]*conj(B[56])
        +A[49]*conj(B[57])
        +A[50]*conj(B[58])
        +A[51]*conj(B[59])
        +A[52]*conj(B[60])
        +A[53]*conj(B[61])
        +A[54]*conj(B[62])
        +A[55]*conj(B[63]);

  C[56] = A[56]*conj(B[0])
        +A[57]*conj(B[1])
        +A[58]*conj(B[2])
        +A[59]*conj(B[3])
        +A[60]*conj(B[4])
        +A[61]*conj(B[5])
        +A[62]*conj(B[6])
        +A[63]*conj(B[7]);

  C[57] = A[56]*conj(B[8])
        +A[57]*conj(B[9])
        +A[58]*conj(B[10])
        +A[59]*conj(B[11])
        +A[60]*conj(B[12])
        +A[61]*conj(B[13])
        +A[62]*conj(B[14])
        +A[63]*conj(B[15]);

  C[58] = A[56]*conj(B[16])
        +A[57]*conj(B[17])
        +A[58]*conj(B[18])
        +A[59]*conj(B[19])
        +A[60]*conj(B[20])
        +A[61]*conj(B[21])
        +A[62]*conj(B[22])
        +A[63]*conj(B[23]);

  C[59] = A[56]*conj(B[24])
        +A[57]*conj(B[25])
        +A[58]*conj(B[26])
        +A[59]*conj(B[27])
        +A[60]*conj(B[28])
        +A[61]*conj(B[29])
        +A[62]*conj(B[30])
        +A[63]*conj(B[31]);

  C[60] = A[56]*conj(B[32])
        +A[57]*conj(B[33])
        +A[58]*conj(B[34])
        +A[59]*conj(B[35])
        +A[60]*conj(B[36])
        +A[61]*conj(B[37])
        +A[62]*conj(B[38])
        +A[63]*conj(B[39]);

  C[61] = A[56]*conj(B[40])
        +A[57]*conj(B[41])
        +A[58]*conj(B[42])
        +A[59]*conj(B[43])
        +A[60]*conj(B[44])
        +A[61]*conj(B[45])
        +A[62]*conj(B[46])
        +A[63]*conj(B[47]);

  C[62] = A[56]*conj(B[48])
        +A[57]*conj(B[49])
        +A[58]*conj(B[50])
        +A[59]*conj(B[51])
        +A[60]*conj(B[52])
        +A[61]*conj(B[53])
        +A[62]*conj(B[54])
        +A[63]*conj(B[55]);

  C[63] = A[56]*conj(B[56])
        +A[57]*conj(B[57])
        +A[58]*conj(B[58])
        +A[59]*conj(B[59])
        +A[60]*conj(B[60])
        +A[61]*conj(B[61])
        +A[62]*conj(B[62])
        +A[63]*conj(B[63]);

};


inline void mult_C_equals_AdaggerB_for_SU8(dc *C, dc *A, dc *B) {

  C[0] = conj(A[0])*B[0]
        +conj(A[8])*B[8]
        +conj(A[16])*B[16]
        +conj(A[24])*B[24]
        +conj(A[32])*B[32]
        +conj(A[40])*B[40]
        +conj(A[48])*B[48]
        +conj(A[56])*B[56];

  C[1] = conj(A[0])*B[1]
        +conj(A[8])*B[9]
        +conj(A[16])*B[17]
        +conj(A[24])*B[25]
        +conj(A[32])*B[33]
        +conj(A[40])*B[41]
        +conj(A[48])*B[49]
        +conj(A[56])*B[57];

  C[2] = conj(A[0])*B[2]
        +conj(A[8])*B[10]
        +conj(A[16])*B[18]
        +conj(A[24])*B[26]
        +conj(A[32])*B[34]
        +conj(A[40])*B[42]
        +conj(A[48])*B[50]
        +conj(A[56])*B[58];

  C[3] = conj(A[0])*B[3]
        +conj(A[8])*B[11]
        +conj(A[16])*B[19]
        +conj(A[24])*B[27]
        +conj(A[32])*B[35]
        +conj(A[40])*B[43]
        +conj(A[48])*B[51]
        +conj(A[56])*B[59];

  C[4] = conj(A[0])*B[4]
        +conj(A[8])*B[12]
        +conj(A[16])*B[20]
        +conj(A[24])*B[28]
        +conj(A[32])*B[36]
        +conj(A[40])*B[44]
        +conj(A[48])*B[52]
        +conj(A[56])*B[60];

  C[5] = conj(A[0])*B[5]
        +conj(A[8])*B[13]
        +conj(A[16])*B[21]
        +conj(A[24])*B[29]
        +conj(A[32])*B[37]
        +conj(A[40])*B[45]
        +conj(A[48])*B[53]
        +conj(A[56])*B[61];

  C[6] = conj(A[0])*B[6]
        +conj(A[8])*B[14]
        +conj(A[16])*B[22]
        +conj(A[24])*B[30]
        +conj(A[32])*B[38]
        +conj(A[40])*B[46]
        +conj(A[48])*B[54]
        +conj(A[56])*B[62];

  C[7] = conj(A[0])*B[7]
        +conj(A[8])*B[15]
        +conj(A[16])*B[23]
        +conj(A[24])*B[31]
        +conj(A[32])*B[39]
        +conj(A[40])*B[47]
        +conj(A[48])*B[55]
        +conj(A[56])*B[63];

  C[8] = conj(A[1])*B[0]
        +conj(A[9])*B[8]
        +conj(A[17])*B[16]
        +conj(A[25])*B[24]
        +conj(A[33])*B[32]
        +conj(A[41])*B[40]
        +conj(A[49])*B[48]
        +conj(A[57])*B[56];

  C[9] = conj(A[1])*B[1]
        +conj(A[9])*B[9]
        +conj(A[17])*B[17]
        +conj(A[25])*B[25]
        +conj(A[33])*B[33]
        +conj(A[41])*B[41]
        +conj(A[49])*B[49]
        +conj(A[57])*B[57];

  C[10] = conj(A[1])*B[2]
        +conj(A[9])*B[10]
        +conj(A[17])*B[18]
        +conj(A[25])*B[26]
        +conj(A[33])*B[34]
        +conj(A[41])*B[42]
        +conj(A[49])*B[50]
        +conj(A[57])*B[58];

  C[11] = conj(A[1])*B[3]
        +conj(A[9])*B[11]
        +conj(A[17])*B[19]
        +conj(A[25])*B[27]
        +conj(A[33])*B[35]
        +conj(A[41])*B[43]
        +conj(A[49])*B[51]
        +conj(A[57])*B[59];

  C[12] = conj(A[1])*B[4]
        +conj(A[9])*B[12]
        +conj(A[17])*B[20]
        +conj(A[25])*B[28]
        +conj(A[33])*B[36]
        +conj(A[41])*B[44]
        +conj(A[49])*B[52]
        +conj(A[57])*B[60];

  C[13] = conj(A[1])*B[5]
        +conj(A[9])*B[13]
        +conj(A[17])*B[21]
        +conj(A[25])*B[29]
        +conj(A[33])*B[37]
        +conj(A[41])*B[45]
        +conj(A[49])*B[53]
        +conj(A[57])*B[61];

  C[14] = conj(A[1])*B[6]
        +conj(A[9])*B[14]
        +conj(A[17])*B[22]
        +conj(A[25])*B[30]
        +conj(A[33])*B[38]
        +conj(A[41])*B[46]
        +conj(A[49])*B[54]
        +conj(A[57])*B[62];

  C[15] = conj(A[1])*B[7]
        +conj(A[9])*B[15]
        +conj(A[17])*B[23]
        +conj(A[25])*B[31]
        +conj(A[33])*B[39]
        +conj(A[41])*B[47]
        +conj(A[49])*B[55]
        +conj(A[57])*B[63];

  C[16] = conj(A[2])*B[0]
        +conj(A[10])*B[8]
        +conj(A[18])*B[16]
        +conj(A[26])*B[24]
        +conj(A[34])*B[32]
        +conj(A[42])*B[40]
        +conj(A[50])*B[48]
        +conj(A[58])*B[56];

  C[17] = conj(A[2])*B[1]
        +conj(A[10])*B[9]
        +conj(A[18])*B[17]
        +conj(A[26])*B[25]
        +conj(A[34])*B[33]
        +conj(A[42])*B[41]
        +conj(A[50])*B[49]
        +conj(A[58])*B[57];

  C[18] = conj(A[2])*B[2]
        +conj(A[10])*B[10]
        +conj(A[18])*B[18]
        +conj(A[26])*B[26]
        +conj(A[34])*B[34]
        +conj(A[42])*B[42]
        +conj(A[50])*B[50]
        +conj(A[58])*B[58];

  C[19] = conj(A[2])*B[3]
        +conj(A[10])*B[11]
        +conj(A[18])*B[19]
        +conj(A[26])*B[27]
        +conj(A[34])*B[35]
        +conj(A[42])*B[43]
        +conj(A[50])*B[51]
        +conj(A[58])*B[59];

  C[20] = conj(A[2])*B[4]
        +conj(A[10])*B[12]
        +conj(A[18])*B[20]
        +conj(A[26])*B[28]
        +conj(A[34])*B[36]
        +conj(A[42])*B[44]
        +conj(A[50])*B[52]
        +conj(A[58])*B[60];

  C[21] = conj(A[2])*B[5]
        +conj(A[10])*B[13]
        +conj(A[18])*B[21]
        +conj(A[26])*B[29]
        +conj(A[34])*B[37]
        +conj(A[42])*B[45]
        +conj(A[50])*B[53]
        +conj(A[58])*B[61];

  C[22] = conj(A[2])*B[6]
        +conj(A[10])*B[14]
        +conj(A[18])*B[22]
        +conj(A[26])*B[30]
        +conj(A[34])*B[38]
        +conj(A[42])*B[46]
        +conj(A[50])*B[54]
        +conj(A[58])*B[62];

  C[23] = conj(A[2])*B[7]
        +conj(A[10])*B[15]
        +conj(A[18])*B[23]
        +conj(A[26])*B[31]
        +conj(A[34])*B[39]
        +conj(A[42])*B[47]
        +conj(A[50])*B[55]
        +conj(A[58])*B[63];

  C[24] = conj(A[3])*B[0]
        +conj(A[11])*B[8]
        +conj(A[19])*B[16]
        +conj(A[27])*B[24]
        +conj(A[35])*B[32]
        +conj(A[43])*B[40]
        +conj(A[51])*B[48]
        +conj(A[59])*B[56];

  C[25] = conj(A[3])*B[1]
        +conj(A[11])*B[9]
        +conj(A[19])*B[17]
        +conj(A[27])*B[25]
        +conj(A[35])*B[33]
        +conj(A[43])*B[41]
        +conj(A[51])*B[49]
        +conj(A[59])*B[57];

  C[26] = conj(A[3])*B[2]
        +conj(A[11])*B[10]
        +conj(A[19])*B[18]
        +conj(A[27])*B[26]
        +conj(A[35])*B[34]
        +conj(A[43])*B[42]
        +conj(A[51])*B[50]
        +conj(A[59])*B[58];

  C[27] = conj(A[3])*B[3]
        +conj(A[11])*B[11]
        +conj(A[19])*B[19]
        +conj(A[27])*B[27]
        +conj(A[35])*B[35]
        +conj(A[43])*B[43]
        +conj(A[51])*B[51]
        +conj(A[59])*B[59];

  C[28] = conj(A[3])*B[4]
        +conj(A[11])*B[12]
        +conj(A[19])*B[20]
        +conj(A[27])*B[28]
        +conj(A[35])*B[36]
        +conj(A[43])*B[44]
        +conj(A[51])*B[52]
        +conj(A[59])*B[60];

  C[29] = conj(A[3])*B[5]
        +conj(A[11])*B[13]
        +conj(A[19])*B[21]
        +conj(A[27])*B[29]
        +conj(A[35])*B[37]
        +conj(A[43])*B[45]
        +conj(A[51])*B[53]
        +conj(A[59])*B[61];

  C[30] = conj(A[3])*B[6]
        +conj(A[11])*B[14]
        +conj(A[19])*B[22]
        +conj(A[27])*B[30]
        +conj(A[35])*B[38]
        +conj(A[43])*B[46]
        +conj(A[51])*B[54]
        +conj(A[59])*B[62];

  C[31] = conj(A[3])*B[7]
        +conj(A[11])*B[15]
        +conj(A[19])*B[23]
        +conj(A[27])*B[31]
        +conj(A[35])*B[39]
        +conj(A[43])*B[47]
        +conj(A[51])*B[55]
        +conj(A[59])*B[63];

  C[32] = conj(A[4])*B[0]
        +conj(A[12])*B[8]
        +conj(A[20])*B[16]
        +conj(A[28])*B[24]
        +conj(A[36])*B[32]
        +conj(A[44])*B[40]
        +conj(A[52])*B[48]
        +conj(A[60])*B[56];

  C[33] = conj(A[4])*B[1]
        +conj(A[12])*B[9]
        +conj(A[20])*B[17]
        +conj(A[28])*B[25]
        +conj(A[36])*B[33]
        +conj(A[44])*B[41]
        +conj(A[52])*B[49]
        +conj(A[60])*B[57];

  C[34] = conj(A[4])*B[2]
        +conj(A[12])*B[10]
        +conj(A[20])*B[18]
        +conj(A[28])*B[26]
        +conj(A[36])*B[34]
        +conj(A[44])*B[42]
        +conj(A[52])*B[50]
        +conj(A[60])*B[58];

  C[35] = conj(A[4])*B[3]
        +conj(A[12])*B[11]
        +conj(A[20])*B[19]
        +conj(A[28])*B[27]
        +conj(A[36])*B[35]
        +conj(A[44])*B[43]
        +conj(A[52])*B[51]
        +conj(A[60])*B[59];

  C[36] = conj(A[4])*B[4]
        +conj(A[12])*B[12]
        +conj(A[20])*B[20]
        +conj(A[28])*B[28]
        +conj(A[36])*B[36]
        +conj(A[44])*B[44]
        +conj(A[52])*B[52]
        +conj(A[60])*B[60];

  C[37] = conj(A[4])*B[5]
        +conj(A[12])*B[13]
        +conj(A[20])*B[21]
        +conj(A[28])*B[29]
        +conj(A[36])*B[37]
        +conj(A[44])*B[45]
        +conj(A[52])*B[53]
        +conj(A[60])*B[61];

  C[38] = conj(A[4])*B[6]
        +conj(A[12])*B[14]
        +conj(A[20])*B[22]
        +conj(A[28])*B[30]
        +conj(A[36])*B[38]
        +conj(A[44])*B[46]
        +conj(A[52])*B[54]
        +conj(A[60])*B[62];

  C[39] = conj(A[4])*B[7]
        +conj(A[12])*B[15]
        +conj(A[20])*B[23]
        +conj(A[28])*B[31]
        +conj(A[36])*B[39]
        +conj(A[44])*B[47]
        +conj(A[52])*B[55]
        +conj(A[60])*B[63];

  C[40] = conj(A[5])*B[0]
        +conj(A[13])*B[8]
        +conj(A[21])*B[16]
        +conj(A[29])*B[24]
        +conj(A[37])*B[32]
        +conj(A[45])*B[40]
        +conj(A[53])*B[48]
        +conj(A[61])*B[56];

  C[41] = conj(A[5])*B[1]
        +conj(A[13])*B[9]
        +conj(A[21])*B[17]
        +conj(A[29])*B[25]
        +conj(A[37])*B[33]
        +conj(A[45])*B[41]
        +conj(A[53])*B[49]
        +conj(A[61])*B[57];

  C[42] = conj(A[5])*B[2]
        +conj(A[13])*B[10]
        +conj(A[21])*B[18]
        +conj(A[29])*B[26]
        +conj(A[37])*B[34]
        +conj(A[45])*B[42]
        +conj(A[53])*B[50]
        +conj(A[61])*B[58];

  C[43] = conj(A[5])*B[3]
        +conj(A[13])*B[11]
        +conj(A[21])*B[19]
        +conj(A[29])*B[27]
        +conj(A[37])*B[35]
        +conj(A[45])*B[43]
        +conj(A[53])*B[51]
        +conj(A[61])*B[59];

  C[44] = conj(A[5])*B[4]
        +conj(A[13])*B[12]
        +conj(A[21])*B[20]
        +conj(A[29])*B[28]
        +conj(A[37])*B[36]
        +conj(A[45])*B[44]
        +conj(A[53])*B[52]
        +conj(A[61])*B[60];

  C[45] = conj(A[5])*B[5]
        +conj(A[13])*B[13]
        +conj(A[21])*B[21]
        +conj(A[29])*B[29]
        +conj(A[37])*B[37]
        +conj(A[45])*B[45]
        +conj(A[53])*B[53]
        +conj(A[61])*B[61];

  C[46] = conj(A[5])*B[6]
        +conj(A[13])*B[14]
        +conj(A[21])*B[22]
        +conj(A[29])*B[30]
        +conj(A[37])*B[38]
        +conj(A[45])*B[46]
        +conj(A[53])*B[54]
        +conj(A[61])*B[62];

  C[47] = conj(A[5])*B[7]
        +conj(A[13])*B[15]
        +conj(A[21])*B[23]
        +conj(A[29])*B[31]
        +conj(A[37])*B[39]
        +conj(A[45])*B[47]
        +conj(A[53])*B[55]
        +conj(A[61])*B[63];

  C[48] = conj(A[6])*B[0]
        +conj(A[14])*B[8]
        +conj(A[22])*B[16]
        +conj(A[30])*B[24]
        +conj(A[38])*B[32]
        +conj(A[46])*B[40]
        +conj(A[54])*B[48]
        +conj(A[62])*B[56];

  C[49] = conj(A[6])*B[1]
        +conj(A[14])*B[9]
        +conj(A[22])*B[17]
        +conj(A[30])*B[25]
        +conj(A[38])*B[33]
        +conj(A[46])*B[41]
        +conj(A[54])*B[49]
        +conj(A[62])*B[57];

  C[50] = conj(A[6])*B[2]
        +conj(A[14])*B[10]
        +conj(A[22])*B[18]
        +conj(A[30])*B[26]
        +conj(A[38])*B[34]
        +conj(A[46])*B[42]
        +conj(A[54])*B[50]
        +conj(A[62])*B[58];

  C[51] = conj(A[6])*B[3]
        +conj(A[14])*B[11]
        +conj(A[22])*B[19]
        +conj(A[30])*B[27]
        +conj(A[38])*B[35]
        +conj(A[46])*B[43]
        +conj(A[54])*B[51]
        +conj(A[62])*B[59];

  C[52] = conj(A[6])*B[4]
        +conj(A[14])*B[12]
        +conj(A[22])*B[20]
        +conj(A[30])*B[28]
        +conj(A[38])*B[36]
        +conj(A[46])*B[44]
        +conj(A[54])*B[52]
        +conj(A[62])*B[60];

  C[53] = conj(A[6])*B[5]
        +conj(A[14])*B[13]
        +conj(A[22])*B[21]
        +conj(A[30])*B[29]
        +conj(A[38])*B[37]
        +conj(A[46])*B[45]
        +conj(A[54])*B[53]
        +conj(A[62])*B[61];

  C[54] = conj(A[6])*B[6]
        +conj(A[14])*B[14]
        +conj(A[22])*B[22]
        +conj(A[30])*B[30]
        +conj(A[38])*B[38]
        +conj(A[46])*B[46]
        +conj(A[54])*B[54]
        +conj(A[62])*B[62];

  C[55] = conj(A[6])*B[7]
        +conj(A[14])*B[15]
        +conj(A[22])*B[23]
        +conj(A[30])*B[31]
        +conj(A[38])*B[39]
        +conj(A[46])*B[47]
        +conj(A[54])*B[55]
        +conj(A[62])*B[63];

  C[56] = conj(A[7])*B[0]
        +conj(A[15])*B[8]
        +conj(A[23])*B[16]
        +conj(A[31])*B[24]
        +conj(A[39])*B[32]
        +conj(A[47])*B[40]
        +conj(A[55])*B[48]
        +conj(A[63])*B[56];

  C[57] = conj(A[7])*B[1]
        +conj(A[15])*B[9]
        +conj(A[23])*B[17]
        +conj(A[31])*B[25]
        +conj(A[39])*B[33]
        +conj(A[47])*B[41]
        +conj(A[55])*B[49]
        +conj(A[63])*B[57];

  C[58] = conj(A[7])*B[2]
        +conj(A[15])*B[10]
        +conj(A[23])*B[18]
        +conj(A[31])*B[26]
        +conj(A[39])*B[34]
        +conj(A[47])*B[42]
        +conj(A[55])*B[50]
        +conj(A[63])*B[58];

  C[59] = conj(A[7])*B[3]
        +conj(A[15])*B[11]
        +conj(A[23])*B[19]
        +conj(A[31])*B[27]
        +conj(A[39])*B[35]
        +conj(A[47])*B[43]
        +conj(A[55])*B[51]
        +conj(A[63])*B[59];

  C[60] = conj(A[7])*B[4]
        +conj(A[15])*B[12]
        +conj(A[23])*B[20]
        +conj(A[31])*B[28]
        +conj(A[39])*B[36]
        +conj(A[47])*B[44]
        +conj(A[55])*B[52]
        +conj(A[63])*B[60];

  C[61] = conj(A[7])*B[5]
        +conj(A[15])*B[13]
        +conj(A[23])*B[21]
        +conj(A[31])*B[29]
        +conj(A[39])*B[37]
        +conj(A[47])*B[45]
        +conj(A[55])*B[53]
        +conj(A[63])*B[61];

  C[62] = conj(A[7])*B[6]
        +conj(A[15])*B[14]
        +conj(A[23])*B[22]
        +conj(A[31])*B[30]
        +conj(A[39])*B[38]
        +conj(A[47])*B[46]
        +conj(A[55])*B[54]
        +conj(A[63])*B[62];

  C[63] = conj(A[7])*B[7]
        +conj(A[15])*B[15]
        +conj(A[23])*B[23]
        +conj(A[31])*B[31]
        +conj(A[39])*B[39]
        +conj(A[47])*B[47]
        +conj(A[55])*B[55]
        +conj(A[63])*B[63];

};


