inline void mult_C_equals_AB_for_SU7(dc *C, dc *A, dc *B) {

  C[0] = A[0]*B[0]
        +A[1]*B[7]
        +A[2]*B[14]
        +A[3]*B[21]
        +A[4]*B[28]
        +A[5]*B[35]
        +A[6]*B[42];

  C[1] = A[0]*B[1]
        +A[1]*B[8]
        +A[2]*B[15]
        +A[3]*B[22]
        +A[4]*B[29]
        +A[5]*B[36]
        +A[6]*B[43];

  C[2] = A[0]*B[2]
        +A[1]*B[9]
        +A[2]*B[16]
        +A[3]*B[23]
        +A[4]*B[30]
        +A[5]*B[37]
        +A[6]*B[44];

  C[3] = A[0]*B[3]
        +A[1]*B[10]
        +A[2]*B[17]
        +A[3]*B[24]
        +A[4]*B[31]
        +A[5]*B[38]
        +A[6]*B[45];

  C[4] = A[0]*B[4]
        +A[1]*B[11]
        +A[2]*B[18]
        +A[3]*B[25]
        +A[4]*B[32]
        +A[5]*B[39]
        +A[6]*B[46];

  C[5] = A[0]*B[5]
        +A[1]*B[12]
        +A[2]*B[19]
        +A[3]*B[26]
        +A[4]*B[33]
        +A[5]*B[40]
        +A[6]*B[47];

  C[6] = A[0]*B[6]
        +A[1]*B[13]
        +A[2]*B[20]
        +A[3]*B[27]
        +A[4]*B[34]
        +A[5]*B[41]
        +A[6]*B[48];

  C[7] = A[7]*B[0]
        +A[8]*B[7]
        +A[9]*B[14]
        +A[10]*B[21]
        +A[11]*B[28]
        +A[12]*B[35]
        +A[13]*B[42];

  C[8] = A[7]*B[1]
        +A[8]*B[8]
        +A[9]*B[15]
        +A[10]*B[22]
        +A[11]*B[29]
        +A[12]*B[36]
        +A[13]*B[43];

  C[9] = A[7]*B[2]
        +A[8]*B[9]
        +A[9]*B[16]
        +A[10]*B[23]
        +A[11]*B[30]
        +A[12]*B[37]
        +A[13]*B[44];

  C[10] = A[7]*B[3]
        +A[8]*B[10]
        +A[9]*B[17]
        +A[10]*B[24]
        +A[11]*B[31]
        +A[12]*B[38]
        +A[13]*B[45];

  C[11] = A[7]*B[4]
        +A[8]*B[11]
        +A[9]*B[18]
        +A[10]*B[25]
        +A[11]*B[32]
        +A[12]*B[39]
        +A[13]*B[46];

  C[12] = A[7]*B[5]
        +A[8]*B[12]
        +A[9]*B[19]
        +A[10]*B[26]
        +A[11]*B[33]
        +A[12]*B[40]
        +A[13]*B[47];

  C[13] = A[7]*B[6]
        +A[8]*B[13]
        +A[9]*B[20]
        +A[10]*B[27]
        +A[11]*B[34]
        +A[12]*B[41]
        +A[13]*B[48];

  C[14] = A[14]*B[0]
        +A[15]*B[7]
        +A[16]*B[14]
        +A[17]*B[21]
        +A[18]*B[28]
        +A[19]*B[35]
        +A[20]*B[42];

  C[15] = A[14]*B[1]
        +A[15]*B[8]
        +A[16]*B[15]
        +A[17]*B[22]
        +A[18]*B[29]
        +A[19]*B[36]
        +A[20]*B[43];

  C[16] = A[14]*B[2]
        +A[15]*B[9]
        +A[16]*B[16]
        +A[17]*B[23]
        +A[18]*B[30]
        +A[19]*B[37]
        +A[20]*B[44];

  C[17] = A[14]*B[3]
        +A[15]*B[10]
        +A[16]*B[17]
        +A[17]*B[24]
        +A[18]*B[31]
        +A[19]*B[38]
        +A[20]*B[45];

  C[18] = A[14]*B[4]
        +A[15]*B[11]
        +A[16]*B[18]
        +A[17]*B[25]
        +A[18]*B[32]
        +A[19]*B[39]
        +A[20]*B[46];

  C[19] = A[14]*B[5]
        +A[15]*B[12]
        +A[16]*B[19]
        +A[17]*B[26]
        +A[18]*B[33]
        +A[19]*B[40]
        +A[20]*B[47];

  C[20] = A[14]*B[6]
        +A[15]*B[13]
        +A[16]*B[20]
        +A[17]*B[27]
        +A[18]*B[34]
        +A[19]*B[41]
        +A[20]*B[48];

  C[21] = A[21]*B[0]
        +A[22]*B[7]
        +A[23]*B[14]
        +A[24]*B[21]
        +A[25]*B[28]
        +A[26]*B[35]
        +A[27]*B[42];

  C[22] = A[21]*B[1]
        +A[22]*B[8]
        +A[23]*B[15]
        +A[24]*B[22]
        +A[25]*B[29]
        +A[26]*B[36]
        +A[27]*B[43];

  C[23] = A[21]*B[2]
        +A[22]*B[9]
        +A[23]*B[16]
        +A[24]*B[23]
        +A[25]*B[30]
        +A[26]*B[37]
        +A[27]*B[44];

  C[24] = A[21]*B[3]
        +A[22]*B[10]
        +A[23]*B[17]
        +A[24]*B[24]
        +A[25]*B[31]
        +A[26]*B[38]
        +A[27]*B[45];

  C[25] = A[21]*B[4]
        +A[22]*B[11]
        +A[23]*B[18]
        +A[24]*B[25]
        +A[25]*B[32]
        +A[26]*B[39]
        +A[27]*B[46];

  C[26] = A[21]*B[5]
        +A[22]*B[12]
        +A[23]*B[19]
        +A[24]*B[26]
        +A[25]*B[33]
        +A[26]*B[40]
        +A[27]*B[47];

  C[27] = A[21]*B[6]
        +A[22]*B[13]
        +A[23]*B[20]
        +A[24]*B[27]
        +A[25]*B[34]
        +A[26]*B[41]
        +A[27]*B[48];

  C[28] = A[28]*B[0]
        +A[29]*B[7]
        +A[30]*B[14]
        +A[31]*B[21]
        +A[32]*B[28]
        +A[33]*B[35]
        +A[34]*B[42];

  C[29] = A[28]*B[1]
        +A[29]*B[8]
        +A[30]*B[15]
        +A[31]*B[22]
        +A[32]*B[29]
        +A[33]*B[36]
        +A[34]*B[43];

  C[30] = A[28]*B[2]
        +A[29]*B[9]
        +A[30]*B[16]
        +A[31]*B[23]
        +A[32]*B[30]
        +A[33]*B[37]
        +A[34]*B[44];

  C[31] = A[28]*B[3]
        +A[29]*B[10]
        +A[30]*B[17]
        +A[31]*B[24]
        +A[32]*B[31]
        +A[33]*B[38]
        +A[34]*B[45];

  C[32] = A[28]*B[4]
        +A[29]*B[11]
        +A[30]*B[18]
        +A[31]*B[25]
        +A[32]*B[32]
        +A[33]*B[39]
        +A[34]*B[46];

  C[33] = A[28]*B[5]
        +A[29]*B[12]
        +A[30]*B[19]
        +A[31]*B[26]
        +A[32]*B[33]
        +A[33]*B[40]
        +A[34]*B[47];

  C[34] = A[28]*B[6]
        +A[29]*B[13]
        +A[30]*B[20]
        +A[31]*B[27]
        +A[32]*B[34]
        +A[33]*B[41]
        +A[34]*B[48];

  C[35] = A[35]*B[0]
        +A[36]*B[7]
        +A[37]*B[14]
        +A[38]*B[21]
        +A[39]*B[28]
        +A[40]*B[35]
        +A[41]*B[42];

  C[36] = A[35]*B[1]
        +A[36]*B[8]
        +A[37]*B[15]
        +A[38]*B[22]
        +A[39]*B[29]
        +A[40]*B[36]
        +A[41]*B[43];

  C[37] = A[35]*B[2]
        +A[36]*B[9]
        +A[37]*B[16]
        +A[38]*B[23]
        +A[39]*B[30]
        +A[40]*B[37]
        +A[41]*B[44];

  C[38] = A[35]*B[3]
        +A[36]*B[10]
        +A[37]*B[17]
        +A[38]*B[24]
        +A[39]*B[31]
        +A[40]*B[38]
        +A[41]*B[45];

  C[39] = A[35]*B[4]
        +A[36]*B[11]
        +A[37]*B[18]
        +A[38]*B[25]
        +A[39]*B[32]
        +A[40]*B[39]
        +A[41]*B[46];

  C[40] = A[35]*B[5]
        +A[36]*B[12]
        +A[37]*B[19]
        +A[38]*B[26]
        +A[39]*B[33]
        +A[40]*B[40]
        +A[41]*B[47];

  C[41] = A[35]*B[6]
        +A[36]*B[13]
        +A[37]*B[20]
        +A[38]*B[27]
        +A[39]*B[34]
        +A[40]*B[41]
        +A[41]*B[48];

  C[42] = A[42]*B[0]
        +A[43]*B[7]
        +A[44]*B[14]
        +A[45]*B[21]
        +A[46]*B[28]
        +A[47]*B[35]
        +A[48]*B[42];

  C[43] = A[42]*B[1]
        +A[43]*B[8]
        +A[44]*B[15]
        +A[45]*B[22]
        +A[46]*B[29]
        +A[47]*B[36]
        +A[48]*B[43];

  C[44] = A[42]*B[2]
        +A[43]*B[9]
        +A[44]*B[16]
        +A[45]*B[23]
        +A[46]*B[30]
        +A[47]*B[37]
        +A[48]*B[44];

  C[45] = A[42]*B[3]
        +A[43]*B[10]
        +A[44]*B[17]
        +A[45]*B[24]
        +A[46]*B[31]
        +A[47]*B[38]
        +A[48]*B[45];

  C[46] = A[42]*B[4]
        +A[43]*B[11]
        +A[44]*B[18]
        +A[45]*B[25]
        +A[46]*B[32]
        +A[47]*B[39]
        +A[48]*B[46];

  C[47] = A[42]*B[5]
        +A[43]*B[12]
        +A[44]*B[19]
        +A[45]*B[26]
        +A[46]*B[33]
        +A[47]*B[40]
        +A[48]*B[47];

  C[48] = A[42]*B[6]
        +A[43]*B[13]
        +A[44]*B[20]
        +A[45]*B[27]
        +A[46]*B[34]
        +A[47]*B[41]
        +A[48]*B[48];

};


inline void mult_C_equals_ABdagger_for_SU7(dc *C, dc *A, dc *B) {

  C[0] = A[0]*conj(B[0])
        +A[1]*conj(B[1])
        +A[2]*conj(B[2])
        +A[3]*conj(B[3])
        +A[4]*conj(B[4])
        +A[5]*conj(B[5])
        +A[6]*conj(B[6]);

  C[1] = A[0]*conj(B[7])
        +A[1]*conj(B[8])
        +A[2]*conj(B[9])
        +A[3]*conj(B[10])
        +A[4]*conj(B[11])
        +A[5]*conj(B[12])
        +A[6]*conj(B[13]);

  C[2] = A[0]*conj(B[14])
        +A[1]*conj(B[15])
        +A[2]*conj(B[16])
        +A[3]*conj(B[17])
        +A[4]*conj(B[18])
        +A[5]*conj(B[19])
        +A[6]*conj(B[20]);

  C[3] = A[0]*conj(B[21])
        +A[1]*conj(B[22])
        +A[2]*conj(B[23])
        +A[3]*conj(B[24])
        +A[4]*conj(B[25])
        +A[5]*conj(B[26])
        +A[6]*conj(B[27]);

  C[4] = A[0]*conj(B[28])
        +A[1]*conj(B[29])
        +A[2]*conj(B[30])
        +A[3]*conj(B[31])
        +A[4]*conj(B[32])
        +A[5]*conj(B[33])
        +A[6]*conj(B[34]);

  C[5] = A[0]*conj(B[35])
        +A[1]*conj(B[36])
        +A[2]*conj(B[37])
        +A[3]*conj(B[38])
        +A[4]*conj(B[39])
        +A[5]*conj(B[40])
        +A[6]*conj(B[41]);

  C[6] = A[0]*conj(B[42])
        +A[1]*conj(B[43])
        +A[2]*conj(B[44])
        +A[3]*conj(B[45])
        +A[4]*conj(B[46])
        +A[5]*conj(B[47])
        +A[6]*conj(B[48]);

  C[7] = A[7]*conj(B[0])
        +A[8]*conj(B[1])
        +A[9]*conj(B[2])
        +A[10]*conj(B[3])
        +A[11]*conj(B[4])
        +A[12]*conj(B[5])
        +A[13]*conj(B[6]);

  C[8] = A[7]*conj(B[7])
        +A[8]*conj(B[8])
        +A[9]*conj(B[9])
        +A[10]*conj(B[10])
        +A[11]*conj(B[11])
        +A[12]*conj(B[12])
        +A[13]*conj(B[13]);

  C[9] = A[7]*conj(B[14])
        +A[8]*conj(B[15])
        +A[9]*conj(B[16])
        +A[10]*conj(B[17])
        +A[11]*conj(B[18])
        +A[12]*conj(B[19])
        +A[13]*conj(B[20]);

  C[10] = A[7]*conj(B[21])
        +A[8]*conj(B[22])
        +A[9]*conj(B[23])
        +A[10]*conj(B[24])
        +A[11]*conj(B[25])
        +A[12]*conj(B[26])
        +A[13]*conj(B[27]);

  C[11] = A[7]*conj(B[28])
        +A[8]*conj(B[29])
        +A[9]*conj(B[30])
        +A[10]*conj(B[31])
        +A[11]*conj(B[32])
        +A[12]*conj(B[33])
        +A[13]*conj(B[34]);

  C[12] = A[7]*conj(B[35])
        +A[8]*conj(B[36])
        +A[9]*conj(B[37])
        +A[10]*conj(B[38])
        +A[11]*conj(B[39])
        +A[12]*conj(B[40])
        +A[13]*conj(B[41]);

  C[13] = A[7]*conj(B[42])
        +A[8]*conj(B[43])
        +A[9]*conj(B[44])
        +A[10]*conj(B[45])
        +A[11]*conj(B[46])
        +A[12]*conj(B[47])
        +A[13]*conj(B[48]);

  C[14] = A[14]*conj(B[0])
        +A[15]*conj(B[1])
        +A[16]*conj(B[2])
        +A[17]*conj(B[3])
        +A[18]*conj(B[4])
        +A[19]*conj(B[5])
        +A[20]*conj(B[6]);

  C[15] = A[14]*conj(B[7])
        +A[15]*conj(B[8])
        +A[16]*conj(B[9])
        +A[17]*conj(B[10])
        +A[18]*conj(B[11])
        +A[19]*conj(B[12])
        +A[20]*conj(B[13]);

  C[16] = A[14]*conj(B[14])
        +A[15]*conj(B[15])
        +A[16]*conj(B[16])
        +A[17]*conj(B[17])
        +A[18]*conj(B[18])
        +A[19]*conj(B[19])
        +A[20]*conj(B[20]);

  C[17] = A[14]*conj(B[21])
        +A[15]*conj(B[22])
        +A[16]*conj(B[23])
        +A[17]*conj(B[24])
        +A[18]*conj(B[25])
        +A[19]*conj(B[26])
        +A[20]*conj(B[27]);

  C[18] = A[14]*conj(B[28])
        +A[15]*conj(B[29])
        +A[16]*conj(B[30])
        +A[17]*conj(B[31])
        +A[18]*conj(B[32])
        +A[19]*conj(B[33])
        +A[20]*conj(B[34]);

  C[19] = A[14]*conj(B[35])
        +A[15]*conj(B[36])
        +A[16]*conj(B[37])
        +A[17]*conj(B[38])
        +A[18]*conj(B[39])
        +A[19]*conj(B[40])
        +A[20]*conj(B[41]);

  C[20] = A[14]*conj(B[42])
        +A[15]*conj(B[43])
        +A[16]*conj(B[44])
        +A[17]*conj(B[45])
        +A[18]*conj(B[46])
        +A[19]*conj(B[47])
        +A[20]*conj(B[48]);

  C[21] = A[21]*conj(B[0])
        +A[22]*conj(B[1])
        +A[23]*conj(B[2])
        +A[24]*conj(B[3])
        +A[25]*conj(B[4])
        +A[26]*conj(B[5])
        +A[27]*conj(B[6]);

  C[22] = A[21]*conj(B[7])
        +A[22]*conj(B[8])
        +A[23]*conj(B[9])
        +A[24]*conj(B[10])
        +A[25]*conj(B[11])
        +A[26]*conj(B[12])
        +A[27]*conj(B[13]);

  C[23] = A[21]*conj(B[14])
        +A[22]*conj(B[15])
        +A[23]*conj(B[16])
        +A[24]*conj(B[17])
        +A[25]*conj(B[18])
        +A[26]*conj(B[19])
        +A[27]*conj(B[20]);

  C[24] = A[21]*conj(B[21])
        +A[22]*conj(B[22])
        +A[23]*conj(B[23])
        +A[24]*conj(B[24])
        +A[25]*conj(B[25])
        +A[26]*conj(B[26])
        +A[27]*conj(B[27]);

  C[25] = A[21]*conj(B[28])
        +A[22]*conj(B[29])
        +A[23]*conj(B[30])
        +A[24]*conj(B[31])
        +A[25]*conj(B[32])
        +A[26]*conj(B[33])
        +A[27]*conj(B[34]);

  C[26] = A[21]*conj(B[35])
        +A[22]*conj(B[36])
        +A[23]*conj(B[37])
        +A[24]*conj(B[38])
        +A[25]*conj(B[39])
        +A[26]*conj(B[40])
        +A[27]*conj(B[41]);

  C[27] = A[21]*conj(B[42])
        +A[22]*conj(B[43])
        +A[23]*conj(B[44])
        +A[24]*conj(B[45])
        +A[25]*conj(B[46])
        +A[26]*conj(B[47])
        +A[27]*conj(B[48]);

  C[28] = A[28]*conj(B[0])
        +A[29]*conj(B[1])
        +A[30]*conj(B[2])
        +A[31]*conj(B[3])
        +A[32]*conj(B[4])
        +A[33]*conj(B[5])
        +A[34]*conj(B[6]);

  C[29] = A[28]*conj(B[7])
        +A[29]*conj(B[8])
        +A[30]*conj(B[9])
        +A[31]*conj(B[10])
        +A[32]*conj(B[11])
        +A[33]*conj(B[12])
        +A[34]*conj(B[13]);

  C[30] = A[28]*conj(B[14])
        +A[29]*conj(B[15])
        +A[30]*conj(B[16])
        +A[31]*conj(B[17])
        +A[32]*conj(B[18])
        +A[33]*conj(B[19])
        +A[34]*conj(B[20]);

  C[31] = A[28]*conj(B[21])
        +A[29]*conj(B[22])
        +A[30]*conj(B[23])
        +A[31]*conj(B[24])
        +A[32]*conj(B[25])
        +A[33]*conj(B[26])
        +A[34]*conj(B[27]);

  C[32] = A[28]*conj(B[28])
        +A[29]*conj(B[29])
        +A[30]*conj(B[30])
        +A[31]*conj(B[31])
        +A[32]*conj(B[32])
        +A[33]*conj(B[33])
        +A[34]*conj(B[34]);

  C[33] = A[28]*conj(B[35])
        +A[29]*conj(B[36])
        +A[30]*conj(B[37])
        +A[31]*conj(B[38])
        +A[32]*conj(B[39])
        +A[33]*conj(B[40])
        +A[34]*conj(B[41]);

  C[34] = A[28]*conj(B[42])
        +A[29]*conj(B[43])
        +A[30]*conj(B[44])
        +A[31]*conj(B[45])
        +A[32]*conj(B[46])
        +A[33]*conj(B[47])
        +A[34]*conj(B[48]);

  C[35] = A[35]*conj(B[0])
        +A[36]*conj(B[1])
        +A[37]*conj(B[2])
        +A[38]*conj(B[3])
        +A[39]*conj(B[4])
        +A[40]*conj(B[5])
        +A[41]*conj(B[6]);

  C[36] = A[35]*conj(B[7])
        +A[36]*conj(B[8])
        +A[37]*conj(B[9])
        +A[38]*conj(B[10])
        +A[39]*conj(B[11])
        +A[40]*conj(B[12])
        +A[41]*conj(B[13]);

  C[37] = A[35]*conj(B[14])
        +A[36]*conj(B[15])
        +A[37]*conj(B[16])
        +A[38]*conj(B[17])
        +A[39]*conj(B[18])
        +A[40]*conj(B[19])
        +A[41]*conj(B[20]);

  C[38] = A[35]*conj(B[21])
        +A[36]*conj(B[22])
        +A[37]*conj(B[23])
        +A[38]*conj(B[24])
        +A[39]*conj(B[25])
        +A[40]*conj(B[26])
        +A[41]*conj(B[27]);

  C[39] = A[35]*conj(B[28])
        +A[36]*conj(B[29])
        +A[37]*conj(B[30])
        +A[38]*conj(B[31])
        +A[39]*conj(B[32])
        +A[40]*conj(B[33])
        +A[41]*conj(B[34]);

  C[40] = A[35]*conj(B[35])
        +A[36]*conj(B[36])
        +A[37]*conj(B[37])
        +A[38]*conj(B[38])
        +A[39]*conj(B[39])
        +A[40]*conj(B[40])
        +A[41]*conj(B[41]);

  C[41] = A[35]*conj(B[42])
        +A[36]*conj(B[43])
        +A[37]*conj(B[44])
        +A[38]*conj(B[45])
        +A[39]*conj(B[46])
        +A[40]*conj(B[47])
        +A[41]*conj(B[48]);

  C[42] = A[42]*conj(B[0])
        +A[43]*conj(B[1])
        +A[44]*conj(B[2])
        +A[45]*conj(B[3])
        +A[46]*conj(B[4])
        +A[47]*conj(B[5])
        +A[48]*conj(B[6]);

  C[43] = A[42]*conj(B[7])
        +A[43]*conj(B[8])
        +A[44]*conj(B[9])
        +A[45]*conj(B[10])
        +A[46]*conj(B[11])
        +A[47]*conj(B[12])
        +A[48]*conj(B[13]);

  C[44] = A[42]*conj(B[14])
        +A[43]*conj(B[15])
        +A[44]*conj(B[16])
        +A[45]*conj(B[17])
        +A[46]*conj(B[18])
        +A[47]*conj(B[19])
        +A[48]*conj(B[20]);

  C[45] = A[42]*conj(B[21])
        +A[43]*conj(B[22])
        +A[44]*conj(B[23])
        +A[45]*conj(B[24])
        +A[46]*conj(B[25])
        +A[47]*conj(B[26])
        +A[48]*conj(B[27]);

  C[46] = A[42]*conj(B[28])
        +A[43]*conj(B[29])
        +A[44]*conj(B[30])
        +A[45]*conj(B[31])
        +A[46]*conj(B[32])
        +A[47]*conj(B[33])
        +A[48]*conj(B[34]);

  C[47] = A[42]*conj(B[35])
        +A[43]*conj(B[36])
        +A[44]*conj(B[37])
        +A[45]*conj(B[38])
        +A[46]*conj(B[39])
        +A[47]*conj(B[40])
        +A[48]*conj(B[41]);

  C[48] = A[42]*conj(B[42])
        +A[43]*conj(B[43])
        +A[44]*conj(B[44])
        +A[45]*conj(B[45])
        +A[46]*conj(B[46])
        +A[47]*conj(B[47])
        +A[48]*conj(B[48]);

};


inline void mult_C_equals_AdaggerB_for_SU7(dc *C, dc *A, dc *B) {

  C[0] = conj(A[0])*B[0]
        +conj(A[7])*B[7]
        +conj(A[14])*B[14]
        +conj(A[21])*B[21]
        +conj(A[28])*B[28]
        +conj(A[35])*B[35]
        +conj(A[42])*B[42];

  C[1] = conj(A[0])*B[1]
        +conj(A[7])*B[8]
        +conj(A[14])*B[15]
        +conj(A[21])*B[22]
        +conj(A[28])*B[29]
        +conj(A[35])*B[36]
        +conj(A[42])*B[43];

  C[2] = conj(A[0])*B[2]
        +conj(A[7])*B[9]
        +conj(A[14])*B[16]
        +conj(A[21])*B[23]
        +conj(A[28])*B[30]
        +conj(A[35])*B[37]
        +conj(A[42])*B[44];

  C[3] = conj(A[0])*B[3]
        +conj(A[7])*B[10]
        +conj(A[14])*B[17]
        +conj(A[21])*B[24]
        +conj(A[28])*B[31]
        +conj(A[35])*B[38]
        +conj(A[42])*B[45];

  C[4] = conj(A[0])*B[4]
        +conj(A[7])*B[11]
        +conj(A[14])*B[18]
        +conj(A[21])*B[25]
        +conj(A[28])*B[32]
        +conj(A[35])*B[39]
        +conj(A[42])*B[46];

  C[5] = conj(A[0])*B[5]
        +conj(A[7])*B[12]
        +conj(A[14])*B[19]
        +conj(A[21])*B[26]
        +conj(A[28])*B[33]
        +conj(A[35])*B[40]
        +conj(A[42])*B[47];

  C[6] = conj(A[0])*B[6]
        +conj(A[7])*B[13]
        +conj(A[14])*B[20]
        +conj(A[21])*B[27]
        +conj(A[28])*B[34]
        +conj(A[35])*B[41]
        +conj(A[42])*B[48];

  C[7] = conj(A[1])*B[0]
        +conj(A[8])*B[7]
        +conj(A[15])*B[14]
        +conj(A[22])*B[21]
        +conj(A[29])*B[28]
        +conj(A[36])*B[35]
        +conj(A[43])*B[42];

  C[8] = conj(A[1])*B[1]
        +conj(A[8])*B[8]
        +conj(A[15])*B[15]
        +conj(A[22])*B[22]
        +conj(A[29])*B[29]
        +conj(A[36])*B[36]
        +conj(A[43])*B[43];

  C[9] = conj(A[1])*B[2]
        +conj(A[8])*B[9]
        +conj(A[15])*B[16]
        +conj(A[22])*B[23]
        +conj(A[29])*B[30]
        +conj(A[36])*B[37]
        +conj(A[43])*B[44];

  C[10] = conj(A[1])*B[3]
        +conj(A[8])*B[10]
        +conj(A[15])*B[17]
        +conj(A[22])*B[24]
        +conj(A[29])*B[31]
        +conj(A[36])*B[38]
        +conj(A[43])*B[45];

  C[11] = conj(A[1])*B[4]
        +conj(A[8])*B[11]
        +conj(A[15])*B[18]
        +conj(A[22])*B[25]
        +conj(A[29])*B[32]
        +conj(A[36])*B[39]
        +conj(A[43])*B[46];

  C[12] = conj(A[1])*B[5]
        +conj(A[8])*B[12]
        +conj(A[15])*B[19]
        +conj(A[22])*B[26]
        +conj(A[29])*B[33]
        +conj(A[36])*B[40]
        +conj(A[43])*B[47];

  C[13] = conj(A[1])*B[6]
        +conj(A[8])*B[13]
        +conj(A[15])*B[20]
        +conj(A[22])*B[27]
        +conj(A[29])*B[34]
        +conj(A[36])*B[41]
        +conj(A[43])*B[48];

  C[14] = conj(A[2])*B[0]
        +conj(A[9])*B[7]
        +conj(A[16])*B[14]
        +conj(A[23])*B[21]
        +conj(A[30])*B[28]
        +conj(A[37])*B[35]
        +conj(A[44])*B[42];

  C[15] = conj(A[2])*B[1]
        +conj(A[9])*B[8]
        +conj(A[16])*B[15]
        +conj(A[23])*B[22]
        +conj(A[30])*B[29]
        +conj(A[37])*B[36]
        +conj(A[44])*B[43];

  C[16] = conj(A[2])*B[2]
        +conj(A[9])*B[9]
        +conj(A[16])*B[16]
        +conj(A[23])*B[23]
        +conj(A[30])*B[30]
        +conj(A[37])*B[37]
        +conj(A[44])*B[44];

  C[17] = conj(A[2])*B[3]
        +conj(A[9])*B[10]
        +conj(A[16])*B[17]
        +conj(A[23])*B[24]
        +conj(A[30])*B[31]
        +conj(A[37])*B[38]
        +conj(A[44])*B[45];

  C[18] = conj(A[2])*B[4]
        +conj(A[9])*B[11]
        +conj(A[16])*B[18]
        +conj(A[23])*B[25]
        +conj(A[30])*B[32]
        +conj(A[37])*B[39]
        +conj(A[44])*B[46];

  C[19] = conj(A[2])*B[5]
        +conj(A[9])*B[12]
        +conj(A[16])*B[19]
        +conj(A[23])*B[26]
        +conj(A[30])*B[33]
        +conj(A[37])*B[40]
        +conj(A[44])*B[47];

  C[20] = conj(A[2])*B[6]
        +conj(A[9])*B[13]
        +conj(A[16])*B[20]
        +conj(A[23])*B[27]
        +conj(A[30])*B[34]
        +conj(A[37])*B[41]
        +conj(A[44])*B[48];

  C[21] = conj(A[3])*B[0]
        +conj(A[10])*B[7]
        +conj(A[17])*B[14]
        +conj(A[24])*B[21]
        +conj(A[31])*B[28]
        +conj(A[38])*B[35]
        +conj(A[45])*B[42];

  C[22] = conj(A[3])*B[1]
        +conj(A[10])*B[8]
        +conj(A[17])*B[15]
        +conj(A[24])*B[22]
        +conj(A[31])*B[29]
        +conj(A[38])*B[36]
        +conj(A[45])*B[43];

  C[23] = conj(A[3])*B[2]
        +conj(A[10])*B[9]
        +conj(A[17])*B[16]
        +conj(A[24])*B[23]
        +conj(A[31])*B[30]
        +conj(A[38])*B[37]
        +conj(A[45])*B[44];

  C[24] = conj(A[3])*B[3]
        +conj(A[10])*B[10]
        +conj(A[17])*B[17]
        +conj(A[24])*B[24]
        +conj(A[31])*B[31]
        +conj(A[38])*B[38]
        +conj(A[45])*B[45];

  C[25] = conj(A[3])*B[4]
        +conj(A[10])*B[11]
        +conj(A[17])*B[18]
        +conj(A[24])*B[25]
        +conj(A[31])*B[32]
        +conj(A[38])*B[39]
        +conj(A[45])*B[46];

  C[26] = conj(A[3])*B[5]
        +conj(A[10])*B[12]
        +conj(A[17])*B[19]
        +conj(A[24])*B[26]
        +conj(A[31])*B[33]
        +conj(A[38])*B[40]
        +conj(A[45])*B[47];

  C[27] = conj(A[3])*B[6]
        +conj(A[10])*B[13]
        +conj(A[17])*B[20]
        +conj(A[24])*B[27]
        +conj(A[31])*B[34]
        +conj(A[38])*B[41]
        +conj(A[45])*B[48];

  C[28] = conj(A[4])*B[0]
        +conj(A[11])*B[7]
        +conj(A[18])*B[14]
        +conj(A[25])*B[21]
        +conj(A[32])*B[28]
        +conj(A[39])*B[35]
        +conj(A[46])*B[42];

  C[29] = conj(A[4])*B[1]
        +conj(A[11])*B[8]
        +conj(A[18])*B[15]
        +conj(A[25])*B[22]
        +conj(A[32])*B[29]
        +conj(A[39])*B[36]
        +conj(A[46])*B[43];

  C[30] = conj(A[4])*B[2]
        +conj(A[11])*B[9]
        +conj(A[18])*B[16]
        +conj(A[25])*B[23]
        +conj(A[32])*B[30]
        +conj(A[39])*B[37]
        +conj(A[46])*B[44];

  C[31] = conj(A[4])*B[3]
        +conj(A[11])*B[10]
        +conj(A[18])*B[17]
        +conj(A[25])*B[24]
        +conj(A[32])*B[31]
        +conj(A[39])*B[38]
        +conj(A[46])*B[45];

  C[32] = conj(A[4])*B[4]
        +conj(A[11])*B[11]
        +conj(A[18])*B[18]
        +conj(A[25])*B[25]
        +conj(A[32])*B[32]
        +conj(A[39])*B[39]
        +conj(A[46])*B[46];

  C[33] = conj(A[4])*B[5]
        +conj(A[11])*B[12]
        +conj(A[18])*B[19]
        +conj(A[25])*B[26]
        +conj(A[32])*B[33]
        +conj(A[39])*B[40]
        +conj(A[46])*B[47];

  C[34] = conj(A[4])*B[6]
        +conj(A[11])*B[13]
        +conj(A[18])*B[20]
        +conj(A[25])*B[27]
        +conj(A[32])*B[34]
        +conj(A[39])*B[41]
        +conj(A[46])*B[48];

  C[35] = conj(A[5])*B[0]
        +conj(A[12])*B[7]
        +conj(A[19])*B[14]
        +conj(A[26])*B[21]
        +conj(A[33])*B[28]
        +conj(A[40])*B[35]
        +conj(A[47])*B[42];

  C[36] = conj(A[5])*B[1]
        +conj(A[12])*B[8]
        +conj(A[19])*B[15]
        +conj(A[26])*B[22]
        +conj(A[33])*B[29]
        +conj(A[40])*B[36]
        +conj(A[47])*B[43];

  C[37] = conj(A[5])*B[2]
        +conj(A[12])*B[9]
        +conj(A[19])*B[16]
        +conj(A[26])*B[23]
        +conj(A[33])*B[30]
        +conj(A[40])*B[37]
        +conj(A[47])*B[44];

  C[38] = conj(A[5])*B[3]
        +conj(A[12])*B[10]
        +conj(A[19])*B[17]
        +conj(A[26])*B[24]
        +conj(A[33])*B[31]
        +conj(A[40])*B[38]
        +conj(A[47])*B[45];

  C[39] = conj(A[5])*B[4]
        +conj(A[12])*B[11]
        +conj(A[19])*B[18]
        +conj(A[26])*B[25]
        +conj(A[33])*B[32]
        +conj(A[40])*B[39]
        +conj(A[47])*B[46];

  C[40] = conj(A[5])*B[5]
        +conj(A[12])*B[12]
        +conj(A[19])*B[19]
        +conj(A[26])*B[26]
        +conj(A[33])*B[33]
        +conj(A[40])*B[40]
        +conj(A[47])*B[47];

  C[41] = conj(A[5])*B[6]
        +conj(A[12])*B[13]
        +conj(A[19])*B[20]
        +conj(A[26])*B[27]
        +conj(A[33])*B[34]
        +conj(A[40])*B[41]
        +conj(A[47])*B[48];

  C[42] = conj(A[6])*B[0]
        +conj(A[13])*B[7]
        +conj(A[20])*B[14]
        +conj(A[27])*B[21]
        +conj(A[34])*B[28]
        +conj(A[41])*B[35]
        +conj(A[48])*B[42];

  C[43] = conj(A[6])*B[1]
        +conj(A[13])*B[8]
        +conj(A[20])*B[15]
        +conj(A[27])*B[22]
        +conj(A[34])*B[29]
        +conj(A[41])*B[36]
        +conj(A[48])*B[43];

  C[44] = conj(A[6])*B[2]
        +conj(A[13])*B[9]
        +conj(A[20])*B[16]
        +conj(A[27])*B[23]
        +conj(A[34])*B[30]
        +conj(A[41])*B[37]
        +conj(A[48])*B[44];

  C[45] = conj(A[6])*B[3]
        +conj(A[13])*B[10]
        +conj(A[20])*B[17]
        +conj(A[27])*B[24]
        +conj(A[34])*B[31]
        +conj(A[41])*B[38]
        +conj(A[48])*B[45];

  C[46] = conj(A[6])*B[4]
        +conj(A[13])*B[11]
        +conj(A[20])*B[18]
        +conj(A[27])*B[25]
        +conj(A[34])*B[32]
        +conj(A[41])*B[39]
        +conj(A[48])*B[46];

  C[47] = conj(A[6])*B[5]
        +conj(A[13])*B[12]
        +conj(A[20])*B[19]
        +conj(A[27])*B[26]
        +conj(A[34])*B[33]
        +conj(A[41])*B[40]
        +conj(A[48])*B[47];

  C[48] = conj(A[6])*B[6]
        +conj(A[13])*B[13]
        +conj(A[20])*B[20]
        +conj(A[27])*B[27]
        +conj(A[34])*B[34]
        +conj(A[41])*B[41]
        +conj(A[48])*B[48];

};

