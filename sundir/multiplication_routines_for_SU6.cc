inline void mult_C_equals_AB_for_SU6(dc *C, dc *A, dc *B) {

  C[0] = A[0]*B[0]
        +A[1]*B[6]
        +A[2]*B[12]
        +A[3]*B[18]
        +A[4]*B[24]
        +A[5]*B[30];

  C[1] = A[0]*B[1]
        +A[1]*B[7]
        +A[2]*B[13]
        +A[3]*B[19]
        +A[4]*B[25]
        +A[5]*B[31];

  C[2] = A[0]*B[2]
        +A[1]*B[8]
        +A[2]*B[14]
        +A[3]*B[20]
        +A[4]*B[26]
        +A[5]*B[32];

  C[3] = A[0]*B[3]
        +A[1]*B[9]
        +A[2]*B[15]
        +A[3]*B[21]
        +A[4]*B[27]
        +A[5]*B[33];

  C[4] = A[0]*B[4]
        +A[1]*B[10]
        +A[2]*B[16]
        +A[3]*B[22]
        +A[4]*B[28]
        +A[5]*B[34];

  C[5] = A[0]*B[5]
        +A[1]*B[11]
        +A[2]*B[17]
        +A[3]*B[23]
        +A[4]*B[29]
        +A[5]*B[35];

  C[6] = A[6]*B[0]
        +A[7]*B[6]
        +A[8]*B[12]
        +A[9]*B[18]
        +A[10]*B[24]
        +A[11]*B[30];

  C[7] = A[6]*B[1]
        +A[7]*B[7]
        +A[8]*B[13]
        +A[9]*B[19]
        +A[10]*B[25]
        +A[11]*B[31];

  C[8] = A[6]*B[2]
        +A[7]*B[8]
        +A[8]*B[14]
        +A[9]*B[20]
        +A[10]*B[26]
        +A[11]*B[32];

  C[9] = A[6]*B[3]
        +A[7]*B[9]
        +A[8]*B[15]
        +A[9]*B[21]
        +A[10]*B[27]
        +A[11]*B[33];

  C[10] = A[6]*B[4]
        +A[7]*B[10]
        +A[8]*B[16]
        +A[9]*B[22]
        +A[10]*B[28]
        +A[11]*B[34];

  C[11] = A[6]*B[5]
        +A[7]*B[11]
        +A[8]*B[17]
        +A[9]*B[23]
        +A[10]*B[29]
        +A[11]*B[35];

  C[12] = A[12]*B[0]
        +A[13]*B[6]
        +A[14]*B[12]
        +A[15]*B[18]
        +A[16]*B[24]
        +A[17]*B[30];

  C[13] = A[12]*B[1]
        +A[13]*B[7]
        +A[14]*B[13]
        +A[15]*B[19]
        +A[16]*B[25]
        +A[17]*B[31];

  C[14] = A[12]*B[2]
        +A[13]*B[8]
        +A[14]*B[14]
        +A[15]*B[20]
        +A[16]*B[26]
        +A[17]*B[32];

  C[15] = A[12]*B[3]
        +A[13]*B[9]
        +A[14]*B[15]
        +A[15]*B[21]
        +A[16]*B[27]
        +A[17]*B[33];

  C[16] = A[12]*B[4]
        +A[13]*B[10]
        +A[14]*B[16]
        +A[15]*B[22]
        +A[16]*B[28]
        +A[17]*B[34];

  C[17] = A[12]*B[5]
        +A[13]*B[11]
        +A[14]*B[17]
        +A[15]*B[23]
        +A[16]*B[29]
        +A[17]*B[35];

  C[18] = A[18]*B[0]
        +A[19]*B[6]
        +A[20]*B[12]
        +A[21]*B[18]
        +A[22]*B[24]
        +A[23]*B[30];

  C[19] = A[18]*B[1]
        +A[19]*B[7]
        +A[20]*B[13]
        +A[21]*B[19]
        +A[22]*B[25]
        +A[23]*B[31];

  C[20] = A[18]*B[2]
        +A[19]*B[8]
        +A[20]*B[14]
        +A[21]*B[20]
        +A[22]*B[26]
        +A[23]*B[32];

  C[21] = A[18]*B[3]
        +A[19]*B[9]
        +A[20]*B[15]
        +A[21]*B[21]
        +A[22]*B[27]
        +A[23]*B[33];

  C[22] = A[18]*B[4]
        +A[19]*B[10]
        +A[20]*B[16]
        +A[21]*B[22]
        +A[22]*B[28]
        +A[23]*B[34];

  C[23] = A[18]*B[5]
        +A[19]*B[11]
        +A[20]*B[17]
        +A[21]*B[23]
        +A[22]*B[29]
        +A[23]*B[35];

  C[24] = A[24]*B[0]
        +A[25]*B[6]
        +A[26]*B[12]
        +A[27]*B[18]
        +A[28]*B[24]
        +A[29]*B[30];

  C[25] = A[24]*B[1]
        +A[25]*B[7]
        +A[26]*B[13]
        +A[27]*B[19]
        +A[28]*B[25]
        +A[29]*B[31];

  C[26] = A[24]*B[2]
        +A[25]*B[8]
        +A[26]*B[14]
        +A[27]*B[20]
        +A[28]*B[26]
        +A[29]*B[32];

  C[27] = A[24]*B[3]
        +A[25]*B[9]
        +A[26]*B[15]
        +A[27]*B[21]
        +A[28]*B[27]
        +A[29]*B[33];

  C[28] = A[24]*B[4]
        +A[25]*B[10]
        +A[26]*B[16]
        +A[27]*B[22]
        +A[28]*B[28]
        +A[29]*B[34];

  C[29] = A[24]*B[5]
        +A[25]*B[11]
        +A[26]*B[17]
        +A[27]*B[23]
        +A[28]*B[29]
        +A[29]*B[35];

  C[30] = A[30]*B[0]
        +A[31]*B[6]
        +A[32]*B[12]
        +A[33]*B[18]
        +A[34]*B[24]
        +A[35]*B[30];

  C[31] = A[30]*B[1]
        +A[31]*B[7]
        +A[32]*B[13]
        +A[33]*B[19]
        +A[34]*B[25]
        +A[35]*B[31];

  C[32] = A[30]*B[2]
        +A[31]*B[8]
        +A[32]*B[14]
        +A[33]*B[20]
        +A[34]*B[26]
        +A[35]*B[32];

  C[33] = A[30]*B[3]
        +A[31]*B[9]
        +A[32]*B[15]
        +A[33]*B[21]
        +A[34]*B[27]
        +A[35]*B[33];

  C[34] = A[30]*B[4]
        +A[31]*B[10]
        +A[32]*B[16]
        +A[33]*B[22]
        +A[34]*B[28]
        +A[35]*B[34];

  C[35] = A[30]*B[5]
        +A[31]*B[11]
        +A[32]*B[17]
        +A[33]*B[23]
        +A[34]*B[29]
        +A[35]*B[35];

};


inline void mult_C_equals_ABdagger_for_SU6(dc *C, dc *A, dc *B) {

  C[0] = A[0]*conj(B[0])
        +A[1]*conj(B[1])
        +A[2]*conj(B[2])
        +A[3]*conj(B[3])
        +A[4]*conj(B[4])
        +A[5]*conj(B[5]);

  C[1] = A[0]*conj(B[6])
        +A[1]*conj(B[7])
        +A[2]*conj(B[8])
        +A[3]*conj(B[9])
        +A[4]*conj(B[10])
        +A[5]*conj(B[11]);

  C[2] = A[0]*conj(B[12])
        +A[1]*conj(B[13])
        +A[2]*conj(B[14])
        +A[3]*conj(B[15])
        +A[4]*conj(B[16])
        +A[5]*conj(B[17]);

  C[3] = A[0]*conj(B[18])
        +A[1]*conj(B[19])
        +A[2]*conj(B[20])
        +A[3]*conj(B[21])
        +A[4]*conj(B[22])
        +A[5]*conj(B[23]);

  C[4] = A[0]*conj(B[24])
        +A[1]*conj(B[25])
        +A[2]*conj(B[26])
        +A[3]*conj(B[27])
        +A[4]*conj(B[28])
        +A[5]*conj(B[29]);

  C[5] = A[0]*conj(B[30])
        +A[1]*conj(B[31])
        +A[2]*conj(B[32])
        +A[3]*conj(B[33])
        +A[4]*conj(B[34])
        +A[5]*conj(B[35]);

  C[6] = A[6]*conj(B[0])
        +A[7]*conj(B[1])
        +A[8]*conj(B[2])
        +A[9]*conj(B[3])
        +A[10]*conj(B[4])
        +A[11]*conj(B[5]);

  C[7] = A[6]*conj(B[6])
        +A[7]*conj(B[7])
        +A[8]*conj(B[8])
        +A[9]*conj(B[9])
        +A[10]*conj(B[10])
        +A[11]*conj(B[11]);

  C[8] = A[6]*conj(B[12])
        +A[7]*conj(B[13])
        +A[8]*conj(B[14])
        +A[9]*conj(B[15])
        +A[10]*conj(B[16])
        +A[11]*conj(B[17]);

  C[9] = A[6]*conj(B[18])
        +A[7]*conj(B[19])
        +A[8]*conj(B[20])
        +A[9]*conj(B[21])
        +A[10]*conj(B[22])
        +A[11]*conj(B[23]);

  C[10] = A[6]*conj(B[24])
        +A[7]*conj(B[25])
        +A[8]*conj(B[26])
        +A[9]*conj(B[27])
        +A[10]*conj(B[28])
        +A[11]*conj(B[29]);

  C[11] = A[6]*conj(B[30])
        +A[7]*conj(B[31])
        +A[8]*conj(B[32])
        +A[9]*conj(B[33])
        +A[10]*conj(B[34])
        +A[11]*conj(B[35]);

  C[12] = A[12]*conj(B[0])
        +A[13]*conj(B[1])
        +A[14]*conj(B[2])
        +A[15]*conj(B[3])
        +A[16]*conj(B[4])
        +A[17]*conj(B[5]);

  C[13] = A[12]*conj(B[6])
        +A[13]*conj(B[7])
        +A[14]*conj(B[8])
        +A[15]*conj(B[9])
        +A[16]*conj(B[10])
        +A[17]*conj(B[11]);

  C[14] = A[12]*conj(B[12])
        +A[13]*conj(B[13])
        +A[14]*conj(B[14])
        +A[15]*conj(B[15])
        +A[16]*conj(B[16])
        +A[17]*conj(B[17]);

  C[15] = A[12]*conj(B[18])
        +A[13]*conj(B[19])
        +A[14]*conj(B[20])
        +A[15]*conj(B[21])
        +A[16]*conj(B[22])
        +A[17]*conj(B[23]);

  C[16] = A[12]*conj(B[24])
        +A[13]*conj(B[25])
        +A[14]*conj(B[26])
        +A[15]*conj(B[27])
        +A[16]*conj(B[28])
        +A[17]*conj(B[29]);

  C[17] = A[12]*conj(B[30])
        +A[13]*conj(B[31])
        +A[14]*conj(B[32])
        +A[15]*conj(B[33])
        +A[16]*conj(B[34])
        +A[17]*conj(B[35]);

  C[18] = A[18]*conj(B[0])
        +A[19]*conj(B[1])
        +A[20]*conj(B[2])
        +A[21]*conj(B[3])
        +A[22]*conj(B[4])
        +A[23]*conj(B[5]);

  C[19] = A[18]*conj(B[6])
        +A[19]*conj(B[7])
        +A[20]*conj(B[8])
        +A[21]*conj(B[9])
        +A[22]*conj(B[10])
        +A[23]*conj(B[11]);

  C[20] = A[18]*conj(B[12])
        +A[19]*conj(B[13])
        +A[20]*conj(B[14])
        +A[21]*conj(B[15])
        +A[22]*conj(B[16])
        +A[23]*conj(B[17]);

  C[21] = A[18]*conj(B[18])
        +A[19]*conj(B[19])
        +A[20]*conj(B[20])
        +A[21]*conj(B[21])
        +A[22]*conj(B[22])
        +A[23]*conj(B[23]);

  C[22] = A[18]*conj(B[24])
        +A[19]*conj(B[25])
        +A[20]*conj(B[26])
        +A[21]*conj(B[27])
        +A[22]*conj(B[28])
        +A[23]*conj(B[29]);

  C[23] = A[18]*conj(B[30])
        +A[19]*conj(B[31])
        +A[20]*conj(B[32])
        +A[21]*conj(B[33])
        +A[22]*conj(B[34])
        +A[23]*conj(B[35]);

  C[24] = A[24]*conj(B[0])
        +A[25]*conj(B[1])
        +A[26]*conj(B[2])
        +A[27]*conj(B[3])
        +A[28]*conj(B[4])
        +A[29]*conj(B[5]);

  C[25] = A[24]*conj(B[6])
        +A[25]*conj(B[7])
        +A[26]*conj(B[8])
        +A[27]*conj(B[9])
        +A[28]*conj(B[10])
        +A[29]*conj(B[11]);

  C[26] = A[24]*conj(B[12])
        +A[25]*conj(B[13])
        +A[26]*conj(B[14])
        +A[27]*conj(B[15])
        +A[28]*conj(B[16])
        +A[29]*conj(B[17]);

  C[27] = A[24]*conj(B[18])
        +A[25]*conj(B[19])
        +A[26]*conj(B[20])
        +A[27]*conj(B[21])
        +A[28]*conj(B[22])
        +A[29]*conj(B[23]);

  C[28] = A[24]*conj(B[24])
        +A[25]*conj(B[25])
        +A[26]*conj(B[26])
        +A[27]*conj(B[27])
        +A[28]*conj(B[28])
        +A[29]*conj(B[29]);

  C[29] = A[24]*conj(B[30])
        +A[25]*conj(B[31])
        +A[26]*conj(B[32])
        +A[27]*conj(B[33])
        +A[28]*conj(B[34])
        +A[29]*conj(B[35]);

  C[30] = A[30]*conj(B[0])
        +A[31]*conj(B[1])
        +A[32]*conj(B[2])
        +A[33]*conj(B[3])
        +A[34]*conj(B[4])
        +A[35]*conj(B[5]);

  C[31] = A[30]*conj(B[6])
        +A[31]*conj(B[7])
        +A[32]*conj(B[8])
        +A[33]*conj(B[9])
        +A[34]*conj(B[10])
        +A[35]*conj(B[11]);

  C[32] = A[30]*conj(B[12])
        +A[31]*conj(B[13])
        +A[32]*conj(B[14])
        +A[33]*conj(B[15])
        +A[34]*conj(B[16])
        +A[35]*conj(B[17]);

  C[33] = A[30]*conj(B[18])
        +A[31]*conj(B[19])
        +A[32]*conj(B[20])
        +A[33]*conj(B[21])
        +A[34]*conj(B[22])
        +A[35]*conj(B[23]);

  C[34] = A[30]*conj(B[24])
        +A[31]*conj(B[25])
        +A[32]*conj(B[26])
        +A[33]*conj(B[27])
        +A[34]*conj(B[28])
        +A[35]*conj(B[29]);

  C[35] = A[30]*conj(B[30])
        +A[31]*conj(B[31])
        +A[32]*conj(B[32])
        +A[33]*conj(B[33])
        +A[34]*conj(B[34])
        +A[35]*conj(B[35]);

};


inline void mult_C_equals_AdaggerB_for_SU6(dc *C, dc *A, dc *B) {

  C[0] = conj(A[0])*B[0]
        +conj(A[6])*B[6]
        +conj(A[12])*B[12]
        +conj(A[18])*B[18]
        +conj(A[24])*B[24]
        +conj(A[30])*B[30];

  C[1] = conj(A[0])*B[1]
        +conj(A[6])*B[7]
        +conj(A[12])*B[13]
        +conj(A[18])*B[19]
        +conj(A[24])*B[25]
        +conj(A[30])*B[31];

  C[2] = conj(A[0])*B[2]
        +conj(A[6])*B[8]
        +conj(A[12])*B[14]
        +conj(A[18])*B[20]
        +conj(A[24])*B[26]
        +conj(A[30])*B[32];

  C[3] = conj(A[0])*B[3]
        +conj(A[6])*B[9]
        +conj(A[12])*B[15]
        +conj(A[18])*B[21]
        +conj(A[24])*B[27]
        +conj(A[30])*B[33];

  C[4] = conj(A[0])*B[4]
        +conj(A[6])*B[10]
        +conj(A[12])*B[16]
        +conj(A[18])*B[22]
        +conj(A[24])*B[28]
        +conj(A[30])*B[34];

  C[5] = conj(A[0])*B[5]
        +conj(A[6])*B[11]
        +conj(A[12])*B[17]
        +conj(A[18])*B[23]
        +conj(A[24])*B[29]
        +conj(A[30])*B[35];

  C[6] = conj(A[1])*B[0]
        +conj(A[7])*B[6]
        +conj(A[13])*B[12]
        +conj(A[19])*B[18]
        +conj(A[25])*B[24]
        +conj(A[31])*B[30];

  C[7] = conj(A[1])*B[1]
        +conj(A[7])*B[7]
        +conj(A[13])*B[13]
        +conj(A[19])*B[19]
        +conj(A[25])*B[25]
        +conj(A[31])*B[31];

  C[8] = conj(A[1])*B[2]
        +conj(A[7])*B[8]
        +conj(A[13])*B[14]
        +conj(A[19])*B[20]
        +conj(A[25])*B[26]
        +conj(A[31])*B[32];

  C[9] = conj(A[1])*B[3]
        +conj(A[7])*B[9]
        +conj(A[13])*B[15]
        +conj(A[19])*B[21]
        +conj(A[25])*B[27]
        +conj(A[31])*B[33];

  C[10] = conj(A[1])*B[4]
        +conj(A[7])*B[10]
        +conj(A[13])*B[16]
        +conj(A[19])*B[22]
        +conj(A[25])*B[28]
        +conj(A[31])*B[34];

  C[11] = conj(A[1])*B[5]
        +conj(A[7])*B[11]
        +conj(A[13])*B[17]
        +conj(A[19])*B[23]
        +conj(A[25])*B[29]
        +conj(A[31])*B[35];

  C[12] = conj(A[2])*B[0]
        +conj(A[8])*B[6]
        +conj(A[14])*B[12]
        +conj(A[20])*B[18]
        +conj(A[26])*B[24]
        +conj(A[32])*B[30];

  C[13] = conj(A[2])*B[1]
        +conj(A[8])*B[7]
        +conj(A[14])*B[13]
        +conj(A[20])*B[19]
        +conj(A[26])*B[25]
        +conj(A[32])*B[31];

  C[14] = conj(A[2])*B[2]
        +conj(A[8])*B[8]
        +conj(A[14])*B[14]
        +conj(A[20])*B[20]
        +conj(A[26])*B[26]
        +conj(A[32])*B[32];

  C[15] = conj(A[2])*B[3]
        +conj(A[8])*B[9]
        +conj(A[14])*B[15]
        +conj(A[20])*B[21]
        +conj(A[26])*B[27]
        +conj(A[32])*B[33];

  C[16] = conj(A[2])*B[4]
        +conj(A[8])*B[10]
        +conj(A[14])*B[16]
        +conj(A[20])*B[22]
        +conj(A[26])*B[28]
        +conj(A[32])*B[34];

  C[17] = conj(A[2])*B[5]
        +conj(A[8])*B[11]
        +conj(A[14])*B[17]
        +conj(A[20])*B[23]
        +conj(A[26])*B[29]
        +conj(A[32])*B[35];

  C[18] = conj(A[3])*B[0]
        +conj(A[9])*B[6]
        +conj(A[15])*B[12]
        +conj(A[21])*B[18]
        +conj(A[27])*B[24]
        +conj(A[33])*B[30];

  C[19] = conj(A[3])*B[1]
        +conj(A[9])*B[7]
        +conj(A[15])*B[13]
        +conj(A[21])*B[19]
        +conj(A[27])*B[25]
        +conj(A[33])*B[31];

  C[20] = conj(A[3])*B[2]
        +conj(A[9])*B[8]
        +conj(A[15])*B[14]
        +conj(A[21])*B[20]
        +conj(A[27])*B[26]
        +conj(A[33])*B[32];

  C[21] = conj(A[3])*B[3]
        +conj(A[9])*B[9]
        +conj(A[15])*B[15]
        +conj(A[21])*B[21]
        +conj(A[27])*B[27]
        +conj(A[33])*B[33];

  C[22] = conj(A[3])*B[4]
        +conj(A[9])*B[10]
        +conj(A[15])*B[16]
        +conj(A[21])*B[22]
        +conj(A[27])*B[28]
        +conj(A[33])*B[34];

  C[23] = conj(A[3])*B[5]
        +conj(A[9])*B[11]
        +conj(A[15])*B[17]
        +conj(A[21])*B[23]
        +conj(A[27])*B[29]
        +conj(A[33])*B[35];

  C[24] = conj(A[4])*B[0]
        +conj(A[10])*B[6]
        +conj(A[16])*B[12]
        +conj(A[22])*B[18]
        +conj(A[28])*B[24]
        +conj(A[34])*B[30];

  C[25] = conj(A[4])*B[1]
        +conj(A[10])*B[7]
        +conj(A[16])*B[13]
        +conj(A[22])*B[19]
        +conj(A[28])*B[25]
        +conj(A[34])*B[31];

  C[26] = conj(A[4])*B[2]
        +conj(A[10])*B[8]
        +conj(A[16])*B[14]
        +conj(A[22])*B[20]
        +conj(A[28])*B[26]
        +conj(A[34])*B[32];

  C[27] = conj(A[4])*B[3]
        +conj(A[10])*B[9]
        +conj(A[16])*B[15]
        +conj(A[22])*B[21]
        +conj(A[28])*B[27]
        +conj(A[34])*B[33];

  C[28] = conj(A[4])*B[4]
        +conj(A[10])*B[10]
        +conj(A[16])*B[16]
        +conj(A[22])*B[22]
        +conj(A[28])*B[28]
        +conj(A[34])*B[34];

  C[29] = conj(A[4])*B[5]
        +conj(A[10])*B[11]
        +conj(A[16])*B[17]
        +conj(A[22])*B[23]
        +conj(A[28])*B[29]
        +conj(A[34])*B[35];

  C[30] = conj(A[5])*B[0]
        +conj(A[11])*B[6]
        +conj(A[17])*B[12]
        +conj(A[23])*B[18]
        +conj(A[29])*B[24]
        +conj(A[35])*B[30];

  C[31] = conj(A[5])*B[1]
        +conj(A[11])*B[7]
        +conj(A[17])*B[13]
        +conj(A[23])*B[19]
        +conj(A[29])*B[25]
        +conj(A[35])*B[31];

  C[32] = conj(A[5])*B[2]
        +conj(A[11])*B[8]
        +conj(A[17])*B[14]
        +conj(A[23])*B[20]
        +conj(A[29])*B[26]
        +conj(A[35])*B[32];

  C[33] = conj(A[5])*B[3]
        +conj(A[11])*B[9]
        +conj(A[17])*B[15]
        +conj(A[23])*B[21]
        +conj(A[29])*B[27]
        +conj(A[35])*B[33];

  C[34] = conj(A[5])*B[4]
        +conj(A[11])*B[10]
        +conj(A[17])*B[16]
        +conj(A[23])*B[22]
        +conj(A[29])*B[28]
        +conj(A[35])*B[34];

  C[35] = conj(A[5])*B[5]
        +conj(A[11])*B[11]
        +conj(A[17])*B[17]
        +conj(A[23])*B[23]
        +conj(A[29])*B[29]
        +conj(A[35])*B[35];

};


