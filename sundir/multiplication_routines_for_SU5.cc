inline void mult_C_equals_AB_for_SU5(dc *C, dc *A, dc *B) {

  C[0] = A[0]*B[0]
        +A[1]*B[5]
        +A[2]*B[10]
        +A[3]*B[15]
        +A[4]*B[20];

  C[1] = A[0]*B[1]
        +A[1]*B[6]
        +A[2]*B[11]
        +A[3]*B[16]
        +A[4]*B[21];

  C[2] = A[0]*B[2]
        +A[1]*B[7]
        +A[2]*B[12]
        +A[3]*B[17]
        +A[4]*B[22];

  C[3] = A[0]*B[3]
        +A[1]*B[8]
        +A[2]*B[13]
        +A[3]*B[18]
        +A[4]*B[23];

  C[4] = A[0]*B[4]
        +A[1]*B[9]
        +A[2]*B[14]
        +A[3]*B[19]
        +A[4]*B[24];

  C[5] = A[5]*B[0]
        +A[6]*B[5]
        +A[7]*B[10]
        +A[8]*B[15]
        +A[9]*B[20];

  C[6] = A[5]*B[1]
        +A[6]*B[6]
        +A[7]*B[11]
        +A[8]*B[16]
        +A[9]*B[21];

  C[7] = A[5]*B[2]
        +A[6]*B[7]
        +A[7]*B[12]
        +A[8]*B[17]
        +A[9]*B[22];

  C[8] = A[5]*B[3]
        +A[6]*B[8]
        +A[7]*B[13]
        +A[8]*B[18]
        +A[9]*B[23];

  C[9] = A[5]*B[4]
        +A[6]*B[9]
        +A[7]*B[14]
        +A[8]*B[19]
        +A[9]*B[24];

  C[10] = A[10]*B[0]
        +A[11]*B[5]
        +A[12]*B[10]
        +A[13]*B[15]
        +A[14]*B[20];

  C[11] = A[10]*B[1]
        +A[11]*B[6]
        +A[12]*B[11]
        +A[13]*B[16]
        +A[14]*B[21];

  C[12] = A[10]*B[2]
        +A[11]*B[7]
        +A[12]*B[12]
        +A[13]*B[17]
        +A[14]*B[22];

  C[13] = A[10]*B[3]
        +A[11]*B[8]
        +A[12]*B[13]
        +A[13]*B[18]
        +A[14]*B[23];

  C[14] = A[10]*B[4]
        +A[11]*B[9]
        +A[12]*B[14]
        +A[13]*B[19]
        +A[14]*B[24];

  C[15] = A[15]*B[0]
        +A[16]*B[5]
        +A[17]*B[10]
        +A[18]*B[15]
        +A[19]*B[20];

  C[16] = A[15]*B[1]
        +A[16]*B[6]
        +A[17]*B[11]
        +A[18]*B[16]
        +A[19]*B[21];

  C[17] = A[15]*B[2]
        +A[16]*B[7]
        +A[17]*B[12]
        +A[18]*B[17]
        +A[19]*B[22];

  C[18] = A[15]*B[3]
        +A[16]*B[8]
        +A[17]*B[13]
        +A[18]*B[18]
        +A[19]*B[23];

  C[19] = A[15]*B[4]
        +A[16]*B[9]
        +A[17]*B[14]
        +A[18]*B[19]
        +A[19]*B[24];

  C[20] = A[20]*B[0]
        +A[21]*B[5]
        +A[22]*B[10]
        +A[23]*B[15]
        +A[24]*B[20];

  C[21] = A[20]*B[1]
        +A[21]*B[6]
        +A[22]*B[11]
        +A[23]*B[16]
        +A[24]*B[21];

  C[22] = A[20]*B[2]
        +A[21]*B[7]
        +A[22]*B[12]
        +A[23]*B[17]
        +A[24]*B[22];

  C[23] = A[20]*B[3]
        +A[21]*B[8]
        +A[22]*B[13]
        +A[23]*B[18]
        +A[24]*B[23];

  C[24] = A[20]*B[4]
        +A[21]*B[9]
        +A[22]*B[14]
        +A[23]*B[19]
        +A[24]*B[24];

};


inline void mult_C_equals_ABdagger_for_SU5(dc *C, dc *A, dc *B) {

  C[0] = A[0]*conj(B[0])
        +A[1]*conj(B[1])
        +A[2]*conj(B[2])
        +A[3]*conj(B[3])
        +A[4]*conj(B[4]);

  C[1] = A[0]*conj(B[5])
        +A[1]*conj(B[6])
        +A[2]*conj(B[7])
        +A[3]*conj(B[8])
        +A[4]*conj(B[9]);

  C[2] = A[0]*conj(B[10])
        +A[1]*conj(B[11])
        +A[2]*conj(B[12])
        +A[3]*conj(B[13])
        +A[4]*conj(B[14]);

  C[3] = A[0]*conj(B[15])
        +A[1]*conj(B[16])
        +A[2]*conj(B[17])
        +A[3]*conj(B[18])
        +A[4]*conj(B[19]);

  C[4] = A[0]*conj(B[20])
        +A[1]*conj(B[21])
        +A[2]*conj(B[22])
        +A[3]*conj(B[23])
        +A[4]*conj(B[24]);

  C[5] = A[5]*conj(B[0])
        +A[6]*conj(B[1])
        +A[7]*conj(B[2])
        +A[8]*conj(B[3])
        +A[9]*conj(B[4]);

  C[6] = A[5]*conj(B[5])
        +A[6]*conj(B[6])
        +A[7]*conj(B[7])
        +A[8]*conj(B[8])
        +A[9]*conj(B[9]);

  C[7] = A[5]*conj(B[10])
        +A[6]*conj(B[11])
        +A[7]*conj(B[12])
        +A[8]*conj(B[13])
        +A[9]*conj(B[14]);

  C[8] = A[5]*conj(B[15])
        +A[6]*conj(B[16])
        +A[7]*conj(B[17])
        +A[8]*conj(B[18])
        +A[9]*conj(B[19]);

  C[9] = A[5]*conj(B[20])
        +A[6]*conj(B[21])
        +A[7]*conj(B[22])
        +A[8]*conj(B[23])
        +A[9]*conj(B[24]);

  C[10] = A[10]*conj(B[0])
        +A[11]*conj(B[1])
        +A[12]*conj(B[2])
        +A[13]*conj(B[3])
        +A[14]*conj(B[4]);

  C[11] = A[10]*conj(B[5])
        +A[11]*conj(B[6])
        +A[12]*conj(B[7])
        +A[13]*conj(B[8])
        +A[14]*conj(B[9]);

  C[12] = A[10]*conj(B[10])
        +A[11]*conj(B[11])
        +A[12]*conj(B[12])
        +A[13]*conj(B[13])
        +A[14]*conj(B[14]);

  C[13] = A[10]*conj(B[15])
        +A[11]*conj(B[16])
        +A[12]*conj(B[17])
        +A[13]*conj(B[18])
        +A[14]*conj(B[19]);

  C[14] = A[10]*conj(B[20])
        +A[11]*conj(B[21])
        +A[12]*conj(B[22])
        +A[13]*conj(B[23])
        +A[14]*conj(B[24]);

  C[15] = A[15]*conj(B[0])
        +A[16]*conj(B[1])
        +A[17]*conj(B[2])
        +A[18]*conj(B[3])
        +A[19]*conj(B[4]);

  C[16] = A[15]*conj(B[5])
        +A[16]*conj(B[6])
        +A[17]*conj(B[7])
        +A[18]*conj(B[8])
        +A[19]*conj(B[9]);

  C[17] = A[15]*conj(B[10])
        +A[16]*conj(B[11])
        +A[17]*conj(B[12])
        +A[18]*conj(B[13])
        +A[19]*conj(B[14]);

  C[18] = A[15]*conj(B[15])
        +A[16]*conj(B[16])
        +A[17]*conj(B[17])
        +A[18]*conj(B[18])
        +A[19]*conj(B[19]);

  C[19] = A[15]*conj(B[20])
        +A[16]*conj(B[21])
        +A[17]*conj(B[22])
        +A[18]*conj(B[23])
        +A[19]*conj(B[24]);

  C[20] = A[20]*conj(B[0])
        +A[21]*conj(B[1])
        +A[22]*conj(B[2])
        +A[23]*conj(B[3])
        +A[24]*conj(B[4]);

  C[21] = A[20]*conj(B[5])
        +A[21]*conj(B[6])
        +A[22]*conj(B[7])
        +A[23]*conj(B[8])
        +A[24]*conj(B[9]);

  C[22] = A[20]*conj(B[10])
        +A[21]*conj(B[11])
        +A[22]*conj(B[12])
        +A[23]*conj(B[13])
        +A[24]*conj(B[14]);

  C[23] = A[20]*conj(B[15])
        +A[21]*conj(B[16])
        +A[22]*conj(B[17])
        +A[23]*conj(B[18])
        +A[24]*conj(B[19]);

  C[24] = A[20]*conj(B[20])
        +A[21]*conj(B[21])
        +A[22]*conj(B[22])
        +A[23]*conj(B[23])
        +A[24]*conj(B[24]);

};


inline void mult_C_equals_AdaggerB_for_SU5(dc *C, dc *A, dc *B) {

  C[0] = conj(A[0])*B[0]
        +conj(A[5])*B[5]
        +conj(A[10])*B[10]
        +conj(A[15])*B[15]
        +conj(A[20])*B[20];

  C[1] = conj(A[0])*B[1]
        +conj(A[5])*B[6]
        +conj(A[10])*B[11]
        +conj(A[15])*B[16]
        +conj(A[20])*B[21];

  C[2] = conj(A[0])*B[2]
        +conj(A[5])*B[7]
        +conj(A[10])*B[12]
        +conj(A[15])*B[17]
        +conj(A[20])*B[22];

  C[3] = conj(A[0])*B[3]
        +conj(A[5])*B[8]
        +conj(A[10])*B[13]
        +conj(A[15])*B[18]
        +conj(A[20])*B[23];

  C[4] = conj(A[0])*B[4]
        +conj(A[5])*B[9]
        +conj(A[10])*B[14]
        +conj(A[15])*B[19]
        +conj(A[20])*B[24];

  C[5] = conj(A[1])*B[0]
        +conj(A[6])*B[5]
        +conj(A[11])*B[10]
        +conj(A[16])*B[15]
        +conj(A[21])*B[20];

  C[6] = conj(A[1])*B[1]
        +conj(A[6])*B[6]
        +conj(A[11])*B[11]
        +conj(A[16])*B[16]
        +conj(A[21])*B[21];

  C[7] = conj(A[1])*B[2]
        +conj(A[6])*B[7]
        +conj(A[11])*B[12]
        +conj(A[16])*B[17]
        +conj(A[21])*B[22];

  C[8] = conj(A[1])*B[3]
        +conj(A[6])*B[8]
        +conj(A[11])*B[13]
        +conj(A[16])*B[18]
        +conj(A[21])*B[23];

  C[9] = conj(A[1])*B[4]
        +conj(A[6])*B[9]
        +conj(A[11])*B[14]
        +conj(A[16])*B[19]
        +conj(A[21])*B[24];

  C[10] = conj(A[2])*B[0]
        +conj(A[7])*B[5]
        +conj(A[12])*B[10]
        +conj(A[17])*B[15]
        +conj(A[22])*B[20];

  C[11] = conj(A[2])*B[1]
        +conj(A[7])*B[6]
        +conj(A[12])*B[11]
        +conj(A[17])*B[16]
        +conj(A[22])*B[21];

  C[12] = conj(A[2])*B[2]
        +conj(A[7])*B[7]
        +conj(A[12])*B[12]
        +conj(A[17])*B[17]
        +conj(A[22])*B[22];

  C[13] = conj(A[2])*B[3]
        +conj(A[7])*B[8]
        +conj(A[12])*B[13]
        +conj(A[17])*B[18]
        +conj(A[22])*B[23];

  C[14] = conj(A[2])*B[4]
        +conj(A[7])*B[9]
        +conj(A[12])*B[14]
        +conj(A[17])*B[19]
        +conj(A[22])*B[24];

  C[15] = conj(A[3])*B[0]
        +conj(A[8])*B[5]
        +conj(A[13])*B[10]
        +conj(A[18])*B[15]
        +conj(A[23])*B[20];

  C[16] = conj(A[3])*B[1]
        +conj(A[8])*B[6]
        +conj(A[13])*B[11]
        +conj(A[18])*B[16]
        +conj(A[23])*B[21];

  C[17] = conj(A[3])*B[2]
        +conj(A[8])*B[7]
        +conj(A[13])*B[12]
        +conj(A[18])*B[17]
        +conj(A[23])*B[22];

  C[18] = conj(A[3])*B[3]
        +conj(A[8])*B[8]
        +conj(A[13])*B[13]
        +conj(A[18])*B[18]
        +conj(A[23])*B[23];

  C[19] = conj(A[3])*B[4]
        +conj(A[8])*B[9]
        +conj(A[13])*B[14]
        +conj(A[18])*B[19]
        +conj(A[23])*B[24];

  C[20] = conj(A[4])*B[0]
        +conj(A[9])*B[5]
        +conj(A[14])*B[10]
        +conj(A[19])*B[15]
        +conj(A[24])*B[20];

  C[21] = conj(A[4])*B[1]
        +conj(A[9])*B[6]
        +conj(A[14])*B[11]
        +conj(A[19])*B[16]
        +conj(A[24])*B[21];

  C[22] = conj(A[4])*B[2]
        +conj(A[9])*B[7]
        +conj(A[14])*B[12]
        +conj(A[19])*B[17]
        +conj(A[24])*B[22];

  C[23] = conj(A[4])*B[3]
        +conj(A[9])*B[8]
        +conj(A[14])*B[13]
        +conj(A[19])*B[18]
        +conj(A[24])*B[23];

  C[24] = conj(A[4])*B[4]
        +conj(A[9])*B[9]
        +conj(A[14])*B[14]
        +conj(A[19])*B[19]
        +conj(A[24])*B[24];

};


