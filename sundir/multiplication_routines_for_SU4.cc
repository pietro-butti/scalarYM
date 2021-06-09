inline void mult_C_equals_AB_for_SU4(dc *C, dc *A, dc *B) {

  C[0] = A[0]*B[0]
        +A[1]*B[4]
        +A[2]*B[8]
        +A[3]*B[12];

  C[1] = A[0]*B[1]
        +A[1]*B[5]
        +A[2]*B[9]
        +A[3]*B[13];

  C[2] = A[0]*B[2]
        +A[1]*B[6]
        +A[2]*B[10]
        +A[3]*B[14];

  C[3] = A[0]*B[3]
        +A[1]*B[7]
        +A[2]*B[11]
        +A[3]*B[15];

  C[4] = A[4]*B[0]
        +A[5]*B[4]
        +A[6]*B[8]
        +A[7]*B[12];

  C[5] = A[4]*B[1]
        +A[5]*B[5]
        +A[6]*B[9]
        +A[7]*B[13];

  C[6] = A[4]*B[2]
        +A[5]*B[6]
        +A[6]*B[10]
        +A[7]*B[14];

  C[7] = A[4]*B[3]
        +A[5]*B[7]
        +A[6]*B[11]
        +A[7]*B[15];

  C[8] = A[8]*B[0]
        +A[9]*B[4]
        +A[10]*B[8]
        +A[11]*B[12];

  C[9] = A[8]*B[1]
        +A[9]*B[5]
        +A[10]*B[9]
        +A[11]*B[13];

  C[10] = A[8]*B[2]
        +A[9]*B[6]
        +A[10]*B[10]
        +A[11]*B[14];

  C[11] = A[8]*B[3]
        +A[9]*B[7]
        +A[10]*B[11]
        +A[11]*B[15];

  C[12] = A[12]*B[0]
        +A[13]*B[4]
        +A[14]*B[8]
        +A[15]*B[12];

  C[13] = A[12]*B[1]
        +A[13]*B[5]
        +A[14]*B[9]
        +A[15]*B[13];

  C[14] = A[12]*B[2]
        +A[13]*B[6]
        +A[14]*B[10]
        +A[15]*B[14];

  C[15] = A[12]*B[3]
        +A[13]*B[7]
        +A[14]*B[11]
        +A[15]*B[15];

};


inline void mult_C_equals_ABdagger_for_SU4(dc *C, dc *A, dc *B) {

  C[0] = A[0]*conj(B[0])
        +A[1]*conj(B[1])
        +A[2]*conj(B[2])
        +A[3]*conj(B[3]);

  C[1] = A[0]*conj(B[4])
        +A[1]*conj(B[5])
        +A[2]*conj(B[6])
        +A[3]*conj(B[7]);

  C[2] = A[0]*conj(B[8])
        +A[1]*conj(B[9])
        +A[2]*conj(B[10])
        +A[3]*conj(B[11]);

  C[3] = A[0]*conj(B[12])
        +A[1]*conj(B[13])
        +A[2]*conj(B[14])
        +A[3]*conj(B[15]);

  C[4] = A[4]*conj(B[0])
        +A[5]*conj(B[1])
        +A[6]*conj(B[2])
        +A[7]*conj(B[3]);

  C[5] = A[4]*conj(B[4])
        +A[5]*conj(B[5])
        +A[6]*conj(B[6])
        +A[7]*conj(B[7]);

  C[6] = A[4]*conj(B[8])
        +A[5]*conj(B[9])
        +A[6]*conj(B[10])
        +A[7]*conj(B[11]);

  C[7] = A[4]*conj(B[12])
        +A[5]*conj(B[13])
        +A[6]*conj(B[14])
        +A[7]*conj(B[15]);

  C[8] = A[8]*conj(B[0])
        +A[9]*conj(B[1])
        +A[10]*conj(B[2])
        +A[11]*conj(B[3]);

  C[9] = A[8]*conj(B[4])
        +A[9]*conj(B[5])
        +A[10]*conj(B[6])
        +A[11]*conj(B[7]);

  C[10] = A[8]*conj(B[8])
        +A[9]*conj(B[9])
        +A[10]*conj(B[10])
        +A[11]*conj(B[11]);

  C[11] = A[8]*conj(B[12])
        +A[9]*conj(B[13])
        +A[10]*conj(B[14])
        +A[11]*conj(B[15]);

  C[12] = A[12]*conj(B[0])
        +A[13]*conj(B[1])
        +A[14]*conj(B[2])
        +A[15]*conj(B[3]);

  C[13] = A[12]*conj(B[4])
        +A[13]*conj(B[5])
        +A[14]*conj(B[6])
        +A[15]*conj(B[7]);

  C[14] = A[12]*conj(B[8])
        +A[13]*conj(B[9])
        +A[14]*conj(B[10])
        +A[15]*conj(B[11]);

  C[15] = A[12]*conj(B[12])
        +A[13]*conj(B[13])
        +A[14]*conj(B[14])
        +A[15]*conj(B[15]);

};


inline void mult_C_equals_AdaggerB_for_SU4(dc *C, dc *A, dc *B) {

  C[0] = conj(A[0])*B[0]
        +conj(A[4])*B[4]
        +conj(A[8])*B[8]
        +conj(A[12])*B[12];

  C[1] = conj(A[0])*B[1]
        +conj(A[4])*B[5]
        +conj(A[8])*B[9]
        +conj(A[12])*B[13];

  C[2] = conj(A[0])*B[2]
        +conj(A[4])*B[6]
        +conj(A[8])*B[10]
        +conj(A[12])*B[14];

  C[3] = conj(A[0])*B[3]
        +conj(A[4])*B[7]
        +conj(A[8])*B[11]
        +conj(A[12])*B[15];

  C[4] = conj(A[1])*B[0]
        +conj(A[5])*B[4]
        +conj(A[9])*B[8]
        +conj(A[13])*B[12];

  C[5] = conj(A[1])*B[1]
        +conj(A[5])*B[5]
        +conj(A[9])*B[9]
        +conj(A[13])*B[13];

  C[6] = conj(A[1])*B[2]
        +conj(A[5])*B[6]
        +conj(A[9])*B[10]
        +conj(A[13])*B[14];

  C[7] = conj(A[1])*B[3]
        +conj(A[5])*B[7]
        +conj(A[9])*B[11]
        +conj(A[13])*B[15];

  C[8] = conj(A[2])*B[0]
        +conj(A[6])*B[4]
        +conj(A[10])*B[8]
        +conj(A[14])*B[12];

  C[9] = conj(A[2])*B[1]
        +conj(A[6])*B[5]
        +conj(A[10])*B[9]
        +conj(A[14])*B[13];

  C[10] = conj(A[2])*B[2]
        +conj(A[6])*B[6]
        +conj(A[10])*B[10]
        +conj(A[14])*B[14];

  C[11] = conj(A[2])*B[3]
        +conj(A[6])*B[7]
        +conj(A[10])*B[11]
        +conj(A[14])*B[15];

  C[12] = conj(A[3])*B[0]
        +conj(A[7])*B[4]
        +conj(A[11])*B[8]
        +conj(A[15])*B[12];

  C[13] = conj(A[3])*B[1]
        +conj(A[7])*B[5]
        +conj(A[11])*B[9]
        +conj(A[15])*B[13];

  C[14] = conj(A[3])*B[2]
        +conj(A[7])*B[6]
        +conj(A[11])*B[10]
        +conj(A[15])*B[14];

  C[15] = conj(A[3])*B[3]
        +conj(A[7])*B[7]
        +conj(A[11])*B[11]
        +conj(A[15])*B[15];

};


