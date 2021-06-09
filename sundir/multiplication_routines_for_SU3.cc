inline void mult_C_equals_AB_for_SU3(dc *C, dc *A, dc *B) {

  C[0] = A[0]*B[0]
        +A[1]*B[3]
        +A[2]*B[6];

  C[1] = A[0]*B[1]
        +A[1]*B[4]
        +A[2]*B[7];

  C[2] = A[0]*B[2]
        +A[1]*B[5]
        +A[2]*B[8];

  C[3] = A[3]*B[0]
        +A[4]*B[3]
        +A[5]*B[6];

  C[4] = A[3]*B[1]
        +A[4]*B[4]
        +A[5]*B[7];

  C[5] = A[3]*B[2]
        +A[4]*B[5]
        +A[5]*B[8];

  C[6] = A[6]*B[0]
        +A[7]*B[3]
        +A[8]*B[6];

  C[7] = A[6]*B[1]
        +A[7]*B[4]
        +A[8]*B[7];

  C[8] = A[6]*B[2]
        +A[7]*B[5]
        +A[8]*B[8];

};


inline void mult_C_equals_ABdagger_for_SU3(dc *C, dc *A, dc *B) {

  C[0] = A[0]*conj(B[0])
        +A[1]*conj(B[1])
        +A[2]*conj(B[2]);

  C[1] = A[0]*conj(B[3])
        +A[1]*conj(B[4])
        +A[2]*conj(B[5]);

  C[2] = A[0]*conj(B[6])
        +A[1]*conj(B[7])
        +A[2]*conj(B[8]);

  C[3] = A[3]*conj(B[0])
        +A[4]*conj(B[1])
        +A[5]*conj(B[2]);

  C[4] = A[3]*conj(B[3])
        +A[4]*conj(B[4])
        +A[5]*conj(B[5]);

  C[5] = A[3]*conj(B[6])
        +A[4]*conj(B[7])
        +A[5]*conj(B[8]);

  C[6] = A[6]*conj(B[0])
        +A[7]*conj(B[1])
        +A[8]*conj(B[2]);

  C[7] = A[6]*conj(B[3])
        +A[7]*conj(B[4])
        +A[8]*conj(B[5]);

  C[8] = A[6]*conj(B[6])
        +A[7]*conj(B[7])
        +A[8]*conj(B[8]);

};


inline void mult_C_equals_AdaggerB_for_SU3(dc *C, dc *A, dc *B) {

  C[0] = conj(A[0])*B[0]
        +conj(A[3])*B[3]
        +conj(A[6])*B[6];

  C[1] = conj(A[0])*B[1]
        +conj(A[3])*B[4]
        +conj(A[6])*B[7];

  C[2] = conj(A[0])*B[2]
        +conj(A[3])*B[5]
        +conj(A[6])*B[8];

  C[3] = conj(A[1])*B[0]
        +conj(A[4])*B[3]
        +conj(A[7])*B[6];

  C[4] = conj(A[1])*B[1]
        +conj(A[4])*B[4]
        +conj(A[7])*B[7];

  C[5] = conj(A[1])*B[2]
        +conj(A[4])*B[5]
        +conj(A[7])*B[8];

  C[6] = conj(A[2])*B[0]
        +conj(A[5])*B[3]
        +conj(A[8])*B[6];

  C[7] = conj(A[2])*B[1]
        +conj(A[5])*B[4]
        +conj(A[8])*B[7];

  C[8] = conj(A[2])*B[2]
        +conj(A[5])*B[5]
        +conj(A[8])*B[8];

};


