inline void mult_C_equals_AB_for_SU2(dc *C, dc *A, dc *B) {

  C[0] = A[0]*B[0]
        +A[1]*B[2];

  C[1] = A[0]*B[1]
        +A[1]*B[3];

  C[2] = A[2]*B[0]
        +A[3]*B[2];

  C[3] = A[2]*B[1]
        +A[3]*B[3];

};


inline void mult_C_equals_ABdagger_for_SU2(dc *C, dc *A, dc *B) {

  C[0] = A[0]*conj(B[0])
        +A[1]*conj(B[1]);

  C[1] = A[0]*conj(B[2])
        +A[1]*conj(B[3]);

  C[2] = A[2]*conj(B[0])
        +A[3]*conj(B[1]);

  C[3] = A[2]*conj(B[2])
        +A[3]*conj(B[3]);

};


inline void mult_C_equals_AdaggerB_for_SU2(dc *C, dc *A, dc *B) {

  C[0] = conj(A[0])*B[0]
        +conj(A[2])*B[2];

  C[1] = conj(A[0])*B[1]
        +conj(A[2])*B[3];

  C[2] = conj(A[1])*B[0]
        +conj(A[3])*B[2];

  C[3] = conj(A[1])*B[1]
        +conj(A[3])*B[3];

};


