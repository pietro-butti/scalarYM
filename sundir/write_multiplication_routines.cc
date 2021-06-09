#include <iostream>
#include <fstream>
#include <sstream>
#include <string> 
#include <cstdio>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

using namespace std;

void write_SUN_multiplication(int Ncol) {

  int i, j, k;
  
  printf("inline void mult_C_equals_AB_for_SU%d(dc *C, dc *A, dc *B) {\n\n", Ncol);
  for (i=0;i<Ncol;i++)
  for (j=0;j<Ncol;j++) {
    printf("  C[%d] = A[%d]*B[%d]",Ncol*i+j,Ncol*i,j);
    for (k=1;k<Ncol;k++) {
      printf("\n        +A[%d]*B[%d]",Ncol*i+k,Ncol*k+j);
    }
    printf(";\n\n");
  }
  printf("};\n\n\n");
 
  
  printf("inline void mult_C_equals_ABdagger_for_SU%d(dc *C, dc *A, dc *B) {\n\n", Ncol);
  for (i=0;i<Ncol;i++)
  for (j=0;j<Ncol;j++) {
    printf("  C[%d] = A[%d]*conj(B[%d])",Ncol*i+j,Ncol*i,Ncol*j);
    for (k=1;k<Ncol;k++) {
      printf("\n        +A[%d]*conj(B[%d])",Ncol*i+k,Ncol*j+k);
    }
    printf(";\n\n");
  }
  printf("};\n\n\n");
  
  printf("inline void mult_C_equals_AdaggerB_for_SU%d(dc *C, dc *A, dc *B) {\n\n", Ncol);
  for (i=0;i<Ncol;i++)
  for (j=0;j<Ncol;j++) {
    printf("  C[%d] = conj(A[%d])*B[%d]",Ncol*i+j,i,j);
    for (k=1;k<Ncol;k++) {
      printf("\n        +conj(A[%d])*B[%d]",Ncol*k+i,Ncol*k+j);
    }
    printf(";\n\n");
  }
  printf("};\n\n\n");
 
}

int main(int argc, char **argv) {

  int Ncol=atoi(argv[1]);
  write_SUN_multiplication(Ncol);

}
