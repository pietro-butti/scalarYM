inline void hermitian_conjugate(dc udag[n][n], dc u[n][n]) {

  int i, j;

  for (i=0;i<n;i++)
  for (j=0;j<n;j++) {
    udag[i][j]=conjg(udag[j][i]);
  }

}
