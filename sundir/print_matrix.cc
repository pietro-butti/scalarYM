void print_matrix(dc tempor[Ncol][Ncol]) {

  int i, j;

  cout << "{";
  for (i=0;i<Ncol;i++) {
    cout << " {";
    for (j=0;j<Ncol;j++) {
      cout << real(tempor[i][j]) << " +I*(" <<imag(tempor[i][j]) << ")";
      if (j==Ncol-1) cout<<"} ";
      else cout <<", ";
    }
    if (i==Ncol-1) cout<<"} ";
    else cout <<", ";
  }
  cout <<endl;

}
