void print_2x2_matrix(dc tempor[2][2]) {

  int i, j;

  cout << "{";
  for (i=0;i<2;i++) {
    cout << " {";
    for (j=0;j<2;j++) {
      cout << real(tempor[i][j]) << " +I*(" <<imag(tempor[i][j]) << ")";
      if (j==1) cout<<"} ";
      else cout <<", ";
    }
    if (i==1) cout<<"} ";
    else cout <<", ";
  }
  cout <<endl;

}
