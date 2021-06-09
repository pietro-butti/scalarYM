void print_diagonal_matrix(dc tempor[Ncol]) {

  int i, j;
  double re, im;

  cout << "{";
  for (i=0;i<Ncol;i++) {
    cout << " {";
    for (j=0;j<Ncol;j++) {
      re=0.0;
      im=0.0;
      if (i==j) {
        re=real(tempor[i]);
        im=imag(tempor[i]);
      }
      cout << re << " +I*(" << im << ")";
      if (j==Ncol-1) cout<<"} ";
      else cout <<", ";
    }
    if (i==Ncol-1) cout<<"} ";
    else cout <<", ";
  }
  cout <<endl;

}

void print_diagonal_matrix(double tempor[Ncol]) {

  int i,j;
  double re;

  cout << "{";
  for (i=0;i<Ncol;i++) {
    cout << " {";
    for (j=0;j<Ncol;j++) {
      re=0.0;
      if (i==j) {
        re=tempor[i];
      }
      cout << re;
      if (j==Ncol-1) cout<<"} ";
      else cout <<", ";
    }
    if (i==Ncol-1) cout<<"} ";
    else cout <<", ";
  }
  cout <<endl;

}
