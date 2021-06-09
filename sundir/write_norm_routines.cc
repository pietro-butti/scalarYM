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

void write_SUN_norm(int Ncol) {

  cout << "#ifndef __norm_su5_h__" << endl;
  cout << "#define __norm_su5_h__" << endl << endl;
  cout << "void norm_su"<< Ncol <<"(dc *vsu"<< Ncol <<" ) {" << endl;
  cout << endl;
  cout << "// Takes the first "<<Ncol-1<<" rows of a "<<Ncol<<"x"<<Ncol<<" complex matrix," << endl;
  cout << "// and makes an SU("<<Ncol<<") matrix out of them." << endl;
  cout << endl;
  cout << "  dc determinant_conjugate_phase=dc(1.,0.);" << endl;
  cout << "  struct complex_double b[Ncol], DUMMY[1][1], WORK[2*Ncol];" << endl;
  cout << "  double AT[2*Ncolsquare];" << endl;
  cout << "  int i, j, ok, c1, c2, c3;" << endl;
  cout << "  char c4;" << endl;
  cout << endl;
  cout << "  dc scalprod=dc(0.,0.);" << endl;
  cout << "  double xn;" << endl;
  cout << endl;
  cout << "// Normalize the row `0' of vsu" << Ncol << endl << endl;
  cout << "  xn = 1./sqrt( norm(vsu" << Ncol <<"[0])";
  for (int i=1; i<Ncol; i++) {
    cout << "\n    + norm(vsu"<< Ncol <<"["<< i <<"])"; 
  }
  cout << " );\n\n";
  for (int i=0; i<Ncol; i++) {
    cout << "  vsu"<< Ncol <<"["<< i <<"]*=xn;\n";
  }
  
  for (int newrow=1; newrow<Ncol; newrow++) {
    for (int oldrow=0; oldrow<newrow; oldrow++) {
     cout << "\n// Make the row `"<< newrow <<"' orthogonal to the row `"<<oldrow<<"':\n\n";
     cout <<"  scalprod = conj( vsu"<< Ncol <<"["<<Ncol*oldrow<<"] )*vsu"<<Ncol<<"["<<Ncol*newrow<<"]";
     for (int i=1;i<Ncol;i++) { 
       cout << "\n    + conj( vsu"<< Ncol<<"["<<Ncol*oldrow+i<<"] )*vsu"<< Ncol<<"["<<Ncol*newrow+i<<"]";
     }
     cout <<";\n\n";
     for (int i=0;i<Ncol;i++) { 
       cout << "  vsu"<< Ncol <<"["<<Ncol*newrow+i<<"] -= scalprod*vsu"<< Ncol <<"["<<Ncol*oldrow+i<<"];\n";
     }
     cout <<"\n";

    }
    
    cout << endl;
    cout << "// Normalize the row `"<< newrow <<"'" << endl << endl;
    cout << "  xn = 1./sqrt( norm(vsu" << Ncol <<"["<< newrow*Ncol << "])";
    for (int i=1; i<Ncol; i++) {
      cout << "\n    + norm(vsu"<< Ncol <<"["<< newrow*Ncol+i <<"])"; 
    }
    cout << " );\n\n";
    for (int i=0; i<Ncol; i++) {
      cout << "    vsu"<< Ncol <<"["<< newrow*Ncol+i <<"]*=xn;\n";
    }
    cout << endl;
  }
  
  cout << "\n// Impose unimodularity by arranging the phase of the last row\n";
  cout << "// to compensate for the phase of the determinant\n";
  cout << "\n";
  cout << "  for (i=0; i<Ncol; i++)\n";
  cout << "  for (j=0; j<Ncol; j++) {\n";
  cout << "    AT[2*(j+Ncol*i)]=real(vsu"<< Ncol <<"[j*Ncol+i]);\n";
  cout << "    AT[2*(j+Ncol*i)+1]=imag(vsu"<< Ncol <<"[j*Ncol+i]);\n";
  cout << "  }\n";
  cout << "  c1=Ncol;\n";
  cout << "  c2=2*Ncol;\n";
  cout << "  c3=1;\n";
  cout << "  c4='N';\n\n";
  cout << "  zgeev_(&c4, &c4,&c1, AT, &c1, b, DUMMY, &c3, DUMMY, &c3, WORK, &c2, WORK, &ok);\n\n";
  cout << "  dc determinante=dc(1.,0.);\n";
  cout << "  if (ok==0) {\n";
  cout << "    for (i=0; i<Ncol; i++) {\n";
  cout << "      determinante*=dc(b[i].re, b[i].im);\n";
  cout << "    }\n";
  cout << "  }\n";
  cout << "  else { printf(\"An error occured\"); exit(0); }\n\n";
  cout << "  determinant_conjugate_phase=conj(determinante);\n\n";

  for (int i=0; i<Ncol; i++) {
    cout << "  vsu"<< Ncol <<"["<< (Ncol-1)*Ncol+i <<"] *= determinant_conjugate_phase;\n";
  }
  cout << "\n\n}\n\n#endif\n";



  
}

int main(int argc, char **argv) {

  int Ncol=atoi(argv[1]);
  write_SUN_norm(Ncol);

}
