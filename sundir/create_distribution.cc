#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iterator>
#include <stdio.h>

using namespace std;

bool myfunction (double i, double j) { return (i<j); }

struct myclass {
  bool operator() (double i, double j) { return (i<j);}
} myobject;

char *file_name = new char[2000];
// char Jarzynski_direction;
bool direct_transformation;
int number_of_data, number_of_columns, dummy_int;
double dummy_double;
double min_value, max_value, col_width;

int main(int argc, char **argv) {

//   if (argc!=4) {
//     cout << "Usage: ./create_distribution.x <file_name> <number_of_data> <type>" << endl; 
//     exit(0);
//   };
  
  if (argc!=2) {
    cout << "Usage: ./create_distribution.x <file_name>" << endl; 
    exit(0);
  };

  sprintf(file_name,"%s",argv[1]);
  
  if ((strstr(file_name, "direct")==NULL)&&(strstr(file_name, "reverse")==NULL)) {
    cout << "Error in the input file name: the direction (\"direct\" or \"reverse\") is not specified" << endl;
      exit(0);
  }
  else {
    if (strstr(file_name, "reverse")==NULL) {
      direct_transformation=true;        
    }
    else{
      direct_transformation=false;        
    }
  }
  
// if (direct_transformation==true) {
//   cout << "direct transformation" << endl;     
// }
// else {
//   cout << "reverse transformation" << endl;     
// }
  
//   number_of_data=atoi(argv[2]);
//   Jarzynski_direction=*(argv[3]);
//   switch (Jarzynski_direction) {
//     case 'd':
//     case 'D': {
//       direct_transformation=true;
//       break;
//     }
//     case 'r':
//     case 'R': {
//       direct_transformation=false;
//       break;
//     }
//     default: {
//       cout << "Error while reading the input file: the direction must be 'd', 'D', 'r' or 'R'" << endl;
//       exit(0);
//       break;
//     }
//   }

  ifstream inputfile;
  inputfile.open(file_name);
  if (!inputfile.is_open()) {
    cout << "Error! " << file_name << " does not exists" << endl;
    cout << "Execution aborted\n";
    exit(0);
  }
  
  // new lines will be skipped unless we stop it from happening:    
  inputfile.unsetf(ios_base::skipws);

  // count the newlines with an algorithm specialized for counting:
  number_of_data = count(
        istream_iterator<char>(inputfile),
        istream_iterator<char>(),
        '\n');

  inputfile.clear();
  inputfile.seekg(0, ios::beg);
  inputfile.setf(ios_base::skipws);
  
  double *data = new double[number_of_data];
  number_of_columns=floor(sqrt(number_of_data));
// number_of_columns/=10;
  int *columns = new int[number_of_columns];

  inputfile >> dummy_int;
  if (dummy_int!=0) {
    cout << "Error! Unexpected value for the measurement index here" << endl;
    cout << "Execution aborted\n";
    exit(0);
  }
  inputfile >> dummy_double;
  if (direct_transformation==false) {
    dummy_double=-dummy_double;
  }
  max_value=min_value=data[0]=dummy_double;
  inputfile >> dummy_double;

  for (int meas_index=1; meas_index<number_of_data; meas_index++) {
    inputfile >> dummy_int;
    if (dummy_int!=meas_index) {
      cout << "Error! Unexpected value for the measurement index" << endl;
      cout << dummy_int << endl;
      cout << "Execution aborted\n";
      exit(0);
    }
    inputfile >> dummy_double;
    if (direct_transformation==false) {
      dummy_double=-dummy_double;
    }
    data[meas_index]=dummy_double;
    if (dummy_double<min_value) {
      min_value=dummy_double;
    }
    if (dummy_double>max_value) {
      max_value=dummy_double;
    }
    inputfile >> dummy_double;
  }
// cout << "min_value=" << min_value << "  max_value=" << max_value << endl;
  col_width=(max_value-min_value)/(number_of_columns-1.);
  for (int col_index=0; col_index<number_of_columns; col_index++) {
    columns[col_index]=0.;
  }

  vector<double> myvector(data, data+number_of_data);               // 32 71 12 45 26 80 53 33

  // using default comparison (operator <):
//   sort(myvector.begin(), myvector.begin()+4);           //(12 32 45 71)26 80 53 33

  // using function as comp
//   sort(myvector.begin()+4, myvector.end(), myfunction); // 12 32 45 71(26 33 53 80)

  // using object as comp
  sort(myvector.begin(), myvector.end());     //(12 26 32 33 45 53 71 80)

  // print out content:
//   cout << "myvector contains:";

  dummy_double=min_value+0.5*col_width;
  dummy_int=0;
  for (vector<double>::iterator it=myvector.begin(); it!=myvector.end(); ++it) {
//     cout << ' ' << *it << endl;
    while ((*it)>dummy_double) {
      dummy_int++;
      dummy_double+=col_width;
    }
    columns[dummy_int]++;
  }
// cout << "dummy_int="<<dummy_int<<endl;exit(0);
  dummy_double=1./(number_of_data*col_width);
  
  char *pch;
  pch=strstr(file_name,"Jarzynski_SF");
  strncpy(pch,"distribution",12);
//   cout << "Creating ";
//   puts(file_name);
  
  
  ofstream outputfile;
  outputfile.open(file_name);
  
  for (int col_index=0; col_index<number_of_columns; col_index++) {
//     cout << min_value << " " << dummy_double*columns[col_index] << endl;
    outputfile << min_value << " " << dummy_double*columns[col_index] << endl;
    min_value+=col_width;
  }
// max_value=0.;
// for (int col_index=0; col_index<number_of_columns; col_index++) {
// max_value+=col_width*dummy_double*columns[col_index];
// }
// cout << "Final sum="<<max_value<< endl;

  outputfile.close();
  return 0;

}
