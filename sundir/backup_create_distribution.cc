#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdio.h>

using namespace std;

bool myfunction (double i, double j) { return (i<j); }

struct myclass {
  bool operator() (double i, double j) { return (i<j);}
} myobject;

char *input_file_name = new char[2000];
char Jarzynski_direction;
bool direct_transformation;
int number_of_data, number_of_columns, dummy_int;
double dummy_double;
double min_value, max_value, col_width;

int main(int argc, char **argv) {

  if (argc!=4) {
    cout << "Usage for scalar mode: ./create_distribution.x <input_file_name> <number_of_data> <type>" << endl; 
    exit(0);
  };

  sprintf(input_file_name,"%s",argv[1]);
  number_of_data=atoi(argv[2]);
  Jarzynski_direction=*(argv[3]);
  switch (Jarzynski_direction) {
    case 'd':
    case 'D': {
      direct_transformation=true;
      break;
    }
    case 'r':
    case 'R': {
      direct_transformation=false;
      break;
    }
    default: {
      cout << "Error while reading the input file: the direction must be 'd', 'D', 'r' or 'R'" << endl;
      exit(0);
      break;
    }
  }
  double *data = new double[number_of_data];
  number_of_columns=floor(sqrt(number_of_data));
// number_of_columns/=10;
  int *columns = new int[number_of_columns];

  ifstream inputfile;
  inputfile.open(input_file_name);
  if (!inputfile.is_open()) {
    cout << "Error! " << input_file_name << " does not exists" << endl;
    cout << "Execution aborted\n";
    exit(0);
  }

  inputfile >> dummy_int;
  if (dummy_int!=0) {
    cout << "Error! Unexpected value for the measurement index" << endl;
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
  for (int col_index=0; col_index<number_of_columns; col_index++) {
    cout << min_value << " " << dummy_double*columns[col_index] << endl;
    min_value+=col_width;
  }
// max_value=0.;
// for (int col_index=0; col_index<number_of_columns; col_index++) {
// max_value+=col_width*dummy_double*columns[col_index];
// }
// cout << "Final sum="<<max_value<< endl;

  return 0;

}
