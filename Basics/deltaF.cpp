using namespace std;

#include <iostream>
#include <fstream> // io file lib
#include <iomanip> // text formatting
#include <cmath>
#include <stdio.h>
#include <string>
const double pi = 4.0 * atan(1.0);

// Jackson kernel to smooth the Cheb polynomials
// to avoid the Gibbs oscillation
inline double j_kernel(int Max, int i){
  double aux = ( (double) 1/(Max + 1) ) * ( ( Max - i + 1 ) *	\
	       cos( (double) i * pi / ( Max + 1 ) ) + \
	       sin( (double) i * pi / ( Max + 1 ) ) / tan( pi / (Max+1) ));
  return aux;}


int main(){
  

  const int N_num = 1000; // number of Chbychev terms
  double x, x_a, delta; // x , and a value to expand the delta-function
  double chx, cha[N_num];
  double x_i = -0.2, x_f = 0.2;
  int steps = 500;   // steps for x
  string junk;
  
  // Reading inital parameters
  ifstream in_file("input.dat", ios :: in);
  in_file >> x_i >> x_f >> junk;
  in_file >> steps >> junk;
  in_file.close();
  cout << "x_i = " << x_i << "   x_f =  " << x_f << endl;
  cout << "steps = " << steps << junk << endl;
  double dx = (x_f - x_i)/steps;

  ofstream out_file("k2_delta.dat", ios::out);
  // kind 2
  x_a = 0.0;
  for (int i = 0; i < N_num; i++){ // generating U_n(a)
    cha[i] = sin( (i+1) * acos(x_a) ) / sin( acos(x_a) );
  }
  
  for (int j = 0; j < steps; j++ ){
    x = x_i + j * dx;
    delta = 0.0;
    for (int i = 0; i < N_num; i++) {
      chx = sin( (i+1) * acos(x) ) / sin( acos(x) );
      delta += cha[i] * chx * j_kernel(N_num,i);
    }
    delta *= ((double) 2/pi);
    out_file <<  setw(8) << setprecision(4) <<  x << "   " ;
    out_file <<  setw(8) << fixed << setprecision(5) << delta << endl;
  }
  
  out_file.close();

  out_file.open("k1_delta.dat", ios :: out);
  // kind 1
  for (int i = 0; i < N_num; i++) { // generating T_n(x)
    cha[i] = cos( i * acos(x_a) ) ;
  }

  for (int j = 0; j < steps; j++){
    x = x_i + j * dx;
    delta = 0.0;
    for(int i = 1; i < N_num; i++){ // second kind m = 0, m != 0
      chx = cos( i * acos(x) );
      delta += cha[i] * chx * j_kernel(N_num,i);
    }
    delta = (2.0 * delta + 1.0 * j_kernel(N_num,0) ) / pi;
    out_file <<  setw(8) << setprecision(4) <<  x << "   ";
    out_file <<  setw(8) << fixed << setprecision(5) << delta << endl;
  }
  out_file.close();

  out_file.open("kR_delta.dat", ios :: out);
  // kind one recursively
  for (int i = 0; i < N_num; i++) { // generating T_n(x)
    cha[i] = cos( i * acos(x_a) ) ;
  }
  double chx_p, chx_n;
  for (int j = 0; j < steps; j++){
    x = x_i + j * dx;
    delta = 0.0; chx_p = 1.0; chx = x;
    for(int i = 1; i < N_num; i++){ // second kind m = 0, m != 0
      delta += cha[i] * chx * j_kernel(N_num,i);
      chx_n = 2 * x * chx - chx_p;
      chx_p = chx;
      chx = chx_n;
    }
    delta = (2.0 * delta + 1.0 * j_kernel(N_num,0) ) / pi;
    out_file <<  setw(8) << setprecision(4) <<  x << "   ";
    out_file <<  setw(8) << fixed << setprecision(5) << delta << endl;
  }
  out_file.close();
  return 0;
  
}
