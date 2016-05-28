// Class of sparse matrix. The matrix is specified
// with the following inputs: (limited to square matrix)
// dimension of matrix :                          dim
// nubmer of nonzero elements :                   nnZero
// array of row indecies of nonzero elements :    indx[nnZero]
// array of column indecies of nonzero elements : indx[nnZero]
// nonzero matrix elements mtx[indx,jndx]       : mtx[nnZero]
// ===========================================================
// Author : A. Ahmadi
// Date   : May 27 2016
//============================================================
#include <iostream>
using namespace std;

class SprsMtx
{
public:
  SprsMtx(int sdim, int snnz);         // constructor
  ~SprsMtx();                          // destructor
  // need a method to set the value of indx,jndx and mtx
  // return the value
  void init(int *in_vec, int *jn_vec, double *mtx, int snnz);
  // then define operator overloading for sparse mtrx-mtrx multiplication
  
  void sget() const;

private:
  int dim, nnZero;
  int *indx, *jndx;
  double *mtx_elmt;
};

SprsMtx::SprsMtx(int sdim, int snnz){
  if (snnz<=(sdim*sdim)){
    dim = sdim;
    nnZero = snnz;
    indx = new int [nnZero];
    jndx = new int [nnZero];
    mtx_elmt = new double [nnZero];
  }else{
    cout << "\nnumber of non-zero elemetns exceeded"
	 << " number matrix elements\n";
  }
}

SprsMtx::~SprsMtx(){
  delete[] indx;
  delete[] jndx;
  delete[] mtx_elmt;
  indx = 0; jndx = 0; mtx_elmt = 0;
}

void SprsMtx::init(int *in_vec, int *jn_vec, double *mtx, int snnz){
  if (snnz == nnZero){
    for(int i = 0; i < nnZero; i++){
      indx[i] = in_vec[i];
      jndx[i] = jn_vec[i];
      mtx_elmt[i] = mtx[i];
    }//endfor
  }
  else{
    cout << "Error : non-zero elements does not match "
	 << " the vectors' dimension \n";
  }//endelse
}

void SprsMtx:: sget() const {
  for(int i = 0; i < nnZero; i++){
    cout << indx[i] << "\t";}
  cout << endl;
  for(int i = 0; i < nnZero; i++){
    cout << jndx[i] << "\t";}
  cout << endl;
  for(int i = 0; i < nnZero; i++){
    cout << mtx_elmt[i] << "\t";}
  cout << endl;
}


int main(){
  
  SprsMtx mtx1(10,5);

  // example of indx and jndx
  int indx[5] = {1,2,3,4,5};
  int jndx[5] = {1,2,3,4,5};
  double mtx[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

  mtx1.init(indx, jndx, mtx,5);
  
  mtx1.sget();
  return 0;
}
