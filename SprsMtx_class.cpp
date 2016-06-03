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
#include <iomanip>
#include "sparse_matrix.hpp"

using namespace std;

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

void SprsMtx::cnvt2D()const{
  // allocate a 2D vector with dim x dim dimension
  // initialize it with zero
  double * mtx2D;
  mtx2D = new double [dim*dim];
  for (int i = 0; i < dim; ++i){
    for (int j = 0; j < dim; ++j){
      mtx2D[i*dim +j] = 0.d;
    }
  }
  // fill with non-zero elemetns
  for (int i = 0; i < nnZero; ++i){
    mtx2D[indx[i]*dim + jndx[i]] = mtx_elmt[i];
  }

  for (int i = 0; i < dim; ++i){
    for (int j = 0; j < dim; ++j){
      cout << setw(3) << mtx2D[i*dim +j] ;
    }
    cout << endl;
  }

  delete [] mtx2D;
}
