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
  void cnvt2D() const;
  //  void mtxPrint() const;

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

// void mtxPrint() const{
//   for (int i = 0; i < dim; ++i){
//     for (int j = 0; j < dim; ++j){
      

int main(){
  int dim = 10;
  int nnz = 5;
  SprsMtx smtx1(dim,nnz);

  // example of indx and jndx
  int indx[5] = {1,2,3,4,5};
  int jndx[5] = {1,2,3,4,5};
  double mtx[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

  smtx1.init(indx, jndx, mtx,nnz);
  
  smtx1.sget();
  
  // double *mtx1 = new double [dim*dim];
  smtx1.cnvt2D();


return 0;
}
