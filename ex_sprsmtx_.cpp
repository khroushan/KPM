#include "sparse_matrix.hpp"
#include <iostream>
#include <iomanip>

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
