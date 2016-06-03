#ifndef SPRSMTXHEAD
#define SPRSMTXHEAD

#include <iostream>
#include <iomanip>

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

#endif
