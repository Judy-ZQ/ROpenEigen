#include "ROpenEigen_types.h"


template<int SizeAtCompileTime, int MaxSizeAtCompileTime>
SEXP eigenPermutationMatrix_toDenseMatrix(Eigen::PermutationMatrix<SizeAtCompileTime, MaxSizeAtCompileTime>* obj){
  try{
    Eigen::Matrix<int,Eigen::Dynamic, Eigen::Dynamic> res = obj->toDenseMatrix();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

RCPP_MODULE(EigenPermutationMatrixModule){
  class_<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>>("EigenPermutationMatrixd")
  // .default_constructor()
  .constructor<Eigen::Transpositions<Eigen::Dynamic, Eigen::Dynamic>>()
  // .constructor<int>()
  .method("toDenseMatrix",&eigenPermutationMatrix_toDenseMatrix<Eigen::Dynamic, Eigen::Dynamic>)
  ;

}
