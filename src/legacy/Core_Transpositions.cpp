#include "ROpenEigen_types.h"


template<int SizeAtCompileTime, int MaxSizeAtCompileTime>
void eigenTranspositions_OutToConsole(Eigen::Transpositions<SizeAtCompileTime, MaxSizeAtCompileTime>* obj){
  Rcout << "Transpositions (indices of permuted rows):" << endl;
  for (int i = 0; i < obj->size(); ++i) {
    Rcout << (*obj).indices()(i) << " ";
  }
  Rcout << endl;
}

RCPP_MODULE(EigenTranspositionsModule){
  class_<Eigen::Transpositions<Eigen::Dynamic, Eigen::Dynamic>>("EigenTranspositionsd")
  // .default_constructor()
  // .constructor<myMatrixXd>()
  // .constructor<int>()
  .method("OutToConsole",&eigenTranspositions_OutToConsole<Eigen::Dynamic, Eigen::Dynamic>)
  ;

}
