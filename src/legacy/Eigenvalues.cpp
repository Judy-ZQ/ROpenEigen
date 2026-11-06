#include "ROpenEigen_types.h"

template<typename _MatrixType>
SEXP eigenSelfAdjointEigenSolver_eigenvalues(Eigen::SelfAdjointEigenSolver<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> result = obj->eigenvalues();
    return Rcpp::wrap(result);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenSelfAdjointEigenSolver_eigenvectors(Eigen::SelfAdjointEigenSolver<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> result = obj->eigenvectors();
    return Rcpp::wrap(result);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}


template<typename _MatrixType>
SEXP eigenEigenSolver_eigenvalues(Eigen::EigenSolver<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> result = obj->eigenvalues();
    return Rcpp::wrap(result);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenEigenSolver_eigenvectors(Eigen::EigenSolver<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> result = obj->eigenvectors();
    return Rcpp::wrap(result);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
_MatrixType eigenEigenSolver_pseudoEigenvalueMatrix(Eigen::EigenSolver<_MatrixType>* obj){
  try{
    _MatrixType result = obj->pseudoEigenvalueMatrix();
    return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
_MatrixType eigenEigenSolver_pseudoEigenvectors(Eigen::EigenSolver<_MatrixType>* obj){
  try{
    _MatrixType result = obj->pseudoEigenvectors();
    return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenEigenSolver_getMaxIterations(Eigen::EigenSolver<_MatrixType>* obj){
  try{
    return obj->getMaxIterations();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
void eigenEigenSolver_setMaxIterations(Eigen::EigenSolver<_MatrixType>* obj,int maxIters){
  try{
    obj->setMaxIterations(maxIters);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

typedef Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> mySelfAdjointEigenSolver4MatrixXd;
typedef Eigen::SelfAdjointEigenSolver<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> mySelfAdjointEigenSolver4MatrixXcd;
typedef Eigen::EigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myEigenSolver4MatrixXd;



RCPP_MODULE(EigenEigenvaluesModule){
  class_<mySelfAdjointEigenSolver4MatrixXd>("SelfAdjointEigenSolver4MatrixXd")
  .default_constructor()
  .constructor<myMatrixXd>()
  .constructor<int>()
  .method("eigenvalues",&eigenSelfAdjointEigenSolver_eigenvalues<myMatrixXd>)
  .method("eigenvectors",&eigenSelfAdjointEigenSolver_eigenvectors<myMatrixXd>)
  ;

  class_<myEigenSolver4MatrixXd>("EigenSolver4MatrixXd")
    .default_constructor()
    .constructor<myMatrixXd,bool>()
    .constructor<int>()
    .method("eigenvalues",&eigenEigenSolver_eigenvalues<myMatrixXd>)
    .method("eigenvectors",&eigenEigenSolver_eigenvectors<myMatrixXd>)
    .method("getMaxIterations",&eigenEigenSolver_getMaxIterations<myMatrixXd>)
    .method("pseudoEigenvalueMatrix",&eigenEigenSolver_pseudoEigenvalueMatrix<myMatrixXd>)
    .method("pseudoEigenvectors",&eigenEigenSolver_pseudoEigenvectors<myMatrixXd>)
    .method("setMaxIterations",&eigenEigenSolver_setMaxIterations<myMatrixXd>)
  ;
}
