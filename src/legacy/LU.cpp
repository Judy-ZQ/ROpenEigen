#include "ROpenEigen_types.h"


template<typename _MatrixType>
void eigenFullPivLU_compute(Eigen::FullPivLU<_MatrixType>* obj,SEXP b){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    SEXP result;
    if(isItMatrix(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(b);
      obj->compute(_b);
    }else if(isItVector(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(b);
      obj->compute(_b);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenFullPivLU_determinant(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    double res = obj->determinant();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenFullPivLU_dimensionOfKernel(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    int res = obj->dimensionOfKernel();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenFullPivLU_image(Eigen::FullPivLU<_MatrixType>* obj,SEXP b){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    SEXP result;
    if(isItMatrix(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(b);
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->image(_b);
      return Rcpp::wrap(res);
    }

  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenFullPivLU_inverse(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->inverse();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenFullPivLU_isInjective(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    bool res = obj->isInjective();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenFullPivLU_isInvertible(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    bool res = obj->isInvertible();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenFullPivLU_isSurjective(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    bool res = obj->isSurjective();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenFullPivLU_kernel(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->kernel();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenFullPivLU_matrixLU(Eigen::FullPivLU<_matrixType>* obj){
  try{
    return obj->matrixLU();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
double eigenFullPivLU_maxPivot(Eigen::FullPivLU<_matrixType>* obj){
  try{
    return obj->maxPivot();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
int eigenFullPivLU_nonzeroPivots(Eigen::FullPivLU<_matrixType>* obj){
  try{
    return obj->nonzeroPivots();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
SEXP eigenFullPivLU_permutationP(Eigen::FullPivLU<_matrixType>* obj){
  try{
    return Rcpp::wrap(obj->permutationP());
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
SEXP eigenFullPivLU_permutationQ(Eigen::FullPivLU<_matrixType>* obj){
  try{
    return Rcpp::wrap(obj->permutationQ());
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenFullPivLU_rank(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    return obj->rank();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenFullPivLU_rcond(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    double res = obj->rcond();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenFullPivLU_reconstructedMatrix(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    // Eigen::LDLT<_MatrixType,_UpLo> decomposer(*obj);
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->reconstructedMatrix();
    return Rcpp::wrap(res);
    // return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
void eigenFullPivLU_setThreshold(Eigen::FullPivLU<_MatrixType>* obj,double b){
  try{
    obj->setThreshold(b);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenFullPivLU_Solve(Eigen::FullPivLU<_MatrixType>* obj,SEXP b){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    SEXP result;
    if(isItMatrix(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(b);
      Eigen::FullPivLU<_MatrixType> decomposer(*obj);
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = decomposer.solve(_b);
      return Rcpp::wrap(res);
    }else if(isItVector(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(b);
      Eigen::FullPivLU<_MatrixType> decomposer(*obj);
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = decomposer.solve(_b);
      return Rcpp::wrap(res);
    }
    // return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenFullPivLU_threshold(Eigen::FullPivLU<_MatrixType>* obj){
  try{
    return obj->threshold();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}


typedef Eigen::FullPivLU<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myFullPivLU4MatrixXd;
typedef Eigen::FullPivLU<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> myFullPivLU4MatrixXcd;




RCPP_MODULE(EigenLUModule){
  class_<myFullPivLU4MatrixXd>("FullPivLU4MatrixXd")
  .default_constructor()
  .constructor<myMatrixXd>()
  .constructor<int,int>()
  .method("compute",&eigenFullPivLU_compute<myMatrixXd>)
  .method("determinant",&eigenFullPivLU_determinant<myMatrixXd>)
  .method("dimensionOfKernel",&eigenFullPivLU_dimensionOfKernel<myMatrixXd>)
  .method("image",&eigenFullPivLU_image<myMatrixXd>)
  .method("inverse",&eigenFullPivLU_inverse<myMatrixXd>)
  .method("isInjective",&eigenFullPivLU_isInjective<myMatrixXd>)
  .method("isInvertible",&eigenFullPivLU_isInvertible<myMatrixXd>)
  .method("isSurjective",&eigenFullPivLU_isSurjective<myMatrixXd>)
  .method("kernel",&eigenFullPivLU_kernel<myMatrixXd>)
  .method("matrixLU",&eigenFullPivLU_matrixLU<myMatrixXd>)
  .method("maxPivot",&eigenFullPivLU_maxPivot<myMatrixXd>)
  .method("nonzeroPivots",&eigenFullPivLU_nonzeroPivots<myMatrixXd>)
  .method("permutationP",&eigenFullPivLU_permutationP<myMatrixXd>)
  .method("permutationQ",&eigenFullPivLU_permutationQ<myMatrixXd>)
  .method("rank",&eigenFullPivLU_rank<myMatrixXd>)
  .method("rcond",&eigenFullPivLU_rcond<myMatrixXd>)
  .method("reconstructedMatrix",&eigenFullPivLU_reconstructedMatrix<myMatrixXd>)
  .method("setThreshold",&eigenFullPivLU_setThreshold<myMatrixXd>)
  .method("solve",&eigenFullPivLU_Solve<myMatrixXd>)
  .method("threshold",&eigenFullPivLU_threshold<myMatrixXd>)
  ;

  // class_<myHouseholderQR4MatrixXcd>("HouseholderQR4MatrixXcd")
  //   .default_constructor()
  //   .constructor<myMatrixXcd>()
  //   .constructor<int,int>()
  //   .method("householderQ",&eigenHouseholderQR_householderQ<myMatrixXcd>)
  // // .method("matrixQR",&eigenHouseholderQR_matrixQR<myMatrixXd>)
  //    .method("matrixQR",&eigenHouseholderQR_matrixQRUpper<myMatrixXcd>)
  // ;
}


