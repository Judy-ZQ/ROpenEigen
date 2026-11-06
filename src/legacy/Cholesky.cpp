#include "ROpenEigen_types.h"

template<typename _MatrixType,int _UpLo>
void eigenLDLT_compute(Eigen::LDLT<_MatrixType,_UpLo>* obj,SEXP b){
  try{
    // typedef _MatrixType::RealScalar Scalar;
    typedef typename _MatrixType::RealScalar Scalar;
    // SEXP result;
    if(isItMatrix(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(b);
      // Eigen::LDLT<_MatrixType,Eigen::Lower> decomposer(*obj);
      // Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = decomposer.solve(_b);
      obj->compute(_b);
    }
    // return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType,int _UpLo>
bool eigenLDLT_isNegative(Eigen::LDLT<_MatrixType,_UpLo>* obj){
  try{
    return obj->isNegative();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType,int _UpLo>
bool eigenLDLT_isPositive(Eigen::LDLT<_MatrixType,_UpLo>* obj){
  try{
    return obj->isPositive();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _MatrixType,int _UpLo>
// SEXP eigenLDLT_inverse(Eigen::LDLT<_MatrixType,_UpLo>* obj){
//   try{
//     typedef typename _MatrixType::RealScalar Scalar;
//     Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->inverse();
//     return Rcpp::wrap(res);
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _MatrixType,int _UpLo>
_MatrixType eigenLDLT_matrixL(Eigen::LDLT<_MatrixType,_UpLo>* obj){
  try{
    // typedef _MatrixType::RealScalar Scalar;
    // typedef typename _MatrixType::RealScalar Scalar;
    // SEXP result;
    // if(isItMatrix(b)){
    //   Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(b);
    //   Eigen::LDLT<_MatrixType,_UpLo> decomposer(*obj);
    //   Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = decomposer.solve(_b);
    //   return Rcpp::wrap(res);
    // }
    // return result;
    // Eigen::LDLT<_MatrixType,_UpLo> decomposer(*obj);
    return obj->matrixL();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType,int _UpLo>
_MatrixType eigenLDLT_matrixLDLT(Eigen::LDLT<_MatrixType,_UpLo>* obj){
  try{
    // Eigen::LDLT<_MatrixType,_UpLo> decomposer(*obj);
    // return decomposer.matrixLDLT();
    return obj->matrixLDLT();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType,int _UpLo>
_MatrixType eigenLDLT_matrixU(Eigen::LDLT<_MatrixType,_UpLo>* obj){
  try{
    // Eigen::LDLT<_MatrixType,_UpLo> decomposer(*obj);
    return obj->matrixU();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType,int _UpLo>
void eigenLDLT_rankUpdate(Eigen::LDLT<_MatrixType,_UpLo>* obj,SEXP u,double sigma){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    if(isItVector(u)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> _u=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(u);
      // Eigen::LDLT<_MatrixType,_UpLo> decomposer(*obj);
      obj->rankUpdate(_u, sigma);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType,int _UpLo>
double eigenLDLT_rcond(Eigen::LDLT<_MatrixType,_UpLo>* obj){
  try{
    double res = obj->rcond();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType,int _UpLo>
SEXP eigenLDLT_reconstructedMatrix(Eigen::LDLT<_MatrixType,_UpLo>* obj){
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

template<typename _MatrixType,int _UpLo>
void eigenLDLT_setZero(Eigen::LDLT<_MatrixType,_UpLo>* obj){
  try{
    obj->setZero();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType,int _UpLo>
SEXP eigenLDLT_Solve(Eigen::LDLT<_MatrixType,_UpLo>* obj,SEXP b){
  try{
    // typedef _MatrixType::RealScalar Scalar;
    typedef typename _MatrixType::RealScalar Scalar;
    // SEXP result;
    if(isItMatrix(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(b);
      // Eigen::LDLT<_MatrixType,_UpLo> decomposer(*obj);
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->solve(_b);
      return Rcpp::wrap(res);
    }
    // return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType,int _UpLo>
SEXP eigenLDLT_vectorD(Eigen::LDLT<_MatrixType,_UpLo>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    // Eigen::LDLT<_MatrixType,_UpLo> decomposer(*obj);
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = obj->vectorD();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType,int _UpLo>
SEXP eigenLDLT_transpositionsP(Eigen::LDLT<_MatrixType,_UpLo>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    // Eigen::LDLT<_MatrixType,_UpLo> decomposer(*obj);
    Eigen::Transpositions<Eigen::Dynamic, Eigen::Dynamic> res = obj->transpositionsP();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

typedef Eigen::LDLT<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,Eigen::Lower> myLDLT4MatrixXd;
typedef Eigen::LDLT<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>,Eigen::Lower> myLDLT4MatrixXcd;




RCPP_MODULE(EigenCholeskyModule){
  class_<myLDLT4MatrixXd>("LDLT4MatrixXd")
  .default_constructor()
  .constructor<myMatrixXd>()
  .constructor<int>()
  .method("compute",&eigenLDLT_compute<myMatrixXd,Eigen::Lower>)
  .method("isNegative",&eigenLDLT_isNegative<myMatrixXd,Eigen::Lower>)
  .method("isPositive",&eigenLDLT_isPositive<myMatrixXd,Eigen::Lower>)
  .method("rankUpdate",&eigenLDLT_rankUpdate<myMatrixXd,Eigen::Lower>)
  .method("reconstructedMatrix",&eigenLDLT_reconstructedMatrix<myMatrixXd,Eigen::Lower>)
  .method("solve",&eigenLDLT_Solve<myMatrixXd,Eigen::Lower>)
  .method("matrixL",&eigenLDLT_matrixL<myMatrixXd,Eigen::Lower>)
  .method("matrixLDLT",&eigenLDLT_matrixLDLT<myMatrixXd,Eigen::Lower>)
  .method("matrixU",&eigenLDLT_matrixU<myMatrixXd,Eigen::Lower>)
  .method("vectorD",&eigenLDLT_vectorD<myMatrixXd,Eigen::Lower>)
  .method("transpositionsP",&eigenLDLT_transpositionsP<myMatrixXd,Eigen::Lower>)
  ;
}
