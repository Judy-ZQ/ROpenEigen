#include "ROpenEigen_types.h"

template<typename _MatrixType>
void eigenBDCSVD_compute(Eigen::BDCSVD<_MatrixType>* obj,SEXP p){
  try{
    Rcpp::List param(p);
    typedef typename _MatrixType::RealScalar Scalar;
    if(param.length()==1){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(param[0]);
      obj->compute(matrix);
    }else if(param.length()==2){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(param[0]);
      unsigned int computationOptions = param[1];
      obj->compute(matrix,computationOptions);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenBDCSVD_computeU(Eigen::BDCSVD<_MatrixType>* obj){
  try{
      return obj->computeU();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenBDCSVD_computeV(Eigen::BDCSVD<_MatrixType>* obj){
  try{
    return obj->computeV();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenBDCSVD_matrixU(Eigen::BDCSVD<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->matrixU();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenBDCSVD_matrixV(Eigen::BDCSVD<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->matrixV();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenBDCSVD_nonzeroSingularValues(Eigen::BDCSVD<_MatrixType>* obj){
  try{
    return obj->nonzeroSingularValues();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenBDCSVD_rank(Eigen::BDCSVD<_MatrixType>* obj){
  try{
    return obj->rank();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
void eigenBDCSVD_setThreshold(Eigen::BDCSVD<_MatrixType>* obj,double threshold){
  try{
    obj->setThreshold(threshold);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenBDCSVD_singularValues(Eigen::BDCSVD<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<double, Eigen::Dynamic, 1> result = obj->singularValues();
    return Rcpp::wrap(result);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}


template<typename _MatrixType>
SEXP eigenBDCSVD_Solve(Eigen::BDCSVD<_MatrixType>* obj,SEXP b){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    SEXP result;
    if(isItMatrix(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(b);
      Eigen::BDCSVD<_MatrixType> decomposer(*obj);
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = decomposer.solve(_b);
      return Rcpp::wrap(res);
    }else if(isItVector(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(b);
      Eigen::BDCSVD<_MatrixType> decomposer(*obj);
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = decomposer.solve(_b);
      return Rcpp::wrap(res);
    }
    // return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenBDCSVD_threshold(Eigen::BDCSVD<_MatrixType>* obj){
  try{
    return obj->threshold();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
void eigenJacobiSVD_compute(Eigen::JacobiSVD<_MatrixType>* obj,SEXP p){
  try{
    Rcpp::List param(p);
    typedef typename _MatrixType::RealScalar Scalar;
    if(param.length()==1){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(param[0]);
      obj->compute(matrix);
    }else if(param.length()==2){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(param[0]);
      unsigned int computationOptions = param[1];
      obj->compute(matrix,computationOptions);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenJacobiSVD_computeU(Eigen::JacobiSVD<_MatrixType>* obj){
  try{
    return obj->computeU();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenJacobiSVD_computeV(Eigen::JacobiSVD<_MatrixType>* obj){
  try{
    return obj->computeV();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenJacobiSVD_matrixU(Eigen::JacobiSVD<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->matrixU();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenJacobiSVD_matrixV(Eigen::JacobiSVD<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->matrixV();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenJacobiSVD_nonzeroSingularValues(Eigen::JacobiSVD<_MatrixType>* obj){
  try{
    return obj->nonzeroSingularValues();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenJacobiSVD_rank(Eigen::JacobiSVD<_MatrixType>* obj){
  try{
    return obj->rank();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
void eigenJacobiSVD_setThreshold(Eigen::JacobiSVD<_MatrixType>* obj,double threshold){
  try{
    obj->setThreshold(threshold);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenJacobiSVD_singularValues(Eigen::JacobiSVD<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<double, Eigen::Dynamic, 1> result = obj->singularValues();
    return Rcpp::wrap(result);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}


template<typename _MatrixType>
SEXP eigenJacobiSVD_Solve(Eigen::JacobiSVD<_MatrixType>* obj,SEXP b){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    SEXP result;
    if(isItMatrix(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(b);
      Eigen::JacobiSVD<_MatrixType> decomposer(*obj);
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = decomposer.solve(_b);
      return Rcpp::wrap(res);
    }else if(isItVector(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(b);
      Eigen::JacobiSVD<_MatrixType> decomposer(*obj);
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = decomposer.solve(_b);
      return Rcpp::wrap(res);
    }
    // return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenJacobiSVD_threshold(Eigen::JacobiSVD<_MatrixType>* obj){
  try{
    return obj->threshold();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

typedef Eigen::BDCSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myBDCSVD4MatrixXd;
typedef Eigen::BDCSVD<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> myBDCSVD4MatrixXcd;

typedef Eigen::JacobiSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myJacobiSVD4MatrixXd;
typedef Eigen::JacobiSVD<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> myJacobiSVD4MatrixXcd;


RCPP_MODULE(EigenSVDModule){
  class_<myBDCSVD4MatrixXd>("BDCSVD4MatrixXd")
  .default_constructor()
  .constructor<myMatrixXd,int>()
  .constructor<int,int,int>()
  .method("compute_",&eigenBDCSVD_compute<myMatrixXd>)
  .method("computeU",&eigenBDCSVD_computeU<myMatrixXd>)
  .method("computeV",&eigenBDCSVD_computeV<myMatrixXd>)
  .method("matrixU",&eigenBDCSVD_matrixU<myMatrixXd>)
  .method("matrixV",&eigenBDCSVD_matrixV<myMatrixXd>)
  .method("nonzeroSingularValues",&eigenBDCSVD_nonzeroSingularValues<myMatrixXd>)
  .method("rank",&eigenBDCSVD_rank<myMatrixXd>)
  .method("setThreshold",&eigenBDCSVD_setThreshold<myMatrixXd>)
  .method("singularValues",&eigenBDCSVD_singularValues<myMatrixXd>)
  .method("solve",&eigenBDCSVD_Solve<myMatrixXd>)
  .method("threshold",&eigenBDCSVD_threshold<myMatrixXd>)
  ;

  class_<myJacobiSVD4MatrixXd>("JacobiSVD4MatrixXd")
    .default_constructor()
    .constructor<myMatrixXd,int>()
    .constructor<int,int,int>()
    .method("compute_",&eigenJacobiSVD_compute<myMatrixXd>)
    .method("computeU",&eigenJacobiSVD_computeU<myMatrixXd>)
    .method("computeV",&eigenJacobiSVD_computeV<myMatrixXd>)
    .method("matrixU",&eigenJacobiSVD_matrixU<myMatrixXd>)
    .method("matrixV",&eigenJacobiSVD_matrixV<myMatrixXd>)
    .method("nonzeroSingularValues",&eigenJacobiSVD_nonzeroSingularValues<myMatrixXd>)
    .method("rank",&eigenJacobiSVD_rank<myMatrixXd>)
    .method("setThreshold",&eigenJacobiSVD_setThreshold<myMatrixXd>)
    .method("singularValues",&eigenJacobiSVD_singularValues<myMatrixXd>)
    .method("solve",&eigenJacobiSVD_Solve<myMatrixXd>)
    .method("threshold",&eigenJacobiSVD_threshold<myMatrixXd>)
  ;
}
