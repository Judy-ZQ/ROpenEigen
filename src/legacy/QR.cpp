#include "ROpenEigen_types.h"

// template<typename _MatrixType>
// class EigenHouseholderQR : public Eigen::HouseholderQR<_MatrixType>{
// public:
//   EigenHouseholderQR():Eigen::HouseholderQR<_MatrixType>(){}
//   EigenHouseholderQR(_MatrixType &m):Eigen::HouseholderQR<_MatrixType>(m){}
//   EigenHouseholderQR(int rows, int cols):Eigen::HouseholderQR<_MatrixType>(rows,cols){}
//
//   _MatrixType householderQ(){
//     return Eigen::HouseholderQR<_MatrixType>::householderQ();
//   }
// };

// typedef EigenHouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> myMatrixXd;
// typedef Eigen::EigenHouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myMatrixXd;
// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> myMatrixXd;
// typedef Eigen::HouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myHouseholderQR4MatrixXd;
// typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> myMatrixXcd;
// typedef Eigen::HouseholderQR<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> myHouseholderQR4MatrixXcd;


typedef Eigen::HouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myHouseholderQR4MatrixXd;
typedef Eigen::HouseholderQR<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> myHouseholderQR4MatrixXcd;

typedef Eigen::CompleteOrthogonalDecomposition<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myCompleteOrthogonalDecomposition4MatrixXd;
typedef Eigen::CompleteOrthogonalDecomposition<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> myCompleteOrthogonalDecomposition4MatrixXcd;

typedef Eigen::ColPivHouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myColPivHouseholderQR4MatrixXd;
typedef Eigen::ColPivHouseholderQR<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> myColPivHouseholderQR4MatrixXcd;

typedef Eigen::FullPivHouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myFullPivHouseholderQR4MatrixXd;
typedef Eigen::FullPivHouseholderQR<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> myFullPivHouseholderQR4MatrixXcd;

// RCPP_MODULE(EigenHouseholderQRModule){
//   class_<myMatrixXd>("HouseholderQR4MatrixXd")
//   .default_constructor()
//   .constructor<myMatrixXd>()
//   .constructor<int,int>()
//   .method("householderQ",&Eigen::EigenHouseholderQR<myMatrixXd>::householderQ)
//   ;
// }

// template<typename _MatrixType>
// SEXP eigenHouseholderQR_compute(Eigen::EigenHouseholderQR<_MatrixType>* obj,SEXP b){
//   try{
//     typedef typename _MatrixType::RealScalar Scalar;
//     SEXP result;
//     if(isItMatrix(b)){
//       Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(b);
//       return Rcpp::wrap(obj->compute(_b));
//     }
//     // return result;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// template<typename _MatrixType>
// SEXP eigenHouseholderQR_compute(Eigen::EigenHouseholderQR<_MatrixType>* obj,_MatrixType b){
//   try{
//     // Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(b);
//     return Rcpp::wrap(obj->compute(b));
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _matrixType>
_matrixType eigenHouseholderQR_householderQ(Eigen::HouseholderQR<_matrixType>* obj){
  try{
    return obj->householderQ();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenHouseholderQR_matrixQR(Eigen::HouseholderQR<_matrixType>* obj){
  try{
    return obj->matrixQR();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenHouseholderQR_matrixQRUpper(Eigen::HouseholderQR<_matrixType>* obj){
  try{
    return obj->matrixQR().template triangularView<Eigen::Upper>();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenColPivHouseholderQR_absDeterminant(Eigen::ColPivHouseholderQR<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    double res = obj->absDeterminant();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
void eigenColPivHouseholderQR_compute(Eigen::ColPivHouseholderQR<_MatrixType>* obj,SEXP b){
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
int eigenColPivHouseholderQR_dimensionOfKernel(Eigen::ColPivHouseholderQR<_MatrixType>* obj){
  try{
    int res = obj->dimensionOfKernel();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenColPivHouseholderQR_hCoeffs(Eigen::ColPivHouseholderQR<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = obj->hCoeffs();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenColPivHouseholderQR_householderQ(Eigen::ColPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->householderQ();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenColPivHouseholderQR_inverse(Eigen::ColPivHouseholderQR<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->inverse();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenColPivHouseholderQR_isInjective(Eigen::ColPivHouseholderQR<_MatrixType>* obj){
  try{
    bool res = obj->isInjective();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenColPivHouseholderQR_isInvertible(Eigen::ColPivHouseholderQR<_MatrixType>* obj){
  try{
    bool res = obj->isInvertible();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenColPivHouseholderQR_isSurjective(Eigen::ColPivHouseholderQR<_MatrixType>* obj){
  try{
    bool res = obj->isSurjective();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenColPivHouseholderQR_logAbsDeterminant(Eigen::ColPivHouseholderQR<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    double res = obj->logAbsDeterminant();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenColPivHouseholderQR_matrixQR(Eigen::ColPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->matrixQR();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenColPivHouseholderQR_matrixR(Eigen::ColPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->matrixR();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenColPivHouseholderQR_matrixQRUpper(Eigen::ColPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->matrixQR().template triangularView<Eigen::Upper>();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
double eigenColPivHouseholderQR_maxPivot(Eigen::ColPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->maxPivot();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
int eigenColPivHouseholderQR_nonzeroPivots(Eigen::ColPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->nonzeroPivots();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenColPivHouseholderQR_rank(Eigen::ColPivHouseholderQR<_MatrixType>* obj){
  try{
    int res = obj->rank();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
void eigenColPivHouseholderQR_setThreshold(Eigen::ColPivHouseholderQR<_MatrixType>* obj,double b){
  try{
    obj->setThreshold(b);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenColPivHouseholderQR_Solve(Eigen::ColPivHouseholderQR<_MatrixType>* obj,SEXP b){
  try{
    // typedef _MatrixType::RealScalar Scalar;
    typedef typename _MatrixType::RealScalar Scalar;
    SEXP result;
    if(isItVector(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(b);
      Eigen::ColPivHouseholderQR<_MatrixType> decomposer(*obj);
      // Eigen::Solve<Eigen::ColPivHouseholderQR, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > res = obj->solve(_b);
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = decomposer.solve(_b);
      return Rcpp::wrap(res);
    }
    // return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenColPivHouseholderQR_threshold(Eigen::ColPivHouseholderQR<_MatrixType>* obj){
  try{
    return obj->threshold();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenCompleteOrthogonalDecomposition_absDeterminant(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    double res = obj->absDeterminant();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
SEXP eigenCompleteOrthogonalDecomposition_colsPermutation(Eigen::CompleteOrthogonalDecomposition<_matrixType>* obj){
  try{
    return Rcpp::wrap(obj->colsPermutation());
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenCompleteOrthogonalDecomposition_dimensionOfKernel(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj){
  try{
    int res = obj->dimensionOfKernel();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenCompleteOrthogonalDecomposition_hCoeffs(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = obj->hCoeffs();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenCompleteOrthogonalDecomposition_householderQ(Eigen::CompleteOrthogonalDecomposition<_matrixType>* obj){
  try{
    return obj->householderQ();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenCompleteOrthogonalDecomposition_isInjective(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj){
  try{
    bool res = obj->isInjective();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenCompleteOrthogonalDecomposition_isInvertible(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj){
  try{
    bool res = obj->isInvertible();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenCompleteOrthogonalDecomposition_isSurjective(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj){
  try{
    bool res = obj->isSurjective();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenCompleteOrthogonalDecomposition_logAbsDeterminant(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    double res = obj->logAbsDeterminant();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenCompleteOrthogonalDecomposition_matrixQTZ(Eigen::CompleteOrthogonalDecomposition<_matrixType>* obj){
  try{
    return obj->matrixQTZ();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenCompleteOrthogonalDecomposition_matrixT(Eigen::CompleteOrthogonalDecomposition<_matrixType>* obj){
  try{
    return obj->matrixT();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenCompleteOrthogonalDecomposition_matrixZ(Eigen::CompleteOrthogonalDecomposition<_matrixType>* obj){
  try{
    return obj->matrixZ();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
double eigenCompleteOrthogonalDecomposition_maxPivot(Eigen::CompleteOrthogonalDecomposition<_matrixType>* obj){
  try{
    return obj->maxPivot();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
int eigenCompleteOrthogonalDecomposition_nonzeroPivots(Eigen::CompleteOrthogonalDecomposition<_matrixType>* obj){
  try{
    return obj->nonzeroPivots();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenCompleteOrthogonalDecomposition_rank(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj){
  try{
    int res = obj->rank();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenCompleteOrthogonalDecomposition_setThreshold(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj,double b){
  try{
    return Rcpp::wrap(obj->setThreshold(b));
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenCompleteOrthogonalDecomposition_Solve(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj,SEXP b){
  try{
    // typedef _MatrixType::RealScalar Scalar;
    typedef typename _MatrixType::RealScalar Scalar;
    SEXP result;
    if(isItVector(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(b);
      Eigen::CompleteOrthogonalDecomposition<_MatrixType> decomposer(*obj);
      // Eigen::Solve<Eigen::ColPivHouseholderQR, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > res = obj->solve(_b);
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = decomposer.solve(_b);
      return Rcpp::wrap(res);
    }
    // return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenCompleteOrthogonalDecomposition_threshold(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj){
  try{
    return obj->threshold();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenCompleteOrthogonalDecomposition_zCoeffs(Eigen::CompleteOrthogonalDecomposition<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = obj->zCoeffs();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenFullPivHouseholderQR_absDeterminant(Eigen::FullPivHouseholderQR<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    double res = obj->absDeterminant();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
SEXP eigenFullPivHouseholderQR_colsPermutation(Eigen::FullPivHouseholderQR<_matrixType>* obj){
  try{
    return Rcpp::wrap(obj->colsPermutation());
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
void eigenFullPivHouseholderQR_compute(Eigen::FullPivHouseholderQR<_MatrixType>* obj,SEXP b){
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
int eigenFullPivHouseholderQR_dimensionOfKernel(Eigen::FullPivHouseholderQR<_MatrixType>* obj){
  try{
    int res = obj->dimensionOfKernel();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenFullPivHouseholderQR_hCoeffs(Eigen::FullPivHouseholderQR<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = obj->hCoeffs();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenFullPivHouseholderQR_inverse(Eigen::FullPivHouseholderQR<_MatrixType>* obj){
  try{
    typedef typename _MatrixType::RealScalar Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->inverse();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenFullPivHouseholderQR_isInjective(Eigen::FullPivHouseholderQR<_MatrixType>* obj){
  try{
    bool res = obj->isInjective();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenFullPivHouseholderQR_isInvertible(Eigen::FullPivHouseholderQR<_MatrixType>* obj){
  try{
    bool res = obj->isInvertible();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
bool eigenFullPivHouseholderQR_isSurjective(Eigen::FullPivHouseholderQR<_MatrixType>* obj){
  try{
    bool res = obj->isSurjective();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenFullPivHouseholderQR_logAbsDeterminant(Eigen::FullPivHouseholderQR<_MatrixType>* obj){
  try{
    // typedef typename _MatrixType::RealScalar Scalar;
    double res = obj->logAbsDeterminant();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenFullPivHouseholderQR_matrixQ(Eigen::FullPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->matrixQ();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
_matrixType eigenFullPivHouseholderQR_matrixQR(Eigen::FullPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->matrixQR();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}


template<typename _matrixType>
_matrixType eigenFullPivHouseholderQR_matrixQRUpper(Eigen::FullPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->matrixQR().template triangularView<Eigen::Upper>();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
double eigenFullPivHouseholderQR_maxPivot(Eigen::FullPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->maxPivot();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _matrixType>
int eigenFullPivHouseholderQR_nonzeroPivots(Eigen::FullPivHouseholderQR<_matrixType>* obj){
  try{
    return obj->nonzeroPivots();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
int eigenFullPivHouseholderQR_rank(Eigen::FullPivHouseholderQR<_MatrixType>* obj){
  try{
    int res = obj->rank();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _matrixType>
// SEXP eigenFullPivHouseholderQR_rowsTranspositions(Eigen::FullPivHouseholderQR<_matrixType>* obj){
//   try{
//     typedef typename _matrixType::RealScalar Scalar;
//     // Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> res = obj->rowsTranspositions();
//     Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->rowsTranspositions();
//     return Rcpp::wrap(res);
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _MatrixType>
void eigenFullPivHouseholderQR_setThreshold(Eigen::FullPivHouseholderQR<_MatrixType>* obj,double b){
  try{
    obj->setThreshold(b);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
SEXP eigenFullPivHouseholderQR_Solve(Eigen::FullPivHouseholderQR<_MatrixType>* obj,SEXP b){
  try{
    // typedef _MatrixType::RealScalar Scalar;
    typedef typename _MatrixType::RealScalar Scalar;
    SEXP result;
    if(isItVector(b)){
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(b);
      Eigen::FullPivHouseholderQR<_MatrixType> decomposer(*obj);
      // Eigen::Solve<Eigen::FullPivHouseholderQR, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > res = obj->solve(_b);
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = decomposer.solve(_b);
      return Rcpp::wrap(res);
    }
    // return result;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _MatrixType>
double eigenFullPivHouseholderQR_threshold(Eigen::FullPivHouseholderQR<_MatrixType>* obj){
  try{
    return obj->threshold();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}


RCPP_MODULE(EigenQRModule){
  class_<myHouseholderQR4MatrixXd>("HouseholderQR4MatrixXd")
  .default_constructor()
  .constructor<myMatrixXd>()
  .constructor<int,int>()
  .method("householderQ",&eigenHouseholderQR_householderQ<myMatrixXd>)
  .method("matrixQR",&eigenHouseholderQR_matrixQRUpper<myMatrixXd>)
  // .method("compute",&eigenHouseholderQR_compute<myMatrixXd>)
  ;

  class_<myHouseholderQR4MatrixXcd>("HouseholderQR4MatrixXcd")
    .default_constructor()
    .constructor<myMatrixXcd>()
    .constructor<int,int>()
    .method("householderQ",&eigenHouseholderQR_householderQ<myMatrixXcd>)
    // .method("matrixQR",&eigenHouseholderQR_matrixQR<myMatrixXd>)
    .method("matrixQR",&eigenHouseholderQR_matrixQRUpper<myMatrixXcd>)
  ;

  class_<myColPivHouseholderQR4MatrixXd>("ColPivHouseholderQR4MatrixXd")
    .default_constructor()
    .constructor<myMatrixXd>()
    .constructor<int,int>()
    .method("absDeterminant",&eigenColPivHouseholderQR_absDeterminant<myMatrixXd>)
    .method("compute",&eigenColPivHouseholderQR_compute<myMatrixXd>)
    .method("dimensionOfKernel",&eigenColPivHouseholderQR_dimensionOfKernel<myMatrixXd>)
    .method("hCoeffs",&eigenColPivHouseholderQR_hCoeffs<myMatrixXd>)
    .method("householderQ",&eigenColPivHouseholderQR_householderQ<myMatrixXd>)
    .method("inverse",&eigenColPivHouseholderQR_inverse<myMatrixXd>)
    .method("isInjective",&eigenColPivHouseholderQR_isInjective<myMatrixXd>)
    .method("isInvertible",&eigenColPivHouseholderQR_isInvertible<myMatrixXd>)
    .method("isSurjective",&eigenColPivHouseholderQR_isSurjective<myMatrixXd>)
    .method("matrixQR",&eigenColPivHouseholderQR_matrixQR<myMatrixXd>)
    .method("matrixR",&eigenColPivHouseholderQR_matrixR<myMatrixXd>)
    .method("matrixQRUpper",&eigenColPivHouseholderQR_matrixQRUpper<myMatrixXd>)
    .method("maxPivot",&eigenColPivHouseholderQR_maxPivot<myMatrixXd>)
    .method("nonzeroPivots",&eigenColPivHouseholderQR_nonzeroPivots<myMatrixXd>)
    .method("logAbsDeterminant",&eigenColPivHouseholderQR_logAbsDeterminant<myMatrixXd>)
    .method("solve",&eigenColPivHouseholderQR_Solve<myMatrixXd>)
    .method("rank",&eigenColPivHouseholderQR_rank<myMatrixXd>)
    .method("setThreshold",&eigenColPivHouseholderQR_setThreshold<myMatrixXd>)
    .method("threshold",&eigenColPivHouseholderQR_threshold<myMatrixXd>)
  ;

  class_<myCompleteOrthogonalDecomposition4MatrixXd>("CompleteOrthogonalDecomposition4MatrixXd")
    .default_constructor()
    .constructor<myMatrixXd>()
    .constructor<int,int>()
    .method("absDeterminant",&eigenCompleteOrthogonalDecomposition_absDeterminant<myMatrixXd>)
    .method("colsPermutation",&eigenCompleteOrthogonalDecomposition_colsPermutation<myMatrixXd>)
    .method("dimensionOfKernel",&eigenCompleteOrthogonalDecomposition_dimensionOfKernel<myMatrixXd>)
    .method("hCoeffs",&eigenCompleteOrthogonalDecomposition_hCoeffs<myMatrixXd>)
    .method("householderQ",&eigenCompleteOrthogonalDecomposition_householderQ<myMatrixXd>)
    .method("isInjective",&eigenCompleteOrthogonalDecomposition_isInjective<myMatrixXd>)
    .method("isInvertible",&eigenCompleteOrthogonalDecomposition_isInvertible<myMatrixXd>)
    .method("isSurjective",&eigenCompleteOrthogonalDecomposition_isSurjective<myMatrixXd>)
    .method("logAbsDeterminant",&eigenCompleteOrthogonalDecomposition_logAbsDeterminant<myMatrixXd>)
    .method("matrixQTZ",&eigenCompleteOrthogonalDecomposition_matrixQTZ<myMatrixXd>)
    .method("matrixT",&eigenCompleteOrthogonalDecomposition_matrixT<myMatrixXd>)
    .method("matrixZ",&eigenCompleteOrthogonalDecomposition_matrixZ<myMatrixXd>)
    .method("maxPivot",&eigenCompleteOrthogonalDecomposition_maxPivot<myMatrixXd>)
    .method("nonzeroPivots",&eigenCompleteOrthogonalDecomposition_nonzeroPivots<myMatrixXd>)
    .method("rank",&eigenCompleteOrthogonalDecomposition_rank<myMatrixXd>)
    .method("setThreshold",&eigenCompleteOrthogonalDecomposition_setThreshold<myMatrixXd>)
    .method("solve",&eigenCompleteOrthogonalDecomposition_Solve<myMatrixXd>)
    .method("threshold",&eigenCompleteOrthogonalDecomposition_threshold<myMatrixXd>)
  ;

  class_<myFullPivHouseholderQR4MatrixXd>("FullPivHouseholderQR4MatrixXd")
    .default_constructor()
    .constructor<myMatrixXd>()
    .constructor<int,int>()
    .method("absDeterminant",&eigenFullPivHouseholderQR_absDeterminant<myMatrixXd>)
    .method("compute",&eigenFullPivHouseholderQR_compute<myMatrixXd>)
    .method("colsPermutation",&eigenFullPivHouseholderQR_colsPermutation<myMatrixXd>)
    .method("dimensionOfKernel",&eigenFullPivHouseholderQR_dimensionOfKernel<myMatrixXd>)
    .method("hCoeffs",&eigenFullPivHouseholderQR_hCoeffs<myMatrixXd>)
    .method("inverse",&eigenFullPivHouseholderQR_inverse<myMatrixXd>)
    .method("isInjective",&eigenFullPivHouseholderQR_isInjective<myMatrixXd>)
    .method("isInvertible",&eigenFullPivHouseholderQR_isInvertible<myMatrixXd>)
    .method("isSurjective",&eigenFullPivHouseholderQR_isSurjective<myMatrixXd>)
    .method("logAbsDeterminant",&eigenFullPivHouseholderQR_logAbsDeterminant<myMatrixXd>)
    .method("matrixQ",&eigenFullPivHouseholderQR_matrixQ<myMatrixXd>)
    .method("matrixQR",&eigenFullPivHouseholderQR_matrixQR<myMatrixXd>)
    .method("matrixQRUpper",&eigenFullPivHouseholderQR_matrixQRUpper<myMatrixXd>)
    .method("maxPivot",&eigenFullPivHouseholderQR_maxPivot<myMatrixXd>)
    .method("nonzeroPivots",&eigenFullPivHouseholderQR_nonzeroPivots<myMatrixXd>)
    .method("solve",&eigenFullPivHouseholderQR_Solve<myMatrixXd>)
    .method("rank",&eigenFullPivHouseholderQR_rank<myMatrixXd>)
    // .method("rowsTranspositions",&eigenFullPivHouseholderQR_rowsTranspositions<myMatrixXd>)
    .method("setThreshold",&eigenFullPivHouseholderQR_setThreshold<myMatrixXd>)
    .method("threshold",&eigenFullPivHouseholderQR_threshold<myMatrixXd>)
  ;
}
