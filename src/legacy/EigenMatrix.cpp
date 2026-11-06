// #ifndef ROpenEigen_TYPE_H_
// #include "ROpenEigen_types.h"
// #endif
#include "ROpenEigen_types.h"


// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> eigenArrayXXdConstant(int rows,int cols,double value){
  return Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>::Constant(rows,cols,value);
}

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> eigenArrayXXdOnes(int rows,int cols){
  return Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>::Ones(rows,cols);
}

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> eigenArrayXXdRandom(int rows,int cols){
  return Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>::Random(rows,cols);
}

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> eigenArrayXXdZero(int rows,int cols){
  return Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows,cols);
}

// // [[Rcpp::export]]
// Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> eigenArrayXXdIdentity(int rows,int cols){
//   return Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>::Identity(rows,cols);
// }

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, 1> eigenArrayXdConstant(int size,double value){
  return Eigen::Array<double, Eigen::Dynamic, 1>::Constant(size,value);
}

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, 1> eigenArrayXdLinSpaced(int size,double low,double high){
  return Eigen::Array<double, Eigen::Dynamic, 1>::LinSpaced(size,low,high);
}

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, 1> eigenArrayXdOnes(int size){
  return Eigen::Array<double, Eigen::Dynamic, 1>::Ones(size);
}

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, 1> eigenArrayXdRandom(int size){
  return Eigen::Array<double, Eigen::Dynamic, 1>::Random(size);
}

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, 1> eigenArrayXdZero(int size){
  return Eigen::Array<double, Eigen::Dynamic, 1>::Zero(size);
}

// // [[Rcpp::export]]
// Eigen::Array<double, Eigen::Dynamic, 1> eigenArrayXdUnit(int index,int i){
//   return Eigen::Array<double, Eigen::Dynamic, 1>::Unit(index,i);
// }

// [[Rcpp::export]]
Eigen::Array<double, 1, Eigen::Dynamic> eigenRowArrayXdConstant(int size,double value){
  return Eigen::Array<double, 1, Eigen::Dynamic>::Constant(size,value);
}

// [[Rcpp::export]]
Eigen::Array<double, 1, Eigen::Dynamic> eigenRowArrayXdLinSpaced(int size,double low,double high){
  return Eigen::Array<double, 1, Eigen::Dynamic>::LinSpaced(size,low,high);
}

// [[Rcpp::export]]
Eigen::Array<double, 1, Eigen::Dynamic> eigenRowArrayXdOnes(int size){
  return Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(size);
}

// [[Rcpp::export]]
Eigen::Array<double, 1, Eigen::Dynamic> eigenRowArrayXdRandom(int size){
  return Eigen::Array<double, 1, Eigen::Dynamic>::Random(size);
}

// [[Rcpp::export]]
Eigen::Array<double, 1, Eigen::Dynamic> eigenRowArrayXdZero(int size){
  return Eigen::Array<double, 1, Eigen::Dynamic>::Zero(size);
}



// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> eigenArrayXXdMul_matAndmat(Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> a,Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> b){
  return a*b;
}

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_matAndcmat(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a,Eigen::Matrix<double, Eigen::Dynamic, 1> b){
//   return a*b;
// }

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_rmatAndmat(Eigen::Matrix<double, 1, Eigen::Dynamic> a,Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b){
//   return a*b;
// }

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_cmatAndrmat(Eigen::Matrix<double, Eigen::Dynamic, 1> a,Eigen::Matrix<double, 1, Eigen::Dynamic> b){
//   // Rcout << a*b << endl;
//   return a*b;
// }


// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_rmatAndcmat(Eigen::Matrix<double, 1, Eigen::Dynamic> a,Eigen::Matrix<double, Eigen::Dynamic, 1> b){
//   // Rcout << a*b << endl;
//   return a*b;
// }

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_sAndmat(double a,Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b){
//   Rcout << a*b << endl;
//   return a*b;
// }

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> eigenMArrayXXdMul_sAndmat(Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> a,double b){
  // Rcout << a*b << endl;
  return a*b;
}

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> eigenArrayXXdDiv_sAndmat(Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> a,double b){
  // Rcout << a*b << endl;
  return a/b;
}

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> eigenArrayXXdAdd_matAndmat(Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> a,Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> b){
  return a+b;
}

// [[Rcpp::export]]
Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> eigenArrayXXdAdd_sAndmat(Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> a,double b){
  return a+b;
}

// [[Rcpp::export]]
Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXcfConstant(int rows,int cols,double value){
  return Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>::Constant(rows,cols,value);
}

// [[Rcpp::export]]
Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXcfOnes(int rows,int cols){
  return Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>::Ones(rows,cols);
}

// [[Rcpp::export]]
Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXcfRandom(int rows,int cols){
  return Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>::Random(rows,cols);
}

// [[Rcpp::export]]
Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXcfZero(int rows,int cols){
  return Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows,cols);
}

// [[Rcpp::export]]
Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXcfIdentity(int rows,int cols){
  return Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>::Identity(rows,cols);
}

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdConstant(int rows,int cols,double value){
  return Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Constant(rows,cols,value);
}

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdOnes(int rows,int cols){
  return Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Ones(rows,cols);
}

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdRandom(int rows,int cols){
  return Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Random(rows,cols);
}

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdZero(int rows,int cols){
  return Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(rows,cols);
}

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdIdentity(int rows,int cols){
  return Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Identity(rows,cols);
}

// [[Rcpp::export]]
SEXP eigenMatrixXdMul_matAndmat(SEXP a,SEXP b){
  if(isItMatrix(a) && isItMatrix(b)){
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _b = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(b);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res = _a * _b;
    return Rcpp::wrap(res);
  }else if(isItMatrix(a) && isItVector(b)){
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, 1> _b = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(b);
    Eigen::Matrix<double, Eigen::Dynamic, 1> res = _a * _b;
    return Rcpp::wrap(res);
  }else if(isItRowVector(a) && isItMatrix(b)){
    Eigen::Matrix<double, 1, Eigen::Dynamic> _a = Rcpp::as<Eigen::Matrix<double, 1, Eigen::Dynamic>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _b = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(b);
    Eigen::Matrix<double, 1, Eigen::Dynamic> res = _a * _b;
    // return Rcpp::wrap(_a + _b); //undefined reference to `SEXPREC* Rcpp::RcppEigen::eigen_wrap<Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op
    return Rcpp::wrap(res);
  }else if(isItRowVector(a) && isItVector(b)){
    Eigen::Matrix<double, 1, Eigen::Dynamic> _a = Rcpp::as<Eigen::Matrix<double, 1, Eigen::Dynamic>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, 1> _b = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(b);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res = _a * _b;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res(res1);
    return Rcpp::wrap(res);
  }else if(isItVector(a) && isItRowVector(b)){
    Eigen::Matrix<double, Eigen::Dynamic, 1> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(a);
    Eigen::Matrix<double, 1, Eigen::Dynamic> _b = Rcpp::as<Eigen::Matrix<double, 1, Eigen::Dynamic>>(b);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res = _a * _b;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res(res1);
    return Rcpp::wrap(res);
  }
  Rcpp::stop("parameter error.");
}

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_matAndmat(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a,Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b){
//   // Rcout << a << endl;
//   // Rcout << b << endl;
//   // Rcout << a*b << endl;
//   return a*b;
// }

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_matAndcmat(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a,Eigen::Matrix<double, Eigen::Dynamic, 1> b){
//   return a*b;
// }
//
// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_rmatAndmat(Eigen::Matrix<double, 1, Eigen::Dynamic> a,Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b){
//   return a*b;
// }
//
// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_cmatAndrmat(Eigen::Matrix<double, Eigen::Dynamic, 1> a,Eigen::Matrix<double, 1, Eigen::Dynamic> b){
//   // Rcout << a*b << endl;
//   return a*b;
// }
//
//
// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_rmatAndcmat(Eigen::Matrix<double, 1, Eigen::Dynamic> a,Eigen::Matrix<double, Eigen::Dynamic, 1> b){
//   // Rcout << a*b << endl;
//   return a*b;
// }

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_sAndmat(double a,Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b){
//   Rcout << a*b << endl;
//   return a*b;
// }

// [[Rcpp::export]]
SEXP eigenMatrixXdMul_sAndmat(SEXP a,double b){
  if(isItMatrix(a)){
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res = _a * b;
    return Rcpp::wrap(res);
  }else if(isItVector(a)){
    Eigen::Matrix<double, Eigen::Dynamic, 1> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, 1> res = _a * b;
    return Rcpp::wrap(res);
  }
  Rcpp::stop("parameter error.");
}

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdMul_sAndmat(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a,double b){
//   // Rcout << a*b << endl;
//   return a*b;
// }

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdMul_sAndmat(Eigen::Matrix<double, Eigen::Dynamic, 1> a,double b){
//   // Rcout << a*b << endl;
//   return a*b;
// }
//
// // [[Rcpp::export]]
// Eigen::Matrix<double, 1, Eigen::Dynamic> eigenRowVectorXdMul_sAndmat(Eigen::Matrix<double, 1, Eigen::Dynamic> a,double b){
//   // Rcout << a*b << endl;
//   return a*b;
// }

// [[Rcpp::export]]
SEXP eigenMatrixXdDiv_sAndmat(SEXP a,double b){
  if(isItMatrix(a)){
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res = _a / b;
    return Rcpp::wrap(res);
  }else if(isItVector(a)){
    Eigen::Matrix<double, Eigen::Dynamic, 1> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, 1> res = _a / b;
    return Rcpp::wrap(res);
  }
  Rcpp::stop("parameter error.");
}

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdDiv_sAndmat(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a,double b){
//   // Rcout << a*b << endl;
//   return a/b;
// }

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdDiv_sAndmat(Eigen::Matrix<double, Eigen::Dynamic, 1> a,double b){
//   // Rcout << a*b << endl;
//   return a/b;
// }
//
// // [[Rcpp::export]]
// Eigen::Matrix<double, 1, Eigen::Dynamic> eigenRowVectorXdDiv_sAndmat(Eigen::Matrix<double, 1, Eigen::Dynamic> a,double b){
//   // Rcout << a*b << endl;
//   return a/b;
// }

// [[Rcpp::export]]
SEXP eigenMatrixXdAdd_matAndmat(SEXP a,SEXP b){
  if(isItMatrix(a) && isItMatrix(b)){
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _b = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(b);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res = _a + _b;
    return Rcpp::wrap(res);
  }else if(isItMatrix(a) && isItVector(b)){
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, 1> _b = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(b);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res = _a + _b;
    return Rcpp::wrap(res);
  }else if(isItVector(a) && isItMatrix(b)){
    Eigen::Matrix<double, Eigen::Dynamic, 1> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> _b = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>(b);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res = _a + _b;
    // return Rcpp::wrap(_a + _b); //undefined reference to `SEXPREC* Rcpp::RcppEigen::eigen_wrap<Eigen::CwiseBinaryOp<Eigen::internal::scalar_sum_op
    return Rcpp::wrap(res);
  }else if(isItVector(a) && isItVector(b)){
    Eigen::Matrix<double, Eigen::Dynamic, 1> _a = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(a);
    Eigen::Matrix<double, Eigen::Dynamic, 1> _b = Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(b);
    Eigen::Matrix<double, Eigen::Dynamic, 1> res = _a + _b;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> res(res1);
    return Rcpp::wrap(res);
  }
  Rcpp::stop("parameter error.");
}

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXdAdd_matAndmat(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a,Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b){
//   return a+b;
// }

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdAdd_matAndmat(Eigen::Matrix<double, Eigen::Dynamic, 1> a,Eigen::Matrix<double, Eigen::Dynamic, 1> b){
//   return a+b;
// }
//
// // [[Rcpp::export]]
// Eigen::Matrix<double, 1, Eigen::Dynamic> eigenRowVectorXdAdd_matAndmat(Eigen::Matrix<double, 1, Eigen::Dynamic> a,Eigen::Matrix<double, 1, Eigen::Dynamic> b){
//   return a+b;
// }

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdConstant(int size,double value){
  return Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(size,value);
}

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdLinSpaced(int size,double low,double high){
  return Eigen::Matrix<double, Eigen::Dynamic, 1>::LinSpaced(size,low,high);
}

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdOnes(int size){
  return Eigen::Matrix<double, Eigen::Dynamic, 1>::Ones(size);
}

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdRandom(int size){
  return Eigen::Matrix<double, Eigen::Dynamic, 1>::Random(size);
}

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdZero(int size){
  return Eigen::Matrix<double, Eigen::Dynamic, 1>::Zero(size);
}

// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdUnit(int index,int i){
  return Eigen::Matrix<double, Eigen::Dynamic, 1>::Unit(index,i);
}

// [[Rcpp::export]]
Eigen::Matrix<double, 1, Eigen::Dynamic> eigenRowVectorXdConstant(int size,double value){
  return Eigen::Matrix<double, 1, Eigen::Dynamic>::Constant(size,value);
}

// [[Rcpp::export]]
Eigen::Matrix<double, 1, Eigen::Dynamic> eigenRowVectorXdLinSpaced(int size,double low,double high){
  return Eigen::Matrix<double, 1, Eigen::Dynamic>::LinSpaced(size,low,high);
}

// [[Rcpp::export]]
Eigen::Matrix<double, 1, Eigen::Dynamic> eigenRowVectorXdOnes(int size){
  return Eigen::Matrix<double, 1, Eigen::Dynamic>::Ones(size);
}

// [[Rcpp::export]]
Eigen::Matrix<double, 1, Eigen::Dynamic> eigenRowVectorXdRandom(int size){
  return Eigen::Matrix<double, 1, Eigen::Dynamic>::Random(size);
}

// [[Rcpp::export]]
Eigen::Matrix<double, 1, Eigen::Dynamic> eigenRowVectorXdZero(int size){
  return Eigen::Matrix<double, 1, Eigen::Dynamic>::Zero(size);
}

// [[Rcpp::export]]
Eigen::Matrix<double, 1, Eigen::Dynamic> eigenRowVectorXdUnit(int index,int i){
  return Eigen::Matrix<double, 1, Eigen::Dynamic>::Unit(index,i);
}


// void std_sort_colVec(Eigen::Matrix<double, Eigen::Dynamic, 1>& v){
//   std::sort(v.begin(),v.end());
// }
//
//
// void std_sort_rowVec(Eigen::Matrix<double, 1, Eigen::Dynamic>& v){
//   std::sort(v.begin(),v.end());
// }

// // [[Rcpp::export]]
// void std_sort(SEXP v){
//   try{
//     if(isItVector(v)){
//       // Eigen::Matrix<double, Eigen::Dynamic, 1> *v1=new Eigen::Matrix<double, Eigen::Dynamic, 1>(Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(v));
//       Eigen::Matrix<double, Eigen::Dynamic, 1> v1(Rcpp::as<Eigen::Matrix<double, Eigen::Dynamic, 1>>(v));
//       // Eigen::Matrix<double, Eigen::Dynamic, 1> *v1=&_v;
//       // SEXP2EigenMatrix<Eigen::Matrix<double, Eigen::Dynamic, 1>,double>(v,*v1);
//       // Rcout << *v1 << endl;
//       // std::sort(v1->begin(),v1->end());
//       std::sort(v1.begin(),v1.end());
//       // (*v1)(0) = -2.0;
//       // Rcout << *v1 << endl;
//     }else if(isItRowVector(v)){
//       Eigen::Matrix<double, 1, Eigen::Dynamic> *v1=new Eigen::Matrix<double, 1, Eigen::Dynamic>(Rcpp::as<Eigen::Matrix<double, 1, Eigen::Dynamic>>(v));
//       // Eigen::Matrix<double, 1, Eigen::Dynamic> *v1=&(Rcpp::as<Eigen::Matrix<double, 1, Eigen::Dynamic>>(v));
//       // SEXP2EigenMatrix<Eigen::Matrix<double, 1, Eigen::Dynamic>,double>(v,*v1);
//       std::sort(v1->begin(),v1->end());
//     }else{
//       throw std::runtime_error("Failed to sort vector.");
//     }
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdAdd_matAndmat(Eigen::Matrix<double, Eigen::Dynamic, 1> a,Eigen::Matrix<double, Eigen::Dynamic, 1> b){
//   Rcout << "here\n";
//   Rcout << a << endl;
//   Rcout << b << endl;
//   Rcout << a+b << endl;
//   return a+b;
// }
//
// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdMul_sAndmat(Eigen::Matrix<double, Eigen::Dynamic, 1> a,double b){
//   Rcout << a*b << endl;
//   return a*b;
// }
//
// // [[Rcpp::export]]
// Eigen::Matrix<double, Eigen::Dynamic, 1> eigenVectorXdDiv_sAndmat(Eigen::Matrix<double, Eigen::Dynamic, 1> a,double b){
//   // Rcout << a*b << endl;
//   return a/b;
// }

// RCPP_EXPOSED_AS(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>)

RCPP_MODULE(EigenArrayModule){
  class_<Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>>("EigenArrayXXd")
  .factory<SEXP>(newEigenArray<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("cv2r",&eigenArray_cv2r<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("r2cv",&eigenArray_r2cv<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("coeff",&eigenArray2d_coeff<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("coeffRef",&eigenArray2d_coeffRef<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("conservativeResize",&eigenArray2d_conservativeResize<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("setConstant",&eigenArray2d_setConstant<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("setZero",&eigenArray2d_setZero<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("setOnes",&eigenArray2d_setOnes<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("setRandom",&eigenArray2d_setRandom<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("abs",&eigenArray_abs<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("abs2",&eigenArray_abs2<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("acos",&eigenArray_acos<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("acosh",&eigenArray_acosh<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("arg",&eigenArray_arg<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("asin",&eigenArray_asin<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("asinh",&eigenArray_asinh<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("atan",&eigenArray_atan<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("atanh",&eigenArray_atanh<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("ceil",&eigenArray_ceil<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("conj",&eigenArray_conj<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("cos",&eigenArray_cos<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("cosh",&eigenArray_cosh<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("cube",&eigenArray_cube<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("digamma",&eigenArray_digamma<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("erf",&eigenArray_erf<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("erfc",&eigenArray_erfc<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("exp",&eigenArray_exp<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("expm1",&eigenArray_expm1<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("floor",&eigenArray_floor<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("imag",&eigenArray_imag<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("inverse",&eigenArray_inverse<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("isfinite",&eigenArray_isfinite<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("isinf",&eigenArray_isinf<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("isnan",&eigenArray_isnan<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("lgamma",&eigenArray_lgamma<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("log",&eigenArray_log<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("log10",&eigenArray_log10<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("log1p",&eigenArray_log1p<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("log2",&eigenArray_log2<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("logistic",&eigenArray_logistic<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("ndtri",&eigenArray_ndtri<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("real",&eigenArray_real<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("rint",&eigenArray_rint<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("round",&eigenArray_round<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("rsqrt",&eigenArray_rsqrt<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("sign",&eigenArray_sign<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("sin",&eigenArray_sin<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("sinh",&eigenArray_sinh<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("sqrt",&eigenArray_sqrt<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("square",&eigenArray_square<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("tan",&eigenArray_tan<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("tanh",&eigenArray_tanh<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("matrix",&eigenArray_matrix<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("adjoint",&eigenMatrix_adjoint<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("adjointInPlace",&eigenMatrix_adjointInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("maxCoeff_",&eigenMatrix_maxCoeff<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("minCoeff_",&eigenMatrix_minCoeff<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("maxCoeff03_",&eigenMatrix_maxCoeff03<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("applyOnTheLeft",&eigenMatrix_applyOnTheLeft<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("applyOnTheRight",&eigenMatrix_applyOnTheRight<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("transpose",&eigenMatrix_transpose<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("transposeInPlace",&eigenMatrix_transposeInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("OutToConsole",&eigenArray_OutToConsole<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("resize",&eigenArray_resize<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("at", &eigenArray_at<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("set", &eigenArray_set<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("rows",&eigenArray_rows<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("cols",&eigenArray_cols<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("size",&eigenArray_size<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("col",&eigenArray_col<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("row",&eigenArray_row<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("nonZeros",&eigenMatrix_nonZeros<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("outerSize",&eigenMatrix_outerSize<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("replicate",&eigenMatrix_replicate<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("reverse",&eigenMatrix_reverse<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("reverseInPlace",&eigenMatrix_reverseInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("colwise",&eigenMatrix_colwise<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("rowSums",&eigenMatrix_rowwise_Sum<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("mean",&eigenMatrix_mean<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("prod",&eigenMatrix_prod<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("sum",&eigenMatrix_sum<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("trace",&eigenMatrix_trace<double, Eigen::Dynamic, Eigen::Dynamic>)
        // .method("greaterthan",&eigenArray_greaterthan<double, Eigen::Dynamic, Eigen::Dynamic>)
  ;

  class_<Eigen::Array<double, Eigen::Dynamic, 1>>("EigenArrayXd")
    .factory<SEXP>(newEigenArray1d<double, Eigen::Dynamic>)
    .method("cv2r",&eigenArray_cv2r<double, Eigen::Dynamic, 1>)
    .method("r2cv",&eigenArray_r2cv<double, Eigen::Dynamic, 1>)
    .method("coeff",&eigenArray2d_coeff<double, Eigen::Dynamic, 1>)
    .method("coeffRef",&eigenArray2d_coeffRef<double, Eigen::Dynamic, 1>)
    .method("conservativeResize",&eigenArray2d_conservativeResize<double, Eigen::Dynamic, 1>)
    .method("setConstant",&eigenArray2d_setConstant<double, Eigen::Dynamic, 1>)
    .method("setZero",&eigenArray2d_setZero<double, Eigen::Dynamic, 1>)
    .method("setOnes",&eigenArray2d_setOnes<double, Eigen::Dynamic, 1>)
    .method("setRandom",&eigenArray2d_setRandom<double, Eigen::Dynamic, 1>)
    .method("abs",&eigenArray_abs<double, Eigen::Dynamic, 1>)
    .method("abs2",&eigenArray_abs2<double, Eigen::Dynamic, 1>)
    .method("acos",&eigenArray_acos<double, Eigen::Dynamic, 1>)
    .method("acosh",&eigenArray_acosh<double, Eigen::Dynamic, 1>)
    .method("arg",&eigenArray_arg<double, Eigen::Dynamic, 1>)
    .method("asin",&eigenArray_asin<double, Eigen::Dynamic, 1>)
    .method("asinh",&eigenArray_asinh<double, Eigen::Dynamic, 1>)
    .method("atan",&eigenArray_atan<double, Eigen::Dynamic, 1>)
    .method("atanh",&eigenArray_atanh<double, Eigen::Dynamic, 1>)
    .method("ceil",&eigenArray_ceil<double, Eigen::Dynamic, 1>)
  // .method("conj",&eigenArray_conj<double, Eigen::Dynamic, 1>)
     .method("cos",&eigenArray_cos<double, Eigen::Dynamic, 1>)
     .method("cosh",&eigenArray_cosh<double, Eigen::Dynamic, 1>)
     .method("cube",&eigenArray_cube<double, Eigen::Dynamic, 1>)
  // .method("digamma",&eigenArray_digamma<double, Eigen::Dynamic, 1>)
  // .method("erf",&eigenArray_erf<double, Eigen::Dynamic, 1>)
  // .method("erfc",&eigenArray_erfc<double, Eigen::Dynamic, 1>)
     .method("exp",&eigenArray_exp<double, Eigen::Dynamic, 1>)
     .method("expm1",&eigenArray_expm1<double, Eigen::Dynamic, 1>)
     .method("floor",&eigenArray_floor<double, Eigen::Dynamic, 1>)
  // .method("imag",&eigenArray_imag<double, Eigen::Dynamic, 1>)
     .method("inverse",&eigenArray_inverse<double, Eigen::Dynamic, 1>)
  // .method("isfinite",&eigenArray_isfinite<double, Eigen::Dynamic, 1>)
  // .method("isinf",&eigenArray_isinf<double, Eigen::Dynamic, 1>)
  // .method("isnan",&eigenArray_isnan<double, Eigen::Dynamic, 1>)
  // .method("lgamma",&eigenArray_lgamma<double, Eigen::Dynamic, 1>)
     .method("log",&eigenArray_log<double, Eigen::Dynamic, 1>)
     .method("log10",&eigenArray_log10<double, Eigen::Dynamic, 1>)
     .method("log1p",&eigenArray_log1p<double, Eigen::Dynamic, 1>)
     .method("log2",&eigenArray_log2<double, Eigen::Dynamic, 1>)
     .method("logistic",&eigenArray_logistic<double, Eigen::Dynamic, 1>)
  // .method("ndtri",&eigenArray_ndtri<double, Eigen::Dynamic, 1>)
     .method("real",&eigenArray_real<double, Eigen::Dynamic, 1>)
     .method("rint",&eigenArray_rint<double, Eigen::Dynamic, 1>)
     .method("round",&eigenArray_round<double, Eigen::Dynamic, 1>)
     .method("rsqrt",&eigenArray_rsqrt<double, Eigen::Dynamic, 1>)
     .method("sign",&eigenArray_sign<double, Eigen::Dynamic, 1>)
     .method("sin",&eigenArray_sin<double, Eigen::Dynamic, 1>)
     .method("sinh",&eigenArray_sinh<double, Eigen::Dynamic, 1>)
     .method("sqrt",&eigenArray_sqrt<double, Eigen::Dynamic, 1>)
     .method("square",&eigenArray_square<double, Eigen::Dynamic, 1>)
     .method("tan",&eigenArray_tan<double, Eigen::Dynamic, 1>)
     .method("tanh",&eigenArray_tanh<double, Eigen::Dynamic, 1>)
     .method("matrix",&eigenArray_matrix<double, Eigen::Dynamic, 1>)
     .method("head",&eigenArray1d_head<double, Eigen::Dynamic, 1>)
     .method("tail",&eigenArray1d_tail<double, Eigen::Dynamic, 1>)
     .method("segment",&eigenArray1d_segment<double, Eigen::Dynamic, 1>)
  // .method("adjoint",&eigenMatrix_adjoint<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("adjointInPlace",&eigenMatrix_adjointInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("maxCoeff_",&eigenMatrix_maxCoeff<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("minCoeff_",&eigenMatrix_minCoeff<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("maxCoeff03_",&eigenMatrix_maxCoeff03<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("applyOnTheLeft",&eigenMatrix_applyOnTheLeft<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("applyOnTheRight",&eigenMatrix_applyOnTheRight<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("transpose",&eigenMatrix_transpose<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("transposeInPlace",&eigenMatrix_transposeInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
        .method("OutToConsole",&eigenArray_OutToConsole<double, Eigen::Dynamic, 1>)
        .method("resize",&eigenArray_resize<double, Eigen::Dynamic, 1>)
        .method("at", &eigenArray_at<double, Eigen::Dynamic, 1>)
        .method("set", &eigenArray_set<double, Eigen::Dynamic, 1>)
        .method("rows",&eigenArray_rows<double, Eigen::Dynamic, 1>)
        .method("cols",&eigenArray_cols<double, Eigen::Dynamic, 1>)
        .method("size",&eigenArray_size<double, Eigen::Dynamic, 1>)
        .method("sort",&eigenArray_sort<double, Eigen::Dynamic, 1>)
  //    .method("nonZeros",&eigenMatrix_nonZeros<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("outerSize",&eigenMatrix_outerSize<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("replicate",&eigenMatrix_replicate<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("reverse",&eigenMatrix_reverse<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("reverseInPlace",&eigenMatrix_reverseInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("colwise",&eigenMatrix_colwise<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("rowSums",&eigenMatrix_rowwise_Sum<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("mean",&eigenMatrix_mean<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("prod",&eigenMatrix_prod<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("sum",&eigenMatrix_sum<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("trace",&eigenMatrix_trace<double, Eigen::Dynamic, Eigen::Dynamic>)
  ;

}

RCPP_MODULE(EigenMatrixModule){
  // class_<Eigen::Matrix<double, 4, 4>>("EigenMatrix4d")
  // .factory<SEXP>(newEigenMatrix<double, 4, 4>)
  // .method("cv2r",&eigenMatrix_cv2r<double, 4, 4>)
  // .method("transpose",&eigenMatrix_transpose<double, 4, 4>)
  // ;

  class_<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>>("EigenMatrixXi")
  .factory<SEXP>(newEigenMatrix<int, Eigen::Dynamic, Eigen::Dynamic>)
  .method("cv2r",&eigenMatrix_cv2r<int, Eigen::Dynamic, Eigen::Dynamic>)
  .method("transpose",&eigenMatrix_transpose<int, Eigen::Dynamic, Eigen::Dynamic>)
  ;

  class_<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>("EigenMatrixXd")
  .factory<SEXP>(newEigenMatrix<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("cv2r",&eigenMatrix_cv2r<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("r2cv",&eigenMatrix_r2cv<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("coeff",&eigenMatrix_coeff<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("coeffRef",&eigenMatrix_coeffRef<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("conservativeResize",&eigenMatrix_conservativeResize<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("setConstant",&eigenMatrix_setConstant<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("setZero",&eigenMatrix_setZero<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("setOnes",&eigenMatrix_setOnes<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("setRandom",&eigenMatrix_setRandom<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("setIdentity",&eigenMatrix_setIdentity<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("acosh",&eigenMatrix_acosh<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("exp",&eigenMatrix_exp<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("adjoint",&eigenMatrix_adjoint<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("adjointInPlace",&eigenMatrix_adjointInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("applyOnTheLeft_",&eigenMatrix_applyOnTheLeft<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("applyOnTheRight_",&eigenMatrix_applyOnTheRight<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("array",&eigenMatrix_array<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("asDiagonal",&eigenMatrix_asDiagonal<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("asinh",&eigenMatrix_asinh<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("atanh",&eigenMatrix_atanh<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("bdcSvd",&eigenMatrix_bdcSvd<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("blueNorm",&eigenMatrix_blueNorm<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("colPivHouseholderQr",&eigenMatrix_colPivHouseholderQr<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("completeOrthogonalDecomposition",&eigenMatrix_completeOrthogonalDecomposition<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("computeInverseAndDetWithCheck_",&eigenMatrix_computeInverseAndDetWithCheck<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("computeInverseWithCheck_",&eigenMatrix_computeInverseWithCheck<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("cos",&eigenMatrix_cos<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("cosh",&eigenMatrix_cosh<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("determinant",&eigenMatrix_determinant<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("eigenvalues",&eigenMatrix_eigenvalues<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("hnormalized",&eigenMatrix_hnormalized<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("homogeneous",&eigenMatrix_homogeneous<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("fullPivLu",&eigenMatrix_fullPivLu<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("householderQr",&eigenMatrix_householderQr<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("ldlt",&eigenMatrix_ldlt<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("hypotNorm",&eigenMatrix_hypotNorm<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("lpNorm1",&eigenMatrix_lpNorm<double, Eigen::Dynamic, Eigen::Dynamic,1>)
  .method("lpNorm2",&eigenMatrix_lpNorm<double, Eigen::Dynamic, Eigen::Dynamic,2>)
  .method("lpNormInf",&eigenMatrix_lpNorm<double, Eigen::Dynamic, Eigen::Dynamic,Eigen::Infinity>)
  .method("norm",&eigenMatrix_norm<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("inverse",&eigenMatrix_inverse<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("maxCoeff_",&eigenMatrix_maxCoeff<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("minCoeff_",&eigenMatrix_minCoeff<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("maxCoeff03_",&eigenMatrix_maxCoeff03<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("applyOnTheLeft",&eigenMatrix_applyOnTheLeft<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("applyOnTheRight",&eigenMatrix_applyOnTheRight<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("transpose",&eigenMatrix_transpose<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("transposeInPlace",&eigenMatrix_transposeInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("OutToConsole",&eigenMatrix_OutToConsole<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("resize",&eigenMatrix_resize<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("at", &eigenMatrix_at<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("set", &eigenMatrix_set<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("reshaped_",&eigenMatrix_reshaped<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("rows",&eigenMatrix_rows<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("cols",&eigenMatrix_cols<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("size",&eigenMatrix_size<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("nonZeros",&eigenMatrix_nonZeros<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("outerSize",&eigenMatrix_outerSize<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("replicate",&eigenMatrix_replicate<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("reverse",&eigenMatrix_reverse<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("reverseInPlace",&eigenMatrix_reverseInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("colwise",&eigenMatrixXd_colwise)
  .method("colAdd",&eigenMatrix_colwise_Add<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("rowAdd",&eigenMatrix_rowwise_Add<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("colSums",&eigenMatrix_colwise_Sum<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("rowSums",&eigenMatrix_rowwise_Sum<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("mean",&eigenMatrix_mean<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("prod",&eigenMatrix_prod<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("sum",&eigenMatrix_sum<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("squaredNorm",&eigenMatrix_squaredNorm<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("trace",&eigenMatrix_trace<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("block",&eigenMatrix_block<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("col",&eigenMatrix_col<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("row",&eigenMatrix_row<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("assign_",&eigenMatrix_assign<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("topLeftCorner",&eigenMatrix_topLeftCorner<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("bottomLeftCorner",&eigenMatrix_bottomLeftCorner<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("topRightCorner",&eigenMatrix_topRightCorner<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("bottomRightCorner",&eigenMatrix_bottomRightCorner<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("topRows",&eigenMatrix_topRows<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("bottomRows",&eigenMatrix_bottomRows<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("leftCols",&eigenMatrix_leftCols<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("rightCols",&eigenMatrix_rightCols<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("middleCols",&eigenMatrix_middleCols<double, Eigen::Dynamic, Eigen::Dynamic>)
  .method("middleRows",&eigenMatrix_middleRows<double, Eigen::Dynamic, Eigen::Dynamic>)

  ;

  class_<Eigen::Matrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>>("EigenMatrixXcf")
    .factory<SEXP>(newEigenMatrix<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
    .method("cv2r",&eigenMatrix_cv2r<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
    .method("r2cv",&eigenMatrix_r2cv<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("coeff",&eigenMatrix_coeff<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("coeffRef",&eigenMatrix_coeffRef<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("conservativeResize",&eigenMatrix_conservativeResize<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("setConstant",&eigenMatrix_setConstant<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("setZero",&eigenMatrix_setZero<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("setOnes",&eigenMatrix_setOnes<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("setRandom",&eigenMatrix_setRandom<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("acosh",&eigenMatrix_acosh<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("exp",&eigenMatrix_exp<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
    .method("adjoint",&eigenMatrix_adjoint<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
    .method("adjointInPlace",&eigenMatrix_adjointInPlace<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
    .method("inverse",&eigenMatrix_inverse<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("applyOnTheLeft",&eigenMatrix_applyOnTheLeft<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("applyOnTheRight",&eigenMatrix_applyOnTheRight<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("transpose",&eigenMatrix_transpose<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
     .method("transposeInPlace",&eigenMatrix_transposeInPlace<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
     .method("array",&eigenMatrix_array<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
     .method("OutToConsole",&eigenMatrix_OutToConsole<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
     .method("resize",&eigenMatrix_resize<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
     .method("at", &eigenMatrix_at<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
     .method("set", &eigenMatrixXc_set<float, Eigen::Dynamic, Eigen::Dynamic>)
     .method("conjugate", &eigenMatrixXc_conjugate<float, Eigen::Dynamic, Eigen::Dynamic>)
     .method("rows",&eigenMatrix_rows<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
     .method("cols",&eigenMatrix_cols<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
     .method("size",&eigenMatrix_size<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("nonZeros",&eigenMatrix_nonZeros<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("outerSize",&eigenMatrix_outerSize<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("replicate",&eigenMatrix_replicate<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("reverse",&eigenMatrix_reverse<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("reverseInPlace",&eigenMatrix_reverseInPlace<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("colwise",&eigenMatrix_colwise<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("rowSums",&eigenMatrix_rowwise_Sum<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("sum",&eigenMatrix_sum<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic>)
  ;

  class_<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>>("EigenMatrixXcd")
    .factory<SEXP>(newEigenMatrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
    .method("cv2r",&eigenMatrix_cv2r<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
    .method("r2cv",&eigenMatrix_r2cv<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("coeff",&eigenMatrix_coeff<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("coeffRef",&eigenMatrix_coeffRef<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("conservativeResize",&eigenMatrix_conservativeResize<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("setConstant",&eigenMatrix_setConstant<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("setZero",&eigenMatrix_setZero<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("setOnes",&eigenMatrix_setOnes<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("setRandom",&eigenMatrix_setRandom<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("acosh",&eigenMatrix_acosh<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //   .method("exp",&eigenMatrix_exp<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
       .method("adjoint",&eigenMatrix_adjoint<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
       .method("adjointInPlace",&eigenMatrix_adjointInPlace<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
       .method("householderQr",&eigenMatrix_householderQr<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
       .method("inverse",&eigenMatrix_inverse<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("applyOnTheLeft",&eigenMatrix_applyOnTheLeft<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("applyOnTheRight",&eigenMatrix_applyOnTheRight<double, Eigen::Dynamic, Eigen::Dynamic>)
        .method("transpose",&eigenMatrix_transpose<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("transposeInPlace",&eigenMatrix_transposeInPlace<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("array",&eigenMatrix_array<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("OutToConsole",&eigenMatrix_OutToConsole<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("resize",&eigenMatrix_resize<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("at", &eigenMatrix_at<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("set", &eigenMatrixXc_set<double, Eigen::Dynamic, Eigen::Dynamic>)
        .method("conjugate", &eigenMatrixXc_conjugate<double, Eigen::Dynamic, Eigen::Dynamic>)
        .method("rows",&eigenMatrix_rows<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("cols",&eigenMatrix_cols<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("size",&eigenMatrix_size<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("nonZeros",&eigenMatrix_nonZeros<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("outerSize",&eigenMatrix_outerSize<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("replicate",&eigenMatrix_replicate<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("reverse",&eigenMatrix_reverse<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("reverseInPlace",&eigenMatrix_reverseInPlace<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("colwise",&eigenMatrix_colwise<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("rowSums",&eigenMatrix_rowwise_Sum<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("sum",&eigenMatrix_sum<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("col",&eigenMatrix_col<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("row",&eigenMatrix_row<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("real",&eigenMatrix_real<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        .method("imag",&eigenMatrix_imag<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)

  ;
}

RCPP_MODULE(EigenVectorModule){
  class_<Eigen::Matrix<std::complex<double>, Eigen::Dynamic,1>>("EigenVectorXcd")
  .factory<SEXP>(newEigenVector<std::complex<double>, Eigen::Dynamic>)
  .method("cv2r",&eigenMatrix_cv2r<std::complex<double>, Eigen::Dynamic, 1>)
  .method("r2cv",&eigenMatrix_r2cv<std::complex<double>, Eigen::Dynamic, 1>)
  .method("coeff",&eigenVector_coeff<std::complex<double>, Eigen::Dynamic, 1>)
  .method("coeffRef",&eigenVector_coeffRef<std::complex<double>, Eigen::Dynamic, 1>)
  .method("conservativeResize",&eigenVector_conservativeResize<std::complex<double>, Eigen::Dynamic, 1>)
  .method("setConstant",&eigenVector_setConstant<std::complex<double>, Eigen::Dynamic, 1>)
  .method("setZero",&eigenVector_setZero<std::complex<double>, Eigen::Dynamic, 1>)
  .method("setOnes",&eigenVector_setOnes<std::complex<double>, Eigen::Dynamic, 1>)
  .method("setRandom",&eigenVector_setRandom<std::complex<double>, Eigen::Dynamic, 1>)
  // .method("acosh",&eigenMatrix_acosh<std::complex<double>, Eigen::Dynamic, 1>)
  .method("adjoint",&eigenVector_adjoint<std::complex<double>, Eigen::Dynamic>)
  .method("array",&eigenMatrix_array<std::complex<double>, Eigen::Dynamic, 1>)
  .method("asDiagonal",&eigenMatrixCPL_asDiagonal<double, Eigen::Dynamic>)
  // .method("cross",&eigenVector_cross<std::complex<double>, Eigen::Dynamic>)
     // .method("maxCoeff_",&eigenMatrix_maxCoeff1<std::complex<double>, Eigen::Dynamic, 1>)
     // .method("minCoeff_",&eigenMatrix_minCoeff1<std::complex<double>, Eigen::Dynamic, 1>)
     .method("dot",&eigenVector_dot<std::complex<double>, Eigen::Dynamic>)
     .method("transpose",&eigenMatrix_transpose<std::complex<double>, Eigen::Dynamic,1>)
  // .method("transposeInPlace",&eigenMatrix_transposeInPlace<std::complex<double>, Eigen::Dynamic,1>)
     .method("OutToConsole",&eigenMatrix_OutToConsole<std::complex<double>, Eigen::Dynamic, 1>)
     .method("resize",&eigenVector_resize<std::complex<double>, Eigen::Dynamic, 1>)
     .method("at", &eigenMatrix_at<std::complex<double>, Eigen::Dynamic, 1>)
     .method("set", &eigenMatrix_set<std::complex<double>, Eigen::Dynamic, 1>)
     .method("rows",&eigenMatrix_rows<std::complex<double>, Eigen::Dynamic, 1>)
     .method("cols",&eigenMatrix_cols<std::complex<double>, Eigen::Dynamic, 1>)
     .method("size",&eigenMatrix_size<std::complex<double>, Eigen::Dynamic, 1>)
  // .method("nonZeros",&eigenMatrix_nonZeros<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("outerSize",&eigenMatrix_outerSize<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("replicate",&eigenMatrix_replicate<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("reverse",&eigenMatrix_reverse<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("reverseInPlace",&eigenMatrix_reverseInPlace<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("colwise",&eigenMatrix_colwise<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("rowSums",&eigenMatrix_rowwise_Sum<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>)
        // .method("mean",&eigenMatrix_sum<std::complex<double>, Eigen::Dynamic, 1>)
        .method("prod",&eigenMatrix_prod<std::complex<double>, Eigen::Dynamic, 1>)
        .method("sum",&eigenMatrix_sum<std::complex<double>, Eigen::Dynamic, 1>)
        .method("real",&eigenMatrix_real<std::complex<double>, Eigen::Dynamic, 1>)
        .method("imag",&eigenMatrix_imag<std::complex<double>, Eigen::Dynamic, 1>)
  ;

  class_<Eigen::Matrix<double, 3,1>>("EigenVector3d")
    .factory<SEXP>(newEigenVector<double, 3>)
    .method("cv2r",&eigenMatrix_cv2r<double, 3, 1>)
    .method("r2cv",&eigenMatrix_r2cv<double, 3, 1>)
    .method("coeff",&eigenVector_coeff<double, 3, 1>)
    .method("coeffRef",&eigenVector_coeffRef<double, 3, 1>)
    // .method("conservativeResize",&eigenVector_conservativeResize<double, 3, 1>)
    .method("setConstant",&eigenVector_setConstant<double, 3, 1>)
    .method("setZero",&eigenVector_setZero<double, 3, 1>)
    .method("setOnes",&eigenVector_setOnes<double, 3, 1>)
    .method("setRandom",&eigenVector_setRandom<double, 3, 1>)
    .method("acosh",&eigenMatrix_acosh<double, 3, 1>)
    .method("adjoint",&eigenVector_adjoint<double, 3>)
    .method("applyOnTheLeft_",&eigenMatrix_applyOnTheLeft<double, 3, 1>)
    .method("applyOnTheRight_",&eigenMatrix_applyOnTheRight<double, 3, 1>)
    .method("array",&eigenMatrix_array<double, 3, 1>)
    .method("asDiagonal",&eigenMatrix_asDiagonal<double, 3, 1>)
    .method("norm",&eigenMatrix_norm<double, 3, 1>)
    .method("squaredNorm",&eigenMatrix_squaredNorm<double, 3, 1>)
    .method("maxCoeff_",&eigenMatrix_maxCoeff1<double, 3, 1>)
    .method("minCoeff_",&eigenMatrix_minCoeff1<double, 3, 1>)
    .method("cross",&eigenVector_cross<double, 3>)
    .method("dot",&eigenVector_dot<double, 3>)
    .method("transpose",&eigenMatrix_transpose<double, 3,1>)
    .method("OutToConsole",&eigenMatrix_OutToConsole<double, 3, 1>)
    // .method("resize",&eigenVector_resize<double, 3, 1>)
    .method("at", &eigenMatrix_at<double, 3, 1>)
    .method("set", &eigenMatrix_set<double, 3, 1>)
    .method("rows",&eigenMatrix_rows<double, 3, 1>)
    .method("cols",&eigenMatrix_cols<double, 3, 1>)
    .method("size",&eigenMatrix_size<double, 3, 1>)
    .method("mean",&eigenMatrix_sum<double, 3, 1>)
    .method("prod",&eigenMatrix_prod<double, 3, 1>)
    .method("sum",&eigenMatrix_sum<double, 3, 1>)
    .method("sort",&eigenMatrix_sort<double, 3, 1>)
  ;

  class_<Eigen::Matrix<double, Eigen::Dynamic,1>>("EigenVectorXd")
  .factory<SEXP>(newEigenVector<double, Eigen::Dynamic>)
  .method("cv2r",&eigenMatrix_cv2r<double, Eigen::Dynamic, 1>)
  .method("r2cv",&eigenMatrix_r2cv<double, Eigen::Dynamic, 1>)
  .method("coeff",&eigenVector_coeff<double, Eigen::Dynamic, 1>)
  .method("coeffRef",&eigenVector_coeffRef<double, Eigen::Dynamic, 1>)
  .method("conservativeResize",&eigenVector_conservativeResize<double, Eigen::Dynamic, 1>)
  .method("setConstant",&eigenVector_setConstant<double, Eigen::Dynamic, 1>)
  .method("setZero",&eigenVector_setZero<double, Eigen::Dynamic, 1>)
  .method("setOnes",&eigenVector_setOnes<double, Eigen::Dynamic, 1>)
  .method("setRandom",&eigenVector_setRandom<double, Eigen::Dynamic, 1>)
  .method("acosh",&eigenMatrix_acosh<double, Eigen::Dynamic, 1>)
  .method("adjoint",&eigenVector_adjoint<double, Eigen::Dynamic>)
  .method("applyOnTheLeft_",&eigenMatrix_applyOnTheLeft<double, Eigen::Dynamic, 1>)
  .method("applyOnTheRight_",&eigenMatrix_applyOnTheRight<double, Eigen::Dynamic, 1>)
  .method("array",&eigenMatrix_array<double, Eigen::Dynamic, 1>)
  .method("asDiagonal",&eigenMatrix_asDiagonal<double, Eigen::Dynamic, 1>)
  .method("norm",&eigenMatrix_norm<double, Eigen::Dynamic, 1>)
  .method("squaredNorm",&eigenMatrix_squaredNorm<double, Eigen::Dynamic, 1>)
  // .method("cross",&eigenVector_cross<double, Eigen::Dynamic>)
     .method("maxCoeff_",&eigenMatrix_maxCoeff1<double, Eigen::Dynamic, 1>)
     .method("minCoeff_",&eigenMatrix_minCoeff1<double, Eigen::Dynamic, 1>)
     // .method("cross",&eigenVector_cross<double, 3>)
  .method("dot",&eigenVector_dot<double, Eigen::Dynamic>)
  .method("transpose",&eigenMatrix_transpose<double, Eigen::Dynamic,1>)
  // .method("transposeInPlace",&eigenMatrix_transposeInPlace<double, Eigen::Dynamic,1>)
  .method("OutToConsole",&eigenMatrix_OutToConsole<double, Eigen::Dynamic, 1>)
  .method("resize",&eigenVector_resize<double, Eigen::Dynamic, 1>)
  .method("at", &eigenMatrix_at<double, Eigen::Dynamic, 1>)
  .method("set", &eigenMatrix_set<double, Eigen::Dynamic, 1>)
  .method("rows",&eigenMatrix_rows<double, Eigen::Dynamic, 1>)
  .method("cols",&eigenMatrix_cols<double, Eigen::Dynamic, 1>)
  .method("size",&eigenMatrix_size<double, Eigen::Dynamic, 1>)
  // .method("nonZeros",&eigenMatrix_nonZeros<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("outerSize",&eigenMatrix_outerSize<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("replicate",&eigenMatrix_replicate<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("reverse",&eigenMatrix_reverse<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("reverseInPlace",&eigenMatrix_reverseInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("colwise",&eigenMatrix_colwise<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("rowSums",&eigenMatrix_rowwise_Sum<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("mean",&eigenMatrix_sum<double, Eigen::Dynamic, 1>)
     .method("prod",&eigenMatrix_prod<double, Eigen::Dynamic, 1>)
     .method("sum",&eigenMatrix_sum<double, Eigen::Dynamic, 1>)
     .method("sort",&eigenMatrix_sort<double, Eigen::Dynamic, 1>)
  ;

  class_<Eigen::Matrix<double, 1, Eigen::Dynamic>>("EigenRowVectorXd")
    .factory<SEXP>(newEigenRowVector<double, Eigen::Dynamic>)
    .method("cv2r",&eigenMatrix_cv2r<double, 1, Eigen::Dynamic>)
    .method("r2cv",&eigenMatrix_r2cv<double, 1, Eigen::Dynamic>)
    .method("coeff",&eigenVector_coeff<double, 1, Eigen::Dynamic>)
    .method("coeffRef",&eigenVector_coeffRef<double, 1, Eigen::Dynamic>)
    .method("conservativeResize",&eigenVector_conservativeResize<double, 1, Eigen::Dynamic>)
    .method("setConstant",&eigenVector_setConstant<double, 1, Eigen::Dynamic>)
    .method("setZero",&eigenVector_setZero<double, 1, Eigen::Dynamic>)
    .method("setOnes",&eigenVector_setOnes<double, 1, Eigen::Dynamic>)
    .method("setRandom",&eigenVector_setRandom<double, 1, Eigen::Dynamic>)
    .method("acosh",&eigenMatrix_acosh<double, 1, Eigen::Dynamic>)
    .method("adjoint",&eigenRowVector_adjoint<double, Eigen::Dynamic>)
    .method("array",&eigenMatrix_array<double, 1, Eigen::Dynamic>)
    // .method("cross",&eigenRowVector_cross<double, Eigen::Dynamic>)
       .method("maxCoeff_",&eigenMatrix_maxCoeff1<double, 1, Eigen::Dynamic>)
       .method("minCoeff_",&eigenMatrix_minCoeff1<double, 1, Eigen::Dynamic>)
    .method("dot",&eigenRowVector_dot<double, Eigen::Dynamic>)
    .method("transpose",&eigenMatrix_transpose<double, 1, Eigen::Dynamic>)
    // .method("transposeInPlace",&eigenMatrix_transposeInPlace<double, 1, Eigen::Dynamic>)
    .method("OutToConsole",&eigenMatrix_OutToConsole<double, 1, Eigen::Dynamic>)
    .method("resize",&eigenVector_resize<double, 1, Eigen::Dynamic>)
    .method("at", &eigenMatrix_at<double, 1, Eigen::Dynamic>)
    .method("set", &eigenMatrix_set<double, 1, Eigen::Dynamic>)
    .method("rows",&eigenMatrix_rows<double, 1, Eigen::Dynamic>)
    .method("cols",&eigenMatrix_cols<double, 1, Eigen::Dynamic>)
    .method("size",&eigenMatrix_size<double, 1, Eigen::Dynamic>)
  // .method("nonZeros",&eigenMatrix_nonZeros<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("outerSize",&eigenMatrix_outerSize<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("replicate",&eigenMatrix_replicate<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("reverse",&eigenMatrix_reverse<double, Eigen::Dynamic, Eigen::Dynamic>)
  // .method("reverseInPlace",&eigenMatrix_reverseInPlace<double, Eigen::Dynamic, Eigen::Dynamic>)
  // // .method("colwise",&eigenMatrix_colwise<double, Eigen::Dynamic, Eigen::Dynamic>)
  //    .method("rowSums",&eigenMatrix_rowwise_Sum<double, Eigen::Dynamic, Eigen::Dynamic>)
     .method("mean",&eigenMatrix_sum<double, 1, Eigen::Dynamic>)
     .method("prod",&eigenMatrix_prod<double, 1, Eigen::Dynamic>)
     .method("sum",&eigenMatrix_sum<double, 1, Eigen::Dynamic>)
     .method("sort",&eigenMatrix_sort<double, 1, Eigen::Dynamic>)
  ;
}


// RCPP_EXPOSED_AS(EigenVectorXd)


// RCPP_MODULE(EigenVectorwiseOpModule){
//   class_<Eigen::VectorwiseOp<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Horizontal>>("EigenVectorRowwiseOp4MatrixXd")
//   .method("sum",&eigenVectorwiseOp_sum<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::Horizontal>)
//   ;
// }
