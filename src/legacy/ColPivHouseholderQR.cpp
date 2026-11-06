// #include "ROpenEigen_types.h"
//
// template<typename _MatrixType>
// SEXP eigenColPivHouseholderQR_Solve(Eigen::ColPivHouseholderQR<_MatrixType>* obj,SEXP b){
//   try{
//     // typedef _MatrixType::RealScalar Scalar;
//     typedef typename _MatrixType::RealScalar Scalar;
//     SEXP result;
//     if(isItVector(b)){
//       Eigen::Matrix<Scalar, Eigen::Dynamic, 1> _b=Rcpp::as<Eigen::Matrix<Scalar, Eigen::Dynamic, 1>>(b);
//       Eigen::ColPivHouseholderQR<_MatrixType> decomposer(*obj);
//       // Eigen::Solve<Eigen::ColPivHouseholderQR, Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > res = obj->solve(_b);
//       Eigen::Matrix<Scalar, Eigen::Dynamic, 1> res = decomposer.solve(_b);
//       return Rcpp::wrap(res);
//     }
//     // return result;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }
//
// typedef Eigen::ColPivHouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> myColPivHouseholderQR4MatrixXd;
// typedef Eigen::ColPivHouseholderQR<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> myColPivHouseholderQR4MatrixXcd;
//
//
//
//
// RCPP_MODULE(EigenColPivHouseholderQRModule){
//   class_<myColPivHouseholderQR4MatrixXd>("ColPivHouseholderQR4MatrixXd")
//   .default_constructor()
//   .constructor<myMatrixXd>()
//   .constructor<int,int>()
//   .method("solve",&eigenColPivHouseholderQR_Solve<myMatrixXd>)
//   // .method("matrixQR",&eigenHouseholderQR_matrixQR<myMatrixXd>)
//      // .method("matrixQR",&eigenHouseholderQR_matrixQRUpper<myMatrixXd>)
//   ;
//
//   // class_<myHouseholderQR4MatrixXcd>("HouseholderQR4MatrixXcd")
//   //   .default_constructor()
//   //   .constructor<myMatrixXcd>()
//   //   .constructor<int,int>()
//   //   .method("householderQ",&eigenHouseholderQR_householderQ<myMatrixXcd>)
//   // // .method("matrixQR",&eigenHouseholderQR_matrixQR<myMatrixXd>)
//   //    .method("matrixQR",&eigenHouseholderQR_matrixQRUpper<myMatrixXcd>)
//   // ;
// }
