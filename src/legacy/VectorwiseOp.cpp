#include "ROpenEigen_types.h"

// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> myMatrixXd;
// typedef Eigen::VectorwiseOp<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>,Eigen::Vertical> myVectorwiseOpCol4MatrixXd;
// // typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> myMatrixXcd;
// // typedef Eigen::HouseholderQR<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>> myHouseholderQR4MatrixXcd;
//
//
// template<typename ExpressionType, int Direction>
// ExpressionType VectorwiseOp_ExpressionType(Eigen::VectorwiseOp<ExpressionType,Direction>* obj){
//   try{
//     return obj->_expression();
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }
//
// RCPP_MODULE(VectorwiseOpModule){
//   class_<myVectorwiseOpCol4MatrixXd>("VectorwiseOpCol4MatrixXd")
//   // .default_constructor()
//   .constructor<myMatrixXd>()
//   // .constructor<int,int>()
//   .method("ExpressionType",&VectorwiseOp_ExpressionType<myMatrixXd,Eigen::Vertical>)
//   // .method("matrixQR",&eigenHouseholderQR_matrixQR<myMatrixXd>)
//      // .method("matrixQR",&eigenHouseholderQR_matrixQRUpper<myMatrixXd>)
//   ;
//
// }
