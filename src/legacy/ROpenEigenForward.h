// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppEigenForward.h: Rcpp/Eigen glue
//
// Copyright (C)      2011 Douglas Bates, Dirk Eddelbuettel and Romain Francois
//
// This file is part of RcppEigen.
//
// RcppEigen is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppEigen is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppEigen.  If not, see <http://www.gnu.org/licenses/>.

#ifndef ROpenEigen__ROpenEigenForward__h
#define ROpenEigen__ROpenEigenForward__h

#include <iterator>
// #include <RcppCommon.h>
#include <Rconfig.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Cholesky>
#include <Eigen/SVD>
#include <Eigen/LU>
#include <Eigen/Geometry>
#include <Eigen/Jacobi>
//#include <unsupported/Eigen/AutoDiff>  // causes problems redefining sign
#include <unsupported/Eigen/IterativeSolvers>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <unsupported/Eigen/Polynomials>
#include <unsupported/Eigen/SparseExtra> // also includes Eigen/Sparse
#include <unsupported/Eigen/Splines>

// namespace Eigen {
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
// }

//自定义派生类会出现很多internal traits访问错误

/* forward declarations */
namespace Rcpp {

template<typename _Scalar, int _Rows, int _Cols> SEXP wrap(const Eigen::Matrix<_Scalar, _Rows, _Cols>&);
// template<typename _Scalar, int _Rows> SEXP wrap(const Eigen::Matrix<_Scalar, _Rows, 1>&);
// template<typename _Scalar, int _Cols> SEXP wrap(const Eigen::Matrix<_Scalar, 1, _Cols>&);
template<typename _Scalar, int _Rows, int _Cols> SEXP wrap(const Eigen::Array<_Scalar, _Rows, _Cols>&);
// template<typename _MatrixType> SEXP wrap(const Eigen::HouseholderQR<_MatrixType>&);
// template<typename _Scalar> SEXP wrap(const Eigen::HouseholderQR<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>&);
// template<typename _Scalar, int _Rows, int _Cols> SEXP wrap(const Eigen::HouseholderQR<Eigen::Matrix<_Scalar, _Rows, _Cols>>&);
// SEXP wrap(const Eigen::HouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>&);
// template<typename _MatrixType> SEXP wrap(const Eigen::EigenHouseholderQR<_MatrixType>&);
template<typename _MatrixType> SEXP wrap(const Eigen::HouseholderQR<_MatrixType>&);
template<typename _MatrixType> SEXP wrap(const Eigen::ColPivHouseholderQR<_MatrixType>&);
template<typename _MatrixType> SEXP wrap(const Eigen::CompleteOrthogonalDecomposition<_MatrixType>&);
template<typename _MatrixType, int _UpLo> SEXP wrap(const Eigen::LDLT<_MatrixType,_UpLo>&);
template<typename _MatrixType> SEXP wrap(const Eigen::BDCSVD<_MatrixType>&);
template<typename _MatrixType> SEXP wrap(const Eigen::FullPivLU<_MatrixType>&);
template<typename _MatrixType> SEXP wrap(const Eigen::SelfAdjointEigenSolver<_MatrixType >&);
template<typename _MatrixType> SEXP wrap(const Eigen::EigenSolver<_MatrixType>&);
template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename _StorageIndex> SEXP wrap(const Eigen::Transpositions<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex>&);
template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename _StorageIndex> SEXP wrap(const Eigen::PermutationMatrix<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex>&);

template<typename _Scalar, int _AmbientDim> SEXP wrap(const Eigen::AlignedBox< _Scalar, _AmbientDim >&);
template<typename _Scalar, int _Dim, int _Mode, int _Options> SEXP wrap(const Eigen::Transform< _Scalar, _Dim, _Mode, _Options >&);
template<typename _Scalar, int _Dim> SEXP wrap(const Eigen::Translation< _Scalar, _Dim>&);
template<typename _Scalar> SEXP wrap(const Eigen::AngleAxis< _Scalar >&);
template<typename _Scalar, int _Options> SEXP wrap(const Eigen::Quaternion< _Scalar, _Options >&);
template<typename _Scalar> SEXP wrap(const Eigen::Rotation2D< _Scalar >&);
template<typename _Scalar> SEXP wrap(const Eigen::JacobiRotation< _Scalar >&);
// template<typename _MatrixType1> SEXP wrap(Eigen::internal::traits<Eigen::HouseholderQR<_MatrixType1>>&);
// template<typename ExpressionType, int Direction> SEXP wrap(const Eigen::VectorwiseOp<ExpressionType,Direction>&);
// template<typename MatrixType, typename MemberOp, int Direction> SEXP wrap(const Eigen::PartialReduxExpr< MatrixType, MemberOp, Direction >&);

    namespace traits {

	/* support for as */
	template<typename T> class Exporter< Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> > >;
	template<typename T> class Exporter< Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > >;
	template<typename T> class Exporter< Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1> > >;
	template<typename T> class Exporter< Eigen::Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> > >;
	template<typename T> class Exporter< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >;
	template<typename T> class Exporter< Eigen::Matrix<T, Eigen::Dynamic, 1> >;
	template<typename T> class Exporter< Eigen::Matrix<T, 1, Eigen::Dynamic> >;
	template<typename T> class Exporter< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> >;
	template<typename T> class Exporter< Eigen::Array<T, Eigen::Dynamic, 1> >;
	template<typename T> class Exporter< Eigen::Array<T, 1, Eigen::Dynamic> >;
    template<typename T> class Exporter< Eigen::Map<Eigen::SparseMatrix<T> > >;
	template<typename T> class Exporter< Eigen::MappedSparseMatrix<T> >;  // Deprecated
	template<typename T> class Exporter< Eigen::SparseMatrix<T> >;
    template<typename T> class Exporter< Eigen::Map<Eigen::SparseMatrix<T, Eigen::RowMajor> > >;
	template<typename T> class Exporter< Eigen::MappedSparseMatrix<T, Eigen::RowMajor> >;  // Deprecated
	template<typename T> class Exporter< Eigen::SparseMatrix<T, Eigen::RowMajor> >;

    } // namespace traits

}

#endif
