// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-
//
// RcppEigenWrap.h: Rcpp wrap methods for Eigen matrices, vectors and arrays
//
// Copyright (C) 2011 - 2022   Douglas Bates, Dirk Eddelbuettel and Romain Francois
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

#ifndef ROpenEigen__ROpenEigenWrap__h
#define ROpenEigen__ROpenEigenWrap__h


namespace Rcpp{

//     namespace ROpenEigen{
//
//         // helper trait to identify if T is a plain object type
//         // TODO: perhaps move this to its own file
//         template <typename T> struct is_plain : Rcpp::traits::same_type<T,typename T::PlainObject>{} ;
//
//         // helper trait to identify if the object has dense storage
//         template <typename T> struct is_dense : Rcpp::traits::same_type<typename T::StorageKind,Eigen::Dense>{} ;
//
//         // for plain dense objects
//         template <typename T>
//         SEXP eigen_wrap_plain_dense( const T& obj, Rcpp::traits::true_type ) {
//             bool needs_dim = T::ColsAtCompileTime != 1;
//             R_xlen_t m = obj.rows(), n = obj.cols();
//             if (needs_dim && (m > INT_MAX || n > INT_MAX)) {
//                 Rcpp::stop("array dimensions cannot exceed INT_MAX");
//             }
//             R_xlen_t size = m * n;
//             // typename Eigen::internal::conditional<
//             //     T::IsRowMajor,
//             //     Eigen::Matrix<typename T::Scalar,
//             //                   T::RowsAtCompileTime,
//             //                   T::ColsAtCompileTime>,
//             //     const T&>::type objCopy(obj);
//             // SEXP ans = PROTECT(::Rcpp::wrap(objCopy.data(),
//             //                                 objCopy.data() + size));
//             // if (needs_dim) {
//             //     SEXP dd = PROTECT(::Rf_allocVector(INTSXP, 2));
//             //     int *d = INTEGER(dd);
//             //     d[0] = m;
//             //     d[1] = n;
//             //     ::Rf_setAttrib(ans, R_DimSymbol, dd);
//             //     UNPROTECT(1);
//             // }
//             // UNPROTECT(1);
//             // return ans;
//
//             typename Eigen::internal::conditional<
//               T::IsRowMajor,
//               Eigen::Matrix<typename T::Scalar,
//                             T::RowsAtCompileTime,
//                             T::ColsAtCompileTime>,
//                             const T&>::type objCopy(obj);
//
//             return Rcpp::internal::make_new_object(objCopy.ptr);
//         }
//
//         // for plain sparse objects
//         template <typename T>
//         SEXP eigen_wrap_plain_dense( const T& object, Rcpp::traits::false_type ){
//             typedef typename T::Scalar     Scalar;
//             const int  RTYPE = Rcpp::traits::r_sexptype_traits<Scalar>::rtype;
//             std::string klass;
//             switch(RTYPE) {
//             case REALSXP: klass = T::IsRowMajor ? "dgRMatrix" : "dgCMatrix";
//                 break;
// //          case INTSXP: klass = T::IsRowMajor ? "igRMatrix" : "igCMatrix";  // classes not exported
// //              break;
//             default:
//                 throw std::invalid_argument("RTYPE not matched in conversion to sparse matrix");
//             }
//             S4           ans(klass);
//             const int    nnz = object.nonZeros();
//             ans.slot("Dim")  = Dimension(object.rows(), object.cols());
//             ans.slot(T::IsRowMajor ? "j" : "i") =
//                 IntegerVector(object.innerIndexPtr(), object.innerIndexPtr() + nnz);
//             ans.slot("p")    = IntegerVector(object.outerIndexPtr(),
//                                              object.outerIndexPtr() + object.outerSize() + 1);
//             ans.slot("x")    = Vector<RTYPE>(object.valuePtr(), object.valuePtr() + nnz);
//             return  ans;
//         }
//
//         // plain object, so we can assume data() and size()
//         template <typename T>
//         inline SEXP eigen_wrap_is_plain( const T& obj, ::Rcpp::traits::true_type ){
//             return eigen_wrap_plain_dense( obj, typename is_dense<T>::type() ) ;
//         }
//
//         // when the object is not plain, we need to eval()uate it
//         template <typename T>
//         inline SEXP eigen_wrap_is_plain( const T& obj, ::Rcpp::traits::false_type ){
//             return eigen_wrap_is_plain( obj.eval(), Rcpp::traits::true_type() ) ;
//         }
//
//
//         // at that point we know that T derives from EigenBase
//         // so it is either a plain object (Matrix, etc ...) or an expression
//         // that eval()uates into a plain object
//         //
//         // so the first thing we need to do is to find out so that we don't evaluate if we don't need to
//         template <typename T>
//         inline SEXP eigen_wrap( const T& obj ){
//             return eigen_wrap_is_plain( obj,
//                 typename is_plain<T>::type()
//                 ) ;
//         }
//
//     } /* namespace RcppEigen */

    // template<typename _MatrixType> SEXP wrap(const Eigen::HouseholderQR<_MatrixType>& data){
    //   const Eigen::HouseholderQR<_MatrixType>* xp = new Eigen::HouseholderQR<_MatrixType>();
    //   xp = &data;
    //   return Rcpp::internal::make_new_object(xp);
    // }

    // template<typename _MatrixType> SEXP wrap(const Eigen::EigenHouseholderQR<_MatrixType>& data){
    //   const Eigen::EigenHouseholderQR<_MatrixType>* xp = new Eigen::EigenHouseholderQR<_MatrixType>();
    //   xp = &data;
    //   return Rcpp::internal::make_new_object(xp);
    // }

    // template<typename ExpressionType, int Direction> SEXP wrap(const Eigen::VectorwiseOp<ExpressionType,Direction>& data){
    //   // const Eigen::HouseholderQR<_MatrixType>* xp = new Eigen::HouseholderQR<_MatrixType>();
    //   // xp = &data;
    //   const Eigen::VectorwiseOp<ExpressionType,Direction>* xp = new Eigen::VectorwiseOp<ExpressionType,Direction>(data);
    //   return Rcpp::internal::make_new_object(xp);
    // }

    // template<typename MatrixType, typename MemberOp, int Direction> SEXP wrap(const Eigen::PartialReduxExpr< MatrixType, MemberOp, Direction >& data){
    //   const Eigen::PartialReduxExpr< MatrixType, MemberOp, Direction >* xp = new Eigen::PartialReduxExpr< MatrixType, MemberOp, Direction>(data);
    //   return Rcpp::internal::make_new_object(xp);
    // }

    template<typename _MatrixType> SEXP wrap(const Eigen::HouseholderQR<_MatrixType>& data){
      // const Eigen::HouseholderQR<_MatrixType>* xp = new Eigen::HouseholderQR<_MatrixType>();
      // xp = &data;
      const Eigen::HouseholderQR<_MatrixType>* xp = new Eigen::HouseholderQR<_MatrixType>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _MatrixType> SEXP wrap(const Eigen::ColPivHouseholderQR<_MatrixType>& data){
      const Eigen::ColPivHouseholderQR<_MatrixType>* xp = new Eigen::ColPivHouseholderQR<_MatrixType>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _MatrixType> SEXP wrap(const Eigen::CompleteOrthogonalDecomposition<_MatrixType>&data){
      const Eigen::CompleteOrthogonalDecomposition<_MatrixType>* xp = new Eigen::CompleteOrthogonalDecomposition<_MatrixType>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _MatrixType, int _UpLo> SEXP wrap(const Eigen::LDLT< _MatrixType,_UpLo>& data){
      const Eigen::LDLT< _MatrixType,_UpLo>* xp = new Eigen::LDLT< _MatrixType,_UpLo>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _MatrixType> SEXP wrap(const Eigen::BDCSVD<_MatrixType>& data){
      const Eigen::BDCSVD<_MatrixType>* xp = new Eigen::BDCSVD<_MatrixType>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _MatrixType> SEXP wrap(const Eigen::FullPivLU<_MatrixType>& data){
      const Eigen::FullPivLU<_MatrixType>* xp = new Eigen::FullPivLU<_MatrixType>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _MatrixType> SEXP wrap(const Eigen::SelfAdjointEigenSolver<_MatrixType>& data){
      const Eigen::SelfAdjointEigenSolver<_MatrixType>* xp = new Eigen::SelfAdjointEigenSolver<_MatrixType>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _MatrixType> SEXP wrap(const Eigen::EigenSolver<_MatrixType >& data){
      const Eigen::EigenSolver<_MatrixType>* xp = new Eigen::EigenSolver<_MatrixType>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename _StorageIndex> SEXP wrap(const Eigen::Transpositions<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex>& data){
      const Eigen::Transpositions<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex>* xp = new Eigen::Transpositions<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename _StorageIndex> SEXP wrap(const Eigen::PermutationMatrix<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex>& data){
      const Eigen::PermutationMatrix<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex>* xp = new Eigen::PermutationMatrix<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    // template<typename _MatrixType1> SEXP wrap(Eigen::internal::traits<Eigen::HouseholderQR<_MatrixType1>& data){
    //   Eigen::internal::traits<Eigen::HouseholderQR<_MatrixType1>>* xp = new Eigen::internal::traits<Eigen::HouseholderQR<_MatrixType1>>();
    //   xp = &data;
    //   return Rcpp::internal::make_new_object(xp);
    // }

    // template<typename _Scalar> SEXP wrap(const Eigen::HouseholderQR<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>& data){
    //   const Eigen::HouseholderQR<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>* xp = new Eigen::HouseholderQR<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>();
    //   xp = &data;
    //   return Rcpp::internal::make_new_object(xp);
    // }

    // template<typename _Scalar, int _Rows, int _Cols> SEXP wrap(const Eigen::HouseholderQR<Eigen::Matrix<_Scalar, _Rows, _Cols>>& data){
    //   const Eigen::HouseholderQR<Eigen::Matrix<_Scalar, _Rows, _Cols>>* xp = new Eigen::HouseholderQR<Eigen::Matrix<_Scalar, _Rows, _Cols>>(data);
    //   return Rcpp::internal::make_new_object(xp);
    // }

    // SEXP wrap(const Eigen::HouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>& data){
    //   const Eigen::HouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>* xp = new Eigen::HouseholderQR<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>();
    //   xp = &data;
    //   return Rcpp::internal::make_new_object(xp);
    // }



    template<typename _Scalar, int _Rows, int _Cols> SEXP wrap(const Eigen::Array<_Scalar, _Rows, _Cols>& data){
      const Eigen::Array<_Scalar, _Rows, _Cols>* xp = new Eigen::Array<_Scalar, _Rows, _Cols>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename Derived> SEXP wrap(const Eigen::EigenBase<Derived>& data){
      const Eigen::EigenBase<Derived>* xp = new Eigen::EigenBase<Derived>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _Scalar, int _Rows, int _Cols> SEXP wrap(const Eigen::Matrix<_Scalar, _Rows, _Cols>& data){
      const Eigen::Matrix<_Scalar, _Rows, _Cols>* xp = new Eigen::Matrix<_Scalar, _Rows, _Cols>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    // template<typename _Scalar, int _Rows> SEXP wrap(const Eigen::Matrix<_Scalar, _Rows, 1>& data){
    //   const Eigen::Matrix<_Scalar, _Rows, 1>* xp = new Eigen::Matrix<_Scalar, _Rows, 1>(data);
    //   return Rcpp::internal::make_new_object(xp);
    // }
    //
    // template<typename _Scalar, int _Cols> SEXP wrap(const Eigen::Matrix<_Scalar, 1, _Cols>& data){
    //   const Eigen::Matrix<_Scalar, 1, _Cols>* xp = new Eigen::Matrix<_Scalar, 1, _Cols>(data);
    //   return Rcpp::internal::make_new_object(xp);
    // }



    // template<typename _MatrixType> SEXP wrap(const Eigen::FullPivHouseholderQR<_MatrixType>& data){
    //   const Eigen::Matrix<_Scalar, _Rows, _Cols>* xp = new Eigen::Matrix<_Scalar, _Rows, _Cols>(data);
    //   return Rcpp::internal::make_new_object(xp);
    // }

    // template<typename ExpressionType, int Direction> SEXP wrap(const Eigen::VectorwiseOp< ExpressionType, Direction >& data){
    //   const Eigen::VectorwiseOp< ExpressionType, Direction >* xp = new Eigen::VectorwiseOp< ExpressionType, Direction >(data._expression());
    //   return Rcpp::internal::make_new_object(xp);
    // }

    // template<typename MatrixType, typename MemberOp, int Direction> SEXP wrap(const Eigen::PartialReduxExpr< MatrixType, MemberOp, Direction >& data){
    //   const Eigen::PartialReduxExpr< MatrixType, MemberOp, Direction >* xp = new Eigen::PartialReduxExpr< MatrixType, MemberOp, Direction >(data.nestedExpression(),data.functor());
    //   return Rcpp::internal::make_new_object(xp);
    // }

    // template<typename XprType, int BlockRows, int BlockCols, bool InnerPanel> SEXP wrap(Eigen::Block< XprType, BlockRows, BlockCols, InnerPanel >& data){
    //   const Eigen::Block< XprType, BlockRows, BlockCols, InnerPanel >* xp=new Eigen::Block< XprType, BlockRows, BlockCols, InnerPanel >(data.nestedExpression(),data.startRow(),data.startCol());
    //   return Rcpp::internal::make_new_object(xp);
    // }

    template<int BlockRows, int BlockCols, bool InnerPanel> SEXP wrap(Eigen::Block< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, BlockRows, BlockCols, InnerPanel >& data){
      const Eigen::Block< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, BlockRows, BlockCols, InnerPanel >* xp=new Eigen::Block< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, BlockRows, BlockCols, InnerPanel >(data.nestedExpression(),data.startRow(),data.startCol());
      return Rcpp::internal::make_new_object(xp);
    }

    // template<int BlockRows, int BlockCols, bool InnerPanel> SEXP wrap(Eigen::Block< Eigen::Matrix<double, -1, -1, 0, -1, -1>, BlockRows, BlockCols, InnerPanel >& data){
    //   const Eigen::Block< Eigen::Matrix<double, -1, -1, 0, -1, -1>, BlockRows, BlockCols, InnerPanel >* xp=new Eigen::Block< Eigen::Matrix<double, -1, -1, 0, -1, -1>, BlockRows, BlockCols, InnerPanel >(data.nestedExpression(),data.startRow(),data.startCol());
    //   return Rcpp::internal::make_new_object(xp);
    // }

    template<typename _Scalar, int _AmbientDim> SEXP wrap(const Eigen::AlignedBox< _Scalar, _AmbientDim >& data){
      const Eigen::AlignedBox< _Scalar, _AmbientDim >* xp=new Eigen::AlignedBox< _Scalar, _AmbientDim >(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _Scalar, int _Dim, int _Mode, int _Options> SEXP wrap(const Eigen::Transform< _Scalar, _Dim, _Mode, _Options >& data){
      const Eigen::Transform< _Scalar, _Dim, _Mode, _Options >* xp=new Eigen::Transform< _Scalar, _Dim, _Mode, _Options >(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _Scalar, int _Dim> SEXP wrap(const Eigen::Translation< _Scalar, _Dim>& data){
      const Eigen::Translation< _Scalar, _Dim>* xp=new Eigen::Translation< _Scalar, _Dim>(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _Scalar> SEXP wrap(const Eigen::AngleAxis< _Scalar >&data){
      const  Eigen::AngleAxis< _Scalar >* xp=new Eigen::AngleAxis< _Scalar >(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _Scalar, int _Options> SEXP wrap(const Eigen::Quaternion< _Scalar, _Options >& data){
      const Eigen::Quaternion< _Scalar, _Options >* xp=new Eigen::Quaternion< _Scalar, _Options >(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _Scalar> SEXP wrap(const Eigen::Rotation2D< _Scalar >&data){
      const  Eigen::Rotation2D< _Scalar >* xp=new Eigen::Rotation2D< _Scalar >(data);
      return Rcpp::internal::make_new_object(xp);
    }

    template<typename _Scalar> SEXP wrap(const Eigen::JacobiRotation< _Scalar >&data){
      const  Eigen::JacobiRotation< _Scalar >* xp=new Eigen::JacobiRotation< _Scalar >(data);
      return Rcpp::internal::make_new_object(xp);
    }

    namespace traits {
        /* support for Rcpp::as */

        template<int SizeAtCompileTime, int MaxSizeAtCompileTime, typename _StorageIndex> class Exporter< Eigen::Transpositions<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex> >{
        public:
          typedef Eigen::Transpositions<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex> trans_t;
          // enum { RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype };
          // typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
          Exporter( SEXP x ) : sp(x){}
          trans_t get(){
            Eigen::Transpositions<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex> result;
            bool isTrans = false;
            // Rcout << sp << endl;
            if(TYPEOF(sp) == S4SXP){
              // Rcout << "ha, environment" << endl;
              Rcpp::Environment ev(sp);
              SEXP objptr = ev[".pointer"];
              SEXP clsxptr = ev.get(".cppclass");
              if (TYPEOF(clsxptr) != EXTPTRSXP){
                std::string msg =
                  "Passing parameter should be a pointer to rcpp clas\n";
                Rcpp::stop(msg);
              }else{
                XPtr<class_Base> xp(clsxptr);
                // Rcout << "clasname: " << xp->name << endl;
                std::string clsname = xp->name;
                if(clsname.find("EigenTranspositions")==string::npos){
                  std::string msg =
                    "Error converting object to eigen::Transpositions\n";
                  Rcpp::stop(msg);
                }else{
                  Eigen::Transpositions<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex>* ptr = (Eigen::Transpositions<SizeAtCompileTime, MaxSizeAtCompileTime, _StorageIndex>*) R_ExternalPtrAddr(objptr);
                  result = *ptr;
                  isTrans = true;
                }
              }

            }
            if(!isTrans){
              throw std::runtime_error("Parameter maybe wrong, error in get a Transpositions.");
            }
            return result;
          }
        private:
          SEXP sp;
        };

        template<typename _Scalar, int _AmbientDim> class Exporter< Eigen::AlignedBox< _Scalar, _AmbientDim > >{
        public:
          typedef Eigen::AlignedBox< _Scalar, _AmbientDim > box_t;
          // enum { RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype };
          // typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
          Exporter( SEXP x ) : sp(x){}
          box_t get(){
            Eigen::AlignedBox< _Scalar, _AmbientDim > result;
            bool isAlignedBox = false;
            if(TYPEOF(sp) == S4SXP){
              Rcpp::Environment ev(sp);
              SEXP objptr = ev[".pointer"];
              SEXP clsxptr = ev.get(".cppclass");
              if (TYPEOF(clsxptr) != EXTPTRSXP){
                std::string msg =
                  "Passing parameter should be a pointer to rcpp clas\n";
                Rcpp::stop(msg);
              }else{
                XPtr<class_Base> xp(clsxptr);
                std::string clsname = xp->name;
                if(clsname.find("EigenAlignedBox")==string::npos){
                  std::string msg =
                    "Error converting object to eigen::AlignedBox\n";
                  Rcpp::stop(msg);
                }else{
                  Eigen::AlignedBox< _Scalar, _AmbientDim >* ptr = (Eigen::AlignedBox< _Scalar, _AmbientDim >*) R_ExternalPtrAddr(objptr);
                  result = *ptr;
                  isAlignedBox = true;
                }
              }

            }
            if(!isAlignedBox){
              throw std::runtime_error("Parameter maybe wrong, error in get a Transpositions.");
            }
            return result;
          }
        private:
          SEXP sp;
        };

        template<typename _Scalar, int _Dim, int _Mode, int _Options> class Exporter< Eigen::Transform< _Scalar, _Dim, _Mode, _Options > >{
        public:
          typedef Eigen::Transform< _Scalar, _Dim, _Mode, _Options > trans_t;
          // enum { RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype };
          // typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
          Exporter( SEXP x ) : sp(x){}
          trans_t get(){
            Eigen::Transform< _Scalar, _Dim, _Mode, _Options > result;
            bool isTransform = false;
            if(TYPEOF(sp) == S4SXP){
              Rcpp::Environment ev(sp);
              SEXP objptr = ev[".pointer"];
              SEXP clsxptr = ev.get(".cppclass");
              if (TYPEOF(clsxptr) != EXTPTRSXP){
                std::string msg =
                  "Passing parameter should be a pointer to rcpp clas\n";
                Rcpp::stop(msg);
              }else{
                XPtr<class_Base> xp(clsxptr);
                std::string clsname = xp->name;
                if(clsname.find("EigenTransform")==string::npos){
                  std::string msg =
                    "Error converting object to eigen::Transform\n";
                  Rcpp::stop(msg);
                }else{
                  Eigen::Transform< _Scalar, _Dim, _Mode, _Options >* ptr = (Eigen::Transform< _Scalar, _Dim, _Mode, _Options >*) R_ExternalPtrAddr(objptr);
                  result = *ptr;
                  isTransform = true;
                }
              }

            }
            if(!isTransform){
              throw std::runtime_error("Parameter maybe wrong, error in get a Transform.");
            }
            return result;
          }
        private:
          SEXP sp;
        };

        template<typename _Scalar, int _Dim> class Exporter< Eigen::Translation< _Scalar, _Dim> >{
        public:
          typedef Eigen::Translation< _Scalar, _Dim > trans_t;
          // enum { RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype };
          // typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
          Exporter( SEXP x ) : sp(x){}
          trans_t get(){
            Eigen::Translation< _Scalar, _Dim > result;
            bool isTranslation = false;
            if(TYPEOF(sp) == S4SXP){
              Rcpp::Environment ev(sp);
              SEXP objptr = ev[".pointer"];
              SEXP clsxptr = ev.get(".cppclass");
              if (TYPEOF(clsxptr) != EXTPTRSXP){
                std::string msg =
                  "Passing parameter should be a pointer to rcpp clas\n";
                Rcpp::stop(msg);
              }else{
                XPtr<class_Base> xp(clsxptr);
                std::string clsname = xp->name;
                if(clsname.find("EigenTranslation")==string::npos){
                  std::string msg =
                    "Error converting object to eigen::Translation\n";
                  Rcpp::stop(msg);
                }else{
                  Eigen::Translation< _Scalar, _Dim >* ptr = (Eigen::Translation< _Scalar, _Dim >*) R_ExternalPtrAddr(objptr);
                  result = *ptr;
                  isTranslation = true;
                }
              }

            }
            if(!isTranslation){
              throw std::runtime_error("Parameter maybe wrong, error in get a Translation.");
            }
            return result;
          }
        private:
          SEXP sp;
        };

        template<typename _Scalar> class Exporter< Eigen::AngleAxis< _Scalar> >{
        public:
          typedef Eigen::AngleAxis< _Scalar > trans_t;
          // enum { RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype };
          // typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
          Exporter( SEXP x ) : sp(x){}
          trans_t get(){
            Eigen::AngleAxis< _Scalar > result;
            bool isAngleAxis = false;
            if(TYPEOF(sp) == S4SXP){
              Rcpp::Environment ev(sp);
              SEXP objptr = ev[".pointer"];
              SEXP clsxptr = ev.get(".cppclass");
              if (TYPEOF(clsxptr) != EXTPTRSXP){
                std::string msg =
                  "Passing parameter should be a pointer to rcpp clas\n";
                Rcpp::stop(msg);
              }else{
                XPtr<class_Base> xp(clsxptr);
                std::string clsname = xp->name;
                if(clsname.find("EigenAngleAxis")==string::npos){
                  std::string msg =
                    "Error converting object to eigen::AngleAxis\n";
                  Rcpp::stop(msg);
                }else{
                  Eigen::AngleAxis< _Scalar >* ptr = (Eigen::AngleAxis< _Scalar >*) R_ExternalPtrAddr(objptr);
                  result = *ptr;
                  isAngleAxis = true;
                }
              }

            }
            if(!isAngleAxis){
              throw std::runtime_error("Parameter maybe wrong, error in get an AngleAxis.");
            }
            return result;
          }
        private:
          SEXP sp;
        };

        template<typename _Scalar, int _Options> class Exporter< Eigen::Quaternion< _Scalar, _Options> >{
        public:
          typedef Eigen::Quaternion< _Scalar, _Options > trans_t;
          // enum { RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype };
          // typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
          Exporter( SEXP x ) : sp(x){}
          trans_t get(){
            Eigen::Quaternion< _Scalar, _Options > result;
            bool isQuaternion = false;
            if(TYPEOF(sp) == S4SXP){
              Rcpp::Environment ev(sp);
              SEXP objptr = ev[".pointer"];
              SEXP clsxptr = ev.get(".cppclass");
              if (TYPEOF(clsxptr) != EXTPTRSXP){
                std::string msg =
                  "Passing parameter should be a pointer to rcpp clas\n";
                Rcpp::stop(msg);
              }else{
                XPtr<class_Base> xp(clsxptr);
                std::string clsname = xp->name;
                if(clsname.find("EigenQuaternion")==string::npos){
                  std::string msg =
                    "Error converting object to eigen::Quaternion\n";
                  Rcpp::stop(msg);
                }else{
                  Eigen::Quaternion< _Scalar, _Options >* ptr = (Eigen::Quaternion< _Scalar, _Options >*) R_ExternalPtrAddr(objptr);
                  result = *ptr;
                  isQuaternion = true;
                }
              }

            }
            if(!isQuaternion){
              throw std::runtime_error("Parameter maybe wrong, error in get a Quaternion.");
            }
            return result;
          }
        private:
          SEXP sp;
        };

        template<typename _Scalar> class Exporter< Eigen::Rotation2D< _Scalar> >{
        public:
          typedef Eigen::Rotation2D< _Scalar > trans_t;
          // enum { RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype };
          // typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
          Exporter( SEXP x ) : sp(x){}
          trans_t get(){
            Eigen::Rotation2D< _Scalar > result;
            bool isRotation2D = false;
            if(TYPEOF(sp) == S4SXP){
              Rcpp::Environment ev(sp);
              SEXP objptr = ev[".pointer"];
              SEXP clsxptr = ev.get(".cppclass");
              if (TYPEOF(clsxptr) != EXTPTRSXP){
                std::string msg =
                  "Passing parameter should be a pointer to rcpp clas\n";
                Rcpp::stop(msg);
              }else{
                XPtr<class_Base> xp(clsxptr);
                std::string clsname = xp->name;
                if(clsname.find("EigenRotation2D")==string::npos){
                  std::string msg =
                    "Error converting object to eigen::Rotation2D\n";
                  Rcpp::stop(msg);
                }else{
                  Eigen::Rotation2D< _Scalar >* ptr = (Eigen::Rotation2D< _Scalar >*) R_ExternalPtrAddr(objptr);
                  result = *ptr;
                  isRotation2D = true;
                }
              }

            }
            if(!isRotation2D){
              throw std::runtime_error("Parameter maybe wrong, error in get an Rotation2D.");
            }
            return result;
          }
        private:
          SEXP sp;
        };

        template<typename _Scalar> class Exporter< Eigen::JacobiRotation< _Scalar> >{
        public:
          typedef Eigen::JacobiRotation< _Scalar > trans_t;
          // enum { RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype };
          // typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
          Exporter( SEXP x ) : sp(x){}
          trans_t get(){
            Eigen::JacobiRotation< _Scalar > result;
            bool isJacobiRotation = false;
            if(TYPEOF(sp) == S4SXP){
              Rcpp::Environment ev(sp);
              SEXP objptr = ev[".pointer"];
              SEXP clsxptr = ev.get(".cppclass");
              if (TYPEOF(clsxptr) != EXTPTRSXP){
                std::string msg =
                  "Passing parameter should be a pointer to rcpp clas\n";
                Rcpp::stop(msg);
              }else{
                XPtr<class_Base> xp(clsxptr);
                std::string clsname = xp->name;
                if(clsname.find("JacobiRotation")==string::npos){
                  std::string msg =
                    "Error converting object to eigen::JacobiRotation\n";
                  Rcpp::stop(msg);
                }else{
                  Eigen::JacobiRotation< _Scalar >* ptr = (Eigen::JacobiRotation< _Scalar >*) R_ExternalPtrAddr(objptr);
                  result = *ptr;
                  isJacobiRotation = true;
                }
              }

            }
            if(!isJacobiRotation){
              throw std::runtime_error("Parameter maybe wrong, error in get an JacobiRotation.");
            }
            return result;
          }
        private:
          SEXP sp;
        };

        template<typename Derived> class Exporter<Eigen::EigenBase<Derived>>{
        public:
          Exporter( SEXP x ) : sp(x){}
          inline Eigen::EigenBase<Derived> get(){
            Eigen::EigenBase<Derived> result;
            bool isEigenBase = false;
            if(TYPEOF(sp) == S4SXP){
              Rcpp::Environment ev(sp);
              SEXP objptr = ev[".pointer"];
              SEXP clsxptr = ev.get(".cppclass");
              if (TYPEOF(clsxptr) != EXTPTRSXP){
                std::string msg =
                  "Passing parameter should be a pointer to rcpp class\n";
                Rcpp::stop(msg);
              }else{
                XPtr<class_Base> xp(clsxptr);
                std::string clsname = xp->name;
                if(clsname.find("EigenMatrixX")!=string::npos){
                  Eigen::EigenBase<Derived>* ptr = (Eigen::EigenBase<Derived>*) R_ExternalPtrAddr(objptr);
                  result = *ptr;
                  isEigenBase = true;
                }else{
                  std::string msg =
                    "Error converting object to EigenBase.\n";
                  Rcpp::stop(msg);
                }
              }
            }
            if(!isEigenBase){
              throw std::runtime_error("Parameter maybe wrong, error in get a EigenBase.");
            }
            return result;
          }
        private:
          SEXP sp;
        };

        // template <typename _Scalar, int _Rows, int _Cols>
        // class Exporter< Eigen::Matrix<_Scalar, _Rows, _Cols> > {
        // public:
        //   Exporter(SEXP x) : sp(x) {}
        //
        //   Eigen::Matrix<_Scalar, _Rows, _Cols> get() {
        //     Eigen::Matrix<_Scalar, _Rows, _Cols> result;
        //     bool isMat=false;
        //     if(TYPEOF(sp) == S4SXP){
        //       // Rcout << sp->get_typeinfo_name() << endl;
        //       // Rcout << "ha, environment" << endl;
        //       Rcpp::Environment ev(sp);
        //       // if(Rcpp::as<std::string>(ev.attr("class")).find("cvSize")==string::npos){
        //       //   std::string msg =
        //       //     "Error converting object to cv::Size:\n"
        //       //     "Input parameter should be cvSize.\n";
        //       //   Rcpp::stop(msg);
        //       // }
        //       SEXP objptr = ev[".pointer"];
        //       SEXP clsxptr = ev.get(".cppclass");
        //       if (TYPEOF(clsxptr) != EXTPTRSXP){
        //         std::string msg =
        //           "Passing parameter should be a pointer to rcpp class\n";
        //         Rcpp::stop(msg);
        //       }else{
        //         XPtr<class_Base> xp(clsxptr);
        //         // Rcout << "Here, clasname: " << xp->name << endl;
        //         std::string clsname = xp->name;
        //         if(clsname.find("EigenMatrixXd")!=string::npos){
        //           Eigen::Matrix<_Scalar, _Rows, _Cols>* ptr = (Eigen::Matrix<_Scalar, _Rows, _Cols>*) R_ExternalPtrAddr(objptr);
        //           result = *ptr;
        //           isMat = true;
        //         }else if(xp.inherits("eigen-mat")){
        //           Rcpp::XPtr< Eigen::Matrix<_Scalar, _Rows, _Cols> > p1(sp) ;
        //           if(p1.get() == NULL)
        //             throw std::runtime_error("data has been destroyed");
        //           result = *p1;
        //           isMat = true;
        //         }else{
        //           std::string msg =
        //             "Error converting object to Eigen::Matrix\n";
        //           Rcpp::stop(msg);
        //         }
        //       }
        //     }
        //     if(!isMat){
        //       throw std::runtime_error("Parameter maybe wrong, error in get a mat.");
        //     }
        //     return result;
        //   }
        // private:
        //   // Rcpp::NumericVector vec;
        //   SEXP sp;
        // };

        template <typename T, int RTYPE>
        class Eigen_Matrix_Exporter {
            public:
            Eigen_Matrix_Exporter(SEXP x) : vec(x), d_ncol(1), d_nrow(Rf_xlength(x)),sp(x) {
                // if (TYPEOF(x) != RTYPE)
                //     throw std::invalid_argument("Wrong R type for mapped vector");
                if (::Rf_isMatrix(x)) {
                    int *dims = vec.dims() ;
                    d_nrow = dims[0];
                    d_ncol = dims[1];
                }
            }
            T get() {return T(vec.begin(), d_nrow, d_ncol );}
            // T get(){
            //   T result;
            //   bool isMat=false;
            //   if(TYPEOF(sp) == S4SXP){
            //     // Rcout << sp->get_typeinfo_name() << endl;
            //     // Rcout << "ha, environment" << endl;
            //     Rcpp::Environment ev(sp);
            //     // if(Rcpp::as<std::string>(ev.attr("class")).find("cvSize")==string::npos){
            //     //   std::string msg =
            //     //     "Error converting object to cv::Size:\n"
            //     //     "Input parameter should be cvSize.\n";
            //     //   Rcpp::stop(msg);
            //     // }
            //     SEXP objptr = ev[".pointer"];
            //     SEXP clsxptr = ev.get(".cppclass");
            //     if (TYPEOF(clsxptr) != EXTPTRSXP){
            //       std::string msg =
            //         "Passing parameter should be a pointer to rcpp class\n";
            //       Rcpp::stop(msg);
            //     }else{
            //       XPtr<class_Base> xp(clsxptr);
            //       // Rcout << "Here, clasname: " << xp->name << endl;
            //       std::string clsname = xp->name;
            //       if(clsname.find("EigenMatrixXd")!=string::npos){
            //         T* ptr = (T*) R_ExternalPtrAddr(objptr);
            //         result = *ptr;
            //         isMat = true;
            //       }else if(xp.inherits("eigen-mat")){
            //         Rcpp::XPtr< T > p1(sp) ;
            //         if(p1.get() == NULL)
            //           throw std::runtime_error("data has been destroyed");
            //         result = *p1;
            //         isMat = true;
            //       }else{
            //         std::string msg =
            //           "Error converting object to Eigen::Matrix\n";
            //         Rcpp::stop(msg);
            //       }
            //     }
            //   }
            //   if(!isMat){
            //     throw std::runtime_error("Parameter maybe wrong, error in get a mat.");
            //   }
            //   return result;
            // }
            private:
                Rcpp::Vector<RTYPE> vec ;
                int d_nrow, d_ncol ;
                SEXP sp;
        } ;

        // modified from "class MatrixExporter" in <Rcpp>/include/Rcpp/internal/Exporter.h
        // MatrixExporter uses brackets [] to index the destination matrix,
        // which is not supported by MatrixXd.
        // Here we copy data to MatrixXd.data() rather than MatrixXd
        template <typename T, typename value_type>
        class MatrixExporterForEigen {
        public:
            typedef value_type r_export_type;

            MatrixExporterForEigen(SEXP x) : object(x){}
            ~MatrixExporterForEigen(){}

            T get() {
                // Shield<SEXP> dims( ::Rf_getAttrib( object, R_DimSymbol ) );
                // if( Rf_isNull(dims) || ::Rf_length(dims) != 2 ){
                //     throw ::Rcpp::not_a_matrix();
                // }
                // int* dims_ = INTEGER(dims);
                // T result( dims_[0], dims_[1] );
                // value_type *data = result.data();
                // ::Rcpp::internal::export_indexing<value_type*, value_type>( object, data );
                // return result ;
                  T result;
                  bool isMat=false;
                  if(TYPEOF(object) == S4SXP){
                    Rcpp::Environment ev(object);
                    SEXP objptr = ev[".pointer"];
                    SEXP clsxptr = ev.get(".cppclass");
                    if (TYPEOF(clsxptr) != EXTPTRSXP){
                      std::string msg =
                        "Passing parameter should be a pointer to rcpp class\n";
                      Rcpp::stop(msg);
                    }else{
                      XPtr<class_Base> xp(clsxptr);
                      // Rcout << "Here, clasname: " << xp->name << endl;
                      std::string clsname = xp->name;
                      if(clsname.find("EigenMatrix")!=string::npos){
                        T* ptr = (T*) R_ExternalPtrAddr(objptr);
                        result = *ptr;
                        // Rcout << result << endl;
                        isMat = true;
                      }else if(clsname.find("EigenVector")!=string::npos){
                        // Eigen::Matrix<value_type, Eigen::Dynamic, 1>* ptr = (Eigen::Matrix<value_type, Eigen::Dynamic, 1>*) R_ExternalPtrAddr(objptr);
                        T* ptr = (T*) R_ExternalPtrAddr(objptr);
                        // Rcout << "here\n";
                        // Rcout << *ptr << endl;
                        result = *ptr;
                        // Rcout << result << endl;
                        isMat = true;
                      }else if(clsname.find("EigenRowVector")!=string::npos){
                        // Eigen::Matrix<value_type, 1, Eigen::Dynamic>* ptr = (Eigen::Matrix<value_type, 1, Eigen::Dynamic>*) R_ExternalPtrAddr(objptr);
                        T* ptr = (T*) R_ExternalPtrAddr(objptr);
                        result = *ptr;
                        // Rcout << result << endl;
                        isMat = true;
                      }else if(clsname.find("EigenArrayXX")!=string::npos){
                        T* ptr = (T*) R_ExternalPtrAddr(objptr);
                        result = *ptr;
                        // Rcout << result << endl;
                        isMat = true;
                      }else if(clsname.find("EigenArray")!=string::npos){
                        // Eigen::Array<value_type, Eigen::Dynamic, 1>* ptr = (Eigen::Array<value_type, Eigen::Dynamic, 1>*) R_ExternalPtrAddr(objptr);
                        T* ptr = (T*) R_ExternalPtrAddr(objptr);
                        result = *ptr;
                        isMat = true;
                      }else if(clsname.find("EigenRowArray")!=string::npos){
                        // Eigen::Array<value_type, 1, Eigen::Dynamic>* ptr = (Eigen::Array<value_type, 1, Eigen::Dynamic>*) R_ExternalPtrAddr(objptr);
                        T* ptr = (T*) R_ExternalPtrAddr(objptr);
                        result = *ptr;
                        isMat = true;
                      }
                    }
                  }else if(TYPEOF(object) == EXTPTRSXP){
                    XPtr<class_Base> xp(object);
                    if(xp.inherits("eigen-mat")){
                      Rcpp::XPtr< T > p1(object) ;
                      if(p1.get() == NULL)
                        throw std::runtime_error("data has been destroyed");
                      result = *p1;
                      isMat = true;
                    }
                  }else if(TYPEOF(object) == REALSXP || TYPEOF(object) == INTSXP){
                    Shield<SEXP> dims( ::Rf_getAttrib( object, R_DimSymbol ) );
                    // if( Rf_isNull(dims) || ::Rf_length(dims) != 2 ){
                    //     throw ::Rcpp::not_a_matrix();
                    // }
                    if(::Rf_length(dims)==2){
                      int* dims_ = INTEGER(dims);
                      T result1( dims_[0], dims_[1] );
                      // Rcout << "here\n";
                      value_type *data = result1.data();
                      ::Rcpp::internal::export_indexing<value_type*, value_type>( object, data );
                      // Rcout << result1 << endl;
                      result = result1 ;
                      isMat = true;
                    }else if( Rf_isNull(dims)){//is an R vector
                      // int _rows = ::Rf_length(object);
                      //
                      // const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype;
                      // Rcpp::Vector<double> vec(object);
                      // Rcpp::NumericVector vec(object);
                      // Eigen::Matrix<value_type, Eigen::Dynamic, 1> result1(vec.size(),1);
                      // Eigen::Matrix<value_type, Eigen::Dynamic, 1> result1(Rf_xlength(object),1);

                      T result1(Rf_xlength(object),1);
                      value_type *data = result1.data();
                      ::Rcpp::internal::export_indexing<value_type*, value_type>( object, data );
                      // Rcout << result1 << endl;
                      // Eigen::Matrix<value_type, Eigen::Dynamic, 1> result1(vec.begin(), Rf_xlength(object), 1 );
                      // for(size_t i=0;i<vec.size();i++){
                      //   result1(i) = vec[i];
                      // }
                      result = result1 ;

                      isMat = true;
                    }

                  }else{
                    std::string msg =
                      "Error converting object to Eigen::Matrix\n";
                    Rcpp::stop(msg);
                  }
                  if(!isMat){
                    throw std::runtime_error("Parameter maybe wrong, error in get a mat.");
                  }
                  return result;
            }

        private:
            SEXP object;
        };

        // Provides only Map::VectorX<t> export
        template<typename T>
        class Exporter<Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> > > {
            typedef typename Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> > OUT ;
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
            Rcpp::Vector<RTYPE> vec ;

        public:
            Exporter(SEXP x) : vec(x) {
                if (TYPEOF(x) != RTYPE)
                    throw std::invalid_argument("Wrong R type for mapped vector"); // #nocov
            }
            OUT get() {return OUT(vec.begin(), vec.size());}
        } ;

        // Provides only Map::RowVectorX<t> export
        template<typename T>
        class Exporter<Eigen::Map<Eigen::Matrix<T, 1, Eigen::Dynamic> > > {
            typedef typename Eigen::Map<Eigen::Matrix<T, 1, Eigen::Dynamic> > OUT ;
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
            Rcpp::Vector<RTYPE> vec ;

        public:
            Exporter(SEXP x) : vec(x) {
                if (TYPEOF(x) != RTYPE)
                    throw std::invalid_argument("Wrong R type for mapped rowvector");
            }
            OUT get() {return OUT(vec.begin(), vec.size());}
        };

        template<typename T>
        class Exporter< Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1> > > {
            typedef typename Eigen::Map<Eigen::Array<T, Eigen::Dynamic, 1> > OUT ;
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
            Rcpp::Vector<RTYPE> vec ;

        public:
            Exporter(SEXP x) : vec(x) {
                if (TYPEOF(x) != RTYPE)
                    throw std::invalid_argument("Wrong R type for mapped 1D array");
            }
            OUT get() {return OUT(vec.begin(), vec.size());}
        } ;

        template<typename T>
        class Exporter<Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > > {
            typedef typename Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > OUT ;
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
            Rcpp::Vector<RTYPE> vec ;
            int d_ncol, d_nrow ;

            public:
            Exporter(SEXP x) : vec(x), d_ncol(1), d_nrow(Rf_xlength(x)) {
                if (TYPEOF(x) != RTYPE)
                    throw std::invalid_argument("Wrong R type for mapped matrix");	// #nocov
                if (::Rf_isMatrix(x)) {
                    int *dims = INTEGER( ::Rf_getAttrib( x, R_DimSymbol ) ) ;
                    d_nrow = dims[0];
                    d_ncol = dims[1];
                }
            }
            OUT get() {return OUT(vec.begin(), d_nrow, d_ncol );}
        } ;

        template<typename T>
        class Exporter<Eigen::Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> > > {
            typedef typename Eigen::Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> > OUT ;
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
            Rcpp::Vector<RTYPE> vec ;
            int d_ncol, d_nrow ;

            public:
            Exporter(SEXP x) : vec(x), d_ncol(1), d_nrow(Rf_xlength(x)) {
                if (TYPEOF(x) != RTYPE)
                    throw std::invalid_argument("Wrong R type for mapped 2D array");
                if (::Rf_isMatrix(x)) {
                    int *dims = INTEGER( ::Rf_getAttrib( x, R_DimSymbol ) ) ;
                    d_nrow = dims[0];
                    d_ncol = dims[1];
                }
            }
            OUT get() {return OUT(vec.begin(), d_nrow, d_ncol );}
        } ;

        // template <typename T>
        // class Exporter<Eigen::Matrix<T, Eigen::Dynamic, 1> >
        //     : public IndexingExporter<Eigen::Matrix<T, Eigen::Dynamic, 1>, T> {
        // public:
        //     Exporter(SEXP x) : IndexingExporter<Eigen::Matrix<T, Eigen::Dynamic, 1>, T >(x){}
        // };

        // template <typename T>
        // class Exporter<Eigen::Array<T, Eigen::Dynamic, 1> >
        //     : public IndexingExporter<Eigen::Array<T, Eigen::Dynamic, 1>, T> {
        // public:
        //     Exporter(SEXP x) : IndexingExporter<Eigen::Array<T, Eigen::Dynamic, 1>, T >(x){}
        // };

        // template <typename T>
        // class Exporter< Eigen::Matrix<T, 1, Eigen::Dynamic> >
        //     : public IndexingExporter< Eigen::Matrix<T, 1, Eigen::Dynamic>, T > {
        // public:
        //     Exporter(SEXP x) : IndexingExporter< Eigen::Matrix<T, 1, Eigen::Dynamic>, T >(x){}
        // };



        // template <typename T>
        // class Exporter< Eigen::Array<T, 1, Eigen::Dynamic> >
        //     : public IndexingExporter< Eigen::Array<T, 1, Eigen::Dynamic>, T > {
        // public:
        //     Exporter(SEXP x) : IndexingExporter< Eigen::Array<T, 1, Eigen::Dynamic>, T >(x){}
        // };

        template <typename T>
        class Exporter< Eigen::Matrix<T, Eigen::Dynamic, 1> >
          : public MatrixExporterForEigen< Eigen::Matrix<T, Eigen::Dynamic, 1>, T > {
        public:
          Exporter(SEXP x) :
          MatrixExporterForEigen< Eigen::Matrix<T, Eigen::Dynamic, 1>, T >(x){}
        };

        template <typename T>
        class Exporter< Eigen::Matrix<T, 1, Eigen::Dynamic> >
          : public MatrixExporterForEigen< Eigen::Matrix<T, 1, Eigen::Dynamic>, T > {
        public:
          Exporter(SEXP x) :
          MatrixExporterForEigen< Eigen::Matrix<T, 1, Eigen::Dynamic>, T >(x){}
        };

        template <typename T>
        class Exporter< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >
            : public MatrixExporterForEigen< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, T > {
        public:
            Exporter(SEXP x) :
                MatrixExporterForEigen< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>, T >(x){}
        };

        template <typename T, int _Rows, int _Cols>
        class Exporter< Eigen::Matrix<T, _Rows, _Cols> >
          : public MatrixExporterForEigen< Eigen::Matrix<T, _Rows, _Cols>, T > {
        public:
          Exporter(SEXP x) :
          MatrixExporterForEigen< Eigen::Matrix<T, _Rows, _Cols>, T >(x){}
        };

        template <typename T,int _Rows>
        class Exporter< Eigen::Matrix<T, _Rows, 1> >
          : public MatrixExporterForEigen< Eigen::Matrix<T, _Rows, 1>, T > {
        public:
          Exporter(SEXP x) :
          MatrixExporterForEigen< Eigen::Matrix<T, _Rows, 1>, T >(x){}
        };

        template <typename T,int _Cols>
        class Exporter< Eigen::Matrix<T, 1, _Cols> >
          : public MatrixExporterForEigen< Eigen::Matrix<T, 1, _Cols>, T > {
        public:
          Exporter(SEXP x) :
          MatrixExporterForEigen< Eigen::Matrix<T, 1, _Cols>, T >(x){}
        };


        template <typename T>
        class Exporter< Eigen::Array<T, Eigen::Dynamic, 1> >
          : public MatrixExporterForEigen< Eigen::Array<T, Eigen::Dynamic, 1>, T > {
        public:
          Exporter(SEXP x) :
          MatrixExporterForEigen< Eigen::Array<T, Eigen::Dynamic, 1>, T >(x){}
        };

        template <typename T,int _Rows>
        class Exporter< Eigen::Array<T, _Rows, 1> >
          : public MatrixExporterForEigen< Eigen::Array<T, _Rows, 1>, T > {
        public:
          Exporter(SEXP x) :
          MatrixExporterForEigen< Eigen::Array<T, _Rows, 1>, T >(x){}
        };

        template <typename T>
        class Exporter< Eigen::Array<T, 1, Eigen::Dynamic> >
          : public MatrixExporterForEigen< Eigen::Array<T, 1, Eigen::Dynamic>, T > {
        public:
          Exporter(SEXP x) :
          MatrixExporterForEigen< Eigen::Array<T, 1, Eigen::Dynamic>, T >(x){}
        };

        template <typename T,int _Cols>
        class Exporter< Eigen::Array<T, 1, _Cols> >
          : public MatrixExporterForEigen< Eigen::Array<T, 1, _Cols>, T > {
        public:
          Exporter(SEXP x) :
          MatrixExporterForEigen< Eigen::Array<T, 1, _Cols>, T >(x){}
        };

        template <typename T>
        class Exporter< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic> >
            : public MatrixExporterForEigen< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>, T > {
        public:
            Exporter(SEXP x) :
                MatrixExporterForEigen< Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>, T >(x){}
        };

        // Starting from Eigen 3.3 MappedSparseMatrix was deprecated.
        // The new type is Map<SparseMatrix>.
        template<typename T>
        class Exporter<Eigen::Map<Eigen::SparseMatrix<T> > > {
        public:
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
            Exporter(SEXP x) : d_x(x), d_dims(d_x.slot("Dim")), d_i(d_x.slot("i")), d_p(d_x.slot("p")), xx( d_x.slot("x") ) {
                if (!d_x.is("dgCMatrix"))
                    throw std::invalid_argument("Need S4 class dgCMatrix for a mapped sparse matrix");
            }
            Eigen::Map<Eigen::SparseMatrix<T> > get() {
                return Eigen::Map<Eigen::SparseMatrix<T> >(d_dims[0], d_dims[1], d_p[d_dims[1]],
                                                           d_p.begin(), d_i.begin(), xx.begin() );
            }
        protected:
            S4            d_x;
            IntegerVector d_dims, d_i, d_p;
            Vector<RTYPE> xx ;
        };
        // Deprecated
        template<typename T>
        class Exporter<Eigen::MappedSparseMatrix<T> > {
        public:
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
            Exporter(SEXP x) : d_x(x), d_dims(d_x.slot("Dim")), d_i(d_x.slot("i")), d_p(d_x.slot("p")), xx( d_x.slot("x") ) {
                if (!d_x.is("dgCMatrix"))
                    throw std::invalid_argument("Need S4 class dgCMatrix for a mapped sparse matrix");
            }
            Eigen::MappedSparseMatrix<T> get() {
                return Eigen::MappedSparseMatrix<T>(d_dims[0], d_dims[1], d_p[d_dims[1]],
                                                    d_p.begin(), d_i.begin(), xx.begin() );
            }
        protected:
            S4            d_x;
            IntegerVector d_dims, d_i, d_p;
            Vector<RTYPE> xx ;
        };

        // Starting from Eigen 3.3 MappedSparseMatrix was deprecated.
        // The new type is Map<SparseMatrix>.
        template<typename T>
        class Exporter<Eigen::Map<Eigen::SparseMatrix<T, Eigen::RowMajor> > > {
        public:
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
            Exporter(SEXP x) : d_x(x), d_dims(d_x.slot("Dim")), d_j(d_x.slot("j")), d_p(d_x.slot("p")), xx( d_x.slot("x") ) {
                if (!d_x.is("dgRMatrix"))
                    throw std::invalid_argument("Need S4 class dgRMatrix for a mapped sparse matrix");
            }
            Eigen::Map<Eigen::SparseMatrix<T, Eigen::RowMajor> > get() {
                return Eigen::Map<Eigen::SparseMatrix<T, Eigen::RowMajor> >(d_dims[0], d_dims[1], d_p[d_dims[1]],
                                                                            d_p.begin(), d_j.begin(), xx.begin() );
            }
        protected:
            S4            d_x;
            IntegerVector d_dims, d_j, d_p;
            Vector<RTYPE> xx ;
        };
        // Deprecated
        template<typename T>
        class Exporter<Eigen::MappedSparseMatrix<T, Eigen::RowMajor> > {
        public:
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
            Exporter(SEXP x) : d_x(x), d_dims(d_x.slot("Dim")), d_j(d_x.slot("j")), d_p(d_x.slot("p")), xx( d_x.slot("x") ) {
                if (!d_x.is("dgRMatrix"))
                    throw std::invalid_argument("Need S4 class dgRMatrix for a mapped sparse matrix");
            }
            Eigen::MappedSparseMatrix<T, Eigen::RowMajor> get() {
                return Eigen::MappedSparseMatrix<T, Eigen::RowMajor>(d_dims[0], d_dims[1], d_p[d_dims[1]],
                                                                     d_p.begin(), d_j.begin(), xx.begin() );
            }
        protected:
            S4            d_x;
            IntegerVector d_dims, d_j, d_p;
            Vector<RTYPE> xx ;
        };

        template<typename T>
        class Exporter<Eigen::SparseMatrix<T> > {
        public:
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;
            Exporter(SEXP x) : d_x(x), d_dims(d_x.slot("Dim")), d_i(d_x.slot("i")), d_p(d_x.slot("p")), xx(d_x.slot("x")) {
                if (!d_x.is("dgCMatrix"))
                    throw std::invalid_argument("Need S4 class dgCMatrix for a sparse matrix");
            }
            Eigen::SparseMatrix<T> get() {
                Eigen::SparseMatrix<T>  ans(d_dims[0], d_dims[1]);
                ans.reserve(d_p[d_dims[1]]);
                for(int j = 0; j < d_dims[1]; ++j) {
                    ans.startVec(j);
                    for (int k = d_p[j]; k < d_p[j + 1]; ++k) ans.insertBack(d_i[k], j) = xx[k];
                }
                ans.finalize();
                return ans;
            }
        protected:
            S4            d_x;
            IntegerVector d_dims, d_i, d_p;
            Vector<RTYPE> xx ;
        };

        template<typename T>
        class Exporter<Eigen::SparseMatrix<T, Eigen::RowMajor> > {
        public:
            const static int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype ;

            Exporter(SEXP x) : d_x(x), d_dims(d_x.slot("Dim")), d_j(d_x.slot("j")), d_p(d_x.slot("p")), xx(d_x.slot("x")) {
                if (!d_x.is("dgRMatrix"))
                    throw std::invalid_argument("Need S4 class dgRMatrix for a sparse matrix");
            }
            Eigen::SparseMatrix<T, Eigen::RowMajor> get() {
                Eigen::SparseMatrix<T, Eigen::RowMajor>  ans(d_dims[0], d_dims[1]);
                ans.reserve(d_p[d_dims[0]]);
                for(int i = 0; i < d_dims[0]; ++i) {
                    ans.startVec(i);
                    for (int k = d_p[i]; k < d_p[i + 1]; ++k) ans.insertBack(i, d_j[k]) = xx[k];
                }
                ans.finalize();
                return ans;
            }
        protected:
            S4            d_x;
            IntegerVector d_dims, d_j, d_p;
            Vector<RTYPE> xx ;
        };

    } // namespace traits
}

#endif
