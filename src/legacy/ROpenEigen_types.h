#pragma once

#define R_NO_REMAP
#define STRICT_R_HEADERS

#ifndef ROpenEigen_TYPE_H_
#define ROpenEigen_TYPE_H_

#include "RcppCommon.h"
#include <cstdint>
#include <cstddef>
#include <string>

template<typename T>
class myStdVector : public std::vector<T>{
public:
  myStdVector():std::vector<T>(){
  }
  myStdVector(int n):std::vector<T>(n){
  }
  myStdVector(const myStdVector& vec):std::vector<T>(vec){
  }
  myStdVector(const std::vector<T>& vec):std::vector<T>(vec){
  }

  size_t size() const{
    // return this->size();
    return std::vector<T>::size();
  }

  size_t max_size() const{
    // return this->max_size();
    return std::vector<T>::max_size();
  }

  size_t capacity() const{
    // return this->capacity();
    return std::vector<T>::capacity();
  }

  bool empty() const{
    // return this->empty();
    return std::vector<T>::empty();
  }

  void reserve(size_t new_cap){
    // this->reserve(new_cap);
    std::vector<T>::reserve(new_cap);
  }

  void pop_back(){
    // this->pop_back();
    std::vector<T>::pop_back();
  }

  void clear(){
    // this->clear();
    std::vector<T>::clear();
  }

  const T& back() const{
    // return this->back();
    return std::vector<T>::back();
  }

  const T& front() const{
    // return this->front();
    return std::vector<T>::front();
  }

  const T& at1(size_t pos) const{
    try{
      // return this->at(pos);
      return std::vector<T>::at(pos);
    }catch(std::exception &e){
      std::string msg =
        "Unknown error.\n";
      Rcpp::stop(msg);
    }
  }

  T& at2(size_t pos){
    try{
      // return this->at(pos);
      return std::vector<T>::at(pos);
    }catch(std::exception &e){
      std::string msg =
        "Unknown error.\n";
      Rcpp::stop(msg);
    }
  }

  // SEXP getXPtr(std::string typ){
  //   Rcpp::XPtr<myStdVector<T>> ptr(this,true);
  //   ptr.attr("class") = Rcpp::CharacterVector::create("xptr-stdVectorOf"+typ);
  //   return ptr;
  // }

  // myStdVector<T>* getThis(){
  //   return this;
  // }
};

namespace Rcpp{
template<typename T> SEXP wrap(const myStdVector<T>&);
template<typename T> SEXP wrap(const std::vector<T>&);
}

#include <ROpenEigenForward.h>

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

#include <ROpenEigenWrap.h>



#pragma region stdVector

namespace Rcpp{
template<typename T> SEXP wrap(const std::vector<T>& data){
  const std::vector<T>* xp = new std::vector<T>(data);
  return Rcpp::internal::make_new_object(xp);
}
}

namespace Rcpp{
template<typename T> SEXP wrap(const myStdVector<T>& data){
  const myStdVector<T>* xp = new myStdVector<T>(data);
  return Rcpp::internal::make_new_object(xp);
}
}

namespace Rcpp {
namespace traits {
template <typename T> class Exporter<myStdVector<T> >{
public:
  typedef myStdVector<T> stdvec_t;
  // enum { RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype };
  // typedef typename Rcpp::traits::storage_type<RTYPE>::type value_t;
  Exporter( SEXP x ) : sp(x){}
  inline stdvec_t get(){
    myStdVector<T> result;
    bool isStdVector = false;
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
        // if(clsname.find("stdVectorOfVectorOf")!=string::npos){//有三重向量？
        //   // myStdVector<myStdVector<T>>* ptr = (myStdVector<myStdVector<T>>*) R_ExternalPtrAddr(objptr);
        //   // result = *ptr;
        //   isStdVector = false;
        //   // return *ptr;
        // }else
        if(clsname.find("stdVectorOf")!=string::npos){
          myStdVector<T>* ptr = (myStdVector<T>*) R_ExternalPtrAddr(objptr);
          result = *ptr;
          isStdVector = true;
          // return *ptr;
        }else{
          std::string msg =
            "Error converting object to std::vector.\n";
          Rcpp::stop(msg);
        }
      }
    }
    if(!isStdVector){
      throw std::runtime_error("Parameter maybe wrong, error in get a std::vector.");
    }
    return result;
  }
private:
  // Rcpp::Vector<RTYPE> vec;
  SEXP sp;
};
}
}


template <typename T>
bool isItStdVector(SEXP sp){
  bool isStdVector=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    // SEXP objptr = ev[".pointer"];
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isStdVector = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      // Rcout << "clasname: " << xp->name << endl;
      std::string clsname = xp->name;
      if(clsname.find("stdVectorOf")!=string::npos){
        // cv::Size_<T>* ptr = (cv::Size_<T>*) R_ExternalPtrAddr(objptr);
        // result = *ptr;
        isStdVector = true;
      }else{
        isStdVector = false;
      }
    }
  }
  return isStdVector;
}

template <typename T>
bool isItStdVector(SEXP sp,string typ){
  bool isStdVector=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    // SEXP objptr = ev[".pointer"];
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isStdVector = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      // Rcout << "clasname: " << xp->name << endl;
      std::string clsname = xp->name;
      if(clsname.find("stdVectorOf"+typ)!=string::npos){
        // cv::Size_<T>* ptr = (cv::Size_<T>*) R_ExternalPtrAddr(objptr);
        // result = *ptr;
        isStdVector = true;
      }else{
        isStdVector = false;
      }
    }
  }
  return isStdVector;
}

bool isItStdVector_notemp(SEXP sp,string typ);

template <typename T>
bool isItStdVectorOfStdVector(SEXP sp){
  bool isStdVector=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    // SEXP objptr = ev[".pointer"];
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isStdVector = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      // Rcout << "clasname: " << xp->name << endl;
      std::string clsname = xp->name;
      if(clsname.find("stdVectorOfVectorOf")!=string::npos){
        // cv::Size_<T>* ptr = (cv::Size_<T>*) R_ExternalPtrAddr(objptr);
        // result = *ptr;
        isStdVector = true;
      }else{
        isStdVector = false;
      }
    }
  }
  return isStdVector;
}

template <typename T>
bool isItStdVectorOfStdVector(SEXP sp,string typ){
  bool isStdVector=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    // SEXP objptr = ev[".pointer"];
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isStdVector = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      // Rcout << "clasname: " << xp->name << endl;
      std::string clsname = xp->name;
      if(clsname.find("stdVectorOfVectorOf"+typ)!=string::npos){
        // cv::Size_<T>* ptr = (cv::Size_<T>*) R_ExternalPtrAddr(objptr);
        // result = *ptr;
        isStdVector = true;
      }else{
        isStdVector = false;
      }
    }
  }
  return isStdVector;
}

template <typename _Tp>
std::vector<_Tp> RcppNumericVector2StlVector(Rcpp::NumericVector nv){
  std::vector<_Tp> res;
  for(int i=0;i<nv.size();i++){
    res.push_back(nv[i]);
  }
  return res;
}


template <typename T>
bool SEXP2stdVectorOfVector(SEXP sp,std::vector<std::vector<T>>& result){
  bool isVector=false;
  // cv::Size_<T> result;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP objptr = ev[".pointer"];
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isVector = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      // Rcout << "clasname: " << xp->name << endl;
      std::string clsname = xp->name;
      if(clsname.find("stdVectorOfVectorOf")!=string::npos){
        std::vector<std::vector<T>>* ptr = (std::vector<std::vector<T>>*) R_ExternalPtrAddr(objptr);
        result = *ptr;
        isVector = true;
      }
    }
  }
  return isVector;
}

template <typename T>
bool SEXP2stdVector(SEXP sp,std::vector<T>& result){
  bool isVector=false;
  // cv::Size_<T> result;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP objptr = ev[".pointer"];
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isVector = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      // Rcout << "clasname: " << xp->name << endl;
      std::string clsname = xp->name;
      if(clsname.find("stdVectorOf")!=string::npos){
        std::vector<T>* ptr = (std::vector<T>*) R_ExternalPtrAddr(objptr);
        result = *ptr;
        isVector = true;
      }else{
        isVector = false;
      }
    }
  }
  return isVector;
}

template <typename T>
void finalize_stdVector(myStdVector<T> *cap){
  delete cap;
  cap = NULL;
}

template <typename T>
Rcpp::XPtr<myStdVector<T>> stdVector_getXPtr(myStdVector<T> *obj,string typ){
  // Rcpp::XPtr<myStdVector<T>,Rcpp::PreserveStorage, finalize_stdVector, true> ptr(obj);
  // XPtrStdVector_ ptr(obj);
  Rcpp::XPtr<myStdVector<T>> ptr(obj);
  ptr.attr("class") = Rcpp::CharacterVector::create("xptr-stdVectorOf"+typ);
  return ptr;
}

bool isItXPtr2Anything(SEXP sp,string typ);

// bool isItStdVectorOfStdVectorOfPoint(SEXP sp){
//   bool isStdVector=false;
//   if(TYPEOF(sp) == S4SXP){
//     Rcpp::Environment ev(sp);
//     // SEXP objptr = ev[".pointer"];
//     SEXP clsxptr = ev.get(".cppclass");
//     if (TYPEOF(clsxptr) != EXTPTRSXP){
//       isStdVector = false;
//     }else{
//       XPtr<class_Base> xp(clsxptr);
//       // Rcout << "clasname: " << xp->name << endl;
//       std::string clsname = xp->name;
//       if(clsname.find("stdVectorOfVectorOf")!=string::npos){
//         // cv::Size_<T>* ptr = (cv::Size_<T>*) R_ExternalPtrAddr(objptr);
//         // result = *ptr;
//         isStdVector = true;
//       }else{
//         isStdVector = false;
//       }
//     }
//   }
//   return isStdVector;
// }

template <typename T>
myStdVector<T>* newstdVector(SEXP p){
  myStdVector<T>* result=new myStdVector<T>();
  try{
    Rcpp::List param(p);
    // Rcout << "what:" << param.length() <<  endl;
    if(param.length()<=0){
      // result = new cv::Point_<T>();
    }else if(param.length()==1){
      bool isStdVector = isItStdVector<T>(param[0]);
      if(isStdVector){
        myStdVector<T> p = param[0];
        result = new myStdVector<T>(p);
      }else{
        int n = param[0];
        result = new myStdVector<T>(n);
      }
      // Rcout << mat << endl;
    }
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }

  return result;
}

template <typename T>
void stdVectorOutToConsole(myStdVector<T>* obj){
  // Rcout << *obj << endl;
  for(size_t i=0;i<obj->size();i++){
    Rcout << obj->at1(i) << " ";
  }
  Rcout << endl;
}

template <typename T>
void stdVector_assign(myStdVector<T>* obj, SEXP data) {
  // const int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype ;
  // Rcpp::Vector<RTYPE> vec(data);
  // Rcpp::Vector<T> vec(data);
  Rcpp::ListOf<T> vec(data);
  obj->assign(vec.begin(), vec.end());
  // T vec(data);
  // std::vector<T> vec(data);
  // if(isItStdVectorOfStdVector<T>(data)){
  //   myStdVector<myStdVector<T>>* vec=new myStdVector<myStdVector<T>>();
  //   SEXP2stdVectorOfVector(data,*vec);
  //   obj->assign(vec->begin(), vec->end());
  // }else
  //   if(isItStdVector<T>(data)){
  //   myStdVector<T>* vec=new myStdVector<T>();
  //   SEXP2stdVector<T>(data,*vec);
  //   obj->assign(vec->begin(), vec->end());
  // }else{
  //   std::string msg =
  //     "Error assigning vector\n";
  //   Rcpp::stop(msg);
  // }

}

template <typename T>
void stdVector_insert(myStdVector<T>* obj, int position, SEXP data) {
  // std::vector<T>::iterator it = obj->begin() + position;
  // const int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype ;
  // Rcpp::Vector<RTYPE> vec(data);
  Rcpp::ListOf<T> vec(data);
  obj->insert(obj->begin()+position,vec.begin(), vec.end());
  // obj->insert(obj->begin()+position, vec.begin(), vec.end());
  // if(isItStdVectorOfStdVector<T>(data)){
  //   myStdVector<myStdVector<T>>* vec=new myStdVector<myStdVector<T>>();
  //   SEXP2stdVectorOfVector(data,*vec);
  //   obj->insert(obj->begin()+position,vec->begin(), vec->end());
  // }else if(isItStdVector<T>(data)){
  //   myStdVector<T>* vec=new myStdVector<T>();
  //   SEXP2stdVector<T>(data,*vec);
  //   obj->insert(obj->begin()+position,vec->begin(), vec->end());
  // }else{
  //   std::string msg =
  //     "Error insertng vector\n";
  //   Rcpp::stop(msg);
  // }
}

template <typename T>
SEXP stdVector_std2r( myStdVector<T>* obj ) {
  // return Rcpp::wrap( *obj );
  Rcpp::List nv(obj->size());
  for(size_t i=0;i<obj->size();i++){
    nv[i] = obj->at1(i);
  }
  return nv;
}

template <typename T>
void stdVector_set(myStdVector<T>* obj, int i, T value) {
  obj->at2( i ) = value;
}

// // Fix for C++11, where we cannot directly expose
// // member functions vec::resize and vec::push_back
template <typename T>
void stdVector_resize (myStdVector<T>* obj, int n) {
  obj->resize(n);
}


template <typename T>
void stdVector_push_back (myStdVector<T>* obj, const T& value) {
  obj->push_back(value);
}

#pragma endregion

#pragma region EigenMatrix
// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Matrix<_Scalar, _Rows, _Cols>* newEigenMatrix(SEXP p){
//   try{
//     Rcpp::List param(p);
//     if(param.length()<=0){
//       Eigen::Matrix<_Scalar, _Rows, _Cols>* result=new Eigen::Matrix<_Scalar, _Rows, _Cols>();
//       return result;
//     }else if(param.length()==1){
//       if(Rcpp::internal::is__simple<int>(param[0]) || Rcpp::internal::is__simple<float>(param[0]) || Rcpp::internal::is__simple<double>(param[0])){
//         Eigen::Matrix<_Scalar, _Rows, _Cols>* result=new Eigen::Matrix<_Scalar, _Rows, _Cols>(param[0]);
//         return result;
//       }
//     }
//   }catch(std::exception &e){
//     std::string msg =
//       "Unknown error.\n";
//     Rcpp::stop(msg);
//   }
// }


bool isItMatrix(SEXP sp);
bool isItVector(SEXP sp);
bool isItRowVector(SEXP sp);

bool isItArray(SEXP sp);
bool isItColArray(SEXP sp);
bool isItRowArray(SEXP sp);

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols>* newEigenArray(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::Array<_Scalar, _Rows, _Cols>* result=new Eigen::Array<_Scalar, _Rows, _Cols>();
      return result;
    }else if(param.length()==1){
      if(isItArray(param[0])){
        Eigen::Array<_Scalar, _Rows, _Cols> m = param[0];
        Eigen::Array<_Scalar, _Rows, _Cols>* result=new Eigen::Array<_Scalar, _Rows, _Cols>(m);
        return result;
      }
    }else if(param.length()==2){
      if(Rcpp::internal::is__simple<int>(param[0]) || Rcpp::internal::is__simple<int>(param[1]) || Rcpp::internal::is__simple<double>(param[0]) || Rcpp::internal::is__simple<double>(param[1])){
        int rows = param[0];
        int cols = param[1];
        Eigen::Array<_Scalar, _Rows, _Cols>* result=new Eigen::Array<_Scalar, _Rows, _Cols>(rows,cols);
        return result;
      }
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows>
Eigen::Array<_Scalar, _Rows, 1>* newEigenArray1d(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::Array<_Scalar, _Rows, 1>* result=new Eigen::Array<_Scalar, _Rows, 1>();
      return result;
    }else if(param.length()==1){
      if(isItColArray(param[0])){
        Eigen::Array<_Scalar, _Rows, 1> m = param[0];
        Eigen::Array<_Scalar, _Rows, 1>* result=new Eigen::Array<_Scalar, _Rows, 1>(m);
        return result;
      }else if(Rcpp::internal::is__simple<int>(param[0]) || Rcpp::internal::is__simple<double>(param[0])){
        int index = param[0];
        Eigen::Array<_Scalar, _Rows, 1>* result=new Eigen::Array<_Scalar, _Rows, 1>(index);
        return result;
      }
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Cols>
Eigen::Array<_Scalar, 1, _Cols>* newEigenRowArray(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::Array<_Scalar, 1, _Cols>* result=new Eigen::Array<_Scalar, 1, _Cols>();
      return result;
    }else if(param.length()==1){
      if(isItRowArray(param[0])){
        // Rcout << "here\n";
        Eigen::Array<_Scalar, 1, _Cols> m = param[0];
        // Rcout << "there: " << m << endl;
        Eigen::Array<_Scalar, 1, _Cols>* result=new Eigen::Array<_Scalar, 1, _Cols>(m);
        return result;
      }else if(Rcpp::internal::is__simple<int>(param[0]) || Rcpp::internal::is__simple<double>(param[0]) ){
        int index = param[0];
        Eigen::Array<_Scalar, 1, _Cols>* result=new Eigen::Array<_Scalar, 1, _Cols>(index);
        return result;
      }
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}


template<typename _Scalar, int _Rows, int _Cols>
SEXP eigenArray_cv2r(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  bool needs_dim = obj->ColsAtCompileTime != 1;
  R_xlen_t m = obj->rows(), n = obj->cols();
  if (needs_dim && (m > INT_MAX || n > INT_MAX)) {
    Rcpp::stop("array dimensions cannot exceed INT_MAX");
  }
  R_xlen_t size = m * n;

  SEXP ans = PROTECT(::Rcpp::wrap(obj->data(),
                                  obj->data() + size));
  if (needs_dim) {
    SEXP dd = PROTECT(::Rf_allocVector(INTSXP, 2));
    int *d = INTEGER(dd);
    d[0] = m;
    d[1] = n;
    ::Rf_setAttrib(ans, R_DimSymbol, dd);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return ans;
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenArray_r2cv(Eigen::Array<_Scalar, _Rows, _Cols>* obj,Rcpp::NumericMatrix mv){
  try{
    R_xlen_t m = obj->rows(), n = obj->cols();
    if (mv.nrow()!=m || mv.ncol()!=n) {
      Rcpp::stop("dimensions are not equal");
    }


    for(R_xlen_t rpos=0;rpos<m;rpos++){
      for(R_xlen_t cpos=0;cpos<n;cpos++){
        (*obj)( rpos,cpos ) = (_Scalar)mv(rpos,cpos);
      }
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }

}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenArray2d_coeff(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int rowId,int colId){
  try{
    return obj->coeff(rowId,colId);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenArray1d_coeff(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    return obj->coeff(index);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenArray2d_coeffRef(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int rowId,int colId){
  try{
    return obj->coeffRef(rowId,colId);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenArray1d_coeffRef(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    return obj->coeffRef(index);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenArray2d_conservativeResize(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int rowId,int colId){
  try{
    obj->conservativeResize(rowId,colId);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenArray1d_conservativeResize(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    obj->conservativeResize(index);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray2d_setConstant(Eigen::Array<_Scalar, _Rows, _Cols>* obj,_Scalar val){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->setConstant(val);
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray1d_setConstant(Eigen::Array<_Scalar, _Rows, _Cols>* obj,_Scalar val){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->setConstant(val);
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray2d_setZero(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->setZero();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray1d_setZero(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->setZero();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray2d_setOnes(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->setOnes();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray1d_setOnes(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->setOnes();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray2d_setRandom(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->setRandom();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray1d_setRandom(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->setRandom();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}


//inherited from EigenBase

template<typename _Scalar, int _Rows, int _Cols>
void eigenArray_set(Eigen::Array<_Scalar, _Rows, _Cols>* obj, Rcpp::NumericVector rng,Rcpp::NumericVector cng,Rcpp::NumericVector p2) {
  try{
    if(rng.length()<0 || cng.length()<0){
      throw std::runtime_error("the length of row and col parameter should be greater than zero.");
    }
    int valind=0;

    for(R_xlen_t rpos=0;rpos<rng.length();rpos++){
      for(R_xlen_t cpos=0;cpos<cng.length();cpos++){
        (*obj)( (int)rng[rpos],(int)cng[cpos] ) = (_Scalar)p2[valind];
        valind++;
      }
    }

  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_at(Eigen::Array<_Scalar, _Rows, _Cols>* obj,Rcpp::NumericVector rng,Rcpp::NumericVector cng){
  try{
    if(rng.length()<0 || cng.length()<0){
      throw std::runtime_error("the length of row and col parameter should be greater than zero.");
    }
    Eigen::Array<_Scalar, _Rows, _Cols> res(rng.length(),cng.length());
    // Eigen::Matrix<_Scalar, _Rows, _Cols> ori(*obj);
    for(R_xlen_t rpos=0;rpos<rng.length();rpos++){
      for(R_xlen_t cpos=0;cpos<cng.length();cpos++){
        res(rpos,cpos) = (*obj)((int)rng[rpos],(int)cng[cpos]);
      }
    }
    return res;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
int eigenArray_rows(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->rows();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
int eigenArray_cols(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->cols();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
int eigenArray_size(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->size();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenArray_OutToConsole(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  Rcout << *obj << endl;
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenArray_resize(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int rows,int cols){
  try{
    obj->resize(rows,cols);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_abs(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->abs();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_abs2(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->abs2();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_acos(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->acos();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_acosh(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->acosh();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_arg(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->arg();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_asin(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->asin();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_asinh(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->asinh();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_atan(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->atan();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_atanh(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->atanh();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_ceil(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->ceil();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_conj(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Array<_Scalar, _Rows, _Cols> m = obj->conj();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_cos(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->cos();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_cosh(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->cosh();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_cube(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->cube();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_digamma(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Array<_Scalar, _Rows, _Cols> m = obj->digamma();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_erf(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Array<_Scalar, _Rows, _Cols> m = obj->erf();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }
//
// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_erfc(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Array<_Scalar, _Rows, _Cols> m = obj->erfc();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_exp(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->exp();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_expm1(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->expm1();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_floor(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->floor();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_imag(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Array<_Scalar, _Rows, _Cols> m = obj->imag();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_inverse(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->inverse();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_isfinite(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Array<_Scalar, _Rows, _Cols> m = obj->isfinite();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_isinf(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Array<_Scalar, _Rows, _Cols> m = obj->isinf();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_isnan(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Array<_Scalar, _Rows, _Cols> m = obj->isnan();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_lgamma(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Array<_Scalar, _Rows, _Cols> m = obj->lgamma();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_log(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->log();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_log10(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->log10();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_log1p(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->log1p();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_log2(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->log2();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_logistic(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->logistic();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_ndtri(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Array<_Scalar, _Rows, _Cols> m = obj->ndtri();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_real(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->real();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_rint(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->rint();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_round(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->round();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_rsqrt(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->rsqrt();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_sign(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->sign();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_sin(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->sin();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_sinh(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->sinh();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_sqrt(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->sqrt();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_square(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->square();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_tan(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->tan();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_tanh(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->tanh();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

//inherited from arraybase
template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenArray_matrix(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->matrix();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

//block operation for vector(array?)
template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray1d_head(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int n){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->head(n);
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray1d_tail(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int n){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->tail(n);
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray1d_segment(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int start,int n){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->segment(start,n);
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_col(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    return obj->col(index);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_row(Eigen::Array<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    return obj->row(index);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenArray_sort(Eigen::Array<_Scalar, _Rows, _Cols>* obj){
  try{
    std::sort(obj->begin(),obj->end());
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Array<_Scalar, _Rows, _Cols> eigenArray_greaterthan(Eigen::Array<_Scalar, _Rows, _Cols>* obj,_Scalar a){
//   try{
//     return (*obj)>a;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }


template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols>* newEigenMatrix(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::Matrix<_Scalar, _Rows, _Cols>* result=new Eigen::Matrix<_Scalar, _Rows, _Cols>();
      return result;
    }else if(param.length()==1){
      // if(isItMatrix<Eigen::Matrix<_Scalar, _Rows, _Cols>>(param[0])){
      if(isItMatrix(param[0])){
        // Rcout << "here\n";
        Eigen::Matrix<_Scalar, _Rows, _Cols> m = param[0];
        // Rcout << "there: " << m << endl;
        Eigen::Matrix<_Scalar, _Rows, _Cols>* result=new Eigen::Matrix<_Scalar, _Rows, _Cols>(m);
        return result;
      }
      // else if(Rcpp::internal::is__simple<int>(param[0]) || Rcpp::internal::is__simple<double>(param[0]) ){
      //   int index = param[0];
      //   Eigen::Matrix<_Scalar, _Rows, _Cols>* result=new Eigen::Matrix<_Scalar, _Rows, _Cols>(index);
      //   return result;
      // }
    }else if(param.length()==2){
      if(Rcpp::internal::is__simple<int>(param[0]) || Rcpp::internal::is__simple<int>(param[1]) || Rcpp::internal::is__simple<double>(param[0]) || Rcpp::internal::is__simple<double>(param[1])){
        int rows = param[0];
        int cols = param[1];
        Eigen::Matrix<_Scalar, _Rows, _Cols>* result=new Eigen::Matrix<_Scalar, _Rows, _Cols>(rows,cols);
        return result;
      }
    }
  }catch(std::exception &e){
    // std::string msg =
    //   "Unknown error.\n";
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows>
Eigen::Matrix<_Scalar, _Rows, 1>* newEigenVector(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::Matrix<_Scalar, _Rows, 1>* result=new Eigen::Matrix<_Scalar, _Rows, 1>();
      return result;
    }else if(param.length()==1){
      if(isItVector(param[0])){
        // Rcout << "here\n";
        Eigen::Matrix<_Scalar, _Rows, 1> m = param[0];
        // Rcout << "there: " << m << endl;
        Eigen::Matrix<_Scalar, _Rows, 1>* result=new Eigen::Matrix<_Scalar, _Rows, 1>(m);
        return result;
      }else if(Rcpp::internal::is__simple<int>(param[0]) || Rcpp::internal::is__simple<double>(param[0])){
        int index = param[0];
        Eigen::Matrix<_Scalar, _Rows, 1>* result=new Eigen::Matrix<_Scalar, _Rows, 1>(index);
        return result;
      }
      // else if(Rcpp::internal::is__simple<IntegerVector>(param[0]) || Rcpp::internal::is__simple<NumericVector>(param[0])){
      //   Eigen::Matrix<_Scalar, _Rows, 1> m = param[0];
      //   Eigen::Matrix<_Scalar, _Rows, 1>* result=new Eigen::Matrix<_Scalar, _Rows, 1>(m);
      //   return result;
      // }
    }
  }catch(std::exception &e){
    // std::string msg =
    //   "Unknown error.\n";
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Cols>
Eigen::Matrix<_Scalar, 1, _Cols>* newEigenRowVector(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::Matrix<_Scalar, 1, _Cols>* result=new Eigen::Matrix<_Scalar, 1, _Cols>();
      return result;
    }else if(param.length()==1){
      if(isItRowVector(param[0])){
        // Rcout << "here\n";
        Eigen::Matrix<_Scalar, 1, _Cols> m = param[0];
        // Rcout << "there: " << m << endl;
        Eigen::Matrix<_Scalar, 1, _Cols>* result=new Eigen::Matrix<_Scalar, 1, _Cols>(m);
        return result;
      }else if(Rcpp::internal::is__simple<int>(param[0]) || Rcpp::internal::is__simple<double>(param[0]) ){
        int index = param[0];
        Eigen::Matrix<_Scalar, 1, _Cols>* result=new Eigen::Matrix<_Scalar, 1, _Cols>(index);
        return result;
      }
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_coeff(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int rowId,int colId){
  try{
    return obj->coeff(rowId,colId);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenVector_coeff(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    return obj->coeff(index);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

//inherited from Eigen::PlainObjectBase

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_coeffRef(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int rowId,int colId){
  try{
    return obj->coeffRef(rowId,colId);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenVector_coeffRef(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    return obj->coeffRef(index);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_conservativeResize(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int rowId,int colId){
  try{
    obj->conservativeResize(rowId,colId);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenVector_conservativeResize(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    obj->conservativeResize(index);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// SEXP eigenMatrix_packet(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int rowId,int colId){
//   try{
//     return obj->packet(rowId,colId);
//   }catch(std::exception &e){
//     std::string msg =
//       "Unknown error.\n";
//     Rcpp::stop(msg);
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_setConstant(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,_Scalar val){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->setConstant(val);
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenVector_setConstant(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,_Scalar val){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->setConstant(val);
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_setZero(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int rows,int cols){
//   try{
//     Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->setZero(rows,cols);
//     return m;
//   }catch(std::exception &e){
//     std::string msg =
//       "Unknown error.\n";
//     Rcpp::stop(msg);
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_setZero(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->setZero();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenVector_setZero(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->setZero();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_setOnes(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->setOnes();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenVector_setOnes(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->setOnes();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_setRandom(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->setRandom();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenVector_setRandom(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->setRandom();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_setIdentity(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->setIdentity();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

//inherited from Eigen::MatrixBase
template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_acosh(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->acosh();
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->unaryExpr([](double x) { return std::acosh(x); });
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_acosh1(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
//   try{
//     Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->acosh();
//     return m;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_adjoint(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->acosh();
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->adjoint();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows>
Eigen::Matrix<_Scalar, 1, _Rows> eigenVector_adjoint(Eigen::Matrix<_Scalar, _Rows, 1>* obj){
  try{
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->acosh();
    Eigen::Matrix<_Scalar, 1, _Rows> m = obj->adjoint();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Cols>
Eigen::Matrix<_Scalar, _Cols, 1> eigenRowVector_adjoint(Eigen::Matrix<_Scalar, 1, _Cols>* obj){
  try{
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->acosh();
    Eigen::Matrix<_Scalar, _Cols, 1> m = obj->adjoint();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_adjointInPlace(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    obj->adjointInPlace();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_applyOnTheLeft(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()==3){
      int p(param[0]);
      int q(param[1]);
      Eigen::JacobiRotation< _Scalar > J = Rcpp::as<Eigen::JacobiRotation< _Scalar >>(param[2]);
      obj->applyOnTheLeft(p,q,J);
    }else if(param.length()==1){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> m = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>(param[0]);
      obj->applyOnTheLeft(m);
    }
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_applyOnTheRight(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()==3){
      int p(param[0]);
      int q(param[1]);
      Eigen::JacobiRotation< _Scalar > J = Rcpp::as<Eigen::JacobiRotation< _Scalar >>(param[2]);
      obj->applyOnTheRight(p,q,J);
    }else if(param.length()==1){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> m = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>(param[0]);
      obj->applyOnTheRight(m);
    }
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// void eigenMatrix_applyHouseholderOnTheLeft(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,Eigen::Matrix<_Scalar, _Rows, _Cols> B){
//   try{
//     obj->applyOnTheLeft(B);
//   }catch(std::exception &e){
//     std::string msg =
//       "Unknown error.\n";
//     Rcpp::stop(msg);
//   }
// }

// template<typename _Scalar, int _Rows, int _Cols>
// void eigenMatrix_applyOnTheLeft(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,Eigen::Matrix<_Scalar, _Rows, _Cols> B){
//   try{
//     obj->applyOnTheLeft(B);
//   }catch(std::exception &e){
//     std::string msg =
//       "Unknown error.\n";
//     Rcpp::stop(msg);
//   }
// }
//
// template<typename _Scalar, int _Rows, int _Cols>
// void eigenMatrix_applyOnTheRight(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,Eigen::Matrix<_Scalar, _Rows, _Cols> B){
//   try{
//     obj->applyOnTheRight(B);
//   }catch(std::exception &e){
//     std::string msg =
//       "Unknown error.\n";
//     Rcpp::stop(msg);
//   }
// }

// template<typename _Scalar, int _Rows>
// Eigen::Matrix<_Scalar, _Rows, 1> eigenVector_cross(Eigen::Matrix<_Scalar, _Rows, 1>* obj,SEXP other){
//   try{
//     Eigen::Matrix<_Scalar, _Rows, 1> res;
//     if(isItVector(other)){
//       Eigen::Matrix<_Scalar, _Rows, 1> _other = Rcpp::as<Eigen::Matrix<_Scalar, _Rows, 1>>(other);
//       if(_other.rows()!=obj->rows()){
//         throw std::runtime_error("dimension not matched");
//       }
//       res = obj->cross(_other);
//     }else if(isItRowVector(other)){
//       Eigen::Matrix<_Scalar, 1, _Rows> _other = Rcpp::as<Eigen::Matrix<_Scalar, 1, _Rows>>(other);
//       if(_other.cols()!=obj->rows()){
//         throw std::runtime_error("dimension not matched");
//       }
//       res = obj->cross(_other);
//     }
//     return res;
//   }catch(std::exception &e){
//     // std::string msg =
//     //   "Unknown error.\n";
//     Rcpp::stop(e.what());
//   }
// }
//
// template<typename _Scalar, int _Cols>
// Eigen::Matrix<_Scalar, 1, _Cols> eigenRowVector_cross(Eigen::Matrix<_Scalar, 1, _Cols>* obj,SEXP other){
//   try{
//     Eigen::Matrix<_Scalar, 1, _Cols> res;
//     if(isItVector(other)){
//       Eigen::Matrix<_Scalar, _Cols, 1> _other = Rcpp::as<Eigen::Matrix<_Scalar, _Cols, 1>>(other);
//       if(_other.rows()!=obj->cols()){
//         throw std::runtime_error("dimension not matched");
//       }
//       res = obj->cross(_other);
//     }else if(isItRowVector(other)){
//       Eigen::Matrix<_Scalar, 1, _Cols> _other = Rcpp::as<Eigen::Matrix<_Scalar, 1, _Cols>>(other);
//       if(_other.cols()!=obj->cols()){
//         throw std::runtime_error("dimension not matched");
//       }
//       res = obj->cross(_other);
//     }
//     return res;
//   }catch(std::exception &e){
//     // std::string msg =
//     //   "Unknown error.\n";
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Array<_Scalar, _Rows, _Cols> eigenMatrix_array(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Array<_Scalar, _Rows, _Cols> m = obj->array();
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
SEXP eigenMatrix_asDiagonal(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> m = obj->asDiagonal();
    return Rcpp::wrap(m);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows>
Eigen::Matrix<std::complex<_Scalar>, _Rows, _Rows> eigenMatrixCPL_asDiagonal(Eigen::Matrix<std::complex<_Scalar>, _Rows, 1>* obj){
  try{
    Eigen::Matrix<std::complex<_Scalar>, _Rows, _Rows> m = obj->asDiagonal();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_asinh(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->unaryExpr([](double x) { return std::asinh(x); });
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->asinh();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_atanh(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->unaryExpr([](double x) { return std::atanh(x); });
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->asinh();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
double eigenMatrix_blueNorm(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    double m = obj->blueNorm();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::ColPivHouseholderQR<Eigen::Matrix<_Scalar, _Rows, _Cols>> eigenMatrix_colPivHouseholderQr(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->colPivHouseholderQr();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::CompleteOrthogonalDecomposition<Eigen::Matrix<_Scalar, _Rows, _Cols>> eigenMatrix_completeOrthogonalDecomposition(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->completeOrthogonalDecomposition();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_computeInverseAndDetWithCheck01(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,SEXP inverse,SEXP determinant,SEXP invertible){
  try{
    if(!Rcpp::is<S4>(inverse)){
      throw std::runtime_error("parameter inverse should be a reference class");
    }
    Rcpp::Reference tmp_inverse(inverse);
    Rcpp::List tmpls_inverse = tmp_inverse.field("name");
    RObject p_inverse(tmpls_inverse[0]);
    Eigen::Matrix<_Scalar, _Rows, _Cols> _inverse;

    if(!Rcpp::is<S4>(determinant)){
      throw std::runtime_error("parameter determinant should be a reference class");
    }
    Rcpp::Reference tmp_determinant(determinant);
    Rcpp::List tmpls_determinant = tmp_determinant.field("name");
    RObject p_determinant(tmpls_determinant[0]);
    double _determinant = 0;

    if(!Rcpp::is<S4>(invertible)){
      throw std::runtime_error("parameter invertible should be a reference class");
    }
    Rcpp::Reference tmp_invertible(invertible);
    Rcpp::List tmpls_invertible = tmp_invertible.field("name");
    RObject p_invertible(tmpls_invertible[0]);
    bool _invertible;

    obj->computeInverseAndDetWithCheck(_inverse,_determinant,_invertible);
    tmpls_inverse[0]=_inverse;
    tmpls_determinant[0]=_determinant;
    tmpls_invertible[0]=_invertible;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_computeInverseAndDetWithCheck02(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,SEXP inverse,SEXP determinant,SEXP invertible,double epsilon){
  try{
    if(!Rcpp::is<S4>(inverse)){
      throw std::runtime_error("parameter inverse should be a reference class");
    }
    Rcpp::Reference tmp_inverse(inverse);
    Rcpp::List tmpls_inverse = tmp_inverse.field("name");
    RObject p_inverse(tmpls_inverse[0]);
    Eigen::Matrix<_Scalar, _Rows, _Cols> _inverse;

    if(!Rcpp::is<S4>(determinant)){
      throw std::runtime_error("parameter determinant should be a reference class");
    }
    Rcpp::Reference tmp_determinant(determinant);
    Rcpp::List tmpls_determinant = tmp_determinant.field("name");
    RObject p_determinant(tmpls_determinant[0]);
    double _determinant = 0;

    if(!Rcpp::is<S4>(invertible)){
      throw std::runtime_error("parameter invertible should be a reference class");
    }
    Rcpp::Reference tmp_invertible(invertible);
    Rcpp::List tmpls_invertible = tmp_invertible.field("name");
    RObject p_invertible(tmpls_invertible[0]);
    bool _invertible;

    obj->computeInverseAndDetWithCheck(_inverse,_determinant,_invertible,epsilon);
    tmpls_inverse[0]=_inverse;
    tmpls_determinant[0]=_determinant;
    tmpls_invertible[0]=_invertible;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_computeInverseAndDetWithCheck(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()==3){
      eigenMatrix_computeInverseAndDetWithCheck01(obj,param[0],param[1],param[2]);
    }else if(param.length()==4){
      return eigenMatrix_computeInverseAndDetWithCheck02(obj,param[0],param[1],param[2],param[3]);
    }else{
      throw std::runtime_error("parameter number error");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }

}
template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_computeInverseWithCheck01(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,SEXP inverse,SEXP invertible){
  try{
    if(!Rcpp::is<S4>(inverse)){
      throw std::runtime_error("parameter inverse should be a reference class");
    }
    Rcpp::Reference tmp_inverse(inverse);
    Rcpp::List tmpls_inverse = tmp_inverse.field("name");
    RObject p_inverse(tmpls_inverse[0]);
    Eigen::Matrix<_Scalar, _Rows, _Cols> _inverse;

    if(!Rcpp::is<S4>(invertible)){
      throw std::runtime_error("parameter invertible should be a reference class");
    }
    Rcpp::Reference tmp_invertible(invertible);
    Rcpp::List tmpls_invertible = tmp_invertible.field("name");
    RObject p_invertible(tmpls_invertible[0]);
    bool _invertible;

    obj->computeInverseWithCheck(_inverse,_invertible);
    tmpls_inverse[0]=_inverse;
    tmpls_invertible[0]=_invertible;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_computeInverseWithCheck02(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,SEXP inverse,SEXP invertible,double epsilon){
  try{
    if(!Rcpp::is<S4>(inverse)){
      throw std::runtime_error("parameter inverse should be a reference class");
    }
    Rcpp::Reference tmp_inverse(inverse);
    Rcpp::List tmpls_inverse = tmp_inverse.field("name");
    RObject p_inverse(tmpls_inverse[0]);
    Eigen::Matrix<_Scalar, _Rows, _Cols> _inverse;

    if(!Rcpp::is<S4>(invertible)){
      throw std::runtime_error("parameter invertible should be a reference class");
    }
    Rcpp::Reference tmp_invertible(invertible);
    Rcpp::List tmpls_invertible = tmp_invertible.field("name");
    RObject p_invertible(tmpls_invertible[0]);
    bool _invertible;

    obj->computeInverseWithCheck(_inverse,_invertible,epsilon);
    tmpls_inverse[0]=_inverse;
    tmpls_invertible[0]=_invertible;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_computeInverseWithCheck(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()==2){
      eigenMatrix_computeInverseWithCheck01(obj,param[0],param[1]);
    }else if(param.length()==3){
      return eigenMatrix_computeInverseWithCheck02(obj,param[0],param[1],param[2]);
    }else{
      throw std::runtime_error("parameter number error");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }

}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_cos(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->unaryExpr([](double x) { return std::cos(x); });
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->asinh();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_cosh(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->unaryExpr([](double x) { return std::cosh(x); });
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->asinh();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows>
Eigen::Matrix<_Scalar, _Rows, 1> eigenVector_cross(Eigen::Matrix<_Scalar, _Rows, 1>* obj,SEXP other){
  try{
    Eigen::Matrix<_Scalar, _Rows, 1> res;
    if(isItVector(other)){
      Eigen::Matrix<_Scalar, _Rows, 1> _other = Rcpp::as<Eigen::Matrix<_Scalar, _Rows, 1>>(other);
      if(_other.rows()!=obj->rows() && _other.rows()!=3){
        throw std::runtime_error("dimension not matched");
      }
      res = obj->cross(_other);
    }else if(isItRowVector(other)){
      Eigen::Matrix<_Scalar, 1, _Rows> _other = Rcpp::as<Eigen::Matrix<_Scalar, 1, _Rows>>(other);
      if(_other.cols()!=obj->rows()){
        throw std::runtime_error("dimension not matched");
      }
      res = obj->cross(_other);
    }
    return res;
  }catch(std::exception &e){
    // std::string msg =
    //   "Unknown error.\n";
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Cols>
Eigen::Matrix<_Scalar, 1, _Cols> eigenRowVector_cross(Eigen::Matrix<_Scalar, 1, _Cols>* obj,SEXP other){
  try{
    Eigen::Matrix<_Scalar, 1, _Cols> res;
    if(isItVector(other)){
      Eigen::Matrix<_Scalar, _Cols, 1> _other = Rcpp::as<Eigen::Matrix<_Scalar, _Cols, 1>>(other);
      if(_other.rows()!=obj->cols()){
        throw std::runtime_error("dimension not matched");
      }
      res = obj->cross(_other);
    }else if(isItRowVector(other)){
      Eigen::Matrix<_Scalar, 1, _Cols> _other = Rcpp::as<Eigen::Matrix<_Scalar, 1, _Cols>>(other);
      if(_other.cols()!=obj->cols()){
        throw std::runtime_error("dimension not matched");
      }
      res = obj->cross(_other);
    }
    return res;
  }catch(std::exception &e){
    // std::string msg =
    //   "Unknown error.\n";
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_determinant(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->unaryExpr([](double x) { return std::cosh(x); });
    _Scalar m = obj->determinant();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_diagonal(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->unaryExpr([](double x) { return std::cosh(x); });
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->diagonal();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows>
_Scalar eigenVector_dot(Eigen::Matrix<_Scalar, _Rows, 1>* obj,SEXP other){
  try{
    _Scalar res;
    if(isItVector(other)){
      Eigen::Matrix<_Scalar, _Rows, 1> _other = Rcpp::as<Eigen::Matrix<_Scalar, _Rows, 1>>(other);
      if(_other.rows()!=obj->rows()){
        throw std::runtime_error("dimension not matched");
      }
      res = obj->dot(_other);
    }else if(isItRowVector(other)){
      Eigen::Matrix<_Scalar, 1, _Rows> _other = Rcpp::as<Eigen::Matrix<_Scalar, 1, _Rows>>(other);
      if(_other.cols()!=obj->rows()){
        throw std::runtime_error("dimension not matched");
      }
      res = obj->dot(_other);
    }
    return res;
  }catch(std::exception &e){
    // std::string msg =
    //   "Unknown error.\n";
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Cols>
_Scalar eigenRowVector_dot(Eigen::Matrix<_Scalar, 1, _Cols>* obj,SEXP other){
  try{
    _Scalar res;
    if(isItVector(other)){
      Eigen::Matrix<_Scalar, _Cols, 1> _other = Rcpp::as<Eigen::Matrix<_Scalar, _Cols, 1>>(other);
      if(_other.rows()!=obj->cols()){
        throw std::runtime_error("dimension not matched");
      }
      res = obj->dot(_other);
    }else if(isItRowVector(other)){
      Eigen::Matrix<_Scalar, 1, _Cols> _other = Rcpp::as<Eigen::Matrix<_Scalar, 1, _Cols>>(other);
      if(_other.cols()!=obj->cols()){
        throw std::runtime_error("dimension not matched");
      }
      res = obj->dot(_other);
    }
    return res;
  }catch(std::exception &e){
    // std::string msg =
    //   "Unknown error.\n";
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_exp(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    // Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->exp();
    Eigen::Matrix<_Scalar, _Rows, _Cols> m = obj->unaryExpr([](double x) { return std::exp(x); });
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<std::complex<double>, _Rows, 1> eigenMatrix_eigenvalues(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    // Eigen::SelfAdjointEigenSolver<Eigen::Matrix<_Scalar, _Rows, _Cols>> eigen_solver(*obj);
    // if (eigen_solver.info() != Eigen::Success) {
    //   throw std::runtime_error("Failed to compute eigenvalues.");
    // }
    // return obj->eigen_solver.eigenvalues();
    return obj->eigenvalues();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_hnormalized(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
//   try{
//     return obj->hnormalized();
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_homogeneous(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
//   try{
//     return obj->homogeneous();
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::BDCSVD<Eigen::Matrix<_Scalar, _Rows, _Cols>> eigenMatrix_bdcSvd(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,unsigned int computationOptions){
  try{
    return obj->bdcSvd(computationOptions);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::FullPivLU<Eigen::Matrix<_Scalar, _Rows, _Cols>> eigenMatrix_fullPivLu(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->fullPivLu();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::HouseholderQR<Eigen::Matrix<_Scalar, _Rows, _Cols>> eigenMatrix_householderQr(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->householderQr();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
double eigenMatrix_hypotNorm(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    double m = obj->hypotNorm();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_inverse(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->inverse();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::LDLT<Eigen::Matrix<_Scalar, _Rows, _Cols>,Eigen::Lower> eigenMatrix_ldlt(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->ldlt();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols,int p>
double eigenMatrix_lpNorm(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    // double m=-1;
    // switch(p){
    // case 1:
    //   m = obj->lpNorm<1>();
    //   break;
    // case 2:
    //   m = obj->lpNorm<2>();
    //   break;
    // default:
    //   m = obj->lpNorm<Eigen::Infinity>();
    // }
    // using lpnorm =  (*obj).lpNorm<p>;
    // double m = lpnorm();
    double m = obj->template lpNorm<p>();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
double eigenMatrix_norm(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    double m = obj->norm();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
double eigenMatrix_squaredNorm(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    double m = obj->squaredNorm();
    return m;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_trace(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->trace();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}




//inheirted from densebase
// typename DenseBase<Derived>::ColwiseReturnType;

// template<typename _Scalar, int _Rows, int _Cols>
// typename Eigen::Matrix<_Scalar, _Rows, _Cols>::ConstColwiseReturnType
// eigenMatrix_colwise(const Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
//   try{
//     return obj->colwise();
//     // return ConstColwiseReturnType(*obj);
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> myMatrixXd1;

// Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> eigenMatrixXd_colwise(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>* obj){
//   try{
//     return obj->colwise();
//     // return ConstColwiseReturnType(*obj);
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// typename Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::ConstColwiseReturnType
// Eigen::DenseBase<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>::ConstColwiseReturnType eigenMatrixXd_colwise(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>* obj){
//   try{
//     return obj->colwise();
//     // return ConstColwiseReturnType(*obj);
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_mean(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->mean();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_prod(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->prod();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Cols, _Rows> eigenMatrix_transpose(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Cols, _Rows> m = obj->transpose();
    // Rcout << m << endl;
    return m;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::Matrix<_Scalar, 1, _Rows> eigenVector_transpose(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
//   try{
//     // Rcout << "here\n";
//     Eigen::Matrix<_Scalar, 1, _Rows> m = obj->transpose();
//     // Rcout << m << endl;
//     return m;
//   }catch(std::exception &e){
//     std::string msg =
//       "Unknown error.\n";
//     Rcpp::stop(msg);
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_transposeInPlace(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    obj->transposeInPlace();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// void eigenVector_transposeInPlace(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
//   try{
//     obj->transposeInPlace();
//   }catch(std::exception &e){
//     std::string msg =
//       "Unknown error.\n";
//     Rcpp::stop(msg);
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_maxCoeff01(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->maxCoeff();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_maxCoeff02(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, SEXP index){
  try{
    if(!Rcpp::is<S4>(index)){
      throw std::runtime_error("parameter index should be a reference class");
    }
    Rcpp::Reference tmp_index(index);
    Rcpp::List tmpls_index = tmp_index.field("name");
    RObject p_index(tmpls_index[0]);
    int _index = 0;

    _Scalar res = obj->maxCoeff(&_index);
    tmpls_index[0]=_index;
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_maxCoeff03(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, SEXP rowId,SEXP colId){
  try{
    if(!Rcpp::is<S4>(rowId)){
      throw std::runtime_error("parameter rowId should be a reference class");
    }
    Rcpp::Reference tmp_rowId(rowId);
    Rcpp::List tmpls_rowId = tmp_rowId.field("name");
    RObject p_rowId(tmpls_rowId[0]);
    int _rowId = 0;

    if(!Rcpp::is<S4>(colId)){
      throw std::runtime_error("parameter colId should be a reference class");
    }
    Rcpp::Reference tmp_colId(colId);
    Rcpp::List tmpls_colId = tmp_colId.field("name");
    RObject p_colId(tmpls_colId[0]);
    int _colId = 0;

    _Scalar res = obj->maxCoeff(&_rowId,&_colId);
    tmpls_rowId[0]=_rowId;
    tmpls_colId[0]=_colId;
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_maxCoeff(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      return eigenMatrix_maxCoeff01(obj);
    }else if(param.length()==2){
      return eigenMatrix_maxCoeff03(obj,param[0],param[1]);
    }else{
      throw std::runtime_error("parameter number error");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }

}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_maxCoeff1(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      return eigenMatrix_maxCoeff01(obj);
    }else if(param.length()==1){
      return eigenMatrix_maxCoeff02(obj,param[0]);
    }else{
      throw std::runtime_error("parameter number error");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }

}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_minCoeff01(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->minCoeff();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_minCoeff02(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, SEXP index){
  try{
    if(!Rcpp::is<S4>(index)){
      throw std::runtime_error("parameter index should be a reference class");
    }
    Rcpp::Reference tmp_index(index);
    Rcpp::List tmpls_index = tmp_index.field("name");
    RObject p_index(tmpls_index[0]);
    int _index = 0;

    _Scalar res = obj->minCoeff(&_index);
    tmpls_index[0]=_index;
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_minCoeff03(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, SEXP rowId,SEXP colId){
  try{
    if(!Rcpp::is<S4>(rowId)){
      throw std::runtime_error("parameter rowId should be a reference class");
    }
    Rcpp::Reference tmp_rowId(rowId);
    Rcpp::List tmpls_rowId = tmp_rowId.field("name");
    RObject p_rowId(tmpls_rowId[0]);
    int _rowId = 0;

    if(!Rcpp::is<S4>(colId)){
      throw std::runtime_error("parameter colId should be a reference class");
    }
    Rcpp::Reference tmp_colId(colId);
    Rcpp::List tmpls_colId = tmp_colId.field("name");
    RObject p_colId(tmpls_colId[0]);
    int _colId = 0;

    _Scalar res = obj->minCoeff(&_rowId,&_colId);
    tmpls_rowId[0]=_rowId;
    tmpls_colId[0]=_colId;
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_minCoeff(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      return eigenMatrix_minCoeff01(obj);
    }else if(param.length()==2){
      return eigenMatrix_minCoeff03(obj,param[0],param[1]);
    }else{
      throw std::runtime_error("parameter number error");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }

}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_minCoeff1(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      return eigenMatrix_minCoeff01(obj);
    }else if(param.length()==1){
      return eigenMatrix_minCoeff02(obj,param[0]);
    }else{
      throw std::runtime_error("parameter number error");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }

}


template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> eigenMatrix_reshaped(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      return obj->reshaped();
    }else if(param.length()==2){
      int nr=param[0];
      int nc=param[1];
      return obj->reshaped(nr,nc);
    }
    // else if(param.length()==3){
    //   int nr=param[0];
    //   int nc=param[1];
    //   std::string rMajor=param[2];
    //   if(rMajor=="RowMajor"){
    //     return (*obj).reshaped<Eigen::RowMajor>(nr,nc);
    //   }
    // }
    // return obj->reshaped(nr,nc);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

//inherited from EigenBase

template<typename _Scalar, int _Rows, int _Cols>
int eigenMatrix_rows(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->rows();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
int eigenMatrix_cols(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->cols();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
int eigenMatrix_size(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->size();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_OutToConsole(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  Rcout << *obj << endl;
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_resize(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int rows,int cols){
  try{
    obj->resize(rows,cols);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenVector_resize(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    obj->resize(index);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}



// template<typename _Scalar, int _Rows, int _Cols>
// void eigenMatrix_set(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, int row, int col, _Scalar value) {
//   try{
//     obj->operator()( row,col ) = value;
//   }catch(std::exception &e){
//     std::string msg =
//       "Unknown error.\n";
//     Rcpp::stop(msg);
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_set(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, Rcpp::NumericVector rng,Rcpp::NumericVector cng,Rcpp::NumericVector p2) {
  try{
    if(rng.length()<0 || cng.length()<0){
      throw std::runtime_error("the length of row and col parameter should be greater than zero.");
    }
    int valind=0;

    for(R_xlen_t rpos=0;rpos<rng.length();rpos++){
      for(R_xlen_t cpos=0;cpos<cng.length();cpos++){
        (*obj)( (int)rng[rpos],(int)cng[cpos] ) = (_Scalar)p2[valind];
        valind++;
      }
    }
    // Rcpp::NumericVector nv(p2);
    // if(Rcpp::internal::is__simple<IntegerVector>(p2) || Rcpp::internal::is__simple<NumericVector>(p2)){
    //   Rcpp::NumericVector nv(p2);
    //   for(R_xlen_t rpos=0;rpos<rng.length();rpos++){
    //     for(R_xlen_t cpos=0;cpos<cng.length();cpos++){
    //       (*obj)( (int)rng[rpos],(int)cng[cpos] ) = (_Scalar)nv[valind];
    //       valind++;
    //     }
    //   }
    // }else if(Rcpp::internal::is__simple<ComplexVector>(p2)){
    //   Rcpp::ComplexVector nv(p2);
    //   for(R_xlen_t rpos=0;rpos<rng.length();rpos++){
    //     for(R_xlen_t cpos=0;cpos<cng.length();cpos++){
    //       (*obj)( (int)rng[rpos],(int)cng[cpos] ) = std::complex<_Scalar>(1, 2);
    //       valind++;
    //     }
    //   }
    // }

  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrixXc_set(Eigen::Matrix<std::complex<_Scalar>, _Rows, _Cols>* obj, Rcpp::NumericVector rng,Rcpp::NumericVector cng,ComplexVector p2) {
  try{
    if(rng.length()<0 || cng.length()<0){
      throw std::runtime_error("the length of row and col parameter should be greater than zero.");
    }
    int valind=0;
    Rcpp::NumericVector rv = Re(p2);
    Rcpp::NumericVector iv = Im(p2);
    for(R_xlen_t rpos=0;rpos<rng.length();rpos++){
      for(R_xlen_t cpos=0;cpos<cng.length();cpos++){
        (*obj)( (int)rng[rpos],(int)cng[cpos] ) = std::complex<_Scalar>((_Scalar)rv[valind],(_Scalar)iv[valind]);
        valind++;
      }
    }

  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<std::complex<_Scalar>, _Rows, _Cols> eigenMatrixXc_conjugate(Eigen::Matrix<std::complex<_Scalar>, _Rows, _Cols>* obj) {
  try{
    return obj->conjugate();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// SEXP eigenMatrix_at(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj, int row, int col) {
//   try{
//     return Rcpp::wrap(obj->operator()( row,col ));
//   }catch(std::exception &e){
//     std::string msg =
//       "Unknown error.\n";
//     Rcpp::stop(msg);
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_at(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,Rcpp::NumericVector rng,Rcpp::NumericVector cng){
  try{
    if(rng.length()<0 || cng.length()<0){
      throw std::runtime_error("the length of row and col parameter should be greater than zero.");
    }
    Eigen::Matrix<_Scalar, _Rows, _Cols> res(rng.length(),cng.length());
    // Eigen::Matrix<_Scalar, _Rows, _Cols> ori(*obj);
    for(R_xlen_t rpos=0;rpos<rng.length();rpos++){
      for(R_xlen_t cpos=0;cpos<cng.length();cpos++){
        res(rpos,cpos) = (*obj)((int)rng[rpos],(int)cng[cpos]);
      }
    }
    return res;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Rows, int _Cols>
SEXP eigenMatrix_cv2r(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  bool needs_dim = obj->ColsAtCompileTime != 1;
  R_xlen_t m = obj->rows(), n = obj->cols();
  if (needs_dim && (m > INT_MAX || n > INT_MAX)) {
    Rcpp::stop("array dimensions cannot exceed INT_MAX");
  }
  R_xlen_t size = m * n;

  SEXP ans = PROTECT(::Rcpp::wrap(obj->data(),
                                  obj->data() + size));
  if (needs_dim) {
    SEXP dd = PROTECT(::Rf_allocVector(INTSXP, 2));
    int *d = INTEGER(dd);
    d[0] = m;
    d[1] = n;
    ::Rf_setAttrib(ans, R_DimSymbol, dd);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return ans;
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_r2cv(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,Rcpp::NumericMatrix mv){
  try{
    R_xlen_t m = obj->rows(), n = obj->cols();
    if (mv.nrow()!=m || mv.ncol()!=n) {
      Rcpp::stop("dimensions are not equal");
    }

    // Rcout << m <<"," << n << "," << mv.nrow() <<"," << mv.ncol() << endl;

    for(R_xlen_t rpos=0;rpos<m;rpos++){
      for(R_xlen_t cpos=0;cpos<n;cpos++){
        // Rcout << rpos <<"," << cpos << endl;
        // Rcout << (*obj)(rpos,cpos) << "," << (_Scalar)mv(rpos,cpos) << endl;
        (*obj)( rpos,cpos ) = (_Scalar)mv(rpos,cpos);
      }
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }

}


template<typename _Scalar, int _Rows, int _Cols>
int eigenMatrix_nonZeros(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->nonZeros();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
int eigenMatrix_outerSize(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->outerSize();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_replicate(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int RowFactor, int ColFactor){
  try{
    return obj->replicate(RowFactor,ColFactor);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_reverse(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->reverse();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_reverseInPlace(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    obj->reverseInPlace();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::VectorwiseOp<Eigen::Matrix<_Scalar, _Rows, _Cols>, Eigen::Horizontal> eigenMatrix_rowwise(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
//   try{
//     return obj->rowwise();
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }


// template<typename _Scalar, int _Rows, int _Cols>
// SEXP eigenMatrix_rowwise(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
//   try{
//     return obj->rowwise();
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// template<typename _Scalar, int _Rows, int _Cols>
// Eigen::VectorwiseOp<Eigen::Matrix<_Scalar, _Rows, _Cols>, Eigen::Vertical>  eigenMatrix_colwise(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
//   try{
//     return obj->colwise();
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, 1, _Cols> eigenMatrix_colwise_Sum(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, 1, _Cols> res = obj->colwise().sum();
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows,1> eigenMatrix_rowwise_Sum(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    Eigen::Matrix<_Scalar, _Rows, 1> res = obj->rowwise().sum();
    // Rcout << res << endl;

    // bool needs_dim = res.ColsAtCompileTime != 1;
    // R_xlen_t m = res.rows(), n = res.cols();
    // if (needs_dim && (m > INT_MAX || n > INT_MAX)) {
    //   Rcpp::stop("array dimensions cannot exceed INT_MAX");
    // }
    // R_xlen_t size = m * n;
    // SEXP ans = PROTECT(::Rcpp::wrap(res.data(),
    //                                 res.data() + size));
    // if (needs_dim) {
    //   SEXP dd = PROTECT(::Rf_allocVector(INTSXP, 2));
    //   int *d = INTEGER(dd);
    //   d[0] = m;
    //   d[1] = n;
    //   ::Rf_setAttrib(ans, R_DimSymbol, dd);
    //   UNPROTECT(1);
    // }
    // UNPROTECT(1);
    // return ans;
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_colwise_Add(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,Eigen::Matrix<_Scalar, _Rows, 1> v){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> res = obj->colwise()+v;
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_rowwise_Add(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,Eigen::Matrix<_Scalar, 1, _Cols> v){
  try{
    Eigen::Matrix<_Scalar, _Rows, _Cols> res = obj->rowwise()+v;
    return res;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
_Scalar eigenMatrix_sum(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->sum();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

//block
template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_block(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int startRow,int startCol,int blockRows,int blockCols){
  try{
    return obj->block(startRow,startCol,blockRows,blockCols);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_col(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    return obj->col(index);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_row(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int index){
  try{
    return obj->row(index);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _Scalar, int _Rows, int _Cols>
// void eigenMatrix_assign(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,Eigen::Matrix<_Scalar, _Rows, _Cols> v){
//   try{
//     // return obj->row(index);
//     *obj = v;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_topLeftCorner(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int p,int q){
  try{
    return obj->topLeftCorner(p,q);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_bottomLeftCorner(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int p,int q){
  try{
    return obj->bottomLeftCorner(p,q);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_topRightCorner(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int p,int q){
  try{
    return obj->topRightCorner(p,q);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_bottomRightCorner(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int p,int q){
  try{
    return obj->bottomRightCorner(p,q);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_topRows(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int q){
  try{
    return obj->topRows(q);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_bottomRows(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int q){
  try{
    return obj->bottomRows(q);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_leftCols(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int q){
  try{
    return obj->leftCols(q);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_rightCols(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int q){
  try{
    return obj->rightCols(q);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_middleCols(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int i,int q){
  try{
    return obj->middleCols(i,q);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<_Scalar, _Rows, _Cols> eigenMatrix_middleRows(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj,int i,int q){
  try{
    return obj->middleRows(i,q);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<double, _Rows, _Cols> eigenMatrix_real(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->real();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
Eigen::Matrix<double, _Rows, _Cols> eigenMatrix_imag(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    return obj->imag();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Rows, int _Cols>
void eigenMatrix_sort(Eigen::Matrix<_Scalar, _Rows, _Cols>* obj){
  try{
    std::sort(obj->begin(),obj->end());
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}


template <typename T, typename value_type>
bool SEXP2EigenMatrix(SEXP object,T& result){
  // T result;
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
      if(clsname.find("EigenMatrixX")!=string::npos){
        T* ptr = (T*) R_ExternalPtrAddr(objptr);
        result = *ptr;
        // Rcout << result << endl;
        isMat = true;
      }else if(clsname.find("EigenVectorX")!=string::npos || clsname.find("EigenRowVectorXd")!=string::npos){
        Eigen::Matrix<value_type, Eigen::Dynamic, 1>* ptr = (Eigen::Matrix<value_type, Eigen::Dynamic, 1>*) R_ExternalPtrAddr(objptr);
        // Rcout << "here\n";
        // Rcout << *ptr << endl;
        result = *ptr;
        // Rcout << result << endl;
        isMat = true;
      }else if(clsname.find("EigenRowVectorX")!=string::npos){
        Eigen::Matrix<value_type, 1, Eigen::Dynamic>* ptr = (Eigen::Matrix<value_type, 1, Eigen::Dynamic>*) R_ExternalPtrAddr(objptr);
        result = *ptr;
        // Rcout << result << endl;
        isMat = true;
      }else if(clsname.find("EigenArrayXX")!=string::npos){
        T* ptr = (T*) R_ExternalPtrAddr(objptr);
        result = *ptr;
        // Rcout << result << endl;
        isMat = true;
      }else if(clsname.find("EigenArrayX")!=string::npos){
        Eigen::Array<value_type, Eigen::Dynamic, 1>* ptr = (Eigen::Array<value_type, Eigen::Dynamic, 1>*) R_ExternalPtrAddr(objptr);
        result = *ptr;
        isMat = true;
      }else if(clsname.find("EigenRowArrayX")!=string::npos){
        Eigen::Array<value_type, 1, Eigen::Dynamic>* ptr = (Eigen::Array<value_type, 1, Eigen::Dynamic>*) R_ExternalPtrAddr(objptr);
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
  return isMat;
}


#pragma endregion

#pragma region Vector

#pragma endregion


#pragma region EigenBase
// template<typename Derived>
// Derived eigenBase_Constant4Matrix(int rows,int cols,double value){
//   try{
//     typedef typename Eigen::internal::traits<Derived>::Scalar Scalar;
//     Scalar _val(value);
//     Derived res = Eigen::EigenBase<Derived>::Constant(rows,cols,_val);
//     return res;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }
#pragma endregion

#pragma region VectorwiseOp
// template<typename ExpressionType, int Direction>
// void eigenVectorwiseOp_sum(Eigen::VectorwiseOp< ExpressionType, Direction >* obj){
//   try{
//     auto res = obj->sum();
//     Rcout << res << endl;
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }
#pragma endregion

#pragma region Geometry
bool isItAlignedBox(SEXP sp);
bool isItTransform(SEXP sp);
bool isItTranslation(SEXP sp);
bool isItQuaternion(SEXP sp);
bool isItAngleAxis(SEXP sp);
bool isItRotation2D(SEXP sp);
#pragma endregion

#pragma region JacobiRotation
bool isItJacobiRotation(SEXP sp);
#pragma endregion


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> myMatrixXd;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> myMatrixXcd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> myVectorXd;

#endif
