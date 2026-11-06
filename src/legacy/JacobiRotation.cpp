#include "ROpenEigen_types.h"

template<typename _Scalar>
Eigen::JacobiRotation< _Scalar >* newJacobiRotation(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::JacobiRotation< _Scalar>* result=new Eigen::JacobiRotation< _Scalar>();
      return result;
    }else if(param.length()==2){
      _Scalar _c(param[0]);
      _Scalar _s(param[1]);
      Eigen::JacobiRotation< _Scalar>* result=new Eigen::JacobiRotation< _Scalar>(_c,_s);
      return result;
    }else{
      throw std::runtime_error("parameter error.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar>
SEXP eigenJacobiRotation_adjoint(Eigen::JacobiRotation< _Scalar>* obj){
  try{
    Eigen::JacobiRotation< _Scalar > result=obj->adjoint();
    return Rcpp::wrap(result);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
void eigenJacobiRotation_makeGivens01(Eigen::JacobiRotation< _Scalar>* obj, _Scalar p, _Scalar q){
  try{
    // if(!Rcpp::is<S4>(r)){
    //   throw std::runtime_error("parameter index should be a reference class");
    // }
    // Rcpp::Reference tmp_index(r);
    // Rcpp::List tmpls_r = tmp_index.field("name");
    // RObject p_index(tmpls_r[0]);
    // _Scalar _r = 0;

    obj->makeGivens(p,q);
    // tmpls_r[0]=_r;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar>
void eigenJacobiRotation_makeGivens02(Eigen::JacobiRotation< _Scalar>* obj, _Scalar p, _Scalar q, SEXP r){
  try{
    if(!Rcpp::is<S4>(r)){
      throw std::runtime_error("parameter index should be a reference class");
    }
    Rcpp::Reference tmp_index(r);
    Rcpp::List tmpls_r = tmp_index.field("name");
    RObject p_index(tmpls_r[0]);
    _Scalar _r = 0;

    obj->makeGivens(p,q,&_r);
    tmpls_r[0]=_r;
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar>
void eigenJacobiRotation_makeGivens(Eigen::JacobiRotation< _Scalar>* obj, SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()==2){
      _Scalar _p(param[0]);
      _Scalar _q(param[1]);
      eigenJacobiRotation_makeGivens01(obj,_p,_q);
    }else if(param.length()==3){
      _Scalar _p(param[0]);
      _Scalar _q(param[1]);
      eigenJacobiRotation_makeGivens02(obj,_p,_q,param[2]);
    }else{
      throw std::runtime_error("parameter number error");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }

}

template<typename _Scalar>
bool eigenJacobiRotation_makeJacobi(Eigen::JacobiRotation< _Scalar>* obj,SEXP p){
  try{
    Rcpp::List param(p);
    if(isItMatrix(param[0])){
      Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic> m1 = Rcpp::as<Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic>>(param[0]);
      Eigen::Matrix<_Scalar,2,2> m = m1;
      int p(param[1]);
      int q(param[2]);
      bool result=obj->makeJacobi(m,p,q);
      return result;
    }else{
      // typedef typename NumTraits<_Scalar>::Real RealScalar;
      // typedef NumTraits<_Scalar>::Real RealScalar;
      // using RealScalar = obj->template RealScalar;
      // obj->template RealScalar x(param[0]);
      using RealScalar = typename Eigen::JacobiRotation< _Scalar>::RealScalar;
      RealScalar x(param[0]);
      _Scalar y(param[1]);
      _Scalar z(param[2]);
      bool result=obj->makeJacobi(x,y,z);
      return result;
    }
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
SEXP eigenJacobiRotation_rightproduct(Eigen::JacobiRotation< _Scalar >* obj,SEXP p){
  try{
    if(isItJacobiRotation(p)){
      Eigen::JacobiRotation< _Scalar > other = Rcpp::as<Eigen::JacobiRotation< _Scalar >>(p);
      Eigen::JacobiRotation< _Scalar > res = (*obj)*other;
      return Rcpp::wrap(res);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar>
SEXP eigenJacobiRotation_transpose(Eigen::JacobiRotation< _Scalar>* obj){
  try{
    Eigen::JacobiRotation< _Scalar > result=obj->transpose();
    return Rcpp::wrap(result);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

typedef Eigen::JacobiRotation<double> myJacobiRotationd;

RCPP_MODULE(JacobiModule){
  class_<myJacobiRotationd>("EigenJacobiRotationd")
  .factory<SEXP>(newJacobiRotation<double>)
  .method("adjoint", &eigenJacobiRotation_adjoint<double>)
  .method("makeGivens_", &eigenJacobiRotation_makeGivens<double>)
  .method("makeJacobi_", &eigenJacobiRotation_makeJacobi<double>)
  .method("rightproduct", &eigenJacobiRotation_rightproduct<double>)
  .method("transpose", &eigenJacobiRotation_transpose<double>)
  ;
}
