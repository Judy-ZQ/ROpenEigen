#include "ROpenEigen_types.h"

template<typename _Scalar, int _AmbientDim>
Eigen::AlignedBox<_Scalar, _AmbientDim>* newAlignedBox(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::AlignedBox<_Scalar, _AmbientDim>* result=new Eigen::AlignedBox<_Scalar, _AmbientDim>();
      return result;
    }else if(param.length()==1){
      if(isItAlignedBox(param[0])){
        Eigen::AlignedBox<_Scalar, _AmbientDim> m = param[0];
        Eigen::AlignedBox<_Scalar, _AmbientDim>* result=new Eigen::AlignedBox<_Scalar, _AmbientDim>(m);
        return result;
      }else if(Rcpp::internal::is__simple<int>(param[0]) || Rcpp::internal::is__simple<double>(param[0])){
        int _dim = param[0];
        Eigen::AlignedBox<_Scalar, _AmbientDim>* result=new Eigen::AlignedBox<_Scalar, _AmbientDim>(_dim);
        return result;
      }
    }else if(param.length()==2){
      if(isItVector(param[0]) && isItVector(param[1])){
        Eigen::Matrix<_Scalar,Eigen::Dynamic,1> _min = Rcpp::as<Eigen::Matrix<_Scalar,Eigen::Dynamic,1>>(param[0]);
        Eigen::Matrix<_Scalar,Eigen::Dynamic,1> _max = Rcpp::as<Eigen::Matrix<_Scalar,Eigen::Dynamic,1>>(param[1]);
        Eigen::AlignedBox<_Scalar, _AmbientDim>* result=new Eigen::AlignedBox<_Scalar, _AmbientDim>(_min,_max);
        return result;
      }
    }else{
      throw std::runtime_error("parameter error.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
SEXP eigenAlignedBox_center(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> res = obj->center();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
Eigen::AlignedBox<_Scalar,_AmbientDim> eigenAlignedBox_clamp(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj,SEXP p){
  try{
    if(isItAlignedBox(p)){
      Eigen::AlignedBox<_Scalar,Eigen::Dynamic> _p = Rcpp::as<Eigen::AlignedBox<_Scalar,Eigen::Dynamic>>(p);
      Eigen::AlignedBox<_Scalar,_AmbientDim> res = obj->clamp(_p);
      return res;
    }else{
      throw std::runtime_error("parameter should be an AlignedBox.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
bool eigenAlignedBox_contains(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj,SEXP p){
  try{
    if(isItMatrix(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> _p = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>(p);
      return obj->contains(_p);
    }else if(isItVector(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> _p = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>>(p);
      return obj->contains(_p);
    }else if(isItAlignedBox(p)){
      Eigen::AlignedBox<_Scalar,Eigen::Dynamic> _p = Rcpp::as<Eigen::AlignedBox<_Scalar,Eigen::Dynamic>>(p);
      return obj->contains(_p);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _Scalar, int _AmbientDim>
// SEXP eigenAlignedBox_corner(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj,int cn){
//   try{
//     // typedef typename Eigen::AlignedBox<_Scalar,_AmbientDim>::CornerType CT;
//     // CT _corner = static_cast<CT>(cn);
//     // Eigen::AlignedBox<_Scalar,_AmbientDim>::CornerType _corner=static_cast<Eigen::AlignedBox<_Scalar,_AmbientDim>::CornerType>(cn);
//     // Eigen::Matrix<_Scalar, _AmbientDim, 1> res = obj->corner(_corner);
//     // Eigen::Matrix<_Scalar, _AmbientDim, 1> res = obj->corner(static_cast<Eigen::AlignedBox<_Scalar, _AmbientDim>::CornerType>(cn));
//     // using CornerType = typename Eigen::AlignedBox<_Scalar, _AmbientDim>::CornerType;
//     int d = obj->dim();
//     if(d==1){
//       // Eigen::AlignedBox<_Scalar,1> obj1(1);
//       using CornerType = typename Eigen::AlignedBox<_Scalar,1>::CornerType;
//       Eigen::Matrix<_Scalar, 1, 1> res(obj->corner(static_cast<CornerType>(cn)));
//       return Rcpp::wrap(res);
//     }else if(d==2){
//       // Eigen::AlignedBox<_Scalar,2> obj1(*obj);
//       using CornerType = typename Eigen::AlignedBox<_Scalar,2>::CornerType;
//       Eigen::Matrix<_Scalar, 2, 1> res(obj->corner(static_cast<CornerType>(cn)));
//       return Rcpp::wrap(res);
//     }else if(d==3){
//       // Eigen::AlignedBox<_Scalar,3> obj1(*obj);
//       using CornerType = typename Eigen::AlignedBox<_Scalar,3>::CornerType;
//       Eigen::Matrix<_Scalar, 3, 1> res(obj->corner(static_cast<CornerType>(cn)));
//       return Rcpp::wrap(res);
//     }else{
//       throw std::runtime_error("THIS_METHOD_IS_ONLY_FOR_VECTORS_OF_A_SPECIFIC_SIZE(1-3).");
//     }
//     // using CornerType = typename Eigen::AlignedBox<_Scalar, _AmbientDim>::CornerType;
//     // Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> res(obj->corner(static_cast<CornerType>(cn)));
//     // return Rcpp::wrap(res);
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _AmbientDim>
SEXP eigenAlignedBox_diagonal(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> res = obj->diagonal();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
int eigenAlignedBox_dim(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    return obj->dim();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
void eigenAlignedBox_extend(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj,SEXP p){
  try{
    if(isItMatrix(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> _p = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>(p);
      obj->extend(_p);
    }else if(isItVector(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> _p = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>>(p);
      obj->extend(_p);
    }else if(isItAlignedBox(p)){
      Eigen::AlignedBox<_Scalar,Eigen::Dynamic> _p = Rcpp::as<Eigen::AlignedBox<_Scalar,Eigen::Dynamic>>(p);
      obj->extend(_p);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
double eigenAlignedBox_exteriorDistance(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj,SEXP p){
  try{
    // if(isItMatrix(p)){
    //   Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> _p = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>(p);
    //   return obj->exteriorDistance(_p);
    // }else
    if(isItVector(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> _p = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>>(p);
      return obj->exteriorDistance(_p);
    }else if(isItAlignedBox(p)){
      Eigen::AlignedBox<_Scalar,Eigen::Dynamic> _p = Rcpp::as<Eigen::AlignedBox<_Scalar,Eigen::Dynamic>>(p);
      return obj->exteriorDistance(_p);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
Eigen::AlignedBox<_Scalar,_AmbientDim> eigenAlignedBox_intersection(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj,SEXP p){
  try{
    if(isItAlignedBox(p)){
      Eigen::AlignedBox<_Scalar,Eigen::Dynamic> _p = Rcpp::as<Eigen::AlignedBox<_Scalar,Eigen::Dynamic>>(p);
      Eigen::AlignedBox<_Scalar,_AmbientDim> res = obj->intersection(_p);
      return res;
    }else{
      throw std::runtime_error("parameter should be an AlignedBox.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
bool eigenAlignedBox_intersects(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj,SEXP p){
  try{
    if(isItAlignedBox(p)){
      Eigen::AlignedBox<_Scalar,Eigen::Dynamic> _p = Rcpp::as<Eigen::AlignedBox<_Scalar,Eigen::Dynamic>>(p);
      bool res = obj->intersects(_p);
      return res;
    }else{
      throw std::runtime_error("parameter should be an AlignedBox.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
bool eigenAlignedBox_isApprox(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj,SEXP p){
  try{
    Rcpp::List param(p);
    if(isItAlignedBox(param[0])){
      Eigen::AlignedBox<_Scalar,Eigen::Dynamic> _p = Rcpp::as<Eigen::AlignedBox<_Scalar,Eigen::Dynamic>>(param[0]);
      bool res;
      if(param.length()==1){
        res = obj->isApprox(_p);
      }else if(param.length()==2){
        double epsilon = param[1];
        res = obj->isApprox(_p,epsilon);
      }
      return res;
    }else{
      throw std::runtime_error("parameter should be an AlignedBox and a (optinal) numeric data.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
bool eigenAlignedBox_isEmpty(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    return obj->isEmpty();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
bool eigenAlignedBox_isNull(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    return obj->isNull();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
SEXP eigenAlignedBox_max(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> res = obj->max();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
SEXP eigenAlignedBox_min(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> res = obj->min();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
SEXP eigenAlignedBox_sample(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> res = obj->sample();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
void eigenAlignedBox_setEmpty(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    obj->setEmpty();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
void eigenAlignedBox_setNull(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    obj->setNull();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
SEXP eigenAlignedBox_sizes(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> res = obj->sizes();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
_Scalar eigenAlignedBox_squaredExteriorDistance(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj,SEXP p){
  try{
    // if(isItMatrix(p)){
    //   Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> _p = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>(p);
    //   return obj->squaredExteriorDistance(_p);
    // }else
    if(isItVector(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> _p = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>>(p);
      return obj->squaredExteriorDistance(_p);
    }else if(isItAlignedBox(p)){
      Eigen::AlignedBox<_Scalar,Eigen::Dynamic> _p = Rcpp::as<Eigen::AlignedBox<_Scalar,Eigen::Dynamic>>(p);
      return obj->squaredExteriorDistance(_p);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _AmbientDim>
_Scalar eigenAlignedBox_volume(Eigen::AlignedBox<_Scalar,_AmbientDim>* obj){
  try{
    return obj->volume();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// template<typename _Scalar, int _Dim, int _Mode, int _Options>
// Eigen::Transform< _Scalar, _Dim, _Mode, _Options >* newTransform(SEXP p){
//   try{
//     Rcpp::List param(p);
//     if(param.length()<=0){
//       Eigen::Transform< _Scalar, _Dim, _Mode, _Options >* result=new Eigen::Transform< _Scalar, _Dim, _Mode, _Options >();
//       return result;
//     }else if(param.length()==1){
//       if(isItTransform(param[0])){
//         Eigen::Transform< _Scalar, _Dim, _Mode, _Options > m = param[0];
//         Eigen::Transform< _Scalar, _Dim, _Mode, _Options >* result=new Eigen::Transform< _Scalar, _Dim, _Mode, _Options >(m);
//         return result;
//       }else if(isItMatrix(param[0])){
//         Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic> m = param[0];
//         Eigen::Transform< _Scalar, _Dim, _Mode, _Options >* result=new Eigen::Transform< _Scalar, _Dim, _Mode, _Options >(m);
//         return result;
//       }
//     }else{
//       throw std::runtime_error("parameter error.");
//     }
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

// template<typename _Scalar, int _Dim, int _Mode, int _Options>
// SEXP eigenTransform_matrix(Eigen::Transform< _Scalar, _Dim, _Mode, _Options >* obj){
//   try{
//     Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->matrix();
//     return Rcpp::wrap(res);
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Dim, int _Mode>
Eigen::Transform< _Scalar, _Dim, _Mode >* newTransform(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      // Rcout << "here\n";
      Eigen::Transform< _Scalar, _Dim, _Mode >* result=new Eigen::Transform< _Scalar, _Dim, _Mode >();
      // Rcout << "there\n";
      return result;
    }else if(param.length()==1){
      if(isItTransform(param[0])){
        Eigen::Transform< _Scalar, _Dim, _Mode > m = param[0];
        Eigen::Transform< _Scalar, _Dim, _Mode >* result=new Eigen::Transform< _Scalar, _Dim, _Mode >(m);
        return result;
      }else if(isItMatrix(param[0])){
        // if(_mode==Eigen::Affine && _Dim==3){
          // Eigen::Matrix<_Scalar,4,4> m = param[0];
          // Eigen::Transform< _Scalar, _Dim, _Mode >* result=new Eigen::Transform< _Scalar, _Dim, _Mode >(m);
          // return result;
        // }
        Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic> m1 = param[0];
        // Rcout << m1 << endl;
        // Eigen::Matrix<_Scalar,_Dim,_Dim+1> m = param[0]; //why can not be Dynamic?
        Eigen::Matrix<_Scalar,_Dim,_Dim+1> m = m1;
        Eigen::Transform< _Scalar, _Dim, _Mode >* result=new Eigen::Transform< _Scalar, _Dim, _Mode >(m);
        // Rcout << m << endl;
        // Rcout << result->matrix() << endl;
        return result;
      }
    }else{
      throw std::runtime_error("parameter error.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}


// Eigen::Affine3d* newAffine3d(SEXP p){
//   try{
//     Rcpp::List param(p);
//     if(param.length()<=0){
//       Rcout << "here\n";
//       Eigen::Affine3d* result=new Eigen::Affine3d();
//       Rcout << "there\n";
//       return result;
//     }else if(param.length()==1){
//       if(isItTransform(param[0])){
//         Eigen::Affine3d m = param[0];
//         Eigen::Affine3d* result=new Eigen::Affine3d(m);
//         return result;
//       }else if(isItMatrix(param[0])){
//         // Eigen::Matrix4d m = Rcpp::as<Eigen::Matrix4d>(param[0]);
//         // Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> m = Rcpp::as<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>>(param[0]);
//         Eigen::Matrix<double,4,4> m = Rcpp::as<Eigen::Matrix<double,4,4>>(param[0]);
//         // Eigen::Matrix4d m = Eigen::Matrix4d::Identity();
//         Eigen::Affine3d* result=new Eigen::Affine3d(m);
//         return result;
//       }
//     }else{
//       throw std::runtime_error("parameter error.");
//     }
//   }catch(std::exception &e){
//     Rcpp::stop(e.what());
//   }
// }

template<typename _Scalar, int _Dim, int _Mode>
SEXP eigenTransform_affine(Eigen::Transform< _Scalar, _Dim, _Mode >* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->affine();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
SEXP eigenTransform_at(Eigen::Transform< _Scalar, _Dim, _Mode >* obj,Rcpp::NumericVector rng,Rcpp::NumericVector cng){
  try{
    if(rng.length()<0 || cng.length()<0){
      throw std::runtime_error("the length of row and col parameter should be greater than zero.");
    }
    Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> m(obj->matrix());
    Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> res(rng.length(),cng.length());
    // Eigen::Matrix<_Scalar, _Rows, _Cols> ori(*obj);
    for(R_xlen_t rpos=0;rpos<rng.length();rpos++){
      for(R_xlen_t cpos=0;cpos<cng.length();cpos++){
        res(rpos,cpos) = m((int)rng[rpos],(int)cng[cpos]);
      }
    }
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
int eigenTransform_cols(Eigen::Transform< _Scalar, _Dim, _Mode >* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> res(obj->matrix());
    return res.cols();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Dim, int _Mode>
SEXP eigenTransform_cv2r(Eigen::Transform< _Scalar, _Dim, _Mode >* obj){
  try{
    Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic> res(obj->matrix());

    bool needs_dim = res.ColsAtCompileTime != 1;
    R_xlen_t m = res.rows(), n = res.cols();
    if (needs_dim && (m > INT_MAX || n > INT_MAX)) {
      Rcpp::stop("array dimensions cannot exceed INT_MAX");
    }
    R_xlen_t size = m * n;

    SEXP ans = PROTECT(::Rcpp::wrap(res.data(),
                                    res.data() + size));
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
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
SEXP eigenTransform_getlinear(Eigen::Transform< _Scalar, _Dim, _Mode >* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->linear();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
void eigenTransform_setlinear(Eigen::Transform< _Scalar, _Dim, _Mode >* obj,SEXP p){
  try{
    if(isItMatrix(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> m1 = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>(p);
      Eigen::Matrix<_Scalar, _Dim, _Dim> m = m1;
      obj->linear() = m;
    }else{
      throw std::runtime_error("patameter error.");
    }

  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
void eigenTransform_makeAffine(Eigen::Transform< _Scalar, _Dim, _Mode >* obj){
  try{
    obj->makeAffine();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
SEXP eigenTransform_matrix(Eigen::Transform< _Scalar, _Dim, _Mode >* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->matrix();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
void eigenTransform_prerotate(Eigen::Transform< _Scalar, _Dim, _Mode >* obj,SEXP p){
  try{
    if(isItMatrix(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> rotm1 = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>(p);
      Eigen::Matrix<_Scalar, _Dim, _Dim> rot = rotm1;
      obj->prerotate(rot);
    }else{
      throw std::runtime_error("parameter error.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
void eigenTransform_prescale(Eigen::Transform< _Scalar, _Dim, _Mode >* obj,SEXP p){
  try{
    if(isItVector(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> sca1 = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>>(p);
      Eigen::Matrix<_Scalar, _Dim, 1> sca = sca1;
      obj->prescale(sca);
    }else{
      // throw std::runtime_error("parameter error.");
      _Scalar sca = Rcpp::as<_Scalar>(p);
      obj->prescale(sca);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
void eigenTransform_pretranslate(Eigen::Transform< _Scalar, _Dim, _Mode >* obj,SEXP p){
  try{
    if(isItTranslation(p)){
      // Eigen::Translation< _Scalar, _Dim > trans = Rcpp::as<Eigen::Translation< _Scalar, _Dim >>(p);
      // obj->pretranslate(trans);
    }else if(isItVector(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> trans1 = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>>(p);
      Eigen::Matrix<_Scalar, _Dim, 1> trans = trans1;
      obj->pretranslate(trans);
    }else{
      throw std::runtime_error("parameter error.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
SEXP eigenTransform_rightproduct(Eigen::Transform< _Scalar, _Dim, _Mode >* obj,SEXP p){
  try{
    if(isItTransform(p)){
      // Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> _p = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>>(p);
      // return Rcpp::wrap((*obj)*_p);
      Eigen::Transform< _Scalar, _Dim, _Mode > other = Rcpp::as<Eigen::Transform< _Scalar, _Dim, _Mode >>(p);
      Eigen::Transform< _Scalar, _Dim, _Mode > res = (*obj)*other;
      return Rcpp::wrap(res);
    }else if(isItTranslation(p)){
      Eigen::Translation< _Scalar, _Dim > other = Rcpp::as<Eigen::Translation< _Scalar, _Dim >>(p);
      Eigen::Transform< _Scalar, _Dim, _Mode > res = (*obj)*other;
      return Rcpp::wrap(res);
    }else if(isItVector(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> _p1 = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>>(p);
      Eigen::Matrix<_Scalar, _Dim, 1> _p = _p1;
      // Eigen::Matrix<_Scalar, _Dim, 1> res = (*obj)*(_p);
      Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> res = (*obj)*(_p);
      // Rcout << typeid(res) == typeid(Eigen::Matrix<_Scalar, _Dim, 1>) << endl;
      // Rcout << res << endl;
      return Rcpp::wrap(res);
    }
    // else if(isItRotation2D(p)){
    //   Eigen::Rotation2D< _Scalar > other = Rcpp::as<Eigen::Rotation2D< _Scalar >>(p);
    //   Eigen::Transform< _Scalar, _Dim, _Mode > res = (*obj)*other;
    //   return Rcpp::wrap(res);
    // }
    // else if(isItAngleAxis(p)){
    //   Eigen::AngleAxis< _Scalar > other = Rcpp::as<Eigen::AngleAxis< _Scalar >>(p);
    //   Eigen::Transform< _Scalar, _Dim, _Mode > res = (*obj)*other;
    //   return Rcpp::wrap(res);
    // }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
SEXP eigenTransform_rotation(Eigen::Transform< _Scalar, _Dim, _Mode >* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> res = obj->rotation();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// bool isInstanceOf(Base* obj, const std::type_info& type) {
//   return typeid(*obj) == type;
// }

template<typename _Scalar, int _Dim, int _Mode>
int eigenTransform_rows(Eigen::Transform< _Scalar, _Dim, _Mode >* obj){
  try{
    Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> res(obj->matrix());
    return res.rows();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Dim, int _Mode>
void eigenTransform_set(Eigen::Transform< _Scalar, _Dim, _Mode >* obj,Rcpp::NumericVector rng,Rcpp::NumericVector cng,Rcpp::NumericVector p2){
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
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
void eigenTransform_setIdentity(Eigen::Transform< _Scalar, _Dim, _Mode >* obj){
  try{
    obj->setIdentity();
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim, int _Mode>
void eigenTransform_shear(Eigen::Transform< _Scalar, _Dim, _Mode >* obj, double sx, double sy){
  try{
    _Scalar _sx = (_Scalar)sx;
    _Scalar _sy = (_Scalar)sy;
    obj->shear(_sx,_sy);
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

// [[Rcpp::export]]
Eigen::Transform<double, 2, Eigen::Affine > EigenTransformAffine2d_Identity(){
  return Eigen::Transform<double, 2, Eigen::Affine>::Identity();
}

// [[Rcpp::export]]
Eigen::Transform<double, 3, Eigen::Affine > EigenTransformAffine3d_Identity(){
  return Eigen::Transform<double, 3, Eigen::Affine>::Identity();
}


template<typename _Scalar, int _Dim>
Eigen::Translation< _Scalar, _Dim >* newTranslation(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::Translation< _Scalar, _Dim>* result=new Eigen::Translation< _Scalar, _Dim>();
      return result;
    }else if(param.length()==1){
      if(isItTranslation(param[0])){
        Eigen::Translation< _Scalar, _Dim > m = param[0];
        Eigen::Translation< _Scalar, _Dim >* result=new Eigen::Translation< _Scalar, _Dim >(m);
        return result;
      }else if(isItVector(param[0])){
        Eigen::Matrix<_Scalar,Eigen::Dynamic,1> m1 = param[0];
        Eigen::Matrix<_Scalar,_Dim,1> m = m1;
        Eigen::Translation< _Scalar, _Dim >* result=new Eigen::Translation< _Scalar, _Dim >(m);
        return result;
      }else if(isItRowVector(param[0])){
        Eigen::Matrix<_Scalar,1,Eigen::Dynamic> m1 = param[0];
        Eigen::Matrix<_Scalar,1,_Dim> m = m1;
        Eigen::Translation< _Scalar, _Dim >* result=new Eigen::Translation< _Scalar, _Dim >(m);
        return result;
      }
    }else{
      throw std::runtime_error("parameter error.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Dim>
_Scalar eigenTranslation_getx(Eigen::Translation< _Scalar, _Dim>* obj){
  try{
    return obj->x();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Dim>
void eigenTranslation_setx(Eigen::Translation< _Scalar, _Dim>* obj,_Scalar value){
  try{
    obj->x() = value;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Dim>
_Scalar eigenTranslation_gety(Eigen::Translation< _Scalar, _Dim>* obj){
  try{
    return obj->y();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Dim>
void eigenTranslation_sety(Eigen::Translation< _Scalar, _Dim>* obj,_Scalar value){
  try{
    obj->y() = value;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Dim>
_Scalar eigenTranslation_getz(Eigen::Translation< _Scalar, _Dim>* obj){
  try{
    return obj->z();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Dim>
void eigenTranslation_setz(Eigen::Translation< _Scalar, _Dim>* obj,_Scalar value){
  try{
    obj->z() = value;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar, int _Options>
Eigen::Quaternion< _Scalar, _Options >* newQuaternion(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::Quaternion< _Scalar, _Options>* result=new Eigen::Quaternion< _Scalar, _Options >();
      return result;
    }else if(param.length()==1){
      if(isItQuaternion(param[0])){
        Eigen::Quaternion< _Scalar, _Options > m = param[0];
        Eigen::Quaternion< _Scalar, _Options >* result=new Eigen::Quaternion< _Scalar, _Options >(m);
        return result;
      }else if(isItVector(param[0])){
        Eigen::Matrix<_Scalar,Eigen::Dynamic,1> m1 = param[0];
        Eigen::Matrix<_Scalar,4,1> m = m1;
        Eigen::Quaternion< _Scalar, _Options >* result=new Eigen::Quaternion< _Scalar, _Options >(m);
        return result;
      }else if(isItRowVector(param[0])){
        // Eigen::Matrix<_Scalar,1,Eigen::Dynamic> m1 = param[0];
        // Eigen::Matrix<_Scalar,1,4> m = m1;
        // Eigen::Quaternion< _Scalar, _Options >* result=new Eigen::Quaternion< _Scalar, _Options >(m);
        // return result;
      }else if(isItAngleAxis(param[0])){
        Eigen::AngleAxis< _Scalar > m = param[0];
        Eigen::Quaternion< _Scalar, _Options >* result=new Eigen::Quaternion< _Scalar, _Options >(m);
        return result;
      }
    }else if(param.length()==4){
      _Scalar _w = param[0];
      _Scalar _x = param[1];
      _Scalar _y = param[2];
      _Scalar _z = param[3];
      Eigen::Quaternion< _Scalar, _Options >* result=new Eigen::Quaternion< _Scalar, _Options >(_w,_x,_y,_z);
      return result;
    }else{
      throw std::runtime_error("parameter error.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar, int _Options >
_Scalar eigenQuaternion_norm(Eigen::Quaternion< _Scalar, _Options >* obj){
  try{
    return obj->norm();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
Eigen::AngleAxis< _Scalar >* newAngleAxis(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::AngleAxis< _Scalar>* result=new Eigen::AngleAxis< _Scalar>();
      return result;
    }else if(param.length()==1){
      if(isItAngleAxis(param[0])){
        Eigen::AngleAxis< _Scalar > m = param[0];
        Eigen::AngleAxis< _Scalar >* result=new Eigen::AngleAxis< _Scalar >(m);
        return result;
      }else if(isItVector(param[0])){
        Eigen::Matrix<_Scalar,Eigen::Dynamic,1> m1 = param[0];
        Eigen::Matrix<_Scalar,4,1> m = m1;
        Eigen::AngleAxis< _Scalar >* result=new Eigen::AngleAxis< _Scalar >(m);
        return result;
      }
      // else if(isItRowVector(param[0])){
      //   Eigen::Matrix<_Scalar,1,Eigen::Dynamic> m1 = param[0];
      //   Eigen::Matrix<_Scalar,1,4> m = m1;
      //   Eigen::AngleAxis< _Scalar >* result=new Eigen::AngleAxis< _Scalar >(m);
      //   return result;
      // }
    }else if(param.length()==2){
      _Scalar _angle = param[0];
      if(isItVector(param[1])){
        Eigen::Matrix<_Scalar,Eigen::Dynamic,1> m1 = param[1];
        Eigen::Matrix<_Scalar,3,1> m = m1;
        Eigen::AngleAxis< _Scalar >* result=new Eigen::AngleAxis< _Scalar >(_angle,m);
        return result;
      }else if(isItRowVector(param[1])){
        Eigen::Matrix<_Scalar,1,Eigen::Dynamic> m1 = param[1];
        Eigen::Matrix<_Scalar,1,3> m = m1;
        Eigen::AngleAxis< _Scalar >* result=new Eigen::AngleAxis< _Scalar >(_angle,m);
        return result;
      }
    }else{
      throw std::runtime_error("parameter error.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar>
_Scalar eigenAngleAxis_getAngle(Eigen::AngleAxis< _Scalar>* obj){
  try{
    return obj->angle();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
void eigenAngleAxis_setAngle(Eigen::AngleAxis< _Scalar>* obj,_Scalar angle){
  try{
    obj->angle() = angle;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
SEXP eigenAngleAxis_getAxis(Eigen::AngleAxis< _Scalar>* obj){
  try{
    Eigen::Matrix<_Scalar,3,1> res = obj->axis();
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
void eigenAngleAxis_setAxis(Eigen::AngleAxis< _Scalar>* obj,SEXP axis){
  try{
    if(isItVector(axis)){
      Eigen::Matrix<_Scalar,Eigen::Dynamic,1> m1 = Rcpp::as<Eigen::Matrix<_Scalar,Eigen::Dynamic,1>>(axis);
      Eigen::Matrix<_Scalar,3,1> m = m1;
      obj->axis() = m;
    }
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
SEXP eigenAngleAxis_toRotationMatrix(Eigen::AngleAxis< _Scalar>* obj){
  try{
    // Eigen::Matrix<_Scalar,2,2> res1 = obj->toRotationMatrix();
    Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic> res = obj->toRotationMatrix();
    // Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic> res = res1;c
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
Eigen::Rotation2D< _Scalar >* newRotation2D(SEXP p){
  try{
    Rcpp::List param(p);
    if(param.length()<=0){
      Eigen::Rotation2D< _Scalar>* result=new Eigen::Rotation2D< _Scalar>();
      return result;
    }else if(param.length()==1){
      if(isItRotation2D(param[0])){
        Eigen::Rotation2D< _Scalar > m = param[0];
        Eigen::Rotation2D< _Scalar >* result=new Eigen::Rotation2D< _Scalar >(m);
        return result;
      }else if(isItMatrix(param[0])){
        Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic> m1 = param[0];
        Eigen::Matrix<_Scalar,2,2> m = m1;
        Eigen::Rotation2D< _Scalar >* result=new Eigen::Rotation2D< _Scalar >(m);
        return result;
      }else{
        _Scalar a(param[0]);
        Eigen::Rotation2D< _Scalar >* result=new Eigen::Rotation2D< _Scalar >(a);
        return result;
      }
    }else{
      throw std::runtime_error("parameter error.");
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar>
SEXP eigenRotation2D_fromRotationMatrix(Eigen::Rotation2D< _Scalar>* obj,SEXP p){
  try{
    if(isItMatrix(p)){
      Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic> m1 = Rcpp::as<Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic>>(p);
      Eigen::Matrix<_Scalar,2,2> m = m1;
      Eigen::Rotation2D< _Scalar > result=obj->fromRotationMatrix(m);
      return Rcpp::wrap(result);
    }else{
      throw std::runtime_error("parameter error.");
    }

  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
_Scalar eigenRotation2D_getAngle(Eigen::Rotation2D< _Scalar>* obj){
  try{
    return obj->angle();
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
void eigenRotation2D_setAngle(Eigen::Rotation2D< _Scalar>* obj,_Scalar angle){
  try{
    obj->angle() = angle;
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

template<typename _Scalar>
SEXP eigenRotation2D_rightproduct(Eigen::Rotation2D< _Scalar >* obj,SEXP p){
  try{
    if(isItRotation2D(p)){
      Eigen::Rotation2D< _Scalar > other = Rcpp::as<Eigen::Rotation2D< _Scalar >>(p);
      Eigen::Rotation2D< _Scalar > res = (*obj)*other;
      return Rcpp::wrap(res);
    }else if(isItVector(p)){
      Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> _p1 = Rcpp::as<Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>>(p);
      Eigen::Matrix<_Scalar, 2, 1> _p = _p1;
      Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> res = (*obj)*(_p);
      return Rcpp::wrap(res);
    }
  }catch(std::exception &e){
    Rcpp::stop(e.what());
  }
}

template<typename _Scalar>
SEXP eigenRotation2D_toRotationMatrix(Eigen::Rotation2D< _Scalar>* obj){
  try{
    Eigen::Matrix<_Scalar,2,2> res1 = obj->toRotationMatrix();
    Eigen::Matrix<_Scalar,Eigen::Dynamic,Eigen::Dynamic> res = res1;
    return Rcpp::wrap(res);
  }catch(std::exception &e){
    std::string msg =
      "Unknown error.\n";
    Rcpp::stop(msg);
  }
}

// typedef Eigen::Transform<double, Eigen::Dynamic, Eigen::Affine> myEigenTransformAffineXd;
typedef Eigen::Transform<double, 3, Eigen::Affine> myEigenTransformAffine3d;
typedef Eigen::Transform<double, 2, Eigen::Affine> myEigenTransformAffine2d;
typedef Eigen::Translation<double, Eigen::Dynamic> myEigenTranslationXd;
typedef Eigen::Translation<double, 3> myEigenTranslation3d;
typedef Eigen::AngleAxis<double> myEigenAngleAxisd;
typedef Eigen::Quaternion<double> myEigenQuaterniond;
typedef Eigen::Rotation2D<double> myEigenRotation2Dd;

RCPP_MODULE(GeometryModule){
  class_<Eigen::AlignedBox<double, Eigen::Dynamic>>("EigenAlignedBoxXd")
  .factory<SEXP>(newAlignedBox<double, Eigen::Dynamic>)
  .method("center",&eigenAlignedBox_center<double, Eigen::Dynamic>)
  .method("clamp",&eigenAlignedBox_clamp<double, Eigen::Dynamic>)
  .method("contains",&eigenAlignedBox_contains<double, Eigen::Dynamic>)
  // .method("corner",&eigenAlignedBox_corner<double, Eigen::Dynamic>)
  .method("diagonal",&eigenAlignedBox_diagonal<double, Eigen::Dynamic>)
  .method("dim",&eigenAlignedBox_dim<double, Eigen::Dynamic>)
  .method("extend",&eigenAlignedBox_extend<double, Eigen::Dynamic>)
  .method("exteriorDistance",&eigenAlignedBox_exteriorDistance<double, Eigen::Dynamic>)
  .method("intersection",&eigenAlignedBox_intersection<double, Eigen::Dynamic>)
  .method("intersects",&eigenAlignedBox_intersects<double, Eigen::Dynamic>)
  .method("isApprox",&eigenAlignedBox_isApprox<double, Eigen::Dynamic>)
  .method("isEmpty",&eigenAlignedBox_isEmpty<double, Eigen::Dynamic>)
  .method("isNull",&eigenAlignedBox_isNull<double, Eigen::Dynamic>)
  .method("max",&eigenAlignedBox_max<double, Eigen::Dynamic>)
  .method("min",&eigenAlignedBox_min<double, Eigen::Dynamic>)
  .method("squaredExteriorDistance",&eigenAlignedBox_squaredExteriorDistance<double, Eigen::Dynamic>)
  .method("volume",&eigenAlignedBox_volume<double, Eigen::Dynamic>)

  ;

  // class_<Eigen::Affine3d>("EigenAffine3d")
  //   .factory<SEXP>(newAffine3d)
  // // .method("matrix",&eigenTransform_matrix<double, Eigen::Dynamic, Eigen::Affine>)
  // ;
  // class_<myEigenTransformAffineXd>("EigenTransformAffineXd")
  //   .factory<SEXP>(newTransform<double, Eigen::Dynamic, Eigen::Affine>)
  //   .method("matrix",&eigenTransform_matrix<double, Eigen::Dynamic, Eigen::Affine>)
  //   ;
  class_<myEigenTransformAffine3d>("EigenTransformAffine3d")
    .factory<SEXP>(newTransform<double, 3, Eigen::Affine>)
    .method("affine",&eigenTransform_affine<double, 3, Eigen::Affine>)
    .method("cv2r",&eigenTransform_cv2r<double, 3, Eigen::Affine>)
    .method("cols",&eigenTransform_cols<double, 3, Eigen::Affine>)
    .method("getlinear",&eigenTransform_getlinear<double, 3, Eigen::Affine>)
    .method("setlinear",&eigenTransform_setlinear<double, 3, Eigen::Affine>)
    .method("makeAffine",&eigenTransform_makeAffine<double, 3, Eigen::Affine>)
    .method("matrix",&eigenTransform_matrix<double, 3, Eigen::Affine>)
    .method("prerotate",&eigenTransform_prerotate<double, 3, Eigen::Affine>)
    .method("prescale",&eigenTransform_prescale<double, 3, Eigen::Affine>)
    .method("pretranslate",&eigenTransform_pretranslate<double, 3, Eigen::Affine>)
    .method("rightproduct",&eigenTransform_rightproduct<double, 3, Eigen::Affine>)
    .method("rotation",&eigenTransform_rotation<double, 3, Eigen::Affine>)
    .method("rows",&eigenTransform_rows<double, 3, Eigen::Affine>)
    .method("setIdentity",&eigenTransform_setIdentity<double, 3, Eigen::Affine>)
    // .method("shear",&eigenTransform_shear<double, 3, Eigen::Affine>)
       .method("at", &eigenTransform_at<double, 3, Eigen::Affine>)
       .method("set", &eigenTransform_set<double, 3, Eigen::Affine>)
  ;

  class_<myEigenTransformAffine2d>("EigenTransformAffine2d")
    .factory<SEXP>(newTransform<double, 2, Eigen::Affine>)
    .method("affine",&eigenTransform_affine<double, 2, Eigen::Affine>)
    .method("cols",&eigenTransform_cols<double, 2, Eigen::Affine>)
    .method("cv2r",&eigenTransform_cv2r<double, 2, Eigen::Affine>)
    .method("getlinear",&eigenTransform_getlinear<double, 2, Eigen::Affine>)
    .method("setlinear",&eigenTransform_setlinear<double, 2, Eigen::Affine>)
    .method("makeAffine",&eigenTransform_makeAffine<double, 2, Eigen::Affine>)
    .method("matrix",&eigenTransform_matrix<double, 2, Eigen::Affine>)
    .method("prerotate",&eigenTransform_prerotate<double, 2, Eigen::Affine>)
    .method("prescale",&eigenTransform_prescale<double, 2, Eigen::Affine>)
    .method("pretranslate",&eigenTransform_pretranslate<double, 2, Eigen::Affine>)
    .method("rightproduct",&eigenTransform_rightproduct<double, 2, Eigen::Affine>)
    .method("rotation",&eigenTransform_rotation<double, 2, Eigen::Affine>)
    .method("rows",&eigenTransform_rows<double, 2, Eigen::Affine>)
    .method("setIdentity",&eigenTransform_setIdentity<double, 2, Eigen::Affine>)
    .method("shear",&eigenTransform_shear<double, 2, Eigen::Affine>)
    .method("at", &eigenTransform_at<double, 2, Eigen::Affine>)
    .method("set", &eigenTransform_set<double, 2, Eigen::Affine>)
  ;

  class_<myEigenTranslationXd>("EigenTranslationXd")
    .factory<SEXP>(newTranslation<double, Eigen::Dynamic>)
    .method("getX", &eigenTranslation_getx<double, Eigen::Dynamic>)
    .method("getY", &eigenTranslation_gety<double, Eigen::Dynamic>)
    .method("getZ", &eigenTranslation_getz<double, Eigen::Dynamic>)
    .method("setX", &eigenTranslation_setx<double, Eigen::Dynamic>)
    .method("setY", &eigenTranslation_sety<double, Eigen::Dynamic>)
    .method("setZ", &eigenTranslation_setz<double, Eigen::Dynamic>)
  ;
  class_<myEigenTranslation3d>("EigenTranslation3d")
    .factory<SEXP>(newTranslation<double, 3>)
    .method("getX", &eigenTranslation_getx<double, 3>)
    .method("getY", &eigenTranslation_gety<double, 3>)
    .method("getZ", &eigenTranslation_getz<double, 3>)
    .method("setX", &eigenTranslation_setx<double, 3>)
    .method("setY", &eigenTranslation_sety<double, 3>)
    .method("setZ", &eigenTranslation_setz<double, 3>)
  ;
  class_<myEigenAngleAxisd>("EigenAngleAxisd")
    .factory<SEXP>(newAngleAxis<double>)
    .method("getAngle", &eigenAngleAxis_getAngle<double>)
    .method("getAxis", &eigenAngleAxis_getAxis<double>)
    .method("setAngle", &eigenAngleAxis_setAngle<double>)
    .method("setAxis", &eigenAngleAxis_setAxis<double>)
    .method("toRotationMatrix", &eigenAngleAxis_toRotationMatrix<double>)
  ;
  class_<myEigenQuaterniond>("EigenQuaterniond")
    .factory<SEXP>(newQuaternion<double,Eigen::AutoAlign>)
    .method("norm", &eigenQuaternion_norm<double,Eigen::AutoAlign>)
  ;
  class_<myEigenRotation2Dd>("EigenRotation2Dd")
    .factory<SEXP>(newRotation2D<double>)
    .method("fromRotationMatrix", &eigenRotation2D_fromRotationMatrix<double>)
    .method("getAngle", &eigenRotation2D_getAngle<double>)
    .method("setAngle", &eigenRotation2D_setAngle<double>)
    .method("rightproduct",&eigenRotation2D_rightproduct<double>)
    .method("toRotationMatrix", &eigenRotation2D_toRotationMatrix<double>)
  ;
}
