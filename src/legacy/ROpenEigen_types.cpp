#include "ROpenEigen_types.h"

bool isItStdVector_notemp(SEXP sp,string typ){
  bool isStdVector=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isStdVector = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("stdVectorOf"+typ)!=string::npos){
        isStdVector = true;
      }else{
        isStdVector = false;
      }
    }
  }
  return isStdVector;
}



// [[Rcpp::export]]
bool isItXPtr2Anything(SEXP sp,string typ){
  bool isxptr=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isxptr = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      // Rcout << clsname << "\n";
      if(clsname.find("xptr-"+typ)!=string::npos){
        isxptr = true;
      }else{
        isxptr = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    XPtr<class_Base> xp(sp);
    std::string clsname = "xptr-"+typ;
    const char* cls = clsname.c_str();
    if(xp.inherits(cls)){
      isxptr = true;
    }else{
      isxptr = false;
    }
    // std::string clsname = xp->name;
    // // Rcout << clsname << "\n";
    // if(clsname.find("xptr-"+typ)!=string::npos){
    //   isxptr = true;
    // }else{
    //   isxptr = false;
    // }
  }
  return isxptr;
}

bool isItScalar_notemp(SEXP sp,string typ=""){
  bool isScalar=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    // SEXP objptr = ev[".pointer"];
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isScalar = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      // Rcout << "clasname: " << xp->name << endl;
      std::string clsname = xp->name;
      if(clsname.find("cvScalar"+typ)!=string::npos){
        // cv::Scalar_<T>* ptr = (cv::Scalar_<T>*) R_ExternalPtrAddr(objptr);
        // result = *ptr;
        isScalar = true;
      }else{
        isScalar = false;
      }
    }
  }
  return isScalar;
}

bool isItMatrix(SEXP sp){
  bool isMat=false;
  // cv::Mat result;
  // Rcout << "TYPE: " << TYPEOF(sp) << endl;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    // SEXP objptr = ev[".pointer"];
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isMat = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      // Rcout << "clasname: " << xp->name << endl;
      std::string clsname = xp->name;
      if(clsname.find("EigenMatrixX")!=string::npos){
        isMat = true;
      }else{
        isMat = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-mat")){
      isMat = false;
    }else{
      isMat = true;
    }
  }else if(TYPEOF(sp) == REALSXP || TYPEOF(sp) == INTSXP){
    Shield<SEXP> dims( ::Rf_getAttrib( sp, R_DimSymbol ) );
    if( Rf_isNull(dims) || ::Rf_length(dims) != 2 ){
      isMat = false;
    }else{
      isMat = true;
    }
  }
  return isMat;
}


bool isItVector(SEXP sp){
  bool isMat=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isMat = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("EigenVector")!=string::npos){
        isMat = true;
      }else{
        isMat = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-vector")){
      isMat = false;
    }else{
      isMat = true;
    }
  }else if(TYPEOF(sp) == REALSXP || TYPEOF(sp) == INTSXP){
    Shield<SEXP> dims( ::Rf_getAttrib( sp, R_DimSymbol ) );
    // if( Rf_isNull(dims) || ::Rf_length(dims) != 2 ){
    if( Rf_isNull(dims) && Rf_xlength(sp)>1){//is an R vector
      isMat = true;
    }else{
      isMat = false;
    }
  }
  return isMat;
}


bool isItRowVector(SEXP sp){
  bool isMat=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isMat = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("EigenRowVector")!=string::npos){
        isMat = true;
      }else{
        isMat = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-rowvector")){
      isMat = false;
    }else{
      isMat = true;
    }
  }else if(TYPEOF(sp) == REALSXP || TYPEOF(sp) == INTSXP){
    Shield<SEXP> dims( ::Rf_getAttrib( sp, R_DimSymbol ) );
    if( Rf_isNull(dims) || ::Rf_length(dims) != 2 ){
      isMat = false;
    }else{
      isMat = true;
    }
  }
  return isMat;
}

bool isItArray(SEXP sp){
  bool isMat=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isMat = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("EigenArrayXX")!=string::npos){
        isMat = true;
      }else{
        isMat = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-arrayxx")){
      isMat = false;
    }else{
      isMat = true;
    }
  }else if(TYPEOF(sp) == REALSXP || TYPEOF(sp) == INTSXP){
    Shield<SEXP> dims( ::Rf_getAttrib( sp, R_DimSymbol ) );
    if( Rf_isNull(dims) || ::Rf_length(dims) != 2 ){
      isMat = false;
    }else{
      isMat = true;
    }
  }
  return isMat;
}


bool isItColArray(SEXP sp){
  bool isMat=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isMat = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("EigenArrayX")!=string::npos){
        isMat = true;
      }else{
        isMat = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-arrayx")){
      isMat = false;
    }else{
      isMat = true;
    }
  }else if(TYPEOF(sp) == REALSXP || TYPEOF(sp) == INTSXP){
    Shield<SEXP> dims( ::Rf_getAttrib( sp, R_DimSymbol ) );
    // if( Rf_isNull(dims) || ::Rf_length(dims) != 2 ){
    if( Rf_isNull(dims) && Rf_xlength(sp)>1){//is an R vector
      isMat = true;
    }else{
      isMat = false;
    }
  }
  return isMat;
}


bool isItRowArray(SEXP sp){
  bool isMat=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isMat = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("EigenRowArrayX")!=string::npos){
        isMat = true;
      }else{
        isMat = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-rowarrayx")){
      isMat = false;
    }else{
      isMat = true;
    }
  }else if(TYPEOF(sp) == REALSXP || TYPEOF(sp) == INTSXP){
    Shield<SEXP> dims( ::Rf_getAttrib( sp, R_DimSymbol ) );
    if( Rf_isNull(dims) || ::Rf_length(dims) != 2 ){
      isMat = false;
    }else{
      isMat = true;
    }
  }
  return isMat;
}

#pragma region Geometry
bool isItAlignedBox(SEXP sp){
  bool isAlignedBox=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isAlignedBox = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("AlignedBoxX")!=string::npos){
        isAlignedBox = true;
      }else{
        isAlignedBox = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-alignedbox")){
      isAlignedBox = false;
    }else{
      isAlignedBox = true;
    }
  }
  return isAlignedBox;
}

bool isItTransform(SEXP sp){
  bool isTransform=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isTransform = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("Transform")!=string::npos){
        isTransform = true;
      }else{
        isTransform = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-transform")){
      isTransform = false;
    }else{
      isTransform = true;
    }
  }
  return isTransform;
}

bool isItTranslation(SEXP sp){
  bool isTranslation=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isTranslation = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("Translation")!=string::npos){
        isTranslation = true;
      }else{
        isTranslation = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-Translation")){
      isTranslation = false;
    }else{
      isTranslation = true;
    }
  }
  return isTranslation;
}

bool isItQuaternion(SEXP sp){
  bool isQuaternion=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isQuaternion = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("Quaternion")!=string::npos){
        isQuaternion = true;
      }else{
        isQuaternion = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-Quaternion")){
      isQuaternion = false;
    }else{
      isQuaternion = true;
    }
  }
  return isQuaternion;
}

bool isItAngleAxis(SEXP sp){
  bool isAngleAxis=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isAngleAxis = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("AngleAxis")!=string::npos){
        isAngleAxis = true;
      }else{
        isAngleAxis = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-AngleAxis")){
      isAngleAxis = false;
    }else{
      isAngleAxis = true;
    }
  }
  return isAngleAxis;
}

bool isItRotation2D(SEXP sp){
  bool isRotation2D=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isRotation2D = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("Rotation2D")!=string::npos){
        isRotation2D = true;
      }else{
        isRotation2D = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-Rotation2D")){
      isRotation2D = false;
    }else{
      isRotation2D = true;
    }
  }
  return isRotation2D;
}
#pragma endregion

#pragma region JacobiRotation
bool isItJacobiRotation(SEXP sp){
  bool isJacobiRotation=false;
  if(TYPEOF(sp) == S4SXP){
    Rcpp::Environment ev(sp);
    SEXP clsxptr = ev.get(".cppclass");
    if (TYPEOF(clsxptr) != EXTPTRSXP){
      isJacobiRotation = false;
    }else{
      XPtr<class_Base> xp(clsxptr);
      std::string clsname = xp->name;
      if(clsname.find("JacobiRotation")!=string::npos){
        isJacobiRotation = true;
      }else{
        isJacobiRotation = false;
      }
    }
  }else if(TYPEOF(sp) == EXTPTRSXP){
    //是mat的xptr指针
    XPtr<class_Base> xp(sp);
    if(!xp.inherits("eigen-JacobiRotation")){
      isJacobiRotation = false;
    }else{
      isJacobiRotation = true;
    }
  }
  return isJacobiRotation;
}

#pragma endregion

