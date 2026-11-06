
#' @export
ArrayXXd_Constant = function(rows,cols,value){
  eigenArrayXXdConstant(rows,cols,value)
}

#' @export
ArrayXXd_Ones = function(rows,cols){
  eigenArrayXXdOnes(rows,cols)
}

#' @export
ArrayXXd_Random = function(rows,cols){
  eigenArrayXXdRandom(rows,cols)
}

#' @export
ArrayXXd_Zero = function(rows,cols){
  eigenArrayXXdZero(rows,cols)
}



#' #' @export
#' ArrayXXd_Identity = function(rows,cols){
#'   eigenArrayXXdIdentity(rows,cols)
#' }

#' @export
ArrayXd_Constant = function(size,value){
  eigenArrayXdConstant(size,value)
}

#' @export
ArrayXd_Ones = function(size){
  eigenArrayXdOnes(rsize)
}

#' @export
ArrayXd_Random = function(size){
  eigenArrayXdRandom(size)
}

#' @export
ArrayXd_Zero = function(size){
  eigenArrayXdZero(size)
}

#' @export
ArrayXd_LinSpaced = function(size,low,high){
  eigenArrayXdLinSpaced(size,low,high)
}

#' #' @export
#' ArrayXd_Unit = function(index,i){
#'   eigenArrayXdUnit(index,i)
#' }

#' @export
MatrixXcf_Constant = function(rows,cols,value){
  eigenMatrixXcfConstant(rows,cols,value)
}

#' @export
MatrixXcf_Ones = function(rows,cols){
  eigenMatrixXcfOnes(rows,cols)
}

#' @export
MatrixXcf_Random = function(rows,cols){
  eigenMatrixXcfRandom(rows,cols)
}

#' @export
MatrixXcf_Zero = function(rows,cols){
  eigenMatrixXcfZero(rows,cols)
}

#' @export
MatrixXcf_Identity = function(rows,cols){
  eigenMatrixXcfIdentity(rows,cols)
}

#' @export
MatrixXd_Constant = function(rows,cols,value){
  eigenMatrixXdConstant(rows,cols,value)
}

#' @export
MatrixXd_Ones = function(rows,cols){
  eigenMatrixXdOnes(rows,cols)
}

#' @export
MatrixXd_Random = function(rows,cols){
  eigenMatrixXdRandom(rows,cols)
}

#' @export
MatrixXd_Zero = function(rows,cols){
  eigenMatrixXdZero(rows,cols)
}

#' @export
MatrixXd_Identity = function(rows,cols){
  eigenMatrixXdIdentity(rows,cols)
}

#' @export
VectorXd_Constant = function(size,value){
  eigenVectorXdConstant(size,value)
}

#' @export
VectorXd_LinSpaced = function(size,low,high){
  eigenVectorXdLinSpaced(size,low,high)
}

#' @export
VectorXd_Ones = function(size){
  eigenVectorXdOnes(size)
}

#' @export
VectorXd_Random = function(size){
  eigenVectorXdRandom(size)
}

#' @export
VectorXd_Zero = function(size){
  eigenVectorXdZero(size)
}

#' @export
VectorXd_Unit = function(index,i){
  eigenVectorXdUnit(index,i)
}

#' @export
RowVectorXd_Constant = function(size,value){
  eigenRowVectorXdConstant(size,value)
}

#' @export
RowVectorXd_LinSpaced = function(size,low,high){
  eigenRowVectorXdLinSpaced(size,low,high)
}

#' @export
RowVectorXd_Ones = function(size){
  eigenRowVectorXdOnes(size)
}

#' @export
RowVectorXd_Random = function(size){
  eigenRowVectorXdRandom(size)
}

#' @export
RowVectorXd_Zero = function(size){
  eigenRowVectorXdZero(size)
}

#' @export
RowVectorXd_Unit = function(index,i){
  eigenRowVectorXdUnit(index,i)
}

#' @export
ArrayXXd = function(...){
  param = list(...)
  if(length(param)<=0){
    obj = new(EigenArrayXXd)
  }else{
    obj = new(EigenArrayXXd,param)
  }
  attr(obj,"cpptype")="EigenArrayXXd"
  obj
}

#' @export
ArrayXd = function(...){
  param = list(...)
  if(length(param)<=0){
    obj = new(EigenArrayXd)
  }else{
    obj = new(EigenArrayXd,param)
  }
  attr(obj,"cpptype")="EigenArrayXd"
  obj
}

#' @export
MatrixXcf = function(...){
  param = list(...)
  if(length(param)<=0){
    obj = new(EigenMatrixXcf)
  }else{
    obj = new(EigenMatrixXcf,param)
  }
  attr(obj,"cpptype")="EigenMatrixXcf"
  obj
}

#' @export
Matrix4d = function(...){
  param = list(...)
  if(length(param)<=0){
    obj = new(EigenMatrix4d)
  }else{
    obj = new(EigenMatrix4d,param)
  }
  attr(obj,"cpptype")="EigenMatrix4d"
  obj
}

#' @export
MatrixXd = function(...){
  param = list(...)
  if(length(param)<=0){
    obj = new(EigenMatrixXd)
  }else{
    obj = new(EigenMatrixXd,param)
  }
  attr(obj,"cpptype")="EigenMatrixXd"
  obj
}

#' @export
VectorXd = function(...){
  param = list(...)
  if(length(param)<=0){
    obj = new(EigenVectorXd)
  }else{
    obj = new(EigenVectorXd,param)
  }
  attr(obj,"cpptype")="EigenVectorXd"
  obj
}

#' @export
Vector3d = function(...){
  param = list(...)
  if(length(param)<=0){
    obj = new(EigenVector3d)
  }else{
    obj = new(EigenVector3d,param)
  }
  attr(obj,"cpptype")="EigenVector3d"
  obj
}


#' @export
RowVectorXd = function(...){
  param = list(...)
  if(length(param)<=0){
    obj = new(EigenRowVectorXd)
  }else{
    obj = new(EigenRowVectorXd,param)
  }
  attr(obj,"cpptype")="EigenRowVectorXd"
  obj
}
