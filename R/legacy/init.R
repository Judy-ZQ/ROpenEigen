
# `=` <- function(e1,e2) {
#   UseMethod(`=`)
# }

.onAttach <- function(libname,pkgname){

}

.onLoad <- function(libname,pkgname){
  loadModule("EigenArrayModule",TRUE)
  loadModule("EigenMatrixModule",TRUE)
  loadModule("EigenVectorModule",TRUE)
  loadModule("EigenQRModule",TRUE)
  # loadModule("EigenHouseholderQRModule",TRUE)
  # loadModule("EigenColPivHouseholderQRModule",TRUE)
  loadModule("EigenCholeskyModule",TRUE)
  loadModule("EigenSVDModule",TRUE)
  loadModule("EigenLUModule",TRUE)
  loadModule("EigenEigenvaluesModule",TRUE)
  loadModule("EigenTranspositionsModule")
  loadModule("EigenPermutationMatrixModule",TRUE)
  loadModule("GeometryModule",TRUE)
  loadModule("JacobiModule",TRUE)


  EigenJacobiRotationd$methods(
    makeGivens = function(p,q,r=0){
      if(missing(r)){
        obj = .self$makeGivens_(list(p,q))
      }else{
        r_q = enquo(r)

        vn_r = get_expr(r_q)
        if(!is.symbol(vn_r)){
          stop(paste0("Error: ",vn_r," is wrong, needs to be a name."))
        }else if(!exists(vn_r,where = get_env(r_q))){
          stop(paste0("Error: ",vn_r," needs to be defined first."))
        }
        ev_r = get_env(r_q)

        oldaddr_r = lobstr::obj_addr(r)

        tmp_r = encloseRef(r)
        obj = .self$makeGivens_(list(p,q,tmp_r))
        newaddr_r= lobstr::obj_addr(tmp_r$name[[1]])

        if(oldaddr_r!=newaddr_r){
          tmp_r1 = uncloseRef(tmp_r)
          if(!is.null(attr(r,"cpptype"))){
            attr(tmp_r1,"cpptype") = attr(r,"cpptype")
          }
          ev_r[[as_string(vn_r)]] = tmp_r1
        }
      }

    },
    makeJacobi = function(...){
      .self$makeJacobi_(list(...))
    }
  )

  BDCSVD4MatrixXd$methods(
    compute = function(...){
      .self$compute_(list(...))
    }
  )

  JacobiSVD4MatrixXd$methods(
    compute = function(...){
      .self$compute_(list(...))
    }
  )

  EigenMatrixXd$methods(
    applyOnTheLeft = function(...){
      .self$applyOnTheLeft_(list(...))
    },
    applyOnTheRight = function(...){
      .self$applyOnTheRight_(list(...))
    },
    # computeInverseAndDetWithCheck = function(inverse,determinant,invertible,precision=1e-6){
    #   inverse_q = enquo(inverse)
    #   determinant_q = enquo(determinant)
    #   invertible_q = enquo(invertible)
    #
    #   vn_inverse = get_expr(inverse_q)
    #   if(!is.symbol(vn_inverse)){
    #     stop(paste0("Error: ",vn_inverse," is wrong, needs to be a name."))
    #   }else if(!exists(vn_inverse,where = get_env(inverse_q))){
    #     stop(paste0("Error: ",vn_inverse," needs to be defined first."))
    #   }
    #   ev_inverse = get_env(inverse_q)
    #
    #   vn_determinant = get_expr(determinant_q)
    #   if(!is.symbol(vn_determinant)){
    #     stop(paste0("Error: ",vn_determinant," is wrong, needs to be a name."))
    #   }else if(!exists(vn_determinant,where = get_env(determinant_q))){
    #     stop(paste0("Error: ",vn_determinant," needs to be defined first."))
    #   }
    #   ev_determinant = get_env(determinant_q)
    #
    #   vn_invertible = get_expr(invertible_q)
    #   if(!is.symbol(vn_invertible)){
    #     stop(paste0("Error: ",vn_invertible," is wrong, needs to be a name."))
    #   }else if(!exists(vn_invertible,where = get_env(invertible_q))){
    #     stop(paste0("Error: ",vn_invertible," needs to be defined first."))
    #   }
    #   ev_invertible = get_env(invertible_q)
    #
    #   oldaddr_inverse = lobstr::obj_addr(inverse)
    #   oldaddr_determinant = lobstr::obj_addr(determinant)
    #   oldaddr_invertible = lobstr::obj_addr(invertible)
    #
    #   tmp_inverse = encloseRef(inverse)
    #   tmp_determinant = encloseRef(determinant)
    #   tmp_invertible = encloseRef(invertible)
    #   # return(0)
    #   if(missing(precision)){
    #     .self$computeInverseAndDetWithCheck_(list(tmp_inverse,tmp_determinant,tmp_invertible))
    #   }else{
    #     .self$computeInverseAndDetWithCheck_(list(tmp_inverse,tmp_determinant,tmp_invertible,precision))
    #   }
    #   newaddr_inverse = lobstr::obj_addr(tmp_inverse$name[[1]])
    #   newaddr_determinant = lobstr::obj_addr(tmp_determinant$name[[1]])
    #   newaddr_invertible = lobstr::obj_addr(tmp_invertible$name[[1]])
    #
    #   if(oldaddr_inverse!=newaddr_inverse){
    #     tmp_inverse1 = uncloseRef(tmp_inverse)
    #     if(!is.null(attr(inverse,"cpptype"))){
    #       attr(tmp_inverse1,"cpptype") = attr(inverse,"cpptype")
    #     }
    #     ev_inverse[[as_string(vn_inverse)]] = tmp_inverse1
    #   }
    #
    #   if(oldaddr_determinant!=newaddr_determinant){
    #     tmp_determinant1 = uncloseRef(tmp_determinant)
    #     if(!is.null(attr(determinant,"cpptype"))){
    #       attr(tmp_determinant1,"cpptype") = attr(determinant,"cpptype")
    #     }
    #     ev_determinant[[as_string(vn_determinant)]] = tmp_determinant1
    #   }
    #
    #   if(oldaddr_invertible!=newaddr_invertible){
    #     tmp_invertible1 = uncloseRef(tmp_invertible)
    #     if(!is.null(attr(invertible,"cpptype"))){
    #       attr(tmp_invertible1,"cpptype") = attr(invertible,"cpptype")
    #     }
    #     ev_invertible[[as_string(vn_invertible)]] = tmp_invertible1
    #   }
    #
    # },
    maxCoeff = function(rowId=-1,colId=-1){
      if(missing(rowId) || missing(colId)){
        obj = .self$maxCoeff_(list())
      }else{
        rowId_q = enquo(rowId)
        colId_q = enquo(colId)

        vn_rowId = get_expr(rowId_q)
        if(!is.symbol(vn_rowId)){
          stop(paste0("Error: ",vn_rowId," is wrong, needs to be a name."))
        }else if(!exists(vn_rowId,where = get_env(rowId_q))){
          stop(paste0("Error: ",vn_rowId," needs to be defined first."))
        }
        ev_rowId = get_env(rowId_q)

        vn_colId = get_expr(colId_q)
        if(!is.symbol(vn_colId)){
          stop(paste0("Error: ",vn_colId," is wrong, needs to be a name."))
        }else if(!exists(vn_colId,where = get_env(colId_q))){
          stop(paste0("Error: ",vn_colId," needs to be defined first."))
        }
        ev_colId = get_env(colId_q)

        oldaddr_rowId = lobstr::obj_addr(rowId)
        oldaddr_colId = lobstr::obj_addr(colId)

        tmp_rowId = encloseRef(rowId)
        tmp_colId = encloseRef(colId)
        # return(0)
        obj = .self$maxCoeff_(list(tmp_rowId,tmp_colId))
        newaddr_rowId = lobstr::obj_addr(tmp_rowId$name[[1]])
        newaddr_colId = lobstr::obj_addr(tmp_colId$name[[1]])

        if(oldaddr_rowId!=newaddr_rowId){
          tmp_rowId1 = uncloseRef(tmp_rowId)
          if(!is.null(attr(rowId,"cpptype"))){
            attr(tmp_rowId1,"cpptype") = attr(rowId,"cpptype")
          }
          ev_rowId[[as_string(vn_rowId)]] = tmp_rowId1
        }

        if(oldaddr_colId!=newaddr_colId){
          tmp_colId1 = uncloseRef(tmp_colId)
          if(!is.null(attr(colId,"cpptype"))){
            attr(tmp_colId1,"cpptype") = attr(colId,"cpptype")
          }
          ev_colId[[as_string(vn_colId)]] = tmp_colId1
        }
      }
      # obj = .self$maxCoeff01_()
      obj
    },
    minCoeff = function(rowId=-1,colId=-1){
      if(missing(rowId) || missing(colId)){
        obj = .self$minCoeff_(list())
      }else{
        rowId_q = enquo(rowId)
        colId_q = enquo(colId)

        vn_rowId = get_expr(rowId_q)
        if(!is.symbol(vn_rowId)){
          stop(paste0("Error: ",vn_rowId," is wrong, needs to be a name."))
        }else if(!exists(vn_rowId,where = get_env(rowId_q))){
          stop(paste0("Error: ",vn_rowId," needs to be defined first."))
        }
        ev_rowId = get_env(rowId_q)

        vn_colId = get_expr(colId_q)
        if(!is.symbol(vn_colId)){
          stop(paste0("Error: ",vn_colId," is wrong, needs to be a name."))
        }else if(!exists(vn_colId,where = get_env(colId_q))){
          stop(paste0("Error: ",vn_colId," needs to be defined first."))
        }
        ev_colId = get_env(colId_q)

        oldaddr_rowId = lobstr::obj_addr(rowId)
        oldaddr_colId = lobstr::obj_addr(colId)

        tmp_rowId = encloseRef(rowId)
        tmp_colId = encloseRef(colId)
        # return(0)
        obj = .self$minCoeff_(list(tmp_rowId,tmp_colId))
        newaddr_rowId = lobstr::obj_addr(tmp_rowId$name[[1]])
        newaddr_colId = lobstr::obj_addr(tmp_colId$name[[1]])

        if(oldaddr_rowId!=newaddr_rowId){
          tmp_rowId1 = uncloseRef(tmp_rowId)
          if(!is.null(attr(rowId,"cpptype"))){
            attr(tmp_rowId1,"cpptype") = attr(rowId,"cpptype")
          }
          ev_rowId[[as_string(vn_rowId)]] = tmp_rowId1
        }

        if(oldaddr_colId!=newaddr_colId){
          tmp_colId1 = uncloseRef(tmp_colId)
          if(!is.null(attr(colId,"cpptype"))){
            attr(tmp_colId1,"cpptype") = attr(colId,"cpptype")
          }
          ev_colId[[as_string(vn_colId)]] = tmp_colId1
        }
      }
      # obj = .self$maxCoeff01_()
      obj
    },
    reshaped = function(...){
      .self$reshaped_(list(...))
    }
  )

  EigenVectorXd$methods(
    applyOnTheLeft = function(...){
      .self$applyOnTheLeft_(list(...))
    },
    applyOnTheRight = function(...){
      .self$applyOnTheRight_(list(...))
    },
    maxCoeff = function(index=-1){
      if(missing(index)){
        obj = .self$maxCoeff_(list())
      }else{
        index_q = enquo(index)

        vn_index = get_expr(index_q)
        if(!is.symbol(vn_index)){
          stop(paste0("Error: ",vn_index," is wrong, needs to be a name."))
        }else if(!exists(vn_index,where = get_env(index_q))){
          stop(paste0("Error: ",vn_index," needs to be defined first."))
        }
        ev_index = get_env(index_q)

        oldaddr_index = lobstr::obj_addr(index)

        tmp_index = encloseRef(index)
        obj = .self$maxCoeff_(list(tmp_index))
        newaddr_index = lobstr::obj_addr(tmp_index$name[[1]])

        if(oldaddr_index!=newaddr_index){
          tmp_index1 = uncloseRef(tmp_index)
          if(!is.null(attr(index,"cpptype"))){
            attr(tmp_index1,"cpptype") = attr(index,"cpptype")
          }
          ev_index[[as_string(vn_index)]] = tmp_index1
        }
      }
      obj
    },
    minCoeff = function(index=-1){
      if(missing(index)){
        obj = .self$minCoeff_(list())
      }else{
        index_q = enquo(index)

        vn_index = get_expr(index_q)
        if(!is.symbol(vn_index)){
          stop(paste0("Error: ",vn_index," is wrong, needs to be a name."))
        }else if(!exists(vn_index,where = get_env(index_q))){
          stop(paste0("Error: ",vn_index," needs to be defined first."))
        }
        ev_index = get_env(index_q)

        oldaddr_index = lobstr::obj_addr(index)

        tmp_index = encloseRef(index)
        obj = .self$minCoeff_(list(tmp_index))
        newaddr_index = lobstr::obj_addr(tmp_index$name[[1]])

        if(oldaddr_index!=newaddr_index){
          tmp_index1 = uncloseRef(tmp_index)
          if(!is.null(attr(index,"cpptype"))){
            attr(tmp_index1,"cpptype") = attr(index,"cpptype")
          }
          ev_index[[as_string(vn_index)]] = tmp_index1
        }
      }
      obj
    }
  )

  EigenRowVectorXd$methods(
    maxCoeff = function(index=-1){
      if(missing(index)){
        obj = .self$maxCoeff_(list())
      }else{
        index_q = enquo(index)

        vn_index = get_expr(index_q)
        if(!is.symbol(vn_index)){
          stop(paste0("Error: ",vn_index," is wrong, needs to be a name."))
        }else if(!exists(vn_index,where = get_env(index_q))){
          stop(paste0("Error: ",vn_index," needs to be defined first."))
        }
        ev_index = get_env(index_q)

        oldaddr_index = lobstr::obj_addr(index)

        tmp_index = encloseRef(index)
        obj = .self$maxCoeff_(list(tmp_index))
        newaddr_index = lobstr::obj_addr(tmp_index$name[[1]])

        if(oldaddr_index!=newaddr_index){
          tmp_index1 = uncloseRef(tmp_index)
          if(!is.null(attr(index,"cpptype"))){
            attr(tmp_index1,"cpptype") = attr(index,"cpptype")
          }
          ev_index[[as_string(vn_index)]] = tmp_index1
        }
      }
      obj
    },
    minCoeff = function(index=-1){
      if(missing(index)){
        obj = .self$minCoeff_(list())
      }else{
        index_q = enquo(index)

        vn_index = get_expr(index_q)
        if(!is.symbol(vn_index)){
          stop(paste0("Error: ",vn_index," is wrong, needs to be a name."))
        }else if(!exists(vn_index,where = get_env(index_q))){
          stop(paste0("Error: ",vn_index," needs to be defined first."))
        }
        ev_index = get_env(index_q)

        oldaddr_index = lobstr::obj_addr(index)

        tmp_index = encloseRef(index)
        obj = .self$minCoeff_(list(tmp_index))
        newaddr_index = lobstr::obj_addr(tmp_index$name[[1]])

        if(oldaddr_index!=newaddr_index){
          tmp_index1 = uncloseRef(tmp_index)
          if(!is.null(attr(index,"cpptype"))){
            attr(tmp_index1,"cpptype") = attr(index,"cpptype")
          }
          ev_index[[as_string(vn_index)]] = tmp_index1
        }
      }
      obj
    }
  )

  setMethod(`[`,
            signature(x = "Rcpp_EigenTransformAffine2d"),
            function (x, i=0,j=0,...)
            {
              if(missing(i)){
                i=1:x$rows()-1
              }
              if(missing(j)){
                j=1:x$cols()-1
              }
              m = x$at(i,j)
              m_r = m$cv2r()
              if(!is.null(dim(m_r)) && length(dim(m_r))>1){
                m_r = drop(m_r)
              }
              m_r
            }
  )

  setMethod(`[<-`,
            signature(x = "Rcpp_EigenTransformAffine2d"),
            function (x, i=0,j=0,...,value)
            {
              if(missing(i)){
                i = 1:x$rows()-1
              }
              if(missing(j)){
                j = 1:x$cols()-1
              }
              if(!is.null(attr(value,"class"))){
                vcls = attr(value,"class")
                if(grepl('Rcpp_EigenArray',vcls) || grepl('Rcpp_EigenRowArray',vcls)){
                  value = value[]
                }
              }
              if(length(i)*length(j)<length(value)){
                stop("The length of the replaced item is less than the length of the replacement item.")
              }
              if((length(i)*length(j)) %% length(value)!=0){
                stop("The replaced item is not a multiple of the length of the replacement value.")
              }
              if(length(i)*length(j)>length(value)){
                value = rep(value,length(i)*length(j)/length(value))
              }
              x$set(i,j,value)
              x
            }
  )

  setMethod(`[`,
            signature(x = "Rcpp_EigenTransformAffine3d"),
            function (x, i=0,j=0,...)
            {
              if(missing(i)){
                i=1:x$rows()-1
              }
              if(missing(j)){
                j=1:x$cols()-1
              }
              m = x$at(i,j)
              m_r = m$cv2r()
              if(!is.null(dim(m_r)) && length(dim(m_r))>1){
                m_r = drop(m_r)
              }
              m_r
            }
  )

  setMethod(`[<-`,
            signature(x = "Rcpp_EigenTransformAffine3d"),
            function (x, i=0,j=0,...,value)
            {
              if(missing(i)){
                i = 1:x$rows()-1
              }
              if(missing(j)){
                j = 1:x$cols()-1
              }
              if(!is.null(attr(value,"class"))){
                vcls = attr(value,"class")
                if(grepl('Rcpp_EigenArray',vcls) || grepl('Rcpp_EigenRowArray',vcls)){
                  value = value[]
                }
              }
              if(length(i)*length(j)<length(value)){
                stop("The length of the replaced item is less than the length of the replacement item.")
              }
              if((length(i)*length(j)) %% length(value)!=0){
                stop("The replaced item is not a multiple of the length of the replacement value.")
              }
              if(length(i)*length(j)>length(value)){
                value = rep(value,length(i)*length(j)/length(value))
              }
              x$set(i,j,value)
              x
            }
  )

  setMethod("print",
            signature(x = "Rcpp_EigenArrayXXd"),
            function(x){
              x[]
            }
  )

  setMethod(`[`,
            signature(x = "Rcpp_EigenArrayXXd"),
            function (x, i=0,j=0,...)
            {
              if(missing(i)){
                i=1:x$rows()-1
              }
              if(missing(j)){
                j=1:x$cols()-1
              }
              m = x$at(i,j)
              m_r = m$cv2r()
              if(!is.null(dim(m_r)) && length(dim(m_r))>1){
                m_r = drop(m_r)
              }
              m_r
            }
  )

  setMethod(`[<-`,
            signature(x = "Rcpp_EigenArrayXXd"),
            function (x, i=0,j=0,...,value)
            {
              if(missing(i)){
                i = 1:x$rows()-1
              }
              if(missing(j)){
                j = 1:x$cols()-1
              }
              if(!is.null(attr(value,"class"))){
                vcls = attr(value,"class")
                if(grepl('Rcpp_EigenArray',vcls) || grepl('Rcpp_EigenRowArray',vcls)){
                  value = value[]
                }
              }
              if(length(i)*length(j)<length(value)){
                stop("The length of the replaced item is less than the length of the replacement item.")
              }
              if((length(i)*length(j)) %% length(value)!=0){
                stop("The replaced item is not a multiple of the length of the replacement value.")
              }
              if(length(i)*length(j)>length(value)){
                value = rep(value,length(i)*length(j)/length(value))
              }
              x$set(i,j,value)
              x
            }
  )

  setMethod("print",
            signature(x = "Rcpp_EigenArrayXd"),
            function(x){
              x[]
            }
  )

  setMethod(`[`,
            signature(x = "Rcpp_EigenArrayXd"),
            function (x, i=0,...)
            {
              if(missing(i)){
                i=1:x$rows()-1
              }
              j = 0
              m = x$at(i,j)
              m_r = m$cv2r()
              if(!is.null(dim(m_r)) && length(dim(m_r))>1){
                m_r = drop(m_r)
              }
              m_r
            }
  )

  setMethod(`[<-`,
            signature(x = "Rcpp_EigenArrayXd"),
            function (x, i=0,...,value)
            {
              if(missing(i)){
                i = 1:x$rows()-1
              }
              j = 0
              if(length(i)*length(j)<length(value)){
                stop("The length of the replaced item is less than the length of the replacement item.")
              }
              if((length(i)*length(j)) %% length(value)!=0){
                stop("The replaced item is not a multiple of the length of the replacement value.")
              }
              if(length(i)*length(j)>length(value)){
                value = rep(value,length(i)*length(j)/length(value))
              }
              x$set(i,j,value)
              x
            }
  )

  setMethod(`[`,
            signature(x = "Rcpp_EigenMatrixXcf"),
            function (x, i=0,j=0,...)
            {
              if(missing(i)){
                i=1:x$rows()-1
              }
              if(missing(j)){
                j=1:x$cols()-1
              }
              m = x$at(i,j)
              m_r = m$cv2r()
              if(!is.null(dim(m_r)) && length(dim(m_r))>1){
                m_r = drop(m_r)
              }
              m_r
            }
  )

  setMethod(`[<-`,
            signature(x = "Rcpp_EigenMatrixXcf"),
            function (x, i=0,j=0,...,value)
            {
              if(missing(i)){
                i = 1:x$rows()-1
              }
              if(missing(j)){
                j = 1:x$cols()-1
              }
              if(length(i)*length(j)<length(value)){
                stop("The length of the replaced item is less than the length of the replacement item.")
              }
              if((length(i)*length(j)) %% length(value)!=0){
                stop("The replaced item is not a multiple of the length of the replacement value.")
              }
              if(length(i)*length(j)>length(value)){
                value = rep(value,length(i)*length(j)/length(value))
              }
              x$set(i,j,value)
              x
            }
  )

  setMethod(`[`,
            signature(x = "Rcpp_EigenMatrixXd"),
            function (x, i=0,j=0,...)
            {
              if(missing(i)){
                i=1:x$rows()-1
              }
              if(missing(j)){
                j=1:x$cols()-1
              }
              m = x$at(i,j)
              m_r = m$cv2r()
              if(!is.null(dim(m_r)) && length(dim(m_r))>1){
                m_r = drop(m_r)
              }
              m_r
            }
  )

  setMethod(`[<-`,
            signature(x = "Rcpp_EigenMatrixXd"),
            function (x, i=0,j=0,...,value)
            {
              if(missing(i)){
                i = 1:x$rows()-1
              }
              if(missing(j)){
                j = 1:x$cols()-1
              }
              if(!is.null(attr(value,"class"))){
                vcls = attr(value,"class")
                if(grepl('Rcpp_EigenMatrix',vcls) || grepl('Rcpp_EigenVector',vcls) || grepl('Rcpp_EigenArray',vcls)){
                  value = value[]
                }
              }

              if(length(i)*length(j)<length(value)){
                stop("The length of the replaced item is less than the length of the replacement item.")
              }
              if((length(i)*length(j)) %% length(value)!=0){
                stop("The replaced item is not a multiple of the length of the replacement value.")
              }
              if(length(i)*length(j)>length(value)){
                value = rep(value,length(i)*length(j)/length(value))
              }
              x$set(i,j,value)
              x
            }
  )

  setMethod(`[`,
            signature(x = "Rcpp_EigenVectorXd"),
            function (x, i=0,...)
            {
              if(missing(i)){
                i=1:x$rows()-1
              }
              j = 0
              m = x$at(i,j)
              m_r = m$cv2r()
              if(!is.null(dim(m_r)) && length(dim(m_r))>1){
                m_r = drop(m_r)
              }
              m_r
            }
  )

  setMethod(`[<-`,
            signature(x = "Rcpp_EigenVectorXd"),
            function (x, i=0,...,value)
            {
              if(missing(i)){
                i = 1:x$rows()-1
              }
              j = 0
              if(length(i)*length(j)<length(value)){
                stop("The length of the replaced item is less than the length of the replacement item.")
              }
              if((length(i)*length(j)) %% length(value)!=0){
                stop("The replaced item is not a multiple of the length of the replacement value.")
              }
              if(length(i)*length(j)>length(value)){
                value = rep(value,length(i)*length(j)/length(value))
              }
              x$set(i,j,value)
              x
            }
  )

  setMethod("*", signature(e1="Rcpp_EigenArrayXXd", e2="Rcpp_EigenArrayXXd"),
            function(e1, e2) {
              res = eigenArrayXXdMul_matAndmat(e1,e2)
              if(res$rows()==1 && res$cols()==1){
                return(res[0,0])
              }else{
                return(res)
              }
            })
  setMethod("*", signature(e1="numeric", e2="Rcpp_EigenArrayXXd"),
            function(e1, e2) {
              return(eigenArrayXXdMul_sAndmat(e2,e1))
            })

  setMethod("*", signature(e1="Rcpp_EigenArrayXXd", e2="numeric"),
            function(e1, e2) {
              return(eigenArrayXXdMul_sAndmat(e1,e2))
            })
  setMethod("-", signature(e1="Rcpp_EigenArrayXXd"),
            function(e1) {
              return(eigenArrayXXdMul_sAndmat(e1,-1))
            })

  setMethod("/", signature(e1="Rcpp_EigenArrayXXd", e2="numeric"),
            function(e1, e2) {
              return(eigenArrayXXdDiv_sAndmat(e1,e2))
            })

  setMethod("+", signature(e1="Rcpp_EigenArrayXXd", e2="Rcpp_EigenArrayXXd"),
            function(e1, e2) {
              return(eigenArrayXXdAdd_matAndmat(e1,e2))
            })
  setMethod("+", signature(e1="Rcpp_EigenArrayXXd", e2="numeric"),
            function(e1, e2) {
              return(eigenArrayXXdAdd_sAndmat(e1,e2))
            })
  setMethod("+", signature(e1="numeric", e2="Rcpp_EigenArrayXXd"),
            function(e1, e2) {
              return(eigenArrayXXdAdd_sAndmat(e2,e1))
            })
  setMethod("-", signature(e1="Rcpp_EigenArrayXXd", e2="Rcpp_EigenArrayXXd"),
            function(e1, e2) {
              return(eigenArrayXXdAdd_matAndmat(e1,eigeArrayXXdMul_sAndmat(e2,-1)))
            })


  # setMethod("=", signature(e1 = "Rcpp_EigenVectorXd", e2="Rcpp_EigenArrayXXd"),
  #           function(e1, e2) {
  #             e1$assign_(e2)
  #           })


  setMethod("*", signature(e1="Rcpp_EigenMatrixXd", e2="Rcpp_EigenMatrixXd"),
            function(e1, e2) {
              res = eigenMatrixXdMul_matAndmat(e1,e2)
              # if(res$rows()==1 && res$cols()==1){
              #   return(res[0,0])
              # }else{
              #   return(res)
              # }
              return(res)
            })

  setMethod("*", signature(e1="numeric", e2="Rcpp_EigenMatrixXd"),
            function(e1, e2) {
              return(eigenMatrixXdMul_sAndmat(e2,e1))
            })

  setMethod("*", signature(e1="Rcpp_EigenMatrixXd", e2="numeric"),
            function(e1, e2) {
              return(eigenMatrixXdMul_sAndmat(e1,e2))
            })
  setMethod("-", signature(e1="Rcpp_EigenMatrixXd"),
            function(e1) {
              return(eigenMatrixXdMul_sAndmat(e1,-1))
            })

  setMethod("/", signature(e1="Rcpp_EigenMatrixXd", e2="numeric"),
            function(e1, e2) {
              return(eigenMatrixXdDiv_sAndmat(e1,e2))
            })

  setMethod("/", signature(e1="Rcpp_EigenVectorXd", e2="numeric"),
            function(e1, e2) {
              return(eigenMatrixXdDiv_sAndmat(e1,e2))
            })

  setMethod("/", signature(e1="Rcpp_EigenRowVectorXd", e2="numeric"),
            function(e1, e2) {
              return(eigenMatrixXdDiv_sAndmat(e1,e2))
            })

  setMethod("+", signature(e1="Rcpp_EigenMatrixXd", e2="Rcpp_EigenMatrixXd"),
            function(e1, e2) {
              return(eigenMatrixXdAdd_matAndmat(e1,e2))
            })
  setMethod("-", signature(e1="Rcpp_EigenMatrixXd", e2="Rcpp_EigenMatrixXd"),
            function(e1, e2) {
              return(eigenMatrixXdAdd_matAndmat(e1,eigenMatrixXdMul_sAndmat(e2,-1)))
            })

  setMethod("*", signature(e1="Rcpp_EigenMatrixXd", e2="Rcpp_EigenVectorXd"),
            function(e1, e2) {
              res = eigenMatrixXdMul_matAndmat(e1,e2)
              # if(res$rows()==1 && res$cols()==1){
              #   return(res[0,0])
              # }else{
              #   return(res)
              # }
              return(res)
            })

  setMethod("*", signature(e1="Rcpp_EigenRowVectorXd", e2="Rcpp_EigenMatrixXd"),
            function(e1, e2) {
              res = eigenMatrixXdMul_matAndmat(e1,e2)
              # if(res$rows()==1 && res$cols()==1){
              #   return(res[0,0])
              # }else{
              #   return(res)
              # }
              return(res)
            })

  setMethod("*", signature(e1="Rcpp_EigenVectorXd", e2="Rcpp_EigenRowVectorXd"),
            function(e1, e2) {
              res = eigenMatrixXdMul_matAndmat(e1,e2)
              # if(res$rows()==1 && res$cols()==1){
              #   return(res[0,0])
              # }else{
              #   return(res)
              # }
              return(res)
            })

  setMethod("*", signature(e1="Rcpp_EigenRowVectorXd", e2="Rcpp_EigenVectorXd"),
            function(e1, e2) {
              res = eigenMatrixXdMul_matAndmat(e1,e2)
              # if(res$rows()==1 && res$cols()==1){
              #   return(res[0,0])
              # }else{
              #   return(res)
              # }
              return(res)
            })

  setMethod("*", signature(e1="Rcpp_EigenTransformAffine2d", e2="Rcpp_EigenVectorXd"),
            function(e1, e2) {
              return(e1$rightproduct(e2))
            })

  setMethod("*", signature(e1="Rcpp_EigenTransformAffine3d", e2="Rcpp_EigenVectorXd"),
            function(e1, e2) {
              return(e1$rightproduct(e2))
            })

  # setMethod("*", signature(e1="Rcpp_EigenTransformAffine3d", e2="Rcpp_EigenTranslationXd"),
  #           function(e1, e2) {
  #             return(e1$rightproduct(e2))
  #           })

  setMethod("*", signature(e1="Rcpp_EigenTransformAffine3d", e2="Rcpp_EigenAngleAxisd"),
            function(e1, e2) {
              return(e1$rightproduct(e2))
            })

  setMethod("*", signature(e1="Rcpp_EigenRotation2Dd", e2="Rcpp_EigenVectorXd"),
            function(e1, e2) {
              return(e1$rightproduct(e2))
            })

  setMethod("*", signature(e1="numeric", e2="Rcpp_EigenVectorXd"),
            function(e1, e2) {
              return(eigenMatrixXdMul_sAndmat(e2,e1))
            })

  setMethod("*", signature(e1="Rcpp_EigenVectorXd", e2="numeric"),
            function(e1, e2) {
              return(eigenMatrixXdMul_sAndmat(e1,e2))
            })

  setMethod("-", signature(e1="Rcpp_EigenVectorXd"),
            function(e1) {
              return(eigenMatrixXdMul_sAndmat(e1,-1))
            })

  setMethod("-", signature(e1="Rcpp_EigenRowVectorXd"),
            function(e1) {
              return(eigenMatrixXdMul_sAndmat(e1,-1))
            })

  setMethod("+", signature(e1="Rcpp_EigenVectorXd", e2="Rcpp_EigenVectorXd"),
            function(e1, e2) {
              return(eigenMatrixXdAdd_matAndmat(e1,e2))
            })
  setMethod("-", signature(e1="Rcpp_EigenVectorXd", e2="Rcpp_EigenVectorXd"),
            function(e1, e2) {
              return(eigenMatrixXdAdd_matAndmat(e1,eigenMatrixXdMul_sAndmat(e2,-1)))
            })

  setMethod("+", signature(e1="Rcpp_EigenMatrixXd", e2="Rcpp_EigenVectorXd"),
            function(e1, e2) {
              return(eigenMatrixXdAdd_matAndmat(e1,e2))
            })
  setMethod("+", signature(e1="Rcpp_EigenVectorXd", e2="Rcpp_EigenMatrixXd"),
            function(e1, e2) {
              return(eigenMatrixXdAdd_matAndmat(e1,e2))
            })
  setMethod("-", signature(e1="Rcpp_EigenMatrixXd", e2="Rcpp_EigenVectorXd"),
            function(e1, e2) {
              return(eigenMatrixXdAdd_matAndmat(e1,eigenMatrixXdMul_sAndmat(e2,-1)))
            })
  setMethod("-", signature(e1="Rcpp_EigenVectorXd", e2="Rcpp_EigenMatrixXd"),
            function(e1, e2) {
              return(eigenMatrixXdAdd_matAndmat(e1,eigenMatrixXdMul_sAndmat(e2,-1)))
            })
  setMethod("*", signature(e1="Rcpp_EigenJacobiRotationd", e2="Rcpp_EigenJacobiRotationd"),
            function(e1, e2) {
              return(e1$rightproduct(e2))
            })

  # print.Rcpp_EigenArrayXd = function(x){
  #   x[]
  # }
  #
  # print.Rcpp_EigenArrayXXd = function(x){
  #   x[]
  # }

}
