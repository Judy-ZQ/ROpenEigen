
AlignedBoxXd = function(...){
  new(EigenAlignedBoxXd,list(...))
}

# TransformAffineXd = function(...){
#   new(EigenTransformAffineXd,list(...))
# }

Affine3d = function(...){
  new(EigenTransformAffine3d,list(...))
}

#' @export
Affine3d_Identity = function(){
  EigenTransformAffine3d_Identity()
}

Affine2d = function(...){
  new(EigenTransformAffine2d,list(...))
}

#' @export
Affine2d_Identity = function(){
  EigenTransformAffine2d_Identity()
}

# Affine3d = function(...){
#   new(EigenAffine3d,list(...))
# }

TranslationXd = function(...){
  new(EigenTranslationXd,list(...))
}

Translation3d = function(...){
  new(EigenTranslation3d,list(...))
}

AngleAxisd = function(...){
  new(EigenAngleAxisd,list(...))
}

Quaterniond = function(...){
  new(EigenQuaterniond,list(...))
}

Rotation2Dd = function(...){
  new(EigenRotation2Dd,list(...))
}
