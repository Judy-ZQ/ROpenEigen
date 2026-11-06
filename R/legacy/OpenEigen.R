setRefClass("Ref",
            fields=list(
              name = "list"
            )
)

encloseRef = function(x){
  tmp = new("Ref",name=list(x))
  tmp
}

uncloseRef = function(x){
  if(!isS4(x)){
    stop("parameter error")
  }
  x$name[[1]]
}

anyEval = function(flags,env=caller_env()){
  needbreak=FALSE
  while(TRUE){
    tryCatch({
      flags = eval_tidy(flags,env=env)
      if(!(is.language(flags) || is.symbol(flags))){
        break
      }
    },error=function(e){
      needbreak <<- TRUE
    })
    # print(needbreak)
    if(needbreak){
      break
    }
  }
  if(is_quosure(flags)){
    get_expr(flags)
  }else{
    flags
  }
}
