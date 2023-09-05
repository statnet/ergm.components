.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("ergm.components", c("statnet"), FALSE)
  if(!is.null(sm)) packageStartupMessage(sm)
}



