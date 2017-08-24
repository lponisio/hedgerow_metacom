## inverse logit
inv.logit <- function(a){
  exp(a)/(exp(a) + 1)
}
