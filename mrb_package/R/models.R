#' func_name()
#'
#' Insert description
#' @param param_1 Insert description
#' @keywords mrRevBayes
#' @export

jc_model = function(ix) {
  return(0)
}


#' func_name()
#'
#' Insert description
#' @param param_1 Insert description
#' @keywords mrRevBayes
#' @export

k80_model = function(ix) {
  kappa_1 = "kappa ~ dnExp(1)"
  kappa_2 = "moves[++move_index] = mvScale(kappa, weight=1.0)"
  
  kappa_1 = gsub("kappa", paste0("kappa_",ix), kappa_1)
  kappa_2 = gsub("kappa", paste0("kappa_",ix), kappa_2)
  
  parameters = c(kappa_1, kappa_2)
  return(parameters)
}

#' func_name()
#'
#' Insert description
#' @param param_1 Insert description
#' @keywords mrRevBayes
#' @export

hky_model = function(ix) {
  pi_1    = "pi ~ dnDirichlet(v(1,1,1,1))"
  pi_2    = "moves[++move_index] = mvBetaSimplex(pi, weight=1.0)"
  kappa_1 = "kappa ~ dnLognormal(0.0, 1.25)"
  kappa_2 = "moves[++move_index] = mvScale(kappa)"
  
  pi_1 = gsub("pi", paste0("pi_",ix), pi_1)
  pi_2 = gsub("pi", paste0("pi_",ix), pi_2)
  kappa_1 = gsub("kappa", paste0("kappa_",ix), kappa_1)
  kappa_2 = gsub("kappa", paste0("kappa_",ix), kappa_2)
  
  parameters = c(pi_1, pi_2, kappa_1, kappa_2)
  return(parameters)
}


#' func_name()
#'
#' Insert description
#' @param param_1 Insert description
#' @keywords mrRevBayes
#' @export

sym_model = function(ix) {
  pi_1 = "pi <- simplex(1,1,1,1)"
  er_1 = "er ~ dnDirichlet(v(1,1,1,1,1,1))"
  er_2 = "moves[++move_index] = mvBetaSimplex(er, weight=1.0)"
  
  pi_1 = gsub("pi", paste0("pi", ix), pi_1)
  er_1 = gsub("er", paste0("er", ix), er_1)
  er_2 = gsub("er", paste0("er", ix), er_2)
  
  parameters = c(pi_1, er_1, er_2)
  return(parameters)
}


#' func_name()
#'
#' Insert description
#' @param param_1 Insert description
#' @keywords mrRevBayes
#' @export

gtr_model = function(ix) {
  pi_1 = "pi ~ dnDirichlet(v(1,1,1,1))"
  pi_2 = "moves[++move_index] = mvBetaSimplex(pi, weight=1.0)"
  er_1 = "er ~ dnDirichlet(v(1,1,1,1,1,1))"
  er_2 = "moves[++move_index] = mvBetaSimplex(er, weight=1.0)"
  
  pi_1 = gsub("pi", paste0("pi_",ix), pi_1)
  pi_2 = gsub("pi", paste0("pi_",ix), pi_2)
  er_1 = gsub("er", paste0("er_",ix), er_1)
  er_2 = gsub("er", paste0("er_",ix), er_2)
  
  parameters = c(pi_1, pi_2, er_1, er_2)
  return(parameters)
}