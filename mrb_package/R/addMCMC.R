#' func_name()
#'
#' Insert description
#' @param param_1 Insert description
#' @keywords mrRevBayes
#' @export

addMCMC = function(model_script, analysis_name, output_folder, generations=10000, burnin=250, tuning=100) {
  
  # * 1.
  # * STARTERS
  # source_line = paste0("source(\"", model_script, "\")")
  name_line = paste0("name = \"", analysis_name, "\"")
  output_line = paste0("output = \"", output_folder, "/\"")
  
  # * 2.
  # * MONITORS
  mon_1 = paste0("monitors[1] = mnModel(filename = output + name + \"/\" + name + \"_posterior_samples.log\", printgen=10, separator=TAB)")
  mon_2 = paste0("monitors[2] = mnFile(filename = output + name + \"/\" + name + \"_tree_samples.trees\", printgen=10, separator=TAB, phylogeny)")
  mon_3 = "monitors[3] = mnScreen(printgen=100, TL)"
  mon_4 = paste0("monitors[4] = mnStochasticVariable(filename = output + name + \"/\" + name + \"_posterior_stoch.var\", separator=TAB)")
  mon_lines = c(mon_1, mon_2, mon_3, mon_4)
  
  # * 3.
  # * MCMC
  an_1 = "analysis = mcmc(my_model, monitors, moves)"
  burn_line = paste0("analysis.burnin(generations=", burnin, ", tuningInterval=", tuning,")")
  an_2 = "analysis.operatorSummary()"
  gen_line = paste0("analysis.run(generations=", generations, ")")
  mcmc_lines = c(an_1, burn_line, an_2, gen_line)
  
  # * 4.
  # * MAP TREE
  trace   = paste0("treetrace = readTreeTrace(output + name + \"/\" + name + \"_tree_samples.trees\", treetype=\"non-clock\")")
  maptree = paste0("map_tree = mapTree(treetrace, output + name + \"/\" + name + \"_map_tree.tre\")")
  
  quit = "q()"
  
  
  # * Z.
  # * SCRIPT
  model_script = scan(model_script, what="character", sep="\n", blank.lines.skip=FALSE)
  script = c(model_script, name_line, output_line, "\n\n",  mon_lines, "\n\n", mcmc_lines, "\n\n", trace, maptree,"\n\n", quit)
  
  return(script)
}
