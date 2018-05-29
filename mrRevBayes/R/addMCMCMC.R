#' addMCMCMC()
#'
#' Takes a model script and adds a Metropolis-Coupled Markov Chain Monte Carlo analyses
#' @param model_script Relative path to the script
#' @param analysis_name str - Choice of name to save file as
#' @param output_folder Relative path to where the RevBayes output should go to
#' @param nchains=4
#' @param deltaheat=0.2
#' @param generations=10000
#' @param burnin=250
#' @param tuning=100
#' @keywords mrRevBayes
#' @export

addMCMCMC = function(model_script, analysis_name, output_folder, nchains=4, deltaheat=0.2, generations=10000, burnin=250, tuning=100) {

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
  an_1 = paste0("analysis = mcmcmc(my_model, monitors, moves, nchains=", nchains,", deltaHeat=",deltaheat,")")
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
  model_script = scan(model_script, what="character", sep="\n")
  script = c(model_script, name_line, output_line, mon_lines, mcmc_lines, trace, maptree, quit)

  return(script)
}
