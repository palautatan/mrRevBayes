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
  model_script = scan(model_script, what="character", sep="\n")
  script = c(model_script, name_line, output_line,  mon_lines, mcmc_lines, trace, maptree, quit)

  return(script)
}


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



likelihood_function = function(gamma_pinv, subset_name=FALSE) {
  gamma_pinv = as.numeric(gamma_pinv)

  lf_1 = "seq ~ dnPhyloCTMC(tree=phylogeny, Q=Q, type=\"DNA\")"
  lf_2 = "seq.clamp(data)"

  if (subset_name != FALSE) {
    lf_1 = gsub("=Q", paste0("=Q_",subset_name), lf_1)
    lf_1 = gsub("seq", paste0("seq_",subset_name), lf_1)

    lf_2 = gsub("data", paste0("data_",subset_name), lf_2)
    lf_2 = gsub("seq", paste0("seq_",subset_name), lf_2)
  }

  if (gamma_pinv[1] != 0) {
    siterates_line = ", siteRates=site_rates"
    lf_1 = gsub(")", paste0(siterates_line,")"), lf_1)

    if (subset_name != FALSE) {
      lf_1 = gsub("site_rates", paste0("site_rates_", subset_name), lf_1)
    }
  }

  if (gamma_pinv[2] != 0) {
    pinv_line = ", pInv=pinvar"
    lf_1 = gsub(")", paste0(pinv_line,")"), lf_1)

    if (subset_name != FALSE) {
      lf_1 = gsub("pinvar", paste0("pinvar_", subset_name), lf_1)
    }
  }

  return(c(lf_1, lf_2))
}


# gamma_pinv = scheme[4:5,(column+1)]

gtr_model = function(gamma_pinv, subset_name="") {
  gamma_pinv = as.numeric(gamma_pinv)

  # * 1. Q MATRIX
  pi_1 = "pi ~ dnDirichlet(v(1,1,1,1))"
  pi_2 = "moves[++move_index] = mvBetaSimplex(pi, weight=1.0)"
  er_1 = "er ~ dnDirichlet(v(1,1,1,1,1,1))"
  er_2 = "moves[++move_index] = mvBetaSimplex(er, weight=1.0)"
  q    = "Q := fnGTR(er, pi)"

  if (nchar(subset_name)>0) {
    pi_1 = gsub("pi", paste0("pi_",subset_name), pi_1)
    pi_2 = gsub("pi", paste0("pi_",subset_name), pi_2)
    er_1 = gsub("er", paste0("er_",subset_name), er_1)
    er_2 = gsub("er", paste0("er_",subset_name), er_2)
    q    = gsub("pi", paste0("pi_",subset_name), q)
    q    = gsub("er", paste0("er_",subset_name), q)
    q    = gsub("Q",  paste0("Q_", subset_name), q)
  }

  q_mat_lines = c(pi_1, pi_2, er_1, er_2, q)


  # * 2. GAMMA
  if (gamma_pinv[1]) {
    gam_1 = "alpha ~ dnExponential(1)"
    gam_2 = "moves[++move_index] = mvScale(alpha, weight=1.0)"
    gam_3 = "site_rates := fnDiscretizeGamma(alpha, alpha, 4)"
    
    if (nchar(subset_name)>0) {
      gam_1 = gsub("alpha", paste0("alpha_",subset_name), gam_1)
      gam_2 = gsub("alpha", paste0("alpha_",subset_name), gam_2)
      gam_3 = gsub("alpha", paste0("alpha_",subset_name), gam_3)
      gam_3 = gsub("site_rates", paste0("site_rates_", subset_name), gam_3)
    }
  }

  if (gamma_pinv[1]) {
    gam_lines = c(gam_1, gam_2, gam_3)
  } else {
    gam_lines = ""
  }

  # * 3. PROPORTION OF INVARIANT SITES
  if (gamma_pinv[2]) {
    pinv_1 = "pinvar ~ dnBeta(1,1)"
    pinv_2 = "moves[++move_index] = mvSlide(pinvar)"
    
    if (nchar(subset_name)>0) {
      pinv_1 = gsub("pinvar", paste0("pinvar_",subset_name), pinv_1)
      pinv_2 = gsub("pinvar", paste0("pinvar_",subset_name), pinv_2)
    }
  }

  if (gamma_pinv[2]) {
    pinv_lines = c(pinv_1, pinv_2)
  } else {
    pinv_lines = ""
  }


  # * 4. CODE
  model_lines = c(q_mat_lines, gam_lines, pinv_lines)

  return(model_lines)
}


sym_model = function(gamma_pinv, subset_name="") {
  gamma_pinv = as.numeric(gamma_pinv)
  
  # * 1. Q MATRIX
  pi_1 = "pi <- simplex(1,1,1,1)"
  er_1 = "er ~ dnDirichlet(v(1,1,1,1,1,1))"
  er_2 = "moves[++move_index] = mvBetaSimplex(er, weight=1.0)"
  q    = "Q := fnGTR(er, pi)"
  
  if (nchar(subset_name)>0) {
    pi_1 = gsub("pi", paste0("pi_",subset_name), pi_1)
    er_1 = gsub("er", paste0("er_",subset_name), er_1)
    er_2 = gsub("er", paste0("er_",subset_name), er_2)
    q    = gsub("pi", paste0("pi_",subset_name), q)
    q    = gsub("er", paste0("er_",subset_name), q)
    q    = gsub("Q",  paste0("Q_", subset_name), q)
  }
  
  q_mat_lines = c(pi_1, er_1, er_2, q)
  
  
  # * 2. GAMMA
  if (gamma_pinv[1]) {
    gam_1 = "alpha ~ dnExponential(1)"
    gam_2 = "moves[++move_index] = mvScale(alpha, weight=1.0)"
    gam_3 = "site_rates := fnDiscretizeGamma(alpha, alpha, 4)"
    
    if (nchar(subset_name)>0) {
      gam_1 = gsub("alpha", paste0("alpha_",subset_name), gam_1)
      gam_2 = gsub("alpha", paste0("alpha_",subset_name), gam_2)
      gam_3 = gsub("alpha", paste0("alpha_",subset_name), gam_3)
      gam_3 = gsub("site_rates", paste0("site_rates_", subset_name), gam_3)
    }
  }
  
  if (gamma_pinv[1]) {
    gam_lines = c(gam_1, gam_2, gam_3)
  } else {
    gam_lines = ""
  }
  
  # * 3. PROPORTION OF INVARIANT SITES
  if (gamma_pinv[2]) {
    pinv_1 = "pinvar ~ dnBeta(1,1)"
    pinv_2 = "moves[++move_index] = mvSlide(pinvar)"
    
    if (nchar(subset_name)>0) {
      pinv_1 = gsub("pinvar", paste0("pinvar_",subset_name), pinv_1)
      pinv_2 = gsub("pinvar", paste0("pinvar_",subset_name), pinv_2)
    }
  }
  
  if (gamma_pinv[2]) {
    pinv_lines = c(pinv_1, pinv_2)
  } else {
    pinv_lines = ""
  }
  
  
  # * 4. CODE
  model_lines = c(q_mat_lines, gam_lines, pinv_lines)
  
  return(model_lines)
}


hky_model = function(gamma_pinv, subset_name="") {
    gamma_pinv = as.numeric(gamma_pinv)

    # * 1. Q MATRIX
    pi_1    = "pi ~ dnDirichlet(v(1,1,1,1))"
    pi_2    = "moves[++move_index] = mvBetaSimplex(pi, weight=1.0)"
    kappa_1 = "kappa ~ dnLognormal(0.0, 1.25)"
    kappa_2 = "moves[++move_index] = mvScale(kappa)"
    q       = "Q := fnHKY(kappa, pi)"

    if (nchar(subset_name)>0) {
      pi_1    = gsub("pi", paste0("pi_",subset_name), pi_1)
      pi_2    = gsub("pi", paste0("pi_",subset_name), pi_2)
      kappa_1 = gsub("kappa", paste0("kappa_",subset_name), kappa_1)
      kappa_2 = gsub("kappa", paste0("kappa_",subset_name), kappa_2)
      q       = gsub("pi", paste0("pi_",subset_name), q)
      q       = gsub("kappa", paste0("kappa_",subset_name), q)
      q       = gsub("Q",  paste0("Q_", subset_name), q)
    }

    q_mat_lines = c(pi_1, pi_2, kappa_1, kappa_2, q)


    # * 2. GAMMA
    if (gamma_pinv[1]) {
      gam_1 = "alpha ~ dnExponential(1)"
      gam_2 = "moves[++move_index] = mvScale(alpha, weight=1.0)"
      gam_3 = "site_rates := fnDiscretizeGamma(alpha, alpha, 4)"
      
      if (nchar(subset_name)>0) {
        gam_1 = gsub("alpha", paste0("alpha_",subset_name), gam_1)
        gam_2 = gsub("alpha", paste0("alpha_",subset_name), gam_2)
        gam_3 = gsub("alpha", paste0("alpha_",subset_name), gam_3)
        gam_3 = gsub("site_rates", paste0("site_rates_",subset_name), gam_3)
      }
    }


    if (gamma_pinv[1]) {
      gam_lines = c(gam_1, gam_2, gam_3)
    } else {
      gam_lines = ""
    }


    # * 3. PROPORTION OF INVARIANT SITES
    if (gamma_pinv[2]) {
      pinv_1 = "pinvar ~ dnBeta(1,1)"
      pinv_2 = "moves[++move_index] = mvSlide(pinvar)"
      
      if (nchar(subset_name)>0) {
        pinv_1 = gsub("pinvar", paste0("pinvar_",subset_name), pinv_1)
        pinv_2 = gsub("pinvar", paste0("pinvar_",subset_name), pinv_2)
      }
    }

    if (gamma_pinv[2]) {
      pinv_lines = c(pinv_1, pinv_2)
    } else {
      pinv_lines = ""
    }


    # * 4. CODE
    model_lines = c(q_mat_lines, gam_lines, pinv_lines)


    return(model_lines)
}


jc_model = function(gamma_pinv, subset_name="") {
  gamma_pinv = as.numeric(gamma_pinv)

    # * 1. Q MATRIX
    q = "Q <- fnJC(4)"
    
    if (nchar(subset_name)>0) {
      q = gsub("Q", paste0("Q_",subset_name), q)
    }

    # * 2. GAMMA
    if (gamma_pinv[1]) {
      gam_1 = "alpha ~ dnExponential(1)"
      gam_2 = "moves[++move_index] = mvScale(alpha, weight=1.0)"
      gam_3 = "site_rates := fnDiscretizeGamma(alpha, alpha, 4)"
      
      if (nchar(subset_name)>0) {
        gam_1 = gsub("alpha", paste0("alpha_",subset_name), gam_1)
        gam_2 = gsub("alpha", paste0("alpha_",subset_name), gam_2)
        gam_3 = gsub("alpha", paste0("alpha_",subset_name), gam_3)
        gam_3 = gsub("site_rates", paste0("site_rates_",subset_name), gam_3)
      }
    }


    if (gamma_pinv[1]) {
      gam_lines = c(gam_1, gam_2, gam_3)
    } else {
      gam_lines = ""
    }


    # * 3. PROPORTION OF INVARIANT SITES
    if (gamma_pinv[2]) {
      pinv_1 = "pinvar ~ dnBeta(1,1)"
      pinv_2 = "moves[++move_index] = mvSlide(pinvar)"
      
      if (nchar(subset_name)>0) {
        pinv_1 = gsub("pinvar", paste0("pinvar_",subset_name), pinv_1)
        pinv_2 = gsub("pinvar", paste0("pinvar_",subset_name), pinv_2)
      }
    }


    if (gamma_pinv[2]) {
      pinv_lines = c(pinv_1, pinv_2)
    } else {
      pinv_lines = ""
    }


    # * 4. CODE
    model_lines = c(q, gam_lines, pinv_lines)


    return(model_lines)
}


# * UNIFORM ANALYSIS
uniform_model = function(scheme_csv, nexus_file, name) {

  # * 1.
  # * READ IN SCHEME CSV
  scheme = read.csv(scheme_csv, as.is=TRUE)


  # * 2.
  # * READ IN NEXUS FILE
  read_in_data = "data = readDiscreteCharacterData("
  read_in_data = paste0(read_in_data,"\"",nexus_file,"\")\n\n")

  taxa_line  = paste0("taxa = data.taxa()")
  ntax_line  = paste0("num_taxa = data.ntaxa()")
  branches   = "num_branches = 2 * num_taxa - 3\n\n"
  tree_stats = c(taxa_line, ntax_line, branches)


  # * 3.
  # * MOVES
  move_start = "move_index = 0"

  # * 4.
  # * SUBSTITUTION MODEL
  substitution_model = scheme[1,2]

  if (substitution_model=="GTR") {
    substitution_model_lines = gtr_model(scheme[4:5,2])
  }
  
  if (substitution_model=="SYM") {
    substitution_model_lines = sym_model(scheme[4:5,2])
  }

  if (substitution_model=="HKY") {
    substitution_model_lines = hky_model(scheme[4:5,2])
  }

  if (substitution_model=="JC") {
    substitution_model_lines = jc_model(scheme[4:5,2])
  }

  # * 5.
  # * TOPOLOGY AND BRANCH LENGTHS

  top_1 = "\n\ntopology ~ dnUniformTopology(taxa)"
  top_2 = "moves[++move_index] = mvNNI(topology, weight=10.0)"
  top_3 = "moves[++move_index] = mvSPR(topology, weight=10.0)"


  br_lens_1 = "\n\nfor(i in 1:num_branches) {"
  br_lens_2 = "  br_lens[i] ~ dnExponential(10.0)"
  br_lens_3 = "  moves[++move_index] = mvScale(br_lens[i], weight=1.0)"
  br_lens_4 = "}"
  br_lens_5 = "TL := sum(br_lens)"


  phylogeny = "\n\nphylogeny := treeAssembly(topology, br_lens)"

  phylo_lines = c(top_1, top_2, top_3, br_lens_1, br_lens_2, br_lens_3, br_lens_4, br_lens_5, phylogeny)


  # * 6.
  # * LIKElIHOOD FUNCTIONS

  lf_lines = likelihood_function(scheme[4:5,2])


  # * 7.
  # * MODEL
  model = "\n\nmy_model = model(phylogeny)"

  # * 8.
  # * SCRIPT
  script = c(read_in_data, tree_stats, move_start, substitution_model_lines, phylo_lines, lf_lines, model)
  return(script)
}



# * PARTITIONED ANALYSIS
partitioned_model = function(scheme_csv, subsets_folder, name) {

    # * 1.
    # * READ IN SCHEME CSV
    scheme = read.csv(scheme_csv, as.is=TRUE)
    num_partitions = length(scheme[1,])-1


    # * 2.
    # * READ IN NEXUS FILE
    subsets = colnames(scheme)[2:ncol(scheme)]
    subsets = gsub("^X", "", subsets)
    full_paths = paste0(subsets_folder, "/", subsets)


    data_reading_lines = c()
    for (ix in 1:length(subsets)) {
      read_in_data = "data = readDiscreteCharacterData("
      read_in_data = gsub("data", paste0("data_", subsets[ix]), read_in_data)
      read_in_data = paste0(read_in_data,"\"", full_paths[ix], ".nex\")")
      data_reading_lines = c(data_reading_lines, read_in_data)
    }

    var_1      = strsplit(data_reading_lines[1], " =")[[1]][1]
    taxa_line  = paste0("taxa = ", var_1, ".taxa()")
    ntax_line  = paste0("num_taxa = ", var_1, ".ntaxa()")
    branches   = "num_branches = 2 * num_taxa - 3"
    tree_stats = c(taxa_line, ntax_line, branches)


    # * 3.
    # * MOVES
    move_start = "move_index = 0"


    # * 4.
    # * PARTITION VERIFICATION
    if (nrow(scheme) > 2) {
      partition = 1
    } else {
      stop("This is not a partitioned model.")
    }


    # * 5.
    # * SUBSTITUTION MODELS
    substitution_model_lines = c()

    for (column in 1:(ncol(scheme)-1)) {
      subset_name = subsets[column]
      substitution_model = scheme[1, (column+1)]

      if (substitution_model=="GTR") {
        substitution_model_lines = c(substitution_model_lines, gtr_model(scheme[4:5,(column+1)], subset_name))
      }
      
      if (substitution_model=="SYM") {
        substitution_model_lines = c(substitution_model_lines, sym_model(scheme[4:5,(column+1)], subset_name))
      }

      if (substitution_model=="HKY") {
        substitution_model_lines = c(substitution_model_lines, hky_model(scheme[4:5,(column+1)], subset_name))
      }

      if (substitution_model=="JC") {
        substitution_model_lines = c(substitution_model_lines, jc_model(scheme[4:5,(column+1)], subset_name))
      }
    }


    # * 6.
    # * TOPOLOGY AND BRANCH LENGTHS

    top_1 = "topology ~ dnUniformTopology(taxa)"
    top_2 = "moves[++move_index] = mvNNI(topology, weight=10.0)"
    top_3 = "moves[++move_index] = mvSPR(topology, weight=10.0)"


    br_lens_1 = "for(i in 1:num_branches){"
    br_lens_2 = "  br_lens[i] ~ dnExponential(10.0)"
    br_lens_3 = "  moves[++move_index] = mvScale(br_lens[i], weight=1.0)"
    br_lens_4 = "}"
    br_lens_5 = "TL := sum(br_lens)"


    phylogeny = "phylogeny := treeAssembly(topology, br_lens)"

    phylo_lines = c(top_1, top_2, top_3, br_lens_1, br_lens_2, br_lens_3, br_lens_4, br_lens_5, phylogeny)


    # * 7.
    # * LIKELIHOOD FUNCTIONS

    lf_lines = c()

    for (column in 1:(ncol(scheme)-1)) {
      subset_name = subsets[column]

      lf_lines = c(lf_lines, likelihood_function(scheme[4:5,(column+1)], subset_name))

    }

    # * 8.
    # * MODEL
    model = "my_model = model(phylogeny)"

    # * 9.
    # * SCRIPT
    script = c(data_reading_lines, tree_stats, move_start, substitution_model_lines, phylo_lines, lf_lines, model)
    return(script)



}
