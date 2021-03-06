setwd("/Users/EDIE/Downloads/S1115")

scheme_csv = "S1115.csv"
subsets_folder = "data/subsets/S1115"

# * PARTITIONED ANALYSIS
partitioned_model = function(scheme_csv, subsets_folder, name) {
  
  # * 0.
  # * BEGIN SCRIPT
  script = c()
  
  # * 1A.
  # * READ IN SCHEME CSV
  scheme = read.csv(scheme_csv, as.is=TRUE)
  num_subsets = length(scheme[1,])-1

  
  # * 2.
  # * READ IN NEXUS FILES
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
  script = c(script, data_reading_lines)
  
  var_1      = strsplit(data_reading_lines[1], " =")[[1]][1]
  taxa_line  = paste0("taxa = ", var_1, ".taxa()")
  ntax_line  = paste0("num_taxa = ", var_1, ".ntaxa()")
  branches   = "num_branches = 2 * num_taxa - 3"
  tree_stats = c(taxa_line, ntax_line, branches)
  script     = c(script, tree_stats)
  
  
  # * 3.
  # * MOVES
  move_start = "move_index = 0"
  script     = c(script, move_start)
  
  
  # * 4.
  # * PARTITION VERIFICATION
  if (ncol(scheme) > 2) {
    partition = 1
  } else {
    stop("This is not a partitioned model.")
  }
  
  
  # * 5.
  # * SUBSTITUTION MODELS
  
  # * GO ROW WISE
  models = as.character(scheme[1,2:ncol(scheme)])
  
  # * IF THE RELATIVE RATES AND STATE FREQS ARE IN THE SAME PARTITIONS
  if (all(as.numeric(scheme[2,2:ncol(scheme)]) == as.numeric(scheme[3,2:ncol(scheme)]))) {
    
    # * GET THE PARTITIONS
    mod_partitions = unique(as.numeric(scheme[2,2:ncol(scheme)]))
    mod_partitions = mod_partitions[!mod_partitions==0]
    
    # * FOR EACH PARTITION,
    for (ix in mod_partitions) {
      # * GET THE SUBSETS THAT ARE IN THE PARTITION
      this_partition = which(scheme[2,]==ix)
      
      # * MAKE SURE THEY'RE ALL THE SAME MODEL
      same_mods = length(unique(models[this_partition]))==1
      
      # * ASSIGN THE MODELS
      if (same_mods) {
        this_mod = models[this_partition[1]-1]
        
        if (this_mod == "GTR") {
          pi_1 = "pi ~ dnDirichlet(v(1,1,1,1))"
          pi_2 = "moves[++move_index] = mvBetaSimplex(pi, weight=1.0)"
          er_1 = "er ~ dnDirichlet(v(1,1,1,1,1,1))"
          er_2 = "moves[++move_index] = mvBetaSimplex(er, weight=1.0)"
          
          pi_1 = gsub("pi", paste0("pi_",ix), pi_1)
          pi_2 = gsub("pi", paste0("pi_",ix), pi_2)
          er_1 = gsub("er", paste0("er_",ix), er_1)
          er_2 = gsub("er", paste0("er_",ix), er_2)
          script = c(script, (c(pi_1, pi_2, er_1, er_2)))
        }
        
      }
      
      if (this_mod != "GTR") {
        print(paste0("Need to write code to accomodate another model."))
      }

    } 
  }

  # * DO SIMILAR STEPS FOR GAMMA PARTITIONS
  gam_partitions = unique(as.numeric(scheme[4,2:ncol(scheme)]))
  gam_partitions = gam_partitions[!gam_partitions==0]
  
  for (ix in gam_partitions){
    # this_partition = which(scheme[2,]==ix)
    
    # scheme[4,this_partition]
    gam_1 = "alpha ~ dnExponential(1)"
    gam_2 = "moves[++move_index] = mvScale(alpha, weight=1.0)"
    gam_3 = "site_rates := fnDiscretizeGamma(alpha, alpha, 4)"
    
    gam_1 = gsub("alpha", paste0("alpha_",ix), gam_1)
    gam_2 = gsub("alpha", paste0("alpha_",ix), gam_2)
    gam_3 = gsub("alpha", paste0("alpha_",ix), gam_3)
    gam_3 = gsub("site_rates", paste0("site_rates_",ix), gam_3)
    
    script = c(script, c(gam_1, gam_2, gam_3))
  }
  
  # * ALSO FOR PROPORTION OF INVARIANT SITES
  pinv_partitions = unique(as.numeric(scheme[2,2:ncol(scheme)]))
  pinv_partitions = pinv_partitions[!pinv_partitions==0]
  
  for (ix in pinv_partitions){
    # this_partition = which(scheme[2,]==ix)
    
    # scheme[4,this_partition]
    pinv_1 = "pinvar ~ dnBeta(1,1)"
    pinv_2 = "moves[++move_index] = mvSlide(pinvar)"
    
    pinv_1 = gsub("pinvar", paste0("pinvar_",ix), pinv_1)
    pinv_2 = gsub("pinvar", paste0("pinvar_",ix), pinv_2)
    
    script = c(script, c(pinv_1, pinv_2))
  }
  
  
  # * 5.5
  # * MAKE THE Q-MATRICES
  subsets = colnames(scheme)[2:ncol(scheme)]
  
  for (ix in 1:length(subsets)) {
    # subsets[ix]
    
    # * CHECK IF THEY'RE THE SAME
    if (var(scheme[2:3,ix+1]) == 0) {
      part_num = scheme[2:3,ix+1][1]
      
      if (models[1] == "GTR") {
        script = c(script, (paste0("Q_", subsets[ix], ":= fnGTR(er_", part_num, ", pi_", part_num,")")))
      }
    }
  }

  
  # * 6.
  # * BRANCH RATES
  br_rate_lvs = unique(as.numeric(scheme[6,2:ncol(scheme)]))
  
  if (var(as.numeric(scheme[6,2:ncol(scheme)]))!=0) {
    for (ix in br_rate_lvs) {
      this_partition = which(scheme[6,]==ix)
      
      sum = paste(paste0("num_sites[", ix, "] = data_", subsets[this_partition-1], ".nchar()"), collapse=" + ")
      script = c(script, sum)
    }
    
    br_1 = paste0("relative_rates ~ dnDirichlet(v(", paste(rep(1, length(this_partition)), collapse=","), "))")
    br_2 = "moves[++move_index] = mvBetaSimplex(relative_rates, weight=1.0)"
    br_3 = "subset_rates := relative_rates * sum(num_sites) / num_sites"
    
    script = c(script, c(br_1, br_2, br_3))
  }
  
 
  
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
  script = c(script, phylo_lines)
  
  
  # * 7.
  # * LIKELIHOOD FUNCTIONS
  
  for (ix in 2:ncol(scheme)) {
    hi = paste0("seq_", subsets[ix-1], " ~ dnPhyloCTMC(tree=phylogeny, Q=Q_", subsets[ix-1], ", siteRates=site_rates_", scheme[4,ix], ", pInv=pinv_", scheme[5,ix], " branchRates=subset_rates[", scheme[6,ix], "])")
    bye = paste0("seq_", subsets[ix-1], ".clamp(data_", subsets[ix-1], ")")
    print(c(hi, bye))
    script = c(script, c(hi, bye))
  }
  
  

  
  # * 8.
  # * MODEL
  model = "my_model = model(phylogeny)"
  script = c(script, model)
  
  # * 9.
  # * SCRIPT
  return(script)
  
  
  
}
