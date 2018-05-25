
# *                             #
#     READS IN BAYES BLOCK      #
#                             * #

getBayesBlock = function(nexus_file) {
  # * READ IN FILE
  nexus_bayes_block = scan(nexus_file, what="character", sep="\n")
  nexus_bayes_block = tolower(nexus_bayes_block)
  
  # * ISOLATE BLOCK
  begin_mrbayes = which(grepl("begin mrbayes;", nexus_bayes_block)==TRUE)
  end_mrbayes   = begin_mrbayes + which(grepl("end;", nexus_bayes_block[begin_mrbayes:length(nexus_bayes_block)])==TRUE)
  bayes_block   = nexus_bayes_block[begin_mrbayes:end_mrbayes]
  
  # * REMOVE GUNK + REFORMAT
  bayes_block = gsub("\\[(.)+\\]|(\t)+", "", bayes_block) # * REMOVE EXTRA COMMENTS OR TABS
  bayes_block = gsub("( )?=( )?", "=", bayes_block) # * REMOVE SPACE AROUND EQUAL SIGNS
  bayes_block = gsub("^( )", "", bayes_block) # * REMOVE ANY LEADING WHITE SPACE
  bayes_block = gsub(" ;", ";", bayes_block) # * REMOVE SEMICOLON WHITE SPACE
  
  # * REMOVE MULTI-LINE COMMENTS
  if (any(grepl("\\[", bayes_block))) {
    begin_mc = which(grepl("\\[", bayes_block))
    end_mc = which(grepl("\\]", bayes_block))
    bayes_block = bayes_block[!bayes_block %in% bayes_block[begin_mc:end_mc]]
  }
  
  return(bayes_block)
  
}



# *                      #
#     READS SUBSETS      #
#                      * #

getSubsets = function(bayes_block) {
  if (any(grepl("^( )?partition ", bayes_block))==TRUE) {
    begin_partition = which(grepl("^( )?partition", bayes_block))[1]
    end_partition = which(grepl("^set partition", bayes_block)) - 1
    partition_line = paste(bayes_block[begin_partition:end_partition], collapse="")
    
    # * DOUBLE CHECK THE PARTITION LINE END
    partition_line = strsplit(partition_line, ";")[[1]][1]
    
    data_subsets = gsub(";", "", strsplit(partition_line, "\\:( |\t)?( )?")[[1]][2])
    data_subsets = strsplit(data_subsets,"\\,")[[1]][which(nchar(strsplit(data_subsets,"\\,")[[1]])>1)]
  } else {
    data_subsets = "unpartitioned"
  }
  
  return(data_subsets)
}




# *                  #
#     READS NST      #
#                  * #

nst_name = function(nst_num) {
  if (is.na(nst_num)) {
    model = "JC"
    return(model)}
  if (nst_num == 1)    {model = "F81"}
  if (nst_num == 2)    {model = "HKY"}
  if (nst_num == 6)    {model = "GTR"}
  return(model)
}




# *                         #
#     READS LSET RATES      #
#                         * #

rate_logical = function(rates) {
  if (is.na(rates)) {
    asrv_info = c(FALSE, FALSE)
    return(asrv_info)
  }
  # * GAMMA FIRST
  if (rates == "equal")    {asrv_info = c(FALSE, FALSE)}
  if (rates == "gamma")    {asrv_info = c(TRUE, FALSE)}
  if (rates == "propinv")  {asrv_info = c(FALSE, TRUE)}
  if (rates == "invgamma") {asrv_info = c(TRUE, TRUE)}
  return(asrv_info)
}




# *                                 #
#     FIXES NST FOR FIXED FREQ      #
#                                 * #

fixed_freq_model = function(current_model) {
  if (current_model=="F81") {
    new_model = "JC"
  }
  
  # * HKY CORRESPONDS TO NST=2
  if (current_model=="HKY") {
    new_model = "K80"
  }
  
  # * GTR CORRESPONDS TO NST=6
  if (current_model=="GTR") {
    new_model = "SYM"
  }
  
  return(new_model)
}



# *                            #
#     APPLYTO # CONVERTER      #
#                            * #

# * FUNCTION TO CONVERT STRING OF NUMBERS
num_fix = function(string_of_nums) {
  if(grepl("-", string_of_nums)) {
    set = as.numeric(unlist(strsplit(x=string_of_nums, split="-")))
    return(set[1]:set[2])
  }
  return(unlist(lapply(strsplit(x=unlist(string_of_nums), split="\\,"), as.numeric)))
}




# *                                     #
#     MRBAYES/REVBAYES TRANSLATION      #
#                                     * #
dict = list(revmat="relative rates", statefreq="stationary frequencies", 
            ratepr="rate multiplier", brlenspr="branch lengths",
            pinvar="invariant sites", shape="gamma",
            revmatpr="relative rates", statefreqpr="stationary frequencies",
            brlens="branch lengths")


mrRevBayes = function(nexus_file) {

  # *                      #
  #     READ IN BLOCK      #
  #                      * #
  
  # bayes_block = getBayesBlock("data/S10036/S10036.nex")
  bayes_block = getBayesBlock(nexus_file)
  bayes_block
  
  
  
  # *                    #
  #     GET SUBSETS      #
  #                    * #
  
  data_subsets = getSubsets(bayes_block)
  data_subsets = gsub(" ", "", data_subsets)
  data_subsets
  
  num_subsets  = length(data_subsets)
  
  
  
  
  # *                      #
  #     OUTPUT MATRIX      #
  #                      * #
  
  config_mat = matrix(nrow=8, ncol=num_subsets+1)
  config_mat[,1] = c("dataset", "model", "relative rates", "stationary frequencies", "gamma", "invariant sites", "rate multiplier", "branch lengths")
  config_mat[1,2:(num_subsets+1)] = gsub("^( )|( )$", "", data_subsets)
  config_mat[which(is.na(config_mat), arr.ind=TRUE)] = "MB" # (FOR CONVENIENCE)
  config_mat
  
  
  
  # *                     #
  #     COLLECT INFO      #
  #                     * #
  
  # * CREATE THE INFORMATION MATRIX
  require(stringr)
  
  ix = grepl("lset", bayes_block, ignore.case=TRUE) & grepl("nst", bayes_block, ignore.case=TRUE)
  ix   = ix + grepl("prset|unlink|link", bayes_block, ignore.case=TRUE)
  info = bayes_block[which(ix==1)]
  info
  
  
  
  # *                   #
  #     READ LINES      #
  #                   * #
  
  # * SPLIT ON WHITE SPACE
  info = strsplit(info, "( )+")
  info
  
  
  
  
  # *                   #
  #     LINE TYPE       #
  #                   * #
  
  # * WHAT IS THIS LINE FOR?
  line_type = c()
  for (x in 1:length(info)) {line_type = c(line_type, info[[x]][1])}
  line_type
  
  
  
  
  # *                     #
  #     DATA SUBSETS      #
  #                     * #
  
  # * WHICH PARTITIONS DOES IT AFFECT?
  applyto = str_extract(string=info, pattern="(?<=applyto=\\()([0-9]{1,2},)*[0-9]{1,2}|(?<=applyto=\\()(all)")
  applyto
  
  # * CONVERT ALL VALUES IN INFO_MAT
  # * APPLYTO VALUES IN THE INFO_MAT TABLE
  table_applyto = lapply(applyto, function(x) {
    if (is.na(x))   {return(0)}
    if (x != "all") {unlist(num_fix(x))} 
    else {1:num_subsets}
  })
  
  
  
  
  # *                   #
  #     PARAMETERS      #
  #                   * #
  
  parameters = str_extract_all(string=info, pattern="([A-z]+( )?=( )?[A-z0-9]+)|([A-z]+( )?=( )?\\(all\\))|([A-z]+( )?=( )?(.+))")
  parameters
  
  # *                          #
  #     BIND ALL TOGETHER      #
  #                          * #
  
  info_mat = cbind(line_type, table_applyto, parameters)
  info_mat
  
  
  
  # *                      #
  #     ROW ITERATION      #
  #                      * #
  
  # for (row_ix in 1:4) {
  for (row_ix in 1:nrow(info_mat)) {
    # for (row_ix in 1:11) {
    # * APPLY FUNCTIONS BASED ON LINE TYPE
    this_line_type = info_mat[row_ix,1]
    
    
    
    # *                 #
    #     LSET ONLY     #
    #                 * #
    
    if (this_line_type == "lset") {
      
      # * MODEL NAME
      nst_num    = str_extract(string=info_mat[row_ix,], pattern="(?<=nst=)[0-9]")
      nst_num    = nst_num[!is.na(nst_num)]
      this_model = nst_name(nst_num)
  
      # * CHECK STATIONARY FREQUENCIES
      applyto_check = sapply(1:length(table_applyto), function(x) {
        identical(as.numeric(unlist(info_mat[row_ix,2])), table_applyto[[x]])
      })
      
      if (any(applyto_check)) {
        matching_row_ix = which(lapply(applyto_check, function(x) TRUE %in% x)==TRUE)
        matching_row_ix = setdiff(matching_row_ix, row_ix)
        matching_row    = unlist(info_mat[matching_row_ix,])
        
        # * IF STATIONARY FREQUENCIES ARE FIXED
        # * THEN THE MODEL CHANGES
        if (any(grepl("statefreqpr", matching_row))) {
          if (any(grepl("fixed", matching_row))) {
            this_model = fixed_freq_model(this_model)
          }
        }
      }
      
      # * AMONG SITE RATE VARIATION
      # * IF THIS LINE IS ASSIGNING USING "NST",
      # * IT IS ALSO ASSIGNING AMONG SITE RATE VARIATION
      asrv      = str_extract(string=info_mat[row_ix,], pattern="(?<=rates=)[A-z]+")
      asrv      = asrv[!is.na(asrv)]
      gamma_inv = rate_logical(asrv)
      
      
      # * FILL IN THE CONFIG FILE
      these_cols = as.numeric(unlist(info_mat[row_ix,2])) + 1
      config_mat[2, these_cols] = this_model                     # * FILL IN THE MODEL
      if (gamma_inv[1]) {                                        # * FILL IN THE GAMMA ASRV
        element_num = length(unique(config_mat[5,2:ncol(config_mat)]))
        config_mat[5, these_cols] = element_num
      } else {
        config_mat[5, these_cols] = 0
      }
      if (gamma_inv[2]) {                                        # * FILL IN THE INVARIANT SITES
        element_num = length(unique(config_mat[6,2:ncol(config_mat)]))
        config_mat[6, these_cols] = element_num
      } else {
        config_mat[6, these_cols] = 0
      }
      
      
      
      
      
      # *                            #
      #     PRSET, LINK, UNLINK      #
      #                            * #
    } else {
      
      # * READ PARAMETERS
      parameters_only = str_extract_all(string=info_mat[row_ix,], pattern="[A-z]+(?=\\=)")
      parameters_only = unlist(Filter(Negate(function(x) length(x)<1), parameters_only))
      # parameters_only = unlist(Filter(Negate(function(x) x=="applyto"), parameters_only))
      parameters_only = unlist(Filter(Negate(function(x) !any(grepl(x, names(dict)))), parameters_only))
      parameters_only
      
      # these_rows = which(config_mat[,1]==as.character(unlist(dict[parameters_only]))) # * TRANSLATE MRBAYES TO REVBAYES
      these_rows = lapply(config_mat[,1], function(x) which(x==as.character(unlist(dict[parameters_only]))))
      these_rows = unlist(Filter(Negate(function(x) length(x)<1), these_rows))
      these_rows
      
      # * FILL IN THE CONFIG FILE FOR PRSET
      if (this_line_type == "prset") {
        these_cols = as.numeric(unlist(info_mat[row_ix,2])) + 1
        
        # * ASSIGN A NUMBER
        if (length(these_cols) != num_subsets) {
          # element_num = sapply(these_rows, function(x) length(unique(config_mat[x,2:ncol(config_mat)])))
          element_num = sapply(these_rows+2, function(x) length(unique(config_mat[x,2:ncol(config_mat)])))
        } else {
          element_num = 1
        }
        
        config_mat[these_rows+2,these_cols] = element_num
      }
      
      
      # * FILL IN THE CONFIG FILE FOR LINK/UNLINK
      if (grepl("link|unlink", this_line_type)) {
        # parameters_only = str_extract_all(string=info_mat[row_ix,], pattern="[A-z]+(?=\\=)")
        # parameters_only = unlist(Filter(Negate(function(x) length(x)<1), parameters_only))
        
        # these_rows = which(config_mat[,1] %in% as.character(unlist(dict[parameters_only]))) # * TRANSLATE MRBAYES TO REVBAYES
        # these_rows = lapply(config_mat[,1], function(x) which(x==as.character(unlist(dict[parameters_only]))))
        # these_rows = unlist(Filter(Negate(function(x) length(x)<1), these_rows))
        
        
        applytos_only   = str_extract_all(string=gsub("\\(|\\)", "", info_mat[row_ix,]), pattern="(?<=[a-z]{5,11}=)[A-z]+|(?<=[a-z]{5,11}=)([0-9]{1,2},)?[0-9]{1,2}+")
        applytos_only   = unlist(Filter(Negate(function(x) length(x)<1), applytos_only))
        
        applyto_link_table = rbind(parameters_only, applytos_only)
        applyto_link_table = applyto_link_table[,applyto_link_table[1,] %in% names(dict)]
        
        clean_applyto = lapply(applyto_link_table[2,], function(x) {
          if (is.na(x))   {return(0)}
          if (x != "all") {unlist(num_fix(x))} 
          else {1:num_subsets}
        })
        
        applyto_link_table = rbind(applyto_link_table, clean_applyto)
        
        # * LINKING LINES
        if (this_line_type == "link") {
          
          for (col_ix in 1:nrow(applyto_link_table)) {
            these_cols = as.numeric(unlist(applyto_link_table[3,col_ix])) + 1
            
            # * ASSIGN A NUMBER
            if (length(these_cols) != num_subsets) {
              element_num = sapply(these_rows+2, function(x) length(unique(config_mat[x,2:ncol(config_mat)]))) 
            } else {
              element_num = 1
            }
            
            config_mat[these_rows+2,these_cols] = element_num
          }
          
        
        # * UNLINKING LINES  
        } else {
          for (col_ix in 1:nrow(applyto_link_table)) {
            these_cols = as.numeric(unlist(applyto_link_table[3,col_ix])) + 1
            
            # * ASSIGN A NUMBER
            if (length(these_cols) != num_subsets) {
              starts = sapply(these_rows+2, function(x) length(unique(config_mat[x,2:ncol(config_mat)])))
              ends   = starts+length(these_cols)-1
              
              for (ix in 1:length(these_rows)) {
                config_mat[these_rows+2,these_cols] = starts[ix]:ends[ix]
              }
            
              
            print(paste("Check row", as.character(these_rows), "for accuracy."))
              
            } else {
              element_num = 1:num_subsets
              # config_mat[these_rows,these_cols] = t(replicate(nrow(config_mat[these_rows,these_cols]), element_num))
              config_mat[these_rows+2,these_cols] = t(replicate(nrow(config_mat[these_rows,these_cols]), element_num))
            }
          }
        }
      }
    }
    # * REMOVE THIS
    # cat ("Press [enter] to continue")
    # line = readline()
  }
  
  
  
  
  # *                          #
  #     LINKED BY DEFAULT      #
  #                          * #
  
  for (row_ix in 2:nrow(config_mat)) {
    if (any(config_mat[row_ix,]=="MB")) {
      these_cols = which(config_mat[row_ix,]=="MB")
      
      # * ASSIGN A NUMBER
      if (length(these_cols) != num_subsets) {
        element_num = length(unique(config_mat[row_ix,2:ncol(config_mat)]))
      } else {
        element_num = 1
      }
      
      config_mat[row_ix,these_cols] = element_num
  
    }
  }
  
  
  
  return(config_mat)
}





