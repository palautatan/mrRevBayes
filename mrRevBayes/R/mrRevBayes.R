#' configSubsetCheck()
#'
#' Cleans up configuration files by separating columns that are oddly together and deleting duplicate columns
#' @param config_mat matrix - NOT A RELATIVE PATH
#' @keywords mrRevBayes
#' @export

configSubsetCheck = function(config_mat) {
  # CHECK COLUMN NAMES
  dataset_cols = lapply(config_mat[1,], function(x) strsplit(x, " ")[[1]])
  multiple_in_col = which(lapply(dataset_cols, length)>1)

  if (length(multiple_in_col)>1) {
    new_cols = c()
    for (each_multiple in multiple_in_col) {
      original_col = config_mat[,each_multiple]
      for (each_col_set in dataset_cols[each_multiple])
        for (each_item in each_col_set) {
          gen_col = c(each_item, original_col[2:length(original_col)])
          new_cols = cbind(new_cols, gen_col)
        }
    }

    # NULL THE ORIGINAL DUPLICATE COLUMNS
    config_mat = config_mat[,-multiple_in_col]
    config_mat = cbind2(config_mat, new_cols)
  }


  return(config_mat)
}


#' getBayesBlock()
#'
#' Returns the Bayes block from a nexus file
#' @param nexus_file Relative path to a nexus file with a Bayes block in it
#' @keywords mrRevBayes
#' @export

# READS IN BAYES BLOCK
getBayesBlock = function(nexus_file) {
  # READ IN FILE
  nexus_bayes_block = scan(nexus_file, what="character", sep="\n")
  orig_bayes_block  = nexus_bayes_block
  nexus_bayes_block = tolower(nexus_bayes_block)

  # ISOLATE BLOCK
  begin_mrbayes = which(grepl("begin mrbayes;", nexus_bayes_block)==TRUE)
  end_mrbayes   = begin_mrbayes + which(grepl("end;", nexus_bayes_block[begin_mrbayes:length(nexus_bayes_block)])==TRUE)
  bayes_block   = orig_bayes_block[begin_mrbayes:end_mrbayes]

  # REMOVE GUNK + REFORMAT
  bayes_block = gsub("\\[(.)+\\]|(\t)+", "", bayes_block) # REMOVE EXTRA COMMENTS OR TABS
  bayes_block = gsub("( )?=( )?", "=", bayes_block) # REMOVE SPACE AROUND EQUAL SIGNS
  bayes_block = gsub("^( )", "", bayes_block) # REMOVE ANY LEADING WHITE SPACE
  bayes_block = gsub(" ;", ";", bayes_block) # REMOVE SEMICOLON WHITE SPACE

  # REMOVE MULTI-LINE COMMENTS
  if (any(grepl("\\[", bayes_block))) {
    begin_mc = which(grepl("\\[", bayes_block))
    end_mc = which(grepl("\\]", bayes_block))
    bayes_block = bayes_block[!bayes_block %in% bayes_block[begin_mc:end_mc]]
  }

  return(bayes_block)

}


#' getSubsets()
#'
#' Returns the subsets specified by a Bayes block
#' @param bayes_block str - Takes the output from getBayesBlock()
#' @keywords mrRevBayes
#' @export

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
    data_subsets = "uniform"
  }

  return(data_subsets)
}

#' nst_name()
#'
#' Translates nst to a model name
#' @param nst_num int/str - Values include NA, 1, 2, or 6
#' @keywords mrRevBayes
#' @export

# MRBAYES LSET NST
nst_name = function(nst_num) {
  if (is.na(nst_num)) {
    model = "JC"
    return(model)}
  if (nst_num == 1)    {model = "F81"}
  if (nst_num == 2)    {model = "HKY"}
  if (nst_num == 6)    {model = "GTR"}
  return(model)
}

#' rate_logical()
#'
#' Translates the ASRV information
#' @param rates str - Values include "equal", "gamma", "propinv", "invgamma"
#' @keywords mrRevBayes
#' @export

# MRBAYES LSET RATES
rate_logical = function(rates) {
  if (is.na(rates)) {
    asrv_info = c(FALSE, FALSE)
    return(asrv_info)
  }
  # GAMMA FIRST
  if (rates == "equal")    {asrv_info = c(FALSE, FALSE)}
  if (rates == "gamma")    {asrv_info = c(TRUE, FALSE)}
  if (rates == "propinv")  {asrv_info = c(FALSE, TRUE)}
  if (rates == "invgamma") {asrv_info = c(TRUE, TRUE)}
  return(asrv_info)
}

#' linkMatrix()
#'
#' Returns a matrix of linkages in Bayes block -- description to be updated later
#' @param link_string TBD
#' @param bayes_links TBD
#' @param headers boolean - the first two lines will be items and config_loc
#' @param main boolean - whether this matrix is standalone; otherwise, the matrix will be appended to another one
#' @param blank boolean - blank matrix with just the headers will be returned
#' @keywords mrRevBayes
#' @export

# READ THE LINKS INTO A MATRIX
linkMatrix = function(link_string, bayes_links, headers=TRUE, main=TRUE, blank=FALSE) {
  # HEADERS: headers=TRUE denotes that the first two lines will be items and config_loc
  # MAIN: main=TRUE denotes that this matrix is standalone; otherwise, the matrix will be appended to another one
  # BLANK: blank=TRUE denotes that a blank matrix with just the headers will be returned

  items = c("revmat","statefreq", "tratio", "shape", "pinvar", "ratemultiplier", "brlenspr")
  config_mat_pos = c(3, 4, NA, 5, 6, 7, 8)

  if (blank) {
    link_mat = matrix(nrow=3, ncol=length(items))
    link_mat[1,] = items
    link_mat[2,] = config_mat_pos
    link_mat[3,] = rep("all", ncol(link_mat))
    return(link_mat)
  }

  if (main) {
    link_mat = matrix(nrow=length(bayes_links)+2, ncol=length(items))
  } else {
    link_mat = matrix(nrow=3, ncol=length(items))
  }
  link_mat[1,] = items
  link_mat[2,] = config_mat_pos

  i = 1
  for (item in items) {
    pattern = paste0("(?<=", item, "=\\()([A-z0-9](,)?)+([0-9])?")
    link_mat[3,i] = str_extract(link_string, pattern)
    i = i+1
  }

  if (headers) {return(link_mat)}
  else {return(link_mat[3:length(link_mat[,1]),])}
}


#' link_params()
#'
#' Classifies whether or not their exists linking lines
#' @param link_string str - Takes in a linking line from a Bayes block
#' @keywords mrRevBayes
#' @export

# LINK OR UNLINKED
link_params = function(link_string) {
  if (grepl("^(\t)?unlink ", link_string)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' readLinks()
#'
#' Reads in links within mrRevBayes()
#' @param bayes_block Output from getBayesBlock()
#' @keywords mrRevBayes
#' @export

# READ LINKS
readLinks = function(bayes_block) {
  links = which(grepl("^( )?(un)?link ", bayes_block)==TRUE)
  link_string = bayes_block[links]

  # 5a. IF THERE EXIST LINKAGES, CREATE LINKAGE TABLE
  if (length(links)==1) { # ONE LINK LINE
    link_mat = linkMatrix(link_string, links)
    if (grepl("^( |\t)?unlink ", link_string)==TRUE)  {linker=FALSE}  else  {linker=TRUE}

  } else { # MORE THAN ONE LINK LINE
    linker = c()
    link_mat = linkMatrix(link_string[1], links, main=FALSE)
    linker[1] = link_params(bayes_block[links[1]])
    for (i in 2:length(links)) {
      link_mat = rbind2(link_mat, linkMatrix(link_string[i], links, headers=FALSE, main=FALSE))
      linker[i] = link_params(bayes_block[links[i]])
    }
  }

  # NAME THE MATRIX ROWS TO NOT BE CONFUSED
  matrix_names = c()
  for (i in 1:length(links)) {
    matrix_names = c(matrix_names, paste0("line_", i))
  }
  rownames(link_mat) = c("parameter", "config_loc", matrix_names)

  return(list(link_mat, linker))
}

#' readLsets()
#'
#' Reads in lsets within mrRevBayes()
#' @param bayes_block Output from getBayesBlock()
#' @param num_subsets int/str - How many subsets specified in the Bayes block
#' @keywords mrRevBayes
#' @export

readLsets = function(bayes_block, num_subsets){
  lsets = which(grepl("^lset ", tolower(bayes_block))==TRUE)
  lset_string = gsub("\\= ", "\\=", bayes_block[lsets])

  # 4a. GET THE INFORMATION FROM LSETS
  # IF THERE ARE NO SUBSETS
  if (num_subsets == 1) {
    lset_applyto = "all"
    lset_nst = str_extract(string=lset_string, pattern="(?<=nst=)[0-9]")
    lset_rates = str_extract(string=lset_string, pattern="(?<=rates=)[A-z]+")
    return(rbind(lset_applyto, lset_nst, lset_rates))

  } else {
    # THERE ARE MULTIPLE SUBSETS
    lset_applyto = str_extract(string=lset_string, pattern="(?<=applyto=\\()([0-9],)*[0-9]{1,2}|(?<=applyto=\\()(all)")
    if (any(grepl(",", lset_applyto))==TRUE) {
      lset_applyto = strsplit(lset_applyto, ",")
    }
    lset_nst = str_extract(string=lset_string, pattern="(?<=nst=)[0-9]")
    lset_rates = str_extract(string=lset_string, pattern="(?<=rates=)[A-z]+")
    return(rbind(lset_applyto, lset_nst, lset_rates))
  }
}


#' readPriors()
#'
#' Reads in priors within mrRevBayes()
#' @param bayes_block Output from getBayesBlock()
#' @keywords mrRevBayes
#' @export

readPriors = function (bayes_block) {

  # CHECK THAT A PRIOR STRING EXISTS
  if (any(grepl("^prset", bayes_block)==TRUE)) {
    priors = which(grepl("^prset", bayes_block)==TRUE)
  } else {
    return(0)
  }

  prior_string = gsub("\\= ", "\\=", bayes_block[priors])

  # APPLYTO
  prset_applyto = str_extract(string=prior_string, pattern="(?<=applyto=\\()([0-9],)*[0-9]{1,2}|(?<=applyto=\\()(all)")

  # WRITTEN TO BE GENERALIZED FOR MORE ITEMS
  items = c("ratepr", "statefreqpr")

  # INITIALISE MATRIX
  pr_mat = matrix(nrow=2+length(prior_string), ncol=1+length(items))
  pr_mat[1,] = c("applyto", items)
  pr_mat[2,2:ncol(pr_mat)] = c("7", NA) # ADD CONFIG LOC
  pr_mat[3:(3+length(prset_applyto)-1),1] = prset_applyto

  # pattern = paste0("(?<=", item, "=)([a-z])+")
  ratepr = paste0("(?<=ratepr=)([a-z])+|[a-z]+\\((([0-9](\\.[0-9])?),)+([0-9](\\.[0-9])?)\\)")
  statefreqpr = "(?<=statefreqpr=)[a-z]+\\((([0-9](\\.[0-9])?),)+([0-9](\\.[0-9])?)\\)|[a-z]+\\([a-z]+\\)"
  ratepr_params = str_extract(prior_string, ratepr)
  ratepr_params = ifelse(!is.na(ratepr_params), ratepr_params, "fixed")
  statefreqpr_params = str_extract(prior_string, statefreqpr)
  statefreqpr_params = ifelse(!is.na(statefreqpr_params), statefreqpr_params, "dirichlet")


  # ADD PARAMETER VALUES
  pr_mat[3:nrow(pr_mat),2:ncol(pr_mat)] = cbind(ratepr_params, statefreqpr_params)

  # NAME
  # NAME THE MATRIX ROWS TO NOT BE CONFUSED
  matrix_names = c()
  for (i in 1:length(priors)) {
    matrix_names = c(matrix_names, paste0("line_", i))
  }
  rownames(pr_mat) = c("parameter", "config_loc", matrix_names)

  return(pr_mat)
}


#' simplifyRows()
#'
#' Reads in priors within mrRevBayes()
#' @param this_row Row from the configuration file scheme
#' @keywords mrRevBayes
#' @export

simplifyRows = function(this_row) {
  unique_vals = unique(this_row)
  no_zeros = setdiff(unique_vals,0)

  new_row = c()
  i=1
  for (value in this_row) {
    if (value==0) {
      new_row[i] = 0
      i = i+1
    } else {
      new_row[i] = which(no_zeros==value)
      i = i+1
    }
  }

  return(new_row)
}


#' mrRevBayes()
#'
#' Takes a nexus file and makes a scheme
#' @param nexus_file Relative path to nexus file with a Bayes block in it
#' @keywords mrRevBayes
#' @export

mrRevBayes = function(nexus_file) {
  # 1. getBayesBlock
  # 2. getSubsets
  # 3. START MATRIX
  # 4. readLsets
  # 5. LINKS
  # 6. readPriors
  # 7. RATE MULTIPLIERS
  # 8. SIMPLIFY ROWS
  # 9. NA REPLACEMENT


  # LIBRARIES
  require(stringr)


  # 1. READ IN NEXUS WITH BAYES BLOCK
  # bayes_block = getBayesBlock(nexus_file)
  bayes_block = getBayesBlock(nexus_file)



  # 1b. TURN ON/OFF SWITCHES
  if (any(grepl("^( |\t)?(un)?link ", bayes_block)))   {got_links = TRUE}   else  {got_links=FALSE}
  if (any(grepl("^prset", bayes_block)))               {got_priors = TRUE}  else  {got_priors=FALSE}


  # 2. FIND PARTITIONS / LACK OF PARTITIONS
  data_subsets = getSubsets(bayes_block)
  num_subsets  = length(data_subsets)



  # 3a. INITIALISE MATRIX
  config_mat = matrix(nrow=8, ncol=num_subsets+1)
  config_mat[,1] = c("dataset", "model", "relative rates", "stationary frequencies", "gamma", "invariant sites", "rate multiplier", "branch lengths")
  config_mat[1,2:(num_subsets+1)] = gsub("^( )|( )$", "", data_subsets)

  # 3b. KILL THE NA'S
  config_mat[which(is.na(config_mat), arr.ind=TRUE)] = "MB" # (FOR CONVENIENCE)



  # 4. READ LSET INFORMATION
  lset_info = readLsets(bayes_block, num_subsets)


  # 4b. LOOP OVER ALL LSETS TO TRANSLATE MRBAYES INTO REV CONFIG
  link_num_now = 1
  for (i in 1:ncol(lset_info)) {
    # A1. IF (ALL) == TRUE
    if (lset_info[1,i] == "all") {

      # A2. MODELS AND ASRV
      for (x in 1:num_subsets) {
        config_mat[2,2:(x+1)] = nst_name(lset_info[2,i]) # REMOVED UNLIST

        if (any(rate_logical(lset_info[3,i])) == TRUE) {
          true_val = which(rate_logical(lset_info[3,i])==TRUE)
          config_mat[(true_val)+4,(x+1)] = link_num_now
          link_num_now = link_num_now + 1
        }
      }

    } else { # B1. IF (ALL) == FALSE
      for (j in 1:length(unlist(lset_info[1,i]))) {
        dataset_integers = as.integer(unlist(lset_info[1,i])[j])
        config_mat[2,(dataset_integers+1)] = nst_name(lset_info[2,i])

        if (any(rate_logical(lset_info[3,i]))==TRUE) {
          true_val = which(rate_logical(lset_info[3,i])==TRUE)
          config_mat[(true_val)+4,(dataset_integers+1)] = link_num_now
        }
      }
      link_num_now = link_num_now + 1
    }
  }

  # 4c. DEFAULT MODEL
  na_models = which(config_mat[2,]=="MB")
  for (i in na_models) {
    config_mat[2,i] = "JC"
  }

  # 4d. DEFAULT ASRV
  na_g = which(config_mat[5,]=="MB")
  for (i in na_g) {
    config_mat[5,i] = FALSE
  }

  na_i = which(config_mat[6,]=="MB")
  for (i in na_i){
    config_mat[6,i] = FALSE
  }



  # --- 5. FIND PRIORS
  prior_mat = readPriors(bayes_block)


  # --- CHECK IF PRIOR INFO EXISTS
  if (class(prior_mat)=="matrix") {

    # * LOOP OVER THE INFO AVAILABLE
    for (i in 3:nrow(prior_mat)) { # ROWS


      # --- CHECK WHERE A MODEL CHANGE IS NEEDED
      # * IF STATE FREQUENCIES ARE EQUAL
      if (grepl("fixed(equal)", prior_mat[i,3], fixed=TRUE)) {

        # --- COLLECT THE MODELS
        applyto = ifelse(unlist(prior_mat[i,1]=="all"), 1:num_subsets, as.numeric(unlist(prior_mat[i,1])))
        line_models = config_mat[2,(applyto+1)]

        # --- CHANGE THE MODELS
        for (current_model in line_models) {

          # * HKY CORRESPONDS TO NST=1
          if (current_model=="F81") {
            config_mat[2, (applyto+1)] = "JC"
          }

          # * HKY CORRESPONDS TO NST=2
          if (current_model=="HKY") {
            config_mat[2, (applyto+1)] = "K80"
          }

          # * GTR CORRESPONDS TO NST=6
          if (current_model=="GTR") {
            config_mat[2, (applyto+1)] = "SYM"
          }
        }
      }
    }
  }


  # --- RATE MULTIPLIERS
  rate_mult_info = prior_mat[which(prior_mat[,2]=="variable"),1:2]

  # --- SEVERAL LINES DICTATE RATE MULTIPLIERS
  if (class(rate_mult_info)=="matrix") {
    for (ix in 1:nrow(rate_mult_info)) {
      if (rate_mult_info[ix,1] == "all") {
        print("This shit says all")
      } else{
        print("This got a numbr")
      }
    }

  # --- THERE IS ONLY ONE LINE TALKING ABOUT RATE MULTIPLIERS
  } else {
    if (rate_mult_info[1] == "all") {
      config_mat[7,2:ncol(config_mat)] = 1:num_subsets
    }
  }





  # 6. FIND LINKAGES
  # LINKED BY DEFAULT (MANUAL p.40)

  # --- 6a. CHECK IF THERE IS LINKING INFORMATION
  # * SAMPLE MATRIX:
  # *             [,1]     [,2]        [,3]     [,4]    [,5]     [,6]
  # *  parameter  "revmat" "statefreq" "tratio" "shape" "pinvar" "ratemultiplier"
  # *  config_loc "3"      "4"         NA       "5"     "6"      "7"
  # *  line_1     "all"    "all"       "all"    "all"   "all"    NA


  # --- READ IN THE LINKING INFORMATION FROM THE BAYES BLOCK
  if (got_links) { # * YES, PRIOR INFORMATION SPECIFIED

    link_info = readLinks(bayes_block)
    link_mat = link_info[[1]]
    linker = link_info[[2]]

  } else {         # * NO, PRIOR INFORMATION NOT SPECIFIED

    link_mat = linkMatrix(none, none, blank=TRUE) # none IS DUMMY
    linker = TRUE
  }


  # --- 6b. REMOVE TRATIO
  # * FOR NOW, REMOVE TRATIO
  # * REMOVES COLUMN
  link_mat = link_mat[,-3]

  # --- 6bii. REMOVE RATE MULTIPLIER
  link_mat = link_mat[,-5]



  # --- 6c. LOOP THROUGH THE LINKING LINES
  # * IN MOST CASES, THERE WILL ONLY BE ONE LINKING LINE

  partition_number_status = 1
  for (i in 3:nrow(link_mat)) { # * FOR EACH OF THE "LINKING" LINES

    # --- 6d. (i). IS THIS LINE FOR LINKING OR UNLINKING?
    # * UNLINKING LINES
    if (linker[i-2]==FALSE) {

      # --- APPLY THE LINK/UNLINK TO THE PARAMETERS
      # * FOR EACH OF THE "LINKING" PARAMETERS
      for (j in 1:ncol(link_mat)) {

        # --- NO NA SUBSETS
        if (!is.na(link_mat[i,j])) {

          these_partitions = unlist(strsplit(link_mat[i,j],","))

          # --- SOME SUBSETS
          # * INTEGERS SPECIFIED
          if (length(these_partitions)>1) {
            config_mat_ix = as.numeric(unlist(link_mat[2,j]))           # THE J'TH COLUMN'S CONFIG INDEX
            these_partitions = as.numeric(these_partitions)        # SOME

            # --- IF THE CELL SAYS 'FALSE', SKIP
            # * NOT IMPORTANT TO THE MODEL
            if (any(config_mat[config_mat_ix,(these_partitions+1)] == FALSE)) {
              falses = which(config_mat[config_mat_ix,(these_partitions+1)] == FALSE)
              # config_mat[config_mat_ix,(these_partitions[falses]+1)] = 0
              these_partitions = setdiff(these_partitions,these_partitions[falses])
            }


            config_mat[config_mat_ix,(these_partitions+1)] = seq(1,length(these_partitions))


            # --- ALL SUBSETS
            # * "ALL" IS SPECIFIED
            # * NOTHING IS SPECIFIED
          } else {

            # * THIS PARTICULAR COLUMN CORRESPONDS TO A CERTAIN ROW
            # * IN THE REVBAYES CONFIGURATION MATRIX
            config_mat_ix = as.numeric(unlist(link_mat[2,j]))           # THE J'TH COLUMN'S CONFIG INDEX
            these_partitions = seq(1,num_subsets)   # ALL PARTITIONS


            # --- IF THE CELL SAYS 'FALSE', SKIP
            # * NOT IMPORTANT TO THE MODEL
            if (any(config_mat[config_mat_ix,(these_partitions+1)] == FALSE)) {
              falses = which(config_mat[config_mat_ix,(these_partitions+1)] == FALSE)
              # config_mat[config_mat_ix,(these_partitions[falses]+1)] = 0
              these_partitions = setdiff(these_partitions,these_partitions[falses])
            }

            # * GIVE ALL OF THE SUBSETS A DIFFERENT NUMBER
            config_mat[config_mat_ix,(these_partitions+1)] = these_partitions
          }
        }
      }



      # --- 6d. (ii). IS THIS LINE FOR LINKING OR UNLINKING?
      # * LINKING LINES
      # * DEFAULT
    } else {

      # --- APPLY THE LINK/UNLINK TO THE PARAMETERS
      # * FOR EACH OF THE "LINKING" PARAMETERS
      for (j in 1:ncol(link_mat)) {

        # --- NO NA SUBSETS
        if (!is.na(link_mat[i,j])) {

          these_partitions = unlist(strsplit(link_mat[i,j],","))

          # --- SOME SUBSETS
          # * INTEGERS SPECIFIED
          if (length(these_partitions)>1) {
            config_mat_ix = as.numeric(unlist(link_mat[2,j]))           # THE J'TH COLUMN'S CONFIG INDEX
            these_partitions = as.numeric(these_partitions)        # SOME

            # --- IF THE CELL SAYS 'FALSE', SKIP
            # * NOT IMPORTANT TO THE MODEL
            if (any(config_mat[config_mat_ix,(these_partitions+1)] == FALSE)) {
              falses = which(config_mat[config_mat_ix,(these_partitions+1)] == FALSE)
              # config_mat[config_mat_ix,(these_partitions[falses]+1)] = 0
              these_partitions = setdiff(these_partitions,these_partitions[falses])
            }

            config_mat[config_mat_ix,(these_partitions+1)] = partition_number_status
            partition_number_status = partition_number_status + 1


            # --- ALL SUBSETS
            # * "ALL" IS SPECIFIED
            # * NOTHING IS SPECIFIED
          } else {
            # * THIS PARTICULAR COLUMN CORRESPONDS TO A CERTAIN ROW
            # * IN THE REVBAYES CONFIGURATION MATRIX
            config_mat_ix = as.numeric(unlist(link_mat[2,j]))            # THE J'TH COLUMN'S CONFIG INDEX
            these_partitions = 1:num_subsets                        # ALL PARTITIONS

            # --- IF THE CELL SAYS 'FALSE', SKIP
            # * NOT IMPORTANT TO THE MODEL
            if (any(config_mat[config_mat_ix,(these_partitions+1)] == FALSE)) {
              falses = which(config_mat[config_mat_ix,(these_partitions+1)] == FALSE)
              # config_mat[config_mat_ix,(these_partitions[falses]+1)] = 0
              these_partitions = setdiff(these_partitions,these_partitions[falses])
            }

            # * GIVE ALL OF THE SUBSETS A DIFFERENT NUMBER
            config_mat[config_mat_ix,(these_partitions+1)] = partition_number_status
            partition_number_status = partition_number_status + 1
          }
        }
      }
    }
  }


  # 6C. KILL THE FALSES
  false_ix = which(config_mat==FALSE, arr.ind=TRUE)
  config_mat[false_ix] = 0




  # 7. LINKAGE DEFAULTS
  # --- SET THE DEFAULT LINKS
  for (i in 3:nrow(config_mat)) {
    if (any(config_mat[i,]=="MB")) {

      # * GET INDICES OF "MB" AND OTHER
      mb_ix = which(config_mat[i,]=="MB")
      other_ix = which(config_mat[i,]!="MB")
      other_ix = setdiff(other_ix,1)

      # * LOOP OVER OTHER, IF IT EXISTS, TO REMOVE LOGICALS
      if (length(other_ix)>0) {
        remove = c()

        for (j in other_ix) {
          if (config_mat[i,j]==TRUE | config_mat[i,j]==FALSE) {
            remove = c(remove, j)
          }
        }

        mb_ix = c(mb_ix, remove)
        other_ix = other_ix[!other_ix %in% other_ix[remove]]
        new_link_num = max(other_ix)+1


      } else {
        new_link_num = 1
      }

      config_mat[i,mb_ix] = new_link_num
    }
  }



  # 8. REMOVE MORPHOLOGICAL DATA
  if(any(grepl(x=tolower(config_mat[1,2:ncol(config_mat)]), pattern="morph"))) {
    morph_ix = which(grepl(x=tolower(config_mat[1,2:ncol(config_mat)]), pattern="morph"))
    config_mat = config_mat[,-(morph_ix+1)]
  }



  # 9. SIMPLIFY ROWS
  for (i in 3:nrow(config_mat)) {
    if (any((config_mat[i,2:ncol(config_mat)]=="MB"))) {
      print("Not all rows were simplified due to an MB/NA")
    } else {
      config_mat[i,2:ncol(config_mat)] = simplifyRows(config_mat[i,2:ncol(config_mat)])
    }
  }


  # RETURN VALUE
  return(config_mat)
}
