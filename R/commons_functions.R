#
# Commonly used functions to analyse Proline outputs
#


list.flatten = function(columns) {
  all_columns = c()
  for (groupname in names(columns)) {
    all_columns = c(all_columns, columns[[groupname]])
  }
  return(all_columns)
}


#' 
#' Clean up Proline ions data by removing:
#'  - duplicates (sequence, modifications) 
#'  - columns not matching the supplied abundance_regexp or included in (peptide_id|sequence|modification|charge|peptide_ion_moz|elution|accession|psm)
#' 
#' @param ions The ions table read from a Proline export
#' 
 
cleanup_ions_dataset = function(ions, abundance_regexp = "raw_abundance") {
  # remove duplicated (sequence, modifications, charge)
  cleaned.ions = ions[!duplicated(ions[, c("sequence","modifications", "master_quant_peptide_ion_charge")]), ]
  # filter columns 
  cleaned.ions = dplyr::select(cleaned.ions, matches(paste0('peptide_id|sequence|modification|charge|peptide_ion_moz|elution|accession|psm|',abundance_regexp)))
  return (cleaned.ions)
}

#' 
#' Clean up Proline peptides data by removing:
#'  - duplicates (sequence, modifications) 
#'  - columns not matching the supplied abundance_regexp or included in (peptide_id|sequence|modification|charge|experimental_moz|elution|accession|psm)
#' 
#' @param peptides The peptides table read from a Proline export
#' 

cleanup_peptides_dataset = function(peptides, abundance_regexp = "raw_abundance") {
  # remove duplicated (sequence, modifications, charge)
  cleaned.peptides = peptides[!duplicated(peptides[, c("sequence","modifications", "charge")]), ]
  # filter columns 
  cleaned.peptides = dplyr::select(cleaned.peptides, matches(paste0('peptide_id|sequence|modification|charge|experimental_moz|elution|accession|psm|',abundance_regexp)))
  return (cleaned.peptides)
}

#' Find columns matching the specified regular expression and group them by conditions. Columns
#' are equaly partioned by conditions, using the same order (first n columns belong to the first group, 
#' the n following to the second group, etc)   
#'
#' @param df 
#' @param conditions 
#' @param regexp 
#'
#' @return
#' @export
#'
#' @examples
#' d = read.csv(text="abundance_F094661,abundance_F094674,abundance_F094684,abundance_F094660,abundance_F094673,abundance_F094683")
#' group_columns_by_conditions(d, c("condition A", "condition B"), "^abundance")

group_columns_by_conditions = function(df, conditions, regexp = "raw") {
  columns_names = grep(regexp, names(df), value=TRUE)
  splitter = as.vector(sapply(conditions, function(x) rep(x,each = length(columns_names)/length(conditions))))
  list = split(columns_names, splitter[1:length(columns_names)])
  # this is mandatory to keep the same order in list and in conditions :
  result = list()
  for (c in conditions) {
    result[[c]] = list[[c]]
  }
  return(result)
}



#' Compute log2(ratios) of columns grouped by conditions for each row. All ratios are based on 
#' the last columns group. This function returns the initial dataframe combined with
#' ratios columns 
#'
#' @param df 
#' @param columns 
#'
#' @return
#' @export
#'
#' @examples

ratios_analysis = function(df, columns) {
  all_columns = list.flatten(columns)
  n = length(columns)
  ref = apply(df[columns[[n]]], 1, function(x) { mean(x, na.rm = TRUE) })
  for (i in 1:(n-1)) {
    df[paste0('ratio_', i)] = apply(df[columns[[i]]], 1, function(x) {
      return (mean(x, na.rm = TRUE)/1.0)
    } )
  }
  
  ratios = select(df, starts_with("ratio"))
  df = select(df, -starts_with("ratio"))
  ratios = -log2(ratios/ref)
  df = (cbind(df, ratios))
  return(df)
}


#' Flag each row depending on the regular expression matching the accession column
#' A flag column is added indicating the group matching the accession. If the accession
#' does not match any group, the flag 'unknown' is used.
#'
#' @param df 
#' @param list.regexp 
#'
#' @return
#' @export
#'
#' @examples
#' df = data.frame(accession = c("WWW_COLI", "XXX_HUMAN", "YYY_YEAST", "ZZZ_YEAS" ))
#' list = list(human = "_HUMAN", yeast = "_YEAS?")
#' flag_by_accession(df, list)

flag_by_accession = function(df, list.regexp) {
  df$flag = 'unknown'
  for (label in names(list.regexp)) {
    df[which(grepl(list.regexp[[label]], df$accession)), "flag"] = label
  }
  return(df)
}


#' Analyze the dataset to extract the frequency of missing values with and without cross assignment
#' and CVs. To identifiy cross assigned values, "psm_count" columns must be present in the dataset  
#'
#' @param ions 
#' @param accession_regexp_filter 
#' @param abundances_column_regexp 
#'
#' @return
#' @export
#'
#' @examples

ions_analysis = function(ions, accession_regexp_filter = NULL, abundances_column_regexp = 'raw_abundance') {
  
  if (!is.null(accession_regexp_filter)) {
    filtered.ions = dplyr::filter(ions, grepl(accession_regexp_filter, accession))
  } else {
    filtered.ions = ions
  }
  
  abundances = dplyr::select(filtered.ions, matches(abundances_column_regexp))
  abundances = replace_zero_by_na(abundances)
  na_count = sum(is.na(abundances))
  values_count = ncol(abundances)*nrow(abundances)
  
  psms = dplyr::select(filtered.ions, matches('psm_count'))
  nocross = apply(psms, c(1,2), function(x) as.numeric(x > 0))
  missing.values = sort(apply(abundances, 1, function(x) sum(is.na(x))))
  missing.values.freq = as.data.frame(cumsum(table(missing.values)))
  missing.values.freq = cbind(c(dim(missing.values.freq)[[1]]:1), missing.values.freq)
  
  #
  # Cell by cell product : if cross assigned, abundance value is replaced by 0 and then NA
  # Since done, missing values are counted again
  #
  nocross.abundances = hadamard.prod(as.matrix(nocross), as.matrix(abundances))
  nocross.abundances[ nocross.abundances == 0 ] = NA
  missing.values = sort(apply(nocross.abundances, 1, function(x) sum(is.na(x))))
  missing.values.freq = cbind(missing.values.freq, cumsum(table(missing.values)))
  colnames(missing.values.freq) = c("runs_count", "with_cross", "without_cross")
  
  #
  # Compute CVs
  #
  cvs = as.data.frame(apply(abundances, 1, function(x) 100*( sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))))
  cvs = cbind(cvs, apply(nocross.abundances, 1, function(x) 100*( sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))))
  colnames(cvs) = c("with_cross", "without_cross")
  
  return( list("missing" = missing.values.freq, "abundances" = abundances, "nocross.abundances" = nocross.abundances, "na.count" = na_count, "CVs" = cvs))
}

#' Analyze the dataset to extract the frequency of missing values and cvs.
#'
#' @param values 
#' @param accession_regexp_filter 
#' @param abundances_column_regexp 
#'
#' @return
#' @export
#'
#' @examples
abundance_values_analysis = function(values, accession_regexp_filter = NULL, abundances_column_regexp = 'raw_abundance') {
  
  if (!is.null(accession_regexp_filter)) {
    filtered_abundances = dplyr::filter(values, grepl(accession_regexp_filter, accession))
  } else {
    filtered_abundances = values
  }
  
  abundances = dplyr::select(filtered_abundances, matches(abundances_column_regexp))
  abundances = replace_zero_by_na(abundances)
  na_count = sum(is.na(abundances))
  values_count = ncol(abundances)*nrow(abundances)
  missing.values = sort(apply(abundances, 1, function(x) sum(is.na(x))))
  missing.values.freq = as.data.frame(cumsum(table(missing.values)))
  missing.values.freq = cbind(c(dim(missing.values.freq)[[1]]:1), missing.values.freq)
  
  colnames(missing.values.freq) = c("runs_count", "count")
  
  cvs = apply(abundances, 1, function(x) 100*( sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)))
  
  return( list("missing" = missing.values.freq, "abundances" = abundances, "na.count" = na_count, "CVs" = cvs))
}

#' Build a summary of identified and quantified (fully or partially) rows in the supplied dataset
#'
#' @param df The dataset
#' @param abundances_column_regexp The regular expression used to indentify abundance columns
#'
#' @return
#' @export
#'
#' @examples

quantified_summary = function(df, abundances_column_regexp = 'raw_abundance') {
  
  abundances = dplyr::select(df, matches(abundances_column_regexp))
  result = abundances
  result$quantified = apply(abundances, 1, function(x) !all(is.na(x)))
  result$has_missing = apply(abundances, 1, function(x) any(is.na(x)))
  result$missing = apply(abundances, 1, function(x) sum(is.na(x)))
  result$partially_quantified = result$quantified & result$has_missing 
  result$fully_quantified = result$quantified & !result$has_missing 
  result$flag = df$flag
  
  l = unique(result$flag)
  
  d <- data.frame(matrix(ncol = length(l), nrow = 0))
  colnames(d) = l
  d = cbind(name = character(0), d)
  
  m = result %>% count(flag) 
  m = data.frame(m) %>% select(flag, n) %>% tidyr::spread(flag, n)
  m$name = "all"
  d = bind_rows(d, m)
  
  m = result %>% filter(quantified == FALSE) %>% count(flag)
  m = data.frame(m) %>% select(flag, n) %>% tidyr::spread(flag, n)
  if (nrow(m) > 0) {
    m$name = "not_quantified"
    d = bind_rows(d, m)
  }
  
  m = result %>% filter(quantified == TRUE) %>% count(flag) 
  m = data.frame(m) %>% select(flag, n) %>% tidyr::spread(flag, n)
  if (nrow(m) > 0) {
    m$name = "quantified"
    d = bind_rows(d, m)
  }
  
  m = result %>% filter(fully_quantified == TRUE) %>% count(flag)
  m = data.frame(m) %>% select(flag, n)  %>% tidyr::spread(flag, n)
  if (nrow(m) > 0) {
    m$name = "fully_quantified"
    d = bind_rows(d, m)
  }

  m = result %>% filter(partially_quantified == TRUE) %>% count(flag) 
  m = data.frame(m) %>% select(flag, n) %>% tidyr::spread(flag, n)
  if (nrow(m) > 0) {
    m$name = "partially_quantified"
    d = bind_rows(d, m)
  }
  
  m = result %>% filter(quantified == TRUE) %>% group_by(flag) %>% summarise(n = sum(missing))
  m = data.frame(m) %>% select(flag, n) %>% tidyr::spread(flag, n)
  t = result %>% filter(quantified == TRUE) %>% group_by(flag) %>% summarize(n = n() * ncol(abundances))
  t = data.frame(t) %>% select(flag, n) %>% tidyr::spread(flag, n)
  if (nrow(m) > 0) {
    m$name = "missing_values"
    d = bind_rows(d, m)
    t$name = "expected_values"
    d = bind_rows(d, t)
  }
  
  d$total = apply(dplyr::select(d, -name), 1, function(x) sum(x, na.rm = TRUE))
  

  rates = 100 * d[nrow(d)-1, -1] / d[nrow(d), -1]
  name = "missing_value_rate"
  rates = cbind(name = as.character(name), rates,stringsAsFactors = FALSE)
  rownames(rates) = NULL
  
  return(list("counts" = d, "rates" = rates)) 
}

flags_summary = function(df) {
  sm = df %>% count(flag)
  sm = data.frame(sm) %>% select(flag, n) %>% tidyr::spread(flag, n)
  return(sm)
}

#' Perform DAPAR differential analysis
#'
#' @param r 
#' @param condition 
#' @param abundance.column.regexp 
#' @param normalize 
#' @param imputation 
#' @param quantile.qval 
#' @param test.method 
#'
#' @return
#' @export
#'
#' @examples

dapar_differential_protein_analysis = function(df, condition, abundance.column.regexp = "abundance", normalize = FALSE, imputation = "quantile", quantile.qval = 0.025, test.method = "Welch") {
  
  dataset = df
  columns = group_columns_by_conditions(dataset, strsplit(condition, "-")[[1]], abundance.column.regexp)
  column.abundance.indexes = which(colnames(dataset) %in% list.flatten(columns))
  column.other.indexes = which(!(colnames(dataset) %in% list.flatten(columns)))
  
  msnset = DAPAR::createMSnset(file = dataset,
                               metadata = .dapar_build_metadata(columns),
                               indExpData = column.abundance.indexes, 
                               indFData = column.other.indexes, 
                               indiceID = which(colnames(dataset) == "accession"),
                               logData = TRUE)
  
  filteredMsnset = DAPAR::mvFilter(msnset, "atLeastOneCond", length(columns[[1]]))
  
  if(normalize) {
    filteredMsnset = DAPAR::wrapper.normalizeD(filteredMsnset, "QuantileCentering", "overall", quantile = 0.5)
  }
  
  if (imputation == "quantile") {
    filteredMsnset = DAPAR::wrapper.impute.detQuant(filteredMsnset, qval = quantile.qval, factor = 1.0)
  } else {
    mec.indexes = DAPAR::findMECBlock(filteredMsnset)
    filteredMsnset = DAPAR::wrapper.impute.slsa(filteredMsnset)
    filteredMsnset = DAPAR::reIntroduceMEC(filteredMsnset, mec.indexes)
    filteredMsnset = DAPAR::wrapper.impute.detQuant(filteredMsnset, qval = quantile.qval, factor = 1.0)
  }
  
  if (test.method == "Limma") {
    ttest = DAPAR::limmaCompleteTest(qData = Biobase::exprs(filteredMsnset), sTab = Biobase::pData(filteredMsnset), comp.type = "OnevsOne")
  }  else {
    ttest = DAPAR::compute.t.tests(qData = Biobase::exprs(filteredMsnset), Conditions = Biobase::pData(filteredMsnset)[,"Condition"], type = "Welch")
  }
  
  result = data.frame(as.character(rownames(filteredMsnset)), ttest$logFC, ttest$P_Value)
  colnames(result) = c("accession", "logFC", "pvalue")
  
  dataset = dataset %>% left_join(result, by = "accession") 
  
  dataset = .rename_columns(dataset,columns)
  dataset$dataset = condition
  
  return(dataset)
}


.dapar_build_metadata = function(columns) {
  all_columns = list.flatten(columns)
  groupCount = 1
  samples = c()
  conditions = c()
  for (groupname in names(columns)) {
    samples = c(samples, columns[[groupname]])
    conditions = c(conditions, rep(groupname, length(columns[[groupname]])))
    groupCount = groupCount + 1
  }
  replicates = c(1:length(all_columns))
  metadata = data.frame(samples, conditions, replicates)
  colnames(metadata) = c("Sample.name", "Condition", "Bio.Rep")
  return(metadata)
}

.rename_columns = function(df, columns) {
  all_columns = list.flatten(columns)
  count = 1
  to = c()
  for (groupname in names(columns)) {
    to = c(to, paste0("C", count, "_R", 1:length(columns[[groupname]])))
    count = count + 1
  }
  renamed = plyr::rename(df, setNames(to, all_columns))
  return(renamed)
}


#' Compute ROC curve based on the pvalue from the supplied dataset. pvalue and flag columns must have been filled before
#' calling this function.  
#'
#' @param df The dataset
#' @param flag The flag value identifying expected differentially expressed proteins  
#' @param expected The theoritically expected number of differential proteins
#'
#' @return
#' @export
#'
#' @examples

compute_roc_curve = function(df, flag = "HUMAN", expected = 144) {
  data = df[order(df$pvalue),]
  result = data.frame(fdp=double(nrow(data)-1),tpr=double(nrow(data)-1),pvalue = double(nrow(data)-1), stringsAsFactors = F)
  tp = nrow(data[data$flag == flag,])
  fp = nrow(data[data$flag != flag,])
  for (k in nrow(data):2) {
    if (data$flag[k] == flag) {
      tp = tp - 1
    } else {
      fp = fp - 1
    }
    
    fdp = (fp/(fp+tp))*100
    tpr = (tp/expected)*100 
    
    result$fdp[nrow(data) - (k-1)] = fdp
    result$tpr[nrow(data) - (k-1)] = tpr
    
    result$pvalue[nrow(data) - (k-1)] = data$pvalue[k]
  }
  
  lastValue = result$fdp[1]
  for (i in 2:(nrow(result))) {
    if (!is.na(result$fdp[i]) & (result$fdp[i] > lastValue)) { 
      result$fdp[i] = NA
    } else {
      lastValue = result$fdp[i]
    }  
  }
  result = result[!is.na(result$fdp),]
  result = result[order(result$tpr, result$fdp, decreasing = TRUE),]
  return(result)
}


#' Replace 0 by NA in columns matching the supplied regexp or in all columns if none
#' supplied 
#'
#' @param values 
#' @param columns_regexp 
#'
#' @return
#' @export
#'
#' @examples
#' values = data.frame(a = c(1,0,3), b = c(4,5,NA))
#' replace_zero_by_na(values)
#' replace_zero_by_na(values, "^a")
#' 
replace_zero_by_na = function(values, columns_regexp) {
  if (!missing(columns_regexp) && !is.null(columns_regexp)) {
    columns = colnames(values)[grep(columns_regexp, colnames(values))]
    is.na(values[, columns]) <- !values[, columns]
  } else {
    is.na(values[]) <- !values[]
  }
  return(values)
}

# ------------------------------------------------------
# Legacy functions, no more used in this code
# ------------------------------------------------------


#' Compute the expected ratio corresponding to the supplied conditions
#' Each item of the conditions list is splitted into words
#' and the first word of each item is converted to numeric value to compute 
#' the ratio. Returned ratio is always > 1.
#' 
#' @param conditions The list of conditions (list of string)
#' @return A ratio value extracted from the supplied conditions labels.
#' 

expected_ratios = function(conditions) {
  values = as.numeric(unlist(lapply(strsplit(conditions, " "), function(x) x[1])))
  if (values[1] > values[2]) {
    ratio = values[1]/values[2]
  } else {
    ratio = values[2]/values[1]
  }
  return(ratio)
}


#' Title
#'
#' @param df 
#' @param columns The list of abundance columns to be considered
#' @param id_column The name of the column defined as row id. if set, the duplicated rows are removed first
#'
#' @return
#' @export
#'

missing_values_analysis = function(df, columns, id_column) {
  if (!missing(id_column)) {
    df = df[!duplicated(df[[id_column]]),]
  }
  for (groupname in names(columns)) {
    mv_colname = paste0("mv.",groupname)
    print(paste("condition : ",groupname, " = ", paste(columns[[groupname]], collapse = ',')))
    df[,mv_colname] = apply(df[columns[[groupname]]], 1, function(x)
      sum(is.na(x)))
    print("missing values frequency : ")
    print(knitr::kable(plyr::count(df[, mv_colname])))
  }
  return (df)
}

# this function suppose that missing_values_analysis have been perfom first
# since only rows without missing values are considered
# df: the abundance matrix, including mv.* columns
# columns: the columns to be considered
# name : the name of the condition to consider 

cv_analysis = function(df, columns, name) {
  mv_colname = paste0("mv.", name)
  complete = df[df[mv_colname] == 0,columns[[name]]]
  complete$cv = apply(complete, MARGIN = 1, FUN = function(x)sd(x)/mean(x))
  complete$name = name
  print(paste0("condition ", name, "(", nrow(complete)," rows)", " cv analysis : "))
  qtl = data.frame(quantile(complete$cv, c(.10, .20, .25, 0.30, .50, .75, .90, .95, .99)))
  names(qtl) = c("percentile value")
  print(knitr::kable(qtl))
  return(complete)
}


filter_missing_values = function(df, columns) {
  
  mv_columns = apply(data.frame(names(columns)), 1, function(x)paste0("mv.",x))
  for (groupname in names(columns)) {
    mv_colname = paste0("mv.",groupname)
    df[,mv_colname] = apply(df[columns[[groupname]]], 1, function(x) sum(is.na(x)))
  }
  df$min_mv = apply(df[mv_columns], 1, min)
  filtered = df[(df$min_mv == 0), ]
  
  return(filtered)
}




