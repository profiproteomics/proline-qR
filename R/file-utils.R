

#' Generates the reports associated to the specified dataset
#'
#' @param dataset_name the dataset name
#' @param folder the folder containing the data to analyse
#' @param params the parameters passed to the report  
#' @param output_format the report format
#' @param output_root_folderthe root folder used to generate the report
#'
#' @return
#' @export
#'

generate_report = function(dataset_name, folder, report_params, output_format, output_root_folder = 'output') {
  
  output_directory = paste0(output_root_folder, '/', dataset_name, '/', folder)
  dataset_directory = paste0('data', '/', dataset_name, '/', folder) 
    
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
#  output_report_file = paste0(output_directory, '/', folder, '_Report', get_outputfile_extension(output_format))

  report_params$input_directory = dataset_directory
  report_params$title = paste0(folder, ' (', dataset_name, ')')

  #
  # warning: using output_file = output_report_file instead of ouput_dir = output_directory fails ! 
  #
  rmarkdown::render(
    paste0('reports', '/', dataset_name ,'.Rmd'), 
    output_format = output_format,
    params = report_params,
    knit_root_dir = rprojroot::find_rstudio_root_file(),
    output_dir = output_directory
  )  
}


#' Returns the list of all subfolders containing one or more xlsx file 
#'
#' @param folder 
#'
#' @return the list of all subfolders containing one or more xlsx file
#'
#' @examples
#' list_dataset_folders('data/PME12')
#' 
list_dataset_folders = function(folder) {
    folders = list.dirs(path = folder, recursive = TRUE, full.names = FALSE)
    result = c()
    for (f in folders) {
      if (length(list.files(path = paste0(folder, "/", f), recursive = FALSE, pattern = "*\\.xlsx")) > 0 ) {
        result = c(result, f)        
      }
    }
    return(result)
}

absolute.file.path = function(files) {
  file.path(normalizePath(dirname(files)), files)
}

absolute.dir.path = function(files) {
  if (dir.exists(files)) {
    return(normalizePath(files)) 
  } else { 
      return("")
  }
}

#' Convert knit_params type to a data frame of string that can be reported as table in a document
#'
#' @param params 
#'
#' @return
#' @export
#'
#'@example

knit_params_as_dataframe = function(params) {
  df = NULL
  names = names(params)
  for (name in names) {
    value = params[[name]]
    if (is.list(value) && (length(names(value)) == length(value))) {
      value = paste(names(value),value,sep="=",collapse="; " )
    } else {
      value = toString(value)
    }
    df = rbind(df, c(name, value))
  }
  df = data.frame(df)
  colnames(df) = c('variable', 'value')
  df = df %>% arrange(variable)
  
  return(df)
}

# ------------------------------------------------------
# Legacy functions, no more used in this code
# ------------------------------------------------------


#'
#' Get the right file extension depending on the outputFormat
#'
#' @param output_format
#'
#' @return a string with the extension associated to the specified format
#' @export
#'
#'

get_outputfile_extension = function(output_format) {
  if (startsWith(output_format, "word")) {
    output_extension = ".docx"
  } else if (startsWith(output_format, "pdf")) {
    output_extension = ".pdf"
  } else {
    output_extension = ".html"
  }
  return(output_extension)
  
}



#' Search dataset files in the specified directory matching the supplied
#' list of conditions. Each item of the conditions list is splitted into words
#' and the first word of each item is used to compose a string. For example
#' c("10 fmol", "100 fmol") is transformed to "10_100" or "100_10" pattern and the
#' search file founction returns files of the directory matching one on this string.
#'
#' @param directory The directory in which files are searched
#' @param conditions The list of conditions (list of string)
#' @return A list of tsv or xlsx filenames matching the supplied conditions
#'
#'

find_dataset_files = function(directory, conditions) {
  regexp = paste0(unlist(lapply(strsplit(conditions, " "), function(x)
    x[1])), collapse = "-")
  #  regexp = paste0(regexp, "_.*\\.(tsv|xlsx)$")
  regexp = paste0(regexp, "\\.(tsv|xlsx)$")
  files = list.files(path = directory, pattern = regexp)
  if (length(files) == 0) {
    rev_regexp = paste0(unlist(lapply(strsplit(
      rev(conditions), " "
    ), function(x)
      x[1])), collapse = "-")
    rev_regexp = paste0(rev_regexp, "\\.(tsv|xlsx)$")
    files = list.files(path = directory, pattern = rev_regexp)
  }
  if (length(files) > 1) {
    filtered_files = subset(files, grepl("\\.tsv$", files))
    if (length(filtered_files) == 0) {
      filtered_files = subset(files, grepl("\\.xslx$", files))
    }
    files = filtered_files
  }
  
  if (length(files) != 1) {
    stop(
      paste(
        "cannot find a single file corresponding to pattern ",
        regexp,
        " or ",
        rev_regexp,
        " in ",
        directory
      )
    )
  }
  return(files)
}

#' Convert now date and time to a string simplified representation
#'
#' @return a String from the current date time
#' @export
#'

.get_timestamp = function() {
  return (as.character(Sys.time(),format='%Y%m%d_%H%M'))
}

