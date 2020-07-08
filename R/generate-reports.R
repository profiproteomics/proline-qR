#' 
#' Reports generation 
#' 

library(knitr)
library(rmarkdown)
library(dplyr)
library(ggplot2)
library(openxlsx)

source('R/file-utils.R')
source('R/commons_functions.R')


#' Generate the report for the specified dataset_name by processing only the data found 
#' in the specified folder. The folder path could be relative to the 'data' folder of this
#' project. In this case the report will be generated in the 'output' folder. The folder path
#' could also be an absolute path, in this case the report is generate into this folder.  
#' 
#' 
#' Allowed dataset_name are "UPS_4Concentrations", "UPS_10Concentrations", "ABRF" and "PME12"
#'
#' @param dataset_name the dataset name
#' @param folder the folder containing the data
#' @param output_folder the output format of the report. default value is rmdformats::readthedown.
#' Alternative formats are 'html_document', 'word_document', 'rmdformats::readthedown', 'pdf_document'
#' 
#' @return 
#' @export
#'
#' @examples
#' generate_dataset_report("PME12", "20191216/spec_bestion_sum")

generate_dataset_report = function(dataset_name, folder, output_format = "rmdformats::readthedown") {
  if (absolute.dir.path(folder) == gsub( "/", "\\\\", folder))  {
    # the folder is an absolute path
    report_params = build_report_parameters(dataset_name, folder)
    # use only the last path component as folder to build parameters
    
    output_directory = folder
    
    if (!dir.exists(output_directory)) {
      dir.create(output_directory, recursive = TRUE)
    }
    
    #  output_report_file = paste0(output_directory, '/', folder, '_Report', get_outputfile_extension(output_format))
    
    report_params$input_directory = folder
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
    
  } else {
    # the folder is a relative path 
    report_params = build_report_parameters(dataset_name, folder)
    generate_report(dataset_name, folder, report_params, output_format)
  }
}

#' Generate all reports of the specified dataset. The data processed are the subfolders found in
#'  the folder 'data/<dataset_name>'.  
#'
#' @param dataset_name the name of the dataset to process.
#' @param output_folder the output format of the report. default value is rmdformats::readthedown.
#' Alternative formats are 'html_document', 'word_document', 'rmdformats::readthedown', 'pdf_document'
#' 
#' @return 
#' @export
#'
#' @examples
#' generate_dataset_reports("PME12")

generate_dataset_reports = function(dataset_name, output_format = "rmdformats::readthedown") {
  
  folders = list_dataset_folders(paste0('data/', dataset_name))
  
  for (folder in folders) {
    report_params = build_report_parameters(dataset_name, folder)
    generate_report(dataset_name, folder, report_params, output_format)
  }
}


#' build the appropriate set of report parameters depending on the dataset and the folder name.
#' Allowed dataset_name are "UPS_4Concentrations", "UPS_10Concentrations", "ABRF" and "PME12".
#' If the folder name does contain the sub string "norm_" then no additional normalization step  
#' is applied. Conversely if the folder name does not contain "norm_" then a DAPAR quantile centering
#' normalization will be applied.  
#'
#' @param dataset_name the dataset name
#' @param folder the folder containing the xlsx reports to process.
#'
#' @return a set of parameters 
#' @export
#'
#' @examples
#' build_report_parameters("PME12", "20191216/spec_bestion_sum")

build_report_parameters = function(dataset_name, folder) {
  switch (dataset_name,
    PME12 = { list(
      abundance_column_regexp = "^abundance",
      flags_accession_regexp = list(human = '_HUMAN', yeast = '_YEAS|_yeas|_SCV?|_SACPS', ecoli = '_ECO'),
      dapar_normalize = !grepl("norm_", tolower(folder)),
      dapar_imputation = "standard",
      dapar_quantile_qval = 0.025,
      dapar_test_method = "Limma"
    ) },
    ABRF = { list(
      abundance_column_regexp = "^abundance",
      flags_accession_regexp = list(ups = 'HUMAN_UPS', background = 'YEAS|yeas'),
      dapar_normalize = !grepl("norm_", tolower(folder)),
      dapar_imputation = "quantile",
      dapar_quantile_qval = 0.025,
      dapar_test_method = "Limma"
    ) },
    UPS_10Concentrations = { list(
      abundance_column_regexp = "^abundance",
      flags_accession_regexp = list(ups = 'HUMAN_UPS', background = 'YEAS|yeas')
    ) },
    UPS_4Concentrations = { list(
      proteins_output_file = paste0('output', '/', folder, '_Proteins.xlsx'),
      abundance_column_regexp = "^abundance",
      flags_accession_regexp = list(ups = 'HUMAN_UPS', background = 'YEAS|yeas'),
      dapar_normalize = !grepl("norm_", tolower(folder)),
      dapar_imputation = "quantile",
      dapar_quantile_qval = 0.025,
      dapar_test_method = "Limma"
    ) }
  )
}