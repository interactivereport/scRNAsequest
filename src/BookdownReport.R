#!/usr/bin/env Rscript
#BookdownReport.R
#This script is part of the scReport program in the scRNASequest pipeline
#By running BookdownReport.R, a bookdown report will be generated in the working directory

args = commandArgs(trailingOnly=T)

if(length(args)<1){
  message("'scAnalyzer' can be used to create a config file for an single-cell RNA-seq project")
  stop("The config.yml file is required!")
}
message("Loading resources ...")
config <- sapply(yaml::read_yaml(args[1]),unlist)

#Initiate a bookdown document
message("Rendering book ...")
suppressMessages(invisible(capture.output(
  bookdown::render_book(args[2]), 
  file = NULL)))