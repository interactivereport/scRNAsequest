
args = commandArgs(trailingOnly=T)
if(length(args)<2){
    #message("'EAinit' can be used to create a config file for an RNAseq project")
    stop("config yaml file is required!")
}
message("loading resource ...")
#suppressMessages(source(paste0(args[1],"utility.R"),chdir=T))
config <- sapply(yaml::read_yaml(args[2]),unlist)
sysConfig <- yaml::read_yaml(paste0(args[1],"sys.yml"))

## submit to CellDepot ------------
if (!("Is_Private" %in% names(config))) {config$Is_Private=0}
message("preparing information for CellDepot")
CellDepotCMD <- paste0("curl -s -k -X POST -d 'data={",
                      paste(paste0('"',names(config),'": "',config,'"'),collapse = ", "),
                      "}' '",sysConfig$CellDepotApp,"api_add_project.php?api_key=",sysConfig$CellDepotAppKey,"'")

message("submitting to CellDepot manager ...")
res <- system(CellDepotCMD,intern=T)
shinyMsg <- tryCatch({
  rjson::fromJSON(res)
},error=function(eMsg){
  stop(paste0(paste(res,collapse="\n"),
              "\nPlease contact ", sysConfig$powerby))
})
if(!shinyMsg$Status){
  stop(paste0(paste(paste(names(shinyMsg),shinyMsg,sep=":"),collapse="\n"),
              "\nPlease contact ",config$powerby))
}

message(paste0("CellDepot access: ", sysConfig$CellDepotApp,"app_project_review.php?ID=",shinyMsg$ID))

##save log
log_file=paste0(args[2], ".scPub.log")
cat(date(), "\n",
    "CellDepot access: ", sysConfig$CellDepotApp,"app_project_review.php?ID=",shinyMsg$ID, "\n\n",
    sep="",  file=log_file )
cat("\n\nPowered by:",sysConfig$powerby, capture.output(sessionInfo()), sep="\n", file=log_file, append=T)

