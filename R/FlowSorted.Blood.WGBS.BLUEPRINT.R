#' @name FlowSorted.Blood.WGBS.BLUEPRINT
#' @title Whole genome bisulfite-sequencing data on sorted 
#' blood cell populations from BLUEPRINT
#'
#' @description 
#' This package measures the DNA methylation of
#' sorted cell types from whole blood (venous and cord) 
#' from the BLUEPRINT consortia using the whole genome 
#' bisulfite-sequencing platform (WGBS). BigWig files
#' were downloaded and converted to HDF5 files. Data 
#' are available on ExperimentHub as a data package.
#' 
#' The HDF5 files were created using 
#' the files: 
#'      /inst/scripts/bp-01-create-data-object.R
#'      /inst/scripts/bp-02-download-data.sh
#'      /inst/scripts/bp-03-create-hdf5-bsseq-object.R
#' and the object downloaded from ExperimentHub.
#' 
#' @docType data
#' @format A BSseq object with 44 WGBS samples and 
#' approximately 29 million CpGs (or 29039352 CpGs
#' to be exact).
#' 
#' @import bsseq
#' @import ExperimentHub
#' @import S4Vectors 
#' @import HDF5Array 
#' @import SummarizedExperiment
#' @import bsseq
#' 
#' @export
#' 
#' @rdname FlowSorted.Blood.WGBS.BLUEPRINT
#' 
#' @examples 
#' library(ExperimentHub)
#' library(FlowSorted.Blood.WGBS.BLUEPRINT)
#' fs_wgbs <- FlowSorted.Blood.WGBS.BLUEPRINT()
#' dim(fs_wgbs)
#' 
FlowSorted.Blood.WGBS.BLUEPRINT <- function()
{
  hub <- ExperimentHub()
  version <- "v1.0.0"
  base <- file.path("FlowSorted.Blood.WGBS.BLUEPRINT", version, 
                    "blueprint_blood")
  # myfiles <- query(eh, "FlowSorted.Blood.WGBS.BLUEPRINT")

  rdatapath <- paste0(base, "_gr.RDS")
  gr_object <- readRDS(rdatapath)
  # gr_object <- query(hub, rdatapath)[[1]]
  # gr_object <- eh[["some_EH_ID"]]
  
  rdatapath <- paste0(base, "_colData.RDS")
  col_data <- readRDS(rdatapath)
  # colData <- query(hub, rdatapath)[[1]]
  # colData <- eh[["another_EH_ID"]] 
  
  rdatapath <- paste0(base, ".h5")
  # h5file <- query(hub, rdatapath)[[1]]
  # h5file <-eh[["some_other_EH_ID"]]
  h5file <- rdatapath

  hdf5_cov <- HDF5Array(filepath = h5file, name = "cov")
  hdf5_meth <- HDF5Array(filepath =  h5file, name = "meth")

  se <- SummarizedExperiment(list(M=hdf5_meth, Cov=hdf5_cov),
                             rowRanges=gr_object,
                             colData=DataFrame(row.names=col_data$sample_name))
  bs <- new2("BSseq", se, check=FALSE)
  colData(bs) <- DataFrame(col_data)
  bs
}
