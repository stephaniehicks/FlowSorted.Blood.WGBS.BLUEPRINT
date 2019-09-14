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
#' The BSseq object was created using 
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
#' 
#' @rdname FlowSorted.Blood.WGBS.BLUEPRINT
#' 
#' @examples 
#' library(ExperimentHub)
#' fs_wgbs <- FlowSorted.Blood.WGBS.BLUEPRINT()
#' dim(fs_wgbs)
#' 
FlowSorted.Blood.WGBS.BLUEPRINT <- function(preloaded = TRUE)
  {
    if(preloaded){
      ## Download a serialized BSseq object from ExperimentHub
      hub <- ExperimentHub()
      version <- "v1.0.0"
      base <- file.path("FlowSorted.Blood.WGBS.BLUEPRINT", version, "files_bsseq_hdf5_col")
      bs <- loadHDF5SummarizedExperiment(base)
    } 
    
    if(!preloaded){
      ## Download HDF5 files and RDS (GRanges and colData) objects, 
      ## and compose into a BSseq object. 
      hub <- ExperimentHub()
      version <- "v1.0.0"
      base <- file.path("FlowSorted.Blood.WGBS.BLUEPRINT", version, 
                        "blueprint_blood")

      ## GRanges object and column data
      rdatapath <- paste0(base, "_gr.rds")
      gr_object <- query(hub, rdatapath)[[1]]
    
      suppressMessages({
        rdatapath <- paste0(base, "_colData.rds")
        colData <- query(hub, rdatapath)[[1]]
      
        ## HDF5, from ExperimentHub:
        rdatapath <- paste0(base, ".h5")
        h5file <- query(hub, rdatapath)[[1]]
      })
      hdf5_cov <- HDF5Array(filepath = h5file, name = "cov")
      hdf5_meth <- HDF5Array(filepath =  h5file, name = "meth")
    
      bs <- BSseq(gr = gr_object, 
                  M = hdf5_meth, 
                  Cov = hdf5_cov, 
                  sampleNames = colData$sample_name)
      pData(bs) <- DataFrame(colData)
    }
  bs
}
