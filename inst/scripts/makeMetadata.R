## validate with `ExperimentHubData::makeExperimentHubMetadata()`
## (above pkg directory)
main.data <- data.frame(
  Title = c("FlowSorted.Blood.WGBS.BLUEPRINT", 
            "FlowSorted.Blood.WGBS.BLUEPRINT (GRanges)", 
            "FlowSorted.Blood.WGBS.BLUEPRINT (col annotation)"),
  Description = rep("Whole genome bisulfite-sequencing data on sorted blood cell populations from BLUEPRINT",3),
  BiocVersion = "3.10",
  Genome = rep("hg19", 3),
  SourceType = "BigWig",
  SourceUrl = rep("ftp://ftp.ebi.ac.uk/pub/databases/blueprint/data/homo_sapiens/GRCh38/",3),
  SourceVersion = rep("BLUEPRINT",3),
  Species="Homo sapiens",
  TaxonomyId="9606",
  Coordinate_1_based=TRUE,
  DataProvider="BLUEPRINT",
  Maintainer="Stephanie Hicks <shicks19@jhu.edu>",
  RDataClass="character",
  DispatchClass=c("H5File", "RDS", "RDS"),
  stringsAsFactors = FALSE,
  RDataPath = as.vector(t(outer(
    paste0("FlowSorted.Blood.WGBS.BLUEPRINT/v1.0.0/"),
    c("blueprint_blood.h5", "blueprint_blood_gr.RDS",
      "blueprint_blood_colData.RDS"), FUN = paste0)))
)

write.csv(file="inst/extdata/metadata.csv", 
          main.data, row.names=FALSE)

# #### validated with `ExperimentHubData::makeExperimentHubMetadata()`
# ExperimentHubData::makeExperimentHubMetadata(
#     pathToPackage = "/users/shicks1/myRpkgs/FlowSorted.Blood.WGBS.BLUEPRINT",
#    fileName = "metadata.csv")
