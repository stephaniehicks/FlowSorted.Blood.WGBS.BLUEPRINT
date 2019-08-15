library(dplyr)
library(DeepBlueR)
library(foreach)
library(data.table)
library(readr)

dataPath <- "/users/shicks1/data/DNAm/blueprint_ihec"

# Test installation and connectivity by saying hello to the DeepBlue server:
deepblue_info("me")

# ok this works
deepblue_list_genomes()
deepblue_list_projects() # "BLUEPRINT Epigenome", "DEEP (IHEC)"
deepblue_list_techniques() # "WGBS", "RRBS", "BisulfiteSeq"
deepblue_list_epigenetic_marks() # "DNA Methylation"
deepblue_list_biosources() # e.g. blood, muscle, etc
# deepblue_list_experiments() 

# next we search experiments
keep_biosource <- c("CD14-positive, CD16-negative classical monocyte", 
                    "CD8-positive, alpha-beta T cell", "CD4-positive, alpha-beta T cell", 
                    "CD38-negative naive B cell", 
                    "cytotoxic CD56-dim natural killer cell", 
                    "mature neutrophil", "mature eosinophil")

blueprint_DNA_meth <- deepblue_list_experiments(genome = "GRCh38",
                                                epigenetic_mark = "DNA Methylation",
                                                technique = "BisulfiteSeq",
                                                biosource = keep_biosource,
                                                project = "BLUEPRINT Epigenome")

# Then we remove `.bed` files and only lookg at `.wig` files
blueprint_DNA_meth <- 
  blueprint_DNA_meth %>% 
  filter(!grepl(".bed", name)) %>%
  data.table()

# To get more information about one experiment, use `deepblue_info()`
deepblue_info("e93346")


## Extract meta-data

# Using the experiment IDs, extract meta data about each sample, 
# including the biosource, etc. 
custom_table = do.call("rbind", apply(blueprint_DNA_meth, 1, function(experiment){
  experiment_id = experiment[1]

  # Obtain the information about the experiment_id
  info = deepblue_info(experiment_id)

  # Print the experiment name, project, biosource, and epigenetic mark.
  with(info, { data.frame(id = `_id`, name = name, project = project,
    technique = technique, epigenetic_mark = epigenetic_mark,
    biosource = sample_info$biosource_name,
    tissue_type = sample_info$TISSUE_TYPE,
    disease_status = extra_metadata$DISEASE,
    donor_id = sample_info$DONOR_ID,
    donor_age = extra_metadata$DONOR_AGE,
    donor_sex = extra_metadata$DONOR_SEX,
    experiment_id = extra_metadata$EXPERIMENT_ID,
    sample_id = sample_id,
    sample_name = sample_info$SAMPLE_NAME,
    ample_barcode = extra_metadata$SAMPLE_BARCODE,
    sample_description = extra_metadata$SAMPLE_DESCRIPTION,
    sample_source = sample_info$source,
    file_path = extra_metadata$FILE,
    first_submission_date = extra_metadata$FIRST_SUBMISSION_DATE,
    instrument_model = extra_metadata$INSTRUMENT_MODEL)
      })
}))
saveRDS(custom_table, file = file.path(dataPath,"blueprint_blood_custom_table.RDS"))

dim(custom_table)
head(custom_table)

# we also write a file with the paths to the bigwigs to download directly
write_csv(data.frame(paste0("ftp://ftp.ebi.ac.uk/pub/databases/", custom_table$file_path)), 
          file.path(dataPath,"blueprint_blood_ftp_paths.csv"), 
          col_names = FALSE)
