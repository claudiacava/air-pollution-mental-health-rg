## Step 1. Candidate phenotypes, download and quality control of Neale Lab GWAS summary statistics
##
## Output: sumstats_raw/ (Neale Lab .tsv.bgz files), sumstats_qc/ (<code>_qc.rds),
##         data/manifest.tsv, data/candidate_mental_health_phenotypes.csv

source("config.R")
source("R/functions.R")
library(data.table)
options(timeout = 3600)

## Manifest
manifest <- fread(MANIFEST_URL, sep = "\t", quote = "")
fwrite(manifest, file.path(DATA_DIR, "manifest.tsv"), sep = "\t")
both <- manifest[Sex == "both_sexes" & !is.na(`AWS File`) & `AWS File` != "NA"]

## Keyword screening of mental health and wellbeing phenotypes (manual review follows,
## see data/trait_classification.csv)
mh_regex <- "happiness|happy|well[- ]?being|wellbeing|life satisfaction|mood|depression|anxiety|nervous|tired|stress"
candidates <- both[grepl(mh_regex, tolower(`Phenotype Description`)),
                   .(code = `Phenotype Code`, description = `Phenotype Description`)]
fwrite(candidates, file.path(DATA_DIR, "candidate_mental_health_phenotypes.csv"))
message("Candidate mental health phenotypes: ", nrow(candidates))

## Phenotypes to download: all classified traits (mental health candidates and
## air pollution measures) plus the socioeconomic traits
classif <- fread(file.path(DATA_DIR, "trait_classification.csv"), colClasses = list(character = "trait"))
to_get <- unique(c(classif$trait, names(SES)))

for (ph in to_get) {
  qc_file <- file.path(QC_DIR, paste0(ph, "_qc.rds"))
  if (file.exists(qc_file)) next
  raw <- list.files(RAW_DIR, pattern = paste0("^", ph, "\\.gwas"), full.names = TRUE)
  if (!length(raw)) {
    raw <- download_neale(manifest, ph, "both_sexes", RAW_DIR)
    if (is.na(raw)) next
    Sys.sleep(1)
  }
  saveRDS(qc_neale(read_neale(raw[1])), qc_file)
  message("QC done: ", ph)
  gc()
}

done <- sub("_qc\\.rds$", "", list.files(QC_DIR, "_qc\\.rds$"))
message("Phenotypes downloaded and QC'd: ", length(intersect(done, to_get)), " of ", length(to_get))
