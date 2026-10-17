## Step 2. Format QC'd summary statistics as LD Score Regression input
##
## Variants are annotated with rsIDs from the Neale Lab variant file (variants.tsv.bgz);
## A1 = alternate (effect) allele, A2 = reference allele, Z = beta / SE.
## Output: sumstats_ldsc/<code>.txt.gz

source("config.R")
source("R/functions.R")
library(data.table)

map <- fread(file.path(LD_DIR, "variants.tsv.bgz"), select = c("variant", "rsid"))
setnames(map, c("variant", "rsid"), c("variant_id", "rsid"))
map <- unique(map[!is.na(rsid) & rsid != ""], by = "variant_id")

qc_files <- list.files(QC_DIR, "_qc\\.rds$", full.names = TRUE)
for (f in qc_files) {
  out <- file.path(LDSC_IN_DIR, sub("_qc\\.rds$", ".txt.gz", basename(f)))
  if (file.exists(out)) next
  tryCatch(make_ldsc_input(f, map, LDSC_IN_DIR),
           error = function(e) message("Skipped ", basename(f), ": ", conditionMessage(e)))
  gc()
}
message("LDSC input files: ", length(list.files(LDSC_IN_DIR, "\\.txt\\.gz$")))
