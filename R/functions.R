## Helper functions: download, quality control and formatting of Neale Lab summary statistics

library(data.table)
library(httr)

## Download the both-sexes GWAS file of one phenotype from the Neale Lab bucket
download_neale <- function(manifest, phenotype_code, sex = "both_sexes", out_dir = ".") {
  row <- manifest[`Phenotype Code` == phenotype_code & Sex == sex & !is.na(`AWS File`) & `AWS File` != "NA"]
  if (nrow(row) == 0) {
    message("No manifest entry for ", phenotype_code, " (", sex, ")")
    return(NA_character_)
  }
  aws_path <- row$`AWS File`[1]
  file_url <- if (startsWith(aws_path, "http")) aws_path else
    paste0("https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/", aws_path)
  local_file <- file.path(out_dir, basename(aws_path))
  if (file.exists(local_file)) file.remove(local_file)   # avoid truncated files from previous attempts

  res <- tryCatch(GET(file_url, write_disk(local_file, overwrite = TRUE)), error = function(e) NULL)
  if (is.null(res) || http_error(res)) {
    message("Download failed: ", phenotype_code)
    if (file.exists(local_file)) file.remove(local_file)
    return(NA_character_)
  }
  local_file
}

## Read a Neale Lab file
read_neale <- function(path) {
  dt <- fread(path)
  setnames(dt, c("variant", "minor_allele", "minor_AF"), c("variant_id", "a1", "af"), skip_absent = TRUE)
  dt
}

## Quality control:
##  - variants flagged as low confidence by the Neale Lab (minor_AF < 0.001 or
##    expected minor allele count < 25 in the smaller group)
##  - missing beta or standard error
##  - fewer than n_min non-missing samples
qc_neale <- function(dat, n_min = 10000) {
  dat <- as.data.table(dat)
  if (!"pval" %in% names(dat) && "tstat" %in% names(dat)) dat[, pval := 2 * pnorm(-abs(tstat))]
  if ("low_confidence_variant" %in% names(dat)) dat <- dat[low_confidence_variant == FALSE]
  dat <- dat[!is.na(beta) & !is.na(se)]
  if ("n_complete_samples" %in% names(dat)) dat <- dat[n_complete_samples >= n_min]
  dat[]
}

## Write an LDSC input file (SNP, A1, A2, N, P, Z) for one QC'd phenotype.
## Neale Lab effect sizes refer to the alternate allele, so A1 = ALT and A2 = REF.
## map: data.table with columns variant_id (chr:pos:ref:alt) and rsid
make_ldsc_input <- function(qc_file, map, out_dir) {
  x <- as.data.table(readRDS(qc_file))
  x <- merge(x, map, by = "variant_id")
  parts <- tstrsplit(x$variant_id, ":", fixed = TRUE, fill = NA_character_)
  if (length(parts) < 4) stop("Unexpected variant_id format in ", qc_file)
  x[, `:=`(REF = parts[[3]], ALT = parts[[4]])]
  x <- x[!is.na(rsid) & !is.na(REF) & !is.na(ALT) & REF != "" & ALT != ""]

  out <- x[, .(SNP = rsid, A1 = ALT, A2 = REF, N = n_complete_samples, P = pval, Z = beta / se)]
  ph <- sub("_qc\\.rds$", "", basename(qc_file))
  out_path <- file.path(out_dir, paste0(ph, ".txt.gz"))
  fwrite(out, out_path, sep = "\t")
  message(ph, ": ", nrow(out), " variants")
  out_path
}

## Parse per-trait results (h2, intercept, attenuation ratio) from a GenomicSEM ldsc() log.
## Traits appear in the log in the order in which they were passed to ldsc().
parse_ldsc_log <- function(log_file, traits) {
  lg <- readLines(log_file)
  hs <- grep("Heritability Results for trait", lg, fixed = TRUE)
  pick <- function(block, key) { x <- block[grepl(key, block, fixed = TRUE)]; if (length(x)) x[1] else NA_character_ }
  nums <- function(x) {
    if (is.na(x)) return(c(NA, NA))
    s <- sub("^[^:]*:", "", x)
    v <- suppressWarnings(as.numeric(regmatches(s, gregexpr("-?[0-9]+\\.?[0-9]*(e-?[0-9]+)?", s))[[1]]))
    c(v, NA, NA)[1:2]
  }
  rbindlist(lapply(seq_along(hs), function(k) {
    b <- lg[hs[k]:min(hs[k] + 8, length(lg))]
    data.table(trait        = traits[k],
               mean_chi2    = nums(pick(b, "Mean Chi^2"))[1],
               intercept    = nums(pick(b, "Intercept:"))[1],
               intercept_se = nums(pick(b, "Intercept:"))[2],
               ratio        = nums(pick(b, "Ratio:"))[1],
               ratio_se     = nums(pick(b, "Ratio:"))[2],
               h2           = nums(pick(b, "Total Observed Scale h2"))[1],
               h2_se        = nums(pick(b, "Total Observed Scale h2"))[2],
               h2_Z         = nums(pick(b, "h2 Z:"))[1])
  }))
}
