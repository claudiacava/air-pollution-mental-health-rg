## Step 4. Main analyses
##
##  - joint LDSC (stand = TRUE) on the selected air pollution, mental health,
##    socioeconomic and sensitivity traits
##  - T1: SNP heritability, intercept and attenuation ratio per trait
##  - T2: genetic correlations air pollution x mental health (Supplementary Table S3)
##  - T3: genetic correlations with socioeconomic traits (Supplementary Table S5)
##  - T4: sensitivity set, help-seeking / treatment phenotypes (Supplementary Table S4)
##  - T5: partial genetic correlations conditional on socioeconomic traits (Supplementary Table S6)
##  - T6: mean genetic correlation by mental health subdomain
##  - T7: latent causal variable (LCV) analysis (Supplementary Table S7)
## Each step saves its result, so the script can be restarted.
## Output: results/main/

source("config.R")
source("R/functions.R")
suppressPackageStartupMessages({ library(data.table); library(GenomicSEM); library(MASS) })

N_DRAWS <- 2000
set.seed(2026)
OUT <- MAIN_DIR
capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))

## 0. Traits ------------------------------------------------------------------
classif <- fread(file.path(DATA_DIR, "trait_classification.csv"), colClasses = list(character = "trait"))
main <- classif[include == TRUE]
sens <- classif[sensitivity == TRUE]
poll <- main[domain == "air_pollution", trait]
mh   <- main[domain != "air_pollution", trait]
message("Main set: ", length(poll), " air pollution and ", length(mh), " mental health traits; ",
        "sensitivity set: ", nrow(sens), "; socioeconomic: ", length(SES))

## 1. LDSC --------------------------------------------------------------------
traits <- c(poll, mh, names(SES), sens$trait)
files  <- file.path(LDSC_IN_DIR, paste0(traits, ".txt.gz"))
stopifnot(all(file.exists(files)))
K <- length(traits)
ldsc_rds <- file.path(OUT, "ldsc_main.rds")
if (!file.exists(ldsc_rds)) {
  old <- setwd(OUT)
  ld <- ldsc(traits = files, sample.prev = rep(NA, K), population.prev = rep(NA, K),
             ld = LD_DIR, wld = LD_DIR, trait.names = traits, stand = TRUE, ldsc.log = "ldsc_main")
  setwd(old)
  saveRDS(ld, ldsc_rds)
} else ld <- readRDS(ldsc_rds)

R  <- ld$S_Stand; dimnames(R) <- list(traits, traits)
VS <- ld$V_Stand
vi <- function(i, j) { if (i < j) { t <- i; i <- j; j <- t }; (j - 1) * K - (j - 1) * (j - 2) / 2 + (i - j + 1) }
se_rg <- function(i, j) sqrt(VS[vi(i, j), vi(i, j)])
idx <- setNames(seq_len(K), traits)

## 2. T1: per-trait estimates --------------------------------------------------
per_trait <- parse_ldsc_log(file.path(OUT, "ldsc_main_ldsc.log"), traits)
per_trait <- merge(per_trait, classif[, .(trait, description, domain, subdomain)], by = "trait", all.x = TRUE)
per_trait[trait %in% names(SES), `:=`(description = SES[trait], domain = "socioeconomic")]
fwrite(per_trait, file.path(OUT, "T1_trait_h2_intercept_ratio.csv"))

## 3. T2-T4: genetic correlations ----------------------------------------------
pair_tab <- function(A, B) rbindlist(lapply(A, function(a) rbindlist(lapply(B, function(b) {
  r <- R[a, b]; s <- se_rg(idx[a], idx[b])
  data.table(trait1 = a, trait2 = b, rg = r, se = s, z = r / s, p = 2 * pnorm(-abs(r / s)))
}))))

rg_main <- pair_tab(poll, mh)
rg_main[, FDR := p.adjust(p, "BH")]
rg_main <- merge(rg_main, classif[, .(trait2 = trait, description2 = description, subdomain)], by = "trait2")
setorder(rg_main, p)
fwrite(rg_main, file.path(OUT, "T2_rg_pollution_mental_health.csv"))

rg_ses <- pair_tab(names(SES), c(poll, mh))
rg_ses[, FDR := p.adjust(p, "BH")]
fwrite(rg_ses, file.path(OUT, "T3_rg_socioeconomic.csv"))

rg_sens <- pair_tab(poll, sens$trait)
rg_sens[, FDR := p.adjust(p, "BH")]
fwrite(rg_sens, file.path(OUT, "T4_rg_sensitivity_help_seeking.csv"))
message("Air pollution x mental health pairs with FDR < 0.05: ", sum(rg_main$FDR < 0.05), " of ", nrow(rg_main))

## 4. Draws from the sampling distribution of the standardized estimates ------
vech_R <- R[lower.tri(R, diag = TRUE)]
draws <- MASS::mvrnorm(N_DRAWS, mu = vech_R, Sigma = VS)
to_mat <- function(v) {
  M <- matrix(0, K, K); M[lower.tri(M, diag = TRUE)] <- v
  M <- M + t(M) - diag(diag(M)); dimnames(M) <- list(traits, traits); M
}

## 5. T5: partial genetic correlations ----------------------------------------
partial <- function(M, x, y, z) {
  s <- M[c(x, y, z), c(x, y, z)]
  P <- tryCatch(solve(s), error = function(e) NULL)
  if (is.null(P)) return(NA_real_)
  -P[1, 2] / sqrt(P[1, 1] * P[2, 2])
}
sets <- list(Townsend = "189_irnt", Degree = "6138_1", Townsend_Degree = c("189_irnt", "6138_1"),
             All_SES = names(SES))
part_rds <- file.path(OUT, "partial_rg.rds")
if (!file.exists(part_rds)) {
  Ms <- lapply(seq_len(N_DRAWS), function(d) to_mat(draws[d, ]))
  part <- rbindlist(lapply(names(sets), function(sn) rbindlist(lapply(poll, function(a) rbindlist(lapply(mh, function(b) {
    est <- partial(R, a, b, sets[[sn]])
    bs  <- vapply(Ms, partial, numeric(1), x = a, y = b, z = sets[[sn]])
    data.table(adjustment = sn, pollutant = a, trait = b, rg = R[a, b], partial_rg = est,
               ci_low = quantile(bs, 0.025, na.rm = TRUE), ci_high = quantile(bs, 0.975, na.rm = TRUE),
               p_approx = 2 * min(mean(bs > 0, na.rm = TRUE), mean(bs < 0, na.rm = TRUE)))
  }))))))
  part[, attenuation := 1 - partial_rg / rg]
  saveRDS(part, part_rds); rm(Ms); gc()
} else part <- readRDS(part_rds)
fwrite(part, file.path(OUT, "T5_partial_rg_socioeconomic.csv"))

## 6. T6: subdomain means -----------------------------------------------------
subd <- classif[trait %in% mh, .(trait, subdomain)]
sub_summary <- rbindlist(lapply(poll, function(a) rbindlist(lapply(unique(subd$subdomain), function(sd) {
  tr <- subd[subdomain == sd, trait]
  ii <- vapply(tr, function(b) vi(idx[a], idx[b]), numeric(1))
  bs <- if (length(ii) > 1) rowMeans(draws[, ii, drop = FALSE]) else draws[, ii]
  data.table(pollutant = a, subdomain = sd, n_traits = length(tr), mean_rg = mean(R[a, tr]),
             ci_low = quantile(bs, 0.025), ci_high = quantile(bs, 0.975))
}))))
fwrite(sub_summary, file.path(OUT, "T6_rg_by_subdomain.csv"))

## 7. T7: LCV for traits with Z(h2) >= 7 ---------------------------------------
lcv_poll <- main[lcv_eligible == TRUE & domain == "air_pollution", trait]
lcv_mh   <- main[lcv_eligible == TRUE & domain != "air_pollution", trait]
source(file.path(LCV_DIR, "MomentFunctions.R")); source(file.path(LCV_DIR, "RunLCV.R"))
ld_files <- list.files(LD_DIR, pattern = "l2\\.ldscore(\\.gz)?$", full.names = TRUE)
ld_dt <- rbindlist(lapply(ld_files, fread), fill = TRUE)[, .(SNP, CHR, BP, L2)]
setorder(ld_dt, CHR, BP)
read_z <- function(t) {
  d <- fread(file.path(LDSC_IN_DIR, paste0(t, ".txt.gz")), showProgress = FALSE, select = c("SNP", "A1", "A2", "Z"))
  d[SNP %in% ld_dt$SNP]
}
lcv_rds <- file.path(OUT, "lcv_partial.rds")
res <- if (file.exists(lcv_rds)) readRDS(lcv_rds) else list()
cache <- list()
for (a in lcv_poll) {
  cache[[a]] <- read_z(a)
  for (b in lcv_mh) {
    key <- paste(a, b, sep = "|"); if (!is.null(res[[key]])) next
    if (is.null(cache[[b]])) cache[[b]] <- read_z(b)
    dat <- merge(merge(ld_dt, cache[[a]], by = "SNP"), cache[[b]], by = "SNP", suffixes = c(".x", ".y"))
    setorder(dat, CHR, BP)
    mm <- which(dat$A1.x != dat$A1.y); if (length(mm)) dat[mm, Z.y := -Z.y]
    out <- tryCatch(RunLCV(dat$L2, dat$Z.x, dat$Z.y), error = function(e) NULL)
    res[[key]] <- data.table(pollutant = a, trait = b,
      gcp    = if (is.null(out)) NA else out$gcp.pm,  gcp_se = if (is.null(out)) NA else out$gcp.pse,
      rho    = if (is.null(out)) NA else out$rho.est, rho_se = if (is.null(out)) NA else out$rho.err,
      p_gcp0 = if (is.null(out)) NA else out$pval.gcpzero.2tailed)
    saveRDS(res, lcv_rds)
  }
}
lcv <- rbindlist(res)[!is.na(p_gcp0)]
lcv[, FDR := p.adjust(p_gcp0, "BH")]
lcv <- merge(lcv, rg_main[, .(pollutant = trait1, trait = trait2, rg_ldsc = rg, rg_FDR = FDR)],
             by = c("pollutant", "trait"), all.x = TRUE)
setorder(lcv, p_gcp0)
fwrite(lcv, file.path(OUT, "T7_lcv.csv"))
message("LCV tests: ", nrow(lcv), "; FDR < 0.05: ", sum(lcv$FDR < 0.05))
