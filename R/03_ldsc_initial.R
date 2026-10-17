## Step 3. Initial LDSC on all candidate traits (used for trait selection)
##
## - joint LDSC on all classified traits with an input file: SNP heritability, Z(h2),
##   intercept and attenuation ratio per trait (Supplementary Table S1)
## - genetic correlations among the air pollution measures, raw and IRNT (Supplementary Table S2)
## The selection itself (include / lcv_eligible / sensitivity) is recorded in data/trait_classification.csv:
## Z(h2) >= 4 for genetic correlations, Z(h2) >= 7 for LCV, one version and one measure per exposure.
## Output: results/initial/

source("config.R")
source("R/functions.R")
library(data.table)
library(GenomicSEM)

classif <- fread(file.path(DATA_DIR, "trait_classification.csv"), colClasses = list(character = "trait"))
traits <- classif$trait[file.exists(file.path(LDSC_IN_DIR, paste0(classif$trait, ".txt.gz")))]
files  <- file.path(LDSC_IN_DIR, paste0(traits, ".txt.gz"))
K <- length(traits)
message("Initial LDSC on ", K, " traits")

## 1. All candidate traits: per-trait estimates
old <- setwd(INIT_DIR)
ld_all <- ldsc(traits = files, sample.prev = rep(NA, K), population.prev = rep(NA, K),
               ld = LD_DIR, wld = LD_DIR, trait.names = traits, ldsc.log = "ldsc_initial")
setwd(old)
saveRDS(ld_all, file.path(INIT_DIR, "ldsc_initial.rds"))

per_trait <- parse_ldsc_log(file.path(INIT_DIR, "ldsc_initial_ldsc.log"), traits)
n_gwas <- as.numeric(ld_all$N)[(seq_len(K) - 1) * K - (seq_len(K) - 1) * (seq_len(K) - 2) / 2 + 1]
per_trait[, N := round(n_gwas)]
per_trait <- merge(classif, per_trait, by = "trait", all.x = TRUE)
fwrite(per_trait, file.path(INIT_DIR, "S1_candidate_traits_ldsc.csv"))

## 2. Air pollution measures: genetic correlations with standard errors
poll <- traits[grepl("^240", traits)]
P <- length(poll)
old <- setwd(INIT_DIR)
ld_poll <- ldsc(traits = file.path(LDSC_IN_DIR, paste0(poll, ".txt.gz")),
                sample.prev = rep(NA, P), population.prev = rep(NA, P),
                ld = LD_DIR, wld = LD_DIR, trait.names = poll, stand = TRUE, ldsc.log = "ldsc_pollution")
setwd(old)
R  <- ld_poll$S_Stand; dimnames(R) <- list(poll, poll)
se <- matrix(NA_real_, P, P, dimnames = list(poll, poll))
k <- 0
for (j in seq_len(P)) for (i in j:P) { k <- k + 1; se[i, j] <- se[j, i] <- sqrt(ld_poll$V_Stand[k, k]) }

rg_long <- rbindlist(lapply(seq_len(P - 1), function(i) rbindlist(lapply((i + 1):P, function(j)
  data.table(trait1 = poll[i], trait2 = poll[j], rg = R[i, j], se = se[i, j])))))
fwrite(rg_long, file.path(INIT_DIR, "S2_rg_pollution_measures.csv"))
message("Done. See ", INIT_DIR)
