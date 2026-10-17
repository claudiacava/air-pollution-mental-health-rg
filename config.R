## Paths used by all scripts. Edit before running.

BASE_DIR <- "~/air-pollution-mental-health-rg"  # project folder (outputs are written here)
LD_DIR   <- "~/ldsc/eur_w_ld_chr"                # eur_w_ld_chr LD scores and variants.tsv.bgz
LCV_DIR  <- "~/LCV/R"                            # MomentFunctions.R and RunLCV.R (github.com/lukejoconnor/LCV)

## Folders created inside BASE_DIR
DATA_DIR    <- file.path(BASE_DIR, "data")
RAW_DIR     <- file.path(BASE_DIR, "sumstats_raw")
QC_DIR      <- file.path(BASE_DIR, "sumstats_qc")
LDSC_IN_DIR <- file.path(BASE_DIR, "sumstats_ldsc")
INIT_DIR    <- file.path(BASE_DIR, "results", "initial")
MAIN_DIR    <- file.path(BASE_DIR, "results", "main")
for (d in c(RAW_DIR, QC_DIR, LDSC_IN_DIR, INIT_DIR, MAIN_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## Neale Lab round 2 manifest ("Manifest 201807" sheet)
MANIFEST_URL <- "https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/export?gid=178908679&format=tsv"

## Socioeconomic traits
SES <- c("189_irnt" = "Townsend deprivation index",
         "6138_1"   = "College or university degree",
         "738"      = "Average total household income before tax",
         "20118_5"  = "Urban residence (England/Wales, urban - less sparse)")
