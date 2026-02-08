#' FAERS Data Preprocessing
#' @description Merges quarterly data and performs initial cleaning

source("R/utils.R")

PROCESSED_DIR <- "data/processed"
dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)

# Define file patterns for each domain
DOMAIN_PATTERNS <- list(
  DEMO = "^DEMO.*\\.txt$",      # Demographics
  DRUG = "^DRUG.*\\.txt$",      # Drug information  
  REAC = "^REAC.*\\.txt$",      # Reactions (ADEs)
  OUTC = "^OUTC.*\\.txt$",      # Outcomes
  INDI = "^INDI.*\\.txt$",      # Indications
  THER = "^THER.*\\.txt$",      # Therapy dates
  RPSR = "^RPSR.*\\.txt$"       # Report sources
)

#' Load and combine all quarters for a domain
load_domain <- function(pattern, raw_dir = "data/raw") {
  quarters <- list.dirs(raw_dir, full.names = TRUE, recursive = FALSE)
  
  all_data <- lapply(quarters, function(q_dir) {
    files <- list.files(q_dir, pattern = pattern, ignore.case = TRUE)
    if (length(files) == 0) return(NULL)
    
    # Read first matching file
    file_path <- file.path(q_dir, files[1])
    dt <- parse_faers_ascii(file_path)
    
    # Add quarter identifier
    dt[, quarter := basename(q_dir)]
    return(dt)
  })
  
  # Combine all quarters
  combined <- rbindlist(all_data, use.names = TRUE, fill = TRUE)
  return(combined)
}

message("Loading DEMO data...")
DEMO <- load_domain(DOMAIN_PATTERNS$DEMO)

message("Loading DRUG data...")
DRUG <- load_domain(DOMAIN_PATTERNS$DRUG)

message("Loading REAC data...")
REAC <- load_domain(DOMAIN_PATTERNS$REAC)

message("Loading OUTC data...")
OUTC <- load_domain(DOMAIN_PATTERNS$OUTC)

# Data cleaning pipeline
message("Cleaning DEMO...")
DEMO <- DEMO %>%
  .[, `:=`(
    primaryid = as.character(PRIMARYID),
    caseid = as.character(CASEID),
    age = suppressWarnings(as.numeric(AGE)),
    age_code = toupper(AGE_COD),
    gender = toupper(GNDR_COD),
    reporter_country = toupper(REPORTER_COUNTRY),
    event_date = parse_date_time(EVENT_DT, orders = c("Ymd", "Ym", "Y")),
    receipt_date = parse_date_time(RECEIPTDATE, orders = c("Ymd")),
    occr_country = toupper(OCCR_COUNTRY)
  )] %>%
  .[!is.na(primaryid)]

message("Cleaning DRUG...")
DRUG <- DRUG %>%
  .[, `:=`(
    primaryid = as.character(PRIMARYID),
    drug_seq = as.integer(ROLE_COD),
    drug_name = toupper(DRUGNAME),
    prod_ai = toupper(PROD_AI),
    route = toupper(ROUTE),
    dose_vbm = toupper(DOSE_VBM),
    drug_recurrence = toupper(DRUG_REC_ACT),
    role_code = toupper(ROLE_COD)  # PS=Primary Suspect, SS=Secondary Suspect, C=Concomitant, I=Interacting
  )] %>%
  .[!is.na(primaryid) & !is.na(drug_name)]

message("Cleaning REAC...")
REAC <- REAC %>%
  .[, `:=`(
    primaryid = as.character(PRIMARYID),
    pt = toupper(PT),  # Preferred Term (MedDRA)
    drug_rec_act = toupper(DRUG_REC_ACT)
  )] %>%
  .[!is.na(primaryid) & !is.na(pt)]

# Save processed files
message("Saving processed data...")
saveRDS(DEMO, file.path(PROCESSED_DIR, "DEMO_combined.rds"))
saveRDS(DRUG, file.path(PROCESSED_DIR, "DRUG_combined.rds"))
saveRDS(REAC, file.path(PROCESSED_DIR, "REAC_combined.rds"))
saveRDS(OUTC, file.path(PROCESSED_DIR, "OUTC_combined.rds"))

# Create lookup tables
message("Creating reference tables...")
DRUG_NAMES <- unique(DRUG[, .(drug_name, prod_ai)])
PT_NAMES <- unique(REAC[, .(pt)])

fwrite(DRUG_NAMES, file.path(PROCESSED_DIR, "drug_dictionary.csv"))
fwrite(PT_NAMES, file.path(PROCESSED_DIR, "event_dictionary.csv"))

message("Preprocessing complete!")
