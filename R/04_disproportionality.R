#' Disproportionality Analysis for FAERS
#' @description Implements ROR, PRR, and BCPNN algorithms

source("R/utils.R")

PROCESSED_DIR <- "data/processed"
OUTPUT_DIR <- "output/tables"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load processed data
message("Loading data...")
DEMO <- readRDS(file.path(PROCESSED_DIR, "DEMO_combined.rds"))
DRUG <- readRDS(file.path(PROCESSED_DIR, "DRUG_combined.rds"))
REAC <- readRDS(file.path(PROCESSED_DIR, "REAC_combined.rds"))

# Focus on Primary Suspect drugs only for signal detection
DRUG_PS <- DRUG[role_code %in% c("PS", "SS")]

#' Build contingency table for specific drug-event pair
build_contingency <- function(target_drug, target_event, 
                             drug_data = DRUG_PS, 
                             event_data = REAC) {
  
  # Standardize inputs
  target_drug <- toupper(target_drug)
  target_event <- toupper(target_event)
  
  # Find cases with target drug
  drug_cases <- drug_data[grepl(target_drug, drug_name, ignore.case = TRUE), 
                         unique(primaryid)]
  
  # Find cases with target event
  event_cases <- event_data[pt == target_event, unique(primaryid)]
  
  # Contingency table
  a <- length(intersect(drug_cases, event_cases))  # Drug + Event
  b <- length(setdiff(drug_cases, event_cases))    # Drug only
  c <- length(setdiff(event_cases, drug_cases))    # Event only
  d <- length(setdiff(
    unique(drug_data$primaryid), 
    union(drug_cases, event_cases)
  ))                                               # Neither
  
  return(list(a = a, b = b, c = c, d = d))
}

#' Screen all events for a specific drug
screen_drug <- function(drug_name, min_cases = 3) {
  message(sprintf("Screening: %s", drug_name))
  
  # Get all cases for this drug
  drug_cases <- DRUG_PS[grepl(drug_name, drug_name, ignore.case = TRUE), 
                       unique(primaryid)]
  
  if (length(drug_cases) < min_cases) {
    message(sprintf("Insufficient cases for %s (%d)", drug_name, length(drug_cases)))
    return(NULL)
  }
  
  # Get all events for these cases
  events <- REAC[primaryid %in% drug_cases, .N, by = pt][order(-N)]
  
  # Calculate disproportionality for each event
  results <- lapply(events$pt, function(evt) {
    cont <- build_contingency(drug_name, evt)
    
    if (cont$a < min_cases) return(NULL)
    
    ror_res <- calc_ror(cont$a, cont$b, cont$c, cont$d)
    prr_res <- calc_prr(cont$a, cont$b, cont$c, cont$d)
    
    data.table(
      drug = drug_name,
      event = evt,
      case_count = cont$a,
      ror = ror_res$ror,
      ror_ci_lower = ror_res$ci_lower,
      ror_ci_upper = ror_res$ci_upper,
      ror_significant = ror_res$significant,
      prr = prr_res$prr,
      prr_ci_lower = prr_res$ci_lower,
      prr_ci_upper = prr_res$ci_upper,
      prr_chi_sq = prr_res$chi_sq
    )
  })
  
  rbindlist(results)
}

# Example: Analyze a list of drugs of interest
DRUGS_OF_INTEREST <- c("METFORMIN", "LISINOPRIL", "ATORVASTATIN")

all_signals <- lapply(DRUGS_OF_INTEREST, screen_drug)
combined_signals <- rbindlist(all_signals)

# Filter significant signals (ROR > 1, lower CI > 1, at least 3 cases)
significant_signals <- combined_signals[
  ror_significant == TRUE & 
  case_count >= 3 &
  ror_ci_lower > 2  # Conservative threshold
][order(-ror)]

# Save results
fwrite(combined_signals, file.path(OUTPUT_DIR, "disproportionality_all.csv"))
fwrite(significant_signals, file.path(OUTPUT_DIR, "significant_signals.csv"))

# Generate summary statistics
summary_stats <- data.table(
  total_drugs_analyzed = length(DRUGS_OF_INTEREST),
  total_drug_event_pairs = nrow(combined_signals),
  significant_signals = nrow(significant_signals),
  median_ror = median(combined_signals$ror, na.rm = TRUE)
)

print(summary_stats)
