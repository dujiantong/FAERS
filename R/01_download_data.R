#' FAERS Data Download Script
#' @description Downloads and extracts FAERS quarterly data files
#' @note Requires ~10GB free space per quarter

library(RCurl)
library(zip)

# Configuration
BASE_URL <- "https://fis.fda.gov/content/Exports/"
QUARTERS <- list(
  c(2023, 1), c(2023, 2), c(2023, 3), c(2023, 4),
  c(2024, 1), c(2024, 2)
)

DATA_DIR <- "data/raw"
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

#' Download single quarter
download_quarter <- function(year, quarter) {
  q_str <- sprintf("%dq%d", year, quarter)
  url <- sprintf("%sfaers_ascii_%s.zip", BASE_URL, q_str)
  dest_file <- file.path(DATA_DIR, basename(url))
  
  if (file.exists(dest_file)) {
    message(sprintf("Already exists: %s", dest_file))
    return(dest_file)
  }
  
  message(sprintf("Downloading %s...", q_str))
  tryCatch({
    download.file(url, dest_file, mode = "wb", timeout = 600)
    message(sprintf("Downloaded: %s", dest_file))
    return(dest_file)
  }, error = function(e) {
    warning(sprintf("Failed to download %s: %s", q_str, e$message))
    return(NULL)
  })
}

#' Extract and organize files
extract_quarter <- function(zip_file) {
  if (is.null(zip_file)) return(NULL)
  
  # Extract to temp then organize
  temp_dir <- tempdir()
  unzip(zip_file, exdir = temp_dir)
  
  # Find ASCII files (case insensitive)
  ascii_files <- list.files(temp_dir, pattern = "\\.txt$", 
                           ignore.case = TRUE, full.names = TRUE)
  
  # Move to structured directory
  quarter_name <- tools::file_path_sans_ext(basename(zip_file))
  target_dir <- file.path(DATA_DIR, quarter_name)
  dir.create(target_dir, showWarnings = FALSE)
  
  file.copy(ascii_files, file.path(target_dir, basename(ascii_files)))
  unlink(temp_dir, recursive = TRUE)
  
  message(sprintf("Extracted %d files to %s", length(ascii_files), target_dir))
  return(target_dir)
}

# Execute download
for (q in QUARTERS) {
  zip_file <- download_quarter(q[[1]], q[[2]])
  extract_quarter(zip_file)
  Sys.sleep(2)  # Be nice to FDA servers
}
