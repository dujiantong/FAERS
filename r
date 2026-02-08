#' FAERS Data Mining Utilities
#' @description Utility functions for FAERS database analysis
#' @author Your Name
#' @date 2026-02-08

library(data.table)
library(tidyverse)
library(lubridate)
library(openxlsx)

#' Check and install required packages
check_packages <- function() {
  required <- c("data.table", "tidyverse", "lubridate", "openxlsx", 
                "PharmacoVigilance", " disproportionality", "meta", "forestplot")
  installed <- rownames(installed.packages())
  
  for (pkg in required) {
    if (!pkg %in% installed) {
      install.packages(pkg)
    }
  }
}

#' Standardize FAERS quarter format
#' @param year Year (e.g., 2023)
#' @param quarter Quarter (1-4)
#' @return Character string in format "23Q1"
standardize_quarter <- function(year, quarter) {
  sprintf("%02dQ%d", year %% 100, quarter)
}

#' Parse FAERS ASCII files efficiently
#' @param file_path Path to ASCII file
#' @param select_cols Columns to keep (NULL for all)
#' @return data.table
parse_faers_ascii <- function(file_path, select_cols = NULL) {
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  
  # Use fread for speed, handle encoding issues
  dt <- fread(file_path, select = select_cols, 
              encoding = "Latin-1", 
              showProgress = TRUE,
              na.strings = c("", "NA", "NULL"))
  
  # Standardize column names (uppercase, remove special chars)
  setnames(dt, names(dt), toupper(gsub("[^[:alnum:]]", "_", names(dt))))
  
  return(dt)
}

#' Calculate reporting odds ratio (ROR)
#' @param a Drug+Event count
#' @param b Drug count without event  
#' @param c Event count without drug
#' @param d Neither count
#' @return List with ROR, CI, and significance
calc_ror <- function(a, b, c, d) {
  # Add 0.5 continuity correction for zero cells
  a <- ifelse(a == 0, 0.5, a)
  b <- ifelse(b == 0, 0.5, b)
  c <- ifelse(c == 0, 0.5, c)
  d <- ifelse(d == 0, 0.5, d)
  
  ror <- (a * d) / (b * c)
  se_log_ror <- sqrt(1/a + 1/b + 1/c + 1/d)
  ci_lower <- exp(log(ror) - 1.96 * se_log_ror)
  ci_upper <- exp(log(ror) + 1.96 * se_log_ror)
  
  # Chi-square test
  chi_sq <- (abs(a*d - b*c) - 0.5)^2 / ((a+b)*(c+d)*(a+c)*(b+d)) * (a+b+c+d)
  p_value <- 1 - pchisq(chi_sq, 1)
  
  return(list(
    ror = ror,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    significant = ci_lower > 1 & p_value < 0.05
  ))
}

#' Calculate Proportional Reporting Ratio (PRR)
#' @param a Drug+Event count
#' @param b Drug count without event
#' @param c Event count without drug  
#' @param d Neither count
calc_prr <- function(a, b, c, d) {
  prr <- (a / (a + b)) / (c / (c + d))
  se_log_prr <- sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d))
  ci_lower <- exp(log(prr) - 1.96 * se_log_prr)
  ci_upper <- exp(log(prr) + 1.96 * se_log_prr)
  
  # Chi-square with Yates correction
  chi_sq <- ifelse(a > 3 & (a+b) > 100 & (a+c) > 100,
                   ((a*d - b*c)^2 * (a+b+c+d)) / ((a+b)*(c+d)*(a+c)*(b+d)),
                   0)
  
  return(list(
    prr = prr,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    chi_sq = chi_sq
  ))
}
