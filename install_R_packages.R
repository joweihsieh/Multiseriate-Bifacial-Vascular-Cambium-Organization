#!/usr/bin/env Rscript

# Install the R packages used in the
# Multiseriate-Bifacial-Vascular-Cambium-Organization repository.
#
# Usage:
#   Rscript install_R_packages.R
#
# Optional environment variables:
#   CRAN_REPO=https://cloud.r-project.org
#   R_PACKAGE_LIBRARY=/path/to/writable/R/library
#   R_INSTALL_NCPUS=4
#
# Notes:
# - The package versions below are the versions recorded for the analyses.
# - grid is distributed with R and is therefore checked, not installed.
# - BiocManager and remotes are installation helpers. Their versions were not
#   specified in the original package inventory, so compatible versions are
#   installed from CRAN when needed.

options(stringsAsFactors = FALSE)

cran_repo <- Sys.getenv("CRAN_REPO", unset = "https://cloud.r-project.org")
options(repos = c(CRAN = cran_repo))

ncpus_text <- Sys.getenv("R_INSTALL_NCPUS", unset = "")
if (nzchar(ncpus_text)) {
  ncpus <- suppressWarnings(as.integer(ncpus_text))
} else {
  detected <- suppressWarnings(parallel::detectCores(logical = FALSE))
  if (is.na(detected) || detected < 2L) {
    ncpus <- 1L
  } else {
    ncpus <- min(8L, detected - 1L)
  }
}
if (is.na(ncpus) || ncpus < 1L) ncpus <- 1L
options(Ncpus = ncpus)

# Use a user-specified library when requested. Otherwise, use the active R
# library. If that library is not writable, fall back to a local R_library.
requested_library <- Sys.getenv("R_PACKAGE_LIBRARY", unset = "")

if (nzchar(requested_library)) {
  target_library <- path.expand(requested_library)
  dir.create(target_library, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(target_library, .libPaths()))
} else {
  active_library <- .libPaths()[1L]
  if (!dir.exists(active_library)) {
    dir.create(active_library, recursive = TRUE, showWarnings = FALSE)
  }

  if (file.access(active_library, 2L) == 0L) {
    target_library <- active_library
  } else {
    target_library <- file.path(getwd(), "R_library")
    dir.create(target_library, recursive = TRUE, showWarnings = FALSE)
    .libPaths(c(target_library, .libPaths()))
    warning(
      "The active R library was not writable. Packages will be installed in: ",
      normalizePath(target_library, mustWork = FALSE),
      call. = FALSE
    )
  }
}

expected_r_version <- "4.3.3"

# install_version is the version string used by CRAN archives.
# expected_version is the normalized value returned by packageVersion().
package_specification <- data.frame(
  package = c(
    "RColorBrewer",
    "magrittr",
    "Matrix",
    "MASS",
    "data.table",
    "readxl",
    "scales",
    "ggplot2",
    "pheatmap",
    "optparse",
    "irlba",
    "RANN",
    "dplyr",
    "igraph",
    "mclust",
    "openxlsx",
    "Seurat"
  ),
  install_version = c(
    "1.1-3",
    "2.0.5",
    "1.6-5",
    "7.3-60.0.1",
    "1.18.4",
    "1.5.0",
    "1.4.0",
    "4.0.3",
    "1.0.13",
    "1.8.2",
    "2.3.7",
    "2.6.2",
    "1.2.1",
    "2.3.2",
    "6.1.2",
    "4.2.8.1",
    "5.5.0"
  ),
  expected_version = c(
    "1.1.3",
    "2.0.5",
    "1.6.5",
    "7.3.60.0.1",
    "1.18.4",
    "1.5.0",
    "1.4.0",
    "4.0.3",
    "1.0.13",
    "1.8.2",
    "2.3.7",
    "2.6.2",
    "1.2.1",
    "2.3.2",
    "6.1.2",
    "4.2.8.1",
    "5.5.0"
  ),
  stringsAsFactors = FALSE
)

base_package_specification <- data.frame(
  package = "grid",
  expected_version = "4.3.3",
  note = "Bundled with R; not installed separately",
  stringsAsFactors = FALSE
)

normalize_version <- function(x) {
  gsub("-", ".", as.character(x), fixed = TRUE)
}

version_matches <- function(installed_version, expected_version) {
  if (is.na(installed_version) || !nzchar(installed_version)) return(FALSE)
  identical(
    normalize_version(installed_version),
    normalize_version(expected_version)
  )
}

get_installed_version <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(package))
}

install_helper <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    message("Installing installation helper: ", package)
    utils::install.packages(
      package,
      lib = target_library,
      repos = cran_repo,
      dependencies = c("Depends", "Imports", "LinkingTo")
    )
  }

  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Failed to install required helper package: ", package, call. = FALSE)
  }
}

cat("Installation started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R version:", R.version.string, "\n")
cat("Expected R version:", expected_r_version, "\n")
cat("Target R library:", normalizePath(target_library, mustWork = FALSE), "\n")
cat("CRAN repository:", cran_repo, "\n")
cat("Parallel installation jobs:", ncpus, "\n\n")

if (!identical(as.character(getRversion()), expected_r_version)) {
  warning(
    "This installation is being run with R ", as.character(getRversion()),
    ", whereas the recorded grid version indicates R ", expected_r_version,
    ". Exact reproduction may require R ", expected_r_version, ".",
    call. = FALSE
  )
}

install_helper("BiocManager")
install_helper("remotes")

installation_records <- vector("list", nrow(package_specification))

for (i in seq_len(nrow(package_specification))) {
  package_name <- package_specification$package[[i]]
  install_version <- package_specification$install_version[[i]]
  expected_version <- package_specification$expected_version[[i]]
  before_version <- get_installed_version(package_name)

  if (version_matches(before_version, expected_version)) {
    message(
      "Already installed: ", package_name, " ", before_version
    )
    install_status <- "already_exact"
    install_error <- NA_character_
  } else {
    message(
      "Installing ", package_name, " ", install_version,
      if (!is.na(before_version)) paste0(" (current: ", before_version, ")") else ""
    )

    install_error <- NA_character_
    install_status <- tryCatch(
      {
        remotes::install_version(
          package = package_name,
          version = install_version,
          repos = cran_repo,
          lib = target_library,
          dependencies = c("Depends", "Imports", "LinkingTo"),
          upgrade = "never",
          force = TRUE,
          quiet = FALSE
        )
        "install_command_completed"
      },
      error = function(e) {
        install_error <<- conditionMessage(e)
        "installation_failed"
      }
    )
  }

  after_version <- get_installed_version(package_name)
  exact_match <- version_matches(after_version, expected_version)

  installation_records[[i]] <- data.frame(
    package = package_name,
    requested_install_version = install_version,
    expected_reported_version = expected_version,
    version_before = before_version,
    version_after = after_version,
    exact_match = exact_match,
    install_status = install_status,
    error = install_error,
    stringsAsFactors = FALSE
  )
}

installation_report <- do.call(rbind, installation_records)

# Check the base/recommended package that is supplied with R.
grid_version <- get_installed_version("grid")
grid_report <- data.frame(
  package = "grid",
  requested_install_version = NA_character_,
  expected_reported_version = base_package_specification$expected_version,
  version_before = grid_version,
  version_after = grid_version,
  exact_match = version_matches(
    grid_version,
    base_package_specification$expected_version
  ),
  install_status = "bundled_with_R_not_installed",
  error = NA_character_,
  stringsAsFactors = FALSE
)

installation_report <- rbind(installation_report, grid_report)
installation_report <- installation_report[
  match(
    c(package_specification$package, "grid"),
    installation_report$package
  ),
  ,
  drop = FALSE
]

report_file <- file.path(getwd(), "R_package_installation_report.tsv")
utils::write.table(
  installation_report,
  file = report_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

session_file <- file.path(getwd(), "R_sessionInfo_after_install.txt")
session_lines <- c(
  paste("Installation completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R version:", R.version.string),
  paste("Expected R version:", expected_r_version),
  paste("Target R library:", normalizePath(target_library, mustWork = FALSE)),
  paste("CRAN repository:", cran_repo),
  paste(
    "BiocManager version:",
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      as.character(utils::packageVersion("BiocManager"))
    } else {
      "NOT INSTALLED"
    }
  ),
  paste(
    "Bioconductor release:",
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      as.character(BiocManager::version())
    } else {
      "UNKNOWN"
    }
  ),
  "",
  "Pinned package verification:",
  capture.output(print(installation_report, row.names = FALSE)),
  "",
  "Full session information:",
  capture.output(utils::sessionInfo())
)
writeLines(session_lines, session_file)

cat("\nPackage verification:\n")
print(installation_report, row.names = FALSE)
cat("\nInstallation report:", normalizePath(report_file, mustWork = FALSE), "\n")
cat("Session information:", normalizePath(session_file, mustWork = FALSE), "\n")

failed_rows <- installation_report[
  !installation_report$exact_match |
    installation_report$install_status == "installation_failed",
  ,
  drop = FALSE
]

if (nrow(failed_rows) > 0L) {
  cat("\nThe following packages did not match the recorded versions:\n")
  print(failed_rows, row.names = FALSE)
  stop(
    "Installation finished with one or more version mismatches. ",
    "See R_package_installation_report.tsv for details.",
    call. = FALSE
  )
}

cat("\nAll recorded R package versions were verified successfully.\n")
