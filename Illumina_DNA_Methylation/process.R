# =========================================================================
# Methylation array processing pipeline (processed-matrix branch)
# -------------------------------------------------------------------------
# Purpose:
#   Many GEO series with methylation array data do NOT include raw IDAT
#   files -- instead they provide already-extracted signal intensities
#   (Methylated/Unmethylated, or equivalently SignalA/SignalB) and
#   detection p-values as plain CSV tables. This script builds a
#   minfi-based QC and normalization pipeline that works directly from
#   those tables, so the same code can be reused across many GEO series
#   that each name their files (and even use different array platforms)
#   a bit differently.
#
#   Three Illumina array platforms are supported: EPIC, HumanMethylation450
#   ("450k"), and HumanMethylation27 ("27k"). You choose which one applies
#   for a given dataset with the array_type argument each time you run the
#   pipeline -- only that platform's manifest/annotation packages need to
#   be installed and get loaded.
#
# How to use it:
#   1. Call download_geo_supplementary() once per GEO accession, to fetch
#      that series' supplementary files.
#   2. Call inspect_supp_file() to peek at the top of a downloaded file
#      before writing patterns for it -- column-naming conventions differ
#      between GEO series, so it's worth a quick look each time you use
#      this on a new dataset. This is also how you confirm whether a
#      dataset uses "Meth"/"Unmeth" columns or "SignalA"/"SignalB" columns.
#   3. Call run_methylation_pipeline(), giving it the accession, the array
#      platform, and two file-NAME patterns (which file to use), plus, if
#      needed, the column-NAME suffixes used inside that file. The rest of
#      the pipeline (QC, filtering, normalization, annotation) runs the
#      same way regardless of the dataset or platform.
# =========================================================================

# -------------------------------------------------------------------------
# One-time setup: install required packages
# -------------------------------------------------------------------------
# Uncomment and run this block once per machine (not needed on every run).
# Most of these packages come from Bioconductor rather than CRAN, so they're
# installed through BiocManager rather than install.packages() directly.
# You only need the manifest/annotation packages for the array platform(s)
# you actually work with -- feel free to comment out the ones you don't need.
#
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install(c(
#   "GEOquery",    # download files from GEO
#   "minfi",       # methylation array data structures, QC, normalization
#   "minfiData",   # example datasets used in minfi's own documentation/testing
# 
#   # -- EPIC array --
#   "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
#   "IlluminaHumanMethylationEPICmanifest",
# 
#   # -- HumanMethylation450 ("450k") array --
#   "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#   "IlluminaHumanMethylation450kmanifest",
# 
#   # -- HumanMethylation27 ("27k") array --
#   "IlluminaHumanMethylation27kanno.ilmn12.hg19",
#   "IlluminaHumanMethylation27kmanifest"
# ))
# 
# # maxprobes isn't on CRAN or Bioconductor -- it's installed from GitHub
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("markgene/maxprobes")
#
# # tidyverse is on CRAN
# install.packages("tidyverse")

library(tidyverse)      # data wrangling and plotting
library(GEOquery)       # for downloading files from GEO
library(minfi)          # methylation array data structures, QC, normalization
library(maxprobes)      # list of known cross-reactive ("blacklisted") probes
# Array-specific manifest/annotation packages are loaded on demand by
# load_array_packages(), based on the array_type you specify at run time,
# rather than all being loaded up front here.

# A colorblind-friendly palette (Okabe & Ito 2008), used for every plot below
okabe_ito <- c(
  orange   = "#E69F00",
  sky_blue = "#56B4E9",
  green    = "#009E73",
  yellow   = "#F0E442",
  blue     = "#0072B2",
  red      = "#D55E00",
  pink     = "#CC79A7",
  black    = "#000000"
)

# -------------------------------------------------------------------------
# Array platform configuration
# -------------------------------------------------------------------------
# Each supported platform needs a slightly different manifest package,
# annotation package, and internal minfi array name. This function is the
# single place that knows those details, so the rest of the pipeline can
# just ask for "EPIC", "450k", or "27k" and not worry about package names.
get_array_config <- function(array_type) {
  array_type <- tolower(array_type)
  
  config <- switch(array_type,
                   "epic" = list(
                     array_name          = "IlluminaHumanMethylationEPIC",
                     annotation_version   = "ilm10b4.hg19",
                     manifest_package     = "IlluminaHumanMethylationEPICmanifest",
                     annotation_package   = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
                     maxprobes_array_type = "EPIC"
                   ),
                   "450k" = list(
                     array_name          = "IlluminaHumanMethylation450k",
                     annotation_version   = "ilmn12.hg19",
                     manifest_package     = "IlluminaHumanMethylation450kmanifest",
                     annotation_package   = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                     maxprobes_array_type = "450K"
                   ),
                   "27k" = list(
                     array_name          = "IlluminaHumanMethylation27k",
                     annotation_version   = "ilmn12.hg19",
                     manifest_package     = "IlluminaHumanMethylation27kmanifest",
                     annotation_package   = "IlluminaHumanMethylation27kanno.ilmn12.hg19",
                     # maxprobes doesn't publish a cross-reactive probe list for the 27k
                     # array (it only covers 450k/EPIC), so that filtering step is
                     # skipped for this platform -- NA is used as a flag for that below.
                     maxprobes_array_type = NA_character_
                   ),
                   stop("array_type must be one of 'EPIC', '450k', or '27k' (got '", array_type, "').")
  )
  
  config
}

# Loads the manifest/annotation packages for one array platform, and
# errors out early with an install hint if they aren't present, rather
# than failing deep inside a minfi call with a less helpful message.
load_array_packages <- function(config) {
  for (pkg in c(config$manifest_package, config$annotation_package)) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for this array type but isn't installed.\n",
           "Install it with: BiocManager::install('", pkg, "')")
    }
    library(pkg, character.only = TRUE)
  }
}

# -------------------------------------------------------------------------
# Shared helper: read a delimited text file with a configurable delimiter
# -------------------------------------------------------------------------
# GEO supplementary files show up as both tab-delimited and comma-delimited
# tables depending on the series, so every function that reads one takes a
# `delim` argument (default "\t") instead of hard-coding read_tsv/read_csv.
# Pass delim = "," for comma-separated files, delim = "\t" (the default)
# for tab-separated ones, or any other single-character delimiter.
# -------------------------------------------------------------------------
# Shared helper: coerce a matrix to numeric, with a clear error if it can't
# -------------------------------------------------------------------------
# If the wrong delimiter is used to read a file, readr often ends up
# cramming an entire row into one character column instead of splitting it
# into the expected numeric columns -- and that failure doesn't surface
# until much later, as a cryptic error deep inside a minfi/matrixStats
# call. This catches it right at the point the data is read, with a
# message that points at the actual cause.
to_numeric_matrix <- function(x, context) {
  if (is.numeric(x)) {
    return(x)
  }
  numeric_x <- suppressWarnings(apply(x, 2, as.numeric))
  dimnames(numeric_x) <- dimnames(x)
  na_fraction <- mean(is.na(numeric_x))
  if (na_fraction > 0.5) {
    # Pull out a handful of actual values that failed to convert, so the
    # real problem (stray text, "1,234" thousands separators, a shifted
    # header row, etc.) is visible instead of just the fact that *some*
    # values failed.
    failed_mask <- is.na(numeric_x) & !is.na(x)
    offending_values <- unique(x[failed_mask])
    sample_values <- utils::head(offending_values, 10)
    
    stop("More than half of the values in ", context, " could not be read as numbers.\n",
         "A sample of the actual values that failed to parse:\n  ",
         paste(sample_values, collapse = ", "), "\n",
         "Common causes: thousands-separator commas in numbers (e.g. \"1,234.5\"), ",
         "an extra header/description row before the real column headers, or a ",
         "genuinely non-numeric column that the suffix pattern matched by mistake -- ",
         "run inspect_supp_file() and check the raw column names/values directly.")
  }
  numeric_x
}

read_table_file <- function(path, delim = "\t", n_max = Inf) {
  readr::read_delim(path, delim = delim, n_max = n_max, show_col_types = FALSE)
}

# -------------------------------------------------------------------------
# Step 1: download supplementary files for a GEO series
# -------------------------------------------------------------------------
download_geo_supplementary <- function(geo_accession, download_dir = NULL) {
  # Default to a fresh temporary directory (via tempfile()) rather than a
  # relative "data" folder in the working directory -- this makes it
  # obvious where files landed, avoids collisions between datasets/runs,
  # and gets cleaned up automatically when the R session ends.
  if (is.null(download_dir)) {
    download_dir <- tempfile(pattern = "geo_downloads_")
  }
  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
  
  # getGEOSuppFiles pulls down every supplementary file listed for the
  # accession (e.g. the Meth/Unmeth matrix, the detection p-value table)
  GEOquery::getGEOSuppFiles(
    GEO = geo_accession,
    baseDir = download_dir,
    makeDirectory = TRUE
  )
  
  # Files land in <download_dir>/<geo_accession>/ -- normalizePath()
  # resolves this to a full absolute path (e.g. "/tmp/geo_downloads_.../
  # GSE140344") regardless of whether download_dir itself was relative
  # or absolute, so it's easy to go find the files directly if needed.
  normalizePath(file.path(download_dir, geo_accession), mustWork = FALSE)
}

# -------------------------------------------------------------------------
# Step 2: peek at a file's structure before writing parsing rules for it
# -------------------------------------------------------------------------
inspect_supp_file <- function(path, n_rows = 5, n_cols = 8, delim = "\t") {
  # Reads just a few rows/columns so you can see the column-naming
  # convention (which varies between GEO series) without loading the
  # whole file -- these are often hundreds of MB -- into memory
  preview <- read_table_file(path, delim = delim, n_max = n_rows)
  preview[, seq_len(min(n_cols, ncol(preview)))]
}

# -------------------------------------------------------------------------
# Step 3: find a downloaded file matching a naming pattern
# -------------------------------------------------------------------------
locate_supp_file <- function(dir_path, file_pattern) {
  matches <- list.files(dir_path, pattern = file_pattern, full.names = TRUE, recursive = TRUE)
  if (length(matches) == 0) {
    stop("No file in ", dir_path, " matched the pattern '", file_pattern, "'.")
  }
  if (length(matches) > 1) {
    warning("Multiple files matched '", file_pattern, "'; using the first one:\n", matches[1])
  }
  matches[1]
}

# -------------------------------------------------------------------------
# Shared helper: read a table whose header row may be one field short
# -------------------------------------------------------------------------
# When a header row has fewer fields than the data rows (because the ID
# column was exported without a column label), readers handle the
# mismatch inconsistently -- some add a placeholder column, some silently
# truncate/misalign the extra data field -- and neither behavior reliably
# preserves probe identity. So this bypasses that automatic handling
# entirely: it reads the header line as plain text, reads the data with
# no header assumptions at all, and aligns the two itself based on their
# actual field counts.
read_table_with_id_column <- function(path, delim = "\t", probe_id_col = "ID_REF") {
  header_line <- readr::read_lines(path, n_max = 1)
  header_tokens <- strsplit(header_line, delim, fixed = TRUE)[[1]]
  
  raw <- readr::read_delim(path, delim = delim, col_names = FALSE, skip = 1, show_col_types = FALSE)
  n_data_cols <- ncol(raw)
  
  if (length(header_tokens) == n_data_cols - 1) {
    # The common case: the ID column was exported with no header of its own
    message("Header row has one fewer field than the data columns (the ID column ",
            "has no header of its own); using column 1 as '", probe_id_col, "'.")
    names(raw) <- c(probe_id_col, header_tokens)
  } else if (length(header_tokens) == n_data_cols) {
    # Header already accounts for every column, including the ID column
    names(raw) <- header_tokens
  } else {
    stop("Header row has ", length(header_tokens), " fields but each data row has ",
         n_data_cols, " fields in ", path, " -- these should differ by at most 1.\n",
         "Open the file directly (e.g. in a text editor) to see what's actually going on.")
  }
  
  raw
}

# -------------------------------------------------------------------------
# Step 4: read the pair of raw intensity columns (Meth/Unmeth or SignalA/SignalB)
# -------------------------------------------------------------------------
# Assumes one row per probe, with a probe ID column plus, for each sample,
# a pair of columns whose names contain that sample's ID along with two
# suffixes -- e.g. "GSM12345.Meth"/"GSM12345.Unmeth", or, for datasets that
# instead label the two color channels directly, "GSM12345.SignalA"/
# "GSM12345.SignalB".
#
# SignalA/SignalB is just a different name for the same pair of raw
# intensities minfi expects as Meth/Unmeth -- but which one (A or B) is
# the methylated channel isn't a fixed rule; it depends on how that
# dataset's processing pipeline assigned the labels (Illumina's own
# GenomeStudio default is SignalA = Unmethylated, SignalB = Methylated,
# but this can vary, so check the series' documentation, or inspect the
# resulting beta value distribution -- it should be roughly bimodal
# between 0 and 1 -- and flip signal_1_is_methylated if it looks inverted).
#
# signal_suffix_1 / signal_suffix_2 identify which two columns to use;
# signal_1_is_methylated says which of the two is the methylated signal.
read_signal_matrix <- function(path,
                               probe_id_col = "ID_REF",
                               signal_suffix_1 = "Meth",
                               signal_suffix_2 = "Unmeth",
                               signal_1_is_methylated = TRUE,
                               delim = "\t") {
  raw <- read_table_with_id_column(path, delim = delim, probe_id_col = probe_id_col)
  
  if (!probe_id_col %in% names(raw)) {
    stop("Column '", probe_id_col, "' (probe_id_col) was not found in ", path, ".\n",
         "This column is used as the probe ID / row name -- if it's missing or NULL, ",
         "everything downstream silently loses track of probe identity instead of erroring here.\n",
         "The file's actual columns include: ", paste(utils::head(names(raw), 10), collapse = ", "), "\n",
         "Set probe_id_col to match the real ID column name for this dataset ",
         "(common alternatives: \"TargetID\", \"ProbeID\", \"Probe_ID\").")
  }
  
  probe_ids <- raw[[probe_id_col]]
  
  cols_1 <- grep(signal_suffix_1, names(raw), value = TRUE)
  cols_2 <- grep(signal_suffix_2, names(raw), value = TRUE)
  
  # Guard against one suffix accidentally matching inside the other
  # (e.g. "Unmeth" contains "Meth"): keep only columns unique to each side
  cols_1 <- setdiff(cols_1, cols_2)
  
  if (length(cols_1) == 0 || length(cols_2) == 0) {
    stop("Could not find columns using suffixes '",
         signal_suffix_1, "' / '", signal_suffix_2,
         "'. Run inspect_supp_file() to check the actual column names.")
  }
  
  matrix_1 <- as.matrix(raw[, cols_1])
  matrix_2 <- as.matrix(raw[, cols_2])
  
  rownames(matrix_1) <- probe_ids
  rownames(matrix_2) <- probe_ids
  
  matrix_1 <- to_numeric_matrix(matrix_1, paste0("the '", signal_suffix_1, "' columns of ", path))
  matrix_2 <- to_numeric_matrix(matrix_2, paste0("the '", signal_suffix_2, "' columns of ", path))
  
  # Strip the suffix so sample names match between the two matrices
  colnames(matrix_1) <- sub(paste0("\\.?", signal_suffix_1, "$"), "", cols_1)
  colnames(matrix_2) <- sub(paste0("\\.?", signal_suffix_2, "$"), "", cols_2)
  
  # Put samples in the same order in both matrices
  common_samples <- intersect(colnames(matrix_1), colnames(matrix_2))
  matrix_1 <- matrix_1[, common_samples, drop = FALSE]
  matrix_2 <- matrix_2[, common_samples, drop = FALSE]
  
  # Assign to Meth/Unmeth based on which signal is the methylated one
  if (signal_1_is_methylated) {
    list(Meth = matrix_1, Unmeth = matrix_2)
  } else {
    list(Meth = matrix_2, Unmeth = matrix_1)
  }
}

# -------------------------------------------------------------------------
# Step 5: read the detection p-value table
# -------------------------------------------------------------------------
read_detection_pvals <- function(path,
                                 probe_id_col = "ID_REF",
                                 detp_suffix = "Detection Pval",
                                 delim = "\t") {
  raw <- read_table_with_id_column(path, delim = delim, probe_id_col = probe_id_col)
  
  if (!probe_id_col %in% names(raw)) {
    stop("Column '", probe_id_col, "' (probe_id_col) was not found in ", path, ".\n",
         "The file's actual columns include: ", paste(utils::head(names(raw), 10), collapse = ", "), "\n",
         "Set probe_id_col to match the real ID column name for this dataset ",
         "(common alternatives: \"TargetID\", \"ProbeID\", \"Probe_ID\").")
  }
  
  probe_ids <- raw[[probe_id_col]]
  detp_cols <- grep(detp_suffix, names(raw), value = TRUE, fixed = TRUE)
  
  if (length(detp_cols) == 0) {
    stop("Could not find detection p-value columns using suffix '", detp_suffix,
         "'. Run inspect_supp_file() to check the actual column names.")
  }
  
  detp_matrix <- as.matrix(raw[, detp_cols])
  rownames(detp_matrix) <- probe_ids
  detp_matrix <- to_numeric_matrix(detp_matrix, paste0("the '", detp_suffix, "' columns of ", path))
  colnames(detp_matrix) <- sub(paste0("\\.?", detp_suffix, "$"), "", detp_cols)
  
  detp_matrix
}

# -------------------------------------------------------------------------
# Step 6: assemble a minfi MethylSet from the intensity matrices
# -------------------------------------------------------------------------
# array_type is one of "EPIC", "450k", or "27k" -- get_array_config() looks
# up the right manifest/annotation package names, and load_array_packages()
# loads just those (so you don't need every platform's packages installed).
build_methyl_set <- function(meth, unmeth, array_type) {
  config <- get_array_config(array_type)
  load_array_packages(config)
  
  minfi::MethylSet(
    Meth = meth,
    Unmeth = unmeth,
    annotation = c(array = config$array_name, annotation = config$annotation_version)
  )
}

# -------------------------------------------------------------------------
# Step 7: quality control
# -------------------------------------------------------------------------
# Runs the QC checks that ARE possible without IDATs/control probes:
#   - overall signal quality (median methylated vs. unmethylated intensity)
#   - per-sample and per-probe detection p-value failure rates
#   - predicted sex, so you can flag any sample/metadata mismatches
run_qc <- function(mset, detp,
                   sample_fail_fraction = 0.1,
                   probe_fail_fraction = 0.1,
                   detection_p_cutoff = 0.01,
                   make_plots = TRUE) {
  
  qc <- minfi::getQC(mset)
  
  # A probe "fails" in a sample if its detection p-value is above the cutoff
  failed <- detp > detection_p_cutoff
  
  # Fraction of probes that failed, per sample -- flags whole bad samples
  sample_fail_rate <- colMeans(failed, na.rm = TRUE)
  failed_samples <- names(sample_fail_rate)[sample_fail_rate > sample_fail_fraction]
  
  # Fraction of samples that failed, per probe -- flags unreliable probes
  probe_fail_rate <- rowMeans(failed, na.rm = TRUE)
  failed_probes <- rownames(failed)[probe_fail_rate > probe_fail_fraction]
  
  # Predicted sex, for comparison against the sample metadata. Wrapped in
  # tryCatch because some platforms -- the 27k array in particular -- don't
  # include enough sex-chromosome probes for this to work reliably.
  gmset <- minfi::mapToGenome(mset)
  predicted_sex <- tryCatch(
    minfi::getSex(gmset),
    error = function(e) {
      warning("Could not predict sex for this array type/dataset: ", conditionMessage(e))
      NULL
    }
  )
  
  if (make_plots) {
    qc_df <- as_tibble(as.data.frame(qc), rownames = "sample")
    
    p_qc <- ggplot(qc_df, aes(mMed, uMed)) +
      geom_point(color = okabe_ito["blue"], size = 2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = okabe_ito["red"]) +
      labs(
        title = "Median methylated vs. unmethylated signal per sample",
        subtitle = "Points far below the dashed line indicate poor overall signal",
        x = "Median methylated intensity (log2)",
        y = "Median unmethylated intensity (log2)"
      ) +
      theme_minimal()
    print(p_qc)
    
    fail_df <- tibble(sample = names(sample_fail_rate), fail_rate = sample_fail_rate)
    p_fail <- ggplot(fail_df, aes(x = reorder(sample, fail_rate), y = fail_rate)) +
      geom_col(fill = okabe_ito["orange"]) +
      geom_hline(yintercept = sample_fail_fraction, linetype = "dashed", color = okabe_ito["red"]) +
      coord_flip() +
      labs(
        title = "Fraction of probes failing detection p-value, per sample",
        subtitle = paste0("Dashed line = ", sample_fail_fraction, " cutoff used to drop samples"),
        x = NULL, y = "Fraction of probes failed"
      ) +
      theme_minimal()
    print(p_fail)
  }
  
  list(
    qc_summary = qc,
    failed_samples = failed_samples,
    failed_probes = failed_probes,
    predicted_sex = predicted_sex
  )
}

# -------------------------------------------------------------------------
# Step 8: filter out unreliable probes
# -------------------------------------------------------------------------
filter_probes <- function(gr_set,
                          failed_probes = character(0),
                          remove_sex_chr = TRUE,
                          remove_snp_probes = TRUE,
                          remove_cross_reactive = TRUE,
                          maxprobes_array_type = NA_character_) {
  
  keep <- rownames(gr_set)
  
  # Drop probes that failed detection p-value QC
  keep <- setdiff(keep, failed_probes)
  
  # Drop probes affected by common SNPs at the CpG or single-base extension site
  if (remove_snp_probes) {
    gr_set <- minfi::dropLociWithSnps(gr_set)
    keep <- intersect(keep, rownames(gr_set))
  }
  
  # Drop probes known to cross-react with multiple genomic locations.
  # maxprobes only publishes a blacklist for 450k and EPIC, so this step
  # is skipped (with a warning) for platforms it doesn't cover, like 27k.
  if (remove_cross_reactive) {
    if (is.na(maxprobes_array_type)) {
      warning("Cross-reactive probe filtering isn't available for this array ",
              "type via maxprobes; skipping this filtering step.")
    } else {
      cross_reactive <- maxprobes::xreactive_probes(array_type = maxprobes_array_type)
      keep <- setdiff(keep, cross_reactive)
    }
  }
  
  # Drop sex-chromosome probes (only relevant if sex isn't a variable of interest)
  if (remove_sex_chr) {
    probe_annotation <- minfi::getAnnotation(gr_set)
    autosomal_probes <- rownames(probe_annotation)[!probe_annotation$chr %in% c("chrX", "chrY")]
    keep <- intersect(keep, autosomal_probes)
  }
  
  gr_set[rownames(gr_set) %in% keep, ]
}

# -------------------------------------------------------------------------
# Step 9: normalize
# -------------------------------------------------------------------------
# preprocessFunnorm/preprocessNoob both need control-probe intensities from
# raw IDATs, which aren't available here -- so preprocessQuantile is the
# normalization method for this processed-matrix pipeline.
normalize_methyl_set <- function(mset) {
  # preprocessQuantile warns that it's "only been tested with preprocessRaw"
  # whenever it can't verify the MethylSet's processing history -- which is
  # always true here, since we build the MethylSet manually from a CSV
  # rather than from preprocessRaw() on IDATs. The underlying data is the
  # same shape (raw, un-normalized Meth/Unmeth intensities) that
  # preprocessRaw() would have produced, so the warning doesn't indicate an
  # actual problem in this pipeline -- it's suppressed here so it doesn't
  # repeat on every dataset run through this function.
  suppressWarnings(minfi::preprocessQuantile(mset))
}

# -------------------------------------------------------------------------
# Step 10: annotate each probe (gene, CpG island context, genomic feature)
# -------------------------------------------------------------------------
# Deliberately kept at the probe level rather than collapsed to one value
# per gene -- a promoter CpG and a gene-body CpG can behave very
# differently, so averaging them together would hide that distinction.
# CpG island context (Island/Shore/Shelf/OpenSea) is included because it's
# one of the strongest predictors of whether a methylation change is
# likely to be functionally meaningful.
annotate_probes <- function(gr_set) {
  beta_values <- minfi::getBeta(gr_set)
  
  probe_annotation <- minfi::getAnnotation(gr_set) %>%
    as.data.frame() %>%
    as_tibble(rownames = "probe_id") %>%
    select(probe_id, chr, pos, UCSC_RefGene_Name, UCSC_RefGene_Group, Relation_to_Island)
  
  beta_tbl <- as_tibble(beta_values, rownames = "probe_id")
  
  probe_annotation %>%
    left_join(beta_tbl, by = "probe_id")
}

# -------------------------------------------------------------------------
# Master wrapper: runs every step above in order
# -------------------------------------------------------------------------
run_methylation_pipeline <- function(geo_accession,
                                     array_type,
                                     meth_unmeth_filename_pattern,
                                     detection_pval_filename_pattern,
                                     download_dir = NULL,
                                     probe_id_col = "ID_REF",
                                     signal_suffix_1 = "Meth",
                                     signal_suffix_2 = "Unmeth",
                                     signal_1_is_methylated = TRUE,
                                     meth_unmeth_delim = "\t",
                                     detp_suffix = "Detection Pval",
                                     detp_delim = "\t",
                                     detection_p_cutoff = 0.01,
                                     sample_fail_fraction = 0.1,
                                     probe_fail_fraction = 0.1,
                                     remove_sex_chr = TRUE,
                                     remove_snp_probes = TRUE,
                                     remove_cross_reactive = TRUE,
                                     make_plots = TRUE) {
  
  message("Step 1/6: downloading supplementary files for ", geo_accession)
  series_dir <- download_geo_supplementary(geo_accession, download_dir)
  message("  files stored at: ", series_dir)
  
  message("Step 2/6: reading signal intensities and detection p-values")
  meth_unmeth_path <- locate_supp_file(series_dir, meth_unmeth_filename_pattern)
  detp_path <- locate_supp_file(series_dir, detection_pval_filename_pattern)
  
  intensities <- read_signal_matrix(meth_unmeth_path, probe_id_col,
                                    signal_suffix_1, signal_suffix_2, signal_1_is_methylated,
                                    delim = meth_unmeth_delim)
  detp <- read_detection_pvals(detp_path, probe_id_col, detp_suffix, delim = detp_delim)
  
  message("Step 3/6: building MethylSet and running QC")
  mset <- build_methyl_set(intensities$Meth, intensities$Unmeth, array_type)
  qc_results <- run_qc(mset, detp, sample_fail_fraction, probe_fail_fraction, detection_p_cutoff, make_plots)
  
  if (length(qc_results$failed_samples) > 0) {
    message("Dropping ", length(qc_results$failed_samples), " sample(s) that failed QC: ",
            paste(qc_results$failed_samples, collapse = ", "))
    keep_samples <- setdiff(colnames(mset), qc_results$failed_samples)
    mset <- mset[, keep_samples]
  }
  
  message("Step 4/6: normalizing (quantile normalization)")
  gr_set <- normalize_methyl_set(mset)
  
  message("Step 5/6: filtering unreliable probes")
  array_config <- get_array_config(array_type)
  gr_set <- filter_probes(gr_set, qc_results$failed_probes, remove_sex_chr, remove_snp_probes,
                          remove_cross_reactive, array_config$maxprobes_array_type)
  
  message("Step 6/6: annotating probes (gene, CpG island context, genomic feature)")
  annotated_data <- annotate_probes(gr_set)
  
  list(
    genomic_ratio_set = gr_set,
    annotated_data = annotated_data,
    qc_results = qc_results
  )
}

# -------------------------------------------------------------------------
# Example usage
# -------------------------------------------------------------------------
# First, take a look at the actual column names, since they vary by series:
#   series_dir <- download_geo_supplementary("GSE140344")
#   inspect_supp_file(locate_supp_file(series_dir, "Meth_UnMethSamples"))
#   inspect_supp_file(locate_supp_file(series_dir, "Beta_DetectionPvals"))
#
# Example 1 -- EPIC array, standard Meth/Unmeth columns (GSE140344):
#
# result <- run_methylation_pipeline(
#   geo_accession = "GSE140344",
#   array_type = "EPIC",
#   meth_unmeth_filename_pattern = "Meth_UnMethSamples",
#   detection_pval_filename_pattern = "Beta_DetectionPvals"
# )
#
# Example 2 -- 450k array, SignalA/SignalB columns instead of Meth/Unmeth:
# Confirm which channel is methylated for this series before running --
# see the note above read_signal_matrix() -- then set signal_1_is_methylated
# accordingly (FALSE here means SignalB is the methylated channel).
#
result <- run_methylation_pipeline(
  geo_accession = "GSE50586",
  array_type = "450k",
  meth_unmeth_filename_pattern = "signal_intensities",
  detection_pval_filename_pattern = "detection_pvalues",
  signal_suffix_1 = "SignalA",
  signal_suffix_2 = "SignalB",
  signal_1_is_methylated = FALSE
)
#
# Example 3 -- 27k array:
#
# result <- run_methylation_pipeline(
#   geo_accession = "GSEYYYYYY",
#   array_type = "27k",
#   meth_unmeth_filename_pattern = "Meth_UnMethSamples",
#   detection_pval_filename_pattern = "Detection_Pvals"
# )
#
# Example 4 -- a series whose files are comma-separated (.csv) rather than
# tab-separated -- override delim on whichever file(s) need it:
#
# result <- run_methylation_pipeline(
#   geo_accession = "GSEZZZZZZ",
#   array_type = "EPIC",
#   meth_unmeth_filename_pattern = "signal_intensities",
#   detection_pval_filename_pattern = "detection_pvalues",
#   meth_unmeth_delim = ",",
#   detp_delim = ","
# )