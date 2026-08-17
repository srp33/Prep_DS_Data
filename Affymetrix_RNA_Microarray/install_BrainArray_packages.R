# BrainArray custom CDF probe packages (Ensembl gene / version 25).
# Source: https://zenodo.org/records/17808641 (archived copy of the Michigan files).
# Example package file: mogene10stmmensgprobe_25.0.0.tar.gz

library(tidyverse)
library(readxl)

annotationSource = "ensg"
version = "25.0.0"

# Split cells that list multiple BrainArray platforms (semicolon-separated),
# then collect the unique platform values.
brainarrayPlatforms = read_xlsx("/project/datasets.xlsx") %>%
  filter(!is.na(brainarray), brainarray != "") %>%
  separate_rows(brainarray, sep = ";") %>%
  mutate(brainarray = str_trim(brainarray)) %>%
  filter(brainarray != "") %>%
  distinct(brainarray) %>%
  pull(brainarray)

# Package names follow BrainArray's {platform}{source}probe convention.
brainarrayPackageNames = paste0(brainarrayPlatforms, annotationSource, "probe")
missingPackageNames = setdiff(brainarrayPackageNames, rownames(installed.packages()))

if (length(missingPackageNames) == 0) {
  message("All required BrainArray packages are already installed")
} else {
  message(
    "Need to install ", length(missingPackageNames),
    " BrainArray package(s): ", paste(missingPackageNames, collapse = ", ")
  )

  # Isolated temp dir so the 700 MB zip and extracted tarballs can be removed together.
  tempDir = tempfile("brainarray_")
  dir.create(tempDir)
  zipFilePath = file.path(tempDir, "BrainArray.zip")

  # Zenodo is used instead of the Michigan HTTP server, which stalls on large files.
  # 30-minute timeout: the archive is ~700 MB.
  options(timeout = max(1800, getOption("timeout")))
  zenodoZipUrl = "https://zenodo.org/records/17808641/files/BrainArray.zip?download=1"
  message("Downloading BrainArray.zip from Zenodo")
  download.file(zenodoZipUrl, zipFilePath, mode = "wb")
  unzip(zipFilePath, exdir = tempDir)

  # The zip contains BrainArray/*.tar.gz (plus macOS sidecar files under __MACOSX/).
  for (packageName in missingPackageNames) {
    packageFileName = paste0(packageName, "_", version, ".tar.gz")
    packageFilePath = file.path(tempDir, "BrainArray", packageFileName)

    if (!file.exists(packageFilePath)) {
      stop("Package file not found in BrainArray.zip: ", packageFileName)
    }

    message("Installing ", packageName, " from ", packageFileName)
    install.packages(packageFilePath, repos = NULL, type = "source")
  }

  # Remove the zip and extracted files only after every missing package installed.
  unlink(tempDir, recursive = TRUE)
  message("Removed temporary BrainArray download files")
}
