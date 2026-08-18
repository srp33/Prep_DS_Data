library(tidyverse)
library(readxl)
library(GEOquery)

options(timeout = 10000) # Allows very large tar files to still download.

datasets = read_xlsx("/project/datasets.xlsx")

geo_ids = datasets %>%
  filter(repository == "GEO") %>%
  pull(series_id)

for (geo_id in geo_ids) {
  celDirPath = "/project/GEO"
  if (!dir.exists(celDirPath)){
    dir.create(celDirPath, recursive = TRUE)
  }

  geoDirPath = paste0(celDirPath, "/", geo_id)
  tarFilePath = paste0(geoDirPath, "/", geo_id, "_RAW.tar")

  # Skip series that were already downloaded and unpacked (failed runs may
  # leave an empty directory or only the RAW.tar).
  existingCelFiles = list.files(
    geoDirPath,
    pattern = "\\.CEL(\\.gz)?$",
    ignore.case = TRUE
  )
  if (length(existingCelFiles) > 0) {
    message("Skipping ", geo_id, ": CEL files already present")
    next
  }

  message("Downloading CEL files for ", geo_id, " (this can take a while for large series)")
  getGEOSuppFiles(geo_id, makeDirectory = TRUE, baseDir = celDirPath)
  untar(tarFilePath, exdir = geoDirPath)
  file.remove(tarFilePath)

  # GEO supplementary files often include CHP or other non-CEL data; keep CEL only.
  allFiles = list.files(geoDirPath, full.names = TRUE)
  allFiles = allFiles[!dir.exists(allFiles)]
  isCelFile = grepl("\\.CEL(\\.gz)?$", basename(allFiles), ignore.case = TRUE)
  nonCelFiles = allFiles[!isCelFile]
  if (length(nonCelFiles) > 0) {
    message(
      "Removing ", length(nonCelFiles), " non-CEL file(s) from ", geo_id, ": ",
      paste(basename(nonCelFiles), collapse = ", ")
    )
    file.remove(nonCelFiles)
  }

  celFiles = allFiles[isCelFile]
  if (length(celFiles) == 0) {
    stop("No CEL files found for ", geo_id, " after removing non-CEL files")
  }

  # Rename CEL files to GSM######.CEL.gz (or .CEL if uncompressed).
  for (celFile in celFiles) {
    gsmMatch = str_match(basename(celFile), regex("^GSM(\\d+)", ignore_case = TRUE))
    if (is.na(gsmMatch[1, 1])) {
      next
    }

    celExtension = if (grepl("\\.gz$", celFile, ignore.case = TRUE)) ".CEL.gz" else ".CEL"
    newCelPath = file.path(geoDirPath, paste0("GSM", gsmMatch[1, 2], celExtension))
    if (!identical(celFile, newCelPath)) {
      file.rename(celFile, newCelPath)
    }
  }
  message("Files for ", geo_id, " downloaded and stored in ", geoDirPath)
}