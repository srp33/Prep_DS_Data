library(tidyverse)
library(readxl)
library(GEOquery)

options(timeout = 10000) # Allows very large tar files to still download.

datasets = read_xlsx("/project/datasets.xlsx")

geo_ids = datasets %>%
  filter(repository == "GEO") %>%
  pull(series_id)

for (geo_id in geo_ids) {
  celDirPath = "/project/GEO/"
  if (!dir.exists(celDirPath)){
    dir.create(celDirPath, recursive = TRUE)
  }

  tarFilePath = paste0(celDirPath, "/", geo_id, "/", geo_id, "_RAW.tar")

  getGEOSuppFiles(geo_id, makeDirectory = TRUE, baseDir = celDirPath)
  untar(tarFilePath, exdir = paste0(celDirPath, "/", geo_id))
  file.remove(tarFilePath)

  # Rename CEL files to GSM######.CEL.gz (case-insensitive match; always uppercase GSM/CEL).
  geoDirPath = paste0(celDirPath, "/", geo_id)
  celFiles = list.files(geoDirPath, pattern = "\\.CEL\\.gz$", ignore.case = TRUE, full.names = TRUE)

  for (celFile in celFiles) {
    gsmMatch = str_match(basename(celFile), regex("^GSM(\\d+)", ignore_case = TRUE))
    if (is.na(gsmMatch[1, 1])) {
      next
    }

    newCelPath = file.path(geoDirPath, paste0("GSM", gsmMatch[1, 2], ".CEL.gz"))
    if (!identical(celFile, newCelPath)) {
      file.rename(celFile, newCelPath)
    }
  }
}