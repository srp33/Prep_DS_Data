# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/entrezg.asp
# Example: http://brainarray.mbni.med.umich.edu/customcdf/25.0.0/entrezg.download/hgu133plus2hsentrezgprobe_25.0.0.tar.gz

library(tidyverse)
library(readxl)

#annotationSource = "entrezg"
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

for (platform in brainarrayPlatforms) {
          # "http://brainarray.mbni.med.umich.edu/customcdf/", version,
  packageUrl = paste0(
    "http://mbni.org/customcdf/", version,
    "/", annotationSource, ".download/",
    platform, annotationSource, "probe_", version, ".tar.gz"
  )
  # https://web.archive.org/web/20240710114225/http://mbni.org/customcdf/25.0.0/ensg.download/mogene10stmmensgprobe_25.0.0.tar.gz
  # https://web.archive.org/web/20240710102334/http://mbni.org/customcdf/25.0.0/ensg.download/mogene10stmmensgprobe_25.0.0.tar.gz

  tempPackageFileName = basename(packageUrl)
  tempPackageFilePath = paste0(tempdir(), "/", tempPackageFileName)

  message("Attempting to download ", tempPackageFileName, " from BrainArray with URL: ", packageUrl)
  download.file(packageUrl, tempPackageFilePath)
  message("Downloaded ", tempPackageFileName, " from BrainArray")
  install.packages(tempPackageFilePath, repos = NULL, type = "source")
}