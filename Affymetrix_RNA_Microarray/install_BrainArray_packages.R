# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/entrezg.asp
# Example: http://brainarray.mbni.med.umich.edu/customcdf/25.0.0/entrezg.download/hgu133plus2hsentrezgprobe_25.0.0.tar.gz

library(tidyverse)
library(readxl)

annotationSource = "entrezg"
version = "25.0.0"

# Split cells that list multiple BrainArray platforms (semicolon-separated),
# then collect the unique platform values.
brainarrayPlatforms = read_xlsx("/project/datasets.xlsx") %>%
  filter(!is.na(brainarray), brainarray != "") %>%
  # filter(brainarray != "hta20hs") %>% # This is a temporary exclusion until we can download this package from BrainArray.
  separate_rows(brainarray, sep = ";") %>%
  mutate(brainarray = str_trim(brainarray)) %>%
  filter(brainarray != "") %>%
  distinct(brainarray) %>%
  pull(brainarray)

brainarrayPlatforms = c("hta20hs")

for (platform in brainarrayPlatforms) {
  packageUrl = paste0(
    "http://brainarray.mbni.med.umich.edu/customcdf/", version,
    "/", annotationSource, ".download/",
    platform, annotationSource, "probe_", version, ".tar.gz"
  )

  tempPackageFileName = basename(packageUrl)
  tempPackageFilePath = paste0(tempdir(), "/", tempPackageFileName)

print(packageUrl)
  download.file(packageUrl, tempPackageFilePath)

  install.packages(tempPackageFilePath, repos = NULL, type = "source")
}