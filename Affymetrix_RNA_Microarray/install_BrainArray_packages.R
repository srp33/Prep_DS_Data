# http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/entrezg.asp
# Example: http://brainarray.mbni.med.umich.edu/customcdf/25.0.0/entrezg.download/hgu133plus2hsentrezgprobe_25.0.0.tar.gz

library(tidyverse)
library(readxl)

#annotationSource = "entrezg"
annotationSource = "ensg"
version = "25.0.0"

# R's download.file()/libcurl often hangs on BrainArray after the progress
# bar fills: the old HTTP server sends the bytes but does not close the
# keep-alive connection. curl with Connection: close, stall detection, and
# retries is more reliable for the large probe packages.
download_brainarray_file = function(url, destfile) {
  if (file.exists(destfile)) {
    unlink(destfile)
  }

  status = system2(
    "curl",
    args = c(
      "--fail",
      "--location",
      "--show-error",
      "--retry", "20",
      "--retry-delay", "5",
      "--retry-all-errors",
      "--connect-timeout", "30",
      "--max-time", "1800",
      "--speed-time", "30",
      "--speed-limit", "1000",
      "--http1.0",
      "--header=Connection:close",
      "-C-",
      "-o", destfile,
      url
    )
  )

  if (!identical(status, 0L) || !file.exists(destfile) || file.info(destfile)$size == 0) {
    stop("Failed to download ", url)
  }
}

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
  packageUrl = paste0(
    "http://brainarray.mbni.med.umich.edu/customcdf/", version,
    "/", annotationSource, ".download/",
    platform, annotationSource, "probe_", version, ".tar.gz"
  )

  tempPackageFileName = basename(packageUrl)
  tempPackageFilePath = paste0(tempdir(), "/", tempPackageFileName)

  message("Attempting to download ", tempPackageFileName, " from BrainArray with URL: ", packageUrl)
  download_brainarray_file(packageUrl, tempPackageFilePath)
  message("Downloaded ", tempPackageFileName, " from BrainArray")
  install.packages(tempPackageFilePath, repos = NULL, type = "source")
}