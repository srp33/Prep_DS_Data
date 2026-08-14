library(GEOquery)
library(tidyverse)
# library(janitor)

get_metadata = function(geoID) {
  metadata = getGEO(geoID)[[1]] 
  
  metadata = as_tibble(pData(metadata))
  # metadata = clean_names(metadata)
  # metadata = fix_bespoke_issues(geoID, metadata) # debug
  return (metadata)
}

# fixes specific issues in the data for certain datasets
fix_bespoke_issues = function(geoID, metadata){
  if(geoID == "GSE1789"){
    metadata = mutate(metadata, age_ch1 = ifelse(is.na(age_ch1), substr(characteristics_ch1_1, start = 5, stop = nchar(characteristics_ch1_1) ), age_ch1)) %>%
      mutate(tissue_ch1 = ifelse(is.na(tissue_ch1), tissue_ch1[1], tissue_ch1)) %>%
      mutate(label_ch1 = ifelse(label_ch1 == "", label_ch1[2], label_ch1)) %>%
      mutate(extract_protocol_ch1 = ifelse(extract_protocol_ch1 == "Total RNA was extract by Trizol", extract_protocol_ch1[1], extract_protocol_ch1)) %>%
      mutate(description_1 = ifelse(startsWith(description_1, "Image data were analyzed using the Affymetrix  expression"), description_1[3], description_1)) %>%
      mutate(tissue_ch1 = ifelse(startsWith(tissue_ch1,"H"), tissue_ch1[1], tissue_ch1)) %>%
      mutate(source_name_ch1 = ifelse(endsWith(source_name_ch1, "week of gestation"), paste0(substr(source_name_ch1, start = 1, stop = 22), "s of gestation"), source_name_ch1))
  }else{
    return(metadata)
  }
}


# returns version of metadata having removed columns with keywords (only called on dataset metadata)
drop_cols = function(metadata) {
  keywords = c("contact", "library", "processing", "relation",
               "platform", "instrument", "file", "date", "row",
               "status", "characteristics", "time", "channel", "taxid", 
               "email", "web.link", "sample.id", "geo.accession") # TODO: figure out what to do for characteristics column.
                                                                        # in some cases, characteristics should be kept. check to see if they got extracted. Add arg?
  cols_to_drop = names(metadata)[grepl(paste(keywords, collapse = "|"), names(metadata))]
  
  # drop columns that are not needed based on above criteria
  metadata_filtered = metadata %>% 
  dplyr::select(-cols_to_drop)
  
  return(metadata_filtered)
}


# goes through each row in the attributes tibble and searches the input for the attributes desired, creating a column for each in the output
standardize_tibble = function(geo_id, input_tbl, attr_tbl) {
  #TODO: refactor so it doesn't count matches (unneeded since we're separating same and diff)
  n = nrow(input_tbl)
  geo_col = rep(as.character(geo_id), n)
  
  # For each column, count how many regex pattern matches
  score_column = function(col, value_dict) {
    sum(map_int(value_dict, function(regex_pattern) {
      any(str_detect(col, regex_pattern), na.rm = TRUE)
    }))
  }
  
  # Iterate over rows of attr_tbl and return one tibble column per attribute
  cols_list = map(seq_len(nrow(attr_tbl)), function(i) { # 1:nrow(attr_tbl)
    attr = attr_tbl[i, ]
    attr_name = as.character(attr$attr_name)
    match_type = as.character(attr$match_type)
    
    if (match_type == "column") {
      pattern = as.character(attr$col_regex)
      col_match = names(input_tbl)[str_detect(names(input_tbl), pattern)]
      if (length(col_match) > 0) {
        # return the first column it matches the title of
        return(tibble(!!attr_name := input_tbl[[col_match[1]]]))
      } else {
        # error handling: if no column title matches, return column of all NA
        return(tibble(!!attr_name := rep(NA_character_, n)))
      }
    }
    
    else if (match_type == "value") {
      value_dict = attr$value_dict[[1]]
      
      # Keep only character or factor columns
      kept_names = names(input_tbl)[map_lgl(names(input_tbl), function(.x){is.character(input_tbl[[.x]]) || is.factor(input_tbl[[.x]])})]
      if (length(kept_names) == 0) {
        return(tibble(!!attr_name := rep(NA_character_, n)))
      }
      
      # create named list with each column stored as a vector
      kept_cols = map(kept_names, function(.x){input_tbl[[.x]]})
      names(kept_cols) = kept_names
      
      # Score each kept column
      scores = map_int(kept_cols, function(.x){score_column(.x, value_dict)})

      # If no column matched any pattern, return NA column
      if (all(scores == 0)) {
        return(tibble(!!attr_name := rep(NA_character_, n), "Original_Value" := rep(NA_character_, n)))
      }
      
      # Choose best column (first with max score)
      best_index = which.max(scores)
      best_col_name = kept_names[best_index]
      best_col_values = kept_cols[[best_index]]
      
      # Build standardized vector (first match wins)
      standardized = rep(NA_character_, length(best_col_values))
      for (key in names(value_dict)) {
        pattern = value_dict[[key]]
        matches = str_detect(best_col_values, pattern)
        # Assign key only where standardized is still NA (first match wins)
        standardized[is.na(standardized) & matches] = key
      }
      
      return(tibble(!!attr_name := standardized, "Original_Value" := best_col_values))
    }
    else{
      stop("Unknown match type: ", match_type)
    }
  })
  
  # Combine all columns into one tibble
  result = bind_cols(tibble(Dataset_ID = geo_col), bind_cols(cols_list))
  rotated_result = pivot_longer(result, !c("Dataset_ID", "ID", "Original_Value"), names_to = "Attribute", values_to = "Value")
  return(rotated_result)
}


get_gse_metadata = function(gse_id) {
  # Get GEO dataset (GSEMatrix = FALSE returns the main metadata structure)
  gse = getGEO(gse_id, GSEMatrix = FALSE)
  
  # Extract metadata list
  meta = Meta(gse)
  
  # Collapse any vector elements into single strings
  meta_collapsed = map(meta, function(x) {
    if (length(x) > 1) paste(x, collapse = "; ") else x
  })
  
  # Create a single-row tibble
  tibble::as_tibble_row(meta_collapsed)
}


append_to_file = function(data, file){
  if (!file.exists(file)) {
    write_tsv(data, file)
  } else {
    write_tsv(data, file, append = TRUE, col_names = FALSE)
  }
}

geoDirs = list.dirs("/project/GEO", full.names = TRUE, recursive = FALSE)

for (geoDir in geoDirs) {
  geoID = basename(geoDir)

  metadataFilePath = paste0(geoDir, "/Metadata.tsv")

  if (!file.exists(metadataFilePath)){
    metadata = get_metadata(geoID)
    write_tsv(metadata, metadataFilePath)
  }
}
stop("got here")

#--------------process_metadata-------------

file_location = "Data/Metadata/"

if (!dir.exists(file_location)){
  dir.create(file_location, recursive = TRUE)
}


# loop through all series IDs
for (geo_id in pull(Datasets, Name)) {
  
  # read metadata into a variable, drop unneeded columns, split into same and diff
  metadata = get_metadata(geo_id)
  diff_metadata = dplyr::select(metadata, where(~n_distinct(.) > 1))
  same_metadata = dplyr::select(metadata, where(~n_distinct(.) == 1))
  
  
  # get sample metadata
  result = standardize_tibble(geo_id, diff_metadata, attr_tbl) %>%
    relocate(Original_Value, .after = Value)

  append_to_file(result, "Data/Metadata/SampleMetadata.tsv")
  

  #get dataset metadata
  dataset_result = same_metadata[1,] %>%
    drop_cols() %>%
    mutate(Dataset_ID = geo_id) %>%
    rename_with(~"data_type", any_of("type"))

  # add attributes from website
  website_metadata = get_gse_metadata(geo_id)
  dataset_result = bind_cols(dataset_result, website_metadata) %>%
    drop_cols() %>%
    rename_with(~"platform_type", any_of("type")) %>%
    rename_with(~"repository", any_of("name"))
  
  # rotate and remove _ch1 from attributes
  rotated_result = pivot_longer(dataset_result, !"Dataset_ID", names_to = "Attribute", values_to = "Value") %>%
    mutate(Attribute = ifelse(endsWith(Attribute, "_ch1"), str_sub(Attribute, 1, -5), Attribute))
  
  append_to_file(rotated_result, "Data/Metadata/DatasetMetadata.tsv")
}