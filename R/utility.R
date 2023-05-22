#' @importFrom dplyr %>%
NULL

########## MLST Cleanup and Verification ##########

#' Determine FastMLST species 'codename' R implementation of
#' best_guess_codename(text) from FastMLST v0.0.15. Used to match the "scheme"
#' listing of FastMLST CSV to PubMLST schema
#'
#' See: https://github.com/EnzoAndree/FastMLST
#'
#' @param text species name to parse
fastmlst_scheme <- function(text){
  text <- stringr::str_trim(text)
  text <- stringr::str_to_lower(text)
  text <- stringr::str_replace_all(text, c("\\." = "","candidatus " = "","/" = "","\\(" = "","\\)" = ""))

  if(length(unlist(stringr::str_split(text, " "))) == 1){
    #"Genus" scheme
    if(grepl("#", text)){
      # Multiple schemes
      split_num <- unlist(stringr::str_split(text, "#"))
      genus <- split_num[1]
      scheme_number <- split_num[length(split_num)]
      scheme <- paste0(genus, "#", scheme_number)
    } else{
      scheme <- text
    }
  }else if(grepl("spp", text)){
    #"Genus spp." scheme
    if(grepl("#", text)){
      split <- unlist(stringr::str_split(text, " "))
      split_num <- unlist(stringr::str_split(text, "#"))
      genus <- split[1]
      scheme_number <- split_num[length(split_num)]
      scheme <- paste0(genus, "#", scheme_number)
    } else{
      split <- unlist(stringr::str_split(text, " "))
      genus <- split[1]
      scheme <- genus
    }
  }else if(length(unlist(stringr::str_split(text, " "))) == 2){
    #"Genus species" scheme
    if(grepl("#", text)){
      split <- unlist(stringr::str_split(text, " "))
      split_num <- unlist(stringr::str_split(text, "#"))
      genus <- str_sub(split[1],1,1)
      species <- unlist(stringr::str_split(split[2], "#"))[1]
      scheme_number <- split_num[length(split_num)]
      scheme <- paste0(genus, species, "#", scheme_number)
    } else{
      split <- unlist(stringr::str_split(text, " "))
      genus <- str_sub(split[1],1,1)
      species <- unlist(stringr::str_split(split[2], "#"))[1]
      scheme <- paste0(genus, species)
    }
  }else if(length(unlist(str_split(text, " "))) > 2){
    #"Genus species extra" scheme
    if(grepl("#", text)){
      split <- unlist(stringr::str_split(text, " "))
      split_num <- unlist(stringr::str_split(text, "#"))
      genus <- stringr::str_sub(split[1],1,1)
      species <- unlist(stringr::str_split(split[2], "#"))[1]
      extra <- paste0(split[2:length(split) - 1],
                      unlist(stringr::str_split(split[length(split)], "#"))[1], collapse = "_")
      scheme_number <- split_num[length(split_num)]
      scheme <- paste0(genus, species,"_", extra, "#", scheme_number)
    } else{
      split <- unlist(stringr::str_split(text, " "))
      genus <- stringr::str_sub(split[1],1,1)
      species <- unlist(stringr::str_split(split[2], "#"))[1]
      extra <- paste0(split[2:length(split) - 1], collapse = "_")
      scheme <- paste0(genus, species, "_", extra)
    }
  }
  return(scheme)
}


#' Download list of all PubMLST MLST schema
#'
#' @param url (string) URL of PubMLST "dbases.xml" (this should only be needed
#'   if the site changes unexpectedly)

pubmlst_listing <- function(url = "https://pubmlst.org/static/data/dbases.xml"){

  xml_data <- xml2::read_xml(url)
  species_list <- xml_data %>%
    xml2::xml_find_all("/data/species")

  species <- vector(mode = "character", length = length(species_list))
  profile <- vector(mode = "character", length = length(species_list))
  url <- vector(mode = "character", length = length(species_list))

  for(n in 1:length(species_list)){
    species[n] <- species_list[[n]] %>%
      xml2::xml_find_first("./text()") %>%
      xml2::xml_text(trim = TRUE)

    profile[n] <- species_list[[n]] %>%
      xml2::xml_find_all("./mlst/database/loci/locus/text()") %>%
      xml2::xml_text(trim = TRUE) %>%
      stringi::stri_remove_empty_na() %>%
      paste0(collapse = ";")

    url[n] <- species_list[[n]] %>%
      xml2::xml_find_all("./mlst/database/profiles/url") %>%
      xml2::xml_text(trim = TRUE) %>%
      stringi::stri_remove_empty_na() %>%
      paste0(collapse = ";")
  }

  return(tibble(species = species, profile = profile, url = url))
}

#' Determine MLST based on allele calls and profiles
#'
#' @param .calls a datatable or tibble of FastMLST allele calls
#' @param .profiles a datatable or tibble of PubMLST MLST profiles
parse_mlst_profiles <- function(.calls, .profiles){

  call_names <- .calls %>%
    dplyr::select(-Genome, -Scheme, -ST, - dplyr::any_of(c("clonal_complex"))) %>%
    names()
  profile_names <- .profiles %>%
    dplyr::select(- dplyr::any_of(c("ST", "clonal_complex", "species"))) %>%
    names()

  if(!setequal(call_names, profile_names))
    warning("Allele names do not match. Ensure correct scheme is being used.")

  .complete <- .calls %>%
    dplyr::select(Genome, all_of(profile_names)) %>%
    separate_rows_sequential(-Genome, sep = "\\|") %>%
    dplyr::distinct() %>%
    dplyr::left_join(mutate(.profiles, match = TRUE), by = profile_names) %>%
    tidyr::unite(col = combined, dplyr::all_of(profile_names), sep = ",", na.rm = TRUE) %>%
    dplyr::mutate(error = dplyr::case_when(
      match ~ NA_character_,
      grepl("~", combined) ~ "new_allele",
      stringr::str_count(combined, ",") < length(profile_names) - 1 ~ "missing_allele",
      TRUE ~ "new_st"
    )) %>%
    dplyr::group_by(Genome) %>%
    dplyr::mutate(ST = dplyr::na_if(paste0(stringr::str_sort(stats::na.omit(unique(ST)), numeric = TRUE), collapse = ";"), "")) %>%
    dplyr::mutate(errors = dplyr::na_if(paste0(stats::na.omit(unique(error)), collapse = ";"), "")) %>%
    dplyr::mutate(row = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::select(Genome, ST, errors) %>%
    dplyr::mutate(assembly = as.character(stringr::str_match(Genome,"GCA_\\d+\\.\\d+"))) %>%
    dplyr::distinct()

  return(.complete)
}

#' Download and save MLST schema
#'
#' @param path path to save schema to
#' @param url path to download schema from
download_mlst_profile <- function(path, url){
  ensure_directory(dirname(path))
  raw_data <- RCurl::getURL(url = url, ssl.verifypeer = FALSE)
  write(raw_data, path)

  message(paste("Saved", tools::file_path_sans_ext(basename(path)), "MLST Schema to", path))
}

#' Automatically re-parse FastMLST calls using profiles downloaded from PubMLST
#'
#' Download correct schema from PubMLST and use to re-parse allele calls from
#' FastMLST. Optionally save the re-parsed file for future use.
#'
#' This is primarily useful to handle situations where multiple alleles present
#' in a single strain caused failed MLST typing in the original file.
#'
#' @param path location of FastMLST csv output file
#' @param save should the re-parsed file be saved?
reparse_mlst <- function(path, save = TRUE){
  raw_mlst <- utils::read.csv(path)
  mlst <- process_mlst(raw_mlst)

  pubmlst <- pubmlst_listing() %>%
    dplyr::mutate(Scheme = purrr::map_chr(species, fastmlst_scheme))

  current_schema <- raw_mlst %>%
    dplyr::left_join(pubmlst, by = "Scheme") %>%
    dplyr::distinct(Scheme, url) %>%
    tidyr::drop_na()

  Scheme <- current_schema %>%
    dplyr::pull(Scheme)

  url <- current_schema %>%
    dplyr::pull(url)

  profile_path <- paste0(dirname(path), "/", stringr::str_replace_all(Scheme, "#", ""),".tsv")
  reparsed_path <- paste0(tools::file_path_sans_ext(path), "_reparsed",".csv")

  if(!file.exists(profile_path))
    download_mlst_profile(profile_path, url)

  profiles <- utils::read.delim(profile_path, na.strings = "", colClasses = "character")

  fixed <- raw_mlst %>%
    parse_mlst_profiles(profiles) %>%
    dplyr::mutate(ST = dplyr::case_when(is.na(ST) ~ NA_character_, TRUE ~ paste0("ST", ST)))

  if(save){
    utils::write.csv(fixed, file = reparsed_path, row.names = FALSE)
    message(paste0("Saved re-parsed ",
                   tools::file_path_sans_ext(path)," to `",
                   reparsed_path, "`"))
  }

  return(fixed)
}

#' Decode and simplify FastMLST ST calls
#'
#' @param mlst tibble or data frame of raw FastMLST csv input
#' @param save should the re-parsed file be saved?
process_mlst <- function(mlst){
  .complete <- mlst %>%
    dplyr::mutate(mlst_note = dplyr::case_when(ST == "-" ~ "Missing Allele(s)",
                                 ST == "new_ST" ~ "New ST",
                                 ST == "new_alleles" ~ "New Allele(s)",
                                 TRUE ~ NA_character_)) %>%
    dplyr::mutate(mlst = dplyr::recode(ST, "new_ST" = NA_character_,
                         "new_alleles" = NA_character_,
                         "-" = NA_character_)) %>%
    dplyr::mutate(assembly = as.character(stringr::str_match(Genome,"GCA_\\d+\\.\\d+"))) %>%
    dplyr::select(assembly, mlst, mlst_note)

  return(.complete)
}

#' Download protein sequences from NCBI
#'
#' Easily download a set of protein sequences by accession number using the
#' Eutilities web interface via the REntrez package
#'
#' @param id character vector of NCBI protein accession numbers
#' @param n how many sequences should be downloaded in a single request?
download_protein_sequences <- function(id, n = 100){
  data <- tibble::tibble(protein = character(), protein_sequence = character())
  progress <- progress::progress_bar$new(format = "Downloaded :current of :total sequence sets :bar :percent (Elapsed: :elapsedfull, Remaining: :eta)",
                                         total = ceiling(length(id) / n), show_after = 0)
  progress$tick(0)

  for(start in seq(1, length(id), n)){
    end <- min(start + n - 1, length(id))

    entrez_xml <- rentrez::entrez_fetch(db = "protein", id = id[start:end],
                                        rettype = "fasta", retmode = "xml") %>%
      xml2::read_xml()

    accn <- entrez_xml %>%
      xml2::xml_find_all("//TSeqSet/TSeq/TSeq_accver") %>%
      xml2::xml_text()

    seq <- entrez_xml %>%
      xml2::xml_find_all("//TSeqSet/TSeq/TSeq_sequence") %>%
      xml2::xml_text()

    data <- dplyr::add_row(data, protein = accn, protein_sequence = seq)

    progress$tick()
  }
  return(data)
}

#' Download (and save) all Reference Gene Catalog protein sequences from NCBI
#'
#' Easily download a set of protein sequences by accession number using the
#' Eutilities web interface via the REntrez package
#'
#' @param path path to an NCBI ReferenceGeneCatalog.txt file
#' @param file filename to save the accession/sequence pairs
download_reference_gene_catalog_proteins <- function(path = refgene_path,
                                                     file = "reference_proteins.tsv"){
  accessions <- import_reference_gene_catalog(path) %>%
    tidyr::drop_na(protein) %>%
    dplyr::pull(protein)

  .sequences <- download_protein_sequences(accessions)

  readr::write_tsv(.sequences, file = paste0(dirname(path),"/",file))

  return(.sequences)
}
