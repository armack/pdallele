
#' @importFrom dplyr %>%
#' @importFrom data.table %like%
NULL

########## Metadata and Clusters ##########

#' Download listing of available organisms
#'
#' Downloads listing of organisms currently available on the NCBI
#' Pathogen Detection FTP server.
#'
#' The FTP directory can be accessed via browser at
#' [https://ftp.ncbi.nlm.nih.gov/pathogen/Results/].
#'
#' @param url Path to the NCBI Pathogen Detection FTP Results directory
#' @returns A tibble with columns `organism` containing organism names and
#'   `updated` containing the date of the most recent update.
#' @export
get_organisms <- function(url = "ftp://ftp.ncbi.nlm.nih.gov/pathogen/Results/"){
  raw_rows <- RCurl::getURL(url = , ftp.use.epsv = FALSE, customrequest = "MLSD") %>%
    stringr::str_split("\n") %>%
    unlist() %>%
    stringr::str_subset("type=dir")

  organism <- raw_rows %>%
    stringr::str_extract("(?<=; ).+$")

  updated <- raw_rows %>%
    stringr::str_extract("(?<=modify\\=)\\d+(?=;)") %>%
    lubridate::ymd_hms(tz = "UTC")

  organism_list <- tibble::tibble(organism, updated) %>%
    dplyr::arrange(organism)

  return(organism_list)
}

#' Download listing of available PDG accession.version releases for a given organism
#'
#' Downloads listing of PDG (Pathogen Detection Organism Group)
#' accession.version releases currently available for `organism` on the NCBI
#' Pathogen Detection FTP server.
#'
#' @inheritParams get_organisms
#' @param organism FTP directory name of an organism group included in the NCBI
#'   Pathogen Detection Project
#' @returns A tibble with columns `version` containing available PDG
#'   accession.version releases `updated` containing the date of the most recent
#'   update.
#' @export
get_versions <- function(url = "ftp://ftp.ncbi.nlm.nih.gov/pathogen/Results/", organism){
  organism_url <- paste0(url, organism, "/")
  raw_rows <- RCurl::getURL(organism_url, ftp.use.epsv = FALSE, customrequest = "MLSD") %>%
    stringr::str_split("\n") %>%
    unlist() %>%
    stringr::str_subset("type=dir")

  version <- raw_rows %>%
    stringr::str_extract("(?<=; ).+$")
  updated <- raw_rows %>%
    stringr::str_extract("(?<=modify\\=)\\d+(?=;)") %>%
    lubridate::ymd_hms(tz = "UTC")
  version_list <- tibble::tibble(version, updated) %>%
    dplyr::arrange(stringr::str_rank(version, numeric = TRUE))

  return(version_list)
}

#' Verify a PDG accession.version release is complete
#'
#' Verifies that the specified PDG accession.version release is complete (i.e.
#' includes AMR Metadata, SNP Clusters, and SNP Trees). Note that SNP Trees are
#' part of the definition of a 'complete' release, but are not downloaded or
#' analyzed by `pdallele`
#'
#' @inheritParams get_versions
#' @param version A PDG accession.version corresponding to `organism`
#' @export
check_complete_pdg <- function(url = "ftp://ftp.ncbi.nlm.nih.gov/pathogen/Results/", organism, version){
  version_url <- paste0(url, organism, "/", version, "/", version, ".final.descriptor.xml")

  .complete <- xml2::read_xml(version_url) %>%
    xml2::xml_find_first("/PathogenGroupDescriptor/PathogenGroup") %>%
    xml2::xml_attr("complete") %>%
    as.integer() %>%
    as.logical()

  return(.complete)
}

#' Download Metadata and Clusters by organism and PDG accession.version
#'
#' @description Download Metadata and SNP Cluster files for the specified
#' `organism` / `version` combination. Verifies the release is complete before
#' downloading.
#'
#' The downloaded files are saved as `asm.metadata.tsv` and `cluster_list.tsv`
#' in `path`.
#'
#' @inheritParams check_complete_pdg
#' @param path Destination directory
#' @export
download_metadata_clusters <- function(url = "ftp://ftp.ncbi.nlm.nih.gov/pathogen/Results/", organism, version, path){
  url_base <- paste0(url, organism, "/", version, "/")
  url_descriptor <- paste0(url_base, version, ".final.descriptor.xml")
  url_metadata <- paste0(url_base, "AMR/", version, ".amr.metadata.tsv")
  url_cluster <- paste0(url_base, "Clusters/", version, ".reference_target.cluster_list.tsv")

  complete <- check_complete_pdg(organism = organism, version = version)

  if(!complete){
    stop("`", organism, "` version `", version, "` is not a complete dataset. Please try a different version.")
  }

  ensure_directory(path)

  metadata_path <- file.path(path, "amr.metadata.tsv")
  cluster_path <- file.path(path, "cluster_list.tsv")

  if(!file.exists(metadata_path)){
    utils::download.file(url = url_metadata, destfile = metadata_path)
    log <- file(file.path(path,"DatasetDetails.txt"), open="a")
    write(paste0(Sys.time(), " - Downloaded ", url_metadata, " as ",
                 basename(metadata_path)), file = log, append = TRUE)
    close(log)
  }else{
    message("Skipped `", basename(metadata_path), "`. File already exists.")
  }

  if(!file.exists(cluster_path)){
    utils::download.file(url = url_cluster, destfile = cluster_path)
    log <- file(file.path(path,"DatasetDetails.txt"), open="a")
    write(paste0(Sys.time(), " - Downloaded ", url_cluster, " as ",
                 basename(cluster_path)), file = log, append = TRUE)
    close(log)
  }else{
    message("Skipped `", basename(cluster_path), "`. File already exists.")
  }
}

########## Reference Gene Catalog ##########

#' Download listing of available Reference Gene Catalog releases
#'
#' Downloads listing of Reference Gene Catalog releases currently available on
#' the NCBI Pathogen Detection FTP server.
#'
#' @param url Path to the NCBI Pathogen Detection FTP Reference Gene Catalog
#'   directory
#' @returns A tibble with columns `version` containing available Reference Gene
#'   Catalog releases and `updated` containing the date of the most recent
#'   update.
#' @export

get_reference_gene_catalogs <- function(url = "ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/"){
  raw_rows <- RCurl::getURL(url, ftp.use.epsv = FALSE, customrequest = "MLSD") %>%
    stringr::str_split("\n") %>%
    unlist() %>%
    stringr::str_subset("type=dir")

  version <- raw_rows %>%
    stringr::str_extract("(?<=; ).+$")
  updated <- raw_rows %>%
    stringr::str_extract("(?<=modify\\=)\\d+(?=;)") %>%
    lubridate::ymd_hms(tz = "UTC")
  catalog_list <- tibble::tibble(version, updated) %>%
    dplyr::arrange(version)

  return(catalog_list)
}

#' Download Reference Gene Catalog
#'
#' @description
#' Downloads a specified `version` of the Reference Gene Catalog from the NCBI
#' Pathogen Detection FTP server.
#'
#' The downloaded file is saved as 'ReferenceGeneCatalog.tsv' in `path`.
#'
#' @inheritParams get_reference_gene_catalogs
#' @param version A Reference Gene Catalog version number
#' @param path Destination directory
#' @export
download_reference_gene_catalog <- function(url = "ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/Data/", version, path){
  url_complete <- file.path(url, version, "ReferenceGeneCatalog.tsv")

  ensure_directory(path)
  version_path <- file.path(path, "refgene.txt")

  if(!file.exists(version_path)){
    utils::download.file(url = url_complete,
                         destfile = version_path)
    log <- file(file.path(path,"DatasetDetails.txt"), open="a")
    write(paste0(Sys.time(), " - Downloaded ", url_complete, " as ",
                 basename(version_path)), file = log, append = TRUE)
    close(log)
  }else{
    message("Skipped `", basename(version_path), "`. File already exists.")
  }
}

########## Identical Protein Groups ##########

#' Download Identical Protein Groups data using NCBI Entrez
#'
#' @description
#' Download Identical Protein Group (IPG) details for a list of protein
#' accessions using the NCBI Entrez API via the rentrez package. The information
#' for `n` accessions is downloaded in each batch and the remaining accessions
#' filtered to reduce duplicate downloads.
#'
#' The downloaded data is saved as 'ipg.tsv' in `path`.
#'
#' A lower value for `n` will require more queries to NCBI servers but will
#' reduce redundant downloads if an IPG is associated with a large number of
#' accession numbers in the dataset.
#'
#' @param accessions A character vector of protein accession numbers
#' @param n Number of accessions to include in each request
#' @param path Destination directory
#' @export
download_identical_protein_groups <- function(accessions, n = 10, path){
  .ipg_data <- tibble(ipg = integer(), ipg_accession = character(),
                      ipg_name = character(), protein = character())
  remaining <- accessions

  progress <- progress::progress_bar$new(format = "Downloaded :current of :total accession numbers :bar :percent (Elapsed: :elapsedfull, Remaining: :eta)",
                                         total = length(accessions), show_after = 0)
  progress$tick(0)

  while(length(remaining) > 0){
    id <- remaining[1:min(n - 1, length(remaining))]

    xml_ipg_report <- xml2::read_xml(rentrez::entrez_fetch(db = "ipg", id = id, retmode = "docsum", rettype = "xml")) %>%
      xml2::xml_find_all("/IPGReportSet/IPGReport")
    if(length(xml_ipg_report) > 0){
      ipg_additional <- xml_ipg_report %>%
        purrr::map(
          ~tibble::tibble(
            ipg = as.integer(xml2::xml_attr(.x, "ipg")),
            ipg_accession = xml2::xml_find_first(.x, "Product") %>%
              xml2::xml_attr("accver"),
            ipg_name = xml2::xml_find_first(.x, "Product") %>%
              xml2::xml_attr("name"),
            protein = c(xml2::xml_attr(.x, "product_acc"), xml2::xml_find_all(.x, "ProteinList/Protein") %>%
                          xml2::xml_attr("accver"))
          )) %>%
        purrr::list_rbind() %>%
        dplyr::filter(.data$protein %in% accessions) %>%
        dplyr::distinct()

      .ipg_data <- dplyr::bind_rows(.ipg_data, ipg_additional)
      remaining <- dplyr::setdiff(remaining, unlist(ipg_additional[["protein"]]))
    }
    remaining <- setdiff(remaining, id)
    progress$update((length(accessions) - length(remaining)) / length(accessions))
  }

  table_path <- file.path(path, "ipg.tsv")
  readr::write_tsv(.ipg_data, file = table_path)

  log <- file(file.path(path,"DatasetDetails.txt"), open="a")
  write(paste0(Sys.time(), ' - Queried data on ', nrow(.ipg_data),
               ' Identical Protein Groups covering ', length(accessions),
               ' protein accession numbers and saved the results as ',
               basename(table_path)), file = log, append = TRUE)
  close(log)
}

########## MicroBIGG-E ##########

#' Download MicroBIGG-E Data from the Google Cloud Platform
#'
#' See bigrquery help for more information about configuration and authentication
#'
#' The downloaded data is saved as 'microbigge.tsv' in `path`.
#'
#' Note that taxgroup, sciname, and element all use SQL syntax for queries. Use
#' '%' as a wildcard rather than '*'.
#'
#' @param billing GCP project identifier for billing/usage tracking
#' @param path directory to save downloaded files to
#' @param taxgroup an NCBI Pathogen Detection organism group to filter by
#' @param sciname scientific name of the organism of interest to filter by
#' @param element name of a gene/allele to filter by
#' @param debug should the query be previewed rather than executed?
download_microbigge_bq <- function(billing, path, taxgroup, sciname, element, debug = FALSE){
  . <- NULL # Workaround to suppress `no visible binding for global variable`

  bq_con <- bigrquery::dbConnect(
    bigrquery::bigquery(),
    project = "ncbi-pathogen-detect",
    dataset = "pdbrowser",
    billing
  )

  taxgroup_missing <- missing(taxgroup)
  element_missing <- missing(element)
  sciname_missing <- missing(sciname)

  table_path <- file.path(path, "microbigge.tsv")

  mbe_data <- bq_con %>%
    tbl("microbigge") %>%
    {if(!taxgroup_missing) filter(., .data$taxgroup_name %like% taxgroup) else . } %>%
    {if(!sciname_missing) filter(., .data$scientific_name %like% sciname) else . } %>%
    {if(!element_missing) filter(., .data$element_symbol %like% element) else . } %>%
    {if (debug) show_query(.) else collect(.) }

  sql_query <- bq_con %>%
    tbl("microbigge") %>%
    {if(!taxgroup_missing) dplyr::filter(., .data$taxgroup_name %like% taxgroup) else . } %>%
    {if(!sciname_missing) dplyr::filter(., .data$scientific_name %like% sciname) else . } %>%
    {if(!element_missing) dplyr::filter(., .data$element_symbol %like% element) else . } %>%
    dbplyr::sql_render() %>%
    stringr::str_replace_all("\n", " ") %>%
    stringr::str_replace_all("   ", " ") %>%
    stringr::str_replace("microbigge", "ncbi-pathogen-detect.pdbrowser.microbigge")

  if(!debug){
    readr::write_tsv(mbe_data, table_path)
    log <- file(file.path(path,"DatasetDetails.txt"), open="a")
    write(paste0(Sys.time(), ' - Ran Query "', sql_query,
                 '" on Google Cloud Platform and saved the results as  ',
                 basename(table_path)), file = log, append = TRUE)
    close(log)
    message(paste0("Successfully downloaded ", nrow(mbe_data), " rows covering ",
                   nrow(dplyr::distinct(mbe_data, .data$biosample_acc)), " isolates."))
  }
}

########## Interactive Download Helper Functions ##########

#' Interactively select an organism from the NCBI Pathogen Detection FTP directory
#'
#' @param url path to NCBI Pathogen Detection FTP "Results" directory
select_ftp_organism <- function(url){
  organism_list <- get_organisms(url = url) %>%
    dplyr::pull(organism)
  organism <- organism_list[utils::menu(stringr::str_replace_all(organism_list, "_", " "),
                                 title = "Select an Organism Group")]

  return(organism)
}

#' Interactively select an accession.version for a given organism from the NCBI
#' Pathogen Detection FTP directory
#'
#' @param url path to NCBI Pathogen Detection FTP "Results" directory
#' @param organism name of a valid organism/organism group
select_ftp_version <- function(url, organism){
  version_list <- get_versions(url = url, organism = organism) %>%
    dplyr::pull(version)
  version <- version_list[utils::menu(stringr::str_replace_all(version_list, "_", " "),
                               title = "Select a version. Higher numbers are newer.")]

  return(version)
}

#' Interactively select a ReferenceGeneCatalog version from the NCBI Pathogen
#' Detection FTP directory
#'
#' @param url path to NCBI Pathogen Detection FTP Reference Gene Catalog directory
select_ftp_refgene <- function(url){
  version_list <- get_refgene_catalogs(url = url) %>%
    dplyr::pull(version)
  version <- version_list[utils::menu(stringr::str_replace_all(version_list, "_", " "),
                               title = "Select a version. The newest is reccomended.")]

  return(version)
}
