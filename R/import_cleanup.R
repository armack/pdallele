#' @importFrom dplyr %>%
NULL

########## configuration ##########
unclear_species_string <- "Unclear or Unknown"

########## Isolates Browser Import and Processing ##########

#' Import, Select, and Trim NCBI Pathogen Detection Isolates Browser Metadata
#'
#' Imports Isolated Browser metadataa, selects and renames desired columns, and
#' cleans up problematic quotation marks
#'
#' @param path Path to an NCBI Isolates Browser
#'             PDG_accession.version.metadata.tsv file

import_isolates_browser_metadata <- function(path) {
  data <- readr::read_tsv(path, quote="", na = character(), show_col_types = FALSE)

  .complete <- data %>%
    dplyr::select(biosample = "biosample_acc", assembly = "asm_acc", pdt = "target_acc",
           bioproject = "bioproject_acc", "bioproject_center", "collection_date",
           sample_type = "epi_type", location = "geo_loc_name", "host", "host_disease",
           "source_type", "isolation_source", "lat_lon", "scientific_name",
           "species_taxid", organism_taxid = "taxid", amr = "AMR_genotypes",
           stress = "stress_genotypes", vir = "virulence_genotypes",
           "refgene_db_version" ) %>%
    dplyr::mutate(across(everything(), ~stringr::str_replace_all(.x, '"', "")))

    return(.complete)
}

#' Replace character vector values with NA across all character columns
#'
#' @param data A dataframe or tibble to modify
#' @param terms A character vector to replace values with NA_character_

na_if_tibble_chr <- function(data, terms) {
  .complete <- data %>%
    dplyr::mutate(across(
      dplyr::where(is.character),
      ~ dplyr::if_else(stringr::str_to_lower(.) %in% terms, NA_character_, .)
    ))

  return(.complete)
}

#' Parse scientific_name to genus and species
#'
#' Parses scientific_name into genus and species columns and their formatted
#' equivalents. Unpublished "<Genus> sp. <strain_identifier>" (See PMC7408187)
#' are treated as genus = <Genus> and species = sp.
#'
#' See "How Many Species" in PMC3245000 for info in "sp." strains
#' Also PMC7408187 (Curation of prokaryotes) and
#' https://www.ncbi.nlm.nih.gov/biosample/docs/organism/
#' And https://ena-docs.readthedocs.io/en/latest/faq/taxonomy_requests.html
#'
#' @param data A dataframe or tibble
#' @param path Path to an NCBI Isolates Browser
#'             PDG_accession.version.metadata.tsv file
parse_genus_species <- function(data, path) {

  .complete <- data %>%
    dplyr::mutate(species = stringr::str_extract(.data$scientific_name, "^[^ ]* [^ ]*")) %>%
    dplyr::mutate(species_markdown = paste0("*", .data$species ,"*")) %>%
    dplyr::mutate(species_math = paste0("italic(", .data$species ,")")) %>%
    dplyr::mutate(genus = stringr::str_extract(.data$scientific_name, "^[^ ]*")) %>%
    dplyr::mutate(genus_markdown = paste0("*", .data$genus ,"*")) %>%
    dplyr::mutate(genus_math = paste0("italic(", .data$genus ,")"))

  return(.complete)
}


#' Reverse geocode "lat_lon" to countries and use to fill in missing "location"
#'
#' @param data A dataframe or tibble to modify
reverse_geocode <- function(data){

  .partial <- data %>%
    dplyr::mutate(location = dplyr::if_else(grepl("/",location), NA_character_, location)) %>%
    tidyr::extract("lat_lon", into=c("lat","lat_dir","lon", "lon_dir"), regex="(\\d*\\.?\\d*) ([[:alpha:]]) (\\d*\\.?\\d*) ([[:alpha:]])",remove=FALSE, convert = TRUE) %>%
    dplyr::mutate(lat = dplyr::case_when(lat_dir == "N" ~ .data$lat, lat_dir == "S" ~ -.data$lat)) %>%
    dplyr::mutate(lon = dplyr::case_when(lon_dir == "E" ~ .data$lon, lon_dir == "W" ~ -.data$lon)) %>%
    dplyr::select(-"lat_lon", -"lat_dir", -"lon_dir")

  geocode_countries <- .partial %>%
    dplyr::distinct(.data$lat, .data$lon) %>%
    tidyr::drop_na() %>%
    dplyr::mutate(geocode_loc = maps::map.where(database = "world", x = .data$lon, y = .data$lat))

  .complete <- .partial %>%
    dplyr::left_join(geocode_countries, by = c("lat", "lon")) %>%
    dplyr::mutate(location = dplyr::coalesce("location", "geocode_loc")) %>%
    dplyr::select(-"geocode_loc")

  return(.complete)
}


#' Splits location into location_broad (country) and location_detail (city/area)
#'
#' Splits INSDC/BioSample geo_loc_name style country:detailed_location format
#' into location_broad and location_detail.
#'
#' @param data A dataframe or tibble
#'
split_location <- function(data) {

  .complete <- data %>%
    tidyr::separate(location, into = c("location_broad", "location_detail"), sep = ":",
             extra = "merge", fill = "right", remove = FALSE) %>%
    dplyr::mutate(location_detail = stringr::str_trim(.data$location_detail, side = "left"))

  return(.complete)
}

#' Adds Region information based on a file and the location_broad column
#'
#' Requires a csv/tsv file with a 'country' and 'region' column. Countries are
#' matched to location_broad using countrycode::countryname to account for names
#' given in as many languages as possible.
#'
#' This function will also append any other columns in the file to the dataset
#' (using the column name from the file). This allows the inclusion of multiple
#' regions or including a region_markdown or similar column
#'
#' @param data A dataframe or tibble
#' @param path path to a csv/tsv file containing columns 'country' and 'region'
import_regions <- function(data, path) {
  if(is.null(path)) return(data)

  regions <- readr::read_csv(file = path, show_col_types = FALSE)

  country_list <- dplyr::distinct(data, .data$location_broad)

  standardized_countries <- country_list %>%
    dplyr::mutate(std_country = countrycode::countryname(.data$location_broad,
                                                  destination = "country.name.en")) %>%
    dplyr::left_join(regions, by = c("std_country" = "country")) %>%
    tidyr::drop_na()

  .complete <- data %>%
    dplyr::left_join(standardized_countries, by = "location_broad")

  return(.complete)
}


#' Import and Add NCBI Pathogen Detection SNP Cluster Data
#'
#' @param data A dataframe or tibble
#' @param path Path to an NCBI Pathogen Detection
#' PDG_accession.version.reference_target.cluster_list.tsv file

import_cluster_list <- function(data, path) {
  if(missing(path) | is.null(path)) return(data)

  .clusters <- readr::read_tsv(file = path, show_col_types = FALSE) %>%
    dplyr::select(pds = "PDS_acc", biosample = "biosample_acc")

  .complete <- data %>%
    dplyr::left_join(.clusters, by = "biosample", multiple = "any")

  return(.complete)
}

#' Import and Add MLST Data
#'
#' @param data A dataframe or tibble
#' @param path Path to a (reprocessed) MLST CSV file

import_mlst <- function(data, path) {
  if(missing(path) | is.null(path)) return(data)

  .mlst <- readr::read_csv(file = path, show_col_types = FALSE) %>%
    dplyr::select(mlst = "ST", mlst_errors = "errors", "assembly")

  .complete <- data %>%
    dplyr::left_join(.mlst, by = "assembly", multiple = "any")

  return(.complete)
}

#' Separate genotypes (amr, stress, vir) by row to allele and allele_type
#'
#' Uses separate_rows to separate genotype columns (amr, stress, vir) originally
#' AMR_genotype, stress_genotype, and virulance_genotype in isolates browser.
#'
#' allele_type -> amr, stress, or vir (virulence)
#' allele_quality -> codes correcponding to reasons for non-"complete" alleles,
#'                   see description of clean_filter_alleles() for details
#'
#' @param data A dataframe or tibble with genotype columns (amr, stress, vir)
#' @param include Which types of alleles (`amr`, `stress`, `vir`) should be kept?
#' This is a regex - separate multiple valeus with `|`

separate_genotypes <- function(data, include = "amr"){

  .complete <- data %>%
    tidyr::pivot_longer(any_of(c("amr", "stress","vir")), names_to = "allele_type",
                 values_to = "allele") %>%
    dplyr::filter(grepl(paste(include, collapse = "|"), .data$allele_type)) %>%
    tidyr::drop_na("allele") %>%
    tidyr::separate_rows("allele", sep = ",") %>%
    tidyr::separate("allele", into = c("allele","allele_quality"), sep = "=", fill = "right")

  return(.complete)
}

#' Clean and filter allele listings
#'
#' allele_quality is a column referring to alleles that are less than perfect
#' matches, which include HMM derived alleles, end of contig issues, and point
#' mutations. See possible values and explanations at:
#' https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#genotype-categories
#'
#' This function also removes any duplicate allele/biosample combinations, which
#' could complicate later analysis. Distinct alleles are maintained.
#'
#' For alleles without a gene listed (i.e. unassigned alleles and lower quality
#' calls not matching an assigned allele), `allele` is coalesced into `gene`
#'
#' @param data A dataframe or tibble with an allele column
#' @param filter A character vector of gene names to keep. Uses grepl(). Partial
#'   matches are kept. Vector is collapsed with "|" (or operator).
#' @param remove A character vector of allele_quality values to remove

clean_filter_alleles <- function(data, filter, remove) {

  .complete <- data %>%
    dplyr::filter(grepl(paste(filter, collapse = "|"), .data$gene)) %>%
    dplyr::filter(!grepl(paste(remove, collapse = "|"), .data$allele_quality)) %>%
    dplyr::distinct(.data$biosample, .data$allele, .keep_all = TRUE)

  return(.complete)
}

#' Import Reference Gene Catalog to tibble
#'
#' Imports Reference Gene Catalog, renames columns, and removes unneeded data
#'
#' @param path Path to an NCBI ReferenceGeneCatalog.txt file
import_reference_gene_catalog <- function(path){
  .rgc <- readr::read_tsv(path, show_col_types = FALSE) %>%
    dplyr::select("allele", gene = "gene_family", name = "product_name", "type", "subtype",
           "class", "subclass", protein = "refseq_protein_accession",
           nucleotide = "refseq_nucleotide_accession", "blacklisted_taxa",
           "whitelisted_taxa")

  return(.rgc)
}


#' Add Reference Gene Catalog information to `data` by allele
#'
#' Filters after joining to remove invalid matches when a row contains
#' whitelisted_taxa. Coalesces empty 'gene' column with 'allele' column
#' to deal with missing Reference Gene genes
#'
#' @param data A dataframe or tibble to add to
#' @param path Path to an NCBI ReferenceGeneCatalog.txt file

add_reference_gene_catalog <- function(data, path) {
  if(is.null(path)) return(data)

  .rgc <- import_reference_gene_catalog(path) %>%
    dplyr::mutate(whitelisted_taxa = stringr::str_replace(.data$whitelisted_taxa, "_", " "))

  .complete <- dplyr::left_join(data, .rgc, by = "allele", multiple = "all") %>%
    dplyr::filter(dplyr::case_when(is.na(.data$whitelisted_taxa) ~ TRUE,
                     whitelisted_taxa == genus ~ TRUE,
                     whitelisted_taxa == species ~ TRUE,
                     TRUE ~ FALSE)) %>%
    dplyr::select(-"whitelisted_taxa") %>%
    dplyr::mutate(gene = dplyr::coalesce("gene", "allele"))

  return(.complete)
}

#' Properly format  bla gene/allele names per ASM Journals guidelines
#'
#' Add formatting for proper display of bla allele names in both markdown (used
#' by ggtext) and 'math' formatting (used by ggplot with expression()).
#'
#' The format is "bla" in italic and gene/allele name in subscript, e.g. using
#' HTML tags: <i>bla</i><sub>PDC-3</sub>. See PMID 35380458 "Consensus on
#' β-Lactamase Nomenclature" by Bush et al. for details
#'
#' @param data A dataframe or tibble with an 'allele' and 'gene' column
parse_bla_formatting <- function(data) {

  .complete <- data %>%
    dplyr::mutate(allele_markdown = sub("^bla(.*?)$", "<i>bla</i><sub>\\1</sub>", .data$allele, perl = TRUE) ) %>%
    dplyr::mutate(allele_math = sub("^bla(.*?)$", "italic\\(bla\\)\\[\\1\\]", .data$allele, perl = TRUE) ) %>%
    dplyr::mutate(gene_markdown = sub("^bla(.*?)$", "<i>bla</i><sub>\\1</sub>", .data$gene, perl = TRUE) ) %>%
    dplyr::mutate(gene_math = sub("^bla(.*?)$", "italic\\(bla\\)\\[\\1\\]", .data$gene, perl = TRUE) )

  return(.complete)
}

#' Add oxa_family parsed from name (originally from Reference Gene Catalog file)
#'
#' Also adds markdown and math formatted oxa_family columns
#'
#' Returns the original table if no OXA alleles are present
#'
#' @param data A dataframe or tibble

parse_ib_oxa_family <- function(data) {
  if(!any(grepl("OXA", data$allele))){return(data)}

  .oxa_families <- data %>%
    dplyr::distinct(.data$allele, .data$name) %>%
    dplyr::filter(grepl("OXA-", .data$allele, ignore.case = TRUE)) %>%
    dplyr::mutate(oxa_family = dplyr::case_when(
      grepl('family', .data$name) ~ paste0("bla",stringr::str_extract(.data$name, "OXA-\\d+?(?= family)")),
      TRUE ~ paste0("bla",stringr::str_extract(.data$name, "OXA-\\d+"))
    )) %>%
    dplyr::mutate(oxa_family_markdown = sub("^bla(.*?)$", "<i>bla</i><sub>\\1</sub>",
                                     .data$oxa_family, perl = TRUE) ) %>%
    dplyr::mutate(oxa_family_math = sub("^bla(.*?)$", "italic\\(bla\\)\\[\\1\\]",
                                 .data$oxa_family, perl = TRUE) ) %>%
    dplyr::select(-"name")

  .complete <- dplyr::left_join(data, .oxa_families, by = "allele")

  return(.complete)
}

#' Add year parsed from collection_date
#'
#' @param data A dataframe or tibble

parse_year <- function(data) {
  .complete <- data %>%
    dplyr::mutate(year = sub("\\d{4}/(\\d{4})", "\\1", .data$collection_date)) %>%
    dplyr::mutate(year = sub("^.*?(\\d{4}).*?$", "\\1", .data$year)) %>%
    dplyr::mutate(across("year", as.integer))

  return(.complete)
}

#' Filter one biosample for each pds/genotype combination (while also keeping
#' non-clustered biosamples)
#'
#' Filters data to keep only a single biosample from each pds/genotype (in this
#' context, the combination of alleles present in a biosample) combination while
#' also keeping all biosamples without a pds cluster assigned.
#'
#' This may allow one to focus more on broader diversity by removing duplicates
#' of successful clones from the dataset
#'
#' @param data A dataframe or tibble
#' @param ignore_genotype should only a single biosample per cluster be kept,
#'                        even when multiple genotypes exist?

filter_pds <- function(data, ignore_genotype = FALSE){
  . <- NULL # Workaround to suppress `no visible binding for global variable`

  cluster_biosamples <- data %>%
    dplyr::select("biosample", "pds", "allele") %>%
    tidyr::drop_na("pds") %>%
    dplyr::group_by("biosample") %>%
    dplyr::mutate(combo = paste(.data$allele, collapse = ",")) %>%
    dplyr::group_by("pds") %>%
    {if(!ignore_genotype) dplyr::group_by(., "combo", .add = TRUE) else .} %>%
    dplyr::arrange("biosample") %>%
    dplyr::slice_head() %>%
    dplyr::pull("biosample")

  noncluster_biosamples <- data %>%
    dplyr::filter(is.na("pds")) %>%
    dplyr::pull("biosample") %>%
    unique()

  .complete <- data %>%
    dplyr::filter(.data$biosample %in% c(cluster_biosamples,noncluster_biosamples))

  return(.complete)
}

########## MicroBIGG-E Import and Processing ##########

#' Import and Select Columns from NCBI Pathogen Detection MicroBIGG-E File
#'
#' Imports NCBI Pathogen DetectionMicrobial Browser for Identification of
#' Genetic and Genomic Elements (MicroBIGG-E) data and selects and renames
#' desired columns for downstream processing steps
#'
#' As of October 2022, up to 100,000 rows of data can be downloaded via the web
#' interface at: https://www.ncbi.nlm.nih.gov/pathogens/microbigge/
#'
#' Larger sets can be downloaded using Google Big Query. See instructions from
#' NCBI at https://www.ncbi.nlm.nih.gov/pathogens/docs/microbigge_gcp/
#'
#' @param path Path to an NCBI MicroBIGG-E csv/tsv file.

import_microbigge <- function(path){
  .complete <- readr::read_tsv(file = path, na = "", show_col_types = FALSE) %>%
    dplyr::select(protein = "protein_acc", biosample = "biosample_acc", assembly = "asm_acc",
           allele = "element_symbol", long_name = "element_name", method = "amr_method",
           coverage = "pct_ref_coverage", identity = "pct_ref_identity",
           nucleotide = "contig_acc", nuc_start = "start_on_contig",
           nuc_stop = "end_on_contig", strand = "strand"
    )

  return(.complete)
}

#' Filter imported MBE data to remove less "reliable" protein calls
#'
#' Also allows for the customization of coverage and identity thresholds. The
#' default identity of 90% is chosen to match the Isolates Browser criteria
#' as stated at:
#' https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#genotype-categories
#'
#' ********** THIS REALLY NEEDS TO BE REVISITED **********
#'
#' @param data a data tabel or tibble
#' @param min_coverage integer minimum percent coverage to keep
#' @param min_identity integer minimum percent identity to keep
#' @param remove regex string of `method` to remove
filter_microbigge <- function(data, min_coverage = 100, min_identity = 90, remove) {
  . <- NULL # Workaround to suppress `no visible binding for global variable`
  remove_missing <- missing(remove)
  .complete <- data %>%
    tidyr::drop_na("protein") %>%
    {if(!remove_missing) dplyr::filter(., !grepl(paste(remove, collapse = "|"), .data$method)) else .} %>%
    dplyr::filter("coverage" >= min_coverage) %>%
    dplyr::filter("coverage" >= min_identity)

  return(.complete)
}

#' Determines OXA family from MBE data and adds a 'family' column
#'
#' @param data A dataframe or tibble
#' @param by passed to left_join as a join by command; the column(s) to match in the join

parse_mbe_oxa_family <- function(data, by = "allele") {
  .complete <- data %>%
    dplyr::mutate(mbe_family = dplyr::case_when(
      grepl('family', long_name) ~ stringr::str_extract(long_name, "OXA-\\d+?(?= family)"),
      TRUE ~ stringr::str_extract(long_name, "OXA-\\d+")
    )) %>%
    dplyr::mutate(ipg_family = dplyr::case_when(
      grepl('family', ipg_name) ~ stringr::str_extract(ipg_name, "OXA-\\d+?(?= family)"),
      TRUE ~ stringr::str_extract(ipg_name, "OXA-\\d+")
    )) %>%
    dplyr::mutate(family = dplyr::coalesce("ipg_family", "mbe_family"))

  return(.complete)
}


#' Import NCBI Identical Protein Groups data and add accession numbers and names
#' to MBE data
#'
#' This data can be obtained by downloading via the web interface, although it
#' seems rather inefficient to download the entire listing for all bla alleles
#' present in an organism.
#'
#' @param data a data table or tibble
#' @param path path to an NCBI IPG csv/tsv file
import_ipg <- function(data, path){
  if(is.null(path)) return(data)

  .ipg <- readr::read_tsv(file = path, show_col_types = FALSE)

  mbe_proteins <- data %>%
    dplyr::pull("protein") %>%
    unique()

  mbe_to_ipg <- .ipg %>%
    dplyr::select(ipg_uid = "Id", accession = "Protein") %>%
    dplyr::filter(.data$accession %in% mbe_proteins) %>%
    dplyr::distinct(.data$accession, .keep_all = TRUE)

  ipg_uid_accession_name <- .ipg %>%
    dplyr::group_by("Id") %>%
    dplyr::arrange(!grepl("REFSEQ", .data$Source)) %>%
    dplyr::select(ipg_uid = "Id", ipg_accession = "Protein", ipg_name = "Protein Name") %>%
    dplyr::slice_head()

  .complete <- data %>%
    dplyr::left_join(mbe_to_ipg, by = c("protein" = "accession")) %>%
    dplyr::left_join(ipg_uid_accession_name, by = "ipg_uid")

  return(.complete)
}

#' Add reference protein sequences
#'
#' Add column "protein_sequence" to a tibble by dplyr::left_join() on "protein"
#' (containing accession numbers)
#'
#' The source data can be generated from a ReferenceGeneCatalog.txt file using
#' download_reference_gene_catalog_proteins()
#'
#' @param data a data table or tibble with a "protein" column of accession numbers
#' @param path path to a TSV containing accession number/sequence pairs
add_ref_protein_sequences <- function(data, path){
  if(is.null(path)) return(data)

  .complete <- data %>%
    dplyr::left_join(readr::read_tsv(path, show_col_types = FALSE), by = "protein")

  return(.complete)
}
