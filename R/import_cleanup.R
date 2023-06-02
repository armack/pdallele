#' @importFrom dplyr %>%
NULL

########## Isolates Browser Metadata Import and Processing ##########

#' Import, select, and trim NCBI Pathogen Detection Isolates Browser Metadata
#'
#' Imports Isolates Browser Metadata, selects and renames relevant columns, and
#' removes quotation marks ()
#'
#' @details
#' The following columns are included in the output:
#'
#' |**Source File**     |**Output Tibble** |
#' |--------------------|-------------------|
#' |biosample_acc       |biosample          |
#' |asm_acc             |assembly           |
#' |pdt                 |target_acc         |
#' |bioproject_acc      |bioproject         |
#' |bioproject_center   |bioproject_center  |
#' |collection_date     |collection_date    |
#' |epi_type            |sample_type        |
#' |geo_loc_name        |location           |
#' |host                |host               |
#' |host_disease        |host_disease       |
#' |source_type         |source_type        |
#' |isolation_source    |isolation_source   |
#' |lat_lon             |lat_lon            |
#' |scientific_name     |scientific_name    |
#' |species_taxid       |species_taxid      |
#' |taxid               |organism_taxid     |
#' |AMR_genotypes       |amr                |
#' |stress_genotypes    |stress             |
#' |virulence_genotypes |virulence          |
#' |refgene_db_version  |refgene_db_version |
#'
#' @param path Path to an NCBI Pathogen Detection Project 'amr.metadata.tsv'
#'   file
#' @returns A [tibble()]. Values are described in Details.
#' @export

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

#' Replace with `NA` across character columns
#'
#' @description Replace values in all `character()` columns in `data` matching a
#'   value in `terms` with `NA`. Ignores case and requires a complete match.
#'
#'   For example, "null" would replace "NULL" with `NA` but leave "nullify"
#'   intact.
#'
#' @section Suggested `terms` Values: Unfortunately, due to the presence of many fields
#'   in BioSample metadata allowing free text entry and no enforced standard
#'   for how to indicate data is unavailable, a wide variety of strings need to
#'   be treated as `NA` in order to maximize the utility of the metadata.
#'
#'   The following suggestions come from examining more than 40,000
#'   *Acinetobacter* and *Pseudomonas* isolates present in the Pathogen
#'   Detection Project databases. It is by no means a complete list, but should
#'   serve as a good starting point to clean your data sets.
#'
#' ```
#' c("na", "n/a", "missing", "none", "unknown", "unknow",
#' "unspecified", "not known", "not applicable", "not_applicable", "not
#' provided", "not available", "not determined", "not specified", "not
#' recorded", "no data", "notfound", "null", "not collected", "not available: to
#' be reported later", "-")
#' ```
#'
#' @param data A dataframe or tibble to modify
#' @param terms A character vector of values to replace with `NA_character_`
#' @return `data` with matching values matching `terms` replaced by
#'   `NA_character_`
#' @export

na_if_tibble_chr <- function(data, terms) {
  .complete <- data %>%
    dplyr::mutate(across(
      dplyr::where(is.character),
      ~ dplyr::if_else(stringr::str_to_lower(.) %in% stringr::str_to_lower(terms), NA_character_, .)
    ))

  return(.complete)
}

#' Parse scientific_name to genus and species
#'
#' @description Parses scientific_name into genus and species columns and their
#'   formatted equivalents. Unpublished "<Genus> sp. <strain_identifier>" (See
#'   PMC7408187) are treated as genus = <Genus> and species = sp.
#'
#' @details See "How Many Species" in PMC3245000 for info in "sp." strains Also
#'   PMC7408187 (Curation of prokaryotes) and
#'   <https://www.ncbi.nlm.nih.gov/biosample/docs/organism/> and
#'   <https://ena-docs.readthedocs.io/en/latest/faq/taxonomy_requests.html>
#'
#' @param data A dataframe or tibble with a `scientific_name` column
#' @returns `data` with columns `species`, `species_markdown`, `species_math`,
#'   `genus`, `genus_markdown`, and `genus_math` added
#' @export
parse_genus_species <- function(data) {

  .complete <- data %>%
    dplyr::mutate(species = stringr::str_extract(.data$scientific_name, "^[^ ]* [^ ]*")) %>%
    dplyr::mutate(species_markdown = paste0("*", .data$species ,"*")) %>%
    dplyr::mutate(species_math = paste0("italic(", .data$species ,")")) %>%
    dplyr::mutate(genus = stringr::str_extract(.data$scientific_name, "^[^ ]*")) %>%
    dplyr::mutate(genus_markdown = paste0("*", .data$genus ,"*")) %>%
    dplyr::mutate(genus_math = paste0("italic(", .data$genus ,")"))

  return(.complete)
}

#' Parse `lat_lon` data and reverse geocode to `location`
#'
#' @description Separates `lat_lon` and convert 'N/S/E/W' labels to '+/-'
#'   format. Uses these values to reverse geocode and determine country when
#'   `location` is not provided in the original metadata.
#'
#'   Also removes (and reverse geocodes if possible) any `location` values
#'   containing a forward slash. These appear to sometimes be used in metadata
#'   to indicate samples tied to multiple locations (e.g. a hospital in one
#'   country treating a patient who acquired an infection in another country),
#'   but is not part of the INSDC standard and is not consistently reliable.
#'
#' @param data A dataframe or tibble to containing a `lat_lon` column
#' @returns `data` with column `lat_lon` replaced by columns `lat` and `lon` and
#'  missing values for `location` added where possible
#' @export
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
    dplyr::mutate(location = dplyr::coalesce(.data$location, .data$geocode_loc)) %>%
    dplyr::select(-"geocode_loc")

  return(.complete)
}

#' Split location into location_broad and location_detail
#'
#' Splits INSDC/BioSample geo_loc_name style country:detailed_location format
#' into location_broad (country) and location_detail (city/area).
#'
#' @param data A dataframe or tibble
#' @returns `data` with columns `location_broad` and `location_detail`
#'   added
#' @export
split_location <- function(data) {

  .complete <- data %>%
    tidyr::separate(location, into = c("location_broad", "location_detail"), sep = ":",
             extra = "merge", fill = "right", remove = FALSE) %>%
    dplyr::mutate(location_detail = stringr::str_trim(.data$location_detail, side = "left"))

  return(.complete)
}

#' Adds `region` column based on `location_broad` and a file
#'
#' Requires a csv file with a 'country' and 'region' column. Countries are
#' matched to location_broad using countrycode::countryname to account for names
#' given in as many languages as possible.
#'
#' This function will also append any other columns in the file to the dataset
#' (using the column name from the file). This allows the inclusion of multiple
#' regions or including a region_markdown or similar column
#'
#' @param data A dataframe or tibble
#' @param path Path to a csv file containing at least columns `country` and
#'   `region`
#' @param warn Should a warning message with unmatched countries be shown?
#' @returns `data` with column `region` (and possibly others) added
#' @export
import_regions <- function(data, path, warn = FALSE) {
  if (is.null(path))
    return(data)

  regions <- readr::read_csv(file = path, show_col_types = FALSE)

  country_list <- dplyr::distinct(data, .data$location_broad)

  standardized_countries <- country_list %>%
    dplyr::mutate(
      std_country = countrycode::countryname(
        .data$location_broad,
        destination = "country.name.en",
        warn = warn
      )
    ) %>%
    dplyr::left_join(regions, by = c("std_country" = "country")) %>%
    tidyr::drop_na()

  .complete <- data %>%
    dplyr::left_join(standardized_countries, by = "location_broad")

  return(.complete)
}


#' Import and add NCBI Pathogen Detection SNP Cluster data
#'
#' Import NCBI **P**athogen **D**etection **S**NP Cluster (PDS) data from `path`
#' and add a `pds` column containing a PDS accession numbers to `data`.
#'
#' @param data A dataframe or tibble containing a `biosample` column
#' @param path Path to an NCBI Pathogen Detection cluster_list.tsv file'
#' @returns `data` with column `pds` added
#' @export

import_cluster_list <- function(data, path) {
  if(missing(path) | is.null(path)) return(data)

  .clusters <- readr::read_tsv(file = path, show_col_types = FALSE, na = "NULL") %>%
    dplyr::select(pds = "PDS_acc", biosample = "biosample_acc")

  .complete <- data %>%
    dplyr::left_join(.clusters, by = "biosample", multiple = "any")

  return(.complete)
}

#' Import and add MLST data
#'
#' Imports MLST calls from a FastMLST
#'
#' To ensure compatibility (particularly around errors), the original FastMLST
#' output file should be reprocessed using <<FUNCTION NAME HERE>> before running
#' this function.
#'
#' For more information on FastMLST, see: <>
#'
#' @param data A dataframe or tibble containing an `assembly` column
#' @param path Path to a (reprocessed) FastMLST CSV file
#' @returns `data` with columns `mlst` and `mlst_errors` added

import_mlst <- function(data, path) {
  if(missing(path) | is.null(path)) return(data)

  .mlst <- readr::read_csv(file = path, show_col_types = FALSE) %>%
    dplyr::select(mlst = "ST", mlst_errors = "errors", "assembly")

  .complete <- data %>%
    dplyr::left_join(.mlst, by = "assembly", multiple = "any")

  return(.complete)
}

#' Separate genotypes by row
#'
#' @description Uses separate_rows to separate genotype columns (`amr`,
#'   `stress`, `vir`) into rows, with each row corresponding to a single allele.
#'   Three columns of allele details are added:
#'
#' * `allele` corresponding to the allele name
#' * `allele_type` corresponding to the source column (`amr`, `stress`, or `vir`)
#' * `allele_call` corresponding to codes used by the NCBI Pathogen Detection Project to describe the method used for calling a given allele. See [clean_filter_alleles()] for details.
#'
#' @param data A dataframe or tibble with at least one `amr`, `stress`, or `vir`
#'   column
#' @param include Which types of alleles (`amr`, `stress`, `vir`) should be
#'   kept? Separate values with '|'.
#' @returns `data` with columns `allele`, `allele_type`, and `allele_call` added
#' @export

separate_genotypes <- function(data, include = "amr"){

  .complete <- data %>%
    tidyr::pivot_longer(any_of(c("amr", "stress","vir")), names_to = "allele_type",
                 values_to = "allele") %>%
    dplyr::filter(grepl(paste(include, collapse = "|"), .data$allele_type)) %>%
    tidyr::drop_na("allele") %>%
    tidyr::separate_rows("allele", sep = ",") %>%
    tidyr::separate("allele", into = c("allele","allele_call"), sep = "=", fill = "right")

  return(.complete)
}

#' Clean and filter allele listings
#'
#' @description The column `allele_call` contains codes representing to the
#'   calling method used for a given allele. Some of these codes correspond to
#'   "less than perfect" matches and may warrant removal depending on the scope
#'   and goals of the analysis. Detailed explanations are available at:
#'   <https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#genotype-categories>
#'
#'   This function removes rows with `allele_call` matching a value in `remove`.
#'
#'   This function also (optionally) keeps only rows with `allele` matching a
#'   value in `filter`, such as when only a specific subset of `allele` is
#'   needed in the dataset.
#'
#'   This function also removes any duplicate allele/biosample combinations,
#'   which could complicate later analysis as some isolates list the same allele
#'   multiple times across contigs. Distinct alleles are maintained.
#'
#' @param data A dataframe or tibble with an `allele` column
#' @param filter A character vector of `allele` names to keep. Uses grepl().
#'   Partial matches are kept. Vector is collapsed with "|" (or operator).
#' @param remove A character vector of allele_call values to remove
#' @returns `data` with rows removed according to the selected parameters
#' @export

clean_filter_alleles <- function(data, filter, remove) {
  missing_filter <- missing(filter)
  missing_remove <- missing(remove)

  .complete <- data %>%
    dplyr::filter(grepl(paste(filter, collapse = "|"), .data$allele)) %>%
    dplyr::filter(!grepl(paste(remove, collapse = "|"), .data$allele_call)) %>%
    dplyr::distinct(.data$biosample, .data$allele, .keep_all = TRUE)

  return(.complete)
}

#' Import Reference Gene Catalog to tibble
#'
#' Imports Reference Gene Catalog, renames columns, and drops data not needed
#' for downstream processing
#'
#' @details
#' The following columns are included in the output:
#'
#' |**Source File**          |**Output Tibble**|
#' |-------------------------|-----------------|
#' |allele                   |allele           |
#' |gene_family              |gene             |
#' |product_name             |name             |
#' |type                     |type             |
#' |subtype                  |subtype          |
#' |class                    |class            |
#' |subclass                 |subclass         |
#' |refseq_protein_accession |protein          |
#' |blacklisted_taxa         |blacklisted_taxa |
#' |whitelisted_taxa         |whitelisted_taxa |
#'
#' @param path Path to an NCBI ReferenceGeneCatalog.txt file
#' @returns A [tibble()].  Values are described in Details.
#' @export
import_reference_gene_catalog <- function(path){
  .rgc <- readr::read_tsv(path, show_col_types = FALSE) %>%
    dplyr::select("allele", gene = "gene_family", name = "product_name", "type", "subtype",
           "class", "subclass", protein = "refseq_protein_accession",
           nucleotide = "refseq_nucleotide_accession", "blacklisted_taxa",
           "whitelisted_taxa")

  return(.rgc)
}


#' Join Reference Gene Catalog metadata to `data` by allele
#'
#' @description Joins Reference Gene Catalog metadata from `path` with `data` by
#' `allele`
#'
#' Filters after joining to remove invalid matches when a row contains
#' whitelisted_taxa.
#'
#' Coalesces empty 'gene' column with 'allele' column to deal with missing
#' Reference Gene genes
#'
#' @param data A dataframe or tibble containing `allele` column
#' @param path Path to an NCBI ReferenceGeneCatalog.txt file
#' @returns `data` with columns `gene`, `name`, `type`, `subtype`, `class`,
#'   `subclass`, and `protein` added

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
    dplyr::mutate(gene = dplyr::coalesce(.data$gene, .data$allele))

  return(.complete)
}

#' Properly format  *bla* gene/allele names
#'
#' Add formatting of formatted strings for proper display of *bla* allele names
#' in both html (used by [ggtext::element_markdown()]) and 'math' formatting
#' (used by [ggplot2::element_text()] when wrapped in `expression()`).
#'
#' @details This formatting is part of the standardized beta-lactamase
#'   nomenclature adopted by ASM Journals. It was developed by a group of
#'   leading beta-lactamase researchers in conjunction with the NCBI Pathogen
#'   Detection Project, which is responsible for assigning most beta-lactamase
#'   allele designations.
#'
#'   The format is "bla" in italic and gene/allele name in subscript, e.g.
#'   <i>bla</i><sub>PDC-3</sub> from `<i>bla</i><sub>PDC-3</sub>`.
#'
#'   See PMID 35380458 "Consensus on β-Lactamase Nomenclature" by Bush et al.
#'   for further details
#'
#' @param data A dataframe or tibble containing`allele` and `gene` columns
#' @returns `data` with columns `allele_markdown`, `allele_math`,
#'   `gene_markdown`, and `gene_math` added
#' @export
parse_bla_formatting <- function(data) {

  .complete <- data %>%
    dplyr::mutate(allele_markdown = sub("^bla(.*?)$", "<i>bla</i><sub>\\1</sub>", .data$allele, perl = TRUE) ) %>%
    dplyr::mutate(allele_math = sub("^bla(.*?)$", "italic\\(bla\\)\\[\\1\\]", .data$allele, perl = TRUE) ) %>%
    dplyr::mutate(gene_markdown = sub("^bla(.*?)$", "<i>bla</i><sub>\\1</sub>", .data$gene, perl = TRUE) ) %>%
    dplyr::mutate(gene_math = sub("^bla(.*?)$", "italic\\(bla\\)\\[\\1\\]", .data$gene, perl = TRUE) )

  return(.complete)
}

#' Parse `oxa_family` from Reference Gene Catalog names
#'
#' @description Uses regular expressions to determine OXA families based on
#'   <i>bla</i><sub>OXA</sub> Reference Gene Catalog names/descriptions and adds
#'   columns `oxa_family`, `oxa_family_markdown`, and `oxa_family_math`. Other
#'   alleles produce `NA` values.
#'
#'   Also adds HTML (`oxa_family_markdown`) and "math" (`oxa_family_math`)
#'   formatted versions for plots. See [parse_bla_formatting()] for details.
#'
#'   Returns `data` unchanged if no <i>bla</i><sub>OXA</sub> alleles are
#'   present.
#'
#'   Two variants are available:
#'
#'   * `parse_ib_oxa_family()` for Isolates Browser data
#'   * `parse_mbe_oxa_family()` for MicroBIGG-E data
#'
#' @param data A dataframe or tibble containing `allele` and `name` columns
#' @returns `data` with columns `oxa_family`, `oxa_family_markdown`, and
#'   `oxa_family_math` added
#' @export

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

#' Parse `year` from `collection_date`
#'
#' @param data A dataframe or tibble with `collection_date` column
#' @returns `data` with column `year` added
#' @export

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
    dplyr::group_by(.data$biosample) %>%
    dplyr::mutate(combo = paste(.data$allele, collapse = ",")) %>%
    dplyr::group_by(.data$pds) %>%
    {if(!ignore_genotype) dplyr::group_by(., .data$combo, .add = TRUE) else .} %>%
    dplyr::arrange(.data$biosample) %>%
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

#' Import and Select Columns from NCBI Pathogen Detection MicroBIGG-E file
#'
#' @description
#' Imports NCBI Pathogen Detection Microbial Browser for Identification of
#' Genetic and Genomic Elements (MicroBIGG-E) data and selects and renames
#' desired columns for downstream processing steps
#'
#' Two variants are available:
#' * `import_microbigge_ncbi()` imports files downloded from the NCBI Website
#' * `import_microbigge_gcp()` imports files downloaded from Google Cloud Platform
#'
#' Slightly different column names are used between NCBI and GCP, resulting in
#' two slightly different import functions.
#'
#' @details
#' The following columns are included in the output:
#'
#' |**Source File**  |**Output Tibble**|
#' |-----------------|-----------------|
#' |protein_acc      |protein          |
#' |biosample_acc    |biosample        |
#' |asm_acc          |assembly         |
#' |element_symbol   |allele           |
#' |element_name     |name             |
#' |amr_method       |method           |
#' |pct_ref_coverage |coverage         |
#' |pct_ref_identity |identity         |
#' |contig_acc       |nucleotide       |
#' |start_on_contig  |nuc_start        |
#' |end_on_contig    |nuc_stop         |
#' |strand           |strand           |
#'
#' @section Obtaining Data: Up to 100,000 rows of data (potentially with an
#'   additional 50 MB cap) can be downloaded via the web interface at:
#'   [https://www.ncbi.nlm.nih.gov/pathogens/microbigge/]
#'
#'   Larger sets can be downloaded using Google Big Query. See instructions from
#'   NCBI at [https://www.ncbi.nlm.nih.gov/pathogens/docs/microbigge_gcp/]
#'
#' @param path Path to a MicroBIGG-E tsv file
#' @returns A [tibble()]. Values are described in Details.
#' @export

import_microbigge_gcp <- function(path){
  .complete <- readr::read_tsv(file = path, show_col_types = FALSE) %>%
    dplyr::select(protein = "protein_acc", biosample = "biosample_acc", assembly = "asm_acc",
           allele = "element_symbol", name = "element_name", method = "amr_method",
           coverage = "pct_ref_coverage", identity = "pct_ref_identity",
           nucleotide = "contig_acc", nuc_start = "start_on_contig",
           nuc_stop = "end_on_contig", strand = "strand"
    )

  return(.complete)
}

#' Filter imported MicroBigg-E data
#'
#' @description Removes allele calls that do not meet the minimum `coverage` and
#'   `identity` thresholds. These values may vary based on the needs of a
#'   particular analysis. Set both to `100` to allow only exact allele calls.
#'
#'   Specific allele calling methods can be excluded from the dataset by name in
#'   `remove`. See <<>> for more information about allele calling methods used
#'   in MicroBIGG-E data.
#'
#' @details Isolates Browser (and thus amr.metadata.tsv files on the NCBI FTP
#'   server) uses a 90% `identity` threshold:
#'   https://www.ncbi.nlm.nih.gov/pathogens/pathogens_help/#genotype-categories
#'
#' @param data A dataframe or tibble
#' @param coverage Minimum percent coverage to keep
#' @param identity Minimum percent identity to keep
#' @param remove Character vector of `method` values to remove alleles
#' @returns `data` with rows removed according to the selected parameters
#' @export
filter_microbigge <- function(data, coverage = 100L, identity = 90L, remove) {
  . <- NULL # Workaround to suppress `no visible binding for global variable`
  remove_missing <- missing(remove)
  .complete <- data %>%
    tidyr::drop_na("protein") %>%
    {if(!remove_missing) dplyr::filter(., !grepl(paste(remove, collapse = "|"), .data$method)) else .} %>%
    dplyr::filter("coverage" >= coverage) %>%
    dplyr::filter("identity" >= identity)

  return(.complete)
}

#' Parse `oxa_family` from MBE & IPG names
#' @rdname parse_ib_oxa_family
#' @export

parse_mbe_oxa_family <- function(data) {
  .complete <- data %>%
    dplyr::mutate(mbe_oxa_family = dplyr::case_when(
      grepl('family', name) ~ stringr::str_extract(name, "OXA-\\d+?(?= family)"),
      TRUE ~ stringr::str_extract(name, "OXA-\\d+")
    )) %>%
    dplyr::mutate(ipg_oxa_family = dplyr::case_when(
      grepl('family', ipg_name) ~ stringr::str_extract(ipg_name, "OXA-\\d+?(?= family)"),
      TRUE ~ stringr::str_extract(ipg_name, "OXA-\\d+")
    )) %>%
    dplyr::mutate(oxa_family = dplyr::coalesce(.data$ipg_oxa_family, .data$mbe_oxa_family)) %>%
    dplyr::mutate(oxa_family_markdown = sub("^bla(.*?)$", "<i>bla</i><sub>\\1</sub>",
                                            .data$oxa_family, perl = TRUE) ) %>%
    dplyr::mutate(oxa_family_math = sub("^bla(.*?)$", "italic\\(bla\\)\\[\\1\\]",
                                        .data$oxa_family, perl = TRUE) )

  return(.complete)
}

########## Identical Protein Group Import and Processing ##########

#' Import Identical Protein Groups data
#'
#' @description Imports Identical Protein Groups data from `path` and adds
#' columns corresponding to the IPG UID (`ipg`), IPG Accession Number
#' (`ipg_accession`), and name of the IPG protein (`iph_name`) to `data`.
#'
#' @param data A dataframe or tibble
#' @param path path to an NCBI IPG tsv file
#' @returns `data` with columns `ipg`, `ipg_accession`, and `ipg_name` added
#' @export
import_ipg <- function(data, path){
  if(is.null(path)) return(data)

  .ipg <- readr::read_tsv(file = path, show_col_types = FALSE) %>%
    dplyr::distinct(.data$protein, .keep_all = TRUE)

  .complete <- data %>%
    dplyr::left_join(.ipg, by = "protein")

  return(.complete)
}
