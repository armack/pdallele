library(dplyr)

# Prepare raw files for `extdata` -----------------------------------------

## Prepare `amr.metadata.tsv` -----
metadata <- readr::read_tsv(file.path(getwd(), "data-raw/amr.metadata.tsv")) %>%
  slice_sample(n = 50)

metadata_biosamples <- metadata %>%
  distinct(biosample_acc) %>%
  pull(biosample_acc)

readr::write_tsv(x = metadata, file = file.path(getwd(), "inst/extdata/amr.metadata.tsv"))

## Prepare `microbigge.tsv` -----
microbigge <- readr::read_tsv(file.path(getwd(), "data-raw/microbigge.tsv")) %>%
  filter(biosample_acc %in% metadata_biosamples)

microbigge_element_symbols <- microbigge %>%
  distinct(element_symbol) %>%
  pull(element_symbol)

microbigge_proteins <- microbigge %>%
  distinct(protein_acc) %>%
  pull(protein_acc)

readr::write_tsv(x = microbigge, file = file.path(getwd(), "inst/extdata/microbigge.tsv"))

## Prepare `cluster_list.tsv` -----
cluster_list <- readr::read_tsv(file.path(getwd(), "data-raw/cluster_list.tsv")) %>%
  filter(biosample_acc %in% metadata_biosamples)

readr::write_tsv(x = cluster_list, file = file.path(getwd(), "inst/extdata/cluster_list.tsv"))

## Prepare `refgene.txt` -----
refgene <- readr::read_tsv(file.path(getwd(), "data-raw/refgene.tsv")) %>%
  filter(allele %in% microbigge_element_symbols)

readr::write_tsv(x = refgene, file = file.path(getwd(), "inst/extdata/refgene.tsv"))

## Prepare `ipg.tsv` -----
ipg <- readr::read_tsv(file.path(getwd(), "data-raw/ipg.tsv")) %>%
  filter(protein %in% microbigge_proteins)

readr::write_tsv(x = ipg, file = file.path(getwd(), "inst/extdata/ipg.tsv"))

# Prepare sample datasets -------------------------------------------------

## File Paths to "extdata"
isolates_path <- system.file("extdata/amr.metadata.tsv", package = "pdallele")
cluster_path <- system.file("extdata/cluster_list.tsv", package = "pdallele")
refgene_path <- system.file("extdata/refgene.tsv", package = "pdallele")
microbigge_path <- system.file("extdata/microbigge.tsv", package = "pdallele")
ipg_path <- system.file("extdata/ipg.tsv", package = "pdallele")
regions_path <- system.file("extdata/regions.csv", package = "pdallele")

## Isolates Browser -----
isolates_raw <- import_isolates_browser_metadata(isolates_path)

isolates_unfiltered <- isolates_raw %>%
  na_if_tibble_chr(terms = na_strings) %>%
  parse_genus_species() %>%
  reverse_geocode() %>%
  split_location() %>%
  import_regions(path = regions_path) %>%
  import_cluster_list(cluster_path) %>%
  separate_genotypes(include = "amr") %>%
  add_reference_gene_catalog(refgene_path) %>%
  parse_bla_formatting() %>%
  parse_ib_oxa_family() %>%
  parse_year()

isolates <- isolates_unfiltered %>%
  clean_filter_alleles(filter = "bla", remove = remove_allele_types)

## Microbial Browser for Identification of Genetic and Genomic Elements (MicroBIGG-E) -----
mbe_raw <- import_microbigge_gcp(microbigge_path)

mbe <- mbe_raw %>%
  na_if_tibble_chr(terms = na_strings) %>%
  filter_microbigge(coverage = 100, identity = 100, remove = remove_mbe_allele_types) %>%
  import_ipg(ipg_path) %>%
  import_cluster_list(cluster_path) %>%
  #import_mlst(mlst_path) %>%
  parse_mbe_oxa_family()
