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

## Prepare `refgene.txt` -----
refgene <- readr::read_tsv(file.path(getwd(), "data-raw/refgene.tsv")) %>%
  filter(allele %in% microbigge_element_symbols)

readr::write_tsv(x = refgene, file = file.path(getwd(), "inst/extdata/refgene.tsv"))

## Prepare `ipg.tsv` -----
ipg <- readr::read_tsv(file.path(getwd(), "data-raw/ipg.tsv")) %>%
  filter(protein %in% microbigge_proteins)

readr::write_tsv(x = ipg, file = file.path(getwd(), "inst/extdata/ipg.tsv"))

# Prepare sample datasets -------------------------------------------------
