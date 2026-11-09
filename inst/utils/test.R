spp_file      <- "data_files/metadata/TaxID2sppCounts.tsv"
taxonomy_path <- "data_files/taxdump/"
ids_file      <- "data_files/metadata/TaxID2SeqID.txt"
seq_file      <- "data_files/fasta/mit_vertebrata.fasta"

## Run once to save species counts to file spp_file
# t0 <- Sys.time()
# spp_df <- get_taxID_spp_counts(taxonomy_path = taxonomy_path,
#                                start_from_species = FALSE,
#                                out_file = NULL,
#                                what = "spp_counts")
# dt <- Sys.time() - t0
# message("\nElapsed time: ", dt, " ", units(dt))

counts <- get_counts(taxonomy_path = taxonomy_path,
                     ids_file      = ids_file,
                     # spp_df        = spp_df,
                     spp_file      = spp_file)
ts_out <- run_TS(taxlist          = counts,
                 taxon            = 40674,
                 m                = 500,
                 seq_file         = seq_file,
                 out_file         = NULL,
                 method           = "diversity",
                 randomize        = "after_first_round",
                 replacement      = FALSE,
                 ignoreIDs        = NULL,
                 requireIDs       = NULL,
                 sampling         = "agnostic")

# summary(ts_out)
# plot(ts_out)
# extract_taxlevel(taxIDs = ts_out$outputIDs,
#                  taxlevel = "order",
#                  nodes = ts_out$nodes)
