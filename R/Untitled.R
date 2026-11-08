library(GenomicRanges)
library(Biostrings)
library(GenomicFeatures)
#library(MotifDb)
#library(seqLogo)
#library(JASPAR2020)
#library(TFBSTools)
#library(motifmatchr)
#library(seqinr)

#working dir
setwd("~/projects/CALANGO/genome_modeling/")

rm(list=ls())

#path to fasta and gff files
fasta_file <- "GCF_029289425.1_NHGRI_mPanPan1-v1.1-0.1.freeze_pri_genomic.fna"
gff_file <- "GCF_029289425.1_NHGRI_mPanPan1-v1.1-0.1.freeze_pri_genomic.gff"

#fasta_file <- "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
#gff_file <- "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff"

#read fasta
genome <- readBStringSet(filepath=fasta_file, format="fasta")

#sanitizing fasta headers
names <- names(genome)
names <- sapply(strsplit(names, " "), "[[", 1)
names <- gsub(",", "", names)

names(genome) <- names

#write.fasta(sequences = genome, names = names(genome), file.out = "tmp.fna")
#importing txdb
txdb <- makeTxDbFromGFF(file=gff_file)

#creating GRanges objects

#transcripts
all_transcripts <- exonsBy(txdb, by="tx", use.names=TRUE)

#genes
all_genes <- cdsBy(txdb, by="gene")

#cds
all_cds <- cdsBy(txdb, by="tx", use.names=TRUE)

tx_length <- data.frame(extractTranscriptSeqs(genome, all_cds)@ranges@NAMES, extractTranscriptSeqs(genome, all_cds)@ranges@width)
colnames(tx_length) <- c("TXNAME", "WIDTH")

gene_db <- select(txdb, keys = keys(txdb), columns=c("TXNAME", "CDSID", "CDSNAME", "CDSSTRAND"), keytype="GENEID")
gene_db <- unique(gene_db[,c("GENEID", "CDSNAME", "TXNAME", "CDSSTRAND")])
gene_db <- gene_db[complete.cases(gene_db),]
gene_db <- merge(gene_db, tx_length, by="TXNAME")
#sort isoforms per gene by length
gene_db <- gene_db[order(gene_db$GENEID, -abs(gene_db$WIDTH)), ]
# Keep isoform with the longest CDS
gene_db <- gene_db[ !duplicated(gene_db$GENEID), ]

#getting genomic coordinates for longest transcripts
longest_transcripts_db <- select(txdb, keys = gene_db$TXNAME, keytype = "TXNAME", columns = c("TXNAME", "TXCHROM", "TXSTART", "TXEND", "TXSTRAND", "GENEID"))

#getting genomic coordinates for the coding exons of longest transcripts
coding_exons_longest_transcripts_db <- select(txdb, keys = gene_db$GENEID, keytype = "GENEID", columns = c("TXID", "TXNAME", "CDSNAME", "CDSCHROM", "CDSSTART", "CDSEND", "CDSSTRAND", "GENEID"))

#subsetting by longest transcript IDs
coding_exons_longest_transcripts_db <- coding_exons_longest_transcripts_db[coding_exons_longest_transcripts_db$TXNAME %in% longest_transcripts_db$TXNAME, ]
coding_exons_longest_transcripts_db <- coding_exons_longest_transcripts_db[complete.cases(coding_exons_longest_transcripts_db),]

tmp_db <- gene_db
tmp_db$STARTCODONBEGIN <- NA
#tmp_db$STARTCODONEND <- NA
i = 0
for (ID in unique(coding_exons_longest_transcripts_db$GENEID)) {
  print(paste(i, "-", ID))
  strand <- unique(coding_exons_longest_transcripts_db$CDSSTRAND[coding_exons_longest_transcripts_db$GENEID == ID])
  strand <- strand[complete.cases(strand)]
  if (strand == "-") { #see A4GALT/XP_003822530.1
    start_codon_begin <- max(coding_exons_longest_transcripts_db$CDSEND[coding_exons_longest_transcripts_db$GENEID == ID])
    #    start_codon_begin <- max(coding_exons_longest_transcripts_db$CDSEND[coding_exons_longest_transcripts_db$GENEID == ID]) + 1
    #    start_codon_end <- start_codon_begin - 3
    tmp_db$STARTCODONBEGIN[tmp_db$GENEID == ID] <- start_codon_begin
    #    tmp_db$STARTCODONEND[tmp_db$GENEID == ID] <- start_codon_end
    #    print("Minus")
  }
  if (strand == "+") { #see A2ML1/XP_008959086.4
    start_codon_begin <- min(coding_exons_longest_transcripts_db$CDSSTART[coding_exons_longest_transcripts_db$GENEID == ID])
    #    start_codon_begin <- min(coding_exons_longest_transcripts_db$CDSSTART[coding_exons_longest_transcripts_db$GENEID == ID]) - 1
    #    start_codon_end <- start_codon_begin + 3
    tmp_db$STARTCODONBEGIN[tmp_db$GENEID == ID] <- start_codon_begin
    #    tmp_db$STARTCODONEND[tmp_db$GENEID == ID] <- start_codon_end
    #    print("Plus")
  }
  i <- i + 1
}


#promoters
#all_promoters <- promoters(txdb, upstream=2000, downstream=50, use.names=TRUE)

#if beginning of promoter is negative, set it to one
#start(all_promoters) <- pmax(start(all_promoters), 1)

#TODO - need to fix the same for end of promoters as well
#end(promoters) <- pmax(end(promoters), 1)

fiveUTR <- fiveUTRsByTranscript(txdb, use.names=TRUE)

#longest_cds <- all_cds[names(all_cds) %in% coding_exons_longest_transcripts$TXNAME,]
#longest_cds_seqs <- extractTranscriptSeqs(genome, longest_cds)

#longest_transcripts <- all_transcripts[names(all_transcripts) %in% longest_transcripts_db$TXNAME,]
#longest_transcript_seqs <- extractTranscriptSeqs(genome, longest_transcripts)

promoter_longest_transcripts_db <- select(txdb, keys = gene_db$TXNAME, keytype = "TXNAME", columns = c("TXNAME", "TXCHROM", "TXSTART", "TXEND", "TXSTRAND", "GENEID"))
promoter_longest_transcripts_db <- promoter_longest_transcripts_db[promoter_longest_transcripts_db$TXNAME %in% longest_transcripts_db$TXNAME, ]
promoter_longest_transcripts <- all_promoters[names(all_promoters) %in% coding_exons_longest_transcripts$TXNAME,]
promoter_longest_tx_seqs <- getSeq(genome, promoter_longest_transcripts[1:(length(promoter_longest_transcripts)-3)])


opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "insects"
motifList <- getMatrixSet(JASPAR2020, opts)

seqLogo::seqLogo(motifList[[1]], ic.scale = FALSE)

start_time <- Sys.time()

motif_positions <- matchMotifs(motifList, promoter_seqs[1:35600], out = "positions", p.cutoff = 0.01)

end_time <- Sys.time()







#sequence <- "CC"


dinuc_count <- function (seq) {

  seq <- seq

  nucs <- c("A", "C", "T", "G")

  #creating data structures
  dinuc_count <- data.frame(matrix(ncol = 16, nrow = 1))
  nuc_count <- data.frame(matrix(ncol = 4, nrow = 1))
  colnames_dinuc <- vector(length = 16)
  colnames_nuc <- vector(length = 4)

  i = 1

  for (i2 in 1:length(nucs)) {
    colnames_nuc[i2] <- nucs[i2]
    for (i3 in 1:length(nucs)) {
      dinuc <- paste0(nucs[i2], nucs[i3])
      colnames_dinuc[i] <- dinuc
      i <- i + 1
    }
  }

  #adding col names
  colnames(dinuc_count) <- colnames_dinuc
  colnames(nuc_count) <- colnames_nuc

  #populating data
  for (i in colnames(dinuc_count)) {
    dinuc_count[i] <- 0
  }
  for (i in colnames(nuc_count)) {
    nuc_count[i] <- 0
  }

  #splitting data
  vec <- unlist(strsplit(seq, split = ""))

  #counting nts
  for (i in 1:(length(vec))) {
    nuc_count[vec[i]] <- nuc_count[vec[i]] + 1
  #  print(paste(i, dinuc, vec[i], vec[i+1]))
  }

  #counting dinucs
  for (i in 2:(length(vec))-1) {
    dinuc <- paste0(vec[i], vec[i+1])
    dinuc_count[dinuc] <- dinuc_count[dinuc] + 1
#  print(paste(i, dinuc, vec[i], vec[i+1]))
  }

  return(dinuc_count)

}


dinuc_odds_ratio <- function (dinuc_count, nuc_count, strandedness) {

}
