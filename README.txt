                       README for Taxon Sampling
             (last updated 08/01/2020 - mm/dd/yyyy format)

AUTHORS
-=-=-=-

Jorge Augusto Hongo (jorgeahongo@gmail.com)
Alison Menezes (alisonpam@gmail.com)
Maycon Oliveira (flayner5@gmail.com)
Francisco Pereira Lobo (franciscolobo@gmail.com, francisco.lobo@embrapa.br)


1 - DESCRIPTION
-=-=-=-=-=-=-=-

Taxon Sampling takes as input NCBI Taxonomy and a database of taxon IDs and
their associated sequences and samples a given amount of sequence records 
while maximizing taxonomic diversity by attempting to cover the most taxa
in each taxonomic level.

2 - HOW TO USE
-=-=-=-=-=-=-=-

2.1 - INSTALL
-=-=-=-=-=-=-

Taxon Sampling requires the CHNOSZ, ape and seqinr R packages. To install them:

  install.packages("CHNOSZ")
  install.packages("ape")
  install.packages("seqinr")

Additionally, Taxon Sampling requires the NCBI Taxonomy files. To obtain them,
execute install.sh (with writing permissions in TaxonSampling root directory).


2.2 - PREPARING INPUT FILES 
-=-=-=-=-=-=-=-=-=-=-=-=-=-

To run Taxon Sampling, prepare the following files:

  - multifasta file with the input sequences to sample
  - tabular file with the taxonomy IDs (NCBI) and the sequence ID of the
    multifasta file, with the columns in this order

Example of the tabular file:

-------------------------------------------------------------------------------
100182	gi|597501527|gb|KF917426.1| Lepus granatensis
1001917	gi|666686397|gb|KM031740.1| Microphysogobio alticorpus voucher
-------------------------------------------------------------------------------

and a corresponding multifasta file:

-------------------------------------------------------------------------------
>gi|597501527|gb|KF917426.1| Lepus granatensis
GTACCTTGTCACTATTGACAAAAATTCCCCTTAACGCTATGTAATTCGTGCATTAGTGCCTTTCCCCATC
AACATGTATCCATACTATCATTTCATAATCAACATTAGACCATTCTATGTTTAATCGTACATTAAAATCT
TATCCCCATGCATATAAGCCAGTACATCCCTACTTACAGGACATAGCACATTCCCTGCAAACTCACAAAC
CCCTATCATCAACACGAATATCCACAACCCAATTACCCACCTTAATCAACATCCAAACATTCATTCCTTA
ATTGGTCATAGACCATCCAAGTCAAATCTTTTCTCGTCCATATGACTATCCCCTCCCCTTAGGAATCCCT
TGATCTACCATCCTCCGTGAAACCAGCAACCCGCCCACCCCGTATCCCTCTTCTCGCTCCGGGCCCATTA
CACTTGGGGGTTTCTAGCGTGAAACT

>gi|666686397|gb|KM031740.1| Microphysogobio alticorpus voucher
ATAGTAGCATGTATGTACTAGTACATATTATGCATAATATTACATAATGTATTAGTACATATATGTATTA
TCACCATTCATTTATTTTAACCCCAAAGCAAGTACTAACGTCTAAGAAGACCATAAACCAAATTATTAAT
ATTCAAAATTAATTTATTTTAACATGATAATATATATTTATCCACTAAAATTTATTTTAATATATATTCT
AGTGAGAGACCACCAACTGGTTGATATAATTGCATATCATGAATGATAGAATCAGGGACAAAACACTAAG
AGACGGTATATTATGAATTATTCCTTGTATCTGGTTATCAATCTCACGTACTTTACTTATTTACTCCACA
CACTGTCATTGAATTGGCATATGGTTAAATGATGTGAGTACATACTCCTCGTTACCCAACATGCCGAGCG
TTCTTTTATATGCATAGGGTTCTCCTTTTTTGGTTTCCTTTCATCATACATCTCAGAGTGCAGGCACAAA
TAATATATCAAGGTGGAACTATTCCTTGCACGAATTAACGTAGGTTAAATTATTAAATGACATAACTCAA
GAATTACATATTAATATCTCAGGTGCATAAGGTATATATTACTTCTTCAACTTACCCCTATATATATGCC
CCCCCTTTTGGCTTACGCGCGACAAACCCCCCTACCCCCTTCGCTCAGCGAATCCTGTTATCCTTGTCAA
ACCC
-------------------------------------------------------------------------------



2.3 - SETTING UP THE PARAMETERS 
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

Set the following parameters in the beginning of the runTS.R file to analyze
your data:

  idsFile: path to the tabular file

  multifasta: path to multifasta file

  taxondir: path to the NCBI taxonomy directory (default is "taxondump/")
  inside root TaxonSampling directory.

  outFile: name of the output file

and the parameters for the sampling algorithm:

  m: number of ID samples to obtain

  method: either "diversity" or "balance". Diversity will try to sample
          sequences as diverse as possible across the phylogeny, while
          balance will try to sample the same amount of sequences fore all
          child nodes. Please notice that balance mode may sample fewer
          than "m" parameter.

  randomize: optional randomization when sampling within a branch

  replacement: optional replacement when sampling within a branch (useful to
  sample more sequences in balance mode in the case where a node contains too
  few sequences, therefore constraining the number of sequences in nodes with
  
  ignoreIDs: which IDs, if any, should be ignored. Ignoring non-leaf nodes here 
             will also ignore its children/offspring

  requireIDs: which terminal IDs, if any, must appear in the result

  ignoreNonLeafID: which IDs, if any, should be ignored. Ignoring non-leaf 
                   nodes here won't prevent its children/offspring from being
                   sampled


The allowed values for the algorithm are:

-------------------------------------------------------------------------------
  m - any positive integer

  method - "diversity" or "balance"

  randomize - "yes" or "no"

  replacement - "yes" or "no"

  ignoreIDs - NULL, or a vector of NCBI taxon IDs for species. Ex: c("1001947",
  "10029")

  requireIDs - NULL, or a vector of NCBI taxon IDs for species. Ex: c("1001947",
  "10029")

  ignoreNonLeafID - NULL, or a vector of IDs. Ex: c("1001947", "10029")

  outFile - any valid character string

-------------------------------------------------------------------------------


2.4 - RUNNING Taxon Sampling
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

After setting the parameters and input files, run:

  source("runTS.R")

Your results should be printed in the path specified by outFile. Note that R 
is case-sensitive and requires quote marks (' or ") to treat your input as a 
character, rather than an existing variable.
