#!/usr/bin/env Rscript

# load packages 
library(Rsubread)
library(tools)
library(readxl)
library(readr)

args=commandArgs(trailingOnly=TRUE)


# Determine counts 
mywtcounts <- featureCounts(args[1],
  annot.ext = args[2] , 
  isGTFAnnotationFile = TRUE, 
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  isPairedEnd = TRUE,
  strandSpecific = 2)

# 
Counts <- as.data.frame(mywtcounts$counts)

mybam=args[1]
myprefix=basename(mybam)
my_counts_pref <- file_path_sans_ext(myprefix)


write_tsv(Counts, file=paste0(my_counts_pref, ".tsv"))