# run the program
library(fs)
library(tidyverse)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(BSgenome.Hsapiens.UCSC.hg38)
library(sigtracker)




get_local_signature_fr <- function(sig_file){
  fr <- suppressMessages(read_tsv(sig_file) %>% dplyr::select(-contains("X")))
  names(fr) %<>% stringr::str_replace_all("\\s","_") %>% tolower
  # get the contexts in the correct form
  fr <- fr %>% dplyr::mutate(alteration = str_replace(substitution_type,">",""),
                             context = paste0(str_sub(trinucleotide,1L,1L),".",str_sub(trinucleotide,3L,3L)),
                             trinucleotide_context = paste0(alteration,"_",context))
  fr <- fr %>% dplyr::select(trinucleotide_context,contains("signature_"))
  fr
}



test_file <- fs::fs_path('~/Desktop/projects/signature_evolution/PD11461_primary.bed')
sig_file <- fs::fs_path("./signatures.txt")
sig_fr <- get_local_signature_fr(sig_file) 
this_sig_fr <- sig_fr %>% dplyr::select(-trinucleotide_context)
fr <- run_tracksig(test_file,BSgenome.Hsapiens.1000genomes.hs37d5,bin_size=200,file_type="bed")

tracksig_matrix <- as.matrix(fr %>% dplyr::select(-time))
sig_matrix <- t(as.matrix(sig_fr %>% dplyr::select(colnames(fr %>% dplyr::select(-time)))))
final_signatures <- tracksig_matrix %*% sig_matrix
colnames(final_signatures) <- sig_fr$trinucleotide_context
final_signatures_fr <- as_tibble(final_signatures)
final_fr <- fr %>% bind_cols(final_signatures_fr)
final_fr %>% readr::write_tsv('./data/pd11461_primary_tracksig_binsize200.tsv')
