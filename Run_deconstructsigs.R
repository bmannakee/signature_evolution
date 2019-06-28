library(fs)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg38)
library(deconstructSigs)

test_file <- fs::fs_path('./PD8948e_primary2.bed')
col_spec <- readr::cols('c','i','i','d','c','c')
col_names <- c("seqnames","start","end","vaf","ref","alt")


sample_name<-"Test2"

this_size <- 10000

fr <- readr::read_tsv(test_file,col_types = col_spec, col_names = col_names)

# Add times here, but name the column "Sample" rather than "time"
order_sample <- function(df, size = this_size){
  max_sample <- ceiling(nrow(df)/size)
  df$Sample <- NA
  for(i in 1:nrow(df)){
    df$Sample[i]<-ceiling(i/size)
  }
  df
}

fr <- order_sample(fr)
# code to drop any sample with less than 50 variant
fr [, -which(colMeans(is.na(fr)) > 50)]
fr <- 
fr <- fr %>% dplyr::select(Sample,"chr"="seqnames","pos"="start",ref,alt)
ds_input <- deconstructSigs::mut.to.sigs.input(as.data.frame(fr), sample.id = "Sample",bsg = BSgenome.Hsapiens.1000genomes.hs37d5)
#for loop with sample ID 
all_sigs <- list()
for( i in unique(fr$Sample)){ 
  print(paste("sample.id=", i))
  sigs <- deconstructSigs::whichSignatures(as.data.frame(ds_input),
                                           contexts.needed = T,
                                           sample.id = i,
                                           signatures.ref=deconstructSigs::signatures.cosmic,
                                           signature.cutoff = 0)
  all_sigs[[i]] <- as_tibble(sigs$weights)
  } 
final_fr <- bind_rows(all_sigs)
final_fr %>% readr::write_tsv('./data/PD8948e_primary2_deconstructsigs_binsize10000.tsv')
