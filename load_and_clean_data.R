library(tidyverse)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(SomaticSignatures)

file1 <- "data/PD11461_dirichlet_subs_plus_timing_v0.1.txt"
file2 <- "data/PD8948_dirichlet_subs_plus_timing_v0.1.txt"


col_spec <- cols(
  ID = col_character(),
  Chrom = col_character(),
  Pos = col_double(),
  Ref = col_character(),
  Alt = col_character(),
  Gene = col_character(),
  Transcript = col_character(),
  RNA = col_character(),
  CDS = col_character(),
  Protein = col_character(),
  Type = col_character(),
  Effect = col_character(),
  SNP = col_character(),
  VD = col_character(),
  VW = col_character(),
  dis_depth_PD11461a = col_double(),
  dis_MtAll_PD11461a = col_double(),
  dis_MtAllPct_PD11461a = col_double(),
  dis_pbinom_PD11461a = col_double(),
  dis_CLASS_PD11461a = col_character(),
  val_depth_PD11461a = col_double(),
  val_MtAll_PD11461a = col_double(),
  val_MtAllPct_PD11461a = col_double(),
  val_pbinom_PD11461a = col_double(),
  val_CLASS_PD11461a = col_character(),
  concordance_PD11461a = col_character(),
  dis_depth_PD11461c = col_double(),
  dis_MtAll_PD11461c = col_double(),
  dis_MtAllPct_PD11461c = col_double(),
  dis_pbinom_PD11461c = col_double(),
  dis_CLASS_PD11461c = col_character(),
  val_depth_PD11461c = col_double(),
  val_MtAll_PD11461c = col_double(),
  val_MtAllPct_PD11461c = col_double(),
  val_pbinom_PD11461c = col_double(),
  val_CLASS_PD11461c = col_character(),
  concordance_PD11461c = col_character(),
  dis_depth_PD11461b = col_double(),
  dis_MtAll_PD11461b = col_double(),
  dis_MtAllPct_PD11461b = col_double(),
  dis_pbinom_PD11461b = col_double(),
  dis_CLASS_PD11461b = col_character(),
  val_depth_PD11461b = col_double(),
  val_MtAll_PD11461b = col_double(),
  val_MtAllPct_PD11461b = col_double(),
  val_pbinom_PD11461b = col_double(),
  val_CLASS_PD11461b = col_character(),
  concordance_PD11461b = col_character(),
  manual_exclusion_preVal = col_character(),
  validation_call = col_character(),
  validation_platform = col_character(),
  INCLUDE_POST_VAL = col_character(),
  visual_inspection = col_character(),
  in.discovery = col_character(),
  in.validation = col_character(),
  CHR = col_character(),
  POS = col_character(),
  dis_MtAllPct_PD11461a.1 = col_double(),
  dis_MtAllPct_PD11461c.1 = col_double(),
  dis_depth_PD11461a.1 = col_double(),
  dis_depth_PD11461c.1 = col_double(),
  Gene.1 = col_character(),
  subclonal.CN.a = col_double(),
  subclonal.CN.c = col_double(),
  nMaj1.a = col_double(),
  nMaj1.c = col_double(),
  nMin1.a = col_double(),
  nMin1.c = col_double(),
  frac1.a = col_double(),
  frac1.c = col_double(),
  nMaj2.a = col_character(),
  nMaj2.c = col_double(),
  nMin2.a = col_character(),
  nMin2.c = col_double(),
  frac2.a = col_character(),
  frac2.c = col_double(),
  mutation.copy.number.a = col_double(),
  subclonal.fraction.a = col_double(),
  no.chrs.bearing.mut.a = col_double(),
  mutation.copy.number.c = col_double(),
  subclonal.fraction.c = col_double(),
  no.chrs.bearing.mut.c = col_double(),
  CHR.1 = col_character(),
  POS.1 = col_character(),
  val_MtAllPct_PD11461a.1 = col_double(),
  val_MtAllPct_PD11461c.1 = col_double(),
  val_depth_PD11461a.1 = col_double(),
  val_depth_PD11461c.1 = col_double(),
  Gene.2 = col_character(),
  subclonal.CN.a.1 = col_double(),
  subclonal.CN.c.1 = col_double(),
  nMaj1.a.1 = col_double(),
  nMaj1.c.1 = col_double(),
  nMin1.a.1 = col_double(),
  nMin1.c.1 = col_double(),
  frac1.a.1 = col_double(),
  frac1.c.1 = col_double(),
  nMaj2.a.1 = col_character(),
  nMaj2.c.1 = col_character(),
  nMin2.a.1 = col_character(),
  nMin2.c.1 = col_character(),
  frac2.a.1 = col_character(),
  frac2.c.1 = col_character(),
  mutation.copy.number.a.1 = col_double(),
  subclonal.fraction.a.1 = col_double(),
  no.chrs.bearing.mut.a.1 = col_double(),
  mutation.copy.number.c.1 = col_double(),
  subclonal.fraction.c.1 = col_double(),
  no.chrs.bearing.mut.c.1 = col_double(),
  discovery.DP.cluster = col_double(),
  validation.DP.cluster = col_double(),
  samples = col_character(),
  chr = col_character(),
  start = col_double(),
  stop = col_double(),
  gene = col_character(),
  Driver.confidence = col_character(),
  TYPE = col_character(),
  change = col_character(),
  timing = col_character()
)
fr1 <- read_tsv(file1, col_types = col_spec)
fr1_primary <- fr1 %>% dplyr::select(ID,Chrom,Pos,Ref,Alt,ends_with("PD11461c")) %>%
  dplyr::mutate(freq = dis_MtAllPct_PD11461c/100.)

# Generate VRanges object for SomaticSignatures
reference <- BSgenome.Hsapiens.1000genomes.hs37d5
vr <- VariantAnnotation::VRanges(seqnames = fr1_primary$Chrom,
                                 ranges = IRanges(start = fr1_primary$Pos, width = rep(1,nrow(fr1_primary))),
                                 ref = fr1_primary$Ref,
                                 alt = fr1_primary$Alt,
                                 freq = fr1_primary$freq)
#VariantAnnotation::refDepth(vr) <- fr$t_ref_count
#VariantAnnotation::altDepth(vr) <- fr$t_alt_count
GenomeInfoDb::seqlevels(vr) <- GenomicAlignments::seqlevelsInUse(vr)
GenomeInfoDb::genome(vr) <- GenomeInfoDb::genome(reference)[1:length(GenomeInfoDb::genome(vr))]
vr <- SomaticSignatures::mutationContext(vr,reference)
vars <- tibble::as_tibble(vr)
vars <- vars %>%  dplyr::select(seqnames, start, end, freq, alteration, context)
vars$trinucleotide_context <- paste0(vars$alteration,"_",vars$context)
vars
vars %>% readr::write_tsv("./PD11461_primary.tsv")


# Note that the "b" sample here is the normal, "a" is the local relapse, and "c" is the primary
fr1_relapse <- fr1 %>% dplyr::select(ID,Chrom,Pos,Ref,Alt,ends_with("PD11461a")) %>%
  dplyr::mutate(freq = dis_MtAllPct_PD11461a/100.)

# Generate VRanges object for SomaticSignatures
reference <- BSgenome.Hsapiens.1000genomes.hs37d5
vr <- VariantAnnotation::VRanges(seqnames = fr1_relapse$Chrom,
                                 ranges = IRanges(start = fr1_relapse$Pos, width = rep(1,nrow(fr1_relapse))),
                                 ref = fr1_relapse$Ref,
                                 alt = fr1_relapse$Alt,
                                 freq = fr1_relapse$freq)
#VariantAnnotation::refDepth(vr) <- fr$t_ref_count
#VariantAnnotation::altDepth(vr) <- fr$t_alt_count
GenomeInfoDb::seqlevels(vr) <- GenomicAlignments::seqlevelsInUse(vr)
GenomeInfoDb::genome(vr) <- GenomeInfoDb::genome(reference)[1:length(GenomeInfoDb::genome(vr))]
vr <- SomaticSignatures::mutationContext(vr,reference)
vars <- tibble::as_tibble(vr)
vars <- vars %>%  dplyr::select(seqnames, start, end, freq, alteration, context)
vars$trinucleotide_context <- paste0(vars$alteration,"_",vars$context)
vars
vars %>% readr::write_tsv("./PD11461_relapse.tsv")
