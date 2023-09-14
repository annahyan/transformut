library(here)
library(tidyverse)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SomaticSigantures)
library(MutationalPatterns)

variants.dir = here("../mutsigints/data/raw/PCAWG/vcfs/consensus_snv_indel")

snv.dir = file.path(variants.dir, "snv_mnv")
indel.dir = file.path(variants.dir, "indel")

output.dir = here("mutation_sequences")
if (!file.exists(output.dir)) {
  dir.create(output.dir)
}

snv.files = list.files(snv.dir)
indel.files = list.files(indel.dir)

get_sample_name = function(x) {
  gsub("(.*).consensus.*", "\\1", x )
}

snv.samples = sapply(snv.files, get_sample_name)
indel.samples = sapply(indel.files, get_sample_name)

common.samples = intersect(snv.samples, indel.samples)

soi = common.samples[1:300]

snv.soi = grep("gz$", names(snv.samples[which(snv.samples %in% soi)]), value = TRUE)
indel.soi = grep("gz$", names(indel.samples[which(indel.samples %in% soi)]), value = TRUE)

snv.granges = read_vcfs_as_granges(file.path(snv.dir, snv.soi), snv.samples[snv.soi],
                                   "BSgenome.Hsapiens.UCSC.hg19")

indel.granges = read_vcfs_as_granges(file.path(indel.dir, indel.soi), indel.samples[indel.soi],
                                   "BSgenome.Hsapiens.UCSC.hg19", type = "indel")


for (sname in soi) {
  sname.out = file.path(output.dir, paste0(sname, ".sequence.mut"))

  snv.mut.context = mut_context(snv.granges[[sname]],
                                  "BSgenome.Hsapiens.UCSC.hg19", extension = 15)
  indel.mut.context = mut_context(indel.granges[[sname]],
                                "BSgenome.Hsapiens.UCSC.hg19", extension = 15)

  write.table(c(snv.mut.context, indel.mut.context), sname.out,
              quote = F, row.names = F, col.names = F, append = TRUE)
}
