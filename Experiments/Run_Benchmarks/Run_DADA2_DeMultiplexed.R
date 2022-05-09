.libPaths("/fs/cbcb-scratch/tluan/dada2")
library(dada2); packageVersion("dada2")

path <-"/fs/cbcb-lab/mpop/projects/SCRAPT/Datasets/Lupus-Microbiome-Published/combined_samples"
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
sample.names <-sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1)

filtFs <- file.path(path, "filtered_combined", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

out <- filterAndTrim(fwd = fnFs, filt = filtFs,  maxN=0, maxEE=Inf, truncQ=2, rm.phix=FALSE,compress=TRUE, multithread=8) 
errF <- learnErrors(filtFs, multithread=8)
dadaFs <- dada(filtFs, err=errF, multithread=8)
seqtab <- makeSequenceTable(dadaFs)

options(max.print=214748364)
dim(seqtab)

capture.output(seqtab, file = "/fs/cbcb-scratch/tluan/corr_exp/dada2_exp/dada2_luo_by_sample_v1_new.txt")

sink("/fs/cbcb-scratch/tluan/corr_exp/dada2_exp/dada2_luo_by_sample_new.txt")

colSums(seqtab)

sink()