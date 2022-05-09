.libPaths("/fs/cbcb-scratch/tluan/dada2")
library(dada2); packageVersion("dada2")

path <- "/fs/cbcb-scratch/tluan/exp/dada2"
filtFs <- file.path(path, "filtered", "a.fastq.gz")
sample.names <- "F3D1"
names(filtFs) <- sample.names
fnFs <- sort(list.files(path, pattern="seqs.fastq", full.names = TRUE))
out <- filterAndTrim(fwd = fnFs, filt = filtFs,  maxN=0, maxEE=Inf, truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE) 
errF <- learnErrors(filtFs, multithread=TRUE)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
write.csv2(dadaFs$denoised, file = "results_luo.csv")
