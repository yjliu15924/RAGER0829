library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg38)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
args <- commandArgs(trailingOnly = TRUE)
#bamfile <- "/home/yjliu/RAGER/case2_test/datasets/ATACseq/bowtie2file/SRR28263042/SRR28263042.bam"
bamfile <- args[1]
bamfile.labels <- gsub(".bam", "", basename(bamfile))
#pdf("/home/yjliu/RAGER/case2_test/datasets/ATACseq/bowtie2file/SRR28263042/1.pdf")
pdf(args[2])
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()