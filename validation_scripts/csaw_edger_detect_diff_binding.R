
#~ Idea: Slide small windows (e.g. 150 bp) across the genome, count fragments
#~ per sample, normalise with edgeR, then merge significant adjacent windows.

#~ Why it often beats DiffBind

#~     - No reliance on peak callers → detects subtle shifts within broad regions.
#~     - edgeR’s tagwise dispersion is robust with two replicates per group.


library(csaw)
library(edgeR)
library(rtracklayer)
library(BiocParallel)

NTHREADS <- 8                       # choose cores (<= nproc)

register(MulticoreParam(NTHREADS))

# Inputs
bam.files <- c(        # filtered BAMs
  "/mnt/data/home/aviv/alignment_replicates_diffbind/MYC-MST1_S28.dedup.filtered.bam",
  "/mnt/data/home/aviv/alignment_replicates_diffbind/MYC-MST2_S29.dedup.filtered.bam",
  "/mnt/data/home/aviv/alignment_replicates_diffbind/MYC-MSC1_S26.dedup.filtered.bam",
  "/mnt/data/home/aviv/alignment_replicates_diffbind/MYC-MSC2_S27.dedup.filtered.bam")
groups <- factor(c("trt","trt","ctl","ctl"))

# 1 count sliding windows
param    <- readParam(minq = 10, pe = "both")
windows <- windowCounts(bam.files,
                        width = 150, ext = 200,
                        param = readParam(minq = 10, pe = "both"),
                        BPPARAM = bpparam())      # ← uses NTHREADS cores

# 2 edgeR QL test
y       <- asDGEList(windows)
y       <- calcNormFactors(y)
design  <- model.matrix(~ groups)
y       <- estimateDisp(y, design)
fit     <- glmQLFit(y, design)
res     <- glmQLFTest(fit, coef = 2)      # trt vs ctl

# 3 merge & combine
merged   <- mergeWindows(rowRanges(windows), tol = 100)
combined <- combineTests(merged$id, res$table)

# check / coerce representative logFC
if ("rep.logFC" %in% colnames(combined)) {
    combined$rep.logFC <- as.numeric(as.character(combined$rep.logFC))
} else {
    stop("combineTests did not return 'rep.logFC'; check column names.")
}

# 4 filter by FDR and |log2FC|
keep <- combined$FDR < 0.05 & abs(combined$rep.logFC) > 0.5
sig.reg <- merged$region[keep]
stats    <- combined[keep ,]

# 5 export
library(rtracklayer)
export(sig.reg, "csaw_diffPeaks.bed")
write.csv(as.data.frame(stats), "csaw_diffPeaks.csv")
