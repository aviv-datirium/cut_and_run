
#~ Idea: Slide small windows (e.g. 150 bp) across the genome, count fragments
#~ per sample, normalise with edgeR, then merge significant adjacent windows.

#~ Why it often beats DiffBind

#~     - No reliance on peak callers → detects subtle shifts within broad regions.
#~     - edgeR’s tagwise dispersion is robust with two replicates per group.


library(csaw)
library(edgeR)

bam.files <- c(        # filtered BAMs
  "/mnt/data/home/aviv/alignment_replicates_diffbind/MYC-MST1_S28.dedup.filtered.bam",
  "/mnt/data/home/aviv/alignment_replicates_diffbind/MYC-MST2_S29.dedup.filtered.bam",
  "/mnt/data/home/aviv/alignment_replicates_diffbind/MYC-MSC1_S26.dedup.filtered.bam",
  "/mnt/data/home/aviv/alignment_replicates_diffbind/MYC-MSC2_S27.dedup.filtered.bam")

groups <- factor(c("trt","trt","ctl","ctl"))

# 1. Count reads in 150-bp windows
param <- readParam(minq=10, pe="both")      # paired-end, MAPQ≥10
windows <- windowCounts(bam.files, width=150, ext=200, param=param)

# 2. EdgeR differential test
y <- asDGEList(windows)
y <- calcNormFactors(y)
design <- model.matrix(~groups)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
res <- glmQLFTest(fit, coef=2)              # treatment vs control

# 3. Merge contiguous sig windows (FDR<0.05)
tab <- res$table
sig  <- tab$FDR < 0.05
merged <- mergeWindows(rowRanges(windows)[sig], tol=100)
combined <- combineTests(merged$id, res$table)

# 4. Export BED of significant regions
library(rtracklayer)
bed <- rowRanges(merged$region)
export(bed[combined$FDR<0.05], "csaw_diffPeaks.bed")
