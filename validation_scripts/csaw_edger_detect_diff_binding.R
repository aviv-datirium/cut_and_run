
#~ Idea: Slide small windows (e.g. 150 bp) across the genome, count fragments
#~ per sample, normalise with edgeR, then merge significant adjacent windows.

#~ Why it often beats DiffBind

#~     - No reliance on peak callers → detects subtle shifts within broad regions.
#~     - edgeR’s tagwise dispersion is robust with two replicates per group.

library(csaw)
library(edgeR)
library(rtracklayer)

# --- inputs ---
bam.files <- c("...S28.bam","...S29.bam","...S26.bam","...S27.bam")
groups    <- factor(c("trt","trt","ctl","ctl"))

# 1 count sliding windows
param    <- readParam(minq = 10, pe = "both")
windows  <- windowCounts(bam.files, width = 150, ext = 200, param = param)

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

keep     <- combined$FDR < 0.05 & abs(combined$logFC) > 0.5
sig.reg  <- merged$region[keep]
final    <- combined[keep ,]

# 4 export
export(sig.reg, "csaw_diffPeaks.bed")
write.csv(final,  "csaw_diffPeaks.csv")
