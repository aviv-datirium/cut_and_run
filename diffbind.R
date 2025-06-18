#!/usr/bin/env Rscript
# diffbind.R   <sampleSheet.csv>  <outdir>
# Produces:    diffbind_DESeq2_report.csv
#              DiffBind_PCA.pdf
#              DiffBind_corr_heatmap.pdf
#              DiffBind_binding_heatmap.pdf
#              DiffBind_MA.pdf

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2)
  stop("Usage: diffbind.R sampleSheet.csv outDir", call. = FALSE)

sheet  <- args[1]
outdir <- args[2]

suppressPackageStartupMessages({
  library(DiffBind)
})

## 1. Load sample sheet, count reads, build contrast
db <- dba(sampleSheet = sheet)
db <- dba.count(db, bUseSummarizeOverlaps = TRUE)
db <- dba.contrast(db, categories = DBA_CONDITION, minMembers = 2)

## 2. Differential analysis (DESeq2)
db <- dba.analyze(db, method = DBA_DESEQ2)

## 3. Save DESeq2 report
rep <- dba.report(db, method = DBA_DESEQ2, th = 0.05)   # FDR ≤ 0.05
write.csv(as.data.frame(rep),
          file = file.path(outdir, "diffbind_DESeq2_report.csv"),
          row.names = FALSE)

## 4. Plots -------------------------------------------------------------------
pdf(file.path(outdir, "DiffBind_PCA.pdf"))
dba.plotPCA(db, DBA_CONDITION, label = DBA_ID)
dev.off()

pdf(file.path(outdir, "DiffBind_corr_heatmap.pdf"))
dba.plotHeatmap(db, correlations = TRUE, scale = "none")
dev.off()

pdf(file.path(outdir, "DiffBind_binding_heatmap.pdf"), width = 8, height = 10)
dba.plotHeatmap(db, contrast = 1, correlations = FALSE, scale = "row")
dev.off()

## 5. MA plot (log2 fold-change vs mean abundance) ----------------------------
pdf(file.path(outdir, "DiffBind_MA.pdf"))
dba.plotMA(db, method = DBA_DESEQ2, bLoess = TRUE, th = 0.05)
dev.off()
