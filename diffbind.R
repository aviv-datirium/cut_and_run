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
  library(BiocParallel)
})

## --------------------------------------------------------------------------
## command-line args
## --------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
csv  <- args[1]        # sample sheet
out  <- args[2]        # output dir
cores <- as.integer(args[3])   # number of cores to use

## --------------------------------------------------------------------------
## parallel back-end
## --------------------------------------------------------------------------
bp <- MulticoreParam(workers = cores)   # on Windows use SnowParam
register(bp)                            # makes bpparam() use these cores

## 1. Load sample sheet, count reads, build contrast
db <- dba(sampleSheet = sheet)
db <- dba.count(db,
                summits = 250,          # keep if you use summit mode
                bUseSummarizeOverlaps = TRUE,
                BPPARAM = bp)           # <- parallel counting
db <- dba.contrast(db, categories = DBA_CONDITION)

## 2. Differential analysis (DESeq2)
db <- dba.analyze(db,
                  method = DBA_DESEQ2,
                  bParallel = TRUE,     # <- parallel DESeq2
                  BPPARAM = bp)
                  
## 3. Save DESeq2 report
report <- dba.report(db, method = DBA_DESEQ2, th = 0.1, fold = 0.5)   # normally: FDR 0.05, |log₂FC| ≥ 1 FDR ≤ 0.05
if(length(report) == 0L){
  message("No sites above threshold – skipping plots")
  quit(save = "no", status = 0)
}

write.csv(as.data.frame(rep),
          file = file.path(outdir, "diffbind_DESeq2_report.csv"),
          row.names = FALSE)

## 4. Plots -------------------------------------------------------------------
dba.plotPCA(db, attributes = DBA_CONDITION,
            file = file.path(out,"DiffBind_PCA.pdf"))
dba.plotHeatmap(db, file = file.path(out,"DiffBind_corr_heatmap.pdf"))
dba.plotMA(db, contrast = 1, method = DBA_DESEQ2,
           file = file.path(out,"DiffBind_MA.pdf"))
dev.off()

## export
write.csv(as.data.frame(report),
          file = file.path(out,"DiffBind_report.csv"))
