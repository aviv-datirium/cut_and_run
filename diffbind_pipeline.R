#!/usr/bin/env Rscript
## diffbind_pipeline.R
## Usage: Rscript diffbind_pipeline.R sample_sheet.csv output_dir

suppressPackageStartupMessages({
  library(DiffBind)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2)
  stop("Usage: diffbind_pipeline.R sample_sheet.csv output_dir")

sheet  <- args[1]
outdir <- args[2]

samples <- read.csv(sheet, stringsAsFactors = FALSE)
if (length(unique(samples$Condition)) < 2) {
  message("Only one condition in sample sheet – nothing to compare.")
  quit(save = "no", status = 0)
}

db <- dba(sampleSheet = samples)
db <- dba.count(db)
db <- dba.contrast(db, categories = DBA_CONDITION)
db <- dba.analyze(db, method = DBA_DESEQ2)

res <- dba.report(db, method = DBA_DESEQ2, th = 0.05)
write.csv(as.data.frame(res), file = file.path(outdir, "diffbind_deseq2.csv"))

png(file.path(outdir, "diffbind_MA.png"), 800, 600)
plot(res, type = "MA")
dev.off()

png(file.path(outdir, "diffbind_heatmap.png"), 800, 600)
dba.plotHeatmap(db, correlations = FALSE, contrast = 1)
dev.off()
