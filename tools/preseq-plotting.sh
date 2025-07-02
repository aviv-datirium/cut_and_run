###############################################################################
#  PRESEQ PLOTTING                                                            #
###############################################################################

plot_preseq_curves () {
  local preseq_dir="$OUTPUT_DIR/preseq"
  local plot_pdf="$OUTPUT_DIR/preseq/preseq_complexity_curves.pdf"
  local r_script="$OUTPUT_DIR/preseq/plot_preseq.R"

  log Preseq Plotting "start"

  # Generate R script
  cat > "$r_script" <<'EOF'
library(ggplot2)
library(data.table)

plot_file <- function(file) {
  df <- fread(file, skip=1)
  df[, sample := gsub("_complexity\\.txt$", "", basename(file))]
  setnames(df, c("total_reads", "expected_unique", "sample"))
  return(df)
}

# Read all preseq outputs
files <- list.files("preseq", pattern = "_complexity\\.txt$", full.names = TRUE)
all_data <- rbindlist(lapply(files, plot_file))

# Plot
p <- ggplot(all_data, aes(x = total_reads, y = expected_unique, color = sample)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(title = "Preseq Library Complexity Projection",
       x = "Total Reads Sequenced",
       y = "Expected Unique Reads") +
  theme(legend.title = element_blank())

ggsave("preseq/preseq_complexity_curves.pdf", plot = p, width = 8, height = 6)
EOF

  # Run the R script
  Rscript "$r_script" >> "$LOG_DIR/preseq_plot.log" 2>&1 \
    && log Preseq Plotting "done: $plot_pdf" \
    || log Preseq Plotting FAIL
}

plot_preseq_curves
