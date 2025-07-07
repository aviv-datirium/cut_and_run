#!/bin/bash
set -euo pipefail

CONFIG=$1   # e.g. config_for_seacr_docker.json

# 0. Load config
ALIGNMENT_DIR=$(jq -r .alignment_dir      "$CONFIG")
OUTPUT_DIR  =$(jq -r .output_dir          "$CONFIG")
LOG_DIR     =$(jq -r .log_dir             "$CONFIG")
BW_DIR      =$(jq -r .output_dir "$CONFIG")/bigwig_bedgraphs
PEAK_DIR    =$(jq -r .output_dir "$CONFIG")/seacr_peaks
SEACR_THRESH=$(jq -r .seacr.threshold     "$CONFIG")
SEACR_NORM  =$(jq -r .seacr.norm          "$CONFIG")
SEACR_STRICT=$(jq -r .seacr.stringency    "$CONFIG")

# 1. Build sample lists
TREAT_R1=( $(jq -r '.samples.treatment[]?.r1'   "$CONFIG") )
TREAT_R2=( $(jq -r '.samples.treatment[]?.r2'   "$CONFIG") )
CTRL_R1=(  $(jq -r '.samples.control[]?.r1 // empty' "$CONFIG") )
CTRL_R2=(  $(jq -r '.samples.control[]?.r2 // empty' "$CONFIG") )

get_basename(){ 
  b=$(basename "$1"); 
  echo "${b%.fastq.gz}" | sed -E 's/_R[12]_?0*//g'; 
}

for i in "${!TREAT_R1[@]}"; do TREAT_NAMES[i]=$(get_basename "${TREAT_R1[i]}"); done
for i in "${!CTRL_R1[@]}";  do CTRL_NAMES[i]=$(get_basename "${CTRL_R1[i]}");   done

# 2. Helpers
make_bg () { bedtools genomecov -bg -ibam "$2" -pc \
                | sort -k1,1 -k2,2n > "$1"; }

seacr_call() (
  # ensure writes go somewhere writable
  TMPDIR=$(mktemp -d)
  cd "$TMPDIR"
  seacr "$@"
)

# 3. Prep directories
mkdir -p "$BW_DIR" "$PEAK_DIR"/{replicate,pooled,merged} "$LOG_DIR"

# 4. Create bedgraphs for all samples & merged control
for n in "${TREAT_NAMES[@]}" "${CTRL_NAMES[@]}"; do
  make_bg "$BW_DIR/${n}.bedgraph" "$ALIGNMENT_DIR/${n}.dedup.filtered.bam"
done

if (( ${#CTRL_NAMES[@]} )); then
  make_bg "$BW_DIR/control_merged.bedgraph" "$ALIGNMENT_DIR/control_merged.bam"
  POOLED_C_BG="$BW_DIR/control_merged.bedgraph"
fi
make_bg "$BW_DIR/treatment_merged.bedgraph" "$ALIGNMENT_DIR/treatment_merged.bam"

# 5. SEACR calls

# A) Per-replicate
for n in "${TREAT_NAMES[@]}"; do
  BG="$BW_DIR/${n}.bedgraph"
  OUT="$PEAK_DIR/replicate/${n}_seacr.bed"
  if [[ -n "${POOLED_C_BG:-}" ]]; then
    seacr_call "$BG" "$POOLED_C_BG" "$SEACR_NORM" "$SEACR_STRICT" "${OUT%.bed}" \
      >>"$LOG_DIR/seacr_${n}.log" 2>&1
    mv "${OUT%.bed}.${SEACR_STRICT}.bed" "$OUT"
  else
    seacr_call "$BG" "$SEACR_THRESH" "$SEACR_NORM" "$SEACR_STRICT" \
      >>"$LOG_DIR/seacr_${n}.log" 2>&1
    mv "${BG%.bedgraph}.${SEACR_STRICT}.bed" "$OUT"
  fi
done

# B) Merged treatment vs control
if [[ -n "${POOLED_C_BG:-}" ]]; then
  seacr_call "$BW_DIR/treatment_merged.bedgraph" "$POOLED_C_BG" \
    "$SEACR_NORM" "$SEACR_STRICT" \
    "${PEAK_DIR}/merged/treatmentMerged_vs_controlMerged_seacr" \
    >>"$LOG_DIR/seacr_merged.log" 2>&1
  mv "${PEAK_DIR}/merged/treatmentMerged_vs_controlMerged_seacr.${SEACR_STRICT}.bed" \
     "${PEAK_DIR}/merged/treatmentMerged_vs_controlMerged_seacr.bed"
fi

# C) Replicate vs pooled control
if [[ -n "${POOLED_C_BG:-}" ]]; then
  for n in "${TREAT_NAMES[@]}"; do
    seacr_call "$BW_DIR/${n}.bedgraph" "$POOLED_C_BG" \
      "$SEACR_NORM" "$SEACR_STRICT" \
      "${PEAK_DIR}/pooled/${n}_vs_ctrlPooled_seacr" \
      >>"$LOG_DIR/seacr_${n}.log" 2>&1
    mv "${PEAK_DIR}/pooled/${n}_vs_ctrlPooled_seacr.${SEACR_STRICT}.bed" \
       "${PEAK_DIR}/pooled/${n}_vs_ctrlPooled_seacr.bed"
  done
fi
