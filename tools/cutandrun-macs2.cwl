cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest

  # Stage all inputs under /tmp, matching your JSON paths
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json
      - entry: $(inputs.fastq_dir)
        entryname: fastq
      - entry: $(inputs.star_indices)
        entryname: star_indices
      - entry: $(inputs.chrom_sizes)
        entryname: chrom/hg38.chrom.sizes
      - entry: $(inputs.annotation_genes)
        entryname: annotation/hg38.refGene.gtf

inputs:
  config_json:
    type: File
  fastq_dir:
    type: Directory
  star_indices:
    type: Directory
  chrom_sizes:
    type: File
  annotation_genes:
    type: File

# Bypass the image’s login‐shell ENTRYPOINT and run exactly what worked:
baseCommand:
  - bash
  - -lc
  - |
    source /opt/conda/etc/profile.d/conda.sh && \
    conda activate cutrun && \
    exec /usr/local/bin/cutrun.sh /tmp/config_for_docker.json

outputs:
  alignment_replicates:
    type: Directory
    outputBinding:
      glob: alignment_replicates
  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates
  cutrun_stdout:
    type: stdout
  cutrun_stderr:
    type: stderr
