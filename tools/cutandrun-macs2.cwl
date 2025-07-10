cwlVersion: v1.2
class: CommandLineTool

# -------------------------------------------------------------------
# Tell CWL to use the container and to stage in ONLY your inputs.
# We will call the image's own entrypoint-script 'cutrun.sh', not a
# host file called run.sh.
# -------------------------------------------------------------------

requirements:
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest

  InitialWorkDirRequirement:
    listing:
      # 1) your config JSON, renamed to what the container expects
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

      # 2) the FASTQ directory
      - entry: $(inputs.fastq_dir)
        entryname: fastq

      # 3) your combined index directory (hg38/ + ecoli_canonical/)
      - entry: $(inputs.star_indices)
        entryname: star_indices

      # 4) chromosome sizes file
      - entry: $(inputs.chrom_sizes)
        entryname: chrom/hg38.chrom.sizes

      # 5) annotation GTF
      - entry: $(inputs.annotation_genes)
        entryname: annotation/hg38.refGene.gtf

# -------------------------------------------------------------------
# Instead of 'bash run.sh', invoke the container's own script:
# -------------------------------------------------------------------

baseCommand:
  - bash
  - -l           # login shell so conda auto-activate kicks in
  - -c
  - "/usr/local/bin/cutrun.sh config_for_docker.json"

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

outputs:
  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates

  alignment_replicates:
    type: Directory
    outputBinding:
      glob:
        - alignment_replicates
        - alignment_replicates/**

  cutrun_stdout:
    type: stdout

  cutrun_stderr:
    type: stderr
