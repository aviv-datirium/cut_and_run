cwlVersion: v1.2
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

  - class: DockerRequirement
    dockerPull: "cutrun-macs2-core:latest"

  - class: InitialWorkDirRequirement
    listing:
      # 1) our launcher script — makes the two output dirs, then runs cutrun.sh
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          mkdir -p output_replicates alignment_replicates
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) copy in your JSON config
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

      # 3) the FASTQ directory
      - entry: $(inputs.fastq_dir.path)
        entryname: fastq

      # 4) genome indices
      - entry: $(inputs.reference_genome_dir.path)
        entryname: star_indices/hg38

      - entry: $(inputs.ecoli_index_dir.path)
        entryname: star_indices/ecoli_canonical

      # 5) chromosome sizes
      - entry: $(inputs.chrom_sizes.path)
        entryname: chrom/hg38.chrom.sizes

      # 6) gene annotation
      - entry: $(inputs.annotation_genes.path)
        entryname: annotation/hg38.refGene.gtf

baseCommand:
  - bash
  - run.sh

inputs:
  config_json:
    type: File

  fastq_dir:
    type: Directory

  reference_genome_dir:
    type: Directory

  ecoli_index_dir:
    type: Directory

  chrom_sizes:
    type: File

  annotation_genes:
    type: File

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

outputs:
  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates

  alignment_replicates:
    type: Directory
    outputBinding:
      glob: alignment_replicates

  cutrun_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log

  cutrun_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
