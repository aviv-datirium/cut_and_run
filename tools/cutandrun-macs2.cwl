cwlVersion: v1.2
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

  - class: DockerRequirement
    dockerPull: "cutrun-macs2-core:latest"

  - class: InitialWorkDirRequirement
    listing:
      # 1) launcher script: make the two dirs, then call cutrun.sh
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          mkdir -p output_replicates alignment_replicates
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) config JSON (full File object!)
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

      # 3) fastq directory
      - entry: $(inputs.fastq_dir)
        entryname: fastq

      # 4) genome indices
      - entry: $(inputs.reference_genome_dir)
        entryname: star_indices/hg38

      - entry: $(inputs.ecoli_index_dir)
        entryname: star_indices/ecoli_canonical

      # 5) chromosome sizes
      - entry: $(inputs.chrom_sizes)
        entryname: chrom/hg38.chrom.sizes

      # 6) gene annotation GTF
      - entry: $(inputs.annotation_genes)
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
