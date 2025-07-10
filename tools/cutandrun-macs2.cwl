cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest

  InitialWorkDirRequirement:
    listing:
      # 1) Your launcher script (static)
      - class: Dirent
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) Config JSON (dynamic, expression)
      - class: ExpressionDirent
        entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

      # 3) FASTQ directory
      - class: ExpressionDirent
        entry: $(inputs.fastq_dir.path)
        entryname: fastq

      # 4) STAR hg38 index
      - class: ExpressionDirent
        entry: $(inputs.reference_genome_dir.path)
        entryname: star_indices/hg38

      # 5) STAR E. coli index
      - class: ExpressionDirent
        entry: $(inputs.ecoli_index_dir.path)
        entryname: star_indices/ecoli_canonical

      # 6) Chrom sizes file
      - class: ExpressionDirent
        entry: $(inputs.chrom_sizes.path)
        entryname: chrom/hg38.chrom.sizes

      # 7) Gene annotation GTF
      - class: ExpressionDirent
        entry: $(inputs.annotation_genes.path)
        entryname: annotation/hg38.refGene.gtf

baseCommand: [ bash, run.sh ]

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
