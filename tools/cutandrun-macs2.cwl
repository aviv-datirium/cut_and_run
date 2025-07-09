cwlVersion: v1.2
class: CommandLineTool

requirements:
  # allow JS expressions everywhere
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}

  # pull your image
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest

  # stage exactly what the container needs in its CWD
  InitialWorkDirRequirement:
    listing:

      # 1) launcher script
      - class: File
        basename: run.sh
        contents: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        writable: true

      # 2) the JSON config
      - class: ExpressionDirent
        name: config_for_docker.json
        entry: $(inputs.config_json)

      # 3) fastq data dir
      - class: ExpressionDirent
        name: fastq
        entry: $(inputs.fastq_dir)

      # 4) host‐genome STAR index
      - class: ExpressionDirent
        name: star_indices/hg38
        entry: $(inputs.reference_genome_dir)

      # 5) spike‐in STAR index
      - class: ExpressionDirent
        name: star_indices/ecoli_canonical
        entry: $(inputs.ecoli_index_dir)

      # 6) chrom sizes file
      - class: ExpressionDirent
        name: chrom/hg38.chrom.sizes
        entry: $(inputs.chrom_sizes)

      # 7) annotation GTF
      - class: ExpressionDirent
        name: annotation/hg38.refGene.gtf
        entry: $(inputs.annotation_genes)

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

  cutrun_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log

  cutrun_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
