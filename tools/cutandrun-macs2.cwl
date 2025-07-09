cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest
  InitialWorkDirRequirement:
    listing:
      # (1) launcher script
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        basename: run.sh
        writable: true

      # (2) your config JSON
      - entry: $(inputs.config_json)
        basename: config_for_docker.json

      # (3) fastq dir
      - entry: $(inputs.fastq_dir)
        basename: fastq

      # (4) host STAR index
      - entry: $(inputs.reference_genome_dir)
        basename: star_indices/hg38

      # (5) spike-in STAR index
      - entry: $(inputs.ecoli_index_dir)
        basename: star_indices/ecoli_canonical

      # (6) chrom sizes
      - entry: $(inputs.chrom_sizes)
        basename: chrom/hg38.chrom.sizes

      # (7) annotation GTF
      - entry: $(inputs.annotation_genes)
        basename: annotation/hg38.refGene.gtf

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
