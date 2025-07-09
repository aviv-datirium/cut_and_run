cwlVersion: v1.2
class: CommandLineTool

requirements:
  # 1) Allow $(...) in InitialWorkDirRequirement
  InlineJavascriptRequirement: {}

  # 2) Pull your freshly-built image
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"

  # 3) Stage only exactly what the pipeline expects under CWD
  InitialWorkDirRequirement:
    listing:
      # (1) tiny launcher script
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$$(pwd)"
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # (2) your JSON config
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json
        class: File

      # (3) the FASTQs tree
      - entry: $(inputs.fastq_dir.path)
        entryname: fastq
        class: Directory

      # (4) host‐genome STAR index
      - entry: $(inputs.reference_genome_dir.path)
        entryname: star_indices/hg38
        class: Directory

      # (5) spike‐in STAR index
      - entry: $(inputs.ecoli_index_dir.path)
        entryname: star_indices/ecoli_canonical
        class: Directory

      # (6) chrom.sizes
      - entry: $(inputs.chrom_sizes.path)
        entryname: chrom/hg38.chrom.sizes
        class: File

      # (7) GTF
      - entry: $(inputs.annotation_genes.path)
        entryname: annotation/hg38.refGene.gtf
        class: File

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
