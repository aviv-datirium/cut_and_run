cwlVersion: v1.2
class: CommandLineTool

requirements:
  # for our $( ) expressions
  InlineJavascriptRequirement: {}

  # pull your image
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"

  # stage exactly what we need into the container CWD
  InitialWorkDirRequirement:
    listing:

      # (1) the launcher script (runs cutrun.sh)
      - class: Dirent
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # (2) your JSON config
      - class: File
        path: $(inputs.config_json.path)
        basename: config_for_docker.json

      # (3) data dirs & files laid out exactly under CWD
      - class: Directory
        path: $(inputs.fastq_dir.path)
        basename: fastq

      - class: Directory
        path: $(inputs.reference_genome_dir.path)
        basename: star_indices/hg38

      - class: Directory
        path: $(inputs.ecoli_index_dir.path)
        basename: star_indices/ecoli_canonical

      - class: File
        path: $(inputs.chrom_sizes.path)
        basename: chrom/hg38.chrom.sizes

      - class: File
        path: $(inputs.annotation_genes.path)
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
