cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"

  InitialWorkDirRequirement:
    listing:
      # 1) launcher script
      - class: Dirent
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) config JSON, under this exact name
      - class: File
        location: $(inputs.config_json.location)
        basename: config_for_docker.json

      # 3) data dirs & files exactly how cutrun.sh expects them:
      - class: Directory
        location: $(inputs.fastq_dir.location)
        basename: fastq

      - class: Directory
        location: $(inputs.reference_genome_dir.location)
        basename: star_indices/hg38

      - class: Directory
        location: $(inputs.ecoli_index_dir.location)
        basename: star_indices/ecoli_canonical

      - class: File
        location: $(inputs.chrom_sizes.location)
        basename: chrom/hg38.chrom.sizes

      - class: File
        location: $(inputs.annotation_genes.location)
        basename: annotation/hg38.refGene.gtf

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

  cutrun_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log

  cutrun_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
