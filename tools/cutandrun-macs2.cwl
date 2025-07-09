cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"

  InitialWorkDirRequirement:
    listing:
      - class: Dirent
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      - class: Dirent
        entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

      - class: Dirent
        entry: $(inputs.fastq_dir.path)
        entryname: fastq

      - class: Dirent
        entry: $(inputs.reference_genome_dir.path)
        entryname: star_indices/hg38

      - class: Dirent
        entry: $(inputs.ecoli_index_dir.path)
        entryname: star_indices/ecoli_canonical

      - class: Dirent
        entry: $(inputs.chrom_sizes.path)
        entryname: chrom/hg38.chrom.sizes

      - class: Dirent
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

  cutrun_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log

  cutrun_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
