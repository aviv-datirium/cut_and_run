cwlVersion: v1.2
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest

  InitialWorkDirRequirement:
    listing:
      # 1) Write out our launcher script
      - class: Dirent
        entryname: run.sh
        writable: true
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          bash /usr/local/bin/cutrun.sh config_for_docker.json

      # 2) Drop in the JSON config
      - class: Dirent
        entryname: config_for_docker.json
        writable: false
        entry: $(inputs.config_json.path)

      # 3) And stage all the data directories/files exactly where cutrun expects them:
      - class: Dirent
        entryname: fastq
        writable: false
        entry: $(inputs.fastq_dir.path)

      - class: Dirent
        entryname: star_indices/hg38
        writable: false
        entry: $(inputs.reference_genome_dir.path)

      - class: Dirent
        entryname: star_indices/ecoli_canonical
        writable: false
        entry: $(inputs.ecoli_index_dir.path)

      - class: Dirent
        entryname: chrom/hg38.chrom.sizes
        writable: false
        entry: $(inputs.chrom_sizes.path)

      - class: Dirent
        entryname: annotation/hg38.refGene.gtf
        writable: false
        entry: $(inputs.annotation_genes.path)

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
