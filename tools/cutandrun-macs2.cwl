cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest
  InitialWorkDirRequirement:
    listing:

      # 1) launcher script
      - class: Dirent
        entry:
          class: File
          basename: run.sh
          contents: |
            #!/usr/bin/env bash
            set -euo pipefail
            cd "$(pwd)"
            bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) config JSON
      - class: Dirent
        entry:
          class: File
          location: $(inputs.config_json.path)
        entryname: config_for_docker.json

      # 3) FASTQ directory
      - class: Dirent
        entry:
          class: Directory
          location: $(inputs.fastq_dir.path)
        entryname: fastq

      # 4) host-genome STAR index
      - class: Dirent
        entry:
          class: Directory
          location: $(inputs.reference_genome_dir.path)
        dirname: star_indices
        entryname: hg38

      # 5) E. coli STAR index
      - class: Dirent
        entry:
          class: Directory
          location: $(inputs.ecoli_index_dir.path)
        dirname: star_indices
        entryname: ecoli_canonical

      # 6) chromosome sizes
      - class: Dirent
        entry:
          class: File
          location: $(inputs.chrom_sizes.path)
        dirname: chrom
        entryname: hg38.chrom.sizes

      # 7) gene annotation
      - class: Dirent
        entry:
          class: File
          location: $(inputs.annotation_genes.path)
        dirname: annotation
        entryname: hg38.refGene.gtf

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
