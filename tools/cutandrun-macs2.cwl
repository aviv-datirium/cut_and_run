cwlVersion: v1.2
class: CommandLineTool

baseCommand:
  - bash
  - run.sh
  - config_for_docker.json

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"
  ResourceRequirement:
    coresMin: 1
    ramMin: 256
  InitialWorkDirRequirement:
    listing:
      # 1) wrapper script
      - class: Dirent
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) our JSON config
      - class: Dirent
        entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

      # 3) fastqs
      - class: Dirent
        entry: $(inputs.fastq_dir.path)
        entryname: fastq

      # 4) hg38 STAR index
      - class: Dirent
        entry: $(inputs.reference_genome_dir.path)
        entryname: star_indices/hg38

      # 5) E.coli STAR index
      - class: Dirent
        entry: $(inputs.ecoli_index_dir.path)
        entryname: star_indices/ecoli_canonical

      # 6) chrom sizes
      - class: Dirent
        entry: $(inputs.chrom_sizes.path)
        entryname: chrom/hg38.chrom.sizes

      # 7) gene annotation
      - class: Dirent
        entry: $(inputs.annotation_genes.path)
        entryname: annotation/hg38.refGene.gtf

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

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

outputs:
  cutrun_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log

  cutrun_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log

  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates

  alignment_replicates:
    type: Directory
    outputBinding:
      glob: alignment_replicates
