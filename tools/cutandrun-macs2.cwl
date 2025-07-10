cwlVersion: v1.2
class: CommandLineTool

baseCommand:
  - bash
  - run.sh
  - config_for_docker.json

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest
  ResourceRequirement:
    coresMin: 1
    ramMin: 256
  InitialWorkDirRequirement:
    listing:
      # 1) our wrapper script
      - class: File
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) the JSON
      - class: File
        entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

      # 3) fastqs
      - class: Directory
        entry: $(inputs.fastq_dir.path)
        entryname: fastq

      # 4) hg38 STAR index
      - class: Directory
        entry: $(inputs.reference_genome_dir.path)
        entryname: star_indices/hg38

      # 5) E. coli index
      - class: Directory
        entry: $(inputs.ecoli_index_dir.path)
        entryname: star_indices/ecoli_canonical

      # 6) chrom sizes
      - class: File
        entry: $(inputs.chrom_sizes.path)
        entryname: chrom/hg38.chrom.sizes

      # 7) gene annotations
      - class: File
        entry: $(inputs.annotation_genes.path)
        entryname: annotation/hg38.refGene.gtf

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

inputs:
  config_json:      File
  fastq_dir:        Directory
  reference_genome_dir: Directory
  ecoli_index_dir:  Directory
  chrom_sizes:      File
  annotation_genes: File

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
