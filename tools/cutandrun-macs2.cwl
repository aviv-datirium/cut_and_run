cwlVersion: v1.2
class: CommandLineTool

requirements:
  # allow JS in entry: and elsewhere
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

      # 2) JSON config
      - class: Dirent
        entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

      # 3) fastq dir
      - class: Dirent
        entry: $(inputs.fastq_dir.path)
        entryname: fastq

      # 4) host STAR index
      - class: Dirent
        entry: $(inputs.reference_genome_dir.path)
        entryname: star_indices/hg38

      # 5) spike STAR index
      - class: Dirent
        entry: $(inputs.ecoli_index_dir.path)
        entryname: star_indices/ecoli_canonical

      # 6) chrom sizes
      - class: Dirent
        entry: $(inputs.chrom_sizes.path)
        entryname: chrom/hg38.chrom.sizes

      # 7) annotation GTF
      - class: Dirent
        entry: $(inputs.annotation_genes.path)
        entryname: annotation/hg38.refGene.gtf

baseCommand: [ bash, run.sh ]

inputs:
  config_json:          File
  fastq_dir:            Directory
  reference_genome_dir: Directory
  ecoli_index_dir:      Directory
  chrom_sizes:          File
  annotation_genes:     File

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
