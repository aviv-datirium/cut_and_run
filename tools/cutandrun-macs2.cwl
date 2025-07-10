cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"
  InitialWorkDirRequirement:
    listing:
      # 1) launcher script
      - class: File
        entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        basename: run.sh

      # 2) copy in the config JSON
      - class: File
        location: $(inputs.config_json.path)
        basename: config_for_docker.json

      # 3) copy in all of your data directories/files
      - class: Directory
        location: $(inputs.fastq_dir.path)
        basename: fastq

      - class: Directory
        location: $(inputs.reference_genome_dir.path)
        basename: star_indices/hg38

      - class: Directory
        location: $(inputs.ecoli_index_dir.path)
        basename: star_indices/ecoli_canonical

      - class: File
        location: $(inputs.chrom_sizes.path)
        basename: chrom/hg38.chrom.sizes

      - class: File
        location: $(inputs.annotation_genes.path)
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
