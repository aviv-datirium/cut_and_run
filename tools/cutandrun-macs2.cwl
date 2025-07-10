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
          bash /usr/local/bin/cutrun.sh '' config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) config JSON
      - class: File
        path: $(inputs.config_json.path)
        basename: config_for_docker.json

      # 3) mount in your inputs exactly as the pipeline expects
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
