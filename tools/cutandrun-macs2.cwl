cwlVersion: v1.2
class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"
  InitialWorkDirRequirement:
    listing:
      # 1) launcher
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          # pass an empty first arg, then your JSON
          bash /usr/local/bin/cutrun.sh '' config_for_docker.json
        entryname: run.sh
        writable: true

      # 2) your config JSON
      - entry: $(inputs.config_json.path)
        entryname: config_for_docker.json

      # 3) data dirs & files as your config expects them
      - entry: $(inputs.fastq_dir.path)
        entryname: fastq
      - entry: $(inputs.reference_genome_dir.path)
        entryname: star_indices/hg38
      - entry: $(inputs.ecoli_index_dir.path)
        entryname: star_indices/ecoli_canonical
      - entry: $(inputs.chrom_sizes.path)
        entryname: chrom/hg38.chrom.sizes
      - entry: $(inputs.annotation_genes.path)
        entryname: annotation/hg38.refGene.gtf

baseCommand: [ bash, run.sh ]

inputs:
  config_json:       File
  fastq_dir:         Directory
  reference_genome_dir: Directory
  ecoli_index_dir:   Directory
  chrom_sizes:       File
  annotation_genes:  File

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
