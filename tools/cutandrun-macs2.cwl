cwlVersion: v1.2
class: CommandLineTool

baseCommand: [ bash, run.sh ]

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: "biowardrobe2/cutrun-macs2-core:latest"
  - class: InitialWorkDirRequirement
    listing:
      # 1) the wrapper that invokes the container’s pipeline script
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          cd "$(pwd)"
          exec bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true
      # 2) your config JSON
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

inputs:
  config_json:
    type: File
    doc: "Your config_for_macs2_cwl.json"
  fastq_dir:
    type: Directory
    doc: "Directory of FASTQ files"
  reference_genome_dir:
    type: Directory
    doc: "STAR genome index (hg38)"
  ecoli_index_dir:
    type: Directory
    doc: "STAR E. coli index"
  chrom_sizes:
    type: File
    doc: "chrom/hg38.chrom.sizes"
  annotation_genes:
    type: File
    doc: "annotation/hg38.refGene.gtf"

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: output_replicates
  log_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log
  log_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log
