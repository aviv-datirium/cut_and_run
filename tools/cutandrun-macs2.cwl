cwlVersion: v1.2
class: CommandLineTool

baseCommand:
  - bash
  - /usr/local/bin/cutrun.sh

requirements:
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"

inputs:
  config_json:
    type: File
    inputBinding:
      position: 1
      prefix: ""

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
