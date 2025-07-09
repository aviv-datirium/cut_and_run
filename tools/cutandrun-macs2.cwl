cwlVersion: v1.2
class: CommandLineTool

baseCommand:
  - bash
  - /usr/local/bin/cutrun.sh

requirements:
  DockerRequirement:
    dockerPull: "cutrun-macs2-core:latest"

inputs:
  # only the config file goes on the command‐line
  config_json:
    type: File
    inputBinding:
      position: 0
      prefix: ""

  # all other inputs are just mounted, not passed
  fastq_dir:            Directory
  reference_genome_dir: Directory
  ecoli_index_dir:      Directory
  chrom_sizes:          File
  annotation_genes:     File

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
