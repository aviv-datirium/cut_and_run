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
      - entry: |
          #!/usr/bin/env bash
          set -euo pipefail
          bash /usr/local/bin/cutrun.sh config_for_docker.json
        entryname: run.sh
        writable: true

      - entry: $(inputs.config_json.location)
        entryname: config_for_docker.json

      - entry:
          class: Directory
          location: $(inputs.fastq_dir.location)
        entryname: fastq

      - entry:
          class: Directory
          location: $(inputs.reference_genome_dir.location)
        entryname: star_indices/hg38

      - entry:
          class: Directory
          location: $(inputs.ecoli_index_dir.location)
        entryname: star_indices/ecoli_canonical

      - entry:
          class: File
          location: $(inputs.chrom_sizes.location)
        entryname: chrom/hg38.chrom.sizes

      - entry:
          class: File
          location: $(inputs.annotation_genes.location)
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
