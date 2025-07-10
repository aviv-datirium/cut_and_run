# tools/cutandrun-macs2.cwl
cwlVersion: v1.2
class: CommandLineTool

baseCommand:
  - bash
  - run.sh
  - config_for_docker.json

requirements:
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest

  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

      - entry: $(inputs.fastq_dir)
        entryname: fastq

      - entry: $(inputs.star_indices)
        entryname: star_indices

      - entry: $(inputs.chrom_sizes)
        entryname: chrom/hg38.chrom.sizes

      - entry: $(inputs.annotation_genes)
        entryname: annotation/hg38.refGene.gtf

inputs:
  config_json:
    type: File

  fastq_dir:
    type: Directory

  star_indices:
    type: Directory

  chrom_sizes:
    type: File

  annotation_genes:
    type: File

outputs:
  output_replicates:
    type: Directory
    outputBinding:
      glob: output_replicates

  alignment_replicates:
    type: Directory
    outputBinding:
      glob: alignment_replicates

  cutrun_stdout:
    type: stdout

  cutrun_stderr:
    type: stderr
