cwlVersion: v1.2
class: CommandLineTool

requirements:
  # allow us to use $(inputs.foo) in InitialWorkDirRequirement
  StepInputExpressionRequirement: {}
  # pull your freshly-built image
  DockerRequirement:
    dockerPull: cutrun-macs2-core:latest

  # stage config + data into exactly the layout cutrun.sh expects
  InitialWorkDirRequirement:
    listing:
      # config JSON must be named exactly this
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

      # raw FASTQs
      - entry: $(inputs.fastq_dir)
        entryname: fastq

      # host reference indices
      - entry: $(inputs.reference_genome_dir)
        entryname: star_indices/hg38

      # spike-in indices
      - entry: $(inputs.ecoli_index_dir)
        entryname: star_indices/ecoli_canonical

      # chromosome sizes & annotation
      - entry: $(inputs.chrom_sizes)
        entryname: chrom/hg38.chrom.sizes
      - entry: $(inputs.annotation_genes)
        entryname: annotation/hg38.refGene.gtf

baseCommand:
  - /usr/local/bin/cutrun.sh
arguments:
  - config_for_docker.json

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
