cwlVersion: v1.2
class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: biowardrobe2/cutrun-macs2-core:latest

  - class: InitialWorkDirRequirement
    listing:
      # bring in your JSON under its own name
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json
      # stage each Directory or File exactly where your JSON expects it
      - entry: $(inputs.fastq_dir)
        entryname: fastq
      - entry: $(inputs.reference_genome_dir)
        entryname: star_indices/hg38
      - entry: $(inputs.ecoli_index_dir)
        entryname: star_indices/ecoli_canonical
      - entry: $(inputs.chrom_sizes)
        entryname: chrom/hg38.chrom.sizes
      - entry: $(inputs.annotation_genes)
        entryname: annotation/hg38.refGene.gtf

baseCommand: ["bash", "-c"]
arguments:
  - |
    # run pipeline in the CWL scratch dir
    bash /usr/local/bin/cutrun.sh config_for_docker.json \
      > cutrun_stdout.log 2> cutrun_stderr.log

inputs:
  config_json:
    type: File
    doc: your JSON config for the pipeline

  fastq_dir:
    type: Directory
    doc: fastq/min_msto211h

  reference_genome_dir:
    type: Directory
    doc: star_indices/hg38

  ecoli_index_dir:
    type: Directory
    doc: star_indices/ecoli_canonical

  chrom_sizes:
    type: File
    doc: chrom/hg38.chrom.sizes

  annotation_genes:
    type: File
    doc: annotation/hg38.refGene.gtf

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
