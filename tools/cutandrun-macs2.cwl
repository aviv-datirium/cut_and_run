cwlVersion: v1.2
class: CommandLineTool

baseCommand: [bash, -c]

requirements:
  # Run inside your Cut&Run Docker image
  - class: DockerRequirement
    dockerPull: biowardrobe2/cutrun-macs2-core:latest

  # Stage exactly the files & dirs your JSON refers to, in-place
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json
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

inputs:
  config_json:
    type: File
    doc: |
      Your JSON with only the _relative_ paths (e.g. `"fastq/min_msto211h/…"`).
  fastq_dir:
    type: Directory
    doc: "The directory containing all your FASTQ files"
  reference_genome_dir:
    type: Directory
    doc: "Your STAR hg38 index"
  ecoli_index_dir:
    type: Directory
    doc: "Your STAR E. coli index"
  chrom_sizes:
    type: File
    doc: "chrom/hg38.chrom.sizes"
  annotation_genes:
    type: File
    doc: "annotation/hg38.refGene.gtf"

arguments:
  - |
    # invoke your pipeline inside the image
    bash /usr/local/bin/cutrun.sh config_for_docker.json

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
