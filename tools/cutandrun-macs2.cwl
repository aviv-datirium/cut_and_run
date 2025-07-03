# cutrun.cwl — wraps the Dockerized bash CUT&RUN pipeline as a CWL CommandLineTool

cwlVersion: v1.2
class: CommandLineTool

baseCommand: ["bash", "/usr/local/bin/cutrun.sh"]

inputs:
  config_json:
    type: File
    inputBinding:
      position: 1
    doc: Path to the JSON config file

  fastq_dir:
    type: Directory
    doc: Directory containing FASTQ files
    
  reference_genome_dir:
    type: Directory
    doc: STAR genome index directory
    
  ecoli_index_dir:
    type: Directory
    doc: E. coli STAR index directory
    
  chrom_sizes:
    type: File
    doc: Chromosome sizes file
    
  annotation_genes:
    type: File
    doc: Gene annotation GTF file
    
outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: output_replicates
    doc: The main output directory produced by the pipeline

  log_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log
    doc: Standard output log from the pipeline

  log_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
    doc: Standard error log from the pipeline

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

hints:
  DockerRequirement:
    dockerPull: biowardrobe2/cutrun-macs2-core:v1.1.0  # <-- update tag if rebuilt

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

      - entry: $(inputs.fastq_dir)
        entryname: inputs/fastq/min_msto211h

      - entry: $(inputs.reference_genome_dir)
        entryname: inputs/star_indices/hg38

      - entry: $(inputs.ecoli_index_dir)
        entryname: inputs/star_indices/ecoli_canonical

      - entry: $(inputs.chrom_sizes)
        entryname: inputs/chrom/hg38.chrom.sizes

      - entry: $(inputs.annotation_genes)
        entryname: inputs/annotation/hg38.refGene.gtf
