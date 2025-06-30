# cutrun.cwl — wraps the Dockerized bash CUT&RUN pipeline as a CWL CommandLineTool

cwlVersion: v1.2
class: CommandLineTool

baseCommand: ["bash", "/usr/local/bin/cutrun.sh"]

inputs:
  config_json:
    type: File
    inputBinding:
      position: 1
    doc: Path to the JSON config file with all parameters

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
    dockerPull: biowardrobe2/cutrun-macs2-core:v1.0.0

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

      # FASTQ files
      - entry: /mnt/data/home/aviv/fastq/min_msto211h
        entryname: fastq/min_msto211h
        writable: false

      # STAR indices
      - entry: /mnt/data/home/aviv/star_indices/hg38
        entryname: star_indices/hg38
        writable: false

      - entry: /mnt/data/home/aviv/star_indices/ecoli_canonical
        entryname: star_indices/ecoli_canonical
        writable: false

      # Chrom sizes file
      - entry: /mnt/data/home/aviv/chrom/hg38.chrom.sizes
        entryname: chrom/hg38.chrom.sizes
        writable: false

      # Annotation GTF
      - entry: /mnt/data/home/aviv/annotation/hg38.refGene.gtf
        entryname: annotation/hg38.refGene.gtf
        writable: false
