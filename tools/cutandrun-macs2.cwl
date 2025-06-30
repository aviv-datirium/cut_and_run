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
    dockerPull: biowardrobe2/cutrun-macs2-core:v1.0.5  # <-- update tag if rebuilt

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      # Mount entire directory with all relevant input files
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

      # FASTQ files (relative to config file)
      - entry: $(inputs.config_json.path.replace(/[^\/]+$/, 'fastq/min_msto211h'))
        entryname: fastq/min_msto211h
        writable: false

      # STAR indices
      - entry: $(inputs.config_json.path.replace(/[^\/]+$/, 'star_indices/hg38'))
        entryname: star_indices/hg38
        writable: false
      - entry: $(inputs.config_json.path.replace(/[^\/]+$/, 'star_indices/ecoli_canonical'))
        entryname: star_indices/ecoli_canonical
        writable: false

      # Chrom sizes file
      - entry: $(inputs.config_json.path.replace(/[^\/]+$/, 'chrom/hg38.chrom.sizes'))
        entryname: chrom/hg38.chrom.sizes
        writable: false

      # Annotation GTF
      - entry: $(inputs.config_json.path.replace(/[^\/]+$/, 'annotation/hg38.refGene.gtf'))
        entryname: annotation/hg38.refGene.gtf
        writable: false
