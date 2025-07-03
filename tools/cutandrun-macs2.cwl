cwlVersion: v1.2
class: CommandLineTool

baseCommand: [ bash, -c ]

requirements:
  DockerRequirement:
    dockerPull: biowardrobe2/cutrun-macs2-core:latest

  InitialWorkDirRequirement:
    listing:
      # copy in your JSON (we still need this as a declared File)
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

      # copy in the entire fastq tree
      - entry: fastq/min_msto211h
        entryname: fastq

      # copy in your STAR indices
      - entry: star_indices/hg38
        entryname: star_indices/hg38
      - entry: star_indices/ecoli_canonical
        entryname: star_indices/ecoli_canonical

      # copy in chrom sizes & GTF
      - entry: chrom/hg38.chrom.sizes
        entryname: chrom/hg38.chrom.sizes
      - entry: annotation/hg38.refGene.gtf
        entryname: annotation/hg38.refGene.gtf

arguments:
  - |
    # now everything lives under the workdir exactly as your JSON expects:
    bash /usr/local/bin/cutrun.sh config_for_docker.json

inputs:
  config_json:
    type: File
    doc: your JSON (with _relative_ paths: e.g. "annotation/hg38.refGene.gtf")

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

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log
