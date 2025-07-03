cwlVersion: v1.2
class: CommandLineTool

doc: |
  Runs CUT&RUN pipeline inside Docker, staging the entire project directory under `/project`.

baseCommand:
  - bash
  - -c

requirements:
  - class: DockerRequirement
    dockerPull: biowardrobe2/cutrun-macs2-core:latest

  - class: InitialWorkDirRequirement
    listing:
      # Mount entire project tree
      - entry: $(inputs.project_dir)
        entryname: project
      # Stage the JSON inside project
      - entry: $(inputs.config_json)
        entryname: project/config_for_docker.json

inputs:
  project_dir:
    type: Directory
    doc: |
      Your project root containing directories:
      fastq/, star_indices/, chrom/, annotation/, etc.

  config_json:
    type: File
    doc: |
      JSON config with relative paths (e.g. "annotation/hg38.refGene.gtf").

arguments:
  - |
    # Change into the staged project and invoke the pipeline
    cd project && bash /usr/local/bin/cutrun.sh config_for_docker.json

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: project/output_replicates

  log_stdout:
    type: File
    outputBinding:
      glob: project/cutrun_stdout.log

  log_stderr:
    type: File
    outputBinding:
      glob: project/cutrun_stderr.log
