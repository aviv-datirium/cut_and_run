cwlVersion: v1.2
class: CommandLineTool

baseCommand: ["/bin/bash","-c"]
arguments:
  - |
    cd $(inputs.project_dir.path) && bash /usr/local/bin/cutrun.sh $(inputs.config_json.path)

requirements:
  DockerRequirement:
    dockerPull: biowardrobe2/cutrun-macs2-core:latest
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

inputs:
  project_dir:
    type: Directory
    doc: "Path to project root containing fastq/, star_indices/, dockerfiles/, etc."
  config_json:
    type: File
    inputBinding:
      position: 1
    doc: "JSON config file for pipeline (dockerfiles/config_for_macs2_cwl.json)"
  fastq_dir:
    type: Directory
    doc: "Staged FASTQ directory"
  reference_genome_dir:
    type: Directory
    doc: "STAR host genome index directory"
  ecoli_index_dir:
    type: Directory
    doc: "STAR E. coli index directory"
  chrom_sizes:
    type: File
    doc: "Chromosome sizes file"
  annotation_genes:
    type: File
    doc: "Gene annotation GTF file"

outputs:
  output_dir:
    type: Directory
    outputBinding:
      glob: output_replicates
    doc: "Main output directory (output_replicates/)"
  log_stdout:
    type: File
    outputBinding:
      glob: cutrun_stdout.log
    doc: "Standard output log"
  log_stderr:
    type: File
    outputBinding:
      glob: cutrun_stderr.log
    doc: "Standard error log"

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log
