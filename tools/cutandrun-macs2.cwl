cwlVersion: v1.2
class: CommandLineTool

# we will shell out to bash -c "...", so use bash -c as our base
baseCommand: [ bash, -c ]

requirements:
  # 1) Docker image
  - class: DockerRequirement
    dockerPull: biowardrobe2/cutrun-macs2-core:latest

  # 2) we copy the config.json into the container as config_for_docker.json
  - class: InitialWorkDirRequirement
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

  # 3) allow JS expressions in `$(...)`
  - class: InlineJavascriptRequirement

inputs:
  # 1) mount your whole project (so that relative paths in the JSON work)
  project_dir:
    type: Directory
    inputBinding:
      position: 1
    doc: Root of your repository (contains fastq/, star_indices/, etc.)

  # 2) the config JSON itself
  config_json:
    type: File
    inputBinding:
      position: 2
    doc: Your CWL‐fed config file

  # these are “dummy” inputs—just so CWL knows to mount them—but note we do
  # *not* pass them on the command line, we rely on the JSON paths instead
  fastq_dir:
    type: Directory
    doc: Mounted so cutrun.sh can read FASTQs
  reference_genome_dir:
    type: Directory
    doc: Mounted so cutrun.sh can read STAR index
  ecoli_index_dir:
    type: Directory
    doc: Mounted so cutrun.sh can read spike‐in index
  chrom_sizes:
    type: File
    doc: Mounted so cutrun.sh can read chrom sizes
  annotation_genes:
    type: File
    doc: Mounted so cutrun.sh can read GTF

arguments:
  # cd into your project root, then call cutrun.sh on the staged
  # config (always named config_for_docker.json)
  - |
    cd "$(inputs.project_dir.path)" && \
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

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log
