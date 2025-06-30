# cutrun.cwl — wraps the Dockerized bash CUT&RUN pipeline as a CWL CommandLineTool

cwlVersion: v1.2
class: CommandLineTool

baseCommand: [bash, /usr/local/bin/cutrun.sh]

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
      glob: $(inputs.config_json.basename.split("config_for_docker.json")[0]) + "output_replicates"
    doc: The main output directory produced by the pipeline

hints:
  DockerRequirement:
    dockerPull: biowardrobe2/cutrun-macs2-core:v1.0.0

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.config_json)
        entryname: config_for_docker.json

stdout: cutrun_stdout.log
stderr: cutrun_stderr.log

arguments:
  - valueFrom: |
      mkdir -p $(inputs.config_json.basename.split("config_for_docker.json")[0])output_replicates
    shellQuote: false
