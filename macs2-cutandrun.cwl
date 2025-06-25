baseCommand: [/usr/local/bin/cutrun.sh]
arguments:
  - --config $(inputs.config.path)
inputs:
  config:
    type: File
    inputBinding: {}
  # mount additional directories as CWL Files / Directories
