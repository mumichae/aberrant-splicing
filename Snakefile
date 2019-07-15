import pandas as pd

## Add parser
from config_parser import ConfigHelper
parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile


htmlOutputPath = config["htmlOutputPath"]  
include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables

rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch("Output/all.done")
