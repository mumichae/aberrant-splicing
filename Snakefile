import pandas as pd

## Add parser
from config_parser import ConfigHelper
parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile


htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"

include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"

rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch(htmlOutputPath + "/all.done")
