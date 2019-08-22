### SNAKEFILE ABERRANT SPLICING

import pandas as pd
import os
from config_parser import ConfigHelper

## ADD tmp/ DIR
if not os.path.exists('tmp'):
    os.makedirs('tmp')

parser = ConfigHelper(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile

htmlOutputPath = config["htmlOutputPath"]
include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables

rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch("tmp/aberrant_splicing.done")
