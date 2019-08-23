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
    input: rules.Index.output, htmlOutputPath + "/aberrant_splicing_readme.html"
    output: touch("tmp/aberrant_splicing.done")
    
    
### RULEGRAPH  
### rulegraph only works without print statements. Call <snakemake produce_rulegraph> for producing output

## For rule rulegraph.. copy configfile in tmp file
import oyaml
with open('tmp/config.yaml', 'w') as yaml_file:
    oyaml.dump(config, yaml_file, default_flow_style=False)

rulegraph_filename = htmlOutputPath + "/" + os.path.basename(os.getcwd()) + "_rulegraph"
rule produce_rulegraph:
    input:
        expand(rulegraph_filename + ".{fmt}", fmt=["svg", "png"])

rule create_graph:
    output:
        rulegraph_filename + ".dot"
    shell:
        "snakemake --configfile tmp/config.yaml --rulegraph > {output}"

rule render_dot:
    input:
        "{prefix}.dot"
    output:
        "{prefix}.{fmt,(png|svg)}"
    shell:
        "dot -T{wildcards.fmt} < {input} > {output}"
