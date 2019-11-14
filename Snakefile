### SNAKEFILE ABERRANT SPLICING

import os
import drop

METHOD = 'AS'
SCRIPT_ROOT = os.getcwd()

parser = drop.config(config, METHOD)
config = parser.parse()
include: config['wBuildPath'] + "/wBuild.snakefile"


rule all:
    input: rules.Index.output, config["htmlOutputPath"] + "/aberrant_splicing_readme.html"
    output: touch(drop.getMethodPath(METHOD, type_='final_file'))

### RULEGRAPH
config_file = drop.getConfFile()
rulegraph_filename = f'{config["htmlOutputPath"]}/{METHOD}_rulegraph'

rule produce_rulegraph:
    input:
        expand(rulegraph_filename + ".{fmt}", fmt=["svg", "png"])

rule create_graph:
    output:
        svg = f"{rulegraph_filename}.svg",
        png = f"{rulegraph_filename}.png"
    shell:
        """
        snakemake --configfile {config_file} --rulegraph | dot -Tsvg > {output.svg}
        snakemake --configfile {config_file} --rulegraph | dot -Tpng > {output.png}
        """

