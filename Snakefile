### SNAKEFILE ABERRANT SPLICING
from pathlib import Path
import os
import drop

cfg = drop.config.DropConfig(config)
sa = cfg.sampleAnnotation
config = cfg.config # for legacy

METHOD = 'AS'
SCRIPT_ROOT = drop.getMethodPath(METHOD, type_='workdir', str_=False)
CONF_FILE = drop.getConfFile(METHOD)

include: drop.utils.getWBuildSnakefile()

rule all:
    input: rules.Index.output, config["htmlOutputPath"] + "/aberrant_splicing_readme.html"
    output: touch(drop.getMethodPath(METHOD, type_='final_file'))

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
        snakemake --configfile {CONF_FILE} --rulegraph | dot -Tsvg > {output.svg}
        snakemake --configfile {CONF_FILE} --rulegraph | dot -Tpng > {output.png}
        """

rule unlock:
    output: touch(drop.getMethodPath(METHOD, type_="unlock"))
    shell: "snakemake --unlock --configfile {CONF_FILE}"

