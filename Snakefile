### SNAKEFILE ABERRANT SPLICING

import os
import drop

parser = drop.config(config)
config = parser.parse()
include: config['wBuildPath'] + "/wBuild.snakefile"

METHOD = 'AS'
SCRIPT_ROOT = drop.getMethodPath(METHOD, link_type='workdir')
TMP_DIR = config['tmpdir']

rule all:
    input: rules.Index.output, config["htmlOutputPath"] + "/aberrant_splicing_readme.html"
    output: touch(drop.getMethodPath(METHOD, link_type='final_file', tmp_dir=TMP_DIR))

### RULEGRAPH
config_file = drop.getMethodPath(METHOD, link_type='config_file', tmp_dir=TMP_DIR)
rulegraph_filename = f'{config["htmlOutputPath"]}/{METHOD}_rulegraph'

rule produce_rulegraph:
    input:
        expand(rulegraph_filename + ".{fmt}", fmt=["svg", "png"])

rule create_graph:
    output:
        rulegraph_filename + ".dot"
    shell:
        "snakemake --configfile {config_file} --rulegraph > {output}"

rule render_dot:
    input:
        "{prefix}.dot"
    output:
        "{prefix}.{fmt,(png|svg)}"
    shell:
        "dot -T{wildcards.fmt} < {input} > {output}"
