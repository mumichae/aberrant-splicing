### SNAKEFILE ABERRANT SPLICING

import os
import drop

tmpdir = os.path.join(config["ROOT"], 'tmp')
config["tmpdir"] = tmpdir
if not os.path.exists(tmpdir+'/AberrantSplicing'):
    os.makedirs(tmpdir+'/AberrantSplicing')
    
    
parser = drop.config(config)
config = parser.config # needed if you dont provide the wbuild.yaml as configfile

include: config['wBuildPath'] + "/wBuild.snakefile"

rule all:
    input: rules.Index.output, config["htmlOutputPath"] + "/aberrant_splicing_readme.html"
    output: touch(tmpdir + "/AS.done")

### RULEGRAPH  
### rulegraph only works without print statements

## For rule rulegraph.. copy configfile in tmp file
import oyaml
with open(tmpdir + '/config.yaml', 'w') as yaml_file:
    oyaml.dump(config, yaml_file, default_flow_style=False)

rulegraph_filename = htmlOutputPath + "/AS_rulegraph" # htmlOutputPath + "/" + os.path.basename(os.getcwd()) + "_rulegraph"
rule produce_rulegraph:
    input:
        expand(rulegraph_filename + ".{fmt}", fmt=["svg", "png"])

rule create_graph:
    output:
        rulegraph_filename + ".dot"
    shell:
        "snakemake --configfile " + tmpdir + "/config.yaml --rulegraph > {output}"

rule render_dot:
    input:
        "{prefix}.dot"
    output:
        "{prefix}.{fmt,(png|svg)}"
    shell:
        "dot -T{wildcards.fmt} < {input} > {output}"
