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

