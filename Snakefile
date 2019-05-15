import pandas as pd

configfile: "wbuild.yaml"

htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"

include: ".wBuild/wBuild.snakefile"  # Has to be here in order to update the config with the new variables
htmlOutputPath = config["htmlOutputPath"]  if (config["htmlOutputPath"] != None) else "Output/html"

rule all:
    input: rules.Index.output, htmlOutputPath + "/readme.html"
    output: touch(htmlOutputPath + "/all.done")
