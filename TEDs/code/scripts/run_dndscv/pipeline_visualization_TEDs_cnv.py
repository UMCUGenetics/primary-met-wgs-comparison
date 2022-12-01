import pandas as pd
import os
# Config
output_path=config["o"]
input_path=config["i"]
with open(input_path,'r') as f:
    d_info=json.load(f)

ttypes,treatments=[],[]

for ttype in d_info:
    for treatment in d_info[ttype]:
        # create list of patients
        file_input = f"{output_path}" + f"{ttype}__{treatment}.json"
        with open(file_input, 'w') as f:
            json.dump(d_info[ttype][treatment], f)
        ttypes.append(ttype)
        treatments.append(treatment)

labels=[]
for l,v in zip(ttypes,treatments):
    labels.append(l+"__"+v)



names =["ignore_amp_hfocal", "ignore_amp_focal"]

# Rules
rule all:
    input:
        expand(f"{output_path}" + "/{conf}/{name}.tsv.gz",name=labels,conf=names)

rule run_pos_selection_LOH:
    input:
        input_path = ancient(f"{output_path}" + "/{name}.json")
    threads: 8
    params: cluster_memory = "16G"
    output:
        output=f"{output_path}" + "/{conf}/{name}.tsv.gz"
    run:
        wgd_status,type_analysis,focal=wildcards.conf.split("_")
        flag=""
        p=f"{output_path}" + f"/{wildcards.conf}/"
        if not(os.path.exists(p)):
            os.mkdir(p)
        cohort="HMF"
        if "untreated" in wildcards.name:
            cohort="PCAWG"
        shell('set +eu '
              ' && PS1=dummy '
              ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate global && ' \
              'python /home/cog/fmartinez/scripts/paper_pancancer/positive_selection_CNV_shuffle_treatment.py --samples_file {input.input_path} ' \
              '--tumor_type "{wildcards.name}" --type_analysis {type_analysis} {flag} --focal {focal}  --output_file {output.output} --chrx --dataset {cohort}')

