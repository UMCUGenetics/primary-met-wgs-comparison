
import os,json
import pandas as pd
import numpy as np

# example run

# Config
output_path=config["o"]
input_path=config["i"]

#snakemake --profile slurm --snakefile run_TEDs_pipeline.py --config i=/home/cog/fmartinez/scripts/resistance/data/control_hmf/samples_mechanisms_updated.json o=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/ -drop-metadata --latency-wait 20
#snakemake --profile slurm --snakefile run_TEDs_pipeline.py --config i=/home/cog/fmartinez/scripts/resistance/data/control_hmf/samples_untreated_updated.json o=/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/processed/resistance/ -drop-metadata --latency-wait 20 -np



with open(input_path,'r') as f: # load dictionary with groups of patients
    d_info=json.load(f)

ttypes,treatments=[],[]
for ttype in d_info:
    os.makedirs(f"{output_path}" + f"/dndscv/{ttype}/",exist_ok=True)
    for treatment in d_info[ttype]:
        # create list of patients
        file_input = f"{output_path}" + f"/dndscv/{ttype}/{treatment}.json"
        with open(file_input, 'w') as f:
            json.dump(d_info[ttype][treatment], f)
        ttypes.append(ttype)
        treatments.append(treatment)

# Rules
rule all:
    input:
        expand(f"{output_path}"+"/dndscv/{ttype}/{treatment}.dndscv.results.tsv.gz",zip,ttype=ttypes,treatment=treatments),

rule create_dataset_input_mutations:
    input:
        input_path=ancient(f"{output_path}" + "/dndscv/{ttype}/{treatment}.json")
    output:
        input_dndscv=f"{output_path}"+"/dndscv/{ttype}/{treatment}.dndscv.input.tsv.gz"
    params: cluster_memory = "64G"
    run:
        shell('set +eu '
        ' && PS1=dummy '
        ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate global && python /home/cog/fmartinez/scripts/paper_pancancer/get_variants_samples.py --patients_file {input.input_path} --treatment_name {wildcards.treatment} --output_file {output.input_dndscv}')

rule run_dndscv:
    input:
        input_dndscv=f"{output_path}"+"/dndscv/{ttype}/{treatment}.dndscv.input.tsv.gz",
    output:
        ouput_dndscv=f"{output_path}"+"/dndscv/{ttype}/{treatment}.dndscv.results.tsv.gz",
    params: cluster_memory = "64G"
    run:
        basepath=os.path.dirname(output.ouput_dndscv)
        shell('set +eu '
            ' && PS1=dummy '
            ' && . $(conda info --base)/etc/profile.d/conda.sh && conda activate dndscv && Rscript /home/cog/fmartinez/scripts/paper_pancancer/run_dndscv.R {input.input_dndscv} {basepath} {wildcards.treatment}')


