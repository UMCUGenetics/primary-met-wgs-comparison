import click
import pandas as pd
import os
import sys
from tqdm import tqdm
import json
import numpy as np

def load_patients(file):
    with open(file) as f:
        patients=json.load(f)
    return patients

@click.command()
@click.option('--patients_file',
              type=click.Path(exists=True),
              help="Input .json with patient IDs",
              required=True)
@click.option('--treatment_name',
              type=click.STRING,
              help="name of ttype",
              required=True)
@click.option('--output_file',
              type=click.Path(),
              help="Output file",
              required=True)
              
def process_mutations(patients_file, treatment_name,output_file):
    patients = load_patients(patients_file)
    # patients PCAWG
    pats_pcawg = set([s for s in patients if s[0:2] == "DO"])
    # patients HMF
    pats_hmf= set([s for s in patients if s[0:2] != "DO"])

    # process Hartwig samples
    for sample in tqdm(pats_hmf):
        if not(os.path.exists(f"/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/vaf/{sample}/{sample}.vaf_info.tsv.gz")):
            continue
        df = pd.read_csv(
            f"/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/vaf/{sample}/{sample}.vaf_info.tsv.gz",
            sep="\t")
        df = df[df["CHROM"] != "chrMT"]
        df["sample"] = sample
        df.rename(columns={"sample": "sampleID", "CHROM": "chr", "POS": "pos", "REF": "ref", "ALT": "mut"})[
            ["sampleID", "chr", "pos", "ref", "mut"]].drop_duplicates().to_csv(output_file, sep="\t", index=False, compression="gzip", mode='a')

    # process PCAWG samples
    for sample in tqdm(pats_pcawg):
        if not(os.path.exists(f"/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/vaf/{sample}/{sample}.vaf_info.tsv.gz")):
            continue
        df = pd.read_csv(
            f"/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/vaf/{sample}/{sample}.vaf_info.tsv.gz",
            sep="\t")
        df = df[df["CHROM"] != "chrMT"]
        df["sample"] = sample
        df.rename(columns={"sample": "sampleID", "CHROM": "chr", "POS": "pos", "REF": "ref", "ALT": "mut"})[
            ["sampleID", "chr", "pos", "ref", "mut"]].drop_duplicates().to_csv(output_file, sep="\t", index=False,
                                                                               compression="gzip", mode='a')



if __name__ == '__main__':
    process_mutations()
