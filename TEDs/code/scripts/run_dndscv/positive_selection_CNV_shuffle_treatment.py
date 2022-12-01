import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy import stats
import click
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import pybedtools
import json
import sys
from joblib import Parallel, delayed


path_exclude_full="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/analysis/results/immune_escape/LOH/positive_selection/chrx_and_chry.bed"
path_exclude_chry="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/analysis/results/immune_escape/LOH/positive_selection/chry.bed"


def load_cytoband():
    cytobands=pd.read_csv("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/analysis/results/immune_escape/LOH/positive_selection/cytoband_sizes.tsv",sep="\t")
    cytobands["chr"] = cytobands.apply(lambda row: row["chrom"][3:],axis=1)
    cytobands.set_index(["chr","arm"],inplace=True)
    return cytobands

def load_region_interest(chunks="100kb",chrx=False):
    chr_=""
    if chrx:
        chr_ = "_xchr"
    if chunks == "100kb":
        df_regions = pd.read_csv(f"/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/analysis/results/immune_escape/LOH/positive_selection/regions_kb{chr_}.tsv.gz",sep="\t")
    elif chunks == "10kb":
        df_regions = pd.read_csv(f"/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/analysis/results/immune_escape/LOH/positive_selection/regions_10kb{chr_}.tsv.gz",sep="\t")
    elif chunks == "1kb": # only recommended for really high sample sizes (>500 samples)
        df_regions = pd.read_csv(f"/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/analysis/results/immune_escape/LOH/positive_selection/regions_1kb{chr_}.tsv.gz",sep="\t")
    return df_regions


def count_unique(grp):
    return len(set(list(grp)))


def count_overlap(df_query,df_regions,simulations=False):
    # Given a dataframe with segments and one with regions, intersect them to count number of events
    # observed segments
    query=pybedtools.BedTool.from_dataframe(df_query[["chromosome","start","end","ID"]])

    if simulations: # perform some prepparations
        df_regions["chromosome"] = df_regions["chromosome"].astype(str)
        df_regions["chr"]="chr"+df_regions["chromosome"]
        regions=pybedtools.BedTool.from_dataframe(df_regions[["chr","start","end","name"]].rename(columns={"chr":"chromosome"}))
    else:
        regions = pybedtools.BedTool.from_dataframe(df_regions[["chromosome", "start", "end", "name"]])
    # Perform pybedtools intersection
    i=intersection=query.intersect(regions,wao=True).to_dataframe()
    i.columns=["chr_i","start_i_q","end_i_q","ID","chr_i_r","start_i_r","end_i_r","name","length"] # overlapping segments
    d_counts=intersection.groupby("name").agg({"ID":count_unique}).to_dict()["ID"] # this counts for every i-th simulation the number of simulated events across samples, ID is the sample identifier

    return d_counts



def get_chr_arm(row,cytobands):
    end = row["end"]
    start = row["start"]
    c = row["chromosome"]

    p_lenght = cytobands.loc[(c, "p")]["length"]
    q_lenght = cytobands.loc[(c, "q")]["length"]
    if start < p_lenght and end < p_lenght:
        return p_lenght
    if start < p_lenght and end > p_lenght:
        return np.nanmin([p_lenght, q_lenght])
    elif start > p_lenght:
        return q_lenght

    return np.nanmin([p_lenght, q_lenght])

def is_highly_focal(row):
    l = row["end"] - row["start"]
    if l < 3*10**6: # 3mb
        return True
    return False

def is_focal(row):
    l = row["end"] - row["start"]
    if l < row["chr_arm_lenght"]*0.75: # lower than 75% of arm lenght
        return True
    return False

def is_nonfocal(row):
    l = row["end"] - row["start"]
    if l > row["chr_arm_lenght"]*0.75: # greater than 75% of arm lenght
        return True
    return False

def is_imbalanced(row):
    ratio =  (row["minor_ploidy_integer"]) / (row["major_ploidy_integer"] + 10 ** -9)
    ratio_diff = row["majorAlleleCopyNumber"] - row["minorAlleleCopyNumber"]

    return (row["minor_ploidy_integer"] > 0) & (ratio < 0.9) & (ratio_diff > 0.4)

def simulate(a):

    i,observed,merge=a[0],a[1],a[2]
    if not(merge): # the number of segments is below <1000, operate normal
        output= observed.shuffle(genome="hg19", noOverlapping=True, excl=path_exclude,allowBeyondChromEnd=False ).sort().to_dataframe().rename(columns={"chrom": "chromosome"})
    else: # we create overlapping and then merge... only used for samples with a very high number of events (i.e., >1000) that cannot be randomized in a short time
        output= observed.shuffle(genome="hg19",excl=path_exclude,allowBeyondChromEnd=False, noOverlapping=False ).sort().merge(d=-1).to_dataframe().rename(
            columns={"chrom": "chromosome"})
    return i,output


def perform_randomizations(samples ,source_path="", type_analysis="loh", focal="nonfocal", n_iterations=10, dataset="HMF", chunk_size = "100kb", chrx=False):

    # load lenght of cytoband
    cytobands=load_cytoband()

    d_simulations,l_observed={},[]

    # create and empty simulations dictionary
    d_simulations = {i: [] for i in range(n_iterations)}

    # for each sample, load the somatic cnv variants
    for sample in tqdm(samples):
        if dataset == "HMF":
            try:
                df = pd.read_csv(f"{source_path}/{sample}/purple/{sample}.purple.cnv.somatic.tsv",sep= "\t")
            except FileNotFoundError:
                continue
            mean_ploidy = pd.read_csv(f"{source_path}/{sample}/purple/{sample}.purple.purity.tsv",sep= "\t")["ploidy"].values[0]
        else:
            df = pd.read_csv(f"{source_path}/{sample}-from-jar/purplesoft3.3/{sample}T.purple.cnv.somatic.tsv", sep="\t")
            mean_ploidy = pd.read_csv(f"{source_path}/{sample}-from-jar/purplesoft3.3/{sample}T.purple.purity.tsv", sep="\t")["ploidy"].values[0]

        # remove sexual chromosomes
        if chrx:
            df=df[(df["chromosome"]!="Y")]
        else:
            df=df[(df["chromosome"]!="X")&(df["chromosome"]!="Y")]


        if df.shape[0] ==0: # health check, next sample if there are no events
            continue

        # extract some useful measurements
        df["chr_arm_lenght"] = df.apply(lambda row: get_chr_arm(row,cytobands), axis=1)
        df["minor_ploidy_integer"] = df.apply(lambda row: round(row["minorAlleleCopyNumber"]),axis=1)
        df["major_ploidy_integer"] = df.apply(lambda row: round(row["majorAlleleCopyNumber"]), axis=1)

        # Type analysis
        if type_analysis == "loh":
            df["loh"] = df.apply(lambda row: (row["minorAlleleCopyNumber"] <= .3) & (row["majorAlleleCopyNumber"] >= 0.7), axis=1)
        elif type_analysis == "deepdel":
            df["deepdel"] = df.apply(lambda row: (row["copyNumber"]<0.5), axis=1)
        elif type_analysis == "amp":
            df["amp"] = df.apply(lambda row: row["copyNumber"] > float(mean_ploidy)+2.5, axis=1) # resistance arne
        elif type_analysis == "imbalance":
            df["imbalance"] = df.apply(lambda row: is_imbalanced(row), axis=1)
        else:
            # not a valid id, stop
            sys.exit(1)

        # Genomic scale
        if focal=="focal": #  if segment length < 0.5 chr arm
            df["focal"] = df.apply(lambda row: is_focal(row), axis=1)
            observed = df[(df[type_analysis])&(df["focal"])][["chromosome", "start", "end"]]

        elif focal =="hfocal":  #  if segment length < 3 Mb
            df["hfocal"] = df.apply(lambda row: is_highly_focal(row), axis=1)
            observed = df[(df[type_analysis]) & (df["hfocal"])][["chromosome", "start", "end"]]

        else: # greater than 0.5
            df["nonfocal"] = df.apply(lambda row: is_nonfocal(row), axis=1)
            observed = df[(df[type_analysis])& (df["nonfocal"])][["chromosome", "start", "end"]]

        if observed.shape[0] ==0: # health check, next sample if there are no events
            continue

        # Check the number of segments, for those cases with very high number the simulations would never converge, the need to set up the overlapping = True flag
        merge=False
        if observed.shape[0] > 10000 or (mean_ploidy < 1.2):
            merge=True
        print (sample,observed.shape[0],mean_ploidy)
        # define segment lenghts
        observed["size"] = observed["end"] - observed["start"]
        observed.sort_values("size",ascending=False,inplace=True)

        # pybed object
        pybed_observed = pybedtools.BedTool.from_dataframe(observed)
        # set the sample ID
        observed["ID"] = sample
        l_observed.append(observed)
        # Invoke the multitreadhing function myfun to perform the n_iterations randomizations, it will shuffle the observed and return the randomizations
        results = Parallel(n_jobs=16,  backend="threading")(
            map(delayed(simulate), [(i,pybed_observed,merge) for i in range(n_iterations)]))
        # Returns the simulated segments
        for i,r in results:
            r["ID"] = sample # assigngs the sample id
            d_simulations[i].append(r) # save into the i-th simulation


    if len(l_observed)>0:

        return d_simulations,pd.concat(l_observed)
    else:
        return d_simulations,pd.DataFrame([])



def get_statistical_test(d_simulations,df_observed,df_query_regions):

    # count number of observed events per bin across samples
    d_counts_observed=count_overlap(df_observed,df_query_regions)

    # Regions are the KB bins
    regions = list(set(df_query_regions["name"].unique()))

    # Counts of simulated events per i-th simulation across amples
    d_dist = {}
    for sim_n,l_sim in tqdm((d_simulations.items())): # for every i-th simulation (10 by default)
        df_sim = pd.concat(l_sim)
        d_counts_simulated=count_overlap(df_sim,df_query_regions,simulations=True) # counts number of events simulated per i-th simulation across samples
        for region in regions:
            if not(region in d_dist):
                d_dist[region] = [] # Initiate
            if region in d_counts_simulated:
                d_dist[region].append(d_counts_simulated[region])
            else:
                d_dist[region].append(0)
    # compare observed with expected
    l = []
    mean_global = np.nanmean([np.nanmean(d_dist[region]) for region in d_dist]) # calculate the global average of events across all regions and all simulations

    for region in regions:
        if region in d_counts_observed:
            observed = d_counts_observed[region]
        else:
            observed = 0
        pvalue_emp = sum(np.array(d_dist[region]) > observed) / len(d_dist[region]) # empirical p-value, only valid for a high number of simulations
        odds_ratio = observed / np.nanmean(d_dist[region]) # local odds ratio
        odds_ratio_global = observed / mean_global # global odds ratio
        q1,q2,q3 = np.percentile( d_dist[region],[25,50,75]) # local distribution of simulations
        if mean_global > observed: # we are only interested on the right tail, assign a non-significaint pvalue to the left tail
            pvalue_ana, pvalue_ana_global = np.random.uniform(0.1,1.0,1)[0], np.random.uniform(0.1,1.0,1)[0] # left-tail we are not interested in depletion
        else:
            _, pvalue_ana = stats.power_divergence(f_obs=[observed,  df_observed.shape[0] - observed],
                                               f_exp=[np.nanmean(d_dist[region]), df_observed.shape[0] - np.nanmean(d_dist[region])],
                                               lambda_="log-likelihood") # g-test
            _, pvalue_ana_global = stats.power_divergence(f_obs=[observed, df_observed.shape[0] - observed],
                                                   f_exp=[mean_global,
                                                          df_observed.shape[0] - mean_global],
                                                   lambda_="log-likelihood")
        l.append([region, observed, np.nanmean(d_dist[region]),q1,q2,q3, pvalue_emp, pvalue_ana, odds_ratio, pvalue_ana_global,odds_ratio_global,mean_global])

    o = pd.DataFrame(l, columns=["region", "n_observed", "n_mean_simulated","Q1_simulated","median_simulated","Q3_simulated", "pvalue_emp_local", "pvalue_ana_local",
                                     "odds_ratio_local", "pvalue_ana_global", "odds_ratio_global","mean_global"])
    o["log2_odds_ratio_global"] = np.log2(o["odds_ratio_global"]+.001) # calculate the log2 odds_ratio
    return o


@click.command()
@click.option('--samples_file',
              type=click.Path(exists=True),
              help="Input dir with summary data",
              required=True)
@click.option('--tumor_type',
              help="tumor type",
              required=True)
@click.option('--type_analysis',
              help="loh, amp, deepdel",
              default="loh",
              required=False)
@click.option('--focal' ,
              help="focal or ignore length",
              default="nonfocal",
              required=False)
@click.option('--output_file',
              type=click.Path(),
              help="Output file)",
              required=True)
@click.option('--dataset',
              type=click.STRING,
              help="Name of the dataset [HMF or PCAWG]",
              default="HMF",
              required=False)
@click.option('--chunk_size',
              type=click.STRING,
              help="Size of chunk [100kb, 10kb, 1kb]",
              default="100kb",
              required=False)
@click.option('--chrx',
              is_flag=True,
              help="Whether include the chrx in the analysis, LOH cannot be used including chrx",
              default=False,
              required=False)
def run(samples_file,tumor_type,type_analysis,focal,output_file,dataset,chunk_size,chrx):
    if chrx and type_analysis =="loh":
        print  ("LOH cannot be used with chrx flag")
        sys.exit(1)
    global path_exclude
    path_exclude = path_exclude_full
    if chrx:
        path_exclude = path_exclude_chry


    # set temp dir, important for randomizations
    pybedtools.helpers.set_tempdir("/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/tmp/pybed2/")

    # get somatic alterations
    if dataset=="HMF":
        path_alterations="/hpc/cuppen/shared_resources/HMF_data/DR-104-update4/somatics/"
    else:
        path_alterations="/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/"


    # group of samples, by cancer type
    with open(samples_file,'r') as f:
        samples_ttype=json.load(f)

    # tumor_type

    tumor_type = tumor_type.replace("_"," ")
    # check that the number of samples is greater than 15, otherwise returns an empty dataframe
    samples = list(samples_ttype)
    if len(samples) < 15:
        pd.DataFrame([]).to_csv(output_file, sep="\t", index=False, compression="gzip")
        return

    d_simulations, df_observed = perform_randomizations(samples,
                                                                source_path=path_alterations,
                                                                type_analysis=type_analysis, focal=focal,dataset=dataset, chunk_size=chunk_size, chrx=chrx)

    if df_observed.shape[0] == 0: # if there are no events across the >=15
        df_observed.to_csv(output_file, sep="\t", index=False, compression="gzip")
        return

    # Load the bins, by default are Kilobases of the automsome
    regions_interest=load_region_interest(chunk_size, chrx)

    # Perform the statistical comparison
    df = get_statistical_test(d_simulations, df_observed, regions_interest)

    # Adjust by FDR
    df["q_value_emp"] = fdrcorrection0(df["pvalue_emp_local"].values)[1]
    df.fillna(value={"pvalue_ana_local":0.5,"pvalue_ana_global":0.5},inplace=True)
    df["q_value_ana_local"] = fdrcorrection0([v for v in  df["pvalue_ana_local"].values if np.isfinite(v)])[1]
    df["q_value_ana_global"] = fdrcorrection0([v for v in  df["pvalue_ana_global"].values if np.isfinite(v) ])[1]
    df["cancer_type"] = tumor_type
    # Save the output
    df.to_csv(output_file,sep="\t",index=False,compression="gzip")

if __name__ == '__main__':
    run()