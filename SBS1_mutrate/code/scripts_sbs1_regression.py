import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as st
import pandas as pd
pd.options.mode.chained_assignment = None  # supress warnings
default='warn'
from matplotlib import gridspec
from statannot import add_stat_annotation
pd.options.display.max_columns=100

ttypes=["Breast carcinoma", "Glioblastoma multiforme", "Colorectal carcinoma", "Esophageal carcinoma",
                          "Stomach carcinoma", "Cholangiocarcinoma", "Hepatocellular carcinoma", "Pancreas carcinoma",
                          "Pancreas neuroendocrine", "Cervical carcinoma", 'Ovarian serous adenocarcinoma', "Uterus carcinoma",
                          "Upper respiratory tract carcinoma", "Kidney renal clear cell carcinoma", "Lung adenocarcinoma","Lung squamous cell carcinoma",
                          "Diffuse large B-cell lymphoma", "Prostate carcinoma", "Skin melanoma", "Leiomyosarcoma",
                          "Liposarcoma", "Thyroid carcinoma", "Bladder urothelial carcinoma"]

ttypes_selected=["Breast carcinoma", "Prostate carcinoma","Kidney renal clear cell carcinoma", "Thyroid carcinoma",
                          "Colorectal carcinoma","Ovarian serous adenocarcinoma"]

def get_calculate_mutations_year(slope_prim,intercept_prim,slope_met,intercept_met,years_range): # this function calculate 
    raw_val,fold_change,prop=[],[],[]
    for year in years_range:
        M=(year*slope_met + intercept_met) # calculate the expected SBS1 burden in met
        P=(year*slope_prim + intercept_prim) # calculate the expected SBS1 burden in met
        diff =  M - P
        raw_val.append(diff)
        fold_change.append( (M / P))
        prop.append(diff / P)
    return raw_val,fold_change,prop
    
def get_offset_ages_met(slope_primary,intercept_primary,y_values,x_real,ttype,sbs,name): # calculat the number of required years (offset) to observed the metasttic SBS1 burden with the same primary SBS1 mut. rate
    l=[]
    for y in y_values:
        x_pred = (y - intercept_primary) / slope_primary
        l.append(x_real - x_pred)
    return l
    
def plot_regression(ttype_name, ttype,df_data,lim_max_sbs1=5000,ylim=2000,column="sbs1_count",name="_",sbs="SBS1",title2="",plot_residuals=False):
    fig,ax = plt.subplots(figsize=(3.5,4))
    
    met=df_data[(df_data["cancer_type_code"]==ttype)&(df_data[column]<5000)&(df_data["cohort"]=="Hartwig")&(np.isfinite(df_data[column]))&(np.isfinite(df_data["age"]))]
    primary=df_data[(df_data["cancer_type_code"]==ttype)&(df_data[column]<5000)&(df_data["cohort"]=="PCAWG")&(np.isfinite(df_data[column]))&(np.isfinite(df_data["age"]))]
    ax.set_xlim(0,100)
    max_v=np.max([list(met[column])+list(primary[column])])+100
    ax.set_ylim(0,max_v)
    ax.scatter(x=met["age"],y=met[column],color="#af8dc3",alpha=0.8,lw=0.25,edgecolor="black")
    ax.scatter(x=primary["age"],y=primary[column],color="#fc8d59",alpha=0.8,lw=0.25,edgecolor="black")
    v_met,v_prim=[],[]
    l_met,l_prim=[],[]
    xseq = np.linspace(0, 100, num=100)
    for N in range(100):
        s=met.sample(frac=0.75,random_state=N)
        s1=primary.sample(frac=0.75,random_state=N)
        y_raw,y_raw1=s[column],s1[column]
        slope,intercept,rvalue,pvalue_r,stde=st.linregress(s["age"], y=y_raw)
        slope1,intercept1,rvalue1,pvalue_r1,stde=st.linregress(s1["age"], y=y_raw1)
        l_met.append([slope,intercept,rvalue,pvalue_r])
        l_prim.append([slope1,intercept1,rvalue1,pvalue_r1])
        
        m=[x*slope + intercept for x in xseq]
        p=[x*slope1 + intercept1 for x in xseq]
        v_met.append(m)
        v_prim.append(p)
    # Metastatic
    ## mean
    slope_met,intercept_met,rvalue_res_met,pvalue_res_met=sorted(l_met)[50]
    y_pred=xseq*slope_met + intercept_met
    ax.plot(xseq,y_pred,color="#5e3c99",linewidth=3,alpha=1.0)
    # background 
    y_min,y_max=[],[]
    values_met_T=list(np.array(v_met).T)
    for i,_ in enumerate(xseq):
        tm_=sorted(values_met_T[i])
        y_min.append(tm_[1])
        y_max.append(tm_[99])
    plt.fill_between(xseq,y_min,y_max,color="#5e3c99",alpha=0.25)
    # Primary
    ## mean
    slope_primary,intercept_primary,rvalue_res_prim,pvalue_res_prim=sorted(l_prim)[50]
    y_pred=xseq*slope_primary + intercept_primary
    ax.plot(xseq,y_pred,color="#e66101",linewidth=3,alpha=1.0)
    # background
    y_min,y_max=[],[]
    values_prim_T=list(np.array(v_prim).T)
    for i,_ in enumerate(xseq):
        tm_=sorted(values_prim_T[i])
        y_min.append(tm_[1])
        y_max.append(tm_[99])
    plt.fill_between(xseq,y_min,y_max,color="#e66101",alpha=0.15)
    ax.set_title(ttype_name+"\n"+title2,fontsize=16)
    ax.set_xlabel("Age at biopsy",fontsize=16)
    ax.set_ylabel(f"{sbs} mutations",fontsize=16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=16)
    # Annotate number of mutations per year
    max_age=int(np.max([list(met["age"])+list(primary["age"])]))
    min_age=int(np.min([list(met["age"])+list(primary["age"])]))
    mut_diff,fold_change,prop=get_calculate_mutations_year(slope_primary,intercept_primary,slope_met, intercept_met,xseq[40:80+1])
    # compute differences in residuals distributions
    # get top most mutated, compared to lowly
    x_raw=list(met["age"])
    y_raw=list(met[column])
    residuals_met = [y_raw[i] - (value*slope_primary + intercept_primary) for i,value in enumerate(x_raw)]
    x_raw=list(primary["age"])
    y_raw=list(primary[column])
    residuals_primary = [y_raw[i] - (value*slope_primary + intercept_primary) for i,value in enumerate(x_raw)]
    try:
        pvalue=st.mannwhitneyu(residuals_met,residuals_primary,alternative="greater")
    except:
        return []
    if pvalue[1] < 0.01 and intercept_met > intercept_primary and rvalue_res_prim>0.1 and rvalue_res_met > 0.1:
        ax.annotate(xy=(2,max_v-200),s=f"x{np.nanmean(fold_change):1.2f} fold change (+{np.nanmean(mut_diff):1.0f} ± {np.nanstd(mut_diff):1.0f} muts)",fontsize=16)
        ax.annotate(xy=(25,max_v-400),s=f"P-value={pvalue[1]:.2e}",fontsize=16)
    plt.savefig(f'../results/figures/{sbs}_regressions/{sbs}_prim_vs_met_{ttype}{name.replace("/","__")}.pdf', dpi=800,bbox_inches="tight")
    if plot_residuals:
        plot_residuals(ttype_name,ttype,residuals_met,residuals_primary,name,sbs)
    return slope_primary,intercept_primary, slope_met,  intercept_met, np.nanmean(mut_diff), np.nanstd(mut_diff), np.nanmean(fold_change),np.nanstd(fold_change), len(residuals_met),len(primary), pvalue[1],pvalue_res_met,rvalue_res_met, pvalue_res_prim,rvalue_res_prim

def plot_regression_combined(df_data,lim_max_sbs1=5000,ylim=2000,column="sbs1_count",name="_",sbs="SBS1",ttypes=ttypes,col_lim=5000,title2="", plot_pearsonr=False):
    fig,ax = plt.subplots(figsize=(30,12))
    gs = gridspec.GridSpec(figure=fig, ncols=8, nrows=3)
    gs.update(hspace=0.5, wspace=0.5)
    d_diff={}
    l_data=[]
    for i,ttype in enumerate(ttypes): # for every cancer type
        ax = plt.subplot(gs[i])
        # select patients fom the cancer type, with annotated age
        met=df_data[(df_data["cancer_type"]==ttype)&(df_data[column]<col_lim)&(df_data["cohort"]=="Hartwig")&(np.isfinite(df_data[column]))&(np.isfinite(df_data["age"]))]
        primary=df_data[(df_data["cancer_type"]==ttype)&(df_data[column]<col_lim)&(df_data["cohort"]=="PCAWG")&(np.isfinite(df_data[column]))&(np.isfinite(df_data["age"]))]
        # if there are less than 5 samples, continue with next cancer type
        if met.shape[0] <5 or primary.shape[0] < 5:
            ax.set_title(ttype+"\n"+title2,fontsize=14)
            continue
        # plot the raw data 
        ax.scatter(x=met["age"],y=met[column],color="#af8dc3",alpha=0.8,lw=0.25,edgecolor="black")
        ax.scatter(x=primary["age"],y=primary[column],color="#fc8d59",alpha=0.8,lw=0.25,edgecolor="black")
        ax.set_xlim(0,100)
        max_v=np.max([list(met[column])+list(primary[column])])+100
        ax.set_ylim(0,max_v)
        # now calculate the regression lines 
        v_met,v_prim=[],[]
        l_met,l_prim=[],[]
        xseq = np.linspace(0, 100, num=100) # this is the range of ages that we will explore
        for N in range(100): # bootstraps
            # select 75% of the primary and met cohort
            s=met.sample(frac=0.75,random_state=N)
            s1=primary.sample(frac=0.75,random_state=N)
            y_raw,y_raw1=s[column],s1[column]
            # calculate the regression parameters for this boostrap
            slope,intercept,rvalue,pvalue_r,stde=st.linregress(s["age"], y=y_raw)
            slope1,intercept1,rvalue1,pvalue_r1,stde=st.linregress(s1["age"], y=y_raw1)
            l_met.append([slope,intercept,rvalue,pvalue_r])
            l_prim.append([slope1,intercept1,rvalue1,pvalue_r1])
            # for the range of ages (0-100), predict the expected SBS1 burden
            m=[x*slope + intercept for x in xseq]
            p=[x*slope1 + intercept1 for x in xseq]
            v_met.append(m)
            v_prim.append(p)
        
        # Metastatic
        ## the the median regressoin (based on the slope). This will be the representative
        slope_met,intercept_met,rvalue_res_met,pvalue_res_met=sorted(l_met)[50]
        # predict the expected SBS1 burden based on this regression
        y_pred_met=xseq*slope_met + intercept_met
        ax.plot(xseq,y_pred_met,color="#5e3c99",linewidth=3,alpha=1.0)
        # now plot the confidence internvals
        y_min,y_max=[],[]
        values_met_T=list(np.array(v_met).T)
        for i,_ in enumerate(xseq):
            tm_=sorted(values_met_T[i])
            y_min.append(tm_[1])
            y_max.append(tm_[99])
        plt.fill_between(xseq,y_min,y_max,color="#5e3c99",alpha=0.25)
        # Primary
        ## the the median regressoin (based on the slope). This will be the representative
        slope_primary,intercept_primary,rvalue_res_prim,pvalue_res_prim=sorted(l_prim)[50]
        y_pred_primary=xseq*slope_primary + intercept_primary
        ax.plot(xseq,y_pred_primary,color="#e66101",linewidth=3,alpha=1.0)
        # background
        y_min,y_max=[],[]
        values_prim_T=list(np.array(v_prim).T)
        for i,_ in enumerate(xseq):
            tm_=sorted(values_prim_T[i])
            y_min.append(tm_[1])
            y_max.append(tm_[99])
        plt.fill_between(xseq,y_min,y_max,color="#e66101",alpha=0.15)
        # adjust visualization parameters
        ax.set_title(ttype+"\n"+title2,fontsize=14)
        ax.set_xlabel("Age at biopsy",fontsize=14)
        ax.set_ylabel(f"{sbs} mutations",fontsize=14)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', which='major', labelsize=12)
        max_age=int(np.max([list(met["age"])+list(primary["age"])]))
        min_age=int(np.min([list(met["age"])+list(primary["age"])]))
        # Calculate fold_change and raw differences across the interval of 40 to 80 years old
        mut_diff,fold_change,prop=get_calculate_mutations_year(slope_primary,intercept_primary,slope_met, intercept_met,xseq[40:80+1])
        # calculate offsets (number of years required to reach met. SBS1 burden)
        offsets=get_offset_ages_met(slope_primary,intercept_primary,met[column],met["age"],ttype,sbs,name)
        # Compute significance of different distributions, for that we will predict based on the primary regression the residuals of both metastatic and primary samples based on the primary regressions
        x_raw=list(met["age"])
        y_raw=list(met[column])
        residuals_met = [y_raw[i] - (value*slope_primary + intercept_primary) for i,value in enumerate(x_raw)]
        x_raw=list(primary["age"])
        y_raw=list(primary[column])
        residuals_primary = [y_raw[i] - (value*slope_primary + intercept_primary) for i,value in enumerate(x_raw)]
        try:
            pvalue_distributions=st.mannwhitneyu(residuals_met,residuals_primary,alternative="greater")
        except: # cases with low numbers may fail
            continue
        if pvalue_distributions[1] < 0.01 and intercept_met > intercept_primary and rvalue_res_prim>0.1 and rvalue_res_met > 0.1: # significance as defined in methods
            ax.annotate(xy=(2,max_v-100),s=f"x{np.nanmean(fold_change):1.2f} fold change",fontsize=12)
            ax.annotate(xy=(2,max_v-200),s=f"+{np.nanmean(mut_diff):1.0f} ± {np.nanstd(mut_diff):1.0f} muts",fontsize=12)
            ax.annotate(xy=(25,max_v-300),s=f"P-value={pvalue_distributions[1]:.2e}",fontsize=12)
        if plot_pearsonr:
            ax.annotate(xy=(50,100),s=f"Rmet={rvalue_res_met:1.2f}",fontsize=12)
            ax.annotate(xy=(50,10),s=f"Rparim={rvalue_res_prim:1.2f}",fontsize=12)
            
        d_diff[ttype]=(mut_diff,y_pred_primary,y_pred_met)
        ttype_code=met["cancer_type_code"].values[0]
        l_data.append([ttype,ttype_code,slope_primary,intercept_primary, slope_met,  intercept_met, np.nanmean(mut_diff), np.nanstd(mut_diff), 
                       np.nanmean(fold_change), np.nanstd(fold_change),len(residuals_met), len(primary), 
                       pvalue_distributions[1],pvalue_res_met,rvalue_res_met, pvalue_res_prim,rvalue_res_prim,np.nanmedian(offsets),np.nanstd(offsets)])
    plt.savefig(f'../results/figures/{sbs}_regressions/combined_{sbs}_prim_vs_met_{name}.pdf', dpi=800,bbox_inches="tight")
    df_stats = pd.DataFrame(l_data,columns=["cancer_type","cancer_type_code","slope_primary","intercept_primary","slope_met","intercept_met",
                                  "mean_diff_residual","std_diff_residual","fold_change_met","std_fold_change_met","n_met","n_prim","pvalue",
                                  "pvalue_res_met","rvalue_res_met", "pvalue_res_prim","rvalue_res_prim","median_offset_ages","std_offset_ages"])
    return d_diff,df_stats
 
        
    

# Function that enables the inclusion and visualization of samples from other cohorts
def plot_regression_with_independent(ttype_name,ttype,df_data,primary_df_independent,met_df_independent,lim_max_sbs1=5000,ylim=2000,column="sbs1_count",name="_",sbs="SBS1",title2="",join=True,plot_regression_met=False,plot_regression_primary=False,exome=False,plot_stats=True):
    fig,ax = plt.subplots(figsize=(4,4))
    plot_stats = plot_stats or exome
    met=df_data[(df_data["cancer_type_code"]==ttype)&(df_data[column]<5000)&(df_data["cohort"]=="Hartwig")&(np.isfinite(df_data[column]))&(np.isfinite(df_data["age"]))]
    primary=df_data[(df_data["cancer_type_code"]==ttype)&(df_data[column]<5000)&(df_data["cohort"]=="PCAWG")&(np.isfinite(df_data[column]))&(np.isfinite(df_data["age"]))]
    ax.set_xlim(0,100)
    if not(exome):
        max_v=np.max([list(met[column])+list(primary[column])])+100
    else:
        max_v=np.max([list(met[column])+list(primary[column])])+5
    ax.set_ylim(0,max_v)
    ax.scatter(x=met["age"],y=met[column],color="#af8dc3",alpha=0.8,lw=0.25,edgecolor="black")
    ax.scatter(x=primary["age"],y=primary[column],color="#fc8d59",alpha=0.8,lw=0.25,edgecolor="black")
    v_met,v_prim=[],[]
    l_met,l_prim=[],[]
    xseq = np.linspace(0, 100, num=100)
    for N in range(100):
        s=met.sample(frac=0.75,random_state=N)
        s1=primary.sample(frac=0.75,random_state=N)
        y_raw,y_raw1=s[column],s1[column]
        slope,intercept,rvalue,pvalue_r,stde=st.linregress(s["age"], y=y_raw)
        slope1,intercept1,rvalue1,pvalue_r1,stde=st.linregress(s1["age"], y=y_raw1)
        l_met.append([slope,intercept,rvalue,pvalue_r])
        l_prim.append([slope1,intercept1,rvalue1,pvalue_r1])
        
        m=[x*slope + intercept for x in xseq]
        p=[x*slope1 + intercept1 for x in xseq]
        v_met.append(m)
        v_prim.append(p)
    # Metastatic
    ## mean
    slope_met,intercept_met,rvalue_res_met,pvalue_res_met=sorted(l_met)[50]
    y_pred=xseq*slope_met + intercept_met
    ax.plot(xseq,y_pred,color="#5e3c99",linewidth=3,alpha=1.0)
    # background 
    y_min,y_max=[],[]
    values_met_T=list(np.array(v_met).T)
    for i,_ in enumerate(xseq):
        tm_=sorted(values_met_T[i])
        y_min.append(tm_[1])
        y_max.append(tm_[99])
    plt.fill_between(xseq,y_min,y_max,color="#5e3c99",alpha=0.25)
    # Primary
    ## mean
    slope_primary,intercept_primary,rvalue_res_prim,pvalue_res_prim=sorted(l_prim)[50]
    y_pred=xseq*slope_primary + intercept_primary
    ax.plot(xseq,y_pred,color="#e66101",linewidth=3,alpha=1.0)
    # background
    y_min,y_max=[],[]
    values_prim_T=list(np.array(v_prim).T)
    for i,_ in enumerate(xseq):
        tm_=sorted(values_prim_T[i])
        y_min.append(tm_[1])
        y_max.append(tm_[99])
    plt.fill_between(xseq,y_min,y_max,color="#e66101",alpha=0.15)
    ax.set_title(ttype_name+"\n"+title2,fontsize=16)
    ax.set_xlabel("Age at biopsy",fontsize=16)
    ax.set_ylabel(f"{sbs} mutations",fontsize=16)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=14)
    # Annotate number of mutations per year
    max_age=int(np.max([list(met["age"])+list(primary["age"])]))
    min_age=int(np.min([list(met["age"])+list(primary["age"])]))
    mut_diff,fold_change,prop=get_calculate_mutations_year(slope_primary,intercept_primary,slope_met, intercept_met,xseq[40:80+1])
    
    # compute differences in residuals distributions
    # get top most mutated, compared to lowly
    x_raw=list(met["age"])
    y_raw=list(met[column])
    residuals_met = [y_raw[i] - (value*slope_primary + intercept_primary) for i,value in enumerate(x_raw)]
    x_raw=list(primary["age"])
    y_raw=list(primary[column])
    residuals_primary = [y_raw[i] - (value*slope_primary + intercept_primary) for i,value in enumerate(x_raw)]
    try:
        pvalue=st.mannwhitneyu(residuals_met,residuals_primary,alternative="greater")
    except:
        return []
    if plot_stats and (pvalue[1] < 0.01 and intercept_met > intercept_primary):
        ax.annotate(xy=(2,max_v-5),s=f"x{np.nanmean(fold_change):1.2f} fold change (+{np.nanmean(mut_diff):1.0f} ± {np.nanstd(mut_diff):1.0f} muts)",fontsize=16)
        ax.annotate(xy=(25,max_v-10),s=f"P-value={pvalue[1]:.2e}",fontsize=16)

    # plot the points of independent study   
    ax.scatter(x=met_df_independent["age"],y=met_df_independent[column],color="#af8dc3",alpha=1.,lw=0.75,edgecolor="black",marker="s",s=100)
    ax.scatter(x=primary_df_independent["age"],y=primary_df_independent[column],color="#fc8d59",alpha=1.0,lw=0.75,edgecolor="black",marker="s",s=100)
    if (join):
        xs=list(zip(primary_df_independent["age"],met_df_independent["age"]))
        ys=list(zip(primary_df_independent[column],met_df_independent[column]))
        for i in range(len(xs)):
            ax.plot(xs[i],ys[i],color="black")
    ax.set_xlim(0,100)
    if plot_regression_met:
        sns.regplot(data=met_df_independent,x="age",y=column,color="#af8dc3",ax=ax,line_kws={"lw":2,"alpha":1.0},truncate=False,robust=True)
    if plot_regression_primary:
        sns.regplot(data=primary_df_independent,x="age",y=column,color="#fc8d59",ax=ax,line_kws={"lw":4,"alpha":1.0},truncate=False,robust=True)
    
    ## Calculate residuals independent met with met regression
    x_raw=list(met_df_independent["age"])
    y_raw=list(met_df_independent[column])
    residuals_met_indepenent = [y_raw[i] - (value*slope_met + intercept_met) for i,value in enumerate(x_raw)]
    x_raw=list(met["age"])
    y_raw=list(met[column])
    residuals_met_depenent = [y_raw[i] - (value*slope_met + intercept_met) for i,value in enumerate(x_raw)]
    pvalue_=st.mannwhitneyu(residuals_met_indepenent,residuals_met_depenent,alternative="two-sided")
    print ("Pvalue met ",pvalue_,np.mean(residuals_met_indepenent),np.mean(residuals_met_depenent))
    
    x_raw=list(primary_df_independent["age"])
    y_raw=list(primary_df_independent[column])
    residuals_pr_indepenent = [y_raw[i] - (value*slope_primary + intercept_primary) for i,value in enumerate(x_raw)]
    x_raw=list(primary["age"])
    y_raw=list(primary[column])
    residuals_pr_depenent = [y_raw[i] - (value*slope_primary + intercept_primary) for i,value in enumerate(x_raw)]
    pvalue_=st.mannwhitneyu(residuals_pr_indepenent,residuals_pr_depenent,alternative="two-sided")
    print ("Pvalue primary ",pvalue_,np.mean(residuals_pr_indepenent),np.mean(residuals_pr_depenent))
    ax.set_ylabel(f"{sbs} mutations",fontsize=16)
    plt.savefig(f'../results/figures/{sbs}_regressions/{sbs}_prim_vs_met_{name}.pdf', dpi=800,bbox_inches="tight")

    return slope_primary,intercept_primary, slope_met,  intercept_met, np.nanmean(mut_diff), np.nanstd(mut_diff), np.nanmean(fold_change), np.nanstd(fold_change), len(residuals_met), len(primary), float(pvalue[1]),pvalue_res_met,rvalue_res_met, pvalue_res_prim,rvalue_res_prim

def ci95(grp):
    return np.nanpercentile(grp,95)
def ci5(grp):
    return np.nanpercentile(grp,5)

def cell_division_rate_increase(df_stats,df_data,xerr=True,title="",column_y="sbs1_per_year",column_x="fold_change_met",std_x="std_fold_change_met",xlabel="Mean difference SBS1 burden (met - prim.) \n (40-80 years)",xlim=(-250,400),name="",type_x="raw",pearson=True,thr_res=0.1):

    v=df_data.groupby(["cancer_type","cohort"],as_index=False).agg(median_sbs1_y=(column_y,np.nanmean),sbs1_95=(column_y,ci95),sbs1_5=(column_y,ci5))
    pm=v[v["cohort"]=="PCAWG"]
    pm_s=pm.merge(df_stats[(df_stats["rvalue_res_met"]>thr_res)&(df_stats["rvalue_res_prim"]>thr_res)][["cancer_type","cancer_type_code",column_x,std_x,"pvalue_diff_distributions","intercept_met","intercept_primary"]])
    
    # now calculate the regression, bootstrap 100
    vals,stats=[],[]
    xseq=range(int(pm_s[column_x].min())-1,int(pm_s[column_x].max())+2,1)
    for N in range(100):
        s=pm_s.sample(frac=0.75,random_state=N)
        slope,intercept,rvalue,pvalue_r,_=st.linregress(s[column_x], y=s["median_sbs1_y"])
        rsp,pvalue_sp=st.spearmanr(s[column_x], s["median_sbs1_y"])
        stats.append([slope,intercept,rvalue,pvalue_r,rsp,pvalue_sp])
        vals.append(xseq*slope + intercept)
    
    slope,intercept,rvalue,pvalue,rsp,pvalue_s=sorted(stats)[50]
    y_pred=xseq*slope + intercept
    
    fig,ax = plt.subplots(figsize=(6,5))

    sns.regplot(data=pm_s,x=column_x,y="median_sbs1_y",line_kws={"lw":5},color="black",scatter=False,robust=True)
    for _i,r in pm_s.iterrows():
        yerr=[[r["median_sbs1_y"]-r["sbs1_5"]],[r["sbs1_95"]-r["median_sbs1_y"]]]
        if xerr and type_x=="fold_change":
            xerr=r[std_x]
        
        sig=r["pvalue_diff_distributions"] <0.01 and r["intercept_met"] > r["intercept_primary"]
        
        ax.errorbar(x=r[column_x],y=r["median_sbs1_y"],yerr=yerr,fmt="o",markeredgewidth=sig*2,color="#888888",markeredgecolor="black",
                   markersize=10,xerr=xerr,lw=1)
        ax.annotate(xy=(r[column_x],r["median_sbs1_y"]+2),s=r["cancer_type_code"],fontsize=14)
    
    #ax.set_ylim(-5,55)
    if pearson:
        ax.annotate(xy=(xlim[1]/3,50),s=f"Pearson correlation coef.={rvalue:1.2f}",fontsize=14)
        ax.annotate(xy=(xlim[1]/3,45),s=f"P-value={pvalue:1.4f}",fontsize=14)
    else:
        ax.annotate(xy=(xlim[1]/3,50),s=f"Spearman ρ.={rsp:1.2f}",fontsize=14)
        ax.annotate(xy=(xlim[1]/3,45),s=f"P-value={pvalue_s:1.4f}",fontsize=14)
        
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_ylabel("SBS1 muts. per year (primary)",fontsize=14)
    ax.set_xlabel(xlabel,fontsize=14)
    if xerr and type_x=="raw":
        ax.set_xlim(-250,400)
    if xerr and type_x=="fold_change":
        ax.set_xlim(xlim)
    
        ax.set_ylim(-5,55)
    
    
    ax.set_title(title,fontsize=16)
    
    plt.savefig(f'../results/figures/turnover_met_prim_{name}.pdf', dpi=800,bbox_inches="tight")
    return slope,intercept,rvalue,pvalue,pm_s

def cell_division_rate_increase_external(df_stats,df_data,column_y="SBS1_genome_wide",column_x="mean_diff_reg_lines",std_x="std_diff_reg_lines",xlabel="Estimated difference SBS1 burden (met - prim.)",name="",ylabel="SBS1 muts. per year \n (Alexandrov et al. 2015)",xlim=(-300,400),xerr=False,pearson=True):

    v=df_data.groupby(["cancer_type_code"],as_index=False).agg(median_sbs1_y=(column_y,np.nanmean))
    pm_s=v.merge(df_stats[(df_stats["rvalue_res_met"]>0.1)&(df_stats["rvalue_res_prim"]>0.1)][["cancer_type","cancer_type_code",column_x,std_x,"pvalue_diff_distributions","intercept_met","intercept_primary"]])
    
    # now calculate the regression, bootstrap 100
    vals,stats=[],[]
    xseq=range(int(pm_s[column_x].min())-1,int(pm_s[column_x].max())+2,1)
    for N in range(100):
        s=pm_s.sample(frac=0.75,random_state=N)
        slope,intercept,rvalue,pvalue_r,_=st.linregress(s[column_x], y=s["median_sbs1_y"])
        rsp,pvalue_sp=st.spearmanr(s[column_x], s["median_sbs1_y"])
        stats.append([slope,intercept,rvalue,pvalue_r,rsp,pvalue_sp])
        
        vals.append(xseq*slope + intercept)
    
    slope,intercept,rvalue,pvalue,rsp,pvalue_s=sorted(stats)[50]
    y_pred=xseq*slope + intercept

    fig,ax = plt.subplots(figsize=(6,5))
    sns.regplot(data=pm_s,x=column_x,y="median_sbs1_y",line_kws={"lw":3},color="black",scatter=False)
    for _i,r in pm_s.iterrows():
        x_err=[0]
        if xerr:
            x_err=[r[std_x]]
        sig=r["pvalue_diff_distributions"] <0.01 and r["intercept_met"] > r["intercept_primary"]
        ax.errorbar(x=r[column_x],y=r["median_sbs1_y"],fmt="o",markeredgewidth=sig*2,color="#888888",markeredgecolor="black",
                   markersize=10,xerr=x_err,lw=1)
        ax.annotate(xy=(r[column_x],r["median_sbs1_y"]+2),s=r["cancer_type_code"],fontsize=14)
    
    ax.set_ylim(-5,80)
    if pearson:
        ax.annotate(xy=(xlim[0],50),s=f"Pearson correlation coef.={rvalue:1.2f}",fontsize=14)
        ax.annotate(xy=(xlim[0],45),s=f"P-value={pvalue:1.4f}",fontsize=14)
    else:
        ax.annotate(xy=(xlim[0],50),s=f"Spearman ρ.={rsp:1.2f}",fontsize=14)
        ax.annotate(xy=(xlim[0],45),s=f"P-value={pvalue_s:1.4f}",fontsize=14)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_ylabel(ylabel,fontsize=14)
    ax.set_xlabel(xlabel,fontsize=14)
    ax.set_xlim(xlim)
    plt.savefig(f'../results/figures/turnover_met_prim_{name}.pdf', dpi=800,bbox_inches="tight")
    return slope,intercept,rvalue,pvalue,pm_s