# 'Pan-cancer whole genome comparison of primary and metastatic solid tumors' 
This repository contains all the data and analysis related to the preprint: https://www.biorxiv.org/content/10.1101/2022.06.17.496528v1


## Project Structure

This repository is structured as follows:

```shell
Analysis type
└── data/            # directory where all external and intermediata data is stored
     
└── code/            # directory where all the code is stored
   
└── results/
    ├── data/        # supplementary tables and results
    ├── figures/     # raw figures
    
├── README.md        # this file
└── environment_analysisA.yaml # Anaconda environment YAML file for specific analysis "karyotype", "driver_enrichment_and_actionability" & "hartwig_pipeline_validation"
└── environment_analysisB.yaml # Anaconda environment YAML file for specific analysis "B"
```

## Data access

### access PCAWG data
Somatic variant calls, gene driver lists, copy number profiles and other core data of the PCAWG cohort generated by the Hartwig analytical pipeline are available for download at https://dcc.icgc.org/releases/PCAWG/Hartwig. Researchers will need to apply to the ICGC data access compliance office (https://daco.icgc-argo.org) for the ICGC portion of the dataset. Similarly, users with authorized access can download the TCGA portion of the PCAWG dataset at https://icgc.bionimbus.org/files/5310a3ac-0344-458a-88ce-d55445540120. Additional information on accessing the data, including raw read files, can be found at https://docs.icgc.org/pcawg/data/.

### access Hartwig data
Metastatic WGS data and metadata from the Hartwig Medical Foundation are freely available for academic use through standardized procedures. Request forms can be found at https://www.hartwigmedicalfoundation.nl/en/data/data-acces-request/
