High level information on the bioinformatics pipeline that was used for the analysis as reported in [Tom Luijts, Kerryn Elliott, Joachim Tetteh Siaw, Joris Van de Velde, Elien Beyls, Arne Claeys, Tim Lammens, Erik Larsson, Wouter Willaert, Anne Vral and Jimmy Van den Eynden. A clinically annotated post-mortem approach to study multi-organ somatic mutational clonality in normal tissues. Scientific Reports, 2022.](https://doi.org/10.1038/s41598-022-14240-8)

# Environment
  
Analysis was performed in a Conda environment. See **mutClon.yml** for details. **scripts/Rpacks** describes R packages that were installed independently.

# Data 

## Downloaded data

The following data were downloaded from external sources:
  - Suppl. data from Healthy Skin (Martincorena et al. 2015)
  - Suppl. data from SCC skin cancer (Chang et al. 2021)

## Raw data

Raw sequencing data are deposited under controlled access in the European Genome–Phenome Archive (EGA; https://ega-archive.org; accession nr. EGAD00001008956).

## Downstream data

Processed data are available in data/ folder
  
# Data processing

## Select genepanel for targeted sequencing

```r
source("scripts/select_TS_genepanel.Rmd")
```

## Aligment

```{bash, eval=F}
scripts/other/align_fastq.sh                    
```

## Mutation calling

ShearwaterML: within patients, use other locations as reference

```{r, eval=F}
source("scripts/run_shearwaterML.R")
```

## Generation of maf file

Processing + Manual curation afterwards: <10bps from each other = manual curation in IGV
Variants annotated using annovar, DNVs added from manual curation

```{r, eval=F}
source("scripts/manuscript_create_maf.R")
```

Filter & add sample alias
```{r, eval=F}
source("scripts/manuscript_filter_maf.R")
```

## Process SimSen results

```{r, eval=F}
source("scripts/process_simSen.R")
```

## Process SCC data

```{r, eval=F}
source("scripts/SCC_process.R")
```

## Process IM2015 data

```{r, eval=F}
source("scripts/IM2015_process.R")
```

# Manuscript

The main analysis, as reported in the manuscript

## Donor description: Fig. 1

```{r, eval=F}
source("scripts/manuscript_donor_description.R")
```

## TS basic description: Fig. 2, Fig. S1c, Fig. S2

```{r, eval=F}
source("scripts/manuscript_TS_description.R")
```

## Selection analysis: Fig. 3, Fig. S3

### Calculate dN/dS & PP2

```{r, eval=F}
# Get simulation data
source("scripts/get_sim_PM.R")
### Calculate dN/dS & PP2
source("scripts/manuscript_selection.R")
```

### Plot
```{r, eval=F}
source("scripts/manuscript_plot_selection.R")
```

## Clonal analysis: Fig. 4 & S4

### Fig. 4
```{r, eval=F}
source("scripts/manuscript_clonal_analysis.R")
```

### Fig. S4
```{r, eval=F}
source("scripts/manuscript_clonal_analysis_ind.R")
```

## Manuscript UV analysis: Fig. 5

```{r, eval=F}
source("scripts/manuscript_UV.R")
```

## Summarizing table: table S1

```{r, eval=F}
source("scripts/manuscript_table_S1.R")
```
