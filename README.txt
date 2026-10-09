Data from: Chromosome evolution model reveals hidden variation in fern diversification rates

Dataset DOI: 10.5061/dryad.612jm64hr
Associated Manuscript DOI: 10.1093/sysbio/syag054

Description of the data and file structure

Description of the data and file structure

Repository contains the datasets and scripts needed to reproduce the findings of the article published by Systematic Biology titled "Chromosome Model Reveals Hidden Variation in Karyotype-Driven Speciation of Ferns."

Repository contains three compressed directories containing: data ('data.zip'), output ('output.zip'), and scripts ('scripts.zip'). 

data*&#xNAN;—*

Directory contains the empirical and simulated datasets used for analysis in the associated study. Chromosome count data are in .tsv or .csv file format and phylogenetic trees are in .tree file format. Data are available in both raw and analysis-ready formats.



/rawdata

Raw data were subset from the Chromosome Counts Database. These are chromosome count data for each species relevant to our study before pruning to the phylogeny. Files include the following columns: raw chromosome count distribution data ('distribution'), max count ('max'), median count ('median'), min count ('min'), and the resolved species name ('resolved_name).

/analyzedSubsets

The analyzed subsets represent the three datasets analyzed in the associated article. Trees are in Newick format. For each clade analyzed, two chromosome count datasets are available. One is structured for the ChromoSSE model ('counts') and the other is structured for the ChromoHiSSE model ('counts_expanded'). These datasets, or others, can be regenerated with the included scripts.

/simulatedDatasets

Simulated datasets are numbered to associate the simulated tree and the chromosome count data in both simple and expanded formats (for each of the models tested).

output —

Directory contains the raw trace files from all MCMC analyses in RevBayes for the empirical analysis (empirical_logs) and simulations (simulation_logs). These are used to compute information criteria and summarize posterior distributions.



/empirical_logs

Directory contains the log files for all analyses of the empirical clades assessed in the study. These files contain the unpruned model output necessary to summarize MCMC convergence, posterior distributions of model parameters, and stochastic character maps produced by the RevBayes analysis. Checkpoint files for each MCMC chain are available which will allow analyses to be restarted from their last checkpoint.

/simulation_logs

Directory contains all analyses of simulated datasets. Subdirectories are named to correspond to the model used to simulate the data followed by the model used to analyze the data. For example, datasets generated under ChromoSSE and analyzed under the ChromoHiSSE model, and vice versa, would be named 'csse_chisse' and 'chisse_csse', respectively.



Log files are structured such that the a set of headers at the top are followed by the values for each variable. 

ChromoSSE log file heater variables and their interpretation: 

Iteration (mcmc generation),

posterior (posterior probability),

likelihood (log likelihood),

prior (prior probability),

clado_demipoly (cladogenetic rate of demipolyploidy (events per lineage per million years (E/L/My))), 

clado_fission (cladogenetic rate of descending dysploidy (E/L/My)), 

clado_fusion (cladogenetic rate of ascending dysploidy (E/L/My)), 

clado_no_change (cladogenetic rate of speciation without chromosome number change (E/L/My)), 

clado_polyploid (cladogenetic rate of polyploidy (E/L/My)), 

clado_rates[1-5] (recoded rates for each cladogenetic parameter: 1 -> clado_no_change, 2 -> clado_fission, 3 -> clado_fusion, 4 -> clado_polyploidy, 5 -> clado_demipoly), 

delta (anagenetic rate of descending dysploidy (events per lineage per million years (E/L/My))), 

eta (anagenetic rate of demipolyploidy (events per lineage per million years (E/L/My))), 

gamma (anagenetic rate of ascending dysploidy (events per lineage per million years (E/L/My))), 

relative_clado[1-5] (cladogenetic rate of demipolyploidy (events per lineage per million years (E/L/My))), 

rho (anagenetic rate of polyploidy (events per lineage per million years (E/L/My))), 

total_speciation (total cladogenetic rate (E/L/My)), 

turnover (extinction (E/L/My)) 

ChromoHiSSE contains the same headers along with  a second version for each rate parameter corresponding to each hidden state. These log files are used for summarizing the model estimates and posterior distributions in R. 

Stochastic character mapping log files (under 'stochasticMapping') share the same basic structure as the model logs. A header specifies the MCMC generation and each state in the model. These log files are used for plotting the stochastic character maps in R. 

scripts —

Directory contains the RevBayes (revbayes git commit 9919664) and R (R version 4.2.3) scripts necessary to regenerate the inputs and outputs from the manuscript. Scripts for every step reported in the associated article are available in the subdirectories:

/convergence_diagnostics - calculating convergence diagnostics in R,

/datapreparation - preparing/reformatting/cleaning data in R for analysis with RevBayes models,

/models_Revbayes - fitting models (ChromoSSE and ChromoHiSSE) in RevBayes,

/modelSelection - performing model selection in R,

/parameterTables - building parameter estimate and summary tables in R,

/plotting - plotting results in R,

/simulations - summarizing simulations in R,

/summaryStats - calculating custom summary statistics in R.



/models_Revbayes

Subdirectories contain the Rev scripts for running our ChromoSSE or ChromoHiSSE models in RevBayes. To run the analyses, a control script (e.g., 'chromohisse_modelspec_tb.rev') must be modified to provide the paths to input and output files/directories, set priors, specify MCMC parameters, and execute the model script (e.g., 'chromohisse_tp_mrm_tb.rev'). Both ChromoSSE and ChromoHiSSE models run on a specific branch of RevBayes (mrm_dev branch (github commit: 9919664). Due to computational burdens, these analyses are highly recommended to be run using the TensorPhylo plugin for RevBayes.

Files and variables

.

├── data

│   ├── analyzedSubsets

│   │   ├── aspleniineae

│   │   ├── polypodiineae

│   │   └── pteridaceae

│   ├── rawdata

│   └── simulatedDatasets

│       ├── chisse

│       │   ├── chisse

│       │   └── csse

│       └── csse

│           ├── chisse

│           └── csse

├── output

│   ├── empirical_logs

│   │   ├── modelFitting

│   │   │   ├── chromohisse

│   │   │   │   ├── asple

│   │   │   │   │   └── checkpts

│   │   │   │   ├── polyp

│   │   │   │   │   └── checkpts

│   │   │   │   └── pteri

│   │   │   │       └── checkpts

│   │   │   └── chromosse

│   │   │       ├── asple

│   │   │       │   └── checkpts

│   │   │       ├── polyp

│   │   │       │   └── checkpts

│   │   │       └── pteri

│   │   │           └── checkpts

│   │   └── stochasticMapping

│   │       ├── chromohisse

│   │       │   ├── asple

│   │       │   ├── polyp

│   │       │   └── pteri

│   │       └── chromosse

│   │           ├── asple

│   │           ├── polyp

│   │           └── pteri

│   └── simulation_logs

│       ├── chisse_chisse

│       │   ├── checkpts

│       │   └── peeks

│       ├── chisse_csse

│       │   ├── checkpts

│       │   └── peeks

│       ├── csse_chisse

│       │   ├── checkpts

│       │   └── peeks

│       └── csse_csse

│           ├── checkpts

│           └── peeks

└── scripts

    ├── convergence_diagnostics

    ├── datapreparation

    ├── models_Revbayes

    │   ├── chromoHiSSE

    │   └── chromoSSE

    ├── modelSelection

    ├── parameterTables

    ├── plotting

    │   ├── stochmaps

    │   ├── summaryStats

    │   └── violins

    ├── simulations

    │   ├── analyzeSimulations

    │   │   ├── chisse

    │   │   │   └── templates

    │   │   ├── csse

    │   │   │   └── templates

    │   │   └── makeBatchScripts

    │   ├── simulateData

    │   │   └── src

    │   │       └── sims

    │   └── summarizeSimulations

    └── summaryStats



77 directories
