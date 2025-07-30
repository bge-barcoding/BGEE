# Barcode Gene Extractor & Evaluator (BGEE) Snakemake workflow #
Snakemake workflow for recovering high-quality barcode sequences from genome skim data, built around MitoGeneExtractor and adapted for genome skims of museum specimens.

# Contents # 
 - [Requirements](#Requirements)
 - [Workflow](#Workflow)
 - [Running](#running)
 - [Cluster configuration](#Cluster-configuration-using-profiles/slurm/config.yaml)
 - [Cluster submission](#Cluster-submission)
 - [Results structure](#Results-structure)
 - [Contributing](#Contributing)

# Requirements #
- [MitoGeneExtractor](https://github.com/cmayer/MitoGeneExtractor) version 1.9.6 installed. Clone repository and follow installation instructons. 
- Paired-end reads in .fastq.gz or .fastq format.
- samples_file.csv (generated manually, or as outlined [below](https://github.com/SchistoDan/BGEE?tab=readme-ov-file#2-generate-samplescsv)) if working from BOLD sample metdata.
- sequence_references_file.csv (generated manually, or using [Gene Fetch](https://github.com/bge-barcoding/gene_fetch?tab=readme-ov-file)).
- Activated conda env (see bgee_env.yaml).

# Workflow #
1. Preprocessing mode (both pre-processing modes are run in parallel):
  - 'concat':
    - Initial raw read quality control - Adapter, trimming, quality trimming, poly-g trimming, and deduplication of paired-end reads using [fastp](https://github.com/OpenGene/fastp) (fastp_pe_concat).
    - Concantenation of trimmed PE reads (fastq_concat) and associated log files (aggregate_concat_logs).
    - Compress trimmed reads (gziped_trimmed).
    - Secondary read quality control - Additional quality trimming of concatenated reads using [trim galore](https://github.com/FelixKrueger/TrimGalore) (quality_trim) and concatenation of associated log files (Aggregate_trim_galore_logs).
  - 'merge:
    - Raw read quality control and merging - Adapter, trimming, quality trimming, poly-g trimming, deduplication, and merging of paired-end reads using [fastp](https://github.com/OpenGene/fastp) (fastp_pe_merge).
    - 'Clean' headers of input files, as required by MitoGeneExtractor (clean_headers_merge), and concentation of associated log files (Aggregate_clean_headers_logs).
![image](https://github.com/user-attachments/assets/139b8c7c-b0dc-465c-8c95-e3a58ea1ab96)
2. Refernce-guided barcode recovery using [MitoGeneExtractor](https://github.com/cmayer/MitoGeneExtractor) (MitoGeneExtractor_concat & MitoGeneExtractor_merge).
3. Rename raw consensus sequence headers and concatenation into multi-FASTA (rename_and_combine_cons_concat & rename_and_combine_cons_merge) (uses supplementary [rename_headers.py](https://github.com/SchistoDan/BGEE/blob/main/workflow/scripts/rename_headers.py).
4. Remove any remaining exonerate intermediate (Concatenated_exonerate_input_*) files leftover after MGE (remove_exonerate_intermediates).
5. Create list of MGE (FASTA) alignment files for downstream processing (create_alignment_log).
6. Filter MGE alignment files to remove low-quality, contaminant, or outlier sequences before repeat consensus sequence generation.
   - Remove aligned reads with similarity to human COI (a common contaminant of museum specimens) (uses supplementary [human_cox1_filter.py](https://github.com/SchistoDan/BGEE/blob/main/workflow/scripts/human_cox1_filter.py).
   - Remove aligned reads with AT content over and/or below a specified threshold (high AT content can be indicative of contamination (e.g. fungal)) (uses supplementary [at_content_filter.py](https://github.com/SchistoDan/BGEE/blob/main/workflow/scripts/at_content_filter.py).
   - Remove reads that are statistical outliers compared to the original consensus (uses supplementary [statistical_outliers.py](https://github.com/SchistoDan/BGEE/blob/main/workflow/scripts/statistical_outlier_filter.py).
   - [optional] Remove reads that are statistical outliers compared to the supplied reference sequence(s) (uses supplementary [reference_filter.py](https://github.com/SchistoDan/BGEE/blob/main/workflow/scripts/reference_filter.py).
   - Generation of 'cleaned' consensus sequence (uses supplementary [consensus_generator.py](https://github.com/SchistoDan/BGEE/blob/main/workflow/scripts/consensus_generator.py).
   - Aggregate metrics from each stage of filtering (uses supplementary [aggregate_filter_metrics.py](https://github.com/SchistoDan/BGEE/blob/main/workflow/scripts/aggregate_filter_metrics.py).
   - Remove intermediate files and unecessary logs generated during the consensus cleaning process (remove_fasta_cleaner_files).
8. Evaluate barcode consensus sequence quality based on various metrics (length, ambiguous base content, etc.), and select the 'best' sequences according to specific ranking criteria: Rank1 = No ambiguous bases & longest stretch ≥ 650, Rank3 = No ambiguous bases & longest stretch ≥ 500, Rank3 = No ambiguous bases & 300 ≤ longest stretch ≤ 499, Rank4 = No ambiguous bases &1 ≤ longest stretch ≤ 299, Rank5 = Has ambiguous bases. (uses supplementary [fasta_compare](https://github.com/SchistoDan/BGEE/blob/main/workflow/scripts/fasta_compare.py).
9. Compile statistics from read QC, MGE, and consensus cleaning metrics into a CSV report for both 'concat' and 'merge' modes (uses supplementary [mge_stats.py](https://github.com/SchistoDan/BGEE/blob/main/workflow/scripts/mge_stats.py) and combine_stats_files).
10. Remove temporary files, sample-specific logs once aggregated, etc. (cleanup_files).



# Workflow execution: #
## Clone BGEE github repository and set up conda environment ##
- [Install miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions).
```bash
git clone https://github.com/bge-barcoding/BGEE.git [path/to/desired/install/location/]
cd BGEE/installation/dir/
conda env create -f /workflow/scripts/bgee_env.yaml
git status
```

## Generate samples.csv ###
- Can be created manually, or via [sample-processing](https://github.com/bge-barcoding/sample-processing) workflow.
- Must contain ID, forward (read paths), reverse (read paths), and taxid columns (see below for example). Column 1 can be named 'ID', 'process_id', 'Process ID', 'process id', 'Process id', 'PROCESS ID', 'sample', 'SAMPLE', or 'Sample'.
- Due to regex matching and statistics aggregation, the sample's 'ID' will be considered the string before the first underscore across output files. It is therefore necessary that sample names do not use '_' characters. E.g. BSNHM002-24 instead of BSNHM002_24, or P3-1-A10-2-G1 instead of P3_1_A10_2_G1.
- Taxid's can be found manually by searching the expected species/genus/family of each sample in the [NCBI taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy), or retrieved from the sample_metadata.csv file output by the [sample-processing](https://github.com/bge-barcoding/sample-processing) workflow.
  
**samples.csv example**
| ID | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177658 |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177627 |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 3084599 |

## Gathering sample-specific pseudo-references ##
- This can be created manually, or using [Gene-fetch](https://github.com/bge-barcoding/gene_fetch). gene-fetch will retrieve the necessary protein pseudo-references for each sample using the samples taxonomic identifier (taxid) or taxonomic lineages for each sample.
- Must contain 'process_id', 'reference'name' and 'protein_reference_path' at a minimum.
- Similarly to `samples.csv`, `proces_id` and `reference_name` must not use '_' characters. E.g. BSNHM002-24 instead of BSNHM002_24, or ref-1 instead of ref_1.

**samples.csv example**
| process_id | reference_name | protein_reference_path | 
| --- | --- | --- |
| BSNHM002-24  | ref-1 | path/to/BSNHM002-24.fasta |
| BSNHM038-24 | ref-2 | path/to/BSNHM038-24.fasta |
| BSNHM046-24 | ref-3 | path/toBSNHM046-24.fasta |

## Customising snakemake configuration file ##
- Update config.yaml with neccessary paths and variables.
  - See [MitoGeneExtractor README.md](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/README.md) for explanation of Exonernate ('r' and 's') run paramters.
  - See [config.yaml](https://github.com/SchistoDan/BGEE/blob/main/config/config.yaml) for information on filtering variables and thresholds (default is likely suitable in most cases).
  - See [Gene Fetch repository](https://github.com/bge-barcoding/gene_fetch) for guidance on creating [sequence_reference_file.csv](https://github.com/bge-barcoding/gene_fetch?tab=readme-ov-file#normal-mode).
- Rule-specific resources in the [config.yaml](https://github.com/SchistoDan/BGEE/blob/main/config/config.yaml) - Each rule can specify the necessary number of threads and memory resources (in Mb) for every job (e.g. specifying 4 threads and 4G memory for fastp_pe_merge would allocate those resources for every 'fastp_pe_merge' job).
```
# config/config.yaml
# MGE run name identifier
run_name: "RUN_1"

# Path to MGE installation (MitoGeneExtractor-vX.X.X file)
mge_path: "/path/to/MitoGeneExtractor-1.9.6beta3/MitoGeneExtractor-v1.9.6beta3"

# Path to samples.csv
samples_file: "/path/to/input/samples/file/samples.csv"

# Path to references.csv
sequence_reference_file: "/path/to/input/reference/sequence/csv/sequence_references.csv"

# Path to output directory. Will make final dir if it does not exist already
output_dir: "/path/to/output/directory"

## MGE running parameters
# Exonerate relative score threshold parameter
r:
  - 1
  - 1.3
  - 1.5

# Exonerate minimum score threshold parameter
s:
  - 50
  - 100
# Number of bp to extend beyond the Exonerate alignment
n: 0
# Genetic code (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
C: 5
# Consensus threshold (e.g. 0.5 = 50%)
t: 0.5

## 'Cleaning' of recovered barcode sequences
fasta_cleaner:
  consensus_threshold: 0.5  # 
  human_threshold: 0.95     # Threshold for human COI similarity (e.g. >95% similarity = read removed)
  at_difference: 0.1        # Maximum allowed difference in AT content between each read and the consensus before filtering (e.g. if consensus AT content = 0.65, and a read has <0.55 or >0.75 AT content, discard it (if in "absolute" mode)
  at_mode: "absolute"       # absolute (remove reads if AT content differs from consensus by more than threshold in either direction), higher (remove reads with AT content above consensus + threshold), lower (remove reads with AT content below consensus - threshold)
  outlier_percentile: 90.0  # Threshold for 'deviant' read removal (e.g. the bottom 90% of sequences (most similar to consensus) are kept)
  disable_human: false      # If false, 'human_cox1_filter.py' will run
  disable_at: false         # If false, 'at_content_filter.py' will run
  disable_outliers: false   # If false, 'statistical_outliers.py' will run
  reference_dir: null       # To skip reference filtering enter "null"/"None". To enable reference filtering, provide path to references


## Comparison of consensus sequences produced directly by MGE and those with additional cleaning (see scripts/fasta_compare.py)
# Options: cox1/COI/CO1, rbcl/RBCL, matk/MATK
run_fasta_compare: true  # Set false to skip fasta_compare step

fasta_compare:
  target: "cox1"
  verbose: false

## Resource allocation for each rule. Memory shown in Mb (mem_mb), e.g. 2048 = 2G memory.
# Rules have dynamic memory scaling upon retry (mem_mb * retry #).
rules:
  fastp_pe_merge:
    mem_mb: 8192
    threads: 4
   
  clean_headers_merge:
    mem_mb: 4096
    threads: 2
    
  aggregate_clean_headers_logs_merge:
    mem_mb: 2048
    threads: 2
    
  fastp_pe_concat:
    mem_mb: 8192
    threads: 4
    
  fastq_concat:
    mem_mb: 4096
    threads: 4
    
  gzip_trimmed:
    mem_mb: 4096
    threads: 4
    
  aggregate_concat_logs:
    mem_mb: 2048
    threads: 2
    
  quality_trim:
    mem_mb: 8192
    threads: 4
    
  aggregate_trim_galore_logs:
    mem_mb: 2048
    threads: 2
    
  MitoGeneExtractor_merge:
    mem_mb: 20480
    threads: 6
    partition: medium
    
  MitoGeneExtractor_concat:
    mem_mb: 20480
    threads: 6
    partition: medium
    
  rename_and_combine_cons_merge:
    mem_mb: 2048
    threads: 4
    
  rename_and_combine_cons_concat:
    mem_mb: 2048
    threads: 4
    
  create_alignment_log_merge:
    mem_mb: 4096
    threads: 2
    
  create_alignment_log_concat:
    mem_mb: 4096
    threads: 2
    
  human_cox1_filter_merge:
    mem_mb: 10240
    threads: 8
    
  at_content_filter_merge:
    mem_mb: 15360
    threads: 8
    
  statistical_outlier_filter_merge:
    mem_mb: 8192
    threads: 8
    
  reference_filter_merge:
    mem_mb: 4096
    threads: 4
    
  consensus_generation_merge:
    mem_mb: 4096
    threads: 4
    
  aggregate_filter_metrics_merge:
    mem_mb: 2048
    threads: 1

  human_cox1_filter_concat:
    mem_mb: 10240
    threads: 8
  
  at_content_filter_concat:
    mem_mb: 15360
    threads: 8
    
  statistical_outlier_filter_concat:
    mem_mb: 4096
    threads: 1
    
  reference_filter_concat:
    mem_mb: 4096
    threads: 4
    
  consensus_generation_concat:
    mem_mb: 4096
    threads: 4
    
  aggregate_filter_metrics_concat:
    mem_mb: 2048
    threads: 1
    
  extract_stats_to_csv_merge:
    mem_mb: 4096
    threads: 4
    
  extract_stats_to_csv_concat:
    mem_mb: 4096
    threads: 4
    
  combine_stats_files:
    mem_mb: 1024
    threads: 1
    
  fasta_compare:
    mem_mb: 4096
    threads: 4
    
  remove_fasta_cleaner_files:
    mem_mb: 4096
    threads: 4
    
  remove_exonerate_intermediates:
    mem_mb: 4096
    threads: 2
    
  cleanup_files:
    mem_mb: 1024
    threads: 1
```

# Cluster configuration using profiles/slurm/config.yaml #
- See `profiles/slurm/config.yaml` below for SLURM cluster submission parameters. You will need to change `slurm_parition` to a suitable parition for your cluster.
```
# profiles/slurm/config.yaml
executor: slurm

default-resources:
  - slurm_partition=medium
  - slurm_extra="--ntasks=1 --signal=USR2@90 --no-requeue"
    # --ntasks=1: Run one instance of the workflow.
    # --signal=USR2@90: Sets up a signal to be sent to the job when it approaches its time limit.
    # --no-requeue: Prevents the job from being automatically requeued if it fails (clashes with snakemake retry memory scaling).
  
jobs: 10                          # Maximum number of concurrent jobs run on cluster. # Remember not to exceed any thread limits. E.g. 20 concurrent jobs × 16 threads per job = up to 320 threads!
latency-wait: 60                  # Wait 60 seconds before checking job status
local-cores: 0                    # Prevents fallback to local execution upon job failure
max-jobs-per-second: 1            # Controls rate Snakemake submits jobs to SLURM. Prevents overwhelming scheduler  
```

# Cluster submission #
- [snakemake_run.sh](https://github.com/SchistoDan/BGEE/blob/main/snakemake_run.sh) handles submission of the snakemake command to the cluster.
- The working directory will first be unlocked (using `--unlock`) and then the snakemake workflow will be run.
- Submit snakemake_run.sh to the cluster with `./snakemake_run.sh` (if in `BGEE/ directory') - This will submit `snakemake_run.sh` to the head/login node of your cluster, and will orchestrator SLURM job submissions for all snakemake workflow steps.
 



# Results structure #
```
output_dir/
├── merge_mode/
│   ├── consensus/
│   │   ├── {sample}_r_{r}_s_{s}_con_{reference_name}.fas      # Individual consensus files for each sample/parameter combination. 
│   │   └── {run_name}_merge.fasta                             # Combined consensus sequences multi-FASTA for merge mode. sequences will have the fasta header '{sample}_r_{r}_s_{s}_con_{reference_name}_merge' after renaming.
│   ├── alignment/                                             # Alignment files containing trimmed reads aligned by MGE/exonerate (i.e reads that went into generating the consensus sequence).
│   ├── trimmed_data/                                          # Intermediate. Removed at end of run.
│   ├── fasta_cleaner/                                         # Consensus 'cleaning' intermediates and outputs. Each 01-05 directory contains metrics related to each cleaning step.
│   │   ├── 01_human_filtered/
│   │   ├── 02_at_filtered/
│   │   ├── 03_outlier_filtered/
│   │   ├── [optional] 04_reference_filtered/
│   │   ├── 05_cleaned_consensus/
│   │   ├── combined_statistics.csv                            # Combined statistics from all cleaning steps.
│   │   └── all_consensus_sequences.fasta                      # Cleaned consensus sequences multi-FASTA.
│   ├── logs/
│   │   ├── clean_headers.log                                  # Aggregated logs for fasta header cleaning for all samples.
│   │   ├── cleaning_complete.txt
│   │   ├── file_cleanup_complete.txt
│   │   ├── rename_fasta.log
│   │   ├── fasta_cleaner/                                     # Contains logs for each consensus sequence cleaning step.
│   │   ├── fastp/                                             # Contains fastp trimming logs for each sample.
│   │   ├── mge/
│   │   │   └── {sample}_r_{r}_s_{s}_{reference_name}/         # MGE vulgar files for each sample
│   │   │   └── alignment_files.log
│   │   │   └── mge_stats.log
│   ├── out/                                                   # MGE output files for each sample.
│   ├── err/                                                   # MGE error/processing logs for each sample.
│   ├── {run_name}_merge-stats.csv                             # Summary statistics file for merge mode. 
│   └── cleanup_complete.txt
│
├── concat_mode/
│   ├── consensus/
│   │   ├── {sample}_r_{r}_s_{s}_con_{reference_name}.fas       # Individual consensus files for each sample/parameter combination.
│   │   └── {run_name}_concat.fasta                            # Combined consensus sequences multi-FASTA for merge mode. sequences will have the fasta header '{sample}_r_{r}_s_{s}_con_{reference_name}' after renaming.
│   ├── alignment/                                             # Alignment files containing trimmed reads aligned by MGE/exonerate (i.e reads that went into generating the consensus sequence).
│   ├── trimmed_data/
│   │   ├── {sample}/                                          # Trimmed fwd and rev reads (fq.gz), fastp JSON and HTML reports, and Trim Galore QC report.
│   ├── fasta_cleaner/                                         # Consensus 'cleaning' intermediates and outputs. Each 01-05 directory contains metrics related to each cleaning step.
│   │   ├── 01_human_filtered/
│   │   ├── 02_at_filtered/
│   │   ├── 03_outlier_filtered/
│   │   ├── [optional] 04_reference_filtered/
│   │   ├── 05_cleaned_consensus/
│   │   ├── combined_statistics.csv                            # Combined statistics from all cleaning steps.
│   │   └── all_consensus_sequences.fasta                      # Cleaned consensus sequences multi-FASTA.
│   ├── logs/
│   │   ├── alignment_files.log
│   │   ├── concat/
│   │   │   └── {sample}.log
│   │   ├── concat_reads.log
│   │   ├── cleaning_complete.txt
│   │   ├── fasta_cleaner.log
│   │   ├── mge/
│   │   │   └── {sample}_vulgar_r_{r}_s_{s}_{reference_name}.txt
│   │   ├── rename_complete.txt
│   │   ├── trim_galore/
│   │   │   └── {sample}.log
│   ├── logs/
│   │   ├── concat_reads.log                                   # Aggregated logs for trimmed PE read concatentation.
│   │   ├── cleaning_complete.txt
│   │   ├── file_cleanup_complete.txt
│   │   ├── rename_fasta.log
│   │   ├── fasta_cleaner/                                     # Contains logs for each consensus sequence cleaning step.
│   │   ├── fastp/                                             # Contains fastp trimming logs for each sample.
│   │   ├── gzip/                                              # Logs for each trimmed read compression.
│   │   ├── mge/
│   │   │   └── {sample}_r_{r}_s_{s}_{reference_name}/         # MGE vulgar files for each sample
│   │   │   └── alignment_files.log
│   │   │   └── mge_stats.log
│   │   └── trim_galore.log
│   ├── out/                                                   # MGE output files for each sample.
│   ├── err/                                                   # MGE error/processing logs for each sample.
│   ├── {run_name}_concat-stats.csv                             # Summary statistics file for concat mode. 
│   └── cleanup_complete.txt
│
├── fasta_compare/
│   ├── {run_name}_fasta_compare.csv                          # Comparison results between modes
│   ├── {run_name}_full_sequences.fasta                       # Full consensus sequences from both modes
│   ├── {run_name}_barcode_sequences.fasta                    # Barcode sequences from both modes
│   └── fasta_compare_{run_name}.log
└── {run_name}_combined_stats.csv                             # Combined statistics from both modes
```


# Contributing #
- This snakemake pipeline was produced by Dan Parsons @ NHMUK for BioDiversity Genomics Europe.
- - Since this snakemake pipeline uses [MitogeneExtractor](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14075) at its core, please cite:
  Brasseur, M.V., Astrin, J.J., Geiger, M.F., Mayer, C., 2023. MitoGeneExtractor: Efficient extraction of mitochondrial genes from next-generation sequencing libraries. Methods in Ecology and Evolution.
- Please feel free to submit issues, fork the repository, and create pull requests for any enhancements or bugs.



  ## To do ##
- Split Snakefile into .smk files
- Integrate Gene fetch into workflow.
- Add pre-MGE rule to subset very large input files (Based on file size or sequence number related)
- Get workflow to generate RO-crates.
 
