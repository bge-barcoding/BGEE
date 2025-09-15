# Barcode Gene Extractor & Evaluator (BGEE) Snakemake workflow #
Snakemake workflow for recovering high-quality barcode sequences from genome skim data, built around MitoGeneExtractor and adapted for genome skims of museum specimens.

# Contents # 
 - [Requirements](#Requirements)
 - [Workflow](#Workflow)
 - [Installation and set up](#Installation-and-set-up)
 - [Cluster configuration](#Cluster-configuration-using-profiles/slurm/config.yaml)
 - [Cluster submission](#Cluster-submission)
 - [Results structure](#Results-structure)
 - [Contributing](#Contributing)

# Requirements #
- [MitoGeneExtractor](https://github.com/cmayer/MitoGeneExtractor) version 1.9.6 installed. Clone repository and follow installation instructons. 
- Paired-end reads in .fastq.gz or .fastq format.
- samples_file.csv (generated manually, or as outlined below if working from BOLD sample metadata).
- sequence_references_file.csv (generated manually, or using [Gene Fetch](https://github.com/bge-barcoding/gene_fetch?tab=readme-ov-file) within the workflow).
- Activated conda env (see bgee_env.yaml).

# Workflow #
1. **Preprocessing mode** (both pre-processing modes are run in parallel):
   - 'concat':
     - Initial raw read quality control - Adapter, trimming, quality trimming, poly-g trimming, and deduplication of paired-end reads using [fastp](https://github.com/OpenGene/fastp) (fastp_pe_concat).
     - Concantenation of trimmed PE reads (fastq_concat) and associated log files (aggregate_concat_logs).
     - Compress trimmed reads (gziped_trimmed).
     - Secondary read quality control - Additional quality trimming of concatenated reads using [trim galore](https://github.com/FelixKrueger/TrimGalore) (quality_trim) and concatenation of associated log files (Aggregate_trim_galore_logs).
   - 'merge:
     - Raw read quality control and merging - Adapter, trimming, quality trimming, poly-g trimming, deduplication, and merging of paired-end reads using [fastp](https://github.com/OpenGene/fastp) (fastp_pe_merge).
     - 'Clean' headers of input files, as required by MitoGeneExtractor (clean_headers_merge), and concentation of associated log files (Aggregate_clean_headers_logs).

<p align="center">
  <img src="https://github.com/user-attachments/assets/139b8c7c-b0dc-465c-8c95-e3a58ea1ab96" width="650"/>
</p>

2. **Sample-specific pseudo-reference retrieval** from GenBank using [Gene-Fetch](https://github.com/bge-barcoding/gene_fetch). (gene_fetch).
3. **Protein reference-guided barcode recovery** using [MitoGeneExtractor](https://github.com/cmayer/MitoGeneExtractor) (MitoGeneExtractor_concat & MitoGeneExtractor_merge).
4. **'Raw' consensus sequence header manipulation and concatenation** into multi-FASTA (rename_and_combine_cons_concat & rename_and_combine_cons_merge) (uses supplementary [rename_headers.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/rename_headers.py).
5. **Remove exonerate intermediates** if remaining after MGE (remove_exonerate_intermediates).
6. **Generate list of MGE alignment files** for downstream processing (create_alignment_log).
7. **Filter MGE alignment files to remove low-quality, contaminant, or outlier sequences before repeat consensus sequence generation.**
   - Remove aligned reads with similarity to human COI (a common contaminant of museum specimens) (uses supplementary [01_human_cox1_filter.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/01_human_cox1_filter.py)).
   - Remove aligned reads with AT content over and/or below a specified threshold (high AT content can be indicative of contamination (e.g. fungal)) (uses supplementary [02_at_content_filter.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/02_at_content_filter.py)).
   - Remove reads that are statistical outliers compared to the original consensus (uses supplementary [03_statistical_outliers.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/03_statistical_outlier_filter.py)).
   - [optional] Remove reads with similarity to supplied reference sequence(s) (uses supplementary [04_reference_filter.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/04_reference_filter.py)).
   - Generation of 'cleaned' consensus sequence (uses supplementary [05_consensus_generator.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/05_consensus_generator.py)).
   - Aggregate metrics from each stage of filtering (uses supplementary [06_aggregate_filter_metrics.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/06_aggregate_filter_metrics.py)).
   - Remove intermediate files and unecessary logs generated during the consensus cleaning process (remove_fasta_cleaner_files).
  
<p align="center">
  <img width="285" height="443" alt="image" src="https://github.com/user-attachments/assets/957d43a7-0c00-40ce-bab8-1827d0e37e1b" />
</p>

8. **Validate and triage barcodes** 
   - Structural validation
   - Local BLAST of structurally validated barcodes
   - Taxonomic validation of BLAST results
10. **Compile statistics** from read QC, MGE, and consensus cleaning metrics into a CSV report for both 'concat' and 'merge' modes (uses supplementary [mge_stats.py](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/workflow/scripts/mge_stats.py) and combine_stats_files).
11. **Clean up** temporary files, sample-specific logs once aggregated, etc. (cleanup_files).



# Installation and set up: #
## Clone BGEE github repository and set up conda environment ##
- [Install miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions).
```bash
git clone https://github.com/bge-barcoding/BGEE.git [path/to/desired/install/location/]
cd installation/dir/BGEE
conda env create -f bgee_env.yaml
git status
```

## Generate sample input file ###
- Can be created manually, or via [sample-processing](https://github.com/bge-barcoding/sample-processing) workflow.
- **Must contain `ID`, `forward` (read paths), `reverse` (read paths), and `taxid` _OR_ `hierarchical taxonomy` columns (see below for examples).**
- Due to regex matching and statistics aggregation, the sample ID will be considered as the string before the first underscore. **It is therefore recommended that sample names do not use '_' characters.** E.g. BSNHM002-24 instead of BSNHM002_24, or P3-1-A10-2-G1 instead of P3_1_A10_2_G1.
- Taxid's can be found manually by searching the expected species/genus/family of each sample in the [NCBI taxonomy database](https://www.ncbi.nlm.nih.gov/taxonomy).
  
**samples.csv example (taxid)**
| ID | forward | reverse | taxid |
| --- | --- | --- | --- |
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177658 |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 177627 |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | 3084599 |

**samples.csv example (hierarchical taxonomy)**
| ID | forward | reverse | phylum | class | order | family | genus | species |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | 
| BSNHM002-24  | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | Arthropoda | Insecta | Hemiptera | Cicadidae | Tibicina | Tibicina tomentosa |
| BSNHM038-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | Tracheophyta | Pinopsida | Pinales | Pinaceae | Abies |  |
| BSNHM046-24 | abs/path/to/R1.fq.gz | abs/path/to/R2.fq.gz | Annelida | Polychaeta | Terebellida | Ampharetidae | Samytha | Samytha sexcirrata |

## Gathering sample-specific pseudo-references ##
- This can be created manually, or using [Gene-fetch](https://github.com/bge-barcoding/gene_fetch) integrated into the workflow (highly recommended). If enabled (in the config.yaml by setting `run_gene_fetch` to 'true'), gene-fetch will retrieve the necessary protein pseudo-references for each sample from NCBI GenBank using the sample's taxonomic identifier (taxid)/taxonomic hierarchy for each sample, a sequence target (e.g. COI or rbcL), and NCBI API credentials (email address & API key - see [guidance](https://support.nlm.nih.gov/kbArticle/?pn=KA-05317) on getting a key).
- If hierarchical taxonomic information is provided (see `samples.csv example (hierarchical taxonomy)` above) for each sample instead of a taxid (see `samples.csv example (taxid)` above), `Gene Fetch` will, starting from the lowest given rank (e.g. species), find the closest valid taxid on NCBI for the provided taxonomy before proceeding with pseudo-reference fetching.
- **Must contain 'process_id', 'reference'name' and 'protein_reference_path' at a minimum.**

**sample_references.csv example**
| process_id | reference_name | protein_reference_path | 
| --- | --- | --- |
| BSNHM002-24  | BSNHM002-24 | path/to/BSNHM002-24.fasta |
| BSNHM038-24 | BSNHM038-24 | path/to/BSNHM038-24.fasta |
| BSNHM046-24 | BSNHM046-24 | path/toBSNHM046-24.fasta |
* **Currently, it is crucial that the sample name (process_id), reference sequence FASTA file, and corresponding reference sequence FASTA header are all identical for correct sample-reference file mapping.**

## Customising snakemake configuration file ##
- Update [config/config.yaml](https://github.com/bge-barcoding/BGEE/blob/main/config/config.yaml) with the neccessary paths and variables.
- **Currently, BGEE's structural_validation only works for COI-5P and rbcL barcodes due to HMM availability (to be updated). BGEE's taxonomic_validation currently only works for COI-5P due to sole use of the BOLDistilled BLAST DB (to be updated).**
- Each rule in the config.yaml can specify the number of requested threads and memory resources (in Mb) for every job (e.g. specifying 4 threads and 4G memory for fastp_pe_merge would allocate those resources for every 'fastp_pe_merge' job).

# Cluster configuration using Snakemake profiles #
- See `profiles/` directory for 'slurm' and 'local' (i.e. non-SLURM) cluster submission parameters. The `jobs` parameter is likely the most important as it dictates how many workflow jobs can be run concurrently.
- The profile (`profiles/local` or `profiles/slurm`) will need to be changed depending on your system (see `$PROFILE` variable in `snakemake_run.sh`).

# Cluster submission #
- [snakemake_run.sh](https://github.com/bge-barcoding/MitoGeneExtractor-BGE/blob/main/snakemake_run.sh) handles submission of the snakemake workflow to the HPC cluster. The working directory will initially be unlocked (using `--unlock`) and then the snakemake workflow will be run. 


### If using `profiles/slurm`, SLURM will orchestrate submission of each step in the workflow as a separate job:
- The safer way is to launch Snakemake inside a persistent screen session. This ensures the workflow keeps running even if you disconnect.
```
# Start a persistent screen session
screen -S [SESSION_NAME]

# Allocate resources for Snakemake to use
salloc --job-name=[SESSION_NAME] \
       --partition=[YOUR_PARTITION] \
       --cpus-per-task=16 \
       --mem=8G
```
- Inside this allocation, you can launch Snakemake within an interactive session with srun:
```
srun ./snakemake_run.sh
```
- To detach from the screen session (disconnect but keep it running): `Ctrl + A + D`
- To reconnect again: `screen -r [SESSION_NAME]`
- You may see a warning such as "You are running snakemake in a SLURM job context. This is not recommended..." - This can generally be ignored because the salloc session is only acting as a submission manager. If you do encounter problems, try running `./snakemake_run.sh` within the screen session without running `salloc`.


### If using `profiles/local`, all workflow steps will be run as a single job:
-  Simply run `./snakemake_run.sh` on your desired cluster compute node. This node will handle all job scheduling and job computation.



# Results structure #
```
output_dir/
├── references/
│   ├── protein/                                               # Contains FASTA pseudo-references for each sample.
│   ├── genbank/                                               # Contains GenBank records corresponding to fetched pseudo-references (if 'genbank: true' set in config.yaml).
│   └── sequence_references.csv                                # Metadata on fetched pseudo-references.
│
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
- Please feel free to submit issues, fork the repository, and create pull requests for any improvements.
- This snakemake pipeline was produced by Dan Parsons @ NHMUK for the BioDiversity Genomics Europe (BGE) consortium.
- Since BGEE uses [MitogeneExtractor](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14075) at its core, please cite:
  Brasseur, M.V., Astrin, J.J., Geiger, M.F., Mayer, C., 2023. MitoGeneExtractor: Efficient extraction of mitochondrial genes from next-generation sequencing libraries. Methods in Ecology and Evolution.


  ## To do ##
- Split Snakefile into modular .smk files.
- Increase flexibility of input CSV headers (e.g. ID column in sample.csv and process_id column in sequence_references.csv could be ID/id/Process ID/PROCESS ID/process_id/sample/sample_id/SAMPLE ID/etc).
- Update 01_human_cox1_filter.py to not just filter against human coi, but the whole human genome/mitogenome.
  
