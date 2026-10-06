# ONT DNA Barcoding: Consensus Generation and NCBI BLAST Analysis

This repository contains scripts for processing Oxford Nanopore Technologies (ONT) DNA barcoding data after demultiplexing with **ONTbarcoder 2.0**. The workflow generates sequence clusters and consensus sequences from demultiplexed ONT reads and subsequently performs similarity searches against the NCBI nucleotide database using BLASTn.

The workflow was developed for DNA barcoding/metabarcoding datasets, particularly COI amplicons, generated using pairs of indexed primers (Srivatsan et al., 2022), but can be applied to other sufficiently uniform amplicon datasets.

The motivation for this script was to use the whole information and recover all barcodes from intrinistically contaminated invertebrate DNA samples, when digested with the Lysis C buffer or Hot Shot buffer DNA sample contain not only DNA from the main study object, but prey, parasites, symbionts, enviromental DNA, which may not be the main goal of sequencing, but informs about ecological ties between species. In addition, universal barcoding primers may not suite well to the study organism or it may have minimal amount of DNA, which makes the barcoding workflow sensitive to background contamination. Detection and dientangling this contamination and simultaneous recovery of minor authentic barcode is possible scenario of barcoding problematic samples.

Additional application of indexed barcoding workflow is using indexed primers for amplification of eDNA for the purpose of metabarcoding. Mixed libraries from different repeats of same sample or different samples may be indexed and multiplexed on the stage of sequencing with later demultiplexing using ONTbarcoder 2.0, clustering of reads originating from the same index combination, and identification of consensus sequences instead of identification of reads.

## Workflow

The analysis consists of two principal steps:

```
\`ONT FASTQ\`  
  
\`   │\`  
  
\`   ▼\`  
  
\`ONTbarcoder 2.0\`  
  
\`demultiplexing\`  
  
\`   │\`  
  
\`   ▼\`  
  
\`Demultiplexed ONTbarcoder files\` (in default folder of ONTbarcoder output)  
  
\`   │\`  
  
\`   │  make\_consensus.sh\`  
  
\`   ▼\`  
  
\`Sequence clustering\`  
  
\`   │\`  
  
\`   ├── representative centroids\`  
  
\`   │\`  
  
\`   └── cluster consensus sequences\`  
  
\`             │\`  
  
\`             ▼\`  
  
\`       all\_consensus.fasta\`  
  
\`             │\`  
  
\`             │  blast\_all.py\`  
  
\`             ▼\`  
  
\`       NCBI BLASTn\`  
  
\`             │\`  
  
\`             ▼\`  
  
\`       Top 10 hits/query\`  
  
\`       blast\_results.csv\`
```

## Input data

The input files are the results of **demultiplexing ONT FASTQ reads using ONTbarcoder 2.0**.

ONTbarcoder produces demultiplexed sequence files corresponding to individual samples/barcodes. The files contain ONT read information together with the nucleotide sequence. In the ONTbarcoder output format, the nucleotide sequence is stored as the final whitespace-separated field of each record.

Example input files may have names such as:

```
\`P1\_A\_1\_all.fa\`  
  
\`P1\_A\_2\_all.fa\`  
  
\`P1\_A\_3\_all.fa\`  
  
\`...\`
```

The input data therefore represent **already basecalled and demultiplexed ONT reads**. Neither script performs basecalling, barcode demultiplexing, adapter trimming, or primary ONT signal processing.

### Sequence filtering

ONTbarcoder output can contain records that do not represent nucleotide sequences, including low-quality reads containing non-IUPAC amino-acid-like characters such as:

```
\`EPFPQEQQEQEEEQEFEFE...\`
```

`make\_consensus.sh` extracts the nucleotide sequence from the ONTbarcoder record and retains sequences composed of valid nucleotide characters (`A`, `C`, `G`, `T`, and `N`) for downstream clustering and consensus generation.

This filtering is performed before consensus generation because non-nucleotide records cannot contribute meaningfully to a nucleotide consensus.

## `make\_consensus.sh`

### Purpose

`make\_consensus.sh` processes the demultiplexed ONTbarcoder files and generates:

1. representative centroid sequences;

2. cluster consensus sequences;

3. multiple-sequence alignments for the clusters;

4. combined FASTA files containing the resulting sequences.

Sequence clustering is performed with **VSEARCH**.

### Clustering

Reads are clustered using:

```
\`vsearch --cluster\_fast \`  
  
\`    INPUT \`  
  
\`    --id 0.75 \`  
  
\`    --strand both \`  
  
\`    --centroids CENTROIDS \`  
  
\`    --sizeout\`
```

The current workflow therefore uses:

- **minimum pairwise identity:** 75%

- **strand:** both

- **clustering method:** VSEARCH `cluster\_fast`

- **cluster abundance:** retained in sequence headers using `size=`

- **minimum cluster size retained for downstream consensus:** 3 reads

The relatively permissive 75% clustering threshold is intended to group reads into broad sequence clusters before generating representative sequences. The threshold should be reconsidered for applications requiring stricter haplotype or species-level clustering.

### Centroids

For each cluster, VSEARCH selects a representative sequence (centroid). These sequences are retained separately from the consensus sequences.

Centroids are useful for:

- checking cluster composition;

- inspecting representative ONT reads;

- identifying problematic clusters;

- comparing centroid and consensus sequences;

- troubleshooting clustering results.

Example centroid header:

```
\`\>95beb4db-31f6-4c81-a624-9d2c57b3b1bc;size=1264\`  
  
\`ATGTCACCACAAACAGAGACT...\`
```

The `size` value indicates the number of reads assigned to the cluster.

### Consensus sequences

For each cluster meeting the minimum size criterion, a consensus sequence is generated from the reads belonging to that cluster.

Consensus sequences are used as the principal input for subsequent database searches because they are intended to reduce the influence of individual ONT sequencing errors and provide a sequence representing the cluster rather than a single raw read.

The consensus output is suitable for downstream:

- BLASTn searches;

- taxonomic identification;

- comparison with reference databases;

- manual sequence validation;

- downstream phylogenetic analyses, subject to appropriate quality control.

### Multiple-sequence alignments

Cluster alignments are retained to allow inspection of the reads contributing to each consensus.

These files are particularly useful for detecting:

- heterogeneous clusters;

- indels;

- strongly divergent reads;

- sequencing artefacts;

- potentially chimeric clusters;

- cases where a consensus may not provide an appropriate representation of the underlying reads.

### Minimum cluster size

Only clusters containing **at least three reads** are retained for consensus analysis.

Singletons and two-read clusters are excluded from the final consensus dataset because a consensus based on very few ONT reads provides limited error correction and is less robust for downstream taxonomic identification.

The threshold can be modified in the script if a different experimental design requires it.

## Outputs of `make\_consensus.sh`

The exact directory structure depends on the script configuration, but the principal outputs are:

```
\`centroids/\`  
  
\`    \<sample\>\_centroids.fa\`  
  
  
\`consensus/\`  
  
\`    \<sample\>\_consensus.fa\`  
  
  
\`msa/\`  
  
\`    \<sample\>\_\*.msa\`  
  
  
\`all\_centroids.fasta\`  
  
\`all\_consensus.fasta\`
```

### `all\_centroids.fasta`

Contains representative centroid sequences from retained clusters.

### `all\_consensus.fasta`

Contains consensus sequences from retained clusters and is the principal input file for `blast\_all.py`.

### Sequence identifiers

Cluster/sample information is retained in the sequence identifiers so that consensus sequences can be traced back to their originating demultiplexed sample and cluster.

This traceability is important because BLAST results alone should not be treated as the complete analytical record; each database identification should remain linked to the original ONT reads and cluster from which the consensus was generated.

# `blast\_all.py`

## Purpose

`blast\_all.py` performs automated similarity searches of the consensus sequences against the **NCBI nucleotide (`nt`) database** using NCBI's BLAST web service.

The script is designed for relatively large batches of barcode sequences and processes queries sequentially rather than submitting the entire dataset as a single BLAST request.

For each consensus sequence, the script:

1. submits a `blastn` search to NCBI;

2. obtains the NCBI Request ID (RID);

3. waits for the search to complete;

4. retrieves the BLAST XML result;

5. parses the result with Biopython;

6. extracts the top 10 database hits;

7. writes the results to a CSV file;

8. stores the raw BLAST XML result for reproducibility.

## BLAST parameters

The searches use:

```
\`Program:      blastn\`  
  
\`Database:     nt\`  
  
\`Maximum hits: 10\`  
  
\`Output:       XML\`
```

The BLAST results are subsequently parsed to obtain information including:

- query identifier;

- query length;

- hit rank;

- NCBI accession;

- hit description/title;

- percentage identity;

- alignment length;

- query coverage;

- number of mismatches;

- number of gaps;

- E-value;

- bit score.

Example output fields:

```
\`query\_id\`  
  
\`query\_length\`  
  
\`rank\`  
  
\`accession\`  
  
\`title\`  
  
\`identity\_percent\`  
  
\`alignment\_length\`  
  
\`query\_coverage\_percent\`  
  
\`mismatches\`  
  
\`gaps\`  
  
\`evalue\`  
  
\`bitscore\`
```

## Output

The principal output is:

```
\`blast\_results.csv\`
```

Each query can produce up to ten BLAST hits.

A separate file is generated for failed searches:

```
\`blast\_failed.csv\`
```

Raw BLAST XML results are stored individually, allowing BLAST results to be re-parsed without submitting the queries again.

Example:

```
\`blast\_xml/\`  
  
\`    P1\_A\_2\_cluster\_1.xml\`  
  
\`    P1\_A\_2\_cluster\_2.xml\`  
  
\`    P1\_A\_2\_cluster\_3.xml\`  
  
\`    ...\`
```

## Resume capability

`blast\_all.py` is designed to support interruption and continuation of long-running BLAST analyses.

Completed query identifiers are read from the existing results file, allowing already processed sequences to be skipped when the script is restarted.

This is important because remote NCBI BLAST searches can take substantially longer than local BLAST searches and may occasionally fail or remain pending.

For a run that was interrupted at a particular query, the script can additionally be configured to start from a specified position.

For example:

```
\`START\_FROM = 318\`
```

The exact value should correspond to the desired position in the input FASTA file.

For a complete restart:

```
\`START\_FROM = 1\`
```

Existing result files should be preserved when continuation is required.

## NCBI BLAST service

The script uses the NCBI BLAST HTTP interface rather than the local `blastn -remote` command.

This approach separates submission and result retrieval:

```
\`submit query\`  
  
\`     │\`  
  
\`     ▼\`  
  
\`NCBI returns RID\`  
  
\`     │\`  
  
\`     ▼\`  
  
\`poll RID until completed\`  
  
\`     │\`  
  
\`     ▼\`  
  
\`retrieve XML\`  
  
\`     │\`  
  
\`     ▼\`  
  
\`parse with Biopython\`  
  
\`     │\`  
  
\`     ▼\`  
  
\`append results to CSV\`
```

A delay between submissions is included to avoid unnecessarily aggressive requests to the NCBI service.

The script should therefore be regarded as a **batch submission client for the NCBI BLAST service**, not as a replacement for a local BLAST installation.

# Prerequisites

## Operating system

The workflow has been developed and tested under:

```
\`Ubuntu Linux 22.04\`
```

Other Linux distributions should work with appropriate versions of the dependencies.

## Hardware

No specialized hardware is required.

The BLAST step is executed by NCBI's remote BLAST service, so GPU acceleration is not used.

The consensus-generation step is CPU-based and can run on ordinary desktop/workstation hardware.

For large datasets, sufficient disk space should be available for:

- input ONTbarcoder files;

- intermediate clustering files;

- multiple-sequence alignments;

- consensus sequences;

- raw BLAST XML results.

## Required software

### 1. ONTbarcoder 2.0

ONTbarcoder 2.0 is required to generate the demultiplexed input files.

The scripts assume that:

- ONT FASTQ data have already been basecalled;

- reads have already been demultiplexed by ONTbarcoder;

- the resulting demultiplexed files are supplied as input.

The scripts do not perform ONTbarcoder processing themselves.

### 2. VSEARCH

`make\_consensus.sh` requires **VSEARCH**.

The workflow was tested with:

```
\`VSEARCH 2.21.1\`
```

VSEARCH must be accessible from the command line:

```
\`vsearch --version\`
```

Expected output should report the installed VSEARCH version.

### 3. Bash

`make\_consensus.sh` is a Bash shell script and requires a standard Unix shell environment.

Check:

```
\`bash --version\`
```

### 4. Python

`blast\_all.py` requires Python 3.

Python 3.10/3.11 is recommended for a reproducible environment.

Check:

```
\`python --version\`
```

### 5. Biopython

`blast\_all.py` uses Biopython for:

- reading FASTA sequences;

- parsing NCBI BLAST XML results.

Install with Conda:

```
\`conda install -c conda-forge biopython\`
```

Verify:

```
\`python -c "from Bio import SeqIO; from Bio.Blast import NCBIXML; print('Biopython OK')"\`
```

### 6. Requests

The script uses Python `requests` for communication with the NCBI BLAST HTTP service.

Install with:

```
\`conda install -c conda-forge requests\`
```

Verify:

```
\`python -c "import requests; print(requests.\_\_version\_\_)"\`
```

# Recommended Conda environment

A dedicated environment is recommended for reproducibility:

```
\`conda create -n ont\_blast python=3.11 biopython requests\`  
  
\`conda activate ont\_blast\`
```

VSEARCH can then be installed in the same environment:

```
\`conda install -c bioconda vsearch\`
```

Verify all dependencies:

```
\`python --version\`  
  
\`python -c "from Bio import SeqIO; from Bio.Blast import NCBIXML; import requests; print('Python dependencies OK')"\`  
  
\`vsearch --version\`
```

A complete dependency check should therefore produce successful results for:

```
\`Python\`  
  
\`Biopython\`  
  
\`Requests\`  
  
\`VSEARCH\`
```

# Running the workflow

## 1. Generate consensus sequences

Place the ONTbarcoder 2.0 demultiplexed files in the input directory and run:

```
\`bash make\_consensus.sh\`
```

The principal downstream file is:

```
\`all\_consensus.fasta\`
```

Before proceeding, it is recommended to inspect the number and approximate length of the resulting sequences:

```
\`grep -c "^\>" all\_consensus.fasta\`
```

and:

```
\`head all\_consensus.fasta\`
```

## 2. Run BLAST

Run:

```
\`python blast\_all.py\`
```

The script will process the consensus sequences sequentially and create:

```
\`blast\_results.csv\`  
  
\`blast\_failed.csv\`  
  
\`blast\_xml/\`
```

For a long analysis, it is recommended to run the script in a persistent terminal session such as `tmux` or `screen`, particularly on a workstation that may otherwise suspend.

Example:

```
\`tmux new -s blast\`  
  
\`python blast\_all.py\`
```

The session can be detached with:

```
\`Ctrl-B D\`
```

and reattached with:

```
\`tmux attach -t blast\`
```

# Reproducibility and interpretation

The workflow deliberately preserves several levels of sequence information:

```
\`raw ONT reads\`  
  
\`      │\`  
  
\`      ▼\`  
  
\`cluster assignments\`  
  
\`      │\`  
  
\`      ├── centroid\`  
  
\`      │\`  
  
\`      └── consensus\`  
  
\`              │\`  
  
\`              ▼\`  
  
\`          BLAST result\`
```

This allows a BLAST identification to be traced back to the sequence cluster and ultimately to the original ONTbarcoder demultiplexed reads.

The BLAST results should not be interpreted solely on the basis of the first database hit. Percentage identity, query coverage, alignment length, E-value, taxonomic consistency among hits, and the quality of the underlying consensus should all be considered.

In particular, a high-scoring BLAST match does not by itself demonstrate correct species identification. Reference database errors, incomplete taxonomic sampling, misidentified sequences, mitochondrial pseudogenes, contamination, chimeric reads, and insufficiently diagnostic barcode regions can all affect the result.

For publication-quality taxonomic identification, BLAST results should therefore be subjected to additional quality control and, where appropriate, phylogenetic or taxonomic validation.

## Exact workflow parameters

### `make\_consensus.sh`

`make\_consensus.sh` processes all ONTbarcoder 2.0 demultiplexed files matching the input-file pattern defined in the script.

For each input file, the script:

1. extracts the nucleotide sequence from each ONTbarcoder record;

2. removes records whose sequence contains characters outside `A`, `C`, `G`, `T`, and `N`;

3. writes the cleaned sequences to an intermediate FASTA file;

4. clusters the sequences with VSEARCH;

5. generates cluster centroids;

6. generates cluster consensus sequences;

7. retains clusters containing at least three reads;

8. combines the retained sequences from all input files into the final FASTA datasets.

The principal VSEARCH clustering operation is:

```
\`vsearch --cluster\_fast INPUT \`  
  
\`    --id 0.75 \`  
  
\`    --strand both \`  
  
\`    --centroids CENTROIDS \`  
  
\`    --consout CONSENSUS \`  
  
\`    --msaout MSA \`  
  
\`    --sizeout\`
```

#### Clustering parameters

| **Parameter** | Value | **Description** |
| :-: | -: | :-: |
| Clustering algorithm | `--cluster\_fast` | Greedy centroid-based clustering |
| Sequence identity | `--id 0.75` | Minimum pairwise identity for cluster assignment |
| Strand | `--strand both` | Both nucleotide orientations are considered |
| Centroids | `--centroids` | Writes one representative sequence per cluster |
| Consensus | `--consout` | Writes one consensus sequence per cluster |
| Alignment | `--msaout` | Writes cluster multiple-sequence alignments and consensus |
| Abundance | `--sizeout` | Adds cluster abundance to FASTA headers |
| Minimum retained cluster | 3 reads | Applied by the shell script after clustering |


`--cluster\_fast` sorts sequences by decreasing length before clustering. Because clustering is centroid-based and sequential, the order in which sequences are processed can influence which sequence becomes the centroid and consequently the resulting clusters.

The 75% identity threshold is intentionally permissive and is used to group related barcode reads prior to consensus generation. It should not be interpreted as a species-level identification threshold.

### Consensus generation

VSEARCH constructs the consensus from a center-star multiple sequence alignment using the cluster centroid as the alignment center. The consensus nucleotide at each alignment position is determined from the majority character, with columns dominated by gaps being removed from the resulting ungapped consensus.

The generated MSA files are therefore retained as an intermediate quality-control product rather than treating the consensus sequence as an unquestionable representation of the cluster.

This is particularly relevant for clusters containing divergent reads or reads with substantial indels. VSEARCH documentation notes that the center-star alignment approach can lose accuracy at low pairwise-identity thresholds.

### Cluster-size filtering

The script retains only clusters with:

```
\`cluster size \>= 3 reads\`
```

This filtering is performed **after VSEARCH clustering**.

Clusters containing one or two reads are excluded from the final consensus dataset because they provide little redundancy for correcting individual ONT sequencing errors.

The original cluster information is retained in the intermediate files, allowing this threshold to be changed without modifying the underlying raw ONTbarcoder data.

# `blast\_all.py`

`blast\_all.py` performs sequential BLASTn searches of the consensus sequences against the NCBI `nt` nucleotide database.

The script uses the **NCBI BLAST Common URL API** rather than the local `blastn -remote` command. NCBI's API uses `CMD=Put` to submit a search and returns a Request ID (RID), which is subsequently used with `CMD=Get` to check the status and retrieve the result.

### BLAST parameters

The relevant search parameters are:

| **Parameter** | **Value** |
| :-: | :-: |
| Program | `blastn` |
| Database | `nt` |
| Maximum reported hits | 10 |
| Result format | XML |
| Query type | nucleotide FASTA |
| Submission mode | NCBI BLAST URL API |
| Processing | sequential, one query at a time |


The NCBI API supports `blastn`, nucleotide databases, configurable hit-list size, and XML result output.

### Submission and polling

For every consensus sequence, the script follows this procedure:

```
\`FASTA sequence\`  
  
\`      │\`  
  
\`      ▼\`  
  
\`CMD=Put\`  
  
\`      │\`  
  
\`      ▼\`  
  
\`NCBI BLAST server\`  
  
\`      │\`  
  
\`      ▼\`  
  
\`Request ID (RID)\`  
  
\`      │\`  
  
\`      ▼\`  
  
\`poll CMD=Get\`  
  
\`      │\`  
  
\`      ▼\`  
  
\`BLAST completed\`  
  
\`      │\`  
  
\`      ▼\`  
  
\`retrieve XML\`  
  
\`      │\`  
  
\`      ▼\`  
  
\`Biopython NCBIXML parser\`  
  
\`      │\`  
  
\`      ▼\`  
  
\`top 10 hits\`  
  
\`      │\`  
  
\`      ▼\`  
  
\`blast\_results.csv\`
```

The script uses a delay between successive submissions and retries failed requests. This is important because NCBI's BLAST service is a shared public resource; NCBI specifically notes that projects involving large numbers of BLAST searches should consider the REST service, cloud providers, or standalone BLAST rather than making excessive use of the public service.

## BLAST result fields

For each query, the script records up to ten hits with the following fields:

```
\`query\_id\`  
  
\`query\_length\`  
  
\`rank\`  
  
\`accession\`  
  
\`title\`  
  
\`identity\_percent\`  
  
\`alignment\_length\`  
  
\`query\_coverage\_percent\`  
  
\`mismatches\`  
  
\`gaps\`  
  
\`evalue\`  
  
\`bitscore\`
```

These fields provide the basic information required for subsequent taxonomic assessment and quality control.

The first BLAST hit is **not automatically treated as a species identification**. The results should be evaluated together with sequence identity, query coverage, alignment length, taxonomic consistency among hits, and the quality of the underlying consensus.

# Resuming interrupted BLAST analyses

`blast\_all.py` is designed for long-running analyses and can resume an interrupted run.

The script checks the existing result file and skips query identifiers that have already been processed.

The recommended procedure is therefore:

```
\`python blast\_all.py\`
```

If the process is interrupted, simply restart the script after preserving:

```
\`blast\_results.csv\`  
  
\`blast\_failed.csv\`  
  
\`blast\_xml/\`
```

Previously completed queries will be skipped.

For manual continuation from a particular position, the `START\_FROM` parameter can be changed in the script:

```
\`START\_FROM = 318\`
```

For example:

```
\`START\_FROM = 1\`
```

starts processing from the first sequence, whereas:

```
\`START\_FROM = 318\`
```

starts from the 318th sequence in the input FASTA.

The resume mechanism based on existing query identifiers should be preferred where possible, because it is independent of the numerical position of a sequence in the input file.

# Reproducibility requirements

For reproducibility, the following information should be retained with each analysis:

```
\`ONTbarcoder version\`  
  
\`VSEARCH version\`  
  
\`Python version\`  
  
\`Biopython version\`  
  
\`Requests version\`  
  
\`make\_consensus.sh\`  
  
\`blast\_all.py\`  
  
\`input FASTA files\`  
  
\`all\_centroids.fasta\`  
  
\`all\_consensus.fasta\`  
  
\`cluster MSAs\`  
  
\`blast\_results.csv\`  
  
\`blast\_failed.csv\`  
  
\`raw BLAST XML files\`
```

The versions used during development were:

```
\`ONTbarcoder      2.0\`  
  
\`VSEARCH          2.21.1\`  
  
\`Python           3.x\`  
  
\`Biopython        current installed version\`  
  
\`Requests         current installed version\`  
  
\`Ubuntu           22.04\`
```

Exact package versions should preferably be recorded at the time of a final analysis:

```
\`python --version\`  
  
\`python -c "import Bio, requests; print('Biopython:', Bio.\_\_version\_\_); print('Requests:', requests.\_\_version\_\_)"\`  
  
\`vsearch --version\`
```

For Conda-based analyses, the complete environment can be exported with:

```
\`conda env export \> environment.yml\`
```

The resulting `environment.yml` should be committed to the repository or archived with the analysis to facilitate reproduction.

# Recommended project structure

A reproducible analysis can be organized as:

```
\`project/\`  
  
\`├── scripts/\`  
  
\`│   ├── make\_consensus.sh\`  
  
\`│   └── blast\_all.py\`  
  
\`│\`  
  
\`├── input/\`  
  
\`│   └── ONTbarcoder\_demultiplexed/\`  
  
\`│       ├── P1\_A\_1\_all.fa\`  
  
\`│       ├── P1\_A\_2\_all.fa\`  
  
\`│       └── ...\`  
  
\`│\`  
  
\`├── intermediate/\`  
  
\`│   ├── cleaned/\`  
  
\`│   ├── centroids/\`  
  
\`│   ├── consensus/\`  
  
\`│   └── msa/\`  
  
\`│\`  
  
\`├── results/\`  
  
\`│   ├── all\_centroids.fasta\`  
  
\`│   ├── all\_consensus.fasta\`  
  
\`│   ├── blast\_results.csv\`  
  
\`│   ├── blast\_failed.csv\`  
  
\`│   └── blast\_xml/\`  
  
\`│\`  
  
\`├── environment.yml\`  
  
\`└── README.md\`
```

This separation keeps the original ONTbarcoder output untouched and makes it possible to reconstruct the analysis from the original demultiplexed reads.

# Important methodological note

The workflow has two distinct analytical stages that should not be conflated:

**Sequence processing and error reduction**

```
\`ONT reads\`  
  
\`   ↓\`  
  
\`clustering\`  
  
\`   ↓\`  
  
\`cluster consensus\`
```

and **taxonomic similarity searching**

```
\`consensus sequence\`  
  
\`   ↓\`  
  
\`NCBI BLASTn\`  
  
\`   ↓\`  
  
\`reference-sequence similarity\`
```

The BLAST result is therefore a database similarity result rather than an independent taxonomic assignment algorithm. Final species-level identification requires appropriate reference sequences and should consider the possibility of erroneous or misidentified database records.

For publication, the clustering threshold, minimum cluster size, consensus-generation method, and BLAST database/search date should be explicitly reported.

# Summary

The complete workflow is:

```
\`ONT FASTQ\`  
  
\`    ↓\`  
  
\`ONTbarcoder 2.0\`  
  
\`demultiplexing\`  
  
\`    ↓\`  
  
\`Demultiplexed ONTbarcoder files\`  
  
\`    ↓\`  
  
\`make\_consensus.sh\`  
  
\`    ↓\`  
  
\`VSEARCH clustering\`  
  
\`    ↓\`  
  
\`Clusters ≥3 reads\`  
  
\`    ├───────────────┐\`  
  
\`    ↓               ↓\`  
  
\`Centroids       Consensus\`  
  
\`    │               │\`  
  
\`    │               ↓\`  
  
\`    │        all\_consensus.fasta\`  
  
\`    │               │\`  
  
\`    │               ↓\`  
  
\`    │         blast\_all.py\`  
  
\`    │               │\`  
  
\`    │               ↓\`  
  
\`    │          NCBI BLASTn\`  
  
\`    │               │\`  
  
\`    │               ↓\`  
  
\`    │        Top 10 hits/query\`  
  
\`    │               │\`  
  
\`    │               ↓\`  
  
\`    │        blast\_results.csv\`  
  
\`    │\`  
  
\`    └── all\_centroids.fasta\`
```

The combination of cluster centroids, consensus sequences, alignments, raw BLAST XML files, and the final tabular BLAST output provides an auditable sequence-identification workflow while retaining the connection between individual ONT reads and their downstream taxonomic assignments.

