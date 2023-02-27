<h1 align="center">
  Single-Molecule Methyladenosine Finder on SLURM cluster
  <br>
</h1>

This program predicts single-molecule 

<p align="center">
  <a href="#key-features">Key Features</a> •
  <a href="#how-to-use">How To Use</a> •
  <a href="#input-file-format">Input File Format</a> •    
  <a href="#output-file-format">Output File Format</a> •  
  <a href="#contact">Contact</a> 
</p>


## Key Features

* Fast analysis pipeline of methyladenosine sites from PacBio subreads on SLURM cluster.
* Consensus N6-methyladenosine sites on circular consensus reads
* BAM format ready for IGV viewer
* BED format for easy downstream analysis

## How To Use
To clone and run this application on slurm cluster, you'll need [SMRTLINK](https://www.pacb.com/support/software-downloads/) deployed on SLURM cluster; [pysam](https://pysam.readthedocs.io/), [Biopython](https://biopython.org/) and [samtools](http://www.htslib.org/) are also required.

  [Installation](#Installation) • [Preprocess](#Preprocess) • [Predict m6A sites](#predict-m6a-sites)


### Installation
```bash
# Clone this repository
$ git clone https://github.com/zhuweix/MethyladenosineFinder.git

# Install dependencies using Conda (Except SMRTLink)
# For SMRTLink installation, please follow the documentation on PacBio website (https://www.pacb.com/support/software-downloads/)

$ conda install -c bioconda pysam biopython samtools

# Go into the repository
$ cd MethyladenosineFinder
# setup the python package
$ python setup.py setup
```
### Preprocess
```bash
# prepare reference index
$ m6apreprocess index \
$     --r [Reference.FASTA] \
$     --output [Reference_index.data]

# Recommended: Filter the Subreads using the PacBio CCS reads
# Get ZMW list by num_pass and AverageQual in PacBio CCS reads
$ m6apreprocess zmwfromccs \
$     --bam [PacBio_CCS.bam] \
$     --output [ZMW_list.txt] \
$     --qual [Min_AveQual, Default=90] \
$     --cov [Min_NumPass, Default=3]
# Filter subreads using the ZMW list (alternatively, bamsieve tool from PacBio SMRT Tools could also be used)
$ m6apreprocess filtersubreads \
$     --bam [PacBio_Subreads.bam] \
$     --zmw [ZMW_list.txt] \
$     --output [Filtered_Subreads.bam]   

# Align the Subreads to Reference using pbmm2 from SMRT Tools
$ pbmm2 index [Reference.FASTA] [Reference.mmi]
$ pbmm2 align [Reference.mmi] [PacBio_Subreads.bam | Filtered_Subreads.bam] \
$     [Output_Aligned_Sorted_Subreads.bam] --preset SUBREAD --sort \
$     -j 16 [Threads for alignment] \
$     -J 8 [Threads for sorting]
# index the aligned bam file using pbindex from SMRT Tools
$ pbindex [Output_Aligned_Sorted_Subreads.bam]
$ samtools index [Output_Aligned_Sorted_Subreads.bam]

```
### Predict m6A sites
```bash
# Generate the batch script
$ m6aworkflow  \
$    --bam [AlignedSubreads.bam | FilteredSubreads.bam] \
$    --sbatchdir [Folder for SBATCH files for SLURM submission] \    
$    --outputdir [Folder for Output BAM files] \
$    --tmpdir [Folder for Temporary Files] \
$    --jobname [JobName For Slurm] \
$    --ref [Reference.fa] \
$    --prefix [Prefix] \ # Prefix of the sbatch files and output files; Default=SMF
$    --cov [MinCoverage] \ # Default=8
$    --batchsize [BatchSize] \# Number of ZMWs per SLURM array job; Default=500
$    --thread [Threads]  \ # Default=16

# IMPORTANT NOTICE: The number of Batch array jobs submitted to SLURM is [Number of ZMWs] / [BatchSize], very small batchsize will result in a large amount of SLURM jobs

# Follow the workflow to submit the batch script
# The commands are stored in [BatchScriptFolder]/[Prefix].workflow.pipeline
$ sbatch [split.sh]
$ sbatch [predict_m6a.sh]
$ sbatch [merge.sh]

# [Optional] Transform BAM file to BED file for m6A sites


```
> **Note**
> The current workflow is tested with SMRTLink/11.1.0.166339 

## Input File Format
* BAM file of PacBio subreads aligned to the reference genome
    - We tested our pipeline with the aligned BAM file using pbalign for SMRTLINK
```bash
# Example command
$ pbmm2 index [refererence.fa] [reference.mmi]

$ pbmm2 align \
    [reference.mmi] \
    [RawSubreads.bam] \
    [AlignedSubreads.bam] \
    --sort -j 32 -J 8 --preset SUBREAD 

$ samtools index -@ 8 [AlignedSubreads.bam]

```
* [Optional] BAM file of PacBio CCS reads
  
  We use the average base quality score to filter the raw subreads. The default setting is to filter ZMWs with >= 90 average base quality and >= 3 passes in the CCS reads. This process is highly recommended to reduce the ZMWs with insufficient quality for methyladenosine prediction.

## Output File Format
* BAM format

This format is ready for IGV viewer. The methyladenosine sites are indicated as insertions.

* BED format

This format is easy for further downstream analyze. Each line represents the consensus methyladenosine sites from one ZMW. The name of each line is {movieName}/{holeNumber}/ which is the same in PacBio BAM file. For subreads within one ZMW that aligned to multiple genome locations, the name of the associated line is {movieName}/{holeNumber}/{Number} where Numer=1,2,3... for each location.







---
## Contact
> [David Clark Lab](https://www.nichd.nih.gov/research/atNICHD/Investigators/clark)
> 
> GitHub [@zhuweix](https://github.com/zhuweix)

