# Castanet fork - Python 3 with self-installation and convenience functions
Rich Mayne 2023
## Installation
1. We assume the user has already downloaded and installed Conda, and has created a fresh environment. For more details see https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html
1. 

# Castanet - Original Readme
Analysis of targeted metagenomics data as described in https://doi.org/10.1101/716902

Dependencies:
* python2.7+
* numpy
* pandas
* matplotlib
* samtools
* optional for mapping - bwa, or mapper of your choice
* optional for adapter removal - trimmomatic
* optional for pretty plots - seaborn
* optional for removing human reads - kraken (v1 or v2) with standard database of human + bacterial/viral RefSeq genomes 

Workflow:

0) Optionally pre-process your raw reads to remove human reads (genomic and mitochondrial) and common contaminants. Assuming you have run kraken on a sample called **SeqName** and saved the output to ${SeqName}.kraken.gz, the following creates a set of filtered fastq files, ${SeqName}_[12]_filt.fastq:
```bash
filter_keep_reads.py -i ${SeqName}_[12].fastq.gz  -k ${SeqName}.kraken.gz --xT Homo,Alteromonas,Achromobacter -x 1969841 --suffix filt
```

1) Trim adapters and poor-quality reads:
```bash
    trimmomatic PE -threads 8 SeqName_1_filt.fastq SeqName_2_filt.fastq SeqName_1_clean.fastq tmp_SeqName_1_trimmings.fq SeqName_2_clean.fastq tmp_SeqName_2_trimmings.fq ILLUMINACLIP:${AdaptersPath}:2:10:7:1:true MINLEN:80 
```

2) Map each sample's reads to the reference containing consensus target sequences with your mapper of choice, and produce a sorted BAM file of mapped reads only, for each sample. Do **not** deduplicate the BAM files.
```bash
    bwa index rmlst_virus_extra_ercc.fasta
    RefStem="rmlst_virus_extra_ercc.fasta"
    bwa mem ${RefStem} SeqName[12]_clean.fastq | samtools view -F4 -Sb -| samtools sort - 1> ${SeqName}.bam 
```

3) Generate a CSV file containing the counts for each uniquely mapped sequence in the entire pool, including improper pairs:
```bash
    for BamFilePath in $(ls \*.bam); do
        count_duplicates_single_bam_improper.sh ${BamFilePath}
    done > PosCounts.csv
```

4) Analyse all organisms detected in the pool, and combine with information on samples (including raw read number) and study participants (eg. clinical diagnosis):
```bash
process_pool_grp.py -h

usage: process_pool_grp.py [-h] -b BATCHNAME -i INFILE [-o OUTDIR]
                           [-p PROBELENGTHS] [-d] --samples SAMPLES
                           [--clin CLIN] [--depth_inf DEPTH_INF]

optional arguments:
  -h, --help            show this help message and exit
  -b BATCHNAME, --batchname BATCHNAME
                        Batch name for these samples. Must be alphanumeric.
  -i INFILE, --infile INFILE
                        Data frame (csv[.gz]) file to process. If gzipped,
                        filename must end in .gz.
  -o OUTDIR, --outdir OUTDIR
                        Output directory. Default: current working directory.
  -p PROBELENGTHS, --probelengths PROBELENGTHS
                        Path to file containing probe lengths. Default:
                        $CASTANET_PATH/scripts/probelengths_rmlst_virus_extra_ercc.csv
  -d, --keepdups        If true, do not reassign duplicates to the sample with
                        the majority in each duplicate cluster (Default:
                        False).
  --samples SAMPLES     Path to CSV file containing information about raw
                        reads (must have at least following fields: sampleid,
                        pt, rawreadnum). Field "pt" must match clinical data.
  --clin CLIN           Path to CSV file containing clinical data (must have
                        at least following fields: pt, clin_int; the field
                        "sampleid" if present will be ignored). Other fields
                        will be ignored.
  --depth_inf DEPTH_INF
                        (For regenerating full CSV with new clinical info):
                        Path to previously generated CSV file of read depth
                        per position for each probe, for all samples in this
                        batch. Must contain the following fields: sampleid,
                        target_id, depth_mean, depth_std, depth_25pc,
                        depth_median, depth_75pc, prop_target_covered,
                        prop_target_covered_mindepth2,
                        prop_target_covered_mindepth5,
                        prop_target_covered_mindepth10, udepth_mean,
                        udepth_std, udepth_25pc, udepth_median, udepth_75pc,
                        uprop_target_covered, uprop_target_covered_mindepth2,
                        uprop_target_covered_mindepth5,
                        uprop_target_covered_mindepth10

Example: process_pool_grp.py -i /path/to/my/dataframe.csv.gz --samples
/path/to/sampleinfo.csv
```
This produces three key outputs: a directory of coverage plots **Depth_BATCHNAME**, and the files **BATCHNAME_depth.csv** and **BATCHNAME_reads_to_drop.csv**.

BATCHNAME_depth.csv contains the number and proportion of reads for each sample and reference genome. Positives may need to be calibrated against a reference set, but in general, the proportion of all clean reads that match the given target (clean_prop_of_reads_on_target) is a good place to start.

BATCHNAME_reads_to_drop.csv can be used to clean BAM files to remove misassigned reads, for downstream processing.

5) If required, you can filter BAM files of interest to remove reads marked as contamination due to index misassignment:
```bash
   samtools -h ${SAMPLENAME}.bam | filter_bam.py SAMPLENAME BATCHNAME_reads_to_drop.csv
```
