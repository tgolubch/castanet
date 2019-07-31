# Castanet
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

0) Optional pre-processing to remove human reads, assuming you have run kraken on a sample called **SeqName** and saved the output to a file called SeqName.kraken.gz:
```bash
filter_keep_reads.py -i SeqName_[12].fastq.gz -k SeqName.kraken.gz --xT Homo,Alteromonas,Achromobacter --suffix filt
```

1) Trim adapters and poor-quality reads:
```bash
    trimmomatic PE -threads 8 SeqName_1_filt.fastq SeqName_2_filt.fastq SeqName_1_clean.fastq tmp_SeqName_1_trimmings.fq SeqName_2_clean.fastq tmp_SeqName_2_trimmings.fq ILLUMINACLIP:${AdaptersPath}:2:10:7:1:true MINLEN:80 
```

2) Map each sample's reads to the reference containing consensus target sequences with your mapper of choice, and produce a sorted BAM file of mapped reads only, for each sample. Do **not** deduplicate the BAM files.
```bash
    bwa index rmlst_virus_extra_ercc.fasta
    RefStem="rmlst_virus_extra_ercc.fasta"
    bwa mem ${RefStem} SeqName[12]_clean.fastq | samtools view -F4 -Sb -| samtools sort - 1> temp_SeqName.bam 
```

3) Generate a CSV file containing the counts for each uniquely mapped sequence in the entire pool:
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
                        /well/resgen2/users/golubchi/CHIMES/PAPER/ChiMES-
                        GAinS/scripts/probelengths_rmlst_virus_extra_ercc.csv
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
