import subprocess as sp
import argparse 

class E2eRunner:
    '''Complete pipeline using functionality of original scripts'''
    def __init__(self) -> None:
        self.a       = self.parse_args()
        self.aliases = self.initiate_aliases()

    def shell(self, args, executable='/bin/bash'):
        sp.run(args, text=True,shell=True, executable=executable)

    def initiate_aliases(self):
        '''Load aliases into shell'''
        aliases = {}
        aliases["trim"] = 'java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar'
        aliases["bwa"]  = './bwa-mem2-2.2.1_x64-linux/bwa-mem2'
        return aliases

    def parse_args(self):
        '''Parse command-line arguments'''
        p = argparse.ArgumentParser(description=__doc__, epilog='Complete pipeline using functionality of original scripts')
        p.add_argument( '-ExpDir',  default="data/",  required=True,   help='Provide directory for input data' )
        p.add_argument( '-SeqName', default=None,     required=True,   help='Provide input sequence (experiment) name.' )
        p.add_argument( '-ExpName', default='myexp',  required=False,  help='Name your experiment batch.' )
        p.add_argument( '-AdaptP',  default="Trimmomatic-0.39/adapters/all.fa", required=False, help='Location of your Trimmomatic adapter sequences - usually in Trimmomatic path.' )
        p.add_argument( '-RefStem', default="",       required=False,  help="Path to mapping file (fasta)")
        p.add_argument( '-PostFilt',default=False,    required=False,  help="Filter BAM file to remove reads marked as contamination")
        p.add_argument( '-Samples',                   required=True,   help="CSV file containing sample information. Absolute path required")
        p.add_argument( '-Probes',  default="",       required=False,  help="CSV file containing probe length mappings. Absolute path required")
        return p.parse_args()
    
    def run_kraken(self):
        '''Call Kraken2 to remove unwanted reads'''
        self.shell(f'kraken2 --db kraken2_human_db/ --threads 4 {self.a.ExpDir}{self.a.SeqName}_1.fastq.gz > {self.a.ExpDir}{self.a.SeqName}_1.kraken')

    def filter_keep_reads(self):
        '''Filter reads using Castanet'''
        self.shell(f"python3 src/filter_keep_reads.py -i {self.a.ExpDir}{self.a.SeqName}_[12].fastq.gz -k {self.a.ExpDir}{self.a.SeqName}_1.kraken --xT Homo,Alteromonas,Achromobacter -x 1969841 --suffix filt")

    def trim(self):
        '''Trim adapters and bad reads'''
        self.shell(f"{self.aliases['trim']} PE -threads 8 {self.a.ExpDir}{self.a.SeqName}_1_filt.fastq {self.a.ExpDir}{self.a.SeqName}_2_filt.fastq {self.a.ExpDir}{self.a.SeqName}_1_clean.fastq {self.a.ExpDir}{self.a.SeqName}_1_trimmings.fq {self.a.ExpDir}{self.a.SeqName}_2_clean.fastq {self.a.ExpDir}{self.a.SeqName}_2_trimmings.fq ILLUMINACLIP:{self.a.AdaptP}:2:10:7:1:true MINLEN:80")

    def map(self):
        '''Use BWA and Samtools to map reads from each sample to targets'''
        self.shell(f"{self.aliases['bwa']} index {self.a.ExpDir}{self.a.RefStem}")
        self.shell(f"{self.aliases['bwa']} mem {self.a.ExpDir}{self.a.RefStem} {self.a.ExpDir}{self.a.SeqName}_[12]_clean.fastq | samtools view -F4 -Sb - | samtools sort - 1> {self.a.ExpDir}${self.a.SeqName}.bam")

    def count_mapped(self):
        '''Generate CSV of counts for each uniq mapped seq'''
        self.shell(f"""for BamFilePath in $(ls data/*.bam); do BamPath=$BamFilePath; BamName=$(basename "BamPath%%.bam"); BamName=$(sed s'/_dedup//' <<< ${{BamName}}); samtools view -F2048 -F4 ${{BamPath}} | python3 src/parse_bam_positions.py {self.a.SeqName} | sort | uniq -c | sed s'/ /,/'g | sed s'/^[,]*//'g; done > {self.a.ExpDir}PosCounts.csv""")

    def analysis(self):
        '''Call main Castanet analysis script'''
        self.shell(f"python3 src/process_pool_grp.py -i {self.a.ExpDir}PosCounts.csv --samples {self.a.Samples} -p {self.a.Probes} -b {self.a.ExpName}")

    def post_filter(self):
        '''OPT: Filter BAM files of interest to remove reads marked as bad'''
        self.shell(f"samtools view -h {self.a.ExpDir}{self.a.SeqName}.bam | python3 src/filter_bam.py {self.a.ExpDir}{self.a.SeqName} {self.a.ExpName}_reads_to_drop.csv")

    def main(self):
        '''Entrypoint'''
        self.initiate_aliases()
        self.run_kraken()
        self.filter_keep_reads()
        self.trim()
        self.map()
        self.count_mapped()
        self.analysis()
        if self.a.PostFilt:
            self.post_filter()

if __name__ == "__main__":
    '''Example input
    python3 app/end_to_end_script.py -ExpDir data/ -SeqName ERR10812890 -RefStem rmlst_virus_extra_ercc_2018.fasta -PostFilt True -Samples data/samples.csv -Probes data/probelengths_rmlst_virus_extra_ercc.csv
    '''
    cls = E2eRunner()
    cls.main()
