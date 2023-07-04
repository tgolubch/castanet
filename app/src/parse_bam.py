from __future__ import print_function
import sys
import re
import os
import pandas as pd

from app.utils.argparsers import parse_args_bam_parse
from app.utils.error_handlers import error_handler_parse_bam_positions
from app.utils.shell_cmds import shell


class Parse_bam_positions:
    '''
    Parse contents of bam file, via reading shell commands passed in.
    Should only get called by another Python script due to requirement for shell input.
    '''

    def __init__(self, argies) -> None:
        self.min_match_length = 40
        self.argies = argies
        self.consensus_reads = {}

    def getmatchsize(self, cigar):
        '''Find matches in cigar string with regex, return count'''
        matches = re.findall(r'([0-9]+)M', cigar)
        if not len(matches):
            return 0
        return sum(int(x) for x in matches)

    def get_gene_orgid(self, target_id):
        '''Find gene and orgid at specific ref; return (gene, orgid).'''
        parts = target_id.split('_')
        return (parts[0], parts[-1] if parts[0].startswith('BACT') else parts[0])

    def build_consensus_db(self, ref, seq, improper_read=False):
        if not ref in self.consensus_reads.keys():
            self.consensus_reads[ref] = [seq]
        else:
            self.consensus_reads[ref].append(seq)

    def parse_bam_position(self, l):
        '''GENERATE COUNTS STAGE: For each line passed in, parse fields of interest; identify matches and print back to stdout.'''
        fields = l.split()
        ref, pos, ref2, tlen, cigar, seq = fields[2], fields[3], fields[6], int(
            fields[8]), fields[5], fields[9]
        if self.argies.Mode == "parse":
            if tlen >= self.min_match_length:
                '''Properly paired and match is of decent mapped length'''
                print(f'{ref},{pos},{tlen},{self.argies.SeqName}')
                self.build_consensus_db(self.get_gene_orgid(ref)[1], seq)

            elif (tlen == 0) and (self.getmatchsize(cigar) >= self.min_match_length) and (self.get_gene_orgid(ref) == self.get_gene_orgid(ref2)):
                '''Improperly paired BUT same gene AND match is of decent mapped length (via CIGAR string lookup) AND RNAME ref organism is same to RNEXT ref org'''
                # RM < TODO De facto we smash all hits from different seqs for same organism together.
                print(f'{ref},{pos},{tlen},{self.argies.SeqName}')
                self.build_consensus_db(self.get_gene_orgid(ref)[1], seq)
        else:
            print(f'{ref},{pos},{tlen},{self.argies.SeqName}')
            self.build_consensus_db(self.get_gene_orgid(ref)[1], seq)

    def filter_bam(self, l, reads_to_drop):
        '''POST FILTER STAGE'''
        if l.startswith('@'):
            sys.stdout.write(l)
            return
        if len(reads_to_drop):
            fields = l.split()
            ref = fields[2]
            pos = int(fields[3])
            tlen = int(fields[8])
            '''output only reads that are not found in the indexed READS_TO_DROP file'''
            if not (ref, pos, tlen) in reads_to_drop.index:
                sys.stdout.write(l)
            else:
                return

    def main(self):
        '''Entrypoint. Multi functional across generate counts and post filter.'''
        error_handler_parse_bam_positions(sys.argv)
        if self.argies.Mode == "filter":
            try:
                reads_to_drop = pd.read_csv(
                    self.argies.ExpDir + self.argies.FilterFile)
            except:
                raise SystemExit(
                    f"Couldn't find reads to drop file: {self.argies.FilterFile}. Did your run generate one?")

        # for l in sys.stdin:
        #     if self.argies.Mode == "parse":
        #         self.parse_bam_position(l, self.argies["SeqNane"])
        #     else:
        #         self.filter_bam(l, reads_to_drop)

        if self.argies.Mode == "parse":
            file = "./bamview.txt"
        elif self.argies.Mode == "reparse":
            file = "./test/remapped_full.sam"
        else:
            raise SystemExit(
                "IM RUNNING IN DEV MODE, CAN'T RUN IN FILTER MODE")

        with open(file) as f:
            for l in f:
                if self.argies.Mode == "parse":
                    self.parse_bam_position(l)
                elif self.argies.Mode == "reparse":
                    self.parse_bam_position(l)
                else:
                    self.filter_bam(l, reads_to_drop)

        if len(self.consensus_reads) == 0:
            return
        else:
            self.save_consensus_db()

    def save_consensus_db(self):
        grp_aln_f = f"experiments/{self.argies.SeqName}/grouped_alignments/"
        remap_f = f"experiments/{self.argies.SeqName}/grouped_alignments/remaps/"
        if not os.path.isdir(grp_aln_f):
            shell(f"mkdir {grp_aln_f}")
            shell(f"mkdir {remap_f}")

        if self.argies.Mode == "parse":
            folder = grp_aln_f
        else:
            folder = remap_f

        for key in self.consensus_reads.keys():
            with open(f"{folder}{key}.fasta", "w") as file:
                [file.write(f">{i+1}\n{self.consensus_reads[key][i]}\n")
                 for i in range(len(self.consensus_reads[key]))]


if __name__ == '__main__':
    cls = Parse_bam_positions(parse_args_bam_parse())
    cls.main()
