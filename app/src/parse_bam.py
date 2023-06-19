from __future__ import print_function
import sys
import re
import pandas as pd

from app.utils.argparsers import parse_args_bam_parse
from app.utils.error_handlers import error_handler_parse_bam_positions


class Parse_bam_positions():
    '''
    Parse contents of bam file, via reading shell commands passed in.
    Should only get called by another Python script due to requirement for shell input.
    '''

    def __init__(self, argies) -> None:
        self.min_match_length = 40
        self.argies = argies

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

    def parse_bam_position(self, l, sampleid):
        '''GENERATE COUNTS STAGE: For each line passed in, parse fields of interest; identify matches and print back to stdout.'''
        fields = l.split()
        ref, pos, ref2, tlen = fields[2], fields[3], fields[6], int(fields[8])
        if tlen >= self.min_match_length:
            '''Properly paired and match is of decent mapped length'''
            print(f'{ref},{pos},{tlen},{sampleid}')
        elif (tlen == 0) and (self.getmatchsize(fields[5]) >= self.min_match_length) and (self.get_gene_orgid(ref) == self.get_gene_orgid(ref2)):
            '''Improperly paired but same gene and match is of decent mapped length'''
            # RM < TODO Check logic of this with @tgolubch, seems to introduce quite a bit of error.
            print(f'{ref},{pos},{tlen},{sampleid}')
        else:
            '''No match'''
            return

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
        sampleid = self.argies.SeqName
        if self.argies.Mode == "filter":
            try:
                reads_to_drop = pd.read_csv(
                    self.argies.ExpDir + self.argies.FilterFile)
            except:
                raise SystemExit(
                    f"Couldn't find reads to drop file: {self.argies.FilterFile}. Did your run generate one?")

        for l in sys.stdin:
            if self.argies.Mode == "parse":
                self.parse_bam_position(l, sampleid)
            else:
                self.filter_bam(l, reads_to_drop)


if __name__ == '__main__':
    cls = Parse_bam_positions(parse_args_bam_parse())
    cls.main()
