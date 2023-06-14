#!/usr/bin/env python

'''
Example:

0 1 K00150:209:HJT2MBBXX:1:1224:26961:10370
1 2 1105
2 3 cytomegallovirus_allrecords_cluster_1
3 4 12002
4 5 10
5 6 17S26M107S
6 7 cardiovirus_allrecords_cluster_5
7 8 227
8 9 0
'''
from __future__ import print_function
import sys
import re

def getmatchsize(cigar):
    matches = re.findall(r'([0-9]+)M', cigar)
    if not len(matches): 
        return 0
    return sum(int(x) for x in matches)

def get_gene_orgid(target_id):
    parts = target_id.split('_')
    gene = parts[0]
    orgid = parts[-1] if gene.startswith('BACT') else gene
    return (gene, orgid)


if __name__=='__main__':
    if len(sys.argv) < 2 or '-h' in sys.argv:
        sys.stderr.write('Usage: samtools view MyBamFile | {} \n\n'.format(sys.argv[0]))
        sys.exit(1)
    min_match_length = 40
    sampleid = sys.argv[1]
    for l in sys.stdin:
        fields = l.split()
        ref = fields[2]
        pos = fields[3]
        ref2 = fields[6]
        tlen = int(fields[8])
        if tlen >= min_match_length:
            # properly paired and match is of decent mapped length
            print(f'{ref},{pos},{tlen},{sampleid}')
        elif (tlen == 0) and (getmatchsize(fields[5]) >= min_match_length) and (get_gene_orgid(ref) == get_gene_orgid(ref2)): 
            # improperly paired but same gene and match is of decent mapped length 
            print(f'{ref},{pos},{tlen},{sampleid}')
        else:
            continue
