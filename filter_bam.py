#!/usr/bin/env python

'''
Example of BAM file fields:

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

import sys
import pandas as pd

if __name__=='__main__':
    if len(sys.argv) < 3 or '-h' in sys.argv:
        sys.stderr.write('Usage: samtools view -h MyBamFile | {} SAMPLENAME READS_TO_DROP.csv\n\n'.format(sys.argv[0])) 
	sys.exit(1)
    sampleid = sys.argv[1]
    # read list of reads to be dropped, header is "sampleid,target_id,startpos,maplen"
    reads_to_drop = pd.read_csv(sys.argv[2])
    reads_to_drop = reads_to_drop[reads_to_drop.sampleid==sampleid].set_index(['target_id','startpos','maplen'])
    
    for l in sys.stdin:
        # output header as is
        if l.startswith('@'):
            sys.stdout.write(l)
            continue
        if len(reads_to_drop):
            fields = l.split()
            ref = fields[2]
            pos = int(fields[3])
            #ref2 = fields[6]
            tlen = int(fields[8])
            # output only reads that are not found in the indexed READS_TO_DROP file
            if not (ref, pos, tlen) in reads_to_drop.index:
                sys.stdout.write(l)
            else:
                continue
