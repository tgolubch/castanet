import argparse, os, sys
from app.utils.argparsemappings import map_args_filter_keep_reads

def parse_args_filter_keep_reads():
    '''
    Parse command-line arguments, if script called from command line. Not used for API entry.
    '''
    parser = argparse.ArgumentParser(description=__doc__, epilog='Example: {progname} -i /path/to/my/fastq.gz -x 9606'.format(progname=os.path.basename(sys.argv[0])))
    parser.add_argument( '-i', required=True, nargs='*', help='Fastq[.gz] files to filter. Files must NOT contain blank lines. If gzipped, filenames must end in .gz. Each file will produce output with the suffix _filt.fastq in the current working directory.' )
    parser.add_argument( '-k', default=None, required=True, help='Path to corresponding Kraken[.gz] output file, or other whitespace-separated text file where the second column is the read name and the third is the TaxID.' )
    parser.add_argument( '-r', default='', help='NCBI TaxID(s) to retain. Use comma-separated list without spaces for multiple TaxIDs). Default is none.' )
    parser.add_argument( '-x', default='', help='NCBI TaxID(s) to exclude. Use comma-separated list without spaces for multiple TaxIDs). Default is none. To remove all reads marked as human set to 9606.' )
    parser.add_argument( '--rT', default='', help='Text names of top-level taxa, to retain this and all taxa below it. Requires "--lineagefile". Eg. "--rT Bacteria,Viruses" to retain all bacterial and viral sequences.' )
    parser.add_argument( '--xT', default='', help='Text names of top-level taxa, to exclude this and all taxa below it. Requires "--lineagefile". Eg. "--xT \"Homo sapiens,Fungi\"" to exclude all human and fungal sequences.' )
    parser.add_argument( '--lineagefile', default='data/ncbi_lineages_2023-06-15.csv.gz', help='Path to CSV file containing lineages of all NCBI taxa. Default is "lineages-2018-03-12.csv.gz".' )
    return map_args_filter_keep_reads(parser.parse_args())