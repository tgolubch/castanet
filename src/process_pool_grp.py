#!/usr/bin/env python

# ============================================================================ #
# Copyright (c) Tanya Golubchik                                                #
# golubchi@bdi.ox.ac.uk                                                        #
# August 2018                                                                  #
# ============================================================================ #

# ============================================================================ #
# Import Modules                                                               #
# ============================================================================ #

from __future__ import division
import os, sys, argparse, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
try:
    # Use seaborn for prettier plots if available
    import seaborn as sns
except ImportError:
    pass

# ============================================================================ #
# GLOBAL VARIABLES                                                             #
# ============================================================================ #
_args = None
_df = None
_progname=os.path.basename(sys.argv[0])

# ============================================================================ #
# LOGGING                                                                      #
# ============================================================================ #
def loginfo(s):
    sys.stderr.write('  Info: {0}\n'.format(s))
def logerr(s):
    sys.stderr.write('  Warning: {0}\n'.format(s))
def stoperr(s, errcode=1):
    errword = 'Finished' if not errcode else 'Error'
    sys.stderr.write('  {0}: {1}\n'.format(errword, s))
    sys.exit(errcode)

# ============================================================================ #
# PROGRAM USAGE                                                                #
# ============================================================================ #
def Initialise():
    '''
    Parse command-line arguments.
    '''
    global _args
    parser = argparse.ArgumentParser(description=__doc__, epilog='Example: {progname} -i /path/to/my/dataframe.csv.gz --samples /path/to/sampleinfo.csv'.format(progname=_progname))
    parser.add_argument( '-b', '--batchname', required=True, help='Batch name for these samples. Must be alphanumeric.' )
    parser.add_argument( '-i', '--infile', required=True, help='Data frame (csv[.gz]) file to process. If gzipped, filename must end in .gz.' )
    parser.add_argument( '-o', '--outdir', default=os.getcwd(), help='Output directory. Default: current working directory.' )
    parser.add_argument( '-p', '--probelengths', default='{}/probelengths_rmlst_virus_extra_ercc.csv'.format(os.path.dirname(sys.argv[0])),\
                                help='Path to file containing probe lengths. Default: {}/probelengths_rmlst_virus_extra_ercc.csv'.format(os.path.dirname(sys.argv[0])) )
    parser.add_argument( '-d', '--keepdups', action='store_true', default=True, help='If true, do not reassign duplicates to the sample with the majority in each duplicate cluster (Default: False).' )
    parser.add_argument( '--samples', required=True, help='Path to CSV file containing information about raw reads (must have at least following fields: sampleid, pt, rawreadnum). Field "pt" must match clinical data.' )
    parser.add_argument( '--clin', required=False, help='Path to CSV file containing clinical data (must have at least following fields: pt, clin_int; the field "sampleid" if present will be ignored). Other fields will be ignored.' )
    parser.add_argument( '--depth_inf', required=False, \
                                help='(For regenerating full CSV with new clinical info): Path to previously generated CSV file of read depth per position for each probe, for all samples in this batch. Must contain the following fields: sampleid, target_id, depth_mean, depth_std, depth_25pc, depth_median, depth_75pc, prop_target_covered, prop_target_covered_mindepth2, prop_target_covered_mindepth5, prop_target_covered_mindepth10, udepth_mean, udepth_std, udepth_25pc, udepth_median, udepth_75pc, uprop_target_covered, uprop_target_covered_mindepth2, uprop_target_covered_mindepth5, uprop_target_covered_mindepth10' )
    _args = parser.parse_args()
    return


def Clean_Commandline():
    '''
    Print errors, and exit if necessary, on bad input data.
    '''
    # Validate main infput file (dataframe from processed BAM files for this pool)
    if not os.path.isfile(_args.infile):
        stoperr('Unable to open input file {0}.'.format(_args.infile))

    # Validate output directory
    if not os.path.isdir(_args.outdir):
        try:
            os.mkdir(_args.outdir)
        except OSError:
            stoperr('Failed to create output directory {}.'.format(_args.outdir))

    # Validate sample info file
    if _args.samples and not os.path.isfile(_args.samples):
        stoperr('Unable to open sample data from input file {0}.'.format(_args.samples))
    with open(_args.samples) as samples_inf:
        samples_header_check = samples_inf.readline()
        if ('sampleid' not in samples_header_check) or ('pt' not in samples_header_check) or ('rawreadnum' not in samples_header_check):
            stoperr('{} must contain at least the following columns: pt, sampleid, rawreadnum'.format(_args.samples))
    if _args.clin:
        if not os.path.isfile(_args.clin):
            stoperr('Unable to open clinical data from input file {0}.'.format(_args.clin))
        # Validate clinical info file
        with open(_args.clin) as clin_inf:
            clin_header_check = clin_inf.readline().split(',')
            if ('pt' not in clin_header_check):
                stoperr('{} must contain at least the following columns: pt, clin_int'.format(_args.clin))

    # Validate batch name
    if not _args.batchname.replace('_', '').isalnum():
        stoperr('{} is not a valid batch name. Must be alphanumeric.'.format(_args.batchname))

    # Open data frame
    try:
        df = pd.read_csv(_args.infile, compression=('gzip' if _args.infile.endswith('.gz') else None), header=None, \
                                names=['n', 'target_id', 'startpos', 'maplen', 'sampleid'])
    except (IOError, TypeError, ParserError) as e:
        stoperr('Failed to read dataframe from {} : {}'.format(_args.df, e))
    if _args.depth_inf and not os.path.isfile(_args.depth_inf):
        stoperr('Unable to open precomputed depth file {0}.'.format(_args.depth_inf))
    return df


# ============================================================================ #
# DATA PROCESSING FUNCTIONS                                                    #
# ============================================================================ #

def Add_Probetype(df):
    ''' For probe aggregation, determine organism and gene ID.'''
    loginfo('Aggregating by organism and gene name.')
    df['genename'] = df.target_id.apply(lambda x: x.replace('_','-').split('-')[0])
    df.loc[df.target_id=='roseolovirus_allrecords_cluster_1', 'genename'] = 'HHV7_roseolovirus_allrecords_cluster_1'
    df.loc[df.target_id=='roseolovirus_allrecords_cluster_2', 'genename'] = 'HHV6_roseolovirus_allrecords_cluster_2'
    # more precise definition for the different virus types
    df.loc[df.genename=='enterovirus', 'genename'] = df.loc[df.genename=='enterovirus'].target_id.apply(lambda x: '_'.join(x.replace('_','-').split('-')[:2]))
    df.loc[df.genename=='Coronaviridae', 'genename'] = df.loc[df.genename=='Coronaviridae'].target_id.apply(lambda x: '_'.join(x.replace('_','-').split('-')[:3]))
    df.loc[df.genename=='Adenoviridae', 'genename'] = df.loc[df.genename=='Adenoviridae'].target_id.apply(lambda x: '_'.join(x.replace('_','-').split('-')[:2]))
    df.loc[df.genename=='Flaviviridae', 'genename'] = df.loc[df.genename=='Flaviviridae'].target_id.apply(lambda x: '_'.join(x.replace('_','-').split('-')[:2]))
    df.loc[df.genename=='Influenza', 'genename'] = df.loc[df.genename=='Influenza'].target_id.apply(lambda x: '_'.join(x.replace('_','-').split('-')[:2]))
    df.loc[df.genename=='Paramyxoviridae', 'genename'] = df.loc[df.genename=='Paramyxoviridae'].target_id.apply(lambda x: '_'.join(x.replace('_','-').split('-')[:2]))
    df.loc[df.genename=='Parvoviridae', 'genename'] = df.loc[df.genename=='Parvoviridae'].target_id.apply(lambda x: '_'.join(x.replace('_','-').split('-')[:2]))
    df['probetype'] = df.genename
    pat = re.compile('BACT[0-9]+_([A-Za-z]+)-[0-9]+[|_]([A-Za-z]+)')
    backup_pat = re.compile('BACT[0-9]+_[0-9]+_[A-Z]+_([A-Za-z_]+)')
    def _pat_search(s, pat=pat, backup_pat=backup_pat):
        ''' Private function to return empty string instead of error when pattern is not matched.'''
        res = pat.findall(s)
        if not res: res = (backup_pat.findall(s),)
        if not res: return ''
        return '_'.join(res[0])
    df.loc[df.genename.str.startswith('BACT'), 'probetype'] = df.loc[df.genename.str.startswith('BACT')].target_id.apply(_pat_search)
    # Streptococcus mitis group (pneumo, mitis, oralis) are cross-mapping,
    # so err on the side of pneumo and classify as S. pneumoniae all targets
    # that are found in S.pneumo at least once in the database
    df.loc[(df.target_id.apply(lambda x: 'pneumoniae' in x)) & (df.probetype.isin(('Streptococcus_pneumoniae', 'Streptococcus_pseudopneumoniae', 'Streptococcus_mitis', 'Streptococcus_oralis'))), 'probetype'] = 'Streptococcus_pneumoniae'
    # Streptococcus_agalactiae and Streptococcus_pyogenes are cross-mapping as well, aggregate them
    df.loc[(df.target_id.apply(lambda x: 'pyogenes' in x or 'agalactiae' in x)) & (df.probetype.isin(('Streptococcus_pyogenes', 'Streptococcus_agalactiae'))), 'probetype'] = 'Streptococcus_agalactiae_pyogenes'
    # Enterobacteriacae are not distinguishable at this level so group them all
    df.loc[df.probetype.apply(lambda x: x.startswith('Escherichia') or x.startswith('Klebsiella') or x.startswith('Enterobacter')), 'probetype'] = 'Enterobacteriaceae'
    loginfo('Organism and gene summary: {} organisms, up to {} genes each.'.format(\
            df.probetype.nunique(), df.groupby('probetype').genename.nunique().max()))
    # loginfo(probelengths.groupby('probetype').genename.nunique().sort_values())
    return df

def Reassign_Dups(df):
    ''' Reassign positional duplicates (reads with same start/end) to the sample with the highest count.'''
    loginfo('Reassaigning duplicates.')
    ttl_prededup_mapped_reads = df.groupby('sampleid').n.sum()
    grouped_prededup_mapped_reads = df.groupby(['sampleid','probetype']).n.sum().unstack().fillna(0)
    df.sort_values(['target_id', 'startpos', 'maplen', 'n'], ascending=True, inplace=True)
    df['n'] = df.groupby(['target_id', 'startpos', 'maplen']).n.cumsum()
    nrows0 = len(df)
    duprows = df.duplicated(['target_id', 'startpos', 'maplen'], keep='last')
    loginfo('Saving {0}_reads_to_drop.csv'.format(_args.batchname))
    df[duprows][['sampleid', 'target_id', 'startpos', 'maplen']].to_csv('{}/{}_reads_to_drop.csv'.format(_args.outdir, _args.batchname), index=False)
    df = df[~duprows]
    # df.drop_duplicates(['target_id', 'startpos', 'maplen'], keep='last', inplace=True)
    loginfo('Kept {} of {} rows ({:.2f}).'.format(len(df), nrows0, float(len(df))/nrows0))
    ttl_dup_mapped_reads = df.groupby('sampleid').n.sum()
    grouped_dup_mapped_reads = df.groupby(['sampleid','probetype']).n.sum().unstack().fillna(0)
    try:
        assert ttl_dup_mapped_reads.sum() == ttl_prededup_mapped_reads.sum()
    except AssertionError:
        stoperr('Total reads after duplicate reassignment does not match starting total.')
    duprate_ttl = ttl_dup_mapped_reads/ttl_prededup_mapped_reads
    duprate_by_probetype = grouped_dup_mapped_reads/grouped_prededup_mapped_reads
    loginfo('Saving {0}_duprate.csv and {0}_duprate_by_probetype.csv.'.format(_args.batchname))
    duprate_ttl.to_csv('{}/{}_duprate.csv'.format(_args.outdir, _args.batchname, header=True))
    duprate_by_probetype.to_csv('{}/{}_duprate_by_probetype.csv'.format(_args.outdir, _args.batchname))
    loginfo('Duplication rate:')
    loginfo(duprate_ttl.describe().to_string())
    return df

def Add_Probelength(df):
    ''' Add length of target_id to each row.'''
    loginfo('Adding probe length information from {}'.format(_args.probelengths))
    try:
        probelengths = pd.read_csv(_args.probelengths, dtype={'target_len':int})
    except:
        stoperr('Failed to read probe information. Is {} a valid CSV file?'.format(_args.probelengths))
    probelengths = Add_Probetype(probelengths)
    df = df.merge(probelengths, left_on='target_id', right_on='target_id', how='left')
    return df, probelengths

def Add_Depth(df, probelengths):
    ''' Calculate read depth per position. '''
    if _args.depth_inf:
        loginfo('Reading read depth information from {}.'.format(_args.depth_inf))
        try:
            depth = pd.read_csv(_args.depth_inf)
            assert len(depth) > 0
        except (IOError, AssertionError, TypeError):
            logerr('Failed to read read depth information from file. Will calculate from data.')
            depth = None
    else:
        loginfo('Calculating read depth information.')
        depth = None
    if depth is None:
        metrics = {}
        odir = '{}/Depth_{}'.format(_args.outdir, _args.batchname)
        if not os.path.isdir(odir):
            try:
                os.mkdir(odir)
            except OSError:
                odir=os.getcwd()
                logerr('Cannot create output directory {} for saving depth information. Proceeding with current working directory.'.format(odir))
        loginfo('Calculating read depth statistics for all probes, for all samples. This is a slow step.')
        for (sampleid, probetype), g in df.groupby(['sampleid', 'probetype']):
            gene_list = g.genename.unique()
            n_genes = len(gene_list)
            target_list = g.target_id.unique()
            n_targets = len(target_list)
            Dd = {}; D1d = {}
            for genename, gg in g.groupby('genename'):
                loginfo('Processing {} - {}'.format(sampleid, genename))
                # Generate two arrays for each probetype-genename group (D = number of occurrences, D = unique start/end positions)
                D = np.zeros(int(gg.target_len.max()), dtype=np.uint32)
                D1 = np.zeros(D.shape, dtype=np.uint32)
                D.fill(0); D1.fill(0) # Clear target arrays
                for target_id, gt in gg.groupby('target_id'):
                    loginfo('..... target: {}'.format(target_id[:100]))
                    for rowname, row in gt.iterrows():
                        D[row.startpos-1:row.startpos-1+row.maplen] += row.n
                        D1[row.startpos-1:row.startpos-1+row.maplen] += 1
                Dd[genename] = D
                D1d[genename] = D1
            #
            # Max possible number of targets for this probetype, aggregating all genes
            nmax_targets = probelengths[probelengths.probetype==probetype].target_id.nunique()
            # Max possible number of genes for this probetype, aggregating all genes
            nmax_genes = probelengths[probelengths.probetype==probetype].genename.nunique()
            # Max possible positions for this probetype, aggregating all genes
            nmax_probetype = probelengths[probelengths.probetype==probetype].groupby('genename').target_len.max().sum()
            #
            # Collapse to a single array for the entire genename group for this probetype in this sample
            D = np.hstack(Dd.values())
            D1 = np.hstack(D1d.values())
            amprate = (D/D1)[~np.isnan(D/D1)]
            # Max possible positions for the genes that were actually in this BAM (accounts for some genes not being captured)
            npos = len(D)
            # Now pad out with zeros to the total number of mappable positions for this probetype (nmax_probetype above)
            D = np.pad(D, (0, nmax_probetype - npos), 'constant', constant_values=0 )
            D1 = np.pad(D1, (0, nmax_probetype - npos), 'constant', constant_values=0 )
            loginfo('Mean depth (all reads) for {}: {:.2f}'.format(probetype, D.mean()))
            loginfo('Mean depth (deduplicated) for {}: {:.2f}'.format(probetype, D1.mean()))
            # Amplification rate calculations use the unpadded (mapped) number of sites as the denominator
            loginfo('Mean amplification ratio for {}: {:.2f}'.format(probetype, amprate.mean()))
            #
            # Save arrays as CSV
            if not os.path.isdir('{}/{}'.format(odir, sampleid)):
                os.mkdir('{}/{}'.format(odir, sampleid))
            with open('{}/{}/{}-{}_depth_by_pos.csv'.format(odir, sampleid, probetype, sampleid), 'a') as o:
                np.savetxt(o, D, fmt='%d', newline=',')
                o.write('\n')
                np.savetxt(o, D1, fmt='%d', newline=',')
                o.write('\n')
            #
            # Save array plots as PDF
            if D1.mean() >= 0.1:
                plt.clf()
                plt.plot(D, label='all reads')
                plt.plot(D1, label='deduplicated reads', c='orange')
                plt.xlabel('position')
                plt.ylabel('number of reads')
                plt.legend()
                plt.title('{}\n{} ({}/{} targets in {}/{} genes)\n'.format(sampleid, probetype, n_targets, nmax_targets, n_genes, nmax_genes))
                try:
                    plt.tight_layout()
                except ValueError:
                    pass
                plt.savefig('{}/{}/{}-{}.pdf'.format(odir, sampleid, probetype, sampleid))
            #
            # Build up dictionary of depth metrics for this sample and probetype
            metrics[sampleid, probetype] = (g.n.sum(), g.n.count(), n_targets, n_genes, nmax_targets, nmax_genes, nmax_probetype, npos,\
                                        amprate.mean(), amprate.std(), np.median(amprate),\
                                        D.mean(), D.std(), np.percentile(D, 25), np.median(D), np.percentile(D, 75),\
                                        (D>0).sum(),\
                                        (D>=2).sum(),\
                                        (D>=5).sum(),\
                                        (D>=10).sum(),\
                                        (D>=100).sum(),\
                                        (D>=1000).sum(),\
                                        D1.mean(), D1.std(), np.percentile(D1, 25), np.median(D), np.percentile(D1, 75),\
                                        (D1>0).sum(),\
                                        (D1>=2).sum(),\
                                        (D1>=5).sum(),\
                                        (D1>=10).sum(),\
                                        (D1>=100).sum(),\
                                        (D1>=1000).sum() )
        #
        # Data frame of all depth metrics
        depth = pd.DataFrame(metrics, index=['n_reads_all', 'n_reads_dedup', 'n_targets', 'n_genes', 'nmax_targets', 'nmax_genes', 'npos_max_probetype', 'npos_cov_probetype',\
                                             'amprate_mean', 'amprate_std', 'amprate_median',\
                                             'depth_mean', 'depth_std', 'depth_25pc', 'depth_median', 'depth_75pc',\
                                             'npos_cov_mindepth1', \
                                             'npos_cov_mindepth2',\
                                             'npos_cov_mindepth5', \
                                             'npos_cov_mindepth10',\
                                             'npos_cov_mindepth100',\
                                             'npos_cov_mindepth1000',\
                                             'udepth_mean', 'udepth_std', 'udepth_25pc', 'udepth_median', 'udepth_75pc',\
                                             'npos_dedup_cov_mindepth1', \
                                             'npos_dedup_cov_mindepth2',\
                                             'npos_dedup_cov_mindepth5', \
                                             'npos_dedup_cov_mindepth10', \
                                             'npos_dedup_cov_mindepth100', \
                                             'npos_dedup_cov_mindepth1000']).T.reset_index()
        depth.rename(columns={'level_0':'sampleid', 'level_1':'probetype'}, inplace=True)
        # Add reads on target
        rot = depth.groupby('sampleid').n_reads_all.sum().reset_index()
        rot.rename(columns={'n_reads_all':'reads_on_target'}, inplace=True)
        depth = depth.merge(rot, on='sampleid', how='left')
        rot_dedup = depth.groupby('sampleid').n_reads_dedup.sum().reset_index()
        rot_dedup.rename(columns={'n_reads_dedup':'reads_on_target_dedup'}, inplace=True)
        depth = depth.merge(rot_dedup, on='sampleid', how='left')
        depth['prop_of_reads_on_target'] = depth.n_reads_all/depth.reads_on_target
        depth['prop_npos_cov1'] = depth.npos_cov_mindepth1/depth.npos_max_probetype
        depth['prop_npos_cov2'] = depth.npos_cov_mindepth2/depth.npos_max_probetype
        depth['prop_npos_cov5'] = depth.npos_cov_mindepth5/depth.npos_max_probetype
        depth['prop_npos_cov10'] = depth.npos_cov_mindepth10/depth.npos_max_probetype
        depth['prop_npos_cov100'] = depth.npos_cov_mindepth100/depth.npos_max_probetype
        depth['prop_npos_cov1000'] = depth.npos_cov_mindepth1000/depth.npos_max_probetype
        depth['prop_ntargets'] = depth.n_targets/depth.nmax_targets
        depth['prop_ngenes'] = depth.n_genes/depth.nmax_genes
        # Add log transforms
        depth['log10_depthmean'] = depth.depth_mean.apply(lambda x: np.log10(x+1))
        depth['log10_udepthmean'] = depth.udepth_mean.apply(lambda x: np.log10(x+1))
        # Duplicated reads only
        depth['clean_n_reads_all'] = depth.n_reads_all-depth.n_reads_dedup
        depth['clean_prop_of_reads_on_target']=(depth.n_reads_all-depth.n_reads_dedup)/depth.reads_on_target
        #
        loginfo('Saving {}/{}_depth.csv.'.format(_args.outdir, _args.batchname))
        depth.to_csv('{}/{}_depth.csv'.format(_args.outdir, _args.batchname), index=False)
        #
        loginfo('Mean read depth per sample:')
        loginfo(depth.groupby('sampleid').depth_mean.mean().to_string())
    return depth


def Add_Clin(df, req_cols_samples = ['sampleid', 'pt', 'rawreadnum'], req_cols_clin = ['pt', 'clin_int']):
    ''' Add raw read numbers and any external categorical/clinical data.
    Samples file must supply at least the following columns: {}.
    Clin file may additionally supply clinical/demographic data,
    with at least the following columns: {}'''.format(req_cols_samples, req_cols_clin)
    # Raw read numbers
    loginfo('Adding sample information and clinical data.')
    samples = pd.read_csv(_args.samples, dtype={'pt':str})
    if set(req_cols_samples) & set(samples.columns) != set(req_cols_samples):
        stoperr('Samples file {} must contain at least the following columns: {}'.format(_args.samples, req_cols_samples))
    if _args.clin:
        # Merge to create simpler metadata for each sample, including patient ID and clinical category
        clin = pd.read_csv(_args.clin, dtype={'pt':str})
        if 'pt' not in clin.columns:
            stoperr('Clin file {} must contain at least the following columns: {}'.format(_args.clin, req_cols_clin))
        if 'sampleid' in clin.columns:
            logerr('Ignoring "sampleid" column from clinical info file {}.'.format(_args.clin))
            clin.drop('sampleid', axis=1, inplace=True)
            # Merge sample information with participant data
            samples = samples.merge(clin, on='pt', how='left')
    df = df.merge(samples, on='sampleid', how='left')
    df['readprop'] = df.n_reads_all/df.rawreadnum
    loginfo('Added the following columns: {}'.format(list(samples.columns)))
    loginfo('Saving {}/{}_depth_with_clin.csv.'.format(_args.outdir, _args.batchname))
    df.to_csv('{}/{}_depth_with_clin.csv'.format(_args.outdir, _args.batchname), index=False)
    return df


def Save_Tophits(depth):
    ''' Simple output of likely best hit for each sample.'''
    loginfo('Target with highest proportion of captured reads:')
    tophits = depth.sort_values(['sampleid', 'clean_prop_of_reads_on_target'], ascending=True).drop_duplicates('sampleid', keep='last')
    loginfo(tophits[['sampleid', 'probetype', 'clean_prop_of_reads_on_target']])
    tophits.to_csv('{}/{}_tophit.csv'.format(_args.outdir, _args.batchname), index=False)
    loginfo('Saved top hits to {}/{}_tophits.csv'.format(_args.outdir, _args.batchname))
    return


# ============================================================================ # 
# MAIN                                                                         #
# ============================================================================ #

if __name__ == '__main__':
    Initialise()
    df = Clean_Commandline()
    df, probelengths = Add_Probelength(df)
    if not _args.keepdups:
        df = Reassign_Dups(df)
    # Depth calculation
    depth = Add_Depth(df, probelengths)
    # Merge in sample info  (including total raw reads) and participant data if specified
    depth = Add_Clin(depth)
    Save_Tophits(depth)
    df = df.merge(depth, on=['sampleid', 'probetype'], how='left')
    df.to_csv('{}/{}_fullDF.csv.gz'.format(_args.outdir, _args.batchname), index=False, compression='gzip')
    loginfo('Finished. Saved final data frame as {}/{}_fullDF.csv.gz'.format(_args.outdir, _args.batchname))

# ============================================================================ #
