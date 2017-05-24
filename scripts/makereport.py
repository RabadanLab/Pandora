#!/usr/bin/env python

import argparse
import sys
import os
import collections

# This script generates the report

# -------------------------------------

def get_arg():
    """
    Get Arguments
    :rtype: object
    """

    # parse arguments
    prog_description = 'Make tsv report'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('-i', '--input', default='blast/top.concat.txt', help='the contig blast file input')
    parser.add_argument('-i2', '--input2', default='discovery/blast/top.concat.txt', help='the ORF blast file input')
    parser.add_argument('-o', '--outputdir', default='report', help='the output directory')
    parser.add_argument('-l', '--logsdir', help='the logs directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--contigreport', default='report.contig.txt', help='name of the contig report (default: report.contig.txt)')
    parser.add_argument('--taxonreport', default='report.taxon.txt', help='name of the taxon report (default: report.taxon.txt)')
    parser.add_argument('--header', default='blast/header', help='the blast header file')
    parser.add_argument('--id2reads', default='assembly/reads2contigs.stats.txt', help='the output file of samtools idxstats mapping ids to #reads')
    parser.add_argument('--taxid2names', help='location of names.dmp file mapping taxid to names')
    parser.add_argument('--blacklist', help='the file of blacklist taxids')
    parser.add_argument('--noclean', type=int, default=0, help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--verbose', type=int, default=0, help='verbose mode: echo commands, etc (default: off)')
    parser.add_argument('--id', required=True, help='sample identifier')
    args = parser.parse_args()

    # add key-value pairs to the args dict
    vars(args)['step'] = 'reporting'

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    global makeHTML
    from helpers import helpers as hp
    from helpers import makeHTML

    # error checking: exit if previous step produced zero output
    for i in [args.input, args.header, args.id2reads]:
        hp.check_file_exists_and_nonzero(i, step=args.step)

    return args

# -------------------------------------

def get_reads_mapped_to_host(args):
    """
    Get the total number of reads properly mapped to host
    :rtype: int (number of reads)
    """

    num_reads_to_host = 0

    for myfile in ['host_separation/mapping_stats.STAR.txt', 'host_separation/mapping_stats.bwt.txt']:
        try:
            with open(myfile, 'r') as f:
                contents = f.read()
            # trying to parse lines that look like this: "115056674 + 0 mapped (80.95%:-nan%)"
            # assume that the first occurrence of the word 'mapped' is the line we want
            num_reads_to_host += int([i for i in contents.split('\n') if 'mapped' in i][0].split()[0])
        except:
            print('[Error] Can\'t find flagstat file')

    return num_reads_to_host

def makerep(args):
    """
    Make report
    """

    # The goal of this function is to filter blast results based on taxid
    # Don't want human sequences or anything in the taxid blacklist,
    # which contains non-pathogens

    hp.echostep(args.step)

    # print args
    print(args)
    print

    # mkdir -p
    hp.mkdirp(args.outputdir)

    # defaults
    # taxid corresponding to homo sapiens
    humantaxid = '9606'
    # desired header
    desiredfields = ['qseqid', 'sseqid', 'qlen','saccver','staxids','evalue', 'bitscore', 'stitle']

    # taxid blacklist
    filterlist = []
    # header
    header = []
    # id 2 #reads dict (output of samtools idxstats)
    idx = {}
    # a dictionary to map taxon to sum of uniq reads, longest contig id, longest contig length, number of contigs
    taxonstats = collections.defaultdict(dict)

    # load blacklist if supplied
    if args.blacklist:
        with open(args.blacklist, 'r') as f:
            # final newline causes empty list elt
            filterlist = f.read().split('\n')[:-1]

    # load idx file
    with open(args.id2reads, 'r') as f:
        for line in f:
            # map id to #reads
            idx[line.split()[0].strip()] = line.split()[2].strip()

    # load header
    with open(args.header, 'r') as f:
        header = map(str.strip, f.read().split())

    # get indicies,fields in file
    myindicies = []
    myfields = []
    for j,k in enumerate(header):
        if k in desiredfields:
            myindicies.append(j)
            myfields.append(k)

    # get index of taxid, qseqid:
    taxidindex = myfields.index('staxids')
    qseqidindex = myfields.index('qseqid')
    qlenindex = myfields.index('qlen')

    # list of input files
    myfiles = [args.input]
    # if there's a valid file produced by the ORF discovery step, examine that file also
    if hp.check_path_bool(args.input2):
        myfiles.append(args.input2)

    # write file keyed on contig
    with open(args.outputdir + '/blast.topfilter.unsort.txt', 'w') as f:
        # print header:
        f.write('sampleid\t' + '\t'.join(myfields) + '\tnum_reads\n')
        # loop through inputs (blast, ORF discovery)
        for myfile in myfiles:
            with open(myfile, 'r') as g:
                for line in g:
                    # don't want lines with predicted (as opposed to real) genes
                    if 'PREDICTED' in line:
                        continue
                    # get desired fields
                    fields = [line.split('\t')[i].strip() for i in myindicies]
                    taxid = fields[taxidindex]
                    qseqid = fields[qseqidindex]
                    qlen = fields[qlenindex]
                    # bypass human taxids or specifically filtered taxids
                    if taxid == humantaxid or taxid in filterlist:
                        continue
                    # multiple taxids not supported
                    if ';' in taxid:
                        print('[WARNING] semicolon detected: multiple taxids not supported')
                    # get read counts
                    readcounts = idx.get(qseqid, '-')
                    # write to file
                    f.write(args.id + '\t' + '\t'.join(fields) + '\t' + readcounts + '\n')
                    # set taxonstats dictionary
                    taxonstats[taxid]['num'] = taxonstats[taxid].get('num', 0) + 1
                    taxonstats[taxid]['longest'] = qseqid
                    taxonstats[taxid]['longestlength'] = max(taxonstats[taxid].get('longestlength', 0), int(qlen))
                    if not readcounts == '-':
                        taxonstats[taxid]['sum'] = taxonstats[taxid].get('sum', 0) + int(readcounts)

    # sort by staxids then qlen (with bash)
    # careful: you're including the header in the file (make sure it's sorted properly)
    cmd = 'sort -k5,5n -k6,6nr {args.outputdir}/blast.topfilter.unsort.txt > {args.outputdir}/{args.contigreport}'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    # get num reads to host
    num_reads_to_host = get_reads_mapped_to_host(args)
    if args.verbose:
        print('Number of reads properly mapped to host: {reads}'.format(reads=num_reads_to_host))
#       print(dict(taxonstats))

    # write file keyed on taxon
    with open(args.outputdir + '/' + args.taxonreport, 'w') as f:
        # print header:
        f.write('sampleid\ttaxid\tsum_reads\tnum_contigs\tid_longest_contig\tlen_longest_contig\tpathogen_reads/(host_reads/10^6)\n')
        for taxid in taxonstats:
            # pathogen reads / host reads
	    sumdivhost = '-'
            # if num_reads_to_host nonzero
            if num_reads_to_host:
                sumdivhost = str(round(1000000*taxonstats[taxid]['sum']/float(num_reads_to_host),2))
            mytaxonattributes = [str(taxonstats[taxid]['sum']),
                                 str(taxonstats[taxid]['num']),
                                 taxonstats[taxid]['longest'],
                                 str(taxonstats[taxid]['longestlength']),
                                 sumdivhost]
            f.write(args.id + '\t' + taxid + '\t' + '\t'.join(mytaxonattributes) + '\n')
    
    # generate and write html report
    if args.taxid2names != 'None':
        makeHTML.generateHTML(
            args.outputdir + '/' + args.taxonreport,
            args.scripts,
            args.taxid2names,
            args.outputdir)
    else: 
        print('[WARNING] missing names.dmp, the file mapping taxids to names. HTML report will not be generated.')
    
#    if args.verbose:
#        print(dict(taxonstats))

    if not args.noclean:
        os.remove(args.outputdir + '/blast.topfilter.unsort.txt')

    hp.echostep(args.step, start=0)

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # make report
    makerep(args)

# -------------------------------------

if __name__ == '__main__':

    main()
