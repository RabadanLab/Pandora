#!/usr/bin/env python

import argparse
import sys
import subprocess
import os
import glob
import shutil
import fileinput

# This script blasts the entries of a fasta file,
# provided they're over the length threshold
# It also concatenates the results and does some minor filtering

humantaxid = '9606'

# -------------------------------------

def get_arg():
    """
    Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'Blast in parallel'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('-i', '--input', default='assembly/contigs_trinity.fasta', help='the input fasta')
    parser.add_argument('-o', '--outputdir', default='blast', help='the output directory')
    parser.add_argument('-l', '--logsdir', default='logs_blast', help='the logs directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--noclean', type=int, default=0, help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--verbose', type=int, default=0, help='verbose mode: echo commands, etc (default: off)')
    parser.add_argument('--threshold', type=int, default=0, help='the length threshold')
    parser.add_argument('--filelength', type=int, default=500, help='the number of rows per split file')
    parser.add_argument('--bmem', help='memory (in G) for qsub of individual blast array job task')
    parser.add_argument('--btime', help='time (in hours) for qsub of individual blast array job task')
    parser.add_argument('--db', help='the database prefix')
    parser.add_argument('--whichblast', default='blastn', choices=['blastn', 'blastp'], help='which blast to use (blastn, blastp)')
    parser.add_argument('--threads', default='1', help='blast -num_threads option')
    parser.add_argument('--nosge', type=int, default=0, help='no SGE bool')
    parser.add_argument('--hpc', type=int, default=0, help='run on the CUMC hpc cluster (add additional qsub flags)')
    parser.add_argument('--id', help='id')
    args = parser.parse_args()

    # blast format string
    fmt="qseqid sseqid saccver staxids pident nident length mismatch gapopen gaps qstart qend qlen qframe qcovs sstart send slen sframe sstrand evalue bitscore stitle"

    # add key-value pairs to the args dict
    vars(args)['fmt'] = fmt
    vars(args)['step'] = 'blast_contigs'

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # error checking: exit if previous step produced zero output
    for i in [args.input]:
        hp.check_file_exists_and_nonzero(i, step=args.step)

    return args

# -------------------------------------

def blast(args):
    """Do blast in parallel"""

    hp.echostep(args.step)

    # print args
    print(args)
    print

    # mkdir -p 
    hp.mkdirp(args.outputdir)

    args.noclean = int(args.noclean)

    # generate header
    with open(args.outputdir + '/header', 'w') as f:
        f.write(args.fmt.replace(' ', '\t') + '\n')

    # if no qsub
    if args.nosge:
        # filter fasta file on contigs above threshold length and hardcode name for blast.py
        # (splitting files doesnt make sense if no cluster)
        (numcontigs, filecount) = hp.fastasplit3(args.input, args.outputdir + '/blast', args.threshold)

        if numcontigs == 0:
            print("No contigs above threshold. Exiting")
            sys.exit(1)
        else:
            print("There are " + str(numcontigs) + " contigs above threshold, and " + str(filecount) + " files to blast.")


        # define log files
        logs_out = args.outputdir + '/' + 'blast_log_' + args.id + '.o'
        logs_err = args.outputdir + '/' + 'blast_log_' + args.id + '.e'
        # define command
        cmd = '{args.scripts}/scripts/blast.py --scripts {args.scripts} --outputdir {args.outputdir} --whichblast {args.whichblast} --db {args.db} --threads {args.threads} --fmt "{args.fmt}" --sgeid 1 > {o} 2> {e}'.format(args=args, o=logs_out, e=logs_err)
        hp.run_cmd(cmd, args.verbose, 0)


        # now filter blast results
        concat(args)

    # if qsub
    else:
        # mkdir -p
        hp.mkdirp(args.logsdir)

        # split fasta file on contigs above threshold length (and return number of contigs, file count)
        (numcontigs, filecount) = hp.fastasplit2(args.input, args.outputdir + '/blast', args.threshold, args.filelength)

        if filecount == 0:
            print("No contigs above threshold. Exiting")
            sys.exit(1)
        else:
            print("There are " + str(numcontigs) + " contigs above threshold, and " + str(filecount) + " files to blast.")

        # qsub part of command (array job)
        qcmd = 'qsub -S {mypython} -N bc_{args.id} -e {args.logsdir} -o {args.logsdir} -l mem={args.bmem}G,time={args.btime}:: -t 1-{filecount} '.format(mypython=sys.executable, args=args, filecount=filecount)
        #if args.hpc:
        #    qcmd += ...
        # regular part of command
        cmd = '{args.scripts}/scripts/blast.py --scripts {args.scripts} --outputdir {args.outputdir} --whichblast {args.whichblast} --db {args.db} --threads {args.threads} --fmt "{args.fmt}"'.format(args=args)
        if args.verbose:
            print(qcmd + cmd)
        message = subprocess.check_output(qcmd + cmd, shell=True)
        print(message)
        # get job id
        jid = hp.getjid(message)

        # hold the script up here, until all the blast jobs finish
        # concat top blast hits; concat log files into one, so as not to clutter the file system
        # qsub part of command
        qcmd = 'qsub -V -b y -cwd -o log.out -e log.err -N wait_{args.id} -hold_jid {jid} -sync y echo wait_here'.format(args=args, jid=jid)
        message = subprocess.check_output(qcmd, shell=True)
        print(message)

        # now concatenate and filter blast results
        concat(args)

    hp.echostep(args.step, start=0)

# -------------------------------------

def concat(args):
    """
    Concatenate blast files and logs, so as not to leave many files messily scattered about.
    Also, implement Ioan's filtering
    """

    print('CONCATENATE START')

    # define commands
    # file of all blast hits
    cmd = 'cat {args.outputdir}/*.result > {args.outputdir}/concat.txt'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)
    # all fasta entries
    cmd = 'cat {args.outputdir}/*.fasta > {args.outputdir}/above_threshold.fa'.format(args=args)
    hp.run_cmd(cmd, args.verbose, 0)

    # set of all input IDs from the concatenated fasta file
    with open(args.outputdir + '/above_threshold.fa', 'r') as g:
        allids = {i[1:] for i in g.read().split('\n') if i and i[0] == '>'}

    # set of IDs seen so far
    seenids = set()

    # file of top blast hits
    tophitsfile = open(args.outputdir + '/top.concat.txt', 'w')
    # file filtered with Ioan's prescription
    ifilterfile = open(args.outputdir + '/ifilter.concat.txt', 'w')

    # Ioan: Top number of BLAST hits to parse through in order to determine whether top hit can be trusted as truly non-human
    topchunk = 10
    # a counter
    minicounter = 0
    # a line representing a top hit
    topline = None
    # a filtering boolean (if true, filter out line)
    filterbool = False

    # glob blast files
    myfiles = glob.glob(args.outputdir + '/*.result')
    f = fileinput.input(files=myfiles)
    for line in f:
        linelist = line.strip().split()
        myid = linelist[0]
        # ID not yet seen (i.e., is top hit)
        if not myid in seenids:
            tophitsfile.write(line)
            seenids.add(myid)

            # this will skip on the loop's first iteration
            if topline:
                # print previous top line
                if not filterbool:
                    ifilterfile.write(topline + '\n')

            # here we're assuming fmt is:
            # qseqid, sseqid, saccver, staxids, pident, nident, length, mismatch, gapopen, gaps, qstart, qend, qlen, qframe, qcovs, sstart, send, slen, sframe, sstrand, evalue, bitscore, stitle
            topbitscore = float(linelist[21])
            topstaxids = linelist[3]
            topline = line.strip()

            # reset counter
            minicounter = 0
            # reset boolean (if tophit human, preemptively filter)
            if topstaxids == humantaxid:
                filterbool = True
            else:
                filterbool = False
            #print(myid + ' ' + str(minicounter) + ' ' + str(filterbool))
        # keep on checking results if filter flag hasn't gone high and #lines < topchunk
        elif (not filterbool) and minicounter < topchunk:
            # here we're assuming fmt is:
            # qseqid, sseqid, saccver, staxids, pident, nident, length, mismatch, gapopen, gaps, qstart, qend, qlen, qframe, qcovs, sstart, send, slen, sframe, sstrand, evalue, bitscore, stitle
            staxids = linelist[3]
            evalue = float(linelist[20])
            bitscore = float(linelist[21])
            filterbool = ioanfilter(staxids, evalue, bitscore, topbitscore)
            minicounter += 1
            #print(myid + ' ' + str(minicounter) + ' ' + str(filterbool))

    # do last entry
    if not filterbool:
        ifilterfile.write(topline)

    f.close()
    tophitsfile.close()
    ifilterfile.close()

    # set of IDs that didn't blast
    # print(allids)
    # print(seenids)
    noblastids = allids - seenids

    # get fasta file of entries that didn't blast
    filecount = hp.fastaidfilter(args.outputdir + '/above_threshold.fa', args.outputdir + '/no_blastn.fa', noblastids)

    if not args.noclean:
        cmd = 'rm {args.outputdir}/*.result {args.outputdir}/*.fasta'.format(args=args)
        hp.run_cmd(cmd, args.verbose, 0)

    print('No blast hits for: ' + ', '.join(list(noblastids)))

    if args.nosge:
        pass
    else:
        # concat blast logs and remove folder
        print('concatenate blast logs')
        cmd = 'head -100 {args.logsdir}/* > {args.outputdir}/log.blast'.format(args=args)
        hp.run_cmd(cmd, args.verbose, 0)

        if not args.noclean:
            shutil.rmtree(args.logsdir)

        print('CONCATENATE END')

# -------------------------------------

def ioanfilter(staxids, evalue, bitscore, top_score):
    """
    Ioan's filter:
    Only include the result in report of top blast hits if it is NOT a suspected high-alignment score human read
    """
    # humantaxid = '9606'
    evalcutoff = 10**(-4)
    return staxids == humantaxid and evalue < evalcutoff and bitscore > 0.95 * top_score

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # blast
    blast(args)
    #concat(args)

# -------------------------------------

if __name__ == '__main__':

    main()
