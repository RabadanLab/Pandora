#!/usr/bin/env python
#$ -V
#$ -cwd

# A simple wrapper for blast array job 

import argparse
import sys
import os

# -------------------------------------

def get_arg():
    """
    Get Arguments

    :rtype: object
    """
    # parse arguments

    # if SGE_TASK_ID is not defined in the shell, this avoids error
    try:
        sgeid_default = os.environ['SGE_TASK_ID']
    except:
        sgeid_default = 'undefined'

    prog_description = 'A wrapper for blast'
    parser = argparse.ArgumentParser(description=prog_description)
    # parser.add_argument('-i', '--input', help='the input fasta')
    parser.add_argument('-o', '--outputdir', default='blast', help='the output directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--sgeid', default=sgeid_default, help='SGE_TASK_ID (set by hand only if qsub is turned off)')
    parser.add_argument('--whichblast', default='blastn', choices=['blastn', 'blastp'], help='which blast to use (blastn, blastp)')
    parser.add_argument('--threads', default='1', help='blast -num_threads option')
    parser.add_argument('--db', help='the database prefix')
    parser.add_argument('--fmt', help='blast format string')
    parser.add_argument('--noclean', type=int, default=0, help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--verbose', type=int, default=0, help='verbose mode: echo commands, etc (default: off)')
    args = parser.parse_args()

    # add key-value pairs to the args dict
    vars(args)['input1'] = args.outputdir + '/blast_1' + '.fasta'
    vars(args)['input2'] = args.outputdir + '/blast_2' + '.fasta'
    vars(args)['step'] = 'blast'

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # error checking: exit if previous step produced zero output
    for i in [args.input]:
        hp.check_file_exists_and_nonzero(i, step=args.step)

    return args

# -------------------------------------

def blastnp(args):
    """Blastn or blastp"""

    hp.echostep(args.step)

    # extra blast flags
    flag = ''

    # if blastp, use blastp-fast
    if args.whichblast == 'blastp':
        flag = '-task blastp-fast'

    # do blastn or blastp
    cmd1 = '{args.whichblast} -outfmt "6 {args.fmt}" -query {args.input1} -db {args.db} -num_threads {args.threads} {flag} > {args.outputdir}/blast_1.result'.format(args=args, flag=flag)
    cmd2 = '{args.whichblast} -outfmt "6 {args.fmt}" -query {args.input2} -db {args.db} -num_threads {args.threads} {flag} > {args.outputdir}/blast_2.result'.format(args=args, flag=flag)

    hp.echostep(args.step, start=0)

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # blast
    blastnp(args)

# -------------------------------------

if __name__ == '__main__':

    main()
