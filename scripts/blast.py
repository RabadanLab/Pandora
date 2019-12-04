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
    vars(args)['prefix'] = args.outputdir + '/blast'
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
    cmd = '{args.scripts}/blast_2.sh {args.whichblast} {args.prefix} {args.db} {args.threads}'.format(args=args, flag=flag)
    hp.run_cmd(cmd, args.verbose, 0)
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
