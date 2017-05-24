#!/usr/bin/env python

import argparse
import sys

# This script performs assembly on the reads leftover after host removal

# -------------------------------------

def get_arg():
    """Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'Trinity assembly'
    parser = argparse.ArgumentParser(description=prog_description)

    ## Trinity run without gunzipping the fastq files
    parser.add_argument('-1', '--mate1', default='host_separation/unmapped_1.fastq.gz', help='mate1')
    parser.add_argument('-2', '--mate2', default='host_separation/unmapped_2.fastq.gz', help='mate2')
    parser.add_argument('--single', default=None)
    parser.add_argument('-o', '--outputdir', default='assembly_trinity', help='the output directory')
    parser.add_argument('--trinitymem', required=True, help='max memory for Trinity')
    parser.add_argument('--trinitycores', required=True, help='number of cores for Trinity')
    parser.add_argument('--trinitythreshold', required=True, help='number of cores for Trinity')
    parser.add_argument('-l', '--logsdir', help='the logs directory')
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('--noclean', help='do not delete temporary intermediate files (default: off)')
    parser.add_argument('--verbose', type=int, default=0, help='verbose mode: echo commands, etc (default: off)')
    args = parser.parse_args()

    # add key-value pairs to the args dict
    vars(args)['step'] = 'assembly'

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    global ahp
    from helpers import helpers as hp
    from helpers import assembly_helpers as ahp

    # error checking: exit if previous step produced zero output

    if (args.single):
        hp.check_file_exists_and_nonzero(args.mate1, step=args.step)
    else:
        for i in [args.mate1, args.mate2]:
            hp.check_file_exists_and_nonzero(i, step=args.step)

    # this silly line casts the string False to the boolean value
    if args.single == 'False':
        args.single = False

    return args

# -------------------------------------

def assembly(args):
    """Do Trinity assembly"""

    hp.echostep(args.step)

    # print args
    print(args)
    print

    # mkdir -p assembly_trinity
    hp.mkdirp(args.outputdir)

    # perform Trinity assembly
    if (args.single):
        cmd = 'Trinity --seqType fq --normalize_reads --min_contig_length={args.trinitythreshold} --max_memory {args.trinitymem}G --CPU {args.trinitycores} --output {args.outputdir} --single {args.mate1}'.format(args=args)
    else:
        cmd = 'Trinity --seqType fq --normalize_reads --min_contig_length={args.trinitythreshold} --max_memory {args.trinitymem}G --CPU {args.trinitycores} --output {args.outputdir} --left {args.mate1} --right {args.mate2}'.format(args=args)
    # use run_long_cmd for programs with verbose output
    hp.run_long_cmd(cmd, args.verbose, 'log.Trinity')

    print('Trinity complete')

    # name of expected output file
    myoutput = args.outputdir + '/Trinity.fasta'

    # exit if no output
    hp.check_file_exists_and_nonzero(myoutput)

    # mkdir -p
    hp.mkdirp('assembly')

    # rename Trinity contigs, join sequence portion of fasta, return number of contigs
    # cat ${outputdir}/Trinity.fasta | awk 'BEGIN{f=0; counter=1}{if ($0~/^>/) {if (f) {printf "\n"; counter++}; print ">contig_"counter; f=1} else printf $0}END{printf "\n"}' > ${output}
    myoutput2 = 'assembly/contigs_trinity.fasta'
    num_contigs = hp.fastajoinlines(myoutput, myoutput2, 'contig')

    # compute simple distribution
    # cat assembly/contigs_trinity.fasta | paste - - | awk '{print length($2)}' | sort -nr | ${d}/scripts/tablecount | awk -v tot=${num_contigs} 'BEGIN{x=0}{x+=$2; print $1"\t"$2"\t"x"/"tot"\t"int(100*x/tot)"%"}' > assembly/contigs.distrib.txt
    ahp.computedistrib(myoutput2, 'assembly/contigs.distrib.txt')

    if not int(args.noclean):
        cmd = 'rm -rf assembly_trinity'
        hp.run_cmd(cmd, args.verbose, 0)

    hp.echostep(args.step, start=0)

    # return the name of assembly file
    return myoutput2

# -------------------------------------

def remap(args, contigs):
    """map contigs back onto assembly"""

    hp.echostep('remap')

    hp.mkdirp('assembly/ref_remap')

    refbowtie="assembly/ref_remap/ref"

    cmd = 'bowtie2-build {} {}'.format(contigs, refbowtie)
    hp.run_cmd(cmd, args.verbose, 0)

    if (args.single):
        cmd = 'bowtie2 -p 4 -x {} -U {} -S {}'.format(refbowtie, args.mate1, 'assembly/reads2contigs.sam')
    else:
        cmd = 'bowtie2 -p 4 -x {} -1 {} -2 {} -S {}'.format(refbowtie, args.mate1, args.mate2, 'assembly/reads2contigs.sam')
    hp.run_cmd(cmd, args.verbose, 0)

    # convert to bam
    ## samtools version compatibility: need .bam extension
    cmd = 'samtools view -bS assembly/reads2contigs.sam | samtools sort -o assembly/reads2contigs.bam'
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'samtools index assembly/reads2contigs.bam'
    hp.run_cmd(cmd, args.verbose, 0)

    cmd = 'rm assembly/reads2contigs.sam'
    hp.run_cmd(cmd, args.verbose, 0)

    # BAM index stats
    cmd = 'samtools idxstats assembly/reads2contigs.bam > assembly/reads2contigs.stats.txt'
    hp.run_cmd(cmd, args.verbose, 0)

    # mpileup
    cmd = 'samtools mpileup -A -B -d 100000 -L 100000 -f assembly/contigs_trinity.fasta assembly/reads2contigs.bam > assembly/reads2contigs.pileup'
    hp.run_cmd(cmd, args.verbose, 0)

    # format pileup file - i.e., add zeros to uncovered positions
    ahp.formatpileup('assembly/reads2contigs.pileup', 'assembly/reads2contigs.stats.txt', 'assembly/reads2contigs.format.pileup', 'assembly/reads2contigs.entropy')

    if not int(args.noclean):
        cmd = 'rm -r assembly/ref_remap'
        hp.run_cmd(cmd, args.verbose, 0)
        cmd = 'rm assembly/reads2contigs.pileup'
        hp.run_cmd(cmd, args.verbose, 0)

    hp.echostep('remap', start=0)

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()
    # assembly
    contigs = assembly(args)
    # remap
    remap(args, contigs)

# -------------------------------------

if __name__ == '__main__':

    main()
