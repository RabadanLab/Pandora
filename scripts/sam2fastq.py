#!/usr/bin/env python

import sys
import os

# convert sam file to fastq
# for this sam file, these reads are unmapped, 
# so no need to worry about complementing or multiple mapping

# argument is output fastq base name
fastqbasename = sys.argv[1]
args_single = sys.argv[2]

# this silly line casts the string False to the boolean value
if args_single == 'False':
    args_single = False

# recapitulate:
# first mate
# samtools view -f 64 $1 | cut -f1,10,11 | awk 'BEGIN{OFS="\t";}{print "@"$1"/1"; print $2; print "+"; print $3}' > ${2}_1.fastq
# second mate
# samtools view -f 128 $1 | cut -f1,10,11 | awk 'BEGIN{OFS="\t";}{print "@"$1"/2"; print $2; print "+"; print $3}' > ${2}_2.fastq

# key: (https://samtools.github.io/hts-specs/SAMv1.pdf)
# 1 1 0x1 template having multiple segments in sequencing
# 2 2 0x2 each segment properly aligned according to the aligner
# 3 4 0x4 segment unmapped
# 4 8 0x8 next segment in the template unmapped
# 5 16 0x10 SEQ being reverse complemented
# 6 32 0x20 SEQ of the next segment in the template being reverse complemented
# 7 64 0x40 the first segment in the template
# 8 128 0x80 the last segment in the template
# 9 256 0x100 secondary alignment
# 10 512 0x200 not passing filters, such as platform/vendor quality controls
# 11 1024 0x400 PCR or optical duplicate
# 12 2048 0x800 supplementary alignment


fastq1 = fastqbasename + '_1.fastq'
fastq2 = fastqbasename + '_2.fastq'

with open(fastq1, 'w') as f, open(fastq2, 'w') as g:
    # loop thro std:in
    for line in sys.stdin:
        id = line.split()[0]
        flag = line.split()[1]
        myread = line.split()[9]
        qual = line.split()[10]

        if (args_single):
            # Ioan: Removing the '/1' read specification before the first '\n' character
            # possibly a source of formatting errors running Trinity with the --single flag
            f.write('@' + id + '\n' + myread + '\n' + '+\n' + qual + '\n')
        else:
            # get the binary value of the flag and bitwise AND with 64 ( i.e, int('1000000',2) )
            # and then shift 6 bits to look at the 7th bit, which tells us if mate 1 (refer to key).
            # I.e., (int('1000000',2) & 64) >> 6 == 1

            # if mate 1 read
            if (( int(flag) & 64 ) >> 6):
                f.write('@' + id + '/1\n' + myread + '\n' + '+\n' + qual + '\n')
            # else if mate 2 read
            elif (( int(flag) & 128 ) >> 7):
                g.write('@' + id + '/2\n' + myread + '\n' + '+\n' + qual + '\n')
