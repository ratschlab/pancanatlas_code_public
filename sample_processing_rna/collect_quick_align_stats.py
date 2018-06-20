import sys
import os
import pysam
import re

if len(sys.argv) < 2:
    print >> sys.stderr, "Usage: %s <align.bam>"
    sys.exit(1)
infile = sys.argv[1]

infile = pysam.Samfile(infile, 'rb')

aln_count = dict()
mm_count = dict()
spliced_count = dict()
total_count = dict()

for read in infile:
    if read.is_secondary:
        continue

    RG = read.get_tag("RG")

    try:
        total_count[RG] += 1
    except KeyError:
        total_count[RG] = 1
        mm_count[RG] = dict()
        aln_count[RG] = dict()
        spliced_count[RG] = 0

    if read.is_unmapped:
        try:
            aln_count[RG][0] += 1
        except KeyError:
            aln_count[RG][0] = 1
        continue
    if 'N' in read.cigarstring:
        spliced_count[RG] += 1

    mm = read.get_tag("NM")
    try:
        mm_count[RG][mm] += 1
    except KeyError:
        mm_count[RG][mm] = 1
    
    nh = read.get_tag("NH")
    try:
        aln_count[RG][nh] += 1
    except:
        aln_count[RG][nh] = 1

OUT = open(re.sub('bam$', '', sys.argv[1]) + 'stats', 'w')
for RG in total_count:
    print >> OUT, "%s\t%i\t%i\t%s\t%s" % (RG, total_count[RG], spliced_count[RG], ','.join(['%s:%s' % (key, value) for key, value in aln_count[RG].iteritems()]), ','.join(['%s:%s' % (key, value) for key, value in mm_count[RG].iteritems()]))
OUT.close()

