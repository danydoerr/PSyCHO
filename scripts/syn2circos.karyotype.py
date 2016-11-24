#!/usr/bin/env python

from sys import stdout, stderr, exit, argv
from os.path import basename 
from Bio import SeqIO
import re


CHR_PAT = re.compile('.*\|chromosome\|([^\|]+)(\|.*|$)')
POS_PAT = re.compile('.*\|(\d+):(\d+)(\|.*|$)')

if __name__ == '__main__':
    
    if len(argv) != 3:
        print '\tusage: %s <FASTA FILE GENOME> <FASTA FILE SEGMENTS>' %argv[0]
        exit(1)

    gid = basename(argv[1]).split('.', 1)[0]
    chr_lengths = dict()
    for rec in SeqIO.parse(argv[1], 'fasta'):
        sid = rec.id.replace('|', '.').lower()
        print 'chr - %s.%s %s.%s %s %s %s.%s' %(gid.lower(), sid, gid, sid, 0,
                len(rec), gid.lower(), sid)
        chr_lengths[sid] = len(rec)

    prev = 0
    prevC = None
    c = 0
    for rec in SeqIO.parse(argv[2], 'fasta'):
        if rec.description.find('[inactivate]') >= 0:
            continue

        chr1 = CHR_PAT.match(rec.description).group(1).lower()

        mpos = POS_PAT.match(rec.description)
        start = int(mpos.group(1))
        end = int(mpos.group(2))+1
        if prevC != None and prevC != chr1:
            if chr_lengths[prevC]-prev > 0:
                bid = '%s.%s.no_seg%s' %(gid.lower(), prevC, c)
                print 'band %s.%s %s %s %s %s %s.no_seg' %(gid.lower(),
                        prevC.lower(), bid, bid, prev, chr_lengths[prevC],
                        gid.lower())
                c += 1
            prev = 0
        prevC = chr1
        if start-prev > 0:
            bid = '%s.%s.no_seg%s' %(gid.lower(), prevC, c)
            print 'band %s.%s %s %s %s %s %s.no_seg' %(gid.lower(),
                    prevC.lower(), bid, bid, prev, start,
                    gid.lower())
            c += 1
        elif start-prev < 0:
            if end-prev <= 0:
                print >> stderr, 'Segment %s is fully contained in previous one. Skipping.' %rec.description
                continue
            else:
                print >> stderr, 'Segment %s is overlapping with previous one. Trimming.' %rec.description
                start = prev


        if end > chr_lengths[chr1]:
            end = chr_lengths[chr1]

        bid = '%s.%s' %(gid.lower(),
                rec.description[:rec.description.find('|')].lower())
        print 'band %s.%s %s %s %s %s %s.seg' %(gid.lower(), chr1, bid,
                bid, start, end, gid.lower())
        prev = end

    if chr_lengths[prevC]-prev > 0:
        bid = '%s.%s.no_seg%s' %(gid.lower(), prevC, c)
        print 'band %s.%s %s %s %s %s %s.no_seg' %(gid.lower(),
                prevC.lower(), bid, bid, prev, chr_lengths[prevC],
                gid.lower())
