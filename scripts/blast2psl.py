#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, stdin, exit
from stat import S_ISFIFO
import os

from Bio import SearchIO
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment

import re

PAT_LENGTH     = re.compile('^\s*Length\s*=\s*(\d+)\s*$')
PAT_HSP_START  = re.compile('^\s*Score\s*=\s*\d')
PAT_HSP_INFO   = re.compile('(\w+)\s*=\s*([^,=]+)\s*(?=,|$)') 
PAT_HSP_SCORE  = re.compile('^([0-9.e-]+) bits \(([0-9.e-]+)\)\s*$')
PAT_HSP_COUNTS = re.compile('^\s*(\d+)/(\d+) \([0-9.]+%\)\s*$')
PAT_HSP_STRAND = re.compile('^\s*(Plus|Minus)\s*/\s*(Plus|Minus)\s*$')
PAT_ALN_START  = re.compile('^\s*Query\s+\d')
PAT_ALN_SEQ    = re.compile('^(\s*)(Query|Sbjct)(\s+)(\d+)(\s+)([A-Za-z-]+)\s+(\d+)\s*$')
PAT_ALN_MATCH  = re.compile('^[ |]+$')


STRAND = {'Plus': 1, 'Minus': -1}

class ParseType:

    QUERY = 0
    HIT = 1
    HIT_ID = 2
    HSP = 3
    HSP_END = 4
    ALN = 5
    ALN_QUERY = 6
    ALN_MATCH = 7
    ALN_SBJCT = 8

def parseData(data, out):

    stack = []

    line = data.readline()
    while line or stack:
        if not stack:
            if line.startswith('Query='):
                qid = line[6:].strip().split(' ', 1)
                query = QueryResult(id = qid[0])
                if len(qid) > 1:
                    query.description = qid[1]

                setattr(query, 'seq_len', 0)
                stack.append((ParseType.QUERY, query))
                stack.append((ParseType.HIT, None))
        else:
            STATE = stack[-1][0]

            if STATE == ParseType.QUERY:
                if not line or line.startswith('Query='):
                    _, query = stack.pop()
                    SearchIO.write(query, out, 'blat-psl')
                    continue

            elif STATE == ParseType.HIT:
                if not line or line.startswith('Query='):
                    _, hit = stack.pop()
                    stack[-1][1].append(hit)
                    continue
                hit_id = None

                # XXX not really clean -- this should be in QUERY state
                m = PAT_LENGTH.match(line)
                if m:
                    stack[-2][1].seq_len = int(m.group(1))
                else:
                    if line.startswith('>'):
                        hit_id = line[1:].strip()
                    elif line.startswith('Subject='):
                        hit_id = line[8:].strip()
                    if hit_id:
                        if stack[-1][1] != None:
                            _, hit = stack.pop()
                            stack[-1][1].append(hit)
                            stack.append((ParseType.HIT, None))
                        stack.append((ParseType.HIT_ID, hit_id))
                    elif PAT_HSP_START.match(line):
                        stack.append((ParseType.HSP, dict()))
                        continue

            elif STATE == ParseType.HIT_ID:
                m = PAT_LENGTH.match(line)
                if m:
                    _, hit_id = stack.pop()
                    hid = hit_id.split(' ', 1)
                    hit = Hit(id = hid[0], query_id = stack[-2][1].id)
                    if len(hid) > 1:
                        hit.description = hid[1]
                    setattr(hit, 'seq_len', int(m.group(1)))
                    stack[-1] = (ParseType.HIT, hit)
                    stack.append((ParseType.HSP, dict()))
                else:
                    stack[-1] = (ParseType.HIT_ID, stack[-1][1] + line[:-1])

            elif STATE == ParseType.HSP and (not line or line.strip()):
                pat_hits = PAT_HSP_INFO.findall(line)
                if pat_hits:
                    info_dict = stack[-1][1] 
                    for key, value in pat_hits:
                        if key == 'Score':
                            info_dict[key] = map(float,
                                    PAT_HSP_SCORE.match(value).groups())
                        elif key == 'Strand':
                            info_dict[key] = map(STRAND.get,
                                    PAT_HSP_STRAND.match(value).groups())
                        elif key == 'Expect':
                            info_dict[key] = float(value) 
                        else:
                            m = PAT_HSP_COUNTS.match(value)
                            if m:
                                info_dict[key] = map(int, m.groups())
                else:
                    if PAT_ALN_START.match(line):
                        aln_len = stack[-1][1].get('Identities', (0, -1))[1]
                        stack.append((ParseType.ALN, list(), aln_len))
                        continue
                    else:
                        raise SyntaxError(('expected alignment, but found ' + \
                                'line \'%s\'') %line[:-1])

            elif STATE == ParseType.HSP_END: 
                _, hsp = stack.pop()
                hit = stack[-1][1]
                hsp.query_id = hit.query_id
                hsp.hit_id = hit.id
                setattr(hsp, 'n_num', sum(map(len, hsp.fragments)))
                setattr(hsp, 'match_num', sum(map(lambda x: x.match_num,
                    hsp.fragments)))
                setattr(hsp, 'mismatch_num', sum(map(lambda x: x.mismatch_num,
                    hsp.fragments)))
                setattr(hsp, 'match_rep_num', 0)
                setattr(hsp, 'query_gap_num', hsp.query_span + 1 - hsp.n_num)
                setattr(hsp, 'hit_gap_num', hsp.hit_span + 1 - hsp.n_num)

                qgapopen = 0
                sgapopen = 0
                for i in xrange(1, len(hsp.fragments)):
                    pf = hsp.fragments[i-1] 
                    f = hsp.fragments[i] 
                    if (f.query_strand == 1 and pf.query_end + 1 != \
                            f.query_start) or (f.query_strand == -1 and \
                            pf.query_start - 1 != f.query_end):
                        qgapopen += 1
                    if (f.hit_strand == 1 and pf.hit_end + 1 != \
                            f.hit_start) or (f.hit_strand == -1 and \
                            pf.hit_start - 1 != pf.hit_end):
                        sgapopen += 1

                setattr(hsp, 'query_gapopen_num', qgapopen)
                setattr(hsp, 'hit_gapopen_num', sgapopen)
                setattr(hsp, '_has_hit_strand', True)

                hit.append(hsp)
                continue
                    
            elif STATE == ParseType.ALN:
                aln_len = stack[-1][2]
                m = PAT_ALN_SEQ.match(line)
                if m:
                    s1, seq_type, s2, start, s3, aln, end = m.groups()
                    if seq_type != 'Query':
                        raise SyntaxError('expected alignment sequence of' + \
                                'Query, but found %s' %seq_type)

                    stack.append((ParseType.ALN_QUERY, aln, int(start),
                        int(end), len(s1) + len(seq_type) + len(s2) + len(start)
                        + len(s3)))
                elif aln_len == 0 or not line or line.strip():
                    if aln_len > 0:
                        raise SyntaxError(('expected %s more character of' + \
                                'alignment') %aln_len)
                    # remove ALN from stack
                    _, frags, _, = stack.pop()
                    # remove HSP from stack
                    _, info_dict = stack.pop()
                    hsp = HSP(frags)
                    if info_dict.has_key('Expect'):
                        hsp.evalue = info_dict['Expect']
                    if info_dict.has_key('Score'):
                        hsp.bitscore = info_dict['Score'][0]
                    stack.append((ParseType.HSP_END, hsp))
                    continue

            elif STATE == ParseType.ALN_QUERY:
                aln = stack[-1][1]
                indent = stack[-1][4]
                m = PAT_ALN_MATCH.match(line[indent:indent+len(aln)])
                if not m:
                    raise SyntaxError('expected alignment match of, but ' + \
                            'found %s' %line[:-1])
                stack.append((ParseType.ALN_MATCH, m.group(0)))

            elif STATE == ParseType.ALN_MATCH:
                m = PAT_ALN_SEQ.match(line)
                if m:
                    _, seq_type, _, sstart, _, sbjct, send = m.groups()
                    sstart, send = int(sstart), int(send)
                    if seq_type != 'Sbjct':
                        raise SyntaxError('expected alignment sequence of' + \
                                'Sbjct, but found %s' %seq_type)

                    # ALN_MATCH
                    _, match = stack.pop()
                    # ALN_QUERY
                    _, query, qstart, qend, _ = stack.pop()
                    # ALN
                    _, pfrags, aln_len = stack[-1]
                    # HSP
                    info_dict = stack[-2][1]

                    frags = list()

                    i, j = 0, 0

                    mismatch_num = 0
                    while i < len(match):
                        j = match.find(' ', j)
                        if j == -1 or query[j] == '-' or sbjct[j] == '-':
                            qstrand, sstrand = info_dict.get('Strand', (0, 0))
                            j = j == -1 and len(match) or j
                            f = HSPFragment(hit=sbjct[i:j], query=query[i:j])
                            f.query_strand = qstrand
                            f.hit_strand = sstrand
                            if qstart < qend:
                                f.query_start = qstart + i
                                f.query_end = qstart + j - 1
                                f.query_strand = 1
                            else:
                                f.query_end = qstart - i
                                f.query_start = qstart - j + 1
                                f.query_strand = -1
                            if sstart < send:
                                f.hit_start = sstart + i
                                f.hit_end = sstart + j - 1
                                f.hit_strand = 1
                            else:
                                f.hit_end = sstart - i
                                f.hit_start = sstart - j + 1
                                f.hit_strand = -1
                            setattr(f, 'mismatch_num', mismatch_num)
                            setattr(f, 'match_num', j-i-mismatch_num)

                            frags.append(f)

                            while j < len(query) and query[j] == '-':
                                j += 1
                                # increase if on plus strand, decrease otherwise
                                qstart -= qstrand
                            while j < len(sbjct) and sbjct[j] == '-':
                                j += 1
                                sstart -= sstrand
                            i = j
                        else:
                            mismatch_num += 1
                        j += 1


                    if pfrags and ((pfrags[-1].query_strand == 1 and \
                            pfrags[-1].query_end+1 == frags[0].query_start) or \
                            (pfrags[-1].query_strand == -1 and  \
                            pfrags[-1].query_start-1 == frags[0].query_end)) \
                            and ((pfrags[-1].hit_strand == 1 and \
                            pfrags[-1].hit_end+1 == frags[0].hit_start) or \
                            (pfrags[-1].hit_strand == -1 and  \
                            pfrags[-1].hit_start-1 == frags[0].hit_end)):

                        query_start = None
                        query_end = None
                        if pfrags[-1].query_strand == 1:
                            query_start = pfrags[-1].query_start
                            query_end = frags[0].query_end 
                        else:
                            query_start = frags[0].query_start
                            query_end = pfrags[-1].query_end 

                        sbjct_start = None
                        sbjct_end = None
                        if pfrags[-1].hit_strand == 1:
                            sbjct_start = pfrags[-1].hit_start 
                            sbjct_end = frags[0].hit_end
                        else:
                            sbjct_start = frags[0].hit_start 
                            sbjct_end = pfrags[-1].hit_end

                        f = HSPFragment(query=pfrags[-1].query + frags[0].query,
                                hit=pfrags[-1].hit + frags[0].hit)
                        f.query_strand = frags[0].query_strand
                        f.query_start = query_start
                        f.query_end = query_end
                        f.hit_strand = frags[0].hit_strand
                        f.hit_start = sbjct_start
                        f.hit_end = sbjct_end
                        setattr(f, 'mismatch_num', frags[0].mismatch_num +
                                pfrags[-1].mismatch_num)
                        setattr(f, 'match_num', frags[0].match_num +
                                pfrags[-1].match_num)
                        pfrags[-1] = f 
                        pfrags.extend(frags[1:])
                    else:
                        pfrags.extend(frags)

                    aln_len -= len(match)
                    stack[-1] = (ParseType.ALN, pfrags, aln_len)
        line = data.readline()


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-i', '--input', type=str, 
            help='Input BLAST file (by default, input is read from stdin')

    args = parser.parse_args()

    data = stdin

    if not S_ISFIFO(os.fstat(0).st_mode) and args.input:
        data = open(args.input)
    elif not S_ISFIFO(os.fstat(0).st_mode):
        print >> stderr, 'FATAL: Input required\n'
        parser.print_usage()
        exit()

    parseData(data, stdout)
