#!/usr/bin/env python

from sys import stdout, stderr, exit
from optparse import OptionParser
import shelve
from psycho import node



if __name__ == '__main__':
    usage = '%prog [options] <SHELVE>'
    parser = OptionParser(usage=usage)

    (options, args) = parser.parse_args()

    # check if input file is given, if so, call function with data
    if len(args) != 1:
        parser.print_help()
        exit(1)


    shObj = shelve.open(args[0], flag='r', protocol=-1)
    ref = shObj['ref']

    #
    # load inclusion tree
    #
    root = shObj['inclusion_tree']
    
    queue = list(map(lambda x: (x, 1), root.children))
    depths = list()
    while queue:
        u, u_depth = queue.pop()
        for child in u.children:
            if len(child.children) > 1:
                queue.append((child, u_depth+1))
            else:
                depths.append(u_depth+1)

    print >> stdout, 'max tree depth: %s' %max(depths)
    print >> stdout, 'avg tree depth: %s' %(sum(depths)/float(len(depths)))
