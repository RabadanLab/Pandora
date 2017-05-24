#!/usr/bin/env python

import sys

# This script takes taxids and traces back up through the taxonomic heirarchy
# E.g., 
# $ echo "9606" | ./getphylo.py nodes.dmp names.dmp
# [('9606', 'Homo sapiens', 'species'), ('9605', 'Homo', 'genus'), ('207598', 'Homininae', 'subfamily'), ('9604', 'Hominidae', 'family'), ('314295', 'Hominoidea', 'superfamily'), ('9526', 'Catarrhini', 'parvorder'), ('314293', 'Simiiformes', 'infraorder'), ('376913', 'Haplorrhini', 'suborder'), ('9443', 'Primates', 'order'), ('314146', 'Euarchontoglires', 'superorder'), ('1437010', 'Boreoeutheria', 'no rank'), ('9347', 'Eutheria', 'no rank'), ('32525', 'Theria', 'no rank'), ('40674', 'Mammalia', 'class'), ('32524', 'Amniota', 'no rank'), ('32523', 'Tetrapoda', 'no rank'), ('1338369', 'Dipnotetrapodomorpha', 'no rank'), ('8287', 'Sarcopterygii', 'no rank'), ('117571', 'Euteleostomi', 'no rank'), ('117570', 'Teleostomi', 'no rank'), ('7776', 'Gnathostomata', 'no rank'), ('7742', 'Vertebrata', 'no rank'), ('89593', 'Craniata', 'subphylum'), ('7711', 'Chordata', 'phylum'), ('33511', 'Deuterostomia', 'no rank'), ('33213', 'Bilateria', 'no rank'), ('6072', 'Eumetazoa', 'no rank'), ('33208', 'Metazoa', 'kingdom'), ('33154', 'Opisthokonta', 'no rank'), ('2759', 'Eukaryota', 'superkingdom'), ('131567', 'cellular organisms', 'no rank')]

# arg1 is the NCBI nodes.dmp file
# arg2 is the NCBI names.dmp file
# download it: wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# nodes.dmp
nodesfile = sys.argv[1]
# names.dmp
namesfile = sys.argv[2]

# taxid to (parent taxid, rank)
id2parentrank = {}
# taxid to name 
id2name = {}

# don't want to loop more than this many iterations through the taxonomic heirarchy
# in case taxid not found 
maxiterations = 50 

# open the nodes.dmp file
with open(nodesfile, 'r') as f:
    for line in f:
        # eliminate tabs and split on |
        l = line.replace('\t', '').split('|')
        # map taxid to (parent taxid, rank)
        id2parentrank[l[0]] = (l[1], l[2])

# open the names.dmp file
with open(namesfile, 'r') as f:
    for line in f:
        l = line.replace('\t', '').split('|')
        if (l[3] == 'scientific name'):
            id2name[l[0]] = l[1]

# print('.dmp files loaded')

for line in sys.stdin:
    # tmp list
    l = []
    # reset iterations
    myiterations = 0
    current_taxid = line.strip()
    while (current_taxid != '1' and myiterations < 50):
        # get current taxid and its level in the taxonomic heirarchy
        # l.append(current_taxid + ':' + id2parentrank[current_taxid][1])
        l.append((current_taxid, id2name.get(current_taxid, 'not found'), id2parentrank.get(current_taxid, ['', 'not found'])[1]))
        # get parent taxid 
        current_taxid = id2parentrank.get(current_taxid, ['1'])[0] 
        myiterations += 1
    print(l)

