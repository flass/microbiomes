#!/usr/bin/python
import os, re, sys
from Bio import Phylo

def get_all_children(clade):
	children = [clade]
	for child in clade.clades:
		children += get_all_children(child)
	return children

resepcadir = sys.argv[1]
pcasca = sys.argv[2]


if pcasca=='scaled_abundances':
	scaling=3
else:
	scaling=1

nfphyloxml = '%s/alledges.33samples.ePCA.%s.4.PC.xml'%(resepcadir, pcasca)
for pcs in ['12', '34']:
	nfcolwid = '%s/alledges.33samples.ePCA-%s.PC%s.colourwidth.csv'%(resepcadir, pcasca, pcs)
	dedgecolwid = {}
	with open(nfcolwid, 'r') as fcolwid:
		for i, line in enumerate(fcolwid):
			dedgecolwid[i] = line.strip('\n').split('\t')
	
	trees = Phylo.parse(nfphyloxml,'phyloxml')
	ltrees = [tree for tree in trees]
	reftree = ltrees[0] # PC1
	
	nfoutrootpath = '%s/alledges.33samples.old2newrootpath.invertededges'%(resepcadir)
	foutrootpath = open(nfoutrootpath, 'w')
	
	bracknumre = re.compile('\{([0-9]+)\}$')
	widthresh = 5e-2
	for clade in get_all_children(reftree.root):
		cladeid = clade.node_id
		if not cladeid:
			cladeid = int(bracknumre.search(clade.name).group(1))
		col = dedgecolwid[cladeid][0]
		wid = float(dedgecolwid[cladeid][1])*scaling
		if wid>widthresh:
			clade.color = col
			clade.width = wid*30
		else:
			clade.color = None
			clade.width = None
	
	#nfphyloxmlout = '%s/alledges.33samples.ePCA.%s.PC%s.xml'%(resepcadir, pcasca, pcs)
	nfphyloxmlout = '%s/alledges.33samples.ePCA.rerooted.%s.PC%s.xml'%(resepcadir, pcasca, pcs)
	Phylo.write(reftree,nfphyloxmlout,'phyloxml')
