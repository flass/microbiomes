#!/usr/bin/python
import os, re, sys
from Bio import Phylo

# Reference tree has to be re-rooted between arbitrary root and the 3-domain root of the common Tree of Life
# guppy-procued edgediff matrix will have to  be modified upon reading in R to change sign of edge mass difference entries
# for edges on the path between old and new root (Matsen and Evans 2013 PLOS One 8(3): e56859. doi:10.1371/journal.pone.0056859)
oldrootlab = 'ENTROBACTER_SP._638[399742]{0}'
newrootlab = '_BACTERIA_{6473}'

resepcadir = os.environ['resepca']
#resepcadir = sys.argv[1]

nfphyloxml = '%s/33samples.xml'%(resepcadir, pcasca)

trees = Phylo.parse(nfphyloxml,'phyloxml')
ltrees = [tree for tree in trees]
reftree = ltrees[0] # PC1

# find path from old to new root
#~ rootpath = [reftree.root]+reftree.root.get_path(newrootlab)
rootpath = reftree.root.get_path(newrootlab)

nfoutrootpath = '%s/alledges.33samples.old2newrootpath.invertededges'%(resepcadir)
with open(nfoutrootpath, 'w')as foutrootpath:
	for clade in get_all_children(reftree.root):
		if clade in rootpath:
			# write out edge index
			foutrootpath.write('%d\t%s\n'%(cladeid, clade.name))

