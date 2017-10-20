#!/usr/bin/python
import os, sys, getopt
ncbitaxodmpdir = '/home/flassall/NCBI/Taxonomy'

fnames = open(os.path.join(ncbitaxodmpdir, 'names.dmp'), 'r')
fnodes = open(os.path.join(ncbitaxodmpdir, 'nodes.dmp'), 'r')
fmerged = open(os.path.join(ncbitaxodmpdir, 'merged.dmp'), 'r')

opts, args = getopt.getopt(sys.argv[1:], 'f:d:', ['filter-tax='])
fin = open(args[0], 'r') if len(args)>0 else sys.stdin
fout = open(args[1], 'w') if len(args)>1 else sys.stdout
dopt = dict(opts)
field2sub = [int(i)-1 for i in dopt.get('-f', '1').split(',')]
delim = dopt.get('-d', '\t')
filtertax = dopt.get('--filter-tax', '').split(',')

dnames = {}
sys.stderr.write("load taxonomy names\n")
for line in fnames:
	lsp = line.rstrip('\t|\n').split('\t|\t')
	if lsp[3]!='scientific name': continue
	dnames[int(lsp[0])] = lsp[1]

fnames.close()

dmerged = {}
sys.stderr.write("load taxonomy merged node pairs\n")
for line in fmerged:
	lsp = line.rstrip('\t|\n').split('\t|\t')
	dmerged[int(lsp[0])] = int(lsp[1])

fmerged.close()

dnodes = {}
sys.stderr.write("load taxonomy node relationships\n")
for line in fnodes:
	lsp = line.rstrip('\t|\n').split('\t|\t')
	p = int(lsp[1])
	dnodes[int(lsp[0])] = p if p!=1 else None

fnodes.close()



def taxostr(taxid, sep='|', top2bottom=True):
	i = dmerged.get(taxid, taxid)
	l = [dnames.get(i, str(i))]
	n = dnodes.get(i)
	while n:
		i = dmerged.get(n, n)
		l.append(dnames[i])
		#~ sys.stderr.write(str(i))
		#~ sys.stderr.flush()
		n = dnodes.get(i)
	if top2bottom: l.reverse()
	s = sep.join(l)
	return s

for line in fin:
	lsp = line.rstrip('\n').split(delim)
	lout = []
	for k, val in enumerate(lsp):
		if (k in field2sub) and (val not in filtertax):
			try:
				tid = int(val)
				ths = taxostr(tid)
			except ValueError:
				ths = val
			lout.append(ths)
		else:
			lout.append(val)
	fout.write(delim.join(lout)+'\n')

fin.close()
fout.close()

