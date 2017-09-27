#!/usr/bin/python
# -*- coding: utf8 -*-

import os
import getopt, sys
import time

from parseKrona import *

def main(nfhtml, outputdir, lftaxa, letaxa=[], abuntag='abundances', concat=False, outmode='a', verbose=False):
	
	starttime = time.time()
	currtime = starttime
	
	cumabuntag = 'cumul_'+abuntag
	
	# create output directory if does not exist
	if not os.path.exists(outputdir):
		os.mkdir(outputdir)
	
	# parsing HTML Krona files
	print "Parsing taxonomy in Krona HTML result file: '%s'"%(nfhtml)
	fhtml = open(nfhtml, 'r')
	kronaparser = KronaHTMLParser() #verbose=verbose
	kronaparser.feed(fhtml.read())
	fhtml.close()
	taxonomy = kronaparser.taxonomy[0]
	if not taxonomy: raise IndexError
	
	excltaxnodes = [taxonomy.get(exclnameortaxid, taxonomy.taxidgetnode(exclnameortaxid)) for exclnameortaxid in letaxa]
	for i in range(len(lftaxa)):
		taxon = lftaxa[i]
		if concat and len(lftaxa)>1:
			outtag = 'all'	# generic output file prefix
			if i == 0: mode = outmode
			else: mode = 'a'	# append to file after first iteration
		else:
			mode = outmode
			outtag = taxon
		nfout = "%s/%s.phylosift_taxon_summary.tab"%(outputdir, outtag)
		taxnode = taxonomy.get(taxon, taxonomy.taxidgetnode(taxon)) # first consider the input as the name, then as a numerical taxid
		rootabun = taxnode.get_cumulative_values(abuntag, exclude=excltaxnodes, filterNull=True, storeAsAttr=cumabuntag, applyFun=sum)
		taxnode.write_attr_to_table(['name', 'taxon', 'rank', abuntag, cumabuntag], nfout, sep='\t')

def usage():
	s =  "python parseKronaGetReadsByTaxid.py [options]\n"
	s += "\transform into a table the data from a Krona HTML file, including abundamces and cumulative abundances\n"
	s += "Mandatory Options:\n"
	s += "\t--krona.html=file.path\n\t\tpath to Krona vizualisation HTML file\n"
	s += "\t--dir.output=directory.path\n\t\tpath to directory where output FASTQ files will be written\n"
	s += "\t--fetch.taxa=taxon1[,taxon2 [, ...]]\n\t\tcomma-separated list of names (enclosed in quotes) or NCBI taxids of taxa for which data must be extracted\n"
	s += "\t--exclude.taxa=taxon1[,taxon2 [, ...]]\n\t\tcomma-separated list of names (enclosed in quotes) or NCBI taxids of taxa for which data must NOT be extracted\n"
	s += "\tA value for at least one option among --fetch.taxa and --exclude.taxa must be specified. If no taxa to be fetched have been specified by --fetch.taxa, all reads but the ones specified by --exclude.taxa will be extracted.\n"
	s += "Facultative Options:\n"
	s += "\t--abundance.field\n\t\tthe name of the field to be taken as abundance in the Lrona HTML file (default to 'abundances')\n"
	s += "\t-k --concatenate\n\t\twhether to pool tables extracted for all querry taxa\n"
	s += "\t-a --append\n\t\twhether ouptut must be appended to existing files (defaults to False)\n"
	s += "\t-p --pickled=file.path\n\t\twether to dump pickled the offsets of reads in FASTQ files for faster running next time or loading the data if they already exist\n"
	s += "\t-v --verbose\n\t\tprint information on processing of reads and taxonomic data (defaults to False)\n"
	s += "\t-h --help\n\t\tprint this help message and quit\n"
	# think of options for rank-specific extraction of reads
	return s
	
	
if __name__=="__main__":

	#~ homedir = os.environ['HOME']
	#~ dirresult = '%s/oral_metagenomes'%(homedir)
	#~ outputdir = '%s/taxonomy-filtered_reads'%(dirresult)
	#~ nfhtml = '%s/run1_AE08.krona.filtered.050.html'%(dirresult)

	try:
		options, args = getopt.getopt(sys.argv[1:], 'kvha:', \
		["krona.html=", "dir.output=", "fetch.taxa=", "exclude.taxa=", "abundance.field=", "concatenate", "verbose", "append", "help"])
	except getopt.GetoptError, err:
		print err
		print usage()
		sys.exit(2)	
	dopt = dict(options)
	#~ print dopt
	
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(1)
	verbose = (('-v' in dopt) or ('--verbose' in dopt))

	nfhtml = dopt['--krona.html']
	outputdir = dopt['--dir.output']
	slftaxa = dopt.get('--fetch.taxa')
	sletaxa = dopt.get('--exclude.taxa')
	abuntag = dopt.get('--abundance.field', 'abundances')
	concat = (('-k' in dopt) or ('--concatenate' in dopt))
	outmode = 'a' if (('-a' in dopt) or ('--append' in dopt)) else 'w'
	
	if slftaxa: lftaxa = [taxon.strip('"\'')for taxon in slftaxa.split(',')]
	elif sletaxa: lftaxa = ["Root"]
	else: print "Must specify a value for at least one option among --fetch.taxa and --exclude.taxa." ; sys.exit(2)
	if sletaxa: letaxa = [taxon.strip('"\'')for taxon in sletaxa.split(',')]
	else: letaxa = []
	
	main(nfhtml, outputdir, lftaxa=lftaxa, letaxa=letaxa, abuntag=abuntag, concat=concat, outmode=outmode, verbose=verbose)
